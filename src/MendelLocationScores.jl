"""
This module orchestrates linkage analysis via location scores.
"""
module MendelLocationScores
#
# Required OpenMendel packages and modules.
#
using MendelBase
# using DataStructures                  # Now in MendelBase.
# using ModelConstruction               # Now in MendelBase.
# using ElstonStewartPreparation        # Now in MendelBase.
# using ElstonStewartEvaluation         # Now in MendelBase.
using Search
using SearchSetup
#
# Required external modules.
#
using DataFrames                        # From package DataFrames.

export LocationScores

"""
This is the wrapper function for the Location Scores analysis option.
"""
function LocationScores(control_file = ""; args...)

  const LOCATION_SCORES_VERSION :: VersionNumber = v"0.1.0"
  #
  # Print the logo. Store the initial directory.
  #
  print(" \n \n")
  println("     Welcome to OpenMendel's")
  println(" Location Scores analysis option")
  println("        version ", LOCATION_SCORES_VERSION)
  print(" \n \n")
  println("Reading the data.\n")
  initial_directory = pwd()
  #
  # The user specifies the analysis to perform via a set of keywords.
  # Start the keywords at their default values.
  #
  keyword = set_keyword_defaults!(Dict{ASCIIString, Any}())
  #
  # Keywords unique to this analysis should be first defined here
  # by setting their default values using the format:
  # keyword["some_keyword_name"] = default_value
  #
  keyword["flanking_distance"] = [0.5, 0.5]
  keyword["flanking_markers"] = 1
  keyword["gender_neutral"] = true
  keyword["lod_score_table"] = "Lod_Score_Frame.txt"
  #
  # Process the run-time user-specified keywords that will control the analysis.
  # This will also initialize the random number generator.
  #
  process_keywords!(keyword, control_file, args)
  #
  # Check that the correct analysis option was specified.
  #
  lc_analysis_option = lowercase(keyword["analysis_option"])
  if (lc_analysis_option != "" &&
      lc_analysis_option != "locationscores")
     throw(ArgumentError(
       "An incorrect analysis option was specified.\n \n"))
  end
  keyword["analysis_option"] = "LocationScores"
  #
  # Read the genetic data from the external files named in the keywords.
  #
  (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
    read_external_data_files(keyword)
  #
  # Execute the specified analysis.
  #
  println(" \nAnalyzing the data.\n")
  execution_error = location_scores_option(pedigree, person, nuclear_family,
    locus, locus_frame, phenotype_frame, pedigree_frame, keyword)
  if execution_error
    println(" \n \nERROR: Mendel terminated prematurely!\n")
  else
    println(" \n \nMendel's analysis is finished.\n")
  end
  #
  # Finish up by closing, and thus flushing, any output files.
  # Return to the initial directory.
  #
  close(keyword["output_unit"])
  cd(initial_directory)
  return nothing
end # function LocationScores

"""
This function maps a trait locus by linkage. In location scores
the trait locus is slid across a marker map. The likelihood of the
data is evaluated in each window with the trait locus positioned
among a fixed number of flanking markers. To the extend possible
the flanking markers occur in balanced numbers to the left and
right of the trait locus.
"""
function location_scores_option(pedigree::Pedigree, person::Person,
  nuclear_family::NuclearFamily, locus::Locus, locus_frame::DataFrame, 
  phenotype_frame::DataFrame, pedigree_frame::DataFrame,
  keyword::Dict{ASCIIString, Any})
  #
  # Problem formulation checks.
  #
  io = keyword["output_unit"]
  if locus.loci <= 1
    println("Error: This option requires at least two loci.")
    println(io, "Error: This option requires at least two loci.")
    return execution_error = true
  end
  trait_place = 0
  for l = 1:locus.loci
    if locus.name[l] == keyword["trait"]
      trait_place = l
      break
    end
  end
  if trait_place == 0
    println("Error: No match to the trait among the possible loci.")
    println(io, "Error:  No match to the trait among the possible loci.")
    return execution_error = true
  end
  flanking_markers = keyword["flanking_markers"]
  if flanking_markers < 1 || flanking_markers > locus.loci
    println("Error: The number of flanking markers is out of range.")
    println(io, "Error: The number of flanking markers is out of range.")
    return execution_error = true
  end
  #
  # Determine the number of parameters and grid points. In grid mode
  # the female and male map distances are maintained in a constant ratio.
  #
  if keyword["travel"] == "grid"
    keyword["parameters"] = 1
    if keyword["points"] == 0
      keyword["points"] = 10
    end
  elseif keyword["travel"] == "search"
    keyword["parameters"] = 2
    if keyword["gender_neutral"]
      keyword["constraints"] = 1
    end
  end
  keyword["goal"] = "maximize"
  #
  # Prepare to eliminate genotypes and lump alleles.
  #
  keyword["eliminate_genotypes"] = true
  keyword["lump_alleles"] = true
  #
  # Prepare a pointer p defining the various marker windows surrounding
  # the trait locus.
  #
  p = collect(1:locus.loci)
  p[1] = trait_place
  for j = 1:trait_place - 1
    p[j+1] = j
  end
  #
  # Initialize the model loci.
  #
  locus.model_loci = flanking_markers + 1
  locus.model_locus = zeros(Int, locus.model_loci)
  locus_list = collect(1:locus.model_loci)
  #
  # Define an empty dataframe to hold the mapping results.
  #
  lodscore_frame = DataFrame(LeftFlank = ASCIIString[],
    Trait = ASCIIString[], RightFlank = ASCIIString[],
    XXCentiMorgans = Float64[], XYCentiMorgans = Float64[],
    LodScore = Float64[])
  #
  # Loop over the map locations of the trait locus, whose position
  # within the model loci is determined by the variable locus.trait.
  #
  locus.trait = 1
  for l = 1:locus.loci
    for j = 1:locus.model_loci
      locus.model_locus[j] = p[locus_list[j]]
    end
    #
    # Prepare to standardize lod scores for the current window by putting
    # the trait locus on the far left of the model loci and at an infinite
    # distance.
    #
    for j = locus.trait:-1:2
      locus.model_locus[j] = locus.model_locus[j - 1]
    end
    locus.model_locus[1] = trait_place
    locus.morgans[:, trait_place] = - Inf
    (saved_locus_trait, locus.trait) = (locus.trait, 1)
    keyword["analysis_option"] = ""
    #
    # Define the parameter data structure.
    #
    parameter = set_parameter_defaults(keyword)
    parameter = initialize_optimization(locus, parameter, keyword)
    #
    # Calculate recombination fractions between adjacent model loci.
    #
    locus.theta = model_recombination_fractions(locus, keyword)
    #
    # Fetch the instructions for conducting the Elston-Stewart algorithm.
    #
    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree, 
      person, nuclear_family, locus, keyword)
    if any(elston_stewart_count .>  keyword["complexity_threshold"])
      println(io, "Marker $locus.name[loc] was skipped because one or more ",
        "pedigrees exceeds the complexity threshold.")
      locus.trait = saved_locus_trait
      (locus.trait, locus_list) = update_location_score_markers(locus_list,
        locus.trait, l, locus.loci, flanking_markers)
      continue
    end
    #
    # Evaluate the standardizing loglikelihood.
    #
    loglikelihood_at_infinity = elston_stewart_loglikelihood(pedigree, person,
      locus, parameter, instruction, keyword)
    #
    # Restore the model loci to their proper order.
    #
    locus.trait = saved_locus_trait
    for j = 2:locus.trait
      locus.model_locus[j - 1] = locus.model_locus[j]
    end
    locus.model_locus[locus.trait] = trait_place
    keyword["analysis_option"] = "LocationScores"
    #
    # Set the optimization title for the current trait window.
    #
    if l == 1
      keyword["title"] =
        "Location scores for " * locus.name[trait_place] *
        " and the closest flanking marker " * locus.name[locus.model_locus[2]]
    elseif l == locus.loci
      keyword["title"] =
        "Location scores for " * locus.name[trait_place] *
        " and the closest flanking marker " *
        locus.name[locus.model_locus[l - 1]]
    else
      keyword["title"] =
        "Location scores for " * locus.name[trait_place] *
        " and the closest flanking markers " *
        locus.name[locus.model_locus[l - 1]] *
        " and " * locus.name[locus.model_locus[l + 1]]
    end
    #
    # Define the parameter data structure.
    #
    parameter = set_parameter_defaults(keyword)
    parameter = initialize_optimization(locus, parameter, keyword)
    #
    # Calculate recombination fractions between adjacent model loci.
    #
    locus.theta = model_recombination_fractions(locus, keyword)
    #
    # Fetch the instructions for conducting the Elston-Stewart algorithm.
    #
    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree, 
      person, nuclear_family, locus, keyword)
    if any(elston_stewart_count .>  keyword["complexity_threshold"])
      println(io, "Marker $locus.name[loc] was skipped because one or more ",
        "pedigrees exceeds the complexity threshold.")
      (locus.trait, locus_list) = update_location_score_markers(locus_list,
        locus.trait, l, locus.loci, flanking_markers)
      continue
    end
    #
    # Pass the variables to optimize for maximum likelihood estimation.
    #
    function fun(par)
      copy!(parameter.par, par)
      f = elston_stewart_loglikelihood(pedigree, person, locus, parameter, 
        instruction, keyword)
      return (f, nothing, nothing)
    end # function fun
    (best_par, best_value) = optimize(fun, parameter)
    #
    # Insert the results in the data frame.
    #
    if locus.trait > 1
      left = locus.name[locus.model_locus[locus.trait - 1]]
    else
      left = ""
    end
    if locus.trait < locus.loci
      right = locus.name[locus.model_locus[locus.trait + 1]]
    else
      right = ""
    end
    if parameter.travel == "search"
      lod = parameter.function_value[end] - loglikelihood_at_infinity
      lod = log10(exp(1.0)) * lod
      push!(lodscore_frame, [left, locus.name[trait_place], right,
        100.*parameter.par[1], 100.0*parameter.par[end], lod])
    else
      for i = 1:parameter.points
        lod = parameter.function_value[i] - loglikelihood_at_infinity
        lod = log10(exp(1.0)) * lod
        push!(lodscore_frame, [left, locus.name[trait_place], right,
          100.0*parameter.grid[i, 1], 100.0*parameter.grid[i, end], lod])
      end
    end
    #
    # Update the marker window.
    #
    (locus.trait, locus_list) = update_location_score_markers(locus_list,
      locus.trait, l, locus.loci, flanking_markers)
  end
  writetable(keyword["lod_score_table"], lodscore_frame)
  show(lodscore_frame)
  return execution_error = false
end # function location_scores_option

"""
This function updates the list of markers involved in the
current location score calculation. The number of flanking
markers equals flanking_markers, the number of loci equals loci,
the position of the trait locus equals i, and the current list
count equals l.
"""
function update_location_score_markers(list::Vector{Int}, i::Int, l::Int,
  loci::Int, flanking_markers::Int)

  if l == loci
    return (i, list)
  elseif l <= div(flanking_markers,2) || l >= loci - div(flanking_markers + 1, 2)
    (list[i], list[i + 1]) = (list[i + 1], list[i])
    i = i + 1
  else
    for j = 1:flanking_markers + 1
      if j != i
        list[j] = list[j] + 1
      end
    end
  end
  return (i, list)
end # function update_location_score_markers

end # module MendelLocationScores

