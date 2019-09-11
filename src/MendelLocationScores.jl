__precompile__()

"""
This module orchestrates linkage analysis via location scores.
"""
module MendelLocationScores
#
# Required OpenMendel packages and modules.
#
using MendelBase
# namely: DataStructures, ModelConstruction,
# ElstonStewartPreparation, ElstonStewartEvaluation
using MendelSearch
#
# Required external modules.
#
using CSV
using DataFrames
using LinearAlgebra

export LocationScores

"""
This is the wrapper function for the Location Scores analysis option.
"""
function LocationScores(control_file = ""; args...)
  #
  # Print the logo. Store the initial directory.
  #
  print(" \n \n")
  println("     Welcome to OpenMendel's")
  println(" Location Scores analysis option")
  print(" \n \n")
  println("Reading the data.\n")
  initial_directory = pwd()
  #
  # The user specifies the analysis to perform via a set of keywords.
  # Start the keywords at their default values.
  #
  keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
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
    locus_frame, phenotype_frame, person_frame, snp_definition_frame) =
    read_external_data_files(keyword)
  #
  # Check if SNP data were read.
  #
  if snpdata.snps != 0
    println(" \n\nERROR: This analysis does not use data from SNP files!\n")
  else
  #
  # Execute the specified analysis.
  #
    println(" \nAnalyzing the data.\n")
    execution_error = location_scores_option(pedigree, person, nuclear_family,
      locus, locus_frame, phenotype_frame, person_frame, keyword)
    if execution_error
      println(" \n \nERROR: Mendel terminated prematurely!\n")
    else
      println(" \n \nMendel's analysis is finished.\n")
    end
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
  phenotype_frame::DataFrame, person_frame::DataFrame,
  keyword::Dict{AbstractString, Any})
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
  lodscore_frame = DataFrame(LeftFlank = AbstractString[],
    Trait = AbstractString[], RightFlank = AbstractString[],
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
    locus.morgans[:, trait_place] .= - Inf
    (saved_locus_trait, locus.trait) = (locus.trait, 1)
    keyword["analysis_option"] = "StandardizingLoglikelihood"
    #
    # Define the parameter data structure.
    #
    parameter = set_parameter_defaults(keyword)
    parameter =
      initialize_optimization_location_score!(locus, parameter, keyword)
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
    loglikelihood_at_infinity =
      elston_stewart_loglikelihood(penetrance_location_score,
      prior_location_score, transmission_location_score,
      pedigree, person, locus, parameter, instruction, person_frame, keyword)
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
    parameter =
      initialize_optimization_location_score!(locus, parameter, keyword)
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
    # Pass the variables to search for maximum likelihood estimation.
    #
    function fun(par)
      copyto!(parameter.par, par)
      f = elston_stewart_loglikelihood(penetrance_location_score,
        prior_location_score, transmission_location_score,
        pedigree, person, locus, parameter, instruction, person_frame, keyword)
      return (f, nothing, nothing)
    end # function fun
    (best_par, best_value) = mendel_search(fun, parameter)
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
        100.0 * parameter.par[1], 100.0 * parameter.par[end], lod])
    else
      for i = 1:parameter.points
        lod = parameter.function_value[i] - loglikelihood_at_infinity
        lod = log10(exp(1.0)) * lod
        push!(lodscore_frame, [left, locus.name[trait_place], right,
          100.0 * parameter.grid[i, 1], 100.0 * parameter.grid[i, end], lod])
      end
    end
    #
    # Update the marker window.
    #
    (locus.trait, locus_list) = update_location_score_markers(locus_list,
      locus.trait, l, locus.loci, flanking_markers)
  end
  lod_table_file = string(keyword["lod_score_table"])
  CSV.write(lod_table_file, lodscore_frame;
    writeheader = true, delim = keyword["output_field_separator"],
    missingstring = keyword["output_missing_value"])
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
  elseif l <= div(flanking_markers,2) ||
         l >= loci - div(flanking_markers + 1, 2)
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

"""
Supply a penetrance for individual i.
"""
function penetrance_location_score(person::Person, locus::Locus,
  multi_genotype::Matrix{Int}, par::Vector{Float64}, person_frame::DataFrame,
  keyword::Dict{AbstractString, Any}, startlocus::Int, endlocus::Int, i::Int)

  pen = 1.0
  for l = startlocus:endlocus
    allele1 = multi_genotype[1, l]
    allele2 = multi_genotype[2, l]
    loc = locus.model_locus[l]
    p = 1.0 # for reduced penetrance let p depend on loc
    pen = p * pen
  end
  return pen
end # function penetrance_location_score

"""
Supply a prior probability for founder i.
"""
function prior_location_score(person::Person, locus::Locus,
  multi_genotype::Matrix{Int}, par::Vector{Float64}, person_frame::DataFrame,
  keyword::Dict{AbstractString, Any}, startlocus::Int, endlocus::Int, i::Int)

  prior_prob = 1.0
  for l = startlocus:endlocus
    loc = locus.model_locus[l]
    allele = multi_genotype[1, l]
    frequency = dot(person.admixture[i, :], locus.frequency[loc][:, allele])
    prior_prob = prior_prob * frequency
    if !locus.xlinked[loc] || !person.male[i]
      allele = multi_genotype[2, l]
      frequency = dot(person.admixture[i, :], locus.frequency[loc][:, allele])
      prior_prob = prior_prob * frequency
    end
  end
  return prior_prob
end # function prior_location_score

"""
Supply the transmission probability that a parent i with a particular
genotype transmits a particular gamete to his or her child j.
"""
function transmission_location_score(person::Person, locus::Locus,
  gamete::Vector{Int}, multi_genotype::Matrix{Int}, par::Vector{Float64},
  person_frame::DataFrame, keyword::Dict{AbstractString, Any},
  startlocus::Int, endlocus::Int, i::Int, j::Int)
  #
  # For male to male inheritance at an x-linked locus,
  # set the transmission probability equal to 1.
  #
  loc = locus.model_locus[startlocus]
  xlinked = locus.xlinked[loc]
  if xlinked && person.male[i] && person.male[j]
    return trans = 1.0
  end
  #
  # Find the map location of the trait locus.
  #
  if keyword["analysis_option"] == "LocationScores"
    trait_place = locus.model_locus[locus.trait]
    locus.morgans[1, trait_place] = par[1]
    locus.morgans[2, trait_place] = par[end]
    if locus.trait > 1
      left = locus.model_locus[locus.trait - 1]
      d = locus.morgans[1, trait_place] - locus.morgans[1, left]
      locus.theta[1, locus.trait - 1] = map_function(d, "Haldane")
      d = locus.morgans[2, trait_place] - locus.morgans[2, left]
      locus.theta[2, locus.trait - 1] = map_function(d, "Haldane")
    end
    if locus.trait < locus.model_loci
      right = locus.model_locus[locus.trait + 1]
      d = locus.morgans[1, right] - locus.morgans[1, trait_place]
      locus.theta[1, locus.trait] = map_function(d, "Haldane")
      d = locus.morgans[2, right] - locus.morgans[2, trait_place]
      locus.theta[2, locus.trait] = map_function(d, "Haldane")
    end
  end
  #
  # Store an indicator of the sex of the parent.
  #
  if person.male[i]
    i_sex = 2
  else
    i_sex = 1
  end
  #
  # Reduce the computations by considering only the heterozygous loci.
  # Use Trow's formula to express the recombination fraction
  # between two heterozygous loci in terms of the recombination
  # fractions between the adjacent loci that separate them.
  # Set the logical variable found to true when the first heterozygous
  # parental locus is found. Phase records the phase of the most
  # recent heterozygous parental locus.
  #
  trans = 1.0
  found = false
  phase = true
  #
  # We use r as 1/2 times the running product in Trow's formula. See equation
  # 7.10 in Mathematical and Statistical Models for Genetic Analysis, 2nd ed.
  #
  r = 0.5
  for l = startlocus:endlocus
    loc = locus.model_locus[l]
    match1 = multi_genotype[1, l] == gamete[l]
    match2 = multi_genotype[2, l] == gamete[l]
    #
    # Check whether either the first or second parental gene at
    # the current locus matches the gamete gene at this locus.
    # If not, then return with 0 for the transmission probability.
    #
    if !match1 && !match2
      return trans = 0.0
    end
    #
    # Check whether the current locus is heterozygous.
    #
    if match1 != match2
      if found
        if phase == match1
          trans = trans * (0.5 + r) # non-recombination
        else
          trans = trans * (0.5 - r) # recombination
        end
      else
        found = true
        if startlocus == 1 || startlocus == endlocus
          trans = 0.5 * trans
##        else
##          trans = 1.0
        end
      end
      phase = match1 # restart Trow's running product
      r = 0.5
    end
    if found && l < endlocus
      r = r * (1.0 - 2.0 * locus.theta[i_sex, l])
    end
  end
##  if !found; trans = 1.0; end
  return trans
end # function transmission_location_score

"""
Initialize the optimization problem.
"""
function initialize_optimization_location_score!(locus::Locus,
  parameter::Parameter, keyword::Dict{AbstractString, Any})
  #
  # Initialize, bound, and name the parameters.
  #
  if parameter.parameters == 1
    parameter.name[1] = "theta"
  else
    parameter.name[1] = "xxtheta"
    parameter.name[2] = "xytheta"
  end
  small = 1e-4
  if locus.trait > 1
    left = locus.model_locus[locus.trait - 1]
    parameter.min[1] = locus.morgans[1, left] + small
    parameter.min[end] = locus.morgans[end, left] + small
  end
  if locus.trait < locus.model_loci
    right = locus.model_locus[locus.trait + 1]
    parameter.max[1] = locus.morgans[1, right] - small
    parameter.max[end] = locus.morgans[end, right] - small
  end
  d = keyword["flanking_distance"]
  if locus.trait == 1
    right = locus.model_locus[2]
    parameter.min[1] = locus.morgans[1, right] - d
    parameter.min[end] = locus.morgans[end, right] - d
  elseif locus.trait == locus.model_loci
    left = locus.model_locus[locus.trait - 1]
    parameter.max[1] = locus.morgans[1, left] + d
    parameter.max[end] = locus.morgans[end, left] + d
  end
  if keyword["travel"] == "search"
    parameter.par = 0.5 * (parameter.min + parameter.max)
    if keyword["gender_neutral"]
      a = (parameter.max[2] - parameter.min[2])
      b = (parameter.max[1] - parameter.min[1])
      c = a / b
      d = parameter.min[2] - c * parameter.min[1]
      parameter.constraint[1, 1] = c
      parameter.constraint[1, 2] = -1.0
      parameter.constraint_level[1] = c * parameter.min[1] - parameter.min[2]
    end
  else
    if parameter.points == 1
      parameter.grid[1, :] = 0.5 * (parameter.min + parameter.max)
    else
      for j = 1:parameter.points
        a = (j - 1.0) / (parameter.points - 1.0)
        parameter.grid[j, :] = a * parameter.max + (1.0 - a) * parameter.min
      end
      parameter.min[1:end] .= -Inf
      parameter.max[1:end] .= Inf
    end
  end
  return parameter
end # function initialize_optimization_location_score!
#
# Method to obtain path to this package's data files
# so they can be used in the documentation and testing routines.
# For example, datadir("Control file.txt") will return
# "/path/to/package/data/Control file.txt"
#
datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module MendelLocationScores
