### Overview
Mendel Location Scores is a component of the umbrella [OpenMendel](https://openmendel.github.io) project. This package performs multipoint linkage analysis. That is, the analysis considers multiple linked markers in combination when mapping the trait.

### Appropriate Problems and Data Sets
Mendel Location Scores analysis can use sib pairs, nuclear families, or large pedigrees, along with marker and trait data, to localize genes.

### Installation
*Note: The three OpenMendel packages (1) [SnpArrays](https://openmendel.github.io/SnpArrays.jl/latest/), (2) [Search](https://openmendel.github.io/Search.jl), and (3) [MendelBase](https://openmendel.github.io/MendelBase.jl) must be installed before any other OpenMendel package will run. It is easiest if these three packages are installed in the above order and before any other OpenMendel package.*

Within Julia, use the package manager to install Mendel:

    Pkg.clone("https://github.com/OpenMendel/MendelLocationScores.jl.git")

This package supports Julia v0.4 and v0.5.

### Input Files
The MendelLocationScores analysis package uses the following input files. Example input files can be found in the [docs]( https://github.com/OpenMendel/MendelLocationScores.jl/tree/master/docs) subfolder of the MendelLocationScores project. (An analysis won't always need every file type below.)

* [Control File](#control-file): Specifies the names of your data input and output files and any optional parameters (*keywords*) for the analysis. (For a list of common keywords, see [Keywords Table](https://openmendel.github.io/MendelBase.jl/#keywords-table)).
* [Locus File]( https://openmendel.github.io/MendelBase.jl/#locus-file): Names and describes the genetic loci in your data.
* [Pedigree File]( https://openmendel.github.io/MendelBase.jl/#pedigree-file): Gives information about your individuals, such as name, sex, family structure, and ancestry.
* [Phenotype File]( https://openmendel.github.io/MendelBase.jl/#phenotype-file): Lists the available phenotypes.
* [SNP Definition File](https://openmendel.github.io/MendelBase.jl/#snp-definition-file): Defines your SNPs with information such as SNP name, chromosome, position, allele names, allele frequencies.
* [SNP Data File](https://openmendel.github.io/MendelBase.jl/#snp-data-file): Holds the genotypes for your data set. Must be a standard binary PLINK BED file in SNP major format. If you have a SNP data file you must have a SNP definition file.

### Control file<a id="control-file"></a>
The Control file is a text file consisting of keywords and their assigned values. The format of the Control file is:

	Keyword = Keyword_Value(s)

Below is an example of a simple Control file to run Location Scores:

	#
	# Input and Output files.
	#
	locus_file = location scores LocusFrame.txt
	pedigree_file = location scores PedigreeFrame.txt
	phenotype_file = location scores PhenotypeFrame.txt
	output_file = location scores Output.txt
	lod-score-table = location scores Table Output.txt
	#
	# Analysis parameters for Location Scores option.
	#
	disease_status = RADIN
	travel = grid
	gender-neutral = false

In the example above, there are eight keywords. Three keywords specify the input files: *location scores LocusFrame.txt*, *location scores PedigreeFrame.txt*, and *location scores PhenotypeFrame.txt*. The next two keywords specify the output files: *location scores Output.txt* the results file - and *location scores Table Output.txt* - a table of the location scores. The last three keywords specify analysis parameters. The text after the '=' are the keyword values.

### Keywords<a id="keywords-table"></a>
This is a list of OpenMendel keywords specific to Location Scores. A list of OpenMendel keywords common to most analysis package can be found [here](https://openmendel.github.io/MendelBase.jl/#keywords-table). The names of keywords are *not* case sensitive. (The keyword values *may* be case sensitive.)

Keyword          |   Default Value    | Allowed Values |  Short Description       
----------------      |  ----------------       |  ----------------      |  ----------------
flanking_distance  |  [0.5, 0.5]  |  Real  |   Marker flanking distance for location scores
  flanking_markers  |  1  |  Integer  |  Number of flanking markers in imputation window
  gender_neutral  |  TRUE  |  TRUE, FALSE  |  Whether different male/female recombination  fractions are considered
lod_score_table  |  "Lod_Score_Frame.txt"  | User defined output text file name  |   Creates a lod score table output file

### Data Files
Location Scores requires a [Control file](https://openmendel.github.io/MendelBase.jl/#control-file), and a [Pedigree file](https://openmendel.github.io/MendelBase.jl/#pedigree-file). Genotype data is provided in a [SNP data file]( https://openmendel.github.io/MendelBase.jl/#snp-data-file), with a [SNP Definition File]( https://openmendel.github.io/MendelBase.jl/#snp-definition-file) describing the SNPs. OpenMendel will also accept [PLINK format](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml) FAM and BIM files. Details on the format and contents of the Control and data files can be found on the [MendelBase](https://openmendel.github.io/MendelBase.jl) documentation page. There are example data files in the Location Scores [docs]( https://github.com/OpenMendel/MendelLocationScores.jl/tree/master/docs) folder.

### Running the Analysis

To run this analysis package, first launch Julia. Then load the package with the command:

     julia> using MendelLocationScores

Next, if necessary, change to the directory containing your files, for example,

     julia> cd("~/path/to/data/files/")

Finally, to run the analysis using the parameters in the control file Control_file.txt use the command:

     julia> LocationScores("Control_file.txt")

*Note: The package is called* MendelLocationScores *but the analysis function is called simply* LocationScores.

<!--- ### Interpreting the results
... --->

### Citation

If you use this analysis package in your research, please cite the following reference in the resulting publications:

*Lange K, Papp JC, Sinsheimer JS, Sripracha R, Zhou H, Sobel EM (2013) Mendel: The Swiss army knife of genetic analysis programs. Bioinformatics 29:1568-1570.*

<!--- ### Contributing
We welcome contributions to this Open Source project. To contribute, follow this procedure ... --->

### Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.