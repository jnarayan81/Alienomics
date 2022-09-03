![Image](alien.jpg) 


# Introduction

Welcome to Alienomics user manual !

Alienomics is a pipeline dedicated to detect both candidate HGTs (cHGTs) and contaminants within an assembled genome, without recquiring additional genome assemblies from closely-related species nor computationaly-heavy phylogenetic inferences. Alienomics combines several sources of information such as read coverage, GC content, expression level, scaffold-scale synteny, sequence similarity (against public databases) and interprete them as proxy for gene origin and its integration within the genome. This allows Alienomics to classify each gene as a "self" gene, a contaminant or a candidate HGT.

**Main authors and contact**: 
Jitendra Narayan <jnarayan81@gmail.com>
Paul Simion <paul.simion@univ-rennes1.fr>
Karine Van Doninck <Karine.Van.Doninck@ulb.be>.

**Citation**
If you used Alienomics in your work, please cite the following article in which a preliminary version of Alienomics was used:
Simion, P., Narayan, J., Houtain, A., Derzelle, A., Baudry, L., Nicolas, E., Arora, R., Cariou, M., Cruaud, C., Gaudray, F.R., Gilbert, C., Guiglielmoni, N., Hespeels, B., Kozlowski, D.K.L., Labadie, K., Limasset, A., Llirós, M., Marbouty, M., Terwagne, M., Virgo, J., Cordaux, R., Danchin, E.G.J., Hallet, B., Koszul, R., Lenormand, T., Flot, J.-F., Van Doninck, K., 2021. *Chromosome-level genome assembly reveals homologous chromosomes and recombination in asexual rotifer* Adineta vaga. Science Advances 7, eabg4216. https://doi.org/10.1126/sciadv.abg4216

We are working actively on the publication of the up-to-date Alienomics software tool.


# Table of content

1. [Quick step-by-step guide](#first-run)
2. [Installing Alienomics](#Installing-Alienomics)
 * [Installing using Conda](#Conda-Installation)
 * [Installing Manually](#Manual-Installation)
3. [Alienomics Usage](#Alienomics-Usage)
4. [Test Run](#Test-run)
5. [Detailed Internal Procedures](#Internal-procedures)
6. [Detailed Options](#Detailed-Options)
7. [Configuration File Parameters](#Configuration-Parameters)
 * [Input Files](#Input-files)
 * [Input Settings](#Input-settings)



##############################
<a name="first-run"></a>
# Quick step-by-step guide: "my first Alienomics run" 

### Step 1: installing Alienomics and its database ###
```bash
conda create -n alienomics
conda activate alienomics
conda install -c jnarayan81 alienomics

alienomics_v.0.3 -dd -ddf BURT -ddl /path/to/database_folder/	# -ddl indicates the location where databases has to be installed.
																# this step might take up to several hours (notably due to the preparation of the uniref50 database)
```
Store the information outputed at the end of this step in order to correctly set the paths in the configuration file in **step 3**.

### Step 2: creating and opening configuration file ###
```bash
alienomics_v.0.3 -cc -ccl /path/to/config_file_folder/			# -cc stands for Create Config file and -ccl stands for Create Config Location
nano /path/to/config_file_folder/config_file 					# use here your prefered text editor
```

### Step 3: preparing configuration file for your alienomics run ###
Several mandatory input/output parameters have to be tailored to your analysis. For more details, see the [Configuration File Parameters](#Configuration-Parameters) section.
```
input_genome_fasta = /path/to/input/genome.fasta 				# location of the genome assembly file (fasta format) of the genome of interest
input_genome_gff = /path/to/input/genome/gene_prediction.gff3 	# location of gff file corresponding to the prediction of genes on the genome of interest
input_mapping_bam = /path/to/input/reads_mapping.bam 			# location of bam file corresponding to the mapping of genomic reads on the genome of interest
RNAseq_location = /path/to/RNA-Seq_folder						# location of the folder containing RNA-Seq reads (fastq format)

output_folder = /path/to/Alienomics_output_folder				# location and name of the folder containing the ouptuts for thei Alienomics run
```
A series of other parameters are **very important to correctly set up**:
- **Self taxonomic group name**. Blast hits within that clade will be considered as "self" origin, while blast hits outside that clade will be considered as "alien". Clade name need to exists in NCBI taxonomy.
```
clade_self = taxonomic_group_name		# example: "metazoa"
```
- **Exclusion of taxa**. All taxonomic ID listed here (seperated by "|" characters) will be discarded from the uniref50 database. This allows to discard closely-related species which might share HGT events with the species under focus. Ancient and shared HGT events might not be detected by Alienomics if closely-related species are still present in unrif50 database.
```
excludeTaxID = taxID1|taxID2|taxID3		# example: "104782|10195|96448"
```
- **Expected GC range** of the organism of interest. GC content outside of this range will be regarded as indication for "alien" origin.
```
expected_GC_range = X:Y           	    # Example 26:38 (two values separated with character ':')
```
- **Mean expected read coverage** of the organism of interest. This value will be used within Alienomics to estimate a coverage range. A gene coverage outside of this range will be regarded as indication for "alien" origin.
```
expected_coverage = N                	# N = expected coverage of genes. Recquire reads using the 'input_mapping_bam' option (use "0" to set this option off)
```
- **Number of threads**. Using more threads will increase Alienomics run speed.
```
max_processors = N                     # N = number of processors to use for parallelized steps
```
- **Paths to databases** need to be set up according to the output of the **step 1**


### Step 4: running Alienomics ###
```bash
alienomics_v.0.3 -c config.file | tee log_file.out
```




##############################
<a name="Installing-Alienomics"></a>
# Installing Alienomics 

<a name="Conda-Installation"></a>
### Conda Installation (highly recommended) 

It is recommended to create an environment to avoid any conflits. 
```bash
conda create -n alienomics
conda activate alienomics
```
Now install alienomics from conda server (and then deactivate conda)
```bash
conda install -c jnarayan81 alienomics
conda deactivate alienomics
```
Conda installation automatically manage all perl module dependencies. 
After activating conda environment, you'll be able to run Alienomics
```bash
conda activate alienomics
```

<a name="Manual-Installation"></a>
### Manual Installation (if using conda is not a possibility)

0. Install perl module dependencies using CPAN.
1. Download or clone the <code>alienomics</code> tool. Consider unzipping in a directory to <code>alienomics/</code>.
2. <code>cd</code> into <code>alienomics/bin/</code> directory.
3. Type <code>perl alienomics_v0.3.pl</code>.

 


##############################
<a name="Test-run"></a>
# Alienomics Test Run

In order to check whether Alienomics has been properly installed, a small set of inputs can be found in the `test_run` folder. At this stage, you should have installed Alienomics and already installed the databases (with the following command : `alienomics_v.0.3 -dd -ddf BURT -ddl /path/to/database_folder/`). At the end of the database install, this option notably outputs the exact database locations to use as parameters in the configuration file.
You now need to modify these paths in the `config/test_run.conf` accordingly and to specify a complete path for the output folder (see the `output_folder` parameter in the configuration file). You can also modify the number of threads used by Alienomics (see the `max_processors` parameter). For more details, see the [Configuration File Parameters](#Configuration-Parameters) section.
```bash
cd config
nano test_run.conf 		# use you preferred text editor here
```
Only after having modified the configuration file, you can procede with Alienomics run test:
```bash
alienomics_v0.3 -c test_run.conf
```
Alienomics should run without encountering issues. When this is the case, you can now copy-paste the `test_run.conf` somewhere else and modify it to suit the analysis of your genome of interest!



##############################
<a name="Alienomics-Usage"></a>
# Alienomics Usage 

Alienomics software tool has several options. The main function is triggered by the `-conf` option (short version: `-c`) and will start a complete Alienomics run corresponding to a specified configuration file:
```bash
alienomics_v0.3 -c config/config.conf 
```

Other functions are available and are mostly needed for the installation of Alienomics and the preparation of your analysis:
- <code>[-dd | --database_download]:</code>		Download and install the third party databases.
- <code>[-cc | --create_config]:</code>			Create a default config file to be modified by the user.
- <code>[-v | --version]:</code>				You can check of the version of the the alienomics program.    
- <code>[-w]:</code>							Know about the developers. 





##############################
<a name="Internal-procedures"></a>
# Detailed internal procedures 

### "Blast score" computing

This score is based on blast hits against uniref50 database. The database actually used within Alienomics might be slightly different than the official uniref50 database according to the potential use of the `excludeTaxID` options ine the configuration file. Diamond is used to compute the blast hits, and up to 50 hits are retained (`-k 50`), with only 1 hit per subject (`--max-hsps 1`). Blast hits are then filtered according to the following criteria: evalue <= 1e-05, length >= 60, 30 <= %identity <= 90. This allows the following step of the procedure to focus on "good" hits, filtering out poor hits but also filtering out hits that are "too good" and moight stem from contaminants in public databases. Taxonomic information is then retrieved for every remaining hits, seperating them into two lists: the hits belonging to the clade defined as "self" (`clade_self` in configuration file) and the hits falling outside this clade, thus interpreted as "alien" hits. The majority list is the one with more hits than the other list ("self" versus "alien"). The best hit, according to bitscore, from the majority list is retained, thus defining whether the "blast score" will be positive (best hit belongs to "self" clade) or negative (best hit is "alien"). The absolute value of the blast score is positively correlated with the bitscore value of the retained best hit, and is computed as follows:

[formula for blast score here]



### Interpreting HGT.gff outputs
0.2 vs 0.5 vs 0.8







##############################
<a name="Detailed-Options"></a>
# Detailed Options 

---

**--conf** (-c)

This option indicates the location of the configuration file for a given Alienomics run. This is the main option, leading to a complete run of Alienomics that will analysis all genes and categorize them as "self", "cHGT" or "contaminant". Minimal example command to run Alienomics (using either `--conf` or `-c`):
```bash
alienomics_v.0.3 -c path/to/configuration_file
```

---

**--database_download** (-dd)

This option is used to download and install databases that will be necessary for future Alienomics runs. this option has to be used in combination with 2 other options: **-ddf** (for Database Download Flags) and **-ddl** (for Database Download Location). Each `-ddf` flags trigger the installation of a given database: **B**usco, **U**niref50, **r**RNA, NCBI **T**axonomy. This allows to re-install only one (or several) of the database when the user wants to update them. The following command download and install all databases needed by Alienomics:
```bash
alienomics_v.0.3 -dd -ddf BURT -ddl /path/to/database_folder/
```

---




##############################
<a name="Configuration-Parameters"></a>
# Configuration file parameters

Alienomics inputs and parameters are set using a configuration file. Below is the detailed description of the variable that you will have to set up to suit the particular requirements of your analysis.

<a name="Input-files"></a>
### Input files

The input files correspond to files that are typically avalaible for all genome assemblies. Their location is indicated to alienomics using the 5 following parameters:
```
input_genome_fasta 
input_genome_gff 
input_mapping_bam 
RNAseq_location 
reference_proteome_file 
```
**input_genome_fasta**: The absolute path to the genome assembly of interest. This file should be in fasta format, and sequence names should not include any `_` characters (this is because Alienomics internal procedures will add such a character when computing gene expression level).

**input_genome_gff**: The absolute path to the gene annotation file corresponding to the genome assembly of interest. This file should be in gff format (or gff3).

**input_mapping_bam**: The absolute path to the bam file corresponding to the mapping of reads onto the genome assembly of interest. This file should correspond to the mapping of genomic reads (i.e. it was not designed with RNA-seq reads in mind), and should be in the bam format. This bam file should be indexed (e.g. with a command line such as `samtools index $input_bam`) and this index should be placed in the same folder as the `input_mapping_bam`, with the `.bai` file extension.

**RNAseq_location**: The absolute path to the folder containing the RNA-Seq reads. These reads should be in fastq format and should use `.fastq.gz` as file extensions. Several files can be placed within the `RNAseq_location` folder,  corresponding either to paired-end reads (possibly multiple pairs of paired-end reads) OR to single-end reads (possibly multiple single-end reads): Alienomics will recognize them and handle them automatically. However, If single-end and paired-end reads are mixed within the `RNAseq_location` folder, Alienomics will only use paired-end data, and ignore the single-end data.

**reference_proteome_file**: The absolute path to a trusted reference proteome. This file should be in fasta format and contain only proteic sequences. This parameter is optionnal, an empty parameter value will set this option off.
*This proteome will be used as a reference for "self" gene*. This means that if a gene from your genome of interest seem homologous to a gene from this reference proteome, this gene will more likely be considered as a "self" gene. Consequently, this gene will be less likely to considered as a candidate HGT or a contaminant. 
*This proteome needs to be of good quality*. Ideally, it should be the proteome of a complete and well-annotated genome assembly, from an organism that belongs to the same clade as your organism of interest. For example, if you want to use alienomics on a worm (annelida, metazoa) you might want to use a good-quality proteome from a fly such as a *Drosophila* species (arthropoda, metazoa).
*This proteome should not share HGT events with your genome of interest*. If the reference proteome itself contains many HGTs or if it is too closely-related to your species of interest, it might share ancient HGT events with your species of interest. In such a case, these shared HGT events might not be detected by Alienomics as candidate HGTs. 


<a name="Input-settings"></a>
### Input settings

The input settings are necessary to tailor the alienomics run according to the characteristics of the genome assembly of interest and the evolutionary scale of the HGT events you would like to detect.
```
clade_self 
excludeTaxID 
expected_GC_range 
expected_coverage 
max_processors 
```
**clade_self**: This crucial parameter is the name of a clade that will define the concept of "**self**" during the analysis. This name as to exist in the NCBI taxonomic database (i.e. `clade_self = metazoa`). This will set the boundary between "self" and "alien". Blast hits within that clade will be regarded as evidence for "self" origin, and blast hits outside of that clade will be regarded for evidence for an "alien" origin. 
*Important note 1*: that HGT events that would have occured between two lineages both included within the `clade_self` can **not** be detected by Alienomics. For example, if `clade_self = metazoa`, then genes potentially originating from HGT events between a worm and a fly will be likely categorized as "self".
*Important note 2*: Both the species of interest and the reference proteome must belong to the clade defined as `clade_self` (see also the `reference_proteome_file` parameter in the [Input Files](#Input-files) section)        

**excludeTaxID**: A list of taxonomic ID of species that need to be ignored during blast steps. This list should contain taxon ID seperatd by `|` characters (e.g. `excludeTaxID = 104782|10195|96448`). It is important to list here all the species that might share HGT events with your genome assembly of interest (i.e. species that are too closely-related to the species of interest and species with genomes containing many HGTs or many contaminants).
*Important note*: Currently, only taxon ID at the level of secies are handled by Alienomics (the taxon ID coresponding to "insecta" does not discard all insects from the database).

**expected_GC_range**: Parameter values setting the boundaries of the expected GC content of the genome assembly of interest. The value for this parameter should correspond to 2 numbers seperated by a `:` character (e.g. `expected_GC_range = 38:45`). A gene GC content value falling *outside* this range will be interpreted as an indication for an "**alien**" origin. If you are unsure about the parameter value to use, we suggest you use the `--exploreGCandCOV` option of alienomics (see the [Alienomics Usage](#Alienomics-Usage) section).

**expected_coverage**: Parameter value setting the expected coverage along the genome assembly of interest. This is the expected coverage (i.e. `expected_coverage = 60` for an expected genome coverage of 60x) corresponding to the input bam file set with the `input_mapping_bam` variable (see the [Input Files](#Input-files) section). A gene coverage value that would deviate too much from this value will be interpreted as an indication for an "**alien**" origin. If you are unsure about the parameter value to use, we suggest you use the `--exploreGCandCOV` option of alienomics (see the [Alienomics Usage](#Alienomics-Usage) section).

**max_processors**: number of threads to be used by Alienomics (i.e. `max_processors = 20`). Although we tried to make Alienomics as fast as possible, blasting steps (especially against uniref50 database) can still take some time, and the more numerous the threads, the better.


### Paths to databases

The parameter values correspond to the path to the databases that you have to install before being able to use Alienomics. Their installation is automated using the `--database_download` option (see the [Alienomics Usage](#Alienomics-Usage) section)
```
uniref50_path 
rRNA_file 
busco_file 
```

The values to set for these 3 parameters is outputed at the end of the Alienomics run when the `--database_download` (or `-dd`) option is used. You can copy-paste these values directly into the configuration file.


### output folder location

```
output_folder 
```
This parameter indicates the location where all Alienomics result files will be placed. User should modify this value between every run of Alienomics, as Alienomics will not overwrite an existing output folder. example: `output_folder = /path/to/folder/`.





























##############################
# Older version of user manual

### config.aom
All the mandatory and optional parameters are needed to be set in config.aom file. 

The primary alienomics configuration file is config.aom. It includes a large number of configuration statements that do not need to be changed for a basic run. In fact, only a few changes to this file are required to get a basic alienomics running. Because the file is quite large, rather than clogging up this article with unnecessary information, I will only show the directives that need to be changed.

First, spend some time reading through the config.aom file to become acquainted with it. One of my favorite aspects of configuration files is the number of comments that describe the various sections and configuration directives. The config.aom file is no exception, as it is well-documented. Use these comments to figure out what the file is doing.

#### Create sample config.aom
To create a sample config file, try -cc or –create_config flags in alienomics:
```sh
alienomics_v.0.3 -cc -ccl /home/download/
```
Here -ccl stands for create config location

#### External programs 
External Program action tasks are used to launch and execute other programs, which can be beneficial when processing your alienomics job file. The external program is run with the parameters specified in the Alienomics configuration file.

In the configuration file, specify the absolute path to all mandatory tools as follows:
Name=Absolute Path

```sh
diamond = /usr/local/bin/diamond
blastn = /usr/bin/blastn
blastx = /usr/bin/blastx
makeblastdb= /usr/bin/makeblastdb
augustus = /usr/local/bin/augustus
bwakit = //usr/bin/bwa
seqtk = /home/urbe/anaconda3/envs/seqtk_paul/bin/seqtk
samtools = /usr/local/bin/samtools
kallisto = /home/urbe/Tools/kallisto/kallisto
```
Because "Name" is a reserved keyword, you can only change the "Absolute path". You can find out the applications path by following command in your terminal

```sh
#The "which" command has had problems getting the proper path (confusion between environment and dot files). For "type", you can get just the path with the -p argument.

type <tool_name> # blast / diamond
which <tool_name>

#The POSIX standard way to do this is command -v git. All UNIX-like systems should support this.
command -v <tool_name>
```
If the tools are missing, you can install them manually with conda or any methods you prefer.

#### Mandatory Inputs
To run the program smoothly, the following files must be provided.
- Genome file
#The genome file should be in multi fasta format, with no special characters in the fasta header.
```sh
Input_genome = genome.fasta 	# Absolute path of genome fasta file
```
- GFF file
GFF file should follow GFF3 format https://learn.gencore.bio.nyu.edu/ngs-file-formats/gff3-format/
```sh
local_gff = genome.gff3		# Absolute path of gene annotation gff file
```
#### Optional Inputs
> Note: `Warning:` all HGTs existing in the reference genome will be considered as 'self' in the genome under analysis

```sh
bam = genome.NDPD.ONT.bam.minimap2.sorted	# sorted and indexed BAM file (for coverage criteria)
kallisto_reads=kallisto_on_AV20		# kallisto reads location input
reference_genome_file = test2.REF.genome.fasta		# fasta file containing reference genomic sequences (reference for self)
reference_proteome_file = test7.REF.protein.fasta		# fasta file containing all reference proteic sequences (reference for self)
```

#### Database Settings
The establishment of an alienomics database has a direct impact on the quality and viability of HGT research, making it a critical component of this tool. Creating a high-quality database from online resources necessitates adhering to a strict data-management procedure.

> Note: `Warning:` do not add a "/" character at the end of the paths

To download well established databases you need to activate -dd or –download_database flags in alienomics script:
```sh
alienomics_v.0.3 -dd -ddf BURT -ddl /home/download/
```
Here download database flags(ddf) are read as follows: **B->Busco, U->Uniref50, R->rRNA, and T->Taxonomy** Databases. Remember these are case sensitive.

```sh
alien_dir = /home/urbe/Tools/Alienomics_v0.2/bin/ #Folder location of alienomics

rRNA_file =  SILVA_128_LSUParc_tax_silva.fasta #rRNA database location
busco_file = BUSCOData/metazoa_odb9/ancestral #BUSCO database location
tnfDB_file = Alienomics/TNFDB/tnf_1975_genomes.longuestseq.tab #In house TNF database

diamonddb_path = Uniref50DB					# or RefSeqDB ? or Uniref50DB (37G) ? or UniprotDB (604M) ?
Uniref50_taxlist = uniref50.taxlist		# location of home-made entryID + taxlist table (see user-manual for setting this up)

#Useful when megablast is used
blastdb_path = blastDB 					# useless now ?
taxdump_path = taxdump					# useless now ?
accession2taxid_path = accession2taxid			# useless now ?
AUGUSTUS_CONFIG_PATH = augustus.2.5.5/config		# useless now ?

localdb = yes					# Check into the local database created by "create_local_db" above -- blast database of you reference genome
create_local_db = yes				# Create a local blast db  ... Let them keep always "yes" for NOW
# Make localDB makeblastdb -in ../WadinetaGapFilled.gapfilled.iteration1.fa -parse_seqids -dbtype nucl -out myDB
# Make localDB makeblastdb -in ../WadinetaGapFilled.gapfilled.iteration1.fa -parse_seqids -dbtype nucl -out myDB -name -myDB

```
#### Additional Settings
For different Alienomics runs, different configurations and settings are required. This section allows you to switch between different sets of settings for each run. You can change the settings to suit your needs. Please read it carefully before selecting the default alienomics profile. This saves your time and confusion. Don't forget to look into the phylogenetic relationships between the species of interest.

You can customize alienomics behavior by adjusting the additional settings. For instance, if you are certain about the feeder species, you may want to adjust to a specific genome and let alienomics think/predict wisely.

```sh
blast_task = blastn 				# when using blast+ on a nucleotidic database (User can choose among blastn/blastn-short/megablast/dc-megablast)
diamond_task = blastx				# when using diamond on a proteic database
mode = gene	                        	# Main analysis mode. Currently ALIEN supports only gene and random models analysis.
verbose = 1      				# Print all detail #1 to print nice log messages telling you what is going on. 0 otherwise
resume = 0					# Resume the last run
ActivateTNF = no				# 'no' by default; if set to 'yes', it will run TNF criteria
ActivateBIGBLAST = on
rejectGi=/home/urbe/Tools/Alienomics_v0.1/Alienomics/config/rejectList.gi 	# GI or Accession list for blast to ignore - neglectGi and rejectGi work in combination
negateGi=0									# ignore the gi list in file: 1 for yes ; 0 for no
taxlevel_species = adineta_vaga			# taxlevel_species you believe their presence make it strong candidate !!!! for example adineta_ricciae
feedlevel_species = escherichia_coli		# Species you feeded with .... currently only one supported !!!!
predict_gene = 0				# 1 to yes / O to NO [default should be "0" !]
blastopt = offline			 	# Blast running mode: Possible values are "online" or "offline" !!!!
pvalue = 0.05                            	# p-values for positive selection detection
qvalue = 0.05                            	# q-values for positive selection detection

```

#### Plot settings
Alienomics automatically assigns different colors, line styles, and markers to plot objects when you plot multiple data sets on the same axes. Colors, line styles, and resolutions can all be changed when using plotting functions.

```sh
resolution = 1000				# Resolution for plotting the gene graph
plot_contig = yes				# "yes" to print the contig plot with score else "no" --- It will increase the computtion time
min_gene_number = 3				# Minumum gene to consider for plotting graph
```

#### Output folder settings
The user must specify the name of the out folder, which will hold all intermediate files. Change the default output path and rename all log files to something else.
```sh
project_dir_path = test1.output.folder	# output folder name
summary = summary.out			# summary file name 
result = scaffoldscore.out		# scaffold info file name
score = genescore.out			# gene info file name
log_file = final_log.log		# Log file name
```



### Update
Version 0.3.1 has been released !
- One of the most significant changes is the third party database download module.
- The sample config.aom is now created automatically with a flag

### Motivation and Reasoning
The Kelly Family, a European-American pop group, wrote the tune "Fell in Love with an Alien." This refers to my fascination with genes that are horizontally transferred between unicellular and/or multicellular organisms rather than vertically transmitted from parent to offspring (reproduction). Many organisms' evolution is influenced by horizontal gene transfer ([HGT](https://scholar.google.co.in/scholar?hl=en&as_sdt=0%2C5&q=hgt&btnG=)). This work was inspired by my love-hate relationship with bdelloid rotifer *Adineta vaga*'s relatively unknown DNA fragments. A surprisingly high proportion of bdelloid genes are homologous to nonmetazoan orthologs, mainly from bacteria but also from fungi and plants, suggesting that the rate of HGT into bdelloid genomes is at least an order of magnitude greater than in other eukaryotes. Interestingly, many genes derived from nonmetazoans via HGT are expressed and known to be functional.

Unfortunately, most published methods so far depend on just a few (if not just one) types of metrics, rendering them susceptible to errors. Plain aligners have limitations and can produce unwanted artifacts. It is important to note that HGT is an evolutionary hypothesis, and that a sequence alignment alone is insufficient to confirm such a hypothesis. For HGT identification and contamination filtration, phylogenetic and parametric methods are commonly used. The phylogenetic approach is based on incongruent relationships among taxa, whereas the parametric method is based on changes in nucleotide composition and homology search. In this study, we used coverage, GC content, sequence similarity, gene expression, and synteny information to categorize genome data as clean, contaminated, or HGT DNA sequences.  



## Frequently Asked Questions

### How to download third party database in alienomics?

You need to activate -dd (or --download_database) flags in alienomics:

```
alienomics_v.0.3 -dd -ddf BURT -ddl /home/download/
```
Here, the download database flags (ddf) are read as follows: B->Busco, U->Uniref50, R->rRNA, and T->Taxonomy Databases. These are case sensitive.

### How to create config files in alienomics

To create a sample config file, try -cc or --create_config flags in alienomics:

```
alienomics_v.0.3 -cc -ccl /home/download/
```
Here -ccl stands for create config location  

## Disclaimer

This software is freely distributed 'as-is.' The programmer is not liable for any damage that may occur. You can share the program's original files with your friends, colleagues, and so on. Report any new releases or bugs on https://github.com/jnarayan81/Alienomics/issues.

## Contribution

Feel free to clone this repository and use it under the licensing terms.

Additionally, as the project is on github, you may submit patches, ticket requests, edit the wiki, send pull requests - anything you like and have the permissions to do. I will encourage any suggestions from followers :)

As always, you can contact the following authors: Jitendra Narayan <jnarayan81@gmail.com>, Paul Simion <paul.simion@univ-rennes1.fr>, Karine Van Doninck <Karine.Van.Doninck@ulb.be>.



