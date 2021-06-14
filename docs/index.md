![Image](alien.jpg) 

## Welcome to Alienomics Homepage !
Thank you for visiting.  
Hello and welcome to the Alienomics's github page. This page contains the following information:

- Motivation and reasoning for creating this tool, step-by-step installation and usage examples
- If you have any questions, remarks, or would like to share your stories, please email me at jnarayan81[at]gmail[dot]com !

### Update
Version 0.3.1 has been released !
- One of the most significant changes is the third party database download module.
- The sample config.aom is now created automatically with a flag

### Motivation and Reasoning
The Kelly Family, a European-American pop group, wrote the tune "Fell in Love with an Alien." This refers to my fascination with genes that are horizontally transferred between unicellular and/or multicellular organisms rather than vertically transmitted from parent to offspring (reproduction). Many organisms' evolution is influenced by horizontal gene transfer ([HGT](https://scholar.google.co.in/scholar?hl=en&as_sdt=0%2C5&q=hgt&btnG=)). This work was inspired by my love-hate relationship with Adineta vaga's relatively unknown DNA fragments. A surprisingly high proportion of bdelloid genes are homologous to nonmetazoan orthologs, mainly from bacteria but also from fungi and plants, suggesting that the rate of HGT into bdelloid genomes is at least an order of magnitude greater than in other eukaryotes. Interestingly, many genes derived from nonmetazoans via HGT are expressed and known to be functional.

Unfortunately, most published methods so far depend on just a few (if not just one) types of metrics, rendering them susceptible to errors. Plain aligners have limitations and can produce unwanted artifacts. It is important to note that HGT is an evolutionary hypothesis, and that a sequence alignment alone is insufficient to confirm such a hypothesis. For HGT identification and contamination filtration, phylogenetic and parametric methods are commonly used. The phylogenetic approach is based on incongruent relationships among taxa, whereas the parametric method is based on changes in nucleotide composition and homology search. In this study, we used coverage, relative GC content, and best similarity match in public databases, expression data, and many other parameters to categorize genome data as clean, contaminated, or HGT DNA sequences.  

## Installation

### Manual Installation
0. Install perl module dependencies using CPAN.
1. Download or clone the <code>alienomics</code> tool. Consider unzipping in a directory to <code>alienomics/</code>.
2. <code>cd</code> into <code>alienomics/bin/</code> directory.
3. Type <code>perl alienomics_v0.3.pl</code>.

### Conda Installation (Recommended)
It is recommended to create an environment to avoid any conflits. 
```
conda create -n alienomics
conda activate alienomics
```
Now install alienomics from conda server
```
conda install -c jnarayan81 alienomics
```
Conda installation automatically manage all perl module dependencies.  

### Alienomics Usage
#### The script alienomics_v0.3.pl / alienomics_v0.3
The alienomics framework, alienomics contains one main script called alienomics_v0.3.* This script coordinates the execution of the whole process of HGT and contaminant detection framework. The following section describes the mandatory and optional flags that need to be used at the time of calling the program.
The script alienomics_v0.3.pl has the following command line:  
```
alienomics_v0.3 -c config/config.aom 
OR 
alienomics_v0.3 --conf config/config.aom
```
The primary flags are as follows:  
- <code>[-c | --conf] <filepath>:</code> Provide the location of alienomics configuration file.    
- <code>[-dd | --database_download]:</code> Download the third party database from server.   
- <code>[-cc | --create_config]:</code> Create the sample config file to test run alienomics.   
- <code>[-w | --version]:</code> Know about the developers.    
- <code>[-v | --version]:</code> You can check of the version of the the alienomics program.    
  
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

## Frequently Asked Questions

### How to download third party database in alienomics?

You need to activate -dd or --download_database flags in alienomics:

```
alienomics_v.0.3 -dd -ddf BURT -ddl /home/download/
```
Here download database flags(ddf) are read as follows: B->Busco, U->Uniref50, R->rRNA, and T->Taxonomy Databases. These are case sensitive.

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

As always, you can contact the authors at <jnarayan81@gmail.com>, Paul Simion <polo.simion@gmail.com>, / Alexandre Mayer / Pr. Karine Van Doninck 

## Citation
Feel free to cite the followings:  

Homologous chromosomes in asexual rotifer Adineta vaga suggest automixis (2020), Paul Simion, Jitendra Narayan, Antoine Houtain, Alessandro Derzelle, Lyam Baudry, Emilien Nicolas, Marie Cariou, Nadège Guiglielmoni, Djampa KL Kozlowski, Florence Rodriguez Gaudray, Matthieu Terwagne, Julie Virgo, Benjamin Noel, Patrick Wincker, Etienne GJ Danchin, Martial Marbouty, Bernard Hallet, Romain Koszul, Antoine Limasset, Jean-François Flot, Karine Van Doninck doi: https://doi.org/10.1101/2020.06.16.155473  
