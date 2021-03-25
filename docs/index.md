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

Unfortunately, most published methods so far rely on only a few (if not only one) type of metrics, making them prone to errors. The use of plain aligners have limits and can result in unwanted artifacts. One must remember that HGT is an evolutionary hypothesis, and a sequence alignment in itself is not enough to validate such hypothesis. Mostly phylogenetic and parametric method are used for HGT detection and contamination filtration. The phylogenetic approach relies on incongruent relationships among taxa, whereas shift in nucleotide composition, and homology search are utilised by parametric method. Here we tried coverage, relative GC content and best similarity match in public databases, expression data and many more parametroc threaholds to divide the genome data into clean, contaminated and HGT DNA sequences.

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
alienomics_v0.3 -c config/test.aom 
OR 
alienomics_v0.3 --conf config/test.aom
```
The primary flags are as follows:  
- <code>[-c | --conf] <filepath>:</code> Provide the location of alienomics configuration file.   
- <code>[-dd | --database_download]:</code> Download the third party database from server.   
- <code>[-cc | --create_config]:</code> Create the sample config file to test run alienomics.   
- <code>[-w | --version]:</code> Know about the developers.    
- <code>[-v | --version]:</code> You can check of the version of the the alienomics program.    
  
All the mandatory and optional parameters are needed to be set in config.aom file. 

#### Mandatory parameters



#### Optional parameters
[Optional] Add more options if required or would like to interfere with the default options. 
 


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

As always, you can contact the authors at <jnarayan81@gmail.com>.

## Citation
Feel free to cite the followings:  

Homologous chromosomes in asexual rotifer Adineta vaga suggest automixis (2020), Paul Simion, Jitendra Narayan, Antoine Houtain, Alessandro Derzelle, Lyam Baudry, Emilien Nicolas, Marie Cariou, Nadège Guiglielmoni, Djampa KL Kozlowski, Florence Rodriguez Gaudray, Matthieu Terwagne, Julie Virgo, Benjamin Noel, Patrick Wincker, Etienne GJ Danchin, Martial Marbouty, Bernard Hallet, Romain Koszul, Antoine Limasset, Jean-François Flot, Karine Van Doninck doi: https://doi.org/10.1101/2020.06.16.155473  
