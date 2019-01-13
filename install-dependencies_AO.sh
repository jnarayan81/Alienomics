#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#!/bin/bash

# Install all dependencies for Alienomics v0.1
# Commenline those programs which is already present elsewhere
# BlastDB will acquire a lot of space (around 120GB)
# Some of the program version and download link might change over time. Please do check if any error...

# Alienomics v0.1 @ Jitendra Narayan

echo -e "\nInstalling perl modules, if not exists ---\n"

allMods=(
Cwd
File::chdir
File::Copy
POSIX
Tie::File
Try::Tiny
Data::Dumper
File::Basename
Bio::SeqIO
FindBin
File::Remove
Capture::Tiny
File::Temp
File::Spec::Functions
Statistics::Multtest
File::Path
Statistics::Distributions
Getopt::Long
Statistics::R
Math::Round
File::Find
Bio::DB::Taxonomy
Pod::Usage
);

for i in "${allMods[@]}"
do
	m=$(perl -M$i -e 1 2>&1)
	rc=$?
	if [[ $rc != 0 ]]
		then
		echo "Installing Perl module $i"
		perl -MCPAN -e "install $i"
	fi
	
	m=$(perl -M$i -e 1 2>&1)
	rc1=$?

	if [[ $rc1 != 0 ]]
		then
		echo ""
		echo "------------============--------------"
		echo "Perl $i modules NOT installed successfully. Please make sure that Perl and CPAN are available on your system."
		echo "For Ubuntu and Debian users: you can use 'apt-get install perl' to get newest Perl and CPAN."
		echo "------------============---------------"
	else
		echo -e "exists \t $i"
	fi

done

echo -e "\nInstalling third party database and tools\n"

set -e

: '
# prinseq
if [ ! -e "$PWD/prinseq/prinseq-lite.pl" ]; then
    if [ ! -d "$PWD/prinseq" ]; then
        mkdir $PWD/prinseq
    fi
    cd /tmp
    wget http://downloads.sourceforge.net/project/prinseq/standalone/prinseq-lite-0.20.4.tar.gz \
    -O /tmp/prinseq-lite-0.20.4.tar.gz;
    tar -xvf prinseq-lite-0.20.4.tar.gz;
    install -v prinseq-lite-0.20.4/prinseq-lite.pl $PWD/prinseq;
else
    echo "prinseq tool exists"
fi
'
# NCBI ncbi-blast-2.8.1+
if [ ! -d "$PWD/ncbi-blast-2.8.1+/bin" ]; then
    cd $PWD;
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/ncbi-blast-2.8.1+-x64-linux.tar.gz;
    tar xzfp ncbi-blast-2.8.1+-x64-linux.tar.gz;
else
   echo -e "exists \t NCBI blast+ tool"
fi


# augustus.2.5.5
if [ ! -d "$PWD/augustus.current" ]; then
    cd $PWD;
    wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus.current.tar.gz
    tar xzfp augustus.current.tar.gz;
else
   echo -e "exists \t augustus.current tool"
fi


# BlastDB download
if [ ! -d "$PWD/blastDB" ]; then
    if [ ! -d "$PWD/blastDB" ]; then
        mkdir $PWD/blastDB
    fi
    cd $PWD/blastDB
    wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz"
    for a in nt.*.tar.gz; do tar xzf $a; done
    #Out of the directory 
    cd ..
else
   echo -e "exists \t blastDB, If not installed sucessfully, delete the folder and re-run it"
fi


# taxdump download
if [ ! -d "$PWD/taxdump" ]; then
    if [ ! -d "$PWD/taxdump" ]; then
        mkdir $PWD/taxdump
    fi
    cd $PWD/taxdump
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar xzfp taxdump.tar.gz;
    cd ..
else
   echo -e "exists \t taxdump database"
fi


# SILVA rRNA database
# Change your database name as per you need
# https://www.arb-silva.de/no_cache/download/archive/release_128/Exports/
# SILVA_128_LSUParc_tax_silva.fasta.gz
# SILVA_128_LSUParc_tax_silva_trunc.fasta.gz	
# SILVA_128_LSURef_tax_silva.fasta.gz
# SILVA_128_LSURef_tax_silva_full_align_trunc.fasta.gz
# SILVA_128_LSURef_tax_silva_trunc.fasta.gz
# SILVA_128_SSUParc_tax_silva.fasta.gz
# SILVA_128_SSUParc_tax_silva_trunc.fasta.gz
# SILVA_128_SSURef_Nr99_tax_silva.fasta.gz
# SILVA_128_SSURef_Nr99_tax_silva_full_align_trunc.fasta.gz
# SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz
# SILVA_128_SSURef_tax_silva.fasta.gz
# SILVA_128_SSURef_tax_silva_full_align_trunc.fasta.gz
# SILVA_128_SSURef_tax_silva_trunc.fasta.gz

if [ ! -d "$PWD/rRNAData" ]; then
    if [ ! -d "$PWD/rRNAData" ]; then
        mkdir $PWD/rRNAData
    fi
    cd $PWD/rRNAData
    wget https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/SILVA_128_LSUParc_tax_silva.fasta.gz
    gzip -d SILVA_128_LSUParc_tax_silva.fasta.gz
    #for a in *.fasta.gz; do gzip $a; done
    cd ..
else
   echo -e "exists \t rRNAData database"
fi


# Download bwakit
if [ ! -d "$PWD/bwa.kit" ]; then
    cd $PWD
    wget https://downloads.sourceforge.net/project/bio-bwa/bwakit/bwakit-0.7.12_x64-linux.tar.bz2
    tar -xjf bwakit-0.7.12_x64-linux.tar.bz2
else
   echo -e "exists \t bwa.kit "
fi


: '
Use database of your interest
For more detail visit file:///home/urbe/Downloads/BUSCO_v2.0_userguide.pdf
# Bacteria
wget http://busco.ezlab.org/v2/datasets/bacteria_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/proteobacteria_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/rhizobiales_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/betaproteobacteria_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/gammaproteobacteria_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/enterobacteriales_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/deltaepsilonsub_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/actinobacteria_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/cyanobacteria_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/firmicutes_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/clostridia_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/lactobacillales_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/bacillales_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/bacteroidetes_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/spirochaetes_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/tenericutes_odb9.tar.gz
        
# Eukaryota
wget http://busco.ezlab.org/v2/datasets/eukaryota_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/fungi_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/microsporidia_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/dikarya_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/ascomycota_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/pezizomycotina_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/eurotiomycetes_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/sordariomyceta_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/saccharomyceta_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/saccharomycetales_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/basidiomycota_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/metazoa_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/nematoda_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/arthropoda_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/insecta_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/endopterygota_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/hymenoptera_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/diptera_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/vertebrata_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/actinopterygii_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/tetrapoda_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/aves_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/mammalia_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/euarchontoglires_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/laurasiatheria_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/embryophyta_odb9.tar.gz
wget http://busco.ezlab.org/v2/datasets/protists_ensembl.tar.gz
wget http://busco.ezlab.org/v2/datasets/alveolata_stramenophiles_ensembl.tar.gz

'

if [ ! -d "$PWD/BUSCOData" ]; then
    if [ ! -d "$PWD/BUSCOData" ]; then
        mkdir $PWD/BUSCOData
    fi
    cd $PWD/BUSCOData
    wget http://busco.ezlab.org/v2/datasets/metazoa_odb9.tar.gz
    tar xzfp metazoa_odb9.tar.gz
    #for a in *.fasta.gz; do gzip $a; done
else
   echo -e "exists \t BUSCOData database"
fi

echo -e "\nAll Done\n"

echo -e "\nInstall R package install.packages('Ckmeans.1d.dp') \n"
