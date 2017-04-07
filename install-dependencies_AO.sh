#!/bin/bash

# Install all dependencies for Alienomics v0.1
# Commenline those programs which is already present elsewhere
# BlastDB will acquire a lot of space (around 120GB)
# Some of the program version and download link might change over time. Please do check if any error...

# Alienomics v0.1 @ Jitendra Narayan

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
# NCBI ncbi-blast-2.6.0+
if [ ! -d "$PWD/ncbi-blast-2.6.0+/bin" ]; then
    cd $PWD;
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz;
    tar xzfp ncbi-blast-2.6.0+-x64-linux.tar.gz;
else
   echo "NCBI blast+ tool exists"
fi


# augustus.2.5.5
if [ ! -d "$PWD/augustus-3.2.3" ]; then
    cd $PWD;
    wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.2.3.tar.gz
    tar xzfp augustus-3.2.3.tar.gz;
else
   echo "augustus.2.5.5 tool exists"
fi

: '
# BlastDB download
if [ ! -d "$PWD/blastDB" ]; then
    if [ ! -d "$PWD/blastDB" ]; then
        mkdir $PWD/blastDB
    fi
    cd $PWD/blastDB
    wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz"
    for a in nt.*.tar.gz; do tar xzf $a; done
else
   echo "blastDB database exists! If not installed sucessfully, delete the folder and re-run it"
fi
'

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
   echo "taxdump database exists"
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
   echo "rRNAData database exists"
fi


# Download bwakit
if [ ! -d "$PWD/bwa.kit" ]; then
    cd $PWD
    wget https://downloads.sourceforge.net/project/bio-bwa/bwakit/bwakit-0.7.12_x64-linux.tar.bz2
    tar -xjf bwakit-0.7.12_x64-linux.tar.bz2
else
   echo "bwa.kit exists"
fi

