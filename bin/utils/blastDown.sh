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
