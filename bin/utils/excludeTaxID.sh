#!/bin/bash

#########
######### Exclude specified taxonomic ID from UNIREF50 fasta file
#########

######## Usage
# bash excludeTaxID.sh -f uniref50_fasta_file -e list_of_taxIDs
#########

# list_of_taxIDs format: a string corresponding to a list of taxonomic IDs seperated by "|" (e.g. "104782|10195|96448|249248|1813166|104781")

# color management
green=`tput setaf 2`
reset=`tput sgr0`

# getting user parameter
while getopts "f:e:" option; do
	case $option in
		f)
			uniref50=$OPTARG;
			#shift; if [ -n "$1" ]; then uniref50=$1; echo -e "uniref50 = $uniref50"; fi; shift;
			;;
		e)
			excludeTaxID=$OPTARG;
			#shift; if [ -n "$1" ]; then excludeTaxID=$1; echo -e "excludeTaxID = $excludeTaxID"; fi; shift;
			;;
	esac
done
#echo -e "${green}Preparing $uniref50 by excluding the following taxon IDs:${reset}"
#echo -e "${green}$excludeTaxID\n${reset}"

outname=taxid_ignored
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' $uniref50 > $uniref50.$outname
echo $excludeTaxID | sed 's/|/\n/g' | while read taxid; do 
	echo -e "   discarding taxon $taxid"
	cat $uniref50.$outname | sed '/^$/d' | paste - - | LC_ALL=C grep -v -F -e "TaxID=$taxid" | sed 's/\t/\n/g' > $uniref50.$outname.tmp
	rm -f $uniref50.$outname; mv $uniref50.$outname.tmp $uniref50.$outname
done

