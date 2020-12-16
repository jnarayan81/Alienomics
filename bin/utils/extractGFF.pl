#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: gffextract.pl
#
#        USAGE: ./gffextract.pl  
#
#  DESCRIPTION: 
#       AUTHOR: Felipe Gajardo (), fgajardoe@gmail.com
#      CREATED: 04/14/2015 05:35:05 PM
#
#===============================================================================

use strict;
use warnings;
use utf8;
use Getopt::Std;
use Data::Dumper;
#use Misc; #qw(gff2genemodel printGene filterGeneTranscripts);
use gffRapeLib;

our($opt_k,$opt_h,$opt_l,$opt_i);
# Defaults
$opt_k='all';


# Usage
my $usage="Usage: perl $0 [-k all|longest|shortest, -l genes_to_extract] -i input.gff
Options:
			-i input_gff	Provide input GFF to extract genes (mandatory).
			-l gene_list	Provide a list of genes to extract (optional).
			-k				Transcript to keep. Arguments are 'longest','shortest',
							'all'(default).
			-h				Print this help

";

# Get options
getopts('hk:l:i:');

# Check if there is input
if(($opt_h)||(!$opt_i)){
	print $usage;
	exit 0;
}

# Main

my $genes=gff2genemodel($opt_i);

my @genelst;
if($opt_l){
	open LST,"<$opt_l" or die "Cannot open $opt_l";
	@genelst=<LST>;
}
else{
	@genelst=keys %$genes;
}

#recorreemos la lista de genes
foreach my $idlst(@genelst){
	chomp $idlst;
	my $gene=$genes->{$idlst};
	my $filtered_gene=$gene;

	if($opt_k ne 'all'){
		$filtered_gene=filterGeneTranscripts($gene,$opt_k);
	}
#	if($filtered_gene!=0){
	if($filtered_gene){ # testeando. cambie la linea de arriba por esta. para evitar un warning que me aparecia
		printGene($filtered_gene);
	}
}

close LST;

