#!/usr/bin/perl


if (@ARGV < 1) {
        print "usage: 

extract_seqs_from_fasta.pl <list of IDs> <input fasta file> <output fasta file>

<list of IDs>        :  list of IDs comma separated
<input fasta file>   :  fasta file name
<output fasta file>  :  fasta file name\n\n";
}


$IDlist = shift or die "Must provide a list of IDs\n";
$infasta = shift or die "Must provide an input fasta file name\n";
$outfasta = shift or die "Must provide an output fasta file name \n";

my @IDs = split(',',$IDlist);

$| = 1;                         # forces immediate prints into files rather than the buffer.

open (INFASTA, "<$infasta") or die "Could not open $infasta";
open (OUTFASTA, ">$outfasta") or die "Could not open $outfasta";

while (<INFASTA>){
  if (/^\>(\S+).*$/){
    $key = "$1";
    $sequences{$key} .= $_;
  } else {
    $sequences{$key} .= $_;
  }
}
close SOURCE;

foreach $key (keys %sequences) {
  foreach $id (@IDs) {
    if($key eq $id) {
      print OUTFASTA "$sequences{$key}\n";
      last;
    }
  }
}
