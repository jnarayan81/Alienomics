package AOM::rRNAHandler;

#Handle the rRNA blast 

#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Text::CSV;
use IO::Handle;
use List::MoreUtils qw( minmax );
use AOM::commonSubs;

sub rRNABlaster {
my ($param, $geneSet)=@_;

#Create local blastDB
#remove the space from header file with "_";
#correctFastaHead("$param->{rRNA_file}", "\s+", "_");

   print "Converting RNA to DNA alphabet ....\n";
   RNA2DNAalphabet ("$param->{rRNA_file}", "$param->{project_dir_path}/intermediate_files/rRNA/rRNA2DNA"); #Formated alphabet to DNA

   print "Creating blast rRNADB for you with updated rRNA file ....";
   system ("$param->{makeblastdb_path} -in $param->{project_dir_path}/intermediate_files/rRNA/rRNA2DNA -parse_seqids -dbtype nucl -out $param->{project_dir_path}/localDB/rRNADB/myrRNADB");

	#For local blast using $param->{blastdb_path} data --- with GENE
	print "Blasting genes onto rRNADB\n";
	system ("$param->{blastn_path} -task $param->{blast_task} -query $param->{project_dir_path}/intermediate_files/bed/gene.fa -db $param->{project_dir_path}/localDB/rRNADB/myrRNADB -evalue $param->{evalue} -num_threads $param->{max_processors} -max_target_seqs 1 -max_hsps 1 -qcov_hsp_perc $param->{qcovper} -outfmt '6 qseqid sallseqid qstart qend sallseqi sstart send evalue length frames qcovs bitscore' -out $param->{project_dir_path}/intermediate_files/rRNA/rRNAGene.megablast");

#Added the detail of the genes in megablast
print "Additng info in rRNA blasthit ....\n";
AOM::commonSubs::addGeneDetail("$param->{project_dir_path}/intermediate_files/rRNA/rRNAGene.megablast", $geneSet, "$param->{project_dir_path}/intermediate_files/rRNA/rRNAGene.megablast.added");

#Store it in hash
my ($rRNAHash_ref, $min, $max)=rRNAHasher("$param->{project_dir_path}/intermediate_files/rRNA/rRNAGene.megablast.added");
return ($rRNAHash_ref, $min, $max);
}

sub addGeneDetail {
my ($file, $geneSet, $out) = @_;
my %gene=%$geneSet;
open(FILE, "<$file") || die "File not found $!";
open(OUT, ">$out") || die "File not found $!";
while (<FILE>) {
 my @values = split /\t/, $_;
 if ($gene{$values[0]}) {
	print OUT "$_\t$gene{$values[0]}\n";
 }
}
close(FILE);
close(OUT);
}


sub correctFastaHead {
my ($file, $fstring, $rstring) = @_;
#open FH, '+<', "$file" or die "Error:$!\n";
#while (<FH>) { if (/^>/) { $_=~ s/\s+/_/g;} print FH; }
#close FH;

open(FILE, "<$file") || die "File not found $!";
my @lines = <FILE>;
close(FILE);

open(FILE, ">$file") || die "File not found $!";
foreach(@lines) {
   chomp;
   s/\s+/_/g if /^>/; #Fixed for now
   print FILE "$_\n";
}
}

#rRNA Handler
sub rRNAHasher {
my $refGeneMBlast = shift;
my %rRNAHash; my @allScore; my $min=0; my $max=0;
open my $fh, "<", $refGeneMBlast or die "$refGeneMBlast: $!";
my $csv = Text::CSV->new ({
      sep_char=>"\t",
      binary    => 1, # Allow special character. Always set this
      auto_diag => 1, # Report irregularities immediately
      });
$fh = IO::Handle->new_from_fd( $fh, 'r');
while ( not $fh->eof ) {
  my $row = $csv->getline( $fh );
  $rRNAHash{$row->[0]} = $row;
  push @allScore, $row->[10];
  #warn Dumper $row; # To see the structure
}
($min, $max) = minmax @allScore;
return (\%rRNAHash, $min, $max);
}

sub RNA2DNAalphabet {
my ($inFile, $outFile) = @_;
open(INFILE, "<$inFile") || die "File not found $!";
open(OUTFILE, ">$outFile") || die "File not found $!";
while (<INFILE>) {
 chomp;
 if ($_ =~ /^>/) {print OUTFILE "$_\n"; next;};
 $_=~tr/U/T/;
 print OUTFILE "$_\n";
}
close(INFILE);
close(OUTFILE);
}

1;
