package AOM::transHandler;

#Handle the transcriptomics data wih kallisto

#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw( minmax );
use AOM::commonSubs;

sub transKallister {
my ($param)=@_;

if ($param->{rna_seq} > 0) {
      print "Creating kallisto index ....\n";
      system ("$param->{kallisto_path} index --make-unique -i $param->{project_dir_path}/intermediate_files/gff/kallisto_index_outfile $param->{project_dir_path}/intermediate_files/gff/local.cds.fasta");
      print "Lets store the rna-seq file names\n";
      my $reads_loc = storeReadsName ("$param->{kallisto_reads}", $param); 
      system ("$param->{kallisto_path} quant -t $param->{max_processors} -i $param->{project_dir_path}/intermediate_files/gff/kallisto_index_outfile -o $param->{project_dir_path}/intermediate_files/kallistoResult $reads_loc");
}
else  {print "Kallisto fail !! Did you forgot to set rna_seq value in config file -- it should be more than zero"; exit; }

#Extract the value in hash
print "Storing trans hits in hash - kallisto based\n";
print "\Filtering the file for duplicated lines -- and keeping the one with best score\n";
#`awk < $param->{project_dir_path}/intermediate_files/kallistoResult/abundance.tsv '{ split($1,f1,/_/) ; printf("%s_%s\t%s\t%s\t%s\t%s\n",f1[1],f1[2],$2,$3,$4,$5) }' > $param->{project_dir_path}/intermediate_files/kallistoResult/abundance_dups.tsv"`;
fileUpdator ("$param->{project_dir_path}/intermediate_files/kallistoResult/abundance.tsv","$param->{project_dir_path}/intermediate_files/kallistoResult/abundance_dups.tsv");
system("sort -k5rn $param->{project_dir_path}/intermediate_files/kallistoResult/abundance_dups.tsv | sort -u -k1,1 > $param->{project_dir_path}/intermediate_files/kallistoResult/abundance_edited.tsv");
my ($transHash_ref)=transHasher("$param->{project_dir_path}/intermediate_files/kallistoResult/abundance_edited.tsv");
return $transHash_ref;
}

#update kallisto out for _*
sub fileUpdator {
my ($kfile,$ofile) = @_;
open my $ifh, "<", $kfile or die "$kfile: $!";
open my $ofh, ">", $ofile or die "$ofile: $!";
my $csv = Text::CSV->new ({
      sep_char=>"\t", # For tab seperation
      binary    => 1, # Allow special character. Always set this
      auto_diag => 1, # Report irregularities immediately
      });
$ifh = IO::Handle->new_from_fd( $ifh, 'r');
while ( not $ifh->eof ) {
  my $row = $csv->getline( $ifh );
  #No underscore should be in fasta headers --- we inserted one at transcriptomics
  my @tmp=split ("_", $row->[0]);
  #fasta header name pattern is crucial here: chr1_23-234 
  #-- ask user to provide genome.fasta header without any underscore or special character 
  print $ofh "$tmp[0]_$tmp[1]\t$row->[1]\t$row->[2]\t$row->[3]\t$row->[4]\n";
};
close $ofh;
}

#Transcriptomics Handler
sub transHasher {
my $kallisto_file = shift;
my %transHash;
open my $fh, "<", $kallisto_file or die "$kallisto_file: $!";
my $csv = Text::CSV->new ({
      sep_char=>"\t", # For tab seperation
      binary    => 1, # Allow special character. Always set this
      auto_diag => 1, # Report irregularities immediately
      });
$fh = IO::Handle->new_from_fd( $fh, 'r');
while ( not $fh->eof ) {
  my $row = $csv->getline( $fh );
#print "$row->[0] ======= $row->[4]\n";
  $transHash{$row->[0]} = $row->[4];
  #warn Dumper $row; # To see the structure
};
return (\%transHash);
}

#Store all the names in order R1 R2 R1 R2
sub storeReadsName {
my ($dir,$loc) = @_;
my @files = glob "$dir/*.fastq.gz";
my @files_sorted = sort { lc($a) cmp lc($b) } @files;
my $allNamesPath=join(' ',@files_sorted);
if (scalar(@files_sorted) % 2 == 0) {print "Seems fine -- even number of files\n";} else { print "You have odd number of files -- most likely not in pair !\n"; exit;}
print "Files used in kallisto \n";
for (0..$#files_sorted){ 
print "$files[$_]\n";
  #$files[$_] =~ s/\.txt$//;
}
return $allNamesPath;
}

1;
