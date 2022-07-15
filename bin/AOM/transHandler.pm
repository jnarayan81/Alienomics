package AOM::transHandler;

#Handle the transcriptomics data wih kallisto

#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw( minmax );
use AOM::commonSubs;

sub transKallister {
my ($param)=@_;

print "RNAseq_location = $param->{RNAseq_location}\n";
if ( length  $param->{RNAseq_location}) {
	print "\nCreating kallisto index\n";
	system ("$param->{kallisto_path} index --make-unique -i $param->{output_folder}/intermediate_files/gff/kallisto_index_outfile $param->{output_folder}/intermediate_files/gff/local.cds.fasta > /dev/null 2>&1");

	my ($reads_loc,$readsformat) = storeReadsName ("$param->{RNAseq_location}", $param); 

	print "\nQuantifying gene expression using kallisto\n";
	if ( $readsformat eq "paired" ) {
		print "$param->{kallisto_path} quant -t $param->{max_processors} -i $param->{output_folder}/intermediate_files/gff/kallisto_index_outfile -o $param->{output_folder}/intermediate_files/kallistoResult $reads_loc > /dev/null 2>&1";
		system ("$param->{kallisto_path} quant -t $param->{max_processors} -i $param->{output_folder}/intermediate_files/gff/kallisto_index_outfile -o $param->{output_folder}/intermediate_files/kallistoResult $reads_loc > /dev/null 2>&1");
	}
	elsif ( $readsformat eq "single" ) {
		print "$param->{kallisto_path} quant --single -l 100 -s 10 -t $param->{max_processors} -i $param->{output_folder}/intermediate_files/gff/kallisto_index_outfile -o $param->{output_folder}/intermediate_files/kallistoResult $reads_loc";    
		system ("$param->{kallisto_path} quant --single -l 100 -s 10 -t $param->{max_processors} -i $param->{output_folder}/intermediate_files/gff/kallisto_index_outfile -o $param->{output_folder}/intermediate_files/kallistoResult $reads_loc"); ;
	}
}
else  {print "Kallisto fail !! Did you forgot to set rna_seq value in config file?\n"; exit; }

#Extract the value in hash
print "Storing and sorting gene expression level\n";
#`awk < $param->{output_folder}/intermediate_files/kallistoResult/abundance.tsv '{ split($1,f1,/_/) ; printf("%s_%s\t%s\t%s\t%s\t%s\n",f1[1],f1[2],$2,$3,$4,$5) }' > $param->{output_folder}/intermediate_files/kallistoResult/abundance_dups.tsv"`;
print "step1\n";
fileUpdator ("$param->{output_folder}/intermediate_files/kallistoResult/abundance.tsv","$param->{output_folder}/intermediate_files/kallistoResult/abundance_dups.tsv");
print "step2\n";
system("sort -k5rn $param->{output_folder}/intermediate_files/kallistoResult/abundance_dups.tsv | sort -u -k1,1 > $param->{output_folder}/intermediate_files/kallistoResult/abundance_edited.tsv");
print "step3\n";
my ($transHash_ref)=transHasher("$param->{output_folder}/intermediate_files/kallistoResult/abundance_edited.tsv");
print "step4\n";
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
  #print "$row";
  #No underscore should be in fasta headers --- we inserted one at transcriptomics  
  #my @tmp=split ("_", $row->[0]);
  ####" trying to not use an underscore
  my @tmp=split ("_1\t", $row->[0]);
  #print "row[0] = $row->[0] ; tmp = @tmp";
  #fasta header name pattern is crucial here: chr1_23-234 
  #-- ask user to provide genome.fasta header without any underscore or special character 
  # trying to use ":" instead of "_"
  print $ofh "$tmp[0]\t$row->[1]\t$row->[2]\t$row->[3]\t$row->[4]\n";
  #print $ofh "$tmp[0]_$tmp[1]\t$row->[1]\t$row->[2]\t$row->[3]\t$row->[4]\n";

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

# Store names of read files (in order: R1 R2 R1 R2)
sub storeReadsName {
my ($dir,$loc) = @_;
my @files = glob "$dir/*.fastq.gz"; #The fq will not be accepted ! only fastq.gz
#print "@files\n";
my @files_sorted = sort { lc($a) cmp lc($b) } @files;
my $allNamesPath=join(' ',@files_sorted);
print "\nallNamesPath: $allNamesPath\n";
# if paired-end
if (scalar(@files_sorted) % 2 == 0) {
	print "\nUsing following paired-end read files:\n";
	for (0..$#files_sorted){
		print "$files[$_]\n";
		#$files[$_] =~ s/\.txt$//;
	}
	my $readsformat="paired";
	return $allNamesPath,$readsformat;
}
# else if single-end
else { 
	my $file = glob "$dir/*.fastq.gz";
	print "\nUsing following single-end read files:\n";
	for (0..$#files_sorted){
                print "$files[$_]\n";
        }
	my $readsformat="single";
        return $allNamesPath,$readsformat;
}

}

1;
