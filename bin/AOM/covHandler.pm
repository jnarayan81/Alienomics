#Handle the mapping expected_coverage
package AOM::covHandler;

use strict;
use Parallel::ForkManager;
use Fcntl qw(:flock);

#extractMappingCov($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]);
#my $covHash=storeInHash($ARGV[3]);
#my %covHash=%$covHash;


#############################
sub extractMappingCov_v2 {

my ($input_mapping_bam, $genebed, $maxCore, $outFile) = @_;
print "Computing coverage for all genes:\nbedfile used: $genebed\noutput coverage: $genebed.cov\n";

# reading coverage info from bam file only once (for all genes)
#previous version from paul : system("samtools depth -aa -b $genebed $input_mapping_bam > $genebed.cov");
# version from Jit : `samtools depth -aa -r $accession_number $input_mapping_bam`
# new test
system("samtools bedcov $genebed $input_mapping_bam > $genebed.cov");

# loading coverage at every positions in hash
my %index;
open(COUNTS, "$genebed.cov") or die ("Could not open $genebed.cov\n");
open(OUTFILE, '>>' ,$outFile);
while ( my $line = <COUNTS> ) {

	# previous code for "samtools depth -aa -b"
	#my @localVal = split(/\t/, $line); # 0=chrom ; 1=pos ; 2=cov 
	#$index{$localVal[0]}{$localVal[1]}=$localVal[2]
	
	# new test for "samtools bedcov"
	my @localVal = split(/\t/, $line); # 0=chrom ; 1=start ; 2=stop ; 7=tot bp coverage
	my $chrom=$localVal[0]; my $start=$localVal[1]; my $stop=$localVal[2]; my $totcoverage=$localVal[7];
	my $genelength=$stop-$start;
	my $averagecov=$totcoverage/$genelength;
	print OUTFILE "$chrom:$start-$stop\t$averagecov\n"; 
}
close(COUNTS);
close(OUTFILE);
}

# previous code for "samtools depth -aa -b" (useless now ? all has been commented)
# reading bed and measuring coverage : 
#open(BED, $genebed) or die ("Could not open $genebed\n");
#open(OUTFILE, '>>' ,$outFile);
#while ( my $line2 = <BED> ) {
#	my @localVal = split(/\t/, $line2);
#	my $chrom=$localVal[0]; my $start=$localVal[1]; my $stop=$localVal[2];
#	my $depthstart=$start + 1; my $depthstop=$stop + 1;
		# collecting coverage values
#	my @tmp_array=();
#	for (my $i = $depthstart; $i < $depthstop; $i++) {
#		push(@tmp_array, $index{$chrom}{$i}); # collect coverage values from hash at all corresponding positions
#	}
		# computing average coverage for current gene
#	my $sum=0; my $total=0; my $average=0;
#	foreach (@tmp_array) {
#		$sum=$sum + $_;
#		++$total;
#	}
#	$average=$sum/$total;
#	print OUTFILE "$chrom:$start-$stop\t$average\n";
#}
#close(OUTFILE);

#}


###########################
sub extractMappingCov {

my ($input_mapping_bam, $genebed, $maxCore, $outFile) = @_;
#print "$input_mapping_bam, $genebed, $maxCore\n";
my $lines_ref=storeLineInHash($genebed);
my %lines=%$lines_ref;
my $max_procs = $maxCore;
my @lCntNo = keys %lines;
# hash to resolve PID's back to child specific information
my $pm =  new Parallel::ForkManager($max_procs);
# Setup a callback for when a child finishes up so we can
# get it's exit code
  $pm->run_on_finish (
    sub { my ($pid, $exit_code, $ident) = @_;
      #print "** $ident just got out of the pool ". "with PID $pid and exit code: $exit_code\n";
    }
  );

  $pm->run_on_start(
    sub { my ($pid,$ident)=@_;
     #print "** $ident started, pid: $pid\n";
    }
  );

  $pm->run_on_wait(
    sub {
      #print "** Have to wait for one children ...\n"
    },
    0.5
  );

  NAMES:
  foreach my $child ( 0 .. $#lCntNo ) {
#print "Workinnnnnnnnnn $child \n";
    my $pid = $pm->start($lCntNo[$child]) and next NAMES;
    storeInFile($input_mapping_bam, "$lines{$lCntNo[$child]}", $outFile); #print "$geneCov-$accession_number\n";
    $pm->finish($child); # pass an exit code to finish
  }
  #print "Waiting for all the annotation jobs to complete...\n";
  $pm->wait_all_children;
  print "Done\n";

}

#Store the expected_coverage information in hash
sub storeInFile {
my ($input_mapping_bam, $line, $outFile)=@_;
my @localVal = split(/\t/, $line);
s{^\s+|\s+$}{}g foreach @localVal; #Delete white space
my $accession_number="$localVal[0]:$localVal[1]-$localVal[2]";

my $geneCov = `samtools depth -aa -r $accession_number $input_mapping_bam | awk '{ c++; s+=\$3; } END {k=(s/c); print k;}'`;
chomp($geneCov);
#print "$geneCov ----- $accession_number ----\n";
collision_free_write($outFile, "$accession_number\t$geneCov"); #Write into outfile
}

#Write into outfile -- collision free because of multicore usesage
sub collision_free_write {
  my($outFile, $msg) = @_;

  open my $ofh, ">>", $outFile or die  "$0 [$$]: open: $!";
  flock $ofh, LOCK_EX      or die  "$0 [$$]: flock: $!";
  print $ofh "$msg\n"      or die  "$0 [$$]: write: $!";
  close $ofh               or warn "$0 [$$]: close: $!";
}

sub storeLineInHash{
my $genebed=shift;
  open my $fh, '<', "$genebed";
  my $lCnt; my %lines;
  while (<$fh>) {
        chomp;
        $lCnt++;
        $lines{$lCnt} = $_;
  }
close $fh;
return \%lines;
}

1;
