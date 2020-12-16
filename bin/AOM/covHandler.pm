#Handle the mapping coverage
package AOM::covHandler;

use strict;
use Parallel::ForkManager;
use Fcntl qw(:flock);

#extractMappingCov($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]);
#my $covHash=storeInHash($ARGV[3]);
#my %covHash=%$covHash;


sub extractMappingCov {
my ($bam, $genebed, $maxCore, $outFile) = @_;
#print "$bam, $genebed, $maxCore\n";
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
    storeInFile($bam, "$lines{$lCntNo[$child]}", $outFile); #print "$geneCov-$accession_number\n";
    $pm->finish($child); # pass an exit code to finish
  }
  print "Waiting for all the annotation jobs to complete...\n";
  $pm->wait_all_children;
  print "Annotation DONE ... Everybody is out of the computation pool!\n\n\n";

}

#Store the coverage information in hash
sub storeInFile {
my ($bam, $line, $outFile)=@_;
my @localVal = split(/\t/, $line);
s{^\s+|\s+$}{}g foreach @localVal; #Delete white space
my $accession_number="$localVal[0]:$localVal[1]-$localVal[2]";

my $geneCov = `samtools depth -aa -r $accession_number $bam | awk '{ c++; s+=\$3; } END {k=(s/c); print k;}'`;
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
