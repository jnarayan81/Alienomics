package AOM::commonSubs;

#Add gene detail
sub addGeneDetail {
my ($file, $gene, $out) = @_;
my %gene=%$gene;
open(FILE, "<$file") || die "File not found $!";
open(OUT, ">$out") || die "File not found $!";
while (<FILE>) {
 chomp;
 my @values = split /\t/, $_;
 if ($gene{$values[0]}) {
	print OUT "$_\t$gene{$values[0]}{len}\t$gene{$values[0]}{gc}\n";
 }
 else { print "Did not find $values[0] in the gene dictionnary (see line below)\n\t$_\n";}
}
close(FILE);
close(OUT);
}

#Check all the  mandatory programs to run Alienomics
sub checkPrograms {
use strict;
use warnings;
use 5.010;
use Term::ANSIColor;

#List all tool name here to check
my @toolNames = ("blastn", "blastx", "diamond", "makeblastdb", "augustus", "bwa", "seqtk", "samtools", "kallisto");
my @notInstalled;

foreach my $name (@toolNames) {
	if (system("which $name")){
		print colored(['red on_bright_yellow'],"Can't find $name in your system", "\n");
		push (@notInstalled, $name);
	} else {
		#print "$name is installed\n";
	}
}
if (!@notInstalled) {print "Done\n";}
else { 
	if (@notInstalled) { print colored(['red on_bright_yellow'],"\nFollowing tools are missing in your system path:", "\n");}
	foreach (@notInstalled) {                        
		  print colored(['red on_bright_yellow'], "==>$_ NOT installed<==", "\n");
	};
exit;
}
}


1;
