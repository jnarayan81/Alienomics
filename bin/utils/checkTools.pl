use strict;
use warnings;
use 5.010;
use Term::ANSIColor;

my @toolNames = ("blastn", "blastx", "diamond", "makeblastdb", "augustus", "bwa", "seqtk", "samtools", "kallisto");
my @notInstalled;
foreach my $name (@toolNames) {
	if (system("which $name")){
		print "Can't find $name in your system\n";
		push (@notInstalled, $name);
	} else {
		print "$name is installed\n";
	}
}
if (!@notInstalled) {print "Congratulations, all necessary programs are installed for Alienomics\n";}
else { 
	foreach (@notInstalled) {                        
		  print colored(['red on_bright_yellow'], "\n==>$_ NOT installed<==", "\n");
	};
exit;
}

