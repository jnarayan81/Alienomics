
package AOM::HGTcategorization;

use strict;
use warnings;

#Turning off BioPerl warnings
$Bio::Root::Root::DEBUG = -1;


#################################################
sub categorize_genes {

my ($Scaff_file, $Gene_file, $scaf_low_th, $gene_low_th, $scaff_up_th, $gene_up_th, $GFF_OUTPUT) = @_;
#my $filename2 = $ARGV[0]; #Scaff scores
#my $filename = $ARGV[1]; #Gene scores
#my $cThresh= $ARGV[2]; #Contig threshold
#my $gThresh = $ARGV[3]; #Gene threshold
#my $cuThresh= $ARGV[4]; #Contig up threshold
#my $guThresh = $ARGV[5]; #Gene up threshold

my %allHash;
open(my $fh2, '<:encoding(UTF-8)', $Scaff_file) or die "Could not open file '$Scaff_file' $!";
while (my $row2 = <$fh2>) {
  chomp $row2;
  #print "$row2\n";
  my @tmpArray2=split ('\t', $row2);
  #next if $tmpArray[0] ne $tmpArray2[0]; #Ignore the rest
  $allHash{$tmpArray2[0]} = $tmpArray2[13];    # look for scaffold_scores columns:     colonne 13 = objet [12] dans perl ; colonne 14 = objet [13] dans perl
}
close $fh2;

#Lets read the gene file
open(my $fh, '<:encoding(UTF-8)', $Gene_file) or die "Could not open file '$Gene_file' $!";
#print $GFF_OUTPUT "Scaffold\tlength\tGC_scaf\tgene\tRef_score\tgeneral_blast_score\tBusco_score\trRNA_score\tfeeder_score\tTNF_score\texpected_coverage\tcov_score\tGC\tGC_score\tgene_scoreFinal_sigScore\tScaffoldScore\tGeneStatus\t[ScaffoldScore(ScafContamFlag,ScafUndefFlag)]\t[GeneScore(GeneAlienFlag,GeneUndefFlag)]\n";
print $GFF_OUTPUT "Scaffold\tSource\tFeature\tstart\tend\tgene_score\tstrand\tframe\tattributes[RefScore;BlastScore;BuscoScore;rRNAScore;feederScore;TNFScore;CoverageScore(expected_coverage);GCScore(GC%);ExpressionScore(TPM);RawGeneScore]\n";

while (my $row = <$fh>) {
	chomp $row;
	#print "$row\n";
	next if $. == 1; #Ignore header
	my @tmpArray=split ('\t', $row);

	my $fScore=$allHash{$tmpArray[0]};	# scaffold score becomes "$fScore"
	my $gScore=$tmpArray[15];		# gene score becomes "$gScore"

	my $scaffname = $tmpArray[0];
	my $source = "Alienomics";
	my @scaf_startstop = split(':',$tmpArray[3]);	# getting "scaff" and "start-stop" fields
	my @start_stop = split('-',$scaf_startstop[1]);	# getting "start" and "stop" fields
	my $start = $start_stop[0];			# getting "start"
	my $end = $start_stop[1];			# getting "end"
	my $EXP = $tmpArray[17];
	my $attributes = "RefS=$tmpArray[4]; BlastS=$tmpArray[5]; BuscoS=$tmpArray[6]; rRNAS=$tmpArray[7]; feederS=$tmpArray[8]; TNFS=$tmpArray[9]; CoverageS=$tmpArray[11]($tmpArray[10]); GCS=$tmpArray[13]($tmpArray[12]); EXPS=$tmpArray[17]($tmpArray[16]); RawGeneS=$tmpArray[14]; GeneSforScaf=$tmpArray[18]";

	### Flaging different conditions
	my $cFlag=0; my $gFlag=0; my $cuFlag=0; my $guFlag=0;
	# contigs/scaffolds
	if ($fScore <= $scaf_low_th) { $cFlag = 1;} else { $cFlag = 0;}  				# => "1" if scaffold is contaminant
	if ($fScore >= $scaf_low_th && $fScore <= $scaff_up_th) { $cuFlag = 1;} else { $cuFlag = 0;}  	# => "1" if scaffold undefined
	# genes
	if ($gScore <= $gene_low_th) { $gFlag = 1;} else  { $gFlag = 0;}  			# => "1" if gene alien
	if ($gScore >= $gene_low_th && $gScore <= $gene_up_th) { $guFlag = 1;} else  { $guFlag = 0;} # => "1" if gene undefined
	# initialising decision
	my $decision = 'NA';

	#contig alien + gene alien = CONTAMINANT
	if ($cFlag == 1 && $gFlag == 1) { $decision = 'contaminant';}
	#contig alien + gene Undef = unknown
	elsif ($cFlag == 1 && $gFlag == 0 && $guFlag == 1) { $decision = 'unknown';}
	#contig alien + gene self = mis-assembly ?
	elsif ($cFlag == 1 && $gFlag == 0 && $guFlag == 0) { $decision = 'mis-assembled self';}

	#contig Undef + gene alien = unknown alien
	elsif ($cFlag == 0 && $cuFlag == 1 && $gFlag == 1) { $decision = 'unknown alien';}
	#contig_Undef + gene Undef = unknown
	elsif ($cFlag == 0 && $cuFlag == 1 && $gFlag == 0 && $guFlag == 1) { $decision = 'unknown';}
	#contig Undef + gene self = unknown self
	elsif ($cFlag == 0 && $cuFlag == 1 && $gFlag == 0 && $guFlag == 0) { $decision = 'unknown self';}

	#contig self + gene alien = HGT
	elsif ($cFlag == 0 && $cuFlag == 0 && $gFlag == 1 && $EXP == -1) { $decision = 'unexpressed HGT';}
	elsif ($cFlag == 0 && $cuFlag == 0 && $gFlag == 1 && $EXP != -1) { $decision = 'expressed HGT';}
	#contig self + gene Undef = unknown
	elsif ($cFlag == 0 && $cuFlag == 0 && $gFlag == 0 && $guFlag == 1) { $decision = 'unknown';}
	#contig self + gene self = SELF
	elsif ($cFlag == 0 && $cuFlag == 0 && $gFlag == 0 && $guFlag == 0) { $decision = 'self';}

	#print $GFF_OUTPUT "$row\t$fScore\t$decision\t[$fScore($cFlag,$cuFlag)]\t[$tmpArray[15]($gFlag,$guFlag)]\n";
	print $GFF_OUTPUT "$scaffname\t$source\t$decision\t$start\t$end\t$gScore\t.\t.\t$attributes\n";

}

}


