use strict;
use warnings;

my $filename2 = $ARGV[0]; #Scaff
my $filename = $ARGV[1]; #Gene
my $cThresh= $ARGV[2]; #Contig threshold
my $gThresh = $ARGV[3]; #Gene threshold
my $cuThresh= $ARGV[4]; #Contig up threshold
my $guThresh = $ARGV[5]; #Gene up threshold

my %allHash;
open(my $fh2, '<:encoding(UTF-8)', $filename2) or die "Could not open file '$filename2' $!";
while (my $row2 = <$fh2>) {
  chomp $row2;
  #print "$row2\n";
  my @tmpArray2=split ('\t', $row2);
  #next if $tmpArray[0] ne $tmpArray2[0]; #Ignore the rest
  $allHash{$tmpArray2[0]} = $tmpArray2[12];
}
close $fh2;

#Lets read the gene file
open(my $fh, '<:encoding(UTF-8)', $filename) or die "Could not open file '$filename' $!";
while (my $row = <$fh>) {
  chomp $row;
  #print "$row\n";
  next if $. == 1; #Ignore header
  my @tmpArray=split ('\t', $row);
  my $fScore=$allHash{$tmpArray[0]}; #Contig score is $fScore

### All conditions here
my $cFlag=0; my $gFlag=0; my $cuFlag=0; my $guFlag=0;
# contigs/scaffolds
if ($fScore <= $cThresh) { $cFlag = 1;} else { $cFlag = 0;}  				# => "1" if scaffold is contaminant
if ($fScore >= $cThresh && $fScore <= $cuThresh) { $cuFlag = 1;} else { $cuFlag = 0;}  	# => "1" if scaffold undefined
# genes
if ($tmpArray[6] <= $gThresh) { $gFlag = 1;} else  { $gFlag = 0;}  			# => "1" if gene alien
if ($tmpArray[6] >= $gThresh && $tmpArray[6] <= $guThresh) { $guFlag = 1;} else  { $guFlag = 0;} # => "1" if gene undefined


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
elsif ($cFlag == 0 && $cuFlag == 0 && $gFlag == 1) { $decision = 'HGT';}
#contig self + gene Undef = unknown
elsif ($cFlag == 0 && $cuFlag == 0 && $gFlag == 0 && $guFlag == 1) { $decision = 'unknown';}
#contig self + gene self = SELF
elsif ($cFlag == 0 && $cuFlag == 0 && $gFlag == 0 && $guFlag == 0) { $decision = 'self';}

  print "$row\t$fScore\t$decision\t[$fScore($cFlag,$cuFlag)]\t[$tmpArray[6]($gFlag,$guFlag)]\n";
}
