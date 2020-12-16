use strict;
use warnings;

#Check opverls between bed files @ Jitendra

####
# Read BED file into an array
####
my $BED_FILE = <<EOBED;
chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
chr7    127472363  127473530  Pos2  0  +  127472363  127473530  255,0,0
chr7    127473530  127474697  Pos3  0  +  127473530  127474697  255,0,0
chr7    127474697  127475864  Pos4  0  +  127474697  127475864  255,0,0
chr7    127475864  127477031  Neg1  0  -  127475864  127477031  0,0,255
chr7    127477031  127478198  Neg2  0  -  127477031  127478198  0,0,255
chr7    127478198  127479365  Neg3  0  -  127478198  127479365  0,0,255
chr7    127479365  127480532  Pos5  0  +  127479365  127480532  255,0,0
chr7    127480532  127481699  Neg4  0  -  127480532  127481699  0,0,255
EOBED

my @ABED;
open my $BFH, '<', \$BED_FILE;
while (my $line = <$BFH>) {
    my ($chrom, $chrBeg, $chrEnd, $name, $score, $strand,
        $thickBeg, $thickEnd, $RGB, $blkCnt, $blkSzs, $blkBegs)
           = split /\s+/, $line;
    push @ABED, [ $chrBeg, $chrEnd, $name ];
}

####
# Compare each coordinate in my list of things to check
# against the items we stored in the BED element array
####
my $COORDS_FILE = <<EOCOORDS;
joe     123456789  123456900   Doesn't touch any of 'em
bill    127473500  127474000   Intersects Pos2 and Pos3
judy    127480550  127480560   Inside Neg4
EOCOORDS

open my $CFH, '<', \$COORDS_FILE;
while (my $line = <$CFH>) {
    my ($item, $beg, $end) = split /\s+/, $line;
    print "Checking item $item [$beg, $end] against BED file\n";
    for my $rBED (@ABED) {
        if (overlaps([$beg, $end], $rBED)) {
            print "Item $item overlaps $rBED->[2]\n";
        }
    }
    print "\n";
}

sub overlaps {
    my ($lBeg, $lEnd) = @{$_[0]};
    my ($rBeg, $rEnd) = @{$_[1]};

    return 1 if $rBeg>$lBeg and $rBeg<$lEnd;
    return 1 if $rEnd>$lBeg and $rEnd<$lEnd;
    return 1 if $lBeg>$rBeg and $lBeg<$rEnd;
    return 1 if $lEnd>$rBeg and $lEnd<$rEnd;
    return 0;
}

