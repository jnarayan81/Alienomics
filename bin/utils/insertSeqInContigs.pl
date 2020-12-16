# sub signature:
#insertSEQintoCONTIGatLOC( SEQ , CONTIG , LOC ) ;

# Fasta file to read
my $file = "$ARGV[0]";

# Sequence var to hold the fasta data
my $sequence = "";

# Open file for reading
open my $FASTAFILE, $file or die "Could not open $file: $!";

$header=<$FASTAFILE>; #Store line 1
print $header;
# Read file line-by-line
while(my $line = <$FASTAFILE>)  {
    chomp;
    # Concatenate everything from the file into a single var
    $sequence .= $line;
}

# Close opened files
close $FASTAFILE;



#$sequence = "ATATGATGATAGATGATAGTAGATAGATAGATAGATAGATAG";
$sequence =~ tr/\r\n//d;
#print $sequence;

use File::Slurp "read_file";
my $toInsert = read_file($ARGV[1]);
$toInsert =~ tr/\r\n//d;

$sequence = insertSEQintoCONTIGatLOC( "$toInsert", $sequence , 4000 ); #2600 #4000

print $sequence;

# OUTPUT:
# s = ATATGATGATAGATGATAGTAGATAGATAGATAGATAGATAG
# insert     = 'CCCC '
#     (now, s =  ATATGATGATAGATGATAGTAGATAGCCCCATAGATAGATAGATAG' )";

sub insertSEQintoCONTIGatLOC{ 
	my ( $SEQ , $CONTIG , $LOC ) = @_;
	substr( $CONTIG , $LOC , -length($CONTIG) ) = $SEQ ;
	return $CONTIG; 
}
