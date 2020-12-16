use strict;
use warnings;
use Bio::DB::GenBank;
use Bio::SeqIO;

my $usage = "\nperl get_seq <file containning genbank IDs, one per line>\n\n";

if (!$ARGV[0]) {
  die ($usage);
}

my $infile = $ARGV[0];

my @ids;

my $i = 0;

open (IN, "$infile");

while (my $line = <IN>) {
  chomp $line;
  $ids[$i] = $line;
  $i++;
}

close IN;

my $db_obj = Bio::DB::GenBank->new;

foreach my $id (@ids) {
  print ("Downloading $id\n");
  my $seq_obj = $db_obj->get_Seq_by_id($id);
  my $seqio_obj = Bio::SeqIO->new(-file => ">$id.gb", -format => 'genbank' );
  $seqio_obj->write_seq($seq_obj);
}
