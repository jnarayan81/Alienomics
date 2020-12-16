package AOM::buscoHandler;

#Handle the busco blast

#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw( minmax );
use AOM::commonSubs;

sub buscoBlaster {
my ($param, $geneSet)=@_;

#Create local blastDB
print "Creating blast buscoDB";
if ($param->{alignment_mode} eq 'blast'){
	system ("$param->{makeblastdb_path} -in $param->{busco_file} -parse_seqids -dbtype prot -out $param->{project_dir_path}/localDB/buscoDB/mybuscoDB");
	}
elsif ($param->{alignment_mode} eq 'diamond'){
        print "$param->{diamond_path} makedb --in $param->{busco_file} --db $param->{project_dir_path}/localDB/buscoDB/mybuscoDB";
        system ("$param->{diamond_path} makedb --in $param->{busco_file} --db $param->{project_dir_path}/localDB/buscoDB/mybuscoDB");
	}
else { print "Did you forgot to provide alignment_mode\n"; exit;}


print "diamond blastx against buscoDB ....\n";
if ($param->{alignment_mode} eq 'blast'){
	system ("$param->{blastx_path} -task blastx-fast -query $param->{project_dir_path}/intermediate_files/bed/gene.fa -db $param->{project_dir_path}/localDB/buscoDB/mybuscoDB -evalue $param->{evalue} -num_threads $param->{max_processors} -max_target_seqs 1 -max_hsps 1 -qcov_hsp_perc $param->{qcovper} -outfmt '6 qseqid qstart qend sstart send evalue length frames qcovs bitscore' -out $param->{project_dir_path}/intermediate_files/busco/buscoGenes.megablast");
}
elsif ($param->{alignment_mode} eq 'diamond'){
        print "$param->{diamond_path} $param->{diamond_task} --more-sensitive --query $param->{project_dir_path}/intermediate_files/bed/gene.fa --db $param->{project_dir_path}/localDB/buscoDB/mybuscoDB -k 1 --max-hsps 1 --evalue $param->{evalue} --threads $param->{max_processors} --query-cover $param->{qcovper} --outfmt 6 qseqid  qstart qend sstart send evalue length qframe qcovhsp bitscore --out $param->{project_dir_path}/intermediate_files/busco/buscoGenes.megablast";
        system ("$param->{diamond_path} $param->{diamond_task} --more-sensitive --query $param->{project_dir_path}/intermediate_files/bed/gene.fa --db $param->{project_dir_path}/localDB/buscoDB/mybuscoDB -k 1 --max-hsps 1 --evalue $param->{evalue} --threads $param->{max_processors} --query-cover $param->{qcovper} --outfmt 6 qseqid qstart qend sstart send evalue length qframe qcovhsp bitscore --out $param->{project_dir_path}/intermediate_files/busco/buscoGenes.megablast");

}
else { print "Did you forgot to provide alignment_mode\n"; exit;}

# warning : there is not equivalent to the blast 'frames' field in diamond output, there is only 'qframe' and 'sframe' fields
# warning : qcovs becomes qcovhsp in diamond
# see details here : http://www.metagenomics.wiki/tools/blast/blastn-output-format-6 and here : diamond --help

# --taxonmap
# --taxonnodes (taxonomy nodes.dmp from NCBI)
# --taxonnames (taxonomy names.dmp from NCBI)
# ===> useless for busco !


#Added the detail of the genes in megablast
print "Adding info into buscoDB blasthit results....\n";
AOM::commonSubs::addGeneDetail("$param->{project_dir_path}/intermediate_files/busco/buscoGenes.megablast", $geneSet, "$param->{project_dir_path}/intermediate_files/busco/buscoGenes.megablast.added");

my ($buscoHash_ref, $min, $max)=buscoHasher("$param->{project_dir_path}/intermediate_files/busco/buscoGenes.megablast.added");
return ($buscoHash_ref, $min, $max);
}

#rRNA Handler
sub buscoHasher {
my $buscoRefMBlast = shift;
my %buscoHash; my @allScore; my $min=0; my $max=0;
open my $fh, "<", $buscoRefMBlast or die "$buscoRefMBlast: $!";
my $csv = Text::CSV->new ({
      sep_char=>"\t",
      binary    => 1, # Allow special character. Always set this
      auto_diag => 1, # Report irregularities immediately
      });
$fh = IO::Handle->new_from_fd( $fh, 'r');
while ( not $fh->eof ) {
  my $row = $csv->getline( $fh );
  $buscoHash{$row->[0]} = $row;
  push @allScore, $row->[10]; #print "$row->[10] ----<<<<<<<<<\n\n";
  #warn Dumper $row; # To see the structure
}
($min, $max) = minmax @allScore;
return (\%buscoHash, $min, $max);
}

1;
