package AOM::refHandler;

#Handle the reference blast

#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw( minmax );
use AOM::commonSubs;

sub refBlaster {
my ($param, $geneSet)=@_;

#Create local blastDB
# Assign permission
#chmod(0660, $parameters_ref->{refDB}) or die "Couldn't chmod $parameters_ref->{refDB}: $!";
#Create local blastDB

if ($param->{create_local_db} eq 'yes') {
   if ($param->{alignment_mode} eq 'blast') { # warning: this will blast against a GENOME
      print "Creating blast DB ....";
      system ("$param->{makeblastdb_path} -in $param->{reference_genome_file} -parse_seqids -dbtype nucl -out $param->{project_dir_path}/localDB/refDB/myDB");
   }
   elsif ($param->{alignment_mode} eq 'diamond') { # warning: this will blast against a PROTEOME
      print "Creating diamond DB for proteic reference ....";
      system ("$param->{diamond_path} makedb --in $param->{reference_proteome_file} --db $param->{project_dir_path}/localDB/refDB/myDB --threads $param->{max_processors}");
   }
}
else  {print "Do database !!"; }

#Blast against localRef database [would be nice to better parallelize this one ?]
if ($param->{alignment_mode} eq 'blast') { # warning: this will blast against a GENOME
   print "Blast against reference genome ....\n";
   system ("$param->{blastn_path} -task $param->{blast_task} -query $param->{project_dir_path}/intermediate_files/bed/gene.fa -db $param->{project_dir_path}/localDB/refDB/myDB -evalue $param->{evalue} -num_threads $param->{max_processors} -max_target_seqs 1 -max_hsps 1 -qcov_hsp_perc $param->{qcovper} -outfmt '6 qseqid staxid qstart qend sallseqi sstart send evalue length frames qcovs bitscore' -out $param->{project_dir_path}/intermediate_files/ref/refgenome.megablast");
}
elsif ($param->{alignment_mode} eq 'diamond') { # warning: this will blast against a PROTEOME
   print "Blast against proteic reference ....\n";
   system ("$param->{diamond_path} $param->{diamond_task} --more-sensitive --query $param->{project_dir_path}/intermediate_files/bed/gene.fa --db $param->{project_dir_path}/localDB/refDB/myDB --max-target-seqs 1 --max-hsps 1 --evalue $param->{evalue} --query-cover $param->{qcovper} --threads $param->{max_processors} --outfmt 6 qseqid sseqid qstart qend sstart send evalue length qframe qcovhsp bitscore --out $param->{project_dir_path}/intermediate_files/ref/refgenome.megablast");
}



#Extract all no_hits_ids
print "Extracting all no_hits_ids - reference based\n";
my $noHitsValues_ref = extractNoHits("$param->{project_dir_path}/intermediate_files/bed/gene.fa", "$param->{project_dir_path}/intermediate_files/ref/refgenome.megablast", "$param->{project_dir_path}/intermediate_files/ref/refgenome.ids.noHits");

print "Extracting fasta sequences - reference based\n";
extractFasta($noHitsValues_ref, "$param->{project_dir_path}/intermediate_files/bed/gene.fa", "$param->{project_dir_path}/intermediate_files/ref/refNoHits.fa");
#extractFasta2("$param->{project_dir_path}/intermediate_files/bed/gene.fa", "$param->{project_dir_path}/intermediate_files/ref/refgenome.ids.noHits", "$param->{project_dir_path}/intermediate_files/ref/refNoHits.fa");

#Added the detail of the genes in megablast
print "Adding info in reference blasthit ....\n";
AOM::commonSubs::addGeneDetail("$param->{project_dir_path}/intermediate_files/ref/refgenome.megablast", $geneSet, "$param->{project_dir_path}/intermediate_files/ref/refgenome.megablast.added");

print "Storing hits in hash - reference based\n";
my ($refHash_ref, $min, $max)=refHasher("$param->{project_dir_path}/intermediate_files/ref/refgenome.megablast.added");
return ($refHash_ref, $min, $max);
}

#rRNA Handler
sub refHasher {
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
  push @allScore, $row->[10];
  #warn Dumper $row; # To see the structure
}
($min, $max) = minmax @allScore;
return (\%buscoHash, $min, $max);
}



############################################################
#Extract HoHits in blastDB Ids
sub extractNoHits {
my ($genefasta, $blastAnnot, $noHitfile) = @_;
my $allFastaIds_ref = extractFastaIds($genefasta);
my $allBlastIds_ref = extractBlastIds($blastAnnot);
 
# Get only NoHitsIds values
my %blastIds = map {$_=>1} @$allBlastIds_ref;
my @noHitsIds = grep { !$blastIds{$_} } @$allFastaIds_ref;

open my $nohitfh, ">", $noHitfile or die "$noHitfile: $!";
#Write in file genome.ids.noHits

foreach my $val (@noHitsIds) {print $nohitfh "$val\n";}
return \@noHitsIds;

}

############################################################
#Extract Fasta Ids
sub extractFastaIds {
my $genefasta = shift;
my @allFastaIds;

open FASTA, "<", $genefasta or die "$genefasta: $!";
while(<FASTA>) {
    chomp($_);
    if ($_ =~  m/^>/ ) {
        $_ =~ s/\>//g;
	my @ids= split /\s+/, $_;
        push @allFastaIds, $ids[0];
    }
}
close FASTA;
return \@allFastaIds;
}

############################################################
#Extract BlastIds
sub extractBlastIds {
my $blastAnnot = shift;
my @blastArr;

open BLAST, "<", $blastAnnot or die "$blastAnnot: $!";
while(<BLAST>) {
        my @tmp = split/\t/;
        push @blastArr,$tmp[0];
}
close BLAST;
my @blastArrUniq = uniq (@blastArr);
return \@blastArrUniq;
}

############################################################
#Extract Fasta sequence
sub extractFasta {
my ($ids_ref, $fasta, $out) = @_;
my %select = map {$_=>1} @$ids_ref;

open OUT, ">$out" or die;
open FASTA, "$fasta" or die;
local $/ = "\n>";  # read by FASTA record
while (<FASTA>) {
    chomp;
    my $seq = $_;
    my ($id) = $seq =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
    if (exists($select{$id})) {
        $seq =~ s/^>*.+\n//;  # remove FASTA header
        $seq =~ s/\n//g;  # remove endlines
        print OUT ">$id\n$seq\n";
    }
}
close FASTA;
close OUT;
}

############################################################
#Get uniq
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}


1;
