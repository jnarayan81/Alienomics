#Handle the scaff/contig blast
package AOM::scaffHandler;

#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw( minmax );
use AOM::commonSubs;
use Parallel::ForkManager;
use File::Copy;
use Fcntl qw/ :flock /;

sub scaffBlaster {
my ($param, $taxid_has_parent_ref, $taxid_has_taxlevel_ref, $taxid_has_name_ref, $geneSet, $fastaSeq, $modetype)= @_;

# Blast against NR database
print "Blast against NR database\n";
if ($param->{negateGi}) {
system ("$param->{blastn_path} -task dc-megablast -query $fastaSeq -negative_gilist $param->{rejectGi} -db $param->{blastdb_path}/nt -evalue $param->{evalue} -num_threads $param->{max_processors} -max_target_seqs 1 -max_hsps 1 -qcov_hsp_perc $param->{qcovper} -outfmt '6 qseqid staxid qstart qend sallseqi sstart send evalue length frames qcovs bitscore' -out $param->{output_folder}/intermediate_files/blast/scaffblastGenome.megablast");
} 
else {
system ("$param->{blastn_path} -task dc-megablast -query $fastaSeq -db $param->{blastdb_path}/nt -evalue $param->{evalue} -num_threads $param->{max_processors} -max_target_seqs 1 -max_hsps 1 -qcov_hsp_perc $param->{qcovper} -outfmt '6 qseqid staxid qstart qend sallseqi sstart send evalue length frames qcovs bitscore' -out $param->{output_folder}/intermediate_files/blast/scaffblastGenome.megablast");
} 

# Annotate the taxonomical information
# If the blast have no hits then create a empty file
print "Annotating the taxonomical information\n";
if (-z "$param->{output_folder}/intermediate_files/blast/scaffblastGenome.megablast") {copy("$param->{output_folder}/intermediate_files/blast/scaffblastGenome.megablast","$param->{output_folder}/intermediate_files/blast/scaffblastGenome.megablast.annotated") or die "Copy failed: $!"; } 
else {
extractBlastCovAnnot($taxid_has_parent_ref, $taxid_has_taxlevel_ref, $taxid_has_name_ref, $param->{tax_list}, "$param->{output_folder}/intermediate_files/blast/scaffblastGenome.megablast", "$param->{output_folder}/intermediate_files/blast/scaffblastGenome.megablast.annotated", $param->{max_processors});
}

#Added the detail of the genes in megablast
print "Additng info in genomeNR blasthit ....\n";
AOM::commonSubs::addGeneDetail("$param->{output_folder}/intermediate_files/blast/scaffblastGenome.megablast.annotated", $geneSet, "$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.annotated.added");

# Keep all the blast hits with annotation
print "Storing all the blast hits with annotation\n";
my ($genomeHash_ref, $min, $max) = genomeHasher("$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.annotated.added");

return ($genomeHash_ref, $min, $max);
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
#rRNA Handler
sub genomeHasher {
my $genomeRefMBlast = shift;
my %genomeHash; my @allScore; my $min=0; my $max=0;
open my $fh, "<", $genomeRefMBlast or die "$genomeRefMBlast: $!";
my $csv = Text::CSV->new ({
      sep_char=>"\t",
      binary    => 1, # Allow special character. Always set this
      auto_diag => 1, # Report irregularities immediately
      });
$fh = IO::Handle->new_from_fd( $fh, 'r');
while ( not $fh->eof ) {
  my $row = $csv->getline( $fh );
  $genomeHash{$row->[0]} = $row;
  push @allScore, $row->[10]; print "$row->[10], $row->[11]------------------------------------>>>>>>";
  #warn Dumper $row; # To see the structure
}
($min, $max) = minmax @allScore;
return (\%genomeHash, $min, $max);
}


############################################################
# min max extractor
sub minMaxExtraction {
my ($genomeRefMBlast, $accession) = @_;
my @allScore; my $min=0; my $max=0;
open my $fh, "<", $genomeRefMBlast or die "$genomeRefMBlast: $!";
my $csv = Text::CSV->new ({
      sep_char=>"\t",
      binary    => 1, # Allow special character. Always set this
      auto_diag => 1, # Report irregularities immediately
      });
$fh = IO::Handle->new_from_fd( $fh, 'r');
while ( not $fh->eof ) {
  my $row = $csv->getline( $fh );
  my @val = split /\:/, $row->[10]; #The name have : seperated accesion number
  if ($val[0] eq $accession) {
	push @allScore, $row->[10];
   }
}
($min, $max) = minmax @allScore;
return ($min, $max);
}

############################################################
sub extractBlastCovAnnot {
  my ($taxid_has_parent_ref, $taxid_has_taxlevel_ref, $taxid_has_name_ref, $tax_list, $fileName, $outFile, $maxCore) = @_;
  open my $fh, '<', "$fileName"; 
  my $lCnt; my %lines; 
  while (<$fh>) {
  	chomp; 
  	$lCnt++;
  	$lines{$lCnt} = $_;
  }

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
    #next if length($blocks{$blkCntNo[$child]}) <= $length;
    my $pid = $pm->start($lCntNo[$child]) and next NAMES;
    extractTax($taxid_has_parent_ref, $taxid_has_taxlevel_ref, $taxid_has_name_ref, "$lines{$lCntNo[$child]}", $tax_list, $outFile);
    $pm->finish($child); # pass an exit code to finish
  }
  print "Waiting for all the annotation jobs to complete...\n";
  $pm->wait_all_children;
  print "Annotation DONE ... Everybody is out of the computation pool!\n";

close $fh;
}

############################################################
#Get uniq
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}


############################################################
sub extractTax {
  my ($taxid_has_parent_ref, $taxid_has_taxlevel_ref, $taxid_has_name_ref, $line, $tax_list, $outFile)=@_;
  my %taxid_has_parent=%$taxid_has_parent_ref; my %taxid_has_taxlevel = %$taxid_has_taxlevel_ref; 
  my %taxid_has_name = %$taxid_has_name_ref;
  my %contig_taxinfo; my @allNames;
  my @tax_list=split(/\,/, $tax_list);
  my %tax_levels;
  my ($nam, $id)= split(/\t/, $line);
  foreach (@tax_list) {$tax_levels{$_} = 1}
	$contig_taxinfo{$nam}= &taxonomy_report($id, \%taxid_has_taxlevel, \%taxid_has_name, \%tax_levels, \%taxid_has_parent);

	for my $tax_level (@tax_list) { push @allNames, ("\t" . (exists(${$contig_taxinfo{$nam}}{$tax_level}) ? ${$contig_taxinfo{$nam}}{$tax_level} : "NA")); }

	#return @allNames;
	foreach my $elem (@allNames) { $elem =~ s{^\s+|\s+$}{}g; $elem =~ s/\s+/_/g; $elem = lc $elem; }
	my $annot;
	if (@allNames) { $annot="$allNames[0]\t$allNames[1]\t$allNames[2]\t$allNames[3]"; } else { $annot="NA\tNA\tNA\tNA";}
	collision_free_write($outFile, "$line\t$annot"); #Write into outfile
	
sub get_parents {
    my (@all) = @_;
    my $current_id = $all[0];
    if (exists $taxid_has_parent{$current_id} and $current_id ne $taxid_has_parent{$current_id}) {
        unshift @all, $taxid_has_parent{$current_id};
        @all = &get_parents(@all);
    }
    return @all;
}

}

#Write into outfile -- collision free because of multicore usesage
sub collision_free_write {
  my($outFile, $msg) = @_;

  open my $ofh, ">>", $outFile or die  "$0 [$$]: open: $!";
  flock $ofh, LOCK_EX      or die  "$0 [$$]: flock: $!";
  print $ofh "$msg\n"      or die  "$0 [$$]: write: $!";
  close $ofh               or warn "$0 [$$]: close: $!";
}


############################################################
#Report the taxonomy using taxdmp
sub taxonomy_report {
    my ($hit_taxid, $taxid_has_taxlevel_ref, $taxid_has_name_ref, $tax_levels_ref, $taxid_has_parent_ref) = @_;
    my %taxid_has_taxlevel=%$taxid_has_taxlevel_ref; my %taxid_has_name=%$taxid_has_name_ref; my %tax_levels=%$tax_levels_ref; my %taxid_has_parent=%$taxid_has_parent_ref;
    my @parents = &get_parents($hit_taxid);
    # convert @parents to tax names:
    my %taxinfo;
    # my $taxonomy_report_string = "";
    for my $parent (@parents) {
        if (exists $taxid_has_taxlevel{$parent} and exists $tax_levels{$taxid_has_taxlevel{$parent}}) {
            $taxinfo{$taxid_has_taxlevel{$parent}} = $taxid_has_name{$parent};
        }
    }
    return \%taxinfo;
}

############################################################
#Load the nodes name from taxdump
sub load_nodes_names {
    my $fh;
    my (%taxid_has_parent, %taxid_has_taxlevel, %taxid_has_name);
    my $nodesfile = shift @_;
    my $namesfile = shift @_;
    $fh = &read_fh($nodesfile);
    while (my $line = <$fh>) {
        # line in nodes.dmp should match the regexp below.
        # Change the regexp if NCBI changes their file format
        next if $line !~ /^(\d+)\s*\|\s*(\d+)\s*\|\s*(.+?)\s*\|/;
        $taxid_has_parent{$1} = $2;
        $taxid_has_taxlevel{$1} = $3;
    }
    close $fh;
    
    $fh = &read_fh($namesfile);
    while (my $line = <$fh>) {
        next unless $line =~ /^(\d+)\s*\|\s*(.+?)\s*\|.+scientific name/;
        $taxid_has_name{$1} = $2;
    }
    return (\%taxid_has_parent, \%taxid_has_taxlevel, \%taxid_has_name);
}

############################################################
#Open and Read a file
sub read_fh {
    my $filename = shift @_;
    my $filehandle;
    if ($filename =~ /gz$/) {
        open $filehandle, "gunzip -dc $filename |" or die $!;
    }
    else {
        open $filehandle, "<$filename" or die $!;
    }
    return $filehandle;
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



1;
