#Handle the genome/gene blast
package AOM::genomeHandler;

#!/usr/bin/env perl
use strict;
#use warnings; 
#-- comment to force avoid Variable "%taxidP" will not stay shared at /home/jit/Downloads/Alienomics_v1.archives-2020-06-01/bin/./AOM/genomeHandler.pm line 383.
use List::MoreUtils qw( minmax );
use AOM::commonSubs;
use Parallel::ForkManager;
use File::Copy;
use Fcntl qw/ :flock /;

sub genomeBlaster {
my ($param, $taxidP_ref, $taxidT_ref, $taxidN_ref, $geneSet, $param_ref, $taxdb, $new_txnDB)= @_;
my %ut;

# Blast against NR database
if ($param->{alignment_mode} eq 'blast') {
	if ($param->{negateGi}) {
		system ("$param->{blastn_path} -task $param->{blast_task} -query $param->{project_dir_path}/intermediate_files/bed/gene.fa -negative_gilist $param->{rejectGi} -db $param->{blastdb_path}/nt -evalue $param->{evalue} -num_threads $param->{max_processors} -max_target_seqs 1 -max_hsps 1 -qcov_hsp_perc $param->{qcovper} -outfmt '6 qseqid staxid qstart qend sstart send evalue length frames qcovs bitscore sscinames' -out $param->{project_dir_path}/intermediate_files/blast/blastGenome.megablast");
	} 
	else {
		system ("$param->{blastn_path} -task $param->{blast_task} -query $param->{project_dir_path}/intermediate_files/bed/gene.fa -db $param->{blastdb_path}/nr -evalue $param->{evalue} -num_threads $param->{max_processors} -max_target_seqs 1 -max_hsps 1 -qcov_hsp_perc $param->{qcovper} -outfmt '6 qseqid staxid qstart qend sstart send evalue length frames qcovs bitscore' -out $param->{project_dir_path}/intermediate_files/blast/blastGenome.megablast");
	} 
}

# Diamond blastx
elsif ($param->{alignment_mode} eq 'diamond') {

my $database="$param->{diamonddb_path}/uniref50.fasta.rotifera2";
	#../RefSeqDB/Refseq.complete.nonredundant.protein.fa
	# UniprotDB/uniprot_sprot.fasta | Uniref50DB/uniref50.fasta     # Uniref50 => PROBLEM TAXONOMY mais +6000 seq avec hits !!!! ===> revenir à UniProt pour le moment ??

print "Preparing diamond database\n";
print "$param->{diamond_path} makedb --in $database --db $database --threads $param->{max_processors}";
system ("$param->{diamond_path} makedb --in $database --db $database --threads $param->{max_processors}");
#system ("$param->{diamond_path} makedb --taxonmap $param->{accession2taxid_path}/prot.accession2taxid --taxonnodes $param_ref->{taxdump_path}/nodes.dmp --taxonnames $param_ref->{taxdump_path}/names.dmp --in $database --db $database --threads $param->{max_processors} > $param->{project_dir_path}/intermediate_files/blast/makediamonddb.log");

print "Running diamond blastx\n";
print "$param->{diamond_path} $param->{diamond_task} --more-sensitive --query $param->{project_dir_path}/intermediate_files/bed/gene.fa  --db $database -k 1 --max-hsps 1 --evalue $param->{evalue} --query-cover $param->{qcovper} --threads $param->{max_processors} --outfmt 6 qseqid qstart qend sstart send evalue length qframe qcovhsp bitscore sseqid --out $param->{project_dir_path}/intermediate_files/blast/blastGenome.megablast > $param->{project_dir_path}/intermediate_files/blast/diamondblastx.log";
system ("$param->{diamond_path} $param->{diamond_task} --more-sensitive --query $param->{project_dir_path}/intermediate_files/bed/gene.fa  --db $database -k 1 --max-hsps 1 --evalue $param->{evalue} --query-cover $param->{qcovper} --threads $param->{max_processors} --outfmt 6 qseqid qstart qend sstart send evalue length qframe qcovhsp bitscore sseqid --out $param->{project_dir_path}/intermediate_files/blast/blastGenome.megablast > $param->{project_dir_path}/intermediate_files/blast/diamondblastx.log");
#system ("$param->{diamond_path} $param->{diamond_task} --more-sensitive --query $param->{project_dir_path}/intermediate_files/bed/gene.fa  --db $database -k 1 --max-hsps 1 --evalue $param->{evalue} --query-cover $param->{qcovper} --threads $param->{max_processors} --outfmt 6 qseqid staxids qstart qend sstart send evalue length qframe qcovhsp bitscore sseqid --out $param->{project_dir_path}/intermediate_files/blast/blastGenome.megablast > $param->{project_dir_path}/intermediate_files/blast/diamondblastx.log");
	# vlock size default + chunks default + 30 CPUs  => ~3h
	# --block-size 10.0 --index-chunks 1 + 30 CPUs   => ~2h30
}
else { print "Did you forgot to provide alignment_mode is config file"; exit; }


# Annotate the taxonomical information
print "Annotating the diamond blast with taxonomic information\n";
# preparing and parsing the taxid hash
#my %ut;
print ("$param->{Uniref50_taxlist}");
open my $TAXLIST, '<', "$param->{Uniref50_taxlist}" or die "Could not open file $param->{Uniref50_taxlist}: $!";
while (my $line = <$TAXLIST>) { 
	if ( $line =~ m/^(\S+)\t(\S+)$/ ) {
		$ut{$1}=$2;
	}
}
# adding the taxid to blast results
print ("$param->{project_dir_path}/intermediate_files/blast/blastGenome.megablast\n");
open my $blastresults, '<:encoding(UTF-8)', "$param->{project_dir_path}/intermediate_files/blast/blastGenome.megablast" or die "Could not open file blastGenome.megablast $!";
print ("$param->{project_dir_path}/intermediate_files/blast/blastGenome.megablast.taxid\n");
open my $BLAST_OUT_TAXID, '>', "$param->{project_dir_path}/intermediate_files/blast/blastGenome.megablast.taxid" or die "Could not create blastGenome.megablast.taxid result file";
while (my $line2 = <$blastresults>) {
	if ( $line2 =~ m/\S+\t(\S+)$/ and exists $ut{$1}) { chomp $line2 ; print $BLAST_OUT_TAXID "$line2\t$ut{$1}\n"; }
	else { chomp $line2 ; print $BLAST_OUT_TAXID "$line2\tnotaxid\n"; } 
}

#Seems much faster normal mode ---
extractTax_and_write_NORMAL($taxdb, "$param->{project_dir_path}/intermediate_files/blast/blastGenome.megablast.taxid","$param->{project_dir_path}/intermediate_files/blast/blastGenome.megablast.taxid.annotated", $param->{max_processors}, $new_txnDB);

#Added the detail of the genes in megablast
print "Adding info into NR blasthit ....\n";
AOM::commonSubs::addGeneDetail("$param->{project_dir_path}/intermediate_files/blast/blastGenome.megablast.taxid.annotated", $geneSet, "$param->{project_dir_path}/intermediate_files/blast/blastGenome.megablast.taxid.annotated.added");

# Keep all the blast hits with annotation
print "Storing all the blast hits with annotation\n";
my ($genomeHash_ref, $min, $max) = genomeHasher("$param->{project_dir_path}/intermediate_files/blast/blastGenome.megablast.taxid.annotated.added");

# Extract all no_hits_ids
print "Extracting all no_hits_ids from NR blast\n";
my $noHitsValues_ref = extractNoHits("$param->{project_dir_path}/intermediate_files/ref/refNoHits.fa", "$param->{project_dir_path}/intermediate_files/blast/blastGenome.megablast.taxid.annotated", "$param->{project_dir_path}/intermediate_files/blast/blastGenome.ids.noHits");

print "Extracting fasta sequences for no_hits_ids\n";
extractFasta($noHitsValues_ref, "$param->{project_dir_path}/intermediate_files/bed/gene.fa", "$param->{project_dir_path}/intermediate_files/blast/blastNoHits.fa");

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

  # pritning for troobleshooting
 # print "genomeHasher subroutine : genomeHash(row[0]) = $genomeHash{$row->[0]}\n\tscore (row[9]) = $row->[9]\n";

  push @allScore, $row->[9]; #print "$row->[9], $row->[10], $row->[11], $row->[12] ------------------------------------>>>>>>\n";
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

#########
sub extractBlastCovAnnot_NEW {
my ($taxdb, $fileName, $outFile, $maxCore) = @_;
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
print "Workinnnnnnnnnn $child \n";
    my $pid = $pm->start($lCntNo[$child]) and next NAMES;
     extractTax_NEW($taxdb, "$lines{$lCntNo[$child]}", $outFile);
    $pm->finish($child); # pass an exit code to finish
  }
  print "Waiting for all the annotation jobs to complete...\n Hold on... assignment takes li`l more time\n";
  $pm->wait_all_children;
  print "Annotation DONE ... Everybody is out of the computation pool!\n\n\n";

close $fh;
}

#############@
sub extractTax_NEW {
  my ($taxdb, $line, $outFile)=@_;
  my ($nam, $id)= split(/\t/, $line);
print "$line, --- $nam, $id <<<<<<\n";
# get a taxon
my ($taxon, $lName)=('NA','NA');
$taxon = $taxdb->get_taxon(-taxonid => $id);
$lName=$taxon->scientific_name; #return name
my $annot="$lName\tNA\tNA\tNA\n";
#print $annot;
collision_free_write($outFile, "$line\t$annot"); #Write into outfile
}

#############@
sub extractTax_and_write_NORMAL {
  my ($taxdb, $infile, $outFile, $corwa, $new_txnDB)=@_;
  open my $fh, '<', "$infile";
  open(my $ofh, '>', "$outFile");
  my $lCnt; my %lines;
  while (<$fh>) {
        chomp;
	my @tmpLine= split(/\t/, $_);
        $lCnt++;
	# get a taxon
	my $id = $tmpLine[11];
	my $lName ='NA';
	my $taxon = $taxdb->get_taxon(-taxonid => $id);
	if ($taxon) {
		$lName=$taxon->scientific_name;  				 # return name
		#print "General blast extractTax\t$taxon\t$lName\t$id\n\t(".join(" ",@tmpLine).")\n";				# print a correctly detected taxon
		print "General blast extractTax\t$taxon\t$lName\t$id\n";							# simpler print of a correctly detected taxon
	} 	
	else {
		#$lName=$new_txnDB->{$id};
		print "General blast extractTax\tWARNING\t$id\tis not a taxonid?\n\t(".join("\t",@tmpLine).")\n";	# print a problematic taxon id and its corresponding line
	}

	#I assume it return 'NA' if not found -- i converted it into lc
	if (lc($lName) eq "na") { $lName=$new_txnDB->{$id};}

	#print $new_txnDB->{$id}; 	#exit;

	my $annot="$_\t$lName\tNA\tNA\tNA\n";
	print "\t\t$id\t=>\t$lName\t=>\t$new_txnDB->{$id}\n";
	print $ofh $annot;
        #$lines{$lCnt} = $_;
  }
  close $ofh;
}

### example of annotated diamond output:
# MTYJ01001406.1:1-964	964	566	142	274	1.1e-69	133	-1	41.4	270.8	UniRef50_A0A1W0WIN8	232323

### exemple of how the new script is missing something: 
# General blast extractTax	WARNING	1977087	is not a taxonid?
#	(scaf_180:1566-2986	438	1403	67	394	1.8e-76	329	3	68.0	293.9	UniRef50_A0A2G6CNT0	1977087)
#	1977087	=>	NA	=>	Proteobacteria_bacterium


############################################################
sub extractBlastCovAnnot {
my ($taxidP_ref, $taxidT_ref, $taxidN_ref, $tax_list, $fileName, $outFile, $maxCore) = @_;
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
    extractTax($taxidP_ref, $taxidT_ref, $taxidN_ref, "$lines{$lCntNo[$child]}", $tax_list, $outFile);
    $pm->finish($child); # pass an exit code to finish
  }
  print "Waiting for all the annotation jobs to complete...\n Hold on... assignment takes li`l more time\n";
  $pm->wait_all_children;
  print "Annotation DONE ... Everybody is out of the computation pool!\n\n\n";

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
my ($taxidP_ref, $taxidT_ref, $taxidN_ref, $line, $tax_list, $outFile)=@_;
my %taxidP=%$taxidP_ref; my %taxidT = %$taxidT_ref; 
my %taxidN = %$taxidN_ref;
my %contig_taxinfo; my @allNames;
my @tax_list=split(/\,/, $tax_list);
my %tax_levels;
my ($nam, $id)= split(/\t/, $line);
foreach (@tax_list) {$tax_levels{$_} = 1}
	$contig_taxinfo{$nam}= &taxonomy_report($id, \%taxidT, \%taxidN, \%tax_levels, \%taxidP);

	for my $tax_level (@tax_list) { push @allNames, ("\t" . (exists(${$contig_taxinfo{$nam}}{$tax_level}) ? ${$contig_taxinfo{$nam}}{$tax_level} : "NA")); }

	#return @allNames;
	foreach my $elem (@allNames) { $elem =~ s{^\s+|\s+$}{}g; $elem =~ s/\s+/_/g; $elem = lc $elem; }
	my $annot;
	#if (@allNames) { $annot="$allNames[0]\tNA\tNA\tNA"; } else { $annot="NA\tNA\tNA\tNA";}
	if (@allNames) { $annot="$allNames[0]"; } else { $annot="NA";}
	collision_free_write($outFile, "$line\t$annot"); #Write into outfile

sub get_parents {
    my (@all) = @_;
    my $current_id = $all[0];
    if (exists $taxidP{$current_id} and $current_id ne $taxidP{$current_id}) {
        unshift @all, $taxidP{$current_id};
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
    my ($hit_taxid, $taxidT_ref, $taxidN_ref, $tax_levels_ref, $taxidP_ref) = @_;
    my %taxidT=%$taxidT_ref; my %taxidN=%$taxidN_ref; my %tax_levels=%$tax_levels_ref; my %taxidP=%$taxidP_ref;
    my @parents = &get_parents($hit_taxid);
    # convert @parents to tax names:
    my %taxinfo;
    # my $taxonomy_report_string = "";
    for my $parent (@parents) {
        if (exists $taxidT{$parent} and exists $tax_levels{$taxidT{$parent}}) {
            $taxinfo{$taxidT{$parent}} = $taxidN{$parent};
        }
    }
    return \%taxinfo;
}

############################################################
#Load the nodes name from taxdump
sub load_nodes_names {
    my $fh;
    my (%taxidP, %taxidT, %taxidN);
    my $nodesfile = shift @_;
    my $namesfile = shift @_;
    $fh = &read_fh($nodesfile);
    while (my $line = <$fh>) {
        # line in nodes.dmp should match the regexp below.
        # Change the regexp if NCBI changes their file format
        next if $line !~ /^(\d+)\s*\|\s*(\d+)\s*\|\s*(.+?)\s*\|/;
        $taxidP{$1} = $2;
        $taxidT{$1} = $3;
    }
    close $fh;
    
    $fh = &read_fh($namesfile);
    while (my $line = <$fh>) {
        next unless $line =~ /^(\d+)\s*\|\s*(.+?)\s*\|.+scientific name/;
        $taxidN{$1} = $2;
    }
    return (\%taxidP, \%taxidT, \%taxidN);
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
