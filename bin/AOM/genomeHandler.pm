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
my ($param, $taxidP_ref, $taxidT_ref, $taxidN_ref, $geneSet, $param_ref, $taxdb, $new_txnDB, $alienomics_location)= @_;
my %ut;

# BLASTN against NR database (highly deprecated)
if ($param->{alignment_mode} eq 'blast') {
	if ($param->{negateGi}) {
		system ("$param->{blastn_path} -task $param->{blast_task} -query $param->{output_folder}/intermediate_files/bed/gene.fa -negative_gilist $param->{rejectGi} -db $param->{blastdb_path}/nt -evalue $param->{evalue} -num_threads $param->{max_processors} -max_target_seqs 1 -max_hsps 1 -qcov_hsp_perc $param->{qcovper} -outfmt '6 qseqid staxid qstart qend sstart send evalue length frames qcovs bitscore sscinames' -out $param->{output_folder}/intermediate_files/blast/blastGenome.megablast");
	} 
	else {
		system ("$param->{blastn_path} -task $param->{blast_task} -query $param->{output_folder}/intermediate_files/bed/gene.fa -db $param->{blastdb_path}/nr -evalue $param->{evalue} -num_threads $param->{max_processors} -max_target_seqs 1 -max_hsps 1 -qcov_hsp_perc $param->{qcovper} -outfmt '6 qseqid staxid qstart qend sstart send evalue length frames qcovs bitscore' -out $param->{output_folder}/intermediate_files/blast/blastGenome.megablast");
	} 
}

# DIAMOND BLASTX (default parameter)
elsif ($param->{alignment_mode} eq 'diamond') {


#print "Removing user-defined list of taxa from Uniref50 database ('excludeTaxID' parameter in config file)\n";
my $uniref50="$param->{uniref50_path}/uniref50.fasta";
my $uniref50oneline="$param->{uniref50_path}/uniref50.fasta.oneline";
my $taxlist="$param->{uniref50_path}/uniref50.taxlist";
my $excludeTaxID="$param->{excludeTaxID}";
my $taxlist_out="$param->{uniref50_path}/uniref50.taxlist.taxid_ignored";
my $database="$param->{uniref50_path}/uniref50.fasta.taxid_ignored";

print "parsing taxID to ignore:\n";
my @exclude_fields = split /\|/, $excludeTaxID;
my %exclude = ();
foreach ( @exclude_fields ) {
        $exclude{$_}=1 ;
        print "\t$_\n"
}

print "transform fasta sequence to oneliners\n";
open(FASTAIN, '<', $uniref50) or die $!;
open(FASTAOUT, '>', $uniref50oneline) or die $!;
while (<FASTAIN>) { $. > 1 and /^>/ ? print FASTAOUT "\n" : chomp; print FASTAOUT }
close(FASTAIN); close (FASTAOUT);

print "filtering taxlist from taxIDs\n";
open(TAXLISTIN, '<', $taxlist) or die $!;
open(TAXLISTOUT, '>', $taxlist_out) or die $!;
my %remaining = (); my $SeqID;
while (my $row = <TAXLISTIN>) {
        #print "\n\nNEW ROW\n";
        #print "$row";
        chomp($row);
        #print "$row\n";
        my @fields = split / /, $row;
        #print @fields;
        my ($SeqID ,$TaxID) = @fields[0, 1];
        #print "SeqID = $SeqID\tTaxID = $TaxID\n";
        #print "checking TaxID in exclusion list: $TaxID\t$exclude{$TaxID}";
        unless ( exists $exclude{$TaxID}) {
                print TAXLISTOUT "$row\n";
                $remaining{$SeqID} = $TaxID;
                #print "\t=> kept\n";
        }
        #else { print "removing: $TaxID\n";}
}
close (TAXLISTIN); close(TAXLISTOUT);

print "writing filtered uniref50 fasta file\n";
open(FASTAONELINE, '<', $uniref50oneline) or die $!;
open(UNIREFOUT, '>', $database) or die $!;
my $name; my $seq;
while (my $row = <FASTAONELINE>) {
        chomp($row);
        if ($row =~ m/>(\S+).*/) {
                $SeqID=$1;
                $name=$row;
        }
        else {
                chomp($row);
                $seq=$row;
                if ( exists $remaining{$SeqID} ) {
                        print UNIREFOUT "$name\n$seq\n";
                }
        }
}

system("rm -f $uniref50oneline");


print "Preparing diamond database\n";
AOM::commonSubs::custom_system ("$param->{diamond_path} makedb --in $database --db $database --threads $param->{max_processors} > /dev/null 2>&1");

print "Running diamond blastx\n";
# SINGLE HIT VERSION
#print "$param->{diamond_path} $param->{diamond_task} --more-sensitive --query $param->{output_folder}/intermediate_files/bed/gene.fa  --db $database -k 1 --max-hsps 1 --evalue $param->{evalue} --query-cover $param->{qcovper} --threads $param->{max_processors} --outfmt 6 qseqid qstart qend sstart send evalue length qframe qcovhsp bitscore sseqid --out $param->{output_folder}/intermediate_files/blast/blastGenome.megablast > $param->{output_folder}/intermediate_files/blast/diamondblastx.log";
#system ("$param->{diamond_path} $param->{diamond_task} --more-sensitive --query $param->{output_folder}/intermediate_files/bed/gene.fa  --db $database -k 1 --max-hsps 1 --evalue $param->{evalue} --query-cover $param->{qcovper} --threads $param->{max_processors} --outfmt 6 qseqid qstart qend sstart send evalue length qframe qcovhsp bitscore sseqid --out $param->{output_folder}/intermediate_files/blast/blastGenome.megablast > $param->{output_folder}/intermediate_files/blast/diamondblastx.log > /dev/null 2>&1");
#system ("$param->{diamond_path} $param->{diamond_task} --more-sensitive --query $param->{output_folder}/intermediate_files/bed/gene.fa  --db $database -k 1 --max-hsps 1 --evalue $param->{evalue} --query-cover $param->{qcovper} --threads $param->{max_processors} --outfmt 6 qseqid staxids qstart qend sstart send evalue length qframe qcovhsp bitscore sseqid --out $param->{output_folder}/intermediate_files/blast/blastGenome.megablast > $param->{output_folder}/intermediate_files/blast/diamondblastx.log");

# MULTI-HITS VERSION (keeping multi-hits + adding pident field !):
AOM::commonSubs::custom_system  ("$param->{diamond_path} $param->{diamond_task} --more-sensitive --query $param->{output_folder}/intermediate_files/bed/gene.fa --db $database -k 50 --max-hsps 1 --evalue $param->{evalue} --query-cover $param->{qcovper} --threads $param->{max_processors} --outfmt 6 qseqid qstart qend sstart send evalue length qframe qcovhsp bitscore pident sseqid --out $param->{output_folder}/intermediate_files/blast/blastGenome.megablast > $param->{output_folder}/intermediate_files/blast/diamondblastx.log > /dev/null 2>&1");
# evalue = 6 ; length = 7 ; bitscore = 10 ; pident = 11

}

#else { 
#	print "Did you forgot to provide alignment_mode in config file ? (default is alignment_mode=diamond)";
#	exit;
#}


# IMPROVMENTS MULTI-HITS (filtering hits)
# filters = evalue < 1e-05 ; length > 60 aa ; 30<pident<90
open BLASTRAW, '<', "$param->{output_folder}/intermediate_files/blast/blastGenome.megablast" or die $!;
open BLASTFILTERED, '>', "$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.filtered" or die $!;
while (my $line = <BLASTRAW>) {
	chomp($line); 
	my @fields = split /\t/, $line;
	if ( ("$fields[5]" <= 1e-05) && ("$fields[6]" >= 60) && (30 <= "$fields[10]") && ( "$fields[10]" < 90) ) {
		print BLASTFILTERED "$line\n";
	}
}
close(BLASTRAW);
close(BLASTFILTERED);



# LOADING TAXONOMIC INFORMATION
print "Annotating blast results with taxonomic information\n";
open my $TAXLIST, '<', "$param->{uniref50_path}/uniref50.taxlist" or die "Could not open file $param->{uniref50_path}/uniref50.taxlist: $!";
while (my $line = <$TAXLIST>) { 
	if ( $line =~ m/^(\S+)\s(\S+)$/ ) {
		$ut{$1}=$2;
		#print "taxlist dictionnary:\t$1\t$ut{$1}\n";
	}
}
# SINGLE HIT VERSION
#print ("$param->{output_folder}/intermediate_files/blast/blastGenome.megablast\n");
#open my $blastresults, '<:encoding(UTF-8)', "$param->{output_folder}/intermediate_files/blast/blastGenome.megablast" or die "Could not open file blastGenome.megablast $!";
#print ("$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.taxid\n");
#open my $BLAST_OUT_TAXID, '>', "$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.taxid" or die "Could not create blastGenome.megablast.taxid result file";
#while (my $line2 = <$blastresults>) {
#	if ( $line2 =~ m/\S+\t(\S+)$/ and exists $ut{$1}) { chomp $line2 ; print $BLAST_OUT_TAXID "$line2\t$ut{$1}\n"; }
#	else { chomp $line2 ; print "no taxonid for: $1\n"; print $BLAST_OUT_TAXID "$line2\tnotaxid\n"; } 
#}

# IMPROVMENTS MULTI-HITS
open my $blastresults, '<:encoding(UTF-8)', "$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.filtered" or die "Could not open file blastGenome.megablast $!";
open my $BLAST_OUT_TAXID, '>', "$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.filtered.taxid" or die "Could not create blastGenome.megablast.taxid result file";
while (my $line2 = <$blastresults>) {
        if ( $line2 =~ m/\S+\t(\S+)$/ and exists $ut{$1}) { chomp $line2 ; print $BLAST_OUT_TAXID "$line2\t$ut{$1}\n"; }
        else { chomp $line2 ; print "no taxonid for: $1\n"; print $BLAST_OUT_TAXID "$line2\tnotaxid\n"; }
}


print "Merging taxonomy info and diamond blast output\n";
# SINGLE HIT VERSION
#extractTax_and_write_NORMAL($taxdb, "$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.taxid","$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.taxid.annotated", $param->{max_processors}, $new_txnDB);
# IMPROVMENTS MULTI-HITS
extractTax_and_write_NORMAL($taxdb, "$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.filtered.taxid","$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.filtered.taxid.annotated", $param->{max_processors}, $new_txnDB);


print "Adding details to annotated blast results\n";
# SINGLE HIT VERSION
#AOM::commonSubs::addGeneDetail("$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.taxid.annotated", $geneSet, "$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.taxid.annotated.added");
# IMPROVMENTS MULTI-HITS
AOM::commonSubs::addGeneDetail("$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.filtered.taxid.annotated", $geneSet, "$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.filtered.taxid.annotated.added");




# IMPROVMENTS MULTI-HITS
# simplifying blast output by keeping a single row: the best representative of the majoritary origin (self vs alien)
print "Detecting best representative hits\n";
my %index = ();
open(BLASTIN, '<', "$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.filtered.taxid.annotated.added") or die $!;
open(BESTBLAST, '>',  "$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.filtered.taxid.annotated.added.best") or die $!;
open (TAXDETECTION, '>', "$param_ref->{output_folder}/intermediate_files/blast/taxonomy_detection_uniref50.log") or die $!;
foreach my $line (<BLASTIN>) {
        chomp($line);
	print TAXDETECTION "\n$line\n";
        my @fields = split(/\t/, $line);
        my $query = $fields[0];
        my $evalue = $fields[5];
        my $bitscore = $fields[9];
        my $pident = $fields[10];
        my $target = $fields[11];
        my $target_taxid = $fields[12];
        my $target_taxname = $fields[13];
        $index{$query}{evalue}=$evalue;
        $index{$query}{bitscore}=$bitscore;
        $index{$query}{pident}=$pident;
        $index{$query}{target}=$target;
        $index{$query}{target_taxid}=$target_taxid;
        $index{$query}{target_taxname}=$target_taxname;

	print TAXDETECTION "investigating $query\n";
	print TAXDETECTION "\tbitscore=$index{$query}{bitscore}\ttarget id=$index{$query}{target_taxid}\n";

        # get taxonomy decision (self vs alien) on this line
        my $taxdecision = parseTaxMultiHits($target_taxid, "$param_ref->{clade_self}");
	#print TAXDETECTION "\t=> $taxdecision\n"	

	if ( $taxdecision eq "self" ) { 
		# if not first self encounter + better self bitscore that previously registered
		if ( (exists $index{$query}{best_self_bitscore}) && ($index{$query}{bitscore} > $index{$query}{best_self_bitscore}) ) {
			$index{$query}{nb_self}=$index{$query}{nb_self} + 1;
                        $index{$query}{best_self_hit}=$line;
                        $index{$query}{best_self_bitscore}=$index{$query}{bitscore};
                        print TAXDETECTION "\tencountering $query (on $target: self) for the $index{$query}{nb_self} time => new best self hit\n";
		}
		# if not first self encounter + not a better hit
		elsif ( (exists $index{$query}{best_self_bitscore}) && ($index{$query}{bitscore} <= $index{$query}{best_self_bitscore}) ) {
			$index{$query}{nb_self}=$index{$query}{nb_self} + 1;
			print TAXDETECTION "\tencountering $query (on $target: self) for the $index{$query}{nb_self} time\n";
		}
		# if first self encounter
		elsif ( ! exists $index{$query}{best_self_bitscore}) {
			$index{$query}{nb_self}=1;
			$index{$query}{best_self_hit}=$line;
                        $index{$query}{best_self_bitscore}=$index{$query}{bitscore};
                        print TAXDETECTION "\tfirst self encounter for $query (on $target: self)\n";
                        print TAXDETECTION "\t$index{$query}{nb_self}";
		}
		# if issue
		else { print TAXDETECTION "\tISSUE: with detecting encounters !!!\n$line\n";}
	}
        elsif ( $taxdecision eq "alien" ) {
		# if not first alien encounter + better alien bitscore that previously registered
                if ( (exists $index{$query}{best_alien_bitscore}) && ($index{$query}{bitscore} > $index{$query}{best_alien_bitscore}) ) {
                        $index{$query}{nb_alien}=$index{$query}{nb_alien} + 1;
                        $index{$query}{best_alien_hit}=$line;
                        $index{$query}{best_alien_bitscore}=$index{$query}{bitscore};
                        print TAXDETECTION "\tencountering $query (on $target: alien) for the $index{$query}{nb_alien} time => new best alien hit\n";
                }
                # if not first alien encounter + not a better hit
                elsif ( (exists $index{$query}{best_alien_bitscore}) && ($index{$query}{bitscore} <= $index{$query}{best_alien_bitscore}) ) {
                        $index{$query}{nb_alien}=$index{$query}{nb_alien} + 1;
                        print TAXDETECTION "\tencountering $query (on $target: alien) for the $index{$query}{nb_alien} time\n";
                }
                # if first alien encounter
                elsif ( ! exists $index{$query}{best_alien_bitscore}) {
                        $index{$query}{nb_alien}=1;
                        $index{$query}{best_alien_hit}=$line;
                        $index{$query}{best_alien_bitscore}=$index{$query}{bitscore};
                        print TAXDETECTION "\tfirst alien encounter for $query (on $target: alien)\n";
                        print TAXDETECTION "\t$index{$query}{nb_alien}";
                }
                # if issue
                else { print TAXDETECTION "\tISSUE: with detecting encounters !!!\n";}
        }
	else {  print TAXDETECTION "\tISSUE: with taxdeceision (taxid not recognized = $target_taxid)\n";}
}

foreach my $query (keys %index) {
	print TAXDETECTION "\nDetecting best representative hit for $query\n";

	if ( (exists $index{$query}{nb_alien}) || (exists  $index{$query}{nb_self}) ) { 

		my $alien_ratio= $index{$query}{nb_alien} / ($index{$query}{nb_alien} + $index{$query}{nb_self});
		print TAXDETECTION "\tnb_self = $index{$query}{nb_self}\tnb_alien = $index{$query}{nb_alien}\talien_ratio = $alien_ratio";


		if ( $alien_ratio > 0.5) { 
			print BESTBLAST $index{$query}{best_alien_hit}."\n";
			print TAXDETECTION "\t=> ALIEN\n";
		}
		if ( $alien_ratio <= 0.5) { 
			print BESTBLAST $index{$query}{best_self_hit}."\n";
			print TAXDETECTION "\t=> SELF\n";
		}
	}
	else {print TAXDETECTION "\t=> UNKNOWN (taxonomy issue ?)\n";}
}

close(BLASTIN);
close(BESTBLAST);


############################################################
sub parseTaxMultiHits {
        #print "\n\t(Starting parseTaxMultiHits subroutine)\n";
        my ($sps,$groupName)= @_;
        my $nam = $taxdb->get_taxon(-taxonid => $sps);
                if(!$nam) {
                print TAXDETECTION "\tNo taxon name for: $sps\n";
                return "unknown";
        }
        elsif($nam eq "NA") {
                print  "\t$sps\tNA\n";
        }
        else {
                print TAXDETECTION "\tsp=$sps\n";
        }
        my $tree_functions = Bio::Tree::Tree->new();
        my $lineage = lc ($tree_functions->get_lineage_string($nam));

        #print TAXDETECTION "\tsp = $sps\n\tnam = $nam\n\tlineage = $lineage\n\tgroupname = ($groupName)\n";

        if (index(lc($lineage), lc($groupName)) != -1) { return "self"; }
        elsif (index(lc($lineage), lc($groupName)) == -1) { return "alien"; }
        else {print TAXDETECTION "\t$sps\tunclassified taxa\n";}       # for troobleshooting
}
########################################################################################################


print "Storing all the blast hits with annotation\n";
# SINGLE HIT VERSION
#my ($genomeHash_ref, $min, $max) = genomeHasher("$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.taxid.annotated.added");
# IMPROVMENTS MULTI-HITS
my ($genomeHash_ref, $min, $max) = genomeHasher("$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.filtered.taxid.annotated.added.best");


print "Extracting all genes with no hits\n";
# SINGLE HIT VERSION
#my $noHitsValues_ref = extractNoHits("$param->{output_folder}/intermediate_files/ref/refNoHits.fa", "$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.taxid.annotated", "$param->{output_folder}/intermediate_files/blast/blastGenome.ids.noHits");
#extractFasta($noHitsValues_ref, "$param->{output_folder}/intermediate_files/bed/gene.fa", "$param->{output_folder}/intermediate_files/blast/blastNoHits.fa");
# IMPROVMENTS MULTI-HITS
my $noHitsValues_ref = extractNoHits("$param->{output_folder}/intermediate_files/ref/refNoHits.fa", "$param->{output_folder}/intermediate_files/blast/blastGenome.megablast.filtered.taxid.annotated", "$param->{output_folder}/intermediate_files/blast/blastGenome.ids.noHits");
extractFasta($noHitsValues_ref, "$param->{output_folder}/intermediate_files/bed/gene.fa", "$param->{output_folder}/intermediate_files/blast/blastNoHits.fa");


return ($genomeHash_ref, $min, $max);

}





###########################################################
### SUBROUTINES
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

###############################################################
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

############################################################
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

############################################################
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
	# SINGLE HIT VERSION
	#my $id = $tmpLine[11];
	# IMPROVMENTS MULTI-HITS
	my $id = $tmpLine[12];
	my $lName ='NA';
	#print("General blast extractTax\tlooking for $id\n");
	my $taxon = $taxdb->get_taxon(-taxonid => $id);
	#print("\tfound: $taxon\n");
	if ($taxon) {
		$lName=$taxon->scientific_name;  				 # return name
		#print "\t$taxon\t$lName\t$id\n";	# simpler print of a correctly detected taxon
	} 	
	else {
		#$lName=$new_txnDB->{$id};
		print "\tWARNING\t$id\tis not a taxonid?\n\t(".join("\t",@tmpLine).")\n";	# print a problematic taxon id and its corresponding line
	}

	#I assume it return 'NA' if not found -- i converted it into lc
	if (lc($lName) eq "NA") { $lName=$new_txnDB->{$id}; 
		print("\tScientific name not found for $id\n");
	}

	#print $new_txnDB->{$id}; 	#exit;

	my $annot="$_\t$lName\tNA\tNA\tNA\n";
	#print "\t\t$id\t=>\t$lName\t=>\t$new_txnDB->{$id}\n";
	print $ofh $annot;
        #$lines{$lCnt} = $_;
  }
  close $ofh;
}


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
#sub extractTax {
#my ($taxidP_ref, $taxidT_ref, $taxidN_ref, $line, $tax_list, $outFile)=@_;
#my %taxidP=%$taxidP_ref; my %taxidT = %$taxidT_ref; 
#my %taxidN = %$taxidN_ref;
#my %contig_taxinfo; my @allNames;
#my @tax_list=split(/\,/, $tax_list);
#my %tax_levels;
#my ($nam, $id)= split(/\t/, $line);
#foreach (@tax_list) {$tax_levels{$_} = 1}
#	$contig_taxinfo{$nam}= &taxonomy_report($id, \%taxidT, \%taxidN, \%tax_levels, \%taxidP);

#	for my $tax_level (@tax_list) { push @allNames, ("\t" . (exists(${$contig_taxinfo{$nam}}{$tax_level}) ? ${$contig_taxinfo{$nam}}{$tax_level} : "NA")); }

	#return @allNames;
#	foreach my $elem (@allNames) { $elem =~ s{^\s+|\s+$}{}g; $elem =~ s/\s+/_/g; $elem = lc $elem; }
#	my $annot;
	#if (@allNames) { $annot="$allNames[0]\tNA\tNA\tNA"; } else { $annot="NA\tNA\tNA\tNA";}
#	if (@allNames) { $annot="$allNames[0]"; } else { $annot="NA";}
#	collision_free_write($outFile, "$line\t$annot"); #Write into outfile

#sub get_parents {
#    my (@all) = @_;
#    my $current_id = $all[0];
#    if (exists $taxidP{$current_id} and $current_id ne $taxidP{$current_id}) {
#        unshift @all, $taxidP{$current_id};
#        @all = &get_parents(@all);
#    }
#    return @all;
#}

#}

############################################################
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
