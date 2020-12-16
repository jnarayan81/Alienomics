package AOM::gffHandler;

# Handle the gff file

#!/usr/bin/env perl
use strict;
#use warnings;
use Data::Dumper;
#use Exporter qw(import);
#our @EXPORT_OK=qw(gff2genemodel);

# re-implementacion de gff2genemodel #########################################################
#sub element{
#	my $(start,end,strand,phase,score,attributes,type,method)=shift;
#	my $e->{start}=$start;
#	my $e->{end}=$end;
#	my $e->{strand}=$strand;
#	my $e->{phase}=$phase;
#	my $e->{score}=$score;
#	my $e->{attributes}=$attributes;
#	my $e->{type}=$type;
#	my $e->{method}=$method;
#
#	
#}

# Title		: gff2genemodel
# Usage		: my $g=gff2genemodel("path/to/gene_model.gff")
#
# Function	: Generate a hash reference from a gff3 file.
# Returns	: Hash reference
# Args		: a string containing the path to gff3 file
# Note		: only gene, mRNA, CDS, exon, intron, start_codon and stop_codon features supported.

# now the script loads all nucleotide sequence files to a hash structure,
# checks their validity and translates them to protein sequence
sub gff2gene{
	my $ingff=shift;
	my $param=shift;

	open GFF, "<$ingff" or die "Cannot open $ingff";
	#while (my $line=<GFF>) {chomp $line; print "$line\n";  exit;}

	my %gene=();
	my %mrna=();
	my %exon=();
	my %intron=();
	my %cds=();
	my %start_codon=();
	my %stop_codon=();
#	my %transposon

	my %count_type=(gene=>0,mRNA=>0,transcript=>0,exon=>0,intron=>0,CDS=>0,start_codon=>0,stop_codon=>0);

	my $dup=0;

	# gff
	#Adineta_ricciae	AUGUSTUS	gene	2	1810	0.28	+	.	g1
	my @lines=<GFF>;	
	foreach my $line(@lines){
		next if(($line=~/^#/g)||(length($line) == 1) || ($line =~ /^\s*$/));
		chomp $line;
		my @col=split /\t/,$line;

		next if(!exists $count_type{$col[2]});

		my @attributes=split /;/,$col[8];
		my %att;
		foreach my $attribute(@attributes){
			my @kv=split /=/,$attribute;
			$att{$kv[0]}=$kv[1];
		}

		if(!$att{ID}){
			my $id=$col[2].$count_type{$col[2]};
			$att{ID}=$id;	
		}
		#Gene dekho
		if($col[2] eq 'gene'){

			$gene{$att{ID}}->{phase}=$col[7];
			$gene{$att{ID}}->{type}=$col[2];
			$gene{$att{ID}}->{seqid}=$col[0];
			$gene{$att{ID}}->{source}=$col[1];
			$gene{$att{ID}}->{start}=$col[3];
			$gene{$att{ID}}->{end}=$col[4];
			$gene{$att{ID}}->{strand}=$col[6];
			$gene{$att{ID}}->{score}=$col[5];
			$gene{$att{ID}}->{attributes}=\%att;
			$count_type{$col[2]}++;
		}
     		#Transcripts ? mRNA ko store karo
		if(($col[2] eq 'mRNA')||($col[2] eq 'transcript')){
			$mrna{$att{ID}}->{type}=$col[2];
			$mrna{$att{ID}}->{seqid}=$col[0];
			$mrna{$att{ID}}->{source}=$col[1];
			$mrna{$att{ID}}->{start}=$col[3];
			$mrna{$att{ID}}->{end}=$col[4];
			$mrna{$att{ID}}->{strand}=$col[6];
			$mrna{$att{ID}}->{score}=$col[5];
			$mrna{$att{ID}}->{phase}=$col[7];
			#my %feats=$mrna{$att{ID}}

			$mrna{$att{ID}}->{attributes}=\%att; 
 
			#$gene{$att{Parent}}->{mrna}->{$att{ID}}=\%mrna;

			$count_type{$col[2]}++;

		}
		#Exon need a special attention
		if($col[2] eq 'exon'){
			
			#my $parent=$att{'Parent'};
			#my $gene_parent=$mrna{$parent}->{attributes}->{Parent};
			$exon{$att{ID}}->{type}=$col[2];
			$exon{$att{ID}}->{seqid}=$col[0];
			$exon{$att{ID}}->{source}=$col[1];
			$exon{$att{ID}}->{start}=$col[3];
			$exon{$att{ID}}->{end}=$col[4];
			$exon{$att{ID}}->{strand}=$col[6];
			$exon{$att{ID}}->{score}=$col[5];
			$exon{$att{ID}}->{attributes}=\%att; 
			$exon{$att{ID}}->{phase}=$col[7];
 			#$gene{$gene_parent}->{mrna}->{$att{Parent}}->{exon}->{$att{ID}}=\%exon;

 			$count_type{$col[2]}=int($count_type{$col[2]})+1;

		}
		#Adineta_ricciae	AUGUSTUS	CDS	30	155	0.98	+	0	transcript_id "g1.t1"; gene_id "g1";
		if($col[2] eq 'CDS'){
			
			if(exists $cds{$att{ID}}){
				#print STDERR "Warning: ".$col[2]." ID(".$att{ID}.") duplicated. Included anyway as ".$att{ID}."_dup".$dup."\n";

				$att{ID}=$att{ID}."_dup".$dup;
				$dup++;
			}
			
			#my $parent=$att{'Parent'};
			#my $gene_parent=$mrna{$parent}->{attributes}->{Parent};
			$cds{$att{ID}}->{phase}=$col[7];
			$cds{$att{ID}}->{type}=$col[2];
			$cds{$att{ID}}->{seqid}=$col[0];
			$cds{$att{ID}}->{source}=$col[1];
			$cds{$att{ID}}->{start}=$col[3];
			$cds{$att{ID}}->{end}=$col[4];
			$cds{$att{ID}}->{strand}=$col[6];
			$cds{$att{ID}}->{score}=$col[5];
			$cds{$att{ID}}->{attributes}=\%att; 
 			#$gene{$gene_parent}->{mrna}->{$att{Parent}}->{exon}->{$att{ID}}=\%exon;

 			$count_type{$col[2]}=int($count_type{$col[2]})+1;

		}
		if($col[2] eq 'intron'){

			
			#my $gene_parent=$mrna{$att{Parent}}->{attributes}->{Parent};
			$intron{$att{ID}}->{type}=$col[2];
			$intron{$att{ID}}->{seqid}=$col[0];
			$intron{$att{ID}}->{source}=$col[1];
			$intron{$att{ID}}->{start}=$col[3];
			$intron{$att{ID}}->{end}=$col[4];
			$intron{$att{ID}}->{strand}=$col[6];
			$intron{$att{ID}}->{score}=$col[5];
			$intron{$att{ID}}->{attributes}=\%att; #
			$intron{$att{ID}}->{phase}=$col[7];
 			#$gene{$gene_parent}->{mrna}->{$att{Parent}}->{exon}->{$att{ID}}=\%intron;
			$count_type{$col[2]}++;
		}

		if($col[2] eq 'start_codon'){

			#my $gene_parent=$mrna{$att{Parent}}->{attributes}->{Parent};
			$start_codon{$att{ID}}->{type}=$col[2];
			$start_codon{$att{ID}}->{seqid}=$col[0];
			$start_codon{$att{ID}}->{source}=$col[1];
			$start_codon{$att{ID}}->{start}=$col[3];
			$start_codon{$att{ID}}->{end}=$col[4];
			$start_codon{$att{ID}}->{strand}=$col[6];
			$start_codon{$att{ID}}->{score}=$col[5];
			$start_codon{$att{ID}}->{attributes}=\%att; 
			$start_codon{$att{ID}}->{phase}=$col[7];
 			#$gene{$gene_parent}->{mrna}->{$att{Parent}}->{exon}->{$att{ID}}=\%start_codon;
			$count_type{$col[2]}++;
		}

		if($col[2] eq 'stop_codon'){

			#my $gene_parent=$mrna{$att{Parent}}->{attributes}->{Parent};
			$stop_codon{$att{ID}}->{type}=$col[2];
			$stop_codon{$att{ID}}->{seqid}=$col[0];
			$stop_codon{$att{ID}}->{source}=$col[1];
			$stop_codon{$att{ID}}->{start}=$col[3];
			$stop_codon{$att{ID}}->{end}=$col[4];
			$stop_codon{$att{ID}}->{strand}=$col[6];
			$stop_codon{$att{ID}}->{score}=$col[5];
			$stop_codon{$att{ID}}->{attributes}=\%att; 
			$stop_codon{$att{ID}}->{phase}=$col[7];
 			#$gene{$gene_parent}->{mrna}->{$att{Parent}}->{exon}->{$att{ID}}=\%stop_codon;
			$count_type{$col[2]}++;
		}

	}# foreach line

	# exons
	foreach my $id(keys %exon){
		my $parent=$exon{$id}->{attributes}->{Parent};
		if(exists $mrna{$parent}){
			$mrna{$parent}->{exon}->{$id}=$exon{$id};
		}
		else{
			print LOG ("Warning: Parent($parent) for $id not found. Skipped.\n") if $param->{verbose};
		}
	}
	# cds
	foreach my $id(keys %cds){
		my $parent=$cds{$id}->{attributes}->{Parent};
		if(exists $mrna{$parent}){
			$mrna{$parent}->{CDS}->{$id}=$cds{$id};
		}
		else{
			print LOG ("Warning: Parent($parent) for $id not found. Skipped.\n") if $param->{verbose};
		}
	}

	# introns
	foreach my $id(keys %intron){
		my $parent=$intron{$id}->{attributes}->{Parent};
		if(exists $mrna{$parent}){
			#$mrna{$parent}->{intron}->{$id}=\$intron{$id}; 
			$mrna{$parent}->{intron}->{$id}=$intron{$id};
		}
		else{
			print LOG ("Warning: Parent($parent) for $id not found. Skipped.\n") if $param->{verbose};
		}
	}
	
	# start
	foreach my $id(keys %start_codon){
		my $parent=$start_codon{$id}->{attributes}->{Parent};
		if(exists $mrna{$parent}){
			$mrna{$parent}->{start_codon}->{$id}=$start_codon{$id};
		}
		else{
			print LOG ("Warning: Parent($parent) for $id not found. Skipped.\n") if $param->{verbose};
		}
	}

	# stop
	foreach my $id(keys %stop_codon){
		my $parent=$stop_codon{$id}->{attributes}->{Parent};
		if(exists $mrna{$parent}){
			$mrna{$parent}->{stop_codon}->{$id}=$stop_codon{$id};
		}
		else{
			print LOG ("Warning: Parent($parent) for $id not found. Skipped.\n") if $param->{verbose};
		}
	}

	# mrna - genes
	foreach my $id(keys %mrna){
		my $parent=$mrna{$id}->{attributes}->{Parent};
		if(exists $gene{$parent}){
			#$gene{$parent}->{mrna}->{$id}=\$mrna{$id};
			$gene{$parent}->{mRNA}->{$id}=$mrna{$id};
		}
		else{
			print LOG ("Warning: Parent($parent) for $id not found. Skipped.\n") if $param->{verbose};
		}
	}
	
	close GFF;
	return \%gene;
}

# print gene and all its childs. this function is RECURSIVE
sub printGene{
	my $gene=shift;
	my $bedFH=shift;
	my @valid_childs=("mRNA","exon","intron","CDS","start_codon","stop_codon","transcript");
		
	my $att=$gene->{attributes};
	my $attributes="";
	foreach my $k(keys %$att){
		$attributes=$attributes.";".$k."=".$$att{$k};
	}
	$attributes=~s/^;//;

	my $seqid = $gene->{seqid};
	my $source = $gene->{source};
	my $type = $gene->{type};
	my $start = $gene->{start};
	my $end = $gene->{end};
	my $score = $gene->{score};
	my $strand = $gene->{strand};
	my $phase = $gene->{phase};

	my $gene_line=join("\t",$seqid,$source,$type,$start,$end,$score,$strand,$phase,$attributes);
	#print $gene_line."\n";


	foreach my $k(keys %$gene){
		foreach my $valid(@valid_childs){
			if($k eq $valid){
				my $childs=$gene->{$k};
				foreach my $child_id(keys %$childs){
					my $child=$childs->{$child_id};
					printGene($child, $bedFH);
				}

			}
		}
	}

}


# print gene and all its childs. this function is RECURSIVE
sub printGeneBed{
	my $gene=shift;
	my $bedFH=shift;
	my @valid_childs=("mRNA","exon","intron","CDS","start_codon","stop_codon","transcript");
		
	my $att=$gene->{attributes};
	my $attributes="";
	foreach my $k(keys %$att){
		$attributes=$attributes.";".$k."=".$$att{$k};
	}
	$attributes=~s/^;//;
	my $seqid = $gene->{seqid};
	#Need to understand and condider ZERO based and ONE based system https://www.biostars.org/p/84686/
	my $start = $gene->{start} - 1;
	my $end = $gene->{end};
	my $phase = $gene->{phase};
	my $score = $gene->{score};
	my $strand = $gene->{strand};

	my $gene_line=join("\t",$seqid,$start,$end,$phase,$score,$strand,$attributes);
	print $bedFH $gene_line."\n";


	foreach my $k(keys %$gene){
		foreach my $valid(@valid_childs){
			if($k eq $valid){
				my $childs=$gene->{$k};
				foreach my $child_id(keys %$childs){
					my $child=$childs->{$child_id};
					printGene($child, $bedFH);
				}

			}
		}
	}

}

# filter transcripts on a gene structure. if returns 0, no transcript associated to the gene
# options are, "longest", that keep only the longest transcript; and "shortest", that keep only shortest transcript. if transcripts are equally length, keep just one. the last that it finds
# usage: $filtered_gene_structure=filterGeneTranscripts($gene_structure, 'longest')
sub filterGeneTranscripts{
	my $gene=shift;
	my $opt=shift;

	my $att=$gene->{attributes};
	my $id=$$att{ID};
	if(!exists $gene->{mRNA}){
		print STDERR "Gene $id not contain mRNA childs. Skipped\n";
		return 0;
	}

	my $mrnas=$gene->{mRNA};

	my $last_len=0;
	my $valid_id;
	foreach my $id(keys %$mrnas){
		my $len=$$mrnas{$id}->{end} - $$mrnas{$id}->{start};
		if($opt eq 'longest'){
			if($len>=$last_len){
				$valid_id=$id;
			}
			else{
				print STDERR "Warning: $id length ($len) < 0?\n";
			}
		}
		if($opt eq 'shortest'){
			if($last_len<=0){
				$last_len=$len;
				$valid_id=$id;
			}
			else{
				if($len<$last_len){
					$last_len=$len;
					$valid_id=$id;
				}
			}
		}
		if(($opt ne 'longest')&&($opt ne 'shortest')){
			die "$opt is not valid option for filterGeneTranscripts";
		}
	}
	my $filtered_gene=$gene;
	my $valid_mrna=$gene->{mRNA}->{$valid_id};
	$gene->{mRNA}=();
	$gene->{mRNA}->{$valid_id}=$valid_mrna;

	return $gene;
}

#Create all mandatory files in intermediate/gff folder 
sub gff2all_files {
my ($file_fasta,$file_gff,$path)=@_;
use strict; 
use warnings; 
use Bio::Seq; 
use Bio::SeqIO; 
use Bio::DB::Fasta;
#add a help message here my $num_args=$#ARGV + 1; if ($num_args != 4) {
#	print "\nUsage: gff2perl Genome.fasta Annotation.gff 
#	OutputPrefix \n\n"; exit;
#}

$| = 1; # Flush output 
my $outfile_cds = Bio::SeqIO->new( -format => 'fasta', -file => ">$path/local.cds.fasta" ); 
my $outfile_pep = Bio::SeqIO->new( -format => 'fasta', -file => ">$path/local.pep.fasta" ); 
my $outfile_cdna = Bio::SeqIO->new( -format => 'fasta', -file => ">$path/local.cdna.fasta" ); 
my $outfile_gene = Bio::SeqIO->new( -format => 'fasta', -file => ">$path/local.gene.fasta" ); 
my $outfile_upstream3000 = Bio::SeqIO->new( -format => 'fasta', -file => ">$path/local.upstream3000.fasta" ); 
my $outfile_exon = Bio::SeqIO->new( -format => 'fasta', -file => ">$path/local.exon.fasta");

###### Output type description ######
#my $file_fasta = $ARGV[0]; 
my $db = Bio::DB::Fasta->new($file_fasta); 
print ("Genome fasta parsed\n");

### Second, parse the GFF3
my %CDS; my %CDNA; my %EXON; my $mRNA_name; my $mRNA_loc; my $frame=''; 
open GFF, "<$file_gff" or die $!; 
while ( my $line = <GFF> ) {
    chomp $line;
    if ($line =~ /^\s*#/) { next; }
    my @array = split( "\t", $line );
    my $type = $array[2]; #Name of types
    if ($type eq 'gene' || $type eq 'mt_gene' ) { #If there is any mito genes
        my @attrs = split( ";", $array[8] );
        $attrs[0] =~ s/ID=//;
        my $gene_name = $attrs[0]; # Name of gene in gff file
        my $gene_loc = "$array[0]:$array[3]-$array[4]"; # Coordinates of the gene -- it is used throught
        my $gene_start = $array[3];
        my $gene_end = $array[4];
        my $gene_seq = $db->seq( $array[0], $gene_start, $gene_end );
        my $output_gene = Bio::Seq->new(
            -seq => $gene_seq,
            -id => $gene_loc,
            -display_id => $gene_loc,
            -alphabet => 'dna',
        );
        # The upstream 3000 - to print the sequences present upstream
        my $upstream_start;
        my $upstream_end;
        if($array[6] eq '+') {
            $upstream_start=$gene_start-3000;
            $upstream_end=$gene_start-1;
        }
        elsif ($array[6] eq '-') {
            $upstream_start=$gene_end+1;
            $upstream_end=$gene_end+3000;
        }
        my $upstream_seq = $db->seq( $array[0], $upstream_start, $upstream_end );
        my $output_upstream3000 = Bio::Seq->new(
            -seq => $upstream_seq,
            -id => $gene_loc."_upstream3000",
            -display_id => $gene_loc."_upstream3000",
            -alphabet => 'dna',
        );
        # Reverse Complement if the frame is minus
        if($array[6] eq '+') {
        }
        elsif ($array[6] eq '-') {
            $output_gene = $output_gene->revcom();
            $output_upstream3000 = $output_upstream3000->revcom();
        }
        else {
            die "Unknown frame! At line $. of the GFF\n";
        }
	# avoid empty fasta headers
        if (length($gene_seq) != 0) {
	$outfile_gene->write_seq($output_gene);
	}
	if (length($upstream_seq) != 0) {
        $outfile_upstream3000->write_seq($output_upstream3000);
	}
    }
#CDS sequences extraction using mRNA
    if ( ( $type eq 'mRNA' || $type eq 'transcript' ) and ( $. > 2 ) ) {
        # CDS: Collect CDSs and extract sequence of the previous mRNA
        my $mergedCDS_seq='';
	# WARNING we must sort by $cds_coord[1]
        foreach my $key (sort {$a <=> $b} keys %CDS) { # Ascending numeric sort of the starting coordinate
            my $coord = $CDS{$key};
            my @cds_coord = split( " ", $coord );
            my $cds_seq = $db->seq( $cds_coord[0], $cds_coord[1], $cds_coord[2] );
            $mergedCDS_seq .= $cds_seq;
        }
        
        my $output_cds = Bio::Seq->new(
            -seq => $mergedCDS_seq,
            -id => $mRNA_loc,
            -display_id => $mRNA_loc,
            -alphabet => 'dna',
        );
        if ($frame eq '-') {
            $output_cds = $output_cds->revcom();
        }
	#translate CDS to peptide for protein sequence
        my $output_pep = $output_cds->translate();
	#write to file
	if (length($mergedCDS_seq) != 0) {
        $outfile_cds->write_seq($output_cds);
	}
	if (length($mergedCDS_seq) != 0) {
        $outfile_pep->write_seq($output_pep);
	}
#exons should be able to add exon output here since exons will be useful 
#in gene models for other organisms can be added in the EVM program
        my $mergedEXON_seq='';
        foreach my $key (sort {$a <=> $b} keys %EXON) { # Ascending numeric sort of the starting coordinatg
            my $coord = $EXON{$key};
            my @exon_coord = split( " ", $coord );
            my $exon_seq = $db->seq( $exon_coord[0], $exon_coord[1], $exon_coord[2] );
            $mergedEXON_seq .= $exon_seq;
        }
       my $output_exon = Bio::Seq->new(
            -seq => $mergedEXON_seq,
            -id => $mRNA_loc,
            -display_id => $mRNA_loc,
            -alphabet => 'dna',
        );
        if ($frame eq '-') {
            $output_exon = $output_exon->revcom();
        }
	#write to file
        if (length($mergedEXON_seq) != 0) {
        $outfile_exon->write_seq($output_exon);
        }
        # CDNA: Collect UTRs and CDSs and extract sequence of the 
        # previous mRNA
        my $mergedCDNA_seq="";
        foreach my $key (sort {$a <=> $b} keys %CDNA) { # Ascending numeric sort of the starting coordinate
            my $coord = $CDNA{$key};
            my @cds_coord = split( " ", $coord );
            my $cds_seq = $db->seq( $cds_coord[0], $cds_coord[1], $cds_coord[2] );
            $mergedCDNA_seq .= $cds_seq;
        }
        my $output_cdna = Bio::Seq->new(
            -seq => $mergedCDNA_seq,
            -id => $mRNA_loc,
            -display_id => $mRNA_loc,
            -alphabet => 'dna',
        );
        if ($frame eq '-') {
            $output_cdna = $output_cdna->revcom();
        }
        if (length($mergedCDNA_seq) != 0) {
	$outfile_cdna->write_seq($output_cdna);
	}
        # Now initialize the next mRNA
        my @attrs = split( ";", $array[8] );
        $attrs[0] =~ s/ID=//;
        $mRNA_name = $attrs[0];
        $mRNA_loc = "$array[0]:$array[3]-$array[4]";
        $frame=$array[6];
        %CDS = (); %CDNA = (); # Empty the chunk arrays
	%EXON = (); %EXON = (); #Empty the EXON chunk arrays
    }
    elsif ( $type eq 'mRNA' ) { # First mRNA
        my @attrs = split( ";", $array[8] );
        $attrs[0] =~ s/ID=//;
        $mRNA_name = $attrs[0];
        $mRNA_loc = "$array[0]:$array[3]-$array[4]";
        $frame=$array[6];
    }
    elsif ( $type eq 'CDS' ) {
        my $cds_coord = $array[0] . " " . $array[3] . " " . $array[4];
        $CDS{$array[3]}=$cds_coord;
        $CDNA{$array[3]}=$cds_coord;
    }
    elsif ($type eq 'UTR' ) {
        my $utr_coord = $array[0] . " " . $array[3] . " " . $array[4];
        $CDNA{$array[3]}=$utr_coord;
    }
    elsif ($type eq 'exon' ) {
	my $exon_coord = $array[0] . " " . $array[3] . " " . $array[4];
	$EXON{$array[3]}=$exon_coord;
    }
}
close GFF;
}

1;

