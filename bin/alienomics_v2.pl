#!/usr/bin/env perl

#TODO 
#Use only one config file - DONE
#Create summary table 
#keep all the message in subs

#Turning off BioPerl warnings
$Bio::Root::Root::DEBUG = -1;
$Bio::DB::IndexedBase::DEBUG = -1;

# Activate this command if working;
# I use cpanm module on my CECI cluster. Hence, in order to use local module we need to acrivate it.
#system ("eval `perl -I ~/perl5/lib/perl5 -Mlocal::lib`"); # Where is ur perl lib

=head1 NAME

alienomics.pl - alienomics script by Jitendra Narayan and Paul Simon

=head1 SYNOPSIS

perl alienomics.pl --conf or -c ../test.conf

=cut

use strict;
use warnings;

use Cwd;
use File::chdir;
use File::Copy;
use POSIX;

use Tie::File;
use Try::Tiny;
use Data::Dumper;
use File::Basename;
use Bio::SeqIO;
use FindBin;
use Math::Trig;
use File::Copy;
use File::Remove;
use Capture::Tiny ':all';
use File::Temp qw(tempfile);
use File::Spec::Functions qw(rel2abs);
use Statistics::Multtest qw(BY qvalue);
use File::Path qw(make_path remove_tree);
use Statistics::Distributions qw(chisqrprob);
use Getopt::Long;
#use Statistics::R;
use Math::Round;
use File::Find;
use Bio::DB::Taxonomy;
use Pod::Usage;
use lib "$FindBin::Bin/.";
use List::Util qw(max);
use List::AllUtils qw(sum);
use List::MoreUtils qw( minmax );
use Fcntl qw/ :flock /;
use Text::Trim;
use Bio::Taxon;
use Bio::Tree::Tree;

use AOM::Module;
use AOM::gffHandler;
use AOM::rRNAHandler;
use AOM::buscoHandler;
use AOM::refHandler;
use AOM::genomeHandler;
use AOM::tnfHandler;
use AOM::plotHandler;
use AOM::covHandler;
use AOM::HGTcategorization;
use AOM::transHandler;

# Basics mandatory alienomics variables
my (
$create_conf, 		# If you wish ALIENOMICS to create an empty configuration file for you
$outfile, 		# Name for ALIENOMICS's main configuration file
$version, 		# ALIENOMICS's version?
$help,			# If you are asking for help
$man,			# Manual for detail help
$Input_genome,		# Path to directory containing input genome
$blastdb_path,		# Path to directory containing blast database
$diamond_path,		# Path to directory containing diamond executables
$diamonddb_path,	# Path to directory containing diamond database

$taxdump_path,		# Path to directory containing taxonomy data from NCIB
$accession2taxid_path,	# Path to directory containing taxonomy accession2taxid from NCBI
$augustusconfig_path,	# Path to augustus config data.
$homology_file_path,	# Path to homology relationship file (OrthoMCL 1.4 format)
$output_dir,		# Path to directory where ALIENOMICS will write results
$genetic_code,		# Genetic code to be used when translating coding sequences
$reference_genome,	# Reference genome ID.
$max_processors,	# Maximum number of processors available. $max_processors=`grep "^processor" /proc/cpuinfo | tail -n 1 | awk '{print $3}'`;
$gc_filter,		# GC thresholds range sepeared with DOT

);

# Default settings here for alienomics
my $current_version = "0.2.0:Camel";	# Alienomics version
my $stringency = "standard";		# Used to automatically create a new configuration file with default values.
my $conf_file = "test.conf";		# The path to a valid alienomics configuration file. This file contains all parameters needed to execute  ALIENOMICS.

print <<'WELCOME';

 /\  |    | |__  |\ | /  \  |\/| | /  ` /__` 
/~~\ |___ | |___ | \| \__/  |  | | \__, .__/ v0.2

Automated "sensitive" contamination filtration and HTG detection tool
Contact: jnarayan81@gmail.com and polo.simion@gmail.com for support

WELCOME

$|++; #OUTPUT_AUTOFLUSH

# Get options for alienomics
GetOptions(
	"create_conf" 		=> \$create_conf,
	"stringency=s" 		=> \$stringency,
	"conf|c=s" 		=> \$conf_file,
	"Input_genome=s" 	=> \$Input_genome,
	"blastdb_path=s" 	=> \$blastdb_path,
	"diamond_path=s" 	=> \$diamond_path,
	"taxdump_path=s" 	=> \$taxdump_path,
	"accession2taxid_path=s" => \$accession2taxid_path,
	"augustusconfig_path=s" => \$augustusconfig_path,
	"homology_file_path=s" 	=> \$homology_file_path,
	"output_dir=s" 		=> \$output_dir,
	"genetic_code=i"	=> \$genetic_code,
	"outfile|o=s" 		=> \$outfile,
	"reference_genome|r=s" 	=> \$reference_genome,
	"max_processors|m=i" 	=> \$max_processors,
	"gc_filter|g=i" 	=> \$gc_filter,
	"version|v" 		=> \$version,
	"help|h" 		=> \$help
) or pod2usage(); # when no command line options are present

#pod2usage(-verbose => 2) if ($man);
#pod2usage(-verbose => 1) if ($help);
#pod2usage(-msg => 'Please supply a valid filename.') unless ($conf_file && -s $conf_file);

# Printing alienomics version
if ($version) {
    print "\n  ALIENOMICS version $current_version\n\n";
    exit;
}

#Check if allready running same script anywhere.
flock(DATA,LOCK_EX|LOCK_NB)
  or  die "This script ($0) is already running.\n Wait the first instances to finish and try again later\n Sorry for inconvenience\n";

#Check if other instances is running
my $cmd = 'ps -eaf | grep -c "\< '."$0".'\>"'; 
chomp(my $instances = `$cmd`);
if($instances > 1) { print "Other instances of your program is running\n Please let them complete first\n"; exit; }


# Creating alienomics configuration file if asked to
if ($create_conf) {
  AOM::Module::create_conf_file($stringency, $Input_genome, $homology_file_path, $output_dir, $outfile, $genetic_code, $max_processors, $gc_filter);
  exit;
}

# Checking if users are asking for help or forgot to provide an essential parameter
if (($help)||(!defined $conf_file)||(!-e $conf_file)) {
  AOM::Module::help($current_version);
  exit;
}

# Used to measure total execution time
my $start_time = time();

my $project_config_file = $conf_file;

# Absolute path of the current working directory
my $ALIENOMICS_path = dirname(rel2abs($0));
print "Path of the dir:\t$ALIENOMICS_path --------\n";
print "Path of the configuration file:\t$conf_file\n";

# Parameters_ref - stores all user-defined parameters, such as file locations and program parameters
my $param_ref = AOM::Module::read_config(\$project_config_file, \$ALIENOMICS_path);  # stores the configuration in a hash reference

# test fun GFF gene extraction
#AOM::Module::extractCoordinates ( "test.fa", "test.gff", "Test");
print "\nChecking manadatory programs for Alienomics\n";
AOM::commonSubs::checkPrograms();

# Check all the parameters for their correctness
print "Checking parameters\n";
AOM::Module::check_parameters($param_ref); #checkin if user setted the parameters right

# Create provided name as folder in requested location
#AOM::createDIR($param_ref, );

# Delete the directory if already exist
if (-e $param_ref->{project_dir_path}) {
	print "Deleting the existing directory named :-$param_ref->{project_dir_path}\n";
	remove_tree( $param_ref->{project_dir_path});
}

# Creating the needed directories if they don't exist
if (!-e $param_ref->{project_dir_path}) {
	 mkdir ($param_ref->{project_dir_path}) || die ("Couldn't create the directory specified as '$param_ref->{project_dir_path}', check if you are trying to create a subdirectory in a non-existent directory or if you don't have permission to create a directory there.\n");
}
else {
  	die("Directory $param_ref->{project_dir_path} already exists.\n");
}

if (!-e "$param_ref->{project_dir_path}/results") { mkdir ("$param_ref->{project_dir_path}/results") || die ("Couldn't create the directory with the results of ALIENOMICS's analysis.\nDetails: $!\n"); }
if (!-e "$param_ref->{project_dir_path}/intermediate_files") { mkdir ("$param_ref->{project_dir_path}/intermediate_files/") || die ("Couldn't create the directory with the steps of ALIENOMICS's analysis.\nDetails: $!\n"); }
if (!-e "$param_ref->{project_dir_path}/intermediate_files/ref") { mkdir ("$param_ref->{project_dir_path}/intermediate_files/ref") || die ("Couldn't create the directory with the steps of ALIENOMICS's analysis.\nDetails: $!\n"); }
if (!-e "$param_ref->{project_dir_path}/intermediate_files/gff") { mkdir ("$param_ref->{project_dir_path}/intermediate_files/gff") || die ("Couldn't create the directory with the steps of ALIENOMICS's analysis.\nDetails: $!\n"); }
if (!-e "$param_ref->{project_dir_path}/intermediate_files/bed") { mkdir ("$param_ref->{project_dir_path}/intermediate_files/bed") || die ("Couldn't create the directory with the steps of ALIENOMICS's analysis.\nDetails: $!\n"); }
if (!-e "$param_ref->{project_dir_path}/intermediate_files/rRNA") { mkdir ("$param_ref->{project_dir_path}/intermediate_files/rRNA") || die ("Couldn't create the directory with the steps of ALIENOMICS's analysis.\nDetails: $!\n"); }
if (!-e "$param_ref->{project_dir_path}/intermediate_files/busco") { mkdir ("$param_ref->{project_dir_path}/intermediate_files/busco") || die ("Couldn't create the directory with the steps of ALIENOMICS's analysis.\nDetails: $!\n"); }
if (!-e "$param_ref->{project_dir_path}/intermediate_files/blast") { mkdir ("$param_ref->{project_dir_path}/intermediate_files/blast") || die ("Couldn't create the directory with the steps of ALIENOMICS's analysis.\nDetails: $!\n"); }
if (!-e "$param_ref->{project_dir_path}/intermediate_files/tnf") { mkdir ("$param_ref->{project_dir_path}/intermediate_files/tnf") || die ("Couldn't create the directory with the steps of ALIENOMICS's analysis.\nDetails: $!\n"); }
if (!-e "$param_ref->{project_dir_path}/localDB") { mkdir ("$param_ref->{project_dir_path}/localDB") || die ("Couldn't create the directory with the steps of ALIENOMICS's analysis.\nDetails: $!\n"); }
if (!-e "$param_ref->{project_dir_path}/localDB/refDB") { mkdir ("$param_ref->{project_dir_path}/localDB/refDB") || die ("Couldn't create the directory with the steps of ALIENOMICS's analysis.\nDetails: $!\n"); }
if (!-e "$param_ref->{project_dir_path}/localDB/rRNADB") { mkdir ("$param_ref->{project_dir_path}/localDB/rRNADB") || die ("Couldn't create the directory with the steps of ALIENOMICS's analysis.\nDetails: $!\n"); }
if (!-e "$param_ref->{project_dir_path}/localDB/buscoDB") { mkdir ("$param_ref->{project_dir_path}/localDB/buscoDB") || die ("Couldn't create the directory with the steps of ALIENOMICS's analysis.\nDetails: $!\n"); }
if (!-e "$param_ref->{project_dir_path}/intermediate_files/stat") { mkdir ("$param_ref->{project_dir_path}/intermediate_files/stat") || die ("Couldn't create the directory with the steps of ALIENOMICS's analysis.\nDetails: $!\n"); }


# Copy file to the locations - just to keep the used config file for future reference
copy($project_config_file, "$param_ref->{project_dir_path}/project_config");

# Write the log files
open (my $LOG, ">", "$param_ref->{project_dir_path}/log") || die ('Could not create log file in ', $param_ref->{project_dir_path}, '. Please check writing permission in your current directory', "\n");
open (my $LOG_ERR, ">", "$param_ref->{project_dir_path}/log.err") || die ('Could not create log.err file in ', $param_ref->{project_dir_path}, '. Please check writing permission in your current directory', "\n");
open (my $SUMMARY, ">", "$param_ref->{project_dir_path}/results/$param_ref->{summary}") || die ('Could not create summary file. Please check writing permission in your current directory', "\n");
open (my $RESULT, ">", "$param_ref->{project_dir_path}/results/$param_ref->{result}") || die ('Could not create result file. Please check writing permission in your current directory', "\n");
open (my $TAXDETECTION, ">", "$param_ref->{project_dir_path}/taxonomy_detection.log") || die ('Could not create result file. Please check writing permission in your current directory', "\n");

# Write the header
print $RESULT "Scaffold\tlength\tGC_scaf\ttotalGenes\trefBlast\tgenomeBlast\tmaxVal,allNames\trRNA\tbusco\ttnfCnt\tallTNF\talienCnt\tscaffoldScore\tnewglobalscore\n";
open (my $SCORE, ">", "$param_ref->{project_dir_path}/results/$param_ref->{score}") || die ('Could not create result file. Please check writing permission in your current directory', "\n");
print $SCORE #"Scaffold\tlength\tGC_scaf\tgene\tGC_gene\tgeneScore_modified\tFinal_sigScore\tgeneClass\tRef_score\tBusco_score\trRNA_score\tfeeder_score\tUniprot_score\tTNF_score\tGC_bonus_score\tgeneScore\tcoverage\tcov_bonus\n";
"Scaffold\tlength\tGC_scaf\tgene\tRef_score\tgeneral_blast_score\tBusco_score\trRNA_score\tfeeder_score\tTNF_score\tcoverage\tcov_score\tGC\tGC_score\tgene_score\tFinal_sigScore\tTPM\tExpression_score\tsigScore_forScaf\n";



##############################################################################
### Parsing input files

# Get one from a NCBI taxonomy database
my $dbh = Bio::DB::Taxonomy->new(-source   => 'flatfile',
                                 -directory=> "$param_ref->{taxdump_path}",
                                 -nodesfile=> "$param_ref->{taxdump_path}/nodes.dmp",
                                 -namesfile=> "$param_ref->{taxdump_path}/names.dmp");
#to check
#parseTax('homo_sapiens', 'metazoa');
#exit;

#$abc = $names{$speciesId};
my (%txn_names, %txn_nodes, @array2d, %finalHash);
my $taxonomy_directory='taxdump'; #Where is your taxdump
#print "Loading the taxonomy $taxonomy_directory -- gonna fill up your RAM :) sorry for that\n";
my ($files,$dirs)=getDirectoryFiles("$param_ref->{taxdump_path}");
process_file($_) for @$files;    # This subrutine create a hash of names and nodes file.
#print "\nProcessing done, yahh ! I guess you have high quality RAM, rich hmmm\n";
#my $abc_nam = $txn_names{9606};
#print $abc_nam; exit;

# Parse the genome file and store in Hash
my $date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
print "\n# Parsing genome file ($date)\n";
my ($sequence_data_ref, $id2tmp_id_ref, $tmp_id2id_ref) = AOM::Module::parse_genome_files($param_ref, $LOG);

# Predict genes with augustus
if ($param_ref->{predict_gene}) {
	AOM::Module::predictGeneAug("$param_ref->{project_dir_path}/intermediate_files/gff/genome.fa", $param_ref, $param_ref->{predict_gene});
}
else {
	print "Using local gff file for analysis\n";
	if (!$param_ref->{local_gff}) { print "Did you forgot to provide me the gff file\n Please check the config file for detail\n";  exit;}
	print "Using local gff file:\t$param_ref->{local_gff}\n";
	copy("$param_ref->{local_gff}","$param_ref->{project_dir_path}/intermediate_files/gff/preGenes.gff") or die "Copying gff file failed: $!";
}

# Store the genes detail
$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
print "\n# Parsing genes ($date)\n";
my $genes=AOM::gffHandler::gff2gene("$param_ref->{project_dir_path}/intermediate_files/gff/preGenes.gff", $param_ref);

#Lets create all files -- the coordinates need to check ;)
AOM::gffHandler::gff2all_files("$param_ref->{project_dir_path}/intermediate_files/gff/genome.fa","$param_ref->{project_dir_path}/intermediate_files/gff/preGenes.gff","$param_ref->{project_dir_path}/intermediate_files/gff");

## my %tmpGenes = %$genes;
## foreach my $aaa (%tmpGenes) {print "$aaa\t $tmpGenes{$aaa}{strand}\n";} exit;

#Lets store the coverage information in Hash
#my $covHash_ref; #To store the hash ref values ---global
#if ($param_ref->{coverage} >=1 ) {
#AOM::covHandler::extractMappingCov($param_ref->{bam},"$param_ref->{project_dir_path}/intermediate_files/bed/gene.bed", $param_ref->{max_processors},"$param_ref->{project_dir_path}/intermediate_files/bed/gene.cov");
#$geneCov = `samtools depth -aa -r $accession_number $param->{bam} | awk '{ c++; s+=\$3; } END {k=(s/c); print k;}'`;
#$covHash_ref=storeInHash("$param_ref->{project_dir_path}/intermediate_files/bed/gene.cov");
#}

# Store transcriptome data 
print "\n# Checking expression ($date)\n";
$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
my $transHash_ref=AOM::transHandler::transKallister($param_ref);
#print $transHash_ref->{'Chrom_1:2155-5638'};

##############################################################################
### Parse genes = GC + coverage

# convert to bed and add gene coverage using BAM file (provided by the user) -- computational expensive step 
$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
print "\n# Convert to bed ($date)\n";
# THIS SHOULD BE PARALLELIZED ?!?! (notably because it includes coverage computation). Maybe with this ? :https://perlmaven.com/speed-up-calculation-by-running-in-parallel
my $geneSet_ref = gffGene2Bed($genes, $param_ref, $transHash_ref);
my %geneSet=%$geneSet_ref;
#my $geneSet_ref;
#my %geneSet;

## extractIds("$param_ref->{project_dir_path}/intermediate_files/bed/gene.bed");

##############################################################################
### Blasting steps

### rRNABlaster
$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
print "\n# Blasting on rRNA database ($date)\n";
my ($rRNAHash_ref, $rRNAmin, $rRNAmax) = AOM::rRNAHandler::rRNABlaster($param_ref, $geneSet_ref);
my %rRNAHash = %$rRNAHash_ref;

### buscoBlaster
$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
print "\n# Blasting on BUSCO database ($date)\n";
my ($buscoHash_ref, $BUSCOmin, $BUSCOmax) = AOM::buscoHandler::buscoBlaster($param_ref, $geneSet_ref);
my %buscoHash = %$buscoHash_ref;

### refBlaster
$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
print "\n# Blasting on genomic REFERENCE (user-provided) ($date)\n";
my ($refHash_ref, $REFmin, $REFmax) = AOM::refHandler::refBlaster($param_ref, $geneSet_ref);
my %refHash = %$refHash_ref;
my %genome=%{$sequence_data_ref}; # foreach my $sss (keys %genome) { print $genome{$sss}{nuc_seq} }
my %finalStats;

# Blast genes against large database and Annotated with taxdb
# Only keeping query with no hits
$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
my ($Gmin, $Gmax); my %genomeHash; # Declare globally
if ($param_ref->{ActivateBIGBLAST} eq 'on') {
        # Lets load the taxonomy database, once loaded used many places -- mostly used in exhaustive blast result
	$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
	print "\n# Loading taxonomical database ($date)\n";
	# This subs is only used while blast ... that why move here (mainly in two modules) 
	#I set it off for setting, keep it on
	my ($taxidP_ref, $taxidT_ref, $taxidN_ref) = AOM::Module::load_nodes_names ("$param_ref->{taxdump_path}/nodes.dmp","$param_ref->{taxdump_path}/names.dmp");

	print "\n# Blasting against exhaustive database set ON ($date)\n";
	my ($genomeHash_ref, $Gmin, $Gmax) = AOM::genomeHandler::genomeBlaster($param_ref, $taxidP_ref, $taxidT_ref, $taxidN_ref, $geneSet_ref, $param_ref, $dbh, \%txn_names);
	%genomeHash = %$genomeHash_ref;
}
else { print "\n# Blasting against exhaustive database set OFF ($date)\n";}

##############################################################################

my %tnfHash;
# trying to render tnf optional:
if ( $param_ref->{ActivateTNF} eq 'on') {
	# Predict the TNF 
	$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
	print "\n# TNF computations ($date)\n";

	# only for remaining noHits sequences :
	#my $coord_ref = AOM::tnfHandler::calculateTNF("$param_ref->{project_dir_path}/intermediate_files/blast/blastNoHits.fa", 0, 4);
	# for all genes :
	my $coord_ref = AOM::tnfHandler::calculateTNF("$param_ref->{project_dir_path}/intermediate_files/bed/gene.fa", 0, 4);

	my %tnfHash;
	foreach my $freq (@$coord_ref) {
		my @tnf=split /\t/, $freq;
		my $gName = shift @tnf; 					# Shift the name of the fasta sequence
		my ($name, $nearest, $min_dist) = AOM::tnfHandler::SearchTNF(\@tnf, "$param_ref->{tnfDB_file}", "$param_ref->{project_dir_path}/intermediate_files/tnf/tnfNoHits.tnf", $gName);
		my @lenCor = split /\-/, (split /\:/, $gName)[-1]; 		#Get the sequence len
		my $len = $lenCor[1]-$lenCor[0];
		$tnfHash{$name}="$nearest-$min_dist-$len";
		#print "$name, $nearest, $min_dist, $len <<<<<<<<<<<<<<\n";
	}
}
else {
	$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
	print "\n# No TNF prediction will be run ($date)\n";
}

##############################################################################
print "\nYou opted for --> $param_ref->{mode} <-- mode of Alienomics\n";
#AOM::Module::pLine("_");

# Everything is ready to roll - Lets start checking
$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
print "\n# Gathering and computing gene scores ($date)\n";
if ($param_ref->{mode} eq 'gene') {
	# Lets check for alien in all the sequence now
	foreach my $accession_number (keys %genome) {
		#AOM::Module::pLine(".");
		#print "$accession_number\n";

		my $busco=0;
		my $rRNA=0;
		my $refBlast=0;
		my $genomeBlast=0;
		my $tnfCnt=0;
		my %namesCnt;
		my %tnfSps;
		my $alienCnt=0;
		my @allNames;
		my @allTNF;
		my @allGenesScoreA;
		my @allGenesScoreB;
		my @allGenesLength;
		my %spsCounts;
		my $division;
		my $roughScore;
		my ($genesLoc_ref) = subsetGFF($genes, $accession_number, "$genome{$accession_number}{len}"); # Gene location detail

		#Total number of genes $accession_number contig
		my $totalGenes = scalar @$genesLoc_ref;
		#my %intronStat = %$intronStatus_ref;

		#Check for all genes
		foreach my $gene (@$genesLoc_ref) {
			# Store gene classification
			my $geneClass="NA";
			# Store rRNA hits detail in one contigs
			$rRNA++ if exists $rRNAHash{$gene};
			# Store buscoDB hits detail in one contigs
			$busco++ if exists $buscoHash{$gene};
			# Store reference hits detail in one contigs
			$refBlast++ if exists $refHash{$gene};
			# if match on REF
			if (exists $refHash{$gene}) { $geneClass="reference";}
			# General blast results count
			$genomeBlast++ if exists $genomeHash{$gene};
			if (exists $genomeHash{$gene}) {
				$namesCnt{$genomeHash{$gene}->[12]}++;
				$spsCounts{$genomeHash{$gene}->[12]}++; #Store only species here
				$geneClass=$genomeHash{$gene}->[12];
			}
			#TNF results count -- one by one
			$tnfCnt++ if exists $tnfHash{$gene};
			my $tnf=0;
			my $tnfLen=0;
			my $min_dist_score=0;
			if (exists $tnfHash{$gene}) {
				my @tnfHits = split /\-/, $tnfHash{$gene}; # The values of the tnf file .. only one value
				$tnfSps{$tnfHits[0]}++;
				#foreach (@tnfHits) { print "$_\n"; }
				#print "$tnfHits[-1]/sqrt($tnfHits[1] >---<\n";
				$min_dist_score=1/(1+$tnfHits[1]); #\C2 0:nearestname/ 1:mindistance/ 2:len
				#$min_dist_score=$tnfHits[-1]/sqrt($tnfHits[1]); # length / sqrt (min_dist)
				my $closestSpsName=lc($tnfHits[0]);
				#print "$tnfHits[0] <<<<<--->>>>> $tnfHits[1]\n";

				#Equal 1 if present $param_ref->{level_upto}
				print $TAXDETECTION "TNF parseTax:\t$gene\t$closestSpsName\t$param_ref->{level_upto})\t";
				my $tnfTaxDecision=parseTax($closestSpsName, $param_ref->{level_upto});
				if ($tnfTaxDecision) { $min_dist_score = ($min_dist_score * 1);
					#print $TAXDETECTION "\tSELF\n";
				}
				else { $min_dist_score = ($min_dist_score * -1); $alienCnt++;
					#print $TAXDETECTION "\tALIEN\n";
				}
				$tnf=1;
				$tnfLen=$tnfHits[-1];
				$geneClass=$closestSpsName;
				print $TAXDETECTION "\t$min_dist_score\t$tnfHits[0]\n";
			}

			#### starting working on gene score

			# GC SCORING
			my @valGC = split /\:/, $param_ref->{gc_filter}; #Provided in config file

			my $GCpenalty=0;
			my $Win="2";    # 
			my $Wout="5";	# 
			if (($valGC[0] <= $geneSet{$gene}{gc}) && ($valGC[1] >= $geneSet{$gene}{gc})) {
				#$GCpenalty=tanh(3*($geneSet{$gene}{gc}-$valGC[0])/$Win)-tanh(3*($geneSet{$gene}{gc}-$valGC[1])/$Win)-1;
				# second formula for positive score => strictly equal 0 (a good coverage can not help decipher between alien and self):
				$GCpenalty=0;
			}
			elsif ( $geneSet{$gene}{gc} < $valGC[0] ) {
				$GCpenalty= max(-(($geneSet{$gene}{gc}-$valGC[0])/$Wout)**2 , -1);
                        }
			elsif ( $geneSet{$gene}{gc} > $valGC[1] ) {
				$GCpenalty= max(-(($geneSet{$gene}{gc}-$valGC[1])/$Wout)**2 , -1);
			}
			else { print "Seems GC value provided in wrong format"; exit;}

			my $GCScore = $GCpenalty; # normal weigh for GC% criteria



			# COVERAGE SCORING

			# adding a condition to avoid computing a score if option is disabled (i.e. coverage=0)
			my $covPenalty=0;
			if ($param_ref->{coverage} != 0) {
				my $cLow=$param_ref->{coverage}/2;      # if cov = 70, $cLow = 35				# 125: 62.5
				my $cUp=$param_ref->{coverage}+$cLow;   # if cov = 70, $cUp = 105				# 125: 187.5
				my $cWin=($cUp-$cLow)/5;                # if cov = 70, $cWin = 21				# 125: 25     (in)
				my $cWout=($cUp-$cLow)/5;               # if cov = 70, $cWout = 35 (equivalent to = $cLow)	# 125: 25   (out)

		                if (($geneSet{$gene}{readcov} >= $cLow) && ( $geneSet{$gene}{readcov} <= $cUp)) {
					# first formula for positive score => from 0 to +1
		                        #$covPenalty=tanh(3*($geneSet{$gene}{readcov}-$cLow)/$cWin)-tanh(3*($geneSet{$gene}{readcov}-$cUp)/$cWin)-1;
					# second formula for positive score => from 0 to +1  (in order to reduce the tendency for true HGT (i.e. alien) to become self)				
					#$covPenalty=(tanh(3*($geneSet{$gene}{readcov}-$cLow)/$cWin)-tanh(3*($geneSet{$gene}{readcov}-$cUp)/$cWin)-1)/10;
					# third formula for positive score => strictly equal 0 (a good coverage can not help decipher between alien and self)
					$covPenalty=0
		                }
		                elsif ( $geneSet{$gene}{readcov} < $cLow ) {
		                        $covPenalty= max(-(($geneSet{$gene}{readcov}-$cLow)/$cWout)**2 , -1);
		                }
		                elsif ( $geneSet{$gene}{readcov} > $cUp ) {
		                        $covPenalty= max(-(($geneSet{$gene}{readcov}-$cUp)/$cWout)**2 , -1);
		                }
		                else { print "Seems something wrong with coverage value"; exit;}
			}

			# EXPRESSION LEVEL SCORING (Kallisto)

			# adding a condition to avoid computing a score if option is disabled (i.e. rna_seq=0)
			my $expressionScore=0;
			if ($param_ref->{rna_seq} > 0) {
				#print $geneSet{'Chrom_1:2155-5638'}{geneexp}; print "\n";
		                if ($geneSet{$gene}{geneexp} == 0) {
					$expressionScore= -1;
		                }
		                elsif (( $geneSet{$gene}{geneexp} > 0 ) and ($geneSet{$gene}{geneexp} < 0.2 )) {
		                        $expressionScore= 0;
		                }
		                elsif (( $geneSet{$gene}{geneexp} >= 0.2 ) and ($geneSet{$gene}{geneexp} < 1 )) {
		                        $expressionScore= 0.5;
		                }
		                elsif ( $geneSet{$gene}{geneexp} >= 1 ) {
		                        $expressionScore= 1;
		                }
		                else { print "Seems something wrong with TPM value"; exit;}
			}


			# finalyzing crriteria scores
			my ($rRNAVal, $buscoVal, $refVal, $speciesName, $orderName, $phylumName, $skName, $feederVal, $acceptVal, $norRNA, $norBUSCO, $norREF, $norGENO);

			# here be carefull of what column is used if blast outputs are modified ine the future...
			if ($genomeHash{$gene}->[12]) { $speciesName=$genomeHash{$gene}->[12];} else {$speciesName="NA"; }
			#if ($genomeHash{$gene}->[13]) { $orderName=$genomeHash{$gene}->[13];} else {$orderName="NA"; }
			#if ($genomeHash{$gene}->[14]) { $phylumName=$genomeHash{$gene}->[14];} else {$phylumName="NA"; }
			#if ($genomeHash{$gene}->[15]) { $skName=$genomeHash{$gene}->[15];} else {$skName="NA"; }

			# Normalization of "blast-based" BITSCORES
			# (check for the exact field used, it could depend on the blast ouput format...)
			if ($rRNAHash{$gene}->[10]) { $norRNA = normalize($rRNAHash{$gene}->[10], $rRNAmin, $rRNAmax); } else {$norRNA=0;}
			if ($buscoHash{$gene}->[10]) { $norBUSCO = normalize($buscoHash{$gene}->[10], $BUSCOmin, $BUSCOmax);} else { $norBUSCO = 0; }
			if ($refHash{$gene}->[10]) { $norREF = normalize ($refHash{$gene}->[10], $REFmin, $REFmax);} else { $norREF = 0;}
			if ($genomeHash{$gene}->[11]) { $norGENO = normalize ($genomeHash{$gene}->[9], $Gmin, $Gmax); } else {$norGENO= 0; }

			#### "Alien" versus "Self" scoring :

			# managing rRNA score (alien)
			if ($rRNAHash{$gene}->[10]) { $rRNAVal=($norRNA * (-1)); } else { $rRNAVal=0; } 
			# managing feeder score (alien)
			if ($speciesName eq $param_ref->{feedlevel_species}) { $feederVal=($norGENO * (-1)); } else { $feederVal=0; }

			# managing busco score (self)
			if ($buscoHash{$gene}->[10]) { $buscoVal=$norBUSCO; } else {$buscoVal=0;}
			# managing reference score (self)
			if ($refHash{$gene}->[10]) { $refVal=$norREF; } else {$refVal=0;}

			# taxonomic decision for general blast score (important step)
			print $TAXDETECTION "GENERAL BLAST parseTax:\t$gene\t$speciesName\t$param_ref->{level_upto}\t";
			my $taxDecision=parseTax($speciesName, $param_ref->{level_upto});
			# managing general bast score (self versus alien)
			if ($taxDecision) {
				$acceptVal = $norGENO; 				# positive score
				print $TAXDETECTION "\t=> SELF\n";
			} 	
			else { 
				$acceptVal=($norGENO * (-1));			# negative score
				$alienCnt++ if not exists $tnfHash{$gene};      # WARNING: For 'NA' it counts alien !!! OVERALL, THE alienCnt VARIABLE SEEMS USELESS IN THE PIPELINE
				print $TAXDETECTION "\t=> ALIEN\n";
			}
			#my $GeneralBlastVal = $acceptVal*2; 	# giving more weigh to general blast criteria
			my $GeneralBlastVal = $acceptVal; 	# giving no additional weigh to general blast criteria

			# COMPUTING THE INFORMATIVENESS OF CRITERIAS
			my $info_threshold = 5;	# expected number of informative criteria to not reduce the absolute score value
			my $informativeness = 2;	# because GC and cov are already informative, even when their value is "0"
			# here "" need to remove, it is used for string not for number 
			if ($refVal != 0) { $informativeness += 1; }
			if ($buscoVal != 0) { $informativeness += 1; }
			if ($rRNAVal != 0) { $informativeness += 1; }
			if ($feederVal != 0) { $informativeness += 1; }
			if ($GeneralBlastVal != 0) { $informativeness += 1; }
			if ($min_dist_score != 0) { $informativeness += 1; }
			#if ($GCScore != "0") { $informativeness += 1; }		# GC il always informative
			#if ($covPenalty != "0") { $informativeness += 1; }		# coverage il always informative
		
			######################
			# COMPUTING GENE SCORE  (set A)
			my $geneScoreA = ($refVal + $buscoVal + $rRNAVal + $feederVal + $GeneralBlastVal + $min_dist_score + $GCScore + $covPenalty);
			my $geneScoreA_modified=$geneScoreA*$informativeness/$info_threshold;	# => genescore absolute value scales with number of informative criterias
			my $geneScoreA_final=AOM::Module::tanhGeneScore($geneScoreA_modified);   	# => gene score from -1 to +1

			my $logLen;
			if ($tnf==0) { $logLen=log10($geneSet{$gene}{len}); } else { $logLen=log10($tnfLen); }
			push @allGenesScoreA, ($geneScoreA_final * $logLen);
			push @allGenesLength, $logLen;
			my $tmpVal=$geneScoreA_final * $logLen;
			my $tmpGC=$GCScore;

			######################
			# COMPUTING GENE SCORE (set B = with expression scores)
			#$roughScore= "$refVal, $buscoVal, $rRNAVal, $feederVal, $acceptVal, $min_dist_score, $GCScore, $expressionScore";
			my $geneScoreB = ($refVal + $buscoVal + $rRNAVal + $feederVal + $GeneralBlastVal + $min_dist_score + $GCScore + $covPenalty + $expressionScore);
			my $geneScoreB_modified=$geneScoreB*($informativeness+1)/$info_threshold;	# => genescore absolute value scales with number of informative criterias
			my $geneScoreB_final=AOM::Module::tanhGeneScore($geneScoreB_modified);   		# => gene score from -1 to +1

			push @allGenesScoreB, ($geneScoreB_final * $logLen);
			my $tmpValB=$geneScoreB_final * $logLen;

			######################
			### OUTPUTS GENE SCORES DETAILS:
			print $SCORE "$accession_number\t$genome{$accession_number}{len}\t$genome{$accession_number}{gc}\t$gene\t$refVal\t$GeneralBlastVal\t$buscoVal\t$rRNAVal\t$feederVal\t$min_dist_score\t$geneSet{$gene}{readcov}\t$covPenalty\t$geneSet{$gene}{gc}\t$GCScore\t$geneScoreA\t$geneScoreA_final\t$geneSet{$gene}{geneexp}\t$expressionScore\t$geneScoreB_final\n";
		}


		# SCAFFOLD SCORING

		my $gSum=AOM::Module::sumArray(@allGenesScoreB); # using gene scores "set B" (with expression level)
		my $lSum=AOM::Module::sumArray(@allGenesLength);
		#print "$gSum\t$lSum ***\n";
		#no warnings 'recursion';
		if ((AOM::Module::sumArray(@allGenesLength)) == 0) {
			push @allGenesLength, 1 ;
		}
		my $rawscaffoldScore=(AOM::Module::sumArray(@allGenesScoreB))*sqrt(scalar(@allGenesScoreB))/(AOM::Module::sumArray(@allGenesLength)); 
			# => results range from -scalar(@allGenesScoreB) to +scalar(@allGenesScoreB)
			# next step: tanh(y*scaffscore/x)  # y=3? ; x=1 ; 
		my $finalscaffoldScore=tanh(2*$rawscaffoldScore);


		##
		# All Hits detail
		foreach my $nam (keys %namesCnt) { push @allNames, "$nam:$namesCnt{$nam}";}
		my $maxVal;
		if (%namesCnt) { $maxVal = max values %namesCnt;} else {$maxVal=0;}
		my $allNames = 'NA';
		$allNames = join ',', @allNames;
		#my $sequence=$genome{$accession_number}{nuc_seq};

		#Deal with TNF
		foreach my $tnfkey (keys %tnfSps) { push @allTNF, "$tnfkey:$tnfSps{$tnfkey}";}
		my $allTNF = 'NA';
		$allTNF = join ',', @allTNF;
		if (!@allTNF){$allTNF='NA'}

		#combine spsCount and tnf species
		my %newSpsHash = (%spsCounts, %tnfSps);
		#contig/chr classification name
		my $spsMax = max values %newSpsHash;
		my $spsClass='Not_Known';
		if (defined $spsMax) {
			my @spsKeys = grep { $newSpsHash{$_} == $spsMax } keys %newSpsHash; # extract the key using val
			#print "$spsKeys[0]\t$spsKeys[1]\t-<->-\n";
			if (scalar (@spsKeys) > 1) { print $LOG "You have multiple hits with similar count for this contig : $accession_number\n"; }
			$spsClass = $spsKeys[0];
		}
		else { print $LOG "Interesting, no classification species for this : $accession_number\n"; }
		#AOM::Module::pLine("*");


		### IF NO GENE ON SCAFFOLD:
		# 1- check TNF
		# 2- check GC
		# 3- check coverage (NOT IMPLEMENTED YET)
                # If ActivateTNF is set 'yes' then it will
		if (($spsClass eq "Not_Known") and ( $param_ref->{ActivateTNF} eq 'on')) { #Remember the case sensitivity 
			print $LOG "\t$spsClass = Not_Known:\t$accession_number\t";

			#print "Dealing with TNF prediction for no-gene contigs/scaffold\n";
			my $tmpSeq=$genome{$accession_number}{nuc_seq};
			my $scaff_fh=AOM::Module::write_fh("$param_ref->{project_dir_path}/intermediate_files/tnf/$accession_number"); #Write the fasta sequence with accession number file
			print $scaff_fh ">$accession_number\n$tmpSeq\n";
			close $scaff_fh;
			my %tnfHash1;
			my $scaffCoord_ref = AOM::tnfHandler::calculateTNF("$param_ref->{project_dir_path}/intermediate_files/tnf/$accession_number", 0, 4);
			foreach my $freq1 (@$scaffCoord_ref) {
				my @tnf1=split /\t/, $freq1;
				my $gName1 = shift @tnf1; #print "$gName1 <<********>>";
				my ($name1, $nearest1, $min_dist1) = AOM::tnfHandler::SearchTNF(\@tnf1, "$param_ref->{tnfDB_file}", "$param_ref->{project_dir_path}/intermediate_files/tnf/tnfNoHits.tnf", $gName1);
				my @lenCor1 = split /\-/, (split /\:/, $gName1)[-1]; #Get the sequence len
				my $len1 = length ($genome{$accession_number}{nuc_seq});
				#my @nearestName1=split /\_/, $nearest1;
				#$tnfHash1{$name1}="$nearestName1[0]-$min_dist1-$len1";
				$tnfCnt++;

				my $newDist=1/(1+$min_dist1); # 1-(1/(1+mindistance))
				my $min_dist_score1=$newDist; 

				print $TAXDETECTION "Secondary TNF parseTax:\t$accession_number\t$param_ref->{level_upto}\t";
				my $tnfTaxDecision1=parseTax($nearest1, $param_ref->{level_upto});
				if ($tnfTaxDecision1) { $min_dist_score1 = ($min_dist_score1 * 1);
					#print $TAXDETECTION "\t=> SELF\n";
				}
				else { $min_dist_score1 = ($min_dist_score1 * -1); $alienCnt++;
					#print $TAXDETECTION "\t=> ALIEN\n";
				}

				# SCAFFOLD GC SCORING
				$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;

				my $scaffGC=$genome{$accession_number}{gc};
				my $scaffabsDist = find_dist($scaffGC , $scaffGC);
				my $scaffwidthDiv = 6.5*$scaffGC/10;
				#my $scaffGCpenalty = 2*(exp(-(($scaffabsDist/$scaffwidthDiv)**4)-0.5));    # what is this ?! (it is cancelled 3 lines later...)

				#print "\n# Scaffold GC scoring: $accession_number ($date)\n";
				#my $W=$valGC[1]-$valGC[0];
				my @valGC = split /\:/, $param_ref->{gc_filter}; #Provided in config file
				my $scaffGCpenalty=0;
				my $Win=3;
				my $Wout=10;

				if (($valGC[0] <= $genome{$accession_number}{gc}) && ($valGC[1] >= $genome{$accession_number}{gc})) {
				       $scaffGCpenalty=tanh(1*($genome{$accession_number}{gc}-$valGC[0])/$Win)-tanh(1*($genome{$accession_number}{gc}-$valGC[1])/$Win)-1;
				}
				elsif ($valGC[0] >= $genome{$accession_number}{gc}) {
				       $scaffGCpenalty= max(-(($genome{$accession_number}{gc}-$valGC[0])/$Wout)**2 , -1);
				}
				elsif ($valGC[1] <= $genome{$accession_number}{gc}) {
				       $scaffGCpenalty= max(-(($genome{$accession_number}{gc}-$valGC[1])/$Wout)**2 , -1);
				}
				else { print "Problem with GC value: ".$genome{$accession_number}; exit;}

				$allTNF=$nearest1;
				my $geneClass='NA';

				
				# SCAFFOLD COVERAGE SCORING
        	                my $covPenalty=0;
	                        # adding a condition to avoid computing a score if option is disabled (i.e. coverage=0)
                	        if ($param_ref->{coverage} != 0) {
                        	        my $cLow=$param_ref->{coverage}/2;      # if cov = 70, $cLow = 35                               # 125: 62.5
                               		my $cUp=$param_ref->{coverage}+$cLow;   # if cov = 70, $cUp = 105                               # 125: 187.5
                                	my $cWin=($cUp-$cLow)/5;                # if cov = 70, $cWin = 21                               # 125: 25     (in)
                                	my $cWout=($cUp-$cLow)/5;               # if cov = 70, $cWout = 35 (equivalent to = $cLow)      # 125: 25   (out)
					my $sLen=$genome{$accession_number}{len};

					my $accession_number_updated="$accession_number:0-$sLen"; # The accession will work file, no need tp provide st and ed
					my $scaffCov = 0; 
					$scaffCov = `samtools depth -aa -r $accession_number_updated $param_ref->{bam} | awk '{ c++; s+=\$3; } END {k=(s/c); print k;}'`;
					chomp $scaffCov;
					print "$accession_number_updated ---- $scaffCov *'''''*\n";

                                	if (($scaffCov >= $cLow) && ( $scaffCov <= $cUp)) {
                                        	$covPenalty=tanh(3*($scaffCov-$cLow)/$cWin)-tanh(3*($scaffCov-$cUp)/$cWin)-1;
                                	}
                                	elsif ( $scaffCov < $cLow ) {
                                        	$covPenalty= max(-(($scaffCov-$cLow)/$cWout)**2 , -1);
                                	}
                                	elsif ( $scaffCov > $cUp ) {
                                        	$covPenalty= max(-(($scaffCov-$cUp)/$cWout)**2 , -1);
                                	}
                               		else { print "Seems something wrong with coverage value"; exit;}
                        	}
				# SCAFFOLD GENERAL SCORE (if no gene)
				#my $scaffoldScore=($scaffGCpenalty+$min_dist_score1+$covPenalty);
				my $rawscaffoldScore=($scaffGCpenalty+$covPenalty); 
				$spsClass=$nearest1;
				my $finalscaffoldScore=tanh(2*($rawscaffoldScore));

				print $LOG "\n\tno-gene scaffold ($accession_number) score:\t$finalscaffoldScore\t(GCscore=$scaffGCpenalty ; Covscore=$covPenalty ; TNFscore=$min_dist_score1)\n";

			}
		}
		### OUTPUT SCAFFOLD RESULTS (various characteristics):
		print $RESULT "$accession_number\t$genome{$accession_number}{len}\t$genome{$accession_number}{gc}\t$totalGenes\t$refBlast\t$genomeBlast\t$maxVal,$allNames\t$rRNA\t$busco\t$tnfCnt\t$allTNF\t$alienCnt\t$rawscaffoldScore\t$finalscaffoldScore\n";

		undef @allGenesScoreA; undef @allGenesScoreB; undef @allGenesLength; undef %newSpsHash;

		# print genomeHash for troobleshooting
		print "ENTIRE GENEHASH:\ngene\tgenomeHash{gene}\n";
		foreach my $gene (keys %genomeHash) {print "$gene\t$genomeHash{$gene}\n";}

	}
}


##############################################################################
#-----------------------------------------
elsif ($param_ref->{mode} eq 'random') {
	print "Thanks for your interest, unfortunatly 'random' is not functional at the moment\n";  exit;
	# Lets check for alien in all the sequence now
	foreach my $key (keys %genome) {
		my $accession_number=$key; # Consider the sequence name here
		#AOM::Module::pLine(".");
		my $sequence=$genome{$key}{nuc_seq};
		#AOM::Module::randomBased($accession_number,$sequence, $key, $param_ref, $taxidP_ref, $taxidT_ref, $taxidN_ref);
	}
}

else { print "Did you forgot to provide 'mode' values:gene or random \n Check your conguration file\n TERMINATING\n"; exit; }


#-------------------------------------------------------------------------------
#Plot the cluster
#AOM::plotHandler::KmeanPlot("$param_ref->{project_dir_path}/results/$param_ref->{result}");
#Use the result/final_results.out file and create a new final_results.out_clustered.csv file with one nice plot
#system ("Rscript utils/plotKmean.R $param_ref->{project_dir_path}/results/$param_ref->{result} $param_ref->{project_dir_path}/results");
#-------------------------------------------------------------------------------

##############################################################################
# Using genescores and scaffold score to discriminate HGTs from contaminant

$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
print "\n# Categoryzing genes (self / HGT / contaminant)  ($date)\n";

# limiting the number of false negatives (i.e. undetected HGT)
my $gene_low_th=-0.2;
my $gene_up_th=0.2;
my $scaf_low_th=-0.2;
my $scaff_up_th=0.2;

open (my $GFF_OUTPUT, ">", "$param_ref->{project_dir_path}/results/HGT.$scaf_low_th\_$gene_low_th\_$scaff_up_th\_$gene_up_th.gff") || die ('Could not create result file. Please check writing permission in your current directory', "\n");
AOM::HGTcategorization::categorize_genes("$param_ref->{project_dir_path}/results/$param_ref->{result}", "$param_ref->{project_dir_path}/results/$param_ref->{score}", $scaf_low_th, $gene_low_th, $scaff_up_th, $gene_up_th, $GFF_OUTPUT);

# balancing false positives and negatives
$gene_low_th=-0.5;
$gene_up_th=0.5;
$scaf_low_th=-0.5;
$scaff_up_th=0.5;

open (my $GFF_OUTPUT2, ">", "$param_ref->{project_dir_path}/results/HGT.$scaf_low_th\_$gene_low_th\_$scaff_up_th\_$gene_up_th.gff") || die ('Could not create result file. Please check writing permission in your current directory', "\n");
AOM::HGTcategorization::categorize_genes("$param_ref->{project_dir_path}/results/$param_ref->{result}", "$param_ref->{project_dir_path}/results/$param_ref->{score}", $scaf_low_th, $gene_low_th, $scaff_up_th, $gene_up_th, $GFF_OUTPUT2);

# limiting the number of false positives (i.e. detection of false HGT)
$gene_low_th=-0.8;
$gene_up_th=0.8;
$scaf_low_th=-0.8;
$scaff_up_th=0.8;

open (my $GFF_OUTPUT3, ">", "$param_ref->{project_dir_path}/results/HGT.$scaf_low_th\_$gene_low_th\_$scaff_up_th\_$gene_up_th.gff") || die ('Could not create result file. Please check writing permission in your current directory', "\n");
AOM::HGTcategorization::categorize_genes("$param_ref->{project_dir_path}/results/$param_ref->{result}", "$param_ref->{project_dir_path}/results/$param_ref->{score}", $scaf_low_th, $gene_low_th, $scaff_up_th, $gene_up_th, $GFF_OUTPUT3);


### Finishing Alienomics run

$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
print "\n# Alienomics run is finished ($date)\n" if $param_ref->{verbose};

close($SUMMARY);
close($RESULT);
close($LOG_ERR);
close($LOG);


###########################################################################################
# sub-routines
###########################################################################################

# Extract the subset of GFF
sub subsetGFF {
my ($genes, $chrName, $len)= @_;
my @genelist=keys %$genes;
my @genesLoc; my %intronStatus;
foreach my $idlist(@genelist) {
	chomp $idlist;
	next if ("$genes->{$idlist}{seqid}" ne "$chrName");
	my $gene=$genes->{$idlist};
	if ($gene) {
		# START base shifted due to ONE based system GFF and ZERO in BED :(
		my $st = $genes->{$idlist}{start};
		my $ed = $genes->{$idlist}{end};
		push @genesLoc, "$genes->{$idlist}{seqid}".":"."$st"."-"."$ed";
		#if ($genes->{$idlist}{intron}) { $intronStatus{$genes->{$idlist}{seqid}}=1; } else {  $intronStatus{$genes->{$idlist}{seqid}}=0; }
	#print "$genes->{$idlist}{seqid}\t$genes->{$idlist}{start}\t$genes->{$idlist}{end}\t$len ---\n";
	#print "$genes->{$idlist}{seqid}\t$genes->{$idlist}{intron} ++++---- \n";
	}
}
return \@genesLoc;
}

# Log base 10 subs
sub log10_10 {
   my $n = shift;
   return log($n)/log(10);
}


########
sub storeInHash { 
my $inFile=shift;
my %gHash; 
open (GL, "$inFile") || die "File not found\n";
     while (<GL>) {
	chomp;
         my @tmpval= split(/\t/, $_);
	#BED and gene.fa header have one basepair differenes !!!!!!!!!! This is due to counting coordinates for O or 1
        #I am FORCE to change the coords here by one base -- it kiil 2 hour of my life :(
	my @nameVal= split(/\:/, $tmpval[0]);
        my @coordVal=split(/\-/, $nameVal[1]);
	my $newCoord1=$coordVal[0]+1; my $newCoord2=$coordVal[1]+0;
	my $updatedName="$nameVal[0]:$newCoord1-$newCoord2";
          $gHash{$updatedName} = $tmpval[1];
     }
close(GL);
return \%gHash;
}

########
sub storeInHash_trans { 
my $inFile=shift;
my %transHash;
open (TL, "$inFile") || die "Transcriptomics/Kallist file not found\n";
     while (<TL>) {
        chomp;
        my @tmpval= split(/\t/, $_); 
        #Input format
        #target_id	length	eff_length	est_counts	tpm
        #Lactococcus_lactis:20916-21373	1924	1746.98	102.328	11129.2 
          $transHash{$tmpval[0]} = $tmpval[4];
        #Stored Lactococcus_lactis:20916-21373	11129.2
     }
close(TL);
return \%transHash;
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
($min, $max) = minmax (@allScore);
if (!$min) { $min = 0; }
if (!$max) { $max = 0; }
return ($min, $max);
}

############################################################
# normalize criteria score
sub normalize {
my ($value, $min, $max)=@_;
my $tanhVal = tanh ( 3 * ($value/$param_ref->{bitscoreCutoff}) )**5 ;
return $tanhVal; #Return tanh value
}

############################################################
sub find_dist {
    my ($x, $v) = @_;
    return abs($x - $v);
}

############################################################
#Convert GFF to BED
sub gffGene2Bed {
	my ($genes, $param, $transHash_ref)= @_;
	# Assume all gene names are unique
	my @genelist=keys %$genes;
	my $bedFH=AOM::Module::write_fh("$param_ref->{project_dir_path}/intermediate_files/bed/gene.bed", 1);
	foreach my $nam(@genelist) {
		chomp $nam;
		my $gene=$genes->{$nam};
		if ($gene) {
			AOM::gffHandler::printGeneBed($gene, $bedFH);
		}
	}
	close $bedFH;

	# Extract all genes using bed coordinates
	$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
	print "\n# Extract genes using bed coordinates ($date)\n"; 

	system ("$param_ref->{seqtk_path} subseq $param_ref->{project_dir_path}/intermediate_files/gff/genome.fa $param_ref->{project_dir_path}/intermediate_files/bed/gene.bed > $param_ref->{project_dir_path}/intermediate_files/bed/gene_raw.fa");
	#Sed all non-AGTC charater
	system ("sed '/^[^>]/ s/[^AGTC]/N/gi' < $param_ref->{project_dir_path}/intermediate_files/bed/gene_raw.fa >$param_ref->{project_dir_path}/intermediate_files/bed/gene.fa");
	#Remove space in header
	#system ("sed '/^>/{s/\s/_/g}' <$param_ref->{project_dir_path}/intermediate_files/bed/gene_raw.fa > $param_ref->{project_dir_path}/intermediate_files/bed/gene.fa");  # \s will also target tabs
	
	
	#Lets store the coverage information in Hash
	$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
	my $covHash_ref;

	if ($param_ref->{coverage} >0 ) {
	print "\n# Compute gene coverage ($date)\n";
	AOM::covHandler::extractMappingCov($param_ref->{bam},"$param_ref->{project_dir_path}/intermediate_files/bed/gene.bed", $param_ref->{max_processors},"$param_ref->{project_dir_path}/intermediate_files/bed/gene.cov");
	#$geneCov = `samtools depth -aa -r $accession_number $param->{bam} | awk '{ c++; s+=\$3; } END {k=(s/c); print k;}'`;
	$covHash_ref=storeInHash("$param_ref->{project_dir_path}/intermediate_files/bed/gene.cov");
	#print "$covHash_ref->{'scaffold_1_Avaga:7035-8356'} -<<<<<<>>>>>>>>>>>**<<<<\n"; #This is not work now, i updates one+ in storeInHash subs
	}
	else {print "\nCoverage criteria will not be used (set to \"0\" in the config file) ($date)\n";}
        #Add transcriptomics data

	#Lets store the transcriptomics information in Hash
	#$date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;

#	print "\n# Compute kallisto score ($date)\n";

#	my $transHash_ref;
#        if ($param_ref->{rna_seq} >0 ) {
#        print "\n# Storing Transcriptomics Data ($date)\n";
#        $transHash_ref=storeInHash_trans("$param_ref->{kallisto_file}");
#        }
        print "\n# Store all score ($date)\n";
	my $geneSet = AOM::Module::parse_gene_files("$param_ref->{project_dir_path}/intermediate_files/bed/gene.fa", $param, $covHash_ref, $transHash_ref);
	return $geneSet;
}

############################################################
# Extract all the IDs from BED ... Point to remeber ZERO based ONE based system
sub extractIds {
my ($bedFile)= @_;
my @allVal;
open(my $fh, '<:encoding(UTF-8)', $bedFile);
while (<$fh>) {
	chomp;
	my @vals=  split /\t/, $_;
	push @allVal, $vals[0];
}
my @allIds = uniq (@allVal);
foreach my $v (@allIds) {print "$v\n";}
}

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}


sub mean { return @_ ? sum(@_) / @_ : 0 }


############################################################
# Parse Taxon name
sub parseTax {
	# print some info for troobleshooting
	print $TAXDETECTION "\n(Starting parseTax subroutine from Alienomics_v2.pl)\n";

	my ($sps,$groupName)= @_;

	# print some info for troobleshooting
	#print $TAXDETECTION "(sps = $sps)\n";
	#print $TAXDETECTION "(groupName = $groupName)\n";

	# SHOULD'T THIS BE "-taxonid" INSTEAD OF "-name" ? (depends on what field is used (the 11th or the 12th) )
	my $nam = $dbh->get_taxon(-name => $sps); 

	# print some info for troobleshooting
	#print $TAXDETECTION "\tnam = $nam\n";

	if(!$nam) {
		return 0;  	# warning here: "return 0" means "alien". But this is only an unknown origin, which is not necessarily "alien".
			 	# We should have a third possibility so that the blast scores become "0" the "parseTax" subroutine fails.
		print $TAXDETECTION "$sps\tNo taxon name !\n";
	}
	elsif($nam eq "NA") {
		print $TAXDETECTION "$sps\tNA\n";
	}	
	else {print $TAXDETECTION "$sps\tNA\n";}
	my $tree_functions = Bio::Tree::Tree->new();
	my $lineage = lc ($tree_functions->get_lineage_string($nam));

	# print some info for troobleshooting
	print $TAXDETECTION "$sps\t$nam\t$lineage\t($groupName)\n";

	if (index(lc($lineage), lc($groupName)) != -1) {
		#print $TAXDETECTION "$sps\tSELF\n";	# for troobleshooting
		return 1;
	}
	elsif (index(lc($lineage), lc($groupName)) == -1) {
		#print $TAXDETECTION "$sps\tALIEN\n";	# for troobleshooting
		return 0;
	}
	else {print $TAXDETECTION "$sps\tunclassified taxa\n";}	# for troobleshooting
}

sub process_file{ # This is your custom subroutine to perform on each file    
    my $f = shift;
    my ($val, $nam) = check_file($f);   
	if ($val == 1 and ($nam eq "names")) {  # print "processing file $f\n";
		%txn_names=file2hash ($f, $nam);
		#return @array;
		}
         elsif ($val == 1 and $nam eq "nodes") { # print "processing file $f\n";
		%txn_nodes=file2hash ($f, $nam);
		}
	}

#-----------------------------------------------------------------------
sub check_file {
    use File::Basename;
    my $filepath = shift; # print $file;
    my $file = basename($filepath);
    my @ff= split /\./, $file;
    if ($ff[0] eq "names" || "nodes" ) 
	{ return 1, $ff[0]; }
}

#-------------------------------------------------------------------------
sub file2hash {
	my ($infile, $n) = @_;
	my %hash;
	open FILE, $infile or die $!;
	while (<FILE>) {
   		chomp;  # s/^\s*(.*)\s*$/$1/;
	        next if (index($_, "scientific") == -1); #{ print "'$string' contains '$substring'\n";}
		my @tmp_array= split /\t/ , $_;
		s{^\s+|\s+$}{}g foreach @tmp_array; # Removing leading and trailing whitespace from array strings.
   		my ($key, $val) = split /\t\|\t/;
   		#Lets add only scientific name -- rest ignore
		if (($infile eq "names.dmp") and  ($tmp_array[6] ne "scientific name")) {next;}  #print "$tmp_array[6]\n";      ## if we want to enter only specific lines.
   		#I edited this line, do not remember why i did add line nu;ber in it
		if($n eq "names") { $key="$key";}    # I make it unique by adding the line number and split later...
		$val =~ s/\s+/_/g;  ## replace the space with underscore ...
		$hash{$key} = $val;  
		# print "$n\t$key\t$val\n";
        } 
	close FILE;
return %hash;      
}

#-------------------------------------------------------------------------
sub printhash {
  	my %hash=%{$_[0]};
	foreach my $key (sort keys %hash) {
     	print "$key : $hash{$key}\n";
	}
}

#-------------------------------------------------------------------------
sub findkey {
        my ($species, $hash) =@_;
	my  %hash=%$hash;  my @all_keys;
	foreach my $key (keys %hash) {
     	if ($hash{$key} =~ m/^$species$/i) { push @all_keys, $key};
	}
s{^\s+|\s+$}{}g foreach @all_keys; # Removing leading and trailing whitespace from array strings.
return @all_keys;
undef @all_keys;
} 

#-------------------------------------------------------------------------
sub findId {
        my ($id, $hash) =@_;
	my  %hash=%$hash;  my @all_values;
	foreach my $key (keys %hash) {
	my @newValue= split(/:/, $hash{$key});       ## What if we have more than two hits for a key !!!!!
     	if ($newValue[0]==$id) { push @all_values, $key};
	}
s{^\s+|\s+$}{}g foreach @all_values; # Removing leading and trailing whitespace from array strings.
my $all_values=join(",",uniq(@all_values));
return $all_values;
undef @all_values;
}


#----------------------------------------------------------------------
sub getDirectoryFiles {          # It get the directory files and return it
     my $taxdir = shift;

     opendir(my $dh, $taxdir) || die "can't opendir $taxdir : $!";
     my @entries = grep {!( /^\.$/ || /^\.\.$/)} readdir($dh);
     @entries =  map { "$taxdir/$_" } @entries; #change to absolute paths
     closedir $dh;

     my @files =  grep( -f $_ , @entries);
     my @dirs = grep(-d $_, @entries);
     return (\@files,\@dirs);     ## return as a reference 

}

############################################################
### POD Documentation

__END__

=head1 NAME

Alienomics.pl	- Script to automated identify contamination and do sensitive filter.

=head1 SYNOPSIS

Usage: alienomics.pl -c or --conf [path to configuration file]\n\n";

Parameters:

  -c :To execute alienomics with the parameters defined in a configuration file
  -v :Alienomics version
  -cc :To create a new configuration file to use with alienomics

Try -h for more detail.

=head1 OPTIONS
