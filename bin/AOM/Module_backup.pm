# $Id$ Module for Module
# Perl module for ALIENOMICS AOM::Module
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

AOM::Module  - DESCRIPTION of Object

=head1 SYNOPSIS
Give standard usage here
=head1 DESCRIPTION
Describe the object here
=cut

=head1 CONTACT
Jitendra <jnarayan81@gmail.com>
=head1 APPENDIX
The rest of the documentation details each of the object methods.

=cut
##-------------------------------------------------------------------------##
## Let the code begin...
##-------------------------------------------------------------------------##
package AOM::Module;

use strict;
use warnings;
# -- comment to force avoid this Variable "%taxidP" will not stay shared at /home/jit/Downloads/Alienomics_v1.archives-2020-06-01/bin/./AOM/Module.pm line 748.
use Bio::SeqIO;
use Cwd;
use File::chdir;
use File::Copy;
use POSIX;
use File::Temp qw(tempfile);
use Statistics::Distributions qw(chisqrprob);
use Statistics::Multtest qw(BY qvalue);
use Tie::File;
use Try::Tiny;
use Data::Dumper;
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use FindBin;
use File::Remove;
use File::Path qw(make_path remove_tree);
use Capture::Tiny ':all';
use Getopt::Long;
#use Statistics::R;
use Math::Round;
use File::Find;
use Bio::DB::Taxonomy;
use Pod::Usage;
use Bio::Root::Root;
use Exporter;
our @EXPORT_OK = "AOMule";

#Turning off BioPerl warnings
$Bio::Root::Root::DEBUG = -1;

sub help {
  my $ver = $_[0];
  print "\n  alienomics $ver\n\n";

  print "Usage: alienomics.pl -c [path to configuration file]\n\n";
  print "       To execute alienomics with the parameters defined in a configuration file\n\n";
  print "                OR\n\n";
  print "       alienomics.pl -cc [additional parameters]\n\n";
  print "       To create a new configuration file to use with alienomics\n\n";
  print "                OR\n\n";
  print "       alienomics.pl  -v\n\n       To print alienomics's version\n\n\n\n";
  print "Currently supported options when creating a new configuration file are:\n\n";
}


############################################################################################
# Read config files in the form element = value #comment --->
sub read_config {
  my $project_config_file = shift;
  my $alien_path = shift;
  my %param;
  open(my $fh_ParamConfig, "<", "$$project_config_file") || die ("Couldn't open the project configuration file: $!\n");

# BASIC PARAMETER FOR LOCATION AND FILE NAME --->
  $param{alien_dir} = read_config_file_line('alien_dir', '', $fh_ParamConfig);

# INPUT FILES --->
  $param{Input_genome} = read_config_file_line('Input_genome', '', $fh_ParamConfig);
  $param{AUGUSTUS_CONFIG_PATH} = read_config_file_line('AUGUSTUS_CONFIG_PATH', $param{alien_dir}, $fh_ParamConfig);

# BLASTDB location --->
  $param{blastdb_path} = read_config_file_line('blastdb_path', $param{alien_dir}, $fh_ParamConfig);
  $param{diamonddb_path} = read_config_file_line('diamonddb_path', $param{alien_dir}, $fh_ParamConfig);

# TAXDUMP location --->
  $param{taxdump_path} = read_config_file_line('taxdump_path', $param{alien_dir}, $fh_ParamConfig);
  $param{accession2taxid_path} = read_config_file_line('accession2taxid_path', $param{alien_dir}, $fh_ParamConfig);

# PROJECT NAME --->
  $param{project_dir_path} = read_config_file_line('project_dir_path', $param{alien_dir}, $fh_ParamConfig);

# PROJECT CONFIGURATION --->
  $param{reference_genome_file} = read_config_file_line('reference_genome_file', '', $fh_ParamConfig);
  $param{reference_proteome_file} = read_config_file_line('reference_proteome_file', '', $fh_ParamConfig);

  $param{coverage} = read_config_file_line('coverage', '', $fh_ParamConfig);
  $param{bam} = read_config_file_line('bam', '', $fh_ParamConfig);
  $param{rna_seq} = read_config_file_line('rna_seq', '', $fh_ParamConfig);
  $param{kallisto_reads} = read_config_file_line('kallisto_reads', '', $fh_ParamConfig);

 $param{rejectGi} = read_config_file_line('rejectGi', '', $fh_ParamConfig);
  $param{negateGi} = read_config_file_line('negateGi', '', $fh_ParamConfig);
  $param{rRNA_file} = read_config_file_line('rRNA_file', '', $fh_ParamConfig);
  $param{busco_file} = read_config_file_line('busco_file', '', $fh_ParamConfig);
  $param{tnfDB_file} = read_config_file_line('tnfDB_file', '', $fh_ParamConfig);
  $param{tax_list} = read_config_file_line('tax_list', '', $fh_ParamConfig);
  $param{taxlevel_species} = read_config_file_line('taxlevel_species', '', $fh_ParamConfig);
  $param{level_upto} = read_config_file_line('level_upto', '', $fh_ParamConfig);
  $param{feedlevel_species} = read_config_file_line('feedlevel_species', '', $fh_ParamConfig);
  $param{gc_filter} = read_config_file_line('gc_filter', '', $fh_ParamConfig);
  $param{verbose} = read_config_file_line('verbose', '', $fh_ParamConfig);

#GENERAL SETTINGS
  $param{mode} = read_config_file_line('mode', '', $fh_ParamConfig);
  $param{alignment_mode} = read_config_file_line('alignment_mode', '', $fh_ParamConfig);
  $param{genescore_penalty} = read_config_file_line('genescore_penalty', '', $fh_ParamConfig);
  $param{bitscoreCutoff} = read_config_file_line('bitscoreCutoff', '', $fh_ParamConfig);
  $param{blast_task} = read_config_file_line('blast_task', '', $fh_ParamConfig);
  $param{diamond_task} = read_config_file_line('diamond_task', '', $fh_ParamConfig);
  $param{predict_gene} = read_config_file_line('predict_gene', '', $fh_ParamConfig);
  $param{random_length} = read_config_file_line('random_length', '', $fh_ParamConfig);
  $param{localdb} = read_config_file_line('localdb', '', $fh_ParamConfig);
  $param{local_gff} = read_config_file_line('local_gff', '', $fh_ParamConfig);
  $param{create_local_db} = read_config_file_line('create_local_db', '', $fh_ParamConfig);
  $param{resolution} = read_config_file_line('resolution', '', $fh_ParamConfig);
  $param{plot_contig} = read_config_file_line('plot_contig', '', $fh_ParamConfig);
  $param{min_gene_number} = read_config_file_line('min_gene_number', '', $fh_ParamConfig);
  $param{ActivateTNF} = read_config_file_line('ActivateTNF', '', $fh_ParamConfig);
  $param{ActivateBIGBLAST} = read_config_file_line('ActivateBIGBLAST', '', $fh_ParamConfig);
  $param{Uniref50_taxlist} = read_config_file_line('Uniref50_taxlist', '', $fh_ParamConfig);

# SCORING --->
  $param{alien_accepted} = read_config_file_line('alien_accepted', '', $fh_ParamConfig);
  $param{doubt_penalty} = read_config_file_line('doubt_penalty', '', $fh_ParamConfig);
  $param{GCup_penalty} = read_config_file_line('GCup_penalty', '', $fh_ParamConfig);
  $param{GCdown_penalty} = read_config_file_line('GCdown_penalty', '', $fh_ParamConfig);
  $param{GCmid_penalty} = read_config_file_line('GCmid_penalty', '', $fh_ParamConfig);
  $param{inLess_penalty} = read_config_file_line('inLess_penalty', '', $fh_ParamConfig);
  $param{blastCov_penalty} = read_config_file_line('blastCov_penalty', '', $fh_ParamConfig);
  $param{feedlevel_penalty} = read_config_file_line('feedlevel_penalty', '', $fh_ParamConfig);

# QUALITY AND PERFORMANCE --->
  $param{max_processors} = read_config_file_line('max_processors', '', $fh_ParamConfig);
  $param{blastopt} = read_config_file_line('blastopt', '', $fh_ParamConfig);
  $param{evalue} = read_config_file_line('evalue', '', $fh_ParamConfig);
  $param{qcovper} = read_config_file_line('qcovper', '', $fh_ParamConfig);

  $param{min_sequence_identity} = read_config_file_line('min_sequence_identity', '', $fh_ParamConfig);
  $param{max_sequence_identity} = read_config_file_line('max_sequence_identity', '', $fh_ParamConfig);
  $param{species_name_augustus} = read_config_file_line('species_name_augustus', '', $fh_ParamConfig);

# PATH TO EXTERNAL PROGRAMS --->
  $param{makeblastdb_path} = read_config_file_line('makeblastdb', $param{alien_dir}, $fh_ParamConfig);

  $param{blastn_path} = read_config_file_line('blastn', $param{alien_dir}, $fh_ParamConfig);
  $param{diamond_path} = read_config_file_line('diamond', $param{alien_dir}, $fh_ParamConfig);

  $param{blastx_path} = read_config_file_line('blastx', $param{alien_dir}, $fh_ParamConfig);
  $param{seqtk_path} = read_config_file_line('seqtk', $param{alien_dir}, $fh_ParamConfig);
  $param{bwakit_path} = read_config_file_line('samtools', $param{alien_dir}, $fh_ParamConfig);
  $param{augustus_path} = read_config_file_line('augustus', $param{alien_dir}, $fh_ParamConfig);
  $param{kallisto_path} = read_config_file_line('kallisto', $param{alien_dir}, $fh_ParamConfig);

# OUTPUT NAMES --->
  $param{summary} = read_config_file_line('summary', '', $fh_ParamConfig);
  $param{result} = read_config_file_line('result', '', $fh_ParamConfig);
  $param{score} = read_config_file_line('score', '', $fh_ParamConfig);

# EXTERNAL ERROR HANDLING --->
  $param{tries} = read_config_file_line('tries', '', $fh_ParamConfig);
  close($fh_ParamConfig);
  return \%param;
}

############################################################################################""
sub read_config_file_line { # file format element = value
  my ($parameter, $alien_dir, $fh_config_file) = @_;

  seek($fh_config_file, 0, 0);    # added to allow any order of parameters in the config files, preventing unfriendly error messages if the user changes the order
  while (my $line = <$fh_config_file>){
    if ($line =~ /^\s*$parameter\s*=\s*([^\s]*)\s*(.*)/) {    # the string to be searched in the file: parameter ($parameter); value ($1); description ($2)

      print $parameter.":\t".$1."\t".$2."\n";
      #print "$parameter\t$1\t$2";

      chomp ($line);
      $line =~ s/^\s*$parameter\s*=\s*//;   # removing what comes before the user input
      $line =~ s/#.*$//;                    # removing what comes after the user input (commentaries)
      $line =~ s/\s*$//;                    # removing what comes after the user input (space caracteres)
      #$line =~ s/\b*$//;                   # removing what comes after the user input (space caracteres)
      #$line =~ s/\r*$//;                   # removing what comes after the user input (space caracteres)
      #$line =~ s/\s*$//;                   # removing what comes after the user input (space caracteres)

      $line =~ s/\$alien_dir/$alien_dir/;     # allows the use of "$alien_dir" in the config file as a reference to the said parameter
      if ($line eq 'undef' || $line eq '') { return; }
      else { return $line; }
    }
  }
  return;
}

############################################################################################""
# function to identify errors in the configuration files and direct the user to the needed adjustments
sub check_parameters { #check for all parameters,
  my $param = shift;

  my $config_path = getcwd();
  $config_path =~ s/\/\w+$/\/config/;

# BASIC PARAMETER FOR LOCATION AND FILE NAME --->
  if (!defined $param->{alien_dir}) { die ("No path to alienomics was specified in SoftConfig at $config_path, please open this file and fill the parameter 'alien_dir'.\n"); }
  if (!-d $param->{alien_dir}) { die ("The path to alienomics isn't a valid directory, please check if the path in 'alien_dir' is correct: $param->{alien_dir}\n"); }
  if (!-w $param->{alien_dir}) { die ("You don't have permission to write in the alienomics directory, please redefine your permissions for this directory.\n"); }

# INPUT FILES --->
  if (!defined $param->{Input_genome}) { die ("No path to the nucleotide files was specified in your project's configuration file, please fill the parameter 'Input_genome'.\n"); }
#  if (!-d $param->{Input_genome}) { die ("The path to your project's nucleotide files isn't a valid directory, please check if the path in 'Input_genome' is correct: $param->{Input_genome}\n"); }
#  if (!-r $param->{Input_genome}) { die ("You don't have permission to read in your project's nucleotide directory, please redefine your permissions.\n"); }

  if (!defined $param->{Uniref50_taxlist}) { die ("No path to the taxlist file was specified in your project's configuration file, please fill the parameter 'Uniref50_taxlist'.\n"); }


# BLASTDB FILES --->
#  if (!defined $param->{blastdb_path}) { die ("No path to the nucleotide files was specified in your project's configuration file, please fill the parameter 'blastdb_path'.\n"); }
#  if (!-d $param->{blastdb_path}) { die ("The path to your project's database isn't a valid directory, please check if the path in 'blastdb_path' is correct: $param->{blastdb_path}\n"); }
#  if (!-r $param->{blastdb_path}) { die ("You don't have permission to read in your blast database, please redefine your permissions.\n"); }

# DIAMONDDB FILES --->
  if (!defined $param->{diamonddb_path}) { die ("No path to database was found in your project's configuration file, please fill the parameter 'diamonddb_path'.\n"); }
  if (!-d $param->{diamonddb_path}) { die ("The path to your project's database isn't a valid directory, please check if the path is correct: $param->{diamonddb_path}\n"); }
  if (!-r $param->{diamonddb_path}) { die ("You don't have permission to read in your diamond protein database, please redefine your permissions.\n"); }


# TAXDUMP FILES --->
  if (!defined $param->{taxdump_path}) { die ("No path to the nucleotide files was specified in your project's configuration file, please fill the parameter 'taxdump_path'.\n"); }
  if (!-d $param->{taxdump_path}) { die ("The path to your project's nucleotide files isn't a valid directory, please check if the path in 'taxdump_path' is correct: $param->{taxdump_path}\n"); }
  if (!-r $param->{taxdump_path}) { die ("You don't have permission to read in your project's nucleotide directory, please redefine your permissions.\n"); }



# AUGUSTUS_CONFIG_PATH FILES --->
  if (!defined $param->{AUGUSTUS_CONFIG_PATH}) { die ("No path to the nucleotide files was specified in your project's configuration file, please fill the parameter 'AUGUSTUS_CONFIG_PATH'.\n"); }
  if (!-d $param->{AUGUSTUS_CONFIG_PATH}) { die ("The path to your project's nucleotide files isn't a valid directory, please check if the path in 'AUGUSTUS_CONFIG_PATH' is correct: $param->{AUGUSTUS_CONFIG_PATH}\n"); }
  if (!-r $param->{AUGUSTUS_CONFIG_PATH}) { die ("You don't have permission to read in your project's nucleotide directory, please redefine your permissions.\n"); }

# PATH TO EXTERNAL PROGRAMS USED BY ALIENOMICS --->

  if (!defined $param->{blastn_path}) { die ("No path to BLASTN was specified in SoftConfig at $config_path, please open this file and fill the parameter 'blastn'.\n"); }
  if (!-s $param->{blastn_path}) { die ("The executable of BLASTN wasn't found in the specified path, please check if the path is correct: $param->{blastn_path}\n"); }
  if (!-x $param->{blastn_path}) { die ("You don't have permission to execute the BLASTN file specified at SoftConfig, please check permissions or replace the file\n"); }

  if (!defined $param->{diamond_path}) { die ("No path to diamond was specified in SoftConfig at $config_path, please fill the parameter 'diamond'.\n"); }
  if (!-s $param->{diamond_path}) { die ("The executable of diamond wasn't found in the specified path, please check if the path to diamond is correct: $param->{diamond_path}\n"); }
  if (!-x $param->{diamond_path}) { die ("You don't have permission to execute the diamond file specified at SoftConfig, please redefine your permissions.\n"); }


if (!defined $param->{seqtk_path}) { die ("No path to BLASTN was specified in SoftConfig at $config_path, please open this file and fill the parameter 'seqtk'.\n"); }
  if (!-s $param->{seqtk_path}) { die ("The executable of BLASTN wasn't found in the specified path, please check if the path is correct: $param->{seqtk_path}\n"); }
  if (!-x $param->{seqtk_path}) { die ("You don't have permission to execute the SEQTK file specified at SoftConfig, please check permissions or replace the file\n"); }

  if (!defined $param->{makeblastdb_path}) { die ("No path to makeblastdb was specified in SoftConfig at $config_path, please open this file and fill the parameter 'makeblastdb'.\n"); }
  if (!-s $param->{makeblastdb_path}) { die ("The executable of makeblastdb wasn't found in the specified path, please check if the path is correct: $param->{makeblastdb_path}\n"); }
  if (!-x $param->{makeblastdb_path}) { die ("You don't have permission to execute the makeblastdb file specified at SoftConfig, please check permissions or replace the file\n"); }

if (!defined $param->{augustus_path}) { die ("No path to augustus was specified in SoftConfig at $config_path, please open this file and fill the parameter 'blastn'.\n"); }
  if (!-s $param->{augustus_path}) { die ("The executable of augustus wasn't found in the specified path, please check if the path is correct: $param->{augustus_path}\n"); }
  if (!-x $param->{augustus_path}) { die ("You don't have permission to execute the AUGUSTUS PATH file specified at SoftConfig, please check permissions or replace the file\n"); }


  if (!defined $param->{project_dir_path}) {die "Project directory not configured. Please set project_dir_path element in configuration file\n";}

  if (!defined $param->{gc_filter}) {die "You must set the parameter \"gc_filter\" with an integer ranging from 0 to 100, please check ALIEN documentation form more instructions about this parameter\n";}

  if (($param->{gc_filter} !~ /[0-9]/)) {die "You must set the parameter \"gc_filter\" with an integer ranging from 0 to 100, please check ALIEN documentation form more instructions about this parameter\n";}


#Score parameters check --->

#if (!defined $param->{alien_penalty}) {die "You must set the parameter \"alien_penalty\" with an integer ranging from 0 to 100, please check ALIEN documentation form more instructions about this parameter\n";}

#  if (($param->{alien_penalty} !~ /^-?[0-9]\d*(\.\d+)?$/)) {die "You must set the parameter \"alien_penalty\" with an integer ranging from 0 to 100, please check ALIEN documentation form more instructions about this parameter\n";}

#default parameters --- set the defaukt here >

  if (!defined $param->{max_processors}) {
    $param->{max_processors} = 2;
    print "max_processors setted to 2, since it was not defined by user\n";
  }
  if (!defined $param->{predict_gene}) {
    $param->{predict_gene} = 0;
    print "predict_gene setted to FALSE(0), since it was not defined by user\n";
  }
}

# now the script loads all nucleotide sequence files to a hash structure,
# checks their validity and translates them to protein sequence

sub parse_genome_files {
  my ($param, $LOG) = @_;
  print ("opening $param->{Input_genome}\n");
  open (my $file, $param->{Input_genome}) || die ("Path to assembly fasta files not found: $!\n");
  #opendir (my $nt_files_dir, $param->{Input_genome}) || die ("Path to assembly fasta files not found: $!\n");
  my (%sequence_data);
  print "Parsing genome file and storing it\n";
  print $LOG ('Parsing assembled genome/contigs/scaffolds files', "\n") if $param->{verbose};
  #while (my $file = $nt_files_dir) {
    if (($file eq '.') || ($file eq '..') || ($file =~ /^\./) || ($file =~ /~$/)) { next; }  # Prevents from reading hidden or backup files
    my $file_content = new Bio::SeqIO(-format => 'fasta',-file => "$param->{Input_genome}");
 
    #Store the genome for augustus
    my $out_content = Bio::SeqIO->newFh(-format => 'fasta', ,-file => ">$param->{project_dir_path}/intermediate_files/gff/genome.fa");
    print $LOG ('Reading file ', $file, "\n") if $param->{verbose};
    while (my $gene_info = $file_content->next_seq()) {
      my $sequence = $gene_info->seq();
      my $accession_number = $gene_info->display_id;
      my $len = $gene_info->length;
      my $GCcount = $sequence =~ tr/GC|gc//;
      my $GCcontent = ($GCcount / $len) * 100;
      $sequence_data{$accession_number}{status} = "OK"; #everybody starts fine
      $sequence_data{$accession_number}{problem_desc} = "-"; #everybody starts fine
      if ($sequence_data{$accession_number}{status} eq "OK") { # Add check points here <<<<<<
        $sequence_data{$accession_number}{nuc_seq} = $sequence;
	$sequence_data{$accession_number}{len} = $len;
	$sequence_data{$accession_number}{gc} = $GCcontent;
	print $out_content $gene_info;
      }
    }
  #}
  #print ('Done', "\n") if $param->{verbose};
  close ($file);
  return (\%sequence_data);
}


#Parse predicted gene location
sub parse_gene_files {
  my $date=`date '+%Y-%m-%d %H:%M:%S'`; chomp $date;
  print "\n# Adding gene GC% and coverage info ($date)\n";
  my ($geneFile, $param, $covHash_ref, $transHash_ref) = @_;
  my (%sequence_data);
    my $file_content = new Bio::SeqIO(-format => 'fasta',-file => "$geneFile");
    while (my $gene_info = $file_content->next_seq()) {
      my $sequence = $gene_info->seq();
      my $accession_number = $gene_info->display_id; #here header is chr:st-ed
         $accession_number =~ s/^\s+|\s+$//g; #del space
      my $len = $gene_info->length;
      my $GCcount = $sequence =~ tr/GC|gc//; 	#For lc as well
      my $GCcontent = ($GCcount / $len) * 100;

      my $geneCov=0; # Default values for coverage
	if ($param->{coverage} > 0 ) {
         $geneCov=$covHash_ref->{$accession_number};
	}
      my $geneTrans=0.0; # Default value for transcriptomics
        if ($param->{rna_seq} > 0 ) {
         $geneTrans=$transHash_ref->{$accession_number} if $transHash_ref->{$accession_number};
#print "$geneTrans : $accession_number -------- $transHash_ref->{'Chrom_1:2155-5638'} ***\n";
        }
      $sequence_data{$accession_number}{status} = "OK";        #everybody starts fine
      $sequence_data{$accession_number}{problem_desc} = "-";   #everybody starts fine

      if ($sequence_data{$accession_number}{status} eq "OK") {       # Add check points here <<<<<<
        $sequence_data{$accession_number}{nuc_seq} = $sequence;
        $sequence_data{$accession_number}{len} = $len;
        $sequence_data{$accession_number}{gc} = $GCcontent;
        $sequence_data{$accession_number}{readcov} = $geneCov;
	$sequence_data{$accession_number}{geneexp} = $geneTrans;
      }
    }
  return (\%sequence_data);
}

sub parse_cluster_id {
  my $line = shift;
  $line =~ s/:$//;
  $line =~ s/\(\s*\d*\s*gene[s]?\s*,\d*\s*tax(a|on)\s*\)$//;
  return $line;
}

sub parse_gene_id {
  my @aux = split (/\(/, $_[0]);
  my $specie = $aux[1];
  $specie =~ s/\)//g;
  return ($aux[0], $specie);  # aux[0] has the if o the gene
}

sub mean {
  my @tmp = @{$_[0]};
  my $soma = 0;
  foreach my $value(@tmp) {
    $soma = $soma + $value;
  }
  my $mean = ($soma/($#tmp+1));
  return $mean;
}

############################################################################################""
#Time check
sub divide_time {
  my $total_time = shift;

  my $hours = POSIX::floor( $$total_time / 3600 );
  my $minutes = POSIX::floor(($$total_time % 3600) / 60);
  if ($minutes < 10) { $minutes = '0' . $minutes; }
  my $seconds = $$total_time % 60;
  if ($seconds < 10) { $seconds = '0' . $seconds; }

  return ($hours, $minutes, $seconds);
}

############################################################################################""!!
#To create the config file need for the alienomics tool
sub create_conf_file {
  my ($mode, $stringency, $fasta_path, $homology_file_path, $output_dir, $outfile, $genetic_code, $max_proc, $hom_filter) = @_;
  my $dummy_file;
  my $conf_file;
  my $radical = ""; #radical for the file name
  if ($mode eq "site") {
    $radical = "site";
  } else {
    die ("Alienomics currently only supports \"site\" for --mode, you used $mode\n.");
  }
  if ($stringency eq "standard") {
    $radical = join ("_", $radical, $stringency);
  } elsif ($stringency eq "soft") {
    $radical = join ("_", $radical, $stringency);
  } elsif ($stringency eq "hard") {
    $radical = join ("_", $radical, $stringency);
  } else {
    die ("Currently supported stringency values are \"standard\", \"soft\" and \"hard\", you provided $stringency\n");
  }

  $conf_file = "alien_main_$radical.conf";

  if ((defined $outfile)&&($outfile ne "")) {
  if ($outfile !~ /\.conf$/) {
    $conf_file = $outfile.".conf";
  } else {
    $conf_file = $outfile;
    }
  }

  if (-e "alien_main_$radical.conf") {
        die ("A configuration file with the parameters you chose (alien_main_$radical.conf) is already present in this directory.\nPlease remove or rename before creating a new alien configuration file\n\n");
  }
  $dummy_file = "../config/dummy_alien_ParamConfig_file_$radical";

  open(FILE, "<$dummy_file") || die ("Configuration file $dummy_file not found, please check you used a valid set of parameters\n");
  my @lines = <FILE>;
  close(FILE);
  for (my $i = 0; $i <= $#lines; $i++) {
    my $line = $lines[$i];
    if (($line =~ /^Input_genome/)&&(defined $fasta_path)) {
      $line = "Input_genome = $fasta_path\n";
      $lines[$i] = $line;
    }
    if (($line =~ /^project_dir_path/)&&(defined $output_dir)) {
      $line = "project_dir_path = $output_dir\n";
      $lines[$i] = $line;
    }
    if (($line =~ /^codon_table/)&&(defined $genetic_code)) {
      $line = "codon_table = $genetic_code\n";
      $lines[$i] = $line;
    }
    if (($line =~ /^homology_file_path\s+/)&&(defined $homology_file_path)) {
      $line = "homology_file_path = $homology_file_path\n";
      $lines[$i] = $line;
    }
    if (($line =~ /^max_processors/)&&(defined $max_proc)) {
      $line = "max_processors = $max_proc\n";
      $lines[$i] = $line;
    }
    if (($line =~ /^homology_filter\s+/)&&(defined $hom_filter)) {
      $line = "homology_filter = $hom_filter\n";
      $lines[$i] = $line;
    }
  }
  open(FILE, ">$conf_file") || die "could not open $conf_file, please check directory permissions\n";
  print FILE @lines;
  close(FILE);
# copy($dummy_file, $conf_file) || (die($!));
}

############################################################################################""!!
#Parse the augustus outfile
sub parseAugustus {
use List::MoreUtils qw{ any };
  my ($param, $seqfile) = @_; #"preGenes.txt"
  my %genCor; my %genIn; my %allGeneInfo;
  open FL, $seqfile;
  local $/ = "\n";  # read by \n
  while (<FL>) {
    chomp;
    next if /^#/;  # discard comments
    my @seqLine = split('\t', $_);
    my @seqL = split('gene_id', $seqLine[8]);
    $seqL[-1] =~ s/"|;//g; $seqL[-1] =~ s/^\s+|\s+$//g;
    $allGeneInfo{$seqL[-1]}{$seqLine[2]} = 1;
    if ($seqLine[2] eq 'intron') { $genIn{$seqLine[8]}=$seqLine[2];} #Store all intron info using their gene ids
    next if ($param->{mode} ne "$seqLine[2]");
    my $region="$seqLine[0]\t$seqLine[3]\t$seqLine[4]\t$seqLine[8]\t$seqLine[6]\t$seqLine[2]";
    $genCor{$region}=$.;
}
return (\%genCor, \%genIn, \%allGeneInfo);
close FL;
}

############################################################################################""!!
#Function to move the files with wild
sub moveFiles {
    my ( $source_ref, $arc_dir ) = @_;
    my @old_files = @$source_ref;
    foreach my $old_file (@old_files)
         {
    #my ($short_file_name) = $old_file =~ m~/(.*?\.dat)$~;
    #my $new_file = $arc_dir . $short_file_name;
    move($old_file, $arc_dir) or die "Could not move $old_file to $arc_dir: $!\n";
   }
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
#Store fasta file to hash
sub fastafile2hash {
    my $fastafile = shift @_;
    my %sequences;
    my $fh = &read_fh($fastafile);
    my $seqid;
    while (my $line = <$fh>) {
        if ($line =~ /^>(\S+)(.*)/) {
            $seqid = $1;
            $sequences{$seqid}{desc} = $2;
        }
        else {
            chomp $line;
            $sequences{$seqid}{seq}     .= $line;
            $sequences{$seqid}{len}     += length $line;
            $sequences{$seqid}{gc}      += ($line =~ tr/gcGC/gcGC/);
            $line =~ s/[^atgc]/N/ig;
            $sequences{$seqid}{nonatgc} += ($line =~ tr/N/N/);
        }
    }
    close $fh;
    return \%sequences;
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


##############################################################
#Process the GCAT sequence
sub processGCAT {
    my $sequence = shift;
    my @letters = split(//, $sequence);
    my $gccount = 0; my $totalcount = 0; my $gccontent = 0;
    my $acount = 0; my $tcount = 0; my $gcount = 0; my $ccount = 0; my $atcontent =0;
    foreach my $i (@letters) {
	if (lc($i) =~ /[a-z]/) { $totalcount++;}
	if (lc($i) eq "g" || lc($i) eq "c") { $gccount++; }
	if (lc($i) eq "a") { $acount++;}
	if (lc($i) eq "t") { $tcount++;}
	if (lc($i) eq "g") { $gcount++;}
	if (lc($i) eq "c") { $ccount++;}
    }
    if ($totalcount > 0) {
	$gccontent = (100 * $gccount) / $totalcount;
    }
    else {
	$gccontent = 0;
    }
    return "$gccontent\t$totalcount\t$gcount\t$ccount\t$acount\t$tcount";
}

############################################################
# Intergenic region extraction
sub intergenicJotter {
  my ($genes_ref, $seq, $acc) = @_;
  my %genes=%$genes_ref; my $start=0; my %igrHash; my @alligrGC;
  for my $key ( sort {$a<=>$b} keys %genes) {
    my $igrSize=$key-$start;
    my $intergenicString = substr $seq, $start, $igrSize;
    my $intergenicStat = &processGCAT($intergenicString);
    $start=$genes{$key};
    my @vals = split(/\t/, $intergenicStat);
    push @alligrGC, $vals[0];
  }
my $GC=meanGCAT(@alligrGC);
return "$GC";
}

############################################################
#Print the hash values
sub print_hash_final {
    my ($href,$fhandler)  = @_;
    while( my( $key, $val ) = each %{$href} ) {
        print $fhandler "$key\n";
	#print $fhandler "$key\t=>$val\n";
    }
}

############################################################
# Returns 1 if present else 0
sub isInList {
   my $needle = shift;
   my @haystack = @_;
   foreach my $hay (@haystack) {
	if ( $needle eq $hay ) {
	return 1;
	}
   }
   return 0;
}

############################################################
#Mean of GCAT
use List::Util qw(sum);
sub meanGCAT { return @_ ? sum(@_) / @_ : 0 }
#sub meanGCAT { return sum(@_)/@_; }


############################################################
#  Sigmoid function
sub sigmoidFun {
my $val = shift;
   my ($h) = @_;
   return 1.0 / ( 1.0 + exp(-$val) );  # produce gene score values from 0 to +1 (sigscore)
}

############################################################
#  tanh alternative to sigmoid function above
sub tanhGeneScore {
my $val = shift;
   my ($h) = @_;
   return tanh ($val) ;            # produce gene score values from -1 to +1 (tanhgenescore)
}

############################################################
sub extractBlastCov {
  my ($taxidP_ref, $taxidT_ref, $taxidN_ref, $tax_list, $fileName)=@_;
  open my $file, '<', "$fileName";
  my $firstLine = <$file>;
  close $file;
  print "$firstLine >><<<<>>\n";
  my @vals = split(/\t/, $firstLine);
  $vals[-1] =~ s/^\s+|\s+$//g;
  my @allAnnot= extractTax($taxidP_ref, $taxidT_ref, $taxidN_ref, $vals[0], $vals[1], $tax_list);
  foreach my $elem (@allAnnot) {
	$elem =~ s{^\s+|\s+$}{}g;
 	$elem =~ s/\s+/_/g;
   	$elem = lc $elem;
  }
print "$allAnnot[0], $allAnnot[1], $allAnnot[2], $allAnnot[3] ==================================\n";
return ($vals[-1], $allAnnot[0], $allAnnot[1], $allAnnot[2], $allAnnot[3]); # the last values of the string

}

############################################################

sub sumArray {
no warnings 'recursion';
    return( defined $_[0] ? $_[0] + sumArray(@_[1..$#_]) : 0 );
}


############################################################
sub extractTax {
	my ($taxidP_ref, $taxidT_ref, $taxidN_ref, $nam, $id, $tax_list)=@_;
	my %taxidP=%$taxidP_ref; my %taxidT = %$taxidT_ref; my %taxidN = %$taxidN_ref;
	my %contig_taxinfo; my @allNames;
	my @tax_list=split(/\,/, $tax_list);
	my %tax_levels;
	foreach (@tax_list) {
		$tax_levels{$_}=1
	}
	$contig_taxinfo{$nam}= &taxonomy_report($id, \%taxidT, \%taxidN, \%tax_levels, \%taxidP);
	for my $tax_level (@tax_list) {
		push @allNames, ("\t" . (exists(${$contig_taxinfo{$nam}}{$tax_level}) ? ${$contig_taxinfo{$nam}}{$tax_level} : "NA"));
	}
	return @allNames;


############################################################
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

############################################################

sub spacer { return (" " x 20000);}

############################################################
#Open and write a file
sub write_fh {
    my ($filename, $append) = @_;
    my $filehandle;
    if ($append) {
        open $filehandle, ">>$filename" or die "Can't open $filename for writing: $!";
    }
    else {
	open $filehandle, ">$filename" or die "Can't open $filename for writing: $!";
    }
    return $filehandle;
}

############################################################
sub extractCoordinates {
my ($gffFile, $fastaFile, $chrName) = @_;

# PURPOSE: for each gene declared in -gff, write its spliced CDS sequence taken from -fasta
# USAGE: program.pl file.fa file.gff

use File::Basename;
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::SeqFeature::Store;
use Bio::SeqFeatureI;

my $outfile = basename($fastaFile, ".fa") . "_cds.fa"; #formats name of infile for outfile
#open (OUT, ">$outfile") or die "Can't open file for writing: $!";
my $wfh = write_fh ($outfile, 0); # Put 1 to append, else 0
#Load GFF
my $db = Bio::DB::SeqFeature::Store->new(   -adaptor => 'memory',
                                            -fasta => $fastaFile,
                                            -gff => $gffFile) or die $!;

my $ss = $db->get_seq_stream(-type => 'gene'); #looks at gff, find genes and their associated information

while (my $gene = $ss->next_seq) #while looking at each gene
{
    next if $gene -> display_name ne $chrName;
    print $wfh ">" , $gene -> display_name , "\n"; #fasta id for CDS
    #print OUT $gene -> start , "\t", $gene -> end , "\t", $gene -> strand, "\n"; #test to see if script was reading gff correctly
    0 && q {
    my @segments = $gene -> segments(-type => 'CDS');
    my %START;
    my %END;
    foreach my $segment(@segments) {
        $START{$segment} = $segment -> start; #create hash of CDS id and start nucleotide
        $END{$segment} = $segment -> end; #create hash of CDS id and end nucleotide
    }
    my @CDS;
    foreach my $segment (sort{$START{$a} <=> $START{$b}} keys %START) #sorts CDS segments by order of starting nucleotide
    {
        my $cds_for_seq = $db->fetch_sequence($gene -> seq_id , $START{$segment}, $END{$segment}); #extracts CDS segment sequence
        push(@CDS, $cds_for_seq); #puts CDS segment sequence to end of array

        #print OUT $segment , "\t" , $START{$segment} , "\t", $END{$segment} ,"\n"; #test to see if CDS segments are ordered correctly
        #print OUT $cds_for_seq , "\n";
    }
    my $CDS_seq = join('', @CDS); #concatenate array of CDS segment sequences

    my $for_obj = Bio::Seq -> new(-seq => $CDS_seq , -alphabet => 'dna'); #load sequence string as sequence object
    my $rev_obj = $for_obj -> revcom; #reverse complement sequence object
    my $rev_seq = $rev_obj -> seq; #convert sequence object to sequence string

    if ($gene -> strand == +1) #check if gene was on plus strand
    {
        print $wfh $CDS_seq , "\n";
    }
    elsif ($gene -> strand == -1) #check if gene was on minus strand
    {
        print $wfh $rev_seq , "\n";
    }
    };
}

}


############################################################

=pod
sub round {

    my ($nr,$decimals) = @_;
    return (-1)*(int(abs($nr)*(10**$decimals) +.5 ) / (10**$decimals)) if $nr<0;
    return int( $nr*(10**$decimals) +.5 ) / (10**$decimals);

}
=cut


############################################################
#Predict the genes and store them in gff folder
sub predictGeneAug {
my ($genome, $param, $decision )=@_;
print "Predicting genes with Augustus: using $param->{species_name_augustus} species\n";
# All parameters set for augustus
system ("$param->{augustus_path} --species=$param->{species_name_augustus} $genome --outfile=$param->{project_dir_path}/intermediate_files/gff/preGenes.gff --AUGUSTUS_CONFIG_PATH=$param->{AUGUSTUS_CONFIG_PATH}");
}

############################################################
#print the lines
sub pLine {
my $msg = shift;
print "$msg" x 80 . "\n";
#print ($msg x 20);
}

1;

__END__

#Make blast to run on all genes at once
#use GFF file to extract all region
#augustus output to gff
