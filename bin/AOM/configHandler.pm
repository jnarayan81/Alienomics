
package AOM::configHandler;

sub create_sample_config {
my ($config_location)=@_;

print AOM::Module::fetchLocalTime()."Generating sample config file for Alienomics\n";

my $str = <<END;

###---CONFIGURATION FILE FOR ALIENOMICS---###

### INPUTS ###
input_genome_fasta = /home/urbe/Tools/Alienomics_v0.3/input_AV20/A_vaga.NDPD.fasta                      # genome fasta file
input_genome_gff = /home/urbe/Tools/Alienomics_v0.3/input_AV20/Adineta_vaga.gff3                        # gene annotation gff file
input_mapping_bam = /home/urbe/Tools/Alienomics_v0.3/input_AV20/A_vaga.NDPD.ONT.bam.minimap2.sorted     # sorted and indexed BAM
RNAseq_location = /home/urbe/Tools/Alienomics_v0.3/RNA-Seq_AV20                                         # RNA-Seq reads location
reference_proteome_file = /home/urbe/Tools/Alienomics_v0.3/input_AV20/REF.proteome.fasta                # fasta file
# (Warning: all HGTs also existing in the reference proteome will likely be considered as 'self' in the genome under analysis)
#reference_genome_file = /home/urbe/Tools/Alienomics_v0.1/Alienomics/test2.REF.paul.fasta               # fasta file

### OUTPUTS ###
output_folder = /home/urbe/Tools/Alienomics_v0.3/alienomics-v0.3-test1.output.AV20                      # output folder name

### IMPORTANT SETTINGS ###
clade_self = metazoa                    # Clade name. Blast hits within that clade = 'self', blast hits outside that clade = 'alien"
excludeTaxID = 104782|10195|96448|249248|1813166|104781         # list of taxonomic IDs that will be removed from uniref50 database
expected_GC_range = 26:38               # Expected GC range of the organism of interest (two values separated with character ':')
expected_coverage = 100                 # Expected expected_coverage of genes. Recquire reads using the 'input_mapping_bam' option
max_processors = 40                     # Number of processors to use for parallelized steps

### BLAST SETTINGS ###
evalue = 1e-01                          # Blast evalue threshold
qcovper = 0                             # Minimum query expected_coverage

### PATH TO DATABASES ###
# (do not add a "/" character at the end of the paths) ###
uniref50_path = /path/to/uniref50/folder                                 
#Uniref50_taxlist = /home/urbe/Tools/Alienomics_v0.1/Alienomics/Uniref50DB/uniref50.taxlist             
rRNA_file =  /path/to/SILVA_128_LSUParc_tax_silva.fasta
busco_file = /path/to/BUSCOData/metazoa_odb9/ancestral
#tnfDB_file = /home/urbe/Tools/Alienomics_v0.1/Alienomics/TNFDB/tnf_1975_genomes.longuestseq.tab

###---ADDITIONAL SETTINGS---###

blast_task = blastn                             # when using blast+ on a nucleotidic database
diamond_task = blastx                           # when using diamond on a proteic database
verbose = 1                                     # Print all detail #1 to print nice log messages telling you what is going on
resume = 0                                      # Resume the last run
ActivateBIGBLAST = yes                          #
ActivateTNF = no                                # 'no' by default; if set to 'yes', it will run TNF criteri

rejectGi=/home/urbe/Tools/Alienomics_v0.1/Alienomics/config/rejectList.gi       # GI or Accession list for blast to ignore
negateGi=0                                                                      # ignore the gi list in file: 1 for yes ; 0 for no

summary = summary.out                           # summary file name
result = scaffoldscore.out                      # per-scaffold info file name
score = genescore.out                           # per-gene info file name
log_file = final_log.log                        # Alienomics log file name

tax_list = species,order,phylum                 # If uncertain, do not change this parameter. Other option: "species,order,phylum,kingdom"

alignment_mode=diamond                          # Using diamond (default). Other option is "blast" (Alienomics run will be much slower)
genescore_penalty=1                             # Keep it >0. It help to sharpen the blob. The bigger, the more it forces gene score values towards 0 and 1
bitscoreCutoff=150                              # bitscore threshold for blast resutlts. smaller values produce scores closer to -1 and 1

##########################################################
########## Useless or Obsolete parameters ################

alien_dir = /home/urbe/Tools/Alienomics_v0.2/bin/                                                       # useless now ?
blastdb_path = /home/urbe/Tools/Alienomics_v0.1/Alienomics/blastDB                                      # useless now ?
taxdump_path = /home/urbe/Tools/Alienomics_v0.1/Alienomics/taxdump                                      # useless now ?
accession2taxid_path = /home/urbe/Tools/Alienomics_v0.1/Alienomics/accession2taxid                      # useless now ?
AUGUSTUS_CONFIG_PATH = /home/urbe/Tools/Alienomics_v0.1/Alienomics/augustus.2.5.5/config                # useless now ?
localdb = yes                                   # Check into the local database created by "create_local_db" above -- blast database of you reference genome
create_local_db = yes                           # Create a local blast db  ... Let them keep always "yes" for NOW
                                                # Make localDB makeblastdb -in ../WadinetaGapFilled.gapfilled.iteration1.fa -parse_seqids -dbtype nucl -out m
                                                # Make localDB makeblastdb -in ../WadinetaGapFilled.gapfilled.iteration1.fa -parse_seqids -dbtype nucl -out m

# Graph Plotting Settings
resolution = 1000                               # Resolution for plotting the gene graph
plot_contig = yes                               # "yes" to print the contig plot with score else "no" --- It will increase the computtion time
min_gene_number = 3                             # Minumum gene to consider for plotting graph

### Sequence Parameters
taxlevel_species = adineta_vaga                 # taxlevel_species you believe their presence make it strong candidate !!!! for example adineta_ricciae
feedlevel_species = escherichia_coli            # Species you feeded with .... currently only one supported !!!!
codon_table = 11
absolute_min_sequence_size = 500                # minimum sequence length cutoff for sequence/group further evaluation !!!
absolute_max_sequence_size = 10000              # maximum sequence length cutoff for sequence/group further evaluation !!!!!
sequence_identity_comparison = nt               # the kind of sequence that will be used when computing aliens !!!!!!
group_identity_comparison = gene                # sequence type that will be used when computing aliens !!!!!!
min_sequence_identity = 70                      # minimum (mean/median) sequence identity cutoff in pairwise sequence alignments
max_sequence_identity = 100                     # maximum (mean/median) sequence identity cutoff in pairwise sequence alignemnts

#random_length = 500                            # Random length of the fragments NOT FUNCTIONAL NOW
#believe=mammal,invertibrate,reference          # Not used, replace with level upto in this update
#remove_identical = yes                         # "yes" to remove, "no" otherwise
predict_gene = 0                                # 1 to yes / O to NO [default should be "0" !]
mode = gene                                     # Main analysis mode. Currently ALIEN supports only gene and random models analysis.
###---Third party tool configuration---###

multiple_alignment = bwa                        # program used for mapping with bowtie, bwa and segemehal !!!!
species_name_augustus = caenorhabditis          # Name of the species for augustus gene predictions !!!
phylogenetic_tree = bwa                         # program used for phylogenetic tree reconstruction. .... NOT USED
blastopt = offline                              # Blast running mode: Possible values are "online" or "offline" !!!!
pvalue = 0.05                                   # p-values for positive selection detection
qvalue = 0.05                                   # q-values for positive selection detection

# path to external programs
diamond = /usr/local/bin/diamond
blastn = /usr/bin/blastn
blastx = /usr/bin/blastx
makeblastdb= /usr/bin/makeblastdb
augustus = /usr/local/bin/augustus
bwakit = //usr/bin/bwa
seqtk = /home/urbe/anaconda3/envs/seqtk_paul/bin/seqtk
samtools = /usr/local/bin/samtools
kallisto = /home/urbe/Tools/kallisto/kallisto

END

my $filename = 'sample.conf';

open(FH, '>', "$config_location/$filename") or die $!;

print FH $str;

close(FH);
#print AOM::Module::fetchLocalTime()."DONE\n";
print AOM::Module::fetchLocalTime()."Alienomics configuration file successfully created!\n";
print AOM::Module::fetchLocalTime()."Find out sample alienomics config file at $config_location/$filename\n";

}

1;
