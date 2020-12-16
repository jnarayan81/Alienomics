These scripts are intended to get and parse sequence files
and make them compatible with ALIENOMICS.

#####################
# get_seq_by_ids.pl #
#####################

Takes as input a single argument, namely a list containing genbank IDs (one per
line). Downloads each genbank file in current directory.

#####################
# get_seq_by_ids.pl #
#####################

Takes as input three arguments: a text file with a list of IDs (one per line),
the path to a fasta file and the name of an output file. Produces a fasta file
containing the sequences in the list of ids.
