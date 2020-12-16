# $Id$ Module for Module
# Perl module for ALIENOMICS AOM::Module
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2020 by Jitendra. All rights reserved.
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
package AOM::tnfHandler;

use strict;
use warnings;
use Getopt::Long;
use Benchmark qw(cmpthese);
use Data::Table;

#Calculate TNF - the main subs
sub calculateTNF {
my ($ffile, $mlength, $kmerl) = @_;
my ($kmerlprobe, $header, $seq, $prevheader, $subseq, @probes, @kmerheader, %kmer, @allFreq);
my $sequence = "AAA";
my $count = -1;
my $prc = 0; #print read count
my $scount = 0;
my $scg = 0; #good sequence counter
my $tkmers = 0;#sum of all kmers
my $output = "Contig";

# CODE HERE
# open the fasta file
open(IN, $ffile) or die("Cannot open $ffile\n");

#Create all possible kmers
$kmerlprobe='';
for (my $kul = 0; $kul < $kmerl; $kul++)  {  $kmerlprobe = $kmerlprobe."N"; }
push (@probes,$kmerlprobe); # store all kmer probes 

foreach my $probe (@probes){
	if ($probe =~ m/N/) { 								
		my $temp1 = $probe;
		$temp1 =~ s/N/A/; push (@probes, "$temp1"); $temp1 = $probe;
		$temp1 =~ s/N/T/; push (@probes, "$temp1"); $temp1 = $probe;
		$temp1 =~ s/N/C/; push (@probes, "$temp1"); $temp1 = $probe;
		$temp1 =~ s/N/G/; push (@probes, "$temp1");								
	}				
}

foreach my $probe (@probes){
	if ($probe =~ m/N/) { 	} #Weired way to skip, chance it later
	else{
		$output = $output."\t$probe";
		push (@kmerheader,$probe);
		$kmer{$probe} = 0;
	}
}
#print OUT "$output\n";

#@kmerheader= permuteG ($kmerl);
#print scalar (@kmerheader); foreach (@kmerheader) { print "$_\n"; }	

while ( my $line = <IN> ) {
	chomp $line;
	$tkmers = 0; $count++;
        #reset kmers here first	
	foreach my $probe1 (@kmerheader) { $kmer{$probe1} = 0; }
	if ($line =~ m/>/) {
		if ($mlength <= length($sequence) and $count > 0) { #print "$line ---->\n";
			for (my $kul = 0; $kul < length($sequence); $kul++)  {
				if (exists($kmer{substr($sequence,$kul,$kmerl)})){#to escape N's in scaffolds
					$kmer{substr($sequence,$kul,$kmerl)}++;				
					my $rc = reverse(substr($sequence,$kul,$kmerl));#Add reverse complement to escape string bias
					$rc =~ tr/ACGT/TGCA/;	
					$kmer{$rc}++;
					$tkmers+= 2;
				}
			}
			if ($tkmers > 0){
				foreach my $probe1 (@kmerheader){				
					my $temp1 = $kmer{$probe1}/$tkmers;	
					my $len=length($sequence);
					#print "$kmer{$probe1}\t$len\t$tkmers\t $temp1 +++\n";			
					$output = $output."\t".sprintf("%.5f",$temp1)*100;
				}
			#print "$output >->\n"; #OUT 
			push @allFreq, $output;
			$scg++;
			}
		}
		$sequence = "";	
		$header = $line;
		$header =~ s/>//;
		$output = $header;									
		if ($prc == 100) {
			$prc = 0;
			print "$scount sequences $scg \>= $mlength bp\n";
		}	
		$scount++;		
		$prc++;		
	}	
	else{		
		$line = uc($line);
		$sequence = $sequence.$line;
	}
}

#Lets count the last sequence frequency
$count++;
$tkmers = 0;
if ($mlength <= length($sequence)) { #weired but work...
	foreach my $probe1 (@kmerheader) { $kmer{$probe1} = 0; }
	for (my $kul = 0; $kul < length($sequence); $kul++)  { # kul hindi
		if (exists($kmer{substr($sequence,$kul,$kmerl)})){                                                    
			$kmer{substr($sequence,$kul,$kmerl)}++;				
			my $rc = reverse(substr($sequence,$kul,$kmerl));                                                 
			$rc =~ tr/ACGT/TGCA/;		
			$kmer{$rc}++;
			$tkmers+= 2;
		}
	}
	$output = $header;
	foreach my $probe1 (@kmerheader){
		if (!$tkmers) {print "Something went wrong .. totalkmers($tkmers) is ZERO in tnfHandler\n";}
		my $temp1 = $kmer{$probe1}/$tkmers;
		my $len=length($sequence);
		#print "$kmer{$probe1}\t$len\t$tkmers\t $temp1 ---\n";
		$output = $output."\t".sprintf("%.5f",$temp1)*100;
	}
	push @allFreq, $output;
}
close IN;
#print scalar(@allFreq); print "+-+-+\n";
return \@allFreq;
}

#Search TNF in the TNF dataset
sub SearchTNF {
my ($TNFrequency_ref, $TNFMatrix, $outfile, $gName) = @_;
open my $ofh, ">>", $outfile or die "$outfile: $!";
# Create a Data::Table with headers (assuming data is tab-delimited)
my $dt = Data::Table::fromTSV( "$TNFMatrix", 1 );
# Get the number of rows in the Data::Table
my $n_rows = $dt->nofRow;

my $query = [ @$TNFrequency_ref ];

my $nearest_name = '';
my $min_dist;
foreach my $i (0..$n_rows - 1){
	my $row_ref = $dt->rowRef($i); # Get row of Data::Table as an ARRAY REF
    	my $name = shift @{$row_ref}; # The name is in the first column
    	my $dist = dist($query, $row_ref);
    	$min_dist = !defined($min_dist) ? $dist
                : $dist < $min_dist ? $dist
                : $min_dist;
    	$nearest_name = $dist <= $min_dist ? $name : $nearest_name;
}
my $queryFreq = join(",", @{$query});
print $ofh "$gName\t$nearest_name\t$queryFreq\n";
#print $ofh "The nearest to:";
#print $ofh join(", ", @{$query});
#print $ofh " is: $nearest_name\n";
close $ofh;
return ("$gName","$nearest_name","$min_dist");
}


# Calculate the Euclidean distance between two vectors
sub dist {
my ($x, $y) = @_;
unless(ref($x) eq 'ARRAY' and ref($y) eq 'ARRAY'){
        die "Vectors must be given as array references"
}
unless (scalar @{$x} == scalar @{$y}) { print "scalar  : {$x} scalar @{$y}\n";
      	die "Vectors are not of equal length";
}

my $sum_sq = 0;
my $len = scalar @{$x};
foreach my $i (0..$len - 1) {
    $sum_sq += ($x->[$i] - $y->[$i])**2;
}
return sqrt($sum_sq);
}

#subs here
#Permutation glob
sub permuteG {
my ($size) = @_;
my  $r     = 'A,T,G,C';
	return glob "{$r}" x $size;
}

#Permutation loop
sub permute_loop {
my ($size) = @_;
my  @a     = qw(A T G C);
while (--$size) {
	@a = map { $_ . 'A', $_ . 'T', $_ . 'G', $_ . 'C'} @a;
}
return @a;
}

1;
