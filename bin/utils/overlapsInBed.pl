#!/usr/local/bin/perl
use strict;
use Switch;

#
# overlap.pl by William Pu
#
# perl overlap.pl afile.bed bfile.bed outfile.txt [-bp bases -name "basename"]
# calculate intersection between a and b, send result to outfile
# does not assume any sort order
# a is read into memory, b is read from disk
# outfile format is a ""modified gff"":
#   chr overlap.pl basename_number start(0base) end overlap_score . . afile.bed a_name chr st end score; bfile.bed b_name chr st end score
#   overlap_score: 1=a; 2=b;3=both
#   for class 1 or 2, the non-overlapping score is randomly assigned a value from -100 to 0.
# -bp bases: min number of bases in the intersection
# -name basename: basename of peak names
my $bed_chr=0;
my $bed_st=1;
my $bed_end=2;
my $bed_name=3;
my $bed_score=4;
my $bed_strand=5;
my $bed_match=6;	# use to record whether or not the line in the a file matched a line in the b file

		
# read parameters, assign output variables
my $afile = $ARGV[0] ;
if ($afile eq "") {
	report_error("invalid input file a");
}
my $bfile = $ARGV[1] ;
if ($bfile eq "") {
	report_error("invalid input file b");
}
my $outfile = $ARGV[2] ;
if ($outfile eq "") {
	report_error("invalid output file");
}

my $min_overlap=1;
my $basename="overlap";
for (my $i=3;$i<@ARGV;$i++) {
	switch($ARGV[$i]) {
		case "-bp"	{ $i++; $min_overlap=$ARGV[$i++] }
		case "-name" { $i++; $basename=$ARGV[$i++] }
		else { report_error("invalid parameter " . $ARGV[$i]) }
	}
}

my @a;
my @b;
my $a_line_entry;
my $b_line_entry;
my $a_row=0;
my $b_row=0;
my $a_lineref;
my @aline;
my @bline;
my $count=0;
my $out_st;
my $out_end;
my $bmatch;
my @amatch;
my $flag;
my $s;
my $achr;
my %atable;
# open files -- read file a into array of arrays a
unless ( open( AFILE, $afile ) ) {
		report_error("Could not open file $afile");
	}

while ($a_line_entry=<AFILE>) {
	chomp($a_line_entry);
	@aline=split /\t/,$a_line_entry;
	push @a,[@aline];
}
close(AFILE);
# now make a hash of array references chr=>(ref1 ref2 ref3...)
for ($a_row=0;$a_row<@a;$a_row++) {
	push @{$atable{$a[$a_row][$bed_chr]} }, \@a[$a_row];
}

unless ( open( BFILE, $bfile ) ) {
		report_error("Could not open file $bfile");
	}

unless ( open( OUTFILE, ">",$outfile ) ) {
		report_error("Could not open file $outfile");
	}

while ($b_line_entry=<BFILE>) {	# read a row of bfile and compare to all a file entries
	chomp($b_line_entry);
	@bline=split(/\t/,$b_line_entry);
	$bmatch=0;
	$achr = \@{$atable{$bline[$bed_chr]}};	# $achr is pointer to the array of arrays from the chromosome matching b
	for ($a_row=0;$a_row<@$achr;$a_row++) {
		$a_lineref=$achr->[$a_row];			#a_lineref is a pointer to an array containing one line of bed entries
		($flag,$out_st,$out_end)=compareab($$a_lineref,\@bline,$min_overlap);	# returns (flag,st,end) where flag=1 if overlap, 0 if no
		if ($flag ) {	# intersect
			# print out the line
			print OUTFILE $$a_lineref->[$bed_chr]."\toverlap.pl\t$basename"."_"."$count\t$out_st\t$out_end\t3\t.\t.\t$afile\t".$$a_lineref->[$bed_name]."\t".$$a_lineref->[$bed_chr]."\t".$$a_lineref->[$bed_st]."\t".$$a_lineref->[$bed_end]."\t".$$a_lineref->[$bed_score]."\t$bfile\t$bline[$bed_name]\t$bline[$bed_chr]\t$bline[$bed_st]\t$bline[$bed_end]\t$bline[$bed_score]\n";
			$count++;
			$bmatch=1;
			$$a_lineref->[$bed_match]=1; 
		}
	}
	# if no match to that line of b, then print it out
	if ($bmatch==0) {
		$out_st=$bline[$bed_st];
		$out_end=$bline[$bed_end];
		print OUTFILE "$bline[$bed_chr]\toverlap.pl\t$basename"."_"."$count\t$out_st\t$out_end\t2\t.\t.\t$afile\t.\t.\t.\t.\t".-1*rand(100)."\t$bfile\t$bline[$bed_name]\t$bline[$bed_chr]\t$bline[$bed_st]\t$bline[$bed_end]\t$bline[$bed_score]\n";
		$count++;
		}
}

# now print out all the rows of a that had no match
for ($a_row=0; $a_row<@a;$a_row++) {
	if ($a[$a_row][$bed_match]!=1) {
		$out_st=$a[$a_row][$bed_st];
		$out_end=$a[$a_row][$bed_end];
		print OUTFILE "$a[$a_row][$bed_chr]\toverlap.pl\t$basename"."_"."$count\t$out_st\t$out_end\t1\t.\t.\t$afile\t$a[$a_row][$bed_name]\t$a[$a_row][$bed_chr]\t$$a[$a_row][$bed_st]\t$a[$a_row][$bed_end]\t$a[$a_row][$bed_score]\t$bfile\t.\t.\t.\t.\t".-1*rand(100)."\n";
		$count++
	}
}
	
close(BFILE);
close(OUTFILE);	
exit;

# compare line a to b for overlap 0<end-st<=min_overlap
# return list (flag,st,end) flag=1 means overlap

sub compareab {
	my ($aline,$bline,$min_overlap) = @_;
	my $left=max($aline->[$bed_st],$bline->[$bed_st]);
	my $right=min($aline->[$bed_end],$bline->[$bed_end]);
	if ($right-$left>=$min_overlap) {
		return(1,$left,$right);
	}
	return(0,0,0);
}

sub min {
	my ($x,$y) = @_;
	if ($x<$y) {
		return($x);
	} else {
		return($y);
	}
}

sub max {
	my ($x,$y) = @_;
	if ($x>$y) {
		return($x);
	} else {
		return($y);
	}
}

sub report_error {
	print "@_\n";
	print "usage: perl overlap.pl afile.bed bfile.bed outfile.txt [-bp bases -name basename]\n";
	print "-bp bases specifies minimum overlap; default = 1";
	print "-name basename specifies name for output peaks; default=overlap\n";
	die;
}

