package AOM::downloadHandler;

use strict;
use warnings;
use ExtUtils::Installed;
use File::Spec;
use File::Basename;
use Archive::Extract;
use LWP::Simple;

#Download databases

sub download_database {
my ($third_party_DB, $flags)=@_;
#Store everthing in third_party_DB folder of your location
if ($third_party_DB) { $third_party_DB = "$third_party_DB/third_party_DB";};

print <<'WELCOME';

   ___                          _                     _                 
  |   \   ___  __ __ ___ _     | |    ___   __ _   __| |   ___     _ _  
  | |) | / _ \ \ V  V | ' \    | |   / _ \ / _` | / _` |  / -_)   | '_| 
  |___/  \___/  \_/\_/|_||_|  _|_|_  \___/ \__,_| \__,_|  \___|  _|_|_  
_|"""""_|"""""_|"""""_|"""""_|"""""_|"""""_|"""""_|"""""_|"""""_|"""""| 
"`-0-0-"`-0-0-"`-0-0-"`-0-0-"`-0-0-"`-0-0-"`-0-0-"`-0-0-"`-0-0-"`-0-0-' v0.0.1

Automated alienomics dependencies downloader >>>

WELCOME

#my $dir = dirname(File::Spec->rel2abs(__FILE__));

#
# First, check if all the required modules have been installed in the system and download the mandatory database
#
print AOM::Module::fetchLocalTime()."Checking alienomics mandatory modules\n";

#BEGIN {
    my @import_modules = (
    'Cwd',
    'File::chdir',
    'File::Copy',
    'POSIX',
    'Tie::File',
    'Try::Tiny',
    'Data::Dumper',
    'File::Basename',
    'Bio::SeqIO',
    'FindBin',
    'File::Remove',
    'Capture::Tiny',
    'File::Temp',
    'File::Spec::Functions',
    'Statistics::Multtest',
    'File::Path',
    'Statistics::Distributions',
    'Getopt::Long',
    'Statistics::R',
    'Math::Round',
    'File::Find',
    'Bio::DB::Taxonomy',
    'Pod::Usage',
    'Archive::Extract',
        );
#Lets commend it for now
=pod
    my ($inst) = ExtUtils::Installed->new();
    my (@installed_modules) = $inst->modules();
    for ( @import_modules ) {
        eval{ $inst->validate($_) };
        if($@) {
            print AOM::Module::fetchLocalTime()."Module $_ Not found! \n";
            #exit 1;
        } # end 'if'
        else { print AOM::Module::fetchLocalTime()."Module $_ OK\n";}
    } # end 'for'
#} # end 'BEGIN' block
=cut

print AOM::Module::fetchLocalTime()."Installing third party databases for Alienomics\n";
print "output folder: $third_party_DB/\n\n";


#---
#Install NCBI taxdump - T for taxonomy
if (index($flags,'T') >= 0) {
print "\n".AOM::Module::fetchLocalTime()."NCBI taxonomic database\n";
print "downloading taxdump\n"; 
chdir('$third_party_DB/taxdump') or die "$!"; # changing working directory
my $url = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz';
my $file = 'taxdump.tar.gz';
print "extracting taxdump 1\n";
my $code = getstore($url, $file); 
if (is_error($code)) { die AOM::Module::fetchLocalTime()."getstore of <$url> failed with $code\n";} 
else { print "$url successfully downloaded\n";}
print "extracting taxdump 2\n";
my $ae = Archive::Extract->new( archive => $file); 
	####print "extracting taxdump 3\n";
	####my $ok = $ae->extract( to => "$third_party_DB/taxdump" );
if ($ae->extract) { print "$file extracted successfully\n";} 
else { warn AOM::Module::fetchLocalTime()."ERROR: " . $ae->error;}
}

#----
#Install silvaDB https://www.arb-silva.de/no_cache/download/archive/release_138_1/Exports/ 
#LSU: Large subunit (23S/28S ribosomal RNAs) 
#SSU: Small subunit (16S/18S ribosomal RNAs)
# R for rRNA/Silva DB
if (index($flags,'R') >= 0) {
print "\n".AOM::Module::fetchLocalTime()."SILVA LSU Parc 138.1 database\n";
chdir('$third_party_DB/silvaDB/') or die "$!"; # changing working directory 
my $silva_url_LSU = 'https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_LSUParc_tax_silva.fasta.gz'; 
my $silva_LSU_fileName = 'SILVA_138.1_LSUParc_tax_silva.fasta.gz';
print "extracting silva 138.1\n";
my $silva_LSU_code = getstore($silva_url_LSU, $silva_LSU_fileName); 
if (is_error($silva_LSU_code)) { die AOM::Module::fetchLocalTime()."getstore of <$silva_url_LSU> failed with $silva_LSU_code\n";}
else { print "$silva_url_LSU successfully downloaded\n";}
print "extracting silva 138.1\n";
my $silva_LSU_ae = Archive::Extract->new( archive => $silva_LSU_fileName); 
	####print "extracting Silva 3\n";
	####my $silva_ok = $silva_LSU_ae->extract( to => "$third_party_DB/silvaDB/SILVA_138.1_LSUParc_tax_silva.fasta" );
if ($silva_LSU_ae->extract) { print "$silva_LSU_fileName extracted successfully\n";}
else { warn AOM::Module::fetchLocalTime()."ERROR: " . $silva_LSU_ae->error;}
}

#----
#Install buscoDB https://busco-archive.ezlab.org/
# B for BUSCO
if (index($flags,'B') >= 0) {
print "\n".AOM::Module::fetchLocalTime()."BUSCO database (metazoa odb9)\n";
chdir('$third_party_DB/') or die "$!"; # changing working directory
my $busco_url_LSU = 'https://busco-archive.ezlab.org/datasets/metazoa_odb9.tar.gz'; 
my $busco_LSU_fileName = 'metazoa_odb9.tar.gz';
print "downloading busco (metazoa odb9)\n";
my $busco_LSU_code = getstore($busco_url_LSU, $busco_LSU_fileName); 
if (is_error($busco_LSU_code)) { die AOM::Module::fetchLocalTime()."getstore of <$busco_url_LSU> failed with $busco_LSU_code\n";}
else { print "$busco_url_LSU successfully downloaded\n";}
print "extracting busco (metazoa odb9)\n";
my $busco_LSU_ae = Archive::Extract->new( archive => $busco_LSU_fileName); 
	####print "extracting busco 3\n";
	####my $busco_ok = $busco_LSU_ae->extract( to => "$third_party_DB/" );
if ($busco_LSU_ae->extract) { print "$busco_LSU_fileName extracted successfully\n";}
else { warn AOM::Module::fetchLocalTime()."ERROR: " . $busco_LSU_ae->error;}
}


#----
#Install Uniref50
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref50/
# U for Uniref50
if (index($flags,'U') >= 0) {


# section to keep !
print "\n".AOM::Module::fetchLocalTime()."Installing Uniref50 database\n\t\t\t(this step will likely take 2 to 3 hours)\n";
chdir('$third_party_DB/') or die "$!"; # changing working directory 
my $uniref_url_fasta = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref50/uniref50.fasta.gz'; 
my $uniref_fileName_fasta = 'uniref50.fasta.gz';
print "downloading uniref50 fasta file\n";
my $uniref_code_fasta = getstore ($uniref_url_fasta, $uniref_fileName_fasta); 
print "extracting uniref50 fasta file\n";
my $uniref_ae_fasta = Archive::Extract -> new (archive => $uniref_fileName_fasta);
	####print "extracting uniref fasta 3\n";
	####my $uniref_ok_fasta = $uniref_ae_fasta -> extract (to => "$third_party_DB/");
if ($uniref_ae_fasta->extract) { print "$uniref_fileName_fasta extracted successfully\n";}
else { warn AOM::Module::fetchLocalTime()."ERROR: " . $uniref_ae_fasta->error;}

my $uniref_url_idmapping = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz';
my $uniref_fileName_idmapping = 'idmapping_selected.tab.gz';
print "downloading uniref50 ID mapping file\n";
my $uniref_code_idmapping = getstore ($uniref_url_idmapping, $uniref_fileName_idmapping);
print "extracting uniref50 ID mapping file\n";
my $uniref_ae_idmapping = Archive::Extract -> new (archive => $uniref_fileName_idmapping);
####print "extracting uniref mapping 3\n";
####my $uniref_ok_idmapping = $uniref_ae_idmapping -> extract (to => "$third_party_DB/");
if ($uniref_ae_idmapping->extract) { print "$uniref_fileName_idmapping extracted successfully\n";}
else { warn AOM::Module::fetchLocalTime()."ERROR: " . $uniref_ae_idmapping->error;}


print "building Uniref50 taxlist\n";
build_uniref50_v2 ("uniref50.fasta", "idmapping_selected.tab", "uniref50.taxlist");
}
#end of uniref50 install



# useless section now ?

#my $uniref_url_xml = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref50/uniref50.xml.gz'; 
#my $uniref_fileName_xml = 'uniref50.xml.gz';
#print "downloading Uniref50 XML\n";
#my $uniref_code_xml = ""; 
#my $uniref_code_xml = getstore ($uniref_url_xml, $uniref_fileName_xml);  
#if (is_error($uniref_code_xml)) { die AOM::Module::fetchLocalTime()."getstore of <$uniref_url_xml> failed with $uniref_code_xml\n";}
#else { print "$uniref_url_xml successfully downloaded\n";}
#print AOM::Module::fetchLocalTime()."Extracting Uniref50 xml file (This will likely take several hours)\n";
#print AOM::Module::fetchLocalTime()."Extraction 1\n";
#my $uniref_ae_xml = Archive::Extract-> new ( archive => $uniref_fileName_xml); 
#print AOM::Module::fetchLocalTime()."Extraction 2\n";
#my $uniref_ok_xml = $uniref_ae_xml -> extract (to => "$third_party_DB/");
#if ($uniref_ae_xml->extract) { print AOM::Module::fetchLocalTime()."$uniref_fileName_xml extracted successfully\n";}
#else { warn AOM::Module::fetchLocalTime()."ERROR: " . $uniref_ae_xml->error;}




# Summary of install script

print "\n".AOM::Module::fetchLocalTime()."All third party databases are now ready to be used within Alienomics\n\n";
print "Use the following parameters in your configuration file for alienomics:\n";
print "rRNA_file =  $third_party_DB/silvaDB/SILVA_138.1_LSUParc_tax_silva.fasta\n" if (index($flags,'R') >= 0);
print "busco_file = $third_party_DB/metazoa_odb9/ancestral\n" if (index($flags,'B') >= 0);
print "uniref50_path = $third_party_DB/\n" if (index($flags,'U') >= 0);
print "taxdump_path = $third_party_DB/taxdump\n" if (index($flags,'T') >= 0);

exit; # end of Alienomics install script


##############################################
# building uniref50 from fasta and idmapping
sub build_uniref50_v2 {
my ($uniref_fasta, $uniref_idmapping, $uniref_out)=@_;
open(UF, '<', $uniref_fasta) or die $!;
open(UI, '<', $uniref_idmapping) or die $!;
open UO, ">$uniref_out" or die $!;

# linking SeqID to RepID
print AOM::Module::fetchLocalTime()."Parsing fasta file: $uniref_fasta\n";
my %repid = (); my %taxid_1 = (); my %seqid = ();
my $SeqID; my $TaxID_1; my $RepID;
while (my $line = <UF>) {
	# example: >UniRef50_A0A5A9P0L4 Peptidylprolyl isomerase n=1 Tax=Triplophysa tibetana TaxID=1572043 RepID=A0A5A9P0L4_9TELE
	#print "$line";
	if ($line =~ m/>(\S+)\s.*TaxID=(\S+)\sRepID=(\S+)/) {
	$SeqID=$1; $TaxID_1=$2; $RepID=$3;
	# print "\tSeqID=$SeqID\tRepID=$RepID\n";
	$repid{$SeqID} = $RepID;
	$taxid_1{$SeqID} = $TaxID_1;
	$seqid{$RepID} = $SeqID;
	}
	#else { print UO "\tNo pattern found\n";}
}

# linking RepID to taxID
print AOM::Module::fetchLocalTime()."Parsing id mapping file: $uniref_idmapping\n";
my $TaxID_2;
while (my $line = <UI>) {
	#print UO "\n\n$line";
	my @fields = split '\t', $line;
	$RepID=$fields[1]; $TaxID_2=$fields[12];
	 	
	if (exists $seqid{$RepID}) {
		#print UO "\tSeqID=$seqid{$RepID}\tTaxID_1=$taxid_1{$seqid{$RepID}}\tRepID=$RepID\tTaxID_2=$TaxID_2\n";
		#print UO "$seqid{$RepID}\t$TaxID_2\n";
		$taxid_1{$seqid{$RepID}} = $TaxID_2;
	}

}

## Outputing taxonomic info
for(keys %taxid_1){
	#print UO "key=$_\tvalue=$taxid_2{$repid{$_}}\n";
	print UO "$_ $taxid_1{$_}\n";
}

} # end of build_uniref50_v2 subroutine








##########################################
#Reformat and build uniref50 from xml file
sub build_uniref50 { 
my ($uniref_in, $uniref_out, $loc)=@_; 
my $entry; my $taxon; my $flag=0; 
open(FH, '<', $uniref_in) or die $!; 
open TX, ">$uniref_out" or die $!; 
while (my $line = <FH>) {
        if ($line =~ /<entry id="(\S+)"/) {
                $entry=$1;
                $taxon="";
                $flag=0;
                #print "entry id $1 $line";
        }
        elsif ($line =~ /"common taxon ID" value="(\S+)"/) {
                $taxon=$1;
                $flag=1;
                #print "common ID $1 $line";
                print TX "$entry $taxon\n";
        }
        elsif ($line =~ /"NCBI taxonomy" value="(\S+)"/ and $flag == 0) {
                $taxon=$1;
                $flag=1;
                #print "NCBI $1 $line";
                print TX "$entry $taxon\n";
        }
}
};

}
1;
