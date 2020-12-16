#!/usr/bin/env perl
use strict;
use warnings;
use Text::CSV_XS;
use autodie;

package AOM::plotHandler;
sub KmeanPlot {
my ($resFile)=shift;
open my $fh, "<", $resFile or die "$resFile: $!";
my $csv = Text::CSV_XS->new( { sep_char => "\t" } );
my @list;
$csv->column_names ($csv->getline (*$fh));
while ( my $hr = $csv->getline_hr(*$fh) ) {
    push @list, $hr->{'globalScore'};
}

=pod
my $R = Statistics::R->new();
$R->set( 'filenames', \@filenames );
$R->set( 'id',        $id );
$R->set( 'par1',      $par1 );
$R->set( 'par2',      $par2 );
$R->run_from_file( "genoplot_by_id.build_df.R" );
$R->run_from_file( "genoplot_by_id.build_plot.R" );
$R->run( qq`setwd("$plot_dir")` );
$R->run(
    qq`ggsave(
      filename = paste("$plot_path", "$plot_format", sep = "."),
      plot     = geno.plot,
      width    = $plot_width,
      height   = $plot_height)`
);
=cut
}
1;
