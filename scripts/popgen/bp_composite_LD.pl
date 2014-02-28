#!/usr/bin/perl
# -*-Perl-*- 

=head1 NAME 

bp_composite_LD -i filename.prettybase.txt --sortbyld E<gt> outfile

=head1 SYNOPSIS

  bp_composite_LD -i filename.prettybase [-o out.LD] [-f prettybase/csv] [--sortbyld] [--noconvertindels]

=head2 DESCRIPTION

This a script which allows an easy way to calculate composite LD.  

=head2 OPTIONS

-i or --in     filename

-f or --format genotype format (prettybase or CSV)

--sortbyld     To see data sorted by LD instead of just all the 
               site1/site2 pair LD values.

-o or --out    output filename, otherwise will print to STDOUT

--noconvert    (applicable for prettybase format file only)
               if specified will NOT attempt to convert indel
               states to 'I' and delete states ('-') to 'D'.

-h or --help   see this documentation

=head2 AUTHOR Jason Stajich, Matthew Hahn

For more information contact:

Matthew Hahn, E<lt>matthew.hahn-at-duke.eduE<gt>
Jason Stajich E<lt>jason-at-bioperl-dot-orgE<gt>


=cut

use strict;
use warnings;

use Bio::PopGen::IO;
use Bio::PopGen::Statistics;
use Bio::PopGen::Population;

use Getopt::Long;

my ($file,$outfile,$sortbyld,$format,$noconvert,$verbose);
$format = 'prettybase'; # default format is prettybase
GetOptions(
	   'i|in:s'       => \$file, # pass the filename as 
	   'o|out:s'      => \$outfile,
	   'f|format:s'   => \$format,
	   'sortbyld'     => \$sortbyld,
	   'noconvert'    => \$noconvert,
	   'v|verbose'    => \$verbose,
	   'h|help'       => sub { system('perldoc', $0);
				   exit; }, 
	   );

if( ! $file ) { 
    $file = shift @ARGV;  # if no -i specified
}

my $io = Bio::PopGen::IO->new(-format => $format,
			      -verbose=> $verbose,
			      -CONVERT_INDEL_STATES => ! $noconvert,
			      -file   => $file);

my $stats = Bio::PopGen::Statistics->new(-verbose => $verbose);
my $pop = $io->next_population;

my %LD_matrix = $stats->composite_LD($pop);

# sites can be ordered by sorting their names

my @sites;
my $out;
if( $outfile ) { 
    open $out, '>', $outfile or die "Could not write file '$outfile': $!\n";
} else { 
    $out = \*STDOUT;
}
foreach my $site1 ( sort keys %LD_matrix ) {
    foreach my $site2 ( sort keys %{ $LD_matrix{$site1} } ) {
	my $LD = $LD_matrix{$site1}->{$site2}; # LD for site1,site2 
	if( $sortbyld ) {
	    push @sites, [ $site1,$site2,@$LD];
	} else { 
	    printf $out "%s,%s - LD=%.4f chisq=%.4f\n",$site1,$site2,@$LD;
	}
    }
}

if( $sortbyld ) {
    foreach my $s ( sort { $b->[3] <=> $a->[3] } @sites ) {
	my ($site1,$site2,$ld,$chisq) = @$s;
	print $out "$site1,$site2 - LD=$ld, chisq=$chisq\n";
    }
}
