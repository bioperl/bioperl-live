#!perl

=head1 NAME

bp_tree2pag - convert Bio::TreeIO parseable format trees to pagel format

=head1 SYNOPSIS

 bp_tree2pag -if nexus -i file.nexus > file.pag

 # OR pipe in through STDIN, and use newick format instead

 cat file.newick | bp_tree2pag -if newick > file.nh

 # OR specify an output and input

 bp_tree2pag -o file.pag -i file.newick

=head1 DESCRIPTION

Convert TreeIO parseable files into Pagel format tree files.  Be
warned that pagel format only really supports a single tree per file
so.  Also Pagel parsing is not yet available in bioperl.

=cut

use strict;
use warnings;
use Bio::TreeIO;
use Getopt::Long;
my ($iformat,$oformat) = ('newick', 'pag');
my ($outfile,$infile);
GetOptions(
	   'if|informat:s'    => \$iformat,
	   'of|outformat:s'   => \$oformat,
	   'i|in:s'           => \$infile,
	   'o|out:s'          => \$outfile,
	   'h|help'           => sub { exec('perldoc', $0);
				       exit(0); },
	   );
my $in;
if( ! $infile ) {
    $in = Bio::TreeIO->new(-format => $iformat,
			   -fh     => \*ARGV);
} else { 
    $in = Bio::TreeIO->new(-format => $iformat,
			   -file   => $infile);
}

my $out;
if( $outfile) {
    $out = Bio::TreeIO->new(-format => $oformat,
			    -file   => ">$outfile");
} else { 
    $out = Bio::TreeIO->new(-format => $oformat); #print to STDOUT instead
}

while( my $t = $in->next_tree ) {
    $out->write_tree($t);
}
