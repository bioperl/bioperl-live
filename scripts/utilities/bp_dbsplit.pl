#!perl
#-*-Perl-*-

=head1 NAME

bp_dbsplit - script to split an input set of database(s) into smaller pieces

=head1 SYNOPSIS

  bp_dbsplit.PLS --size 50 [-i inputfile] [-if inputformat] [-of outputformat]
              [--prefix outputprefix] [ < file1 file 2  OR file1 file2]

=head1 DESCRIPTION

This script will take as input a list of filenames or a single file or
from STDIN a sequence database and split the database into separate
files of X numbers of sequences.  You specify X with the C<--size/-s>
parameter.  The input and output sequence format is any that is
supported by bioperl (fasta,embl,genbank,gcg, swissprot, etc).

You can specify the input data either as a single file with -i
filename, or as a single file as an argument like

  % bp_dbsplit file1 file2

or as a list of sequence data with 

  % cat file1 file2 file3 | bp_dbsplit

You'll want to use the C<--prefix> to specify what the output prefix will
be.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Jason Stajich, jason-at-bioperl-dot-org

=cut

use strict;
use warnings;
use Bio::SeqIO;
use Bio::SeqIO::MultiFile;

use Getopt::Long;
my $dbsize = 100;
my $prefix;
my ($informat,$outformat,$infile) = ( 'fasta', 'fasta');

GetOptions (
	    's|size:s'     => \$dbsize,
	    'if:s'         => \$informat,
	    'of:s'         => \$outformat,
	    'i:s'          => \$infile,
	    'p|prefix:s'   => \$prefix,
	    
);
if( @ARGV == 1 ) {
    $infile = shift @ARGV;
}
$prefix ||= $infile || $ARGV[0] || 'db';

my $in;
if( @ARGV ) {
    $in = new Bio::SeqIO::MultiFile(-files => [@ARGV],
				    -format => $informat || 'fasta');
} elsif( $infile ) {
    $in = new Bio::SeqIO(-file  => $infile,
			 -format=> $informat);
} else { 
    $in = new Bio::SeqIO(-format=> $informat);
}
my $count = 1;
my $out = new Bio::SeqIO(-format => $outformat,
			 -file   => ">$prefix.$count");
my $scount = 0;
while( my $seq = $in->next_seq ) {    
    if( ++$scount > $dbsize && $count ) { 
	$out->close();
	undef($out);
	$count++;
	$out = new Bio::SeqIO(-format => $outformat,
			      -file   => ">$prefix.$count");
	$scount = 1;
    }
    $out->write_seq($seq);
}


__END__
