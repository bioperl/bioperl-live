#!/usr/local/bin/perl


=head1 NAME

bpindex.pl - indexes files for use by bpfetch

=head1 SYNOPSIS

  bpindex.pl index_name file1 file2 etc

=head1 DESCRIPTION

bpindex.pl builds a bioperl index for the sequence files given in the
argument list, under the index name. For example

   bpindex.pl nrdb /data/nrdb/nrdb.fasta

would build an index called 'nrdb' as the index name for the file
nrdb.fasta, and

   bpindex.pl -fmt EMBL swiss /data/swiss/*.dat

would build an index called swiss for all the files in /data/swiss
which end in .dat which are in EMBL format.

The indexes are build using the Bio/Index/* modules, in particular,
Bio::Index::EMBL and the Bio::Index::Fasta modules. Any script which
uses these modules can use the index. A good example script is bpfetch
which fetches sequences and pipes them to stdout, for example

   bpfetch swiss:ROA1_HUMAN

gets the ROA1_HUMAN sequence from the swiss index and writes it as
fasta format on stdout

=head1 OPTIONS

  --fmt <format> - Fasta (default), or EMBL
  --dir <dir>    - directory to find the index files

=head1 ENVIRONMENT

bpindex and bpfetch coordinate where the databases lie using the
enviroment variable BIOPERL_INDEX. This can be overridden using the
--dir option

=head1 USING IT YOURSELF

bpindex.pl is a script that drives the Index modules. If you want to 
use this script heavily in your work, if it is Perl based, it is 
almost certainly better to look at the code in this script and copy
it across (probably you will be more likely to want to use the bpfetch
code).

=head1 EXTENDING IT

bpindex is just a wrapper around James Gilbert's excellent Index modules
found in bioperl

=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

    vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
    vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
    http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via email
or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR

Ewan Birney, birney@sanger.ac.uk

=cut

use strict;
use Getopt::Long;

#
# Dofus catcher for people who are trying this script without
# installing bioperl
#

BEGIN {
    eval {
	require Bio::Index::Fasta;
	require Bio::Index::EMBL;
    };
    if ( $@ ) {
	# one up from here is Bio directory - we hope!
	push(@INC,"..");
	eval {
	    require Bio::Index::Fasta;
	    require Bio::Index::EMBL;
	};
	if ( $@ ) {
	    print STDERR ("\nbpindex cannot find Bio::Index::Fasta and Bio::Index::EMBL\nbpindex needs to have bioperl installed for it to run.\nBioperl is very easy to install\nSee http://bio.perl.org for more information\n\n");
	    exit(1);
	} else {
	    print STDERR ("\nYou are running bpindex.pl without installing bioperl.\nYou have done it from bioperl/scripts, and so we can find the necessary information\nbut it is much better to install bioperl\n\nPlease read the README in the bioperl distribution\n\n");
	}
    }
}

my $name = shift;

if( !$name ) {
    system("perldoc $0");
    exit(1);
}

my @files = @ARGV;
my $dir = $ENV{'BIOPERL_INDEX'};
my $fmt = 'Fasta';
my $ret = GetOptions('dir=s' => \$dir,'fmt=s' => \$fmt );

if( !$dir ) {
    print STDERR "\nNo directory specified for index\nDirectory must be specified by the environment varaible BIOPERL_INDEX or --dir option\ngo bpindex with no arguments for more help\n\n";
    exit(1);
}

#
# Rock and roll...
# 
my $index;
$_ = $fmt;
SWITCH : {
    /Fasta/ && do {
	$index = Bio::Index::Fasta->new("$dir/$name", 'WRITE');
	last;
    };
    /EMBL/ && do {
	$index = Bio::Index::EMBL->new("$dir/$name", 'WRITE');
	last;
    };
    die("No index format called $fmt");
}

$index->make_index(@files);

# finished. Neat eh.

# check things close/deallocate etc.
$index = undef;


	
