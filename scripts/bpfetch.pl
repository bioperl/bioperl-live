#!/usr/local/bin/perl


=head1 NAME

bpfetch.pl - fetches sequences from bioperl indexed databases

=head1 SYNOPSIS 

  bpfetch.pl swiss:ROA1_HUMAN

  bpfetch.pl net::genbank:X47072

  bpfetch.pl net::genpept:ROA1_HUMAN

  bpfetch.pl -fmt GCG swiss:ROA1_HUMAN

=head1 DESCRIPTION

Fetches sequences using the DB access systems in Bioperl. The most
common use of this is to fetch sequences from bioperl indices built
using bpindex.pl, or to fetch sequences from the NCBI website

=head1 OPTIONS

  -fmt  <format> - Output format
                   Fasta (default), EMBL, Raw or GCG
  -acc           - string is an accession number, not an
                   id. 

options only for expert use

  -dir  <dir>    - directory to find the index files
                  (overrides BIOPERL_INDEX environment varaible)
  -type <type>   - type of DBM file to open 
                  (overrides BIOPERL_INDEX_TYPE environment variable)

=head1 ENVIRONMENT

bpindex and bpfetch coordinate where the databases lie using the
enviroment variable BIOPERL_INDEX. This can be overridden using the
-dir option. The index type (SDBM or DB_File or another index file)
is controlled by the BIOPERL_INDEX_TYPE variable. This defaults to 
SDBM_File 

=head1 USING IT YOURSELF

bpfetch is a wrapper around the bioperl modules which support 
the Bio::DB::BioSeqI abstract interface. These include:

  Author          Code

  James Gilbert - Fasta indexer, Abstract indexer
  Aaron Mackay  - GenBank and GenPept DB access
  Ewan Birney   - EMBL .dat indexer
  Many people   - SeqIO code

These modules can be used directly, which is far better than using
this script as a system call or a pipe to read from. Read the
source code for bpfetch to see how it is used.

=head1 EXTENDING IT

bpfetch uses a number of different modules to provide access to
databases. Any module which subscribes to the Bio::DB::BioSeqI
interface can be used here. For flat file indexers, this is
best done by extending Bio::Index::Abstract, as is done in
Bio::Index::EMBL and Bio::Index::Fasta. For access to other
databases you will need to roll your own interface.

For new output formats, you need to add a new SeqIO module. The
easiest thing is to look at Bio::SeqIO::Fasta and figure out
how to hack it for your own format (call it something different
obviously).

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
# installing bioperl. In your own script, you can just go
#
# use Bio::Index::Fasta etc, rather than this
#

BEGIN {
    eval {
	require Bio::Index::Fasta;
	require Bio::Index::EMBL;
        require Bio::DB::GenBank;
        require Bio::DB::GenPept;
        require Bio::SeqIO;
	
    };
    if ( $@ ) {
	# one up from here is Bio directory - we hope!
	push(@INC,"..");
	eval {
	    require Bio::Index::Fasta;
	    require Bio::Index::EMBL;
            require Bio::DB::GenBank;
            require Bio::DB::GenPept;
            require Bio::SeqIO;
	};
	if ( $@ ) {
	    print STDERR ("\nbpindex cannot find Bio::Index::Fasta and Bio::Index::EMBL\nbpindex needs to have bioperl installed for it to run.\nBioperl is very easy to install\nSee http://bio.perl.org for more information\n\n");
	    exit(1);
	} else {
	    print STDERR ("\nYou are running bpindex.pl without installing bioperl.\nYou have done it from bioperl/scripts, and so we can find the necessary information\nbut it is much better to install bioperl\n\nPlease read the README in the bioperl distribution\n\n");
	}
    }
}

#
# Start processing the command line
#

my $dir = $ENV{'BIOPERL_INDEX'};
my $type = $ENV{'BIOPER_INDEX_TYPE'};
my $fmt = 'Fasta';
my $useacc = 0;
my $ret = GetOptions('dir=s' => \$dir,'fmt=s' => \$fmt , 'type=s' => \$type , 'acc!' => \$useacc);

#
# print pod documentation if we have no arguments
#

exec('perldoc',$0) unless @ARGV;

my($isnet,$db,$id,$seq,$seqio,$out);

#
# Reset the type if needed
#

if( $type ) {
   $Bio::Index::Abstract::USE_DBM_TYPE = $type;
}

#
# Build at run time the SeqIO output
#

$out = Bio::SeqIO->new(-fh => \*STDOUT , -format => $fmt);

#
# Main loop over remaining arguments
#

foreach my $arg ( @ARGV ) {
    $_= $arg;

    # strip out net:: if there
    if( /^net::/ ) {
	$isnet = 1;
	s/^net:://;
    } else {
	$isnet = 0;
    }

    # parse to db:id 

    /^(\w+)\:(\S+)$/ || do { print STDERR "$_ is not parsed as db:name\n"; next;};
    $db = $1;
    $id = $2;

    #
    # the eval block catches exceptions if they occur
    # in the code in the block. The exception goes in $@
    #

    eval {
	if( $isnet == 1 ) {
	    if( $db =~ /genbank/ ) {
		$db = Bio::DB::GenBank->new();
	    }
	    if( $db =~ /genpept/ ) {
		$db = Bio::DB::GenPept->new();
	    }
	} else {
	    if( !$dir ) {
		
		print STDERR "\nNo directory specified for index\nDirectory must be specified by the environment varaible BIOPERL_INDEX or --dir option\ngo bpindex with no arguments for more help\n\n";
		exit(1);
	    }

	    #
	    # $db gets re-blessed to the correct index when
	    # it is made from the abstract class. Cute eh?
	    #

	    $db = Bio::Index::Abstract->new("$dir/$db");

	}
    };
    if( $@ ) {
	warn("Database $db in $arg is not loadable. Skipping\n\nError $@");
	next;
    }

    #
    # We expect the databases to adhere to the BioSeqI
    # the sequence index databases and the GenBank/GenPept do already
    #

    if( ! $db->isa('Bio::DB::BioSeqI') ) {
	warn("$db in $arg does not inherit from Bio::DB::BioSeqI, so is not expected to work under the DB guidlines. Going to try it anyway");
    }

    eval {
	if( $useacc == 0 ) {
	    $seq = $db->get_Seq_by_id($id);
	} else {
	    $seq = $db->get_Seq_by_acc($id);
	}

    };
    if( $@ ) {
	warn("Sequence $id in Database $db in $arg is not loadable. Skipping.\n\nError $@");
	next;
    }

    $out->write_seq($seq);
}


	




