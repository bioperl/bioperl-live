#!/usr/local/bin/perl


=head1 NAME

bpfetch.pl - fetches sequences from bioperl indexed databases

=head1 SYNOPSIS 

  bpfetch.pl swiss:ROA1_HUMAN

  bpfetch.pl net::genbank:X47072

  bpfetch.pl net::genpept:ROA1_HUMAN

  bpfetch.pl --fmt GCG swiss:ROA1_HUMAN

=head1 DESCRIPTION

Fetches sequences using the DB access systems in Bioperl. The most
common use of this is to fetch sequences from bioperl indices built
using bpindex.pl

=head1 OPTIONS

  --fmt <format> - Output format
                   Fasta (default), EMBL, raw or GCG
  --dir <dir>    - directory to find the index files

=head1 ENVIRONMENT

bpindex and bpfetch coordinate where the databases lie using the
enviroment variable BIOPERL_INDEX. This can be overridden using the
--dir option

=head1 USING IT YOURSELF

bpfetch is a wrapper around the bioperl modules which support 
the Bio::DB::BioSeqI abstract interface. These include:

  James Gilbert - Fasta indexer
  Ewan Birney   - EMBL .dat indexer
  Aaron Mackay  - GenBank and GenPept DB access

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



my @files = @ARGV;
my $dir = $ENV{'BIOPERL_INDEX'};
my $fmt = 'Fasta';
my $ret = GetOptions('dir=s' => \$dir,'fmt=s' => \$fmt );

if( $#ARGV == -1 ) {
    system("perldoc $0");
    exit(1);
}

my($isnet,$db,$id,$seq,$seqio,$out);

$out = Bio::SeqIO->new(-fh => \*STDOUT , -format => $fmt);

foreach my $arg ( @ARGV ) {
    $_= $arg;
    if( /^net::/ ) {
	$isnet = 1;
	s/^net:://;
    } else {
	$isnet = 0;
    }


    /^(\w+)\:(\S+)$/ || do { print STDERR "$_ is not parsed as db:name\n"; next;};
    $db = $1;
    $id = $2;

    eval {
	if( $isnet == 1 ) {
	    $db =~ /genbank/ && do {
		$db = Bio::DB::GenBank->new();
	    };
	    $db =~ /genpept/ && do {
		$db = Bio::DB::GenPept->new();
	    };
	} else {
	    if( !$dir ) {
		
		print STDERR "\nNo directory specified for index\nDirectory must be specified by the environment varaible BIOPERL_INDEX or --dir option\ngo bpindex with no arguments for more help\n\n";
		exit(1);
	    }

	    $db = Bio::Index::Abstract->new("$dir/$db");

	    #
	    # $db gets re-blessed to the correct index when
	    # it is made from the abstract class. Cute eh?
	    #
	}
    };
    if( $@ ) {
	warn("Database $db in $arg is not loadable. Skipping\n\nError $@");
	next;
    }
    if( ! $db->isa('Bio::DB::BioSeqI') ) {
	warn("$db in $arg does not inheriet from Bio::DB::BioSeqI, so is not expected to work under the DB guidlines. Going to try it anyway");
    }

    eval {
	$seq = $db->get_Seq_by_id($id);
    };
    if( $@ ) {
	warn("Sequence $id in Database $db in $arg is not loadable. Skipping.\n\nError $@");
	next;
    }

    $out->write_seq($seq);
}


	
