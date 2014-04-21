#!/usr/bin/perl

=head1 NAME

bp_fetch.pl - fetches sequences from bioperl indexed databases

=head1 SYNOPSIS

  bp_fetch.pl swiss:ROA1_HUMAN

  bp_fetch.pl net::genbank:JX295726

  bp_fetch.pl net::genpept:ROA1_HUMAN

  bp_fetch.pl ace::myserver.somewhere.edu,21000:X56676

  bp_fetch.pl -fmt GCG swiss:ROA1_HUMAN

=head1 DESCRIPTION

Fetches sequences using the DB access systems in Bioperl. The most
common use of this is to bp_fetch sequences from bioperl indices built
using bpindex.pl, or to fetch sequences from the NCBI website

The format for retrieving sequences is delibrately like the
GCG/EMBOSS format like the following:

  db:name

with the potential of putting in a 'meta' database type, being

  meta::db:name

The meta information can be one of three types

  local - local indexed flat file database
  net   - networked http: based database
  ace   - ACeDB database

This information defaults to 'local' for database names with no meta
db information

=head1 OPTIONS

  -fmt  <format> - Output format
                   Fasta (default), EMBL, Raw, swiss or GCG
  -acc           - string is an accession number, not an
                   id.

options only for expert use

  -dir  <dir>    - directory to find the index files
                  (overrides BIOPERL_INDEX environment varaible)
  -type <type>   - type of DBM file to open
                  (overrides BIOPERL_INDEX_TYPE environment variable)

=head1 ENVIRONMENT

bp_index and bp_fetch coordinate where the databases lie using the
environment variable BIOPERL_INDEX. This can be overridden using the
-dir option. The index type (SDBM or DB_File or another index file)
is controlled by the BIOPERL_INDEX_TYPE variable. This defaults to
SDBM_File

=head1 USING IT YOURSELF

bp_fetch is a wrapper around the bioperl modules which support
the Bio::DB::BioSeqI abstract interface. These include:

  Author          Code

  James Gilbert - Fasta indexer, Abstract indexer
  Aaron Mackay  - GenBank and GenPept DB access
  Ewan Birney   - EMBL .dat indexer
  Many people   - SeqIO code

These modules can be used directly, which is far better than using
this script as a system call or a pipe to read from. Read the
source code for bp_fetch to see how it is used.

=head1 EXTENDING IT

bp_fetch uses a number of different modules to provide access to
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
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Ewan Birney E<lt>birney@ebi.ac.ukE<gt>

=cut

use strict;
use warnings;
use Getopt::Long;

#
# Dofus catcher for people who are trying this script without
# installing bioperl. In your own script, you can just go
#
# use Bio::Index::Fasta etc, rather than this
#
use Bio::SeqIO;
BEGIN {
  eval {
    require Bio::Index::Fasta;
    require Bio::Index::EMBL;
    require Bio::Index::GenBank;
    require Bio::Index::Swissprot;
    require Bio::Index::SwissPfam;
  };
  if ( $@ ) {
    # one up from here is Bio directory - we hope!
    push(@INC,"..");
    eval {
      require Bio::Index::Fasta;
      require Bio::Index::EMBL;
      require Bio::Index::GenBank;
      require Bio::Index::Swissprot;
      require Bio::Index::SwissPfam;
    };
    if ( $@ ) {
      print STDERR ("\nbp_index cannot find Bio::Index::Fasta and Bio::Index::EMBL\nbp_index needs to have bioperl installed for it to run.\nBioperl is very easy to install\nSee http://bio.perl.org for more information\n\n");
      exit(1);
    } else {
      print STDERR ("\nYou are running bp_index.pl without installing bioperl.\nYou have done it from bioperl/scripts, and so we can find the necessary information\nbut it is much better to install bioperl\n\nPlease read the README in the bioperl distribution\n\n");
    }
  }

  eval {
    require Bio::DB::GenBank;
    require Bio::DB::GenPept;
    require Bio::DB::EMBL;
    require Bio::DB::SwissProt;
  };

  if ( $@ ) {
    if ( !exists $ENV{'BIOPERL_SAVVY'} ) {
      print STDERR ("\nbp_fetch cannot find Bio::DB::GenBank and Bio::DB::EMBL modules\nThis is most likely because LWP has not been installed\nThis does not effect local indexing\nset environment variable BIOPERL_SAVVY to supress this message\n\n");
    }
  }
}

#
# Start processing the command line
#

my $dir = $ENV{'BIOPERL_INDEX'};
my $type = $ENV{'BIOPERL_INDEX_TYPE'};
my $fmt = 'Fasta';
my $useacc = 0;
my $ret = GetOptions('d|dir=s' => \$dir,
		     'f|fmt=s' => \$fmt ,
		     't|type=s' => \$type ,
		     'acc!' => \$useacc);

#
# print pod documentation if we have no arguments
#

exec('perldoc',$0) unless @ARGV;

my ($isnet,$db,$dbobj,$id,$seq,$seqio,$out,$meta);

#
# Reset the type if needed
#

if( $type ) {
   $Bio::Index::Abstract::USE_DBM_TYPE = $type;
}

#
# Build at run time the SeqIO output
#
if ( $fmt !~ /swisspfam|pfam/ ) {
  $out = Bio::SeqIO->new(-fh => \*STDOUT , -format => $fmt);
}

#
# Main loop over remaining arguments
#

for my $arg ( @ARGV ) {
  $_= $arg;
  # strip out meta:: if there
  if ( /^(\w+)::/ ) {
    $meta = $1;
    s/^(\w+):://;
  } else {
    $meta = 'local';
  }

  # parse to db:id

  /^(\S+)\:(\S+)$/ || do { warn "$_ is not parsed as db:name\n"; next; };
  ($db,$id) = split/:/,$_,2;
  #
  # the eval block catches exceptions if they occur
  # in the code in the block. The exception goes in $@
  #

  eval {
    SWITCH : {
      $_ = $meta;
      /^net$/ && do {
	if ( $db =~ /genbank/i ) {
	  $dbobj = Bio::DB::GenBank->new(-format => $fmt);
	} elsif ( $db =~ /genpept/i ) {
	  $dbobj = Bio::DB::GenPept->new();
	} elsif ( $db =~ /embl/i ) {
	  $dbobj = Bio::DB::EMBL->new();
	} else {
	  die "Net database $db not available";
	}
	last SWITCH;
      };
      /^ace$/ && do {
	# yank in Bio::DB::Ace at runtime
	eval {
	  require Bio::DB::Ace;
	};
	if ( $@ ) {
	  die "Unable to load Bio::DB::Ace for ace::$db\n\n$@\n";
	}

	# db is server,port
	my ($server,$port);

	$db =~ /(\S+)\,(\d+)/ || die "$db is not server.name,port for acedb database";
	$server = $1;
	$port = $2;
	# print STDERR "Connecting to $server,$port\n";

	$dbobj = Bio::DB::Ace->new(-host => $server, -port => $port);
	last SWITCH;
      };
      /^local$/ && do {
	if ( !$dir ) {
	  die "\nNo directory specified for index\nDirectory must be specified by the environment varaible BIOPERL_INDEX or --dir option\ngo bp_index with no arguments for more help\n\n";
	}

	#
	# $db gets re-blessed to the correct index when
	# it is made from the abstract class. Cute eh?
	#

	$dbobj = Bio::Index::Abstract->new("$dir/$db");
	last SWITCH;
      };
      die "Meta database $meta is not valid";
    }
    };				# end of eval to get db
  if ( $@ ) {
    warn("Database $db in $arg is not loadable. Skipping\n\nError $@");
    next;
  }

  #
  # We expect the databases to adhere to the BioSeqI
  # the sequence index databases and the GenBank/GenPept do already
  #
  if ( $dbobj->isa("Bio::Index::SwissPfam") ) {
    my $seq = $dbobj->fetch($id);
    if ( $seq ) {
      my $started;
      while ( <$seq> ) {
	last if ( /^\s+$/ );
	print;
      }
    } else {
      warn("Cannot find $id\n");
    }
    next;
  }
  if ( ! $dbobj->isa('Bio::DB::RandomAccessI') ) {
    warn("$db in $arg does not inherit from Bio::DB::RandomAccessI, so is not expected to work under the DB guidlines. Going to try it anyway");
  }
  eval {
    if ( $useacc == 0 ) {
      $seq = $dbobj->get_Seq_by_id($id);
    } else {
      $seq = $dbobj->get_Seq_by_acc($id);
    }
  };
  if ( $@ ) {
    warn("Sequence $id in Database $db in $arg is not loadable. Skipping.\n\nError $@");
    next;
  } elsif ( !defined $seq ) {
    warn("Sequence $id in Database $db is not present\n");
    next;
  }
  $out->write_seq($seq);
}
