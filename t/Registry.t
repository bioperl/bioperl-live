##-*-Perl-*-
# $Id$
# test for Bio::DB::Registry
# written by Brian Osborne

use strict;
my $old_search_path;

BEGIN {
   eval { require Test; };
   if( $@ ) {
      use lib 't','..';
   }
   use Test;
   plan tests => 6;
   $old_search_path = $ENV{OBDA_SEARCH_PATH} if defined $ENV{OBDA_SEARCH_PATH};
   $ENV{OBDA_SEARCH_PATH} = "t/data";
}

use Bio::DB::Registry;
use Bio::DB::Flat;
use Bio::Root::IO;

my $tmpdir = "t/tmp";
my $dbname = "testflat";
mkdir($tmpdir,0777);

my $registry = Bio::DB::Registry->new;
ok defined($registry);
my @available_services = $registry->services;
ok $available_services[1] eq "testflat";
my $flat = Bio::DB::Flat->new(-directory  => $tmpdir,
			      -dbname     => $dbname,
			      -index      => "flat",
                              -write_flag => 1,
                              -format     => "fasta" );
my $entries = $flat->build_index("t/data/cysprot.fa");
ok $entries == 7;
my $db = $registry->get_database($dbname);
ok defined($db);
my $seq = $db->get_Seq_by_id("ALEU_HORVU");
ok defined($seq);
ok ($seq->seq) eq "MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLGALGRTRHALRFARFAVRYGKSYESAAEVRRRFRIFSESLEEVRSTNRKGLPYRLGINRFSDMSWEEFQATRLGAAQTCSATLAGNHLMRDAAALPETKDWREDGIVSPVKNQAHCGSCWTFSTTGALEAAYTQATGKNISLSEQQLVDCAGGFNNFGCNGGLPSQAFEYIKYNGGIDTEESYPYKGVNGVCHYKAENAAVQVLDSVNITLNAEDELKNAVGLVRPVSVAFQVIDGFRQYKSGVYTSDHCGTTPDDVNHAVLAVGYGVENGVPYWLIKNSWGADWGDNGYFKMEMGKNMCAIATCASYPVVAA";

sub cleanup {
   eval {
      Bio::Root::IO->rmtree($tmpdir) if (-d $tmpdir);
   };
   $ENV{OBDA_SEARCH_PATH} = $old_search_path;
}

END { cleanup(); }
