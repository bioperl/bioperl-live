#!/usr/bin/perl

# This script will report the names of the tissues which were seen
# in a BLAST/FASTA report against an EST (or cDNA possibly) library.

# this script assumes you have a directory which you have downloaded
# gbestXX.seq.gz from ncbi genbank release.  This will run faster if 
# they are uncompressed already, but if will uncompress the files 
# on demand.  Be sure that there is sufficient space and the uid
# has write permission on the files and in that directory if you
# plan to run this script on compressed files.

# Alternatively you can use this with the -r option and it will 
# use the remote sequence databases either genbank or embl to 
# retrieve the specific EST (so one can attempt to guess the tissue type)
##
# cmd line options are
# -i/--index=indexname
# -d/--dir=dir where gbest data files are located
# -b/--blast=filename blast filename which compared against an EST db
# -f/--format=(blast|blastxml|fasta) - type of search output either 
#                                      from BLAST or FASTA suites 
# -c/--cache=filename cache for accession number to tissue
# -p/pvalue=pvalue pvalue to limit search to
# -r/--remote=[GenBank|EMBL] use remote db for searching

use strict;
use DB_File;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::DB::EMBL;
use Bio::DB::GenBank;
use Bio::Index::GenBank;

use Getopt::Long;
my $GZIP = '/usr/bin/gzip';
my $GUNZIP = '/usr/bin/gunzip';

my $dir = '/home/data/libs/gbest'; # local dir for gbest files
my $index = 'dbest_tissue.idx';    # local index filename
my $cache;      # filename to create cache of accession->tissue
my $VERBOSE = 0;# verbosity option
my $blastfile;  # blastfile to parse
my $pvalue;     # Max P-Value allowed when parsing blastfile
my $remote;     # flag for remote database 
my $db;         # generic database handle
my %accessions; # cache results
my $format = 'blast';

&GetOptions( 'd|dir:s'   => \$dir,
	     'i|index:s' => \$index,
	     'v|verbose' => \$VERBOSE,
	     'b|blast:s' => \$blastfile,
	     'f|format:s' => \$format,
	     'c|cache:s' => \$cache,
	     'p|pvalue:s'       => \$pvalue,
	     'r|remote:s'=> \$remote);

if( $cache && -w $cache ) {
    print "creating cache file\n";
    tie %accessions, "DB_File", $cache,  O_RDWR|O_CREAT,0660, $DB_HASH;
}

if( ! $remote ) {
    opendir(GBEST, $dir) or die("cannot open $dir");
    
    my $indexfile = new Bio::Index::GenBank(-filename   => $index,
					    -write_flag => 'WRITE');
    foreach my $file  ( readdir(GBEST) ) {
#	print "file is $file\n";
	    next unless ( $file =~ /(gbest\d+\.seq)(.gz)?$/ );
	    if( $2 ) {		
		`$GUNZIP $dir/$file`;
	    }
	    $indexfile->make_index("$dir/$1");
    }

    $indexfile = undef;
    $db = new Bio::Index::GenBank(-filename => $index);
    
} else { 
    if( $remote =~ /(ncbi)|(genbank)/i ) {

	$db = new Bio::DB::GenBank;
    } elsif( $remote =~ /embl/i ) {
	$db = new Bio::DB::EMBL;
    } else { 
	die("remote must be either 'NCBI' or 'EMBL'");
    }
    # would need to add code to set proxy info for those who need it
}

if(! $blastfile || ! -r $blastfile ) {
    die("Must specify a valid blastfile");
}

my $parser = new Bio::SearchIO(-format => $format,
			       -file => $blastfile);

my %tissues_seen = ();
my ($result,$hit,$hsp);
while( my $result = $parser->next_result )  {
  HIT: while( my $hit = $result->next_hit ) {
      if( defined $pvalue ) {
	  while( my $hsp = $hit->next_hsp ) {
	      if( $hsp->evalue > $pvalue ) {
		  print "skipping ", $hit->name, " because of low evalue \n";
		  # skip this Subject if it contains a pvalue of > $pvalue
		  next HIT;
	      }
	  }
      }
      my  ($id) = split(/\s+/, $hit->name);
      # get the last value
      my @ids = split(/\|/, $id);
      $id = pop @ids;
      my ($tissuetype) = get_tissue($id);
      if( defined $tissuetype ) {
	  push @{$tissues_seen{$tissuetype}}, $hit->name;
      } else { 
	  print STDERR "could not find tissue for $id\n" if( $VERBOSE);
      }
  }
  print "tissues seen for: ", $result->query_name, "\n";

  foreach my $tissue ( sort keys %tissues_seen ) {
      print "* $tissue\n-----------\n\t", 
      join("\n\t",@{$tissues_seen{$tissue}}), "\n\n";
  }
}

# cleanup -- avoid segfault here
$db = undef;

# subroutines

sub get_tissue {
    my ($id) = @_;
    my $tissue;
    if( $tissue = $accessions{$id} ) {
	return $tissue;
    }

    my $seq = $db->get_Seq_by_acc($id);
    return  unless(  $seq );

    foreach my $feature ( $seq->all_SeqFeatures ) {
	if( $feature->primary_tag eq 'source' ) {
	    foreach my $tag ( sort { $b cmp $a }
			      $feature->all_tags ) {
		if( $tag =~ /tissue/i  || 
		    ( ! $tissue && 
		      $tag =~ /clone_lib/i ) ){
		    ($tissue) = $feature->each_tag_value($tag);
		    $accessions{$seq->display_id} = $tissue;
		    return $tissue;
		}
	    }
	}
    }	    
    return;
}
