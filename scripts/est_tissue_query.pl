#!/usr/local/bin/perl -w

# this script assumes you have a directory which you have downloaded
# gbestXX.seq.gz from ncbi genbank release (leave them compressed)
#
# the script is unix centric for the time being, could use Compress::Gzip to
# uncompress the files in a platform indepedent way
#
# cmd line options are
# -i/--index=indexname
# -d/--dir=dir where gbest data files are located
# -b/--blast=filename blast filename which compared against an EST db
# -p=pvalue pvalue to limit search to
# -r/--remote=[GenBank|EMBL] use remote db for searching

use strict;
use DB_File;
use Bio::SeqIO;
use Bio::Tools::BPlite;
use Bio::DB::EMBL;
use Bio::DB::GenBank;

use Getopt::Long;

my $dir = '/home/data/libs/gbest';
my $index = 'dbest_tissue.idx';
my $VERBOSE = 0;
my $blastfile;
my $pvalue;
my $remote;
my %accessions;

&GetOptions( 'd|dir:s'   => \$dir,
	     'i|index:s' => \$index,
	     'v|verbose' => \$VERBOSE,
	     'b|blast:s' => \$blastfile,
	     'p:s'       => \$pvalue,
	     'r|remote:s'=> \$remote);


if( ! $remote ) {
    opendir(GBEST, $dir) or die("cannot open $dir");


    my $already = ( -r $index); 
    tie %accessions, "DB_File", $index,  O_RDWR|O_CREAT,0640, $DB_HASH;

    unless( $already ) {
	foreach my $file  ( readdir(GBEST) ) {
#	print "file is $file\n";
	    next unless ( $file =~ /gbest\d+\.seq\.gz$/ );
	    open(IN, "gunzip -c $dir/$file |") or 
		die("Cannot open file $dir/$file");
	    my $seqio = new Bio::SeqIO(-format => 'genbank', -fh => \*IN);

	    while( my $seq = $seqio->next_seq ) {
		my $tissue ='';
	      FEATURE: foreach my $feature ( $seq->all_SeqFeatures ) {
		  if( $feature->primary_tag eq 'source' ) {
		      foreach my $tag ( sort { $b cmp $a }
					$feature->all_tags ) {
			  if( $tag =~ /tissue/i  || 
			      ( ! $tissue && 
				$tag =~ /clone_lib/i ) ){
			      ($tissue) = $feature->each_tag_value($tag);
			      $accessions{$seq->display_id} = $tissue;
			      last FEATURE;
			  }
		      }
		  }
	      }
		if( ! defined $accessions{$seq->display_id} ) {
		    print STDERR "accession ", $seq->display_id, " skipped no tissue information\n" if $VERBOSE;
		}
	    }
	} 
    }
    print "there are ", scalar keys %accessions, " accessions in the db\n";
} else { 
    if( $remote =~ /(ncbi)|(genbank)/i ) {

	$remote = new Bio::DB::GenBank;
    } elsif( $remote =~ /embl/i ) {
	$remote = new Bio::DB::EMBL;
    } else { 
	die("remote must be either 'NCBI' or 'EMBL'");
    }
    # would need to add code to set proxy info for those who need it
}

if(! $blastfile || ! -r $blastfile ) {
    die("Must specify a valid blastfile");
}

my $parser = new Bio::Tools::BPlite(-file => $blastfile);
if( $parser->database !~ /est/i ) {
    print "BLAST db was '", $parser->database, "'\n which does not appear to be an EST database.  Proceeding anyways.\n"; 
}


my %tissues_seen;
SUBJECT: while( my $sbjct = $parser->nextSbjct ) {
    if( defined $pvalue ) {
	while( my $hsp = $sbjct->nextHSP ) {
	    if( $hsp->P > $pvalue ) {
		print "skipping ", $sbjct->name, "\n";		
		# skip this Subject if it contains a pvalue of > $pvalue
		next SUBJECT;
	    }
	}	
    }
    my  ($id) = split(/\s+/, $sbjct->name);
    # get the last value
    my @ids = split(/\|/, $id);
    $id = pop @ids;
    my ($tissuetype) = get_tissue($id);
    if( defined $tissuetype ) {
	push @{$tissues_seen{$tissuetype}}, $sbjct->name;
    } else { 
	print STDERR "could not find tissue for $id\n";
    }
}

print "tissues seen for: ", $parser->query, "\n";
foreach my $tissue ( sort keys %tissues_seen ) {
    print "* $tissue\n-----------\n\t", 
    join("\n\t",@{$tissues_seen{$tissue}}), "\n\n";
}

sub get_tissue {
    my ($id) = @_;
    my $tissue;
    if( $tissue = $accessions{$id} ) {
	return $tissue;
    }
    
    if( $remote ) {
	my $seq = $remote->get_Seq_by_acc($id);
	if( ! $seq ) {
	    print STDERR "unable to find seq for id $id\n";
	    return '';
	}
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
    }
    return undef;
}
