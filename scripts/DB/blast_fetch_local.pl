#!/usr/local/bin/perl -w

# This is a working script to extract sequence ids from a blast 
# result file and extract those sequences from a local db.

# This may not work if the blastdb to index is est depending on 
# wha your machine supports

# $Id$

use strict;
use Bio::Tools::BPlite;
use Bio::SeqIO;
use Bio::Index::Fasta;

use Getopt::Long;
use vars qw($USAGE);
$USAGE = 'usage: blast_fetch_local.pl -b report.bls -o outseqs.fa -d /path/to/db -i /path/to/index -p pvalue';

my ($db, $index, $output, $blastfile,$pvalue);

$pvalue = 1e-5;

&GetOptions
    (
     'd|db:s'  => \$db,
     'i|index:s' => \$index,
     'o|output:s' => \$output,
     'b|blast:s' => \$blastfile,
     'p|pvalue:s' => \$pvalue);

die $USAGE if( ! $db || ! -r $db || ! $index || ! $output 
	       || ! $blastfile || ! -r $blastfile );
if( ! -r $index ) {
    # skip if it is already created
    my $indexdb = new Bio::Index::Fasta('-filename' => $index,
					'-write_flag' => 1);
    $indexdb->id_parser( \&parse_ncbi_id );
    $indexdb->make_index($db);
    undef $indexdb;
}

my $seqout = new Bio::SeqIO(-file => ">$output");
my $indexdb = new Bio::Index::Fasta('-filename' => $index);
my $bplite = new Bio::Tools::BPlite( -file => $blastfile );

SBJCT: while( my $sbjct = $bplite->nextSbjct ) {
    HSP: while( my $hsp = $sbjct->nextHSP ) { 
	if( $hsp->P > $pvalue) {	
	    # skip Sbjcts that don't meet the minimum P value
	    print 'skipping Sbjct with a HSP with pvalue=', $hsp->P, "\n";
	    next SBJCT;    
	} 
# could also accept any Subject with an HSP that has appropriate e-value
#	else {
#	    last HSP;
#	}
    }
    my  ($id) = split(/\s+/, $sbjct->name);
    # get the last value
    my @ids = split(/\|/, $id);
    $id = pop @ids;
    my $seq = $indexdb->fetch($id);
    if( ! $seq ) {
	print "could not find id '$id'\n";
    } else {
	$seqout->write_seq($seq);
    }
}

sub parse_ncbi_id {
    my @retvals;
    if( $_[0] =~ /^>(\S+)/ ) {
	my $val = $1;
	my (@elements) = split(/\|/, $val);
	while( @elements ) {
	    my $id = shift @elements;
	    if( $id eq 'gb' || $id eq 'gi' ) {
		$id = $id . '|' . shift @elements;
	    }
	    push @retvals, $id;
	}
    }
    return @retvals;
}
