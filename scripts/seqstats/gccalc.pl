#!/usr/bin/perl -w
# $Id$

# Author Jason Stajich <jason@bioperl.org> 
# based on script code (see bottom) submitted by cckim@stanford.edu
#
# Submitted as part of bioperl script project 2001/08/06

use strict;

use Bio::SeqIO;
use Bio::Tools::SeqStats;
use Getopt::Long;
my $format = 'fasta';
my $file;
my $help =0;
GetOptions( 
	    'f|format:s' => \$format,
	    'i|in:s'     => \$file,
	    'h|help'    => \$help,
	    );

if( $help ) { 
    die("usage: gccalc.pl -i filename -f format\n\t or gccalc.pl -f format < filename");
}
my $seqin;

if( defined $file ) { 
    $seqin = new Bio::SeqIO(-format => $format,
			    -file   => $file);
} else { 
    $seqin = new Bio::SeqIO(-format => $format,
			    -fh     => \*STDIN);
}

while( my $seq = $seqin->next_seq ) {
    next if( $seq->length == 0 );
    my $seq_stats  =  Bio::Tools::SeqStats->new('-seq'=>$seq);
    my $hash_ref = $seq_stats->count_monomers();  # for DNA sequence
    print "Seq: ", $seq->display_id, " ", $seq->desc, 
    " Len:", $seq->length, "\n";
    printf "GC content is %.4f\n", ($hash_ref->{'G'} + $hash_ref->{'C'}) / 
	$seq->length();
    
    foreach my $base (sort keys %{$hash_ref}) {
	print "Number of bases of type ", $base, "= ", $hash_ref->{$base},"\n";
    }
    print "--\n";
}

# alternatively one could use code submitted by
# cckim@stanford.edu

sub calcgc {
    my $seq = $_[0];
    my @seqarray = split('',$seq);
    my $count = 0;
    foreach my $base (@seqarray) {
	$count++ if $base =~ /[G|C]/i;
    }
    
    my $len = $#seqarray+1;
    return $count / $len;
}
