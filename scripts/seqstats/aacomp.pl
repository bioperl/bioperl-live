#!/usr/bin/perl -w
use strict;
use Carp;
# $Id$

# Author Jason Stajich <jason@bioperl.org> 
# based on aacomp.c from EMBOSS
#

use Bio::SeqIO;
use Getopt::Long;
use Bio::Tools::CodonTable;
use Bio::Tools::IUPAC;
my $table = new Bio::Tools::CodonTable(-id => 1);
my @BASES = $table->valid_aa(0);
my %all = $table->valid_aa(2);
my ($file,$format,$help) = ( undef, 'fasta');
GetOptions(
	   'i|in:s'  => \$file,
	   'f|format:s' => \$format,
	   'h|help'  => \$help,
	   );

die("usage: aacomp.pl -i filename [-f format]\n\tdefault format is fasta\n")
    if( $help );

my $seqio = new Bio::SeqIO(-format => $format,
			   -file   => $file);
my %composition;
my $total;
foreach my $base ( @BASES ) {
    $composition{$base} = 0;
}
while ( my $seq = $seqio->next_seq ) {
    if( $seq->alphabet ne 'protein' ) {
	confess("Must only provide amino acid sequences to aacomp...skipping this seq");
	next;
    }
    foreach my $base ( split(//,$seq->seq()) ) {
	$composition{$base}++;
	$total++;
    }
}

printf("%d aa\n",$total); 
printf("%5s %4s\n", 'aa', '#' );
my $ct = 0;
foreach my $base ( @BASES ) {
    printf(" %s %s %3d\n", $base, $all{$base}, $composition{$base} );
    $ct += $composition{$base};
}
printf( "%6s %s\n", '','-'x5);
printf( "%6s %3d\n", '',$ct);
