#!/usr/local/bin/perl -w
use strict;

use Bio::SeqIO;
use Bio::Root::IO;
use Algorithm::Diff qw(diff LCS);
use IO::ScalarArray;
use IO::String;

my %files = ( 
#'test.embl'      => 'embl',
#	      'test.ace'       => 'ace',
	      'test.fasta'     => 'fasta',
#	      'test.game'      => 'game',
	      'test.gcg'       => 'gcg',
#	      'test.genbank'   => 'genbank',
	      'test.raw'       => 'raw',
#	      'test_badlf.gcg' => 'gcg'
	      );

while( my ($file, $type) = each %files ) {
    my $filename = Bio::Root::IO->catfile('t','data',$file);
    print "processing file $filename\n";
    open(FILE, "< $filename") or die("cannot open $filename");
    my @datain = <FILE>;
    my $in = new IO::String(join('', @datain));
    my $seqin = new Bio::SeqIO( -fh => $in,
				-format => $type);
    my $out = new IO::String;
    my $seqout = new Bio::SeqIO( -fh => $out,
				 -format => $type);
    my $seq;
    while( defined($seq = $seqin->next_seq) ) {	
	$seqout->write_seq($seq);
    }
    $seqout->close();
    $seqin->close();
    my $strref = $out->string_ref;
    my @dataout = map { $_."\n"} split(/\n/, $$strref );
    my @diffs = &diff( \@datain, \@dataout);
    foreach my $d ( @diffs ) {
	foreach my $diff ( @$d ) {
	    chomp($diff->[2]);
	    print $diff->[0], $diff->[1], "\n>", $diff->[2], "\n";
	}
    }
    if( @diffs ) {
	print "in is \n", join('', @datain), "\n";
	print "out is \n", join('',@dataout), "\n";	
    }
}
