# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

# test for Bio::Tools::GuessSeqFormat
# written by Heikki Lehvaslaiho

use strict;
my $NUMTESTS;

BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't','..';
    }
    use Test;
    $NUMTESTS = 43;
    plan tests => $NUMTESTS;

}

use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Tools::GuessSeqFormat;
use Data::Dumper;

ok 1;

my $format;
my $verbose =1;
#
# Seqio formats
#

#not tested:  waba 
my @seqformats = qw{ ace embl fasta game gcg 
                     genbank mase  pfam pir raw swiss tab };

my %no_seqio_module = map {$_=>1} qw {gcgblast gcgfasta mase pfam};

foreach $format (@seqformats) {
    my $guessed_format = new Bio::Tools::GuessSeqFormat
        (-file => Bio::Root::IO->catfile("t","data","test.$format"),
         #-verbose=> $verbose;
        )->guess;
    $format =~ s/\..*$//;
    ok $guessed_format, $format;
    next if $no_seqio_module{$format};

    eval {
        my $input = Bio::SeqIO->new
            (-file=>Bio::Root::IO->catfile("t","data","test.$format"));
        ok my $seq = $input->next_seq();
    };
    ok 0, 1, $@ if $@;
}


#
# AlignIO formats
#

@seqformats = qw{ aln:clustalw fasta mase msf nexus pfam phylip
                  prodom stockholm}; # not selex (same as pfam, mainly)

my %no_alignio_module = map {$_=>1} qw {};

foreach my $ext (@seqformats) {
    my $format;
    ($ext, $format) = split /:/, $ext;
    my $guesser = new Bio::Tools::GuessSeqFormat
        (-file => Bio::Root::IO->catfile("t","data","testaln.$ext"));
    $format ||= $ext;
    ok $guesser->guess(), $format;

    next if $no_alignio_module{$format};

    eval {
        my $input = Bio::AlignIO->new
            (-file=>Bio::Root::IO->catfile("t","data","testaln.$ext"));
        ok my $seq = $input->next_aln();
    };
    ok 0, 1, $@ if $@;

}


#
# File handle tests
#

use IO::String;

my $string = ">test1 no comment
agtgctagctagctagctagct
>test2 no comment
gtagttatgc
";

my $stringfh = new IO::String($string);

my $seqio = new Bio::SeqIO(-fh => $stringfh);
while( my $seq = $seqio->next_seq ) {
    ok $seq->id =~ /test/;
}

