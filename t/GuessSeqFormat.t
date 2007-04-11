# -*-Perl-*- mode (to keep my emacs happy)
# $Id$
# test for Bio::Tools::GuessSeqFormat
# written by Heikki Lehvaslaiho

use strict;
my $NUMTESTS;
my $error;

BEGIN {
   eval { require Test::More; };
   if( $@ ) {
      use lib 't/lib';
   }
   use Test::More;
   $error = 0;
   # SeqIO::game needs XML::Writer and XML::Parser::PerlSAX
   eval {require XML::Writer; require XML::Parser::PerlSAX;};
   if ($@) {
      print STDERR "XML::Writer or XML::Parser::PerlSAX not found, skipping game test\n";
      $error = 1;
   }
   $NUMTESTS = ($error == 1) ? 48 : 50;
   plan tests => $NUMTESTS;
   use_ok('Bio::SeqIO');
   use_ok('Bio::AlignIO');
   use_ok('Bio::Tools::GuessSeqFormat');
   use_ok('Data::Dumper');
}

my @seqformats = qw{ ace embl fasta gcg genbank mase
                        pfam pir raw swiss tab };

push @seqformats,"game" if ($error == 0);

my $format;
my $verbose =1;
#
# Seqio formats
#

#not tested:  waba

my %no_seqio_module = map {$_=>1} qw {gcgblast gcgfasta mase pfam};

my $guessed_format = new Bio::Tools::GuessSeqFormat
        (-file => Bio::Root::IO->catfile("t","data","test.waba"))->guess;
is $guessed_format, undef ;

my $seq;

eval {
   my $input = Bio::SeqIO->new
       (-file=>Bio::Root::IO->catfile("t","data","test.waba"));
   $seq = $input->next_seq();
};

ok !$seq;

$@ ? ok 1 : ok 0;

foreach $format (@seqformats) {
    my $guessed_format = new Bio::Tools::GuessSeqFormat
        (-file => Bio::Root::IO->catfile("t","data","test.$format"),
         #-verbose=> $verbose;
        )->guess;
    $format =~ s/\..*$//;
    is $guessed_format, $format, "Guessed:$format";
    next if $no_seqio_module{$format};

    eval {
        my $input = Bio::SeqIO->new
            (-file=>Bio::Root::IO->catfile("t","data","test.$format"));
        $seq = $input->next_seq();
    };
    
    my $implemented = $format eq 'ace' ? 'Bio::PrimarySeqI' : 'Bio::SeqI';
    
    isa_ok $seq, $implemented;
    
    is 0, 1, $@ if $@;
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
        $seq = $input->next_aln();
    };
    
    isa_ok $seq, 'Bio::Align::AlignI';
    #ok 0, 1, $@ if $@;
}


#
# File handle tests
#
if( eval 'require IO::String; 1' ) {

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
    
#
# text guessing
#

    ok new Bio::Tools::GuessSeqFormat( -text => $string )->guess, 'fasta';
} else {
    for (1..3) {
	skip("skipping guessing format from string, IO::String not installed",1);
    }
}
