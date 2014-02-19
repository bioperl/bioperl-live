# -*-Perl-*- Test Harness script for Bioperl
# $Id$


use strict;

BEGIN {
   use lib '.';
   use Bio::Root::Test;
   
   test_begin(-tests => 96);
   
   use_ok 'Bio::Tools::GuessSeqFormat';
   use_ok 'Bio::SeqIO';
   use_ok 'Bio::AlignIO';
}


my $fmt;
my $seq;
my $verbose = test_debug();

ok my $guesser = Bio::Tools::GuessSeqFormat->new;
isa_ok $guesser, 'Bio::Tools::GuessSeqFormat';


#
# Test guesser interfaces
#

# File guess
ok $guesser = Bio::Tools::GuessSeqFormat->new(
    -file => test_input_file('test.fasta'),
), 'File input';
is $guesser->guess, 'fasta';

# String guess
my $string = ">test1 no comment
agtgctagctagctagctagct
>test2 no comment
gtagttatgc
";
ok $guesser = Bio::Tools::GuessSeqFormat->new(
    -text => $string,
), 'String input';
is $guesser->guess, 'fasta';

# Filehandle guess
SKIP: {
    test_skip(-tests => 2, -requires_modules => [qw(IO::String)]);
    require IO::String;
    my $fh = IO::String->new($string);
    ok $guesser = Bio::Tools::GuessSeqFormat->new(
        -text => $string,
    ), 'Filehandle input';
    is $guesser->guess, 'fasta';
}


#
# Test SeqIO formats
#

# waba is not guessed
is $guesser = Bio::Tools::GuessSeqFormat->new(
    -file => test_input_file('test.waba'),
)->guess, undef;

throws_ok {
    Bio::SeqIO->new( -file=>test_input_file('test.waba') );
} qr/Could not guess format/;

# other seq formats
my @fmts = qw{ace embl fasta fastq gcg genbank mase pfam pir raw swiss tab game};
my %skip_module = map {$_=>1} qw {gcgblast gcgfasta mase pfam};

for $fmt (@fmts) {
    SKIP: {
        test_skip(
            -tests => 4,
            -requires_modules => [qw(XML::Writer XML::Parser::PerlSAX)]
        ) if $fmt eq 'game';

        my $guess = Bio::Tools::GuessSeqFormat->new(
            -file => test_input_file("test.$fmt"),
            #-verbose=> $verbose;
        )->guess;
        is $guess, $fmt, "$fmt format";
        next if $skip_module{$fmt};

        ok my $input = Bio::SeqIO->new( -file=>test_input_file("test.$fmt") );
        ok $seq = $input->next_seq();
        isa_ok $seq, 'Bio::PrimarySeqI';
    }
}


#
# Test AlignIO formats
#

@fmts = qw{ aln:clustalw fasta fastq mase msf nexus pfam phylip prodom stockholm};
# not selex (same as pfam, mainly)
%skip_module = map {$_=>1} qw { fastq };

for my $ext (@fmts) {
    ($ext, my $fmt) = split /:/, $ext;
    $fmt ||= $ext;

    my $guess = Bio::Tools::GuessSeqFormat->new(
        -file => test_input_file("testaln.$ext")
    )->guess;
    is $guess, $fmt;
    next if $skip_module{$fmt};

    ok my $input = Bio::AlignIO->new( -file=>test_input_file("testaln.$ext") );
    ok $seq = $input->next_aln();
    isa_ok $seq, 'Bio::Align::AlignI';
}

