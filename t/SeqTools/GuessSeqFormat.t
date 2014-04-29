# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
   
    test_begin(-tests => 105);
   
    use_ok 'Bio::Tools::GuessSeqFormat';
    use_ok 'Bio::SeqIO';
    use_ok 'Bio::AlignIO';
}


ok my $guesser = Bio::Tools::GuessSeqFormat->new;
isa_ok $guesser, 'Bio::Tools::GuessSeqFormat';


#
# Test guesser interfaces
#

# 1/ File guess
ok $guesser = Bio::Tools::GuessSeqFormat->new(
    -file => test_input_file('test.fasta'),
), 'File input';
is $guesser->guess, 'fasta';

# 2/ String guess
my $string = ">test1 no comment
agtgctagctagctagctagct
>test2 no comment
gtagttatgc
";
ok $guesser = Bio::Tools::GuessSeqFormat->new(
    -text => $string,
), 'String input';
is $guesser->guess, 'fasta';

# 3/ Filehandle guess
SKIP: {
    test_skip(-tests => 2, -requires_modules => [qw(IO::String)]);
    require IO::String;
    my $fh = IO::String->new($string);
    ok $guesser = Bio::Tools::GuessSeqFormat->new(
        -fh => $fh,
    ), 'Filehandle input';
    is $guesser->guess, 'fasta';
}


#
# Test behavior with unknown format
#

is $guesser = Bio::Tools::GuessSeqFormat->new(
    -file => test_input_file('test.waba'), # Bio::SearchIO::waba
)->guess, undef, 'Unknown file format';

throws_ok {
    Bio::SeqIO->new( -file=>test_input_file('test.waba') );
} qr/Could not guess format/;


#
# Test SeqIO formats
#

my @fmts = qw{ace embl fasta fastq game gcg genbank pir raw swiss tab};

for my $fmt (@fmts) {
    SKIP: {
        test_skip(
            -tests => 4,
            -requires_modules => [qw(XML::Writer XML::Parser::PerlSAX)]
        ) if $fmt eq 'game';
        test_skip(
            -tests => 4,
            -requires_module  => 'Data::Stag'
        ) if $fmt eq 'swiss';

        my $guess = Bio::Tools::GuessSeqFormat->new(
            -file => test_input_file("test.$fmt"),
        )->guess;
        is $guess, $fmt, "$fmt format";

        ok my $input = Bio::SeqIO->new( -file=>test_input_file("test.$fmt") );
        ok my $seq = $input->next_seq();
        isa_ok $seq, 'Bio::PrimarySeqI';
    }
}


#
# Test AlignIO formats
#

@fmts = qw{clustalw fasta fastq mase mega msf nexus pfam phylip prodom selex stockholm};
my %skip_module = map {$_=>1} qw{ fastq };

for my $fmt (@fmts) {
    my $guess = Bio::Tools::GuessSeqFormat->new(
        -file => test_input_file("testaln.$fmt")
    )->guess;
    is $guess, $fmt, "$fmt format";
    next if $skip_module{$fmt};

    ok my $input = Bio::AlignIO->new( -file=>test_input_file("testaln.$fmt") );
    ok my $seq = $input->next_aln();
    isa_ok $seq, 'Bio::Align::AlignI';
}


#
# Other formats
#

my %fmts = (
   blast    => test_input_file('blastp2215.blast'),
   gcgblast => test_input_file('test.gcgblast'),
   vcf      => test_input_file('example.vcf'),
);

while (my ($fmt, $file) = each %fmts) {
    my $guess = Bio::Tools::GuessSeqFormat->new(
        -file => $file,
    )->guess;
    is $guess, $fmt, "$fmt format";
}
