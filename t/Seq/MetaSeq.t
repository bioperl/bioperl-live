# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 132);
	
	use_ok('Bio::Seq::Meta');
	use_ok('Bio::Seq::Meta::Array');
	use_ok('Bio::SeqIO');
	use_ok('Bio::AlignIO');
	use_ok('Bio::Seq::Quality');
}

my $DEBUG = test_debug();

ok my $seq = Bio::Seq::Meta->new( -seq => "AT-CGATCGA");
is $seq->is_flush, 1;
is $seq->revcom->seq, 'TCGATCG-AT';
is $seq->meta, "";
ok $seq->force_flush(1);
is $seq->meta, "          ";
$seq->seq("AT-CGATCGATT");
is $seq->meta, "            ";
ok not $seq->force_flush(0);

ok $seq = Bio::Seq::Meta::Array->new( -seq => "AT-CGATCGA");
is $seq->is_flush, 1;
is $seq->revcom->seq, 'TCGATCG-AT';
is $seq->meta_text, "";
ok $seq->force_flush(1);
$seq->seq("AT-CGATCGATT");
is $seq->meta_text, "0 0 0 0 0 0 0 0 0 0 0 0";
ok not $seq->force_flush(0);

ok $seq = Bio::Seq::Quality->new( -seq => "AT-CGATCGA");
is $seq->meta_text, "";
ok $seq->force_flush(1);
is $seq->meta_text, "0 0 0 0 0 0 0 0 0 0";
$seq->seq("AT-CGATCGATT");
is $seq->meta_text, "0 0 0 0 0 0 0 0 0 0 0 0";
ok not $seq->force_flush(0);

ok $seq = Bio::Seq::Meta->new
    ( -seq => "",
      -meta => "",
      -alphabet => 'dna',
      -id => 'myid'
    );

# create a sequence object
ok $seq = Bio::Seq::Meta->new( -seq => "AT-CGATCGA",
                               -id => 'test',
                               -verbose => 2,
                               -force_flush => 1
                             );

is $seq->meta, "          ";
is $seq->meta_length, 10;

# Create some random meta values, but gap in the wrong place
my $metastring = "a-abb  bb ";
$seq->meta($metastring);
$seq->verbose(1);

# create some random meta values, but not for the last residue
$metastring = "aa-bb  bb";
ok $seq->meta($metastring), $metastring. " ";

# truncate the sequence by assignment
$seq->force_flush(1);
$seq->seq('AT-CGA');
$seq->alphabet('dna');
is $seq->meta, 'aa-bb ';
is $seq->start, 1;
is $seq->end, 5;
$seq->force_flush(0);

# truncate the sequence with trunc()
is $seq->strand(-1), -1;
ok $seq = $seq->trunc(1,5);
is $seq->start, 2;
is $seq->end, 5;
is $seq->seq, 'AT-CG';
is $seq->meta, 'aa-bb';
is $seq->strand, -1;

# revcom
ok $seq = $seq->revcom;
is $seq->seq, 'CG-AT';
is $seq->meta, 'bb-aa';
is $seq->strand, 1;

# submeta
is $seq->subseq(2,4), 'G-A';
is $seq->submeta(2,4), 'b-a';
is $seq->submeta(2,undef, 'c-c'), 'c-ca';
is $seq->submeta(2,4), 'c-c';
is $seq->meta, 'bc-ca';
is $seq->meta(''), '     ';
is $seq->submeta(2,undef, 'c-c'), 'c-c ';
is $seq->meta, ' c-c ';

# add named meta annotations

my $first = '11-22';
is $seq->named_meta('first', $first), $first;
is $seq->named_meta('first'), $first;

my $second = '[[-]]';
ok $seq->named_meta('second', $second);

# undefined range arguments
is $seq->named_submeta('second', 3, 4), '-]';
is $seq->named_submeta('second', 3), '-]]';
is $seq->named_submeta('second'), '[[-]]';

my @names =  $seq->meta_names;
is @names, 3;
is $names[0], 'DEFAULT';



#
# IO tests
#

sub diff {
    my ($infile, $outfile) = @_;
    my ($in, $out);
    open my $FH_IN, '<', $infile or die "Could not read file '$infile': $!\n";
    $in .= $_ while (<$FH_IN>);
    close $FH_IN;

    open my $FH_OUT, '<', $outfile or die "Could not read file '$outfile': $!\n";
    $out .= $_ while (<$FH_OUT>);
    close $FH_OUT;
    print "|$in||$out|\n" if $DEBUG;
    is $in, $out;
}


# SeqIO
my $str = Bio::SeqIO->new
    ( '-file'=> test_input_file('test.metafasta'),
      '-format' => 'metafasta');
ok  $seq = $str->next_seq;

my $outfile = test_output_file();
my $strout = Bio::SeqIO->new
    ('-file'=> ">". $outfile,
     '-format' => 'metafasta');
ok $strout->write_seq($seq);

diff (test_input_file('test.metafasta'),
      $outfile
     );

# AlignIO

$str = Bio::AlignIO->new
    ( '-file'=> test_input_file('testaln.metafasta'),
      '-format' => 'metafasta');
ok my $aln = $str->next_aln;

$outfile = test_output_file();
$strout = Bio::AlignIO->new
    ('-file'=> ">". $outfile,
     '-format' => 'metafasta');
ok $strout->write_aln($aln);

diff (test_input_file('testaln.metafasta'),
      $outfile
     );

#
##
### tests for Meta::Array
##
#

ok $seq = Bio::Seq::Meta::Array->new
    ( -seq => "",
      -meta => "",
      -alphabet => 'dna',
      -id => 'myid'
    );

# create a sequence object
ok $seq = Bio::Seq::Meta::Array->new( -seq => "AT-CGATCGA",
                                      -id => 'test',
                                      -force_flush => 1,
                                      -verbose => 2
                             );

is $seq->is_flush, 1;
#is $seq->meta_text, "          ";
is $seq->meta_text, '0 0 0 0 0 0 0 0 0 0';

# create some random meta values, but not for the last residue
$metastring = "a a - b b 0 b b 0";
is join (' ',  @{$seq->meta($metastring)}), $metastring. ' 0';
is $seq->meta_text, $metastring. ' 0';

# truncate the sequence by assignment
$seq->seq('AT-CGA');
$seq->alphabet('dna');
is $seq->meta_text, 'a a - b b 0';

# truncate the sequence with trunc()
is $seq->strand(-1), -1;
ok $seq = $seq->trunc(1,5);
is $seq->seq, 'AT-CG';
is $seq->meta_text, 'a a - b b';
is $seq->strand, -1;

#is $seq->length, 5;
#is $seq->meta_length, 6;
#ok $seq->force_flush(1);
#is $seq->meta_length, 5;

# revcom
ok $seq = $seq->revcom;
is $seq->seq, 'CG-AT';
is $seq->meta_text, 'b b - a a';
is $seq->strand, 1;

# submeta

is $seq->subseq(2,4), 'G-A';

is $seq->submeta_text(2,4), 'b - a';
is $seq->submeta_text(2,undef, 'c - c'), 'c - c';
is $seq->submeta_text(2,4), 'c - c';
is $seq->meta_text, 'b c - c a';

is $seq->meta_text(''), '0 0 0 0 0';
is $seq->submeta_text(2,undef, 'c - c'), 'c - c';
is $seq->meta_text, '0 c - c 0';

# add named meta annotations
$first = '1 10 - 222 23';
is $seq->named_meta_text('first', $first), $first;
is $seq->named_meta_text('first'), $first;
$second = '[ [ - ] ]';
ok $seq->named_meta_text('second', $second);

# undefined range arguments
is $seq->named_submeta_text('second', 3, 4), '- ]';
is $seq->named_submeta_text('second', 3), '- ] ]';
is $seq->named_submeta_text('second'), '[ [ - ] ]';

@names =  $seq->meta_names;
is @names, 3;
is $names[0], 'DEFAULT';




#
# testing the forcing of flushed meta values
#




ok $seq = Bio::Seq::Meta->new( -seq =>  "AT-CGATCGA",
                                  -id => 'test',
                                  -verbose => 2
                             );
is $seq->submeta(4, 6, '456'), '456';
is $seq->meta_length, 6;
is $seq->length, 10;

is $seq->meta, "   456";

ok $seq->force_flush(1);
is $seq->meta, "   456    ";
ok $seq->seq('aaatttc');
is $seq->meta, "   456 ";

ok $seq = Bio::Seq::Meta::Array->new( -seq =>  "AT-CGATCGA",
                                  -id => 'test',
                                  -verbose => 2
                             );
is join (' ', @{$seq->submeta(4, 6, '4 5 6')}), '4 5 6';
is $seq->meta_length, 6;
is $seq->length, 10;

is $seq->meta_text, "0 0 0 4 5 6";
ok $seq->force_flush(1);
is $seq->meta_text, "0 0 0 4 5 6 0 0 0 0";

ok $seq->seq('aaatttc');
is $seq->meta_text, "0 0 0 4 5 6 0";
is $seq->meta_length, 7;


ok  $seq = Bio::Seq::Quality->new( -seq =>  "AT-CGATCGA",
                                  -id => 'test',
                                  -verbose => 2
                             );
is join (' ', @{$seq->submeta(4, 6, '4 5 6')}), '4 5 6';
is $seq->meta_length, 6;
is $seq->length, 10;

is $seq->meta_text, "0 0 0 4 5 6";

ok $seq->force_flush(1);

is $seq->meta_text, "0 0 0 4 5 6 0 0 0 0";

ok $seq->seq('aaatttc');
is $seq->meta_text, "0 0 0 4 5 6 0";
is $seq->meta_length, 7;
is $seq->trace_length, 7;
#is $seq->quality_length, 7;

is $seq->is_flush, 1;
is $seq->trace_is_flush, 1;
is $seq->quality_is_flush, 1;

# quality: trace_lengths, trace_is_flush, quality_is_flush
