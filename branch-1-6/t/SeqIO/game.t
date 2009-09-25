# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 24,
		   -requires_modules => [qw(XML::Parser::PerlSAX XML::Writer)]);
    use_ok('Bio::SeqIO::game');
}

my $verbose = test_debug() || -1;
my $str = Bio::SeqIO->new('-file'=> test_input_file('test.game'), 
			  '-format' => 'game',
			  '-verbose' => $verbose);
isa_ok ($str, 'Bio::SeqIO');
my $seq = $str->next_seq();
isa_ok($seq, 'Bio::Seq::RichSeq');

# exercise game parsing
$str = Bio::SeqIO->new(
    -format =>'game',
    -file => test_input_file('test.game')
		      );
$seq = $str->next_seq;
ok(defined $seq);
ok(defined $seq->seq);
is($seq->alphabet, 'dna');
is($seq->display_id, 'L16622');
is($seq->length, 28735);
is($seq->species->binomial, 'Caenorhabditis elegans');
my @feats = $seq->get_SeqFeatures;
is(scalar(@feats), 7);
my $source = grep { $_->primary_tag eq 'source' } @feats;
ok($source);
my @genes = grep { $_->primary_tag eq 'gene' } @feats;
is(scalar(@genes), 3);
ok($genes[0]->has_tag('gene'));
my $gname;
if ( $genes[0]->has_tag('gene') ) {
    ($gname) = $genes[0]->get_tag_values('gene');
}
is($gname, 'C02D5.3');
is($genes[0]->strand, 1);
my $cds   = grep { $_->primary_tag eq 'CDS' } @feats;
is($cds, 3);

# make sure we can read what we write
# test XML-writing
my $testfile = test_output_file();
# map argument is require to write a <map_position> element
my $out = Bio::SeqIO->new(-format => 'game', -file => ">$testfile", -map => 1);
$out->write_seq($seq);
$out->close();

$str = Bio::SeqIO->new(-format =>'game', -file => $testfile);
$seq = $str->next_seq;
ok(defined $seq);
ok(defined $seq->seq);
is($seq->alphabet, 'dna');
is($seq->display_id, 'L16622');
is($seq->length, 28735);
is($seq->species->binomial, 'Caenorhabditis elegans');

my $genes = grep { $_->primary_tag eq 'gene' } @feats;
$cds   = grep { $_->primary_tag eq 'CDS' } @feats;
is($genes, 3);
is($cds, 3);
