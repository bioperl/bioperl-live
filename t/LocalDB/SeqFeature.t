# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;
use constant TEST_COUNT => 59;

BEGIN {
    use lib '/home/lstein/projects/bioperl-live';
    use lib '.';
    use Bio::Root::Test;
	
    test_begin(-tests => TEST_COUNT);
    
    $ENV{ORACLE_HOME} ||= '/home/oracle/Home';
	
    use_ok('Bio::DB::SeqFeature::Store');
    use_ok('Bio::DB::SeqFeature::Store::GFF3Loader');
}

my $gff_file = test_input_file('seqfeaturedb','test.gff3');

my (@f,$f,@s,$s,$seq1,$seq2);

my @args = @ARGV;
@args = (-adaptor => 'memory') unless @args;

SKIP: {
my $db = eval { Bio::DB::SeqFeature::Store->new(@args) };
skip "DB load failed? Skipping all! $@", (TEST_COUNT - 2) if $@;
ok($db);

my $loader = eval { Bio::DB::SeqFeature::Store::GFF3Loader->new(-store=>$db) };
skip "GFF3 loader failed? Skipping all! $@", (TEST_COUNT - 3) if $@;
ok($loader);

# exercise the loader
ok($loader->load($gff_file));

# there should be one gene named 'abc-1'
@f = $db->get_features_by_name('abc-1');
is(@f,1);

$f = $f[0];
# there should be three subfeatures of type "exon" and three of type "CDS"
is($f->get_SeqFeatures('exon'),3);
is($f->get_SeqFeatures('CDS'),3);

# the sequence of feature abc-1 should match the sequence of the first exon at the beginning
$seq1 = $f->seq->seq;
$seq2 = (sort {$a->start<=>$b->start} $f->get_SeqFeatures('exon'))[0]->seq->seq;
is(substr($seq1,0,length $seq2),$seq2);

# sequence lengths should match
is(length $seq1, $f->length);

# if we pull out abc-1 again we should get the same object
($s) = $db->get_features_by_name('abc-1');
is($f, $s);

# we should get two objects when we ask for abc-1 using get_features_by_alias
# this also depends on selective subfeature indexing
@f = $db->get_features_by_alias('abc-1');
is(@f,2);

# the two features should be different
isnt($f[0], $f[1]);

# test that targets are working
($f) = $db->get_features_by_name('match1');
ok(defined $f);
$s = $f->target;
ok(defined $s);
ok($s->seq_id  eq 'CEESC13F');
$seq1 = $s->seq->seq;
is(substr($seq1,0,10), 'ttgcgttcgg');

# can we fetch subfeatures?
# gene3.a has the Index=1 attribute, so we should fetch it
($f) = $db->get_features_by_name('gene3.a');
ok($f);

# gene 3.b doesn't have an index, so we shouldn't get it
($f) = $db->get_features_by_name('gene3.b');
ok(!$f);

# test three-tiered genes
($f) = $db->get_features_by_name('gene3');
ok($f);
my @transcripts = $f->get_SeqFeatures;
is(@transcripts, 2);
is($transcripts[0]->method,'mRNA');
is($transcripts[0]->source,'confirmed');

# test that exon #2 is shared between the two transcripts
my @exons1      = $transcripts[0]->get_SeqFeatures('CDS');
is(@exons1, 3);
my @exons2      = $transcripts[1]->get_SeqFeatures('CDS');
my ($shared1)   = grep {$_->display_name||'' eq 'shared_exon'} @exons1;
my ($shared2)   = grep {$_->display_name||'' eq 'shared_exon'} @exons2;
ok($shared1 && $shared2);
is($shared1, $shared2);
is($shared1->primary_id, $shared2->primary_id);

# test attributes
is($shared1->phase, 0);
is($shared1->strand, +1);
is(($f->attributes('expressed'))[0], 'yes');

# test autoloading
my ($gene3a) = grep { $_->display_name eq 'gene3.a'} @transcripts;
my ($gene3b) = grep { $_->display_name eq 'gene3.b'} @transcripts;
ok($gene3a);
ok($gene3b);
ok($gene3a->Is_expressed);
ok(!$gene3b->Is_expressed);

# the representation of the 3'-UTR in the two transcripts a and b is
# different (not recommended but supported by the GFF3 spec). In the
# first case, there are two 3'UTRs existing as independent
# features. In the second, there is one UTR with a split location.
is($gene3a->Three_prime_UTR, 2);
is($gene3b->Three_prime_UTR, 1);
my ($utr) = $gene3b->Three_prime_UTR;
is($utr->segments, 2);
my $location = $utr->location;
isa_ok($location, 'Bio::Location::Split');
is($location->sub_Location,2);

# ok, test that queries are working properly.
# find all features with the attribute "expressed"
@f = $db->get_features_by_attribute({expressed=>'yes'});
is(@f, 2);

# find all top-level features on Contig3 -- there should be two
@f = $db->get_features_by_location(-seq_id=>'Contig3');
is(@f, 2);

# find all top-level features on Contig3 of type 'assembly_component'
@f = $db->features(-seq_id=>'Contig3',-type=>'assembly_component');
is(@f, 1);

# test iteration
@f = $db->features;
my $feature_count = @f;
cmp_ok($feature_count, '>', 0);

my $i = $db->get_seq_stream;
ok($i);

my $count;
while ($i->next_seq) { $count++ }
is($feature_count,$count);

# regression test on bug in which get_SeqFeatures('type') did not filter inline segments
@f = $db->get_features_by_name('agt830.3');
ok(@f && !$f[0]->get_SeqFeatures('exon'));
ok(@f && $f[0]->get_SeqFeatures('EST_match'));

# regression test on bug in which the load_id disappeared
is(@f && $f[0]->load_id, 'Match2');

# regress on proper handling of multiple ID features
my ($alignment) = $db->get_features_by_name('agt830.5');
ok($alignment);
is($alignment->target->start,1);
is($alignment->target->end, 654);
is($alignment->get_SeqFeatures, 2);
my $gff3 = $alignment->gff3_string(1);
my @lines = split "\n",$gff3;
is (@lines, 2);
ok("@lines" !~ /Parent=/s);
ok("@lines" =~ /ID=/s);

SKIP: {
	test_skip(-tests => 2, -excludes_os => 'mswin');
	
	if (my $child = open(F,"-|")) { # parent reads from child
		cmp_ok(scalar <F>,'>',0);
		close F;
		# The challenge is to make sure that the handle
		# still works in the parent!
		my @f = $db->features();
		cmp_ok(scalar @f,'>',0);
	}
	else { # in child
		$db->clone;
		my @f = $db->features();
		my $feature_count = @f;
		print $feature_count;
		exit 0;
	}
}


# test the -ignore_seqregion flag

# the original should have a single feature named 'Contig1'
my @f   = $db->get_features_by_name('Contig1');
is(scalar @f,1);

$db     = eval { Bio::DB::SeqFeature::Store->new(@args) };
$loader = eval { Bio::DB::SeqFeature::Store::GFF3Loader->new(-store=>$db,
							     -ignore_seqregion=>1)
               };
$loader->load($gff_file);
@f      = $db->get_features_by_name('Contig1');
is(scalar @f,0);

}
