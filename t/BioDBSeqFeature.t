##-*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use Bio::Root::IO;
use FindBin '$Bin';
use constant TEST_COUNT => 52;
use constant GFF_FILE    => Bio::Root::IO->catfile('t','data',
					   'seqfeaturedb','test.gff3');

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan test => TEST_COUNT;
    $ENV{ORACLE_HOME} ||= '/home/oracle/Home';
}

use lib "$Bin/..","$Bin/../blib/lib";
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;

sub bail ($;$) {
  my $count = shift;
  my $explanation = shift;
  for (1..$count) {
    skip($explanation,1);
  }
  exit 0;
}

sub fail ($) {
  my $count = shift;
  for (1..$count) {
    ok(0);
  }
  exit 0;
}

my (@f,$f,@s,$s,$seq1,$seq2);

my @args = @ARGV;
@args = (-adaptor => 'memory') unless @args;

my $db = eval { Bio::DB::SeqFeature::Store->new(@args) };
warn $@ if $@;
ok($db);
fail(TEST_COUNT - 1) unless $db;

my $loader = eval { Bio::DB::SeqFeature::Store::GFF3Loader->new(-store=>$db) };
warn $@ if $@;
ok($loader);
fail(TEST_COUNT - 2) unless $loader;

# exercise the loader
ok($loader->load(GFF_FILE));

# there should be one gene named 'abc-1'
@f = $db->get_features_by_name('abc-1');
ok(@f==1);

$f = $f[0];
# there should be three subfeatures of type "exon" and three of type "CDS"
ok($f->get_SeqFeatures('exon')==3);
ok($f->get_SeqFeatures('CDS')==3);

# the sequence of feature abc-1 should match the sequence of the first exon at the beginning
$seq1 = $f->seq->seq;
$seq2 = (sort {$a->start<=>$b->start} $f->get_SeqFeatures('exon'))[0]->seq->seq;
ok(substr($seq1,0,length $seq2) eq $seq2);

# sequence lengths should match
ok(length $seq1 == $f->length);

# if we pull out abc-1 again we should get the same object
($s) = $db->get_features_by_name('abc-1');
ok($f eq $s);

# we should get two objects when we ask for abc-1 using get_features_by_alias
# this also depends on selective subfeature indexing
@f = $db->get_features_by_alias('abc-1');
ok(@f==2);

# the two features should be different
ok($f[0] ne $f[1]);

# test that targets are working
($f) = $db->get_features_by_name('match1');
ok(defined $f);
$s = $f->target;
ok(defined $s);
ok($s->seq_id  eq 'CEESC13F');
$seq1 = $s->seq->seq;
ok(substr($seq1,0,10) eq 'ttgcgttcgg');

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
ok(@transcripts == 2);
ok($transcripts[0]->method eq 'mRNA');
ok($transcripts[0]->source eq 'confirmed');

# test that exon #2 is shared between the two transcripts
my @exons1      = $transcripts[0]->get_SeqFeatures('CDS');
ok(@exons1 == 3);
my @exons2      = $transcripts[1]->get_SeqFeatures('CDS');
my ($shared1)   = grep {$_->display_name||'' eq 'shared_exon'} @exons1;
my ($shared2)   = grep {$_->display_name||'' eq 'shared_exon'} @exons2;
ok($shared1 && $shared2);
ok($shared1 eq $shared2);
ok($shared1->primary_id eq $shared2->primary_id);

# test attributes
ok($shared1->phase == 0);
ok($shared1->strand eq +1);
ok(($f->attributes('expressed'))[0] eq 'yes');

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
ok($gene3a->Three_prime_UTR == 2);
ok($gene3b->Three_prime_UTR == 1);
my ($utr) = $gene3b->Three_prime_UTR;
ok($utr->segments == 2);
my $location = $utr->location;
ok($location->isa('Bio::Location::Split'));
ok($location->sub_Location == 2);

# ok, test that queries are working properly.
# find all features with the attribute "expressed"
@f = $db->get_features_by_attribute({expressed=>'yes'});
ok(@f == 2);

# find all top-level features on Contig3 -- there should be two
@f = $db->get_features_by_location(-seq_id=>'Contig3');
ok(@f == 2);

# find all top-level features on Contig3 of type 'assembly_component'
@f = $db->features(-seq_id=>'Contig3',-type=>'assembly_component');
ok(@f==1);

# test iteration
@f = $db->features;
my $feature_count = @f;
ok($feature_count > 0);

my $i = $db->get_seq_stream;
ok($i);

my $count;
while ($i->next_seq) { $count++ }
ok($feature_count == $count);

# regression test on bug in which get_SeqFeatures('type') did not filter inline segments
@f = $db->get_features_by_name('agt830.3');
ok(@f && !$f[0]->get_SeqFeatures('exon'));
ok(@f && $f[0]->get_SeqFeatures('EST_match'));

# regression test on bug in which the load_id disappeared
ok(@f && $f[0]->load_id eq 'Match2');

# regress on proper handling of multiple ID features
my ($alignment) = $db->get_features_by_name('agt830.5');
ok($alignment);
ok($alignment->target->start == 1 && $alignment->target->end == 654);
ok($alignment->get_SeqFeatures == 2);
my $gff3 = $alignment->gff3_string(1);
my @lines = split "\n",$gff3;
ok (@lines == 2);
ok ("@lines" !~ /Parent=/s);
ok ("@lines" =~ /ID=/s);

1;

__END__

