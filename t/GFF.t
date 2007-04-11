#-*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;

BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    plan tests => 36;
    use_ok('Bio::Seq');
    use_ok('Bio::Tools::GFF');
    use_ok('Bio::SeqFeatureI');
    use_ok('Bio::SeqFeature::Generic');
}

END {
    unlink("out1.gff", "out2.gff");
}

my $feat = new Bio::SeqFeature::Generic( -start => 10, -end => 100,
				-strand => -1, -primary => 'repeat',
				-source => 'repeatmasker',
				-score  => 1000,
				-tag    => {
				    new => 1,
				    author => 'someone',
				    sillytag => 'this is silly!;breakfast' } );
ok($feat);
my $gff1out = Bio::Tools::GFF->new(-gff_version => 1, -file => ">out1.gff");
ok($gff1out);
my $gff2out = Bio::Tools::GFF->new(-gff_version => 2, -file => ">out2.gff");
ok($gff2out);

$gff1out->write_feature($feat);
$gff2out->write_feature($feat);

$gff1out->close();
$gff2out->close();

my $gff1in = Bio::Tools::GFF->new(-gff_version => 1,  -file => "out1.gff");
ok($gff1in);
my $gff2in = Bio::Tools::GFF->new(-gff_version => 2, -file => "out2.gff");
ok($gff2in);

my $feat1 = $gff1in->next_feature();
ok($feat1);
is($feat1->start, $feat->start);
is($feat1->end, $feat->end);
is($feat1->primary_tag, $feat->primary_tag);
is($feat1->score, $feat->score);

my $feat2 = $gff2in->next_feature();
ok($feat2);
is($feat2->start, $feat->start);
is($feat2->end, $feat->end);
is($feat2->primary_tag, $feat->primary_tag);
is($feat2->score, $feat->score);
is(($feat2->each_tag_value('sillytag'))[0], 'this is silly!;breakfast');

#test sequence-region parsing
$gff2in = Bio::Tools::GFF->new(-gff_version => 2, -file => Bio::Root::IO->catfile("t","data","hg16_chroms.gff"));
is($gff2in->next_feature(),undef);
my $seq = $gff2in->next_segment;
is($seq->display_id, 'chr1');
is($seq->end, 246127941);
is($seq->start, 1);


# GFF3
SKIP: {
    eval { require IO::String };
    skip('cannot verify GFF3 writing tests without IO::String installed',12)
        if $@;
    my $str = IO::String->new;
    my $gffout = new Bio::Tools::GFF(-fh => $str, -gff_version => 3);
    my $feat_test = new Bio::SeqFeature::Generic
    (-primary_tag => 'tag',
     -source_tag  => 'exon',
     -seq_id      => 'testseq',
     -score       => undef,
     -start       => 10,
     -end         => 120,
     -strand      => 1,
     -tag         => { 
         'bungle' => 'jungle;mumble',
         'lion'   => 'snake=tree'
         });
    $feat_test->add_tag_value('giant_squid', 'lakeshore manor');
    $gffout->write_feature($feat_test);
    seek($str,0,0);
    my $in = new Bio::Tools::GFF(-fh          => $str,
                 -gff_version => 3);
    my $f_recon = $in->next_feature;
    is($f_recon->primary_tag, $feat_test->primary_tag);
    is($f_recon->source_tag,  $feat_test->source_tag);
    is($f_recon->score, $feat_test->score);
    is($f_recon->start, $feat_test->start);
    is($f_recon->end, $feat_test->end);
    is($f_recon->strand, $feat_test->strand);
    for my $tag ( $feat_test->get_all_tags ) {
        ok($f_recon->has_tag($tag));
        if( $f_recon->has_tag($tag) ) {
            my @v = $feat_test->get_tag_values($tag);
            my @g = $f_recon->get_tag_values($tag);
            while( @v && @g ) {
               is(shift @v, shift @g);
            }
        }
    }
}


