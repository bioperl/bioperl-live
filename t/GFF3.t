#!/usr/bin/perl -w

# author - chris mungall cjm@fruitfly.org

use strict;

BEGIN {     
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan test => 4;
}

use Bio::Root::IO;
use Bio::Tools::GFF;

# doesn't test writing yet

my $gfffn = Bio::Root::IO->catfile("t","data","flybase.gff3");
my $gffio = Bio::Tools::GFF->new(-file => $gfffn);
my @features = $gffio->features;
$gffio->close();

printf "n top features=%s\n", scalar(@features);
ok(scalar(@features) == 7);
my ($gene) = grep { $_->unique_id && $_->unique_id eq 'gene00001' } @features;
ok($gene->type->name eq 'gene');

my @transcripts = $gene->sub_SeqFeature;
printf "n trs=%s\n", scalar(@transcripts);

ok(grep {$_->type->name eq 'mRNA'} @transcripts);

foreach my $t (@transcripts) {
    my @s = $t->sub_SeqFeature;
    foreach (@s) {
        printf "%s %s PART-OF %s %s\n", 
          $_->unique_id, $_->type->name,
            $t->unique_id, $t->type->name,
            ;
    }
}

my @subf = map {$_->sub_SeqFeature} @transcripts;

ok(grep {$_->type->name eq 'cds' or $_->type->name eq 'exon'} @subf);

