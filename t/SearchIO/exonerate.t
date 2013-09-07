# -*-Perl-*- Test Harness script for Bioperl
# $Id: SearchIO_exonerate.t 11525 2007-06-27 10:16:38Z sendu $

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin( -tests => 52 );

    use_ok('Bio::SearchIO');
}

my ( $searchio, $result, $hit, $hsp );

$searchio = Bio::SearchIO->new(
    -file   => test_input_file('testdat.exonerate'),
    -format => 'exonerate'
);
my @data = (
    [
        qw(ln27 416 Contig124 939
          293 416 -1
          1   124 1

          107 292 -1
          178 363 1

          66 106 -1
          899 939 1
          )
    ],
    [
        qw(ln74 644 Contig275 1296
          601 644 -1
          901 944 1

          436 600 -1
          998 1162    1

          386 435 -1
          1247 1296 1

          )
    ]
);

my $val;

while ( my $r = $searchio->next_result ) {
    my $d = shift @data;
    is( $r->query_name, shift @$d );
  SKIP: {
        $val = shift @$d;
        skip( 'no query length available in default output', 1 );
        is( $r->query_length, $val );
    }

    my $h = $r->next_hit;
    is( $h->name, shift @$d );
  SKIP: {
        $val = shift @$d;
        skip( 'no hit length available in default output', 1 );
        is( $h->length, $val );
    }
    while ( my $hsp = $h->next_hsp ) {
        is( $hsp->query->start,  shift @$d );
        is( $hsp->query->end,    shift @$d );
        is( $hsp->query->strand, shift @$d );

        is( $hsp->hit->start,  shift @$d );
        is( $hsp->hit->end,    shift @$d );
        is( $hsp->hit->strand, shift @$d );
    }
}

# bug 2346

my @cts = (7, 1, 7);

$searchio = Bio::SearchIO->new(
   -format => 'exonerate',
   -file   => test_input_file('exonerate.output.works'),
);
parse($searchio);

$searchio = Bio::SearchIO->new(
   -format => 'exonerate',
   -file   => test_input_file('exonerate.output.dontwork'),
);
parse($searchio);

$searchio = Bio::SearchIO->new(
   -format => 'exonerate',
   -file   => test_input_file('exonerate.whitespace_before_query.works'),
);
parse($searchio);

sub parse {
    my ($searchio) = @_;
    while( my $r = $searchio->next_result ) {
        my $hsp_ct = 0;
        while(my $hit = $r->next_hit){
            while(my $hsp = $hit->next_hsp){
                $hsp_ct++;
            }
        }
        ok($r->query_name, "query_name");
        is($hsp_ct, shift @cts);
    }
}

$searchio = Bio::SearchIO->new(
    -format => 'exonerate',
    -file   => test_input_file('exonerate.output.negativescore.works'),
);
my $r   = $searchio->next_result;
$hit = $r->next_hit;
is( $hit->score, "-3" );
