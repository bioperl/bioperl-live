# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#

use strict;

BEGIN {
    use vars qw($DEBUG);
    $DEBUG = $ENV{'BIOPERLDEBUG'};
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 23;
}

END {
}

use Bio::Matrix::Generic;

my $raw = [ [ 0, 10, 20],
	    [ 2, 17,  4],
	    [ 3,  4,  5] ];

my $matrix = new Bio::Matrix::Generic(-values     => $raw,
				      -matrix_id  => 'fakeid00',
				      -matrix_name=> 'matname',
				      -rownames   => [qw(A B C)],
				      -colnames   => [qw(D E F)]);

ok($matrix->matrix_name, 'matname');
ok($matrix->matrix_id,   'fakeid00');
ok($matrix->entry('A','F'), $raw->[0]->[2]);
my @colE = $matrix->get_column('E');
ok($colE[0], $raw->[0]->[1]);
ok($colE[1], $raw->[1]->[1]);
ok($colE[2], $raw->[2]->[1]);

my @rowC = $matrix->get_row('C');
ok($rowC[0], $raw->[2]->[0]);
ok($rowC[1], $raw->[2]->[1]);
ok($rowC[2], $raw->[2]->[2]);

ok($matrix->row_num_for_name('A'),0);
ok($matrix->column_num_for_name('D'),0);

ok($matrix->row_header(1),'B');
ok($matrix->column_header(0),'D');

ok($matrix->add_row(1, 'b', [qw(21 13 14)]),4);
ok($matrix->add_column(2, 'f', [qw(71 81 14 3)]),4);

ok($matrix->row_num_for_name('B'),2);
ok($matrix->row_num_for_name('b'),1);
ok($matrix->column_num_for_name('D'),0);
ok($matrix->column_num_for_name('F'),3);
ok($matrix->column_num_for_name('f'),2);

ok($matrix->row_header(2),'B');
ok($matrix->column_header(3),'F');


ok($matrix->get_entry('b', 'f'), 81);
