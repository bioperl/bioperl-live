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
    plan tests => 75;
}

use Bio::Matrix::Generic;
use Bio::Matrix::IO;
use Bio::Root::IO;

END {
	unlink(Bio::Root::IO->catfile("t","data","nucmatrix.out"));
}

my $raw = [ [ 0, 10, 20],
	    [ 2, 17,  4],
	    [ 3,  4,  5] ];

my $matrix = new Bio::Matrix::Generic(-values => $raw,
				      -matrix_id  => 'fakeid00',
				      -matrix_name=> 'matname',
				      -rownames   => [qw(A B C)],
				      -colnames   => [qw(D E F)] );

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

ok($matrix->add_row(4, 'c', [qw(22 11 17)]),5);
ok($matrix->remove_row(4),4);

ok($matrix->add_column(4, 'g', [qw(11 10 100 71)]),5);
ok($matrix->remove_column(4),4);

ok($matrix->row_num_for_name('B'),2);
ok($matrix->row_num_for_name('b'),1);

ok($matrix->column_num_for_name('D'),0);
ok($matrix->column_num_for_name('F'),3);
ok($matrix->column_num_for_name('f'),2);

ok($matrix->row_header(2),'B');

ok($matrix->column_header(3),'F');

ok($matrix->get_entry('b', 'f'), 81);


# read in a scoring matrix

my $io = Bio::Matrix::IO->new(-format => 'scoring',
			      -file   => Bio::Root::IO->catfile
			      (qw(t data BLOSUM50)));
my $blosum_matrix = $io->next_matrix;
ok($blosum_matrix->isa('Bio::Matrix::Scoring'));
ok($blosum_matrix->entropy, 0.4808);
ok($blosum_matrix->expected_score, -0.3573);
ok($blosum_matrix->scale, '1/3');
ok($blosum_matrix->get_entry('*','A'), -5);
ok($blosum_matrix->get_entry('V','Y'), -1);
ok($blosum_matrix->get_entry('Y','V'), -1);
ok($blosum_matrix->get_entry('L','I'), 2);
my @diag = $blosum_matrix->get_diagonal;
ok($diag[2],7);
my @row = $blosum_matrix->get_row('D');
ok($row[5], $blosum_matrix->get_entry('D','Q'));
ok($blosum_matrix->num_rows,24);
ok($blosum_matrix->num_columns,24);
 
$io = Bio::Matrix::IO->new(-format => 'scoring',
			   -file   => Bio::Root::IO->catfile
			   (qw(t data PAM250)));
my $pam_matrix = $io->next_matrix;
ok($pam_matrix->isa('Bio::Matrix::Scoring'));
ok($pam_matrix->entropy, 0.354);
ok($pam_matrix->expected_score, -0.844,);
ok($pam_matrix->scale, 'ln(2)/3');
ok($pam_matrix->num_rows,24);
ok($pam_matrix->get_entry('G','*'), -8);
ok($pam_matrix->get_entry('V','Y'), -2);
ok($pam_matrix->get_entry('Y','V'), -2);
ok($pam_matrix->get_entry('L','I'), 2);
@diag = $pam_matrix->get_diagonal;
ok($diag[2],2);
@row = $pam_matrix->get_row('D');
ok($row[5], $pam_matrix->get_entry('D','Q'));

# test Phylip parsing

$io = new Bio::Matrix::IO(-format  => 'phylip',
			  -program => 'phylipdist',
			  -file    => Bio::Root::IO->catfile
			  (qw(t data phylipdist.out)));

my $phy = $io->next_matrix;
ok $phy->program, 'phylipdist';
ok $phy->get_entry('Alpha','Beta'), '4.23419';
ok $phy->get_entry('Gamma','Alpha'),'3.63330';

my @column =  $phy->get_column('Alpha');
ok $column[0], '0.00000';
ok $column[1], '4.23419';
ok $column[2], '3.63330';
ok $column[3], '6.20865';
ok $column[4], '3.45431';

@row    = $phy->get_row('Gamma');
ok $row[0], '3.63330';
ok $row[1], '3.49289';
ok $row[2], '0.00000';
ok $row[3], '3.68733';
ok $row[4], '5.84929';

@diag   = $phy->get_diagonal;

ok $diag[0], '0.00000';
ok $diag[1], '0.00000';
ok $diag[2], '0.00000';
ok $diag[3], '0.00000';
ok $diag[4], '0.00000';


# test mlagan parsing

$io = new Bio::Matrix::IO(-format => 'mlagan',
						  -file   => Bio::Root::IO->catfile(qw(t data nucmatrix.txt)));

my $mlag = $io->next_matrix;
ok $mlag->get_entry('A', 'C'), -150;
ok $mlag->get_entry('.', 'A'), 0;
ok $mlag->gap_open, -300;
ok $mlag->gap_continue, -25;

# test output round-trip
$mlag->entry('A', 'C', -149);
$mlag->gap_open(-150);
$mlag->gap_continue(-5);

my $out = Bio::Root::IO->catfile(qw(t data nucmatrix.out));
$io = new Bio::Matrix::IO(-format => 'mlagan',
						  -file   => ">$out");
$io->write_matrix($mlag);

$io = new Bio::Matrix::IO(-format => 'mlagan',
						  -file   => $out);
$mlag = $io->next_matrix;
ok $mlag->get_entry('A', 'C'), -149;
ok $mlag->gap_open, -150;
ok $mlag->gap_continue, -5;