# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 77);
	
	use_ok('Bio::Matrix::Generic');
	use_ok('Bio::Matrix::IO');
}

my $raw = [ [ 0, 10, 20],
	    [ 2, 17,  4],
	    [ 3,  4,  5] ];

my $matrix = Bio::Matrix::Generic->new(-values => $raw,
				      -matrix_id  => 'fakeid00',
				      -matrix_name=> 'matname',
				      -rownames   => [qw(A B C)],
				      -colnames   => [qw(D E F)] );

is($matrix->matrix_name, 'matname');
is($matrix->matrix_id,   'fakeid00');
is($matrix->entry('A','F'), $raw->[0]->[2]);
my @colE = $matrix->get_column('E');
is($colE[0], $raw->[0]->[1]);
is($colE[1], $raw->[1]->[1]);
is($colE[2], $raw->[2]->[1]);

my @rowC = $matrix->get_row('C');
is($rowC[0], $raw->[2]->[0]);
is($rowC[1], $raw->[2]->[1]);
is($rowC[2], $raw->[2]->[2]);

is($matrix->row_num_for_name('A'),0);
is($matrix->column_num_for_name('D'),0);

is($matrix->row_header(1),'B');
is($matrix->column_header(0),'D');

is($matrix->add_row(1, 'b', [qw(21 13 14)]),4);
is($matrix->add_column(2, 'f', [qw(71 81 14 3)]),4);

is($matrix->add_row(4, 'c', [qw(22 11 17)]),5);
is($matrix->remove_row(4),4);

is($matrix->add_column(4, 'g', [qw(11 10 100 71)]),5);
is($matrix->remove_column(4),4);

is($matrix->row_num_for_name('B'),2);
is($matrix->row_num_for_name('b'),1);

is($matrix->column_num_for_name('D'),0);
is($matrix->column_num_for_name('F'),3);
is($matrix->column_num_for_name('f'),2);

is($matrix->row_header(2),'B');

is($matrix->column_header(3),'F');

is($matrix->get_entry('b', 'f'), 81);


# read in a scoring matrix

my $io = Bio::Matrix::IO->new(-format => 'scoring',
			      -file   => test_input_file('BLOSUM50'));
my $blosum_matrix = $io->next_matrix;
isa_ok($blosum_matrix,'Bio::Matrix::Scoring');
is($blosum_matrix->entropy, 0.4808);
is($blosum_matrix->expected_score, -0.3573);
is($blosum_matrix->scale, '1/3');
is($blosum_matrix->get_entry('*','A'), -5);
is($blosum_matrix->get_entry('V','Y'), -1);
is($blosum_matrix->get_entry('Y','V'), -1);
is($blosum_matrix->get_entry('L','I'), 2);
my @diag = $blosum_matrix->get_diagonal;
is($diag[2],7);
my @row = $blosum_matrix->get_row('D');
is($row[5], $blosum_matrix->get_entry('D','Q'));
is($blosum_matrix->num_rows,24);
is($blosum_matrix->num_columns,24);
 
$io = Bio::Matrix::IO->new(-format => 'scoring',
			   -file   => test_input_file('PAM250'));
my $pam_matrix = $io->next_matrix;
isa_ok($pam_matrix, 'Bio::Matrix::Scoring');
is($pam_matrix->entropy, 0.354);
is($pam_matrix->expected_score, -0.844);
is($pam_matrix->scale, 'ln(2)/3');
is($pam_matrix->num_rows,24);
is($pam_matrix->get_entry('G','*'), -8);
is($pam_matrix->get_entry('V','Y'), -2);
is($pam_matrix->get_entry('Y','V'), -2);
is($pam_matrix->get_entry('L','I'), 2);
@diag = $pam_matrix->get_diagonal;
is($diag[2],2);
@row = $pam_matrix->get_row('D');
is($row[5], $pam_matrix->get_entry('D','Q'));

# test Phylip parsing

$io = Bio::Matrix::IO->new(-format  => 'phylip',
			  -program => 'phylipdist',
			  -file    => test_input_file('phylipdist.out'));

my $phy = $io->next_matrix;
is $phy->program, 'phylipdist';
is $phy->get_entry('Alpha','Beta'), '4.23419';
is $phy->get_entry('Gamma','Alpha'),'3.63330';

my @column =  $phy->get_column('Alpha');
is $column[0], '0.00000';
is $column[1], '4.23419';
is $column[2], '3.63330';
is $column[3], '6.20865';
is $column[4], '3.45431';

@row    = $phy->get_row('Gamma');
is $row[0], '3.63330';
is $row[1], '3.49289';
is $row[2], '0.00000';
is $row[3], '3.68733';
is $row[4], '5.84929';

@diag   = $phy->get_diagonal;

is $diag[0], '0.00000';
is $diag[1], '0.00000';
is $diag[2], '0.00000';
is $diag[3], '0.00000';
is $diag[4], '0.00000';


# test mlagan parsing

$io = Bio::Matrix::IO->new(-format => 'mlagan',
						  -file   => test_input_file('nucmatrix.txt'));

my $mlag = $io->next_matrix;
is $mlag->get_entry('A', 'C'), -150;
is $mlag->get_entry('.', 'A'), 0;
is $mlag->gap_open, -300;
is $mlag->gap_continue, -25;

# test output round-trip
$mlag->entry('A', 'C', -149);
$mlag->gap_open(-150);
$mlag->gap_continue(-5);

my $out = test_output_file();
$io = Bio::Matrix::IO->new(-format => 'mlagan',
						  -file   => ">$out");
$io->write_matrix($mlag);

$io = Bio::Matrix::IO->new(-format => 'mlagan',
						  -file   => $out);
$mlag = $io->next_matrix;
is $mlag->get_entry('A', 'C'), -149;
is $mlag->gap_open, -150;
is $mlag->gap_continue, -5;
