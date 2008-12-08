# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 47);
	
    use_ok('Bio::Tools::Phylo::Phylip::ProtDist');
}

my $inputfilename= test_input_file('phylipdist.out');
my $tool= Bio::Tools::Phylo::Phylip::ProtDist->new(-file => $inputfilename);
my $phy = $tool->next_matrix;
is(@{$phy->names}, 5);
is $phy->get_entry('Alpha','Beta'), '4.23419';
is $phy->get_entry('Gamma','Alpha'),'3.63330';
my @column =  $phy->get_column('Alpha');
is $column[0], '0.00000';
is $column[1], '4.23419';
is $column[2], '3.63330';
is $column[3], '6.20865';
is $column[4], '3.45431';

my @row    = $phy->get_row('Gamma');
is $row[0], '3.63330';
is $row[1], '3.49289';
is $row[2], '0.00000';
is $row[3], '3.68733';
is $row[4], '5.84929';

my @diag   = $phy->get_diagonal;


is $diag[0], '0.00000';
is $diag[1], '0.00000';
is $diag[2], '0.00000';
is $diag[3], '0.00000';
is $diag[4], '0.00000';

my $matrix =<<END;
    5
Alpha          0.00000  4.23419  3.63330  6.20865  3.45431
Beta           4.23419  0.00000  3.49289  3.36540  4.29179
Gamma          3.63330  3.49289  0.00000  3.68733  5.84929
Delta          6.20865  3.36540  3.68733  0.00000  4.43345
Epsilon        3.45431  4.29179  5.84929  4.43345  0.00000
END
;
is $phy->print_matrix , $matrix;

# now parse Phylip 3.6 output

$inputfilename= test_input_file('phylipdist-36.out');
$tool= Bio::Tools::Phylo::Phylip::ProtDist->new(-file => $inputfilename);
$phy = $tool->next_matrix;

is(@{$phy->names}, 39);
is $phy->get_entry('CBG01299','CBG00435'), '4.7793';
is $phy->get_entry('CBG22788','CBG22521'),'5.3195';
is $phy->get_entry('CBG01466', 'CBG01473'), '3.3944';

@row = $phy->get_row('CBG01473');
is(scalar @row, 39);
@column =  $phy->get_column('CBG01300');
is $column[0], '0.0817';
is $column[1], '0.0000';
is $column[2], '0.0950';
is $column[3], '0.3111';
is $column[37], '4.7190';
is $column[38], '4.7592';

@row    = $phy->get_row('CBG17433');
is $row[0], '4.8451';
is $row[1], '4.5982';
is $row[2], '4.0620';
is $row[3], '5.9673';
is $row[4], '4.6224';
is $row[5], '5.1993';
is $row[6], '5.4427';
is $row[7], '4.2783';

@diag   = $phy->get_diagonal;

is $diag[0], '0.0000';
is $diag[1], '0.0000';
is $diag[2], '0.0000';
is $diag[3], '0.0000';
is $diag[4], '0.0000';
is $diag[5], '0.0000';
is $diag[37], '0.0000';
is $diag[38], '0.0000';
