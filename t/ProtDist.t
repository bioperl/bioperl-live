# -*-Perl-*-
## Bioperl Test Harness Script for Modules


use strict;
BEGIN {
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    use vars qw($NTESTS);
    $NTESTS = 46;
    plan tests => $NTESTS;
}
use Bio::Tools::Phylo::Phylip::ProtDist;

END {
    for ( $Test::ntest..$NTESTS ) {
        skip("Error in ProtDist.pm",1);
    }
}

my $inputfilename= Bio::Root::IO->catfile("t","data","phylipdist.out");
my $tool= Bio::Tools::Phylo::Phylip::ProtDist->new(-file => $inputfilename);
my $phy = $tool->next_matrix;
ok(@{$phy->names}, 5);
ok $phy->get_entry('Alpha','Beta'), '4.23419';
ok $phy->get_entry('Gamma','Alpha'),'3.63330';
my @column =  $phy->get_column('Alpha');
ok $column[0] = '0.00000';
ok $column[1] = '4.23419';
ok $column[2] = '3.63330';
ok $column[3] = '6.20865';
ok $column[4] = '3.45431';

my @row    = $phy->get_row('Gamma');
ok $row[0] = '3.63330';
ok $row[1] = '3.49289';
ok $row[2] = '0.00000';
ok $row[3] = '3.68733';
ok $row[4] = '5.84929';

my @diag   = $phy->get_diagonal;


ok $diag[0] = '0.00000';
ok $diag[1] = '0.00000';
ok $diag[2] = '0.00000';
ok $diag[3] = '0.00000';
ok $diag[4] = '0.00000';

my $matrix =<<END;
    5
Alpha          0.00000  4.23419  3.63330  6.20865  3.45431
Beta           4.23419  0.00000  3.49289  3.36540  4.29179
Gamma          3.63330  3.49289  0.00000  3.68733  5.84929
Delta          6.20865  3.36540  3.68733  0.00000  4.43345
Epsilon        3.45431  4.29179  5.84929  4.43345  0.00000
END
;
ok $phy->print_matrix , $matrix;

# now parse Phylip 3.6 output

$inputfilename= Bio::Root::IO->catfile("t","data","phylipdist-36.out");
$tool= Bio::Tools::Phylo::Phylip::ProtDist->new(-file => $inputfilename);
$phy = $tool->next_matrix;

ok(@{$phy->names}, 39);
ok $phy->get_entry('CBG01299','CBG00435'), '4.7793';
ok $phy->get_entry('CBG22788','CBG22521'),'5.3195';
ok $phy->get_entry('CBG01466', 'CBG01473'), '3.3944';

@row = $phy->get_row('CBG01473');
ok(scalar @row, 39);
@column =  $phy->get_column('CBG01300');
ok $column[0] = '0.0817';
ok $column[1] = '0.0000';
ok $column[2] = '0.0950';
ok $column[3] = '0.3111';
ok $column[37] = '4.7190';
ok $column[38] = '4.7592';

@row    = $phy->get_row('CBG17433');
ok $row[0] = '4.8451';
ok $row[1] = '4.5982';
ok $row[2] = '4.0620';
ok $row[3] = '5.9673';
ok $row[4] = '4.6224';
ok $row[5] = '5.1993';
ok $row[6] = '5.4427';
ok $row[7] = '4.2783';

@diag   = $phy->get_diagonal;

ok $diag[0] = '0.00000';
ok $diag[1] = '0.00000';
ok $diag[2] = '0.00000';
ok $diag[3] = '0.00000';
ok $diag[4] = '0.00000';
ok $diag[5] = '0.00000';
ok $diag[37] = '0.00000';
ok $diag[38] = '0.00000';




