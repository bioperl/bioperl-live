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
    $NTESTS = 19;
    plan tests => $NTESTS;
}
use Bio::Matrix::PhylipDist;

END {
    for ( $Test::ntest..$NTESTS ) {
        skip("Error in PhylipDist.pm",1);
    }
}
use Bio::Tools::Phylo::Phylip::ProtDist;

my $inputfilename= Bio::Root::IO->catfile("t","data","phylipdist.out");
my $parser = Bio::Tools::Phylo::Phylip::ProtDist->new(-program => 'phylipdist',
						      -file => $inputfilename);

my $phy = $parser->next_matrix;
ok $phy->program, 'phylipdist';
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




