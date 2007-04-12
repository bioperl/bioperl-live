# -*-Perl-*-
## Bioperl Test Harness Script for Modules


use strict;
BEGIN {
    eval { require Test::More; };
    if( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
    use vars qw($NTESTS);
    $NTESTS = 21;
    plan tests => $NTESTS;
	use_ok('Bio::Matrix::PhylipDist');
	use_ok('Bio::Tools::Phylo::Phylip::ProtDist');
}

my $inputfilename= Bio::Root::IO->catfile("t","data","phylipdist.out");
my $parser = Bio::Tools::Phylo::Phylip::ProtDist->new(-program => 'phylipdist',
						      -file => $inputfilename);

my $phy = $parser->next_matrix;
is $phy->program, 'phylipdist';
is $phy->get_entry('Alpha','Beta'), '4.23419';
is $phy->get_entry('Gamma','Alpha'),'3.63330';
my @column =  $phy->get_column('Alpha');
is $column[0], '0.00000';
is $column[1], '4.23419';
is $column[2], '3.63330';
is $column[3], '6.20865';
is $column[4], '3.45431';

my @row = $phy->get_row('Gamma');
is $row[0], '3.63330';
is $row[1], '3.49289';
is $row[2], '0.00000';
is $row[3], '3.68733';
is $row[4], '5.84929';

my @diag  = $phy->get_diagonal;


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




