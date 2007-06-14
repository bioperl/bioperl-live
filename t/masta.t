# -*-Perl-*-
#Some simple tests for meme and transfac parsers

use strict;

BEGIN {
    eval { require Test::More; };
    if( $@ ) {
	use lib 't/lib';
    }
    use Test::More;

    plan tests => 17;
	use_ok('Bio::Matrix::PSM::IO');
	use_ok('Bio::Root::IO');
}

END {
 unlink(Bio::Root::IO->catfile(qw(t data masta_w.dat)));
}

#Let's try masta formats here
my $mio =  Bio::Matrix::PSM::IO->new(-format=>'masta', 
				      -file=>Bio::Root::IO->catfile(qw(t data masta.dat)));
my $wmio=Bio::Matrix::PSM::IO->new(-format=>'masta', 
				      -file=>">".Bio::Root::IO->catfile(qw(t data masta_w.dat)));
$wmio->_flush_on_write(1);
ok $mio;
ok $wmio;
my @cons;
my $carry;
while (my $site=$mio->next_matrix) {
	ok $site;
	push @cons,$site->consensus;	
	$carry=$site if ($site->id eq 'm1logs');
} 
is $cons[0],'CAGAAAAATNGAATNCCCACCCCCC';
is $cons[1],'CAGAAAAATAGAATCCCCACCCCCC';
is $cons[2],'CAGAAAAATNNAATNCCCACCNCCC';

$wmio->write_psm($carry,'PWM');
$carry->id('m1freq');
$wmio->write_psm($carry,'PFM');
$carry->id('m1seq');
$wmio->write_psm($carry,'SEQ');
$wmio->DESTROY;
my $chio=Bio::Matrix::PSM::IO->new(-format=>'masta', 
				      -file=>Bio::Root::IO->catfile(qw(t data masta_w.dat)));
ok $chio;
my $site=$chio->next_matrix;
is $site->id,'m1logs';
is $site->consensus,'CAGAAAAATAGAATCCCCACCCCCC';	
$site=$chio->next_matrix;
is $site->id,'m1freq';
is $site->consensus,'CAGAAAAATAGAATCCCCACCCCCC';  
$site=$chio->next_matrix;
is $site->id,'m1seq';
is $site->consensus,'CAGAAAAATAGAATCCCCACCCCCC';
