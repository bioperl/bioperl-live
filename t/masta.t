# -*-Perl-*-
#Some simple tests for meme and transfac parsers

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

    plan tests => 16;
}

use Bio::Matrix::PSM::IO;
use Bio::Root::IO;

ok(1);
END {
 unlink(Bio::Root::IO->catfile(qw(t data masta_w.dat)));
}
#Let's try masta formats here
my $mio =  new Bio::Matrix::PSM::IO(-format=>'masta', 
				      -file=>Bio::Root::IO->catfile(qw(t data masta.dat)));
my $wmio=new Bio::Matrix::PSM::IO(-format=>'masta', 
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
ok $cons[0],'CAGAAAAATNGAATNCCCACCCCCC';
ok $cons[1],'CAGAAAAATAGAATCCCCACCCCCC';
ok $cons[2],'CAGAAAAATNNAATNCCCACCNCCC';

$wmio->write_psm($carry,'PWM');
$carry->id('m1freq');
$wmio->write_psm($carry,'PFM');
$carry->id('m1seq');
$wmio->write_psm($carry,'SEQ');
$wmio->DESTROY;
my $chio=new Bio::Matrix::PSM::IO(-format=>'masta', 
				      -file=>Bio::Root::IO->catfile(qw(t data masta_w.dat)));
ok $chio;
my $site=$chio->next_matrix;
ok $site->id,'m1logs';
ok $site->consensus,'CAGAAAAATAGAATCCCCACCCCCC';	
$site=$chio->next_matrix;
ok $site->id,'m1freq';
ok $site->consensus,'CAGAAAAATAGAATCCCCACCCCCC';  
$site=$chio->next_matrix;
ok $site->id,'m1seq';
ok $site->consensus,'CAGAAAAATAGAATCCCCACCCCCC';
