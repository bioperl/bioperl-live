# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

my $error = 0;

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
    plan tests => 23; 

#    eval { require XML::Parser::PerlSAX; };
#    if( $@ ) {
#	print STDERR "XML::Parser::PerlSAX not loaded. This means TreeIO::phyloxml test cannot be executed. Skipping\n";
#	foreach ( 1..43 ) {
#	    skip(1,1);
#	}
#       $error = 1;
#	
#    } 

}

if( $error == 1 ) {
    exit(0);
}

use vars qw($FILE1 $FILE2);

$FILE1= 'testmap.map';

END { unlink $FILE1; }
use Bio::MapIO;
use Bio::Root::IO;
my $verbose = $ENV{'BIOPERLDEBUG'};
ok(1);

my $mapio = new Bio::MapIO('-verbose' => $verbose,
			    '-format' => 'mapmaker',
			    '-file'   => Bio::Root::IO->catfile('t','data', 
			    'mapmaker.out'));

ok($mapio);
my $map = $mapio->next_map;

ok(ref($map) && $map->isa('Bio::Map::MapI'));

ok($map->units, 'cM');
ok($map->type, 'Genetic');
ok($map->name('test map'), 'test map'); # map name is unset for this data type

my $count = 1;
foreach my $marker ( $map->each_element ) {
    ok(($marker->position->each_position)[0],$count++);
}
