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
    plan tests => 14;
}

#END {
#}


#
# Testing physical maps
#


#Bio::MapIO::fpc needs testing, too

use Bio::Map::Clone;
use Bio::Map::Contig;
use Bio::Map::FPCMarker;
use Bio::Map::OrderedPositionWithDistance;


use Bio::Map::Physical;
ok 1;


ok my $phm = new Bio::Map::Physical;
ok $phm->version(2), 2;
ok $phm->version(), 2;
ok $phm->modification_user('me'), 'me';
ok $phm->modification_user(), 'me';

ok $phm->group_type('xx'), 'xx';
ok $phm->group_type(), 'xx';

ok $phm->group_abbr('xx'), 'xx';
ok $phm->group_abbr(), 'xx';

ok $phm->core_exists, undef, 'code holds and returns a string, definition requires a boolean';

ok $phm->core_exists(3), 1, 'code holds and returns a string, definition requires a boolean';

ok $phm->core_exists(1), 1;
ok $phm->core_exists(), 1;


#my @clones = $map->each_cloneid();
#my $cloneobj = $map->get_cloneobj('CLONEA');
#my $cloneobj = $map->get_cloneobj('CLONEA');
#my @markers = $map->each_markerid();

#my $markerobj = $map->get_markerobj('MARKERA');
#my @contigs = $map->each_contigid();
#my $contigobj = $map->get_contigobj('CONTIG1');
#$self->matching_bands('cloneA','cloneB',[$tol]);
#$self->coincidence_score('cloneA','cloneB'[,$tol,$gellen]);
#$map->print_contiglist([showall]); #[Default 0]
#$map->print_markerlist();
#$map->print_gffstyle([style]);
