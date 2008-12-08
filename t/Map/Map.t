# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 150);
    
    use_ok('Bio::Map::SimpleMap');
    use_ok('Bio::Map::Marker');
    use_ok('Bio::Map::Position');
    use_ok('Bio::Map::Relative');
    use_ok('Bio::Map::Mappable');
}

###
# We explicitly test Bio::Map::SimpleMap, Bio::Map::Mappable, Bio::Map::Position,
# Bio::Map::Marker and Bio::Map::Relative.
#
# We implicitly test Bio::Map::MapI, Bio::Map::MappableI, Bio::Map::PositionI,
# and Bio::Map::PositionHandler.
###

# Test map basics
my $map;
{
    ok $map = Bio::Map::SimpleMap->new(-name  => 'my');
    ok $map->type('cyto');
    is $map->type, 'cyto';
    is $map->units, '';
    is $map->length, 0, "Length is ". $map->length;
    is $map->name, 'my';
    is $map->species('human'), 'human';
    is $map->species, 'human';
    is $map->unique_id, '1';
}

# Test marker basics
my $marker;
{
    # make a plane one and add details after
    ok $marker = Bio::Map::Marker->new();
    is $marker->name('gene1'), 'gene1';
    ok $marker->position($map, 100);
    is $marker->position->value, 100;
    is $marker->map, $map;
    
    # make positions a little easier to add by setting a default map first
    ok my $marker2 = Bio::Map::Marker->new(-name => 'gene3');
    ok $map->add_element($marker2); # one way of setting default
    is $marker2->default_map, $map;
    $marker2 = Bio::Map::Marker->new(-name => 'gene3');
    ok $marker2->default_map($map); # the other way of setting default
    is $marker2->default_map, $map;
    ok $marker2->position(300);
    is $marker2->position->value, 300;
    ok my $position = $marker2->position();
    is $position->value, 300;
    
    # make one with details set in new()
    ok my $marker3 = Bio::Map::Marker->new(-name => 'gene2', -position => [$map, 200]);
    is $marker3->default_map, $map;
    is $marker3->position->value, 200;
    
    # make one with multiple positions on multiple maps
    my $map2 = Bio::Map::SimpleMap->new();
    $marker2->positions([[$map, 150], [$map, 200], [$map2, 200], [$map2, 400]]);
    my @p = map($_->numeric, $marker2->each_position);
    is $p[0], 150;
    is $p[1], 200;
    is $p[2], 200;
    is $p[3], 300;
    is $p[4], 400;
    $marker2->purge_positions($map2);
    @p = map($_->numeric, $marker2->each_position);
    is $p[0], 150;
    is $p[1], 200;
    is $p[2], 300;
    
    # make sure we can add positions with 0 value
    my $map3 = Bio::Map::SimpleMap->new();
    $marker->add_position($map3, 0);
    ok my @positions = $marker->get_positions($map3);
    is @positions, 1;
    is $positions[0]->value, 0;
}

# Test position basics
my $pos;
{
    ok $pos = Bio::Map::Position->new();
    ok $pos->map($map);
    is $pos->map(), $map;
    ok $pos->element($marker);
    is $pos->element(), $marker;
    
    ok $pos->value('10');
    is $pos->value(), '10';
    is $pos->numeric, 10;
    is $pos->sortable, 10;
    is $pos->start, 10;
    is $pos->end, 10;
    
    # give a marker a single position with explicit position creation
    ok $pos = Bio::Map::Position->new(-map => $map, -value => 500);
    ok $marker->position($pos);
    ok my $got_pos = $marker->position();
    is $got_pos, $pos;
    is $marker->position->value, 500;
    
    # add a position
    my $map2 = Bio::Map::SimpleMap->new(-name => 'genethon', -type => 'Genetic');
    my $pos2 = Bio::Map::Position->new(-map => $map2, -value => 100);
    $marker->add_position($pos2);
    ok my @positions = $marker->get_positions($map2);
    is @positions, 1;
    is $positions[0]->value, 100;
}

# Test interaction of Markers and Maps via Positions
{
    # markers know what maps they are on
    $marker->purge_positions;
    is $marker->known_maps, 0;
    $pos->element($marker);
    is $marker->known_maps, 1;
    ok $marker->in_map(1);
    ok $marker->in_map($map);
    
    # maps know what markers are on themselves
    $map->purge_positions;
    my @els = $map->get_elements;
    is @els, 0;
    $pos->map($map);
    ok my @elements = $map->get_elements;
    is @elements, 1;
    is $elements[0], $marker;
    
    # positions know what marker they are for and what map they are on
    is $pos->map, $map;
    is $pos->element, $marker;
}

# We can compare Map objects to their own kind
{
    # positions to positions
    {
        ok $pos->equals($pos);
        # these get tested properly when testing Relative, later
    }
    
    # markers to markers
    {
        ok $marker->equals($marker);
        # these get tested properly when testing Mappables, later
    }
    
    # maps to maps
    {
        my $human = Bio::Map::SimpleMap->new();
        my $mouse = Bio::Map::SimpleMap->new();
        my $chicken = Bio::Map::SimpleMap->new();
        my $aardvark = Bio::Map::SimpleMap->new();
        
        # scenario 1: we have information about where some factors bind upstream
        # of a gene in 4 different species. Which factors are found in all the
        # species?
        my $fac1 = Bio::Map::Mappable->new();
        my $pos1 = Bio::Map::Position->new(-map => $human, -element => $fac1);
        my $pos2 = Bio::Map::Position->new(-map => $mouse, -element => $fac1);
        my $pos3 = Bio::Map::Position->new(-map => $chicken, -element => $fac1);
        my $pos4 = Bio::Map::Position->new(-map => $aardvark, -element => $fac1);
        my $fac2 = Bio::Map::Mappable->new();
        my $pos5 = Bio::Map::Position->new(-map => $human, -element => $fac2);
        my $pos6 = Bio::Map::Position->new(-map => $mouse, -element => $fac2);
        my $pos7 = Bio::Map::Position->new(-map => $chicken, -element => $fac2);
        my $fac3 = Bio::Map::Mappable->new();
        my $pos8 = Bio::Map::Position->new(-map => $human, -element => $fac3);
        my $pos9 = Bio::Map::Position->new(-map => $mouse, -element => $fac3);
        
        # scenario 1 answer:
        ok my @factors = $human->common_elements([$mouse, $chicken, $aardvark]);
        is @factors, 1;
        ok @factors = $human->common_elements([$mouse, $chicken, $aardvark], -min_percent => 50);
        is @factors, 3;
        ok @factors = $human->common_elements([$mouse, $chicken, $aardvark], -min_percent => 50, -min_num => 3);
        is @factors, 2;
        ok @factors = $chicken->common_elements([$mouse, $human, $aardvark], -min_percent => 50, -require_self => 1);
        is @factors, 2;
        ok @factors = Bio::Map::SimpleMap->common_elements([$human, $mouse, $human, $aardvark], -min_percent => 50, -required => [$aardvark]);
        is @factors, 1;
    }
}

# Test relative positions
{
    my $map = Bio::Map::SimpleMap->new();
    my $pos1 = Bio::Map::Position->new(-map => $map, -start => 50, -length => 5);
    my $pos2 = Bio::Map::Position->new(-map => $map, -start => 100, -length => 5);
    ok my $relative = Bio::Map::Relative->new(-position => $pos2);
    ok $pos1->relative($relative);
    is $pos1->start, 50;
    is $pos1->absolute(1), 1;
    is $pos1->start, 150;
    is $pos1->absolute(0), 0;
    ok my $relative2 = Bio::Map::Relative->new(-map => 10);
    my $pos3 = Bio::Map::Position->new(-map => $map, -element => $marker, -start => -5, -length => 5);
    $pos3->relative($relative2);
    my $relative3 = Bio::Map::Relative->new(-position => $pos3);
    is $pos1->start($relative3), 145;
    is $pos1->numeric($relative3), 145;
    is $pos1->end($relative3), 149;
    
    # Test the RangeI-related methods on relative positions
    {
        my $pos1 = Bio::Map::Position->new(-map => $map, -start => 50, -length => 10);
        my $pos2 = Bio::Map::Position->new(-map => $map, -start => 100, -length => 10);
        my $pos3 = Bio::Map::Position->new(-map => $map, -start => 45, -length => 1);
        my $pos4 = Bio::Map::Position->new(-map => $map, -start => 200, -length => 1);
        my $relative = Bio::Map::Relative->new(-position => $pos3);
        my $relative2 = Bio::Map::Relative->new(-position => $pos4);
        ok ! $pos1->overlaps($pos2);
        $pos1->relative($relative);
        ok $pos1->overlaps($pos2);
        ok $pos2->overlaps($pos1);
        ok $pos1->overlaps($pos2, undef, $relative2);
        
        # Make sure it works with normal Ranges
        use Bio::Range;
        my $range = Bio::Range->new(-start => 100, -end => 109);
        ok $pos1->overlaps($range);
        ok ! $range->overlaps($pos1);
        $pos1->absolute(1);
        ok $range->overlaps($pos1);
        $pos1->absolute(0);
        
        # Try the other methods briefly
        ok my $i = $pos1->intersection($pos2); # returns a mappable
        ($i) = $i->get_positions; # but we're just interested in the first (and only) position of mappable
        is $i->toString, '100..104';
        ok $i = $pos1->intersection($pos2, undef, $relative2);
        ($i) = $i->get_positions;
        is $i->toString, '-100..-96';
        is $i->map, $map;
        is $i->relative, $relative2;
        $i->absolute(1);
        is $i->toString, '100..104';
        
        ok my $u = $pos1->union($pos2);
        ($u) = $u->get_positions;
        is $u->toString, '95..109';
        ok $u = $pos1->union($pos2, $relative2);
        ($u) = $u->get_positions;
        is $u->toString, '-105..-91';
        is $u->map, $map;
        is $u->relative, $relative2;
        $u->absolute(1);
        is $u->toString, '95..109';
        
        ok ! $pos1->contains($pos2);
        $pos2->end(104);
        ok $pos1->contains($pos2);
        ok $pos1->contains(100);
        
        ok ! $pos1->equals($pos2);
        $pos2->start(95);
        ok $pos1->equals($pos2);
    }
}

# Test Mappables
{
    ok my $gene = Bio::Map::Mappable->new();
    my $human = Bio::Map::SimpleMap->new();
    my $mouse = Bio::Map::SimpleMap->new();
    ok my $pos1 = Bio::Map::Position->new(-map => $human, -element => $gene, -start => 50, -length => 1000);
    my $pos2 = Bio::Map::Position->new(-map => $mouse, -start => 100, -length => 1000);
    $gene->add_position($pos2);
    my $gene_rel = Bio::Map::Relative->new(-element => $gene);
    
    # scenario 1a: we know where a TF binds upstream of a gene in human.
    # we use four different programs to predict the site; how good were they?
    # scenaria 1b: to what extent do the predictions and known agree?
    my $factor = Bio::Map::Mappable->new();
    my $pos3 = Bio::Map::Position->new(-map => $human, -element => $factor, -start => -25, -length => 10, -relative => $gene_rel);
    my $perfect_prediction = Bio::Map::Mappable->new();
    my $pos4 = Bio::Map::Position->new(-map => $human, -element => $perfect_prediction, -start => 25, -length => 10);
    my $good_prediction = Bio::Map::Mappable->new();
    my $pos5 = Bio::Map::Position->new(-map => $human, -element => $good_prediction, -start => 24, -length => 10);
    my $ok_prediction = Bio::Map::Mappable->new();
    my $pos6 = Bio::Map::Position->new(-map => $human, -element => $ok_prediction, -start => 20, -length => 10);
    my $bad_prediction = Bio::Map::Mappable->new();
    my $pos7 = Bio::Map::Position->new(-map => $human, -element => $bad_prediction, -start => 10, -length => 10);
    
    # scenario 2: we have the same program making a prediciton for a site
    # in two different species; is the predicted site conserved in terms of
    # its position relative to the gene?
    my $human_prediction = Bio::Map::Mappable->new();
    my $pos8 = Bio::Map::Position->new(-map => $human, -element => $human_prediction, -start => 25, -length => 10);
    my $mouse_prediction = Bio::Map::Mappable->new();
    my $pos9 = Bio::Map::Position->new(-map => $mouse, -element => $mouse_prediction, -start => 75, -length => 10);
    
    # Test the RangeI-related methods
    {
        # scenario 1a answers:
        ok $perfect_prediction->equals($factor);
        ok $perfect_prediction->contains($factor);
        ok ! $ok_prediction->equals($factor);
        ok $ok_prediction->overlaps($factor);
        ok ! $bad_prediction->overlaps($factor);
        ok $bad_prediction->less_than($factor);
        ok ! $bad_prediction->greater_than($factor);
        ok $factor->greater_than($bad_prediction);
        
        # scenario 1b answer:
        my $predictions = [$perfect_prediction, $good_prediction, $ok_prediction, $bad_prediction];
        ok my @groups = $factor->overlapping_groups($predictions, -relative => $gene_rel);
        is @groups, 2;
        is ${$groups[0]}[0], $pos7;
        is ${$groups[1]}[0], $pos6;
        is ${$groups[1]}[1], $pos5;
        is ${$groups[1]}[2]->toString($gene_rel), $pos4->toString($gene_rel);
        is ${$groups[1]}[3]->toString($gene_rel), $pos3->toString($gene_rel);
        ok my $di = $factor->disconnected_intersections($predictions, -relative => $gene_rel, -min_mappables_num => 3);
        my @di = $di->get_positions;
        is @di, 1;
        is $di[0]->toString, '-25..-21';
        ok my $du = $factor->disconnected_unions($predictions, -relative => $gene_rel, -min_mappables_num => 3);
        my @du = $du->get_positions;
        is @du, 1;
        is $du[0]->toString, '-30..-16';
        
        # test the flags on overlapping_groups a bit more
        @groups = $factor->overlapping_groups($predictions, -relative => $gene_rel, -min_pos_num => 2);
        is @groups, 1;
        @groups = $factor->overlapping_groups($predictions, -relative => $gene_rel, -min_pos_num => 1, -min_mappables_num => 2);
        is @groups, 1;
        @groups = $factor->overlapping_groups($predictions, -relative => $gene_rel, -min_pos_num => 1, -min_mappables_num => 1, -min_mappables_percent => 50);
        is @groups, 1;
        @groups = $factor->overlapping_groups($predictions, -relative => $gene_rel, -min_pos_num => 1, -min_mappables_num => 1, -min_mappables_percent => 5);
        is @groups, 2;
        @groups = $factor->overlapping_groups($predictions, -relative => $gene_rel, -require_self => 1);
        is @groups, 1;
        @groups = $factor->overlapping_groups($predictions, -relative => $gene_rel, -required => [$factor]);
        is @groups, 1;
        
        # scenario 2 answer:
        ok ! $human_prediction->overlaps($mouse_prediction);
        ok $human_prediction->overlaps($mouse_prediction, -relative => $gene_rel);
    }
}
