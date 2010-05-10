# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 267);
    
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

# complex (multi-mappable, multi-map) test of disconnected_*
# and test Bio::Map::GeneMap, Bio::Map::Gene, Bio::Map::TranscriptionFactor,
# Bio::Map::GeneRelative, Bio::Map::GenePosition and Bio::Map::Prediction
use_ok('Bio::Map::Gene');
use_ok('Bio::Map::GeneMap');
use_ok('Bio::Map::TranscriptionFactor');
use_ok('Bio::Map::GeneRelative');
use_ok('Bio::Map::GenePosition');
use_ok('Bio::Map::Prediction');
{
    my @genes;
    my @predictions;
    
    $genes[0] = Bio::Map::Gene->get(-universal_name => 'gene1');
    $genes[1] = Bio::Map::Gene->get(-universal_name => 'gene2');
    $genes[2] = Bio::Map::Gene->get(-universal_name => 'gene3');
    Bio::Map::GeneMap->get(-gene => 'gene1',
                           -species => 'species1',
                           -upstream => 1000);
    Bio::Map::GeneMap->get(-gene => 'gene1',
                           -species => 'species2',
                           -upstream => 2000);
    Bio::Map::GeneMap->get(-gene => 'gene2',
                           -species => 'species1',
                           -upstream => 1000);
    Bio::Map::GeneMap->get(-gene => 'gene2',
                           -species => 'species2',
                           -upstream => 2000);
    Bio::Map::GeneMap->get(-gene => 'gene3',
                           -species => 'species1',
                           -upstream => 1000);
    Bio::Map::GeneMap->get(-gene => 'gene3',
                           -species => 'species2',
                           -upstream => 2000);
    
    $predictions[0] = Bio::Map::Prediction->new(-source => 'meme');
    Bio::Map::Position->new(-element => $predictions[0],
                            -map => Bio::Map::GeneMap->get(-gene => 'gene1', -species => 'species1'),
                            -start => 950,
                            -end => 960);
    Bio::Map::Position->new(-element => $predictions[0],
                            -map => Bio::Map::GeneMap->get(-gene => 'gene1', -species => 'species2'),
                            -start => 1950,
                            -end => 1960);
    Bio::Map::Position->new(-element => $predictions[0],
                            -map => Bio::Map::GeneMap->get(-gene => 'gene2', -species => 'species1'),
                            -start => 955,
                            -end => 965);
    Bio::Map::Position->new(-element => $predictions[0],
                            -map => Bio::Map::GeneMap->get(-gene => 'gene2', -species => 'species2'),
                            -start => 1955,
                            -end => 1965);
    $predictions[1] = Bio::Map::Prediction->new(-source => 'meme');
    Bio::Map::Position->new(-element => $predictions[1],
                            -map => Bio::Map::GeneMap->get(-gene => 'gene1', -species => 'species1'),
                            -start => 950,
                            -end => 960);
    Bio::Map::Position->new(-element => $predictions[1],
                            -map => Bio::Map::GeneMap->get(-gene => 'gene1', -species => 'species2'),
                            -start => 1950,
                            -end => 1960);
    Bio::Map::Position->new(-element => $predictions[1],
                            -map => Bio::Map::GeneMap->get(-gene => 'gene3', -species => 'species1'),
                            -start => 956,
                            -end => 966);
    Bio::Map::Position->new(-element => $predictions[1],
                            -map => Bio::Map::GeneMap->get(-gene => 'gene3', -species => 'species2'),
                            -start => 1956,
                            -end => 1966);
    
    Bio::Map::Position->new(-element => $predictions[0],
                            -map => Bio::Map::GeneMap->get(-gene => 'gene2', -species => 'species2'),
                            -start => 19850,
                            -end => 19860);
    Bio::Map::Position->new(-element => $predictions[1],
                            -map => Bio::Map::GeneMap->get(-gene => 'gene2', -species => 'species2'),
                            -start => 19850,
                            -end => 19860);
    
    my $rel = Bio::Map::GeneRelative->new(-gene => 0);
    my $di = Bio::Map::Mappable->disconnected_intersections(\@predictions,
                                                            -min_mappables_percent => 100,
                                                            -min_map_percent => 100,
                                                            -relative => $rel);
    my @positions = $di->get_positions;
    my $expected;
    $expected->{gene1}->{species1} = '-45..-41';
    $expected->{gene1}->{species2} = '-45..-41';
    $expected->{gene2}->{species1} = '-45..-41';
    $expected->{gene2}->{species2} = '-45..-41';
    $expected->{gene3}->{species1} = '-45..-41';
    $expected->{gene3}->{species2} = '-45..-41';
    foreach my $pos (@positions) {
        my $map = $pos->map;
        my $gname = $map ? $map->gene->universal_name : 'n/a';
        my $species = $map ? $map->species : 'n/a';
        if (defined $expected->{$gname}->{$species}) {
            is $expected->{$gname}->{$species}, $pos->toString;
            delete $expected->{$gname}->{$species};
        }
        unless (keys %{$expected->{$gname}} > 0) {
            delete $expected->{$gname};
        }
    }
    is scalar(keys %{$expected}), 0;
    
    # check we don't have any extraneous positions
    $expected = 8;
    foreach my $map ($genes[0]->known_maps) {
        foreach my $pos ($map->get_positions) {
            $expected--;
        }
    }
    is $expected, 0;
    $expected = 8;
    foreach my $map ($genes[1]->known_maps) {
        foreach my $pos ($map->get_positions) {
            $expected--;
        }
    }
    is $expected, 0;
}

{
	# make some maps that will represent an area around a particular gene in
	# particular species
    ok my $map1 = Bio::Map::GeneMap->get(-gene => 'BRCA2',
                                         -species => 'human',
                                         -description => 'breast cancer 2, early onset');
    ok my $gene = $map1->gene;
	my $map2 = Bio::Map::GeneMap->get(-gene => $gene,
                                      -species => 'mouse',
                                      -upstream => 500);
    is $map1->gene, $map2->gene;
    is $gene->universal_name, 'BRCA2';
    is $gene->description, 'breast cancer 2, early onset';
    is $map1->universal_name, 'BRCA2';
    is $map1->upstream, 1000;
    is $map2->upstream, 500;
    my $map3 = Bio::Map::GeneMap->get(-gene => 'BRCA2', -species => 'human');
    is $map3, $map1;
    is $map3->gene, $gene;
	
	# model a TF that binds 500bp upstream of the BRCA2 gene in humans and
	# 250bp upstream of BRCA2 in mice
	ok my $rel = Bio::Map::GeneRelative->new(-description => "gene start");
    ok my $tf = Bio::Map::TranscriptionFactor->get(-universal_name => 'tf1');
    is $tf->universal_name, 'tf1';
	Bio::Map::Position->new(-map => $map1,
                            -element => $tf,
                            -start => -450,
                            -length => 10,
                            -relative => $rel);
	Bio::Map::Position->new(-map => $map2,
                            -element => $tf,
                            -start => -200,
                            -length => 10,
                            -relative => $rel);
	
	# find out all the things that map near BRCA2 in all species
    my %answers = (human => ['human', 'tf1', -450, 'gene start', 551], mouse => ['mouse', 'tf1', -200, 'gene start', 301]);
	foreach my $map ($gene->known_maps) {
        my @answers = @{$answers{$map->species}};
		foreach my $thing ($map->get_elements) {
            next if $thing eq $gene;
            foreach my $pos ($thing->get_positions($map)) {
                is $map->species, shift @answers;
                is $thing->universal_name, shift @answers;
                is $pos->value, shift @answers;
                is $pos->relative->description, shift @answers;
                $pos->absolute(1);
                is $pos->start, shift @answers;
            }
		}
        is @answers, 0;
        delete $answers{$map->species};
	}
    is keys %answers, 0;
    
    # test getting abolute positions of things relative to things that don't
    # even exist in the map yet: 1st exon of default transcript
    ok $rel = Bio::Map::GeneRelative->new(-exon => [1]);
    my $pos = Bio::Map::Position->new(-map => $map1,
                                      -element => $tf,
                                      -start => 5,
                                      -length => 5,
                                      -relative => $rel);
    is $pos->start, 5;
    is $pos->relative->absolute_conversion($pos), 1006;
    is $pos->start($pos->absolute_relative), 1006;
    $pos->absolute(1);
    is $pos->start, 1006;
    is $pos->end, 1010;
    
    # now actually add some transcripts, exons, introns, coding etc. and retest
    ok my $trans = Bio::Map::GenePosition->new(-start => 0, -length => 700, -map => $map1, -type => 'transcript');
    ok $gene->add_transcript_position($trans);
    ok my $gene_pos = $gene->position($map1);
    is $gene_pos->start, 1001;
    is $gene_pos->end, 1700; # the gene position is big enough to hold the transcript
    
    $trans = Bio::Map::GenePosition->new(-start => 100, -length => 800, -map => $map1, -type => 'transcript');
    $gene->add_transcript_position($trans);
    is $gene_pos->end, 1900;
    is $gene->active_transcript($map1), 2;
    my @t_pos = $gene->get_transcript_positions($map1);
    is @t_pos, 2;
    
    # pos was relative to the active transcript, which has now changed...
    is $pos->start, 1106;
    
    # make a new one relative to an explicit transcript
    $rel = Bio::Map::GeneRelative->new(-exon => [1, 2]);
    my $pos2 = Bio::Map::Position->new(-map => $map1,
                            -element => $tf,
                            -start => 15,
                            -length => 15,
                            -relative => $rel);
    is $pos2->start, 15;
    $pos2->absolute(1);
    is $pos2->start, 1116;
    is $pos2->end, 1130;
    
    # which isn't affected by changing the active
    is $gene->active_transcript($map1, 1), 1;
    is $pos->start, 1006;
    is $pos2->start, 1116;
    
    $map1->get_position_handler->purge_positions($pos2);
    
    # add some exons to the first transcript
    ok my $exon = Bio::Map::GenePosition->new(-start => 0, -length => 100, -map => $map1, -type => 'exon');
    $gene->add_exon_position($exon, 1);
    is $pos->start, 1006;
    $exon->start(30); # not something you'd normally do; just for demo purposes
    is $exon->relative->absolute_conversion($exon), 1031;
    is $pos->start, 1036;
    
    # add another exon before the previous one - this will be considered exon 1
    my $exon1 = Bio::Map::GenePosition->new(-start => 0, -length => 20, -map => $map1, -type => 'exon');
    $gene->add_exon_position($exon1, 1);
    is $gene->get_exon_position($map1, 2), $exon;
    ok my @exons = $gene->get_exon_positions($map1);
    is @exons, 2;
    is $exons[0], $exon1;
    is $exons[1], $exon;
    is $pos->start, 1006;
    
    # add the intervening intron
    ok my $intron = Bio::Map::GenePosition->new(-start => 20, -length => 10, -map => $map1, -type => 'intron');
    ok ! $gene->get_intron_position($map1, 1);
    $gene->add_intron_position($intron, 1);
    is $gene->get_intron_position($map1, 1), $intron;
    ok ! $gene->get_intron_position($map1, 1, 2);
    ok ! $gene->get_intron_position($map1, 2);
    is $gene->get_intron_positions($map1), 1;
    is $intron->relative->absolute_conversion($intron), 1021;
    
    # like for exon 1, we can be relative to the coding region without
    # defining it
    $rel = Bio::Map::GeneRelative->new(-coding => 0);
    my $pos3 = Bio::Map::Position->new(-map => $map1,
                            -element => $tf,
                            -start => -10,
                            -length => 5,
                            -relative => $rel);
    is $pos3->start, -10;
    $pos3->absolute(1);
    is $pos3->start, 991;
    
    # add the coding region for transcript 1
    ok my $coding1a = Bio::Map::GenePosition->new(-start => 50, -length => 20, -map => $map1, -type => 'coding');
    $gene->coding_position($coding1a);
    is $pos3->start, 1041;
    
    # try adding another coding region to the same transcript: we can't, so
    # the existing coding is replaced
    my $coding1b = Bio::Map::GenePosition->new(-start => 60, -length => 20, -map => $map1, -type => 'coding');
    $gene->coding_position($coding1b);
    is $pos3->start, 1051;
    ok ! $coding1a->element;
    # try adding things without using the add_x_position methods of Gene
    #...
    
    # GenePositions can have sequence
    like $exon->seq, qr/N{70}/;
    my $pos4 = Bio::Map::GenePosition->new(-start => 200, -seq => 'ATGCCCAAAG', -map => $map1, -type => 'exon');
    is $pos4->seq, 'ATGCCCAAAG';
    $gene->add_exon_position($pos4, 1);
    is $gene->get_exon_positions($map1), 3;
    is $pos4->length, 10;
    $pos4->absolute(1);
    is $pos4->end, 1210;
    is $pos4->seq('ATGCC'), 'ATGCC';
    is $pos4->length, 5;
    is $pos4->end, 1205;
    
    # so can GeneMaps
    my $map4 = Bio::Map::GeneMap->get(-gene => $gene,
                                      -species => 'chicken',
                                      -seq => 'ATGCCCAAAG');
    like $map4->seq, qr/ATGCCCAAAGN{991}/;
    is $map4->subseq(3,6), 'GCCC';
    is $map4->subseq(9,15), 'AGNNNNN'; # subseq is within map but beyond supplied seq, pads with Ns
    
    # map sequence can be built from the sequence of PositionWithSequences on the map
    my $pos5 = Bio::Map::PositionWithSequence->new(-start => 11, -seq => 'ATG', -map => $map4);
    is $pos5->seq, 'ATG';
    is $map4->subseq(9,15), 'AGATGNN';
    
    SKIP: {
        test_skip(-tests => 19,
                  -requires_modules => [qw(Bio::Tools::Run::Ensembl Bio::EnsEMBL::Registry)],
                  -requires_networking => 1);
        
        # make transcript, coding, exon, intron positions on all maps for gene,
        # purging manually added GenePositions
        my $success = $gene->set_from_db;
        
        skip('Failed to retreive anything from Ensembl; not sure why', 19) unless $success;
        
        is $gene->get_transcript_position($map1)->toString($pos->absolute_relative), '1001..84737';
        is $gene->get_transcript_position($map2)->toString($pos->absolute_relative), '501..47617';
        is $gene->get_transcript_position($map4)->toString($pos->absolute_relative), '1373..37665';
        like $gene->description($map1), qr/Breast cancer/i;
        is $gene->display_id($map1), 'ENSG00000139618';
        is $gene->display_id($map2), 'ENSMUSG00000041147';
        is $gene->display_id($map4), 'ENSGALG00000017073';
        is $gene->display_xref($map4), 'NP_989607.1';
        is $gene->external_name($map1), 'BRCA2';
        is $gene->biotype($map2), 'protein_coding';
        is $gene->source($map4), 'ensembl';
        
        # we can add to a new map and the database info will be automatically there
        my $species = Bio::Species->new(-name => 'dog');
        $species->db_handle(Bio::DB::Taxonomy->new(-source => 'entrez'));
        my $map5 = Bio::Map::GeneMap->get(-gene => $gene, -species => $species);
        is $gene->display_id($map5), 'ENSCAFG00000006383';
        
        # now the gene has a database connection, its maps and positions can get sequence
        ok my $seq = $map1->seq;
        is length($seq), 84737;
        is substr($seq, 0, 20), 'TGTTACAGAACCAACGAATT'; # start of upstream
        is substr($seq, -20, 20), 'CTACAAGTATTATTTTACAA'; # end of gene since no downstream
        is substr($map1->subseq($gene->coding_position($map1)), 0, 3), 'ATG';
        my $exon1_str = 'GGGCTTGTGGCGCGAGCTTCTGAAACTAGGCGGCAGAGGCGGAGCCGCTGTGGCACTGCTGCGCCTCTGCTGCG';
        my $exon1_pos = $gene->get_exon_position($map1, 1);
        is $map1->subseq($exon1_pos), $exon1_str;
        is $exon1_pos->seq, $exon1_str;
    }
    
    # test a gene with multiple transcripts...
    #...
}

# test predictor<->map integration
# this seems to work and give an ok looking result, but the tests are very slow
# to complete and the results as yet unvalidated
if (0) {
    # we respresent all the genes of interest in species of interest as Genes on
    # GeneMaps with attached upstream sequence
    use_ok('Bio::SeqIO');
    my @genes;
    my $seqs_dir = test_input_file('map_hem');
    opendir(SEQSDIR, $seqs_dir) or die "couldn't open seqs dir '$seqs_dir'\n";
    while (my $thing = readdir(SEQSDIR)) {
        next unless $thing =~ /(\w+\d+)\.ups\.fa_$/;
        my $gene_name = $1;
        
        my $gene = Bio::Map::Gene->get(-universal_name => $gene_name);
        push(@genes, $gene);
        
        my $seqin = Bio::SeqIO->new(-file => "$seqs_dir/$thing", -format => "fasta");
        my $seqout = Bio::SeqIO->new(-file => ">$seqs_dir/$thing.revcom", -format => "fasta");
        while (my $seq = $seqin->next_seq) {
            my $id = $seq->id;
            my ($species) = $id =~ /^[a-z]+_([a-z]+)/i;
            
            # skip Scas and Sklu
            next if ($species eq 'Scas' || $species eq 'Sklu');
            $seqout->write_seq($seq->revcom);
            
            # (we don't need to do anything with the return value of this get()
            # since we've already stored the Bio::Map::Gene above, which knows
            # about the maps it is on as soon as they're created)
            Bio::Map::GeneMap->get(-gene => $gene_name,
                                   -species => $species,
                                   -seq => $seq->revcom->seq,
                                   -upstream => $seq->length);
        }
    }
    closedir(SEQSDIR);
    
    # then we supply these Genes to a prediction program, in this case Meme,
    # adding the predictions to the maps (we're doing all pairwise-combos of
    # genes here)
    my @predicitons_of_interest;
    @genes = sort { $a->universal_name cmp $b->universal_name } @genes;

    #use_ok('Bio::Tools::Run::Meme');
    #my %params = (dna => 1, mod => "oops", revcomp => 1, nmotifs => 5, bfile => "$seqs_dir/yeast.nc.1.freq", maxw => 20);
    #my $factory = Bio::Tools::Run::Meme->new(%params);
    #my %predictions; # not used for anything here, but this is how we might store all our predictions
    #
    #for my $i (0..$#genes) {
    #    my $gene_a = $genes[$i];
    #    my $gene_name_a = $gene_a->universal_name;
    #    
    #    my $done = 0;
    #    for my $j (($i+1)..$#genes) {
    #        my $gene_b = $genes[$j];
    #        my $gene_name_b = $gene_b->universal_name;
    #        
    #        my $prediction = $factory->run($gene_a, $gene_b);
    #        $prediction->name("$gene_name_a vs $gene_name_b in yeasts");
    #        print "got prediction for ", $prediction->name, ":\n";
    #        foreach my $pos ($prediction->get_positions) {
    #            print "  pos ", $pos->toString, " on map for gene ", $pos->map->gene->universal_name, " and species ", $pos->map->species, "\n";
    #        }
    #        
    #        $predictions{$gene_name_a}->{$gene_name_b} = $prediction;
    #        push(@predicitons_of_interest, $prediction) if $gene_name_a eq 'HEM1';
    #        
    #        $done++;
    #        #last if $done == 2; # temp for testing
    #    }
    #    
    #    #last; # shortcut for testing, since we're only interested in predictions featuring MEM1 anyway
    #}
    
    # shortcut for testing - avoid actually calling meme, just create the results directly
    # $pred_num = 0; while (<>) { if (/got prediction/) { $pred_num++; print "my \$pred$pred_num = Bio::Map::Prediction->new(-source => \"meme\");\n"; } else { ($start, $end, $gene, $spe) = $_ =~ /(\d+)\.\.(\d+) on map for gene (\S+) and species (\S+)/; print "Bio::Map::Position->new(-element => \$pred$pred_num, -map => Bio::Map::GeneMap->get(-gene => \"$gene\", -species => \"$spe\"), -start => $start, -end => $end);\n"; } }'
    my $pred1 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 51, -end => 70);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 60, -end => 79);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 90, -end => 109);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 96, -end => 115);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 97, -end => 115);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 100, -end => 118);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 100, -end => 118);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 102, -end => 120);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 103, -end => 121);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 117, -end => 136);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 120, -end => 139);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 120, -end => 139);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 122, -end => 141);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 123, -end => 142);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 133, -end => 152);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 173, -end => 192);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 175, -end => 194);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 176, -end => 195);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 178, -end => 197);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 180, -end => 199);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 197, -end => 216);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 199, -end => 218);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 199, -end => 218);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 202, -end => 221);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 204, -end => 223);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 217, -end => 236);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 219, -end => 238);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 275, -end => 293);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 275, -end => 293);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 299, -end => 318);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 299, -end => 318);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 302, -end => 321);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 303, -end => 322);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 303, -end => 322);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 304, -end => 323);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 331, -end => 350);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 335, -end => 354);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 337, -end => 356);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 485, -end => 503);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 489, -end => 507);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 604, -end => 623);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 728, -end => 747);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 796, -end => 815);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 807, -end => 826);
    Bio::Map::Position->new(-element => $pred1, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 810, -end => 829);
    my $pred2 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 108, -end => 127);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 111, -end => 130);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 111, -end => 130);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 113, -end => 132);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 114, -end => 133);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 128, -end => 147);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 131, -end => 150);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 131, -end => 150);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 133, -end => 152);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 134, -end => 153);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 167, -end => 186);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 169, -end => 188);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 170, -end => 189);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 172, -end => 191);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 174, -end => 193);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 196, -end => 215);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 198, -end => 217);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 198, -end => 217);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 201, -end => 220);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 203, -end => 222);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 272, -end => 291);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 277, -end => 296);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 277, -end => 296);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 278, -end => 297);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 282, -end => 301);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 386, -end => 405);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 388, -end => 407);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 388, -end => 407);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 394, -end => 413);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 642, -end => 661);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 643, -end => 662);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 645, -end => 664);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 650, -end => 669);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 669, -end => 688);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 670, -end => 689);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 673, -end => 692);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 676, -end => 695);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 710, -end => 729);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 711, -end => 730);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 717, -end => 736);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 721, -end => 740);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 805, -end => 824);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 807, -end => 826);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 811, -end => 830);
    Bio::Map::Position->new(-element => $pred2, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 815, -end => 834);
    my $pred3 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 15, -end => 34);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 50, -end => 69);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 50, -end => 69);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 65, -end => 84);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 70, -end => 89);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 70, -end => 89);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 72, -end => 91);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 94, -end => 113);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 101, -end => 120);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 104, -end => 123);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 104, -end => 123);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 106, -end => 125);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 107, -end => 126);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 107, -end => 126);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 123, -end => 142);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 126, -end => 145);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 126, -end => 145);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 128, -end => 147);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 129, -end => 148);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 132, -end => 151);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 133, -end => 152);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 138, -end => 157);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 169, -end => 188);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 170, -end => 189);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 175, -end => 194);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 177, -end => 196);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 191, -end => 210);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 192, -end => 211);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 192, -end => 211);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 194, -end => 213);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 194, -end => 213);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 197, -end => 216);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 199, -end => 218);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 199, -end => 218);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 262, -end => 281);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 272, -end => 291);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 277, -end => 296);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 277, -end => 296);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 278, -end => 297);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 282, -end => 301);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 316, -end => 335);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 319, -end => 338);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 320, -end => 339);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 320, -end => 339);
    Bio::Map::Position->new(-element => $pred3, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 321, -end => 340);
    my $pred4 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 106, -end => 125);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 109, -end => 128);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 109, -end => 128);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 111, -end => 130);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 112, -end => 131);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 142, -end => 161);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 147, -end => 166);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 153, -end => 172);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 155, -end => 174);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 169, -end => 188);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 202, -end => 221);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 211, -end => 230);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 212, -end => 231);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 214, -end => 233);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 228, -end => 247);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 261, -end => 280);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 281, -end => 300);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 283, -end => 302);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 286, -end => 305);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 286, -end => 305);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 288, -end => 307);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 306, -end => 325);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 308, -end => 327);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 310, -end => 329);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 311, -end => 330);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 313, -end => 332);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 346, -end => 365);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 347, -end => 366);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 349, -end => 368);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 349, -end => 368);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 352, -end => 371);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 458, -end => 477);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 460, -end => 479);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 461, -end => 480);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 464, -end => 483);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 466, -end => 485);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 505, -end => 524);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 509, -end => 528);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 513, -end => 532);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 515, -end => 534);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 517, -end => 536);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 529, -end => 548);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 549, -end => 568);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 553, -end => 572);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 555, -end => 574);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 558, -end => 577);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 573, -end => 592);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 574, -end => 593);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 576, -end => 595);
    Bio::Map::Position->new(-element => $pred4, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 578, -end => 597);
    my $pred5 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 88, -end => 107);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 104, -end => 114);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 107, -end => 117);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 107, -end => 117);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 109, -end => 119);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 110, -end => 120);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 115, -end => 134);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 118, -end => 137);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 118, -end => 137);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 120, -end => 139);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 121, -end => 140);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 121, -end => 140);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 137, -end => 156);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 140, -end => 159);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 140, -end => 159);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 142, -end => 161);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 144, -end => 163);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 153, -end => 172);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 165, -end => 184);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 167, -end => 186);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 168, -end => 187);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 170, -end => 189);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 172, -end => 191);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 281, -end => 300);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 284, -end => 303);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 286, -end => 305);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 289, -end => 308);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 289, -end => 308);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 291, -end => 310);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 311, -end => 321);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 455, -end => 474);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 460, -end => 479);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 744, -end => 763);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 780, -end => 799);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 804, -end => 823);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 804, -end => 823);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 805, -end => 824);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 805, -end => 824);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 806, -end => 825);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 824, -end => 834);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 824, -end => 834);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 825, -end => 835);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 825, -end => 835);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 842, -end => 861);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 846, -end => 865);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 847, -end => 866);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 847, -end => 866);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 885, -end => 904);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 890, -end => 909);
    Bio::Map::Position->new(-element => $pred5, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 912, -end => 931);
    my $pred6 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 73, -end => 92);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 76, -end => 95);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 77, -end => 96);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 78, -end => 97);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 78, -end => 97);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 118, -end => 137);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 134, -end => 153);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 137, -end => 156);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 137, -end => 156);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 137, -end => 156);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 138, -end => 157);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 139, -end => 158);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 141, -end => 160);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 144, -end => 163);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 144, -end => 163);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 145, -end => 164);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 157, -end => 176);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 164, -end => 183);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 164, -end => 183);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 165, -end => 184);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 195, -end => 214);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 197, -end => 216);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 197, -end => 216);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 200, -end => 219);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 202, -end => 221);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 300, -end => 319);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 302, -end => 321);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 304, -end => 323);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 305, -end => 324);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 307, -end => 326);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 333, -end => 352);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 334, -end => 353);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 336, -end => 355);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 336, -end => 355);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 337, -end => 356);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 362, -end => 381);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 363, -end => 382);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 365, -end => 384);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 365, -end => 384);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 368, -end => 387);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 475, -end => 494);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 476, -end => 495);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 477, -end => 496);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 478, -end => 497);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 479, -end => 498);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 707, -end => 726);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 708, -end => 727);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 709, -end => 728);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 709, -end => 728);
    Bio::Map::Position->new(-element => $pred6, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 709, -end => 728);
    my $pred7 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 91, -end => 110);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 99, -end => 118);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 100, -end => 119);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 104, -end => 123);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 117, -end => 136);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 118, -end => 137);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 121, -end => 140);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 121, -end => 140);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 123, -end => 142);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 124, -end => 143);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 164, -end => 183);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 166, -end => 185);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 167, -end => 186);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 169, -end => 188);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 171, -end => 190);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 186, -end => 205);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 188, -end => 207);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 188, -end => 207);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 191, -end => 210);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 193, -end => 212);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 216, -end => 235);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 218, -end => 237);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 223, -end => 242);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 224, -end => 243);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 227, -end => 246);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 285, -end => 304);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 287, -end => 306);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 292, -end => 311);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 293, -end => 312);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 296, -end => 315);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 296, -end => 315);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 300, -end => 319);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 300, -end => 319);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 305, -end => 324);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 308, -end => 327);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 317, -end => 336);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 320, -end => 339);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 321, -end => 340);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 321, -end => 340);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 322, -end => 341);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 378, -end => 397);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 381, -end => 400);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 384, -end => 403);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 386, -end => 405);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 390, -end => 409);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Skud"), -start => 422, -end => 441);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Smik"), -start => 423, -end => 442);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Scer"), -start => 423, -end => 442);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Spar"), -start => 425, -end => 444);
    Bio::Map::Position->new(-element => $pred7, -map => Bio::Map::GeneMap->get(-gene => "HEM1", -species => "Sbay"), -start => 427, -end => 446);
    my $pred8 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 53, -end => 72);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 62, -end => 81);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 92, -end => 111);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 109, -end => 123);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 118, -end => 132);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 209, -end => 223);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 243, -end => 259);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 247, -end => 263);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 253, -end => 269);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 255, -end => 271);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 325, -end => 344);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 378, -end => 397);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 463, -end => 482);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 465, -end => 484);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 510, -end => 529);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 550, -end => 566);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 642, -end => 660);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 643, -end => 661);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 645, -end => 663);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 650, -end => 668);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 672, -end => 686);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 673, -end => 687);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 676, -end => 690);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 679, -end => 693);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 714, -end => 732);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 714, -end => 732);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 731, -end => 749);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 751, -end => 767);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 756, -end => 772);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 761, -end => 780);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 763, -end => 782);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 769, -end => 788);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 773, -end => 792);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 784, -end => 800);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 814, -end => 833);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 816, -end => 835);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 820, -end => 839);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 824, -end => 843);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 856, -end => 870);
    Bio::Map::Position->new(-element => $pred8, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 912, -end => 930);
    my $pred9 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 10, -end => 28);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 10, -end => 28);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 24, -end => 42);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 49, -end => 63);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 58, -end => 72);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 67, -end => 85);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 96, -end => 110);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 135, -end => 149);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 136, -end => 150);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 140, -end => 154);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 141, -end => 155);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 166, -end => 184);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 167, -end => 185);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 172, -end => 190);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 174, -end => 192);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 188, -end => 202);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 189, -end => 203);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 192, -end => 206);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 196, -end => 210);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 232, -end => 246);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 234, -end => 248);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 234, -end => 248);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 260, -end => 274);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 273, -end => 287);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 273, -end => 287);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 282, -end => 296);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 493, -end => 511);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 540, -end => 554);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 682, -end => 696);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 689, -end => 707);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 691, -end => 709);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 739, -end => 757);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 742, -end => 760);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 742, -end => 756);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 743, -end => 761);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 745, -end => 759);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 774, -end => 788);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 814, -end => 832);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 824, -end => 842);
    Bio::Map::Position->new(-element => $pred9, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 914, -end => 928);
    my $pred10 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 15, -end => 34);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 51, -end => 70);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 56, -end => 75);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 61, -end => 80);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 65, -end => 84);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 76, -end => 95);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 118, -end => 137);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 138, -end => 157);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 143, -end => 162);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 149, -end => 168);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 151, -end => 170);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 165, -end => 184);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 185, -end => 204);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 199, -end => 218);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 208, -end => 227);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 209, -end => 228);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 211, -end => 230);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 212, -end => 231);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 225, -end => 244);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 387, -end => 406);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 389, -end => 408);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 509, -end => 528);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 533, -end => 552);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 545, -end => 564);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 553, -end => 572);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 557, -end => 576);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 559, -end => 578);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 562, -end => 581);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 577, -end => 596);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 578, -end => 597);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 580, -end => 599);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 582, -end => 601);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 610, -end => 629);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 675, -end => 694);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 676, -end => 695);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 694, -end => 713);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 747, -end => 766);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 750, -end => 769);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 779, -end => 798);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 808, -end => 827);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 824, -end => 843);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 831, -end => 850);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 834, -end => 853);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 859, -end => 878);
    Bio::Map::Position->new(-element => $pred10, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 932, -end => 951);
    my $pred11 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 40, -end => 59);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 86, -end => 105);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 100, -end => 119);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 129, -end => 148);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 132, -end => 151);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 218, -end => 231);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 265, -end => 284);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 306, -end => 319);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 366, -end => 385);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 370, -end => 389);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 422, -end => 435);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 425, -end => 444);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 428, -end => 447);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 501, -end => 514);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 508, -end => 521);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 526, -end => 539);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 557, -end => 570);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 697, -end => 710);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 699, -end => 712);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 732, -end => 745);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 747, -end => 760);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 767, -end => 786);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 780, -end => 799);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 781, -end => 800);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 781, -end => 800);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 781, -end => 800);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 782, -end => 801);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 802, -end => 821);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 802, -end => 821);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 802, -end => 821);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 803, -end => 822);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 803, -end => 822);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 804, -end => 823);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 812, -end => 831);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 822, -end => 841);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 822, -end => 841);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 823, -end => 842);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 823, -end => 842);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 842, -end => 855);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 846, -end => 859);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 847, -end => 860);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 847, -end => 860);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 863, -end => 876);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 910, -end => 923);
    Bio::Map::Position->new(-element => $pred11, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 911, -end => 924);
    my $pred12 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 10, -end => 29);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 13, -end => 32);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 107, -end => 126);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 126, -end => 145);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 130, -end => 149);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 133, -end => 152);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 133, -end => 152);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 134, -end => 153);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 149, -end => 168);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 156, -end => 175);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 156, -end => 175);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 157, -end => 176);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 306, -end => 325);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 316, -end => 335);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 367, -end => 386);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 491, -end => 510);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 521, -end => 540);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 522, -end => 541);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 522, -end => 541);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 524, -end => 543);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 525, -end => 544);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 533, -end => 552);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 549, -end => 568);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 552, -end => 571);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 552, -end => 571);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 553, -end => 572);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 554, -end => 573);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 605, -end => 624);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 720, -end => 739);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 721, -end => 740);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 722, -end => 741);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 722, -end => 741);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 722, -end => 741);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 730, -end => 749);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 733, -end => 752);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 761, -end => 780);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 763, -end => 782);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 767, -end => 786);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 808, -end => 827);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 811, -end => 830);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 838, -end => 857);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 868, -end => 887);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 893, -end => 912);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 923, -end => 942);
    Bio::Map::Position->new(-element => $pred12, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 932, -end => 951);
    my $pred13 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 40, -end => 59);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 59, -end => 78);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 97, -end => 116);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 133, -end => 152);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 143, -end => 162);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 149, -end => 168);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 169, -end => 188);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 212, -end => 230);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 222, -end => 241);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 224, -end => 243);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 229, -end => 248);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 230, -end => 249);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 233, -end => 252);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 262, -end => 281);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 262, -end => 281);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 275, -end => 294);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 277, -end => 296);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 282, -end => 301);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 283, -end => 302);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 286, -end => 305);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 331, -end => 349);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 332, -end => 351);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 334, -end => 352);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 336, -end => 355);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 338, -end => 357);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 339, -end => 357);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 339, -end => 357);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 343, -end => 361);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 360, -end => 379);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 363, -end => 382);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 364, -end => 383);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 369, -end => 388);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 371, -end => 390);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 385, -end => 404);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 388, -end => 407);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 391, -end => 410);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 393, -end => 412);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 397, -end => 416);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 427, -end => 445);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 430, -end => 448);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Skud"), -start => 491, -end => 510);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Scer"), -start => 533, -end => 552);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 598, -end => 617);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Sbay"), -start => 717, -end => 735);
    Bio::Map::Position->new(-element => $pred13, -map => Bio::Map::GeneMap->get(-gene => "HEM12", -species => "Spar"), -start => 862, -end => 881);
    my $pred14 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 91, -end => 110);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 105, -end => 122);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 109, -end => 126);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 112, -end => 129);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 113, -end => 130);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 136, -end => 155);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 137, -end => 156);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 141, -end => 160);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 142, -end => 161);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 164, -end => 183);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 168, -end => 187);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 172, -end => 191);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 174, -end => 193);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 175, -end => 189);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 176, -end => 190);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 181, -end => 195);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 183, -end => 197);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 233, -end => 252);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 235, -end => 254);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 235, -end => 254);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 240, -end => 259);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 262, -end => 281);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 264, -end => 283);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 264, -end => 283);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 312, -end => 326);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 315, -end => 329);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 318, -end => 332);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 319, -end => 333);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 636, -end => 653);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 637, -end => 654);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 639, -end => 656);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 644, -end => 661);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 706, -end => 725);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 707, -end => 726);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 713, -end => 732);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 717, -end => 736);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 733, -end => 752);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 735, -end => 754);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 741, -end => 760);
    Bio::Map::Position->new(-element => $pred14, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 745, -end => 764);
    my $pred15 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 89, -end => 108);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 92, -end => 111);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 95, -end => 114);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 96, -end => 115);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 115, -end => 134);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 120, -end => 139);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 126, -end => 145);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 128, -end => 147);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 142, -end => 161);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 193, -end => 212);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 202, -end => 221);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 203, -end => 222);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 205, -end => 224);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 219, -end => 238);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 260, -end => 279);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 264, -end => 283);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 270, -end => 289);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 335, -end => 349);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 337, -end => 351);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 338, -end => 352);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 343, -end => 357);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 503, -end => 522);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 512, -end => 526);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 526, -end => 540);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 528, -end => 542);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 530, -end => 544);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 539, -end => 553);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 546, -end => 560);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 547, -end => 566);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 551, -end => 570);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 553, -end => 572);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 556, -end => 575);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 560, -end => 574);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 562, -end => 576);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 568, -end => 582);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 573, -end => 587);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 671, -end => 685);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 672, -end => 686);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 675, -end => 689);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 678, -end => 692);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 782, -end => 801);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 806, -end => 825);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 808, -end => 827);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 812, -end => 831);
    Bio::Map::Position->new(-element => $pred15, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 816, -end => 835);
    my $pred16 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 141, -end => 160);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 145, -end => 164);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 150, -end => 169);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 151, -end => 170);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 215, -end => 234);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 232, -end => 248);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 232, -end => 248);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 234, -end => 250);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 238, -end => 254);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 265, -end => 281);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 338, -end => 354);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 342, -end => 358);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 447, -end => 466);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 451, -end => 470);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 521, -end => 540);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 539, -end => 558);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 545, -end => 564);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 547, -end => 566);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 547, -end => 566);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 623, -end => 639);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 644, -end => 660);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 713, -end => 732);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 714, -end => 733);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 720, -end => 739);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 724, -end => 743);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 779, -end => 795);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 780, -end => 796);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 780, -end => 796);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 780, -end => 796);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 781, -end => 797);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 805, -end => 824);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 805, -end => 824);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 806, -end => 825);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 806, -end => 825);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 807, -end => 826);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 830, -end => 849);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 862, -end => 878);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 863, -end => 879);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 866, -end => 882);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 870, -end => 886);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 881, -end => 900);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 884, -end => 903);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 884, -end => 903);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 885, -end => 904);
    Bio::Map::Position->new(-element => $pred16, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 886, -end => 905);
    my $pred17 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 107, -end => 126);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 120, -end => 139);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 124, -end => 143);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 126, -end => 145);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 127, -end => 146);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 128, -end => 147);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 131, -end => 150);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 133, -end => 152);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 133, -end => 152);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 134, -end => 153);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 150, -end => 169);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 157, -end => 176);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 157, -end => 176);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 158, -end => 177);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 172, -end => 191);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 187, -end => 206);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 194, -end => 213);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 194, -end => 213);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 195, -end => 214);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 231, -end => 250);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 232, -end => 251);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 232, -end => 251);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 233, -end => 252);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 236, -end => 255);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 243, -end => 262);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 249, -end => 268);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 251, -end => 270);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 252, -end => 271);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 266, -end => 280);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 272, -end => 286);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 275, -end => 289);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 276, -end => 290);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 277, -end => 291);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 307, -end => 321);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 310, -end => 324);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 313, -end => 327);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 314, -end => 328);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 713, -end => 732);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 714, -end => 733);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 720, -end => 739);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 724, -end => 743);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 816, -end => 835);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 818, -end => 837);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 822, -end => 841);
    Bio::Map::Position->new(-element => $pred17, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 826, -end => 845);
    my $pred18 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 92, -end => 111);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 100, -end => 119);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 101, -end => 120);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 105, -end => 124);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 118, -end => 137);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 215, -end => 234);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 216, -end => 235);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 218, -end => 237);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 219, -end => 238);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 223, -end => 242);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 224, -end => 243);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 227, -end => 246);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 240, -end => 259);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 242, -end => 261);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 247, -end => 266);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 248, -end => 267);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 251, -end => 270);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 281, -end => 300);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 281, -end => 300);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 283, -end => 302);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 289, -end => 308);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 313, -end => 332);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 316, -end => 335);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 321, -end => 340);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 321, -end => 340);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 325, -end => 344);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 328, -end => 347);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 331, -end => 350);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 332, -end => 351);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 337, -end => 356);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 339, -end => 358);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 376, -end => 395);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 408, -end => 427);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 692, -end => 711);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 693, -end => 712);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 698, -end => 717);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 699, -end => 718);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 727, -end => 746);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 729, -end => 748);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 735, -end => 754);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 739, -end => 758);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Scer"), -start => 757, -end => 776);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Sbay"), -start => 759, -end => 778);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Spar"), -start => 765, -end => 784);
    Bio::Map::Position->new(-element => $pred18, -map => Bio::Map::GeneMap->get(-gene => "HEM13", -species => "Smik"), -start => 769, -end => 788);
    my $pred19 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 29, -end => 43);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 34, -end => 48);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 36, -end => 50);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 45, -end => 59);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 52, -end => 66);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 57, -end => 71);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 62, -end => 76);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 66, -end => 80);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 77, -end => 91);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 97, -end => 116);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 135, -end => 154);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 137, -end => 151);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 138, -end => 152);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 140, -end => 159);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 142, -end => 156);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 143, -end => 157);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 146, -end => 165);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 148, -end => 167);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 162, -end => 181);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 182, -end => 201);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 183, -end => 202);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 186, -end => 205);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 190, -end => 209);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 199, -end => 218);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 208, -end => 227);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 209, -end => 228);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 211, -end => 230);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 225, -end => 244);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 239, -end => 258);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 243, -end => 262);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 245, -end => 264);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 245, -end => 264);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 272, -end => 291);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 274, -end => 293);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 294, -end => 313);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 504, -end => 523);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 527, -end => 541);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 529, -end => 543);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 548, -end => 567);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 552, -end => 571);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 554, -end => 573);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 557, -end => 576);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 571, -end => 585);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 572, -end => 586);
    Bio::Map::Position->new(-element => $pred19, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 574, -end => 588);
    my $pred20 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 98, -end => 117);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 99, -end => 118);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 100, -end => 119);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 100, -end => 119);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 103, -end => 122);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 135, -end => 154);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 136, -end => 155);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 140, -end => 159);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 141, -end => 160);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 168, -end => 178);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 171, -end => 185);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 172, -end => 186);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 172, -end => 186);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 173, -end => 187);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 173, -end => 187);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 174, -end => 188);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 174, -end => 188);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 177, -end => 191);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 179, -end => 193);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 184, -end => 194);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 188, -end => 198);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 189, -end => 199);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 192, -end => 202);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 196, -end => 206);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 200, -end => 219);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 201, -end => 220);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 201, -end => 220);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 201, -end => 220);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 201, -end => 211);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 202, -end => 221);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 202, -end => 212);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 236, -end => 255);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 238, -end => 257);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 238, -end => 257);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 243, -end => 262);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 272, -end => 282);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 280, -end => 290);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 306, -end => 316);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 312, -end => 322);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 437, -end => 447);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 802, -end => 812);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 802, -end => 812);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 803, -end => 813);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 803, -end => 813);
    Bio::Map::Position->new(-element => $pred20, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 804, -end => 814);
    my $pred21 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 32, -end => 45);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 37, -end => 50);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 39, -end => 52);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 48, -end => 61);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 92, -end => 111);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 111, -end => 130);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 117, -end => 136);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 118, -end => 137);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 118, -end => 137);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 119, -end => 138);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 136, -end => 155);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 137, -end => 156);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 138, -end => 157);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 142, -end => 161);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 143, -end => 162);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 143, -end => 162);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 143, -end => 162);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 144, -end => 163);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 185, -end => 204);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 186, -end => 205);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 189, -end => 208);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 193, -end => 212);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 233, -end => 252);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 235, -end => 254);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 235, -end => 254);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 240, -end => 259);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 285, -end => 304);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 286, -end => 305);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 286, -end => 305);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 290, -end => 309);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 292, -end => 311);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 294, -end => 313);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 295, -end => 314);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 296, -end => 315);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 297, -end => 316);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 513, -end => 532);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 514, -end => 533);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 514, -end => 533);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 516, -end => 535);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 517, -end => 536);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 734, -end => 747);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 735, -end => 748);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 736, -end => 749);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 736, -end => 749);
    Bio::Map::Position->new(-element => $pred21, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 736, -end => 749);
    my $pred22 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 33, -end => 52);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 38, -end => 57);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 40, -end => 59);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 49, -end => 68);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 90, -end => 109);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 98, -end => 117);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 99, -end => 118);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 103, -end => 122);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 116, -end => 135);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 133, -end => 152);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 134, -end => 153);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 138, -end => 157);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 139, -end => 158);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 165, -end => 184);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 166, -end => 185);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 171, -end => 190);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 173, -end => 192);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 188, -end => 207);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 189, -end => 208);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 192, -end => 211);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 196, -end => 215);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 219, -end => 238);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 221, -end => 240);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Sbay"), -start => 222, -end => 240);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 226, -end => 245);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 227, -end => 246);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 230, -end => 249);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Smik"), -start => 253, -end => 271);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Scer"), -start => 255, -end => 273);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM14", -species => "Spar"), -start => 255, -end => 273);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 288, -end => 306);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 290, -end => 308);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 295, -end => 313);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 296, -end => 314);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 299, -end => 317);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 339, -end => 358);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 342, -end => 361);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 343, -end => 362);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 348, -end => 367);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 350, -end => 369);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 380, -end => 399);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 383, -end => 402);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 386, -end => 405);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 388, -end => 407);
    Bio::Map::Position->new(-element => $pred22, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 392, -end => 411);
    my $pred23 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 5, -end => 24);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 24, -end => 43);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 46, -end => 65);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 69, -end => 88);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 71, -end => 90);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 71, -end => 90);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 99, -end => 118);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 146, -end => 160);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 151, -end => 165);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 157, -end => 171);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 159, -end => 173);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 173, -end => 187);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 253, -end => 272);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 256, -end => 275);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 257, -end => 276);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 414, -end => 433);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 475, -end => 494);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 478, -end => 497);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 493, -end => 512);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 493, -end => 512);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 507, -end => 526);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 509, -end => 528);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 510, -end => 529);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 519, -end => 538);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 522, -end => 541);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 527, -end => 546);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 529, -end => 548);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 529, -end => 548);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 531, -end => 550);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 573, -end => 592);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 574, -end => 593);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 576, -end => 595);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 578, -end => 597);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 634, -end => 653);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 645, -end => 664);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 649, -end => 668);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 654, -end => 673);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 659, -end => 678);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 666, -end => 685);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 697, -end => 716);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 779, -end => 793);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 780, -end => 794);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 780, -end => 794);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 780, -end => 794);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 781, -end => 795);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 802, -end => 821);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 802, -end => 821);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 803, -end => 822);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 803, -end => 822);
    Bio::Map::Position->new(-element => $pred23, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 804, -end => 823);
    my $pred24 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 104, -end => 123);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 104, -end => 122);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 109, -end => 128);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 115, -end => 134);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 117, -end => 136);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 123, -end => 141);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 130, -end => 148);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 130, -end => 148);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 131, -end => 149);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 131, -end => 150);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 143, -end => 162);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 155, -end => 169);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 162, -end => 181);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 169, -end => 188);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 169, -end => 188);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 170, -end => 189);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 485, -end => 499);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 485, -end => 503);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 489, -end => 507);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 503, -end => 521);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 503, -end => 521);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 505, -end => 523);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 506, -end => 524);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 515, -end => 533);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 529, -end => 543);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 534, -end => 548);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 536, -end => 550);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 538, -end => 552);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 544, -end => 563);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 554, -end => 572);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 556, -end => 574);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 558, -end => 577);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 559, -end => 577);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 560, -end => 579);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 566, -end => 585);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 571, -end => 590);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 579, -end => 593);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 580, -end => 594);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 583, -end => 597);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 589, -end => 603);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 634, -end => 653);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 635, -end => 654);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 636, -end => 655);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 637, -end => 656);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 637, -end => 656);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 731, -end => 749);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 732, -end => 750);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 733, -end => 751);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 733, -end => 751);
    Bio::Map::Position->new(-element => $pred24, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 733, -end => 751);
    my $pred25 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 152, -end => 171);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 157, -end => 176);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 163, -end => 182);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 165, -end => 184);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 179, -end => 198);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 197, -end => 216);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 206, -end => 225);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 207, -end => 226);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 209, -end => 228);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 217, -end => 231);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 219, -end => 233);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 223, -end => 242);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 224, -end => 238);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 225, -end => 239);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 228, -end => 242);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 247, -end => 266);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 249, -end => 268);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 254, -end => 273);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 255, -end => 274);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 258, -end => 277);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 260, -end => 274);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 277, -end => 296);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 279, -end => 298);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 284, -end => 303);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 285, -end => 304);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 288, -end => 307);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 380, -end => 399);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 383, -end => 402);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 386, -end => 405);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 388, -end => 407);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 392, -end => 411);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 476, -end => 495);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 520, -end => 539);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 525, -end => 544);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 525, -end => 544);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 527, -end => 546);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 529, -end => 548);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 529, -end => 548);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 536, -end => 555);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 540, -end => 559);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 542, -end => 561);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 549, -end => 568);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 573, -end => 592);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 574, -end => 593);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 576, -end => 595);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Smik"), -start => 578, -end => 597);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Scer"), -start => 599, -end => 613);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Spar"), -start => 601, -end => 615);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Skud"), -start => 607, -end => 621);
    Bio::Map::Position->new(-element => $pred25, -map => Bio::Map::GeneMap->get(-gene => "HEM15", -species => "Sbay"), -start => 612, -end => 626);
    my $pred26 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 10, -end => 29);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 42, -end => 61);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 56, -end => 72);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 57, -end => 73);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 57, -end => 73);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 58, -end => 74);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 59, -end => 75);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 108, -end => 127);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 127, -end => 146);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 130, -end => 149);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 134, -end => 153);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 134, -end => 153);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 135, -end => 154);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 149, -end => 168);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 156, -end => 175);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 156, -end => 175);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 157, -end => 176);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 159, -end => 178);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 220, -end => 239);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 290, -end => 309);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 321, -end => 340);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 527, -end => 546);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 581, -end => 600);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 607, -end => 626);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 610, -end => 629);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 611, -end => 630);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 612, -end => 631);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 673, -end => 692);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 721, -end => 740);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 722, -end => 741);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 723, -end => 742);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 723, -end => 742);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 723, -end => 742);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 747, -end => 766);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 748, -end => 767);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 769, -end => 788);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 770, -end => 789);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 770, -end => 789);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 771, -end => 790);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 785, -end => 804);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 804, -end => 823);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 804, -end => 823);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 805, -end => 824);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 805, -end => 824);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 806, -end => 825);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 885, -end => 901);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 888, -end => 904);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 888, -end => 904);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 889, -end => 905);
    Bio::Map::Position->new(-element => $pred26, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 890, -end => 906);
    my $pred27 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 89, -end => 108);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 97, -end => 116);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 98, -end => 117);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 102, -end => 121);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 104, -end => 123);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 111, -end => 130);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 115, -end => 134);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 168, -end => 187);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 184, -end => 203);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 254, -end => 273);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 256, -end => 275);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 261, -end => 280);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 262, -end => 281);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 265, -end => 284);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 272, -end => 291);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 280, -end => 299);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 281, -end => 300);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 283, -end => 302);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 288, -end => 307);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 289, -end => 308);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 292, -end => 311);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 380, -end => 399);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 383, -end => 402);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 383, -end => 402);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 386, -end => 405);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 386, -end => 405);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 387, -end => 406);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 388, -end => 407);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 392, -end => 411);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 392, -end => 411);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 394, -end => 413);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 437, -end => 456);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 543, -end => 562);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 560, -end => 579);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 660, -end => 679);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 661, -end => 680);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 680, -end => 699);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 778, -end => 797);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 779, -end => 798);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 779, -end => 798);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 779, -end => 798);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 780, -end => 799);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 808, -end => 827);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Sbay"), -start => 808, -end => 827);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 809, -end => 828);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 809, -end => 828);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Smik"), -start => 810, -end => 829);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Skud"), -start => 840, -end => 859);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Scer"), -start => 932, -end => 951);
    Bio::Map::Position->new(-element => $pred27, -map => Bio::Map::GeneMap->get(-gene => "HEM2", -species => "Spar"), -start => 932, -end => 951);
    my $pred28 = Bio::Map::Prediction->new(-source => "meme");
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 203, -end => 222);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 205, -end => 224);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 210, -end => 229);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 211, -end => 230);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 214, -end => 233);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 228, -end => 247);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 230, -end => 249);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 235, -end => 254);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 236, -end => 255);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 239, -end => 258);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 272, -end => 291);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 274, -end => 293);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 279, -end => 298);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 280, -end => 299);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 283, -end => 302);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 305, -end => 324);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 308, -end => 327);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 313, -end => 332);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 313, -end => 332);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 317, -end => 336);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Sbay"), -start => 379, -end => 398);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Smik"), -start => 382, -end => 401);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Skud"), -start => 385, -end => 404);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Spar"), -start => 387, -end => 406);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM4", -species => "Scer"), -start => 391, -end => 410);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 466, -end => 485);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 467, -end => 486);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 468, -end => 487);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 469, -end => 488);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 470, -end => 489);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 515, -end => 534);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 516, -end => 535);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 516, -end => 535);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 518, -end => 537);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 519, -end => 538);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 546, -end => 565);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 549, -end => 568);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 549, -end => 568);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 550, -end => 569);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 551, -end => 570);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 712, -end => 731);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 713, -end => 732);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 713, -end => 732);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 720, -end => 739);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Smik"), -start => 724, -end => 743);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 725, -end => 744);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Spar"), -start => 726, -end => 745);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Scer"), -start => 726, -end => 745);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Skud"), -start => 726, -end => 745);
    Bio::Map::Position->new(-element => $pred28, -map => Bio::Map::GeneMap->get(-gene => "HEM3", -species => "Sbay"), -start => 739, -end => 758);
    push(@predicitons_of_interest, $pred1, $pred2, $pred3, $pred4, $pred5, $pred6, $pred7, $pred8, $pred9, $pred10,
                                   $pred11, $pred12, $pred13, $pred14, $pred15, $pred16, $pred17, $pred18, $pred19, $pred20,
                                   $pred21, $pred22, $pred23, $pred24, $pred25, $pred26, $pred27, $pred28);
    
    # our run class has placed the predictions on our maps, lets find all the
    # intersections of predictions shared by all HEM1-HEMx
    my $rel = Bio::Map::GeneRelative->new(-gene => 0);
    my $min_pables_num = 2;
    my $di_able = Bio::Map::Mappable->disconnected_intersections(\@predicitons_of_interest,
                                                                   -min_mappables_num => $min_pables_num,
                                                                   -relative => $rel,
                                                                   -min_overlap_percent => 66);
    
    my @positions = $di_able->get_positions;
    #print "--\nAsked for disconnected intersections of all predictions that agree for at least $min_pables_num meme runs:\n"; 
    #foreach my $pos (sort { $a->start($rel) <=> $b->start($rel) } @positions) {
    #    #*** disconnected_intersections should (?) somehow manage to transfer across the appropriate element to each generated position
    #    print "  pos ", $pos->toString($rel), " on map for gene ", $pos->map->gene->universal_name, " and species ", $pos->map->species, "\n";
    #}
    #print "__\n";
    is @positions, 17; # 223??
    
    $min_pables_num = 3;
    my @overlapping_groups = Bio::Map::Mappable->overlapping_groups(\@predicitons_of_interest,
                                                                   -min_mappables_num => $min_pables_num,
                                                                   -relative => $rel,
                                                                   -min_overlap_percent => 66);
    
    #print "\nAsked for all the overlapping groups of predictions amongst all predictions, but each group needed predictions from at least $min_pables_num meme runs:\n";
    #foreach my $group (@overlapping_groups) {
    #    print "  - got a group with ", scalar(@{$group}), " positions\n";
    #    foreach my $pos (@{$group}) {
    #        print "    * pos ", $pos->toString($rel), " on map for gene ", $pos->map->gene->universal_name, " and species ", $pos->map->species, " from meme run '", $pos->element->name, "'\n";
    #    }
    #}
    #print "__\n";
    
    # not sure what the correct answer is
}
