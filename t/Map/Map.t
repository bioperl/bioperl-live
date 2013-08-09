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
        is ${$groups[1]}[0], $pos7;
        is ${$groups[0]}[0], $pos6;
        is ${$groups[0]}[1], $pos5;
        is ${$groups[0]}[2]->toString($gene_rel), $pos4->toString($gene_rel);
        is ${$groups[0]}[3]->toString($gene_rel), $pos3->toString($gene_rel);
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
                  -requires_modules => [qw(Bio::Tools::Run::Ensembl Bio::EnsEMBL::Registry XML::Twig)],
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
        ok defined($gene->display_xref($map4));
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
        ok substr($map1->subseq($gene->coding_position($map1)), 0, 3);
        #my $exon1_str = 'GGGCTTGTGGCGCGAGCTTCTGAAACTAGGCGGCAGAGGCGGAGCCGCTGTGGCACTGCTGCGCCTCTGCTGCGCCTCGGGTGTCTTTTGCGGCGGTGGGTCGCCGCCGGGAGAAGCGTGAGGGGACAGA';
        my $exon1_pos = $gene->get_exon_position($map1, 1);
        cmp_ok(length($map1->subseq($exon1_pos)), '>', 20);
        is $exon1_pos->seq, $map1->subseq($exon1_pos);
    }

    # test a gene with multiple transcripts...
    #...
}
