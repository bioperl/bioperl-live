#-*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use ExtUtils::MakeMaker;
use Bio::Root::IO;
use constant TEST_COUNT => 1;
use constant TODO_LIST	=> "[1]";

BEGIN {
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan (	test => TEST_COUNT,
    			todo => TODO_LIST );
}

use lib '..','./blib/lib';
use lib "$ENV{HOME}/cvswork/bioperl-live/";
use Bio::Graphics::SimpleRenderer;
use Bio::SeqFeature::SimpleCollection;
use Bio::Graphics::SimpleConfigurator;
use Bio::SeqFeature::Generic;

##
## create a new renderer
##

### NOTE: There is not yet an implementation of Renderer, so this example won't
###       work. Change the following line to match the implementation you
###       want to test. Also replace the "use" line above.
#
 my $rend	= new Bio::Graphics::SimpleRenderer;
#
ok($rend);

##
## prepare for rendering
##

### NOTE: Before any rendering can happen, we have to know what we need to
###       render, how we want it rendered, and optionally a panel. So! We'll
###       use the Bio::Graphics::Configurator method get to get our preferences
###       ready, we'll use the Bio::SeqFeature::Generic module to create some
###       features, and we'll use the Bio::SeqFeature::Collection method 
###       add_features to get our collection of features to render.

my $feat1	= new Bio::SeqFeature::Generic (	-start		=> 10,
															-end			=> 100,
															-strand		=> -1,
															-primary		=> 'repeat',
															-source_tag	=> 'repeatmasker' );

ok( $feat1 );
														
my $feat2	= new Bio::SeqFeature::Generic (	-start		=> 120,
															-end			=> 150,
															-strand		=> -1,
															-primary		=> 'repeat',
															-source_tag	=> 'repeatmasker' );

ok( $feat2 );

my $collection	= new Bio::SeqFeature::SimpleCollection;
ok( $collection );

$collection->add_features($feat1,$feat2);
ok( $features );

my $prefs = new Bio::Graphics::SimpleConfigurator;
ok( $prefs );

##
## exercise the renderer
##

my $number_of_tracks	= $rend->render();
ok( $number_of_tracks, 0 );

$number_of_tracks	= $rend->render( $features, $prefs );
ok( $number_of_tracks, 1 );

my ( $tracks_rendered, $panel ) = $rend->render( $features, $prefs );

ok( $tracks_rendered );
ok( $panel );

my ( $tracks_rendered_pass_2, $panel_pass_2 )	= $rend->render( $features, $prefs, $panel );

ok( $tracks_rendered_pass_2, $tracks_rendered );
ok( $panel_pass_2, $panel );

