#-*-Perl-*-
## Bioperl Test Harness Script for Modules

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

use strict;
use ExtUtils::MakeMaker;
use Bio::Root::IO;
use constant TEST_COUNT => 11;
use constant TODO_LIST	=> "[1,2,3,4,5,6,7,8,9,10,11]";

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
use Bio::Graphics::SimpleConfigurator;

##
## create a new configurator
##

### NOTE: There is not yet an implementation of Renderer, so this example won't
###       work. Change the following line to match the implementation you
###       want to test. Also replace the "use" line above.
#
# my $conf	= new Bio::Graphics::SimpleConfigurator;
#

ok($conf);

##
## exercise the setter
##
ok( $conf->set('EST','fgcolor','chartreuse'), 'chartreuse' );
ok( $conf->set('fgcolor','blue'), 'blue' );
ok( $conf->set('EST','height','sub { return \"Short\"; }'), 'sub { return \"Short\"; }');

##
## exercise the getters after setting
##

# exercise get_sections
my @sections = $conf->get_sections();
ok( $sections[0], 'EST' );

# exercise get_tags
my @tags	= $conf->get_tags();
ok( $tags[0], 'fgcolor' );

@tags	= $conf->get_tags('EST');
@tags	= sort(@tags);
$tags	= join(" ", @tags);
ok( $tags eq 'fgcolor height' );

# exercise get
my @value	= $conf->get('fgcolor');
ok( $value[0], 'blue' );
@value	= $conf->get('EST','fgcolor');
ok( $value[0], 'chartreuse' );
@value	= $conf->get('EST','height');
ok( $value[0], 'sub { return \"Short\"; }');

# exercise get_and_eval
@value	= $conf->get('EST','height');
ok( $value[0], 'Short' );

