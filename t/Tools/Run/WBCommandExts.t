#-*-perl-*-
#$Id$
# testing CommandExts
use strict;
use warnings;
our $home;
BEGIN {
    use Bio::Root::Test;
    use lib '.';
    use lib 't/Tools/Run';
    $home = '../../..'; # set to '.' for Build use, 
                      # '../../..' for debugging from .t file
    unshift @INC, $home;
    test_begin(-tests => 22,
	       -requires_modules => [qw(Bio::Tools::Run::WrapperBase
                                        Bio::Tools::Run::WrapperBase::CommandExts)]);
}

use_ok( 'Dummy::Config' );
use_ok( 'Dummy' );
use_ok('Bio::Tools::Run::WrapperBase');
use_ok('Bio::Tools::Run::WrapperBase::CommandExts');

ok my $fac = Dummy->new( -command => 'goob',
			 -narf => 42.0,
			 -schlurb => 'breb',
			 -freen => 1 ), "make factory";
ok $fac->parameters_changed, "parm changed flag set";
is $fac->program_name, 'flurb', "correct prog name";
ok $fac->is_pseudo, "is pseudo";
is $fac->narf, 42, "correct parm set";
ok !$fac->parameters_changed, "parm flag cleared";
is join(' ',@{$fac->_translate_params}), '--schlurb breb -n 42 -f', "translate opts to command line";

ok $fac->reset_parameters, "parm reset";
ok !$fac->narf, "parm cleared after reset";

is_deeply( [$fac->available_parameters('parameters')], [qw( command narf schlurb )], "avail parms");
is_deeply( [$fac->available_parameters('switches')], ['freen'], "avail switches");
is_deeply( [$fac->available_parameters('commands')], [qw(rpsblast cat goob blorb multiglob)], "avail commands");

ok $fac = Dummy->new( -command => 'multiglob',
		     -g_freen => 1,
		     -b_scroob => 10.5,
		      -trud => 'sklim' ), "make composite cmd factory";

is $fac->trud, 'sklim', "comp cmd parm set";

ok my %facs = $fac->_create_factory_set, "make subfactories";
is $facs{goob}->freen, 1, "subfactory 1 parm correct";
is $facs{blorb}->scroob, 10.5, "subfactory 2 parm correct";

$fac->program_dir('.');
# ok $fac->executables('rpsblast'), "find in program_dir";
ok $fac->executables('cat'), "find in syspath";

1;
