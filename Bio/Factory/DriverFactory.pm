# $Id$
#
# BioPerl module for Bio::Factory::DriverFactory
#
# Cared for by Jason Stajich <jason@chg.mc.duke.edu> and
#              Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Jason Stajich, Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Factory::DriverFactory - Base class for factory classes loading drivers

=head1 SYNOPSIS


=head1 DESCRIPTION

This a base class for factory classes that load drivers. Normally, you don't
instantiate this class directly.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                - General discussion
  http://bio.perl.org/MailList.html    - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email Jason Stajich <jason@chg.mc.duke.edu>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

#'

package Bio::Factory::DriverFactory;
use strict;
use Bio::Root::RootI;
use vars qw(@ISA %DRIVERS);

@ISA = qw(Bio::Root::RootI); 

BEGIN {
    %DRIVERS = ();
}

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    return $self;
}

=head2 register_driver

 Title   : register_driver
 Usage   : $factory->register_driver("genscan", "Bio::Tools::Genscan");
 Function: Registers a driver a factory class should be able to instantiate.
 Example : 
 Returns : 
 Args    : Key of the driver (string) and the module implementing the driver
           (string).

=cut

sub register_driver {
    my ($self, @args) = @_;
    my %drivers = @args;

    foreach my $drv (keys(%drivers)) {
	$DRIVERS{$drv} = $drivers{$drv};
    }
}

=head2 get_driver

 Title   : get_driver
 Usage   : $module = $factory->get_driver("genscan");
 Function: Returns the module implementing a driver registered under the
           given key.
 Example : 
 Returns : A string.
 Args    : Key of the driver (string).

=cut

sub get_driver {
    my ($self, $key) = @_;
    return $DRIVERS{$key} if(exists($DRIVERS{$key}));
    return undef;
}

=head2 _load_module

 Title   : _load_module
 Usage   : $self->_load_module("Bio::Tools::Genscan");
 Function: Loads up (like use) a module at run time on demand.
 Example : 
 Returns : TRUE on success
 Args    :

=cut

sub _load_module {
    my ($self, $name) = @_;
    my ($module, $load, $m);
    
    $module = "_<$name.pm";
    $load = "$name.pm";
    $load =~ s/:/\//g;

    return 1 if $main::{$module};
    eval {
	require $load;
    };
    if ( $@ ) {
	$self->throw("$load: $name cannot be found: ".$@);
    }
    return 1;
}

1;
