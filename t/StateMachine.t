# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
##
## Tests the SearchIO::blast::blast module for parsing traditional
## BLAST reports (non-XML).

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
    plan tests => 2; 
}

if( $error == 1 ) {
    exit(0);
}

use Bio::Tools::StateMachine::IOStateMachine;
use Bio::Tools::StateMachine::AbstractStateMachine;
ok(1);

package TestAbstractStateMachine;

use Bio::Root::RootI;
use Bio::Tools::StateMachine::AbstractStateMachine qw($INITIAL_STATE $FINAL_STATE);
use vars qw( @ISA );

@ISA = qw( Bio::Root::RootI
           Bio::Tools::StateMachine::AbstractStateMachine );

my @state_transitions = (  [ $INITIAL_STATE, 'State1'],
                           [ 'State1', 'State2' ],
                           [ 'State2', $FINAL_STATE]
                        );
sub new {
    my($caller,@args) = @_;
    my $self = $caller->SUPER::new( @args);
    $self->_init_state_machine( -transition_table => \@state_transitions );
    return $self;
}

package main;

my $sm = TestAbstractStateMachine->new();

ok($sm);


