#!/usr/bin/perl

package SimpleStateMachine;

use strict;
use Bio::Root::Root;
use Bio::Tools::StateMachine::AbstractStateMachine qw($INITIAL_STATE $FINAL_STATE);
use vars qw( @ISA );

@ISA = qw( Bio::Root::Root
           Bio::Tools::StateMachine::AbstractStateMachine );

my $data = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
my @state_transitions = (  [ $INITIAL_STATE, 'State1'],
                           [ 'State1', 'State2' ],
                           [ 'State2', 'State3' ],
                           [ 'State3', 'State4' ],
                           [ 'State4', $FINAL_STATE]
                        );
sub new {
    my($caller,@args) = @_;
    my $self = $caller->SUPER::new( @args);
    print STDERR "test machine: initial = $INITIAL_STATE, final = $FINAL_STATE\n";

    $self->_init_state_machine( -transition_table => \@state_transitions );
    return $self;
}


sub check_for_new_state {
    my ($self) = @_;

    my $chunk = $self->_next_data_chunk();
    $self->verbose and print STDERR "Checking for new state ($chunk)\n";

    my $newstate = undef;
    if( $chunk eq 'A' ) {
        $newstate = 'State1';
    }
    elsif ($chunk eq 'G' ) {
        $newstate = 'State2';
    }
    elsif ($chunk eq 'M' ) {
        $newstate = 'State3';
    }
    elsif ($chunk eq 'T' ) {
        $newstate = 'State4';
    }
    elsif ($chunk eq 'Z' ) {
        $newstate = $FINAL_STATE;
    }
    return $newstate;
}

sub _next_data_chunk {
    my $self = shift;
    $self->verbose and print STDERR "Getting next input chunk...\n";
    if( ! defined $self->{'unread_data'}) {
        $self->{'unread_data'} = reverse $self->{'_data'};
    }
    chop $self->{'unread_data'};
}

sub process_data {
    my ($self, $data) = @_;
    $self->{'_data'} = $data;
    $self->run();
}

#---------------------------------
package main;

use Error qw(:try);

try {

    my $sm = SimpleStateMachine->new( -verbose => 1 );

    print "\nRunning state machine...\n";

    $sm->process_data( "ABCDEFGHIJKLMNOPQRSTUVWXYZ" );
}
catch Bio::Root::Exception with {
    my $err = shift;
    print "\nCaught exception:\n\n$err\n";
};

