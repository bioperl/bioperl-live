#-----------------------------------------------------------------
# $Id$
#
# BioPerl module Bio::Tools::StateMachine::AbstractStateMachine
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

=head1 NAME

Bio::Tools::StateMachine::AbstractStateMachine - Abstract state machine object

=head1 SYNOPSIS

Here is a portion of an implementation. For the full example, see
examples/state-machine/state-machine.pl in the Bioperl distribution.

    package SimpleStateMachine;

    use Bio::Root::Root;
    use Bio::Tools::StateMachine::AbstractStateMachine qw($INITIAL_STATE 
    							  $FINAL_STATE);
    use vars qw( @ISA );

    @ISA = qw( Bio::Root::Root
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



=head1 DESCRIPTION

B<AbstractStateMachine> provides a generic framework for representing a
state machine. This is not an event-based framework where you register
handlers to be called when certain events occur. Instead, it provides
a set of methods that define the basic logic of an object that has
state behavior, that logic being:

=over 4

=item 1. Check for whether or not a new state has occurred in the external world.

=item 2. If so, change the state of the machine to the new state.

=item 3. Otherwise, keep checking current conditions for a new state.

=item 4. Stop checking for new states if we reach the final state, or if an error occurs.

=back

A B<state> is just a string representing the name of the state.  A
state machine is initialized with a B<state transition table>
consisting of a set of allowable transitions, where each B<transition>
is a two-element array in which the first element is the B<from
state>, and the second element is the B<to state>.  This table permits
the AbstractStateMachine to determine if a requested transition is
valid.

This module is flexible enough to represent both deterministic and
non-deterministic finite automata (DFAs and NFAs), but it is fairly
new and should be considered experimental.

The key methods in AbstractStateMachine that define this logic of
operation are:

=over 4

=item check_for_new_state().

Does whatever checking is necessary to determine if a state transition
should occur (for example, read a line of input from STDIN). If a
transition should occur, a string is returned containing the name of
the new state. Otherwise, it returns C<undef>.

This method B<must be implemented> as AbstractStateMachine does not
define it (and in fact will throw a NotImplemented exception if you
fail to implement it).

=item change_state( C<new_state> )

Causes the machine to change its state to the state specified in the
argument. change_state() allows you to mapping a state transition to a
particular handler method that does whatever processing is needed to
deal with the state transition.

=item run()

This method keeps calling check_for_new_state() and if that method
returns a defined value (the name of the state to change to), it then
calls change_state( $state ), where $state is the value returned by
check_for_new_state().

Before calling check_for_new_state(), the run() method checks the
current state of the machine and exits the loop if the current state
ever becomes $PAUSE_STATE, $ERROR_STATE, or $FINAL_STATE.

=item append_input_cache( C<data> )

Adds data to a buffer for processing at the next state
transition. check_for_new_state() should call
append_input_cache() passing it any data it receives while checking
for a new state that should be processed later.

=item get_input_cache()

Retrieves the data stored by calling
append_input_cache(). change_state() should call get_input_cache() to
obtain the data to be processed for the current transition.

=back

This object defines an abstract class, meaning that some but not all methods 
have been implemented. Subclasses must define the methods not implemented here.
These include at a minimum:

=over 4

=item check_for_new_state()

=item change_state()

A default simple implementation of change_state() is provided, but
subclasses of AbstractStateMachine most likely will want to override
this method to do something useful according to the particular state
change.

=back

If your state machine needs to cache input while processing, you will
also need to provide implementations of these methods (which are no-op
in AbstractStateMachine):

=over 4

=item append_input_cache()

=item get_input_cache()

=item clear_input_cache()

=back

There are some other nuances provided by AbstractStateMachine, such as
the ability to pause() and resume() the running of the machine.


=head1 EXAMPLES

To get a feel for how to use this, have look at scripts in the
examples/state-machine directory of the Bioperl distribution. Also
have a look at Bio::Tools::StateMachine::IOStateMachine which provides
a Bio::Root::IO-based implementation of
AbstractStateMachine. Bio::SearchIO::psiblast subclasses
IOStateMachine.


=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

    bioperl-l@bioperl.org              - General discussion
    http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via email
or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR - Steve Chervitz

sac@bioperl.org

See the L<FEEDBACK> section for where to send bug reports and comments.

=head1 ACKNOWLEDGEMENTS

I would like to acknowledge my colleagues at Affymetrix for useful
feedback.

=head1 COPYRIGHT

Copyright (c) 2001 Steve Chervitz. All Rights Reserved.

=cut

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Tools::StateMachine::AbstractStateMachine;

use strict;
use Bio::Root::Interface;
use Exporter ();

use vars qw( @ISA @EXPORT_OK $INITIAL_STATE $FINAL_STATE $PAUSE_STATE $ERROR_STATE );
@ISA = qw( Bio::Root::Interface  Exporter );
@EXPORT_OK = qw( $INITIAL_STATE $FINAL_STATE $PAUSE_STATE $ERROR_STATE );

@Bio::Tools::StateMachine::StateException::ISA = qw( Bio::Root::Exception );

$INITIAL_STATE = 'Initial';
$FINAL_STATE = 'Final';
$PAUSE_STATE = 'Pause';
$ERROR_STATE = 'Error';

sub _init_state_machine {
    my  ($self, @args ) = @_;
    my ($transition_table) = $self->_rearrange( [qw(TRANSITION_TABLE)], @args);

    $self->verbose and print STDERR "Initializing State Machine...\n";

    if($transition_table) {
        $self->_set_transition_table( $transition_table );
    }

    $self->add_transition( $INITIAL_STATE, $FINAL_STATE );
    $self->_set_current_state( $INITIAL_STATE );
}

sub reset {
    my $self = shift;
    $self->verbose and print STDERR "Resetting state machine\n";
    $self->_set_current_state( $INITIAL_STATE );
}

sub _set_current_state {
    my ($self, $state) = @_;
    if( defined $state) {
	$self->verbose and print STDERR "  setting current state to $state\n";
	$self->{'_current_state'} = $state;
    }
}

sub current_state { shift->{'_current_state'} }

sub initial_state { $INITIAL_STATE }

sub final_state { $FINAL_STATE }

sub pause_state { $PAUSE_STATE }

sub error_state { $ERROR_STATE }

sub resume_state {
    my ($self, $state) = @_;
    if( $state ) {
      $self->{'_resume_state'} = $state;
    }
    $self->{'_resume_state'};
}

sub _clear_resume_state {
    my $self = shift;
    undef $self->{'_resume_state'};
}

=head2 running

The machine is either running or not running. 
Once the machine has stopped running, it cannot be re-started.
Use pause() to temporarily halt a machine without exiting the run state.

=cut 

sub running { shift->{'_running'} }

sub _set_running {
    my $self = shift;
    $self->{'_running'} = shift;
}

sub run {
    my ($self, @args) = @_;

    my $verbose = $self->verbose;
    my $curr_state = $self->current_state;
    $self->_set_running( 1 );

    while( not ($curr_state eq $PAUSE_STATE ||
                $curr_state eq $ERROR_STATE ||
                $curr_state eq $FINAL_STATE )) {

	$verbose and print STDERR "Current state (run): ${\$self->current_state}\n";

        if( my $state = $self->check_for_new_state()) {
            $self->change_state( $state );
        }

        $curr_state = $self->current_state;
    }

    # Handle EOF situations
    if( not ($curr_state eq $PAUSE_STATE ||
             $curr_state eq $FINAL_STATE )) {

        $self->change_state( $FINAL_STATE );
	$self->_set_running( 0 );
    }

    $verbose and print STDERR "StateMachine Run complete ($curr_state).\n";
}

# The pause() and resume() methods don't go through change_state()
sub pause {
    my ($self) = @_;
#    print "PAUSING...\n";
    $self->resume_state( $self->current_state );
    $self->_set_current_state( $PAUSE_STATE );
#    print "After pause(): Current state: ${\$self->current_state}\n";
}

sub paused {
    my ($self) = @_;
    return $self->current_state eq $PAUSE_STATE;
}

sub throw{
   my ($self,@args) = @_;
   $self->_set_current_state( $ERROR_STATE );
   $self->_set_running( 0 );
   $self->SUPER::throw( @args );
}

sub error {
    my ($self, $err) = @_;
    return $self->current_state eq $ERROR_STATE;
}

sub resume {
    my ($self) = @_;

    # Don't resume if we're done.
    return if $self->current_state eq $FINAL_STATE;

#    print "RESUMING...\n";
    $self->_set_current_state( $self->resume_state );
    $self->_clear_resume_state;
    $self->run();
}

=head2 transition_table

 Arg      : n/a
 Returns  : An array of array references to two-element arrays.
            Each array ref defines a single transition where
            the first element is the name of the "from" state and
            the second element is the name of the "to" state.

 Example  : $sm->transition_table( [ $INITIAL_STATE, 'State1'],
				   [ 'State1', 'State2' ],
				   [ 'State2', 'State3' ],
				   [ 'State3', $FINAL_STATE]
				 );

=cut

sub transition_table {
    my ($self) = @_;

    return @{$self->{'_transition_table'}};
}

sub _set_transition_table {
    my ($self, $table_ref) = @_;

    my $verbose = $self->verbose;
    $verbose and print STDERR "Setting state transition table:\n";

    if( not ref($table_ref) eq 'ARRAY') {
	$self->throw( -class => 'Bio::Root::BadParameter',
                      -text => "Can't set state transition table: Arg wasn't an array reference."
                    );
    }

    foreach my $t (@$table_ref) {
        if( ref($t) and scalar(@$t) == 2 ) {
            push @{$self->{'_transition_table'}->{$t->[0]}}, $t->[1];
            $verbose and print STDERR "  adding: $t->[0] -> $t->[1]\n";
        }
        else {
            $self->throw( -class => 'Bio::Root::BadParameter',
                          -text => "Can't add state transition from table: Not a 2-element array reference ($t)"
                        );
        }
    }
}

=head2 add_transition

 Arg      : Two string arguments where:
            First string = name of the "from" state.
            Second string = name of the "to" state.
 Throws   : A Bio::Root::BadParameter exception if two arguments
            are not provided.

=cut

sub add_transition {
    my ($self, $from, $to) = @_;

    if( defined($from) and defined($to) ) {
	push @{$self->{'_transition_table'}->{$from}}, $to;
    }
    else {
	$self->throw( -class => 'Bio::Root::BadParameter',
                      -text => "Can't add state transition: Insufficient arguments."
                    );
    }
}


=head2 change_state

 Purpose  : To cause the machine to change its state.
 Argument : A String containing the name of the the new state.
 Returns  : n/a
 Throws   : A Bio::Tools::StateMachine::StateException exception if the
            state transition cannot be carried out.

This is a default implementation that simply validates the state change
(by calling  validate_transition) and then calls finalize_state_change()
if the transition is valid.

Subclasses of AbstractStateMachine most likely will want to override this 
method to do something useful according to the particular state change.

=cut

sub change_state {
    my ($self, $new_state) = @_;

    $self->verbose and print STDERR "  changing state to $new_state\n";

    if ( $self->validate_transition( $self->current_state, $new_state, 1 ) ) {
      $self->finalize_state_change( $new_state, 1 );
    }

}


=head2 get_transitions_from

 Purpose  : Returns a list array references that have the indicated state
            in their 'from' slot.

=cut

sub get_transitions_from {
    my ($self, $state) = @_;

    my @trans = ();
    if( ref $self->{'_transition_table'}->{$state}) {
        @trans = @{$self->{'_transition_table'}->{$state}};
    }

    return @trans;
}

=head2 validate_transition

 Purpose  : Determines if the desired state change is defined within 
            the set of registered transitions for this StateMachine.
 Arg      : Two required arguments:
            [0] string defining the name of the "from" state (case sensitive)
            [1] string defining the name of the "to" state (case sensitive)
 Returns  : True if the transition is valid.
            If not valid, throws an exception.
 Throws   : A Bio::Tools::StateMachine::StateException if the desired 
            transition does not exist with the registered transitions
            for this machine.
 Throws   : A Bio::Root::BadParameter if insufficient arguments are given.

=cut

sub validate_transition {
    my ($self, $from_state, $to_state ) = @_;

    #print STDERR "  validating transition $from_state -> $to_state\n";

    if( not( defined($from_state) and defined($to_state))) {
        $self->throw( -class => 'Bio::Root::BadParameter',
                      -text => "Can't validate state transition: Insufficient arguments.");
    }

    my $is_valid = 0;

    foreach my $t ( $self->get_transitions_from( $from_state ) ) {
        if( $t eq $to_state ) {
#        if( $t->[1] eq $to_state ) {
            $is_valid = 1;
            last;
        }
    }

    if( not $is_valid ) {
        $self->throw( -class => 'Bio::Tools::StateMachine::StateException',
                      -text => "The desired state change is not valid for this machine: $from_state -> $to_state");
    }

    #print STDERR "  valid!\n";

    return $to_state;
}

=head2 check_for_new_state

 Purpose : To do whatever checking is necessary to determine if 
            a state transition should occur. 
 Argument : Any necessary data required to determine if the state 
            machine should change to a new state.
 Returns  : A string containing the name of the new state if the 
            state machine should change to a new state. 
            Otherwise returns undef.

This is a virtual method and must be implemented by a subclass to do 
whatever checking is necessary to determine if a state transition should occur.
If not implemented, calling this method will result in a 
Bio::Root::NotImplemented exception.

=cut

sub check_for_new_state {
    my ($self, $data) = @_;
    $self->throw_not_implemented;
}

sub append_input_cache {
    my ($self, $data) = @_;
}

sub get_input_cache {
    my $self = shift;
}

sub clear_input_cache {
    my $self = shift;
}

sub state_change_cache {
    my ($self, $data) = @_;
    if( defined $data ) {
        $self->{'_state_change_cache'} = $data;
    }
    return $self->{'_state_change_cache'};
}

sub clear_state_change_cache {
    my ($self, $data) = @_;
    $self->{'_state_change_cache'} = undef;
}


=head2 finalize_state_change

 Purpose  : Performs routine operations to finish changing state.
            This method should be called at the end of change_state().
 Usage    : finalize_state_change( $new_state, $clear_input_cache );
 Argument : $new_state = the name of the state to change to.
            $clear_input_cache = boolean whether or not to zap whatever 
                                 was in the input cache. Depends on 
                                 the logic of your state machine.

=cut

sub finalize_state_change {
    my ($self, $to_state, $clear_input_cache ) = @_;

    if( $self->paused ) {
        $self->resume_state( $to_state );
    }
    else {
        $self->_set_current_state( $to_state );
    }
    $self->clear_input_cache() if $clear_input_cache;
    $self->append_input_cache( $self->state_change_cache );
    $self->clear_state_change_cache();
}


1;


