#-----------------------------------------------------------------
# $Id$
#
# BioPerl module Bio::Tools::StateMachine::IOStateMachine
#
# Cared for by Steve Chervitz <steve_chervitz@affymetrix.com>
#
# You may distribute this module under the same terms as perl itself
#-----------------------------------------------------------------

=head1 NAME

Bio::Tools::StateMachine::IOStateMachine - IO-based implementation of AbstractStateMachine

=head1 SYNOPSIS

    use Bio::Tools::IOStateMachine;

    # A state machine that reads input from a file
    my $sm = Bio::Tools::IOStateMachine->new( -file => 'data.txt' );

    # A state machine that reads input from a STDIN
    my $sm = Bio::Tools::IOStateMachine->new();

    # A state machine that reads input from a STDIN
    # and times out if input doesn't arrive within 30 seconds.
    my $sm = Bio::Tools::IOStateMachine->new( -timeout_sec => 30 );


=head1 DESCRIPTION

An implementation of AbstractStateMachine that samples an input stream
to determine whether a state change has occurred.

=head1 EXAMPLES

To get a feel for how to use this, have look at
Bio::SearchIO::psiblast which subclasses IOStateMachine.


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

steve_chervitz@affymetrix.com

See the L<FEEDBACK> section for where to send bug reports and comments.

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

package Bio::Tools::StateMachine::IOStateMachine;

use strict;
use vars qw( @ISA @EXPORT_OK );

use Bio::Root::IO;
use Bio::Tools::StateMachine::AbstractStateMachine qw($INITIAL_STATE $FINAL_STATE);

@ISA = qw( Bio::Root::IO
           Bio::Tools::StateMachine::AbstractStateMachine
         );

# Propagating the initial and final states from AbstractStateMachine
@EXPORT_OK = qw( $INITIAL_STATE $FINAL_STATE );

=head2 _init_state_machine()

 Argument : Named parameter -TIMEOUT_SEC => seconds,
            to specify the number of seconds to allow before throwing
            an exception if input fails to arrive within that amount of time. 

=cut

sub _init_state_machine {
    my($self, @args) = @_;

    $self->SUPER::_init_state_machine(@args);

    my ($timeout) = $self->_rearrange( [qw(TIMEOUT_SECS)], @args);

    if( defined $timeout ) {
	if($timeout =~ /^\d+$/ ) {
	    $self->{'_timeout_secs'} = $timeout;
	}
	else {
	    $self->throw(-class =>'Bio::Root::BadParameter',
			 -text => "TIMEOUT_SECS must be a number: $timeout",
			 -value => $timeout
			);
	}
    }
}

=head2 check_for_new_state()
 
 Purpose  : Obtains data from the input stream to be checked 
            for the existence of a new state.
 Usage    : check_for_new_state( [$ignore_blank_lines] );
 Argument : boolean: true if you want to ignore blank lines
 Returns  : the next chunk of input ($/ is not altered)
            If there is no more input, returns undef.

Subclasses should override this method and call it to obtain
the chunk of data for new state testing. 

=cut

sub check_for_new_state {
    my ($self, $ignore_blank_lines) = @_;

    $self->verbose and print STDERR "Checking for new state...\n";

    my $chunk = $self->next_input_chunk();

    # Determine if we're supposed to ignore blanks and if so, loop
    # until we're either out of input or hit a non-blank line.
    if( $ignore_blank_lines and $chunk =~ /^\s*$/ ) {
        while(  $chunk = $self->next_input_chunk()) {
            last unless not $chunk or $chunk =~ /^\s*$/;
        }
    }

    $self->verbose and print STDERR "  Input chunk: " . $chunk, "\n";

    return $chunk;
}

=head2 next_input_chunk()

 Argument : n/a
 Returns  : The next chunk of input data from the IO stream
            To be used in determining what state the machine should be in.

=cut

sub next_input_chunk {
    my $self = shift;

    $self->verbose and print STDERR "Getting next input chunk...\n", ;

    if(not defined $self->{'_alarm_available'}) {
        $self->_check_if_alarm_available();
    }

    $SIG{ALRM} = sub { die "Timed out!"; };

    my $chunk;

    eval {
        if( $self->{'_alarm_available'} and defined $self->{'_timeout_secs'}) {
	    alarm($self->{'_timeout_secs'});
	}

        $chunk = $self->_readline();

    };
    if($@ =~ /Timed out!/) {
	 $self->throw(-class => 'Bio::Root::IOException',
                      -text => "Timed out while waiting for input (timeout=$self->{'_timeout_secs'}s).");
     } elsif($@ =~ /\S/) {
         my $err = $@;
         $self->throw(-class => 'Bio::Root::IOException',
                      -text => "Unexpected error during readline: $err");
    }

    return $chunk;
}



# alarm() not available (ActiveState perl for win32 doesn't have it.
# See jitterbug PR#98)
sub _check_if_alarm_available {
    my $self = shift;
    eval {
        alarm(0);
    };
    if($@) {
        $self->{'_alarm_available'} = 0;
    }
    else {
        $self->{'_alarm_available'} = 1;
    }
}

sub append_input_cache {
    my ($self, $data) = @_;
    push( @{$self->{'_input_cache'}}, $data) if defined $data;
}

sub get_input_cache {
    my $self = shift;
    my @cache =  ();
    if( ref $self->{'_input_cache'} ) {
       @cache = @{$self->{'_input_cache'}};
    }
    return @cache;
}

sub clear_input_cache {
    my $self = shift;
    @{$self->{'_input_cache'}} = ();
}



1;



