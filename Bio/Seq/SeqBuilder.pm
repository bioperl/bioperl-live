# $Id$
#
# BioPerl module for Bio::Seq::SeqBuilder
#
# Cared for by Hilmar Lapp <hlapp at gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

#
# (c) Hilmar Lapp, hlapp at gmx.net, 2002.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2002.
#
# You may distribute this module under the same terms as perl itself.
# Refer to the Perl Artistic License (see the license accompanying this
# software package, or see http://www.perl.com/language/misc/Artistic.html)
# for the terms under which you may use, modify, and redistribute this module.
# 
# THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::SeqBuilder - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::SeqBuilder;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Factory::ObjectBuilderI;

@ISA = qw(Bio::Root::Root Bio::Factory::ObjectBuilderI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Seq::SeqBuilder();
 Function: Builds a new Bio::Seq::SeqBuilder object 
 Returns : an instance of Bio::Seq::SeqBuilder
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  $self->{'wanted_slots'} = [];
  $self->want_all(1);

  return $self;
}

=head1 Methods for implementing L<Bio::Factory::ObjectBuilderI>

=cut

=head2 want_slot

 Title   : want_slot
 Usage   :
 Function: Whether or not the object builder wants to populate the
           specified slot of the object to be built.

           The slot can be specified either as the name of the
           respective method, or the initialization parameter that
           would be otherwise passed to new() of the object to be
           built.

 Example :
 Returns : TRUE if the object builder wants to populate the slot, and
           FALSE otherwise.
 Args    : the name of the slot (a string)


=cut

sub want_slot{
    my ($self,@args) = @_;


}

=head2 add_slot_value

 Title   : add_slot_value
 Usage   :
 Function: Adds one or more values to the specified slot of the object
           to be built.

           Naming the slot is the same as for want_slot().

           The object builder may further filter the content to be
           set, or even completely ignore the request.

           If this method reports failure, the caller should not add
           more values to the same slot. In addition, the caller may
           find it appropriate to abandon the object being built
           altogether.

 Example :
 Returns : TRUE on success, and FALSE otherwise
 Args    : the name of the slot (a string)
           parameters determining the value to be set


=cut

sub add_slot_value{
   my ($self,@args) = @_;


}

=head2 want_object

 Title   : want_object
 Usage   :
 Function: Whether or not the object builder is still interested in
           continuing with the object being built.

           If this method returns FALSE, the caller should not add any
           more values to slots, or otherwise risks that the builder
           throws an exception. In addition, make_object() is likely
           to return undef after this method returned FALSE.

 Example :
 Returns : TRUE if the object builder wants to continue building
           the present object, and FALSE otherwise.
 Args    : none


=cut

sub want_object{
   my ($self,@args) = @_;


}

=head2 make_object

 Title   : make_object
 Usage   :
 Function: Get the built object.

           This method is allowed to return undef if no value has ever
           been added since the last call to make_object(), or if
           want_object() returned FALSE (or would have returned FALSE)
           before calling this method.

           For an implementation that allows consecutive building of
           objects, a caller must call this method once, and only
           once, between subsequent objects to be built. I.e., a call
           to make_object implies 'end_object.'

 Example :
 Returns : the object that was built
 Args    : none


=cut

sub make_object{
   my ($self,@args) = @_;


}

=head1 Implementation specific methods

=cut

=head2 get_wanted_slots

 Title   : get_wanted_slots
 Usage   : $obj->get_wanted_slots($newval)
 Function: 
 Example : 
 Returns : value of wanted_slots (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub get_wanted_slots{
    my $self = shift;

    return @{$self->{'wanted_slots'}};
}

=head2 add_wanted_slots

 Title   : add_wanted_slots
 Usage   :
 Function: Adds the specified slots to the list of wanted slots.
 Example :
 Returns : TRUE
 Args    : an array of slot names (strings)


=cut

sub add_wanted_slots{
    my ($self,@slots) = @_;

    my $myslots = $self->{'wanted_slots'};
    foreach my $slot (@slots) {
	if(! grep { $slot eq $_; } @$myslots) {
	    push(@$myslots, $slot);
	}
    }
    return 1;
}

=head2 remove_wanted_slots

 Title   : remove_wanted_slots
 Usage   :
 Function: Removes all wanted slots added previously through
           add_wanted_slots().
 Example :
 Returns : the previous list of wanted slot names
 Args    : none


=cut

sub remove_wanted_slots{
    my $self = shift;
    my @slots = $self->get_wanted_slots();
    $self->{'wanted_slots'} = [];
    return @slots;
}

=head2 want_none

 Title   : want_none
 Usage   :
 Function: Disables all slots. After calling this method, want_slot()
           will return FALSE regardless of slot name.
 Example :
 Returns : TRUE
 Args    : none


=cut

sub want_none{
    my $self = shift;

    $self->want_all(0);
    $self->remove_wanted_slots();
    return 1;
}

=head2 want_all

 Title   : want_all
 Usage   : $obj->want_all($newval)
 Function: Whether or not this sequence object builder wants to
           populate all slots that the object has.

           This will be ON by default. Call $obj->want_none() to
           disable all slots.

 Example : 
 Returns : TRUE if this builder wants to populate all slots, and
           FALSE otherwise.
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub want_all{
    my $self = shift;

    return $self->{'want_all'} = shift if @_;
    return $self->{'want_all'};
}

=head2 sequence_factory

 Title   : sequence_factory
 Usage   : $obj->sequence_factory($newval)
 Function: Get/set the sequence factory to be used by this object
           builder.
 Example : 
 Returns : the Bio::Factory::SequenceFactoryI implementing object to use
 Args    : on set, new value (a Bio::Factory::SequenceFactoryI
           implementing object or undef, optional)


=cut

sub sequence_factory{
    my $self = shift;

    return $self->{'sequence_factory'} = shift if @_;
    return $self->{'sequence_factory'};
}

1;
