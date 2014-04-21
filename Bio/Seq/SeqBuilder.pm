#
# BioPerl module for Bio::Seq::SeqBuilder
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

Bio::Seq::SeqBuilder - Configurable object builder for sequence stream parsers

=head1 SYNOPSIS

   use Bio::SeqIO;

   # usually you won't instantiate this yourself - a SeqIO object -
   # you will have one already
   my $seqin = Bio::SeqIO->new(-fh => \*STDIN, -format => "genbank");
   my $builder = $seqin->sequence_builder();

   # if you need only sequence, id, and description (e.g. for 
   # conversion to FASTA format):
   $builder->want_none();
   $builder->add_wanted_slot('display_id','desc','seq');

   # if you want everything except the sequence and features
   $builder->want_all(1); # this is the default if it's untouched
   $builder->add_unwanted_slot('seq','features');

   # if you want only human sequences shorter than 5kb and skip all
   # others
   $builder->add_object_condition(sub {
       my $h = shift;
       return 0 if $h->{'-length'} > 5000;
       return 0 if exists($h->{'-species'}) &&
                   ($h->{'-species'}->binomial() ne "Homo sapiens");
       return 1;
   });

   # when you are finished with configuring the builder, just use
   # the SeqIO API as you would normally
   while(my $seq = $seqin->next_seq()) {
       # do something
   }

=head1 DESCRIPTION

This is an implementation of L<Bio::Factory::ObjectBuilderI> used by
parsers of rich sequence streams. It provides for a relatively
easy-to-use configurator of the parsing flow.

Configuring the parsing process may be for you if you need much less
information, or much less sequence, than the stream actually
contains. Configuration can in both cases speed up the parsing time
considerably, because unwanted sections or the rest of unwanted
sequences are skipped over by the parser. This configuration could
also conserve memory if you're running out of available RAM.

See the methods of the class-specific implementation section for
further documentation of what can be configured.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::SeqBuilder;
use strict;

# Object preamble - inherits from Bio::Root::Root


use base qw(Bio::Root::Root Bio::Factory::ObjectBuilderI);

my %slot_param_map = ("add_SeqFeature" => "features",
		      );
my %param_slot_map = ("features"       => "add_SeqFeature",
		      );

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Seq::SeqBuilder->new();
 Function: Builds a new Bio::Seq::SeqBuilder object 
 Returns : an instance of Bio::Seq::SeqBuilder
 Args    :

=cut

sub new {
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);

    $self->{'wanted_slots'} = [];
    $self->{'unwanted_slots'} = [];
    $self->{'object_conds'} = [];
    $self->{'_objhash'} = {};
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

           Note that usually only the parser will call this
           method. Use add_wanted_slots and add_unwanted_slots for
           configuration.

 Example :
 Returns : TRUE if the object builder wants to populate the slot, and
           FALSE otherwise.
 Args    : the name of the slot (a string)


=cut

sub want_slot{
	my ($self,$slot) = @_;
	my $ok = 0;

	$slot = substr($slot,1) if substr($slot,0,1) eq '-';
	if($self->want_all()) {
	foreach ($self->get_unwanted_slots()) {
		# this always overrides in want-all mode
		return 0 if($slot eq $_);
	}
	if(! exists($self->{'_objskel'})) {
		$self->{'_objskel'} = $self->sequence_factory->create_object();
	}
	if(exists($param_slot_map{$slot})) {
		$ok = $self->{'_objskel'}->can($param_slot_map{$slot});
	} else {
		$ok = $self->{'_objskel'}->can($slot);
	}
	return $ok if $ok;
	# even if the object 'cannot' do this slot, it might have been
	# added to the list of wanted slot, so carry on
}
	foreach ($self->get_wanted_slots()) {
		if($slot eq $_) {
			$ok = 1;
			last;
		}
	}
	return $ok;
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

           This implementation will allow the caller to overwrite the
           return value from want_slot(), because the slot is not
           checked against want_slot().

           Note that usually only the parser will call this method,
           but you may call it from anywhere if you know what you are
           doing. A derived class may be used to further manipulate
           the value to be added.

 Example :
 Returns : TRUE on success, and FALSE otherwise
 Args    : the name of the slot (a string)
           parameters determining the value to be set

                 OR

           alternatively, a list of slotname/value pairs in the style
           of named parameters as they would be passed to new(), where
           each element at an even index is the parameter (slot) name
           starting with a dash, and each element at an odd index is
           the value of the preceding name.

=cut

sub add_slot_value{
	my ($self,$slot,@args) = @_;

	my $h = $self->{'_objhash'};
	return unless $h;
	# multiple named parameter variant of calling?
	if((@args > 1) && (@args % 2) && (substr($slot,0,1) eq '-')) {
		unshift(@args, $slot);
		while(@args) {
			my $key = shift(@args);
			$h->{$key} = shift(@args);
		}
	} else {
		if($slot eq 'add_SeqFeature') {
			$slot = '-'.$slot_param_map{$slot};
			$h->{$slot} = [] unless $h->{$slot};
			push(@{$h->{$slot}}, @args);
		} else {
			$slot = '-'.$slot unless substr($slot,0,1) eq '-';
			$h->{$slot} = $args[0];
		}
	}
	return 1;
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

           Note that usually only the parser will call this
           method. Use add_object_condition for configuration.

 Example :
 Returns : TRUE if the object builder wants to continue building
           the present object, and FALSE otherwise.
 Args    : none

=cut

sub want_object{
	my $self = shift;

	my $ok = 1;
	foreach my $cond ($self->get_object_conditions()) {
		$ok = &$cond($self->{'_objhash'});
		last unless $ok;
	}
	delete $self->{'_objhash'} unless $ok;
	return $ok;
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
	my $self = shift;

	my $obj;
	if(exists($self->{'_objhash'}) && %{$self->{'_objhash'}}) {
		$obj = $self->sequence_factory->create_object(%{$self->{'_objhash'}});
	}
	$self->{'_objhash'} = {}; # reset
	return $obj;
}

=head1 Implementation specific methods

These methods allow to conveniently configure this sequence object
builder as to which slots are desired, and under which circumstances a
sequence object should be abandoned altogether. The default mode is
want_all(1), which means the builder will report all slots as wanted
that the object created by the sequence factory supports.

You can add specific slots you want through add_wanted_slots(). In
most cases, you will want to call want_none() before in order to relax
zero acceptance through a list of wanted slots.

Alternatively, you can add specific unwanted slots through
add_unwanted_slots(). In this case, you will usually want to call
want_all(1) before (which is the default if you never touched the
builder) to restrict unrestricted acceptance.

I.e., want_all(1) means want all slots except for the unwanted, and
want_none() means only those explicitly wanted.

If a slot is in both the unwanted and the wanted list, the following
rules hold. In want-all mode, the unwanted list overrules. In
want-none mode, the wanted list overrides the unwanted list. If this
is confusing to you, just try to avoid having slots at the same time
in the wanted and the unwanted lists.

=cut

=head2 get_wanted_slots

 Title   : get_wanted_slots
 Usage   : $obj->get_wanted_slots($newval)
 Function: Get the list of wanted slots
 Example : 
 Returns : a list of strings
 Args    : 


=cut

sub get_wanted_slots{
	my $self = shift;

	return @{$self->{'wanted_slots'}};
}

=head2 add_wanted_slot

 Title   : add_wanted_slot
 Usage   :
 Function: Adds the specified slots to the list of wanted slots.
 Example :
 Returns : TRUE
 Args    : an array of slot names (strings)

=cut

sub add_wanted_slot{
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

=head2 get_unwanted_slots

 Title   : get_unwanted_slots
 Usage   : $obj->get_unwanted_slots($newval)
 Function: Get the list of unwanted slots.
 Example : 
 Returns : a list of strings
 Args    : none

=cut

sub get_unwanted_slots{
	my $self = shift;

	return @{$self->{'unwanted_slots'}};
}

=head2 add_unwanted_slot

 Title   : add_unwanted_slot
 Usage   :
 Function: Adds the specified slots to the list of unwanted slots.
 Example :
 Returns : TRUE
 Args    : an array of slot names (strings)

=cut

sub add_unwanted_slot{
	my ($self,@slots) = @_;

	my $myslots = $self->{'unwanted_slots'};
	foreach my $slot (@slots) {
		if(! grep { $slot eq $_; } @$myslots) {
			push(@$myslots, $slot);
		}
	}
	return 1;
}

=head2 remove_unwanted_slots

 Title   : remove_unwanted_slots
 Usage   :
 Function: Removes the list of unwanted slots added previously through
           add_unwanted_slots().
 Example :
 Returns : the previous list of unwanted slot names
 Args    : none

=cut

sub remove_unwanted_slots{
	my $self = shift;
	my @slots = $self->get_unwanted_slots();
	$self->{'unwanted_slots'} = [];
	return @slots;
}

=head2 want_none

 Title   : want_none
 Usage   :
 Function: Disables all slots. After calling this method, want_slot()
           will return FALSE regardless of slot name.

           This is different from removed_wanted_slots() in that it
           also sets want_all() to FALSE. Note that it also resets the
           list of unwanted slots in order to avoid slots being in
           both lists.

 Example :
 Returns : TRUE
 Args    : none

=cut

sub want_none{
	my $self = shift;

	$self->want_all(0);
	$self->remove_wanted_slots();
	$self->remove_unwanted_slots();
	return 1;
}

=head2 want_all

 Title   : want_all
 Usage   : $obj->want_all($newval)
 Function: Whether or not this sequence object builder wants to
           populate all slots that the object has. Whether an object
           supports a slot is generally determined by what can()
           returns. You can add additional 'virtual' slots by calling
           add_wanted_slot.

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

=head2 get_object_conditions

 Title   : get_object_conditions
 Usage   :
 Function: Get the list of conditions an object must meet in order to
           be 'wanted.' See want_object() for where this is used.

           Conditions in this implementation are closures (anonymous
           functions) which are passed one parameter, a hash reference
           the keys of which are equal to initialization
           paramaters. The closure must return TRUE to make the object
           'wanted.'

           Conditions will be implicitly ANDed.

 Example :
 Returns : a list of closures
 Args    : none

=cut

sub get_object_conditions{
	my $self = shift;

	return @{$self->{'object_conds'}};
}

=head2 add_object_condition

 Title   : add_object_condition
 Usage   :
 Function: Adds a condition an object must meet in order to be 'wanted.'
           See want_object() for where this is used.

           Conditions in this implementation must be closures
           (anonymous functions). These will be passed one parameter,
           which is a hash reference with the sequence object
           initialization parameters being the keys.

           Conditions are implicitly ANDed. If you want other
           operators, perform those tests inside of one closure
           instead of multiple.  This will also be more efficient.

 Example :
 Returns : TRUE
 Args    : the list of conditions

=cut

sub add_object_condition{
	my ($self,@conds) = @_;

	if(grep { ref($_) ne 'CODE'; } @conds) {
		$self->throw("conditions against which to validate an object ".
						 "must be anonymous code blocks");
	}
	push(@{$self->{'object_conds'}}, @conds);
	return 1;
}

=head2 remove_object_conditions

 Title   : remove_object_conditions
 Usage   :
 Function: Removes the conditions an object must meet in order to be
           'wanted.'
 Example :
 Returns : The list of previously set conditions (an array of closures)
 Args    : none

=cut

sub remove_object_conditions{
	my $self = shift;
	my @conds = $self->get_object_conditions();
	$self->{'object_conds'} = [];
	return @conds;
}

=head1 Methods to control what type of object is built

=cut

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

	if(@_) {
		delete $self->{'_objskel'};
		return $self->{'sequence_factory'} = shift;
	}
	return $self->{'sequence_factory'};
}

1;
