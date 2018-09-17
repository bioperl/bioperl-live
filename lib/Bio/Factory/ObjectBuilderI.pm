#
# BioPerl module for Bio::Factory::ObjectBuilderI
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

Bio::Factory::ObjectBuilderI - Interface for an object builder

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

An object builder is different from an object factory in that it
accumulates information for the object and finally, or constantly,
depending on the implementation, builds the object. It also allows for
implementations that can tell the information feed in which kind of
information the builder is interested in which not. In addition, the
implementation may choose to filter, transform, or completely ignore
certain content it is fed for certain slots.

Implementations will hence be mostly used by stream-based parsers to
parse only desired content, and/or skip over undesired entries.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Factory::ObjectBuilderI;
use strict;
use Carp;

use base qw(Bio::Root::RootI);

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
    shift->throw_not_implemented();
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
    shift->throw_not_implemented();
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
    shift->throw_not_implemented();
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
    shift->throw_not_implemented();
}

1;
