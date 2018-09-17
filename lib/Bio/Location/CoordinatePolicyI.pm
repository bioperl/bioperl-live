#
# BioPerl module for Bio::Location::CoordinatePolicyI
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp@gmx.net>
#          and Jason Stajich <jason@bioperl.org>
#
# Copyright Hilmar Lapp, Jason Stajich
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Location::CoordinatePolicyI - Abstract interface for objects implementing
a certain policy of computing integer-valued coordinates of a Location

=head1 SYNOPSIS

    # get a location, e.g., from a SeqFeature
    $location = $feature->location();
    # examine its coordinate computation policy
    print "Location of feature ", $feature->primary_tag(), " employs a ",
          ref($location->coordinate_policy()), 
          " instance for coordinate computation\n";
    # change the policy, e.g. because the user chose to do so
    $location->coordinate_policy(Bio::Location::NarrowestCoordPolicy->new());

=head1 DESCRIPTION

Objects implementing this interface are used by Bio::LocationI
implementing objects to determine integer-valued coordinates when
asked for it. While this may seem trivial for simple locations, there
are different ways to do it for fuzzy or compound (split)
locations. Classes implementing this interface implement a certain
policy, like 'always widest range', 'always smallest range', 'mean for
BETWEEN locations', etc. By installing a different policy object in a
Location object, the behaviour of coordinate computation can be changed
on-the-fly, and with a single line of code client-side.

=head1 FEEDBACK

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Hilmar Lapp, Jason Stajich

Email hlapp@gmx.net, jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::Location::CoordinatePolicyI;
use strict;

use base qw(Bio::Root::RootI);

=head2 start

  Title   : start
  Usage   : $start = $policy->start($location);
  Function: Get the integer-valued start coordinate of the given location as
            computed by this computation policy.
  Returns : A positive integer number.
  Args    : A Bio::LocationI implementing object.

=cut

sub start {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 end

  Title   : end
  Usage   : $end = $policy->end($location);
  Function: Get the integer-valued end coordinate of the given location as
            computed by this computation policy.
  Returns : A positive integer number.
  Args    : A Bio::LocationI implementing object.

=cut

sub end {
    my ($self) = @_;
    $self->throw_not_implemented();
}

1;
