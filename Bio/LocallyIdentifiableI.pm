package Bio::LocallyIdentifiableI;

# $Id $
# Interface for uniquely identifiable objects

=head1 NAME

Bio::LocallyIdentifiableI - interface for objects with unique identifiers

=head1 SYNOPSIS

    # to test this is an identifiable object

    $obj->isa("Bio::LocallyIdentifiableI") || 
      $obj->throw("$obj does not implement the Bio::LocallyIdentifiableI interface");

    $unique_id = $obj->unique_id();

=head1 DESCRIPTION

This interface describes a single method, expected on all uniquely
identifiable bioperl objects.  Unique identification in this context
means that the string returned from the unique_id() method is
sufficient to differentiate the object from any others presently in
memory on the same computer.  See L<Bio::GloballyIdentifiableI> for a
stricter requirement.

Note: An object that implements Bio::LocallyIdentifiable that is
unable to provide a unique identifier B<must> return undef from its
unique_id() method.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                       - General discussion
  http://bio.perl.org/MailList.html           - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR

Paul Edlefsen E<lt>paul@systemsbiology.orgE<gt>.

Copyright (c) 2003 Institute for Systems Biology

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=head1 CONTRIBUTORS

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

use vars qw( @ISA );
use strict;
use Bio::Root::RootI;

@ISA = qw( Bio::Root::RootI );

=head2 unique_id

 Title   : unique_id
 Usage   : $string = $obj->unique_id()
 Function: Differentiate between this object and others
 Returns : A string which is unique to this object, or undef in none exists
 Status  : Public

Note that an object that implements Bio::LocallyIdentifiable that is
unable to provide a unique identifier B<must> return undef from its
unique_id() method.

=cut

sub unique_id {
   my ($self) = @_;
   $self->throw_not_implemented();
}

1;

__END__
