#
# BioPerl module for Bio::DB::ReferenceI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chris Fields <cjfields at bioperl dot org>
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::DB::ReferenceI - A RandomAccessI-like abstract interface for
retrieving Reference data from a sequence database and returning
Bio::Annotation::Reference objects 

=head1 SYNOPSIS

  #
  # get a database object somehow using a concrete class
  #

  $ref = $db->get_Reference_by_id('123456');

  #
  # $ref is a Bio::Annotation::Reference object
  #

=head1 DESCRIPTION

This is a pure interface class - in other words, all this does is define
methods which other (concrete) classes will actually implement. 

The Bio::DB::ReferenceI class defines methods used to retrieve reference data
from a sequence.  This is returned in the form of Bio::Annotation::Reference
objects.

At the moment it is just the ability to make Bio::Annotation::Reference
objects from a name or unique id (id), an accession number (acc), and so on.

=head1 CONTACT

Ewan Birney originally wrote Bio::DB::RandomAccessI, from which this class
is based.

=head2 Mailing Lists

User feedback is an integral part of the 
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@lists.open-bio.org               - General discussion
  http://www.bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via the web.

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR 

Email cjfields at bioperl dot org

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::ReferenceI;

use strict;

=head2 get_Reference_by_id

 Title   : get_Reference_by_id
 Usage   : $ref = $db->get_Reference_by_id('123456')
 Function: Gets a Bio::Annotation::Reference-implementing object by its name (id)
 Returns : a Bio::Annotation::Reference object or undef if not found
 Args    : the id (as a string) of a sequence

=cut

sub get_Reference_by_id{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 get_Reference_by_acc

 Title   : get_Reference_by_acc
 Usage   : $ref = $db->get_Reference_by_acc('X77802');
 Function: Gets a Bio::Annotation::Reference object by accession number
 Returns : A Bio::Annotation::Reference object or undef if not found
 Args    : accession number (as a string)
 Throws  : "more than one sequences correspond to this accession"
            if the accession maps to multiple primary ids and
            method is called in a scalar context

=cut

sub get_Reference_by_acc{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 get_Reference_by_version

 Title   : get_Reference_by_version
 Usage   : $ref = $db->get_Reference_by_version('X77802.1');
 Function: Gets a Bio::Annotation::Reference object by sequence version
 Returns : A Bio::Annotation::Reference object
 Args    : accession.version (as a string)
 Throws  : "acc.version does not exist" exception

=cut

sub get_Reference_by_version{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

## End of Package

1;
