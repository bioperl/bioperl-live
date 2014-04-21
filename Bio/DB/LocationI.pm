#
# BioPerl module for Bio::DB::LocationI
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

Bio::DB::LocationI - A RandomAccessI-like abstract interface for
retrieving location data from a sequence database and returning
Bio::LocationI objects 

=head1 SYNOPSIS

  #
  # get a database object somehow using a concrete class
  #

  $loc = $db->get_Location_by_id('123456');

  #
  # $loc is a Bio::LocationI object
  #

=head1 DESCRIPTION

This is a pure interface class - in other words, all this does is define
methods which other (concrete) classes will actually implement. 

The Bio::DB::LocationI class defines methods used to retrieve location data
from a sequence.  This is returned in the form of Bio::LocationI objects,
which can include:

Bio::Location::Simple
Bio::Location::Fuzzy
Bio::Location::Split

At the moment it is just the ability to make Bio::LocationI objects
from a name or unique id (id), an accession number (acc), and so on.

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

package Bio::DB::LocationI;

use strict;

use Bio::Root::RootI;

use base qw(Bio::Root::Root);

=head2 get_Location_by_id

 Title   : get_Location_by_id
 Usage   : $loc = $db->get_Location_by_id('123456')
 Function: Gets a Bio::LocationI-implementing object by its name (id)
 Returns : a Bio::LocationI object or undef if not found
 Args    : the id (as a string) of a sequence

=cut

sub get_Location_by_id{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 get_Location_by_acc

 Title   : get_Location_by_acc
 Usage   : $loc = $db->get_Location_by_acc('X77802');
 Function: Gets a Bio::LocationI object by accession number
 Returns : A Bio::LocationI object or undef if not found
 Args    : accession number (as a string)
 Throws  : "more than one sequences correspond to this accession"
            if the accession maps to multiple primary ids and
            method is called in a scalar context

=cut

sub get_Location_by_acc{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 get_Location_by_version

 Title   : get_Location_by_version
 Usage   : $loc = $db->get_Location_by_version('X77802.1');
 Function: Gets a Bio::LocationI object by sequence version
 Returns : A Bio::LocationI object
 Args    : accession.version (as a string)
 Throws  : "acc.version does not exist" exception

=cut

sub get_Location_by_version{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

## End of Package

1;
