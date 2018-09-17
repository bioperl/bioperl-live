# POD documentation - main docs before the code
#
#

=head1 NAME

Bio::DB::RandomAccessI - Abstract interface for a sequence database

=head1 SYNOPSIS

  #
  # get a database object somehow using a concrete class
  #

  $seq = $db->get_Seq_by_id('ROA1_HUMAN');

  #
  # $seq is a Bio::Seq object
  #

=head1 DESCRIPTION

This is a pure interface class - in other words, all this does is define
methods which other (concrete) classes will actually implement.

The Bio::DB::RandomAccessI class defines what methods a generic database class
should have. At the moment it is just the ability to make Bio::Seq objects
from a name (id) or an accession number.

=head1 CONTACT

Ewan Birney E<lt>birney@ebi.ac.ukE<gt> originally wrote this class.

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::DB::RandomAccessI;

use strict;

use Bio::Root::RootI;

use base qw(Bio::Root::Root);

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id('ROA1_HUMAN')
 Function: Gets a Bio::Seq object by its name
 Returns : a Bio::Seq object or undef if not found
 Args    : the id (as a string) of a sequence,

=cut

sub get_Seq_by_id{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 get_Seq_by_acc

 Title   : get_Seq_by_acc
 Usage   : $seq = $db->get_Seq_by_acc('X77802');
           $seq = $db->get_Seq_by_acc(Locus => 'X77802');
 Function: Gets a Bio::Seq object by accession number
 Returns : A Bio::Seq object or undef if not found
 Args    : accession number (as a string), or a two
               element list consisting of namespace=>accession
 Throws  : "more than one sequences correspond to this accession"
            if the accession maps to multiple primary ids and
            method is called in a scalar context

NOTE: The two-element form allows you to choose the namespace for the
accession number.

=cut

sub get_Seq_by_acc{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}


=head2 get_Seq_by_version

 Title   : get_Seq_by_version
 Usage   : $seq = $db->get_Seq_by_version('X77802.1');
 Function: Gets a Bio::Seq object by sequence version
 Returns : A Bio::Seq object
 Args    : accession.version (as a string)
 Throws  : "acc.version does not exist" exception

=cut


sub get_Seq_by_version{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

## End of Package

1;

