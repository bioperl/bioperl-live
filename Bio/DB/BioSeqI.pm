# POD documentation - main docs before the code

=head1 NAME

Bio::DB::BioSeqI - Abstract interface for a sequence database

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

The Bio::DB::BioSeqI class defines what methods a generic database class
should have. At the moment it is just the ability to make Bio::Seq objects
from a name (id) or a accession number.

=head1 CONTACT

Ewan Birney originally wrote this class.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the bugs and 
their resolution. Bug reports can be submitted via email or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::DB::BioSeqI;
use vars qw($AUTOLOAD @ISA @EXPORT_OK);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;

@ISA = qw(Bio::Root::Object Exporter);
@EXPORT_OK = qw();

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  # set stuff in self from @args
  return $make; # success - we hope!
}


=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : $seq = $db->get_Seq_by_id('ROA1_HUMAN')
 Function: Gets a Bio::Seq object by its name
 Returns : a Bio::Seq object
 Args    : the id (as a string) of a sequence
 Throws  : "id does not exist" exception


=cut

sub get_Seq_by_id{
   my ($self,@args) = @_;

   $self->throw("Abstract database call of get_Seq_by_id. Your database has not implemented this method!");
}

=head2 get_Seq_by_acc

 Title   : get_Seq_by_acc
 Usage   : $seq = $db->get_Seq_by_acc('X77802');
 Function: Gets a Bio::Seq object by accession number
 Returns : A Bio::Seq object
 Args    : accession number (as a string)
 Throws  : "acc does not exist" exception


=cut

sub get_Seq_by_acc{
   my ($self,@args) = @_;

   $self->throw("Abstract database call of get_Seq_by_acc. Your database has not implemented this method!");

}


## End of Package

1;

__END__

