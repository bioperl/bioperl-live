
#
# BioPerl module for Bio::Index::Abstract
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Index::Abstract - Abstract interface for indexing a flat file

=head1 SYNOPSIS

You should not be using this module directly

=head1 DESCRIPTION

This object provides the basic mechanism to associate positions
in files with names. The position and filenames are stored in DBM
which can then be accessed later on. It is the equivalent of flat
file indexing (eg, SRS or efetch).

This object is the guts to the mechanism, which will be used by the
specific objects inherieting from it.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via
email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email - birney@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Index::Abstract;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

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


=head2 start_record

 Title   : start_record
 Usage   : Abstract - should be overriden by the subclass
 Function: This function is called for every line in the 
           indexing function. It should return 1 when the
           start of a record occurs
 Example :
 Returns : 1 or 0
 Args    :


=cut

sub start_record{
   my ($self,@args) = @_;

   $self->throw("In Bio::Index::Abstract, defaulting to start_record - the sub class should be providing this");
}

=head2 end_record

 Title   : end_record
 Usage   : Internal - should be overriden in subclass
 Function: This function is called for every line in the indexing function
           It should return 1 at the end of a record and 0 otherwise
 Example :
 Returns : 
 Args    :


=cut

sub end_record{
   my ($self,@args) = @_;
   
   $self->throw("In Bio::Index::Abstract - defaulted for end_record. Sub class should provide this");
}

=head2 record_id

 Title   : record_id
 Usage   : Internal - should be overridden in subclass
 Function: should return the id (if wanted) for this entry
           to be indexed under
 Example :
 Returns : scalar or undef
 Args    :


=cut

sub record_id{
   my ($self,@args) = @_;

   $self->throw("In Bio::Index::Abstract - defaulted for record_id. Sub class should provide this");
}


=head2 _index_hash

 Title   : _index_hash
 Usage   : %hash = $obj->_index_hash();
 Function: returns the hash from the dbm index
 Returns : a hash
 Args    : none


=cut

sub _index_hash{
   my ($self,@args) = @_;

}

=head2 _dbm_file

 Title   : _dbm_file
 Usage   : $obj->_dbm_file($newval)
 Function: _dbm_file should be a fully qualified 
           filename for the dbm file
 Returns : value of _dbm_file
 Args    : newvalue (optional)


=cut

sub _dbm_file{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_dbm_file'} = $value;
    }
    return $obj->{'_dbm_file'};
}



