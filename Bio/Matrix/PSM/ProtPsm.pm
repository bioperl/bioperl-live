#---------------------------------------------------------

#ISA ProtMatrix, HAS InstanceSite

=head1 NAME

Bio::Matrix::PSM::ProtPsm - handle combination of site matricies

=head1 SYNOPSIS

  use Bio::Matrix::PSM::IO;

  #To get a ProtPsm object from a file use the Psm parser:
  my $psmIO =  Bio::Matrix::PSM::IO->new(-format=>'psiblast', -file=>$file);

  # Now go through all entities in the file with next_psm, which
  # returns a Psm object see Bio::Matrix::PSM::IO for detailed
  # documentation (matrix predictions or matrix sequence matches or
  # both):

  while (my $psm=$psmIO->next_psm) {
     my %psm_header = $psm->header;
     my $ic    = $psm_header{IC};
     my $sites = $psm_header{sites};
     my $width = $psm_header{width};
     my $score = $psm_header{e_val};
     my $IUPAC = $psm->IUPAC;
     my $instances = $psm->instances;
     foreach my $instance (@{$instances}) {
       my $id = $instance->primary_id;
       #Do something with the id
     }
   }

=head1 DESCRIPTION

To handle a combination of site matrices and/or their corresponding sequence
matches (instances). This object inherits from Bio::Matrix::PSM::ProtMatrix, so
you can methods from that class. It may hold also an array of
Bio::Matrix::PSM::InstanceSite object, but you will have to retrieve these
through Bio::Matrix::PSM::ProtPsm-E<gt>instances method (see below). To some
extent this is an expanded ProtMatrix object, holding data from analysis that
also deal with sequence matches of a particular matrix.


=head2 DESIGN ISSUES

This does not make too much sense to me I am mixing PSM with PSM sequence
matches Though they are very closely related, I am not satisfied by the way
this is implemented here.  Heikki suggested different objects when one has
something like meme But does this mean we have to write a different objects for
mast, meme, transfac, theiresias, etc.?  To me the best way is to return
SiteMatrix object + arrray of InstanceSite objects and then mast will return
undef for SiteMatrix and transfac will return undef for InstanceSite. Probably
I cannot see some other design issues that might arise from such approach, but
it seems more straightforward.  Hilmar does not like this beacause it is an
exception from the general BioPerl rules. Should I leave this as an option?
Also the header rightfully belongs the driver object, and could be retrieved as
hashes.  I do not think it can be done any other way, unless we want to create
even one more object with very unclear content.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

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
the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - James Thompson

Email tex@biosysadmin.com


=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 SEE ALSO

ProtMatrix, meme, transfac, psiblast, InstanceSite

=head1 APPENDIX

=cut


# Let the code begin...
package Bio::Matrix::PSM::ProtPsm;
use Bio::Matrix::PSM::InstanceSite;
use strict;

use base qw(Bio::Matrix::PSM::ProtMatrix Bio::Matrix::PSM::PsmI Bio::Annotation::Collection);

@Bio::Matrix::PSM::Psm::HEADER = qw(e_val sites IC width);

=head2 new

 Title   : new
 Usage   : my $psm = Bio::Matrix::PSM::ProtPsm->new(
              -pS => [ '0', '33', '0', '16', '1', '12', '11', '25' ],
              -pF => [ '0', '0', '2', '0', '3', '0', '0', '0' ],
              -pT => [ '0', '8', '7', '10', '1', '2', '7', '8' ],
              -pN => [ '0', '0', '2', '13', '0', '36', '1', '4' ],
              -pK => [ '0', '5', '0', '13', '1', '15', '0', '2' ],
              -pY => [ '0', '0', '0', '0', '0', '0', '0', '0' ],
              -pE => [ '0', '41', '1', '12', '0', '0', '0', '15' ],
              -pV => [ '0', '3', '9', '0', '2', '0', '3', '1' ],
              -pQ => [ '0', '0', '0', '15', '0', '4', '0', '3' ],
              -pM => [ '100', '0', '66', '0', '2', '0', '0', '0' ],
              -pC => [ '0', '0', '0', '0', '0', '0', '0', '0' ],
              -pL => [ '0', '0', '8', '0', '25', '0', '4', '0' ],
              -pA => [ '0', '10', '1', '9', '2', '0', '22', '16' ],
              -pW => [ '0', '0', '0', '0', '0', '0', '0', '0' ],
              -pP => [ '0', '0', '0', '0', '3', '1', '45', '0' ],
              -pH => [ '0', '0', '0', '0', '0', '0', '1', '0' ],
              -pD => [ '0', '0', '1', '7', '2', '2', '0', '22' ],
              -pR => [ '0', '0', '0', '3', '0', '27', '0', '0' ],
              -pI => [ '0', '0', '3', '0', '59', '1', '2', '3' ],
              -pG => [ '0', '0', '0', '1', '0', '0', '4', '1' ],
              -IC => $ic,
              -sites => $istes,
              -width => $width,
              -e_val => $e_val, 
              -instances => $instances, 
           }

 Function: Creates a new Bio::Matrix::PSM::ProtPsm object
 Throws  :
 Example :
 Returns :  Bio::Matrix::PSM::Psm object
 Args    :  hash


=cut

sub new {
    my ($caller,@args) = @_;
    my $class = ref($caller) || $caller;
    my $self = $class->SUPER::new(@args);
    $self->{'_annotation'} = {};  #Init from Annotation::Collection
    $self->_typemap(Bio::Annotation::TypeManager->new()); #same
    ($self->{instances})=$self->_rearrange(['INSTANCES'], @args);
    return $self;
}


=head2 instances

 Title   : instances
 Usage   :   my @instances=@{$psm->instances};
 Function: Gets/sets the instances (Bio::Matrix::PSM::InstanceSite objects)
        associated with the Psm object
 Throws  :
 Example :
 Returns :  array reference (Bio::Matrix::PSM::InstanceSite objects)
 Args    :  array reference (Bio::Matrix::PSM::InstanceSite objects)

=cut

sub instances {
    my $self = shift;
    my $prev = $self->{instances};
    if (@_) { $self->{instances} = shift; }
    return $prev;
}


=head2 header

 Title   : header
 Usage   :  my %header=$psm->header;
        my $ic=$psm->header('IC');
 Function: Gets the general information, common for most files,
       dealing with PSM such as information content (IC), score
       (e-value, etc.), number of sites (sites) and width. This
       list may expand. The current list should be in
       @Bio::Matrix::PSM::Psm::HEADER. Returns an epty list if an
       argument is supplied that is not in
       @Bio::Matrix::PSM::meme::HEADER.
 Throws  :
 Example :
 Returns :  hash or string
 Args    :  string (IC, e_val...)

=cut

sub header {
    my $self = shift;
    return  if ($self->{end});
    my %header;
    if (@_) {my $key=shift; return $self->{$key}; }
    foreach my $key (@Bio::Matrix::PSM::ProtPsm::HEADER) {
	$header{$key}=$self->{$key};
    }
    return %header;
}


=head2 matrix

 Title   :  matrix
 Usage   :  my $matrix = $psm->matrix;
 Function:  Gets/sets the SiteMatrix related information
 Throws  :
 Example :
 Returns :  Bio::Matrix::PSM::SiteMatrix objects
 Args    :  Bio::Matrix::PSM::SiteMatrix objects

=cut


sub matrix {
   my $self = shift;

   if (@_) {
      my $matrix = shift;
      my @alphabet = $self->alphabet;

      foreach my $char (@alphabet) {
         $self->{"log$char"}  = $matrix->{"log$char"};
         $self->{"prob$char"} = $matrix->{"prob$char"};
      }
      $self->{IC}    = $matrix->IC;
      $self->{e_val} = $matrix->e_val;
      $self->{id}    = $matrix->id;
    }

    return $self;
}

1;
