#---------------------------------------------------------

#ISA SiteMatrix, HAS InstanceSite

=head1 NAME

Bio::Matrix::PSM::Psm - handle combination of site matricies

=head1 SYNOPSIS

  use Bio::Matrix::PSM::IO;

  #To get a Psm object from a file use the Psm parser:
  my $psmIO =  Bio::Matrix::PSM::IO->new(-format=>'meme', -file=>$file);

  # Now go through all entities in the file with next_psm, which
  # returns a Psm object see Bio::Matrix::PSM::IO for detailed
  # documentation (matrix predictions or matrix sequence matches or
  # both):

  while (my $psm=$psmIO->next_psm) {
    my %psm_header=$psm->header;
    my $ic=$psm_header{IC};
    my $sites=$psm_header{sites};
    my $width=$psm_header{width};
    my $score=$psm_header{e_val};
    my $IUPAC=$psm->IUPAC;
    my $instances=$psm->instances;
    foreach my $instance (@{$instances}) {
      my $id=$instance->primary_id;
      #Do something with the id
    }
  }

 #or create from memmory:
  my $psm= Bio::Matrix::PSM::Psm->new( -pA=>\@pA,-pC=>\@pC,-pG=>\@pG,-pT=>\@pT,
       -id=>$id,
       -instances=>$instances, -e_val=>$e_val,
       -IC=>$ic, -width=>$width, -sites=>$sites)

  # where pA through pG are the respective frequencies of the matrix (see also
  # Bio::Matrix::PSM::SiteMatrix), and everything else is self-explenatory, 
  # except for -instances (reference to an array of 
  #  Bio::Matrix::PSM::InstanceSite objects) which is documented bellow.

=head1 DESCRIPTION

To handle a combination of site matrices and/or their corresponding
sequence matches (instances). This object inherits from
Bio::Matrix::PSM::SiteMatrix, so you can use the respective
methods. It may hold also an array of Bio::Matrix::PSM::InstanceSite
object, but you will have to retrieve these through
Bio::Matrix::PSM::Psm-E<gt>instances method (see below). To some extent
this is an expanded SiteMatrix object, holding data from analysis that
also deal with sequence matches of a particular matrix.


=head2 DESIGN ISSUES

This does not make too much sense to me I am mixing PSM with PSM
sequence matches Though they are very closely related, I am not
satisfied by the way this is implemented here.  Heikki suggested
different objects when one has something like meme But does this mean
we have to write a different objects for mast, meme, transfac,
theiresias, etc.?  To me the best way is to return SiteMatrix object +
arrray of InstanceSite objects and then mast will return undef for
SiteMatrix and transfac will return undef for InstanceSite. Probably I
cannot see some other design issues that might arise from such
approach, but it seems more straightforward.  Hilmar does not like
this beacause it is an exception from the general BioPerl rules Should
I leave this as an option?  Also the header rightfully belongs the
driver object, and could be retrieved as hashes.  I do not think it
can be done any other way, unless we want to create even one more
object with very unclear content.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Stefan Kirov

Email skirov@utk.edu


=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 SEE ALSO

SiteMatrix, meme, transfac, InstanceSite

=head1 APPENDIX

=cut


# Let the code begin...
package Bio::Matrix::PSM::Psm;
use Bio::Matrix::PSM::InstanceSite;
use strict;

use base qw(Bio::Matrix::PSM::SiteMatrix Bio::Matrix::PSM::PsmI Bio::Annotation::Collection);

@Bio::Matrix::PSM::Psm::HEADER = qw(e_val sites IC width);

=head2 new

 Title   : new
 Usage   : my $psm= Bio::Matrix::PSM::Psm->new( -pA=>\@pA,-pC=>\@pC,
					       -pG=>\@pG,-pT=>\@pT,-id=>$id,
					       -instances=>$instances, 
					       -e_val=>$e_val,
					       -IC=>$ic, -width=>$width, 
					       -sites=>$sites)
 Function: Creates a new Bio::Matrix::PSM::Psm object
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
           @Bio::Matrix::PSM::Psm::HEADER. Returns undef if an
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
    foreach my $key (@Bio::Matrix::PSM::Psm::HEADER) {
	$header{$key}=$self->{$key};
    }
    return %header;
}


=head2 matrix

 Title   :  matrix
 Usage   :  my $matrix=$psm->matrix;
 Function:  Gets/sets the SiteMatrix related information
 Throws  :
 Example :
 Returns :  Bio::Matrix::PSM::SiteMatrix objects
 Args    :  Bio::Matrix::PSM::SiteMatrix objects

=cut


sub matrix {
    my $self = shift;
    my $prev = Bio::Matrix::PSM::SiteMatrix->new(-pA=>$self->{probA}, 
						-pC=>$self->{probC},
						-pG=>$self->{probG},
						-pT=>$self->{probT},
						-lA=>$self->{logA},
						-lC=>$self->{logC},
						-lG=>$self->{logG},
						-lT=>$self->{logT},
						-IC=>$self->{IC},
						-e_val=>$self->{e_val},
						-id=>$self->{id});
    if (@_) {
	my $matrix=shift;
	$self->{IC} = $matrix->IC;
	$self->{probA}=$matrix->{probA};
	$self->{probC}=$matrix->{probC};
	$self->{probG}=$matrix->{probG};
	$self->{probT}=$matrix->{probT};
	$self->{e_val}=$matrix->e_val;
	$self->{id}=$matrix->id;
    }
    return $prev;
}
 
1;

