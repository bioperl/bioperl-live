package Bio::Matrix::PSM::InstanceSiteI;

=head1 NAME

Bio::Matrix::PSM::InstanceSiteI - InstanceSite interface, holds an instance of a PSM

=head1 SYNOPSIS

use Bio::Matrix::PSM::InstanceSite;
#You can get an InstanceSite object either from a file:
  my ($instances,$matrix)=$SomePSMFile->parse_next;
#or from memory
  my %params=(seq=>'TATAAT',id=>"TATAbox1", accession='ENSG00000122304', mid=>'TB1',
              desc=>'TATA box, experimentally verified in PRM1 gene',relpos=>-35);

=head1 DESCRIPTION

Abstract interface to PSM site occurrence (PSM sequence match). InstanceSite objects
may be used to describe a PSM (See Bio::Matrix::PSM::SiteMatrix) sequence matches.
The usual characteristic of such a match is sequence coordinates, score, sequence and
sequence (gene) identifier- accession number or other id. This object inherits from
Bio::LocatableSeq (which defines the real sequence) and might hold a SiteMatrix object,
used to detect the CRE (cis-regulatory element), or created from this CRE.
While the documentation states that the motif id and gene id (accession) combination
should be unique, this is not entirely true- there might be more than one occurrence
of the same cis-regulatory element in the upstream region of the same gene.
Therefore relpos would be the third element to create a really unique combination.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Stefan Kirov

Email skirov@utk.edu

=head1 APPENDIX

=head1 SEE ALSO

Bio::Matrix::PSM::SiteMatrix, Bio::Matrix::PSM::Psm, Bio::Matrix::PSM::IO
=cut


# Let the code begin...
use Bio::LocatableSeq;
use Bio::Root::RootI;
use vars qw(@ISA);
use strict;

 @ISA=qw(Bio::Root::RootI  Bio::LocatableSeqI);
 
=head2 new

 Title   : new
 Usage   : my $isntance=new Bio::Matrix::PSM::InstanceSite (-seq=>'TATAAT', -id=>"TATAbox1",
                          -accession_numbaer='ENSG00000122304', -mid=>'TB1',
                          -desc=>'TATA box, experimentally verified in PRM1 gene',-relpos=>-35)
 Function: Creates an InstanceSite object from memory.
 Throws  :
 Example :
 Returns : Bio::Matrix::PSM::InstanceSite object
 Args    : hash


=cut

 sub new {
     my $self = shift;
    $self->throw_not_implemented();
 }

=head2 mid

 Title   : mid
 Usage   : my $mid=$instance->mid;
 Function: Get/Set the motif id
 Throws  :
 Example :
 Returns : scalar
 Args    : scalar


=cut
sub mid {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 score

 Title   : score
 Usage   : my $score=$instance->score;
 Function: Get/Set the score (mismatches) between the instance and the attached (or
            initial) PSM
 Throws  :
 Example :
 Returns : real number
 Args    : real number
 
=cut
 
sub score {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 start

 Title   : start
 Usage   : my $start=$instance->start;
 Function: Get/Set the position of the instance on the sequence used
 Throws  :
 Example :
 Returns : integer
 Args    : integer

=cut

sub start {
    my $self = shift;
    $self->throw_not_implemented();
}

=head2 relpos

 Title   : relpos
 Usage   : my $seqpos=$instance->relpos;
 Function: Get/Set the relative position of the instance with respect to the transcription start
            site (if known). Can and usually is negative.
 Throws  :
 Example :
 Returns : integer
 Args    : integer

=cut

sub relpos {
    my $self = shift;
    $self->throw_not_implemented();
}


=head2 minstance

 Title   : minstance
 Usage   : my $minstance=$misntance->score;
 Function: Get/Set the unique identifier- sequence id/motif id, for example PRM1_TATAbox.
          Not necessarily human readable.
 Throws  :
 Example :
 Returns : string
 Args    : string

=cut

sub minstance {
    my $self = shift;
    $self->throw_not_implemented();
}


1;
