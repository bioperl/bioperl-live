package Bio::Matrix::PSM::InstanceSite;
=head1 NAME

Bio::Matrix::PSM::InstanceSite

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

=head2 Description

Bio::Matrix::PSM::InstanceSiteI implementation

=head1 AUTHOR - Stefan Kirov

Email skirov@utk.edu


=head1 APPENDIX

=cut


# Let the code begin...
use Bio::Matrix::PSM::SiteMatrix;
use Bio::Root::Root;
use Bio::Matrix::PSM::InstanceSiteI;
use Bio::LocatableSeq;
use vars qw(@ISA);
use strict;

 @ISA=qw(Bio::LocatableSeq  Bio::Matrix::PSM::InstanceSiteI Bio::Root::Root);
 
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
	(my $class, my @args)=@_;
	my %args=@args; #Too many things to rearrange, and I am creating >1K such objects routinely, so this is a performance issue
	my $end=$args{-start}+length($args{-seq});
	if (!defined($args{-strand})) {
		$args{-strand}=1;
		@args=%args;
	}
my $self = $class->SUPER::new(@args,'-end',$end);
while( @args ) {
	(my $key = shift @args) =~ s/-//gi; #deletes all dashes (only dashes)!
	$args{$key} = shift @args;
}
#should throw exception if seq is null, for now just warn
if (($args{seq} eq '') || (!defined($args{seq}))) {
	$args{seq}="AGCT";
	warn "No sequence?!\n";
}
$self->{mid}=$args{mid};
my $minstance=$args{mid} . "@" . $args{accession};
$self->seq($args{seq});
$self->desc($args{desc});
$self->accession_number($args{accession});
$self->primary_id($minstance); #Since this is a unique key this is the place for it?
$self->{score}=$args{score};
$self->{relpos}=$args{relpos};
return $self;
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
    my $prev = $self->{mid};
    if (@_) { $self->{mid} = shift; }
    return $prev;
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
    my $prev = $self->{score};
    if (@_) { $self->{score} = shift; }
    return $prev;
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
    my $start = $self->{start};
    return $start;
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
    my $prev = $self->{minstance};
    if (@_) { $self->{minstance} = shift; }
    return $prev;
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
    my $prev = $self->{relpos};
    if (@_) { $self->{relpos} = shift; }
    return $prev;
}

1;
