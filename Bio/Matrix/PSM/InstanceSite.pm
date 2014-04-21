
=head1 NAME

Bio::Matrix::PSM::InstanceSite - A PSM site occurance

=head1 SYNOPSIS

  use Bio::Matrix::PSM::InstanceSite;

  #You can get an InstanceSite object either from a file:

  my ($instances,$matrix)=$SomePSMFile->parse_next;

  #or from memory

  my %params=(seq=>'TATAAT',
    id=>"TATAbox1", accession=>'ENSG00000122304', mid=>'TB1',
    desc=>'TATA box, experimentally verified in PRM1 gene',
    -relpos=>-35, -anchor=>'CHR7', -start=>35000921, -end=>35000926);

  #Last 2 arguments are passed to create a Bio::LocatableSeq object
  #Anchor shows the coordinates system for the Bio::LocatableSeq object

=head1 DESCRIPTION

Abstract interface to PSM site occurrence (PSM sequence
match). InstanceSite objects may be used to describe a PSM (See
L<Bio::Matrix::PSM::SiteMatrix>) sequence matches.  The usual
characteristic of such a match is sequence coordinates, score,
sequence and sequence (gene) identifier- accession number or other id.

This object inherits from Bio::LocatableSeq (which defines the real
sequence) and might hold a SiteMatrix object, used to detect the CRE
(cis-regulatory element), or created from this CRE.

While the documentation states that the motif id and gene id
(accession) combination should be unique, this is not entirely true-
there might be more than one occurrence of the same cis-regulatory
element in the upstream region of the same gene.  Therefore relpos
would be the third element to create a really unique combination.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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

=head2 Description

Bio::Matrix::PSM::InstanceSiteI implementation

=head1 AUTHOR - Stefan Kirov

Email skirov@utk.edu


=head1 APPENDIX

=cut


# Let the code begin...
package Bio::Matrix::PSM::InstanceSite;
use strict;

use base qw(Bio::LocatableSeq Bio::Matrix::PSM::InstanceSiteI);

=head2 new

 Title   : new
 Usage   : my $isntance=Bio::Matrix::PSM::InstanceSite->new 
                         (-seq=>'TATAAT', -id=>"TATAbox1",
                          -accession_number='ENSG00000122304', -mid=>'TB1',
                          -desc=>'TATA box, experimentally verified in PRM1 gene',
                          -relpos=>-35, -anchor=>'CHR7', -start=>35000921, -end=>35000926, strand=>1)
 Function: Creates an InstanceSite object from memory.
 Throws  :
 Example :
 Returns : Bio::Matrix::PSM::InstanceSite object
 Args    : hash


=cut

sub new {
    my ($class, @args) = @_;
    my %args = @args; #Too many things to rearrange, and I am creating >1K such objects routinely, so this is a performance issue    
    $args{'-start'} ||= 1;
    my $end = $args{'-start'} + length($args{-seq}) -1;
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
    $self->seq($args{seq});
    $self->desc($args{desc});
    $self->{score}=$args{score};
    $self->{relpos}=$args{relpos};
    $self->{frame}=$args{frame};
    $self->{anchor}=$args{anchor};
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

=head2 anchor

 Title   : anchor
 Usage   : my $anchor=$instance->anchor;
 Function: Get/Set the anchor which shows what coordinate system start/end use
 Throws  :
 Example :
 Returns : string
 Args    : string

=cut

sub anchor {
    my $self = shift;
    my $prev = $self->{anchor};
    if (@_) { $self->{anchor} = shift; }
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


#Provided by LocatableSeq

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

=head2 annotation

 Title   : annotation
 Usage   : $ann = $seq->annotation or $seq->annotation($annotation)
 Function: Gets or sets the annotation
 Returns : L<Bio::AnnotationCollectionI> object
 Args    : None or L<Bio::AnnotationCollectionI> object

See L<Bio::AnnotationCollectionI> and L<Bio::Annotation::Collection>
for more information

=cut

sub annotation {
    my ($obj,$value) = @_;
    if( defined $value ) {
	$obj->throw("object of class ".ref($value)." does not implement ".
		    "Bio::AnnotationCollectionI. Too bad.")
	    unless $value->isa("Bio::AnnotationCollectionI");
	$obj->{'_annotation'} = $value;
    } elsif( ! defined $obj->{'_annotation'}) {
	$obj->{'_annotation'} = Bio::Annotation::Collection->new();
    }
    return $obj->{'_annotation'};
}

=head2 species

 Title   : species
 Usage   : $species = $seq->species() or $seq->species($species)
 Function: Gets or sets the species
 Returns : L<Bio::Species> object
 Args    : None or L<Bio::Species> object

See L<Bio::Species> for more information

=cut

sub species {
    my ($self, $species) = @_;
    if ($species) {
        $self->{'species'} = $species;
    } else {
        return $self->{'species'};
    }
}


=head2 frame

 Title   : frame
 Usage   : my $frane=$instance->frame;
 Function: Get/Set the frame of a DNA instance with respect to a protein motif used.
            Returns undef if the motif was not protein or the DB is protein.
 Throws  :
 Example :
 Returns : integer
 Args    : integer (0, 1, 2)

=cut

sub frame {
    my $self = shift;
    my $prev = $self->{frame};
    if (@_) { $self->{frame} = shift; $self->throw("This is not a legitimate frame") unless (grep(/$self->{frame}/,qw[0 1 2])); }
    return $prev;
}

1;
