#
# BioPerl module for Bio::PopGen::Marker
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PopGen::Marker - A genetic marker which one uses to generate genotypes

=head1 SYNOPSIS

  my $name = $marker->name();            # marker name
  my $description = $marker->description(); # description
  my $type = $marker->type();            # coded type of the marker
  my $unique_id = $marker->unique_id;    # optional unique ID
  my @alleles = $marker->get_Alleles();  # the known alleles
  my %allele_freqs = $marker->get_Allele_Frequencies(); # keys are marker names
                                         # vals are frequencies
                                         # may change to handle multiple populations

=head1 DESCRIPTION

This object will not contain genotype information pertaining to an
individual, but rather population level statistics and descriptive
information about a marker.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 CONTRIBUTORS

Matthew Hahn, matthew.hahn-at-duke.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::PopGen::Marker;
use strict;

# Object preamble - inherits from Bio::Root::Root


use vars qw($UniqueCounter);

$UniqueCounter = 0;

use base qw(Bio::Root::Root Bio::PopGen::MarkerI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::PopGen::Marker->new();
 Function: Builds a new Bio::PopGen::Marker object 
 Returns : an instance of Bio::PopGen::Marker
 Args    : -name          => [string] marker name
           -description   => [string] marker description
           -type          => [string] marker type
           -unique_id     => [string/int] unique id
           -allele_freq   => [hash ref] allele frequencies 

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($name,$desc,$type,$uid,$af) = $self->_rearrange([qw(NAME
							  DESCRIPTION
							  TYPE
							  UNIQUE_ID
							  ALLELE_FREQ)],@args);
  $self->{'_allele_freqs'} = {};
  if( ! defined $uid ) {
      $uid = $UniqueCounter++;
  }
  if( defined $name) {
      $self->name($name);
  } else { 
      $self->throw("Must provide a name when initializing a Marker");
  }
  defined $desc && $self->description($desc);
  defined $type && $self->type($type);
  $self->unique_id($uid);
  if( defined $af) {
      if( ref($af) !~ /HASH/i ) {
	  $self->warn("Must provide valid Hash reference for allele_freq method");
      } else { 
	  foreach my $allele ( keys %$af ) {
	      $self->add_Allele_Frequency($allele, $af->{$allele});
	  }
      }
  }
  return $self;
}

=head2 name

 Title   : name
 Usage   : my $name = $marker->name();
 Function: Get the name of the marker
 Returns : string representing the name of the marker
 Args    : [optional] name


=cut

sub name{
    my $self = shift;

    return $self->{'_name'} = shift if @_;
    return $self->{'_name'};
}


=head2 description

 Title   : description
 Usage   : my $desc = $marker->description
 Function: Get the marker description free text
 Returns : string
 Args    : [optional] string


=cut

sub description{
    my $self = shift;

    return $self->{'_description'} = shift if @_;
    return $self->{'_description'};
}

=head2 type

 Title   : type
 Usage   : my $type = $marker->type;
 Function: Get coded string for marker type
 Returns : string
 Args    : [optional] string


=cut

sub type{
    my $self = shift;

    return $self->{'_type'} = shift if @_;
    return $self->{'_type'};
}


=head2 unique_id

 Title   : unique_id
 Usage   : my $id = $marker->unique_id;
 Function: Get the unique marker ID
 Returns : unique ID string
 Args    : [optional ] string


=cut

sub unique_id{
    my $self = shift;

    return $self->{'_uniqueid'} = shift if @_;
    return $self->{'_uniqueid'};
}


=head2 annotation

 Title   : annotation
 Usage   : my $annotation_collection = $marker->annotation;
 Function: Get/set a Bio::AnnotationCollectionI for this marker
 Returns : Bio::AnnotationCollectionI object
 Args    : [optional set] Bio::AnnotationCollectionI object

=cut

sub annotation{
   my ($self, $arg) = @_;
   return $self->{_annotation} unless $arg;
   $self->throw("Bio::AnnotationCollectionI required for argument") unless
       ref($arg) && $arg->isa('Bio::AnnotationCollectionI');
   return $self->{_annotation} = $arg;
}

=head2 get_Alleles

 Title   : get_Alleles
 Usage   : my @alleles = $marker->get_Alleles();
 Function: Get the available marker alleles
 Returns : Array of strings
 Args    : none

=cut

sub get_Alleles{
    my $self = shift;
    my (@numeric,@alpha);

    for ( keys %{$self->{'_allele_freqs'}} ) {
	if( /[^\d\.\-e]/ ) { push @alpha, $_ }
	else { push @numeric, $_ }
    }
    @numeric = sort { $b <=> $a } @numeric;
    @alpha   = sort { $b cmp $a } @alpha;
    return @numeric,@alpha;
}


=head2 get_Allele_Frequencies

 Title   : get_Allele_Frequencies
 Usage   : my %allele_freqs = $marker->get_Allele_Frequencies;
 Function: Get the alleles and their frequency (set relative to
           a given population - you may want to create different
           markers with the same name for different populations
           with this current implementation
 Returns : Associative array where keys are the names of the alleles
 Args    : none


=cut

sub get_Allele_Frequencies{
   return %{$_[0]->{'_allele_freqs'}};
}

=head2 add_Allele_Frequency

 Title   : add_Allele_Frequency
 Usage   : $marker->add_Allele_Frequency($allele,$freq)
 Function: Adds an allele frequency
 Returns : None
 Args    : $allele - allele name
           $freq   - frequency value


=cut

sub add_Allele_Frequency{
   my ($self,$allele,$freq) = @_;
   $self->{'_allele_freqs'}->{$allele} = $freq;
}

=head2 reset_alleles

 Title   : reset_alleles
 Usage   : $marker->reset_alleles();
 Function: Reset the alleles for a marker
 Returns : None
 Args    : None


=cut

sub reset_alleles{
   my ($self) = @_;
   $self->{'_allele_freqs'} = {};
}

=head2 marker_coverage

 Title   : marker_coverage
 Usage   : $marker->marker_coverage();
 Function: Get marker coverage, that is, the number of 
           individuals where the marker is present 
           excluding missing or ambiguous alleles
 Returns : integer, representing marker coverage
 Args    : 


=cut

sub marker_coverage{
    my ($self) = @_;
 
    return $self->{_marker_coverage};
}

1;
