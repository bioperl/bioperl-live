# $Id$
#
# BioPerl module for Bio::Search::Hit::GenericHit
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::Hit::GenericHit - A generic implementation of the Bio::Search::Hit::HitI interface

=head1 SYNOPSIS

  {
    use Bio::Search::Hit::GenericHit;
    my $hit = new Bio::Search::Hit::GenericHit(-algorithm => 'blastp');

  }

=head1 DESCRIPTION

This object handles the hit data from a Database Sequence Search such
as FASTA or BLAST.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich and Steve Chervitz

Email jason@bioperl.org
Email sac@bioperl.org

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::Hit::GenericHit;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::Search::Hit::HitI;

@ISA = qw(Bio::Root::Root Bio::Search::Hit::HitI );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Search::Hit::GenericHit();
 Function: Builds a new Bio::Search::Hit::GenericHit object 
 Returns : Bio::Search::Hit::GenericHit
 Args    : -name         => Name of Hit (required)
           -description  => Description (optional)
           -accession    => Accession number (optional)
           -length       => Length of the Hit (optional)
           -score        => Raw Score for the Hit (optional)
           -significance => Significance value for the Hit (optional)
           -algorithm    => Algorithm used (BLASTP, FASTX, etc...)
           -hsps         => Array ref of HSPs for this Hit. 

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($hsps, $name,$desc, $acc, $length,
      $score,$algo,$signif) = $self->_rearrange([qw(HSPS NAME DESCRIPTION
						   ACCESSION
						   LENGTH SCORE ALGORITHM 
						   SIGNIFICANCE)], @args);

  if( ! defined $name ) { 
      $self->throw("Must have defined a valid name for Hit");
  } else { 
      $self->name($name);
  }  
  defined $acc    && $self->accession($acc);
  defined $desc   && $self->description($desc);
  defined $length && $self->length($length);
  defined $algo   && $self->algorithm($algo);
  defined $signif && $self->significance($signif);
  defined $score  && $self->raw_score($score);
  
  $self->{'_iterator'} = 0;
  $self->{'_hsps'} = [];
  if( defined $hsps  ) {
      if( ref($hsps) !~ /array/i ) {
	  $self->warn("Did not specify a valid array ref for the param HSPS ($hsps)");
      } else {
	  while( @$hsps ) { 
	      $self->add_hsp(shift @$hsps );
	  }
      }
  } 
  return $self;
}

=head2 add_hsp

 Title   : add_hsp
 Usage   : $hit->add_hsp($hsp)
 Function: Add a HSP to the collection of HSPs for a Hit
 Returns : number of HSPs in the Hit
 Args    : Bio::Search::HSP::HSPI object


=cut

sub add_hsp {
   my ($self,$hsp) = @_;
   if( !defined $hsp || ! $hsp->isa('Bio::Search::HSP::HSPI') ) { 
       $self->warn("Must provide a valid Bio::Search::HSP::HSPI object to object: $self method: add_hsp");
       return undef;
   }
   push @{$self->{'_hsps'}}, $hsp;
   return scalar @{$self->{'_hsps'}};
}



=head2 Bio::Search::Hit::HitI methods

Implementation of Bio::Search::Hit::HitI methods

=head2 name

 Title   : name
 Usage   : $hit_name = $hit->name();
 Function: returns the name of the Hit sequence
 Returns : a scalar string
 Args    : [optional] scalar string to set the name

=cut

sub name {
    my ($self,$value) = @_;
    my $previous = $self->{'_name'};
    if( defined $value || ! defined $previous ) {
	$value = $previous = '' unless defined $value;
	$self->{'_name'} = $value;
    } 
    return $previous;
}

=head2 accession

 Title   : accession
 Usage   : $acc = $hit->accession();
 Function: Retrieve the accession (if available) for the hit
 Returns : a scalar string (empty string if not set)
 Args    : none

=cut

sub accession {
    my ($self,$value) = @_;
    my $previous = $self->{'_accession'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_accession'} = $value;
    } 
    return $previous;
}

=head2 description

 Title   : description
 Usage   : $desc = $hit->description();
 Function: Retrieve the description for the hit
 Returns : a scalar string
 Args    : [optional] scalar string to set the descrition

=cut

sub description {
    my ($self,$value) = @_;
    my $previous = $self->{'_description'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_description'} = $value;
    } 
    return $previous;
}

=head2 length

 Title   : length
 Usage   : my $len = $hit->length
 Function: Returns the length of the hit 
 Returns : integer
 Args    : [optional] integer to set the length

=cut

sub length {
    my ($self,$value) = @_;
    my $previous = $self->{'_length'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = 0 unless defined $value;
	$self->{'_length'} = $value;
    } 
    return $previous;
}


=head2 algorithm

 Title   : algorithm
 Usage   : $alg = $hit->algorithm();
 Function: Gets the algorithm specification that was used to obtain the hit
           For BLAST, the algorithm denotes what type of sequence was aligned 
           against what (BLASTN: dna-dna, BLASTP prt-prt, BLASTX translated 
           dna-prt, TBLASTN prt-translated dna, TBLASTX translated 
           dna-translated dna).
 Returns : a scalar string 
 Args    : [optional] scalar string to set the algorithm

=cut

sub algorithm {
    my ($self,$value) = @_;
    my $previous = $self->{'_algorithm'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_algorithm'} = $value;
    } 
    return $previous;
}

=head2 raw_score

 Title   : raw_score
 Usage   : $score = $hit->raw_score();
 Function: Gets the "raw score" generated by the algorithm.  What
           this score is exactly will vary from algorithm to algorithm,
           returning undef if unavailable.
 Returns : a scalar value
 Args    : [optional] scalar value to set the raw score

=cut

sub raw_score {
    my ($self,$value) = @_;
    my $previous = $self->{'_score'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_score'} = $value;
    } 
    return $previous;
}

=head2 significance

 Title   : significance
 Usage   : $significance = $hit->significance();
 Function: Used to obtain the E or P value of a hit, i.e. the probability that
           this particular hit was obtained purely by random chance.  If
           information is not available (nor calculatable from other
           information sources), return undef.
 Returns : a scalar value or undef if unavailable
 Args    : [optional] scalar value to set the significance

=cut

sub significance {
    my ($self,$value) = @_;
    my $previous = $self->{'_significance'};
    if( defined $value || ! defined $previous ) { 
	$value = $previous = '' unless defined $value;
	$self->{'_significance'} = $value;
    } 
    return $previous;
}

=head2 next_hsp

 Title    : next_hsp
 Usage    : while( $hsp = $obj->next_hsp()) { ... }
 Function : Returns the next available High Scoring Pair
 Example  : 
 Returns  : Bio::Search::HSP::HSPI object or null if finished
 Args     : none

=cut

sub next_hsp {
    my ($self) = @_;
    $self->{'_iterator'} = 0 unless defined $self->{'_iterator'};
    return undef if $self->{'_iterator'} > scalar @{$self->{'_hsps'}};
    return $self->{'_hsps'}->[$self->{'_iterator'}++];    
}

1;
