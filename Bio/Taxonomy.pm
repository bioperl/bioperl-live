# $Id$
#
# BioPerl module for Bio::Taxonomy
#
# Cared for by Dan Kortschak
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Taxonomy - Conversion used bt the Taxonomy classes

=head1 SYNOPSIS

    use Bio::Taxonomy;

=head1 DESCRIPTION

Provides methods for converting classifications into taxonomic
structures.

=head1 CONTACT

Dan Kortschak email B<kortschak@rsbs.anu.edu.au>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# code begins...


package Bio::Taxonomy;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object
use Bio::Root::Root;

@ISA = qw(Bio::Root::Root);


=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Taxonomy();
 Function: Builds a new Bio::Taxonomy object
 Returns : Bio::Taxonomy
 Args    : -method  -> method used to decide classification
                       (none|trust|lookup)
           -ranks   -> what ranks are there

=cut


sub new {
   my ($class,@args) = @_;

   my $self = $class->SUPER::new(@args);

   $self->{'_method'}='none';
   $self->{'_ranks'}=[];
   $self->{'_rank_hash'}={};

   my ($method,$ranks,$order) = $self->_rearrange([qw(METHOD RANKS ORDER)], @args);

   if ($method) {
      $self->method($method);
   }

   if (defined $ranks &&
      (ref($ranks) eq "ARRAY") ) {
      $self->ranks(@$ranks);
   } else {
      # default ranks
      # I think these are in the right order, but not sure:
      # some parvorder|suborder and varietas|subspecies seem
      # to be at the same level - any taxonomists?
      # I don't expect that these will actually be used except as a way
      # to find what ranks there are in taxonomic use 
      $self->ranks(('root',
                    'superkingdom',
                    'kingdom',
                    'superphylum',
                    'phylum',
                    'subphylum',
                    'superclass',
                    'class',
                    'subclass',
                    'infraclass',
                    'superorder',
                    'order',
                    'suborder',
                    'parvorder',
                    'infraorder',
                    'superfamily',
                    'family',
                    'subfamily',
                    'tribe',
                    'subtribe',
                    'genus',
                    'subgenus',
                    'species group',
                    'species subgroup',
                    'species',
                    'subspecies',
                    'varietas',
                    'forma',
                    'no rank'));
   }

   return $self;
}


=head2 method

 Title   : method
 Usage   : $obj = taxonomy->method($method);
 Function: set or return the method used to decide classification
 Returns : $obj
 Args    : $obj

=cut


sub method {
   my ($self,$value) = @_;
   if (defined $value && $value=~/none|trust|lookup/) {
       $self->{'_method'} = $value;
   }
   return $self->{'_method'};
}


=head2 classify

 Title   : classify
 Usage   : @obj[][0-1] = taxonomy->classify($species);
 Function: return a ranked classification
 Returns : @obj of taxa and ranks as word pairs separated by "@"
 Args    : Bio::Species object

=cut


sub classify {
   my ($self,$value) = @_;
   my @ranks;

   if (! $value->isa('Bio::Species') ) {
      $self->throw("Trying to classify $value which is not a Bio::Species object");
   }

   my @classes=reverse($value->classification);

   if ($self->method eq 'none') {
      for (my $i=0; $i < @classes-2; $i++) {
         ($ranks[$i][0],$ranks[$i][1])=($classes[$i],'no rank');
      }
      push @ranks,[$classes[-2],'genus'];
      push @ranks,[$value->binomial,'species'];
   } elsif ($self->method eq 'trust') {
      if (scalar(@classes)==scalar($self->ranks)) {
         for (my $i=0; $i < @classes; $i++) {
            if ($self->rank_of_number($i) eq 'species') {
               push @ranks,[$value->binomial,$self->rank_of_number($i)];
            } else {
               push @ranks,[$classes[$i],$self->rank_of_number($i)];
            }
         }
      } else {
         $self->throw("Species object and taxonomy object cannot be reconciled");
      }
   } elsif ($self->method eq 'lookup') {
      # this will lookup a DB for the rank of a taxon name
      # I imagine that some kind of Bio::DB class will be need to
      # be given to the taxonomy object to act as an DB interface
      # (I'm not sure how useful this is though - if you have a DB of
      # taxonomy - why would you be doing things this way?)
      $self->throw("Not yet implemented");
   }

   return @ranks;
}


=head2 level_of_rank

 Title   : level_of_rank
 Usage   : $obj = taxonomy->level_of_rank($obj);
 Function: returns the level of a rank name
 Returns : $obj
 Args    : $obj

=cut


sub level_of {
   my ($self,$value) = @_;

   return $self->{'_rank_hash'}{$value};
}


=head2 rank_of_number

 Title   : rank_of_number
 Usage   : $obj = taxonomy->rank_of_number($obj);
 Function: returns the rank name of a rank level
 Returns : $obj
 Args    : $obj

=cut


sub rank_of_number {
   my ($self,$value) = @_;

   return ${$self->{'_ranks'}}[$value];
}


=head2 ranks

 Title   : ranks
 Usage   : @obj = taxonomy->ranks(@obj);
 Function: set or return all ranks
 Returns : @obj
 Args    : @obj

=cut


sub ranks {
   my ($self,@value) = @_;

   # currently this makes no uniqueness sanity check (this should be done)
   # I am think that adding a way of converting multiple 'no rank' ranks 
   # to unique 'no rank #' ranks so that the level of a 'no rank' is 
   # abstracted way from the user - I'm not sure of the vlaue of this

   if (defined @value) {
      $self->{'_ranks'}=\@value;
   }

   for (my $i=0; $i <= @{$self->{'_ranks'}}-1; $i++) {
      $self->{'_rank_hash'}{$self->{'_ranks'}[$i]}=$i unless $self->{'_ranks'}[$i] eq 'no rank';
   }
      
   return @{$self->{'_ranks'}};
}


1;
