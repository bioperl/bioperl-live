# $Id$
#
# BioPerl module for Bio::Tools::Phylo::Molphy::Result
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Phylo::Molphy::Result - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

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

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Phylo::Molphy::Result;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;


@ISA = qw(Bio::Root::Root );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::Phylo::Molphy::Result();
 Function: Builds a new Bio::Tools::Phylo::Molphy::Result object 
 Returns : Bio::Tools::Phylo::Molphy::Result
 Args    : 


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($trees,
      $smat,$tmat,$freq,
      $model, $sspace,
      ) = $self->_rearrange([qw(TREES SUBSTITUTION_MATRIX
				TRANSITION_MATRIX FREQUENCIES
				MODEL SEARCH_SPACE)], @args);

  if( $trees ) {
      if(ref($trees) !~ /ARRAY/i ) { 
	  $self->warn("Must have provided a valid array reference to initialize trees");
      } else {
	  foreach my $t ( @$trees ) {
	      $self->add_tree($t);
	  }
      }
  }
  # initialize things through object methods to be a good 
  # little OO programmer
  if( ref($smat) =~ /HASH/i ) {
      $self->substitution_matrix($smat);
  }
  if( ref($tmat) =~ /HASH/i ) { 
      $self->transition_probability_matrix($tmat);
  }
  if( ref($freq) =~ /HASH/i ) {
      $self->residue_frequencies($freq);
  }
  
  $model && $self->model($model); 
  $sspace && $self->search_space($sspace);
  $self->{'_treeiterator'} = 0;

  return $self;
}

=head2 model

 Title   : model
 Usage   : $obj->model($newval)
 Function: 
 Returns : value of model
 Args    : newvalue (optional)


=cut

sub model{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'model'} = $value;
    }
    return $self->{'model'};

}

=head2 substitution_matrix

 Title   : substitution_matrix
 Usage   : my $smat = $result->subsitution_matrix;
 Function: Get the relative substitution matrix calculated in the ML procedure
 Returns : reference to hash of hashes where key is the aa/nt name and value
           is another hash ref which contains keys for all the aa/nt 
           possibilities
 Args    : none


=cut

sub substitution_matrix{
   my ($self,$val) = @_;
   if(defined $val ) { 
       if( ref($val) =~ /HASH/ ) {
	   foreach my $v (values %{$val} ) {
	       if( ref($v) !~ /HASH/i ) { 
		   $self->warn("Must be a valid hashref of hashrefs for substition_matrix");
		   return undef;
	       }
	   }
	   $self->{'_substitution_matrix'} = $val;
       } else { 
	   $self->warn("Must be a valid hashref of hashrefs for substition_matrix");
	   return undef;
       }
   }
   return $self->{'_substitution_matrix'};
}

=head2 transition_probability_matrix

 Title   : transition_probability_matrix
 Usage   : my $matrixref = $molphy->transition_probablity_matrix();
 Function: Gets the observed transition probability matrix
 Returns : hash of hashes of aa/nt transition to each other aa/nt 
 Args    : none


=cut

sub transition_probability_matrix{
   my ($self,$val) = @_;
   if(defined $val ) { 
       if( ref($val) =~ /HASH/ ) {
	   foreach my $v (values %{$val} ) {
	       if( ref($v) !~ /HASH/i ) { 
		   $self->warn("Must be a valid hashref of hashrefs for transition_probability_matrix");
		   return undef;
	       }
	   } 
	   $self->{'_TPM'} = $val;
       } else { 
	   $self->warn("Must be a valid hashref of hashrefs for transition_probablity_matrix");
	   return undef;
       }
   }

   # fix this for nucml where there are 2 values (one is just a transformation
   # of the either, but how to represent?)
   return $self->{'_TPM'};
}

=head2 residue_frequencies

 Title   : residue_frequencies
 Usage   : my %data = $molphy->residue_frequencies()
 Function: Get the modeled and expected frequencies for
           each of the residues in the sequence
 Returns : hash of either aa (protml) or nt (nucml) frequencies
           each key will point to an array reference where
           1st slot is model's expected frequency
           2nd slot is observed frequency in the data
           $hash{'A'}->[0] = 
 Args    : none


=cut

#'

sub residue_frequencies{
   my ($self,$val) = @_;
   if(defined $val ) { 
       if( ref($val) =~ /HASH/ ) {
	   $self->{'_residue_frequencies'} = $val;
       } else { 
	   $self->warn("Must be a valid hashref of hashrefs for residue_frequencies");
       }
   }
   return %{$self->{'_residue_frequencies'}};
}

=head2 next_tree

 Title   : next_tree
 Usage   : my $tree = $factory->next_tree;
 Function: Get the next tree from the factory
 Returns : L<Bio::Tree::TreeI>
 Args    : none

=cut

sub next_tree{
   my ($self,@args) = @_;
   return $self->{'_trees'}->[$self->{'_treeiterator'}++] || undef;
}

=head2 rewind_tree

 Title   : rewind_tree_iterator
 Usage   : $result->rewind_tree()
 Function: Rewinds the tree iterator so that next_tree can be 
           called again from the beginning
 Returns : none
 Args    : none

=cut

sub rewind_tree_iterator {
    shift->{'_treeiterator'} = 0;
}

=head2 add_tree

 Title   : add_tree
 Usage   : $result->add_tree($tree);
 Function: Adds a tree 
 Returns : integer which is the number of trees stored
 Args    : L<Bio::Tree::TreeI>

=cut

sub add_tree{
   my ($self,$tree) = @_;
   if( $tree && ref($tree) && $tree->isa('Bio::Tree::TreeI') ) {
       push @{$self->{'_trees'}},$tree;
   }
   return scalar @{$self->{'_trees'}};
}

=head2 search_space

 Title   : search_space
 Usage   : $obj->search_space($newval)
 Function: 
 Returns : value of search_space
 Args    : newvalue (optional)


=cut

sub search_space{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'search_space'} = $value;
    }
    return $self->{'search_space'};
}

1;
