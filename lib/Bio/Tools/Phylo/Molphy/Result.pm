#
# BioPerl module for Bio::Tools::Phylo::Molphy::Result
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

Bio::Tools::Phylo::Molphy::Result - container for data parsed from a ProtML run

=head1 SYNOPSIS

  # do not use this object directly, you will get it back as part of a 
  # Molphy parser
  use Bio::Tools::Phylo::Molphy;
  my $parser = Bio::Tools::Phylo::Molphy->new(-file => 'output.protml');
  while( my $r = $parser->next_result ) {
    # r is a Bio::Tools::Phylo::Molphy::Result object

    # print the model name
    print $r->model, "\n";

    # get the substitution matrix
    # this is a hash of 3letter aa codes -> 3letter aa codes representing
    # substitution rate
    my $smat = $r->substitution_matrix;
    print "Arg -> Gln substitution rate is %d\n", 
          $smat->{'Arg'}->{'Gln'}, "\n";

    # get the transition probablity matrix
    # this is a hash of 3letter aa codes -> 3letter aa codes representing
    # transition probabilty
    my $tmat = $r->transition_probability_matrix;
    print "Arg -> Gln transition probablity is %.2f\n", 
          $tmat->{'Arg'}->{'Gln'}, "\n";

    # get the frequency for each of the residues
    my $rfreqs = $r->residue_frequencies;

    foreach my $residue ( keys %{$rfreqs} ) {
       printf "residue %s  expected freq: %.2f observed freq: %.2f\n",
              $residue,$rfreqs->{$residue}->[0], $rfreqs->{$residue}->[1];
    }

    my @trees;
    while( my $t = $r->next_tree ) {
        push @trees, $t;
    }

    print "search space is ", $r->search_space, "\n",
          "1st tree score is ", $trees[0]->score, "\n";

    # writing to STDOUT, use -file => '>filename' to specify a file
    my $out = Bio::TreeIO->new(-format => "newick");
    $out->write_tree($trees[0]); # writing only the 1st tree
  }


=head1 DESCRIPTION

A container for data parsed from a ProtML run.


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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Phylo::Molphy::Result;
use strict;

# Object preamble - inherits from Bio::Root::Root



use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Phylo::Molphy::Result->new();
 Function: Builds a new Bio::Tools::Phylo::Molphy::Result object 
 Returns : Bio::Tools::Phylo::Molphy::Result
 Args    : 


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($trees,$smat,$freq,
      $model, $sspace,
      ) = $self->_rearrange([qw(TREES SUBSTITUTION_MATRIX
				FREQUENCIES
				MODEL SEARCH_SPACE)], @args);

  if( $trees ) {
      if(ref($trees) !~ /ARRAY/i ) { 
	  $self->warn("Must provide a valid array reference to initialize trees");
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
		   return;
	       }
	   }
	   $self->{'_substitution_matrix'} = $val;
       } else { 
	   $self->warn("Must be a valid hashref of hashrefs for substition_matrix");
	   return;
       }
   }
   return $self->{'_substitution_matrix'};
}

=head2 transition_probability_matrix

 Title   : transition_probability_matrix
 Usage   : my $matrixref = $molphy->transition_probablity_matrix();
 Function: Gets the observed transition probability matrix
 Returns : hash of hashes of aa/nt transition to each other aa/nt 
 Args    : Transition matrix type, typically
           '1PAM-1.0e05' or '1PAM-1.0e07'


=cut

sub transition_probability_matrix {
   my ($self,$type,$val) = @_;
   $type = '1PAM-1.0e7' unless defined $type;
   if(defined $val ) { 
       if( ref($val) =~ /HASH/ ) {
	   foreach my $v (values %{$val} ) {
	       if( ref($v) !~ /HASH/i ) { 
		   $self->warn("Must be a valid hashref of hashrefs for transition_probability_matrix");
		   return;
	       }
	   } 
	   $self->{'_TPM'}->{$type} = $val;
       } else { 
	   $self->warn("Must be a valid hashref of hashrefs for transition_probablity_matrix");
	   return;
       }
   }

   # fix this for nucml where there are 2 values (one is just a transformation
   # of the either, but how to represent?)
   return $self->{'_TPM'}->{$type};
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

sub residue_frequencies {
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
