# Result.pm,v 1.3 2002/06/20 18:50:39 amackey Exp
#
# BioPerl module for Bio::Tools::Phylo::PAML::Result
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich, Aaron Mackey
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Phylo::PAML::Result - A PAML result set object

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
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Jason Stajich, Aaron Mackey

Email jason@bioperl.org
Email amackey@virginia.edu

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Phylo::PAML::Result;
use vars qw(@ISA);
use strict;


use Bio::Root::Root;
use Bio::AnalysisResultI;
@ISA = qw(Bio::Root::Root Bio::AnalysisResultI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::Phylo::PAML::Result(%data);
 Function: Builds a new Bio::Tools::Phylo::PAML::Result object
 Returns : Bio::Tools::Phylo::PAML::Result
 Args    : -trees => array reference of L<Bio::Tree::TreeI> objects
           -MLmatrix => ML matrix
           .... MORE ARGUMENTS LISTED HERE BY AARON AND JASON 

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($trees,$mlmat,$seqs,$ngmatrix,
      $codonpos,$codonfreq,$version) = $self->_rearrange([qw(TREES MLMATRIX 
						    SEQS NGMATRIX
						    CODONPOS CODONFREQ
						    VERSION)], @args);
  $self->reset_seqs;
  if( $trees ) {
      if(ref($trees) !~ /ARRAY/i ) { 
	  $self->warn("Must have provided a valid array reference to initialize trees");
      } else { 
	  foreach my $t ( @$trees ) {
	      $self->add_tree($t);
	  }
      }
  }
  $self->{'_treeiterator'} = 0;

  if( $mlmat ) {
      if( ref($mlmat) !~ /ARRAY/i ) {
	  $self->warn("Must have provided a valid array reference to initialize MLmatrix");
      } else { 
	  $self->set_MLmatrix($mlmat);
      }
  } 
  if( $seqs ) { 
      if( ref($seqs) !~ /ARRAY/i ) {
	  $self->warn("Must have provided a valid array reference to initialize seqs");
      } else {
	  foreach my $s ( @$seqs ) {
	      $self->add_seq($s);
	  }
      }
  }
  if( $ngmatrix ) {
      if( ref($ngmatrix) !~ /ARRAY/i ) {
	  $self->warn("Must have provided a valid array reference to initialize NGmatrix");
      } else { 
	  $self->set_NGmatrix($ngmatrix);
      }
  } 
  
  if( $codonfreq ) {
      
  
  }

  if( $codonpos ) {
      if( ref($codonpos) !~ /ARRAY/i ) {
	  $self->warn("Must have provided a valid array reference to initialize codonpos");
      } else { 
	  $self->set_codon_pos_basefreq(@$codonpos);
      }
  }

  $self->version($version) if defined $version;

  return $self;
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


=head2 set_MLmatrix

 Title   : set_MLmatrix
 Usage   : $result->set_MLmatrix($mat)
 Function: Set the ML Matrix
 Returns : none
 Args    : Arrayref to MLmatrix (must be arrayref to 2D matrix whic is 
	   lower triangle pairwise)


=cut

sub set_MLmatrix{
   my ($self,$mat) = @_;
   return unless ( defined $mat );
   if( ref($mat) !~ /ARRAY/i ) {
       $self->warn("Did not provide a valid 2D Array reference for set_MLmatrix");
       return;
   }
   $self->{'_mlmatrix'} = $mat;
}

=head2 get_MLmatrix

 Title   : get_MLmatrix
 Usage   : my $mat = $result->get_MLmatrix()
 Function: Get the ML matrix
 Returns : 2D Array reference
 Args    : none


=cut

sub get_MLmatrix{
   my ($self,@args) = @_;
   return $self->{'_mlmatrix'};
}

=head2 set_NGmatrix

 Title   : set_NGmatrix
 Usage   : $result->set_NGmatrix($mat)
 Function: Set the Nei & Gojobori Matrix
 Returns : none
 Args    : Arrayref to NGmatrix (must be arrayref to 2D matrix whic is 
	   lower triangle pairwise)


=cut

sub set_NGmatrix{
   my ($self,$mat) = @_;
   return unless ( defined $mat );
   if( ref($mat) !~ /ARRAY/i ) {
       $self->warn("Did not provide a valid 2D Array reference for set_NGmatrix");
       return;
   }
   $self->{'_ngmatrix'} = $mat;
}

=head2 get_NGmatrix

 Title   : get_NGmatrix
 Usage   : my $mat = $result->get_NGmatrix()
 Function: Get the Nei & Gojobori matrix
 Returns : 2D Array reference
 Args    : none


=cut

sub get_NGmatrix{
   my ($self,@args) = @_;
   return $self->{'_ngmatrix'};
}


=head2 add_seq

 Title   : add_seq
 Usage   : $obj->add_seq($seq)
 Function: Add a Bio::PrimarySeq to the Result
 Returns : none
 Args    : Bio::PrimarySeqI
See also : L<Bio::PrimarySeqI>

=cut

sub add_seq{
   my ($self,$seq) = @_;
   if( $seq ) { 
       unless( $seq->isa("Bio::PrimarySeqI") ) {
	   $self->warn("Must provide a valid Bio::PrimarySeqI to add_seq");
	   return;
       }
       push @{$self->{'_seqs'}},$seq;
   }

}

=head2 reset_seqs

 Title   : reset_seqs
 Usage   : $result->reset_seqs
 Function: Reset the OTU seqs stored
 Returns : none
 Args    : none


=cut

sub reset_seqs{
   my ($self) = @_;
   $self->{'_seqs'} = [];
}

=head2 get_seqs

 Title   : get_seqs
 Usage   : my @otus = $result->get_seqs
 Function: Get the seqs Bio::PrimarySeq (OTU = Operational Taxonomic Unit)
 Returns : Array of Bio::PrimarySeq
 Args    : None
See also : L<Bio::PrimarySeq>

=cut

sub get_seqs{
   my ($self) = @_;
   return @{$self->{'_seqs'}};
}

=head2 set_codon_pos_basefreq

 Title   : set_codon_pos_basefreq
 Usage   : $result->set_codon_pos_basefreq(@freqs)
 Function: Set the codon position base frequencies
 Returns : none
 Args    : Array of length 3 where each slot has a hashref 
           keyed on DNA base


=cut

sub set_codon_pos_basefreq {
    my ($self,@codonpos) = @_;
    if( scalar @codonpos != 3 ) { 
	$self->warn("invalid array to set_codon_pos_basefreq, must be an array of length 3");
	return;
    }
    foreach my $pos ( @codonpos ) { 
	if( ref($pos) !~ /HASH/i ||
	    ! exists $pos->{'A'} ) { 
	    $self->warn("invalid array to set_codon_pos_basefreq, must be an array with hashreferences keyed on DNA bases, C,A,G,T");
	}
    }
    $self->{'_codonposbasefreq'} = [@codonpos];
}

=head2 get_codon_pos_basefreq

 Title   : get_codon_pos_basefreq
 Usage   : my @basepos = $result->get_codon_pos_basefreq;
 Function: Get the codon position base frequencies
 Returns : Array of length 3 (each codon position), each 
           slot is a hashref keyed on DNA bases, the values are
           the frequency of the base at that position for all sequences
 Args    : none
 Note    : The array starts at 0 so position '1' is in position '0' 
           of the array

=cut

sub get_codon_pos_basefreq{
   my ($self) = @_;
   return @{$self->{'_codonposbasefreq'}};
}

=head2 version

 Title   : version
 Usage   : $obj->version($newval)
 Function: Get/Set version
 Returns : value of version
 Args    : newvalue (optional)


=cut

sub version{
   my $self = shift;
   $self->{'_version'} = shift if @_;
   return $self->{'_version'};
}

1;
