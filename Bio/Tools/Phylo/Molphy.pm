# $Id$
#
# BioPerl module for Bio::Tools::Phylo::Molphy
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Phylo::Molphy - DESCRIPTION of Object

=head1 SYNOPSIS

    use Bio::Tools::Phylo::Molphy;
    my $parser = new Bio::Tools::Phylo::Molphy(-file => 'output.protml');
    while( my $result = $parser->next_result ) {

    }

=head1 DESCRIPTION

A parser for Molphy output (protml,dnaml)

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


package Bio::Tools::Phylo::Molphy;
use vars qw(@ISA);
use strict;

use Bio::Tools::Phylo::Molphy::Result;
use Bio::Root::Root;
use Bio::Root::IO;
use Bio::TreeIO;

@ISA = qw(Bio::Root::Root Bio::Root::IO );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::Phylo::Molphy();
 Function: Builds a new Bio::Tools::Phylo::Molphy object 
 Returns : Bio::Tools::Phylo::Molphy
 Args    : -fh/-file => $val, # for initing input, see Bio::Root::IO


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->_initialize_io(@args);

  return $self;
}

=head2 next_result

 Title   : next_result
 Usage   : my $r = $molphy->next_result
 Function: Get the next result set from parser data
 Returns : Bio::Tools::Phylo::Molphy::Result object
 Args    : none


=cut

sub next_result{
   my ($self) = @_;

   # A little statemachine for the parser here
   my ($state,$transition_ct,
       @transition_matrix, %transition_mat, @resloc,) = ( 0,0);
   my ( %subst_matrix, @treelines, @treedata, %frequencies);
   my ( $treenum,$possible_trees, $model);
   while( defined ( $_ = $self->_readline()) ) {
       if( /^Relative Substitution Rate Matrix/ ) {
	   if( %subst_matrix ) { 
	       $self->_pushback($_);
	       last;
	   }
	   $state = 0;
	   my ( @tempdata);
	   @resloc = ();
	   while( defined ($_ = $self->_readline) ) {
	       last if (/^\s+$/);
	       # remove leading/trailing spaces
	       s/^\s+//;
	       s/\s+$//;
	       my @data = split;
	       my $i = 0;
	       for my $l ( @data ) {
		   if( $l =~ /\D+/ ) { 
		       push @resloc, $l;
		   }
		   $i++;
	       }
	       push @tempdata, \@data;
	   }
	   my $i = 0;
	   for my $row ( @tempdata ) {
	       my $j = 0;
	       for my $col ( @$row ) {
		   if( $i == $j ) {
		       # empty string for diagonals
		       $subst_matrix{$resloc[$i]}->{$resloc[$j]} = '';
		   } else {
		       $subst_matrix{$resloc[$i]}->{$resloc[$j]} = $col;
		   }
		   $j++;
	       }
	       $i++;
	   }
       } elsif( /^Transition Probability Matrix/ ) {	   
	   if( /1\.0e7/ ) { 
	       $state = 1;
	       $transition_ct = 0;
	   } else { 
	       $state = 0;
	   }
       } elsif ( /Acid Frequencies/ ) {
	   $state = 0;
	   $self->_readline(); # skip the next line
	   while( defined( $_ = $self->_readline) ) {
	       unless( /^\s+/) {
		   $self->_pushback($_);
		   last;
	       }
	       s/^\s+//;
	       s/\s+$//;
	       my ($index,$res,$model,$data) = split;
	       $frequencies{$res} = [ $model,$data];
	   }
       } elsif( /^(\d+)\s*\/\s*(\d+)\s+(.+)\s+model/ ) {
	   my @save = ($1,$2,$3);
	   # finish processing the transition_matrix
	   my $i =0;
	   foreach my $row ( @transition_matrix ) {
	       my $j = 0;
	       foreach my $col ( @$row ) {
		   $transition_mat{$resloc[$i]}->{$resloc[$j]} = $col;
		   $j++;
	       }
	       $i++;
	   }
	   
	   if( defined $treenum ) { 	       
	       $self->_pushback($_);
	       last;
	   }
	   
	   $state = 2;	   
	   ($treenum,$possible_trees, $model) = @save;
	   $model =~ s/\s+/ /g;
       } elsif( $state == 1 ) {
	   next if( /^\s+$/ );
	   s/^\s+//;
	   s/\s+$//;
	   # because the matrix is split up into 2-10 column sets 
	   push @{$transition_matrix[$transition_ct++]}, split ;
	   $transition_ct = 0 if $transition_ct % 20 == 0;
       } elsif( $state == 2 ) {
	   if( s/^(\d+)\s+(\-?\d+(\.\d+)?)\s+// ) {
	       push @treedata, [ $1,$2];
	   }
	   # save this for the end so that we can 
	   # be efficient and only open one tree parser
	   push @treelines, $_;
       }
   }
   # waiting till the end to do this, is it better
   my @trees;
   if( @treelines ) {
       my $strdat = IO::String->new(join('',@treelines));
       my $treeio = new Bio::TreeIO(-fh => $strdat,
				    -format => 'newick');
       while( my $tree = $treeio->next_tree ) {
	   if( @treedata ) {
	       my $dat = shift @treedata;
	       # set the associated information
	       $tree->id($dat->[0]);
	       $tree->score($dat->[1]);
	   }
	   push @trees, $tree;
       }
   }

   my $result = new Bio::Tools::Phylo::Molphy::Result
       (-trees => \@trees,
	-substitution_matrix => \%subst_matrix,
	-transition_matrix   => \%transition_mat,
	-frequencies         => \%frequencies,
	-model               => $model,
	-search_space        => $possible_trees,
	);

}

1;
