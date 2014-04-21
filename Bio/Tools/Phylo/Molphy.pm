#
# BioPerl module for Bio::Tools::Phylo::Molphy
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Phylo::Molphy - parser for Molphy output

=head1 SYNOPSIS

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

A parser for Molphy output (protml,dnaml)

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


package Bio::Tools::Phylo::Molphy;
use strict;

use Bio::Tools::Phylo::Molphy::Result;
use Bio::TreeIO;
use IO::String;

use base qw(Bio::Root::Root Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Phylo::Molphy->new();
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
   my ($trans_type,$trans_amount);
   my $parsed = 0;
   while( defined ( $_ = $self->_readline()) ) {
       $parsed = 1;
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
	   if( /(1\.0e(5|7))\)\s+(\S+)/ ) {
	       $state = 1;
	       my $newtrans_type = "$3-$1";
	       $trans_amount = $1;
	       if( defined $trans_type ) {
		   # finish processing the transition_matrix
		   my $i =0;
		   foreach my $row ( @transition_matrix ) {
		       my $j = 0;
		       foreach my $col ( @$row ) {
			   $transition_mat{$trans_type}->{$resloc[$i]}->{$resloc[$j]} = $col;
			   $j++;
		       }
		       $i++;
		   }
	       }
	       $trans_type = $newtrans_type;
	       $transition_ct = 0;
	       @transition_matrix = ();
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
		   $transition_mat{$trans_type}->{$resloc[$i]}->{$resloc[$j]} = $col;
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
	   next if( /^\s+$/ || /^\s+Ala/);
	   s/^\s+//;
	   s/\s+$//;
	   if( $trans_type eq '1PAM-1.0e7' ) {
	       # because the matrix is split up into 2-10 column sets 
	       push @{$transition_matrix[$transition_ct++]}, split ;	   
	       $transition_ct = 0 if $transition_ct % 20 == 0;
	   } elsif( $trans_type eq '1PAM-1.0e5' ) {
	       # because the matrix is split up into 2-10 column sets 
	       my ($res,@row) = split;
	       next if $transition_ct >= 20; # skip last 
	       push @{$transition_matrix[$transition_ct++]}, @row;	   	       
	   }
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
       my $treeio = Bio::TreeIO->new(-fh => $strdat,
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
   return unless( $parsed );
   my $result = Bio::Tools::Phylo::Molphy::Result->new
       (-trees => \@trees,
	-substitution_matrix => \%subst_matrix,
	-frequencies         => \%frequencies,
	-model               => $model,
	-search_space        => $possible_trees,
	);
   while( my ($type,$mat) = each %transition_mat ) {
       $result->transition_probability_matrix( $type,$mat);
   }
   $result;
}

1;
