# $Id$
#
# BioPerl module for Bio::Tree::Statistics
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::Statistics - DESCRIPTION of Object

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

=head1 CONTRIBUTORS

Matt Hahn <matthew.hahn@duke.duke>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tree::Statistics;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;

@ISA = qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tree::Statistics();
 Function: Builds a new Bio::Tree::Statistics object 
 Returns : Bio::Tree::Statistics
 Args    :


=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
}


=head2 fu_and_li_D

 Title   : fu_and_li_D
 Usage   : my $D = $statistics->fu_an_li_D($tree,$nummut);
 Function:
           For this we assume that the tree is made up of
           Bio::Tree::AlleleNode\'s which contain markers and alleles
           each marker is a 'mutation' 
 Returns : Fu and Li\'s D statistic for this Tree
 Args    : $tree - Bio::Tree::TreeI which contains Bio::Tree::AlleleNodes

=cut

sub fu_and_li_D{
   my ($self,$tree) = @_;
   
   # for this we assume that the tree is made up of
   # allele nodes which contain markers and alleles
   # each marker is a 'mutation' 
   my @nodes = $tree->get_nodes();
   my $muttotal =0;
   my $tipmutcount = 0;
   my $sampsize = 0;
   foreach my $n ( @nodes ) {
       if ($n->is_Leaf() ) {
	   $sampsize++;
	   $tipmutcount += $n->get_marker_names();
       }
       $muttotal += $n->get_marker_names();
   }

   if( $muttotal <= 0 ) { 
       $self->warn("mutation total was not > 0, cannot calculate a Fu and Li D");
       return 0;
   }
   my $a = 0;
   for(my $k= 1; $k < $sampsize; $k++ ) {
        $a += ( 1 / $k );
    }
   
   my $b = 0;
    for(my $k= 1; $k < $sampsize; $k++ ) {
        $b += ( 1 / $k**2 );
    }
 
    my $c = 2 * ( ( ( $sampsize * $a ) - (2 * ( $sampsize -1 ))) /
                  ( ( $sampsize - 1) * ( $sampsize - 2 ) ) );
 
    my $v = 1 + ( ( $a**2 / ( $b + $a**2 ) ) * ( $c - ( ( $sampsize + 1) /
                                                        ( $sampsize - 1) ) ));
 
    my $u = $a - 1 - $v;
    my $D = ( $muttotal - (  $a * $tipmutcount) ) /
            ( sqrt ( ($u * $muttotal) + ( $v * $muttotal**2) ) );
 
    return $D;
}


1;
