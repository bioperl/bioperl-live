# $Id$
#
# BioPerl module for Bio::Tree::AlleleNode
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::AlleleNode - DESCRIPTION of Object

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


package Bio::Tree::AlleleNode;
use vars qw(@ISA);
use strict;

use Bio::Tree::Node;

@ISA = qw(Bio::Tree::Node );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tree::AlleleNode();
 Function: Builds a new Bio::Tree::AlleleNode object 
 Returns : Bio::Tree::AlleleNode
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  my ($alleles) = $self->_rearrange([qw(ALLELES)], @args);
  $self->{'_data'} = {};
  if( defined $alleles ) { 
      if( ref($alleles) !~ /HASH/i ) {
	  $self->warn("Must specify a valid HASH reference for the -alleles value...Ignoring initializing input");
	  
      } else { 
	  foreach my $mkr ( keys %{$alleles} ) { 
	      $self->add_alleles($mkr,@{$alleles->{$mkr}});
	  }
      }
  }
  return $self;
}

=head2 add_alleles

 Title   : add_alleles
 Usage   : $node->add_alleles($mkr,@alleles);
 Function: Adds allele(s) for $mkr, @alleles can be a single or 
           multiple alleles. If the same marker is added more than one, the 
           previous value will be overwritten with a warning.
 Returns : none
 Args    : $marker => marker name
           @alleles => alleles for the marker


=cut

sub add_alleles{
   my ($self,$marker,@alleles) = @_;
   if( ! defined $marker || $marker eq '' ) { 
       $self->warn("must specify a valid marker name for add_alleles");
       return;
   } 
   if( $self->{'_data'}->{$marker} ) { 
       $self->warn("Overwriting value of $marker");
   }
   $self->{'_data'}->{$marker} = []; # reset the array ref
   foreach my $a ( sort @alleles ) { 
       next if ! defined $a; # skip undef alleles 
       push @{$self->{'_data'}->{$marker}},$a;
   }
}

=head2 get_alleles

 Title   : get_alleles
 Usage   : my @alleles = $node->get_alleles($marker);
 Function: Return the alleles for a marker $marker
 Returns : Array of Alleles for a marker or empty array
 Args    : $marker name

=cut

sub get_alleles{
   my ($self,$marker) = @_;
   if( defined $self->{'_data'}->{$marker} ) { 
       return @{$self->{'_data'}->{$marker}};
   }
   return ();
}

=head2 get_marker_names

 Title   : get_marker_names
 Usage   : my @names =$node->get_marker_names();
 Function: Return the names of the markers that have been added to this node
 Returns : List of Marker Names
 Args    : none

=cut

sub get_marker_names{
   my ($self) = @_;
   return keys %{$self->{'_data'}};
}


1;
