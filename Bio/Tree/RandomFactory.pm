# $Id$
#
# BioPerl module for Bio::Tree::RandomFactory
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tree::RandomFactory - TreeFactory for generating Random Trees

=head1 SYNOPSIS

use Bio::Tree::RandomFactory
my $factory = new Bio::Tree::RandomFactory( -samples => \@taxonnames,
					    -maxcount => 10);
					     
# or for anonymous samples
    
my $factory = new Bio::Tree::RandomFactory( -sample_size => 6, 
					    -maxcount = 50);

=head1 DESCRIPTION

Builds a random tree every time next_tree is called.

This algorithm is based on the make_tree algorithm from Richard Hudson 199? 
XXXX.

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


package Bio::Tree::RandomFactory;
use vars qw(@ISA);
use strict;

use Bio::Factory::TreeFactoryI;
use Bio::Root::Root;
use Bio::TreeIO::TreeEventBuilder;
use Bio::Tree::AlleleNode;

@ISA = qw(Bio::Root::Root Bio::Factory::TreeFactoryI );

=head2 new

 Title   : new
 Usage   : my $factory = new Bio::Tree::RandomFactory(-samples => \@samples,
						      -maxcount=> $N);
 Function: Initializes a Bio::Tree::RandomFactory object
 Returns : Bio::Tree::RandomFactory
 Args    :


=cut

sub new{
   my ($class,@args) = @_;
   my $self = $class->SUPER::new(@args);
   
   $self->{'_eventbuilder'} = new Bio::TreeIO::TreeEventBuilder();
   $self->{'_treecounter'} = 0;
   $self->{'_maxcount'} = 0;
   my ($maxcount, $samps,$samplesize ) = $self->_rearrange([qw(MAXCOUNT
							       SAMPLES
							       SAMPLE_SIZE)],
							   @args);
   my @samples;
   
   if( ! defined $samps ) { 
       if( ! defined $samplesize || $samplesize <= 0 ) { 
	   $self->throw("Must specify a valid samplesize if parameter -SAMPLE is not specified");
       }
       foreach ( 1..$samplesize ) { push @samples, "Samp$_"; }      
   } else { 
       if( ref($samps) =~ /ARRAY/i ) { 
	   $self->throw("Must specify a valid ARRAY reference to the parameter -SAMPLES, did you forget a leading '\\'?");
       }
       @samples = @$samps;
   }
   
   $self->samples(\@samples);
   $self->sample_size(scalar @samples);
   if( defined $maxcount ) { 
       $self->maxcount($maxcount);
   }
   return $self;
}

=head2 next_tree

 Title   : next_tree
 Usage   : my $tree = $factory->next_tree
 Function: Returns a random tree based on the initialized number of nodes
           NOTE: if maxcount is not specified on initialization or
                 set to a valid integer, subsequent calls to next_tree will 
                 continue to return random trees and never return undef
 Returns : Bio::Tree::TreeI object
 Args    : none

=cut

sub next_tree{
   my ($self) = @_;
   my $size = $self->sample_size;
   my $i;
   # adopted from Hudson, 19??
   my @nodes;
   my $start = 2 * $size; 
   my $bl = 0;
   my @list;
   for($i=0;$i<$start; $i++) { 
       $nodes[$i] = new Bio::Tree::AlleleNode(-id => "node$i");
       $list[$i] = $nodes[$i];
   }

   for($i=$start;$i > 1;$i-- ) {
       $bl += -2.0 * log ( 1.0000 - rand()) / ( $i * ($i - 1) );
       $nodes[2*$size - $i]->branch_length(sprintf("%.5f",$bl));       
   }
   for( $i= $size; $i > 1; $i--) {
       my $pick = $self->random($i);
       my $node1 = $nodes[$pick];
       $nodes[2*$size - $i]->add_Descendent($node1);
       $nodes[$pick] = $nodes[$i-1];
       $pick = $self->random($i-1);       
       my $node2 = $nodes[$pick];
       $nodes[2*$size - $i]->add_Descendent($node2);       
       $nodes[$pick] = $nodes[2*$size - $i];
   }
   my $tree = new Bio::Tree::Tree(-root => $nodes[0]);
   return $tree;
}

=head2 maxcount

 Title   : maxcount
 Usage   : $obj->maxcount($newval)
 Function: 
 Example : 
 Returns : value of maxcount
 Args    : newvalue (optional)


=cut

sub maxcount{
   my ($self,$value) = @_;
   if( defined $value) {
       if( $value =~ /^(\d+)/ ) { 
	   $self->{'maxcount'} = $1;
       } else { 
	   $self->warn("Must specify a valid Positive integer to maxcount");
	   $self->{'maxcount'} = 0;
       }
  }
   return $self->{'_maxcount'};
}

=head2 samples

 Title   : samples
 Usage   : $obj->samples($newval)
 Function: 
 Example : 
 Returns : value of samples
 Args    : newvalue (optional)


=cut

sub samples{
   my ($self,$value) = @_;
   if( defined $value) {
       if( ref($value) !~ /ARRAY/i ) { 
	   $self->warn("Must specify a valid array ref to the method 'samples'");
	   $value = [];
       } 
      $self->{'samples'} = $value;
    }
    return $self->{'samples'};

}

=head2 sample_size

 Title   : sample_size
 Usage   : $obj->sample_size($newval)
 Function: 
 Example : 
 Returns : value of sample_size
 Args    : newvalue (optional)


=cut

sub sample_size{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'sample_size'} = $value;
    }
    return $self->{'sample_size'};

}

=head2 attach_EventHandler

 Title   : attach_EventHandler
 Usage   : $parser->attatch_EventHandler($handler)
 Function: Adds an event handler to listen for events
 Returns : none
 Args    : Bio::Event::EventHandlerI

=cut

sub attach_EventHandler{
    my ($self,$handler) = @_;
    return if( ! $handler );
    if( ! $handler->isa('Bio::Event::EventHandlerI') ) {
	$self->warn("Ignoring request to attatch handler ".ref($handler). ' because it is not a Bio::Event::EventHandlerI');
    }
    $self->{'_handler'} = $handler;
    return;
}

=head2 _eventHandler

 Title   : _eventHandler
 Usage   : private
 Function: Get the EventHandler
 Returns : Bio::Event::EventHandlerI
 Args    : none


=cut

sub _eventHandler{
   my ($self) = @_;

   return $self->{'_handler'};
}

=head2 random

 Title   : random
 Usage   : my $rint = $node->random($size)
 Function: Generates a random number between 0..$size
 Returns : Integer
 Args    : $maximum size for random number (defaults to 1)


=cut

sub random{
   my ($self,$max) = @_;
   $max = 2 unless defined $max || $max < 0;
   return int ( $max * rand());
}

1;
