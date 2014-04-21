#
# BioPerl module for Bio::TreeIO::TreeEventBuilder
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

Bio::TreeIO::TreeEventBuilder - Build Bio::Tree::Tree's and 
  Bio::Tree::Node's from Events 

=head1 SYNOPSIS

# internal use only

=head1 DESCRIPTION

This object will take events and build a Bio::Tree::TreeI compliant
object makde up of Bio::Tree::NodeI objects.

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


package Bio::TreeIO::TreeEventBuilder;
use strict;

use Bio::Tree::Tree;
use Bio::Tree::Node;

use base qw(Bio::Root::Root Bio::Event::EventHandlerI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::TreeIO::TreeEventBuilder->new();
 Function: Builds a new Bio::TreeIO::TreeEventBuilder object 
 Returns : Bio::TreeIO::TreeEventBuilder
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($treetype, $nodetype) = $self->_rearrange([qw(TREETYPE 
						    NODETYPE)], @args);
  $treetype ||= 'Bio::Tree::Tree';
  $nodetype ||= 'Bio::Tree::Node';

  eval { 
      $self->_load_module($treetype);
      $self->_load_module($nodetype);
  };

  if( $@ ) {
      $self->throw("Could not load module $treetype or $nodetype. \n$@\n")
  }
  $self->treetype($treetype);
  $self->nodetype($nodetype);
  $self->{'_treelevel'} = 0;
  return $self;
}

=head2 treetype

 Title   : treetype
 Usage   : $obj->treetype($newval)
 Function: 
 Returns : value of treetype
 Args    : newvalue (optional)


=cut

sub treetype{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'treetype'} = $value;
    }
    return $self->{'treetype'};
}

=head2 nodetype

 Title   : nodetype
 Usage   : $obj->nodetype($newval)
 Function: 
 Returns : value of nodetype
 Args    : newvalue (optional)


=cut

sub nodetype{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'nodetype'} = $value;
    }
    return $self->{'nodetype'};
}


=head2 SAX methods

=cut

=head2 start_document

 Title   : start_document
 Usage   : $handler->start_document
 Function: Begins a Tree event cycle
 Returns : none 
 Args    : none

=cut

sub start_document {
   my ($self) = @_;   
   $self->{'_lastitem'} = {};
   $self->{'_currentitems'} = [];
   $self->{'_currentnodes'} = [];
   return;
}

=head2 end_document

 Title   : end_document
 Usage   : my @trees = $parser->end_document
 Function: Finishes a Phylogeny cycle
 Returns : An array  Bio::Tree::TreeI
 Args    : none

=cut

sub end_document {
    my ($self,$label) = @_; 

    my ($root) = @{$self->{'_currentnodes'}};

    $self->debug("Root node is " . $root->to_string()."\n");
    if( $self->verbose > 0 ) { 
	foreach my $node ( $root->get_Descendents  ) {
	    $self->debug("node is ". $node->to_string(). "\n");
	}
    }
    my $tree = $self->treetype->new(-verbose => $self->verbose,
				    -root => $root);
    return $tree;       
}

=head2 start_element

 Title   : start_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    : $data => hashref with key 'Name'

=cut

sub start_element{
   my ($self,$data) =@_;
   $self->{'_lastitem'}->{$data->{'Name'}}++;   

   $self->debug("starting element: $data->{Name}\n");   
   push @{$self->{'_lastitem'}->{'current'}},$data->{'Name'};
   
   my %data;
   
   if( $data->{'Name'} eq 'node' ) {
       push @{$self->{'_currentitems'}}, \%data; 
       $self->{'_treelevel'}++;
   } elsif ( $data->{Name} eq 'tree' ) {
   }
}

=head2 end_element

 Title   : end_element
 Usage   : 
 Function:
 Returns : none
 Args    : $data => hashref with key 'Name'

=cut

sub end_element{
   my ($self,$data) = @_;   

   $self->debug("end of element: $data->{Name}\n");
   # this is the stack where we push/pop items from it
   my $curcount = scalar @{$self->{'_currentnodes'}};
   my $level   = $self->{'_treelevel'};
   my $levelct = $self->{'_nodect'}->[$self->{'_treelevel'}+1] || 0;

   if( $data->{'Name'} eq 'node' ) {
       my $tnode;
       my $node = pop @{$self->{'_currentitems'}};	   

       $tnode = $self->nodetype->new( -verbose => $self->verbose,
				      %{$node});       
       $self->debug( "new node will be ".$tnode->to_string."\n");
       if ( !$node->{'-leaf'} && $levelct > 0) {
	   $self->debug(join(',', map { $_->to_string } 
			     @{$self->{'_currentnodes'}}). "\n");
	   $self->throw("something wrong with event construction treelevel ".
			"$level is recorded as having $levelct nodes  ".
			"but current nodes at this level is $curcount\n")
	       if( $levelct > $curcount);	
	   for ( splice( @{$self->{'_currentnodes'}}, - $levelct)) {
	       $self->debug("adding desc: " . $_->to_string . "\n");
	       $tnode->add_Descendent($_);
	   }
	   $self->{'_nodect'}->[$self->{'_treelevel'}+1] = 0;
       }
       push @{$self->{'_currentnodes'}}, $tnode;
       $self->{'_nodect'}->[$self->{'_treelevel'}]++;
       
       $curcount = scalar @{$self->{'_currentnodes'}};
       $self->debug ("added node: count is now $curcount, treelevel: $level, nodect: $levelct\n");

       $self->{'_treelevel'}--;       
   } elsif(  $data->{'Name'} eq 'tree' ) { 
       $self->debug("end of tree: nodes in stack is $curcount\n");
   }

   $self->{'_lastitem'}->{ $data->{'Name'} }--; 
   pop @{$self->{'_lastitem'}->{'current'}};
}


=head2 in_element

 Title   : in_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub in_element{
   my ($self,$e) = @_;

   return 0 if ! defined $self->{'_lastitem'} || 
       ! defined $self->{'_lastitem'}->{'current'}->[-1];
   return ($e eq $self->{'_lastitem'}->{'current'}->[-1]);

}

=head2 within_element

 Title   : within_element
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub within_element{
   my ($self,$e) = @_;
   return $self->{'_lastitem'}->{$e};
}

=head2 characters

 Title   : characters
 Usage   : $handler->characters($text);
 Function: Processes characters 
 Returns : none
 Args    : text string


=cut

sub characters{
   my ($self,$ch) = @_;
   if( $self->within_element('node') ) {
       my $hash = pop @{$self->{'_currentitems'}};
       if( $self->in_element('bootstrap') ) {
	   # leading/trailing Whitespace-B-Gone
	   $ch =~ s/^\s+//; $ch =~ s/\s+$//;  
	   $hash->{'-bootstrap'} = $ch;
       } elsif( $self->in_element('branch_length') ) {
	   # leading/trailing Whitespace-B-Gone
	   $ch =~ s/^\s+//; $ch =~ s/\s+$//;
	   $hash->{'-branch_length'} = $ch;
       } elsif( $self->in_element('id')  ) {
	   $hash->{'-id'} = $ch;
       } elsif( $self->in_element('description') ) {
	   $hash->{'-desc'} = $ch;
       } elsif ( $self->in_element('tag_name') ) {
	   $hash->{'-NHXtagname'} = $ch;
       } elsif ( $self->in_element('tag_value') ) {
	   $hash->{'-nhx'}->{$hash->{'-NHXtagname'}} = $ch;
	   delete $hash->{'-NHXtagname'};
       } elsif( $self->in_element('leaf') ) {
	   $hash->{'-leaf'} = $ch;
       }
       push @{$self->{'_currentitems'}}, $hash;
   }
   $self->debug("chars: $ch\n");
}


1;
