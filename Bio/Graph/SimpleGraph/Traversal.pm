# Traversal.pm,v 1.10.2.1 2005/10/09 15:16:25 jason Exp
#
# BioPerl module for Bio::Graph::SimpleGraph::Traversal;
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Graph::SimpleGraph::Traversal - graph traversal operations for Bio::Graph::SimpleGraph and Bio::Graph::Protein::Graph objects 

=head1 SYNOPSIS

  use Bio::Graph::SimpleGraph::Traversal;
  use Bio::Graph::SimpleGraph;

  ## get a graph , $g.

  my $traversal = Bio::Graph::SimpleGraph::Traversal->new(-graph=>$g,
                                                          -start=>$start,
                                                          -order=>$order,
                                                          -what =>$what);
 ## cycle through nodes one at a time
 while ($traversal->has_next() ) {
        my $node = $traversal->get_next();
      }
 ## reset traversal to start
  $traversal->reset;

 ## get all nodes
  my @all_nodes = $traversal->get_all();



=head1 DESCRIPTION

This is a helper class for performing graph traversal operations for
Bio::Graph::SimpleGraph objects and Bio::Graph::Protein::Graph
objects. The documentation concerning the use of this class is
described in the "Graph algorithms" section of the
Bio::Graph::SimpleGraph modules. Only the methods are documented here.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Nat Goodman, Richard Adams

Email natg@shore.net, richard.adams@ed.ac.uk

=cut

package Bio::Graph::SimpleGraph::Traversal;
use vars qw(@AUTO_ATTRIBUTES @OTHER_ATTRIBUTES %SYNONYMS %DEFAULTS);
use Bio::Graph::SimpleGraph;
use strict;
use base qw(Class::AutoClass);

@AUTO_ATTRIBUTES=qw(order what graph start is_initialized
		    _past _present _future);
@OTHER_ATTRIBUTES=qw();
%SYNONYMS=();
%DEFAULTS=(order   => 'dfs',
	        what   => 'node',
	       _past   => {},
	      _future  => []);
Class::AutoClass::declare(__PACKAGE__);

sub _init_self {
	my($self,$class,$args)=@_;
	return unless $class eq __PACKAGE__; 
	# to prevent subclasses from re-running this
	$self->graph or $self->graph(new Bio::Graph::SimpleGraph);
	# can't be in DEFAULTS - circular includes!
}

=head2      has_next

 name      : has_next
 usage     : while (my $traversal->has_next() ) {..
 purpose   : returns true if there are more items in traversal, else undef
 arguments : none
 returns   : true or unde;

=cut 

sub has_next {
  my($self)=@_;
  $self->reset unless $self->is_initialized;
  @{$self->_future}>0;
}

=head2      get_next

 name      : get_next
 usage     : my $node =  $traversal->get_next() ;
 purpose   : returns  next item in traversal or undef if traversal is exhausted. 
 arguments : none
 returns   : a node  or undef;

=cut 

sub get_next {
  my($self)= @_;
  $self->reset unless $self->is_initialized;
  my $past   = $self->_past;
  my $future = $self->_future;
  my $present;
  my $graph  = $self->graph;
  while (@$future) {
    $present = shift @$future;
    unless($past->{$present}) {	# this is a new node
      $self->_present($present);
      $past->{$present}=1;
      if ($self->order =~ /^d/i) {
		unshift(@$future,$graph->neighbors($present,$self->what));
      } else {
		push(@$future,$graph->neighbors($present,$self->what));
      }
       return $present;
    }
  }
  $self->_present(undef);
}

=head2      get_all

 name      : get_all
 usage     : my @nodes =  $traversal->get_all() ;
 purpose   : get all remaining items in traversal as ARRAY (in array context)
              or ARRAY ref.
 arguments : none
 returns   : an array, an array reference or undef.

=cut 

sub get_all {
  my($self, $val)   = @_;
  $self->reset unless $self->is_initialized;
  my $past    = $self->_past;
  my $future  = $self->_future;
  my $i = 0; 
  my $present;
  my $graph   = $self->graph;
  my $nodes   = $graph->_nodes;

  my $results =[];
  while (@$future) {
    $present = shift @$future;
     if(!$past->{$present}) {	# this is a new node
         $past->{$present} = 1;

         push(@$results,$present);
		 $i++;
         if ($self->order =~ /^d/i) {
		    unshift(@$future,$graph->neighbors($present,$self->what));
             } else {
			push(@$future,$graph->neighbors($present,$self->what));
           }
        }
  }
  $self->_present(undef);
  wantarray? @$results: $results;
}

=head2      get_this

 name      : get_all
 usage     : my @nodes =  $traversal->get_all() ;
 purpose   : gets current node in traversal 
 arguments : none
 returns   : the current node or undef.

=cut 

sub get_this {
  my($self)=@_;
  $self->reset unless $self->is_initialized;
  $self->_present;
}

=head2      reset

 name      : reset
 usage     : $traversal->reset() ;
 purpose   : restarts traversal from first node
 arguments : none
 returns   : void.

=cut 

sub reset {
  my($self)= @_;
  $self->_past({});
  $self->order('d');
  $self->_present(undef);
  $self->_future([]);
  $self->is_initialized(1);
  my $graph = $self->graph;
  my $start = $self->start;
  my $what  = $self->what || 'node';
  if ($what=~/^n/i) {
    defined $start or $start=$graph->nodes->[0];
  } elsif ($what=~/^e/i) {
    $start=defined $start? $graph->edge($start): $graph->edges->[0];
  } else {
    $self->throw("Unrecognized \$what parameter $what: should be 'node' or 'edge'");
  }
  return unless defined $start;
  $self->_future([$start]);
}

1;
