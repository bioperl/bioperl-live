package Bio::Graph::SimpleGraph::Traversal;
use vars qw(@ISA @AUTO_ATTRIBUTES @OTHER_ATTRIBUTES %SYNONYMS %DEFAULTS);
use Class::AutoClass;
use Bio::Graph::SimpleGraph;
use strict;
@ISA = qw(Class::AutoClass); # AutoClass must be first!!

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
  return unless $class eq __PACKAGE__; # to prevent subclasses from re-running this
  $self->graph or $self->graph(new Bio::Graph::SimpleGraph); # can't be in DEFAULTS - circular includes!
}

sub has_next {
  my($self)=@_;
  $self->reset unless $self->is_initialized;
  @{$self->_future}>0;
}

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
sub get_all {
  my($self, $val)   = @_;
  $self->reset unless $self->is_initialized;
  my $past    = $self->_past;
  my $future  = $self->_future;
  my $i = 0; 
  my $present;
  my $graph   = $self->graph;
  my $nodes   = $graph->_nodes;

 ##temporary fix for dealing with nodes that are references
 my $is_href = 1;
 for my $n (keys %$nodes) {
	if( !ref ($nodes->{$n}) || ref($nodes->{$n}) eq  'ARRAY') {
		$is_href = 0;
		last;
		}
 	 $nodes->{$n}{'_node_id'} = undef;
	}
	
  my $results =[];
  while (@$future) {
    $present = shift @$future;
     if(!$past->{$present}) {	# this is a new node
         $past->{$present} = 1;
		
         push(@$results,$present);
		 $nodes->{$present}{'_node_id'} = $i if $is_href;
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

sub get_this {
  my($self)=@_;
  $self->reset unless $self->is_initialized;
  $self->_present;
}
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
