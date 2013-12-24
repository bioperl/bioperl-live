package Bio::Coordinate::Graph;
use utf8;
use strict;
use warnings;
use parent qw(Bio::Root::Root);

# ABSTRACT: Finds shortest path between nodes in a graph.
# AUTHOR:   Heikki Lehvaslaiho <heikki@bioperl.org>
# OWNER:    Heikki Lehvaslaiho
# LICENSE:  Perl_5

=head1 SYNOPSIS

  # get a hash of hashes representing the graph. E.g.:
  my $hash= {
             '1' => {
                     '2' => 1
                    },
             '2' => {
                     '4' => 1,
                     '3' => 1
                    },
             '3' => undef,
             '4' => {
                     '5' => 1
                    },
             '5' => undef
            };

  # create the object;
  my $graph = Bio::Coordinate::Graph->new(-graph => $hash);

  # find the shortest path between two nodes
  my $a = 1;
  my $b = 6;
  my @path = $graph->shortest_paths($a);
  print join (", ", @path), "\n";

=head1 DESCRIPTION

This class calculates the shortest path between input and output
coordinate systems in a graph that defines the relationships between
them. This class is primarely designed to analyze gene-related
coordinate systems. See L<Bio::Coordinate::GeneMapper>.

Note that this module can not be used to manage graphs.

Technically the graph implemented here is known as Directed Acyclic
Graph (DAG). DAG is composed of vertices (nodes) and edges (with
optional weights) linking them. Nodes of the graph are the coordinate
systems in gene mapper.

The shortest path is found using the Dijkstra's algorithm. This
algorithm is fast and greedy and requires all weights to be
positive. All weights in the gene coordinate system graph are
currently equal (1) making the graph unweighted. That makes the use of
Dijkstra's algorithm an overkill. A simpler and faster breadth-first
would be enough. Luckily the difference for small graphs is not
significant and the implementation is capable of taking weights into
account if needed at some later time.

=head2 Input format

The graph needs to be primed using a hash of hashes where there is a
key for each node. The second keys are the names of the downstream
neighboring nodes and values are the weights for reaching them. Here
is part of the gene coordiante system graph:

    $hash = {
             '6' => undef,
             '3' => {
                     '6' => 1
                    },
             '2' => {
                     '6' => 1,
                     '4' => 1,
                     '3' => 1
                    },
             '1' => {
                     '2' => 1
                    },
             '4' => {
                     '5' => 1
                    },
             '5' => undef
            };

Note that the names need to be positive integers. Root should be '1'
and directness of the graph is taken advantage of to speed
calculations by assuming that downsream nodes always have larger
number as name.

An alternative (shorter) way of describing input is to use hash of
arrays. See L<Bio::Coordinate::Graph::hash_of_arrays>.

=cut

=head2 new
=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my($graph, $hasharray) =
        $self->_rearrange([qw(
                              GRAPH
                              HASHARRAY
                             )],
                         @args);

    $graph  && $self->graph($graph);
    $hasharray  && $self->hasharray($hasharray);

    $self->{'_root'} = undef;

    return $self; # success - we hope!
}

=head2 graph

 Title   : graph
 Usage   : $obj->graph($my_graph)
 Function: Read/write method for the graph structure
 Example :
 Returns : hash of hashes grah structure
 Args    : reference to a hash of hashes

=cut

sub graph {

  my ($self,$value) = @_;

  if ($value) {
      $self->throw("Need a hash of hashes")
          unless  ref($value) eq 'HASH' ;
      $self->{'_dag'} = $value;

      # empty the cache
      $self->{'_root'} = undef;

  }

  return $self->{'_dag'};

}

=head2 hash_of_arrays

 Title   : hash_of_arrays
 Usage   : $obj->hash_of_array(%hasharray)
 Function: An alternative method to read in the graph structure.
           Hash arrays are easier to type. This method converts
           arrays into hashes and assigns equal values "1" to
           weights.

 Example : Here is an example of simple structure containing a graph.

           my $DAG = {
                      6  => [],
                      5  => [],
                      4  => [5],
                      3  => [6],
                      2  => [3, 4, 6],
                      1  => [2]
                     };

 Returns : hash of hashes graph structure
 Args    : reference to a hash of arrays

=cut

sub hash_of_arrays {

  my ($self,$value) = @_;

  # empty the cache
  $self->{'_root'} = undef;

  if ($value) {

      $self->throw("Need a hash of hashes")
          unless  ref($value) eq 'HASH' ;

      #copy the hash of arrays into a hash of hashes;
      my %hash;
      foreach my $start ( keys %{$value}){
          $hash{$start} = undef;
          map { $hash{$start}{$_} = 1 } @{$value->{$start}};
      }

      $self->{'_dag'} = \%hash;
  }

  return $self->{'_dag'};

}

=head2 shortest_path

 Title   : shortest_path
 Usage   : $obj->shortest_path($a, $b);
 Function: Method for retrieving the shortest path between nodes.
           If the start node remains the same, the method is sometimes
           able to use cached results, otherwise it will recalculate
           the paths.
 Example :
 Returns : array of node names, only the start node name if no path
 Args    : name of the start node
         : name of the end node

=cut

sub shortest_path {
    my ($self, $root, $end) = @_;

    $self->throw("Two arguments needed") unless @_ == 3;
    $self->throw("No node name [$root]")
        unless exists $self->{'_dag'}->{$root};
    $self->throw("No node name [$end]")
        unless exists $self->{'_dag'}->{$end};

    my @res;     # results
    my $reverse;

    if ($root > $end) {
        ($root, $end) = ($end, $root );
        $reverse++;
    }

    # try to use cached paths
    $self->dijkstra($root) unless
        defined $self->{'_root'} and $self->{'_root'} eq $root;

    return @res unless $self->{'_paths'} ;

    # create the list
    my $node = $end;
    my $prev = $self->{'_paths'}->{$end}{'prev'};
    while ($prev) {
        unshift @res, $node;
        $node = $self->{'_paths'}->{$node}{'prev'};
        $prev = $self->{'_paths'}->{$node}{'prev'};
    }
    unshift @res, $node;

    $reverse ? return reverse @res : return @res;
}

=head2 dijkstra

 Title   : dijkstra
 Usage   : $graph->dijkstra(1);
 Function: Implements Dijkstra's algorithm.
           Returns or sets a list of mappers. The returned path
           description is always directed down from the root.
           Called from shortest_path().
 Example :
 Returns : Reference to a hash of hashes representing a linked list
           which contains shortest path down to all nodes from the start
           node. E.g.:

            $res = {
                      '2' => {
                               'prev' => '1',
                               'dist' => 1
                             },
                      '1' => {
                               'prev' => undef,
                               'dist' => 0
                             },
                    };

 Args    : name of the start node

=cut

sub dijkstra {
    my ($self,$root) = @_;

    $self->throw("I need the name of the root node input") unless $root;
    $self->throw("No node name [$root]")
        unless exists $self->{'_dag'}->{$root};

    my %est = ();          # estimate hash
    my %res = ();          # result hash
    my $nodes = keys %{$self->{'_dag'}};
    my $maxdist = 1000000;

    # cache the root value
    $self->{'_root'} = $root;

    foreach my $node ( keys %{$self->{'_dag'}} ){
        if ($node eq $root) {
            $est{$node}{'prev'} = undef;
            $est{$node}{'dist'} = 0;
        } else {
            $est{$node}{'prev'} = undef;
            $est{$node}{'dist'} = $maxdist;
        }
    }

    # remove nodes from %est until it is empty
    while (keys %est) {

        #select the node closest to current one, or root node
        my $min_node;
        my $min = $maxdist;
        foreach my $node (reverse sort keys %est) {
            if ( $est{$node}{'dist'} < $min ) {
                $min = $est{$node}{'dist'};
                $min_node = $node;
            }
        }

        # no more links between nodes
        last unless ($min_node);

        # move the node from %est into %res;
        $res{$min_node} = delete $est{$min_node};

        # recompute distances to the neighbours
        my $dist = $res{$min_node}{'dist'};
        foreach my $neighbour ( keys %{$self->{'_dag'}->{$min_node}} ){
            next unless $est{$neighbour}; # might not be there any more
            $est{$neighbour}{'prev'} = $min_node;
            $est{$neighbour}{'dist'} =
                $dist + $self->{'_dag'}{$min_node}{$neighbour}
                if $est{$neighbour}{'dist'} > $dist + 1 ;
        }
    }
    return $self->{'_paths'} = \%res;
}

1;
