#!/bin/perl -w
use strict;
package Bio::Graph::ProteinGraph;
use Bio::Graph::SimpleGraph;
use Clone qw(clone);
use vars  qw(@ISA);
our @ISA = qw(SimpleGraph);

=head1     NAME

Bio::Graph::Protein - a representation of a protein interaction graph.

=head1     SYNOPSIS

     ## read in from file

    my $graphio = Bio::Graph::IO->new(-file=>'myfile.dat', -format=>'dip');
    my $graph   = $graphio->next_network();

    ## get interactors of your favourite protein

    my $node      = $graph->nodes_by_id('NP_023232');
    my @neighbors = $graph->neighbors($node); 
    print "      NP_023232 interacts with ";
    print join " ," map{$_->object_id()} @neighbors;
    print "\n";

    ## annotate your sequences with interaction info

    my @my_seqs; ##array of sequence objects
    for my $seq(@seqs) {
		if ($graph->has_node($seq->accession_number)) {
          my $node      = $graph->nodes_by_id(
                                    $seq->accession_number);
          my @neighbors = $graph->neighbors($node);
          for my $n (@neighbors) {
            my $ft = Bio::SeqFeature::Generic->new(
                      -primary_tag => 'Interactor',
                      -tags        => { id => $n->accession_number
                                     }
                      );
            $seq->add_SeqFeature($ft);
         }
        }
    }

    ## get proteins with > 10 interactors

    my @nodes = $graph->nodes();
    my @hubs;
    for my $node (@nodes) {
		my @interactors = $graph->neighbors($node);
        if ($#interactors > 10) {	
				push @hubs, $node;
		}
	}
	print "the following proteins have > 10 interactors:\n"
    print join "\n", map{$_->object_id()} @hubs;

    ## merge 2 graphs, flag duplicate edges ##

   # get second graph $g2
   $g1->union($g2);
   my @duplicates = $g1->dup_edges();

   print "these interactions exist in $g1 and $g2:\n";
   print join "\n", map{$_->object_id} @duplicates;


=head1          DESCRIPTION

 A Protein graph is a representation of a protein interaction network.
 It derives most of its functionality from Nat Goodman's SimpleGraph module,
 but is adapted to be able to use protein identifiers to identify the nodes.

 At present it is fairly 'lightweight' in that it represents nodes and edges
 but does not contain all the data about experiment ids etc found in the 
 Protein Standards Initiative schema. Hopefully that will be available 
 soon.  

 For developers:

 In this module, nodes are represented by Bio::Seq::RichSeq objects containing 
 all possible database identifiers but no sequence, as parsed from the
 interaction files. However, a node represented by a Bio::PrimarySeq object 
 should work fine too. 

 Edges are represented by Bio::Graph::ProteinEdge objects. IN order to 
 work with SimpleGraph these objects must be array references, with the 
 first 2 elements being references to the 2 nodes. More data can be added
 in $e[2]. etc. Edges should implement the Bio::Graph::ProteinEdgeI interface
 which basically just demands an object_id() method. At present edges only
 have an identifier and a weight() method, to hold confidence data, but
 subclasses of this could hold all the interaction data  held in 
 an XML document.
 
=head1  REQUIREMENTS

To use this code you will need the Clone.pm module availabe from CPAN.
You also need Class::AutoClass availabe from CPAN as well. 

=head1 SEE ALSO

L<Bio::Graph::SimpleGraph>, 
L<Bio::Graph::IO>,
L<Bio::Graph::ProteinEdgeI>
L<Bio::DB::CUTG>


=head1 FEEDBACK


=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Richard Adams 

Email richard.adams@ed.ac.uk

=cut

=head2  has_node

 name      : has_node
 purpose   : Is a protein in the graph?
 usage     : if ($g->has_node('NP_23456')) {....
 returns   : 1 if true, 0 if false
 arguments : A sequence identifier.

=cut

sub has_node {

 my ($self, $arg) = @_;
 if (!$arg) {
		$self->throw ("I need a sequence identifier!");
   }
 my @nodes = $self->nodes_by_id($arg);
 if (@nodes){return 1;}else{return 0};



}


=head2    nodes_by_id

 name      : nodes_by_id
 purpose   : get node memory address from an id
 usage     : my @neighbors= $self->neighbors($self->nodes_by_id('O232322'))
 returns   : a SimpleGraph node representation ( a text representation
               of a node needed fro other graph methods e.g., neighbors()
              , edges()
 arguments : a protein identifier., e.g., its accession number. 

=cut

sub nodes_by_id {

	my $self  = shift;
	my @refs  = $self->_ids(@_);
	my @nodes = $self->nodes(@refs);
	wantarray? @nodes: $nodes[0];

}


=head2    union

 name        :  union
 purpose     : To merge two graphs together, flagging interactions as duplicate.
 usage       : $g1->union($g2), where g1 and g2 are 2 graph objects. 
 returns     : void, $g1 is modified
 arguments   : A Graph object of the same class as the calling object. 
 description : This method merges 2 graphs. The calling graph is modified, the
                 parameter graph ($g2) in usage) is unchanged. To take account
                 of differing IDs identifying the same protein, all ids are
                 compared. The following rules are used to modify $g1.

               First of all both graphs are scanned for nodes that share an id
                in common. 

               1. If 2 nodes(proteins) share an interaction in both graphs,
                     the edge in graph 2 is copied to graph 1 and added as a
                     duplicate edge to graph 1,

               2. If 2 nodes interact in $g2 but not $g1, but both nodes exist
                     in $g1, the attributes of the interaction in $g2 are 
                     used to make a new edge in $g1.

               3. If 2 nodes interact in g2 but not g1, and 1 of them is a new
                     protein, that protein is put in $g1 and a new edge made to
                     it. 

               4. At present, if there is an interaction in $g2 composed of a
                  pair of interactors that are not present in $g1, they are 
                   not copied to $g1. This is rather conservative but prevents
                  the problem of having redundant nodes in $g1 due to the same
                  protein being identified by different ids in the same graph.

               So, for example 

                        Edge   N1  N2 Comment

              Graph 1:  E1     P1  P2
                        E2     P3  P4
                        E3     P1  P4

              Graph 2:  X1     P1  P2 - will be added as duplicate to Graph1
                        X2     P1  X4 - X4 added to Graph 1 and new edge made
                        X3     P2  P3 - new edge links existing proteins in G1
                        X4     Z4  Z5 - not added to Graph1. Are these different
                                       proteins or synonyms for proteins in G1?
                         



=cut

sub union {
	
	my ($self, $other) = @_;
	my $class      = ref($self);
	if (!$other->isa($class)) {
		$self->throw("I need a ". $class . " object, not a [".
						ref($other). "] object");
		}
	my @common_nodes;
	my %detected_common_nodes;

	## for each node see if Ids are in common between the 2 graphs
	## just get1 common id per sequence
	for my $id (keys %{$self->{'_id_map'}}) {
		if (exists($other->{'_id_map'}{$id}) &&
			!exists($detected_common_nodes{$self->{'_id_map'}{$id}})) {
			push @common_nodes, $id;
			$detected_common_nodes{$self->{'_id_map'}{$id}} = undef;
			}
	}			
	
	## now cyle through common nodes..
	for my $common (@common_nodes) {
		
		## get neighbours of common node for self and other
		my @self_ns   = $self->neighbors($self->nodes_by_id($common));
		my @other_ns  = $other->neighbors($other->nodes_by_id($common));
	
		## now get all ids of all neighbours
		my %self_n_ids = $self->_get_ids(@self_ns); # get all ids of neighbors

		##cycle through other neighbors
		for my $other_n(@other_ns){ 
			my %other_n_ids = $self->_get_ids($other_n); # get ids of single other neighbor

            ## case (1) in description
			## do any ids exist in any self neighbors?
			#if yes,  @int_match is defined, is not a new interaction##
			my @int_match = grep{exists($self_n_ids{$_}) } keys %other_n_ids;
			my $is_in_self = 0;
			if (@int_match){
				my $i = 0;
				my $dup_edge;
				while (!$dup_edge && $i <= $#int_match){
					my $other_edge = $other->edge(
									[$other->nodes_by_id($common),
										$other->nodes_by_id($other_n->object_id)]
										);
					my $dup_edge = edge->new(
										-weight=> $other_edge->weight(),
										-id    => $other_edge->object_id(),
										-nodes=>[$self->nodes_by_id($common),
								   				 $self->nodes_by_id($int_match[0])
												]);
					$self->add_dup_edge($dup_edge);

					$i++;
				}
				
				eval{
					$self->add_dup_edge($dup_edge);
					};
				next if $@;
				$is_in_self = 1;
				}
			## this is a new interaction. Does interactor alresay 
			## now, is this new neighbour already in self, but not connected?
			## if it is, clonig not needed, just connect nodes in self using edge
			## attributes of other graph.
			else {
				my @existing_node = grep{exists($self->{'_id_map'}{$_} )  }keys %other_n_ids;
				
				if (@existing_node) {
					my $other_edge = $other->edge($other->nodes_by_id($other_n->object_id()),
												  $other->nodes_by_id($common));
					$self->add_edge(edge->new(
										-weight=> $other_edge->weight(),
										-id    => $other_edge->object_id(),
										-nodes=>[$self->nodes_by_id($common),
								   				 $self->nodes_by_id($existing_node[0])
												]),
									);
					$is_in_self = 1;
					}
				}
			## but if other neighbour is entirely new, clone it and make connection.
			if (!$is_in_self) {
				my $other_edge = $other->edge($other->nodes_by_id($other_n->object_id()),
											  $other->nodes_by_id($common));

				my $new = clone($other_n);
				$self->add_edge(edge->new(
									-weight => $other_edge->weight(),
									-id     => $other_edge->object_id(),
								    -nodes  =>[$new, $self->nodes_by_id($common)],
							    			)
								);

				## add new ids to self graph look up table
				map {$self->{'_id_map'}{$_} = $new} keys %other_n_ids;
	
				}
			
			}
	}
}
	
=head2      _get_ids
 
 name     : _get_ids
 purpose  : gets all ids for a node, assuming its Bio::Seq object
 arguments : A Bio::PrimarySeqI object
 returns  : A hash: Keys are sequence ids, values are undef
 usage    : my %ids = _get_ids($seqobj);

=cut

sub _get_ids {
	my %ids;
	my $dummy_self = shift;
	while (my $n = shift @_ ){  #ref to node, assume is a Bio::Seq
		if (!$n->isa('Bio::PrimarySeqI')) {
			$n->throw("I need a Bio::Seq object, not a [" .ref($n) ."]");
		}
		## get ids
		#map{$ids{$_} = undef}($n->accession_number, $n->primary_id);

		##if BioSeq getdbxref ids as well.
		if ($n->can('annotation')) {
			my $ac = $n->annotation();	
			for my $an($ac->get_Annotations('dblink')) {
				$ids{$an->database()} = $an->primary_id();
			}
		}
	}
	return %ids;
}

sub add_edge {
  my $self = shift;
  my $edges = $self->_edges;
  my $neighbors = $self->_neighbors;
  my $dup_edges = $self->_dup_edges;
	my $edge;
  while (@_) {
    if ( $_[0]->isa('Bio::Graph::Edge') ) {	# it's already an edge
       $edge = shift;
    } 
	elsif(ref($_[0]) eq 'ARRAY' || !ref($_[0])) {
      $self->SUPER::add_edges(@_);
		return;
    }
	else {
		$self->throw(" Invalid edge! - must be an array of nodes, or an edge object");
	}
	my ($m, $n) = $edge->nodes();
    next if $m eq $n;		# no self edges
    last unless defined $m && defined $n;
    ($m,$n)=($n,$m) if "$n" lt "$m";
    unless ($edges->{$m,$n}) {
      $self->add_node($m,$n);
      ($m,$n) = $self->nodes($m,$n);
      $edges->{$m,$n} = $edge;
      push(@{$neighbors->{$m}},$n);
      push(@{$neighbors->{$n}},$m);
    } else {
     $self->add_dup_edge($edge); 
    }
  }
  $self->_is_connected(undef);	# clear cached value

}

=head2      add_dup_edge

 name       : add_dup_edge
 purpose    : to flag an interaction as a duplicate, take advantage of edge ids
 returns    : 1 on successful addition, 0 on there being an existing duplicate. 
 usage      : $gr->add_dup_edge(edge->new (-nodes => [$n1, $n2],
                                           -score => $score
									       -id    => $id);
 arguments  : an edgeI implementing object.
 descripton : 
 );

=cut 

sub add_dup_edge {

	## get the 2 nodes
	my ($self, $edge) = @_;
	## prelimaries
	my $newedge_id   = $edge->object_id();

	## now we have node objects, an edge id.
	## is edge id new?
	my $dup_edges = $self->_dup_edges();
	if(!grep{$_->object_id eq $newedge_id } @$dup_edges) {
		push @$dup_edges, $edge;
		return 1;
		}
	else {
		$self->warn("2nd duplicate edge - $newedge_id");
		return 0;
	}
}


sub _ids {
	my $self = shift;
	my @refs;
	while (@_) {
		push @refs, $self->{'_id_map'}{shift @_ };
	}
	return @refs;
}
		
	

		
		
	

