# $Id$
#
# BioPerl module for Bio::Graph::IO::dip
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Graph::IO::dip - class for parsing interaction data in dip format

=head1 SYNOPSIS

Do not use this module directly, use Bio::Graph::IO, for example:

  my $graph_io = Bio::Graph::IO->new(-format => 'dip',
                                     -file   => 'data.dip');

=head1 METHODS

The naming system is analagous to the SeqIO system, although usually
next_network() will be called only once per file.

=cut

package Bio::Graph::IO::dip;
use vars qw($FAC);

use Bio::Graph::ProteinGraph;
use Bio::Seq::SeqFactory;
use Bio::Annotation::DBLink;
use Bio::Annotation::Collection;
use Bio::Graph::Edge;
use strict;
use base qw(Bio::Graph::IO);

BEGIN{
	$FAC = Bio::Seq::SeqFactory->new(-type=>'Bio::Seq::RichSeq');
}

=head2        next_network

  name        : next_network
  purpose     : parses a graph file and returns a Bio::Graph::ProteinGraph 
                object
  usage       : my $g = $graph_io->next_network();
  arguments   : none
  returns     : a Bio::Graph::ProteinGraph object

=cut

sub next_network {

	my $self = shift;
	my %seen_nodes = ();
	my $graph = new Bio::Graph::ProteinGraph();

	while (my $l = $self->_readline() ) {

		##get line, only gi and node_id always defined
		my ($edge_id, $n1, $s1, $p1, $g1, $n2, $s2, $p2, $g2, $score) =
	 	$l =~/^DIP:(\d+E)\t+
			(DIP\S+)\t+
			(SWP\S+)?\t*
			(PIR\S+)?\t*
			(GI\S+)\t+
			(DIP\S+)\t+
			(SWP\S+)?\t*
			(PIR\S+)?\t*
			(GI\S+)\t*
			(\d\.\d+)? #optional confidence  or weight score 
			/x;

		## skip if score is below threshold
		if (defined($self->{'_th'}) && defined($score)) {
			next unless $score >= $self->{'_th'};
		}

	   ## build node object if is new node
	   my ($node1, $node2);
	 	$n1 =~ s/DIP://;  
	   if(!exists($seen_nodes{$n1}) ) {
			if($g1){ $g1 =~ s/GI://;     }
			if($p1){ $p1 =~ s/PIR://;  }
			if($s1){ $s1 =~ s/SWP://;  }
			my $acc = $s1 || $p1 || $g1;
			my $ac  = $self->_add_db_links($acc, $s1, $p1,  $n1, $g1);
			$node1 = $FAC->create(
								-accession_number => $acc,
								-primary_id       => $g1,
								-display_id		  => $acc,
								-annotation       => $ac,
								);		
			for my $n ($g1, $p1, $s1, $n1) {
				$seen_nodes{$n} = $node1 if $n;
				}
			} else {
			$node1 = $seen_nodes{$n1};
		}
	 	$n2 =~ s/DIP://;  
		if(!exists($seen_nodes{$n2}) ) {
			if($g2){$g2 =~ s/GI://; }
			if($p2){$p2 =~ s/PIR://; }
			if($s2){$s2 =~ s/SWP://; }
			my $acc = $s2 || $p2 || $g2;
			my $ac  = $self->_add_db_links($acc, $s2, $p2,  $n2, $g2);
			$node2  = $FAC->create(
								-accession_number => $acc,
								-primary_id       => $g2,
								-display_id		  => $acc,
								-annotation       => $ac,
								);		
			for my $n ($g2, $p2, $s2, $n2) {
				$seen_nodes{$n} = $node2 if $n;
				}
		  } else {
			$node2 = $seen_nodes{$n2};
		}

		## create new edge object based on node, weight. 
		$graph->add_edge(Bio::Graph::Edge->new( -nodes  => [$node1, $node2],
									-weight => $score,
									-id     => $edge_id),
									);
	}

	## now ensure nodes are accessible by either 1ary or 2ndary ids. 
	$graph->{'_id_map'} = \%seen_nodes;
	return $graph;
}

=head2      write_network

 name     : write_network
 purpose  : write graph out in dip format
 arguments: a Bio::Graph::ProteinGraph object
 returns  : void
 usage    : $out->write_network($gr);

=cut

sub write_network {

my ($self, $gr) = @_;
if (!$gr || !$gr->isa('Bio::Graph::ProteinGraph')) {
	$self->throw("I need a Bio::Graph::ProteinGraph, not a [".
	              ref($gr) . "]");
	  }
my @edges = $gr->edges();

# need to have all ids as annotations with database ids as well
# idea is to be able to round trip, to write it in same way as 

#for each edge	
for my $edge (@edges) {
	my $str  = "DIP:" .$edge->object_id(). "\t"; #output string
	my @nodes = $edge->nodes();

	# add node ids to string in correct order
	for my $n (@nodes){

	    # print out nodes in dip order
		my %ids = $gr->_get_ids_by_db($n); #need to modify this in graph()
		for my $db (qw(DIP SWP PIR GI Ref-Seq RefSeq psixml ens)){
			if (exists($ids{$db})){
				$str .= "$db:$ids{$db}\t";
			} else {
				$str .= "\t";
			}
		}
	}
	# add weights if defined
	$str =~ s/\t$//;
	if(defined($edge->weight)) {
		$str .= "\t" .$edge->weight. "\n";
		}else {
		$str .= "\n";
	}
	$self->_print($str);
 }# next edge
$self->flush();
}

sub _add_db_links {
	my ($self, $acc, $s1, $p1,  $n1, $g1) = @_;
	my %ids;
	$ids{'PIR'} = $p1 if $p1;
	$ids{'SWP'} = $s1 if $s1;
	$ids{'DIP'} = $n1 if $n1;
	$ids{'GI'}  = $g1 if $g1;
	my $ac = Bio::Annotation::Collection->new();
	for my $db (keys %ids) {
		#next if  $ids{$db}  eq $acc;
		my $an = Bio::Annotation::DBLink->new( -database   => $db,
											   -primary_id => $ids{$db},
											);
		$ac->add_Annotation('dblink', $an);
	}
	return $ac;
}

1;
