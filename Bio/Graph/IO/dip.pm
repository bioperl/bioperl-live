#!/bin/perl -w
package Bio::Graph::IO::dip;
use vars qw(@ISA $FAC);
use lib '../';
use Bio::Graph::IO;
use Bio::Graph::ProteinGraph;
use Bio::Seq::SeqFactory;
use Bio::Annotation::DBLink;
use Bio::Annotation::Collection;
use Bio::Graph::Edge;
use strict;
@ISA = qw(Bio::Graph::IO);
BEGIN{
	$FAC = Bio::Seq::SeqFactory->new(-type=>'Bio::Seq::RichSeq');
}
#sub new {
#	my $class = shift;
#	my $self = $class->SUPER::new(@_);
#	return $self;
#}

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

	   ## build node object if is new node
	   my ($node1, $node2);
	   if(!exists($seen_nodes{$n1}) ) {

			if($g1){ $g1 =~ s/GI://;     }
			if($p1){ $p1 =~ s/PIR://;  }
			if($s1){ $s1 =~ s/SWP://;  }
			my $acc = $s1 || $p1 || $g1;
			my $ac  = $self->_add_db_links($acc, $s1, $p1,  $n1);
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
		if(!exists($seen_nodes{$n2}) ) {
			if($g2){$g2 =~ s/GI://; }
			if($p2){$p2 =~ s/PIR://; }
			if($s2){$s2 =~ s/SWP://; }
			my $acc = $s2 || $p2 || $g2;
			my $ac  = $self->_add_db_links($acc, $s2, $p2,  $n2);
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

sub write_network{

}

sub _add_db_links {
	my ($self, $acc, $s1, $p1,  $n1) = @_;
	my %ids;
	$ids{'pir'} = $p1 if $p1;
	$ids{'swp'} = $s1 if $s1;
	$ids{'dip'} = $n1 if $n1;
	my $ac = Bio::Annotation::Collection->new();
	for my $db (keys %ids) {
		next if  $ids{$db}  eq $acc;
		my $an = Bio::Annotation::DBLink->new( -database   => $db,
											   -primary_id => $ids{$db},
											);
		$ac->add_Annotation('dblink', $an);
	}
	return $ac;
}
		
		
	

1;


