#!/bin/perl -w
package NetworkIO::psi_xml;
use strict;
use XML::Twig;
use Bio::Seq::SeqFactory;
use prot_graph;
use edge;
use NetworkIO;
use Bio::Annotation::DBLink;
use Bio::Annotation::Collection;
use Bio::Species;
use vars qw(@ISA  %species @g $c $fac);
@ISA = qw(NetworkIO);

BEGIN{
		 $fac  = Bio::Seq::SeqFactory->new(
							-type => 'Bio::Seq::RichSeq'
						);
$c = 0;
}
#parsing done by XML::Twig, not by RootIO, therefore override usual new
sub new {

my ($class,@args) = @_;
my $self = bless {}, $class;
$self->_initialize(@args);
return $self;

}

sub _initialize  {

  my($self,@args) = @_;
  return unless $self->SUPER::_initialize_io(@args);

}

=head2     next_network

 name       : next_network
 purpose    : to construct a protein interaction graph from xml data
 usage      : my $gr = $netowkio->next_network();
 arguments  : void
 returns    : A Bio::Network::Proteingraph object

=cut

sub next_network {

 my $self = shift;
 push @g, prot_graph->new();
 my $t    = XML::Twig->new
				(  TwigHandlers => {
									 proteinInteractor   => \&_proteinInteractor,
								 	 interaction         => \&_addEdge
									});
 $t->parsefile($self->file);
 return $g[$#g];	 


}

=head2   _proteinInteractor

 name      : _proteinInteractor
 purpose   : parses protein information into Bio::Seq::RichSeq objects
 returns   : void
 usage     : internally called by next_network(), 
 arguments : none.

=cut

sub _proteinInteractor {
	my ($twig, $pi) = @_;
	$c++;
	if ($c % 50 == 0) {
		print  STDERR ".";
	}

	my ($acc, $sp, $desc, $taxid,  $prim_id);

	my $org =  $pi->first_child('organism');
	$taxid  = $org->att('ncbiTaxId');
	
	## just make new species object if doesn't already exist ##
	if (!exists($species{$taxid})) {
		my $full       =  $org->first_child('names')->first_child('fullName')->text;
		my ($gen, $sp) = $full =~ /(\S+)\s+(.+)/;
		my $sp_obj     = Bio::Species->new(-ncbi_taxid    => $taxid,
									    -classification=> [$sp, $gen],
									);
		$species{$taxid} = $sp_obj;
		} 
	
	## next extract sequence id info ##
	my @ids          = $pi->first_child('xref')->children();
	my %ids          = map{$_->att('db'), $_->att('id')} @ids;
	 $ids{'psixml'}  = $pi->att('id');
	
	
	$prim_id = defined ($ids{'GI'})?  $ids{'GI'}:'';
	$acc        = $ids{'RefSeq'} || $ids{'SWP'} || $ids{'PIR'} || $ids{'GI'};
	
	## get description line
	$desc    = $pi->first_child('names')->first_child('fullName')->text;

	## use ids that aren't accession_no or primary_tag to build dbxref Annotations
	my $ac = Bio::Annotation::Collection->new();	
	for my $db (keys %ids) {
		next if $ids{$db} eq $acc or $ids{$db} eq $prim_id;
		my $an = Bio::Annotation::DBLink->new( -database   => $db,
											   -primary_id => $ids{$db},
											);
		$ac->add_Annotation('dblink',$an);
	}

	## now we can make sequence object ##
	my $node = $fac->create(
						-accession_number => $acc,
						-desc             => $desc,
						-display_id       => $acc,
						-primary_id       => $prim_id,
						-species          => $species{$taxid},
						-annotation       => $ac);
	
	## now fill hash with keys = ids and vals = node refs to have lookip
	## hash for nodes by any id.	
	$g[$#g]{'_id_map'}{$ids{'psixml'}}          = $node;
	if (defined($node->primary_id)) {
		$g[$#g]{'_id_map'}{$node->primary_id} = $node;
		}
	if (defined($node->accession_number)) {
		$g[$#g]{'_id_map'}{$node->accession_number} = $node;
		}
	## cycle thru annotations
	 $ac = $node->annotation();
	for my $an ($ac->get_Annotations('dblink')) {
		$g[$#g]{'_id_map'}{$an->primary_id} = $node;
		}
	$twig->purge();
}

=head2      add_edge

 name     : add_edge
 purpose  : adds a new edge to a graph
 usage    : do not call, called by next_network
 returns  : void
 
=cut

sub _addEdge {
	$c++; 
	if ($c % 50 ==0 ) {
		print STDERR ",";
		}
	my ($twig, $i) = @_;
	my @ints = $i->first_child('participantList')->children;
	my @node = map{$_->first_child('proteinInteractorRef')->att('ref')} @ints;
	$g[$#g]->add_edge(edge->new(
					-nodes =>[($g[$#g]{'_id_map'}{$node[0]}, $g[$#g]{'_id_map'}{$node[1]})],
					-id    => undef));
	$twig->purge();
}
1;

