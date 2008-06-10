# $Id$
#
# BioPerl module for Bio::Graph::IO::psi_xml
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Graph::IO::psi_xml - access and manipulate PSI XML graphs

=head1 SYNOPSIS

Do not use this module directly, use Bio::Graph::IO, for example:

  my $graph_io = Bio::Graph::IO->new(-format => 'psi_xml',
                                     -file   => 'data.xml');

=head1 DESCRIPTION

PSI XML is a format to describe protein-protein interactions and 
interaction networks. The following databases support PSI XML:

=over 4

=item BIND  

L<http://www.bind.ca>

=item DIP 

L<http://dip.doe-mbi.ucla.edu/>

=item HPRD    

L<http://www.hprd.org>

=item IntAct  

L<http://www.ebi.ac.uk/intact>

=item MINT    

L<http://cbm.bio.uniroma2.it/mint/>

=back 

Notes on PSI XML from various databases can be found in the Bioperl Wiki
at L<http://bioperl.org/wiki/Module:Bio::Graph::IO::psi_xml>

Documentation for PSI XML can be found at L<http://psidev.sourceforge.net>.

=head1 METHODS

The naming system is analagous to the SeqIO system, although usually
L<next_network()> will be called only once per file.

=cut

package Bio::Graph::IO::psi_xml;
use strict;
use XML::Twig;
use Bio::Seq::SeqFactory;
use Bio::Graph::ProteinGraph;
use Bio::Graph::Edge;
use Bio::Annotation::DBLink;
use Bio::Annotation::Collection;
use Bio::Species;
use vars qw(%species $g $c $fac);
use base qw(Bio::Graph::IO);

BEGIN{
	$fac  = Bio::Seq::SeqFactory->new(-type => 'Bio::Seq::RichSeq');
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
 usage      : my $gr = $io->next_network();
 arguments  : void
 returns    : A Bio::Graph::ProteinGraph object

=cut

sub next_network {

 my $self = shift;
 $g = Bio::Graph::ProteinGraph->new(); ## bugfix, now is reset each time
 my $t    = XML::Twig->new
	(  TwigHandlers => {
							  proteinInteractor   => \&_proteinInteractor,
							  interaction         => \&_addEdge
							 });
  $t->parsefile($self->file);
 return $g;
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

	my ($acc, $sp, $desc, $taxid,  $prim_id);

	my $org =  $pi->first_child('organism');
	$taxid  = $org->att('ncbiTaxId');

	## just make new species object if doesn't already exist ##
	if (!exists($species{$taxid})) {
		my $common     =  $org->first_child('names')->first_child('shortLabel')->text;
		my $full       =  $org->first_child('names')->first_child('fullName')->text;
		my ($gen, $sp) = $full =~ /(\S+)\s+(.+)/;
		my $sp_obj     = Bio::Species->new(-ncbi_taxid     => $taxid,
											-classification => [$sp, $gen],
											-common_name    => $common
										   );
		$sp_obj->name('scientific', $full);
		$species{$taxid} = $sp_obj;
        print "species parse error $@" if $@;
      }
      

	## next extract sequence id info ##
	my @ids          = $pi->first_child('xref')->children();
	my %ids          = map {$_->att('db'), $_->att('id')} @ids;
	$ids{'psixml'}  = $pi->att('id');

	$prim_id = defined ($ids{'GI'}) ?  $ids{'GI'} : '';
	$acc = $ids{'RefSeq'} || 
	       $ids{'SWP'} || 
			 $ids{'Swiss-Prot'} || # db name from HPRD
			 $ids{'Ref-Seq'} ||    # db name from HPRD
			 $ids{'GI'} || 
			 $ids{'PIR'} ||
			 $ids{'intact'} ||     # db name from IntAct
			 $ids{'psi-mi'};       # db name from IntAct

	## get description line - certain files, like PSI XML from HPRD, have
	## "shortLabel" but no "fullName"
	eval {
		$desc = $pi->first_child('names')->first_child('fullName')->text; 
	};
	if ($@) {
		warn("No fullName, use shortLabel for description instead");
		$desc = $pi->first_child('names')->first_child('shortLabel')->text;
	}
	
	# use ids that aren't accession_no or primary_tag to build 
   # dbxref Annotations
	my $ac = Bio::Annotation::Collection->new();	
	for my $db (keys %ids) {
		next if $ids{$db} eq $acc;
		next if $ids{$db} eq $prim_id;
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

	## now fill hash with keys = ids and vals = node refs to have lookup
	## hash for nodes by any id.	
	$g->{'_id_map'}{$ids{'psixml'}} = $node;
	if (defined($node->primary_id)) {
		$g->{'_id_map'}{$node->primary_id} = $node;
	}
	if (defined($node->accession_number)) {
		$g->{'_id_map'}{$node->accession_number} = $node;
	}

	## cycle thru annotations
	$ac = $node->annotation();
	for my $an ($ac->get_Annotations('dblink')) {
		$g->{'_id_map'}{$an->primary_id} = $node;
	}
	$twig->purge();
}

=head2 add_edge

 name     : _addEdge
 purpose  : adds a new edge to a graph
 usage    : do not call, called by next_network
 returns  : void

=cut

sub _addEdge {

	my ($twig, $i) = @_;
	my @ints = $i->first_child('participantList')->children;
	my @node = map {$_->first_child('proteinInteractorRef')->att('ref')} @ints;
	my $edge_id = $i->first_child('xref')->first_child('primaryRef')->att('id');
	$g->add_edge(Bio::Graph::Edge->new(
					-nodes =>[($g->{'_id_map'}{$node[0]}, 
                               $g->{'_id_map'}{$node[1]})],
					-id    => $edge_id));
	$twig->purge();

}

1;
