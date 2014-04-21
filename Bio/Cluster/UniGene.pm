#
# BioPerl module for Bio::Cluster::UniGene.pm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Andrew Macgregor <andrew at cbbc.murdoch.edu.au>
#
# Copyright Andrew Macgregor, Jo-Ann Stanton, David Green
# Molecular Embryology Group, Anatomy & Structural Biology, University of Otago
# http://meg.otago.ac.nz/
#
# You may distribute this module under the same terms as perl itself
#
# _history
# April 17, 2002 - Initial implementation by Andrew Macgregor
# POD documentation - main docs before the code

=head1 NAME

Bio::Cluster::UniGene - UniGene object

=head1 SYNOPSIS

	use Bio::Cluster::UniGene;
	use Bio::ClusterIO;

	$stream  = Bio::ClusterIO->new('-file' => "Hs.data", 
                                       '-format' => "unigene");
	# note: we quote -format to keep older perl's from complaining.

	while ( my $in = $stream->next_cluster() ) {
		print $in->unigene_id() . "\n";
		while ( my $sequence = $in->next_seq() ) {
			print $sequence->accession_number() . "\n";
		}
       }

=head1 DESCRIPTION

This UniGene object implements the L<Bio::Cluster::UniGeneI> interface
for the representation if UniGene clusters in Bioperl. It is returned
by the L<Bio::ClusterIO> parser for unigene format and contains all
the data associated with one UniGene record.

This class implements several interfaces and hence can be used
wherever instances of such interfaces are expected. In particular, the
interfaces are L<Bio::ClusterI> as the base interface for all cluster
representations, and in addition L<Bio::IdentifiableI> and
L<Bio::DescribableI>.

The following lists the UniGene specific methods that are available
(see below for details). Be aware next_XXX iterators take a snapshot
of the array property when they are first called, and this snapshot is
not reset until the iterator is exhausted. Hence, once called you need
to exhaust the iterator to see any changes that have been made to the
property in the meantime. You will usually want to use the
non-iterator equivalents and loop over the elements yourself.

new() - standard new call

unigene_id() - set/get unigene_id

title() - set/get title (description)

gene() - set/get gene

cytoband() - set/get cytoband

mgi() - set/get mgi

locuslink() - set/get locuslink

homol() - set/get homologene

gnm_terminus() - set/get gnm_terminus

scount() - set/get scount

express() - set/get express, currently takes/returns a reference to an
array of expressed tissues

next_express() - returns the next tissue expression from the expressed
tissue array

chromosome() - set/get chromosome, currently takes/returns a reference
to an array of chromosome lines

next_chromosome() - returns the next chromosome line from the array of
chromosome lines

sts() - set/get sts, currently takes/returns a reference to an array
of sts lines

next_sts() - returns the next sts line from the array of sts lines

txmap() - set/get txmap, currently takes/returns a reference to an
array of txmap lines

next_txmap() - returns the next txmap line from the array of txmap
lines

protsim() - set/get protsim, currently takes/returns a reference to an
array of protsim lines

next_protsim() - returns the next protsim line from the array of
protsim lines

sequences() - set/get sequence, currently takes/returns a reference to
an array of references to seq info

next_seq() - returns a Seq object that currently only contains an
accession number


=head1 Implemented Interfaces

This class implementes the following interfaces.

=over 4

=item Bio::Cluster::UniGeneI

This includes implementing Bio::ClusterI.

=item Bio::IdentifiableI

=item Bio::DescribableI

=item Bio::AnnotatableI

=item Bio::Factory::SequenceStreamI

=back

=head1 FEEDBACK


=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Andrew Macgregor

Email andrew at cbbc.murdoch.edu.au

=head1 CONTRIBUTORS

Hilmar Lapp, hlapp at gmx.net

=head1 APPENDIX


The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a "_".

=cut

# Let the code begin...


package Bio::Cluster::UniGene;
use strict;

use Bio::Annotation::Collection;
use Bio::Annotation::DBLink;
use Bio::Annotation::SimpleValue;
use Bio::Species;
use Bio::Seq::SeqFactory;

use base qw(Bio::Root::Root Bio::Cluster::UniGeneI Bio::IdentifiableI Bio::DescribableI Bio::AnnotatableI Bio::Factory::SequenceStreamI);

my %species_map = (
		   'Aga' => "Anopheles gambiae",
		   'Ame' => "Apis mellifera",
		   'At'  => "Arabidopsis thaliana",
		   'Bmo' => "Bombyx mori",
		   'Bt'  => "Bos taurus",
		   'Cel' => "Caenorhabditis elegans",
		   'Cfa' => "Canine familiaris",
		   'Cin' => "Ciona intestinalis",
		   'Cre' => "Chlamydomonas reinhardtii",
		   'Csa' => "Ciona savignyi",
		   'Csi' => "Citrus sinensis",
		   'Ddi' => "Dictyostelium discoideum",
		   'Dr'  => "Danio rerio",
		   'Dm'  => "Drosophila melanogaster",
		   'Gga' => "Gallus gallus",
		   'Gma' => "Glycine max",
		   'Han' => "Helianthus annus",
		   'Hs'  => "Homo sapiens",
		   'Hma' => "Hydra magnipapillata",
		   'Hv'  => "Hordeum vulgare",
		   'Lco' => "Lotus corniculatus",
		   'Les' => "Lycopersicon esculentum",
		   'Lsa' => "Lactuca sativa",
		   'Mdo' => "Malus x domestica",
                   'Mgr' => "Magnaporthe grisea",
		   'Mm'  => "Mus musculus",
		   'Mtr' => "Medicago truncatula",
                   'Ncr' => "Neurospora crassa",
		   'Oar' => "Ovis aries",
		   'Omy' => "Oncorhynchus mykiss",
		   'Os'  => "Oryza sativa",
		   'Ola' => "Oryzias latipes",
		   'Ppa' => "Physcomitrella patens",
		   'Pta' => "Pinus taeda",
		   'Ptp' => "Populus tremula x Populus tremuloides",
		   'Rn'  => "Rattus norvegicus",
		   'Sbi' => "Sorghum bicolor",
		   'Sma' => "Schistosoma mansoni",
		   'Sof' => "Saccharum officinarum",
		   'Spu' => "Strongylocentrotus purpuratus",
		   'Ssa' => "Salmo salar",
		   'Ssc' => "Sus scrofa",
		   'Str' => "Xenopus tropicalis",
		   'Stu' => "Solanum tuberosum",
		   'Ta'  => "Triticum aestivum",
		   'Tgo' => "Toxoplasma gondii",
                   'Tru' => "Takifugu rubripes",
		   'Vvi' => "Vitis vinifera",
		   'Xl'  => "Xenopus laevis",
		   'Zm'  => "Zea mays",
		   );


=head2 new

 Title   : new
 Usage   : used by ClusterIO
 Returns : a new Bio::Cluster::Unigene object

=cut

sub new {
    # standard new call..
    my($caller,@args) = @_;
    my $self = $caller->SUPER::new(@args);

    my ($ugid,$desc,$mems,$size,$species,$dispid,$id,$ns,$auth,$v,$seqfact) =
	$self->_rearrange([qw(UNIGENE_ID
			      DESCRIPTION
			      MEMBERS
			      SIZE
			      SPECIES
			      DISPLAY_ID
			      OBJECT_ID
			      NAMESPACE
			      AUTHORITY
			      VERSION
			      SEQFACTORY
			      )], @args);

    $self->{'_alphabet'} = 'dna';

    $self->unigene_id($ugid) if $ugid;
    $self->description($desc) if $desc;
    $self->sequences($mems) if $mems;
    $self->size($size) if defined($size);
    $self->display_id($dispid) if $dispid; # overwrites ugid
    $self->object_id($id) if $id;          # overwrites dispid
    $self->namespace($ns || 'UniGene');
    $self->authority($auth || 'NCBI');
    $self->version($v) if defined($v);
    if( ! defined $seqfact ) {
	$seqfact = Bio::Seq::SeqFactory->new
	    (-verbose => $self->verbose(), 
	     -type => 'Bio::Seq::RichSeq');
    }
    $self->sequence_factory($seqfact);
    if( (! $species) && (defined $self->unigene_id() && 
			 $self->unigene_id() =~ /^([A-Za-z]+)\.[0-9]/)) {
	# try set a default one depending on the ID
	$species = $species_map{$1};
    }
    $self->species($species);
    return $self;
}


=head1 L<Bio::Cluster::UniGeneI> methods

=cut

=head2 unigene_id

 Title   : unigene_id
 Usage   : unigene_id();
 Function: Returns the unigene_id associated with the object.
 Example : $id = $unigene->unigene_id or $unigene->unigene_id($id)
 Returns : A string
 Args    : None or an id


=cut

sub unigene_id {
	my ($obj,$value) = @_;
	if( defined $value) {
		$obj->{'unigene_id'} = $value;
	}
	return $obj->{'unigene_id'};
}



=head2 title

 Title   : title
 Usage   : title();
 Function: Returns the title associated with the object.
 Example : $title = $unigene->title or $unigene->title($title)
 Returns : A string
 Args    : None or a title


=cut

sub title {
	my ($obj,$value) = @_;
	if( defined $value) {
		$obj->{'title'} = $value;
	}
	return $obj->{'title'};
}


=head2 gene

 Title   : gene
 Usage   : gene();
 Function: Returns the gene associated with the object.
 Example : $gene = $unigene->gene or $unigene->gene($gene)
 Returns : A string
 Args    : None or a gene


=cut

sub gene {
    my $self = shift;
    return $self->_annotation_value('gene_name', @_);
}


=head2 cytoband

 Title   : cytoband
 Usage   : cytoband();
 Function: Returns the cytoband associated with the object.
 Example : $cytoband = $unigene->cytoband or $unigene->cytoband($cytoband)
 Returns : A string
 Args    : None or a cytoband


=cut

sub cytoband {
    my $self = shift;
    return $self->_annotation_value('cyto_band', @_);
}

=head2 mgi

 Title   : mgi
 Usage   : mgi();
 Function: Returns the mgi associated with the object.
 Example : $mgi = $unigene->mgi or $unigene->mgi($mgi)
 Returns : A string
 Args    : None or a mgi


=cut

sub mgi {
    my $self = shift;
    my $acc;

    if(@_) {
	# purge first
	$self->_remove_dblink('dblink','MGI');
	# then add if a valid value is present
	if($acc = shift) {
	    $self->_annotation_dblink('dblink','MGI',$acc);
	}
    } else {
	($acc) = $self->_annotation_dblink('dblink','MGI');
    }
    return $acc;
}


=head2 locuslink

 Title   : locuslink
 Usage   : locuslink();
 Function: Returns or stores a reference to an array containing locuslink data.
 Returns : An array reference
 Args    : None or an array reference

=cut

sub locuslink {
    my ($self,$ll) = @_;
    
    if($ll) {
	# purge first
	$self->_remove_dblink('dblink','LocusLink');
	# then add as many accessions as are present
	foreach my $acc (@$ll) {
	    $self->_annotation_dblink('dblink','LocusLink',$acc);
	}
    } else {
	my @accs = $self->_annotation_dblink('dblink','LocusLink');
	$ll = [@accs];
    }
    return $ll;
}


=head2 homol

 Title   : homol
 Usage   : homol();
 Function: Returns the homol entry associated with the object.
 Example : $homol = $unigene->homol or $unigene->homol($homol)
 Returns : A string
 Args    : None or a homol entry

=cut

sub homol {
    my $self = shift;
    return $self->_annotation_value('homol', @_);
}


=head2 restr_expr

 Title   : restr_expr
 Usage   : restr_expr();
 Function: Returns the restr_expr entry associated with the object.
 Example : $restr_expr = $unigene->restr_expr or $unigene->restr_expr($restr_expr)
 Returns : A string
 Args    : None or a restr_expr entry

=cut

sub restr_expr {
    my $self = shift;
    return $self->_annotation_value('restr_expr', @_);
}


=head2 gnm_terminus

 Title   : gnm_terminus
 Usage   : gnm_terminus();
 Function: Returns the gnm_terminus associated with the object.
 Example : $gnm_terminus = $unigene->gnm_terminus or 
           $unigene->gnm_terminus($gnm_terminus)
 Returns : A string
 Args    : None or a gnm_terminus

=cut

sub gnm_terminus {
    my $self = shift;
    return $self->_annotation_value('gnm_terminus', @_);
}

=head2 scount

 Title   : scount
 Usage   : scount();
 Function: Returns the scount associated with the object.
 Example : $scount = $unigene->scount or $unigene->scount($scount)
 Returns : A string
 Args    : None or a scount

=cut

sub scount {
	my ($obj,$value) = @_;
	if( defined $value) {
	    $obj->{'scount'} = $value;
	} elsif((! defined($obj->{'scount'})) && defined($obj->sequences())) {
	    $obj->{'scount'} = $obj->size();
	}
	return $obj->{'scount'};
}


=head2 express

 Title   : express
 Usage   : express();
 Function: Returns or stores a reference to an array containing 
           tissue expression data
 Returns : An array reference
 Args    : None or an array reference

=cut

sub express {
    my $self = shift;

    return $self->_annotation_value_ary('expressed',@_);
}


=head2 chromosome

 Title   : chromosome
 Usage   : chromosome();
 Function: Returns or stores a reference to an array containing
           chromosome lines
 Returns : An array reference
 Args    : None or an array reference

=cut

sub chromosome {
    my $self = shift;

    return $self->_annotation_value_ary('chromosome',@_);
 }


=head2 sts

 Title   : sts
 Usage   : sts();
 Function: Returns or stores a reference to an array containing sts lines

 Returns : An array reference
 Args    : None or an array reference

=cut

sub sts {
    my $self = shift;

    return $self->_annotation_value_ary('sts',@_);
}


=head2 txmap

 Title   : txmap
 Usage   : txmap();
 Function: Returns or stores a reference to an array containing txmap lines

 Returns : An array reference
 Args    : None or an array reference

=cut

sub txmap {
    my $self = shift;

    return $self->_annotation_value_ary('txmap',@_);
}


=head2 protsim

 Title   : protsim
 Usage   : protsim();
 Function: Returns or stores a reference to an array containing protsim lines
	   This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub protsim {
    my $self = shift;

    return $self->_annotation_value_ary('protsim',@_);
}


=head2 sequences

 Title   : sequences
 Usage   : sequences();
 Function: Returns or stores a reference to an array containing
           sequence data.

           This is mostly reserved for ClusterIO parsers. You should
           use get_members() for get and add_member()/remove_members()
           for set.

 Returns : An array reference, or undef
 Args    : None or an array reference or undef

=cut

sub sequences {
    my $self = shift;

    return $self->{'members'} = shift if @_;
    return $self->{'members'};
}

=head2 species

 Title   : species
 Usage   : $obj->species($newval)
 Function: Get/set the species object for this Unigene cluster.
 Example : 
 Returns : value of species (a L<Bio::Species> object)
 Args    : on set, new value (a L<Bio::Species> object or 
           the binomial name, or undef, optional)


=cut

sub species{
    my $self = shift;

    if(@_) {
	my $species = shift;
	if($species && (! ref($species))) {
	    my @class = reverse(split(' ',$species));
	    $species = Bio::Species->new(-classification => \@class);
	}
	return $self->{'species'} = $species;
    }
    return $self->{'species'};
}


=head1 L<Bio::ClusterI> methods

=cut

=head2 display_id

 Title   : display_id
 Usage   : 
 Function: Get/set the display name or identifier for the cluster

           This is aliased to unigene_id().

 Returns : a string
 Args    : optional, on set the display ID ( a string)

=cut

sub display_id{
    return shift->unigene_id(@_);
}

=head2 description

 Title   : description
 Usage   : Bio::ClusterI->description("POLYUBIQUITIN")
 Function: get/set for the consensus description of the cluster

           This is aliased to title().

 Returns : the description string 
 Args    : Optional the description string 

=cut

sub description{
    return shift->title(@_);
}

=head2 size

 Title   : size
 Usage   : Bio::ClusterI->size();
 Function: get for the size of the family, 
           calculated from the number of members

           This is aliased to scount().

 Returns : the size of the cluster
 Args    : 

=cut

sub size {
    my $self = shift;

    # hard-wiring the size is allowed if there are no sequences
    return $self->scount(@_) unless defined($self->sequences());
    # but we can't change the number of members through this method
    my $n = scalar(@{$self->sequences()});
    if(@_ && ($n != $_[0])) {
	$self->throw("Cannot change cluster size using size() from $n to ".
		     $_[0]);
    }
    return $n;
}

=head2 cluster_score

 Title   : cluster_score
 Usage   : $cluster ->cluster_score(100);
 Function: get/set for cluster_score which
           represent the score in which the clustering
           algorithm assigns to this cluster.

           For UniGene clusters, there really is no cluster score that
           would come with the data. However, we provide an
           implementation here so that you can score UniGene clusters
           if you want to.

 Returns : a number
 Args    : optionally, on set a number

=cut

sub cluster_score{
    my $self = shift;

    return $self->{'cluster_score'} = shift if @_;
    return $self->{'cluster_score'};
}

=head2 get_members

 Title   : get_members
 Usage   : Bio::ClusterI->get_members(($seq1, $seq2));
 Function: retrieve the members of the family by some criteria

           Will return all members if no criteria are provided.

           At this time this implementation does not support
           specifying criteria and will always return all members.

 Returns : the array of members
 Args    : 

=cut

sub get_members {
    my $self = shift;

    my $mems = $self->sequences() || [];
    # already objects?
    if(@$mems && (ref($mems->[0]) eq "HASH")) {
	# nope, we need to build the object list from scratch
	my @memlist = ();
	while(my $seq = $self->next_seq()) {
	    push(@memlist, $seq);
	}
	# we cache this array of objects as the new member list
	$mems = \@memlist;
	$self->sequences($mems);
    }
    # done
    return @$mems;
}


=head1 Annotatable view at the object properties

=cut

=head2 annotation

 Title   : annotation
 Usage   : $obj->annotation($newval)
 Function: Get/set the L<Bio::AnnotationCollectionI> object for
           this UniGene cluster.

           Many attributes of this class are actually stored within
           the annotation collection object as L<Bio::AnnotationI>
           compliant objects, so you can conveniently access them
           through the same interface as you would e.g. access
           L<Bio::SeqI> annotation properties.

           If you call this method in set mode and replace the
           annotation collection with another one you should know
           exactly what you are doing.

 Example : 
 Returns : a L<Bio::AnnotationCollectionI> compliant object
 Args    : on set, new value (a L<Bio::AnnotationCollectionI> 
           compliant object or undef, optional)


=cut

sub annotation{
    my $self = shift;

    if(@_) {
	return $self->{'annotation'} = shift;
    } elsif(! exists($self->{'annotation'})) {
	$self->{'annotation'} = Bio::Annotation::Collection->new();
    }
    return $self->{'annotation'};
}


=head1 Implementation specific methods

 These are mostly for adding/removing to array properties, and for
 methods with special functionality.

=cut

=head2 add_member

 Title   : add_member
 Usage   :
 Function: Adds a member object to the list of members.
 Example :
 Returns : TRUE if the new member was successfuly added, and FALSE
           otherwise.
 Args    : The member to add.


=cut

sub add_member{
    my ($self,@mems) = @_;

    my $memlist = $self->{'members'} || [];
    # this is an object interface; is the member list already objects?
    if(@$memlist && (ref($memlist->[0]) eq "HASH")) {
	# nope, convert to objects
        $memlist = [$self->get_members()];
    }
    # add new member(s)
    push(@$memlist, @mems);
    # store if we created this array ref ourselves
    $self->sequences($memlist);
    # done
    return 1;
}

=head2 remove_members

 Title   : remove_members
 Usage   :
 Function: Remove the list of members for this cluster such that the
           member list is undefined afterwards (as opposed to zero members).
 Example :
 Returns : the previous list of members
 Args    : none


=cut

sub remove_members{
    my $self = shift;

    my @mems = $self->get_members();
    $self->sequences(undef);
    return @mems;
}


=head2 next_locuslink

 Title   : next_locuslink
 Usage   : next_locuslink();
 Function: Returns the next locuslink from an array referred 
           to using $obj->{'locuslink'}

           If you call this iterator again after it returned undef, it
           will re-cycle through the list of elements. Changes in the
           underlying array property while you loop over this iterator
           will not be reflected until you exhaust the iterator.

 Example : 	while ( my $locuslink = $in->next_locuslink() ) {
				print "$locuslink\n";
			}
 Returns : String
 Args    : None

=cut

sub next_locuslink {
    my ($obj) = @_;

    return $obj->_next_element("ll","locuslink");
}

=head2 next_express

 Title   : next_express
 Usage   : next_express();
 Function: Returns the next tissue from an array referred 
           to using $obj->{'express'}

           If you call this iterator again after it returned undef, it
           will re-cycle through the list of elements. Changes in the
           underlying array property while you loop over this iterator
           will not be reflected until you exhaust the iterator.

 Example : 	while ( my $express = $in->next_express() ) {
				print "$express\n";
			}
 Returns : String
 Args    : None

=cut

sub next_express {
    my ($obj) = @_;

    return $obj->_next_element("express","express");
}


=head2 next_chromosome

 Title   : next_chromosome
 Usage   : next_chromosome();
 Function: Returns the next chromosome line from an array referred
           to using $obj->{'chromosome'}

           If you call this iterator again after it returned undef, it
           will re-cycle through the list of elements. Changes in the
           underlying array property while you loop over this iterator
           will not be reflected until you exhaust the iterator.

 Example : 	while ( my $chromosome = $in->next_chromosome() ) {
				print "$chromosome\n";
			}
 Returns : String
 Args    : None

=cut

sub next_chromosome {
    my ($obj) = @_;

    return $obj->_next_element("chr","chromosome");
}


=head2 next_protsim

 Title   : next_protsim
 Usage   : next_protsim();
 Function: Returns the next protsim line from an array referred 
           to using $obj->{'protsim'}

           If you call this iterator again after it returned undef, it
           will re-cycle through the list of elements. Changes in the
           underlying array property while you loop over this iterator
           will not be reflected until you exhaust the iterator.

 Example : 	while ( my $protsim = $in->next_protsim() ) {
				print "$protsim\n";
			}
 Returns : String
 Args    : None

=cut

sub next_protsim {
    my ($obj) = @_;

    return $obj->_next_element("protsim","protsim");
}


=head2 next_sts

 Title   : next_sts
 Usage   : next_sts();
 Function: Returns the next sts line from an array referred 
           to using $obj->{'sts'}

           If you call this iterator again after it returned undef, it
           will re-cycle through the list of elements. Changes in the
           underlying array property while you loop over this iterator
           will not be reflected until you exhaust the iterator.

 Example : 	while ( my $sts = $in->next_sts() ) {
				print "$sts\n";
			}
 Returns : String
 Args    : None

=cut

sub next_sts {
    my ($obj) = @_;

    return $obj->_next_element("sts","sts");
}


=head2 next_txmap

 Title   : next_txmap
 Usage   : next_txmap();
 Function: Returns the next txmap line from an array 
           referred to using $obj->{'txmap'}

           If you call this iterator again after it returned undef, it
           will re-cycle through the list of elements. Changes in the
           underlying array property while you loop over this iterator
           will not be reflected until you exhaust the iterator.

 Example : 	while ( my $tsmap = $in->next_txmap() ) {
				print "$txmap\n";
			}
 Returns : String
 Args    : None

=cut

sub next_txmap {
    my ($obj) = @_;

    return $obj->_next_element("txmap","txmap");
}

###############################
# private method
#
# args: prefix name for the queue
#       name of the method from which to re-fill
# returns: the next element from that queue, or undef if the queue is empty
###############################
sub _next_element{
    my ($self,$queuename,$meth) = @_;

    $queuename = "_".$queuename."_queue";
    if(! exists($self->{$queuename})) {
	# re-initialize from array of sequence data
	$self->{$queuename} = [@{$self->$meth() }];
    }
    my $queue = $self->{$queuename};
    # is queue exhausted (equivalent to end of stream)?
    if(! @$queue) {
	# yes, remove queue and signal to the caller
	delete $self->{$queuename};
	return;
    }
    return shift(@$queue);
}

=head1 L<Bio::IdentifiableI> methods

=cut

=head2 object_id

 Title   : object_id
 Usage   : $string    = $obj->object_id()
 Function: a string which represents the stable primary identifier
           in this namespace of this object. For DNA sequences this
           is its accession_number, similarly for protein sequences

           This is aliased to unigene_id().

 Returns : A scalar


=cut

sub object_id {
    return shift->unigene_id(@_);
}

=head2 version

 Title   : version
 Usage   : $version    = $obj->version()
 Function: a number which differentiates between versions of
           the same object. Higher numbers are considered to be
           later and more relevant, but a single object described
           the same identifier should represent the same concept

           Unigene clusters usually won't have a version, so this
           will be mostly undefined.

 Returns : A number
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub version {
    my $self = shift;

    return $self->{'version'} = shift if @_;
    return $self->{'version'};
}


=head2 authority

 Title   : authority
 Usage   : $authority    = $obj->authority()
 Function: a string which represents the organisation which
           granted the namespace, written as the DNS name for  
           organisation (eg, wormbase.org)

 Returns : A scalar
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub authority {
    my $self = shift;

    return $self->{'authority'} = shift if @_;
    return $self->{'authority'};
}


=head2 namespace

 Title   : namespace
 Usage   : $string    = $obj->namespace()
 Function: A string representing the name space this identifier
           is valid in, often the database name or the name
           describing the collection 

 Returns : A scalar
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub namespace {
    my $self = shift;

    return $self->{'namespace'} = shift if @_;
    return $self->{'namespace'};
}

=head1 L<Bio::DescribableI> methods

=cut

=head2 display_name

 Title   : display_name
 Usage   : $string    = $obj->display_name()
 Function: A string which is what should be displayed to the user
           the string should have no spaces (ideally, though a cautious
           user of this interface would not assumme this) and should be
           less than thirty characters (though again, double checking 
           this is a good idea)

           This is aliased to unigene_id().

 Returns : A scalar
 Status  : Virtual

=cut

sub display_name {
    return shift->unigene_id(@_);
}


=head2 description()

 Title   : description
 Usage   : $string    = $obj->description()
 Function: A text string suitable for displaying to the user a 
           description. This string is likely to have spaces, but
           should not have any newlines or formatting - just plain
           text. The string should not be greater than 255 characters
           and clients can feel justified at truncating strings at 255
           characters for the purposes of display

           This is already demanded by Bio::ClusterI and hence is
           present anyway.

 Returns : A scalar


=cut


=head1 L<Bio::Factory::SequenceStreamI> methods

=cut

=head2 next_seq

 Title   : next_seq
 Usage   : next_seq();
 Function: Returns the next seq as a Seq object as defined by 
           $seq->sequence_factory(), 
           at present an empty Bio::Seq::RichSeq object with 
           just the accession_number() and pid() set

           This iterator will not exhaust the array of member
           sequences. If you call next_seq() again after it returned
           undef, it will re-cycle through the list of member
           sequences.

 Example :  while ( my $sequence = $in->next_seq() ) {
             print $sequence->accession_number() . "\n";
	    }
 Returns : Bio::PrimarySeqI object
 Args    : None

=cut

sub next_seq {
    my ($obj) = @_;

    if(! exists($obj->{'_seq_queue'})) {
	# re-initialize from array of sequence data
	$obj->{'_seq_queue'} = [@{$obj->sequences()}];
    }
    my $queue = $obj->{'_seq_queue'};
    # is queue exhausted (equivalent to end of stream)?
    if(! @$queue) {
	# yes, remove queue and signal to the caller
	delete $obj->{'_seq_queue'};
	return;
    }
    # no, still data in the queue: get the next one from the queue
    my $seq_h = shift(@$queue);
    # if this is not a simple hash ref, it's an object already, and we'll
    # return just that
    return $seq_h if(ref($seq_h) ne 'HASH');
    # nope, we need to assemble this object from scratch
    #
    # assemble the annotation collection
    my $ac = Bio::Annotation::Collection->new();
    foreach my $k (keys %$seq_h) {
	next if $k =~ /acc|pid|nid|version/;
	my $ann = Bio::Annotation::SimpleValue->new(-tagname => $k,
						    -value   => $seq_h->{$k});
	$ac->add_Annotation($ann);
    }
    # assemble the initialization parameters and create object
    my $seqobj = $obj->sequence_factory->create(
	  -accession_number => $seq_h->{acc},
	  -pid              => $seq_h->{pid},
	  # why does NCBI prepend a 'g' to its own identifiers??
	  -primary_id       => $seq_h->{nid} && $seq_h->{nid} =~ /^g\d+$/ ?
				     substr($seq_h->{nid},1) : $seq_h->{nid},
	  -display_id       => $seq_h->{acc},
	  -seq_version	    => $seq_h->{version},
	  -alphabet         => $obj->{'_alphabet'},
	  -namespace        => $seq_h->{acc} =~ /^NM_/ ? 'RefSeq' : 'GenBank',
	  -authority        => $obj->authority(), # default is NCBI
	  -species          => $obj->species(),
	  -annotation       => $ac
	  );
    return $seqobj;
}

=head2 sequence_factory

 Title   : sequence_factory
 Usage   : $seqio->sequence_factory($seqfactory)
 Function: Get/Set the Bio::Factory::SequenceFactoryI
 Returns : Bio::Factory::SequenceFactoryI
 Args    : [optional] Bio::Factory::SequenceFactoryI


=cut

sub sequence_factory {
    my ($self,$obj) = @_;   
    if( defined $obj ) {
	if( ! ref($obj) || ! $obj->isa('Bio::Factory::SequenceFactoryI') ) {
	    $self->throw("Must provide a valid Bio::Factory::SequenceFactoryI object to ".ref($self)." sequence_factory()");
	}
	$self->{'_seqfactory'} = $obj;
    }
    $self->{'_seqfactory'};
}

=head1 Private methods

=cut

=head2 _annotation_value

 Title   : _annotation_value
 Usage   :
 Function: Private method.
 Example :
 Returns : the value (a string)
 Args    : annotation key (a string)
           on set, annotation value (a string)


=cut

sub _annotation_value{
    my $self = shift;
    my $key = shift;

    my ($ann, $val);
    if(@_) {
	$val = shift;
	if(! defined($val)) {
	    ($ann) = $self->annotation->remove_Annotations($key);
	    return $ann ? $ann->value() : undef;
	}
    }
    ($ann) = $self->annotation->get_Annotations($key);
    if(defined $ann && (! $val)) {
	# get mode and exists
	$val = $ann->value();
    } elsif($val) {
	# set mode
	if(!defined $ann) {
	    $ann = Bio::Annotation::SimpleValue->new(-tagname => $key);
	    $self->annotation->add_Annotation($ann);
	}
	$ann->value($val);
    }
    return $val;
}


=head2 _annotation_value_ary

 Title   : _annotation_value_ary
 Usage   :
 Function: Private method.
 Example :
 Returns : reference to the array of values
 Args    : annotation key (a string)
           on set, reference to an array holding the values


=cut

sub _annotation_value_ary{
    my ($self,$key,$arr) = @_;

    my $ac = $self->annotation;
    if($arr) {
	# purge first
	$ac->remove_Annotations($key);
	# then add as many values as are present
	foreach my $val (@$arr) {
	    my $ann = Bio::Annotation::SimpleValue->new(-value => $val,
							-tagname => $key
							);
	    $ac->add_Annotation($ann);
	}
    } else {
	my @vals = map { $_->value(); } $ac->get_Annotations($key);
	$arr = [@vals];
    }
    return $arr;
}


=head2 _annotation_dblink

 Title   : _annotation_dblink
 Usage   :
 Function: Private method.
 Example :
 Returns : array of accessions for the given database (namespace)
 Args    : annotation key (a string)
           dbname (a string) (optional on get, mandatory on set)
           on set, accession or ID (a string), and version


=cut

sub _annotation_dblink{
    my ($self,$key,$dbname,$acc,$version) = @_;

    if($acc) {
	# set mode -- this is adding here
	my $ann = Bio::Annotation::DBLink->new(-tagname    => $key,
					       -primary_id => $acc,
					       -database   => $dbname,
					       -version    => $version);
	$self->annotation->add_Annotation($ann);
	return 1;
    } else {
	# get mode
	my @anns = $self->annotation->get_Annotations($key);
	# filter out those that don't match the requested database
	if($dbname) {
	    @anns = grep { $_->database() eq $dbname; } @anns;
	}
	return map { $_->primary_id(); } @anns;
    }
}

=head2 _remove_dblink

 Title   : _remove_dblink
 Usage   :
 Function: Private method.
 Example :
 Returns : array of accessions for the given database (namespace)
 Args    : annotation key (a string)
           dbname (a string) (optional)


=cut

sub _remove_dblink{
    my ($self,$key,$dbname) = @_;

    my $ac = $self->annotation();
    my @anns = ();
    if($dbname) {
	foreach my $ann ($ac->remove_Annotations($key)) {
	    if($ann->database() eq $dbname) {
		push(@anns, $ann);
	    } else {
		$ac->add_Annotation($ann);
	    }
	}
    } else {
	@anns = $ac->remove_Annotations($key);
    }
    return map { $_->primary_id(); } @anns;
}


#####################################################################
# aliases for naming consistency or other reasons                   #
#####################################################################

*sequence = \&sequences;

1;
