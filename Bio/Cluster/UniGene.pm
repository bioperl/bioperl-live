# $Id$
#
# BioPerl module for Bio::Cluster::UniGene.pm
#
# Cared for by Andrew Macgregor <andrew@anatomy.otago.ac.nz>
#
# Copyright Andrew Macgregor, Jo-Ann Stanton, David Green
# Molecular Embryology Group, Anatomy & Structural Biology, University of Otago
# http://anatomy.otago.ac.nz/meg
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
(see below for details). Be aware that with the exception of next_seq,
all next_XXX methods exhaust the array over which they iterate. You
will usually want to use their non-iterator equivalents and loop over
the elements yourself.

new() - standard new call

unigene_id() - set/get unigene_id

title() - set/get title (description)

gene() - set/get gene

cytoband() - set/get cytoband

mgi() - set/get mgi

locuslink() - set/get locuslink

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

=head1 AUTHOR - Andrew Macgregor

Email andrew@anatomy.otago.ac.nz

=head1 CONTRIBUTORS

Hilmar Lapp, hlapp at gmx.net

=head1 APPENDIX


The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a "_".

=cut

# Let the code begin...


package Bio::Cluster::UniGene;
use vars qw(@ISA $VERSION);
use strict;


use Bio::Root::Root;
use Bio::Annotation::Collection;
use Bio::Seq;
use Bio::Seq::RichSeq;
use Bio::Factory::SequenceStreamI;
use Bio::Seq::SeqFactory;
use Bio::Cluster::UniGeneI;

$VERSION = '1.0';
@ISA = qw(Bio::Root::Root Bio::Cluster::UniGeneI
	Bio::IdentifiableI Bio::DescribableI
	Bio::Factory::SequenceStreamI);


=head2 new

 Title   : new
 Usage   : used by ClusterIO
 Returns : a new Bio::Cluster::Unigene object

=cut

sub new {
    # standard new call..
    my($caller,@args) = @_;
    my $self = $caller->SUPER::new(@args);

    my ($ugid,$desc,$mems,$dispid,$id,$ns,$auth,$v,$seqfact) =
	$self->_rearrange([qw(UNIGENE_ID
			      DESCRIPTION
			      MEMBERS
			      DISPLAY_ID
			      OBJECT_ID
			      NAMESPACE
			      AUTHORITY
			      VERSION
			      SEQFACTORY
			      )], @args);

    $self->{'_alphabet'} = 'dna';
    $self->{'sequences'} = [];

    $self->unigene_id($ugid) if $ugid;
    $self->description($desc) if $desc;
    $self->sequences($mems) if $mems;
    $self->display_id($dispid) if $dispid; # overwrites ugid
    $self->object_id($id) if $id;          # overwrites dispid
    $self->namespace($ns || 'UniGene');
    $self->authority($auth || 'NCBI');
    $self->version($v) if defined($v);
    if( ! defined $seqfact ) {
	$seqfact = new Bio::Seq::SeqFactory
	    (-verbose => $self->verbose(), 
	     -type => 'Bio::Seq::RichSeq');
    }
    $self->sequence_factory($seqfact);
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
    return $self->_annotation_value('mgi', @_);
}


=head2 locuslink

 Title   : locuslink
 Usage   : locuslink();
 Function: Returns or stores a reference to an array containing locuslink data.
 Returns : An array reference
 Args    : None or an array reference

=cut

sub locuslink {
	my ($obj,$value) = @_;
	if( defined $value) {
		$obj->{'locuslink'} = $value;
	}
	return $obj->{'locuslink'};
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
	my ($obj,$value) = @_;
	if( defined $value) {
		$obj->{'express'} = $value;
	}
	return $obj->{'express'};
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
	my ($obj,$value) = @_;
	if( defined $value) {
		$obj->{'chromosome'} = $value;
	}
	return $obj->{'chromosome'};
}


=head2 sts

 Title   : sts
 Usage   : sts();
 Function: Returns or stores a reference to an array containing sts lines
 	   This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub sts {
	my ($obj,$value) = @_;
	if( defined $value) {
		$obj->{'sts'} = $value;
	}
	return $obj->{'sts'};
}


=head2 txmap

 Title   : txmap
 Usage   : txmap();
 Function: Returns or stores a reference to an array containing txmap lines
	   This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub txmap {
	my ($obj,$value) = @_;
	if( defined $value) {
		$obj->{'txmap'} = $value;
	}
	return $obj->{'txmap'};
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
	my ($obj,$value) = @_;
	if( defined $value) {
		$obj->{'protsim'} = $value;
	}
	return $obj->{'protsim'};
}


=head2 sequences

 Title   : sequences
 Usage   : sequences();
 Function: Returns or stores a reference to an array containing sequence data
 	   This should really only be used by ClusterIO, not directly
 Returns : An array reference
 Args    : None or an array reference

=cut

sub sequences {
	my ($obj,$value) = @_;
	if( defined $value) {
		$obj->{'sequences'} = $value;
	}
	return $obj->{'sequences'};
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
    return shift->scount(@_);
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

    my @seqs = ();
    while(my $seq = $self->next_seq()) {
	push(@seqs, $seq);
    }
    return @seqs;
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
           compliant objects. If you call this method in set mode and
           replace the annotation collection with another one you
           should know exactly what you are doing.

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

 These are mostly for adding/removing to array properties, and for methods with special functionality.

=cut

=head2 next_locuslink

 Title   : next_locuslink
 Usage   : next_locuslink();
 Function: Returns the next locuslink from an array referred 
           to using $obj->{'locuslink'}

           Note that this iterator will exhaust the array!

 Example : 	while ( my $locuslink = $in->next_locuslink() ) {
				print "$locuslink\n";
			}
 Returns : String
 Args    : None

=cut

sub next_locuslink {
	my ($obj) = @_;
	shift @{$obj->{'locuslink'}};
}

=head2 next_express

 Title   : next_express
 Usage   : next_express();
 Function: Returns the next tissue from an array referred 
           to using $obj->{'express'}

           Note that this iterator will exhaust the array!

 Example : 	while ( my $express = $in->next_express() ) {
				print "$express\n";
			}
 Returns : String
 Args    : None

=cut

sub next_express {
	my ($obj) = @_;
	shift @{$obj->{'express'}};
}


=head2 next_chromosome

 Title   : next_chromosome
 Usage   : next_chromosome();
 Function: Returns the next chromosome line from an array referred
           to using $obj->{'chromosome'}

           Note that this iterator will exhaust the array!

 Example : 	while ( my $chromosome = $in->next_chromosome() ) {
				print "$chromosome\n";
			}
 Returns : String
 Args    : None

=cut

sub next_chromosome {
	my ($obj) = @_;
	shift @{$obj->{'chromosome'}};
}


=head2 next_protsim

 Title   : next_protsim
 Usage   : next_protsim();
 Function: Returns the next protsim line from an array referred 
           to using $obj->{'protsim'}

           Note that this iterator will exhaust the array!

 Example : 	while ( my $protsim = $in->next_protsim() ) {
				print "$protsim\n";
			}
 Returns : String
 Args    : None

=cut

sub next_protsim {
	my ($obj) = @_;
	shift @{$obj->{'protsim'}};
}


=head2 next_sts

 Title   : next_sts
 Usage   : next_sts();
 Function: Returns the next sts line from an array referred 
           to using $obj->{'sts'}

           Note that this iterator will exhaust the array!

 Example : 	while ( my $sts = $in->next_sts() ) {
				print "$sts\n";
			}
 Returns : String
 Args    : None

=cut

sub next_sts {
	my ($obj) = @_;
	shift @{$obj->{'sts'}};
}


=head2 next_txmap

 Title   : next_txmap
 Usage   : next_txmap();
 Function: Returns the next txmap line from an array 
           referred to using $obj->{'txmap'}

           Note that this iterator will exhaust the array!

 Example : 	while ( my $tsmap = $in->next_txmap() ) {
				print "$txmap\n";
			}
 Returns : String
 Args    : None

=cut

sub next_txmap {
	my ($obj) = @_;
	shift @{$obj->{'txmap'}};
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

           Unigene clusters usually won''t have a version, so this
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


=head2 description

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

    if(! exists($obj->{'_queue'})) {
	# re-initialize from array of sequence data
	$obj->{'_queue'} = [@{$obj->{'sequences'}}];
    }
    my $queue = $obj->{'_queue'};
    # is queue exhausted (equivalent to end of stream)?
    if(! @$queue) {
	# yes, remove queue and signal to the caller
	delete $obj->{'_queue'};
	return undef;
    }
    # no, still data in the queue
    my $seq_h = shift(@$queue);
    my $seqobj = $obj->sequence_factory->create
	( -accession_number => $seq_h->{acc},
	  -pid              => $seq_h->{pid},
	  -display_id       => $seq_h->{acc},
	  -alphabet         => $obj->{'_alphabet'},
	  -namespace        => 'GenBank',
	  -authority        => $obj->authority(), # default is NCBI
	  -desc => join(' ',
			map { uc($_) ."=". $seq_h->{$_}} sort keys %{$seq_h})
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
 Function:
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
    if($ann && (! $val)) {
	# get mode and exists
	$val = $ann->value();
    } elsif($val) {
	# set mode
	if(! $ann) {
	    $ann = Bio::Annotation::SimpleValue->new(-tagname => $key);
	    $self->annotation->add_Annotation($ann);
	}
	$ann->value($val);
    }
    return $val;
}


#####################################################################
# aliases for naming consistency or other reasons                   #
#####################################################################

*sequence = \&sequences;

1;
