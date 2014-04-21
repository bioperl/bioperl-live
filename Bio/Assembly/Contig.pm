#
# BioPerl module for Bio::Assembly::Contig
#   Mostly based on Bio::SimpleAlign by Ewan Birney
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Robson Francisco de Souza <rfsouza@citri.iq.usp.br>
#
# Copyright Robson Francisco de Souza
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Assembly::Contig - Perl module to hold and manipulate
                     sequence assembly contigs.

=head1 SYNOPSIS

    # Module loading
    use Bio::Assembly::IO;

    # Assembly loading methods
    $aio = Bio::Assembly::IO->new(-file=>"test.ace.1",
                                  -format=>'phrap');

    $assembly = $aio->next_assembly;
    foreach $contig ($assembly->all_contigs) {
      # do something
    }

    # OR, if you want to build the contig yourself,

    use Bio::Assembly::Contig;
    $c = Bio::Assembly::Contig->new(-id=>"1");

    $ls  = Bio::LocatableSeq->new(-seq=>"ACCG-T",
                                  -id=>"r1",
                                  -alphabet=>'dna');
    $ls2 = Bio::LocatableSeq->new(-seq=>"ACA-CG-T",
                                  -id=>"r2",
                                  -alphabet=>'dna');

    $ls_coord = Bio::SeqFeature::Generic->new(-start=>3,
                                              -end=>8,
                                              -strand=>1);
    $ls2_coord = Bio::SeqFeature::Generic->new(-start=>1,
                                               -end=>8,
                                               -strand=>1);
    $c->add_seq($ls);
    $c->add_seq($ls2);
    $c->set_seq_coord($ls_coord,$ls);
    $c->set_seq_coord($ls2_coord,$ls2);

    $con = Bio::LocatableSeq->new(-seq=>"ACACCG-T",
                                  -alphabet=>'dna');
    $c->set_consensus_sequence($con);

    $l = $c->change_coord('unaligned r2','ungapped consensus',6);
    print "6 in unaligned r2 => $l in ungapped consensus\n";


=head1 DESCRIPTION

A contig is as a set of sequences, locally aligned to each other, so
that every sequence has overlapping regions with at least one sequence
in the contig, such that a continuous of overlapping sequences is
formed, allowing the deduction of a consensus sequence which may be
longer than any of the sequences from which it was deduced.

In this documentation we refer to the overlapping sequences used to
build the contig as "aligned sequences" and to the sequence deduced
from the overlap of aligned sequences as the "consensus". Methods to
deduce the consensus sequence from aligned sequences were not yet
implemented in this module, but its posssible to add a consensus
sequence deduced by other means, e.g, by the assembly program used to
build the alignment.

All aligned sequences in a Bio::Assembly::Contig must be Bio::Assembly::Locatable
objects and have a unique ID. The unique ID restriction is due to the
nature of the module's internal data structures and is also a request
of some assembly programs. If two sequences with the same ID are added
to a contig, the first sequence added is replaced by the second one.

=head2 Coordinate_systems

There are four base coordinate systems in Bio::Assembly::Contig.  When
you need to access contig elements or data that exists on a certain
range or location, you may be specifying coordinates in relation to
different sequences, which may be either the contig consensus or one
of the aligned sequences that were used to do the assembly.

 =========================================================
          Name           | Referenced sequence
 ---------------------------------------------------------
   "gapped consensus"    | Contig (with gaps)
   "ungapped consensus"  | Contig (without gaps)
   "aligned $seqID"      | sequence $seqID (with gaps)
   "unaligned $seqID"    | sequence $seqID (without gaps)
 =========================================================

"gapped consensus" refers to positions in the aligned consensus
sequence, which is the consensus sequence including the gaps inserted
to align it agains the aligned sequences that were used to assemble
the contig. So, its limits are [ 1, (consensus length + number of gaps
in consensus) ]

"ungapped consensus" is a coordinate system based on the consensus
sequence, but excluding consensus gaps. This is just the coordinate
system that you have when considering the consensus sequence alone,
instead of aligned to other sequences.

"aligned $seqID" refers to locations in the sequence $seqID after
alignment of $seqID against the consensus sequence (reverse
complementing the original sequence, if needed).  Coordinate 1 in
"aligned $seqID" is equivalent to the start location (first base) of
$seqID in the consensus sequence, just like if the aligned sequence
$seqID was a feature of the consensus sequence.

"unaligned $seqID" is equivalent to a location in the isolated
sequence, just like you would have when considering the sequence
alone, out of an alignment.  When changing coordinates from "aligned
$seq2" to "unaligned $seq2", if $seq2 was reverse complemented when
included in the alignment, the output coordinates will be reversed to
fit that fact, i.e. 1 will be changed to length($seq2), 2 will be
length($seq)-1 and so on.

An important note: when you change gap coordinates from a gapped
system ("gapped consensus" or "aligned $seqID") to a system that does
not include gaps ("ungapped consensus" or "unaligned $seqID"), the
position returned will be the first location before all gaps
neighboring the input location.

=head2 Feature_collection

Bio::Assembly::Contig stores much information about a contig in a
Bio::Assembly::SeqFeature::Collection object. Relevant information on the
alignment is accessed by selecting features based on their primary
tags (e.g. all features which have a primary tag of the form
'_aligned_coord:$seqID', where $seqID is an aligned sequence ID, are
coordinates for sequences in the contig alignment) and, by using
methods from Bio::Assembly::SeqFeature::Collection, it's possible to select
features by overlap with other features.

We suggest that you use the primary tags of features as identifiers
for feature classes. By convention, features with primary tags
starting with a '_' are generated by modules that populate the contig
data structure and return the contig object, maybe as part of an
assembly object, e.g.  drivers from the Bio::Assembly::IO set.

Features in the features collection may be associated with particular
aligned sequences. To obtain this, you must attach the sequence to the
feature, using attach() seq from Bio::Assembly::SeqFeatureI, before you add the
feature to the feature collection. We also suggest to add the sequence
id to the primary tag, so that is easy to select feature for a
particular sequence.

There is only one feature class that some methods in
Bio::Assembly::Contig expect to find in the feature collection: features
with primary tags of the form '_aligned_coord:$seqID', where $seqID is
the aligned sequence id (like returned by $seq-E<gt>id()). These features
describe the position (in "gapped consensus" coordinates) of aligned
sequences, and the method set_seq_coord() automatically changes a
feature's primary tag to this form whenever the feature is added to
the collection by this method. Only two methods in Bio::Assembly::Contig
will not work unless there are features from this class:
change_coord() and get_seq_coord().

Other feature classes will be automatically available only when
Bio::Assembly::Contig objects are created by a specific module. Such
feature classes are (or should be) documented in the documentation of
the module which create them, to which the user should refer.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

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

=head1 AUTHOR - Robson Francisco de Souza

rfsouza@citri.iq.usp.br

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
package Bio::Assembly::Contig;

use strict;

use Bio::DB::SeqFeature::Store; # isa Bio::SeqFeature::CollectionI
use Bio::Seq::PrimaryQual;      # isa Bio::Seq::QualI

use Scalar::Util qw(weaken);

use base qw(Bio::Root::Root Bio::Align::AlignI);

=head1 Object creator

=head2 new

 Title     : new
 Usage     : my $contig = Bio::Assembly::Contig->new();
 Function  : Creates a new contig object
 Returns   : Bio::Assembly::Contig
 Args      : -id         => unique contig ID
             -source     => string for the sequence assembly program used
             -collection => Bio::SeqFeature::CollectionI instance

=cut

#-----------
sub new {
#-----------
    my ($class, @args) = @_;

    my $self = $class->SUPER::new(@args);

    my ($src, $id, $collection) = $self->_rearrange([qw(SOURCE ID COLLECTION)], @args);
    $src && $self->source($src);
    ($id && $self->id($id)) || ($self->{'_id'} = 'NoName'); # Alignment (contig) name
    ($id && $self->id($id)) || ($self->{'_source'} = 'Unknown'); # Program used to build the contig
    # we need to set up internal hashes first!

    # Bio::SimpleAlign derived fields (check which ones are needed for AlignI compatibility)
    $self->{'_elem'} = {}; # contig elements: aligned sequence objects (keyed by ID)
    $self->{'_order'} = {}; # store sequence order
    # $self->{'start_end_lists'} = {}; # References to entries in {'_seq'}. Keyed by seq ids.
    # $self->{'_dis_name'} = {}; # Display names for each sequence
    $self->{'_symbols'} = {}; # List of symbols

    #Contig specific slots
    $self->{'_consensus_sequence'} = undef;
    $self->{'_consensus_quality'} = undef;
    $self->{'_nof_residues'} = 0;
    $self->{'_nof_seqs'} = 0;
    # $self->{'_nof_segments'} = 0; # Let's not make it heavier than needed by now...
    
    # for cases where SF::Collection is shared between Bio::Assembly::Contig 
    if ($collection) {
        $self->throw("Collection must implement Bio::SeqFeature::CollectionI") unless $collection->isa('Bio::SeqFeature::CollectionI');
        $self->{'_sfc'} = $collection;
    } else {
        $self->{'_sfc'} = Bio::DB::SeqFeature::Store->new(
            -adaptor           => 'memory',
            -index_subfeatures => 1,
        );
    }

    # Assembly specifics
    $self->{'_assembly'} = undef; # Bio::Assembly::Scaffold the contig belongs to
    $self->{'_strand'} = 0; # Reverse (-1) or forward (1), if contig is in a scaffold. 0 otherwise
    $self->{'_neighbor_start'} = undef; # Neighbor Bio::Assembly::Contig
    $self->{'_neighbor_end'}   = undef; # Neighbor Bio::Assembly::Contig

    return $self; # success - we hope!
}

=head1 Assembly related methods

These methods exist to enable adding information about possible
relations among contigs, e.g. when you already have a scaffold for
your assembly, describing the ordering of contigs in the final
assembly, but no sequences covering the gaps between neighboring
contigs.

=head2 source

 Title     : source
 Usage     : $contig->source($program);
 Function  : Get/Set program used to build this contig
 Returns   : string
 Argument  : [optional] string

=cut

sub source {
    my $self = shift;
    my $source = shift;

    $self->{'_source'} = $source if (defined $source);
    return $self->{'_source'};
}

=head2 assembly

 Title     : assembly
 Usage     : $contig->assembly($assembly);
 Function  : Get/Set assembly object for this contig
 Returns   : a Bio::Assembly::Scaffold object
 Argument  : a Bio::Assembly::Scaffold object

=cut

sub assembly {
    my $self = shift;
    my $assembly = shift;

    $self->throw("Using non Bio::Assembly::Scaffold object when assign contig to assembly")
    if (defined $assembly && ! $assembly->isa("Bio::Assembly::Scaffold"));
    # We create a circular reference to a Scaffold object. It is made weak
    # to prevent memory leaks.
    $self->{'_assembly'} = $assembly if (defined $assembly); 
    weaken($self->{'_assembly'});

    return $self->{'_assembly'};
}

=head2 strand

 Title     : strand
 Usage     : $contig->strand($num);
 Function  : Get/Set contig orientation in a scaffold/assembly.
             Its equivalent to the strand property of sequence
             objects and sets whether the contig consensus should
             be reversed and complemented before being added to a
             scaffold or assembly.
 Returns   : integer
 Argument  : 1 if orientaion is forward, -1 if reverse and
             0 if none

=cut

sub strand {
    my $self = shift;
    my $ori = shift;

    if (defined $ori) {
        $self->throw("Contig strand must be either 1, -1 or 0")
            unless $ori == 1 || $ori == 0 || $ori == -1;
        $self->{'_strand'} = $ori;
    }

    return $self->{'_strand'};
}

=head2 upstream_neighbor

 Title     : upstream_neighbor
 Usage     : $contig->upstream_neighbor($contig);
 Function  : Get/Set a contig neighbor for the current contig when
             building a scaffold. The upstream neighbor is
             located before $contig first base
 Returns   : nothing
 Argument  : Bio::Assembly::Contig

=cut

sub upstream_neighbor {
    my $self = shift;
    my $ref = shift;

    $self->throw("Trying to assign a non Bio::Assembly::Contig object to upstream contig")
        if (defined $ref && ! $ref->isa("Bio::Assembly::Contig"));

    $self->{'_neighbor_start'} = $ref if (defined $ref);
    return $self->{'_neighbor_start'};
}

=head2 downstream_neighbor

 Title     : downstream_neighbor
 Usage     : $contig->downstream_neighbor($num);
 Function  : Get/Set a contig neighbor for the current contig when
             building a scaffold. The downstream neighbor is
             located after $contig last base
 Returns   : nothing
 Argument  : Bio::Assembly::Contig

=cut

sub downstream_neighbor {
    my $self = shift;
    my $ref = shift;

    $self->throw("Trying to assign a non Bio::Assembly::Contig object to downstream contig")
        if (defined $ref && ! $ref->isa("Bio::Assembly::Contig"));
    $self->{'_neighbor_end'} = $ref if (defined $ref);
    return $self->{'_neighbor_end'};
}

=head1 Contig feature collection methods

=head2 add_features

 Title     : add_features
 Usage     : $contig->add_features($feat,$flag)
 Function  :

             Add an array of features to the contig feature
             collection. The consensus sequence may be attached to the
             added feature, if $flag is set to 1. If $flag is 0 and
             the feature attached to one of the contig aligned
             sequences, the feature is registered as an aligned
             sequence feature. If $flag is 0 and the feature is not
             attched to any sequence in the contig, the feature is
             simply added to the feature collection and no attachment
             or registration is made.

             Note: You must attach aligned sequences to their features
             prior to calling add_features, otherwise you won't be
             able to access the feature through get_seq_feat_by_tag()
             method.

 Returns   : number of features added.
 Argument  :
             $feat : A reference to an array of Bio::SeqFeatureI
             $flag : boolean - true if consensus sequence object
                     should be attached to this feature, false if
                     no consensus attachment should be made.
                     Default: false.

=cut

sub add_features {
    my ($self, $args, $flag) = @_;

    # Adding shortcuts for aligned sequence features
    $flag = 0 unless (defined $flag);
    if ($flag && defined $self->{'_consensus_sequence'}) {
        foreach my $feat (@$args) {
            next if (defined $feat->seq);
            $feat->attach_seq($self->{'_consensus_sequence'});
        }
    } elsif (!$flag) { # Register aligned sequence features
        foreach my $feat (@$args) {
            if (my $seq = $feat->entire_seq()) {
                my $seqID = $seq->id() || $seq->display_id || $seq->primary_id;
                $self->warn("Adding contig feature attached to unknown sequence $seqID!")
                    unless (exists $self->{'_elem'}{$seqID});
                my $tag = $feat->primary_tag;
                $self->{'_elem'}{$seqID}{'_feat'}{$tag} = $feat;
            }
        }
    }

    # Add feature to feature collection
    my $nof_added = $self->get_features_collection->add_features($args);

    return $nof_added;
}

=head2 remove_features

 Title     : remove_features
 Usage     : $contig->remove_features(@feat)
 Function  : Remove an array of contig features
 Returns   : true if successful
 Argument  : An array of Bio::SeqFeature::Generic (Bio::SeqFeatureI)

=cut

sub remove_features {
    my ($self, @args) = @_;

    # Removing shortcuts for aligned sequence features
    for my $feat (@args) {
        if (my $seq = $feat->entire_seq()) {
            my $seqID = $seq->id || $seq->display_id || $seq->primary_id;
            my $tag = $feat->primary_tag;
            $tag =~ s/:$seqID$/$1/g;
            delete( $self->{'_elem'}{$seqID}{'_feat'}{$tag} )
                if (exists $self->{'_elem'}{$seqID}{'_feat'}{$tag} &&
                $self->{'_elem'}{$seqID}{'_feat'}{$tag} eq $feat);
        }
    }
   
    # Removing Bio::SeqFeature objects
    return $self->get_features_collection->delete(@args);
}

=head2 get_features_collection

 Title     : get_features_collection
 Usage     : $contig->get_features_collection()
 Function  : Get the collection of all contig features and seqfeatures
 Returns   : Bio::DB::SeqFeature::Store (Bio::SeqFeature::CollectionI)
 Argument  : none

=cut

sub get_features_collection {
    my $self = shift;
    return $self->{'_sfc'};
}

=head2 remove_features_collection

 Title     : remove_features_collection
 Usage     : $contig->remove_features_collection()
 Function  : Remove the collection of all contig features. It is useful
             to save some memory (when contig features are not needed).
 Returns   : none
 Argument  : none

=cut

sub remove_features_collection {
    my $self = shift;
    # Removing shortcuts for aligned sequence features
    for my $seqID (keys %{$self->{'_elem'}}) {
        delete $self->{'_elem'}{$seqID};
    }
    # Removing Bio::SeqFeature::Collection features
    $self->{'_sfc'} = {};
    return;
}

=head1 Coordinate system's related methods

See L<Coordinate_Systems> above.

=head2 change_coord

 Title     : change_coord
 Usage     : $contig->change_coord($in,$out,$query)
 Function  :

             Change coordinate system for $query.  This method
             transforms locations between coordinate systems described
             in section "Coordinate Systems" of this document.

             Note: this method will throw an exception when changing
             coordinates between "ungapped consensus" and other
             systems if consensus sequence was not set. It will also
             throw exceptions when changing coordinates among aligned
             sequence, either with or without gaps, and other systems
             if sequence locations were not set with set_seq_coord().

 Returns   : integer
 Argument  :
             $in    : [string]  input coordinate system
             $out   : [string]  output coordinate system
             $query : [integer] a position in a sequence

=cut

sub change_coord {
    my $self     = shift;
    my $type_in  = shift;
    my $type_out = shift;
    my $query    = shift;

    # Parsing arguments
    # Loading read objects (these calls will throw exceptions whether $read_in or
    # $read_out is not found
    my ($read_in,$read_out) = (undef,undef);
    my $in_ID  = ( split(' ',$type_in)  )[1];
    my $out_ID = ( split(' ',$type_out) )[1];

    if ($in_ID  ne 'consensus') {
        $read_in  = $self->get_seq_coord( $self->get_seq_by_name($in_ID)  );
        $self->throw("Can't change coordinates without sequence location for $in_ID")
            unless (defined $read_in);
    }
    if ($out_ID ne 'consensus') {
        $read_out = $self->get_seq_coord( $self->get_seq_by_name($out_ID) );
        $self->throw("Can't change coordinates without sequence location for $out_ID")
            unless (defined $read_out);
    }

    # Performing transformation between coordinates
    SWITCH1: {

        # Transformations between contig padded and contig unpadded
        (($type_in eq 'gapped consensus') && ($type_out eq 'ungapped consensus')) && do {
            $self->throw("Can't use ungapped consensus coordinates without a consensus sequence")
                unless (defined $self->{'_consensus_sequence'});
            $query = &_padded_unpadded($self->{'_consensus_gaps'}, $query);
            last SWITCH1;
        };
        (($type_in eq 'ungapped consensus') && ($type_out eq 'gapped consensus')) && do {
            $self->throw("Can't use ungapped consensus coordinates without a consensus sequence")
                unless (defined $self->{'_consensus_sequence'});
            $query = &_unpadded_padded($self->{'_consensus_gaps'},$query);
            last SWITCH1;
        };

        # Transformations between contig (padded) and read (padded)
        (($type_in  eq 'gapped consensus') &&
        ($type_out =~ /^aligned /) && defined($read_out)) && do {
            $query = $query - $read_out->start() + 1;
            last SWITCH1;
        };
        (($type_in =~ /^aligned /) && defined($read_in) &&
        ($type_out  eq 'gapped consensus')) && do {
            $query = $query + $read_in->start() - 1;
            last SWITCH1;
        };

        # Transformations between contig (unpadded) and read (padded)
        (($type_in eq 'ungapped consensus') &&
        ($type_out =~ /^aligned /) && defined($read_out)) && do {
            $query = $self->change_coord('ungapped consensus','gapped consensus',$query);
            $query = $self->change_coord('gapped consensus',"aligned $out_ID",$query);
            last SWITCH1;
        };
        (($type_in =~ /^aligned /) && defined($read_in) &&
        ($type_out eq 'ungapped consensus')) && do {
            $query = $self->change_coord("aligned $in_ID",'gapped consensus',$query);
            $query = $self->change_coord('gapped consensus','ungapped consensus',$query);
            last SWITCH1;
        };

        # Transformations between seq $read_in padded and seq $read_out padded
        (defined($read_in)  && ($type_in  =~ /^aligned /)  &&
        defined($read_out) && ($type_out =~ /^aligned /)) && do {
            $query = $self->change_coord("aligned $in_ID",'gapped consensus',$query);
            $query = $self->change_coord('gapped consensus',"aligned $out_ID",$query);
            last SWITCH1;
        };

        # Transformations between seq $read_in padded and seq $read_out unpadded
        (defined($read_in)  && ($type_in  =~ /^aligned /)    &&
        defined($read_out) && ($type_out =~ /^unaligned /)) && do {
            if ($read_in ne $read_out) {
                $query = $self->change_coord("aligned $in_ID",'gapped consensus',$query);
                $query = $self->change_coord('gapped consensus',"aligned $out_ID",$query);
            }
            my $list_out = $self->{'_elem'}{$out_ID}{'_gaps'};
            $query = &_padded_unpadded($list_out,$query);
            # Changing read orientation if read was reverse complemented when aligned
            if ($read_out->strand == -1) {
                my ($length) = $read_out->length();
                $length = $length - &_nof_gaps($list_out,$length);
                $query  = $length - $query + 1;
            }
            last SWITCH1;
        };
        (defined($read_in)  && ($type_in  =~ /^unaligned /) &&
        defined($read_out) && ($type_out =~ /^aligned /))  && do {
            my $list_in = $self->{'_elem'}{$in_ID}{'_gaps'};
            # Changing read orientation if read was reverse complemented when aligned
            if ($read_in->strand == -1) {
                my ($length) = $read_in->length();
                $length = $length - &_nof_gaps($list_in,$length);
                $query  = $length - $query + 1;
            }
            $query = &_unpadded_padded($list_in,$query);
            if ($read_in ne $read_out) {
                $query = $self->change_coord("aligned $in_ID",'gapped consensus',$query);
                $query = $self->change_coord('gapped consensus',"aligned $out_ID",$query);
            }
            last SWITCH1;
        };

        # Transformations between seq $read_in unpadded and seq $read_out unpadded
        (defined($read_in)  && ($type_in  =~ /^unaligned /)    &&
        defined($read_out) && ($type_out =~ /^unaligned /)) && do {
            $query = $self->change_coord("unaligned $in_ID","aligned $out_ID",$query);
            $query = $self->change_coord("aligned $out_ID","unaligned $out_ID",$query);
            last SWITCH1;
        };

        # Transformations between contig (padded) and read (unpadded)
        (($type_in eq 'gapped consensus') &&
        ($type_out =~ /^unaligned /) && defined($read_out)) && do {
            $query = $self->change_coord('gapped consensus',"aligned $out_ID",$query);
            $query = $self->change_coord("aligned $out_ID","unaligned $out_ID",$query);
            last SWITCH1;
        };
        (($type_in =~ /^unaligned /) && defined($read_in) &&
        ($type_out eq 'gapped consensus')) && do {
            $query = $self->change_coord("unaligned $in_ID","aligned $in_ID",$query);
            $query = $self->change_coord("aligned $in_ID",'gapped consensus',$query);
            last SWITCH1;
        };

        # Transformations between contig (unpadded) and read (unpadded)
        (($type_in eq 'ungapped consensus') &&
        ($type_out =~ /^unaligned /) && defined($read_out)) && do {
            $query = $self->change_coord('ungapped consensus','gapped consensus',$query);
            $query = $self->change_coord('gapped consensus',"unaligned $out_ID",$query);
            last SWITCH1;
        };
        (($type_in =~ /^unaligned /) && defined($read_in) &&
        ($type_out eq 'ungapped consensus')) && do {
            $query = $self->change_coord("unaligned $in_ID",'gapped consensus',$query);
            $query = $self->change_coord('gapped consensus','ungapped consensus',$query);
            last SWITCH1;
        };

        $self->throw("Unknow coordinate system. Args: $type_in, $type_out.");
        $query = undef; # If a coordinate systems just requested is unknown
    }

    return $query;
}

=head2 get_seq_coord

 Title     : get_seq_coord
 Usage     : $contig->get_seq_coord($seq);
 Function  : Get "gapped consensus" location for aligned sequence
 Returns   : Bio::SeqFeature::Generic for coordinates or undef.
             A warning is printed if sequence coordinates were not set.
 Argument  : Bio::LocatableSeq object

=cut

sub get_seq_coord {
    my ($self,$seq) = @_;

    if( !ref $seq || ! $seq->isa('Bio::LocatableSeq') ) {
        $self->throw("$seq is not a Bio::LocatableSeq");
    }
    my $seqID = $seq->id() || $seq->display_id || $seq->primary_id;

    unless (exists( $self->{'_elem'}{$seqID} )) {
        $self->warn("No such sequence ($seqID) in contig ".$self->id);
        return;
    }
    unless (exists( $self->{'_elem'}{$seqID}{'_feat'}{"_aligned_coord:$seqID"} )) {
        # $self->warn("Chad. Location not set for sequence ($seqID) in contig ".$self->id);
        return;
    }

    return $self->{'_elem'}{$seqID}{'_feat'}{"_aligned_coord:$seqID"};
}

=head2 set_seq_coord

 Title     : set_seq_coord
 Usage     : $contig->set_seq_coord($feat,$seq);
 Function  :

             Set "gapped consensus" location for an aligned
             sequence. If the sequence was previously added using
             add_seq, its coordinates are changed/set.  Otherwise,
             add_seq is called and the sequence is added to the
             contig.

 Returns   : Bio::SeqFeature::Generic for old coordinates or undef.
 Argument  :
             $feat  : a Bio::SeqFeature::Generic object
                      representing a location for the
                      aligned sequence, in "gapped
                      consensus" coordinates.

             Note: the original feature primary tag will
                   be lost.

             $seq   : a Bio::LocatableSeq object

=cut

sub set_seq_coord {
    my ($self,$feat,$seq) = @_;

    if( !ref $seq || ! $seq->isa('Bio::LocatableSeq') ) {
        $self->throw("Unable to process non locatable sequences [".ref($seq)."]");
    }

    # Complaining about inadequate feature object
    $self->throw("Coordinates must be a Bio::SeqFeature::Generic object!")
        unless ( $feat->isa("Bio::SeqFeature::Generic") );
    $self->throw("Sequence coordinates must have an end!")
        unless (defined $feat->end);
    $self->throw("Sequence coordinates must have a start!")
        unless (defined $feat->start);

    my $seqID = $seq->id() || $seq->display_id || $seq->primary_id;
    if ( exists( $self->{'_elem'}{$seqID} ) &&
         exists( $self->{'_elem'}{$seqID}{'_seq'} ) &&
         defined( $self->{'_elem'}{$seqID}{'_seq'} ) &&
         ($seq ne $self->{'_elem'}{$seqID}{'_seq'}) ) {
        $self->warn("Replacing sequence $seqID\n");
        $self->remove_seq($self->{'_elem'}{$seqID}{'_seq'});
        $self->remove_features($feat);
    }

    # Add new sequence and Bio::Generic::SeqFeature
    $self->add_seq($seq);

    $feat->add_tag_value('contig',$self->id) unless ( $feat->has_tag('contig') );
    $feat->primary_tag("_aligned_coord");
    $feat->source_tag($seqID);
    $feat->attach_seq($seq);

    $self->{'_elem'}{$seqID}{'_feat'}{"_aligned_coord:$seqID"} = $feat;
    $self->add_features([ $feat ]);
}

=head1 Bio::Assembly::Contig consensus methods

=head2 set_consensus_sequence

 Title     : set_consensus_sequence
 Usage     : $contig->set_consensus_sequence($seq)
 Function  : Set the consensus sequence object for this contig
 Returns   : consensus length
 Argument  : Bio::LocatableSeq

=cut

sub set_consensus_sequence {
    my $self = shift;
    my $seq  = shift;

    $self->throw("Consensus sequence must be a Bio::LocatableSeq!")
        unless ($seq->isa("Bio::LocatableSeq"));

    $self->{'_consensus_gaps'} = []; # Consensus Gap registry
    $self->_register_gaps( $seq->seq, $self->{'_consensus_gaps'} );
    $self->{'_consensus_sequence'} = $seq;

    $seq->start(1);
    $seq->end($seq->_ungapped_len);

    my $con_len = $seq->length;

    return $con_len;
}

=head2 set_consensus_quality

 Title     : set_consensus_quality
 Usage     : $contig->set_consensus_quality($qual)
 Function  : Set the quality object for consensus sequence
 Returns   : nothing
 Argument  : Bio::Seq::QualI object

=cut

sub set_consensus_quality {
    my ($self, $qual) = @_;

    $self->throw("Consensus quality must be a Bio::Seq::QualI object!")
        unless ( $qual->isa("Bio::Seq::QualI") );

    $self->throw("Consensus quality can't be added before you set the consensus sequence!")
        unless (defined $self->{'_consensus_sequence'});

    $self->{'_consensus_quality'} = $qual;
}

=head2 get_consensus_length

 Title     : get_consensus_length
 Usage     : $contig->get_consensus_length()
 Function  : Get consensus sequence length
 Returns   : integer
 Argument  : none

=cut

sub get_consensus_length {
    my $self = shift;

    return $self->{'_consensus_sequence'}->length();
}

=head2 get_consensus_sequence

 Title     : get_consensus_sequence
 Usage     : $contig->get_consensus_sequence()
 Function  : Get a reference to the consensus sequence object
             for this contig
 Returns   : Bio::SeqI object
 Argument  : none

=cut

sub get_consensus_sequence {
    my ($self, @args) = @_;

    return $self->{'_consensus_sequence'};
}

=head2 get_consensus_quality

 Title     : get_consensus_quality
 Usage     : $contig->get_consensus_quality()
 Function  : Get a reference to the consensus quality object
             for this contig.
 Returns   : A Bio::Seq::QualI object
 Argument  : none

=cut

sub get_consensus_quality {
    my ($self, @args) = @_;

    return $self->{'_consensus_quality'};
}

=head1 Bio::Assembly::Contig aligned sequences methods

=head2 set_seq_qual

 Title     : set_seq_qual
 Usage     : $contig->set_seq_qual($seq,$qual);
 Function  : Adds quality to an aligned sequence.
 Returns   : nothing
 Argument  : a Bio::LocatableSeq object and
             a Bio::Seq::QualI object

See L<Bio::LocatableSeq> for more information.

=cut

sub set_seq_qual {
    my ($self,$seq,$qual) = @_;

    if( !ref $seq || ! $seq->isa('Bio::LocatableSeq') ) {
        $self->throw("Unable to process non locatable sequences [".ref($seq)."]");
    }
    my $seqID = $seq->id() || $seq->display_id || $seq->primary_id;

    $self->throw("Consensus quality must be a Bio::Seq::QualI object!")
        unless ( $qual->isa("Bio::Seq::QualI") );
    $self->throw("Use add_seq first: aligned sequence qualities can't be added before you load the sequence!")
        unless (exists $self->{'_elem'}{$seqID}{'_seq'});
    $self->throw("Use set_seq_coord first: aligned sequence qualities can't be added before you add coordinates for the sequence!") unless (defined( $self->get_seq_coord($seq) ));

    # Adding gaps to quality object
    my $sequence = $self->{'_elem'}{$seqID}{'_seq'}->seq();
    my $tmp = $qual->qual();
    @{$tmp} = reverse(@{$tmp}) if ($self->get_seq_coord($seq)->strand() == -1);
    my @quality  = ();
    my $previous = 0;
    my $next     = 0;
    my $i = 0; my $j = 0;
    while ($i <= $#{$tmp}) {
        # IF base is a gap, quality is the average for neighbouring sites
        if ($j > $i && substr($sequence,$j,1) eq '-') {
            $previous = $tmp->[$i-1] unless ($i == 0);
            if ($i < $#{$tmp}) {
                $next = $tmp->[$i+1];
            } else {
                $next = 0;
            }
            push(@quality,int( ($previous+$next)/2 ));
        } else {
            push(@quality,$tmp->[$i]);
            $i++;
        }
        $j++;
    }

    $self->{'_elem'}{$seqID}{'_qual'} = Bio::Seq::PrimaryQual->new(
        -qual=>join(" ",@quality), -id=>$seqID );
}

=head2 get_seq_ids

 Title     : get_seq_ids
 Usage     : $contig->get_seq_ids( -start => $start,
                                   -end   => $end,
                                   -type  => "gapped A0QR67B08.b" );
 Function  : Get list of sequence IDs overlapping interval [$start, $end]
             The default interval is [1,$contig->length]
             Default coordinate system is "gapped contig"
 Returns   : An array
 Argument  : A hash with optional elements:
             -start : consensus subsequence start
             -end   : consensus subsequence end
             -type  : the coordinate system type for $start and $end arguments
                      Coordinate system available are:
                      "gapped consensus"   : consensus coordinates with gaps
                      "ungapped consensus" : consensus coordinates without gaps
                      "aligned $ReadID"    : read $ReadID coordinates with gaps
                      "unaligned $ReadID"  : read $ReadID coordinates without gaps


=cut

sub get_seq_ids {
    my ($self, @args) = @_;

    my ($type, $start, $end) = $self->_rearrange([qw(TYPE START END)], @args);

    my @list;
    if (defined($start) && defined($end)) {
        if (defined($type) && ($type ne 'gapped consensus')) {
            $start = $self->change_coord($type,'gapped consensus',$start);
            $end   = $self->change_coord($type,'gapped consensus',$end);
        }
        @list = $self->get_features_collection->features(
           -type         => '_aligned_coord', # primary tag
           -start        => $start,
           -end          => $end,
           #-contain     => 0,
           #-strandmatch => 'ignore',
        );
        @list = map { $_->entire_seq->id } @list;
    } else {
        # Entire aligned sequences list
        @list = map { $self->{'_order'}{$_} } sort { $a cmp $b } keys %{ $self->{'_order'} };
    }

    return @list;
}

=head2 get_seq_feat_by_tag

 Title     : get_seq_feat_by_tag
 Usage     : $seq = $contig->get_seq_feat_by_tag($seq,"_aligned_coord:$seqID")
 Function  : Get a sequence feature based on its primary_tag.
 Returns   : a Bio::SeqFeature object
 Argument  : a Bio::LocatableSeq and a string (feature primary tag)

=cut

sub get_seq_feat_by_tag {
    my ($self,$seq,$tag) = @_;

    if( !ref $seq || ! $seq->isa('Bio::LocatableSeq') ) {
        $self->throw("Unable to process non locatable sequences [".ref($seq)."]");
    }
    my $seqID = $seq->id || $seq->display_id || $seq->primary_id;

    return $self->{'_elem'}{$seqID}{'_feat'}{$tag};
}

=head2 get_seq_by_name

 Title     : get_seq_by_name
 Usage     : $seq = $contig->get_seq_by_name('Seq1')
 Function  : Gets a sequence based on its id.
 Returns   : a Bio::LocatableSeq object
             undef if name is not found
 Argument  : string

=cut

sub get_seq_by_name {
    my $self = shift;
    my ($seqID) = @_;

    unless (exists $self->{'_elem'}{$seqID}{'_seq'}) {
        $self->throw("Could not find sequence $seqID in contig ".$self->id);
    return;
    }

    return $self->{'_elem'}{$seqID}{'_seq'};
}

=head2 get_qual_by_name

 Title     : get_qual_by_name
 Usage     : $seq = $contig->get_qual_by_name('Seq1')
 Function  :

             Gets Bio::Seq::QualI object for a sequence
             through its id ( as given by $qual->id() ).

 Returns   : a Bio::Seq::QualI object.
             undef if name is not found
 Argument  : string

=cut

sub get_qual_by_name {
    my $self = shift;
    my ($seqID) = @_;

    unless (exists $self->{'_elem'}{$seqID}{'_qual'}) {
        $self->warn("Could not find quality for $seqID in contig!");
        return;
    }

    return $self->{'_elem'}{$seqID}{'_qual'};
}

=head1 Bio::Align::AlignI compatible methods

=head2 Modifier methods

These methods modify the MSE by adding, removing or shuffling complete
sequences.

=head2 add_seq

 Title     : add_seq
 Usage     : $contig->add_seq($newseq);
 Function  :

             Adds a sequence to the contig. *Does*
             *not* align it - just adds it to the
             hashes.

 Returns   : nothing
 Argument  : a Bio::LocatableSeq object

See L<Bio::LocatableSeq> for more information.

=cut

sub add_seq {
    my $self = shift;
    my $seq = shift;

    if( !ref $seq || ! $seq->isa('Bio::LocatableSeq') ) {
        $self->throw("Unable to process non locatable sequences [".ref($seq)."]");
    }

    my $seqID = $seq->id() || $seq->display_id || $seq->primary_id;
    $self->{'_elem'}{$seqID} = {} unless (exists $self->{'_elem'}{$seqID});

    if (exists( $self->{'_elem'}{$seqID}{'_seq'} ) &&
    ($seq eq $self->{'_elem'}{$seqID}{'_seq'}) ) {
        $self->warn("Adding sequence $seqID, which has already been added");
    }

    # Our locatable sequences are always considered to be complete sequences
    $seq->start(1);
    $seq->end($seq->_ungapped_len);
    
    my $alphabet = $seq->alphabet;
    
    $alphabet = lc($alphabet) if defined $alphabet;
    
    $self->warn("Adding non-nucleotidic sequence ".$seqID)
        if (!$alphabet || ($alphabet ne 'dna' && $alphabet ne 'rna'));

    # build the symbol list for this sequence,
    # will prune out the gap and missing/match chars
    # when actually asked for the symbol list in the
    # symbol_chars
    if (defined $seq->seq) {
        map { $self->{'_symbols'}->{$_} = 1; } split(//,$seq->seq);
    } else {
        $self->{'_symbols'} = {};
    }

    my $seq_no = ++$self->{'_nof_seqs'};

    if (ref( $self->{'_elem'}{$seqID}{'_seq'} )) {
        $self->warn("Replacing one sequence [$seqID]\n");
    } else {
        #print STDERR "Assigning $seqID to $order\n";
        $self->{'_order'}->{$seq_no} = $seqID;
        # $self->{'_start_end_lists'}->{$id} = []
        # unless(exists $self->{'_start_end_lists'}->{$id});
        # push @{$self->{'_start_end_lists'}->{$id}}, $seq;
    }

    $self->{'_elem'}{$seqID}{'_seq'}  = $seq;
    $self->{'_elem'}{$seqID}{'_feat'} = {};
    $self->{'_elem'}{$seqID}{'_gaps'} = [];
    my $dbref = $self->{'_elem'}{$seqID}{'_gaps'};
    my $nofgaps = $self->_register_gaps($seq->seq,$dbref);

    # Updating residue count
    $self->{'_nof_residues'} += $seq->length - $nofgaps;

    return 1;
}

=head2 remove_seq

 Title     : remove_seq
 Usage     : $contig->remove_seq($seq);
 Function  : Removes a single sequence from a contig
 Returns   : 1 on success, 0 otherwise
 Argument  : a Bio::LocatableSeq object

=cut

sub remove_seq {
    my ($self,$seq) = @_;

    if( !ref $seq || ! $seq->isa('Bio::LocatableSeq') ) {
        $self->throw("Unable to process non locatable sequences [".ref($seq)."]");
    }

    my $seqID = $seq->id() || $seq->display_id || $seq->primary_id;
    unless (exists $self->{'_elem'}{$seqID} ) {
        $self->warn("No sequence named $seqID  [$seq]");
        return 0;
    }

    # Updating residue count
    $self->{'_nof_residues'} -= $seq->length() +
    &_nof_gaps( $self->{'_elem'}{$seqID}{'_gaps'}, $seq->length );
    
    # Update number of sequences
    $self->{'_nof_seqs'}--; 
    
    # Update order of sequences (order starts at 1)
    my $max_order = $self->{'_nof_seqs'} + 1;
    my $target_order = $max_order + 1;
    for (my $order = 1 ; $order <= $max_order ; $order++) {
      if ($self->{'_order'}->{$order} eq $seqID) {
        # Found the wanted sequence order
        $target_order = $order;
      }
      if ($order > $target_order) {
        # Decrement this sequence order by one order
        $self->{'_order'}->{$order-1} = $self->{'_order'}->{$order};
      }
      if ($order == $max_order) {
        # Remove last order
        delete $self->{'_order'}->{$order};
      }
    }

    # Remove all references to features of this sequence
    my @feats = ();
    for my $tag (keys %{ $self->{'_elem'}{$seqID}{'_feat'} }) {
        push(@feats, $self->{'_elem'}{$seqID}{'_feat'}{$tag});
    }
    $self->{'_sfc'}->remove_features(\@feats);
    delete $self->{'_elem'}{$seqID};

    return 1;
}

=head2 purge

 Title   : purge
 Usage   : $contig->purge(0.7);
 Function:

           Removes sequences above whatever %id.

           This function will grind on large alignments. Beware!
           (perhaps not ideally implemented)

 Example :
 Returns : An array of the removed sequences
 Argument:


=cut

sub purge {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 sort_alphabetically

 Title     : sort_alphabetically
 Usage     : $contig->sort_alphabetically
 Function  :

             Changes the order of the alignemnt to alphabetical on name
             followed by numerical by number.

 Returns   :
 Argument  :

=cut

sub sort_alphabetically {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 Sequence selection methods

Methods returning one or more sequences objects.

=head2 each_seq

 Title     : each_seq
 Usage     : foreach $seq ( $contig->each_seq() )
 Function  : Gets an array of Seq objects from the alignment
 Returns   : an array
 Argument  :

=cut

sub each_seq {
    my ($self) = @_;

    my (@arr,$seqID);

    foreach $seqID ( map { $self->{'_order'}{$_} } sort { $a <=> $b } keys %{$self->{'_order'}} ) {
        push(@arr,$self->{'_elem'}{$seqID}{'_seq'});
    }

    return @arr;
}

=head2 each_alphabetically

 Title     : each_alphabetically
 Usage     : foreach $seq ( $contig->each_alphabetically() )
 Function  :

             Returns an array of sequence object sorted alphabetically
             by name and then by start point.
             Does not change the order of the alignment

 Returns   :
 Argument  :

=cut

sub each_alphabetically {
    my($self) = @_;
    $self->throw_not_implemented();
}

=head2 each_seq_with_id

 Title     : each_seq_with_id
 Usage     : foreach $seq ( $contig->each_seq_with_id() )
 Function  :

             Gets an array of Seq objects from the
             alignment, the contents being those sequences
             with the given name (there may be more than one)

 Returns   : an array
 Argument  : a seq name

=cut

sub each_seq_with_id {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 get_seq_by_pos

 Title     : get_seq_by_pos
 Usage     : $seq = $contig->get_seq_by_pos(3)
 Function  :

             Gets a sequence based on its position in the alignment.
             Numbering starts from 1.  Sequence positions larger than
             num_sequences() will thow an error.

 Returns   : a Bio::LocatableSeq object
 Argument  : positive integer for the sequence osition

=cut

sub get_seq_by_pos {
    my $self = shift;
    my ($pos) = @_;

    $self->throw("Sequence position has to be a positive integer, not [$pos]")
        unless $pos =~ /^\d+$/ and $pos > 0;
    $self->throw("No sequence at position [$pos]")
        unless $pos <= $self->num_sequences ;

    my $seqID = $self->{'_order'}->{--$pos};
    return $self->{'_elem'}{$seqID}{'_seq'};
}

=head2 Create new alignments

The result of these methods are horizontal or vertical subsets of the
current MSE.

=head2 select

 Title     : select
 Usage     : $contig2 = $contig->select(1, 3) # three first sequences
 Function  :

             Creates a new alignment from a continuous subset of
             sequences.  Numbering starts from 1.  Sequence positions
             larger than num_sequences() will thow an error.

 Returns   : a Bio::Assembly::Contig object
 Argument  : positive integer for the first sequence
             positive integer for the last sequence to include (optional)

=cut

sub select {
    my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 select_noncont

 Title     : select_noncont
 Usage     : $contig2 = $contig->select_noncont(1, 3) # first and 3rd sequences
 Function  :

             Creates a new alignment from a subset of
             sequences.  Numbering starts from 1.  Sequence positions
             larger than num_sequences() will throw an error.

 Returns   : a Bio::Assembly::Contig object
 Args      : array of integers for the sequences

=cut

sub select_noncont {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 slice

 Title     : slice
 Usage     : $contig2 = $contig->slice(20, 30)
 Function  :

             Creates a slice from the alignment inclusive of start and
             end columns.  Sequences with no residues in the slice are
             excluded from the new alignment and a warning is printed.
             Slice beyond the length of the sequence does not do
             padding.

 Returns   : a Bio::Assembly::Contig object
 Argument  : positive integer for start column
             positive integer for end column

=cut

sub slice {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 Change sequences within the MSE

These methods affect characters in all sequences without changeing the
alignment.


=head2 map_chars

 Title     : map_chars
 Usage     : $contig->map_chars('\.','-')
 Function  :

             Does a s/$arg1/$arg2/ on the sequences. Useful for gap
             characters

             Notice that the from (arg1) is interpretted as a regex,
             so be careful about quoting meta characters (eg
             $contig->map_chars('.','-') wont do what you want)

 Returns   :
 Argument  : 'from' rexexp
             'to' string

=cut

sub map_chars {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 uppercase

 Title     : uppercase()
 Usage     : $contig->uppercase()
 Function  : Sets all the sequences to uppercase
 Returns   :
 Argument  :

=cut

sub uppercase {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 match_line

 Title    : match_line()
 Usage    : $contig->match_line()
 Function : Generates a match line - much like consensus string
            except that a line indicating the '*' for a match.
 Argument : (optional) Match line characters ('*' by default)
            (optional) Strong match char (':' by default)
            (optional) Weak match char ('.' by default)

=cut

sub match_line {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 match

 Title     : match()
 Usage     : $contig->match()
 Function  :

             Goes through all columns and changes residues that are
             identical to residue in first sequence to match '.'
             character. Sets match_char.

             USE WITH CARE: Most MSE formats do not support match
             characters in sequences, so this is mostly for output
             only. NEXUS format (Bio::AlignIO::nexus) can handle
             it.

 Returns   : 1
 Argument  : a match character, optional, defaults to '.'

=cut

sub match {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 unmatch

 Title     : unmatch()
 Usage     : $contig->unmatch()
 Function  :

             Undoes the effect of method match. Unsets match_char.

 Returns   : 1
 Argument  : a match character, optional, defaults to '.'

=cut

sub unmatch {
    my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 MSE attibutes

Methods for setting and reading the MSE attributes.

Note that the methods defining character semantics depend on the user
to set them sensibly.  They are needed only by certain input/output
methods. Unset them by setting to an empty string ('').

=head2 id

 Title     : id
 Usage     : $contig->id("Ig")
 Function  : Gets/sets the id field of the alignment
 Returns   : An id string
 Argument  : An id string (optional)

=cut

sub id {
    my ($self, $contig_name) = @_;

    if (defined( $contig_name )) {
        $self->{'_id'} = $contig_name;
    }

    return $self->{'_id'};
}

=head2 missing_char

 Title     : missing_char
 Usage     : $contig->missing_char("?")
 Function  : Gets/sets the missing_char attribute of the alignment
             It is generally recommended to set it to 'n' or 'N'
             for nucleotides and to 'X' for protein.
 Returns   : An missing_char string,
 Argument  : An missing_char string (optional)

=cut

sub missing_char {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 match_char

 Title     : match_char
 Usage     : $contig->match_char('.')
 Function  : Gets/sets the match_char attribute of the alignment
 Returns   : An match_char string,
 Argument  : An match_char string (optional)

=cut

sub match_char {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 gap_char

 Title     : gap_char
 Usage     : $contig->gap_char('-')
 Function  : Gets/sets the gap_char attribute of the alignment
 Returns   : An gap_char string, defaults to '-'
 Argument  : An gap_char string (optional)

=cut

sub gap_char {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 symbol_chars

 Title   : symbol_chars
 Usage   : my @symbolchars = $contig->symbol_chars;
 Function: Returns all the seen symbols (other than gaps)
 Returns : array of characters that are the seen symbols
 Argument: boolean to include the gap/missing/match characters

=cut

sub symbol_chars{
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 Alignment descriptors

These read only methods describe the MSE in various ways.


=head2 consensus_string

 Title     : consensus_string
 Usage     : $str = $contig->consensus_string($threshold_percent)
 Function  : Makes a strict consensus
 Returns   :
 Argument  : Optional threshold ranging from 0 to 100.
             The consensus residue has to appear at least threshold %
             of the sequences at a given location, otherwise a '?'
             character will be placed at that location.
             (Default value = 0%)

=cut

sub consensus_string {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 consensus_iupac

 Title     : consensus_iupac
 Usage     : $str = $contig->consensus_iupac()
 Function  :

             Makes a consensus using IUPAC ambiguity codes from DNA
             and RNA. The output is in upper case except when gaps in
             a column force output to be in lower case.

             Note that if your alignment sequences contain a lot of
             IUPAC ambiquity codes you often have to manually set
             alphabet.  Bio::PrimarySeq::_guess_type thinks they
             indicate a protein sequence.

 Returns   : consensus string
 Argument  : none
 Throws    : on protein sequences


=cut

sub consensus_iupac {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 is_flush

 Title     : is_flush
 Usage     : if( $contig->is_flush() )
           :
           :
 Function  : Tells you whether the alignment
           : is flush, ie all of the same length
           :
           :
 Returns   : 1 or 0
 Argument  :

=cut

sub is_flush {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 length

 Title     : length()
 Usage     : $len = $contig->length()
 Function  : Returns the maximum length of the alignment.
             To be sure the alignment is a block, use is_flush
 Returns   :
 Argument  :

=cut

sub length {
    my ($self) = @_;

    $self->throw_not_implemented();
}

=head2 maxname_length

 Title     : maxname_length
 Usage     : $contig->maxname_length()
 Function  :

             Gets the maximum length of the displayname in the
             alignment. Used in writing out various MSE formats.

 Returns   : integer
 Argument  :

=cut

sub maxname_length {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 num_residues

 Title     : num_residues
 Usage     : $no = $contig->num_residues
 Function  : number of residues in total in the alignment
 Returns   : integer
 Argument  :
 Note      : replaces no_residues

=cut

sub num_residues {
    my ($self) = @_;

    return $self->{'_nof_residues'};
}

=head2 num_sequences

 Title     : num_sequences
 Usage     : $depth = $contig->num_sequences
 Function  : number of sequence in the sequence alignment
 Returns   : integer
 Argument  : None
 Note      : replaces no_sequences

=cut

sub num_sequences {
    my ($self) = @_;

    return scalar( keys %{ $self->{'_elem'} } );
}

=head2 percentage_identity

 Title   : percentage_identity
 Usage   : $id = $contig->percentage_identity
 Function: The function calculates the percentage identity of the alignment
 Returns : The percentage identity of the alignment (as defined by the
                             implementation)
 Argument: None

=cut

sub percentage_identity{
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 overall_percentage_identity

 Title   : percentage_identity
 Usage   : $id = $contig->percentage_identity
 Function: The function calculates the percentage identity of
           the conserved columns
 Returns : The percentage identity of the conserved columns
 Args    : None

=cut

sub overall_percentage_identity{
    my ($self) = @_;
    $self->throw_not_implemented();
}


=head2 average_percentage_identity

 Title   : average_percentage_identity
 Usage   : $id = $contig->average_percentage_identity
 Function: The function uses a fast method to calculate the average
           percentage identity of the alignment
 Returns : The average percentage identity of the alignment
 Args    : None

=cut

sub average_percentage_identity {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 Alignment positions

Methods to map a sequence position into an alignment column and back.
column_from_residue_number() does the former. The latter is really a
property of the sequence object and can done using
L<Bio::LocatableSeq::location_from_column>:

    # select somehow a sequence from the alignment, e.g.
    my $seq = $contig->get_seq_by_pos(1);
    #$loc is undef or Bio::LocationI object
    my $loc = $seq->location_from_column(5);


=head2 column_from_residue_number

 Title   : column_from_residue_number
 Usage   : $col = $contig->column_from_residue_number( $seqname, $resnumber)
 Function:

           This function gives the position in the alignment
           (i.e. column number) of the given residue number in the
           sequence with the given name. For example, for the
           alignment

           Seq1/91-97 AC..DEF.GH
           Seq2/24-30 ACGG.RTY..
           Seq3/43-51 AC.DDEFGHI

           column_from_residue_number( "Seq1", 94 ) returns 5.
           column_from_residue_number( "Seq2", 25 ) returns 2.
           column_from_residue_number( "Seq3", 50 ) returns 9.

           An exception is thrown if the residue number would lie
           outside the length of the aligment
           (e.g. column_from_residue_number( "Seq2", 22 )

      Note: If the the parent sequence is represented by more than
      one alignment sequence and the residue number is present in
      them, this method finds only the first one.

 Returns : A column number for the position in the alignment of the
           given residue in the given sequence (1 = first column)
 Args    : A sequence id/name (not a name/start-end)
           A residue number in the whole sequence (not just that
           segment of it in the alignment)

=cut

sub column_from_residue_number {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 Sequence names

Methods to manipulate the display name. The default name based on the
sequence id and subsequence positions can be overridden in various
ways.

=head2 displayname

 Title     : displayname
 Usage     : $contig->displayname("Ig", "IgA")
 Function  : Gets/sets the display name of a sequence in the alignment
           :
 Returns   : A display name string
 Argument  : name of the sequence
             displayname of the sequence (optional)

=cut

sub displayname { # Do nothing
}

=head2 set_displayname_count

 Title     : set_displayname_count
 Usage     : $contig->set_displayname_count
 Function  :

             Sets the names to be name_# where # is the number of
             times this name has been used.

 Returns   : None
 Argument  : None

=cut

sub set_displayname_count {
    my ($self) = @_;
    $self->throw_not_implemented();
}

=head2 set_displayname_flat

 Title     : set_displayname_flat
 Usage     : $contig->set_displayname_flat()
 Function  : Makes all the sequences be displayed as just their name,
             not name/start-end
 Returns   : 1
 Argument  : None

=cut

sub set_displayname_flat { # Do nothing!
}

=head2 set_displayname_normal

 Title     : set_displayname_normal
 Usage     : $contig->set_displayname_normal()
 Function  : Makes all the sequences be displayed as name/start-end
 Returns   : None
 Argument  : None

=cut

sub set_displayname_normal { # Do nothing!
}

=head1 Internal Methods

=head2 _binary_search

 Title     : _binary_search
 Usage     : _binary_search($list,$query)
 Function  :

             Find a number in a sorted list of numbers.  Return values
             may be on or two integers. One positive integer or zero
             (>=0) is the index of the element that stores the queried
             value.  Two positive integers (or zero and another
             number) are the indexes of elements among which the
             queried value should be placed. Negative single values
             mean:

             -1: $query is smaller than smallest element in list
             -2: $query is greater than greatest element in list

 Returns   : array of integers
 Argument  :
             $list  : array reference
             $query : integer

=cut

sub _binary_search {
    my $list   = shift;
    my $query  = shift;
    #
    # If there is only one element in list
    if (!$#{$list} && ($query == $list->[0])) { return (0) }
    # If there are others...
    my $start = 0;
    my $end   = $#{$list};
    (&_compare($query,$list->[$start]) == 0) && do { return ($start) };
    (&_compare($query,$list->[$end])   == 0) && do { return ($end) };
    (&_compare($query,$list->[$start])  < 0) && do { return (-1) };
    (&_compare($query,$list->[$end])    > 0) && do { return (-2) };
    my $middle = 0;
    while ($end - $start > 1) {
        $middle = int(($end+$middle)/2);
        (&_compare($query,$list->[$middle]) == 0) && return ($middle);
        (&_compare($query,$list->[$middle]) <  0) && do { $end   = $middle ; $middle = 0; next };
        $start = $middle; # If &_compare() > 0, move region beggining
    }
    return ($start,$end);
}

=head2 _compare

    Title   : _compare
    Usage   : _compare($arg1,$arg2)
    Function: Perform numeric or string comparisons
    Returns : integer (0, 1 or -1)
    Args    : values to be compared

=cut

sub _compare {
    my $arg1 = shift;
    my $arg2 = shift;
    #
    if (($arg1 =~ /^\d+$/) && ($arg2 =~ /^\d+$/)) { return $arg1 <=> $arg2 }
    else { return $arg1 cmp $arg2 }
}

=head2 _nof_gaps

    Title   : _nof_gaps
    Usage   : _nof_gaps($array_ref, $query)
    Function: number of gaps found before position $query
    Returns : integer
    Args    :
              $array_ref : gap registry reference
              $query     : [integer] a position in a sequence

=cut

#' emacs...
sub _nof_gaps {
    my $list  = shift;
    my $query = shift;
    # If there are no gaps in this contig
    return 0 unless (defined($list) && scalar(@{$list}));
    # Locate query index in gap list (if any)
    my @index = &_binary_search($list,$query);
    # If after all alignments, correct using total number of align
    if ($index[0] == -2) { $query = scalar(@{$list}) }
    # If before any alignment, return 0
    elsif ($index[0] == -1) { $query = 0 }
    elsif ($index[0] >= 0) {
    # If query is between alignments, translate coordinates
    if ($#index > 0) { $query = $index[0] + 1 }
    # If query sits upon an alignment, do another correction
    elsif ($#index == 0) { $query = $index[0] }
    }
    #
    return $query;
}

=head2 _padded_unpadded

    Title   : _padded_unpadded
    Usage   : _padded_unpadded($array_ref, $query)
    Function:

              Returns a coordinate corresponding to
              position $query after gaps were
              removed from a sequence.

    Returns : integer
    Args    :
              $array_ref : reference to this gap registry
              $query     : [integer] coordionate to change

=cut

sub _padded_unpadded {
    my $list  = shift;
    my $query = shift;

    my $align = &_nof_gaps($list,$query);
    $query-- if (defined($list->[$align]) && ($list->[$align] == $query));
    $query = $query - $align;
    #
    return $query;
}

=head2 _unpadded_padded

    Title   : _unpadded_padded
    Usage   : _unpadded_padded($array_ref, $query)
    Function:

              Returns the value corresponding to
              ungapped position $query when gaps are
              counted as valid sites in a sequence

    Returns :
    Args    : $array_ref = a reference to this sequence's gap registry
              $query = [integer] location to change

=cut

#'
sub _unpadded_padded {
    my $list  = shift;
    my $query = shift;

    my $align  = &_nof_gaps($list,$query);
    $query = $query + $align;
    my $new_align = &_nof_gaps($list,$query);
    while ($new_align - $align > 0) {
        $query = $query + $new_align - $align;
        $align  = $new_align;
        $new_align = &_nof_gaps($list,$query);
    }
    # If current position is also a align, look for the first upstream base
    while (defined($list->[$align]) && ($list->[$align] == $query)) {
        $query++; $align++;
    }
    #
    return $query;
}

=head2 _register_gaps

    Title   : _register_gaps
    Usage   : $self->_register_gaps($seq, $array_ref)
    Function: stores gap locations for a sequence
    Returns : number of gaps found
    Args    :
              $seq       : sequence string
              $array_ref : a reference to an array,
                           where gap locations will
                           be stored

=cut

sub _register_gaps {
    my $self     = shift;
    my $sequence = shift;
    my $dbref    = shift;

    $self->throw("Not an aligned sequence string to register gaps")
        if (ref($sequence));

    $self->throw("Not an array reference for gap registry")
        unless (ref($dbref) eq 'ARRAY');

    # Registering alignments
    @{$dbref} = (); # Cleaning registry
    if (defined $sequence) {
        my $i = -1;
        while(1) {
            $i = index($sequence,"-",$i+1);
            last if ($i == -1);
            push(@{$dbref},$i+1);
        }
    } else {
        # $self->warn("Found undefined sequence while registering gaps");
        return 0;
    }

    return scalar(@{$dbref});
}

=head1 Deprecated methods

=cut

=head2 no_residues

 Title     : no_residues
 Usage     : $no = $ali->no_residues
 Function  : number of residues in total in the alignment
 Returns   : integer
 Argument  :
 Note      : deprecated in favor of num_residues() 

=cut

sub no_residues {
    my $self = shift;
    $self->deprecated(-warn_version  => 1.0069,
                      -throw_version => 1.0075);
    $self->num_residues(@_);
}

=head2 no_sequences

 Title     : no_sequences
 Usage     : $depth = $ali->no_sequences
 Function  : number of sequence in the sequence alignment
 Returns   : integer
 Argument  :
 Note      : deprecated in favor of num_sequences()

=cut

sub no_sequences {
    my $self = shift;
    $self->deprecated(-warn_version => 1.0069,
                      -throw_version => 1.0075);
    $self->num_sequences(@_);
}

1;
