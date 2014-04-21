#
# BioPerl module for Bio::SeqFeature::Generic
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Generic - Generic SeqFeature

=head1 SYNOPSIS

   $feat = Bio::SeqFeature::Generic->new( 
            -start        => 10, 
            -end          => 100,
            -strand       => -1, 
            -primary      => 'repeat', # -primary_tag is a synonym
            -source_tag   => 'repeatmasker',
            -display_name => 'alu family',
            -score        => 1000,
            -tag          => { new => 1,
                               author => 'someone',
                               sillytag => 'this is silly!' } );

   $feat = Bio::SeqFeature::Generic->new( -gff_string => $string );
   # if you want explicitly GFF1
   $feat = Bio::SeqFeature::Generic->new( -gff1_string => $string );

   # add it to an annotated sequence

   $annseq->add_SeqFeature($feat);

=head1 DESCRIPTION

Bio::SeqFeature::Generic is a generic implementation for the
Bio::SeqFeatureI interface, providing a simple object to provide all
the information for a feature on a sequence.

For many Features, this is all you will need to use (for example, this
is fine for Repeats in DNA sequence or Domains in protein
sequence). For other features, which have more structure, this is a
good base class to extend using inheritence to have new things: this
is what is done in the L<Bio::SeqFeature::Gene>,
L<Bio::SeqFeature::Transcript> and L<Bio::SeqFeature::Exon>, which provide
well coordinated classes to represent genes on DNA sequence (for
example, you can get the protein sequence out from a transcript
class).

For many Features, you want to add some piece of information, for
example a common one is that this feature is 'new' whereas other
features are 'old'.  The tag system, which here is implemented using a
hash can be used here.  You can use the tag system to extend the
L<Bio::SeqFeature::Generic> programmatically: that is, you know that you have
read in more information into the tag 'mytag' which you can then
retrieve. This means you do not need to know how to write inherited
Perl to provide more complex information on a feature, and/or, if you
do know but you do not want to write a new class every time you need
some extra piece of information, you can use the tag system to easily
store and then retrieve information.

The tag system can be written in/out of GFF format, and also into EMBL
format via the L<Bio::SeqIO> system

=head1 Implemented Interfaces

This class implements the following interfaces.

=over 4

=item L<Bio::SeqFeatureI>

Note that this includes implementing Bio::RangeI.

=item L<Bio::AnnotatableI>

=item L<Bio::FeatureHolderI>

Features held by a feature are essentially sub-features.

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via 
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Ewan Birney

Ewan Birney E<lt>birney@sanger.ac.ukE<gt>

=head1 DEVELOPERS

This class has been written with an eye out for inheritance. The fields
the actual object hash are:

   _gsf_tag_hash  = reference to a hash for the tags
   _gsf_sub_array = reference to an array for subfeatures

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Generic;
use strict;

use Bio::Annotation::Collection;
use Bio::Location::Simple;
use Bio::Location::Split;
use Bio::Tools::GFF;
#use Tie::IxHash;

use base qw(Bio::Root::Root Bio::SeqFeatureI Bio::FeatureHolderI Bio::AnnotatableI);

sub new {
    my ( $caller, @args) = @_;
    my ($self) = $caller->SUPER::new(@args);
    $self->_register_for_cleanup(\&cleanup_generic);
    $self->{'_parse_h'}       = {};
    $self->{'_gsf_tag_hash'}  = {};

    # bulk-set attributes
    $self->set_attributes(@args);

    # done - we hope
    return $self;
}

=head2 set_attributes

 Title   : set_attributes
 Usage   :
 Function: Sets a whole array of parameters at once.
 Example :
 Returns : none
 Args    : Named parameters, in the form as they would otherwise be passed
           to new(). Currently recognized are:

                    -start          start position
                    -end            end position
                    -strand         strand
                    -phase          the phase of the feature (0..2)
                    -primary_tag    primary tag 
                    -primary        (synonym for -primary_tag)
                    -source         source tag
                    -frame          frame
                    -score          score value
                    -tag            a reference to a tag/value hash
                    -gff_string     GFF v.2 string to initialize from
                    -gff1_string    GFF v.1 string to initialize from
                    -seq_id         the display name of the sequence
                    -annotation     the AnnotationCollectionI object
                    -location       the LocationI object

=cut

sub set_attributes {
    my ($self,@args) = @_;
    my ($start, $end, $strand, $primary_tag, $source_tag, $primary, 
        $source, $frame, $score, $tag, $gff_string, $gff1_string,
        $seqname, $seqid, $annot, $location,$display_name, $pid,$phase) =
            $self->_rearrange([qw(START
                                  END
                                  STRAND
                                  PRIMARY_TAG
                                  SOURCE_TAG
                                  PRIMARY
                                  SOURCE
                                  FRAME
                                  SCORE
                                  TAG
                                  GFF_STRING
                                  GFF1_STRING
                                  SEQNAME
                                  SEQ_ID
                                  ANNOTATION
                                  LOCATION
                                  DISPLAY_NAME
                                  PRIMARY_ID
                                  PHASE
                                  )], @args);
    $location    && $self->location($location);
    $gff_string  && $self->_from_gff_string($gff_string);
    $gff1_string  && do {
        $self->gff_format(Bio::Tools::GFF->new('-gff_version' => 1));
        $self->_from_gff_stream($gff1_string);
    };
    
    $pid                    && $self->primary_id($pid);
    $primary_tag            && $self->primary_tag($primary_tag);
    $source_tag             && $self->source_tag($source_tag);
    $primary                && $self->primary_tag($primary);
    $source                 && $self->source_tag($source);
    defined $start          && $self->start($start);
    defined $end            && $self->end($end);
    defined $strand         && $self->strand($strand);
    defined $frame          && $self->frame($frame);
    defined $display_name   && $self->display_name($display_name);
    defined $score          && $self->score($score);
    $annot                  && $self->annotation($annot);
    if($seqname) {
        $self->warn("-seqname is deprecated. Please use -seq_id instead.");
        $seqid = $seqname unless $seqid;
    }
    $self->seq_id($seqid) if (defined($seqid));
    $tag            && do {
        foreach my $t ( keys %$tag ) {
            $self->add_tag_value($t, UNIVERSAL::isa($tag->{$t}, "ARRAY") ? @{$tag->{$t}} : $tag->{$t});
        }
    };
    defined $phase          && $self->phase($phase);
}


=head2 direct_new

 Title   : direct_new
 Usage   : my $feat = Bio::SeqFeature::Generic->direct_new;
 Function: create a blessed hash - for performance improvement in 
           object creation
 Returns : Bio::SeqFeature::Generic object
 Args    : none

=cut

sub direct_new {
    my ( $class) = @_;
    my ($self) = {};

    bless $self,$class;

    return $self;
}


=head2 location

 Title   : location
 Usage   : my $location = $feat->location();
 Function: returns a location object suitable for identifying location 
           of feature on sequence or parent feature  
 Returns : Bio::LocationI object
 Args    : [optional] Bio::LocationI object to set the value to.

=cut

sub location {
    my($self, $value ) = @_;  

    if (defined($value)) {
        unless (ref($value) and $value->isa('Bio::LocationI')) {
            $self->throw("object $value pretends to be a location but ".
                         "does not implement Bio::LocationI");
        }
        $self->{'_location'} = $value;
    }
    elsif (! $self->{'_location'}) {
        # guarantees a real location object is returned every time
        $self->{'_location'} = Bio::Location::Simple->new();
    }
    return $self->{'_location'};
}


=head2 start

 Title   : start
 Usage   : my $start = $feat->start;
           $feat->start(20);
 Function: Get/set on the start coordinate of the feature
 Returns : integer
 Args    : none

=cut

sub start {
    my ($self, $value) = @_;
    # Return soon if setting value
    if (defined $value) {
        return $self->location->start($value);
    }

    return $self->location->start() if not defined $self->{'_gsf_seq'};
    # Check circular sequences cut by origin
    my $start;
    if (    $self->{'_gsf_seq'}->is_circular
        and $self->location->isa('Bio::Location::SplitLocationI')
        ) {
        my $primary_seq_length = $self->{'_gsf_seq'}->length;
        my @sublocs = $self->location->sub_Location;

        my $cut_by_origin = 0;
        my ($a_end,   $a_strand) = (0, 0);
        my ($b_start, $b_strand) = (0, 0);
        for (my $i = 1; $i < scalar @sublocs; $i++) {
            $a_end    = $sublocs[$i-1]->end;
            $a_strand = $sublocs[$i-1]->strand;
            $b_start  = $sublocs[$i]->start;
            $b_strand = $sublocs[$i]->strand;
            # cut by origin condition
            if (    $a_end    == $primary_seq_length
                and $b_start  == 1
                and $a_strand == $b_strand
                ) {
                $cut_by_origin = 1;
                last;
            }
        }
        $start = ($cut_by_origin == 1) ? ($sublocs[0]->start) : ($self->location->start);
    }
    else {
        $start = $self->location->start;
    }
    return $start;
}


=head2 end

 Title   : end
 Usage   : my $end = $feat->end;
           $feat->end($end);
 Function: get/set on the end coordinate of the feature
 Returns : integer
 Args    : none

=cut

sub end {
    my ($self, $value) = @_;
    # Return soon if setting value
    if (defined $value) {
        return $self->location->end($value);
    }

    return $self->location->end() if not defined $self->{'_gsf_seq'};
    # Check circular sequences cut by origin
    my $end;
    if (    $self->{'_gsf_seq'}->is_circular
        and $self->location->isa('Bio::Location::SplitLocationI')
        ) {
        my $primary_seq_length = $self->{'_gsf_seq'}->length;
        my @sublocs = $self->location->sub_Location;

        my $cut_by_origin = 0;
        my ($a_end,   $a_strand) = (0, 0);
        my ($b_start, $b_strand) = (0, 0);
        for (my $i = 1; $i < scalar @sublocs; $i++) {
            $a_end    = $sublocs[$i-1]->end;
            $a_strand = $sublocs[$i-1]->strand;
            $b_start  = $sublocs[$i]->start;
            $b_strand = $sublocs[$i]->strand;
            # cut by origin condition
            if (    $a_end    == $primary_seq_length
                and $b_start  == 1
                and $a_strand == $b_strand
                ) {
                $cut_by_origin = 1;
                last;
            }
        }
        $end = ($cut_by_origin == 1) ? ($sublocs[-1]->end) : ($self->location->end);
    }
    else {
        $end = $self->location->end;
    }
    return $end;
}


=head2 length

 Title   : length
 Usage   : my $len = $feat->length;
 Function: Get the feature length computed as:
              $feat->end - $feat->start + 1
 Returns : integer
 Args    : none

=cut

sub length {
    my $self   = shift;
    my $length = $self->end() - $self->start() + 1;

    # In circular sequences cut by origin $start > $end,
    # e.g., join(5075..5386,1..51)), $start = 5075, $end = 51,
    # then adjust using the primary_seq length (5386)
    if ($length < 0 and defined $self->{'_gsf_seq'}) {
        $length += $self->{'_gsf_seq'}->length;
    }
    return $length;
}


=head2 strand

 Title   : strand
 Usage   : my $strand = $feat->strand();
           $feat->strand($strand);
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none

=cut

sub strand {
    my $self = shift;
    return $self->location->strand(@_);
}


=head2 score

 Title   : score
 Usage   : my $score = $feat->score();
           $feat->score($score);
 Function: get/set on score information
 Returns : float
 Args    : none if get, the new value if set

=cut

sub score {
    my $self = shift;

    if (@_) {
        my $value = shift;

        if ( defined $value && $value && $value !~ /^[A-Za-z]+$/ &&
            $value !~ /^[+-]?\d+\.?\d*(e-\d+)?/ and $value != 0) {
            $self->throw(-class=>'Bio::Root::BadParameter',
                         -text=>"'$value' is not a valid score",
                         -value=>$value);
        }
        if ($self->has_tag('score')) {
            $self->warn("Removing score value(s)");
            $self->remove_tag('score');
        }
        $self->add_tag_value('score',$value);
    }
    my ($score) = $self->has_tag('score') ? $self->get_tag_values('score') : undef;
    return $score;
}


=head2 frame

 Title   : frame
 Usage   : my $frame = $feat->frame();
           $feat->frame($frame);
 Function: get/set on frame information
 Returns : 0,1,2, '.'
 Args    : none if get, the new value if set

=cut

sub frame {
    my $self = shift;

    if ( @_ ) {
        my $value = shift;
        if ( defined $value && 
            $value !~ /^[0-2.]$/ ) {
            $self->throw("'$value' is not a valid frame");
        }
        if( defined $value && $value eq '.' ) { $value = '.' } 
        return $self->{'_gsf_frame'} = $value;
    }
    return $self->{'_gsf_frame'};
}


=head2 primary_tag

 Title   : primary_tag
 Usage   : my $tag = $feat->primary_tag();
           $feat->primary_tag('exon');
 Function: get/set on the primary tag for a feature,
           eg 'exon'
 Returns : a string
 Args    : none

=cut

sub primary_tag {
    my $self = shift;
    return $self->{'_primary_tag'} = shift if @_;
    return $self->{'_primary_tag'} || '';
}


=head2 source_tag

 Title   : source_tag
 Usage   : my $tag = $feat->source_tag();
           $feat->source_tag('genscan');
 Function: Returns the source tag for a feature,
           eg, 'genscan'
 Returns : a string
 Args    : none

=cut

sub source_tag {
    my $self = shift;
    return $self->{'_source_tag'} = shift if @_;
    return $self->{'_source_tag'} || '';
}


=head2 has_tag

 Title   : has_tag
 Usage   : my $value = $feat->has_tag('some_tag');
 Function: Tests wether a feature contaings a tag
 Returns : TRUE if the SeqFeature has the tag,
           and FALSE otherwise.
 Args    : The name of a tag

=cut

sub has_tag {
    my ($self, $tag) = @_;
    return exists $_[0]->{'_gsf_tag_hash'}->{$tag};
}


=head2 add_tag_value

 Title   : add_tag_value
 Usage   : $feat->add_tag_value('note',"this is a note");
 Returns : TRUE on success
 Args    : tag (string) and one or more values (any scalar(s))

=cut

sub add_tag_value {
    my $self = shift;
    my $tag = shift;
    $self->{'_gsf_tag_hash'}->{$tag} ||= [];
    push(@{$self->{'_gsf_tag_hash'}->{$tag}},@_);
}


=head2 get_tag_values

 Title   : get_tag_values
 Usage   : my @values = $feat->get_tag_values('note');
 Function: Returns a list of all the values stored
           under a particular tag.
 Returns : A list of scalars
 Args    : The name of the tag

=cut

sub get_tag_values {
    my ($self, $tag) = @_;

    if( ! defined $tag ) { return (); }
    if ( ! exists $self->{'_gsf_tag_hash'}->{$tag} ) {
        $self->throw("asking for tag value that does not exist $tag");
    }
    return @{$self->{'_gsf_tag_hash'}->{$tag}};
}


=head2 get_all_tags

 Title   : get_all_tags
 Usage   : my @tags = $feat->get_all_tags();
 Function: Get a list of all the tags in a feature
 Returns : An array of tag names
 Args    : none

# added a sort so that tags will be returned in a predictable order
# I still think we should be able to specify a sort function
# to the object at some point
# -js

=cut

sub get_all_tags {
    my ($self, @args) = @_;   
    return sort keys %{ $self->{'_gsf_tag_hash'}};
}


=head2 remove_tag

 Title   : remove_tag
 Usage   : $feat->remove_tag('some_tag');
 Function: removes a tag from this feature
 Returns : the array of values for this tag before removing it
 Args    : tag (string)

=cut

sub remove_tag {
    my ($self, $tag) = @_;

    if ( ! exists $self->{'_gsf_tag_hash'}->{$tag} ) {
        $self->throw("trying to remove a tag that does not exist: $tag");
    }
    my @vals = @{$self->{'_gsf_tag_hash'}->{$tag}};
    delete $self->{'_gsf_tag_hash'}->{$tag};
    return @vals;
}


=head2 attach_seq

 Title   : attach_seq
 Usage   : $feat->attach_seq($seq);
 Function: Attaches a Bio::Seq object to this feature. This
           Bio::Seq object is for the *entire* sequence: ie
           from 1 to 10000
 Example :
 Returns : TRUE on success
 Args    : a Bio::PrimarySeqI compliant object

=cut

sub attach_seq {
    my ($self, $seq) = @_;

    if ( ! ($seq && ref($seq) && $seq->isa("Bio::PrimarySeqI")) ) {
        $self->throw("Must attach Bio::PrimarySeqI objects to SeqFeatures but got '".ref($seq)."'");
    }

    $self->{'_gsf_seq'} = $seq;

    # attach to sub features if they want it
    foreach ( $self->sub_SeqFeature() ) {
        $_->attach_seq($seq);
    }
    return 1;
}


=head2 seq

 Title   : seq
 Usage   : my $tseq = $feat->seq();
 Function: returns the truncated sequence (if there) for this
 Example :
 Returns : sub seq (a Bio::PrimarySeqI compliant object) on attached sequence
           bounded by start & end, or undef if there is no sequence attached
 Args    : none

=cut

sub seq {
    my ($self, $arg) = @_;

    if ( defined $arg ) {
        $self->throw("Calling SeqFeature::Generic->seq with an argument. You probably want attach_seq");
    }

    if ( ! exists $self->{'_gsf_seq'} ) {
        return;
    }

    # assumming our seq object is sensible, it should not have to yank
    # the entire sequence out here.

    my $seq = $self->{'_gsf_seq'}->trunc($self->start(), $self->end());


    if ( defined $self->strand &&
        $self->strand == -1 ) {

        # ok. this does not work well (?)
        #print STDERR "Before revcom", $seq->str, "\n";
        $seq = $seq->revcom;
        #print STDERR "After  revcom", $seq->str, "\n";
    }

    return $seq;
}


=head2 entire_seq

 Title   : entire_seq
 Usage   : my $whole_seq = $feat->entire_seq();
 Function: gives the entire sequence that this seqfeature is attached to
 Example :
 Returns : a Bio::PrimarySeqI compliant object, or undef if there is no
           sequence attached
 Args    :

=cut

sub entire_seq {
    return shift->{'_gsf_seq'};
}


=head2 seq_id

 Title   : seq_id
 Usage   : $feat->seq_id($newval)
 Function: There are many cases when you make a feature that you
           do know the sequence name, but do not know its actual
           sequence. This is an attribute such that you can store
           the ID (e.g., display_id) of the sequence.

           This attribute should *not* be used in GFF dumping, as
           that should come from the collection in which the seq
           feature was found.
 Returns : value of seq_id
 Args    : newvalue (optional)

=cut

sub seq_id {
    my $obj = shift;
    return $obj->{'_gsf_seq_id'} = shift if @_;
    return $obj->{'_gsf_seq_id'};
}


=head2 display_name

 Title   : display_name
 Usage   : my $featname = $feat->display_name;
 Function: Implements the display_name() method, which is a human-readable
           name for the feature. 
 Returns : value of display_name (a string)
 Args    : Optionally, on set the new value or undef 

=cut

sub display_name {
    my $self = shift;
    return $self->{'display_name'} = shift if @_;
    return $self->{'display_name'} || '';
}


=head1 Methods for implementing Bio::AnnotatableI

=head2 annotation

 Title   : annotation
 Usage   : $feat->annotation($annot_obj);
 Function: Get/set the annotation collection object for annotating this
           feature.

 Example : 
 Returns : A Bio::AnnotationCollectionI object
 Args    : newvalue (optional)

=cut

sub annotation {
    my ($obj,$value) = @_;

    # we are smart if someone references the object and there hasn't been
    # one set yet
    if(defined $value || ! defined $obj->{'annotation'} ) {
        $value = Bio::Annotation::Collection->new() unless ( defined $value );
        $obj->{'annotation'} = $value;
    }
    return $obj->{'annotation'};
}


=head1 Methods to implement Bio::FeatureHolderI

This includes methods for retrieving, adding, and removing
features. Since this is already a feature, features held by this
feature holder are essentially sub-features.

=head2 get_SeqFeatures

 Title   : get_SeqFeatures
 Usage   : my @feats = $feat->get_SeqFeatures();
 Function: Returns an array of sub Sequence Features
 Returns : An array
 Args    : none

=cut

sub get_SeqFeatures {
    return @{ shift->{'_gsf_sub_array'} || []};    
}


=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   : $feat->add_SeqFeature($subfeat);
           $feat->add_SeqFeature($subfeat,'EXPAND');
 Function: Adds a SeqFeature into the subSeqFeature array.
           With no 'EXPAND' qualifer, subfeat will be tested
           as to whether it lies inside the parent, and throw
           an exception if not.

           If EXPAND is used, the parent's start/end/strand will
           be adjusted so that it grows to accommodate the new
           subFeature

           !IMPORTANT! The coordinates of the subfeature should not be relative
           to the parent feature it is attached to, but relative to the sequence
           the parent feature is located on.

 Returns : nothing
 Args    : An object which has the SeqFeatureI interface

=cut

sub add_SeqFeature {
    my ($self,$feat,$expand) = @_;
    unless( defined $feat ) {
        $self->warn("Called add_SeqFeature with no feature, ignoring");
        return;
    }
    if ( !$feat->isa('Bio::SeqFeatureI') ) {
        $self->warn("$feat does not implement Bio::SeqFeatureI. Will add it anyway, but beware...");
    }

    if($expand && ($expand eq 'EXPAND')) {
        $self->_expand_region($feat);
    } else {
        if ( !$self->contains($feat) ) {
            $self->throw("$feat is not contained within parent feature, and expansion is not valid");
        }
    }

    $self->{'_gsf_sub_array'} = [] unless exists($self->{'_gsf_sub_array'});
    push(@{$self->{'_gsf_sub_array'}},$feat);

}


=head2 remove_SeqFeatures

 Title   : remove_SeqFeatures
 Usage   : $feat->remove_SeqFeatures;
 Function: Removes all SeqFeatures

           If you want to remove only a subset of features then remove that 
           subset from the returned array, and add back the rest.
 Example :
 Returns : The array of Bio::SeqFeatureI implementing features that was
           deleted.
 Args    : none

=cut

sub remove_SeqFeatures {
    my ($self) = @_;
    my @subfeats = @{$self->{'_gsf_sub_array'} || []};
    $self->{'_gsf_sub_array'} = []; # zap the array implicitly.
    return @subfeats;
}


=head1 GFF-related methods

=head2 gff_format

 Title   : gff_format
 Usage   : # get:
           my $gffio = $feat->gff_format();
           # set (change the default version of GFF2):
           $feat->gff_format(Bio::Tools::GFF->new(-gff_version => 1));
 Function: Get/set the GFF format interpreter. This object is supposed to 
           format and parse GFF. See Bio::Tools::GFF for the interface.

           If this method is called as class method, the default for all
           newly created instances will be changed. Otherwise only this
           instance will be affected.
 Example : 
 Returns : a Bio::Tools::GFF compliant object
 Args    : On set, an instance of Bio::Tools::GFF or a derived object.

=cut

sub gff_format {
    my ($self, $gffio) = @_;
    if(defined($gffio)) {
        if(ref($self)) {
            $self->{'_gffio'} = $gffio;
        } else {
            $Bio::SeqFeatureI::static_gff_formatter = $gffio;
        }
    }
    return (ref($self) && exists($self->{'_gffio'}) ?
            $self->{'_gffio'} : $self->_static_gff_formatter);
}


=head2 gff_string

 Title   : gff_string
 Usage   : my $str = $feat->gff_string;
           my $str = $feat->gff_string($gff_formatter);
 Function: Provides the feature information in GFF format.

           We override this here from Bio::SeqFeatureI in order to use the
           formatter returned by gff_format().

 Returns : A string
 Args    : Optionally, an object implementing gff_string().

=cut

sub gff_string {
    my ($self,$formatter) = @_;
    $formatter = $self->gff_format() unless $formatter;
    return $formatter->gff_string($self);
}


=head2 slurp_gff_file

 Title   : slurp_file
 Usage   : my @features = Bio::SeqFeature::Generic::slurp_gff_file(\*FILE);
 Function: Sneaky function to load an entire file as in memory objects.
           Beware of big files.

           This method is deprecated. Use Bio::Tools::GFF instead, which can
           also handle large files.

 Example :
 Returns :
 Args    :

=cut

sub slurp_gff_file {
    my ($f) = @_;
    my @out;
    if ( !defined $f ) {
        Bio::Root::Root->throw("Must have a filehandle");
    }

    Bio::Root::Root->deprecated( -message => "deprecated method slurp_gff_file() called in Bio::SeqFeature::Generic. Use Bio::Tools::GFF instead.",
                                 -warn_version  => '1.005',
                                 -throw_version => '1.007',
                               );

    while(<$f>) {
        my $sf = Bio::SeqFeature::Generic->new('-gff_string' => $_);
        push(@out, $sf);
    }

    return @out;
}


=head2 _from_gff_string

 Title   : _from_gff_string
 Usage   :
 Function: Set feature properties from GFF string. 

           This method uses the object returned by gff_format() for the
           actual interpretation of the string. Set a different GFF format
           interpreter first if you need a specific version, like GFF1. (The
           default is GFF2.)
 Example :
 Returns : 
 Args    : a GFF-formatted string

=cut

sub _from_gff_string {
    my ($self, $string) = @_;
    $self->gff_format()->from_gff_string($self, $string);
}


=head2 _expand_region

 Title   : _expand_region
 Usage   : $feat->_expand_region($feature);
 Function: Expand the total region covered by this feature to
           accommodate for the given feature.

           May be called whenever any kind of subfeature is added to this
           feature. add_SeqFeature() already does this.
 Returns : 
 Args    : A Bio::SeqFeatureI implementing object.

=cut

sub _expand_region {
    my ($self, $feat) = @_;
    if(! $feat->isa('Bio::SeqFeatureI')) {
        $self->warn("$feat does not implement Bio::SeqFeatureI");
    }
    # if this doesn't have start set - forget it!
    # changed to reflect sanity checks for LocationI
    if(!$self->location->valid_Location) {
        $self->start($feat->start);
        $self->end($feat->end);
        $self->strand($feat->strand) unless $self->strand;
    } else {
        my ($start,$end,$strand) = $self->union($feat);
        $self->start($start);
        $self->end($end);
        $self->strand($strand);
    }
}


=head2 _parse

 Title   : _parse
 Usage   :
 Function: Parsing hints
 Example :
 Returns :
 Args    :

=cut

sub _parse {
    my ($self) = @_;
    return $self->{'_parse_h'};
}


=head2 _tag_value

 Title   : _tag_value
 Usage   : 
 Function: For internal use only. Convenience method for those tags that
           may only have a single value.
 Returns : The first value under the given tag as a scalar (string)
 Args    : The tag as a string. Optionally, the value on set.

=cut

sub _tag_value {
    my $self = shift;
    my $tag = shift;

    if(@_ || (! $self->has_tag($tag))) {
        $self->remove_tag($tag) if($self->has_tag($tag));
        $self->add_tag_value($tag, @_);
    }
    return ($self->get_tag_values($tag))[0];
}


#######################################################################
# aliases for methods that changed their names in an attempt to make  #
# bioperl names more consistent                                       #
#######################################################################

sub seqname {
    my $self = shift;
    $self->warn("SeqFeatureI::seqname() is deprecated. Please use seq_id() instead.");
    return $self->seq_id(@_);
}

sub display_id {
    my $self = shift;
    $self->warn("SeqFeatureI::display_id() is deprecated. Please use display_name() instead.");
    return $self->display_name(@_);
}

# this is towards consistent naming
sub each_tag_value { return shift->get_tag_values(@_); }
sub all_tags { return shift->get_all_tags(@_); }

# we revamped the feature containing property to implementing
# Bio::FeatureHolderI
*sub_SeqFeature = \&get_SeqFeatures;
*add_sub_SeqFeature = \&add_SeqFeature;
*flush_sub_SeqFeatures = \&remove_SeqFeatures;
# this one is because of inconsistent naming ...
*flush_sub_SeqFeature = \&remove_SeqFeatures;

sub cleanup_generic {
    my $self = shift;
    foreach my $f ( @{$self->{'_gsf_sub_array'} || []} ) {
        $f = undef;
    }
    $self->{'_gsf_seq'} = undef;
    foreach my $t ( keys %{$self->{'_gsf_tag_hash'} } ) {
        $self->{'_gsf_tag_hash'}->{$t} = undef;
        delete($self->{'_gsf_tag_hash'}->{$t}); # bug 1720 fix
    }
}

1;
