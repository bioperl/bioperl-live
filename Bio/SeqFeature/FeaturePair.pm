#
# BioPerl module for Bio::SeqFeature::FeaturePair
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

Bio::SeqFeature::FeaturePair - hold pair feature information e.g. blast hits

=head1 SYNOPSIS

    my $feat  = Bio::SeqFeature::FeaturePair->new(
        -feature1 => $f1,
        -feature2 => $f2,
    );

    # Bio::SeqFeatureI methods can be used

    my $start = $feat->start;
    my $end   = $feat->end;

    # Bio::FeaturePair methods can be used
    my $hstart = $feat->hstart;
    my $hend   = $feat->hend;

   my $feature1 = $feat->feature1;  # returns feature1 object

=head1 DESCRIPTION

A sequence feature object where the feature is itself a feature on
another sequence - e.g. a blast hit where residues 1-40 of a protein
sequence SW:HBA_HUMAN has hit to bases 100 - 220 on a genomic sequence
HS120G22.  The genomic sequence coordinates are used to create one
sequence feature $f1 and the protein coordinates are used to create
feature $f2.  A FeaturePair object can then be made

    my $fp = Bio::SeqFeature::FeaturePair->new(
        -feature1 => $f1,   # genomic
        -feature2 => $f2,   # protein
    );

This object can be used as a standard Bio::SeqFeatureI in which case

    my $gstart = $fp->start  # returns start coord on feature1 - genomic seq.
    my $gend   = $fp->end    # returns end coord on feature1.

In general standard Bio::SeqFeatureI method calls return information
in feature1.

Data in the feature 2 object are generally obtained using the standard
methods prefixed by h (for hit!)

    my $pstart = $fp->hstart # returns start coord on feature2 = protein seq.
    my $pend   = $fp->hend   # returns end coord on feature2.

If you wish to swap feature1 and feature2 around :

    $feat->invert

so... 

    $feat->start # etc. returns data in $feature2 object


No sub_SeqFeatures or tags can be stored in this object directly.  Any
features or tags are expected to be stored in the contained objects
feature1, and feature2.

=head1 CONTACT

Ewan Birney E<lt>birney@sanger.ac.ukE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::FeaturePair;
use vars qw($AUTOLOAD);
use strict;

use Bio::SeqFeatureI;
use Bio::Factory::ObjectFactory;

use base qw(Bio::SeqFeature::Generic);

=head2 new

 Title   : new
 Usage   :
 Function: Constructor for this module. Accepts the following parameters:

             -feature1   Bio::SeqFeatureI-compliant object
             -feature2   Bio::SeqFeatureI-compliant object
             -feature_factory  Bio::Factory::ObjectFactoryI compliant
                         object to be used when feature1 and/or feature2
                         are accessed without explicitly set before. This
                         is mostly useful for derived classes who want to
                         set their preferred class for feature objects.

 Example :
 Returns : 
 Args    : see above


=cut

sub new {
    my ($class, @args) = @_;

    #
    # We've got a certain problem here that somewhat relates to chicken and
    # eggs. The problem is, we override a lot of SeqFeatureI methods here
    # to delegate them to either feature1 or feature2. If we pass along
    # those attributes right away, we need feature1 or feature2 or the feature
    # factory in place, or there is no way around the dreaded default, which
    # is ugly too (as it necessitates subsequent copying if you wanted a
    # different feature object class).
    #
    # So I decided to go with the lesser of two evils here: we need to assume
    # here that we can set all attributes through set_attributes(), which we
    # assume is no different from setting them through the constructor. This
    # gives us a window to set the feature objects and the factory, such that
    # any derived class doesn't have to worry about this any more.
    #
    # I'm happy to hear a better solution, but I think this one isn't so bad.
    #
    my $self = $class->SUPER::new();
    my ($feature1,$feature2,$featfact) = 
        $self->_rearrange([qw( FEATURE1
                               FEATURE2
                               FEATURE_FACTORY )],@args);
    
    $self->_register_for_cleanup(\&cleanup_fp);
    # initialize the feature object factory if not provided
    if(! $featfact) {
        $featfact = Bio::Factory::ObjectFactory->new(
            -type => "Bio::SeqFeature::Generic",
            -interface => "Bio::SeqFeatureI"
        );
    }
    $self->feature_factory($featfact);
    # Store the features in the object
    $feature1 && $self->feature1($feature1);
    $feature2 && $self->feature2($feature2);
    
    # OK. Now we're setup to store all the attributes, and they'll go right
    # away into the right objects.
    $self->set_attributes(@args);

    # done - we hope
    return $self;
}

=head2 feature1

 Title   : feature1
 Usage   : $f = $featpair->feature1
           $featpair->feature1($feature)
 Function: Get/set for the query feature
 Returns : Bio::SeqFeatureI
 Args    : Bio::SeqFeatureI


=cut

sub feature1 {
    my ($self,$arg) = @_;    
    if ( defined($arg) || !defined $self->{'feature1'} ) {
        $self->throw("internal error: feature factory not set!") 
            unless $self->feature_factory;
        $arg = $self->feature_factory->create_object() unless( defined $arg);
        $self->throw("Argument [$arg] must be a Bio::SeqFeatureI") 
            unless (ref($arg) && $arg->isa("Bio::SeqFeatureI"));
        $self->{'feature1'} = $arg;
    }
    return $self->{'feature1'};
}

=head2 feature2

 Title   : feature2
 Usage   : $f = $featpair->feature2
           $featpair->feature2($feature)
 Function: Get/set for the hit feature
 Returns : Bio::SeqFeatureI
 Args    : Bio::SeqFeatureI


=cut

sub feature2 {
    my ($self,$arg) = @_;

    if ( defined($arg) || ! defined $self->{'feature2'}) {
        $self->throw("internal error: feature factory not set!") 
            unless $self->feature_factory;
        $arg = $self->feature_factory->create_object() unless( defined $arg);
        $self->throw("Argument [$arg] must be a Bio::SeqFeatureI") 
            unless (ref($arg) && $arg->isa("Bio::SeqFeatureI"));
        $self->{'feature2'} = $arg;
    }
    return $self->{'feature2'};
}

=head2 start

 Title   : start
 Usage   : $start = $featpair->start
           $featpair->start(20)
 Function: Get/set on the start coordinate of feature1
 Returns : integer
 Args    : [optional] beginning of feature

=cut

sub start {
    return shift->feature1->start(@_);
}

=head2 end

 Title   : end
 Usage   : $end = $featpair->end
           $featpair->end($end)
 Function: get/set on the end coordinate of feature1
 Returns : integer
 Args    : [optional] ending point of feature


=cut

sub end {
    return shift->feature1->end(@_);    
}

=head2 strand

 Title   : strand
 Usage   : $strand = $feat->strand()
           $feat->strand($strand)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : [optional] strand information to set


=cut

sub strand {
    return shift->feature1->strand(@_);    
}

=head2 location

 Title   : location
 Usage   : $location = $featpair->location
           $featpair->location($location)
 Function: Get/set location object (using feature1)
 Returns : Bio::LocationI object
 Args    : [optional] LocationI to store

=cut

sub location {
    return shift->feature1->location(@_);
}

=head2 score

 Title   : score
 Usage   : $score = $feat->score()
           $feat->score($score)
 Function: get/set on score information
 Returns : float
 Args    : none if get, the new value if set


=cut

sub score {
    return shift->feature1->score(@_);    
}

=head2 frame

 Title   : frame
 Usage   : $frame = $feat->frame()
           $feat->frame($frame)
 Function: get/set on frame information
 Returns : 0,1,2
 Args    : none if get, the new value if set


=cut

sub frame {
    return shift->feature1->frame(@_);    
}

=head2 primary_tag

 Title   : primary_tag
 Usage   : $ptag = $featpair->primary_tag
 Function: get/set on the primary_tag of feature1
 Returns : 0,1,2
 Args    : none if get, the new value if set


=cut

sub primary_tag {
    return shift->feature1->primary_tag(@_);    
}

=head2 source_tag

 Title   : source_tag
 Usage   : $tag = $feat->source_tag()
           $feat->source_tag('genscan');
 Function: Returns the source tag for a feature,
           eg, 'genscan' 
 Returns : a string 
 Args    : none


=cut

sub source_tag {
    return shift->feature1->source_tag(@_);    
}

=head2 seqname

 Title   : seqname
 Usage   : $obj->seq_id($newval)
 Function: There are many cases when you make a feature that you
           do know the sequence name, but do not know its actual
           sequence. This is an attribute such that you can store 
           the seqname.

           This attribute should *not* be used in GFF dumping, as
           that should come from the collection in which the seq
           feature was found.
 Returns : value of seqname
 Args    : newvalue (optional)


=cut

sub seq_id {
    return shift->feature1->seq_id(@_);    
}

=head2 hseqname

 Title   : hseqname
 Usage   : $featpair->hseqname($newval)
 Function: Get/set method for the name of
           feature2.
 Returns : value of $feature2->seq_id
 Args    : newvalue (optional)


=cut

sub hseq_id {
    return shift->feature2->seq_id(@_);
}


=head2 hstart

 Title   : hstart
 Usage   : $start = $featpair->hstart
           $featpair->hstart(20)
 Function: Get/set on the start coordinate of feature2
 Returns : integer
 Args    : none

=cut

sub hstart {
    return shift->feature2->start(@_);    
}

=head2 hend

 Title   : hend
 Usage   : $end = $featpair->hend
           $featpair->hend($end)
 Function: get/set on the end coordinate of feature2
 Returns : integer
 Args    : none


=cut

sub hend {
    return shift->feature2->end(@_);    
}


=head2 hstrand

 Title   : hstrand
 Usage   : $strand = $feat->strand()
           $feat->strand($strand)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none


=cut

sub hstrand {
    return shift->feature2->strand(@_);
}

=head2 hscore

 Title   : hscore
 Usage   : $score = $feat->score()
           $feat->score($score)
 Function: get/set on score information
 Returns : float
 Args    : none if get, the new value if set


=cut

sub hscore {
    return shift->feature2->score(@_);    
}

=head2 hframe

 Title   : hframe
 Usage   : $frame = $feat->frame()
           $feat->frame($frame)
 Function: get/set on frame information
 Returns : 0,1,2
 Args    : none if get, the new value if set


=cut

sub hframe {
    return shift->feature2->frame(@_);    
}

=head2 hprimary_tag

 Title   : hprimary_tag
 Usage   : $ptag = $featpair->hprimary_tag
 Function: Get/set on the primary_tag of feature2
 Returns : 0,1,2
 Args    : none if get, the new value if set


=cut

sub hprimary_tag {
    return shift->feature2->primary_tag(@_);    
}

=head2 hsource_tag

 Title   : hsource_tag
 Usage   : $tag = $feat->hsource_tag()
           $feat->source_tag('genscan');
 Function: Returns the source tag for a feature,
           eg, 'genscan' 
 Returns : a string 
 Args    : none


=cut

sub hsource_tag {
    return shift->feature2->source_tag(@_);
}

=head2 invert

 Title   : invert
 Usage   : $tag = $feat->invert
 Function: Swaps feature1 and feature2 around
 Returns : Nothing
 Args    : none


=cut

sub invert {
    my ($self) = @_;

    my $tmp = $self->feature1;
    
    $self->feature1($self->feature2);
    $self->feature2($tmp);
    return 1;
}

=head2 feature_factory

 Title   : feature_factory
 Usage   : $obj->feature_factory($newval)
 Function: Get/set the feature object factory for this feature pair.

           The feature object factory will be used to create a feature
           object if feature1() or feature2() is called in get mode
           without having been set before.

           The default is an instance of Bio::Factory::ObjectFactory
           and hence allows the type to be changed dynamically at any
           time.

 Example : 
 Returns : The feature object factory in use (a 
           Bio::Factory::ObjectFactoryI compliant object)
 Args    : on set, a Bio::Factory::ObjectFactoryI compliant object


=cut

sub feature_factory {
    my $self = shift;

    return $self->{'feature_factory'} = shift if @_;
    return $self->{'feature_factory'};
}

#################################################################
# aliases for backwards compatibility                           #
#################################################################

# seqname() is already aliased in Generic.pm, and we overwrite seq_id

sub hseqname {
    my $self = shift;
    $self->warn("SeqFeatureI::seqname() is deprecated. Please use seq_id() instead.");
    return $self->hseq_id(@_);
}

sub cleanup_fp {
    my $self = shift;
    $self->{'feature1'} = $self->{'feature2'} = undef;
}
1;
