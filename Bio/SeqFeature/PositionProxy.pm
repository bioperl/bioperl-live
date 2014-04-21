#
# BioPerl module for Bio::SeqFeature::PositionProxy
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::PositionProxy - handle features when truncation/revcom sequences span a feature

=head1 SYNOPSIS

   $proxy = Bio::SeqFeature::PositionProxy->new( -loc => $loc,
                                                 -parent => $basefeature);

   $seq->add_SeqFeature($feat);

=head1 DESCRIPTION

PositionProxy is a Proxy Sequence Feature to handle truncation
and revcomp without duplicating all the data within the sequence features.
It holds a new location for a sequence feature and the original feature
it came from to provide the additional annotation information.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Ewan Birney

Ewan Birney E<lt>birney@sanger.ac.ukE<gt>

=head1 DEVELOPERS

This class has been written with an eye out of inheritence. The fields
the actual object hash are:

   _gsf_tag_hash  = reference to a hash for the tags
   _gsf_sub_array = reference to an array for sub arrays
   _gsf_start     = scalar of the start point
   _gsf_end       = scalar of the end point
   _gsf_strand    = scalar of the strand

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::PositionProxy;
use strict;

use Bio::Tools::GFF;


use base qw(Bio::Root::Root Bio::SeqFeatureI);

sub new {
    my ($caller, @args) = @_;
    my $self = $caller->SUPER::new(@args);

    my ($feature,$location) = $self->_rearrange([qw(PARENT LOC)],@args);

    if( !defined $feature || !ref $feature || !$feature->isa('Bio::SeqFeatureI') ) {
      $self->throw("Must have a parent feature, not a [$feature]");
    }

    if( $feature->isa("Bio::SeqFeature::PositionProxy") ) {
      $feature = $feature->parent();
    }

    if( !defined $location || !ref $location || !$location->isa('Bio::LocationI') ) {
      $self->throw("Must have a location, not a [$location]");
    }


    return $self;
}


=head2 location

 Title   : location
 Usage   : my $location = $seqfeature->location()
 Function: returns a location object suitable for identifying location 
           of feature on sequence or parent feature  
 Returns : Bio::LocationI object
 Args    : none

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


=head2 parent

 Title   : parent
 Usage   : my $sf = $proxy->parent()
 Function: returns the seqfeature parent of this proxy
 Returns : Bio::SeqFeatureI object
 Args    : none

=cut

sub parent {
    my($self, $value ) = @_;  

    if (defined($value)) {
        unless (ref($value) and $value->isa('Bio::SeqFeatureI')) {
            $self->throw("object $value pretends to be a location but ".
                "does not implement Bio::SeqFeatureI");
        }
        $self->{'_parent'} = $value;
    }

    return $self->{'_parent'};
}



=head2 start

 Title   : start
 Usage   : $start = $feat->start
           $feat->start(20)
 Function: Get
 Returns : integer
 Args    : none

=cut

sub start {
   my ($self,$value) = @_;
   return $self->location->start($value);
}


=head2 end

 Title   : end
 Usage   : $end = $feat->end
           $feat->end($end)
 Function: get
 Returns : integer
 Args    : none

=cut

sub end {
   my ($self,$value) = @_;
   return $self->location->end($value);
}


=head2 length

 Title   : length
 Usage   :
 Function:
 Example :
 Returns :
 Args    :

=cut

sub length {
   my ($self) = @_;
   return $self->end - $self->start() + 1;
}


=head2 strand

 Title   : strand
 Usage   : $strand = $feat->strand()
           $feat->strand($strand)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none

=cut

sub strand {
   my ($self,$value) = @_;
   return $self->location->strand($value);
}


=head2 attach_seq

 Title   : attach_seq
 Usage   : $sf->attach_seq($seq)
 Function: Attaches a Bio::Seq object to this feature. This
           Bio::Seq object is for the *entire* sequence: ie
           from 1 to 10000
 Example :
 Returns : TRUE on success
 Args    :

=cut

sub attach_seq {
   my ($self, $seq) = @_;

   if ( !defined $seq || !ref $seq || ! $seq->isa("Bio::PrimarySeqI") ) {
       $self->throw("Must attach Bio::PrimarySeqI objects to SeqFeatures");
   }

   $self->{'_gsf_seq'} = $seq;

   # attach to sub features if they want it

   foreach my $sf ( $self->sub_SeqFeature() ) {
       if ( $sf->can("attach_seq") ) {
           $sf->attach_seq($seq);
       }
   }
   return 1;
}


=head2 seq

 Title   : seq
 Usage   : $tseq = $sf->seq()
 Function: returns the truncated sequence (if there) for this
 Example :
 Returns : sub seq on attached sequence bounded by start & end
 Args    : none

=cut

sub seq {
   my ($self, $arg) = @_;

   if ( defined $arg ) {
       $self->throw("Calling SeqFeature::PositionProxy->seq with an argument. You probably want attach_seq");
   }

   if ( ! exists $self->{'_gsf_seq'} ) {
       return;
   }

   # assumming our seq object is sensible, it should not have to yank
   # the entire sequence out here.

   my $seq = $self->{'_gsf_seq'}->trunc($self->start(), $self->end());


   if ( $self->strand == -1 ) {
       $seq = $seq->revcom;
   }

   return $seq;
}


=head2 entire_seq

 Title   : entire_seq
 Usage   : $whole_seq = $sf->entire_seq()
 Function: gives the entire sequence that this seqfeature is attached to
 Example :
 Returns :
 Args    :

=cut

sub entire_seq {
   my ($self) = @_;

   return unless exists($self->{'_gsf_seq'});
   return $self->{'_gsf_seq'};
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

sub seqname {
    my ($obj,$value) = @_;
    if ( defined $value ) {
        $obj->{'_gsf_seqname'} = $value;
    }
    return $obj->{'_gsf_seqname'};
}


=head2 Proxies

These functions chain back to the parent for all non sequence related stuff.

=cut

=head2 primary_tag

 Title   : primary_tag
 Usage   : $tag = $feat->primary_tag()
 Function: Returns the primary tag for a feature,
           eg 'exon'
 Returns : a string 
 Args    : none

=cut

sub primary_tag {
   my ($self,@args) = @_;

   return $self->parent->primary_tag();
}


=head2 source_tag

 Title   : source_tag
 Usage   : $tag = $feat->source_tag()
 Function: Returns the source tag for a feature,
           eg, 'genscan' 
 Returns : a string 
 Args    : none

=cut

sub source_tag {
   my ($self) = @_;

   return $self->parent->source_tag();
}


=head2 has_tag

 Title   : has_tag
 Usage   : $tag_exists = $self->has_tag('some_tag')
 Function: 
 Returns : TRUE if the specified tag exists, and FALSE otherwise
 Args    :

=cut

sub has_tag {
   my ($self,$tag) = @_;

   return $self->parent->has_tag($tag);
}


=head2 get_tag_values

 Title   : get_tag_values
 Usage   : @values = $self->get_tag_values('some_tag')
 Function: 
 Returns : An array comprising the values of the specified tag.
 Args    :

=cut

*each_tag_value = \&get_tag_values;

sub get_tag_values {
   my ($self,$tag) = @_;

   return $self->parent->get_tag_values($tag);
}


=head2 get_all_tags

 Title   : get_all_tags
 Usage   : @tags = $feat->get_all_tags()
 Function: gives all tags for this feature
 Returns : an array of strings
 Args    : none

=cut

*all_tags = \&get_all_tags;

sub get_all_tags {
   my ($self) = @_;

   return $self->parent->all_tags();
}

1;
