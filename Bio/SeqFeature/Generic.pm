# $Id$
#
# BioPerl module for Bio::SeqFeature::Generic
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

   $feat = new Bio::SeqFeature::Generic ( -start => 10, -end => 100,
				-strand => -1, -primary => 'repeat',
				-source => 'repeatmasker',
				-score  => 1000,
				-tag    => {
				    new => 1,
				    author => 'someone',
				    sillytag => 'this is silly!' } );

   $feat = new Bio::SeqFeature::Generic ( -gff_string => $string );
   # if you want explicitly GFF1
   $feat = new Bio::SeqFeature::Generic ( -gff1_string => $string );

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
is what is done in the Bio::SeqFeature::Gene,
Bio::SeqFeature::Transcript and Bio::SeqFeature::Exon, which provide
well coordinated classes to represent genes on DNA sequence (for
example, you can get the protein sequence out from a transcript
class).

For many Features, you want to add some piece of information, for
example a common one is that this feature is 'new' whereas other
features are 'old'.  The tag system, which here is implemented using a
hash can be used here.  You can use the tag system to extend the
SeqFeature::Generic programmatically: that is, you know that you have
read in more information into the tag 'mytag' which you can the
retrieve. This means you do not need to know how to write inherieted
Perl to provide more complex information on a feature, and/or, if you
do know but you donot want to write a new class every time you need
some extra piece of information, you can use the tag system to easily
store and then retrieve information.

The tag system can be written in/out of GFF format, and also into EMBL
format via AnnSeqIO::EMBL.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Ewan Birney <birney@sanger.ac.uk>

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


package Bio::SeqFeature::Generic;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
use Bio::SeqFeatureI;
use Bio::Annotation;
use Bio::Location::Simple;
use Bio::Tools::GFF;

@ISA = qw(Bio::Root::RootI Bio::SeqFeatureI);

sub new {
    my ( $caller, @args) = @_;   
    my ($self) = $caller->SUPER::new(@args); 

    $self->{'_gsf_tag_hash'} = {};
    $self->{'_gsf_sub_array'} = [];
    $self->{'_parse_h'} = {};
    my ($start, $end, $strand, $primary, $source, $frame, 
	$score, $tag, $gff_string, $gff1_string, $seqname, $annot, $location) =
	    $self->_rearrange([qw(START
				  END
				  STRAND
				  PRIMARY
				  SOURCE
				  FRAME
				  SCORE
				  TAG
				  GFF_STRING
				  GFF1_STRING
				  SEQNAME
				  ANNOTATION
				  LOCATION
				  )], @args);
    $location    && $self->location($location);
    $gff_string  && $self->_from_gff_string($gff_string);
    $gff1_string  && do {
	$self->gff_format(Bio::Tools::GFF->new('-gff_version' => 1));
	$self->_from_gff_stream($gff1_string);
    };
    $primary     && $self->primary_tag($primary);
    $source      && $self->source_tag($source);
    $start       && $self->start($start);
    $end         && $self->end($end);
    $strand      && $self->strand($strand);
    $frame       && $self->frame($frame);
    $score       && $self->score($score);
    $seqname     && $self->seqname($seqname);
    $annot       && $self->annotation($annot);
    $tag         && do {
	foreach my $t ( keys %$tag ) {
	    $self->add_tag_value($t,$tag->{$t});
	}
    };
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
   my ($self,$value) = @_;  

   # guarantees a real location object is returned every time
   if( defined($value) || !exists($self->{'_location'}) ) {
       if(defined($value)) {
	   if(! $value->isa('Bio::LocationI')) {
	       $self->throw("object $value pretends to be a location but ".
			    "does not implement Bio::LocationI");
	   }
       } else {
	   $value = Bio::Location::Simple->new();
       }
       $self->{'_location'} = $value;
   }
   return $self->{'_location'};
}


=head2 start

 Title   : start
 Usage   : $start = $feat->start
           $feat->start(20)
 Function: Get/set on the start coordinate of the feature
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
 Function: get/set on the end coordinate of the feature
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

=head2 score

 Title   : score
 Usage   : $score = $feat->score()
           $feat->score($score)
 Function: get/set on score information
 Returns : float
 Args    : none if get, the new value if set


=cut

sub score {
  my ($self,$value) = @_;

  if (defined($value)) {
       if ( $value !~ /^[+-]?\d+\.?\d*(e-\d+)?/ ) {
	   $self->throw("'$value' is not a valid score");
       }
       $self->{'_gsf_score'} = $value;
  }

  return $self->{'_gsf_score'};
}

=head2 frame

 Title   : frame
 Usage   : $frame = $feat->frame()
           $feat->frame($frame)
 Function: get/set on frame information
 Returns : 0,1,2, '.'
 Args    : none if get, the new value if set


=cut

sub frame {
  my ($self,$value) = @_;

  if ( defined $value ) {
       if ( $value !~ /^[0-2.]$/ ) {
	   $self->throw("'$value' is not a valid frame");
       }
       if( $value eq '.' ) { $value = undef; } 
       $self->{'_gsf_frame'} = $value;
  }
  return $self->{'_gsf_frame'};
}

=head2 sub_SeqFeature

 Title   : sub_SeqFeature
 Usage   : @feats = $feat->sub_SeqFeature();
 Function: Returns an array of sub Sequence Features
 Returns : An array
 Args    : none


=cut

sub sub_SeqFeature {
   my ($self) = @_;   
   return @{$self->{'_gsf_sub_array'}};
}

=head2 add_sub_SeqFeature

 Title   : add_sub_SeqFeature
 Usage   : $feat->add_sub_SeqFeature($subfeat);
           $feat->add_sub_SeqFeature($subfeat,'EXPAND')
 Function: adds a SeqFeature into the subSeqFeature array.
           with no 'EXPAND' qualifer, subfeat will be tested
           as to whether it lies inside the parent, and throw
           an exception if not.

           If EXPAND is used, the parent's start/end/strand will
           be adjusted so that it grows to accommodate the new
           subFeature
 Returns : nothing
 Args    : An object which has the SeqFeatureI interface


=cut
#'

sub add_sub_SeqFeature{
   my ($self,$feat,$expand) = @_;

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

   push(@{$self->{'_gsf_sub_array'}},$feat);

}

=head2 flush_sub_SeqFeature

 Title   : flush_sub_SeqFeature
 Usage   : $sf->flush_sub_SeqFeature
 Function: Removes all sub SeqFeature
           (if you want to remove only a subset, take
	    an array of them all, flush them, and add
            back only the guys you want)
 Example :
 Returns : none
 Args    : none


=cut

sub flush_sub_SeqFeature {
   my ($self) = @_;

   $self->{'_gsf_sub_array'} = []; # zap the array implicitly.
}


=head2 primary_tag

 Title   : primary_tag
 Usage   : $tag = $feat->primary_tag()
           $feat->primary_tag('exon')
 Function: get/set on the primary tag for a feature,
           eg 'exon'
 Returns : a string
 Args    : none


=cut

sub primary_tag {
   my ($self,$value) = @_;
   if ( defined $value ) {
       $self->{'_primary_tag'} = $value;
   }
   return $self->{'_primary_tag'};
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
   my ($self,$value) = @_;

   if( defined $value ) {
       $self->{'_source_tag'} = $value;
   }
   return $self->{'_source_tag'};
}

=head2 has_tag

 Title   : has_tag
 Usage   : $value = $self->has_tag('some_tag')
 Function: Tests wether a feature contaings a tag
 Returns : TRUE if the SeqFeature has the tag,
           and FALSE otherwise.
 Args    : The name of a tag


=cut

sub has_tag {
   my ($self, $tag) = (shift, shift);
   return exists $self->{'_gsf_tag_hash'}->{$tag};
}

=head2 add_tag_value

 Title   : add_tag_value
 Usage   : $self->add_tag_value('note',"this is a note");
 Returns : TRUE on success
 Args    : tag (string) and value (any scalar)


=cut

sub add_tag_value {
   my ($self, $tag, $value) = @_;

   if ( !defined $self->{'_gsf_tag_hash'}->{$tag} ) {
       $self->{'_gsf_tag_hash'}->{$tag} = [];
   }

   push(@{$self->{'_gsf_tag_hash'}->{$tag}},$value);
}


=head2 each_tag_value

 Title   : each_tag_value
 Usage   : @values = $gsf->each_tag_value('note');
 Function: Returns a list of all the values stored
           under a particular tag.
 Returns : A list of scalars
 Args    : The name of the tag


=cut

sub each_tag_value {
   my ($self, $tag) = @_;
   if ( ! exists $self->{'_gsf_tag_hash'}->{$tag} ) {
       $self->throw("asking for tag value that does not exist $tag");
   }

   return @{$self->{'_gsf_tag_hash'}->{$tag}};
}


=head2 all_tags

 Title   : all_tags
 Usage   : @tags = $feat->all_tags()
 Function: Get a list of all the tags in a feature
 Returns : An array of tag names
 Args    : none


=cut

sub all_tags {
   my ($self, @args) = @_;

   return keys %{$self->{'_gsf_tag_hash'}};
}


=head2 remove_tag

 Title   : remove_tag
 Usage   : $feat->remove_tag('some_tag')
 Function: removes a tag from this feature
 Returns : nothing
 Args    : tag (string)


=cut

sub remove_tag {
   my ($self, $tag) = @_;

   if ( ! exists $self->{'_gsf_tag_hash'}->{$tag} ) {
       $self->throw("trying to remove a tag that does not exist: $tag");
   }

   delete $self->{'_gsf_tag_hash'}->{$tag};
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
       $self->throw("Calling SeqFeature::Generic->seq with an argument. You probably want attach_seq");
   }

   if ( ! exists $self->{'_gsf_seq'} ) {
       return undef;
   }

   # assumming our seq object is sensible, it should not have to yank
   # the entire sequence out here.

   my $seq = $self->{'_gsf_seq'}->trunc($self->start(), $self->end());


   if ( $self->strand == -1 ) {

       # ok. this does not work well (?)
       #print STDERR "Before revcom", $seq->str, "\n";
       $seq = $seq->revcom;
       #print STDERR "After  revcom", $seq->str, "\n";
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

   return undef unless exists($self->{'_gsf_seq'});
   return $self->{'_gsf_seq'};
}


=head2 seqname

 Title   : seqname
 Usage   : $obj->seqname($newval)
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

=head2 annotation

 Title   : annotation
 Usage   : $obj->annotation($annot_obj)
 Function: 
 Example : 
 Returns : A Bio::Annotation object
 Args    : newvalue (optional)


=cut

sub annotation {
    my ($obj,$value) = @_;

    # we are smart if someone references the object and there hasn't been
    # one set yet
    if(defined $value || ! defined $obj->{'annotation'} ) {
	$value = new Bio::Annotation unless ( defined $value );
	$obj->{'annotation'} = $value;
    }
    return $obj->{'annotation'};

}

=head2 gff_format

 Title   : gff_format
 Usage   : # get:
           $gffio = $feature->gff_format();
           # set (change the default version of GFF2):
           $feature->gff_format(Bio::Tools::GFF(-gff_version => 1));
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
 Usage   : $str = $feat->gff_string;
           $str = $feat->gff_string($gff_formatter);
 Function: Provides the feature information in GFF format.

           We override this here from Bio::SeqFeatureI in order to use the
           formatter returned by gff_format().

 Returns : A string
 Args    : Optionally, an object implementing gff_string().


=cut

sub gff_string{
   my ($self,$formatter) = @_;

   $formatter = $self->gff_format() unless $formatter;
   return $formatter->gff_string($self);
}

#  =head2 slurp_gff_file
#
#   Title   : slurp_file
#   Usage   : @features = Bio::SeqFeature::Generic::slurp_gff_file(\*FILE);
#   Function: Sneaky function to load an entire file as in memory objects.
#             Beware of big files.
#
#             This method is deprecated. Use Bio::Tools::GFF instead, which can
#             also handle large files.
#
#   Example :
#   Returns :
#   Args    :
#
#  =cut

sub slurp_gff_file {
   my ($f) = @_;
   my @out;
   if ( !defined $f ) {
       die "Must have a filehandle";
   }

   Bio::Root::RootI->warn("deprecated method slurp_gff_file() called in Bio::SeqFeature::Generic. Use Bio::Tools::GFF instead.");
  
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
 Usage   : $self->_expand_region($feature);
 Function: Expand the total region covered by this feature to
           accomodate for the given feature.

           May be called whenever any kind of subfeature is added to this
           feature. add_sub_SeqFeature() already does this.
 Returns : 
 Args    : A Bio::SeqFeatureI implementing object.


=cut

sub _expand_region {
    my ($self, $feat) = @_;

    if(! $feat->isa('Bio::SeqFeatureI')) {
	$self->warn("$feat does not implement Bio::SeqFeatureI");
    }
    # if this doesn't have start/end set - forget it!
    if((! defined($self->start())) && (! defined $self->end())) {
	$self->start($feat->start());
	$self->end($feat->end());
	$self->strand($feat->strand) unless defined($self->strand());
    } else {
	my ($start, $end, $strand) = $self->union($feat);
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
 Returns : 
 Args    : 


=cut

sub _tag_value {
    my ($self, $tag, $value) = @_;

    if(defined($value) || (! $self->has_tag($tag))) {
	$self->remove_tag($tag) if($self->has_tag($tag));
	$self->add_tag_value($tag, $value);
    }
    return ($self->each_tag_value($tag))[0];
}

1;
