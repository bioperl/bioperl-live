
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

   # add it to an annotated sequence

   $annseq->add_SeqFeature($feat);



=head1 DESCRIPTION

Bio::SeqFeature::Generic is a generic implementation for the
Bio::SeqFeatureI interface, providing a simple object to provide
all the information for a feature on a sequence.

For many Features, this is all you will need to use (for example, this
is fine for Repeats in DNA sequence or Domains in protein
sequence). For other features, which have more structure, this is a
good base class to extend using inheritence to have new things: this
is what is done in the Bio::SeqFeature::Gene,
Bio::SeqFeature::Transcript and Bio::SeqFeature::Exon, which provide
well coordinated classes to represent genes on DNA sequence (for
example, you can get the protein sequence out from a transcript class).

For many Features, you want to add some piece of information, for example
a common one is that this feature is 'new' whereas other features are 'old'.
The tag system, which here is implemented using a hash can be used here.
You can use the tag system to extend the SeqFeature::Generic programmatically:
that is, you know that you have read in more information into the tag 
'mytag' which you can the retrieve. This means you do not need to know 
how to write inherieted Perl to provide more complex information on a feature,
and/or, if you do know but you donot want to write a new class every time
you need some extra piece of information, you can use the tag system
to easily store and then retrieve information.

The tag system can be written in/out of GFF format, and also into EMBL
format via AnnSeqIO::EMBL.

=head1 CONTACT

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

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Generic;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::SeqFeatureI;


@ISA = qw(Bio::Root::Object Bio::SeqFeatureI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  my $make = $self->SUPER::_initialize;

  $self->{'_gsf_tag_hash'} = {};
  $self->{'_gsf_sub_array'} = [];
  $self->{'_parse_h'} = {};

  my($start,$end,$strand,$primary,$source,$frame,$score,$tag,$gff_string) = 
      $self->_rearrange([qw(START
			    END
			    STRAND
			    PRIMARY
			    SOURCE
			    FRAME
			    SCORE
			    TAG
			    GFF_STRING
			    )],@args);

  $gff_string && $self->_from_gff_string($gff_string);
  $start && $self->start($start);
  $end   && $self->end($end);
  $strand && $self->strand($strand);
  $primary && $self->primary_tag($primary);
  $source  && $self->source_tag($source);
  $frame   && $self->frame($frame);
  $score   && $self->score($score);
  $tag     && do {
      foreach my $t ( keys %$tag ) {
	  $self->has_tag($t,$tag->{$t});
      }
  };

  # set stuff in self from @args
  return $make; # success - we hope!
}

=head2 start

 Title   : start
 Usage   : $start = $feat->start
           $feat->start(20)
 Function: Get/set on the start coordinate of the feature
 Returns : integer
 Args    : none


=cut

sub start{
   my $self = shift;

   if( @_ ) {
       my $value = shift;
       if( $value !~ /^\-?\d+/ ) {
	   $self->throw("$value is not a valid start");
       }
       $self->{'_gsf_start'} = $value
   } 

   return $self->{'_gsf_start'};
}

=head2 end

 Title   : end
 Usage   : $end = $feat->end
           $feat->end($end)
 Function: get/set on the end coordinate of the feature
 Returns : integer
 Args    : none


=cut

sub end{
   my $self = shift;

   if( @_ ) {
       my $value = shift;
       if( $value !~ /^\-?\d+/ ) {
	   $self->throw("$value is not a valid end");
       }
       $self->{'_gsf_end'} = $value
   } 

   return $self->{'_gsf_end'};
}

=head2 length

 Title   : length
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub length{
   my ($self) = @_;

   return $self->end - $self->start +1;
}


=head2 strand

 Title   : strand
 Usage   : $strand = $feat->strand()
           $feat->strand($strand)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none


=cut

sub strand{
   my $self = shift;

   if( @_ ) {
       my $value = shift;
       if( $value eq '+' ) { $value = 1; }
       if( $value eq '-' ) { $value = -1; }
       if( $value eq '.' ) { $value = 0; }
       
       if( $value != -1 && $value != 1 && $value != 0 ) {
	   $self->throw("$value is not a valid strand info");
       }
       $self->{'_gsf_strand'} = $value
   } 

   return $self->{'_gsf_strand'};
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
  my $self = shift;
  
  if(@_) {
       my $value = shift;
       if( $value !~ /^[+-]?\d+\.?\d*(e-\d+)?/ ) {
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
 Returns : 0,1,2
 Args    : none if get, the new value if set


=cut

sub frame {
  my $self = shift;
  
  if(@_) {
       my $value = shift;
       if( $value != 1 && $value != 2 && $value != 3 ) {
	   $self->throw("'$value' is not a valid frame");
       }
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

sub sub_SeqFeature{
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

sub add_sub_SeqFeature{
   my ($self,$feat,$expand) = @_;

   if( !$feat->isa('Bio::SeqFeatureI') ) {
       $self->warn("$feat does not implement Bio::SeqFeatureI. Will add it anyway, but beware...");
   }

   if( $expand eq 'EXPAND' ) {
       # if this doesn't have start/end set - forget it!
       if( !defined $self->start && !defined $self->end ) {
	   $self->start($feat->start());
	   $self->end($feat->end());
	   $self->strand($feat->strand);
       } else {
	   my ($start,$end,$strand) = $self->union($feat);
	   $self->start($start);
	   $self->end($end);
	   $self->strand($strand);
       }
   } else {
       if( !$self->contains($feat) ) {
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

sub primary_tag{
   my $self = shift;
   if( @_ ) {
       $self->{'_primary_tag'} = shift;
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

sub source_tag{
   my $self = shift;

   if( @_ ) {
       $self->{'_source_tag'} = shift;
   }
   return $self->{'_source_tag'};
}

=head2 has_tag

 Title   : has_tag
 Usage   : $value = $self->has_tag('some_tag')
 Function: Returns the value of the tag (undef if 
           none)
 Returns : 
 Args    :


=cut

sub has_tag{
   my ($self,$tag) = (shift, shift);

   return exists $self->{'_gsf_tag_hash'}->{$tag};
}

=head2 add_tag_value

 Title   : add_tag_value
 Usage   : $self->add_tag_value('note',"this is a note");
 Returns : nothing
 Args    : tag (string) and value (any scalar)


=cut

sub add_tag_value{
   my ($self,$tag,$value) = @_;

   if( !defined $self->{'_gsf_tag_hash'}->{$tag} ) {
       $self->{'_gsf_tag_hash'}->{$tag} = [];
   }

   push(@{$self->{'_gsf_tag_hash'}->{$tag}},$value);
}

=head2 each_tag_value

 Title   : each_tag_value
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_tag_value {
   my ($self,$tag) = @_;

   return @{$self->{'_gsf_tag_hash'}->{$tag}};
}


=head2 all_tags

 Title   : all_tags
 Usage   : @tags = $feat->all_tags()
 Function: gives all tags for this feature
 Returns : an array of strings
 Args    : none


=cut

sub all_tags{
   my ($self,@args) = @_;

   return keys %{$self->{'_gsf_tag_hash'}};
}

=head2 attach_seq

 Title   : attach_seq
 Usage   : $sf->attach_seq($seq)
 Function: Attaches a Bio::Seq object to this feature. This
           Bio::Seq object is for the *entire* sequence: ie
           from 1 to 10000
 Example :
 Returns : 
 Args    :


=cut

sub attach_seq{
   my ($self,$seq) = @_;

   if( !defined $seq  || !ref $seq || ! $seq->isa("Bio::Seq") ) {
       $self->throw("Must attach Bio::Seq objects to SeqFeatures");
   }

   $self->{'_gsf_seq'} = $seq;

   # attach to sub features if they want it

   foreach my $sf ( $self->sub_SeqFeature() ) {
       if( $sf->can("attach_seq") ) {
	   $sf->attach_seq($seq);
       }
   }
}

=head2 seq

 Title   : seq
 Usage   : $tseq = $sf->seq()
 Function: returns the truncated sequence (if there) for this
 Example :
 Returns : 
 Args    :


=cut

sub seq{
   my ($self,$arg) = @_;

   if( defined $arg ) {
       $self->throw("Calling SeqFeature::Generic->seq with an argument. You probably want attach_seq");
   }

   if( ! exists $self->{'_gsf_seq'} ) {
       return undef;
   }

   # assumming our seq object is sensible, it should not have to yank
   # the entire sequence out here.

   my $seq = $self->{'_gsf_seq'}->trunc($self->start(),$self->end());


   if( $self->strand == -1 ) {

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

sub entire_seq{
   my ($self) = @_;

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

sub seqname{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_gsf_seqname'} = $value;
    }
    return $obj->{'_gsf_seqname'};

}

=head2 _from_gff_string

 Title   : _from_gff_string
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _from_gff_string{
   my ($self,$string) = @_;

   my($seqname,$source,$primary,$start,$end,$score,$strand,$frame,@group) = split(/\s+/,$string);
   if( !defined $frame ) {
       $self->throw("[$string] does not look like GFF to me");
   }
   $self->seqname($seqname);
   $self->source_tag($source);
   $self->primary_tag($primary);
   $self->start($start);
   $self->end($end);
   if( $score eq '.' ) {
       $self->score(undef);
   } else {
       $self->score($score);
   }
   if( $strand eq '-' ) { $self->strand(-1); }
   if( $strand eq '+' ) { $self->strand(1); }   
   if( $strand eq '.' ) { $self->strand(0); }
   
}

=head2 _parse

 Title   : _parse
 Usage   :
 Function: Parsing hints
 Example :
 Returns : 
 Args    :


=cut

sub _parse{
   my ($self) = @_;

   return $self->{'_parse_h'};
}



