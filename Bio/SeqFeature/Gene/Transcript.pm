# $Id$
#
# BioPerl module for Bio::SeqFeature::Gene::Transcript
#
# Cared for by Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Gene::Transcript - A feature representing a transcript

=head1 SYNOPSIS

See documentation of methods.

=head1 DESCRIPTION

A feature representing a transcript.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Hilmar Lapp

Email hlapp@gmx.net

Modifications by Paul Edlefsen paul@systemsbiology.org

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::SeqFeature::Gene::Transcript;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::SeqFeature::Gene::TranscriptI;
use Bio::SeqFeature::Generic;
use Bio::PrimarySeq;

@ISA = qw( Bio::SeqFeature::Generic Bio::SeqFeature::Gene::TranscriptI );

sub new {
  my ( $caller, @args ) = @_;
  my $self = $caller->SUPER::new( @args );
  my ( $primary ) = $self->_rearrange( [ qw( PRIMARY ) ], @args );

  $primary |= 'transcript';
  $self->primary_tag( $primary );
  $self->strand( 0 ) unless defined( $self->strand() );
  return $self;
} # new(..)

=head2 promoters

 Title   : promoters()
 Usage   : @proms = $transcript->promoters();
 Function: Get the promoter features/sites of this transcript. 

           Note that OO-modeling of regulatory elements is not stable yet.
           This means that this method might change or even disappear in a
           future release. Be aware of this if you use it.

 Returns : An array of Bio::SeqFeatureI implementing objects representing the
           promoter regions or sites.
 Args    : 


=cut

sub promoters {
    my ($self) = @_;
    return $self->get_feature_type('Bio::SeqFeature::Gene::Promoter');
}

=head2 add_promoter

 Title   : add_promoter()
 Usage   : $transcript->add_promoter($feature);
 Function: Add a promoter feature/site to this transcript.


           Note that OO-modeling of regulatory elements is not stable yet.
           This means that this method might change or even disappear in a
           future release. Be aware of this if you use it.

 Returns : 
 Args    : A Bio::SeqFeatureI implementing object.

=cut

sub add_promoter {
  my $self = shift;
  my ( $fea ) = @_;
  $self->_add( $fea, 'Bio::SeqFeature::Gene::Promoter' );
} # add_promorter(..)

=head2 flush_promoters

 Title   : flush_promoters()
 Usage   : $transcript->flush_promoters();
 Function: Remove all promoter features/sites from this transcript.

           Note that OO-modeling of regulatory elements is not stable yet.
           This means that this method might change or even disappear in a
           future release. Be aware of this if you use it.

 Returns : the removed features as a list
 Args    : none


=cut

sub flush_promoters {
  my $self = shift;
  return $self->_flush( 'Bio::SeqFeature::Gene::Promoter' );
}

=head2 exons

 Title   : exons()
 Usage   : @exons = $gene->exons();
           ($inital_exon) = $gene->exons('Initial');
 Function: Get all exon features or all exons of specified type of this 
           transcript.

           Exon type is treated as a case-insensitive regular expression and 
           is optional. For consistency, use only the following types: 
           initial, internal, terminal.

 Returns : An array of Bio::SeqFeature::Gene::ExonI implementing objects.
 Args    : An optional string specifying the primary_tag of the feature.


=cut

sub exons {
  my $self = shift;
  my ( $type ) = @_;
  return
    $self->get_unordered_feature_type(
      'Bio::SeqFeature::Gene::ExonI', 
      $type
    );
} # exons(..)

=head2 exons_ordered

 Title   : exons_ordered
 Usage   : @exons = $gene->exons_ordered();
           @exons = $gene->exons_ordered("Internal");
 Function: Get an ordered list of all exon features or all exons of specified
           type of this transcript.

           Exon type is treated as a case-insensitive regular expression and 
           is optional. For consistency, use only the following types:

 Returns : An array of Bio::SeqFeature::Gene::ExonI implementing objects.
 Args    : An optional string specifying the primary_tag of the feature.

=cut

sub exons_ordered { 
  my $self = shift;
  my ( $type ) = @_;
  return
    $self->get_feature_type(
      'Bio::SeqFeature::Gene::ExonI', 
      $type
    );

  ## TODO: REMOVE
  #my @exons = 
  #  $self->get_feature_type(
  #    'Bio::SeqFeature::Gene::ExonI', 
  #    $type
  #  );
  #print STDOUT "exons, ordered, of $self are: ( ";
  #for( my $i = 0; $i < scalar( @exons ); $i++ ) {
  #  unless( $i == 0 ) {
  #    print STDOUT ', ';
  #  }
  #  print STDOUT $exons[ $i ]->toRelRangeString( 'both', 'plus' );
  #}
  #print STDOUT " ).\n";
  #return @exons;
}

=head2 add_exon

 Title   : add_exon()
 Usage   : $transcript->add_exon($exon,'initial');
 Function: Add a exon feature to this transcript.

           The second argument denotes the type of exon. Mixing exons with and
           without a type is likely to cause trouble in exons(). Either
           leave out the type for all exons or for none.

           Presently, the following types are known: initial, internal, 
           terminal, utr, utr5prime, and utr3prime (all case-insensitive).
           UTR should better be added through utrs()/add_utr().

           If you wish to use other or additional types, you will almost
           certainly have to call exon_type_sortorder() in order to replace
           the default sort order, or mrna(), cds(), protein(), and exons()
           may yield unexpected results.

 Returns : 
 Args    : A Bio::SeqFeature::Gene::ExonI implementing object.
           A string indicating the type of the exon (optional).


=cut

sub add_exon {
  my $self = shift;
  my ( $fea, $pri ) = @_;
  unless( $fea->isa( 'Bio::SeqFeature::Gene::ExonI' ) ) {
    $self->throw( "$fea does not implement Bio::SeqFeature::Gene::ExonI" );
  }
  $fea->primary_tag( $pri ) if defined( $pri );
  $self->_add( $fea, 'Bio::SeqFeature::Gene::Exon' );
} # add_exon(..)

=head2 flush_exons

 Title   : flush_exons()
 Usage   : $transcript->flush_exons();
           $transcript->flush_exons('terminal');
 Function: Remove all or a certain type of exon features from this transcript.

           See add_exon() for documentation about types.

           Calling without a type will not flush UTRs. Call flush_utrs() for
           this purpose.
 Returns : the deleted features as a list
 Args    : A string indicating the type of the exon (optional).

=cut

sub flush_exons {
  my $self = shift;
  my ( $type ) = @_;
  return $self->_flush( 'Bio::SeqFeature::Gene::Exon', $type );
}

=head2 introns

 Title   : introns()
 Usage   : @introns = $gene->introns();
 Function: Get all intron features this gene structure.

           Note that this implementation generates these features
           on-the-fly, that is, it simply treats all regions between
           exons as introns, assuming that exons do not overlap. A
           consequence is that a consistent correspondence between the
           elements in the returned array and the array that exons()
           returns will exist only if the exons are properly sorted
           within their types (forward for plus- strand and reverse
           for minus-strand transcripts). To ensure correctness the
           elements in the array returned will always be sorted.

 Returns : An array of Bio::SeqFeature::Gene::Intron objects representing
           the intron regions.
 Args    : 


=cut

sub introns {
  my $self = shift;
  my @introns = ();
  my @exons = $self->exons();

  # if there's 1 or fewer exons we're done
  return () unless( $#exons > 0 );

  # record strand and order (a minus-strand transcript is likely to have
  # the exons stacked in reverse order)
  my $strand;
  foreach my $exon ( @exons ) {
    $strand = $exon->strand();
    last if $strand; # we're done if we've got 1 or -1
  }
  my $rev_order = ( ( $exons[ 0 ]->end() < $exons[ 1 ]->start() ) ? 0 : 1 );

  # Make sure exons are sorted. Because we assume they don't overlap, we
  # simply sort by start position.
  if( ( !defined( $strand ) ) || ( $strand != -1 ) || ( !$rev_order ) ) {
    # always sort forward for plus-strand transcripts, and for negative-
    # strand transcripts that appear to be unsorted or forward sorted
    @exons =
      map  { $_->[ 0 ] }
      sort { $a->[ 1 ] <=> $b->[ 1 ] }
      map  { [ $_, $_->start() ] }
      @exons;
  } else {
    # sort in reverse order for transcripts on the negative strand and
    # found to be in reverse order
    @exons =
      map  { $_->[ 0 ] }
      sort { $b->[ 1 ] <=> $a->[ 1 ] }
      map  { [ $_, $_->start() ] }
      @exons;
  }
  # loop over all intervening gaps
  my $seq = $self->entire_seq();
  my ( $start, $end, $intron );
  for( my $i = 0; $i < $#exons; $i++ ) { # all but the last exon.  TODO: Why?
    if( defined( $exons[ $i ]->strand() ) &&
        ( ( $exons[ $i ]->strand() * $strand ) < 0 )
      ) {
      $self->throw( "Transcript mixes plus and minus strand exons. ".
                    "Computing introns makes no sense then." );
    }
    $start = $exons[ $i     + $rev_order ]->end()   + 1;   # $i or $i+1
    $end   = $exons[ $i + 1 - $rev_order ]->start() - 1;   # $i+1 or $i
    $intron = Bio::SeqFeature::Gene::Intron->new(
      '-start'   => $start,
      '-end'     => $end,
      '-strand'  => $strand,
      '-primary' => 'intron',
      '-source'  => ref( $self )
    );
    $intron->seq_id( $self->seq_id() );
    if( $seq ) {
      $intron->attach_seq( $seq );
    }
    push( @introns, $intron );
  }
  return @introns;
} # introns(..)

=head2 poly_A_site

 Title   : poly_A_site()
 Usage   : $polyAsite = $transcript->poly_A_site();
 Function: Get/set the poly-adenylation feature/site of this transcript.
 Returns : A Bio::SeqFeatureI implementing object representing the
           poly-adenylation region.
 Args    : A Bio::SeqFeatureI implementing object on set, or FALSE to flush
           a previously set object.


=cut

sub poly_A_site {
  my $self = shift;
  my ( $fea ) = @_;
  if( $fea ) {
    $self->_add( $fea, 'Bio::SeqFeature::Gene::Poly_A_site' );
  }
  return
    ( $self->get_feature_type( 'Bio::SeqFeature::Gene::Poly_A_site' ) )[ 0 ];
} # poly_A_site(..)

=head2 utrs

 Title   : utrs()
 Usage   : @utr_sites = $transcript->utrs('utr3prime');
           @utr_sites = $transcript->utrs('utr5prime');
           @utr_sites = $transcript->utrs();
 Function: Get the features representing untranslated regions (UTR) of this
           transcript.

           You may provide an argument specifying the type of UTR. Currently
           the following types are recognized: utr5prime utr3prime for UTR on the
           5' and 3' end of the CDS, respectively.

 Returns : An array of Bio::SeqFeature::Gene::UTR objects
           representing the UTR regions or sites.
 Args    : Optionally, either utr3prime, or utr5prime for the the type of UTR
           feature.


=cut

sub utrs {
  my $self = shift;
  my ( $type ) = @_;
  return $self->get_feature_type( 'Bio::SeqFeature::Gene::UTR', $type );
} # utrs(..)

=head2 add_utr

 Title   : add_utr()
 Usage   : $transcript->add_utr($utrobj, 'utr3prime');
           $transcript->add_utr($utrobj);
 Function: Add a UTR feature/site to this transcript.

           The second parameter is optional and denotes the type of the UTR
           feature. Presently recognized types include 'utr5prime' and 'utr3prime'
           for UTR on the 5' and 3' end of a gene, respectively.

           Calling this method is the same as calling 
           add_exon($utrobj, 'utr'.$type). In this sense a UTR object is a
           special exon object, which is transcribed, not spliced out, but
           not translated.

           Note that the object supplied should return FALSE for is_coding().
           Otherwise cds() and friends will become confused.

 Returns : 
 Args    : A Bio::SeqFeature::Gene::UTR implementing object.


=cut

sub add_utr {
  my $self = shift;
  my ( $fea, $pri ) = @_;
  $fea->primary_tag( $pri ) if defined( $pri );
  $self->_add( $fea, 'Bio::SeqFeature::Gene::UTR' );
} # add_utr(..)

=head2 flush_utrs

 Title   : flush_utrs()
 Usage   : $transcript->flush_utrs();
           $transcript->flush_utrs('utr3prime');
 Function: Remove all or a specific type of UTR features/sites from this
           transcript.

           Cf. add_utr() for documentation about recognized types.
 Returns : a list of the removed features
 Args    : Optionally a string denoting the type of UTR feature.


=cut

sub flush_utrs {
  my $self = shift;
  my ( $type ) = @_;
  return $self->_flush( 'Bio::SeqFeature::Gene::UTR', $type );
} # flush_utrs(..)

=head2 sub_SeqFeature

 Title   : sub_SeqFeature
 Usage   : @feats = $transcript->sub_SeqFeature();
 Function: Returns an array of all subfeatures.

           This method is defined in Bio::SeqFeatureI. We override this here
           to include the exon etc features.

 Returns : An array Bio::SeqFeatureI implementing objects.
 Args    : none


=cut

sub sub_SeqFeature {
  my $self = shift;

  # get what the parent already has
  my @features = $self->SUPER::sub_SeqFeature();
  # add the features we have in addition
  push( @features, $self->exons() ); # this includes UTR features
  push( @features, $self->promoters() );
  if( my $poly_A_site = $self->poly_A_site() ) {
    push( @features, $poly_A_site );
  }
  return @features;
} # sub_SeqFeature(..)

=head2 flush_sub_SeqFeature

 Title   : flush_sub_SeqFeature
 Usage   : $transcript->flush_sub_SeqFeature();
           $transcript->flush_sub_SeqFeature(1);
 Function: Removes all subfeatures.

           This method is overridden from Bio::SeqFeature::Generic to flush
           all additional subfeatures like exons, promoters, etc., which is
           almost certainly not what you want. To remove only features added
           through $transcript->add_sub_SeqFeature($feature) pass any
           argument evaluating to TRUE.

 Example :
 Returns : none
 Args    : Optionally, an argument evaluating to TRUE will suppress flushing
           of all transcript-specific subfeatures (exons etc.).

=cut

sub flush_sub_SeqFeature {
  my $self = shift;
  my ( $fea_only ) = @_;

  $self->SUPER::flush_sub_SeqFeature();
  unless( $fea_only ) {
    $self->flush_promoters();
    $self->flush_exons();
    $self->flush_utrs();
    $self->poly_A_site( 0 ); ## TODO: Is passing 0 working to flush it?
  }
} # flush_sub_SeqFeature(..)

=head2 cds

 Title   : cds
 Usage   : $seq = $transcript->cds();
 Function: Returns the CDS (coding sequence) as defined by the exons
           of this transcript and the attached sequence.

           If no sequence is attached this method will return undef.

           Note that the implementation provided here returns a
           concatenation of all coding exons, thereby assuming that
           exons do not overlap.

           Note also that you cannot set the CDS via this method. Set
           a single CDS feature as a single exon, or derive your own
           class if you want to store a predicted CDS.

 Example :
 Returns : A Bio::PrimarySeqI implementing object.
 Args    : 

=cut

sub cds {
  my $self = shift;

  # this is always sorted properly according to strand
  my @exons = $self->exons_ordered();
  return undef unless @exons;

  # record strand (a minus-strand transcript must have the exons sorted in
  # reverse order)
  my $common_strand;
  foreach my $exon ( @exons ) {
    if( !$common_strand && defined( $exon->strand() ) ) {
      $common_strand = $exon->strand();
    }
    if( $exon->strand() &&
        ( ( $exon->strand() * $common_strand ) < 0 )
      ) {
      $self->throw( "Transcript mixes coding exons on plus and minus ".
                    "strand. This makes no sense." );
    }
  }
  my $cds = $self->_make_cds( @exons );
  return undef unless $cds;
  return
    Bio::PrimarySeq->new(
      '-id'       => $self->seq_id(),
      '-seq'      => $cds,
      '-alphabet' => 'dna'
    );
} # cds(..)

=head2 protein

 Title   : protein()
 Usage   : $protein = $transcript->protein();
 Function: Get the protein encoded by the transcript as a sequence object.

           The implementation provided here simply calls translate() on the
           object returned by cds().

 Returns : A Bio::PrimarySeqI implementing object.
 Args    : 


=cut

sub protein {
  my $self = shift;
  my $seq = $self->cds();
  if( $seq ) {
    return $seq->translate();
  } else {
    return undef;
  }
} # protein(..)

=head2 mrna

 Title   : mrna()
 Usage   : $mrna = $transcript->mrna();
 Function: Get the mRNA of the transcript as a sequence object.

           The difference to cds() is that the sequence object returned by
           this methods will also include UTR and the poly-adenylation site,
           but not promoter sequence (TBD).

           HL: do we really need this method?

 Returns : A Bio::PrimarySeqI implementing object.
 Args    : 

=cut

sub mrna {
  my $self = shift;

  # get the coding part
  my $seq = $self->cds();
  unless( $seq ) {
    $seq =
      Bio::PrimarySeq->new(
        '-id' => $self->abs_range(),
        '-alphabet' => "rna",
	'-seq' => ""
      );
  }
  # get and add UTR sequences
  my $mrna = "";
  foreach my $elem ( $self->utrs( 'utr5prime' ) ) {
    $mrna .= $elem->seq()->seq();
  }
  $seq->seq( $mrna . $seq->seq() );
  $mrna = "";
  foreach my $elem ( $self->utrs( 'utr3prime' ) ) {
    $mrna .= $elem->seq()->seq();
  }
  $seq->seq( $seq->seq() . $mrna );
  if( my $poly_A_site = $self->poly_A_site() ) {
    $seq->seq( $seq->seq() . $poly_A_site->seq()->seq() );
  }
  if( $seq->length() == 0 ) {
    return undef;
  } else {
    return $seq;
  }
} # mrna(..)

## This method doesn't seem to be used at all by anyone.  Can we remove it? -PE
sub _get_typed_keys {
  my $self = shift;
  my ( $keyprefix, $type ) = @_;

  # make case-insensitive
  $type = ( $type ? lc( $type ) : '' );

  # pull out all feature types that exist and match
  return grep { /^_$keyprefix$type/i; } ( keys( %{ $self } ) );
} # _get_typed_keys(..)

sub _make_cds {
  my $self = shift;
  my @exons = @_;

  my $cds = "";
  foreach my $exon ( @exons ) {

    my $exon_seq = $exon->seq();
    next unless( defined( $exon_seq ) && $exon->is_coding() );

    my $phase = length( $cds ) % 3;
    # let's check the simple case 
    if( !defined( $exon->frame() ) || ( $phase == $exon->frame() ) ) {
      # this one fits exactly, or frame of the exon is undefined (should
      # we warn about that?); we bypass the $exon->cds() here (hmm,
      # not very clean style, but I don't see where this screws up)
      ## TODO: REMOVE
      #warn "(A) Exon $exon has frame ".$exon->frame().".  The cds so far has phase $phase.  The exon isa ".ref( $exon ).".  Adding its seq, which when translated is\n'".$exon_seq->translate()->seq()."'";
      $cds .= $exon_seq->seq();
    } else {
      $exon_seq = $exon->cds(); # now $seq is guaranteed to be in frame 0
      ## TODO: REMOVE
      #warn "(B) Exon $exon has frame ".$exon->frame().".  The cds so far has phase $phase.  Adding its seq, which when translated is \n'".$exon_seq->translate()->seq()."'";
      # this one is probably from exon shuffling and needs some work
      next unless $exon_seq;

      my $seq = $exon_seq->seq();
      # adjustment needed?
      if( $phase > 0 ) {
        # how many Ns can we chop off the piece to be added?
        my $n_crop = 0;
        if( $seq =~ /^(n+)/i ) {
          $n_crop = length( $1 );
        }
        if( $n_crop >= $phase ) {
          # chop off to match the phase
          $seq = substr( $seq, $phase );
        } else {
          # fill in Ns
          $seq = ( 'n' x ( 3 - $phase ) ) . $seq;
        }
      }
      $cds .= $seq;
    }
  }
  return $cds;
} # _make_cds(..)

=head2 features_ordered

 Title   : features_ordered
 Usage   : my @features=$transcript->features_ordered;
 Function: returns all the features associated with this transcript,
           in order by feature start, according to strand
 Returns : a list of SeqFeatureI implementing objects
 Args    : none


=cut

sub features_ordered {
  my $self = shift;
  return $self->_stranded_sort( $self->features( @_ ) );
}

sub get_unordered_feature_type {
  my $self = shift;
  my ( $type, $pri ) = @_;
  my @list;

  foreach my $feature ( $self->features() ) {
    if( $feature->isa( $type ) ) {
      if( $pri && ( $feature->primary_tag() !~ /$pri/i ) ) {
        next;
      }
      push( @list, $feature );
    }
  }
  return @list;
} # get_unordered_feature_type(..)

sub get_feature_type {
  my $self = shift;
  return $self->_stranded_sort( $self->get_unordered_feature_type( @_ ) );
}

#This was fixed by Gene Cutler - the indexing on the list being reversed
#fixed a bad bug.  Thanks Gene!
sub _flush {
  my $self = shift;
  my ( $type, $pri ) = @_;

  my @features = $self->features();
  my @cut;
  for( reverse ( 0..$#features ) ) {
    if( $features[ $_ ]->isa( $type ) ) {
      if( $pri && ( $features[ $_ ]->primary_tag() !~ /$pri/i ) ) {
        next;
      }
      $self->remove_feature( $features[ $_ ] );
      push( @cut, $features[ $_ ] );
    }
  }
  return reverse @cut;
} # _flush(..)

sub _add {
  my $self = shift;
  my @features_and_types = @_;
  require Bio::SeqFeature::Gene::Promoter;
  require Bio::SeqFeature::Gene::UTR;
  require Bio::SeqFeature::Gene::Exon;
  require Bio::SeqFeature::Gene::Intron;
  require Bio::SeqFeature::Gene::Poly_A_site;
  
  my @features;
  my ( $feature, $type );
  for( my $i = 0; $i < scalar( @features_and_types ); $i += 2 ) {
    $feature = $features_and_types[ $i ];
    $type = $features_and_types[ $i + 1 ];
    unless( $feature->isa( 'Bio::SeqFeature::SegmentI' ) ) {
      $self->throw( "$feature does not implement Bio::SeqFeature::SegmentI" );
    }
    unless( $feature->isa( $type ) ) {
      $feature = $self->_new_of_type( $feature, $type );
    }
    unless( $self->strand() ) {
      $self->strand( $feature->strand() );
    } elsif( ( $self->strand() * $feature->strand() ) == -1 ) {
      $self->throw( "$feature is on opposite strand from \$self" );
    }
    if(
       defined( $self->entire_seq() ) &&
       ( !defined( $feature->entire_seq() ) ) &&
       $feature->can( 'attach_seq' )
      ) {
      $feature->attach_seq( $self->entire_seq() );
    }
    push( @features, $feature );
  }
  
  $self->add_features( @features );
  
  ## TODO: What we really want here is to trigger our grandest
  ## ancestor's adjust_bounds method.
  if( defined $self->parent() ) {
    $self->parent()->adjust_bounds( $self );
  } else {
    $self->adjust_bounds( @features );
  }
  
  return 1;
} # _add(..)

## TODO: Can we replace this with the -sort flag?
sub _stranded_sort {
  my $self = shift;
  my @features = @_;

  my $common_strand;
  foreach my $feature ( @features ) {
    if( $feature->abs_strand() ) {
      # defined and != 0
      unless( $common_strand ) {
        $common_strand = $feature->abs_strand();
      }
      if( ( $feature->abs_strand() * $common_strand ) < 0 ) {
        undef $common_strand;
        last;
      }
    }
  }
  if( defined( $common_strand ) && ( $common_strand < 0 ) ) {  #reverse strand
    ## TODO: REMOVE
    #warn "Hey, it's on the reverse strand.. stranded sort..";
    return map  { $_->[ 0 ] }
           sort { $b->[ 1 ] <=> $a->[ 1 ] }
           map  { [ $_, $_->abs_high( 'plus' ) ] }
           @features;
  } else {                                       #undef or forward strand
    ## TODO: REMOVE
    #warn "Hey, it's on the plus strand.. stranded sort.. common strand is $common_strand";
    return map  { $_->[ 0 ] }
           sort { $a->[ 1 ] <=> $b->[ 1 ] }
           map  { [ $_, $_->abs_start( 'plus' ) ] }
           @features;
  }
} # _stranded_sort;

sub _new_of_type {
  my $self = shift;
  my ( $fea, $type, $pri ) = @_;

  my $primary;
  if( $pri ) {
    # can set new primary tag if desired
    $primary = $pri;
  } else {
    # or else primary is just end of type string
    ( $primary ) = ( $type =~ /.*::(.+)/ );
  }
  bless $fea, $type;
  $fea->primary_tag( $primary );
  return $fea;
} # _new_of_type(..)

1;
