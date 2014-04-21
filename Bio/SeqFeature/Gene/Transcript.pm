#
# BioPerl module for Bio::SeqFeature::Gene::Transcript
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

  # See documentation of methods.

=head1 DESCRIPTION

A feature representing a transcript.

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

=head1 AUTHOR - Hilmar Lapp

Email hlapp@gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::SeqFeature::Gene::Transcript;
use strict;


use Bio::PrimarySeq;

use base qw(Bio::SeqFeature::Generic Bio::SeqFeature::Gene::TranscriptI);

sub new {
    my ($caller, @args) = @_;
    my $self = $caller->SUPER::new(@args);
    $self->_register_for_cleanup(\&transcript_destroy);
    my ($primary) = $self->_rearrange([qw(PRIMARY)],@args);

    $primary = 'transcript' unless $primary;
    $self->primary_tag($primary);
    $self->strand(0) if(! defined($self->strand()));
    return $self;
}


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
    my ($self, $fea) = @_;
    $self->_add($fea,'Bio::SeqFeature::Gene::Promoter');
}

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
    my ($self) = @_;
    return $self->_flush('Bio::SeqFeature::Gene::Promoter');
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
    my ($self, $type) = @_;
    return $self->get_unordered_feature_type('Bio::SeqFeature::Gene::ExonI', 
					     $type);
}

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
    my ($self,$type) = @_;
    return $self->get_feature_type('Bio::SeqFeature::Gene::ExonI', $type);
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
    my ($self, $fea, $type) = @_;
    if(! $fea->isa('Bio::SeqFeature::Gene::ExonI') ) {
	$self->throw("$fea does not implement Bio::SeqFeature::Gene::ExonI");
    }
    $self->_add($fea,'Bio::SeqFeature::Gene::Exon', $type);
}

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
    my ($self, $type) = @_;
    return $self->_flush('Bio::SeqFeature::Gene::Exon',$type);
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
    my ($self) = @_;
    my @introns = ();
    my @exons = $self->exons();
    my ($strand, $rev_order);

    # if there's 1 or less exons we're done
    return () unless($#exons > 0);
    # record strand and order (a minus-strand transcript is likely to have
    # the exons stacked in reverse order)
    foreach my $exon (@exons) {
	$strand = $exon->strand();
	last if $strand; # we're done if we've got 1 or -1
    }
    $rev_order = ($exons[0]->end() < $exons[1]->start() ? 0 : 1);

    # Make sure exons are sorted. Because we assume they don't overlap, we
    # simply sort by start position.
    if((! defined($strand)) || ($strand != -1) || (! $rev_order)) {
	# always sort forward for plus-strand transcripts, and for negative-
	# strand transcripts that appear to be unsorted or forward sorted
        @exons = map { $_->[0] } sort { $a->[1] <=> $b->[1] } 
	         map { [ $_, $_->start * ($_->strand || 1)] } @exons;
    } else {
	# sort in reverse order for transcripts on the negative strand and
	# found to be in reverse order
        @exons = map { $_->[0] } sort { $b->[1] <=> $a->[1] } map { [ $_, $_->start()] } @exons;
    }
    # loop over all intervening gaps
    while ((my $exonA = shift (@exons)) &&(my $exonB = shift(@exons))){
       my $intron = Bio::SeqFeature::Gene::Intron->new(-primary=>'intron');
       $intron->upstream_Exon($exonA);
       $intron->downstream_Exon($exonB);
       $intron->attach_seq($self->entire_seq) if $self->entire_seq;
       unshift(@exons,$exonB);
       push @introns,$intron;
    }
    return @introns;
}

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
    my ($self, $fea) = @_;
    if ($fea) {
	$self->_add($fea,'Bio::SeqFeature::Gene::Poly_A_site');
    }
    return ($self->get_feature_type('Bio::SeqFeature::Gene::Poly_A_site'))[0];
}

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
    my ($self, $type) = @_;
    return $self->get_feature_type('Bio::SeqFeature::Gene::UTR',$type);

}

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
    my ($self, $fea, $type) = @_;
    $self->_add($fea,'Bio::SeqFeature::Gene::UTR',$type);
}

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
    my ($self, $type) = @_;
    return $self->_flush('Bio::SeqFeature::Gene::UTR',$type);
}

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
   my ($self) = @_;   
   my @feas;
   
   # get what the parent already has
   @feas = $self->SUPER::sub_SeqFeature();
   # add the features we have in addition
   push(@feas, $self->exons()); # this includes UTR features
   push(@feas, $self->promoters());
   push(@feas, $self->poly_A_site()) if($self->poly_A_site());
   return @feas;
}

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
   my ($self,$fea_only) = @_;

   $self->SUPER::flush_sub_SeqFeature();
   if(! $fea_only) {
       $self->flush_promoters();
       $self->flush_exons();
       $self->flush_utrs();
       $self->poly_A_site(0);
   }
}


=head2 cds

 Title   : cds
 Usage   : $seq = $transcript->cds();
 Function: Returns the CDS (coding sequence) as defined by the exons
           of this transcript and the attached sequence.

           If no sequence is attached this method will return false.

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
    my ($self) = @_;
    my @exons = $self->exons_ordered();  #this is always sorted properly according to strand
    my $strand;

    return  unless(@exons);
    # record strand (a minus-strand transcript must have the exons sorted in
    # reverse order)
    foreach my $exon (@exons) {
	if(defined($exon->strand()) && (! $strand)) {
	    $strand = $exon->strand();
	}
	if($exon->strand() && (($exon->strand() * $strand) < 0)) {
	    $self->throw("Transcript mixes coding exons on plus and minus ".
			 "strand. This makes no sense.");
	}
    }
    my $cds = $self->_make_cds(@exons);
    return unless $cds;
    return Bio::PrimarySeq->new('-id' => $self->seq_id(),
				'-seq' => $cds,
				'-alphabet' => "dna");
}

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
    my ($self) = @_;
    my $seq;

    $seq = $self->cds();
    return $seq->translate() if $seq;
    return;
}

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
    my ($self) = @_;
    my ($seq, $mrna, $elem);

    # get the coding part
    $seq = $self->cds();
    if(! $seq) {
	$seq = Bio::PrimarySeq->new('-id' => $self->seq_id(),
				    '-alphabet' => "rna",
				    '-seq' => "");
    }
    # get and add UTR sequences
    $mrna = "";
    foreach $elem ($self->utrs('utr5prime')) {
	$mrna .= $elem->seq()->seq();
    }
    $seq->seq($mrna . $seq->seq());
    $mrna = "";
    foreach $elem ($self->utrs('utr3prime')) {
	$mrna .= $elem->seq()->seq();
    }
    $seq->seq($seq->seq() . $mrna);
    if($self->poly_A_site()) {
	$seq->seq($seq->seq() . $self->poly_A_site()->seq()->seq());
    }
    return if($seq->length() == 0);
    return $seq;
}

sub _get_typed_keys {
    my ($self, $keyprefix, $type) = @_;
    my @keys = ();
    my @feas = ();

    # make case-insensitive
    $type = ($type ? lc($type) : "");
    # pull out all feature types that exist and match
    @keys = grep { /^_$keyprefix$type/i; } (keys(%{$self}));
    return @keys;
}

sub _make_cds {
    my ($self,@exons) = @_;
    my $cds = "";

    foreach my $exon (@exons) {
	next if((! defined($exon->seq())) || (! $exon->is_coding()));
	my $phase = length($cds) % 3;
	# let's check the simple case 
	if((! defined($exon->frame())) || ($phase == $exon->frame())) {
	    # this one fits exactly, or frame of the exon is undefined (should
	    # we warn about that?); we bypass the $exon->cds() here (hmm,
	    # not very clean style, but I don't see where this screws up)
	    $cds .= $exon->seq()->seq();
	} else {
	    # this one is probably from exon shuffling and needs some work
	    my $seq = $exon->cds(); # now $seq is guaranteed to be in frame 0
	    next if(! $seq);
	    $seq = $seq->seq();
	    # adjustment needed?
	    if($phase > 0) {
		# how many Ns can we chop off the piece to be added?
		my $n_crop = 0;
		if($seq =~ /^(n+)/i) {
		    $n_crop = length($1);
		}
		if($n_crop >= $phase) {
		    # chop off to match the phase
		    $seq = substr($seq, $phase);
		} else {
		    # fill in Ns
		    $seq = ("n" x (3-$phase)) . $seq;
		}
	    }
	    $cds .= $seq;
	}
    }
    return $cds;
}

=head2 features

 Title   : features
 Usage   : my @features=$transcript->features;
 Function: returns all the features associated with this transcript
 Returns : a list of SeqFeatureI implementing objects
 Args    : none


=cut


sub features {
    my $self = shift;    
    return grep { defined } @{$self->{'_features'} || []};
}

=head2 features_ordered

 Title   : features_ordered
 Usage   : my @features=$transcript->features_ordered;
 Function: returns all the features associated with this transcript,
           in order by feature start, according to strand
 Returns : a list of SeqFeatureI implementing objects
 Args    : none


=cut

sub features_ordered{
   my ($self) = @_;
   return $self->_stranded_sort(@{$self->{'_features'} || []});
}


sub get_unordered_feature_type{
    my ($self, $type, $pri)=@_;
    my @list;
    foreach ( $self->features) {
	if ($_->isa($type)) {
	    if ($pri && $_->primary_tag !~ /$pri/i) {
		next;
	    }
	    push @list,$_;
	}
    }
    return @list;

}

sub get_feature_type {
    my ($self)=shift;
    return $self->_stranded_sort($self->get_unordered_feature_type(@_));
}

#This was fixed by Gene Cutler - the indexing on the list being reversed
#fixed a bad bug.  Thanks Gene!
sub _flush {
     my ($self, $type, $pri)=@_;
     my @list=$self->features;
     my @cut;
     for (reverse (0..$#list)) {
         if (defined $list[$_] &&
	     $list[$_]->isa($type)) {
             if ($pri && $list[$_]->primary_tag !~ /$pri/i) {
                 next;
             }
             push @cut, splice @list, $_, 1;  #remove the element of $type from @list
                                              #and return each of them in @cut
         }
     }
     $self->{'_features'}=\@list;
     return reverse @cut;
}

sub _add {
    my ($self, $fea, $type, $pri)=@_;
    require Bio::SeqFeature::Gene::Promoter;
    require Bio::SeqFeature::Gene::UTR;
    require Bio::SeqFeature::Gene::Exon;
    require Bio::SeqFeature::Gene::Intron;
    require Bio::SeqFeature::Gene::Poly_A_site;

    if(! $fea->isa('Bio::SeqFeatureI') ) {
	$self->throw("$fea does not implement Bio::SeqFeatureI");
    }
    if(! $fea->isa($type) || $pri) {
	$fea=$self->_new_of_type($fea,$type,$pri);
    }
    if (! $self->strand) {
	$self->strand($fea->strand);
    } else {
	if ($self->strand * $fea->strand == -1) {
	    $self->throw("$fea is on opposite strand from $self");
	}
    }

    $self->_expand_region($fea);
    if(defined($self->entire_seq()) && (! defined($fea->entire_seq())) &&
       $fea->can('attach_seq')) {
	$fea->attach_seq($self->entire_seq());
    }
    if (defined $self->parent) {
	$self->parent->_expand_region($fea);
    }
    push(@{$self->{'_features'}}, $fea);
    1;
}

sub _stranded_sort {
    my ($self,@list)=@_;
    my $strand;
    foreach my $fea (@list) {
	if($fea->strand()) {
	    # defined and != 0
	    $strand = $fea->strand() if(! $strand);
	    if(($fea->strand() * $strand) < 0) {
		$strand = undef;
		last;
	    }
	}
    }
    if (defined $strand && $strand == - 1) {  #reverse strand
	return map { $_->[0] } sort {$b->[1] <=> $a->[1]} map { [$_, $_->start] } @list;
    } else {               #undef or forward strand
	return map { $_->[0] } sort {$a->[1] <=> $b->[1]} map { [$_, $_->start] } @list;
    }
}

sub _new_of_type {
    my ($self, $fea, $type, $pri)= @_;
    my $primary;
    if ($pri) {
	$primary = $pri;    #can set new primary tag if desired
    } else {
	($primary) = $type =~ /.*::(.+)/;  #or else primary is just end of type string
    }
    bless $fea,$type;
    $fea->primary_tag($primary);
    return $fea;
}

sub transcript_destroy {
    my $self = shift;    
    # We're going to be really explicit to insure memory leaks 
    # don't occur
    foreach my $f ( $self->features ) {
	$f = undef;
    }
    $self->parent(undef);
}

1;
