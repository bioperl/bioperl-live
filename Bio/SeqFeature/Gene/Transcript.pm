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
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Hilmar Lapp

Email hlapp@gmx.net or hilmar.lapp@pharma.novartis.com

Describe contact details here

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

@ISA = qw(Bio::SeqFeature::Generic Bio::SeqFeature::Gene::TranscriptI);

# Default sort order for exon types propagated to a new instance.
# This is only of relevance to exons(), and therefore mrna(), cds(), and
# protein().
my %exon_sortorder = ("utr5prime" => 0,
		      "initial" => 1,
		      "internal" => 2,
		      "terminal" => 3,
		      "utr3prime" => 4,
		      "polyA" => 5);

sub new {
    my ($caller, @args) = @_;
    my $self = $caller->SUPER::new(@args);
    my ($primary) = $self->_rearrange([qw(PRIMARY)],@args);

    $primary = 'transcript' unless $primary;
    $self->primary_tag($primary);
    $self->strand(0) if(! defined($self->strand()));
    $self->{'_exonsortorder'} = \%exon_sortorder;
    return $self;
}


=head2 promoters

 Title   : promoters()
 Usage   : @proms = $transcript->promoters();
 Function: Get the promoter features/sites of this transcript. 

           Note that mrna() expects that there is only a single promoter
           part of the sequence.

 Returns : An array of Bio::SeqFeatureI implementing objects representing the
           promoter regions or sites.
 Args    : 


=cut

sub promoters {
    my ($self) = @_;

    return () unless exists($self->{'_promoters'});
    return @{$self->{'_promoters'}};
}

=head2 add_promoter

 Title   : add_promoter()
 Usage   : $transcript->add_promoter($feature);
 Function: Add a promoter feature/site to this transcript.

           Note that mrna() expects that there is only a single promotor
           part of the sequence.

 Returns : 
 Args    : A Bio::SeqFeatureI implementing object.


=cut

sub add_promoter {
    my ($self, $fea) = @_;

    if(! $fea->isa('Bio::SeqFeatureI') ) {
	$self->throw("$fea does not implement Bio::SeqFeatureI");
    }
    if(! exists($self->{'_promoters'})) {
	$self->{'_promoters'} = [];
    }
    $self->_expand_region($fea);
    if(defined($self->entire_seq()) && (! defined($fea->entire_seq())) &&
       $fea->can('attach_seq')) {
	$fea->attach_seq($fea->entire_seq());
    }
    push(@{$self->{'_promoters'}}, $fea);
}

=head2 flush_promoters

 Title   : flush_promoters()
 Usage   : $transcript->flush_promoters();
 Function: Remove all promoter features/sites from this transcript.
 Returns : 
 Args    :


=cut

sub flush_promoters {
    my ($self) = @_;

    if(exists($self->{'_promoters'})) {
	delete($self->{'_promoters'});
    }
}

=head2 exons

 Title   : exons()
 Usage   : @exons = $gene->exons();
           ($inital_exon) = $gene->exons('Initial');
 Function: Get all exon features or all exons of specified type of this gene
           structure.

           Exon type is treated as a case-insensitive regular expression and 
           optional. For consistency, use only the following types: 
           initial, internal, terminal, utr, utr5prime, and utr3prime. 
           A special and virtual type is 'coding', which refers to all types
           except utr.

 Returns : An array of Bio::SeqFeature::Gene::ExonI implementing objects.
 Args    : An optional string specifying the primary_tag of the feature.


=cut

sub exons {
    my ($self, $type) = @_;
    my @keys = ();
    my @exons = ();

    # pull out all exon types that exist and match
    @keys = _get_typed_keys("exons_",
			    ($type && (lc($type) eq 'coding') ? "" : $type));
    # if none matched we're done
    return () unless(@keys);
    # bring keys into the right order
    my %order = %{$self->exon_type_sortorder()};
    @keys = sort {$order{substr($a,rindex($a,'_')+1)} <=>
		      $order{substr($b,rindex($b,'_'))};} @keys;
    # gather the individual arrays and flatten out into one
    foreach my $key (@keys) {
	if($type && (lc($type) eq 'coding')) {
	    push(@exons, grep { $_->is_coding(); } @{$self->{$key}});
	} else {
	    push(@exons, @{$self->{$key}});
	}
    }
    return @exons;
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
           may deliver unexpected results.

 Returns : 
 Args    : A Bio::SeqFeature::Gene::ExonI implementing object.
           A string indicating the type of the exon (optional).


=cut

sub add_exon {
    my ($self, $fea, $type) = @_;
    my $key;

    if(! $fea->isa('Bio::SeqFeature::Gene::ExonI') ) {
	$self->throw("$fea does not implement Bio::SeqFeature::Gene::ExonI");
    }
    $type = ($type ? lc($type) : "");
    # treat utr separately
    return $self->add_utr($fea,$type) if($type =~ /utr/);
    # prefix key
    $key = "_exons_$type";
    if(! exists($self->{$key})) {
	$self->{$key} = [];
    }
    $self->_expand_region($fea);
    if(defined($self->entire_seq()) &&
       (! defined($fea->entire_seq())) && $fea->can('attach_seq')) {
	$fea->attach_seq($fea->entire_seq());
    }
    push(@{$self->{$key}}, $fea);
}

=head2 flush_exons

 Title   : flush_exons()
 Usage   : $transcript->flush_exons();
           $transcript->flush_exons('terminal');
 Function: Remove all or a certain type of exon features from this transcript.

           See add_exon() for documentation about types.
 
           Calling without a type will not flush UTRs. Call flush_utrs() for
           this purpose.
 Returns : 
 Args    : A string indicating the type of the exon (optional).


=cut

sub flush_exons {
    my ($self, $type) = @_;

    # pull out all exon types that exist and match
    my @keys = grep { $_ !~ /utr/i; } _get_typed_keys("exons_", $type);
    # delete the keys pulled out
    foreach my $key (@keys) {
	delete($self->{$key});
    }
}

=head2 exon_type_sortorder

 Title   : exon_type_sortorder
 Usage   : $transcript->exon_type_sortorder('type1','type2','type3');
           $tableref = $transcript->exon_type_sortorder();
           Bio::SeqFeature::Gene::Transcript->exon_type_sortorder('type1','type2','type3');
 Function: Install or retrieve the sort-order of exon types.
 
           Ordering exons by type is only of relevance for methods
           constructing a sequence corresponding to the transcript object
           (e.g., mrna(), cds(), etc). For instance, an initial exon would
           have to be put before an internal exon.

           The sort-order installed by default covers only the exon types
           documented in exons(). If you use a different set or additional
           types, you must call this method with the proper ordering before
           you can obtain meaningful results from methods that concatenate
           the sequence of exons. On installing a sort-order, pass all types
           you are using in ascending order.

           You can also retrieve the installed order. Do not modify the
           returned hash table (the method returns a reference) unless you know
           exactly what you are doing. Inspecting the keys of the table will
           tell you which types are currently known.

           This method can also be called as a class method. When called as
           class method it will operate on the default sort-order, which will
           be propagated to every object instantiated thereafter. When called
           as instance method it will operate only on the private copy of the
           object.

 Returns : A reference to a hash table representing the sort-order.
 Args    : On installing an order, an array of exon types in ascending order.


=cut

sub exon_type_sortorder {
    my ($caller, @order) = @_;
    my %sorttable;

    if(@order) {
	my $num = 0;
	foreach my $t (@order) {
	    $sorttable{$t} = $num;
	    $num++;
	}
	if(ref($caller)) {
	    $caller->{'_exonsortorder'} = \%sorttable;
	} else {
	    %exon_sortorder = %sorttable;
	}
    }
    return $caller->{'_exonsortorder'} if(ref($caller));
    return \%exon_sortorder;
}

=head2 introns

 Title   : introns()
 Usage   : @introns = $gene->introns();
 Function: Get all intron features this gene structure.

           Note that this implementation generates these features on-the-fly,
           that is, it simply treats all regions between exons as introns,
           assuming that exons do not overlap. A consequence is that a
           consistent correspondence between the elements in the returned
           array and the array exons() returns will exist only if the
           exons are properly sorted within their types (forward for plus-
           strand and reverse for minus-strand transcripts). To ensure
           correctness the elements in the array returned will always be
           sorted.

 Returns : An array of Bio::SeqFeatureI implementing objects representing the
           intron regions.
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
        @exons = sort { $a->start() <=> $b->start(); } @exons;
    } else {
	# sort in reverse order for transcripts on the negative strand and
	# found to be in reverse order
        @exons = sort { $b->start() <=> $a->start(); } @exons;	
    }
    # loop over all intervening gaps
    for(my $i = 0; $i < $#exons; $i++) {
	my ($start, $end);
	my $intron;

	if(defined($exons[$i]->strand()) &&
	   (($exons[$i]->strand() * $strand) < 0)) {
	    $self->throw("Transcript mixes plus and minus strand. ".
			 "This makes no sense.");
	}
	$start = $exons[$i+$rev_order]->end() + 1;     # $i or $i+1
	$end = $exons[$i+1-$rev_order]->start() - 1;   # $i+1 or $i
	$intron = Bio::SeqFeature::Generic->new(
                                        '-start'   => $start,
                                        '-end'     => $end,
                                        '-strand'  => $strand,
                                        '-primary' => 'intron',
					'-source'  => ref($self));
	my $seq = $self->entire_seq();
	$intron->attach_seq($seq) if $seq;
	$intron->seqname($self->seqname());
	push(@introns, $intron);
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

    if(defined($fea)) {
	if(! $fea) {
	    delete($self->{'_poly_A_site'});
	} else {
	    if(! $fea->isa('Bio::SeqFeatureI') ) {
		$self->throw("$fea does not implement Bio::SeqFeatureI");
	    }
	    $self->_expand_region($fea);
	    if(defined($self->entire_seq()) &&
	       (! defined($fea->entire_seq())) && $fea->can('attach_seq')) {
		$fea->attach_seq($fea->entire_seq());
	    }
	    $self->{'_poly_A_site'} = $fea;
	}
    }
    return $self->{'_poly_A_site'};
}

=head2 utrs

 Title   : utrs()
 Usage   : @utr_sites = $transcript->utrs('3prime');
           @utr_sites = $transcript->utrs('5prime');
           @utr_sites = $transcript->utrs();
 Function: Get the features representing untranslated regions (UTR) of this
           transcript.

           You may provide an argument specifying the type of UTR. Currently
           the following types are recognized: 5prime 3prime for UTR on the
           5' and 3' end of the CDS, respectively.

 Returns : An array of Bio::SeqFeatureI implementing objects representing the
           UTR regions or sites.
 Args    : Optionally, either 3prime, or 5prime for the the type of UTR
           feature.


=cut

sub utrs {
    my ($self, $type) = @_;
    my @utrs = ();
    my @keys;

    # pull out all exon types that exist and match
    @keys = _get_typed_keys("exons_utr", $type);
    # gather the individual arrays and flatten out into one
    foreach my $key (@keys) {
	push(@utrs, @{$self->{$key}});
    }
    return @utrs;
}

=head2 add_utr

 Title   : add_utr()
 Usage   : $transcript->add_utr($feature, '3prime');
           $transcript->add_utr($feature);
 Function: Add a UTR feature/site to this transcript.

           The second parameter is optional and denotes the type of the UTR
           feature. Presently recognized types include '5prime' and '3prime'
           for UTR on the 5' and 3' end of a gene, respectively.

           The feature object must at least implement Bio::SeqFeatureI, but
           for increased consistency you should try to supply an object
           that implements Bio::SeqFeature::Gene::ExonI. In the future, this
           may become a requirement.
 Returns : 
 Args    : A Bio::SeqFeatureI or, better yet, a Bio::SeqFeature::Gene::ExonI
           implementing object.


=cut

sub add_utr {
    my ($self, $fea, $type) = @_;

    if(! $fea->isa('Bio::SeqFeatureI') ) {
	$self->throw("$fea does not implement Bio::SeqFeatureI");
    }
    # prefix key
    $type = ($type ? lc($type) : "");
    $type = "utr".$type if($type !~ /^utr/);
    my $key = "_exons_utr$type";
    if(! exists($self->{$key})) {
	$self->{$key} = [];
    }
    $self->_expand_region($fea);
    if(defined($self->entire_seq()) &&
       (! defined($fea->entire_seq())) && $fea->can('attach_seq')) {
	$fea->attach_seq($fea->entire_seq());
    }
    $self->_expand_region($fea);
    push(@{$self->{$key}}, $fea);
}

=head2 flush_utrs

 Title   : flush_utrs()
 Usage   : $transcript->flush_utrs();
           $transcript->flush_utrs('3prime');
 Function: Remove all or a specific type of UTR features/sites from this
           transcript.

           Cf. add_utr() for documentation about recognized types.
 Returns : 
 Args    : Optionally a string denoting the type of UTR feature.


=cut

sub flush_utrs {
    my ($self, $type) = @_;

    # prefix key
    $type = ($type ? lc($type) : "");
    $type = "utr".$type if($type !~ /^utr/);
    # pull out all types that exist and match
    my @keys = _get_typed_keys("exons_", $type);
    # delete the keys pulled out
    foreach my $key (@keys) {
	delete($self->{$key});
    }
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
   push(@feas, $self->poly_A_site());
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
           of all added subfeatures (exons etc.).


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
 Function: Returns the CDS (coding sequence) as defined by the exons of this
           transcript and the attached sequence.

           If no sequence is attached this method will return undef.

           Note that the implementation provided here returns a concatenation
           of all coding exons, thereby assuming that exons do not overlap. 

           Note also that you cannot set the CDS via this method. Set a single
           CDS feature as a single exon, or derive your own class if you want
           to store a predicted CDS.

 Example :
 Returns : A Bio::PrimarySeqI implementing object.
 Args    : 

=cut

sub cds {
    my ($self) = @_;
    my @exons = $self->exons('coding');
    my $strand;

    return undef unless(@exons);
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
    # Make sure exons are sorted if we know the strand. If we don't know the
    # the strand we go with the order found.
    if($strand) {
	if($strand == 1) {
	    # always sort forward for plus-strand transcripts
	    @exons = sort { $a->start() <=> $b->start(); } @exons;
	} else {
	    # sort in reverse order for transcripts on the negative strand
	    @exons = sort { $b->start() <=> $a->start(); } @exons;
	}
    }
    my $cds = $self->_make_cds(@exons);
    return undef unless $cds;
    return Bio::PrimarySeq->new('-id' => $self->seqname(),
				'-seq' => $cds);
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
    return undef;
}

=head2 mrna

 Title   : mrna()
 Usage   : $mrna = $transcript->mrna();
 Function: Get the mRNA of the transcript as a sequence object.

           The difference to cds() is that the sequence object returned by
           this methods will also include UTR and the poly-adenylation site,
           but not promoter sequence (TBD);

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
	$seq = Bio::PrimarySeq->new('-id' => $self->seqname(),
				    '-moltype' => "rna",
				    '-seq' => "");
    }
    # get and add UTR sequences
    $mrna = "";
    foreach $elem ($self->utrs('5prime')) {
	$mrna .= $elem->seq()->seq();
    }
    $seq->seq($mrna . $seq->seq());
    $mrna = "";
    foreach $elem ($self->utrs('3prime')) {
	$mrna .= $elem->seq()->seq();
    }
    $seq->seq($seq->seq() . $mrna);
    if($self->poly_A_site()) {
	$seq->seq($seq->() . $self->poly_A_site()->seq()->seq());
    }
    return undef if($seq->length() == 0);
    return $seq;
}

sub _get_typed_keys {
    my ($self, $keyprefix, $type) = @_;
    my @keys = ();
    my @feas = ();

    # make case-insensitive
    $type = ($type ? lc($type) : "");
    # pull out all feature types that exist and match
    @keys = grep { $_ =~ /^_$keyprefix$type/i; } (keys(%{$self}));
    return @keys;
}

sub _make_cds {
    my ($self,@exons) = @_;
    my $cds = "";

    foreach my $exon (@exons) {
	next if((! defined($exon->seq())) || (! $exon->is_coding()));
	my $phase = length($cds) % 3;
	# let's check the simple case 
	if(defined($exon->frame()) && ($phase == $exon->frame())) {
	    # this one fits exactly
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

1;
