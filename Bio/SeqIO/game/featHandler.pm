# $Id: 
# 
#
# Helper module for Bio::SeqIO::game::featHandler
#
# Cared for by Sheldon McKay <smckay@bcgsc.bc.ca>
#
# Copyright Sheldon McKay
#
# You may distribute this module under the same terms as perl itself
#

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::game::featHandler -- a class for handling feature elements

=head1 SYNOPSIS

This module is not used directly

=head1 DESCRIPTION

Bio::SeqIO::game::featHandler converts game XML <annotation> elements into
flattened Bio::SeqFeature::Generic objects to be added to the sequence

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
to one of the Bioperl mailing lists.

Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/MailList.shtml      - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution.

_Bug reports can be submitted via email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Sheldon McKay

Email smckay@bcgsc.bc.ca

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::SeqIO::game::featHandler;

use Bio::SeqIO::game::gameSubs;
use Bio::SeqFeature::Generic;
use Bio::Location::Split;
use strict;

use vars qw { @ISA };                                                                                

@ISA = qw { Bio::SeqIO::game::gameSubs };

=head2 new

 Title   : new
 Usage   : my $featHandler = Bio::SeqIO::game::featHandler->new($seq, $seq_h, $ann_l)
 Function: creates an object to deal with sequence features 
 Returns : a handler object
 Args    : $seq   -- a Bio::SeqI compliant object
           $seq_h -- ref. to a hash of other sequences associated 
                     with the main sequence (proteins, etc)
           $ann_l -- ref. to a list of annotations

=cut

sub new {
    my ($caller, $seq, $seq_h, $ann_l  ) = @_;
    my $class = ref($caller) || $caller;

    my $self = bless ({                                                                             
        seq           => $seq,                                                                           
        curr_feats    => [],
	curr_coords   => [],
	seq_h         => $seq_h,
	ann_l         => $ann_l,
    }, $class);

    return $self;
}

=head2 add_source

 Title   : add_source
 Usage   : $featHandler->add_source($seq->length, \%tags);
 Function: creates a source feature
 Returns : a Bio::SeqFeature::Generic object 
 Args    : sequence length and a ref. to a hash of tag/value attributes

=cut

sub add_source {
    my ($self, $length, $tags) = @_;
    my $feat = Bio::SeqFeature::Generic->new( -primary => 'source',
					      -start   => 1,
					      -end     => $length,
					    );
    for ( keys %{$tags} ) {
	for my $val ( @{$tags->{$_}} ) {
	    $feat->add_tag_value( $_ => $val );
	}
    }

    return $feat;
}

=head2 has_gene

 Title   : has_gene
 Usage   : my $gene = $self->_has_gene($gene, $gname, $id)
 Function: method to get/set the current gene feature
 Returns : a Bio::SeqFeature::Generic object (if there is a gene)
 Args    : (optional)
           $gene  -- an XML element for the annotation
           $gname -- gene name
           $id    -- gene ID (not always the same as the name)

=cut

sub has_gene {
    my ($self, $gene, $gname, $id) = @_;
    
    unless ( $gene ) {
	if ( defined $self->{curr_gene} ) {
	    return $self->{curr_gene};
	}
        else {
	    return 0;
        }
    }
    else {
        if ( $id && !$self->{curr_ltag} ) {
	    $self->{curr_ltag} = $id;
	}
	if ( $gname && $gname ne $id && !$self->{curr_gname} ) {
	    $self->{curr_gname} = $gname;
	}
	    
	my $tags  = {};
	
	for my $child ( @{$gene->{Children}} ) {
	    my $name = $child->{Name};

	    if ( $name eq 'dbxref' ) {
	        $tags->{dbxref} ||= [];
		push @{$tags->{dbxref}}, $self->dbxref( $child );
	    }
            elsif ( $name !~ /name/ ){
                $self->complain("Unrecognized element '$name'. I don't " .
			        "know what to do with $name elements");
	    }
	}
	
	my $feat = Bio::SeqFeature::Generic->new( 
	    -primary => 'gene',
	);

	for ( keys %{$tags} ) {
	    for my $val ( @{$tags->{$_}} ) {
		$feat->add_tag_value( $_ => $val );
	    }
	}
	    
	$self->{curr_gene} = $feat;
	return $feat;
    }	
}

=head2 _has_CDS

 Title   : _has_CDS
 Usage   : my $cds = $self->_has_CDS
 Function: internal getter/setter for CDS features
 Returns : a Bio::SeqFeature::Generic transcript object (or nothing)
 Args    : a Bio::SeqFeature::Generic transcript feature

=cut

sub _has_CDS {
    my ($self, $transcript) = @_;

    if ( !$transcript ) {
	if ( defined $self->{curr_cds} ) {
	    return $self->{curr_cds};
        }
	else {
	    return 0;
	}
    }
    else {
	my $tags = $self->{curr_tags};

	unless ( defined $tags->{product} ){
	    return 0;
	}
	
	$self->{curr_cds} = $self->_add_CDS( $transcript, $tags );
    }
}

=head2 add_annotation

 Title   : add_annotation
 Usage   : $featHandler->add_annotation($seq, $type, $id, $tags, $feats)
 Function: converts a containment hierarchy into an ordered list of flat features
 Returns : nothing
 Args    : $seq   -- a Bio::SeqI compliant object
           $type  -- the annotation type
           $id    -- the anotation ID
           $tags  -- ref. to a hash of tag/value attributes
           $feats -- ref to an array of Bio::SeqFeature::Generic objects

=cut

sub add_annotation {
    my ($self, $seq, $type, $id, $tags, $feats) = @_;

    # is this a generic feature?
    unless ( $self->has_gene ) {
	shift;
	$self->_add_generic_annotation(@_);
	return 0;
    }

    my $feat;

    if ( $type eq 'gene' ) {
	$feat = $self->has_gene;
	$feat->add_tag_value( standard_name => $id );
	$feat->add_tag_value( gene => ($self->{curr_gname} || $id) );
        
    }
    else {
	$feat = Bio::SeqFeature::Generic->new;
	$feat->primary_tag($type);
	my $gene = $self->has_gene;
	$gene->add_tag_value( standard_name => $id );
	$gene->add_tag_value( gene => ($self->{curr_gname} || $id) );
    }
    for ( keys %{$tags} ) {
	# or else add simple tag/value pairs
	if ( $_ eq 'name' && $tags->{type}->[0] eq 'gene' ) {
	    unless ( $feat->has_tag( 'gene' ) ) {
		$feat->add_tag_value( gene => $tags->{name}->[0] );
	    }
	    delete $tags->{name};
	}
	else {
	    next if $_ eq 'type' && $tags->{$_}->[0] eq 'gene';
	    next if $_ eq 'gene' && $feat->has_tag( 'gene' );
	    for my $val ( @{$tags->{$_}} ) {
		$feat->add_tag_value( $_ => $val );
	    }
	}
    }


    $feat->strand( $self->{curr_strand} );
    $feat->start( $self->{curr_coords}->[0] );
    $feat->end( $self->{curr_coords}->[1] );

    # create an array of features for the annotation (order matters)
    my @annotations = ( $feat );

    # add the gene feature if the annotation is not a gene
    if ( $self->has_gene && $type ne 'gene') {
	my $gene = $self->has_gene;
	$gene->strand( $self->{curr_strand} );
	$gene->start( $self->{curr_coords}->[0] );
        $gene->end( $self->{curr_coords}->[-1] );
	push @annotations, $gene;
	$self->{curr_gene} = {};
    }

    # add the subfeatures
    for ( @{$feats} ) {
	$self->complain("bad feature $_") unless ref($_) =~ /Bio/;
	push @annotations, $_;
    }
    
    # add the annotation array to the list for this sequence
    my $seqid = $seq->id;
    my $list = $self->{ann_l};
    
    # make sure the feature_sets appear in ascending order
    if ( $list->[0] && $annotations[0]->start < $list->[0]->start ) {
	    unshift @{$list}, @annotations;
       }
    else {
        push @{$list}, @annotations;
    }

    # garbage collection
    $self->{curr_ltag}   = '';
    $self->{curr_gname}  = '';
    $self->{curr_coords} = [];
    $self->{curr_feats}  = [];
    $self->{curr_strand} = 0;
    $self->{ann_seq}     = $seq;    
    $self->flush;
}


=head2 _add_generic_annotation

 Title   : _add_generic_annotation
 Usage   : $self->_add_generic_annotation($seq, $type, $id, $tags, $feats)
 Function: an internal method to handle non-gene annotations
 Returns : nothing
 Args    : $seq   -- a Bio::SeqI compliant object
           $type  -- the annotation type
           $id    -- the anotation ID
           $tags  -- ref. to a hash of tag/value attributes
           $feats -- ref to an array of Bio::SeqFeature::Generic objects

=cut

sub _add_generic_annotation {
    my ($self, $seq, $type, $id, $tags, $feats) = @_;

    push @{$self->{ann_l}}, @$feats;

    $self->{curr_coords} = [];
    $self->{curr_feats}  = [];
    $self->{curr_strand} = 0;
    $self->{ann_seq}     = $seq;
    $self->flush;
}


=head2 feature_set

 Title   : feature_set
 Usage   : push @feats, $featHandler->feature_set($id, $gname, $set, $anntype);
 Function: handles <feature_span> hierarchies (usually a transcript)
 Returns : a list of Bio::SeqFeature::Generic objects
 Args    : $id      -- ID of the feature set
           $gname   -- name of the gene
           $set     -- the <feature_set> object
           $anntype -- type of the parent annotation


=cut


sub feature_set {
    my ($self, $id, $gname, $set, $anntype) = @_;
    my $stype = $set->{_type}->{Characters};
    $self->{curr_loc}      = [];
    $self->{curr_tags}     = {};
    $self->{curr_subfeats} = [];
    $self->{curr_strand}   = 0;
    my @feats = ();
    my $generic;
    my $tags = $self->{curr_tags};
    my $sname = $set->{_name}->{Characters} ||
        $set->{Attributes}->{id};

    if ( $set->{Attributes}->{problem} ) {
        $tags->{problem} = [$set->{Attributes}->{problem}];
    }

    my @fcount = grep { $_->{Name} eq 'feature_span' } @{$set->{Children}};
    
    if ( @fcount == 1 ) {
	$self->_build_feature_set($set, 1);
	my ($feat) = @{$self->{curr_subfeats}};
	$feat->primary_tag('transcript') if $feat->primary_tag eq 'exon';

	for my $tag ( keys %{$tags} ) {
	    for my $val ( @{$tags->{$tag}} ) {
		$feat->add_tag_value( $tag => $val ) if $val;
	    }
	}
	$feat->add_tag_value( gene => ($gname || $id) );
	@feats = ($feat);
    }
    else {
	$self->{curr_ltag}     = $id;
	$self->{curr_cds}      = '';
	$gname = $id if $gname eq 'gene';
	$self->{curr_gname} = $gname;
    
	if ( $self->has_gene ) {
	    unless ( $anntype =~/RNA/i ) {
		$stype =~ s/transcript/mRNA/;
	    }
	}
	
	$self->{curr_feat}  = Bio::SeqFeature::Generic->new(
							    -primary => $stype,
							    -id      => $id,
							    );
	my $feat = $self->{curr_feat};
	$self->_build_feature_set($set);
	
	if ( $self->{curr_ltag} || $gname ) {
	    $self->{curr_feat}->add_tag_value( gene => ($gname || $self->{curr_ltag}) );
	}
	
	# if there is an annotated protein product
	my $cds = $self->_has_CDS( $feat );

	my $gene = $gname || $self->{curr_ltag};
	$feat->add_tag_value( gene => $gene );
	
	if ( $cds ) {
	    my $gene = $gname || $self->{curr_ltag};
	    # we really just want one value here
	    $cds->remove_tag('standard_name') if $cds->has_tag('standard_name');
	    $cds->add_tag_value( standard_name => $sname );
	    $cds->remove_tag('gene') if $cds->has_tag('gene');
	    $cds->add_tag_value( gene => $gene );
	    
	    # catch missing/empty protein ids
            if ( $cds->has_tag('protein_id' )) {
		if ( !$cds->get_tag_values('protein_id') ) {
		    $cds->remove_tag('protein_id');
		    if ( $cds->has_tag('product') ) {
			$cds->add_tag_value($cds->get_tag_values('product'));
		    }
		}
	    }

	    # make sure other subfeats are tied to the transcript
            # via a 'standard_name' qualifier and the gene via a 'gene' qualifier
	    my @subfeats = @{$self->{curr_subfeats}};
            for my $sf ( @ subfeats ) {
                $sf->add_tag_value( standard_name => $sname );
                $sf->add_tag_value( gene => $gene );
            }
	    
	    $feat->add_tag_value( standard_name => $sname );

            @feats = sort { $a->start <=> $b->start } ($cds, @subfeats);
	    unshift @feats, $feat;
	}
	else {
	    if ( @{$self->{curr_loc}} > 1 ) {
		my $loc = Bio::Location::Split->new( -splittype => 'JOIN' );
		
		# sort the exons in ascending start order
		my @loc = sort { $a->start <=> $b->start } @{$self->{curr_loc}};
		
		# then add them to the transcript location
		for ( @loc ) {
		    $loc->add_sub_Location( $_ ) 
		}
		$feat->location( $loc );
	    }
	    else {
		$feat->location( $self->{curr_loc}->[0] );
	    }	
	    
	    # make sure other subfeats are tied to the transcript
	    my @subfeats = @{$self->{curr_subfeats}};
	    for my $sf ( @ subfeats ) {
		$sf->add_tag_value( standard_name => $sname );
		$sf->add_tag_value( gene => $gene );
	    }

	    @feats = ( $feat, @subfeats );
	}
    }
    
    # adjust the maximum extent of the annotated feature 
    # if req'd (ie the <annotation> element)
    $self->{curr_coords}->[0] ||= 1000000000000;
    $self->{curr_coords}->[1] ||= -1000000000000;
    for ( @feats ) {
	if ( $self->{curr_coords}->[0] > $_->start ) {
	    $self->{curr_coords}->[0] = $_->start;
	}
	if ( $self->{curr_coords}->[1] < $_->end ) {
	    $self->{curr_coords}->[1] = $_->end;
	}
    }
    
    $self->flush( $set );

    return @feats;
}


=head2 _build_feature_set

 Title   : _build_feature_set
 Usage   : $self->_build_feature_set($set, 1) # 1 flag means retain the exon as a subfeat
 Function: an internal method to process attributes and subfeats of a feature set
 Returns : nothing
 Args    : $set -- a <feature_set> element
           1    -- optional flag to retain exons as subfeats.  Otherwise, they will
                   be converted to sublocations of a parent CDS feature

=cut


sub _build_feature_set {
    my ($self, $set, $keep_subfeat) = @_;

    for my $child ( @{$set->{Children}} ) {
        my $name = $child->{Name};

        # these elements require special handling
        if ( $name eq 'date' ) {
            $self->date( $child );
        }
        elsif ( $name eq 'comment' ) {
            $self->comment( $child );
        }
        elsif ( $name eq 'evidence' ) {
            $self->evidence( $child );
        }
        elsif ( $name eq 'feature_span' ) {
            $self->_add_feature_span( $child, $keep_subfeat );
	}
        elsif ( $name eq 'property' ) {
            $self->property( $child );
        }

        # need to add the db_xref tags to the gene??
        # otherwise, simple tag/value pairs
        elsif ( $name =~ /synonym|author|description/) {
            $self->{curr_tags}->{$name} = [$child->{Characters}];
        }
        elsif ( $name !~ /name|type|seq/ ){
            $self->complain("Unrecognized element '$name'. I don't " .
                            "know what to do with $name elements");

        }
    }
}

=head2 _add_feature_span
 
 Title   : _add_feature_span
 Usage   : $self->_add_feature_span($el, 1)
 Function: an internal method to process <feature_span> elements
 Returns : nothing
 Args    : $el -- a <feature_span> element
           1   -- an optional flag to retain exons as subfeatures


=cut


sub _add_feature_span {
    my ($self, $el, $keep_subfeat) = @_;

    my $tags  = $self->{curr_tags};
    my $feat  = $self->{curr_feat};
    my $type  = $el->{_type}->{Characters} || $el->{Name};
    my $id    = $el->{Attributes}->{id} || $el->{_name}->{Characters};
    my $seqr  = $el->{_seq_relationship};
    my $start = int $seqr->{_span}->{_start}->{Characters};
    my $end   = int $seqr->{_span}->{_end}->{Characters};
    my $stype = $seqr->{Attributes}->{type}; 
    my $seqid = $seqr->{Attributes}->{seq};

    push @{$self->{seq_l}}, $self->{seq_h}->{$seqid};

    if ( $start > $end ) {
	$self->{curr_strand} = -1;
	($start, $end) = ($end, $start);
    }
    else {
	$self->{curr_strand} = 1;
    }

    # add exons to the transcript
    if ( $type eq 'exon' ) {
	my $sl = Bio::Location::Simple->new( -start  => $start,
                                             -end    => $end,
                                             -strand => $self->{curr_strand} );
        push @{$self->{curr_loc}}, $sl;
    }
    
    # apollo and gadfly use different tags for the same thing 
    if ( $type =~ /start_codon|translate offset/ ) {
        $self->{curr_tags}->{codon_start} = [$start];
    }
    else { 
	if ( $type eq 'exon' ) {
	    return unless $keep_subfeat;
	}
	push @{$self->{curr_subfeats}}, 
	Bio::SeqFeature::Generic->new( -start   => $start,
				       -end     => $end,
				       -strand  => $self->{curr_strand},
				       -primary => $type );
    }

    # identify the translation product                                             
    if ( $el->{Attributes}->{produces_seq} ) {
	my $subseq = $self->{seq_h}->{$el->{Attributes}->{produces_seq}};                       
        $self->{curr_tags}->{product} = [$el->{Attributes}->{produces_seq}];
	$self->{curr_tags}->{translation} = [$subseq->seq];
    }      

    $self->flush( $el );
}

=head2 _add_CDS

 Title   : _add_CDS
 Usage   : my $cds = $self->_add_CDS($transcript, $tags)
 Function: an internal method to create a CDS feature from a transcript feature
 Returns : a Bio::SeqFeature::Generic object
 Args    : $transcript -- a Bio::SeqFeature::Generic object for a transcript
           $tags       -- ref. to a hash of tag/value attributes
           
=cut

sub _add_CDS {
    my ($self, $feat, $tags) = @_;
    my $loc  = {};
    my $single = 0;

    if ( @{$self->{curr_loc}} > 1 ) {        
        $loc = Bio::Location::Split->new;

        # sort the exons in ascending start order
        my @loc = map  { $_->[1] }
                  sort { $a->[0] <=> $b->[0] }
                  map  { [$_->start, $_] } @{$self->{curr_loc}};

        # then add them to the location object
        for ( @loc ) {
            $loc->add_sub_Location( $_ );
        }
    }
    else {
        $loc = $self->{curr_loc}->[0];
        $single++;
    }

    # create a CDS
    my @exons = $single ? $loc : $loc->sub_Location(1);

    $feat->location($loc);

    # try to find a peptide
    my $seq  = $self->{seq_h}->{ $tags->{product}->[0] } ||
	       $self->{seq_h}->{ $tags->{protein_id}->[0] } ||
	       $self->{seq_h}->{ $tags->{gene}->[0] } ||
	       $self->{seq_h}->{ $tags->{standard_name}->[0] };
    
    unless ( $seq ) {
	$self->warn("I did not find a protein sequence for " . 
		     $feat->display_name . ".  I can't figure out the CDS boundaries");
	$feat->location($loc);
	next;
    }	
    my $peptide = $seq->display_id;

    # Can we count on the description format being consistent?
    # Why is critical CDS info saved as unconstrained text not 
    # specified in the DTD? 
    # Anyone have a better idea?    
    my $desc = $seq->description;
    $desc =~ s/,//g;
    $desc =~ s/\]\)/\]\) /g;
    my ($start, $stop) = ();
    
    if ( $desc =~ /cds_boundaries:.+?(\d+)\.\.(\d+)/ ) {
	($start, $stop) = ($1 - $self->{offset}, $2 - $self->{offset});
    }
    else {
	$self->warn("I could not find the CDS boundaries");
	$feat->location($loc);
	next;
    }
    
    delete $tags->{transcript};
    
    # now chop off the UTRs to create a CDS
    my @exons_to_add = ();
    for ( @exons ) {
        my $exon = Bio::Location::Simple->new;

	if ( $_->end < $start || $_->start > $stop ) {
	    next;
	}
	if ( $_->start < $start && $_->end > $start ) {
	    $exon->start( $start );
	}
	if ( $_->end > $stop && $_->start < $stop ) {
	    $exon->end( $stop );
	}

	$exon->start( $_->start ) unless $exon->start;
	$exon->end( $_->end ) unless $exon->end;
	$exon->strand ( $self->{curr_strand} );
	push @exons_to_add, $exon;
    }

    my $cds_loc;
    if ( @exons_to_add > 1 ) {
        $cds_loc = Bio::Location::Split->new( -splittype => 'JOIN'  );
        for ( @exons_to_add ) {
	    $cds_loc->add_sub_Location( $_ );
        }
    }
    else {
	$cds_loc = $exons_to_add[0];	
    }

    my $parent = $self->{curr_gname} || $self->{curr_ltag};
    my $cds_tags = $tags; 
    $tags = {};
    
    #delete $feat->{_gsf_tag_hash};
    for ( $feat->all_tags ) {
	$feat->remove_tag($_);
    }    

    $cds_tags->{product_desc} = [$desc];

    if ( $self->{curr_gname} ) {
        $cds_tags->{gene} = [$self->{curr_gname}];	
    }
    
    my $gene = $self->has_gene;
    
    my $cds = Bio::SeqFeature::Generic->new( 
        -primary  => 'CDS',
	-location => $cds_loc,
    );

    $cds_tags->{translation} = [$seq->seq];
  
    for ( keys %{$cds_tags} ) {
	for my $val (@{$cds_tags->{$_}}) {
	    $cds->add_tag_value( $_ => $val );
	}        
    }
    
    $cds;

}

1;



