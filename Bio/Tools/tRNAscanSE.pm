#
# BioPerl module for Bio::Tools::tRNAscanSE
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::tRNAscanSE - A parser for tRNAscan-SE output

=head1 SYNOPSIS

   use Bio::Tools::tRNAscanSE;

   my $parser = Bio::Tools::tRNAscanSE->new(-file => 'result.tRNAscanSE');

   # parse the results
   while( my $gene = $parser->next_prediction ) {

       @exon_arr = $gene->get_SeqFeatures();

   }

=head1 DESCRIPTION

This script will parse tRNAscan-SE output.  Just the tabular output of
the tRNA locations in the genome for now.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::tRNAscanSE;
use strict;

use Bio::SeqFeature::Generic;

use base qw(Bio::Tools::AnalysisResult);

use vars qw($GeneTag $SrcTag $ExonTag);
($GeneTag,$SrcTag,$ExonTag) = qw(gene tRNAscan-SE exon);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::tRNAscanSE->new();
 Function: Builds a new Bio::Tools::tRNAscanSE object 
 Returns : an instance of Bio::Tools::tRNAscanSE
 Args    : -fh/-file for input filename
           -genetag => primary tag used in gene features (default 'tRNA_gene')
           -exontag => primary tag used in exon features (default 'tRNA_exon')
           -srctag  => source tag used in all features (default 'tRNAscan-SE')


=cut

sub _initialize {
  my($self,@args) = @_;
  $self->SUPER::_initialize(@args);
  my ($genetag,$exontag,$srctag) =  $self->SUPER::_rearrange([qw(GENETAG
								 SRCTAG
								 EXONTAG)],
							      @args);
  $self->gene_tag(defined $genetag ? $genetag : $GeneTag);
  $self->source_tag(defined $srctag ? $srctag : $SrcTag);
  $self->exon_tag(defined $exontag ? $exontag : $ExonTag);
  $self->{'_seen'} = {};
}

=head2 gene_tag

 Title   : gene_tag
 Usage   : $obj->gene_tag($newval)
 Function: Get/Set the value used for the 'gene_tag' of genes
           Default is 'tRNA_gene' as set by the global $GeneTag
 Returns : value of gene_tag (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub gene_tag{
    my $self = shift;

    return $self->{'gene_tag'} = shift if @_;
    return $self->{'gene_tag'};
}

=head2 source_tag

 Title   : source_tag
 Usage   : $obj->source_tag($newval)
 Function: Get/Set the value used for the 'source_tag' of exons and genes
           Default is 'tRNAscan-SE' as set by the global $SrcTag
 Returns : value of source_tag (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub source_tag{
    my $self = shift;

    return $self->{'_source_tag'} = shift if @_;
    return $self->{'_source_tag'};
}

=head2 exon_tag

 Title   : exon_tag
 Usage   : $obj->exon_tag($newval)
 Function: Get/Set the value used for the 'primary_tag' of exons
           Default is 'tRNA_exon' as set by the global $ExonTag
 Returns : value of exon_tag (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub exon_tag{
    my $self = shift;

    return $self->{'_exon_tag'} = shift if @_;
    return $self->{'_exon_tag'};
}


=head2 analysis_method

 Usage     : $genscan->analysis_method();
 Purpose   : Inherited method. Overridden to ensure that the name matches
             /tRNAscan-SE/i.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method { 
#-------------
    my ($self, $method) = @_;  
    if($method && ($method !~ /tRNAscan-SE/i)) {
	$self->throw("method $method not supported in " . ref($self));
    }
    return $self->SUPER::analysis_method($method);
}

=head2 next_feature

 Title   : next_feature
 Usage   : while($gene = $genscan->next_feature()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the Genscan result
           file. Call this method repeatedly until FALSE is returned.

           The returned object is actually a SeqFeatureI implementing object.
           This method is required for classes implementing the
           SeqAnalysisParserI interface, and is merely an alias for 
           next_prediction() at present.

 Example :
 Returns : A Bio::SeqFeature::Generic object.
 Args    :
See also : L<Bio::SeqFeature::Generic>

=cut

sub next_feature {
    my ($self,@args) = @_;
    # even though next_prediction doesn't expect any args (and this method
    # does neither), we pass on args in order to be prepared if this changes
    # ever
    return $self->next_prediction(@args);
}

=head2 next_prediction

 Title   : next_prediction
 Usage   : while($gene = $genscan->next_prediction()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the Genscan result
           file. Call this method repeatedly until FALSE is returned.

 Example :
 Returns : A Bio::SeqFeature::Generic object.
 Args    :
See also : L<Bio::SeqFeature::Generic>

=cut

sub next_prediction {
    my ($self) = @_;
    my ($genetag,$srctag,$exontag) = ( $self->gene_tag,
				       $self->source_tag,
				       $self->exon_tag);

    while( defined($_ = $self->_readline) ) {
	if( m/^(\S+)\s+       # sequence name
	    (\d+)\s+          # tRNA #
	    (\d+)\s+(\d+)\s+  # tRNA start,end
	    (\w{3})\s+        # tRNA type
	    ([CAGT]{3})\s+    # Codon
	    (\d+)\s+(\d+)\s+  # Intron Begin End
	    (\d+\.\d+)/ox     # Cove Score
	    ) {
	    
	    my ($seqid,$tRNAnum,$start,$end,$type,
		$codon,$intron_start,$intron_end,
		$score) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);
	    
	    my $strand = 1;
	    if( $start > $end ) { 
		($start,$end,$strand) = ($end,$start,-1);
	    }
	    if( $self->{'_seen'}->{$type}++ ) {
		$type .= "-".$self->{'_seen'}->{$type};
	    }
	    my $gene = Bio::SeqFeature::Generic->new
		( -seq_id => $seqid,
		  -start  => $start,
		  -end    => $end,
		  -strand => $strand,
		  -score  => $score,
		  -primary_tag => $genetag,
		  -source_tag  => $srctag,
		  -tag     => {
		      'ID'    => "tRNA:$type",
		      'Name'  => "tRNA:$type",
		      'AminoAcid' => $type,
		      'Codon'     => $codon,
		  });
	    if( $intron_start ) {
		if( $intron_start > $intron_end ) {
		    ($intron_start,$intron_end) = ($intron_end,$intron_start);
		}
		$gene->add_SeqFeature(Bio::SeqFeature::Generic->new
				      ( -seq_id=> $seqid,
					-start => $start,
					-end   => $intron_start-1,
					-strand=> $strand,
					-primary_tag => $exontag,
					-source_tag  => $srctag,
					-tag => { 
					    'Parent' => "tRNA:$type",
					    }));
		$gene->add_SeqFeature(Bio::SeqFeature::Generic->new
				      ( -seq_id=> $seqid,
					-start => $intron_end+1,
					-end   => $end,
					-strand=> $strand,
					-primary_tag => $exontag,
					-source_tag  => $srctag,
					-tag => { 
					    'Parent' => "tRNA:$type" 
					    }));
	    } else {
		$gene->add_SeqFeature(Bio::SeqFeature::Generic->new
				      ( -seq_id=> $seqid,
					-start => $start,
					-end   => $end,
					-strand=> $strand,
					-primary_tag => $exontag,
					-source_tag  => $srctag,
					-tag => { 
					     'Parent' => "tRNA:$type" 
					     }));
	    }
	    return $gene;
	} 
    }
}

1;
