#
# BioPerl module for Bio::SeqFeature::SiRNA::Pair
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Donald Jackson, donald.jackson@bms.com
#
# Copyright Donald Jackson
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::SiRNA::Pair - Perl object for small inhibitory RNA
(SiRNA) oligo pairs

=head1 SYNOPSIS

  use Bio::SeqFeature::SiRNA::Pair;
  my $pair = Bio::SeqFeature::SiRNA::Pair->
      new( -sense       => $bio_seqfeature_sirna_oligo, # strand=1
           -antisense	=> $bio_seqfeature_sirna_oligo, # strand= -1
	   -primary	=> 'SiRNA::Pair',
	   -source_tag 	=> 'Bio::Tools::SiRNA',
	   -start	=> 8,
	   -end		=> 31,
	   -rank	=> 1,
	   -fxgc	=> 0.5,
	   -tag		=> { note => 'a note' } );

  $target_sequence->add_SeqFeature($pair);

=head1 DESCRIPTION

Object methods for (complementary) pairs of L<Bio::SeqFeature::SiRNA::Oligo> 
objects - inherits L<Bio::SeqFeature::Generic>. See that package for information
on inherited methods.

Does B<not> include methods for designing SiRNAs -- see L<Bio::Tools::SiRNA>

=head1 SEE ALSO

L<Bio::SeqFeature::Oligo>, L<Bio::Tools::SiRNA>.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Donald Jackson (donald.jackson@bms.com)

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::SeqFeature::SiRNA::Pair;

use strict;
use warnings;

use base qw(Bio::SeqFeature::Generic);

# arguments to new().  Taken from Bio::SeqFeature Generic.
# Omit frame (not relevant), GFF_STRING and GFF1_STRING 
# because I'm not sure how to handle them.  Add RANK, FXGC, SENSE, ANTISENSE
our @ARGNAMES = qw(RANK FXGC SENSE ANTISENSE START END STRAND PRIMARY SOURCE_TAG
		   SCORE TAG SEQNAME ANNOTATION LOCATION);

=head1 METHODS

=head2 new

  Title		: new
  Usage		: my $sirna_pair = Bio::SeqFeature::SiRNA::Pair->new();
  Purpose	: Create a new SiRNA::Pair object
  Returns	: Bio::Tools::SiRNA object
  Args		: -start 	10
                  -end		31
                  -rank		1 #  'Rank' in Tuschl group's rules
                  -fxgc		0.5 # GC fraction for target sequence
		  -primary	'SiRNA::Pair', # default value
		  -source_tag	'Bio::Tools::SiRNA'
		  -tag		{ note => 'A note' }
                  -sense	a Bio::SeqFeature::SiRNA::Oligo object
                                with strand = 1
                  -antisense	a Bio::SeqFeature::SiRNA::Oligo object
                                with strand = -1
);

  Note		: SiRNA::Pair objects are typically created by a design 
                  algorithm such as Bio::Tools::SiRNA

=cut

sub new {
    my ($proto, @args) = @_;

    my $pkg = ref($proto) || $proto;

    my $self = $pkg->SUPER::new();
    my %args;
    @args{@ARGNAMES} = $self->_rearrange(\@ARGNAMES, @args); 
    # default primary tag
    $args{'PRIMARY'} ||= 'SiRNA::Pair';

    $args{'PRIMARY'}		&& $self->primary_tag($args{'PRIMARY'});
    $args{'SOURCE_TAG'}		&& $self->source_tag($args{'SOURCE_TAG'});
    $args{'SEQNAME'}		&& $self->seqname($args{'SEQNAME'});
    $args{'ANNOTATION'}		&& $self->annotation($args{'ANNOTATION'});
    $args{'LOCATION'}		&& $self->location($args{'LOCATION'});
    $args{'SENSE'}		&& $self->sense($args{'SENSE'});
    $args{'ANTISENSE'}		&& $self->antisense($args{'ANTISENSE'});
    defined($args{'START'})	&& $self->start($args{'START'});
    defined($args{'END'})	&& $self->end($args{'END'});
    defined($args{'STRAND'})	&& $self->strand($args{'STRAND'});
    defined($args{'SCORE'})	&& $self->score($args{'SCORE'});
    defined($args{'RANK'})	&& $self->rank($args{'RANK'});
    defined($args{'FXGC'})	&& $self->fxGC($args{'FXGC'});

    if ($args{'TAG'}) {	
	foreach my $t (keys %{$args{'TAG'}}) {
	    $self->add_tag_value($t, $args{'TAG'}->{$t});
	}
    }


    return $self;
}

=head2 rank

  Title		: rank
  Usage		: my $pair_rank = $sirna_pair->rank()
  Purpose	: Get/set the 'quality rank' for this pair.
                  See Bio::Tools::SiRNA for a description of ranks.
  Returns	: scalar
  Args		: scalar (optional) indicating pair rank

=cut

sub rank {
    my ($self, $rank) = @_;

    if (defined $rank) {
	# first clear out old tags
	$self->remove_tag('rank') if ( $self->has_tag('rank') );
	$self->add_tag_value('rank', $rank);
    }
    else {
	if ($self->has_tag('rank')) {
	    my @ranks = $self->get_tag_values('rank');
	    return shift @ranks;
	}
	else {
	    $self->throw("Rank not defined for this Pair\n");
	    return;
	}
    }
}

=head2 fxGC

  Title		: fxGC
  Usage		: my $fxGC = $sirna_pair->fxGC();
  Purpose 	: Get/set the fraction of GC for this pair - based on TARGET sequence, not oligos.
  Returns 	: scalar between 0-1
  Args		: scalar between 0-1 (optional)

=cut


sub fxGC {
    my ($self, $fxGC) = @_;

    if (defined $fxGC) {
	# is this an integer?
	if ($fxGC =~ /[^.\d]/) {
	    $self->throw(  -class => 'Bio::Root::BadParameter',
			   -text  => "Fraction GC must be a number between 0, 1 - NOT <$fxGC>",
			   -value => $fxGC
			   );
	}
	if  ( $fxGC < 0 or $fxGC > 1 ) {
	    $self->throw( -class => 'Bio::Root::BadParameter',
			  -text  => "Fraction GC must be a number between 0, 1 - NOT <$fxGC>",
			   -value => $fxGC
);
	}
	    
	#  clear out old tags
	$self->remove_tag('fxGC') if ( $self->has_tag('fxGC') );
	$self->add_tag_value('fxGC', $fxGC)
	    or $self->throw("Unable to set fxGC");
    }
    else {
	if ($self->has_tag('fxGC')) {
	    my @fxGCs = $self->get_tag_values('fxGC');
	    return shift @fxGCs;
	}
	else {
	    $self->throw("FxGC not defined for this Pair");
	}
    }
}

=head2 sense

  Title		: sense
  Usage		: my $sense_oligo = $sirna_pair->sense()
  Purpose	: Get/set the SiRNA::Oligo object corresponding to the sense strand
  Returns 	: Bio::SeqFeature::SiRNA::Oligo object
  Args		: Bio::SeqFeature::SiRNA::Oligo object

=cut


sub sense {
    my ($self, $soligo) = @_;

    if ($soligo) {
	$self->_add_oligo($soligo, 1) or return;
    }
    else {
	return $self->_get_oligo(1);
    }
}

=head2 antisense

  Title		: antisense
  Usage		: my $antisense_oligo = $sirna_pair->antisense()
  Purpose	: Get/set the SiRNA::Oligo object corresponding to the antisense strand
  Returns 	: Bio::SeqFeature::SiRNA::Oligo object
  Args		: Bio::SeqFeature::SiRNA::Oligo object

=cut

sub antisense {
    my ($self, $asoligo) = @_;

    if ($asoligo) {
	$self->_add_oligo($asoligo, -1) or return;
    }
    else {
	return $self->_get_oligo(-1);
    }
}
	
sub _add_oligo {
    my ($self, $oligo, $strand) = @_;

    unless ($oligo->isa('Bio::SeqFeature::SiRNA::Oligo')) {
	$self->throw( -class => 'Bio::Root::BadParameter',
		      -text  =>  "Oligos must be passed as Bio::SeqFeature::SiRNA::Oligo objects\n");	
    }

    $oligo->strand($strand);
    return $self->add_sub_SeqFeature($oligo, 'EXPAND');
}

sub _get_oligo {
    my ($self, $strand) = @_;
    my $feat;

    my @feats = $self->sub_SeqFeature;

    foreach $feat (@feats) {
	next unless ($feat->primary_tag eq 'SiRNA::Oligo');
	next unless ($feat->strand == $strand);
	return $feat;
    }
    return;
}

1;
