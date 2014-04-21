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

Bio::SeqFeature::SiRNA::Oligo - Perl object for small inhibitory RNAs.

=head1 SYNOPSIS

  use Bio::SeqFeature::SiRNA::Oligo;

  my $oligo = Bio::SeqFeature::SiRNA::Oligo->
      new( -seq		=> 'AUGCCGAUUGCAAGUCAGATT',
	   -start 	=> 10,
	   -end		=> 31,
	   -strand	=> -1,
	   -primary	=> 'SiRNA::Oligo',
	   -source_tag	=> 'Bio::Tools::SiRNA',
	   -tag		=> { note => 'A note' }, );

  # normally two complementary Oligos are combined in an SiRNA::Pair
  # object
  $pair->antisense($oligo);


=head1 DESCRIPTION

Object methods for single SiRNA oligos - inherits
L<Bio::SeqFeature::Generic>.  Does B<not> include methods for designing
SiRNAs - see L<Bio::Tools::SiRNA> for that.

=head1 SEE ALSO

L<Bio::Tools::SiRNA>, L<Bio::SeqFeature::SiRNA::Pair>.

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

package Bio::SeqFeature::SiRNA::Oligo;

use strict;
use warnings;

use base qw(Bio::SeqFeature::Generic);

our @ARGNAMES = qw(SEQ START END STRAND PRIMARY SOURCE_TAG SCORE TAG 
                   SEQ_ID ANNOTATION LOCATION);

=head2 new

  Title		: new
  Usage		: my $sirna_oligo = Bio::SeqFeature::SiRNA::Oligo->new();
  Function	: Create a new SiRNA::Oligo object
  Returns	: Bio::Tools::SiRNA object
  Args    	: -seq		  sequence of the RNAi oligo.  Should be in RNA alphabet
                                  except for the final TT overhang.
                  -start          start position
 	 	  -end            end position
 	 	  -strand         strand
 	 	  -primary        primary tag - defaults to 'SiRNA::Oligo'
 	 	  -source         source tag
 	 	  -score          score value
 	 	  -tag            a reference to a tag/value hash
 	 	  -seq_id         the display name of the sequence
 	 	  -annotation     the AnnotationCollectionI object
 	 	  -location       the LocationI object

Currently passing arguments in gff_string or gff1_string is not
supported.  SiRNA::Oligo objects are typically created by a design
algorithm such as Bio::Tools::SiRNA

=cut

sub new {
    my ($proto, @args) = @_;

    my $pkg = ref($proto) || $proto;

    my (%args);

    my $self = $pkg->SUPER::new();

    @args{@ARGNAMES} = $self->_rearrange(\@ARGNAMES, @args); 
    # default primary tag
    $args{'PRIMARY'} ||= 'SiRNA::Oligo';

    $args{'PRIMARY'}		&& $self->primary_tag($args{'PRIMARY'});
    $args{'SOURCE_TAG'}		&& $self->source_tag($args{'SOURCE_TAG'});
    $args{'SEQNAME'}		&& $self->seqname($args{'SEQNAME'});
    $args{'SEQ'}		&& $self->seq($args{'SEQ'});
    $args{'ANNOTATION'}		&& $self->annotation($args{'ANNOTATION'});
    $args{'LOCATION'}		&& $self->location($args{'LOCATION'});
    defined($args{'START'})	&& $self->start($args{'START'});
    defined($args{'END'})	&& $self->end($args{'END'});
    defined($args{'STRAND'})	&& $self->strand($args{'STRAND'});
    defined($args{'SCORE'})	&& $self->score($args{'SCORE'});

    if ($args{'TAG'}) {	
	foreach my $t ( keys %{ $args{'TAG'} } ) {
	    $self->add_tag_value($t, $args{'TAG'}->{$t});
	}
    }

    return $self;
}

=head2 seq

  Title		: Seq
  Usage		: my $oligo_sequence = $sirna_oligo->seq();
  Purpose	: Get/set the sequence of the RNAi oligo
  Returns 	: Sequence for the RNAi oligo
  Args		: Sequence of the RNAi oligo (optional)
  Note		: Overloads Bio::SeqFeature::Generic seq method - the oligo and 
                  parent sequences are different. 
                  Note that all but the last 2 nucleotides are RNA (per Tuschl and colleagues).
                  SiRNA::Pair objects are typically created by a design algorithm such as
                  Bio::Tools::SiRNA.

=cut

sub seq {
    my ($self, $seq) = @_;
    if ($seq) {
	# check alphabet
	if ($seq =~ /[^ACGTUacgtu]/ ) {
	    warn "Sequence contains illegal characters";
	    return;
	}
	else {
	    $self->{'seq'} = $seq;
	}
    }
    return $self->{'seq'};
}

1;
