package Bio::SeqFeature::SiRNA::Oligo;

require 5.005_62;
use strict;
use warnings;

use Bio::SeqFeature::Generic;

our @ISA = qw(Bio::SeqFeature::Generic);

our $VERSION = '1.0';

our @ARGNAMES = qw(SEQ START END STRAND PRIMARY SOURCE_TAG SCORE TAG SEQ_ID ANNOTATION LOCATION);


1;

=head1 NAME

Bio::SeqFeature::SiRNA::Oligo - Perl object for small inhibitory RNAs.

=head1 SYNOPSIS

  use Bio::SeqFeature::SiRNA::Oligo;

  my $oligo = Bio::SeqFeature::SiRNA::Oligo->new( -seq		=> 'AUGCCGAUUGCAAGUCAGATT',
						  -start 	=> 10,
						  -end		=> 31,
						  -strand	=> -1,
						  -primary	=> 'SiRNA::Oligo',
						  -source_tag	=> 'Bio::Tools::SiRNA',
						  -tag		=> { note => 'A note' }, );

  # normally two complementary Oligos are combined in an SiRNA::Pair object
  $pair->antisense($oligo);

						  
=head1 DESCRIPTION

Object methods for single SiRNA oligos - inherits Bio::SeqFeature::Generic.
DOES NOT include methods for designing SiRNAs -- see Bio::Tools::SiRNA.pm for that.

=head2 EXPORT

None by default.

=head1 METHODS

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

  Note		: Currently passing arguments in gff_string or gff1_string is not supported.
                  SiRNA::Oligo objects are typically created by a design 
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
	    return undef;
	}
	else {
	    $self->{'seq'} = $seq;
	}
    }
    return $self->{'seq'};
}



=head1 AUTHOR

Donald Jackson (donald.jackson@bms.com)

=head1 SEE ALSO

Bio::Tools::SiRNA.pm, Bio::SeqFeature::SiRNA::Pair.pm, perl(1).

=cut
