#
# BioPerl module for Bio::SeqFeature::Similarity
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

Bio::SeqFeature::Similarity - A sequence feature based on similarity

=head1 SYNOPSIS

    # obtain a similarity feature somehow
    print "significance: ", $sim_fea->significance(), "\n";
    print "bit score: ", $sim_fea->bits(), "\n";
    print "score: ", $sim_fea->score(), "\n";
    print "fraction of identical residues: ", $sim_fea->frac_identical(), "\n";

=head1 DESCRIPTION

This module is basically a sequence features based on similarity, and therefore
has support for measures assessing the similarity.

Everything else is inherited from L<Bio::SeqFeature::Generic>.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

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

Email hlapp@gmx.net or hilmar.lapp@pharma.novartis.com

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Similarity;
use strict;


use base qw(Bio::SeqFeature::Generic);

sub new {
    my ( $caller, @args) = @_;
    my ($self) = $caller->SUPER::new(@args);

    my ($primary,$evalue, $bits, $frac,$seqlen,$seqdesc) =
	$self->_rearrange([qw(PRIMARY
			      EXPECT
			      BITS
			      FRAC
			      SEQLENGTH
			      SEQDESC
			      )],@args);

    defined $evalue && $self->significance($evalue);
    defined $bits   && $self->bits($bits);
    defined $frac   && $self->frac_identical($frac);
    defined $seqlen && $self->seqlength($seqlen);
    defined $seqdesc && $self->seqdesc($seqdesc);
    $primary  = 'similarity' unless defined $primary;
    $self->primary_tag($primary) unless( defined $self->primary_tag() );
    $self->strand(0) unless( defined $self->strand() );

    return $self;
}

=head2 significance

 Title   : significance
 Usage   : $evalue = $obj->significance();
           $obj->significance($evalue);
 Function: 
 Returns : 
 Args    : 


=cut

sub significance {
    return shift->_tag_value('signif', @_);
}

=head2 bits

 Title   : bits
 Usage   : $bits = $obj->bits();
           $obj->bits($value);
 Function: 
 Returns : 
 Args    : 


=cut

sub bits {
    return shift->_tag_value('Bits', @_);
}

=head2 frac_identical

 Title   : frac_identical
 Usage   : $fracid = $obj->frac_identical();
           $obj->frac_identical($value);
 Function: 
 Returns : 
 Args    : 


=cut

sub frac_identical {
    return shift->_tag_value('FracId', @_);
}

=head2 seqlength

 Title   : seqlength
 Usage   : $len = $obj->seqlength();
           $obj->seqlength($len);
 Function: 
 Returns : 
 Args    : 


=cut

sub seqlength {
    return shift->_tag_value('SeqLength', @_);
}

=head2 seqdesc

 Title   : seqdesc
 Usage   : $desc = $obj->seqdesc();
           $obj->seqdesc($desc);
 Function: At present this method is a shorthand for 
           $obj->annotation()->description().

           Note that this is not stored in the tag system and hence will
           not be included in the return value of gff_string().
 Returns : 
 Args    : 


=cut

sub seqdesc {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        my $v = Bio::Annotation::SimpleValue->new();
        $v->value($value);
        $self->annotation->add_Annotation( 'description', $v );
    }
    my ($v) = $self->annotation()->get_Annotations('description');
    return defined $v ? $v->value : undef;
}

#
# Everything else is just inherited from SeqFeature::Generic.
#

1;
