# $Id$
#
# BioPerl module for Bio::SeqFeature::Similarity
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


package Bio::SeqFeature::Similarity;
use vars qw(@ISA);
use strict;

use Bio::SeqFeature::Generic;

@ISA = qw(Bio::SeqFeature::Generic);

sub new {
    my ( $caller, @args) = @_;   
    my ($self) = $caller->SUPER::new(@args); 

    my ($primary,$evalue, $bits, $frac,$seqlen,$seqdesc) =
	$self->_rearrange([qw(PRIMARY
			      EXPECT
			      BITS
			      FRAC
			      SEQDESC
			      SEQLENGTH				      
			      )],@args);

    $evalue && $self->significance($evalue);
    $bits   && $self->bits($bits);
    $frac   && $self->frac_identical($frac);
    $seqlen && $self->seqlength($seqlen);
    $seqdesc && $self->seqdesc($seqdesc);
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
    my ($self, $value) = @_;

    return $self->_tag_value('signif', $value);
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
    my ($self, $value) = @_;

    return $self->_tag_value('Bits', $value);
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
    my ($self, $value) = @_;

    return $self->_tag_value('FracId', $value);
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
    my ($self, $value) = @_;

    return $self->_tag_value('SeqLength', $value);
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
    my ($self, $value) = @_;

    if( defined $value ) { 
	my $v = Bio::Annotation::SimpleValue->new();
	$v->value($value);
	$self->annotation->add_Annotation('description',$v);
    }
    my ($v) = $self->annotation()->get_Annotations('description');
    return $v ? $v->value : undef;
}

#
# Everything else is just inherited from SeqFeature::Generic.
#

1;
