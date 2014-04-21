# $Id: PositionWithSequence.pm,v 1.19 2006/09/20 10:20:01 sendu Exp $
#
# BioPerl module for Bio::Map::PositionWithSequence
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::PositionWithSequence - A position with a sequence.

=head1 SYNOPSIS

    use Bio::Map::PositionWithSequence;
    
    my $pos = Bio::Map::PositionWithSequence->new(-map => $map, 
                                -element => $element,
                                -start => 0,
                                -seq => 'ATGC');


=head1 DESCRIPTION

Have a position with a sequence, eg. define what the binding site sequence of
a certain transcription factor binding site is by modelling it as one of these
objects with the -element assigned to a Bio::Map::TranscriptionFactor instance.

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

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::PositionWithSequence;
use strict;

use base qw(Bio::Map::Position Bio::LocatableSeq);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Map::PositionWithSequence->new();
 Function: Builds a new Bio::Map::PositionWithSequence object 
 Returns : Bio::Map::PositionWithSequence
 Args    : -map      => Bio::Map::GeneMap object
           -element  => Bio::Map::Gene object
           -relative => Bio::Map::GeneRelative object
           -seq      => string, length of this string will set the length
                        of this position's range

           * If this position has no range, or if a single value can describe
             the range *
           -value => scalar             : something that describes the single
                                          point position or range of this
                                          Position, most likely an int

           * Or if this position has a range, at least two of *
           -start => int                : value of the start co-ordinate
           -end => int                  : value of the end co-ordinate
           -length => int               : length of the range

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my ($seq) = $self->_rearrange([qw( SEQ )], @args);
    
    $self->seq($seq) if $seq;
    
    return $self;
}

=head2 seq

 Title   : seq
 Usage   : my $string = $obj->seq();
 Function: Get/set the sequence as a string of letters.
 Returns : scalar
 Args    : Optionally on set the new value (a string). An optional second
           argument presets the alphabet (otherwise it will be guessed).

=cut

sub seq {
    my ($self, $str, $alpha) = @_;
    
    # done like this because SUPER will set seq to undef if undef supplied,
    # but GeneMap wants to send undef, undef, 1 to decendants of this method
    
    my $seq;
    if ($str) {
        $alpha ? ($seq = $self->SUPER::seq($str, $alpha)) : ($seq = $self->SUPER::seq($str));
    }
    else {
        $seq = $self->SUPER::seq;
    }
    
    if ($seq) {
        $self->length(length($seq));
        return $seq;
    }
    return;
}

1;
