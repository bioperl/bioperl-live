# $Id: GenePosition.pm,v 1.19 2006/09/20 10:20:01 sendu Exp $
#
# BioPerl module for Bio::Map::GenePosition
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

Bio::Map::GenePosition - A typed position, suitable for modelling the various
                         regions of a gene.

=head1 SYNOPSIS

    use Bio::Map::GenePosition;
    use Bio::Map::GeneMap;

    # say that the first transcript of a particular gene on a particular map
    # (species) is 1000bp long
    my $map = Bio::Map:GeneMap->get(-universal_name => 'BRCA2',
                                    -species => 'human');
    my $gene = $map->gene;
    Bio::Map::GenePosition->new(-map => $map, 
                                -element => $gene,
                                -start => 0,
                                -length => 1000,
                                -type => 'transcript');

    # say that the coding region of the gene starts 30bp into the first
    # transcript
    Bio::Map::GenePosition->new(-map => $map, 
                                -element => $gene,
                                -start => 30,
                                -length => 600,
                                -type => 'coding');

    # A GenePosition isa PositionWithSequence, so can have sequence associated
    # with it
    my $exon = Bio::Map::GenePosition->new(-map => $map, 
                                -element => $gene,
                                -start => 0,
                                -type => 'exon',
                                -seq => 'ATGGGGTGGG');
    my $length = $exon->length; # $length is 10

=head1 DESCRIPTION

Define where various sub-regions (transcripts, exons, introns etc.) of a gene
are. Do this so that you can then go onto to model other mappable elements as
having positions 20bp upstream of transcript 2, or 10bp into intron 3 etc., all
without having to know the absolute position of anything.

See Bio::Map::GeneRelative and t/Map/Map.t for more example usage.

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

package Bio::Map::GenePosition;
use strict;

use Bio::Map::GeneRelative;

use base qw(Bio::Map::PositionWithSequence);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Map::GenePosition->new();
 Function: Builds a new Bio::Map::GenePosition object 
 Returns : Bio::Map::GenePosition
 Args    : -map      => Bio::Map::GeneMap object
           -element  => Bio::Map::Gene object
           -relative => Bio::Map::GeneRelative object
           -type     => 'transcript|coding|exon|intron', REQUIRED
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
    
    my ($type) = $self->_rearrange([qw( TYPE )], @args);
	$type || $self->throw("type must be supplied");
    $self->type($type);
    
    $self->{_relative_not_implicit} = 1;
    
    return $self;
}

=head2 map

 Title   : map
 Usage   : my $map = $position->map();
           $position->map($map);
 Function: Get/set the map the position is in.
 Returns : L<Bio::Map::MapI>
 Args    : none to get
           new L<Bio::Map::MapI> to set

=cut

sub map {
    my ($self, $map) = @_;
    if ($map) {
        $map->isa('Bio::Map::GeneMap') || $self->throw("This is [$map], not a Bio::Map::GeneMap");
    }
    return $self->SUPER::map($map);
}

=head2 element

 Title   : element
 Usage   : my $element = $position->element();
           $position->element($element);
 Function: Get/set the element the position is for.
 Returns : L<Bio::Map::MappableI>
 Args    : none to get
           new L<Bio::Map::MappableI> to set

=cut

sub element {
    my ($self, $element) = @_;
    if ($element) {
        $element->isa('Bio::Map::Gene') || $self->throw("This is [$element], not a Bio::Map::Gene");
    }
    return $self->SUPER::element($element);
}

=head2 type

 Title   : type
 Usage   : my $type = $position->type();
           $position->type($type);
 Function: Get/set the type of this position.
 Returns : string
 Args    : none to get, OR
           string transcript|coding|exon|intron to set

=cut

sub type {
    my $self = shift;
    if (@_) {
        my $type = shift;
        if ($type !~ /transcript|coding|exon|intron/i) {
            $self->throw("type must be supplied and be one of 'transcript', 'coding', 'exon', 'intron'");
        }
        $self->{type} = $type;
    }
    return $self->{type};
}

=head2 relative

  Title   : relative
  Usage   : my $relative = $position->relative();
            $position->relative($relative);
  Function: Get/set the thing this Position's coordinates (numerical(), start(),
            end()) are relative to, as described by a RelativeI object.
  Returns : Bio::Map::GeneRelative. The default GeneRelative returned has a
            meaning that depends on the type() of GenePosition this is:
            'transcript'         : "relative to the start of the gene on the
                                    Position's map"
            'coding|exon|intron' : "relative to the start of the default
                                    transcript of the gene on the Position's map"
  Args    : none to get, OR
            Bio::Map::GeneRelative to set

=cut

sub relative {
    my ($self, $relative) = @_;
    if ($relative) {
        $self->throw("Must supply an object") unless ref($relative);
        $self->throw("This is [$relative], not a Bio::Map::GeneRelative") unless $relative->isa('Bio::Map::GeneRelative');
        $self->{_relative} = $relative;
    }
    return $self->{_relative} || $self->_default_relative;
}

=head2 seq

 Title   : seq
 Usage   : my $string = $position->seq();
 Function: Get/set the sequence as a string of letters. If no sequence is
           manually set by you, the position's map will be asked for the
           sequence, and if available, that will be returned.
 Returns : scalar
 Args    : Optionally on set the new value (a string). An optional second
           argument presets the alphabet (otherwise it will be guessed).

=cut

sub seq {
    # $shortcut is internal-use only by GeneMap
    my ($self, $str, $alpha, $shortcut) = @_;
    
    my $seq = $self->SUPER::seq($str, $alpha);
    
    if ($seq) {
        $self->length(CORE::length($seq));
        return $seq;
    }
    elsif (! $shortcut && defined(my $map = $self->map) && ! defined $self->{_getting_seq}) {
        $self->{_getting_seq} = 1;
        $seq = $map->subseq($self);
        delete $self->{_getting_seq};
        return $seq;
    }
    return;
}

# return a Relative that is suitable for the type
sub _default_relative {
    my $self = shift;
    my $type = $self->type;
    if ($type eq 'transcript') {
        return Bio::Map::GeneRelative->new(-gene => 0, -description => 'start of gene');
    }
    else {
        return Bio::Map::GeneRelative->new(-transcript => 0, -description => 'start of default transcript');
    }
}

1;
