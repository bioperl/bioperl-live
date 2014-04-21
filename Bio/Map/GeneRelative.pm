# $Id: GeneRelative.pm,v 1.6 2006/09/20 11:53:29 sendu Exp $
#
# BioPerl module for Bio::Map::GeneRelative
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

Bio::Map::GeneRelative - Represents being relative to named sub-regions of a
                         gene.

=head1 SYNOPSIS

    use Bio::Map::GeneRelative;

    # say that a somthing will have a position relative to the start of the
    # gene on map
    my $rel = Bio::Map::GeneRelative->new(-gene => 0);

    # or that something will be relative to the third transcript of a gene
    # on a map
    $rel = Bio::Map::GeneRelative->new(-transcript => 3);

    # or to the 5th intron of the default transcript
    $rel = Bio::Map::GeneRelative->new(-intron => [0, 5]);

    # use the $rel as normal; see L<Bio::Map::Relative>

=head1 DESCRIPTION

Be able to say that a given position is relative to some standard part of a
gene.

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

package Bio::Map::GeneRelative;
use strict;

use Scalar::Util qw(looks_like_number);

use base qw(Bio::Map::Relative);

=head2 new

 Title   : new
 Usage   : my $relative = Bio::Map::Relative->new();
 Function: Build a new Bio::Map::Relative object.
 Returns : Bio::Map::Relative object
 Args    : -gene => int       : coordinates are relative to the int'th base
                                downstream of the Position's map's gene
                                [default is gene => 0, ie. relative to the
                                start of the gene],
           -transcript => int : or relative to the start of the int'th
                                transcript of the Position's map's gene,
           -exon => [i, n]    : or relative to the start of the n'th
                                transcript's i'th exon,
           -intron => [i, n]  : or intron,
           -coding => int     : or the start of the int'th transcript's coding
                                region.

           -description => string : Free text description of what this relative
                                    describes

           (To say a Position is relative to something and upstream of it,
            the Position's start() co-ordinate should be set negative)
           In all cases, a transcript number of 0 means the active transcript.

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my ($gene, $transcript, $exon, $intron, $coding) =
        $self->_rearrange([qw( GENE TRANSCRIPT EXON INTRON CODING )], @args);
    
    my $set = (defined $gene) + (defined $transcript) + (defined $exon) + (defined $intron) + (defined $coding);
    if ($set > 1) {
        $self->throw("-gene, -transcript, -exon, -intron and -coding are mutually exclusive");
    }
    if ($exon && (! ref($exon) || ref($exon) ne 'ARRAY')) {
        $self->throw("-exon takes an array ref");
    }
    if ($intron && (! ref($intron) || ref($intron) ne 'ARRAY')) {
        $self->throw("-intron takes an array ref");
    }
    if ($set == 0) {
        # type could have been set already in the call to SUPER::new
        if ($self->type) {
            $self->warn("You set a type of relative not supported by GeneRelative; resetting to type 'gene'");
        }
        $gene = 0;
    }
    
    $self->gene($gene) if defined $gene;
    $self->transcript($transcript) if defined $transcript;
    $self->exon(@{$exon}) if defined $exon;
    $self->intron(@{$intron}) if defined $intron;
    $self->coding($coding) if defined $coding;
    
    return $self;
}

=head2 absolute_conversion

 Title   : absolute_conversion
 Usage   : my $absolute_coord = $relative->absolute_conversion($pos);
 Function: Convert the start co-ordinate of the supplied position into a number
           relative to the start of its map.
 Returns : scalar number
 Args    : Bio::Map::PositionI object

=cut

sub absolute_conversion {
    my ($self, $pos) = @_;
    $self->throw("Must supply an object") unless ref($pos);
    $self->throw("This is [$pos], not a Bio::Map::PositionI") unless $pos->isa('Bio::Map::PositionI');
    
    # get the raw start position of our position
    my $raw = $pos->start($pos->relative);
    $self->throw("Can't convert co-ordinates when start isn't set") unless defined($raw); #*** needed? return undef?
    
    # what are we relative to?
    my $type = $self->type;
    my $value = $self->$type;
    $self->throw("Details not yet set for this Relative, cannot convert") unless defined($value);
    
    # get the absolute start of the thing we're relative to
    if ($type =~ /gene|transcript|exon|intron|coding/) {
        my $map = $pos->map;
        my $throw_desc = $type eq 'gene' ? 'gene' : "gene's transcript";
        $self->throw("Relative to a map's $throw_desc, but the Position has no map") unless $map;
        $self->throw("Relative to a map's $throw_desc, but the Position's map isn't a Bio::Map::GeneMap") unless $map->isa('Bio::Map::GeneMap');
        my $gene = $map->gene;
        
        if ($type eq 'gene') {
            my $gene_pos = $gene->position($map);
            my $rel = $gene_pos->relative;
            my $start = $rel->absolute_conversion($gene_pos);
            $value += $start;
        }
        else {
            my @values = ref($value) ? @{$value} : ($value);
            my $trans = ref($value) ? $values[1] : $value;
            my $throw_txt = $trans == 0 ? 'default/active transcript' : "transcript $trans";
            my $throw_txt2 = ref($value) ? ", or no $type $values[0]" : '';
            my $method = $type eq 'coding' ? 'coding_position' : "get_${type}_position";
            $value = $gene->$method($map, @values) || $self->throw("Relative to $throw_txt of the map's gene, but there is no such transcript$throw_txt2");
        }
    }
    else {
        return $self->SUPER::absolute_conversion($pos);
    }
    if (ref($value)) {
        # psuedo-recurse
        my $rel = $value->relative;
        $value = $rel->absolute_conversion($value);
    }
    
    if (defined($value)) {
        return $value + $raw;
    }
    return;
}

=head2 type

 Title   : type
 Usage   : my $type = $relative->type();
 Function: Get the type of thing we are relative to. The types correspond
           to a method name, so the value of what we are relative to can
           subsequently be found by $value = $relative->$type;

           Note that type is set by the last method that was set, or during
           new().

 Returns : 'gene', 'transcript', 'exon', 'intron' or 'coding'
 Args    : none

=cut

=head2 gene

 Title   : gene
 Usage   : my $int = $relative->gene();
           $relative->gene($int);
 Function: Get/set the distance from the start of the gene that the Position's
           co-ordiantes are relative to.
 Returns : int
 Args    : none to get, OR
           int to set; a value of 0 means relative to the start of the gene.

=cut

sub gene {
    my ($self, $num) = @_;
    if (defined($num)) {
        $self->throw("This is [$num], not a number") unless looks_like_number($num);
        $self->{_use} = 'gene';
        $self->{_gene} = $num;
    }
    return defined($self->{_gene}) ? $self->{_gene} : return;
}

=head2 transcript

 Title   : transcript
 Usage   : my $int = $relative->transcript();
           $relative->transcript($int);
 Function: Get/set which transcript of the Position's map's gene the Position's
           co-ordinates are relative to.
 Returns : int
 Args    : none to get, OR
           int to set; a value of 0 means the active (default) transcript.

=cut

sub transcript {
    my ($self, $num) = @_;
    if (defined($num)) {
        $self->throw("This is [$num], not a number") unless looks_like_number($num);
        $self->{_use} = 'transcript';
        $self->{_transcript} = $num;
    }
    return defined($self->{_transcript}) ? $self->{_transcript} : return;
}

=head2 exon

 Title   : exon
 Usage   : my ($exon_number, $transcript_number) = @{$relative->exon()};
           $relative->exon($exon_number, $transcript_number);
 Function: Get/set which exon of which transcript of the Position's map's gene
           the Position's co-ordinates are relative to.
 Returns : reference to list with two ints, exon number and transcript number
 Args    : none to get, OR
           int (exon number) AND int (transcript number) to set. The second int
           is optional and defaults to 0 (meaning default/active transcript).

=cut

sub exon {
    my ($self, $num, $t_num) = @_;
    if (defined($num)) {
        if (defined($t_num)) {
            $self->throw("This is [$t_num], not a number") unless looks_like_number($t_num);
        }
        $t_num ||= 0;
        $self->throw("This is [$num], not a number") unless looks_like_number($num);
        $self->{_use} = 'exon';
        $self->{_exon} = [$num, $t_num];
    }
    return $self->{_exon} || return;
}

=head2 intron

 Title   : intron
 Usage   : my ($intron_number, $transcript_number) = @{$relative->intron()};
           $relative->intron($intron_number, $transcript_number);
 Function: Get/set which intron of which transcript of the Position's map's gene
           the Position's co-ordinates are relative to.
 Returns : reference to list with two ints, intron number and transcript number
 Args    : none to get, OR
           int (intron number) AND int (transcript number) to set. The second
           int is optional and defaults to 0 (meaning default/active
           transcript).

=cut

sub intron {
    my ($self, $num, $t_num) = @_;
    if (defined($num)) {
        if (defined($t_num)) {
            $self->throw("This is [$t_num], not a number") unless looks_like_number($t_num);
        }
        $t_num ||= 0;
        $self->throw("This is [$num], not a number") unless looks_like_number($num);
        $self->{_use} = 'intron';
        $self->{_intron} = [$num, $t_num];
    }
    return $self->{_intron} || return;
}

=head2 coding

 Title   : coding
 Usage   : my $transcript_number = $relative->coding;
           $relative->coding($transcript_number);
 Function: Get/set which transcript's coding region of the Position's map's gene
           the Position's co-ordinates are relative to.
 Returns : int
 Args    : none to get, OR
           int to set (the transcript number, see transcript())

=cut

sub coding {
    my ($self, $num) = @_;
    if (defined($num)) {
        $self->throw("This is [$num], not a number") unless looks_like_number($num);
        $self->{_use} = 'coding';
        $self->{_coding} = $num;
    }
    return defined($self->{_coding}) ? $self->{_coding} : return;
}

1;
