# $Id: GeneMap.pm,v 1.17 2006/07/17 14:16:53 sendu Exp $
#
# BioPerl module for Bio::Map::GeneMap
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu bala
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::GeneMap - A MapI implementation to represent the area around a gene

=head1 SYNOPSIS

    use Bio::Map::GeneMap;
    use Bio::Map::Gene;
    use Bio::Map::TranscriptionFactor;
    use Bio::Map::GeneRelative;

	# make some maps that will represent an area around a particular gene in
	# particular species (by default, the map represents the area in the genome
    # 1000bp upstream of the gene)
    my $map1 = Bio::Map::GeneMap->get(-gene => 'BRCA2',
                                      -species => 'human',
                                      -description => 'breast cancer 2, early onset');
	my $map2 = Bio::Map::GeneMap->get(-gene => 'BRCA2',
                                      -species => 'mouse');

	# model a TF that binds 500bp upstream of the BRCA2 gene in humans and
	# 250bp upstream of BRCA2 in mice
	my $rel = Bio::Map::GeneRelative->new(-description => "gene start");
    my $tf = Bio::Map::TranscriptionFactor->get(-universal_name => 'tf1');
	Bio::Map::Position->new(-map => $map1,
                            -element => $tf,
                            -start => -500,
                            -length => 10,
                            -relative => $rel);
	Bio::Map::Position->new(-map => $map2,
                            -element => $tf,
                            -start => -250,
                            -length => 10,
                            -relative => $rel);

	# find out all the things that map near BRCA2 in all species
	foreach my $map ($gene->known_maps) {
		foreach my $thing ($map->get_elements) {
            next if $thing eq $gene;
            foreach my $pos ($thing->get_positions($map)) {
                print "In species ", $map->species, ", ",
                      $thing->universal_name, " maps at ", $pos->value,
                      " relative to ", $pos->relative->description, " of gene ",
                      $gene->universal_name, "\n";
            }
		}
	}
    
    # a GeneMap isa PrimarySeq and so can have sequence associated with it
    $map1->seq('ATGC');
    my $subseq = $map1->subseq(2,3); # TG

=head1 DESCRIPTION

Model the abstract notion of the area around a gene - you don't care exactly
where this area is in the genome, you just want to be able to say "something
binds upstream of gene X" and "something else binds 20bp upstream of the first
something" etc.

It's useful for modelling transcription factor bindings sites, letting you find
out which transcription factors bind near a gene of interest, or which genes
are bound by a transcription factor of interest.

See t/Map/Map.t for more example usage.

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

package Bio::Map::GeneMap;
use strict;

use Bio::Map::Gene;
use Bio::Map::Position;

use base qw(Bio::Map::SimpleMap Bio::PrimarySeq);

our $GENEMAPS = {};

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Map::GeneMap->new();
 Function: Builds a new Bio::Map::GeneMap object (that has placed on it a
           mappable element (Bio::Map::Gene) representing a gene).
 Returns : Bio::Map::GeneMap
 Args    : -gene        => string name of the gene this map will be for
                           (in a form common to all species that have the gene,
                           but unique amongst non-orthologous genes) or a
                           Bio::Map::Gene object, REQUIRED
           -species     => Bio::Taxon or string representing species, REQUIRED
           -uid         => string, unique identifier for this map (must be
                           unique amongst all gene/species combinations)
           -description => string, free text description of the gene
           -upstream    => int, the number of bases the map extends before the
                           start of the gene element (default 1000).
           -downstream  => int, the number of bases the map extends beyond the
                           end of the gene element (default 0).
           -seq         => string, the sequence of the map, presumably the
                           genomic sequence -upstream bases of the gene,
                           including the gene, and -downstream bases of the gene

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my ($uid, $gene, $species, $desc, $up, $down, $seq) = $self->_rearrange([qw(UID
                                                    GENE
                                                    SPECIES
                                                    DESCRIPTION
                                                    UPSTREAM
                                                    DOWNSTREAM
                                                    SEQ)], @args);
    
    unless (defined $gene && defined $species) {
        $self->throw("You must supply both -species and -gene");
    }
    
    $self->gene(-gene => $gene, -description => $desc, -upstream => $up, -downstream => $down);
    $self->seq($seq) if $seq;
    
    unless (defined($uid)) {
        # trigger the special behaviour in our unique_id method by supplying it
        # the unique_id we got from our parent class
        $self->unique_id($self->unique_id);
    }
    
    return $self;
}

=head2 get

 Title   : get
 Usage   : my $map = Bio::Map::GeneMap->get();
 Function: Builds a new Bio::Map::GeneMap object (like new()), or gets a
           pre-existing one that corresponds to your arguments.
 Returns : Bio::Map::GeneMap
 Args    : -gene        => string name of the gene this map will be for
                           (in a form common to all species that have the gene,
                           but unique amongst non-orthologous genes) or a
                           Bio::Map::Gene object, REQUIRED
           -species     => Bio::Taxon or string representing species, REQUIRED
           -uid         => string, unique identifier for this map (must be
                           unique amongst all gene/species combinations)
           -description => string, free text description of the gene
           -upstream    => int, the number of bases the map extends before the
                           start of the gene element (default 1000).
           -downstream  => int, the number of bases the map extends beyond the
                           end of the gene element (default 0).
           -seq         => string, the sequence of the map, presumably the
                           genomic sequence -upstream bases of the gene,
                           including the gene, and -downstream bases of the gene

           If you supply a -uid, and a map had previously been created and
           given that uid, that same map object will be returned. Otherwise, the
           combination of -gene and -species will be used to determine
           if the same map had previously been made. If a corresponding map
           hadn't previously been made, a new map object will be created and
           returned.

=cut

sub get {
    my ($class, @args) = @_;
    my ($uid, $gene, $species, $desc, $up, $down, $seq) = Bio::Root::Root->_rearrange([qw(UID
                                                    GENE
                                                    SPECIES
                                                    DESCRIPTION
                                                    UPSTREAM
                                                    DOWNSTREAM
                                                    SEQ)], @args);
    
    my $gene_map;
    if ($uid && defined $GENEMAPS->{by_uid}->{$uid}) {
        $gene_map = $GENEMAPS->{by_uid}->{$uid};
    }
    elsif ($gene && $species) {
        my $name = ref($gene) ? $gene->universal_name : $gene;
        if (defined $GENEMAPS->{by_ns}->{$name}->{$species}) {
            $gene_map = $GENEMAPS->{by_ns}->{$name}->{$species};
        }
    }
    if ($gene_map) {
        $gene_map->gene->description($desc) if $desc;
        $gene_map->upstream($up) if defined($up);
        $gene_map->downstream($down) if defined($down);
        $gene_map->seq($seq) if $seq;
        return $gene_map;
    }
    
    return $class->new(@args);
}

=head2 unique_id

 Title   : unique_id
 Usage   : my $id = $map->unique_id;
 Function: Get/set the unique ID for this map
 Returns : string
 Args    : none to get, OR string to set

=cut

sub unique_id {
    my ($self, $id) = @_;
    if (defined $id) {
        delete $GENEMAPS->{by_uid}->{$self->{'_uid'}};
        $self->{'_uid'} = $id;
        $GENEMAPS->{by_uid}->{$id} = $self;
    }
    return $self->{'_uid'};
}

=head2 species

 Title   : species
 Usage   : my $species = $map->species;
 Function: Get/set Species for a map. It is not recommended to change this once
           set.
 Returns : Bio::Taxon object or string
 Args    : none to get, OR Bio::Taxon or string to set

=cut

sub species {
    my ($self, $value) = @_;
    if ($value) {
        my $old_species = $self->{_species};
        $self->{'_species'} = $value;
        my $name = $self->universal_name || return $value;
        if ($old_species) {
            delete $GENEMAPS->{by_ns}->{$name}->{$old_species};
        }
        $GENEMAPS->{by_ns}->{$name}->{$value} = $self;
    }
    return $self->{'_species'};
}

=head2 type

 Title   : type
 Usage   : my $type = $map->type
 Function: Get Map type
 Returns : string 'gene'
 Args    : none

=cut

sub type {
    return 'gene';
}

=head2 gene

 Title   : gene
 Usage   : my $gene = $map->gene;
           $map->gene(-gene => $gene);
 Function: Get/set the mappable element on this map that represents the gene
           this map is for. Once set, it is not recommended to re-set the gene
           to something else. Behaviour in that case is undefined.
 Returns : Bio::Map::Gene
 Args    : none to get, OR to set:
           -gene        => Bio::Map::Gene or string of the universal name (see
                           Bio::Map::Gene docs), REQUIRED
           -description => string, applied to the Bio::Map::Gene
           -upstream    => int, the number of bases the map extends before the
                           start of the gene element (default 1000).
           -downstream  => int, the number of bases the map extends beyond the
                           end of the gene element (default 0).

=cut

sub gene {
    my ($self, @args) = @_;
    
    if (@args > 0) {
        my ($gene, $desc, $up, $down) = $self->_rearrange([qw(GENE
                                                    DESCRIPTION
                                                    UPSTREAM
                                                    DOWNSTREAM)], @args);
        $self->throw("You must supply -gene") unless $gene;
        
        my $gene_obj = ref($gene) ? $gene : Bio::Map::Gene->get(-universal_name => $gene, -description => $desc);
        if (defined $self->{gene}) {
            if ($self->{gene} ne $gene_obj) {
                $self->warn("Changing the gene that this map is for, which could be bad");
                $self->purge_positions($self->{gene});
                delete $GENEMAPS->{by_ns}->{$self->universal_name}->{$self->species};
                $self->{gene} = $gene_obj;
            }
            
            # change the gene's position on us if necessary
            $self->upstream($up) if defined $up;
            $self->downstream($down) if defined $down;
        }
        else {
            # give the gene object a position on us
            $up ||= 1000;
            $up >= 0 || $self->throw("-upstream must be a positive integer");
            Bio::Map::Position->new(-map => $self, -start => ($up + 1), -element => $gene_obj);
            $self->{gene} = $gene_obj;
            $self->downstream($down || 0);
            
            # set other gene positions from db if already user-requested
            $gene_obj->_set_from_db($self);
        }
        
        $GENEMAPS->{by_ns}->{$self->universal_name}->{$self->species} = $self;
    }
    
    return $self->{gene};
}

=head2 universal_name

 Title   : universal_name
 Usage   : my $name = $map->universal_name
 Function: Get/set the name of Bio::Map::Gene object associated with this map.
           It is not recommended to change this once set.
 Returns : string
 Args    : none to get, OR string to set

=cut

sub universal_name {
    my ($self, $value) = @_;
    $self->gene || return;
    if ($value) {
        my $species = $self->species;
        delete $GENEMAPS->{by_ns}->{$self->gene->universal_name}->{$species};
        $self->gene->universal_name($value);
        $GENEMAPS->{by_ns}->{$value}->{$species} = $self;
    }
    return $self->gene->universal_name;
}

=head2 upstream

 Title   : upstream
 Usage   : my $distance = $map->upstream;
           $map->upstream($distance);
 Function: Get/set how long the map is before the start of the Bio::Map::Gene
           object on this map.
 Returns : int
 Args    : none to get, OR int to set (the number of bases the map extends
           before the start of the gene)

=cut

sub upstream {
    my ($self, $value) = @_;
    
    my $pos = $self->gene->position($self);
    if (defined($value)) {
        $value >= 0 || $self->throw("Supplied value must be a positive integer");
        $pos->start($value + 1);
    }
    
    return $pos->start - 1;
}

=head2 downstream

 Title   : downstream
 Usage   : my $distance = $map->downstream;
           $map->downstream($distance);
 Function: Get/set the nominal end of the map relative to the end of the
           Bio::Map::Gene object on this map.
 Returns : int
 Args    : none to get, OR int to set (the number of bases the map extends
           beyond the end of the gene)

=cut

sub downstream {
    my $self = shift;
    if (@_) { $self->{_downstream} = shift }
    return $self->{_downstream} || 0;
}

=head2 length

 Title   : length
 Usage   : my $length = $map->length();
 Function: Retrieves the length of the map. This is normally the length of the
           upstream region + length of the gene + length of the downstream
           region, but may be longer if positions have been placed on the map
           beyond the end of the nominal downstream region.
 Returns : int
 Args    : none

=cut

sub length {
	my $self = shift;
	my $expected_length = $self->gene->position($self)->length + $self->upstream + $self->downstream;
    my $actual_length = $self->SUPER::length;
    return $actual_length > $expected_length ? $actual_length : $expected_length;
}

=head2 seq

 Title   : seq
 Usage   : $string = $obj->seq()
 Function: Get/set the sequence as a string of letters. When getting, If the
           GeneMap object didn't have sequence attached directly to it for the
           region requested, the map's gene's database will be asked for the
           sequence, and failing that, the map's gene's positions will be asked
           for their sequences. Areas for which no sequence could be found will
           be filled with Ns, unless no sequence was found anywhere, in which
           case undef is returned.
 Returns : string
 Args    : Optionally on set the new value (a string). An optional second
           argument presets the alphabet (otherwise it will be guessed).

=cut

sub seq {
    my ($self, @args) = @_;
    my $seq = $self->SUPER::seq(@args);
    my $expected_length = $self->length;
    if (! $seq || CORE::length($seq) < $expected_length) {
        my @have = split('', $seq || '');
        my @result;
        for (0..($expected_length - 1)) {
            $result[$_] = shift(@have) || 'N';
        }
        
        # build map sequence by asking gene or positions
        my @slice_stuff = $self->gene->_get_slice($self);
        if (@slice_stuff) {
            my ($slice_adaptor, $slice, $strand) = @slice_stuff;
            my ($start, $end, $gene_start) = (CORE::length($seq || '') + 1, $expected_length, $self->upstream + 1);
            
            # convert map coords to genomic coords
            my $adjust = $strand == -1 ? $slice->end : $slice->start;
            my $adjustment = sub { return $strand == -1 ? $adjust - shift() : shift() + $adjust; };
            my $converted_start = &$adjustment($start - $gene_start);
            my $converted_end = &$adjustment($end - $gene_start);
            ($converted_start, $converted_end) = ($converted_end, $converted_start) if $converted_start > $converted_end;
            
            # get sequence from a new slice of desired region
            #*** what happens if desired region starts or ends off end of chromo?...
            my $new_slice = $slice_adaptor->fetch_by_region($slice->coord_system_name, $slice->seq_region_name, $converted_start, $converted_end);
            if ($new_slice && (my $seq_str = $new_slice->seq)) {
                if ($strand == -1) {
                    $seq_str = $self->_revcom($seq_str);
                }
                splice(@result, CORE::length($seq || ''), CORE::length($seq_str), split('', $seq_str));
            }
        }
        else {
            foreach my $pos ($self->get_positions) {
                next unless $pos->can('seq');
                my @pos_seq = split('', $pos->seq(undef, undef, 1) || next);
                for my $i ($pos->start($pos->absolute_relative)..$pos->end($pos->absolute_relative)) {
                    $i--;
                    my $base = shift(@pos_seq);
                    if ($result[$i] eq 'N') {
                        $result[$i] = $base;
                    }
                }
            }
        }
        
        $seq = join('', @result);
    }
    return $seq;
}

=head2 subseq

 Title   : subseq
 Usage   : $substring = $obj->subseq(10, 40);
 Function: Returns the subseq from start to end, where the first base
           is 1 and the number is inclusive, ie 1-2 are the first two
           bases of the sequence. If the GeneMap object didn't have sequence
           attached directly to it for the region requested, the map's gene's
           database will be asked for the sequence, and failing that, the map's
           gene's positions will be asked for their sequences. Areas for which
           no sequence could be found will be filled with Ns, unless no
           sequence was found anywhere, in which case undef is returned. subseq
           requests that extend beyond the end of the map will throw.
 Returns : string
 Args    : integer for start position AND integer for end position
                 OR
           Bio::LocationI location for subseq (strand honored)
                 OR
           Bio::RangeI (eg. a Bio::Map::PositionI)

=cut

sub subseq {
    my ($self, $start, $end) = @_;
    
    if ($start && ref($start) && $start->isa('Bio::RangeI')) {
        my $thing = $start;
        if ($start->isa('Bio::Map::Position')) {
            ($start, $end) = ($thing->start($thing->absolute_relative), $thing->end($thing->absolute_relative));
        }
        else {
            ($start, $end) = ($thing->start, $thing->end);
        }
    }
    
    # *** this implementation potentially wastefull? Should duplicate code
    #     from seq() to do this just for the desired region??
    my $orig_seq = $self->{seq};
    $self->{seq} = $self->seq();
    my $subseq = $self->{seq} ? $self->SUPER::subseq($start, $end) : '';
    $self->{seq} = $orig_seq;
    
    return $subseq;
}

# quick revcom for strings (silly to create a PrimarySeq just to revcom and then
# return a string again)
sub _revcom {
    my ($self, $seq) = @_;
    $seq or return;
    $seq = reverse($seq);
    $seq =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    return $seq;
}

1;
