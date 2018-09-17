# $Id: Gene.pm,v 1.6 2006/07/17 14:16:53 sendu Exp $
#
# BioPerl module for Bio::Map::Gene
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

Bio::Map::Gene - An gene modelled as a mappable element.

=head1 SYNOPSIS

  use Bio::Map::Gene;

  my $gene = Bio::Map::Gene->get(-universal_name => 'BRCA2',
                                 -description => 'breast cancer 2, early onset');

  # Normally you get Gene objects from GeneMaps
  use Bio::Map::GeneMap;

  # Model a gene with its orthologous versions found in different species,
  # but at abstract locations within each genome
  my $map1 = Bio::Map::GeneMap->get(-universal_name => 'BRCA2', -species => $human);
  my $map2 = Bio::Map::GeneMap->get(-universal_name => 'BRCA2', -species => $mouse);

  $gene = $map1->gene;

  # Genes can have special kinds of positions (Bio::Map::GenePosition) that
  # define where various sub-regions of the gene are, relative to one of the
  # normal Positions the gene has placing it on a map.
  my $trans = Bio::Map::GenePosition->new(-start => 0, -length => 700,
                                          -map => $map1, -type => 'transcript');
  $gene->add_transcript_position($trans);
  my $exon = Bio::Map::GenePosition->new(-start => 0, -length => 100,
                                         -map => $map1, -type => 'exon');
  $gene->add_exon_position($exon, 1);
  # (so now the gene has 1 transcript 700bp long which starts at the beginning
  #  of the gene, and we've defined the first of many exons which starts at the
  #  start of the transcript and is 100bp long)

=head1 DESCRIPTION

Model a gene as an abstract mappable element. This is for when you don't care
exactly where a gene is in a genome, but just want to model other things (like
transcription factor binding sites) that are near it so you can answer questions
like "what binds near this gene?", or "which genes does this bind near?".

See t/Map/Map.t for more example usage.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

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

package Bio::Map::Gene;
use strict;

use Bio::Map::GenePosition;

use base qw(Bio::Map::Mappable);

our $USE_ENSEMBL;
our $GENES = {};
our $SET_FROM_DB = 0;

BEGIN {
    # Bio::Tools::Run::Ensembl is in bioperl-run package which may not be
    # installed, but its functionality is only optional here
    eval {require Bio::Tools::Run::Ensembl;};
    $USE_ENSEMBL = ! $@;
}

=head2 new

 Title   : new
 Usage   : my $gene = Bio::Map::Gene->new();
 Function: Builds a new Bio::Map::Gene object
 Returns : Bio::Map::Gene
 Args    : -universal_name => string : name of the gene (in a form common to all
                                       species that have the gene, but unique
                                       amongst non-orthologous genes), REQUIRED
           -description => string    : free text description of the gene

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my ($u_name, $desc) = $self->_rearrange([qw(UNIVERSAL_NAME DESCRIPTION)], @args);
    $u_name || $self->throw("You must supply a -universal_name");
    $self->universal_name($u_name);
    
    defined $desc && $self->description($desc);
    
    return $self;
}

=head2 get

 Title   : get
 Usage   : my $gene = Bio::Map::Gene->get();
 Function: Builds a new Bio::Map::Gene object (like new()), or gets a
           pre-existing one that shares the same universal_name.
 Returns : Bio::Map::Gene
 Args    : -universal_name => string, name of the gene (in a form common to all
                              species that have the gene, but unique amongst
                              non-orthologous genes), REQUIRED
           -description    => string, free text description of the gene

=cut

sub get {
    my ($class, @args) = @_;
    my ($u_name, $desc) = Bio::Root::Root->_rearrange([qw(UNIVERSAL_NAME DESCRIPTION)], @args);
    
    if ($u_name && defined $GENES->{$u_name}) {
        $GENES->{$u_name}->description($desc) if $desc;
        return $GENES->{$u_name};
    }
    
    return $class->new(@args);
}

=head2 universal_name

 Title   : universal_name
 Usage   : my $name = $gene->universal_name
 Function: Get/Set Mappable name, corresponding to the name of the gene in a
           form shared by orthologous versions of the gene in different species,
           but otherwise unique.
 Returns : string
 Args    : none to get, OR string to set

=cut

sub universal_name {
    my ($self, $value) = @_;
    if (defined $value) {
        delete $GENES->{$self->{'_uname'}} if $self->{'_uname'};
        $self->{'_uname'} = $value;
        $GENES->{$value} = $self;
    }
    return $self->{'_uname'};
}

=head2 description

 Title   : description
 Usage   : my $description = $gene->description();
           $gene->description($description);
 Function: Get/set information relating to the gene, in this case the
           description (eg. 'full name of gene')
 Returns : string (empty string if not defined)
 Args    : none to get general version, OR Bio::Map::GeneMap to get map-specific
           version.
           string to set general version, optionally AND Bio::Map::GeneMap to
           set map-specific version

=cut

sub description {
    my $self = shift;
    return $self->_gene_data('description', @_);
}

=head2 display_id

 Title   : display_id
 Usage   : my $display_id = $gene->display_id();
           $gene->display_id($display_id);
 Function: Get/set information relating to the gene, in this case the
           display_id (eg. 'ENSG00000155287')
 Returns : string (empty string if not defined)
 Args    : none to get general version, OR Bio::Map::GeneMap to get map-specific
           version.
           string to set general version, optionally AND Bio::Map::GeneMap to
           set map-specific version

=cut

sub display_id {
    my $self = shift;
    return $self->_gene_data('display_id', @_);
}

=head2 display_xref

 Title   : display_xref
 Usage   : my $display_xref = $gene->display_xref();
           $gene->display_xref($display_xref);
 Function: Get/set information relating to the gene, in this case the
           display_xref (eg. 'HUGO:23472').
 Returns : string (empty string if not defined)
 Args    : none to get general version, OR Bio::Map::GeneMap to get map-specific
           version.
           string to set general version, optionally AND Bio::Map::GeneMap to
           set map-specific version

=cut

sub display_xref {
    my $self = shift;
    return $self->_gene_data('display_xref', @_);
}

=head2 external_db

 Title   : external_db
 Usage   : my $external_db = $gene->external_db();
           $gene->external_db($external_db);
 Function: Get/set information relating to the gene, in this case the
           external_db (eg. 'HUGO').
 Returns : string (empty string if not defined)
 Args    : none to get general version, OR Bio::Map::GeneMap to get map-specific
           version.
           string to set general version, optionally AND Bio::Map::GeneMap to
           set map-specific version

=cut

sub external_db {
    my $self = shift;
    return $self->_gene_data('external_db', @_);
}

=head2 external_name

 Title   : external_name
 Usage   : my $external_name = $gene->external_name();
           $gene->external_name($external_name);
 Function: Get/set information relating to the gene, in this case the (eg.
           'gene_name', probably the same as or similar to what you set
           universal_name() to, but could be a species-specific alternative).
 Returns : string (empty string if not defined)
 Args    : none to get general version, OR Bio::Map::GeneMap to get map-specific
           version.
           string to set general version, optionally AND Bio::Map::GeneMap to
           set map-specific version

=cut

sub external_name {
    my $self = shift;
    return $self->_gene_data('external_name', @_);
}

=head2 biotype

 Title   : biotype
 Usage   : my $biotype = $gene->biotype();
           $gene->biotype($biotype);
 Function: Get/set information relating to the gene, in this case the biotype
           (eg. 'protein_coding').
 Returns : string (empty string if not defined)
 Args    : none to get general version, OR Bio::Map::GeneMap to get map-specific
           version.
           string to set general version, optionally AND Bio::Map::GeneMap to
           set map-specific version

=cut

sub biotype {
    my $self = shift;
    return $self->_gene_data('biotype', @_);
}

=head2 source

 Title   : source
 Usage   : my $source = $gene->source();
           $gene->source($source);
 Function: Get/set information relating to the gene, in this case the source
           (eg. '??').
 Returns : string (empty string if not defined)
 Args    : none to get general version, OR Bio::Map::GeneMap to get map-specific
           version.
           string to set general version, optionally AND Bio::Map::GeneMap to
           set map-specific version

=cut

sub source {
    my $self = shift;
    return $self->_gene_data('source', @_);
}

=head2 position

 Title   : position
 Usage   : my $position = $mappable->position($map);
 Function: Get the main Position of this Mappable on a given map. (A gene may
           have many positions on a map, but all but one of them are
           Bio::Map::GenePosition objects that describe sub-regions of the gene
           which are relative to the 'main' Bio::Map::Position position, which
           is the only one that is directly relative to the map - this is the
           Position returned by this method.)
 Returns : Bio::Map::Position
 Args    : L<Bio::Map::MapI> object.

=cut

sub position {
    my ($self, $map) = @_;
    ($map && $self->in_map($map)) || return;
    
    foreach my $pos ($self->get_positions($map, 1)) {
        next if $pos->isa('Bio::Map::GenePosition');
        return $pos;
        #*** could do sanity checking; there should only be 1 non-GenePosition
        #    object here, and it should have a relative of type 'map', and it
        #    should sort before or equal to all other positions
    }
}

=head2 add_transcript_position

 Title   : add_transcript_position
 Usage   : $gene->add_transcript_position($position);
 Function: Set the bounds of a transcript on a map (that of the supplied
           position). All transcript positions added this way must have
           coordinates relative to the main position of the 'gene' mappable on
           this transcript's map. The first position added using this method
           must have a start of 0. The supplied Position will be given a type of
           'transcript' and relative of (gene => 0). The active_transcript for
           the Position's map will be set to this one.
 Returns : n/a
 Args    : Bio::Map::GenePosition (which must have its map() defined, and be for
           a map this gene is on)

=cut

sub add_transcript_position {
    my ($self, $pos) = @_;
    ($pos && $pos->isa('Bio::Map::GenePosition')) || return;
    
    my $map = $pos->map || $self->throw("Supplied GenePosition has no map");
    $self->in_map($map) || $self->throw("Supplied GenePosition is not on a map that this gene belong to");
    my @transcripts = $self->get_transcript_positions($map);
    if (@transcripts == 0) {
        # first transcript needs start of 0
        if ($pos->start != 0) {
            $self->warn("The first transcript position added to a map needs a start of 0, not adding");
            return;
        }
    }
    
    $pos->type('transcript');
    $pos->relative->gene(0);
    $self->SUPER::add_position($pos);
    
    # need to remember the order these were added, but remember what we store
    # here could become invalid if positions are purged outside of this class
    push(@{$self->{t_order}->{$map}}, $pos);
    
    # adjust main position's length to hold this transcript
    my $main_pos = $self->position($map);
    my $increase = ($pos->length + $pos->start($pos->absolute_relative)) - ($main_pos->end + 1);
    if ($increase > 0) {
        $main_pos->end($main_pos->end + $increase);
    }
    
    # make this new transcript the active one
    $self->active_transcript($map, scalar(@transcripts) + 1);
}

=head2 active_transcript

 Title   : active_transcript
 Usage   : my $active = $gene->active_transcript($map);
           $gene->active_transcript($map, $int);
 Function: Get/set the active transcript number (an int of 1 would mean the 1st
           transcript position added to the object for the given map, ie. would
           correspond to the the 1st Position object in the list returned by
           get_transcript_positions($map)). The active transcript is the one
           considered by other methods and objects when dealing with positions
           relative to 'the' transcript.
 Returns : int, 0 means there were no transcript positions on the given map,
           undef is some other problem
 Args    : Just Bio::Map::GeneMap to get
           Bio::Map::GeneMap AND int to set

=cut

sub active_transcript {
    my ($self, $map, $int) = @_;
    $map or return;
    
    my @transcripts = $self->get_transcript_positions($map);
    if (@transcripts > 0) {
        if (defined($int)) {
            if ($int > 0 && $int <= @transcripts) {
                $self->{active_transcript}->{$map} = $int;
                return $int;
            }
            else {
                $self->warn("Supplied int '$int' not a good number (higher than the number of transcripts on the map?)");
                return;
            }
        }
        else {
            if (defined $self->{active_transcript}->{$map}) {
                return $self->{active_transcript}->{$map};
            }
            else {
                # default to the total number of transcripts on the map, ie. the
                # most recently added
                $self->{active_transcript}->{$map} = @transcripts;
                return $self->{active_transcript}->{$map};
            }
        }
    }
    return 0;
}

=head2 get_transcript_positions

 Title   : get_transcript_positions
 Usage   : my @transcript_positions = $gene->get_transcript_positions($map);
 Function: Get all the transcript positions of this gene on the given map, in
           the order they were added to the map.
 Returns : list of Bio::Map::GenePosition
 Args    : Bio::Map::GeneMap

=cut

sub get_transcript_positions {
    my ($self, $map) = @_;
    $map or return;
    $map->isa('Bio::Map::GeneMap') or return;
    return $self->_get_typed_positions($map, 'transcript');
}

=head2 get_transcript_position

 Title   : get_transcript_position
 Usage   : my $position = $gene->get_transcript_position($map, $int);
 Function: Get the $int'th transcript position added to the map. If no
           transcripts have been added to the map, and the default transcript
           was requested, $gene->position is returned, as that will have the
           same start and end as the first transcript.
 Returns : Bio::Map::GenePosition
 Args    : Bio::Map::GeneMap AND int (if int not supplied, or 0, returns
           the currently active transcript position)

=cut

sub get_transcript_position {
    my ($self, $map, $value) = @_;
    $map or return;
    $value ||= $self->active_transcript($map);
    my @transcripts = $self->get_transcript_positions($map);
    if (@transcripts == 0 && $value == 0) {
        return $self->position($map);
    }
    return $self->_get_list_element($value, @transcripts);
}

=head2 coding_position

 Title   : coding_position
 Usage   : $gene->coding_position($position, $transcript_number);
           $gene->coding_position($map, $transcript_number);
 Function: Get/set the bounds of a coding region of a given transcript on a map
           (that of the supplied position).

           When setting, coordinates must be relative to the transcript start.
           The supplied position will be given a type 'coding' and a relative
           (-transcript => $transcript_number). There can be only one coding
           position per transcript (hence this is a get/set).

           When getting, if a coding region has not been defined for the
           requested transcript, $gene->get_transcript_position($map,
           $transcript_number) is returned, as if assuming the entirety of the
           transcript is coding.

 Returns : Bio::Map::GenePosition
 Args    : Bio::Map::GeneMap AND int (the transcript number) to get, OR to set:
           Bio::Map::GenePosition (which must have its map() defined, and be for
           a map this gene is on) AND int (the transcript number)
           In both cases, if transcript number not supplied or 0 this will be
           resolved to the current active transcript number - there must be at
           least one transcript on the map

=cut

sub coding_position {
    my ($self, $thing, $transcript_num) = @_;
    ref($thing) || return;
    $transcript_num ||= 0;
    
    # deliberate test for PositionI so _add_type_position can do nothing if
    # its not a GenePosition
    if ($thing->isa('Bio::Map::PositionI')) {
        my $map = $thing->map || return;
        my ($existing_pos) = $self->_get_typed_positions($map, 'coding', $transcript_num);
        if ($existing_pos) {
            # purge it
            $self->purge_positions($existing_pos);
        }
        $self->_add_type_position('coding', $thing, $transcript_num);
        $thing = $map;
    }
    
    my ($pos) = $self->_get_typed_positions($thing, 'coding', $transcript_num);
    return $pos || $self->get_transcript_position($thing, $transcript_num);
}

=head2 add_exon_position

 Title   : add_exon_position
 Usage   : $gene->add_exon_position($position, $transcript_number);
 Function: Set the bounds of an exon of a given transcript on a map (that of the
           supplied position). Coordinates must be relative to the transcript
           start. The supplied position will be given a type 'exon' and a
           relative (-transcript => $transcript_number).
 Returns : n/a
 Args    : Bio::Map::GenePosition (which must have its map() defined, and be for
           a map this gene is on) AND int (the transcript number; if not
           supplied or 0 this will be resolved to the current active transcript
           number - there must be at least one transcript on the map)

=cut

sub add_exon_position {
    my $self = shift;
    $self->_add_type_position('exon', @_);
}

=head2 get_exon_positions

 Title   : get_exon_positions
 Usage   : my @positions = $gene->get_exon_positions($map, $int);
 Function: Get all the exon positions that are relative to the $int'th
           transcript position added to the map. Exons are returned sorted by
           their start positions.
 Returns : array of Bio::Map::GenePosition
 Args    : Bio::Map::GeneMap AND int (the transcript number; if second int not
           supplied, or 0, considers the currently active transcript)

=cut

sub get_exon_positions {
    my ($self, $map, $value) = @_;
    $map || return;
    $value ||= 0;
    return $self->_get_typed_positions($map, 'exon', $value);
}

=head2 get_exon_position

 Title   : get_exon_position
 Usage   : my $position = $gene->get_exon_position($map, $exon_num, $int);
 Function: Get the $exon_num'th exon position that is relative to the $int'th
           transcript position added to the map. Exons are numbered in Position
           order, not the order they were added to the map. If no exons have
           been added to the map, and the first exon was requested,
           $gene->get_transcript_position($map, $int) is returned, as that will
           have the same start as the first exon, and could have the same end
           for a single exon gene.
 Returns : Bio::Map::GenePosition
 Args    : Bio::Map::GeneMap AND int (the exon you want) AND int (the transcript
           number; if second int not supplied, or 0, considers the currently
           active transcript)

=cut

sub get_exon_position {
    my ($self, $map, $exon_num, $value) = @_;
    my @exons = $self->get_exon_positions($map, $value);
    if (@exons == 0 && $exon_num == 1) {
        return $self->get_transcript_position($map, $value);
    }
    return $self->_get_list_element($exon_num, @exons);
}

=head2 add_intron_position

 Title   : add_intron_position
 Usage   : $gene->add_intron_position($position, $transcript_number);
 Function: Set the bounds of an intron of a given transcript on a map (that of
           the supplied position). Coordinates must be relative to the
           transcript start. The supplied position will be given a type 'intron'
           and a relative (-transcript => $transcript_number).
 Returns : n/a
 Args    : Bio::Map::GenePosition (which must have its map() defined, and be for
           a map this gene is on) AND int (the transcript number; if not
           supplied or 0 this will be resolved to the current active transcript
           number - there must be at least one transcript on the map)

=cut

sub add_intron_position {
    my $self = shift;
    $self->_add_type_position('intron', @_);
}

=head2 get_intron_positions

 Title   : get_intron_positions
 Usage   : my @positions = $gene->get_intron_positions($map, $int);
 Function: Get all the intron positions that are relative to the $int'th
           transcript position added to the map. Introns are returned sorted by
           their start positions.
 Returns : array of Bio::Map::GenePosition
 Args    : Bio::Map::GeneMap AND int (the transcript number; if second int not
           supplied, or 0, considers the currently active transcript)

=cut

sub get_intron_positions {
    my ($self, $map, $value) = @_;
    $map || return;
    $value ||= 0;
    return $self->_get_typed_positions($map, 'intron', $value);
}

=head2 get_intron_position

 Title   : get_intron_position
 Usage   : my $position = $gene->get_intron_position($map, $intron_num, $int);
 Function: Get the $intron_num'th intron position that is relative to the
           $int'th transcript position added to the map. Introns are numbered in
           Position order, not the order they were added to the map.
 Returns : Bio::Map::GenePosition
 Args    : Bio::Map::GeneMap AND int (the intron you want) AND int (the
           transcript number; if second int not supplied, or 0, considers the
           currently active transcript)

=cut

sub get_intron_position {
    my ($self, $map, $intron_num, $value) = @_;
    my @introns = $self->get_intron_positions($map, $value);
    return $self->_get_list_element($intron_num, @introns);
}

=head2 set_from_db

 Title   : set_from_db
 Usage   : $gene->set_from_db(); # for an instance only
           Bio::Map::Gene->set_from_db(); # decide that all future genes added
                                          # to maps will be set from db
 Function: Creates all the various types of positions (transcripts, coding,
           exons, introns) for this gene on all its maps. The information comes
           from an Ensembl database via Bio::Tools::Run::Ensembl. NB: will
           purge any existing Bio::Map::GenePosition objects that were
           previously on the maps this gene is one.
 Returns : undef on failure, otherwise the number of maps that successfully
           had positions added to them
 Args    : boolean (no argument/undef is treated as 1, ie. do set from db;
           supply 0 to turn off)

           NB: Bio::Tools::Run::Ensembl is available in the bioperl-run package;
           see it for details on setting up a database to use.

           Once set, any new maps (species) this gene is added to will
           automatically also have their positions set_from_db

=cut

sub set_from_db {
    my ($self, $bool) = @_;
    return unless $USE_ENSEMBL;
    return unless Bio::Tools::Run::Ensembl->registry_setup();
    defined($bool) || ($bool = 1);
    
    unless (ref($self)) {
        $SET_FROM_DB = $bool;
        return 0;
    }
    
    $self->{_set_from_db} = $bool;
    
    my $success = 0;
    foreach my $map ($self->known_maps) {
        $success += $self->_set_from_db($map);
    }
    
    return $success;
}

# set from db for a particular map (species)
sub _set_from_db {
    my ($self, $map) = @_;
    my $gene_name = $self->universal_name || return 0;
    $SET_FROM_DB || $self->{_set_from_db} || return;
    
    my $species = $map->species;
    
    my $slice_adaptor = Bio::Tools::Run::Ensembl->get_adaptor($species, 'Slice') || return 0;
    my $gene = Bio::Tools::Run::Ensembl->get_gene_by_name(-species => $species,
                                                          -name => $gene_name,
                                                          -use_orthologues => 'Homo sapiens',
                                                          -use_swiss_lookup => 1,
                                                          -use_entrez_lookup => 1) || return 0;
    
    # attach species(map)-specific gene info to self
    $self->description($gene->description, $map);
    $self->display_id($gene->display_id, $map);
    $self->display_xref($gene->display_xref->display_id, $map);
    $self->external_db($gene->external_db, $map);
    $self->external_name($gene->external_name, $map);
    $self->biotype($gene->biotype, $map);
    $self->source($gene->source, $map);
    
    # get the transcripts for this map
    my $trans_ref = $gene->get_all_Transcripts;
    unless ($trans_ref && @{$trans_ref} > 0) {
        return 0;
    }
    
    # purge all existing GenePositions from the map
    my $handler = $map->get_position_handler();
    foreach my $pos ($map->get_positions) {
        if ($pos->isa('Bio::Map::GenePosition')) {
            $handler->purge_positions($pos);
        }
    }
    
    # assume all transcripts on the same strand, sort them
    my $strand = ${$trans_ref}[0]->strand;
    my @transcripts = sort { $strand == -1 ? ($b->end <=> $a->end) : ($a->start <=> $b->start) } @{$trans_ref};
    
    # store slice of first transcript so we can use it to get seq data, and
    # add chromosome info to our map if not set
    my $primary_slice = $slice_adaptor->fetch_by_transcript_stable_id($transcripts[0]->stable_id, 0);
    my $uid = $map->unique_id;
    @{$self->{_ensembl}->{$uid}} = ($slice_adaptor, $primary_slice, $strand);
    
    #my $cyto = $map->location || Bio::Map::CytoPosition->new();
    #unless ($cyto->chr) {
    #    $cyto->chr($primary_slice->seq_region_name);
    #}
    #$map->location($cyto);
    
    # adjustment needed to make all transcript coords relative to the start of
    # the first transcript which must start at 0
    my $adjust = $strand == -1 ? $transcripts[0]->end : $transcripts[0]->start;
    my $orig_adjust = $adjust;
    my $adjustment = sub { return $strand == -1 ? $adjust - shift() : shift() - $adjust; };
    
    # go through all the transcripts, remembering the longest
    my $longest_trans = 0;
    my $longest = 1;
    my $count = 1;
    foreach my $transcript (@transcripts) {
        # length is the total number of bases the exons cover, not genomic span
        my $length = $transcript->length();
        if ($length > $longest_trans) {
            $longest_trans = $length;
            $longest = $count;
        }
        
        # make positions for this transcript
        my $slice = $slice_adaptor->fetch_by_transcript_stable_id($transcript->stable_id, 0);
        my $start = &$adjustment($slice->start());
        my $end = &$adjustment($slice->end());
        ($start, $end) = ($end, $start) if $start > $end;
        
        my $trans_pos = Bio::Map::GenePosition->new(-map => $map, -start => $start, -end => $end, -type => 'transcript');
        $self->add_transcript_position($trans_pos);
        
        # all subsequent coordinates need to be relative to the start of this
        # transcript
        $adjust = $strand == -1 ? $slice->end : $slice->start;
        
        # there may not be a coding region
        if (defined($transcript->coding_region_start)) {
            my $atg = &$adjustment($transcript->coding_region_start());
            my $stop = &$adjustment($transcript->coding_region_end());
            ($atg, $stop) = ($stop, $atg) if $atg > $stop;
            
            my $cod_pos = Bio::Map::GenePosition->new(-map => $map, -start => $atg, -end => $stop, -type => 'coding');
            $self->coding_position($cod_pos);
        }
        
        # exons
        foreach my $exon (@{$transcript->get_all_Exons}) {
            my $start = &$adjustment($exon->start());
            my $end = &$adjustment($exon->end());
            ($start, $end) = ($end, $start) if $start > $end;
            
            my $throw_species = ref($species) ? $species->scientific_name : $species;
            defined($end) || $self->throw("gene $gene_name in species $throw_species (".$gene->display_id.") had exon $start with no end");
            my $pos = Bio::Map::GenePosition->new(-map => $map, -start => $start, -end => $end, -type => 'exon');
            $self->add_exon_position($pos);
        }
        
        # introns
        foreach my $intron (@{$transcript->get_all_Introns}) {
            my $start = &$adjustment($intron->start());
            my $end = &$adjustment($intron->end());
            ($start, $end) = ($end, $start) if $start > $end;
            
            my $pos = Bio::Map::GenePosition->new(-map => $map, -start => $start, -end => $end, -type => 'intron');
            $self->add_intron_position($pos);
        }
        
        $adjust = $orig_adjust;
    } continue { $count++ };
    
    $self->active_transcript($map, $longest);
    
    return 1;
}

# get safely sorted positions of a certain type
sub _get_typed_positions {
    my ($self, $map, $type, $transcript_number) = @_;
    if (defined $transcript_number && $transcript_number == 0) {
        $transcript_number = $self->active_transcript($map);
    }
    
    my @positions;
    foreach my $pos ($self->get_positions($map, 1)) {
        $pos->isa('Bio::Map::GenePosition') || next;
        $pos->type eq $type || next;
        
        if (defined $transcript_number) {
            my $rel = $pos->relative || next;
            $rel->type eq 'transcript' || next;
            my $rel_transcript_num = $rel->transcript || $self->active_transcript($map);
            $rel_transcript_num == $transcript_number || next;
        }
        
        push(@positions, $pos);
    }
    
    # avoid sorting using $pos->sortable since we would go infinite from the
    # call to absolute_conversion - we don't need absolute_conversion here
    # since we know the raw starts are all relative to the same thing, or in
    # the case of transcripts, we want them sorted in the way they were added
    if (defined $transcript_number) {
        # ensure we get raw start; ask for starts relative to the things
        # the positions are relative to. Precompute answer for efficiency
        my @sort = map { $_->[1] }
                   sort { $a->[0] <=> $b->[0] }
                   map { [$_->start($_->relative), $_] }
                   @positions;
        return @sort;
    }
    else {
        my @known_order = @{$self->{t_order}->{$map} || []};
        @known_order || return;
        
        # transcripts might have been removed, so known_order could be invalid
        return @known_order if @known_order == @positions; #*** dangerous assumption?
        my %exists = map { $_ => $_ } @positions;
        my @new_order;
        foreach my $pos (@known_order) {
            exists $exists{$pos} || next;
            push(@new_order, $pos);
        }
        @{$self->{t_order}->{$map}} = @new_order;
        return @new_order;
    }
}

# get a certain element from an array, checking the array has that element
sub _get_list_element {
    my ($self, $wanted, @list) = @_;
    ($wanted && $wanted > 0) || return;
    @list > 0 || return;
    my $index = $wanted - 1;
    if ($index >= 0 && $index <= $#list) {
        return $list[$index];
    }
    return;
}

# add a certain type of posiiton
sub _add_type_position {
    my ($self, $type, $pos, $transcript_num) = @_;
    ($pos && $pos->isa('Bio::Map::GenePosition')) || return;
    
    my $map = $pos->map || $self->throw("Supplied GenePosition has no map");
    $self->in_map($map) || $self->throw("Supplied GenePosition is not on a map that this gene belong to");
    
    $transcript_num ||= $self->active_transcript($map) || $self->throw("Asked to be relative to the active transcript, but there is no transcript");
    
    # sanity check - must be within the transcript
    my $transcript_pos = $self->get_transcript_position($map, $transcript_num) || $self->throw("Asked to be relative to transcript $transcript_num, but there is no such transcript");
    $transcript_pos->end || ($self->warn("no transcript pos end for pos for gene ".$self->universal_name." and species ".$pos->map->species."!") && exit);
    $pos->end || ($self->warn("no pos end for pos for gene ".$self->universal_name." and species ".$pos->map->species."!") && exit);
    unless ($transcript_pos->contains($pos)) {
        $self->warn("$type coordinates must lie within those of the transcript, not adding $type");
        return;
    }
    
    $pos->type($type);
    $pos->relative->transcript($transcript_num);
    $self->SUPER::add_position($pos);
}

# get/setter for general/map-specific data
sub _gene_data {
    my ($self, $type, $thing, $map) = @_;
    $thing or return ($self->{$type}->{general} || '');
    
    if (ref($thing) && $thing->isa('Bio::Map::GeneMap')) {
        return $self->{$type}->{$thing} || '';
    }
    
    if ($map && $map->isa('Bio::Map::GeneMap')) {
        $self->{$type}->{$map} = $thing;
    }
    else {
        $self->{$type}->{general} = $thing;
    }
    return $thing;
}

# for exclusive use by GeneMap so it can get sequence data
sub _get_slice {
    my ($self, $map) = @_;
    $map || return;
    my $uid = $map->unique_id || return;
    if (defined $self->{_ensembl}->{$uid}) {
        return @{$self->{_ensembl}->{$uid}};
    }
    return;
}

1;
