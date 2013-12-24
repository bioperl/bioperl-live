package Bio::Coordinate::GeneMapper;
use utf8;
use strict;
use warnings;
use Bio::Coordinate::Result;
use Bio::Location::Simple;
use Bio::Coordinate::Graph;
use Bio::Coordinate::Collection;
use Bio::Coordinate::Pair;
use Bio::Coordinate::ExtrapolatingPair;
use parent qw(Bio::Root::Root Bio::Coordinate::MapperI);

# ABSTRACT: Transformations between gene related coordinate systems.
# AUTHOR:   Heikki Lehvaslaiho <heikki@bioperl.org>
# OWNER:    Heikki Lehvaslaiho
# LICENSE:  Perl_5

=head1 SYNOPSIS

  use Bio::Coordinate::GeneMapper;

  # get a Bio::RangeI representing the start, end and strand of the CDS
  # in chromosomal (or entry) coordinates
  my $cds;

  # get a Bio::Location::Split or an array of Bio::LocationI objects
  # holding the start, end and strand of all the exons in chromosomal
  # (or entry) coordinates
  my $exons;

  # create a gene mapper and set it to map from chromosomal to cds coordinates
  my $gene = Bio::Coordinate::GeneMapper->new(-in   =>'chr',
                                              -out  =>'cds',
                                              -cds  =>$cds,
                                              -exons=>$exons
                                             );

  # get a a Bio::Location or sequence feature in input (chr) coordinates
  my $loc;

  # map the location into output coordinates and get a new location object
  $newloc = $gene->map($loc);

=head1 DESCRIPTION

Bio::Coordinate::GeneMapper is a module for simplifying the mappings
of coodinate locations between various gene related locations in human
genetics. It also adds a special human genetics twist to coordinate
systems by making it possible to disable the use of zero
(0). Locations before position one start from -1. See method
L<nozero>.

It understands by name the following coordinate systems and mapping
between them:

                          peptide (peptide length)
                             ^
                             | -peptide_offset
                             |
                    frame  propeptide (propeptide length)
                        ^    ^
                         \   |
             translate    \  |
                           \ |
                            cds  (transcript start and end)
                             ^
      negative_intron        | \
              ^              |  \  transcribe
               \             |   \
              intron        exon  \
               ^   ^         ^     /
      splice    \   \      / |    /
                 \   \    /  |   /
                  \   inex   |  /
                   \    ^    | /
                    \    \   |/
                     ----- gene (gene_length)
                             ^
                             | - gene_offset
                             |
                            chr (or entry)

This structure is kept in the global variable $DAG which is a
representation of a Directed Acyclic Graph. The path calculations
traversing this graph are done in a helper class. See
L<Bio::Coordinate::Graph>.

Of these, two operations are special cases, translate and splice.
Translating and reverse translating are implemented as internal
methods that do the simple 1E<lt>-E<gt>3 conversion. Splicing needs
additional information that is provided by method L<exons> which takes
in an array of Bio::LocationI objects.

Most of the coordinate system names should be selfexplanatory to
anyone familiar with genes. Negative intron coordinate system is
starts counting backwards from -1 as the last nucleotide in the
intron. This used when only exon and a few flanking intron nucleotides
are known.

This class models coordinates within one transcript of a gene, so to
tackle multiple transcripts you need several instances of the
class. It is therefore valid to argue that the name of the class
should be TranscriptMapper. GeneMapper is a catchier name, so it
stuck.

=cut

# first set internal values for all translation tables

our %COORDINATE_SYSTEMS = (
    peptide          => 10,
    propeptide       => 9,
    frame            => 8,
    cds              => 7,
    negative_intron  => 6,
    intron           => 5,
    exon             => 4,
    inex             => 3,
    gene             => 2,
    chr              => 1,
);

our %COORDINATE_INTS = (
    10 => 'peptide',
    9  => 'propeptide',
    8  => 'frame',
    7  => 'cds',
    6  => 'negative_intron',
    5  => 'intron',
    4  => 'exon',
    3  => 'inex',
    2  => 'gene',
    1  => 'chr'
);

our $TRANSLATION = $COORDINATE_SYSTEMS{'cds'}. "-". $COORDINATE_SYSTEMS{'propeptide'};

our $DAG = {
    10 => [],
    9  => [10],
    8  => [],
    7  => [8, 9],
    6  => [],
    5  => [6],
    4  => [7],
    3  => [4, 5],
    2  => [3, 4, 5, 7],
    1  => [2],
};

our $NOZERO_VALUES = {
    0        => 0,
    'in'     => 1,
    'out'    => 2,
    'in&out' => 3,
};

our $NOZERO_KEYS = {
    0 => 0,
    1 => 'in',
    2 => 'out',
    3 => 'in&out',
};

=head2 new
=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    # prime the graph
    my $graph = Bio::Coordinate::Graph->new();
    $graph->hash_of_arrays($DAG);
    $self->graph($graph);

    my($in, $out, $peptide_offset, $exons,
       $cds, $nozero, $strict) =
        $self->_rearrange([qw(IN
                              OUT
                              PEPTIDE_OFFSET
                              EXONS
                              CDS
                              NOZERO
                              STRICT
                             )],
                         @args);

    # direction of mapping when going chr to protein
    $self->{_direction} = 1;

    $in  && $self->in($in);
    $out  && $self->out($out);
    $cds && $self->cds($cds);
    $exons  && ref($exons) =~ /ARRAY/i && $self->exons(@$exons);
    $peptide_offset && $self->peptide_offset($peptide_offset);
    $nozero && $self->nozero($nozero);
    $strict && $self->strict($strict);

    return $self; # success - we hope!
}

=head2 in

 Title   : in
 Usage   : $obj->in('peptide');
 Function: Set and read the input coordinate system.
 Example :
 Returns : value of input system
 Args    : new value (optional)

=cut

sub in {
   my ($self,$value) = @_;
   if( defined $value) {
       $self->throw("Not a valid input coordinate system name [$value]\n".
                    "Valid values are ". join(", ", keys %COORDINATE_SYSTEMS ))
           unless defined $COORDINATE_SYSTEMS{$value};

       $self->{'_in'} = $COORDINATE_SYSTEMS{$value};
   }
   return $COORDINATE_INTS{ $self->{'_in'} };
}

=head2 out

 Title   : out
 Usage   : $obj->out('peptide');
 Function: Set and read the output coordinate system.
 Example :
 Returns : value of output system
 Args    : new value (optional)

=cut

sub out {
   my ($self,$value) = @_;
   if( defined $value) {
       $self->throw("Not a valid input coordinate system name [$value]\n".
                    "Valid values are ". join(", ", keys %COORDINATE_SYSTEMS ))
           unless defined $COORDINATE_SYSTEMS{$value};

       $self->{'_out'} = $COORDINATE_SYSTEMS{$value};
   }
   return $COORDINATE_INTS{ $self->{'_out'} };
}

=head2 strict

 Title   : strict
 Usage   : $obj->strict('peptide');
 Function: Set and read whether strict boundaried of coordinate
           systems are enforced.
           When strict is on, the end of the coordinate range must be defined.
 Example :
 Returns : boolean
 Args    : boolean (optional)

=cut

sub strict {
   my ($self,$value) = @_;
   if( defined $value) {
       $value ? ( $self->{'_strict'} = 1 ) : ( $self->{'_strict'} = 0 );
       ## update in each mapper !!
   }
   return $self->{'_strict'} || 0 ;
}

=head2 nozero

 Title   : nozero
 Usage   : $obj->nozero(1);
 Function: Flag to disable the use of zero in the input,
           output or both coordinate systems. Use of coordinate
           systems without zero is a peculiarity  common in
           human genetics community.
 Example :
 Returns : 0 (default), or 'in', 'out', 'in&out'
 Args    : 0 (default), or 'in', 'out', 'in&out'

=cut

sub nozero {
   my ($self,$value) = @_;

   if (defined $value) {
       $self->throw("Not a valid value for nozero [$value]\n".
                    "Valid values are ". join(", ", keys %{$NOZERO_VALUES} ))
           unless defined $NOZERO_VALUES->{$value};
       $self->{'_nozero'} = $NOZERO_VALUES->{$value};
   }

   my $res = $self->{'_nozero'} || 0;
   return $NOZERO_KEYS->{$res};
}

=head2 graph

 Title   : graph
 Usage   : $obj->graph($new_graph);
 Function: Set and read the graph object representing relationships
           between coordinate systems
 Example :
 Returns : Bio::Coordinate::Graph object
 Args    : new Bio::Coordinate::Graph object (optional)

=cut

sub graph {
   my ($self,$value) = @_;
   if( defined $value) {
       $self->throw("Not a valid graph [$value]\n")
           unless $value->isa('Bio::Coordinate::Graph');
       $self->{'_graph'} = $value;
   }
   return $self->{'_graph'};
}

=head2 peptide

 Title   : peptide
 Usage   : $obj->peptide_offset($peptide_coord);
 Function: Read and write the offset of peptide from the start of propeptide
           and peptide length
 Returns : a Bio::Location::Simple object
 Args    : a Bio::LocationI object

=cut

sub peptide {
   my ($self, $value) = @_;
   if( defined $value) {
       $self->throw("I need a Bio::LocationI, not  [". $value. "]")
           unless $value->isa('Bio::LocationI');

       $self->throw("Peptide start not defined")
           unless defined $value->start;
       $self->{'_peptide_offset'} = $value->start - 1;

       $self->throw("Peptide end not defined")
           unless defined $value->end;
       $self->{'_peptide_length'} = $value->end - $self->{'_peptide_offset'};

       my $a = $self->_create_pair
           ('propeptide', 'peptide', $self->strict,
            $self->{'_peptide_offset'}, $self->{'_peptide_length'} );
       my $mapper =  $COORDINATE_SYSTEMS{'propeptide'}. "-".  $COORDINATE_SYSTEMS{'peptide'};
       $self->{'_mappers'}->{$mapper} = $a;
   }
   return  Bio::Location::Simple->new
       (-seq_id => 'propeptide',
        -start => $self->{'_peptide_offset'} + 1 ,
        -end => $self->{'_peptide_length'} + $self->{'_peptide_offset'},
        -strand => 1,
        -verbose => $self->verbose,
       );
}

=head2 peptide_offset

 Title   : peptide_offset
 Usage   : $obj->peptide_offset(20);
 Function: Set and read the offset of peptide from the start of propeptide
 Returns : set value or 0
 Args    : new value (optional)

=cut

sub peptide_offset {
   my ($self,$offset, $len) = @_;
   if( defined $offset) {
       $self->throw("I need an integer, not [$offset]")
           unless $offset =~ /^[+-]?\d+$/;
       $self->{'_peptide_offset'} = $offset;

       if (defined $len) {
           $self->throw("I need an integer, not [$len]")
               unless $len =~ /^[+-]?\d+$/;
           $self->{'_peptide_length'} = $len;
       }

       my $a = $self->_create_pair
           ('propeptide', 'peptide', $self->strict, $offset, $self->{'_peptide_length'} );
       my $mapper =  $COORDINATE_SYSTEMS{'propeptide'}. "-". $COORDINATE_SYSTEMS{'peptide'};
       $self->{'_mappers'}->{$mapper} = $a;
   }
   return $self->{'_peptide_offset'} || 0;
}

=head2 peptide_length

 Title   : peptide_length
 Usage   : $obj->peptide_length(20);
 Function: Set and read the offset of peptide from the start of propeptide
 Returns : set value or 0
 Args    : new value (optional)

=cut

sub peptide_length {
   my ($self, $len) = @_;
   if( defined $len) {
       $self->throw("I need an integer, not [$len]")
           if defined $len && $len !~ /^[+-]?\d+$/;
       $self->{'_peptide_length'} = $len;
   }
   return $self->{'_peptide_length'};
}

=head2 exons

 Title   : exons
 Usage   : $obj->exons(@exons);
 Function: Set and read the offset of CDS from the start of transcript
           You do not have to sort the exons before calling this method as
           they will be sorted automatically.
           If you have not defined the CDS, is will be set to span all
           exons here.
 Returns : array of Bio::LocationI exons in genome coordinates or 0
 Args    : array of Bio::LocationI exons in genome (or entry) coordinates

=cut

sub exons {
   my ($self,@value) = @_;
   my $cds_mapper =  $COORDINATE_SYSTEMS{'gene'}. "-". $COORDINATE_SYSTEMS{'cds'};
   my $inex_mapper =
       $COORDINATE_SYSTEMS{'gene'}. "-". $COORDINATE_SYSTEMS{'inex'};
   my $exon_mapper =
       $COORDINATE_SYSTEMS{'gene'}. "-". $COORDINATE_SYSTEMS{'exon'};
   my $intron_mapper =
       $COORDINATE_SYSTEMS{'gene'}. "-". $COORDINATE_SYSTEMS{'intron'};
   my $negative_intron_mapper =
       $COORDINATE_SYSTEMS{'intron'}. "-". $COORDINATE_SYSTEMS{'negative_intron'};
   my $exon_cds_mapper =  $COORDINATE_SYSTEMS{'exon'}. "-". $COORDINATE_SYSTEMS{'cds'};

   if(@value) {
       if (ref($value[0]) &&
           $value[0]->isa('Bio::SeqFeatureI') and
           $value[0]->location->isa('Bio::Location::SplitLocationI')) {
           @value = $value[0]->location->each_Location;
       } else {
           $self->throw("I need an array , not [@value]")
               unless ref \@value eq 'ARRAY';
           $self->throw("I need a reference to an array of Bio::LocationIs, not to [".
                        $value[0]. "]")
               unless ref $value[0] and $value[0]->isa('Bio::LocationI');
       }

       #
       # sort the input array
       #
       # and if the used has not defined CDS assume it is the complete exonic range
       if (defined $value[0]->strand &&
           $value[0]->strand == - 1) {  #reverse strand
           @value = map { $_->[0] }
                    sort { $b->[1] <=> $a->[1] }
                    map { [ $_, $_->start] }
                    @value;

           unless ($self->cds) {
               $self->cds(Bio::Location::Simple->new
                          (-start   => $value[-1]->start,
                           -end     => $value[0]->end,
                           -strand  => $value[0]->strand,
                           -seq_id  => $value[0]->seq_id,
                           -verbose => $self->verbose,
                           )
                          );
           }
       } else {               # undef or forward strand
           @value = map { $_->[0] }
                    sort { $a->[1] <=> $b->[1] }
                    map { [ $_, $_->start] }
                    @value;
           unless ($self->cds) {
               $self->cds(Bio::Location::Simple->new
                          (-start   => $value[0]->start,
                           -end     => $value[-1]->end,
                           -strand  => $value[0]->strand,
                           -seq_id  => $value[0]->seq_id,
                           -verbose => $self->verbose,
                           )
                         );
           }

       }

       $self->{'_chr_exons'} = \@value;

       # transform exons from chromosome to gene coordinates
       # but only if gene coordinate system has been set
       my @exons ;
       #my $gene_mapper = $self->$COORDINATE_SYSTEMS{'chr'}. "-". $COORDINATE_SYSTEMS{'gene'};
       my $gene_mapper = "1-2";
       if (defined $self->{'_mappers'}->{$gene_mapper} ) {

           my $tmp_in = $self->{'_in'};
           my $tmp_out = $self->{'_out'};
           my $tmp_verb = $self->verbose;
           $self->verbose(0);

           $self->in('chr');
           $self->out('gene');
           @exons = map {$self->map($_) } @value;

           $self->{'_in'} = ($tmp_in);
           $self->{'_out'} = ($tmp_out);
           $self->verbose($tmp_verb);
       } else {
           @exons = @value;
       }

       my $cds_map = Bio::Coordinate::Collection->new;
       my $inex_map = Bio::Coordinate::Collection->new;
       my $exon_map = Bio::Coordinate::Collection->new;
       my $exon_cds_map = Bio::Coordinate::Collection->new;
       my $intron_map = Bio::Coordinate::Collection->new;
       my $negative_intron_map = Bio::Coordinate::Collection->new;

       my $tr_end = 0;
       my $coffset;
       my $exon_counter;
       my $prev_exon_end;

       for my $exon ( @exons ) {
           $exon_counter++;

           #
           # gene -> cds
           #

           my $match1 = Bio::Location::Simple->new
               (-seq_id =>'gene' ,
                -start  => $exon->start,
                -end    => $exon->end,
                -strand => 1,
                -verbose=> $self->verbose);

           my $match2 = Bio::Location::Simple->new
               (-seq_id => 'cds',
                -start => $tr_end + 1,
                -end => $tr_end + $exon->end - $exon->start +1,
                -strand=>$exon->strand,
                -verbose=>$self->verbose);

           $cds_map->add_mapper(Bio::Coordinate::Pair->new
                                (-in => $match1,
                                 -out => $match2,
                                )
                               );

           if ($exon->start <= 1 and $exon->end >= 1) {
               $coffset = $tr_end - $exon->start + 1;
           }
           $tr_end = $tr_end  + $exon->end - $exon->start + 1;

           #
           # gene -> intron
           #

           if (defined $prev_exon_end) {
               my $match3 = Bio::Location::Simple->new
                   (-seq_id  => 'gene',
                    -start   => $prev_exon_end + 1,
                    -end     => $exon->start -1,
                    -strand  => $exon->strand,
                    -verbose => $self->verbose);

               my $match4 = Bio::Location::Simple->new
                   (-seq_id  => 'intron'. ($exon_counter -1),
                    -start   => 1,
                    -end     => $exon->start - 1 - $prev_exon_end,
                    -strand  =>$exon->strand,
                    -verbose => $self->verbose,);

               # negative intron coordinates
               my $match5 = Bio::Location::Simple->new
                   (-seq_id  => 'intron'. ($exon_counter -1),
                    -start   => -1 * ($exon->start - 2 - $prev_exon_end) -1,
                    -end     => -1,
                    -strand  => $exon->strand,
                    -verbose => $self->verbose);

               $inex_map->add_mapper(Bio::Coordinate::Pair->new
                                     (-in => $match3,
                                      -out => $match4
                                     )
                                    );
               $intron_map->add_mapper(Bio::Coordinate::Pair->new
                                       (-in => $self->_clone_loc($match3),
                                        -out => $self->_clone_loc($match4)
                                       )
                                      );
               $negative_intron_map->add_mapper(Bio::Coordinate::Pair->new
                                                (-in => $self->_clone_loc($match4),
                                                 -out => $match5
                                                ));

           }

           # store the value
           $prev_exon_end = $exon->end;

           #
           # gene -> exon
           #
           my $match6 = Bio::Location::Simple->new
               (-seq_id => 'exon'. $exon_counter,
                -start  => 1,
                -end    => $exon->end - $exon->start +1,
                -strand => $exon->strand,
                -verbose=> $self->verbose,);

           my $pair2 = Bio::Coordinate::Pair->new(-in => $self->_clone_loc($match1),
                                                  -out => $match6
                                                 );
           my $pair3 = Bio::Coordinate::Pair->new(-in => $self->_clone_loc($match6),
                                                  -out => $self->_clone_loc($match2)
                                                 );
           $inex_map->add_mapper(Bio::Coordinate::Pair->new
                                 (-in => $self->_clone_loc($match1),
                                  -out => $match6
                                 )
                                );
           $exon_map->add_mapper(Bio::Coordinate::Pair->new
                                 (-in => $self->_clone_loc($match1),
                                  -out => $self->_clone_loc($match6)
                                 )
                                );
           $exon_cds_map->add_mapper(Bio::Coordinate::Pair->new
                                     (-in => $self->_clone_loc($match6),
                                      -out => $self->_clone_loc($match2)
                                     )
                                    );

       }

       # move coordinate start if exons have negative values
       if ($coffset) {
           foreach my $m ($cds_map->each_mapper) {
               $m->out->start($m->out->start - $coffset);
               $m->out->end($m->out->end - $coffset);
           }

       }

       $self->{'_mappers'}->{$cds_mapper} = $cds_map;
       $self->{'_mappers'}->{$exon_cds_mapper} = $exon_cds_map;
       $self->{'_mappers'}->{$inex_mapper} = $inex_map;
       $self->{'_mappers'}->{$exon_mapper} = $exon_map;
       $self->{'_mappers'}->{$intron_mapper} = $intron_map;
       $self->{'_mappers'}->{$negative_intron_mapper} = $negative_intron_map;
   }
   return  @{$self->{'_chr_exons'}}  || 0;
}

=head2 _clone_loc

 Title   : _clone_loc
 Usage   : $copy_of_loc = $obj->_clone_loc($loc);
 Function: Make a deep copy of a simple location
 Returns : a Bio::Location::Simple object
 Args    : a Bio::Location::Simple object to be cloned

=cut

sub _clone_loc { # clone a simple location
   my ($self,$loc) = @_;

   $self->throw("I need a Bio::Location::Simple , not [". ref $loc. "]")
       unless $loc->isa('Bio::Location::Simple');

   return  Bio::Location::Simple->new
       (-verbose       => $self->verbose,
        -seq_id        => $loc->seq_id,
        -start         => $loc->start,
        -end           => $loc->end,
        -strand        => $loc->strand,
        -location_type => $loc->location_type
       );
}

=head2 cds

 Title   : cds
 Usage   : $obj->cds(20);
 Function: Set and read the offset of CDS from the start of transcipt

           Simple input can be an integer which gives the start of the
           coding region in genomic coordinate. If you want to provide
           the end of the coding region or indicate the use of the
           opposite strand, you have to pass a Bio::RangeI
           (e.g. Bio::Location::Simple or Bio::SegFeature::Generic)
           object to this method.

 Returns : set value or 0
 Args    : new value (optional)

=cut

sub cds {
   my ($self,$value) = @_;
   if( defined $value) {
       if ($value =~ /^[+-]?\d+$/ ) {
           my $loc = Bio::Location::Simple->new(-start=>$value, -end => $value,
                                                -verbose=>$self->verbose);
           $self->{'_cds'} = $loc;
       }
       elsif (ref $value &&  $value->isa('Bio::RangeI') ) {
           $self->{'_cds'} = $value;
       } else {
           $self->throw("I need an integer or Bio::RangeI, not [$value]")
       }
       # strand !!
       my $len;

       $len = $self->{'_cds'}->end - $self->{'_cds'}->start +1
           if defined $self->{'_cds'}->end;

       my $a = $self->_create_pair
           ('chr', 'gene', 0,
            $self->{'_cds'}->start-1,
            $len,
            $self->{'_cds'}->strand);
       my $mapper =  $COORDINATE_SYSTEMS{'chr'}. "-". $COORDINATE_SYSTEMS{'gene'};
       $self->{'_mappers'}->{$mapper} = $a;

       # recalculate exon-based mappers
       if ( defined $self->{'_chr_exons'} ) {
           $self->exons(@{$self->{'_chr_exons'}});
       }

   }
   return $self->{'_cds'} || 0;
}

=head2 map

 Title   : map
 Usage   : $newpos = $obj->map(5);
 Function: Map the location from the input coordinate system
           to a new value in the output coordinate system.
 Example :
 Returns : new value in the output coordiante system
 Args    : a Bio::Location::Simple

=cut

sub map {
   my ($self,$value) = @_;
   my ($res);
   $self->throw("Need to pass me a Bio::Location::Simple or ".
                "Bio::Location::Simple or Bio::SeqFeatureI, not [".
                ref($value). "]")
       unless ref($value) && ($value->isa('Bio::Location::Simple') or
                              $value->isa('Bio::Location::SplitLocationI') or
                              $value->isa('Bio::SeqFeatureI'));
   $self->throw("Input coordinate system not set")
       unless $self->{'_in'};
   $self->throw("Output coordinate system not set")
       unless $self->{'_out'};
   $self->throw("Do not be silly. Input and output coordinate ".
                "systems are the same!")
       unless $self->{'_in'} != $self->{'_out'};

   $self->_check_direction();

   $value = $value->location if $value->isa('Bio::SeqFeatureI');
   $self->debug( "=== Start location: ". $value->start. ",".
                 $value->end. " (". ($value->strand || ''). ")\n");

   # if nozero coordinate system is used in the input values
   if ( defined $self->{'_nozero'} &&
        ( $self->{'_nozero'} == 1 || $self->{'_nozero'} == 3 ) ) {
       $value->start($value->start + 1)
           if defined $value->start && $value->start < 1;
       $value->end($value->end + 1)
           if defined $value->end && $value->end < 1;
   }

   my @steps = $self->_get_path();
   $self->debug( "mapping ". $self->{'_in'}. "->". $self->{'_out'}.
                 "  Mappers: ". join(", ", @steps). "\n");

   foreach my $mapper (@steps) {
       if ($mapper eq $TRANSLATION) {
           if ($self->direction == 1) {

               $value = $self->_translate($value);
               $self->debug( "+   $TRANSLATION cds -> propeptide (translate) \n");
           } else {
               $value = $self->_reverse_translate($value);
               $self->debug("+   $TRANSLATION propeptide -> cds (reverse translate) \n");
           }
       }
       # keep the start and end values, and go on to next iteration
       #  if this mapper is not set
       elsif ( ! defined $self->{'_mappers'}->{$mapper} ) {
           # update mapper name
           $mapper =~ /\d+-(\d+)/;   my ($counter) = $1;
           $value->seq_id($COORDINATE_INTS{$counter});
           $self->debug( "-   $mapper\n");
       } else {
           #
           # the DEFAULT : generic mapping
           #

           $value = $self->{'_mappers'}->{$mapper}->map($value);

           $value->purge_gaps
               if ($value && $value->isa('Bio::Location::SplitLocationI') &&
                   $value->can('gap'));

           $self->debug( "+  $mapper (". $self->direction. "):  start ".
                         $value->start. " end ". $value->end. "\n")
               if $value && $self->verbose > 0;
       }
   }

   # if nozero coordinate system is asked to be used in the output values
   if ( defined $value && defined $self->{'_nozero'} &&
        ( $self->{'_nozero'} == 2 || $self->{'_nozero'} == 3 ) ) {

       $value->start($value->start - 1)
           if defined $value->start && $value->start < 1;
       $value->end($value->end - 1)
           if defined $value->end && $value->end < 1;
   }

   # handle merging of adjacent split locations!

   if (ref $value eq "Bio::Coordinate::Result" && $value->each_match > 1 ) {
       my $prevloc;
       my $merging = 0;
       my $newvalue;
       my @matches;
       foreach my $loc ( $value->each_Location(1) ) {
           unless ($prevloc) {
               $prevloc = $loc;
               push @matches, $prevloc;
               next;
           }
           if ($prevloc->end == ($loc->start - 1) &&
               $prevloc->seq_id eq $loc->seq_id) {
               $prevloc->end($loc->end);
               $merging = 1;
           } else {
               push @matches, $loc;
               $prevloc = $loc;
           }
       }
       if ($merging) {
           if (@matches > 1 ) {
               $newvalue = Bio::Coordinate::Result->new;
               map {$newvalue->add_sub_Location} @matches;
           } else {
               $newvalue = Bio::Coordinate::Result::Match->new
                   (-seq_id   => $matches[0]->seq_id,
                    -start    => $matches[0]->start,
                    -end      => $matches[0]->end,
                    -strand   => $matches[0]->strand,
                    -verbose  => $self->verbose,);
           }
           $value = $newvalue;
       }
   }
   elsif (ref $value eq "Bio::Coordinate::Result" &&
          $value->each_match == 1 ){
       $value = $value->match;
   }

   return $value;
}

=head2 direction

 Title   : direction
 Usage   : $obj->direction('peptide');
 Function: Read-only method for the direction of mapping deduced from
           predefined input and output coordinate names.
 Example :
 Returns : 1 or -1, mapping direction
 Args    : new value (optional)

=cut

sub direction {
   my ($self) = @_;
   return $self->{'_direction'};
}

=head2 swap

 Title   : swap
 Usage   : $obj->swap;
 Function: Swap the direction of transformation
           (input <-> output)
 Example :
 Returns : 1
 Args    :

=cut

sub swap {
   my ($self,$value) = @_;

   ($self->{'_in'}, $self->{'_out'}) = ($self->{'_out'}, $self->{'_in'});
   map { $self->{'_mappers'}->{$_}->swap } keys %{$self->{'_mappers'}};

   # record the changed direction;
   $self->{_direction} *= -1;

   return 1;
}

=head2 to_string

 Title   : to_string
 Usage   : $newpos = $obj->to_string(5);
 Function: Dump the internal mapper values into a human readable format
 Example :
 Returns : string
 Args    :

=cut

sub to_string {
   my ($self) = shift;

   print "-" x 40, "\n";

   # chr-gene
   my $mapper_str = 'chr-gene';
   my $mapper = $self->_mapper_string2code($mapper_str);

   printf "\n     %-12s (%s)\n", $mapper_str, $mapper ;
   if (defined $self->cds) {
       my $end = $self->cds->end -1 if defined $self->cds->end;
       printf "%16s%s: %s (%s)\n", ' ', 'gene offset', $self->cds->start-1 , $end || '';
       printf "%16s%s: %s\n", ' ', 'gene strand', $self->cds->strand || 0;
   }

   # gene-intron
   $mapper_str = 'gene-intron';
   $mapper = $self->_mapper_string2code($mapper_str);
   printf "\n     %-12s (%s)\n", $mapper_str, $mapper ;

   my $i = 1;
   foreach my $pair ( $self->{'_mappers'}->{$mapper}->each_mapper ) {
       printf "%8s :%8s -> %-12s\n", $i, $pair->in->start, $pair->out->start ;
       printf "%8s :%8s -> %-12s\n", '', $pair->in->end, $pair->out->end ;
       $i++;
   }

   # intron-negative_intron
   $mapper_str = 'intron-negative_intron';
   $mapper = $self->_mapper_string2code($mapper_str);
   printf "\n     %-12s (%s)\n", $mapper_str, $mapper ;

   $i = 1;
   foreach my $pair ( $self->{'_mappers'}->{$mapper}->each_mapper ) {
       printf "%8s :%8s -> %-12s\n", $i, $pair->in->start, $pair->out->start ;
       printf "%8s :%8s -> %-12s\n", '', $pair->in->end, $pair->out->end ;
       $i++;
   }

   # gene-exon
   $mapper_str = 'gene-exon';
   $mapper = $self->_mapper_string2code($mapper_str);
   printf "\n     %-12s (%s)\n", $mapper_str, $mapper ;

   $i = 1;
   foreach my $pair ( $self->{'_mappers'}->{$mapper}->each_mapper ) {
       printf "%8s :%8s -> %-12s\n", $i, $pair->in->start, $pair->out->start ;
       printf "%8s :%8s -> %-12s\n", '', $pair->in->end, $pair->out->end ;
       $i++;
   }

   # gene-cds
   $mapper_str = 'gene-cds';
   $mapper = $self->_mapper_string2code($mapper_str);
   printf "\n     %-12s (%s)\n", $mapper_str, $mapper ;

   $i = 1;
   foreach my $pair ( $self->{'_mappers'}->{$mapper}->each_mapper ) {
       printf "%8s :%8s -> %-12s\n", $i, $pair->in->start, $pair->out->start ;
       printf "%8s :%8s -> %-12s\n", '', $pair->in->end, $pair->out->end ;
       $i++;
   }

   # cds-propeptide
   $mapper_str = 'cds-propeptide';
   $mapper = $self->_mapper_string2code($mapper_str);
   printf "\n     %-12s (%s)\n", $mapper_str, $mapper ;
   printf "%9s%-12s\n", "", '"translate"';

   # propeptide-peptide
   $mapper_str = 'propeptide-peptide';
   $mapper = $self->_mapper_string2code($mapper_str);
   printf "\n     %-12s (%s)\n", $mapper_str, $mapper ;
   printf "%16s%s: %s\n", ' ', "peptide offset", $self->peptide_offset;

   print "\nin : ", $self->in, "\n";
   print "out: ", $self->out, "\n";
   my $dir;
   $self->direction ? ($dir='forward') : ($dir='reverse');
   printf "direction: %-8s(%s)\n",  $dir, $self->direction;
   print "\n", "-" x 40, "\n";

   1;
}

=head2 _mapper_code2string
=cut

sub _mapper_code2string {
    my ($self, $code) = @_;
    my ($a, $b) = $code =~ /(\d+)-(\d+)/;
    return $COORDINATE_INTS{$a}. '-'.  $COORDINATE_INTS{$b};

}

=head2 _mapper_string2code
=cut

sub _mapper_string2code {
    my ($self, $string) =@_;
    my ($a, $b) = $string =~ /([^-]+)-(.*)/;
    return $COORDINATE_SYSTEMS{$a}. '-'.  $COORDINATE_SYSTEMS{$b};
}

=head2 _create_pair

 Title   : _create_pair
 Usage   : $mapper = $obj->_create_pair('chr', 'gene', 0, 2555, 10000, -1);
 Function: Internal helper method to create a mapper between
           two coordinate systems
 Returns : a Bio::Coordinate::Pair object
 Args    : string, input coordinate system name,
           string, output coordinate system name,
           boolean, strict mapping
           positive integer, offset
           positive integer, length
           1 || -1 , strand

=cut

sub _create_pair {
   my ($self, $in, $out, $strict, $offset, $length, $strand ) = @_;
   $strict ||= 0;
   $strand ||= 1;
   $length ||= 20;

   my $match1 = Bio::Location::Simple->new
       (-seq_id  => $in,
        -start   => $offset+1,
        -end     => $offset+$length,
        -strand  => 1,
        -verbose => $self->verbose);

   my $match2 = Bio::Location::Simple->new
       (-seq_id  => $out,
        -start   => 1,
        -end     => $length,
        -strand  => $strand,
        -verbose => $self->verbose);

   my $pair = Bio::Coordinate::ExtrapolatingPair->new
       (-in      => $match1,
        -out     => $match2,
        -strict  => $strict,
        -verbose => $self->verbose,
       );

   return $pair;
}

=head2 _translate

 Title   : _translate
 Usage   : $newpos = $obj->_translate($loc);
 Function: Translate the location from the CDS coordinate system
           to a new value in the propeptide coordinate system.
 Example :
 Returns : new location
 Args    : a Bio::Location::Simple or Bio::Location::SplitLocationI

=cut

sub _translate {
   my ($self,$value) = @_;

   $self->throw("Need to pass me a Bio::Location::Simple or ".
                "Bio::Location::SplitLocationI, not [". ref($value). "]")
       unless defined $value &&
           ($value->isa('Bio::Location::Simple') || $value->isa('Bio::Location::SplitLocationI'));

   my $seqid = 'propeptide';

   if ($value->isa("Bio::Location::SplitLocationI") ) {
       my $split = Bio::Location::Split->new(-seq_id=>$seqid);
       foreach my $loc ( $value->each_Location(1) ) {
           my $match = Bio::Location::Simple->new
               (-start   => int ($loc->start / 3 ) +1,
                -end     => int ($loc->end / 3 ) +1,
                -seq_id  => $seqid,
                -strand  => 1,
                -verbose => $self->verbose,
                );
           $split->add_sub_Location($match);
       }
       return $split;

   } else {
       return new Bio::Location::Simple(-start  => int($value->start / 3 )+1,
                                        -end    => int($value->end / 3 )+1,
                                        -seq_id => $seqid,
                                        -strand => 1,
                                        -verbose=> $self->verbose,
                                       );
   }
}

=head2 _frame
=cut

sub _frame {
   my ($self,$value) = @_;

   $self->throw("Need to pass me a Bio::Location::Simple or ".
                "Bio::Location::SplitLocationI, not [". ref($value). "]")
       unless defined $value &&
           ($value->isa('Bio::Location::Simple') || $value->isa('Bio::Location::SplitLocationI'));

   my $seqid = 'propeptide';

   if ($value->isa("Bio::Location::SplitLocationI")) {
       my $split = Bio::Location::Split->new(-seq_id=>$seqid);
       foreach my $loc ( $value->each_Location(1) ) {

           my $match = Bio::Location::Simple->new
               (-start  => ($value->start-1) % 3 +1,
                -end    => ($value->end-1) % 3 +1,
                -seq_id => 'frame',
                -strand => 1,
                -verbose=> $self->verbose);
           $split->add_sub_Location($match);
       }
       return $split;
   } else {
       return new Bio::Location::Simple(-start   => ($value->start-1) % 3 +1,
                                        -end     => ($value->end-1) % 3 +1,
                                        -seq_id  => 'frame',
                                        -strand  => 1,
                                        -verbose => $self->verbose,
                                        );
   }
}

=head2 _reverse_translate

 Title   : _reverse_translate
 Usage   : $newpos = $obj->_reverse_translate(5);
 Function: Reverse translate the location from the propeptide
           coordinate system to a new value in the CSD.
           Note that a single peptide location expands to cover
           the codon triplet
 Example :
 Returns : new location in the CDS coordinate system
 Args    : a Bio::Location::Simple or Bio::Location::SplitLocationI

=cut

sub _reverse_translate {
   my ($self,$value) = @_;

   $self->throw("Need to pass me a Bio::Location::Simple or ".
                "Bio::Location::SplitLocationI, not [". ref($value). "]")
       unless defined $value &&
           ($value->isa('Bio::Location::Simple') || $value->isa('Bio::Location::SplitLocationI'));

   my $seqid = 'cds';

   if ($value->isa("Bio::Location::SplitLocationI")) {
       my $split = Bio::Location::Split->new(-seq_id=>$seqid);
       foreach my $loc ( $value->each_Location(1) ) {

           my $match = Bio::Location::Simple->new
               (-start   => $value->start * 3 - 2,
                -end     => $value->end * 3,
                -seq_id  => $seqid,
                -strand  => 1,
                -verbose => $self->verbose,
                );
           $split->add_sub_Location($match);
       }
       return $split;

   } else {
       return new Bio::Location::Simple(-start   => $value->start * 3 - 2,
                                        -end     => $value->end * 3,
                                        -seq_id  => $seqid,
                                        -strand  => 1,
                                        -verbose => $self->verbose,
                                       );
   }
}

=head2 _check_direction

 Title   : _check_direction
 Usage   : $obj->_check_direction();
 Function: Check and swap when needed the direction the location
           mapping Pairs based on input and output values
 Example :
 Returns : new location
 Args    : a Bio::Location::Simple

=cut

sub _check_direction {
   my ($self) = @_;

   my $new_direction = 1;
   $new_direction = -1 if $self->{'_in'} > $self->{'_out'};

   unless ($new_direction == $self->{_direction} ) {
       map { $self->{'_mappers'}->{$_}->swap } keys %{$self->{'_mappers'}};
       # record the changed direction;
       $self->{_direction} *= -1;
   }
   1;
}

=head2 _get_path

 Title   : _get_path
 Usage   : $obj->_get_path('peptide');
 Function: internal method for finding that shortest path between
           input and output coordinate systems.
           Calculations and caching are handled by the graph class.
           See L<Bio::Coordinate::Graph>.
 Example :
 Returns : array of the mappers
 Args    : none

=cut

sub _get_path {
   my ($self) = @_;

   my $start = $self->{'_in'} || 0;
   my $end = $self->{'_out'} || 0;

   # note the order
   # always go from smaller to bigger: it  makes caching more efficient
   my $reverse;
   if ($start > $end) {
       ($start, $end) = ($end, $start );
       $reverse++;
   }

   my @mappers;
   if (exists $self->{'_previous_path'} and
       $self->{'_previous_path'} eq "$start$end" ) {
       # use cache
       @mappers = @{$self->{'_mapper_path'}};
   } else {
       my $mapper;
       my $prev_node = '';
       @mappers =
           map { $mapper = "$prev_node-$_"; $prev_node = $_; $mapper; }
               $self->{'_graph'}->shortest_path($start, $end);
       shift @mappers;

       $self->{'_previous_path'} = "$start$end";
       $self->{'_mapper_path'} = \@mappers;
   }

   $reverse ? return reverse @mappers : return @mappers;
}

1;
