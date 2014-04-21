#
# BioPerl module for Bio::Assembly::Tools::ContigSpectrum
#
# Copyright by Florent Angly
#
# You may distribute this module under the same terms as Perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Assembly::Tools::ContigSpectrum - create and manipulate contig spectra

=head1 SYNOPSIS

  # Simple contig spectrum creation
  my $csp1 = Bio::Assembly::Tools::ContigSpectrum->new(
    -id       => 'csp1',
    -spectrum => { 1 => 10,
                   2 => 2,
                   3 => 1 } );

  # ...or another way to create a simple contig spectrum
  my $csp2 = Bio::Assembly::Tools::ContigSpectrum->new;
  $csp2->id('csp2');
  $csp2->spectrum({ 1 => 20, 2 => 1, 4 => 1 });

  # Get some information
  print "This is contig spectrum ".$csp->id."\n";
  print "It contains ".$csp->nof_seq." sequences\n";
  print "The largest contig has ".$csp->max_size." sequences\n";
  print "The spectrum is: ".$csp->to_string($csp->spectrum)."\n";

  # Let's add the contig spectra
  my $summed_csp = Bio::Assembly::Tools::ContigSpectrum->new;
  $summed_csp->add($csp1);
  $summed_csp->add($csp2);
  print "The summed contig spectrum is ".$summed_csp->to_string."\n";

  # Make an average
  my $avg_csp = Bio::Assembly::Tools::ContigSpectrum->new;
  $avg_csp = $avg_csp->average([$csp1, $csp2]);
  print "The average contig spectrum is ".$avg_csp->to_string."\n";

  # Get a contig spectrum from an assembly
  my $from_assembly = Bio::Assembly::Tools::ContigSpectrum->new(
    -assembly       => $assembly_object,
    -eff_asm_params => 1);
  print "The contig spectrum from assembly is ".$from_assembly->to_string."\n";

  # Report advanced information (possible because eff_asm_params = 1)
  print "Average sequence length: ".$from_assembly->avg_seq_len." bp\n";
  print "Minimum overlap length: ".$from_assembly->min_overlap." bp\n";
  print "Average overlap length: ".$from_assembly->avg_overlap." bp\n";
  print "Minimum overlap match: ".$from_assembly->min_identity." %\n";
  print "Average overlap match: ".$from_assembly->avg_identity." %\n";

  # Assuming the assembly object contains sequences from several different
  # metagenomes, we have a mixed contig spectrum from which a cross contig
  # spectrum and dissolved contig spectra can be obtained
  my $mixed_csp = $from_assembly;

  # Calculate a dissolved contig spectrum
  my $meta1_dissolved = Bio::Assembly::Tools::ContigSpectrum->new(
    -dissolve => [$mixed_csp, 'metagenome1'] );
  my $meta2_dissolved = Bio::Assembly::Tools::ContigSpectrum->new(
    -dissolve => [$mixed_csp, 'metagenome2'] );
  print "The dissolved contig spectra are:\n".
    $meta1_dissolved->to_string."\n".
    $meta2_dissolved->to_string."\n";

  # Determine a cross contig spectrum
  my $cross_csp = Bio::Assembly::Tools::ContigSpectrum->new(
    -cross => $mixed_csp );
  print "The cross contig spectrum is ".$cross_csp->to_string."\n";

  # Score a contig spectrum (the more abundant the contigs and the larger their
  # size, the larger the score)
  my $csp_score = $csp->score( $csp->nof_seq );

=head1 DESCRIPTION

The Bio::Assembly::Tools::ContigSpectrum Perl module enables to
manually create contig spectra, import them from assemblies,
manipulate them, transform between different types of contig spectra
and output them.

Bio::Assembly::Tools::ContigSpectrum is a module to create, manipulate
and output contig spectra, assembly-derived data used in metagenomics
(community genomics) for diversity estimation.

=head2 Background

A contig spectrum is the count of the number of contigs of different
size in an assembly. For example, the contig spectrum [100 5 1 0 0
...] means that there were 100 singlets (1-contigs), 5 contigs of 2
sequences (2-contigs), 1 contig of 3 sequences (3-contig) and no
larger contigs.

An assembly can be produced from a mixture of sequences from different
metagenomes. The contig obtained from this assembly is a mixed contig
spectrum. The contribution of each metagenome in this mixed contig
spectrum can be obtained by determining a dissolved contig spectrum.

Finally, based on a mixed contig spectrum, a cross contig spectrum can
be determined. In a cross contig spectrum, only contigs containing
sequences from different metagenomes are kept; "pure" contigs are
excluded. Additionally, the total number of singletons (1-contigs)
from each region that assembles with any fragments from other regions
is the number of 1-contigs in the cross contig spectrum.

=head2 Implementation

The simplest representation of a contig spectrum is as a hash
representation where the key is the contig size (number of sequences
making up the contig) and the value the number of contigs of this
size.

In fact, it is useful to have more information associated with the
contig spectrum, hence the Bio::Assembly::Tools::ContigSpectrum module
implements an object containing a contig spectrum hash and additional
information. The get/set methods to access them are:

    id              contig spectrum ID
    nof_rep         number of repetitions (assemblies) used
    max_size        size of (number of sequences in) the largest contig
    spectrum        hash representation of a contig spectrum

    nof_seq         number of sequences
    avg_seq_len     average sequence length

    eff_asm_params  reports effective assembly parameters

    nof_overlaps    number of overlaps (needs eff_asm_params)
    min_overlap     minimum overlap length in a contig (needs eff_asm_params)
    min_identity    minimum sequence identity percentage (needs eff_asm_params)
    avg_overlap     average overlap length (needs eff_asm_params)
    avg_identity    average overlap identity percentage (needs eff_asm_params)

  Operations on the contig spectra:

    to_string       create a string representation of the spectrum
    spectrum        import a hash contig spectrum
    assembly        determine a contig spectrum from an assembly, contig or singlet
    dissolve        calculate a dissolved contig spectrum (depends on assembly)
    cross           produce a cross contig spectrum (depends on assembly)
    add             add a contig spectrum to an existing one
    average         make an average of several contig spectra
    score           score a contig spectrum: the higher the number of contigs
                      and the larger their size, the higher the score.

When using operations that rely on knowing "where" (from what
metagenomes) a sequence came from (i.e. when creating a dissolved or
cross contig spectrum), make sure that the sequences used for the
assembly have a name header, e.g.  E<gt>metagenome1|seq1,
E<gt>metagenome2|seq1, ...

Note: The following operations require the C<Graph::Undirected> module:
   eff_asm_params, cross, dissolve

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

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

Report bugs to the BioPerl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Florent E Angly

Email florent_dot_angly_at_gmail_dot_com

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a "_".

=cut

package Bio::Assembly::Tools::ContigSpectrum;

use strict;

use Bio::Root::Root;
use Bio::Assembly::Scaffold;
use Bio::SimpleAlign;
use Bio::LocatableSeq;

use base 'Bio::Root::Root';


=head2 new

  Title   : new
  Usage   : my $csp = Bio::Assembly::Tools::ContigSpectrum->new();
              or
            my $csp = Bio::Assembly::Tools::ContigSpectrum->new(
              -id => 'some_name',
              -spectrum =>  { 1 => 90 , 2 => 3 , 4 => 1 },
            );
              or
            my $csp = Bio::Assembly::Tools::ContigSpectrum->new(
              -assembly =>  $assembly_obj
            );
  Function: create a new contig spectrum object
  Returns : reference to a contig spectrum object
  Args    : none

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ( $id, $nof_seq, $nof_rep, $max_size, $nof_overlaps, $min_overlap,
    $min_identity, $avg_overlap, $avg_identity, $avg_seq_len, $spectrum,
    $assembly, $eff_asm_params, $dissolve, $cross) = $self->_rearrange(
    [qw(ID NOF_SEQ NOF_REP MAX_SIZE NOF_OVERLAPS MIN_OVERLAP MIN_IDENTITY
    AVG_OVERLAP AVG_IDENTITY AVG_SEQ_LEN SPECTRUM ASSEMBLY EFF_ASM_PARAMS
    DISSOLVE CROSS)], @args );

  # First set up some defauts
  $self->{'_id'}             = 'NoName';
  $self->{'_nof_seq'}        = 0;
  $self->{'_nof_rep'}        = 0;
  $self->{'_max_size'}       = 0;
  $self->{'_nof_overlaps'}   = 0;
  $self->{'_min_overlap'}    = undef;
  $self->{'_min_identity'}   = undef;
  $self->{'_avg_overlap'}    = 0;
  $self->{'_avg_identity'}   = 0;
  $self->{'_avg_seq_len'}    = 0;
  $self->{'_eff_asm_params'} = 0;
  $self->{'_spectrum'}       = {1 => 0}; # contig spectrum hash representation
  $self->{'_assembly'}       = []; # list of assembly, contigs and singlet objects

  # Then, according to user desires, override defaults
  $self->{'_id'}             = $id             if (defined $id);
  $self->{'_nof_seq'}        = $nof_seq        if (defined $nof_seq);
  $self->{'_nof_rep'}        = $nof_rep        if (defined $nof_rep);
  $self->{'_max_size'}       = $max_size       if (defined $max_size);
  $self->{'_nof_overlaps'}   = $nof_overlaps   if (defined $nof_overlaps);
  $self->{'_min_overlap'}    = $min_overlap    if (defined $min_overlap);
  $self->{'_avg_overlap'}    = $avg_overlap    if (defined $avg_overlap);
  $self->{'_min_identity'}   = $min_identity   if (defined $min_identity);
  $self->{'_avg_identity'}   = $avg_identity   if (defined $avg_identity);
  $self->{'_avg_seq_len'}    = $avg_seq_len    if (defined $avg_seq_len);
  $self->{'_eff_asm_params'} = $eff_asm_params if (defined $eff_asm_params);

  # Finally get stuff that can be obtained in an automated way
  $self->_import_spectrum($spectrum) if defined($spectrum);
  $self->_import_assembly($assembly) if defined($assembly);
  $self->_import_cross_csp($cross)   if defined($cross);
  if (defined($dissolve)) {
    my ($mixed_csp, $header) = (@$dissolve[0], @$dissolve[1]);
    $self->_import_dissolved_csp($mixed_csp, $header);
  }

  return $self;
}


=head2 id

  Title   : id
  Usage   : $csp->id
  Function: get/set contig spectrum id
  Returns : string
  Args    : string [optional]

=cut

sub id {
  my ($self, $id) = @_;
  if (defined $id) {
    $self->{'_id'} = $id;
  }
  $id = $self->{'_id'};
  return $id;
}


=head2 nof_seq

  Title   : nof_seq
  Usage   : $csp->nof_seq
  Function: get/set the number of sequences making up the contig spectrum
  Returns : integer
  Args    : integer [optional]

=cut

sub nof_seq {
  my ($self, $nof_seq) = @_;
  if (defined $nof_seq) {
    $self->throw("The number of sequences must be strictly positive. Got ".
      "'$nof_seq'") if $nof_seq < 1;
    $self->{'_nof_seq'} = $nof_seq;
  }
  $nof_seq = $self->{'_nof_seq'};
  return $nof_seq;
}


=head2 nof_rep

  Title   : nof_rep
  Usage   : $csp->nof_rep
  Function: Get/Set the number of repetitions (assemblies) used to create the 
            contig spectrum
  Returns : integer
  Args    : integer [optional]

=cut

sub nof_rep {
  my ($self, $nof_rep) = @_;
  if (defined $nof_rep) {
    $self->throw("The number of repetitions must be strictly positive. Got ".
      "'$nof_rep'") if $nof_rep < 1;
    $self->{'_nof_rep'} = $nof_rep;
  }
  $nof_rep = $self->{'_nof_rep'};
  return $nof_rep;
}


=head2 max_size

  Title   : max_size
  Usage   : $csp->max_size
  Function: get/set the size of (number of sequences in) the largest contig
  Returns : integer
  Args    : integer [optional]

=cut

sub max_size {
  my ($self, $max_size) = @_;
  if (defined $max_size) {
    $self->throw("The contig maximum size must be strictly positive. Got ".
      "'$max_size'") if $max_size < 1;
    $self->{'_max_size'} = $max_size;
  }
  $max_size = $self->{'_max_size'};
  return $max_size;
}


=head2 nof_overlaps

  Title   : nof_overlaps
  Usage   : $csp->nof_overlaps
  Function: Get/Set the number of overlaps in the assembly.
  Returns : integer
  Args    : integer [optional]

=cut

sub nof_overlaps {
  my ($self, $nof_overlaps) = @_;
  if (defined $nof_overlaps) {
    $self->throw("The number of overlaps must be strictly positive. Got ".
      "'$nof_overlaps'") if $nof_overlaps < 1;
    $self->{'_nof_overlaps'} = $nof_overlaps;
  }
  $nof_overlaps = $self->{'_nof_overlaps'};
  return $nof_overlaps;
}


=head2 min_overlap

  Title   : min_overlap
  Usage   : $csp->min_overlap
  Function: get/set the assembly minimum overlap length
  Returns : integer
  Args    : integer [optional]

=cut

sub min_overlap {
  my ($self, $min_overlap) = @_;
  if (defined $min_overlap) {
    $self->throw("The minimum of overlap length must be strictly positive. Got".
      " '$min_overlap'") if $min_overlap < 1;
    $self->{'_min_overlap'} = $min_overlap;
  }
  $min_overlap = $self->{'_min_overlap'};
  return $min_overlap;
}


=head2 avg_overlap

  Title   : avg_overlap
  Usage   : $csp->avg_overlap
  Function: get/set the assembly average overlap length
  Returns : decimal
  Args    : decimal [optional]

=cut

sub avg_overlap {
  my ($self, $avg_overlap) = @_;
  if (defined $avg_overlap) {
    $self->throw("The average overlap length must be strictly positive. Got ".
      "'$avg_overlap'") if $avg_overlap < 1;
    $self->{'_avg_overlap'} = $avg_overlap;
  }
  $avg_overlap = $self->{'_avg_overlap'};
  return $avg_overlap;
}


=head2 min_identity

  Title   : min_identity
  Usage   : $csp->min_identity
  Function: get/set the assembly minimum overlap identity percent
  Returns : 0 < decimal < 100
  Args    : 0 < decimal < 100 [optional]

=cut

sub min_identity {
  my ($self, $min_identity) = @_;
  if (defined $min_identity) {
    $self->throw("The minimum overlap percent identity must be strictly ".
      "positive. Got '$min_identity'") if $min_identity < 1;
    $self->{'_min_identity'} = $min_identity;
  }
  $min_identity = $self->{'_min_identity'};
  return $min_identity;
}


=head2 avg_identity

  Title   : avg_identity
  Usage   : $csp->avg_identity
  Function: get/set the assembly average overlap identity percent
  Returns : 0 < decimal < 100
  Args    : 0 < decimal < 100 [optional]

=cut

sub avg_identity {
  my ($self, $avg_identity) = @_;
  if (defined $avg_identity) {
    $self->throw("The average overlap percent identity must be strictly ".
      "positive. Got '$avg_identity'") if $avg_identity < 1;
    $self->{'_avg_identity'} = $avg_identity;
  }
  $avg_identity = $self->{'_avg_identity'};
  return $avg_identity;
}


=head2 avg_seq_len

  Title   : avg_seq_len
  Usage   : $csp->avg_seq_len
  Function: get/set the assembly average sequence length
  Returns : avg_seq_len
  Args    : real [optional]

=cut

sub avg_seq_len {
  my ($self, $avg_seq_len) = @_;
  if (defined $avg_seq_len) {
    $self->throw("The average sequence length must be strictly positive. Got ".
      "'$avg_seq_len'") if $avg_seq_len < 1;
    $self->{'_avg_seq_len'} = $avg_seq_len;
  }
  $avg_seq_len = $self->{'_avg_seq_len'};
  return $avg_seq_len;
}


=head2 eff_asm_params

  Title   : eff_asm_params
  Usage   : $csp->eff_asm_params(1)
  Function: Get/set the effective assembly parameters option. It defines if the
            effective assembly parameters should be determined when a contig
            spectrum based or derived from an assembly is calculated. The
            effective assembly parameters include avg_seq_length, nof_overlaps,
            min_overlap, avg_overlap, min_identity and avg_identity.
            1 = get them, 0 = don't.
  Returns : integer
  Args    : integer [optional]

=cut

sub eff_asm_params {
  my ($self, $eff_asm_params) = @_;
  if (defined $eff_asm_params) {
    $self->throw("eff_asm_params can only take values 0 or 1. Input value was ".
      "'$eff_asm_params'") unless $eff_asm_params == 0 || $eff_asm_params == 1;
    $self->{'_eff_asm_params'} = $eff_asm_params;
  }
  $eff_asm_params = $self->{'_eff_asm_params'};
  return $eff_asm_params;
}


=head2 spectrum

  Title   : spectrum
  Usage   : my $spectrum = $csp->spectrum({1=>10, 2=>2, 3=>1});
  Function: Get the current contig spectrum represented as a hash / Update a
            contig spectrum object based on a contig spectrum represented as a
            hash
            The hash representation of a contig spectrum is as following:
              key   -> contig size (in number of sequences)
              value -> number of contigs of this size
  Returns : contig spectrum as a hash reference
  Args    : contig spectrum as a hash reference [optional]

=cut

sub spectrum {
  my ($self, $spectrum) = @_;
  if (defined $spectrum) {
    $self->_import_spectrum($spectrum);
  }
  $spectrum = $self->{'_spectrum'};
  return $spectrum;
}


=head2 assembly

  Title   : assembly
  Usage   : my @obj_list = $csp->assembly();
  Function: get/set the contig spectrum object by adding an assembly, contig or
            singlet object to it, or get the list of objects associated with it
  Returns : arrayref of assembly, contig and singlet objects used in the contig
            spectrum object (Bio::Assembly::Scaffold, Bio::Assembly::Contig and
            Bio::Assembly::Singlet objects)
  Args    : Bio::Assembly::Scaffold, Contig or Singlet object

=cut

sub assembly {
  my ($self, $assembly) = @_;
  if (defined $assembly) {
    $self->_import_assembly($assembly);
  }
  my @obj_list = @{$self->{'_assembly'}} if defined $self->{'_assembly'};
  return @obj_list;
}


=head2 drop_assembly

  Title   : drop_assembly
  Usage   : $csp->drop_assembly();
  Function: Remove all assembly objects associated with a contig spectrum.
            Assembly objects can take a lot of memory, which can be freed by
            calling this method. Don't call this method if you need the assembly
            object later on, for example for creating a dissolved or cross
            contig spectrum.
  Returns : 1 for success
  Args    : none

=cut

sub drop_assembly {
  my ($self) = @_;
  $self->{'_assembly'} = [];
  return 1;
}


=head2 dissolve

  Title   : dissolve
  Usage   : $dissolved_csp->dissolve($mixed_csp, $seq_header);
  Function: Dissolve a mixed contig spectrum for the set of sequences that
            contain the specified header, i.e. determine the contribution of
            these sequences to the mixed contig spectrum. The mixed contig
            spectrum object must have one or several assembly object(s). In
            addition, min_overlap, min_identity and eff_asm_params are taken
            from the mixed contig spectrum, unless they are specified manually
            for the dissolved contig spectrum. The dissolved contigs underlying
            the contig spectrum can be obtained by calling the assembly() method.
  Returns : 1 for success
  Args    : Bio::Assembly::Tools::ContigSpectrum reference
            sequence header string

=cut

sub dissolve {
  my ($self, $mixed_csp, $seq_header) = @_;
  $self->_import_dissolved_csp($mixed_csp, $seq_header);
  return 1;
}


=head2 cross

  Title   : cross
  Usage   : $cross_csp->cross($mixed_csp);
  Function: Calculate a cross contig_spectrum based on a mixed contig_spectrum.
            The underlying cross-contigs themselves can be obtained by calling 
            the assembly() method.
  Returns : 1 for success
  Args    : Bio::Assembly::Tools::ContigSpectrum reference

=cut

sub cross {
  my ($self, $mixed_csp) = @_;
  $self->_import_cross_csp($mixed_csp);
  return 1;
}

=head2 to_string

  Title   : to_string
  Usage   : my $csp_string = $csp->to_string;
  Function: Convert the contig spectrum into a string (easy to print!!).
  Returns : string
  Args    : element separator (integer) [optional]
              1 -> space-separated
              2 -> tab-separated
              3 -> newline-separated

=cut

sub to_string {
  my ($self, $element_separator) = @_;
  return 0 if $self->{'_max_size'} == 0;
  $element_separator ||= 1;
  if ($element_separator == 1) {
    $element_separator = ' ';
  } elsif ($element_separator == 2) {
    $element_separator = "\t";
  } elsif ($element_separator == 3) {
    $element_separator = "\n";
  } else {
    $self->throw("Unknown separator type '$element_separator'\n");
  }
  my $str = '';
  for (my $q = 1 ; $q <= $self->{'_max_size'} ; $q++) {
    my $val = 0;
    if (exists $self->{'_spectrum'}{$q}) {
      $val = $self->{'_spectrum'}{$q};
    }
    $str .= $val.$element_separator;
  }
  $str =~ s/\s$//;
  return $str;
}


=head2 add

  Title   : add
  Usage   : $csp->add($additional_csp);
  Function: Add a contig spectrum to an existing one: sums the spectra, update
            the number of sequences, number of repetitions, ...
  Returns : 1 for success
  Args    : Bio::Assembly::Tools::ContigSpectrum object

=cut

sub add {
  my ($self, $csp) = @_;
  # Sanity check
  if( !ref $csp || ! $csp->isa('Bio::Assembly::Tools::ContigSpectrum') ) {
        $self->throw("Unable to process non Bio::Assembly::Tools::ContigSpectrum ".
        "object [".ref($csp)."]");
  }
  # Update overlap statistics
  if ( $self->{'_eff_asm_params'} > 0 ) {

    # Warnings
    if ( $csp->{'_eff_asm_params'} == 0 ) {
      $self->warn("The parent contig spectrum needs effective assembly ".
        "parameters (eff_asm_params = ".$self->{'_eff_asm_params'}.") but the ".
        "child contig spectrum doesn't have them (eff_asm_params = ".
        $csp->{'_eff_asm_params'}."). Skipping them...");
    } elsif ( $csp->{'_eff_asm_params'} != $self->{'_eff_asm_params'} ) {
      $self->warn("The parent contig spectrum needs a different method for ".
        "detecting the effective assembly parameters (eff_asm_params = ".
        $self->{'_eff_asm_params'}.") than the one specified for the child ".
        "contig spectrum (eff_asm_params = ".$csp->{'_eff_asm_params'}."). ".
        "Ignoring the differences...");
    }

    # Update existing stats
    ( $self->{'_avg_overlap'} , $self->{'_avg_identity'}, $self->{'_min_overlap'},
      $self->{'_min_identity'}, $self->{'_nof_overlaps'} ) = $self->_update_overlap_stats(
      $self->{'_avg_overlap'} , $self->{'_avg_identity'}, $self->{'_min_overlap'},
      $self->{'_min_identity'}, $self->{'_nof_overlaps'},
      $csp->{'_avg_overlap'}  , $csp->{'_avg_identity'} , $csp->{'_min_overlap'},
      $csp->{'_min_identity'} , $csp->{'_nof_overlaps'} );

  }

  # Update sequence average length (not number of sequences)
  ( $self->{'_avg_seq_len'} ) = $self->_update_seq_stats(
    $self->{'_avg_seq_len'}, $self->{'_nof_seq'}, $csp->{'_avg_seq_len'},
    $csp->{'_nof_seq'} );

  # Update spectrum (and nof_seq, max_size, and increment nof_rep by 1)
  $self->_import_spectrum($csp->{'_spectrum'});

  # Update nof_rep
  $self->{'_nof_rep'}--;
  $self->{'_nof_rep'} += $csp->{'_nof_rep'};
  # Update list of assembly objects used
  push @{$self->{'_assembly'}}, @{$csp->{'_assembly'}}
    if defined $csp->{'_assembly'};

  return 1;
}


=head2 average

  Title   : average
  Usage   : my $avg_csp = $csp->average([$csp1, $csp2, $csp3]);
  Function: Average one contig spectrum or the sum of several contig spectra by
            the number of repetitions
  Returns : Bio::Assembly::Tools::ContigSpectrum
  Args    : Bio::Assembly::Tools::ContigSpectrum array reference
            eff_asm_params

=cut

sub average {
  my ($self, $list) = @_;
  # Sanity check
  if ( ! ref $list || ! ref $list eq 'ARRAY') {
    $self->throw("Average takes an array reference but got [".ref($list)."]");
  }
  # New average contig spectrum object
  my $avg = Bio::Assembly::Tools::ContigSpectrum->new;
  $avg->{'_eff_asm_params'} = $self->{'_eff_asm_params'};

  # Cycle through contig spectra
  my $tot_nof_rep = 0;
  for my $csp (@$list) {
    # Sanity check
    if (not $csp->isa('Bio::Assembly::Tools::ContigSpectrum')) {
      $csp->throw("Unable to process non Bio::Assembly::Tools::ContigSpectrum ".
        "object [".ref($csp)."]");
    }
    # Import contig spectrum
    $avg->add($csp);
  }

  # Average sum of contig spectra by number of repetitions
  for (my $q = 1 ; $q <= $avg->{'_max_size'} ; $q++) {
    $avg->{'_spectrum'}{$q} /= $avg->{'_nof_rep'}
      if (defined $avg->{'_spectrum'}{$q});
  }
  # Average number of sequences
  $avg->{'_nof_seq'} /= $avg->{'_nof_rep'};
  # Average number of overlaps
  $avg->{'_nof_overlaps'} /= $avg->{'_nof_rep'};

  return $avg;
}


=head2 score

  Title   : score
  Usage   : my $score = $csp->score();
  Function: Score a contig spectrum (or cross-contig spectrum) such that the
             higher the number of contigs (or cross-contigs) and the larger their 
             size, the higher the score.
             Let n   : total number of sequences
                 c_q : number of contigs of size q
                 q   : number of sequence in a contig
             We define: score = n/(n-1) * (X - 1/n)
                  where X = sum ( c_q * q^2 ) / n**2
             The score ranges from 0 (singlets only) to 1 (a single large contig)
             It is possible to specify a value for the number of sequences to
              assume in the contig spectrum.
  Returns : contig score, or undef if there were no sequences in the contig spectrum
  Args    : number of total sequences to assume [optional]

=cut

sub score {
  my ($self, $nof_seqs) = @_;
  # Sanity check
  my $n = $self->nof_seq;
  return undef if ($n <= 0);
  # Calculate X
  my $score = 0;
  my $q_max = $self->max_size;
  my $spec = $self->spectrum;
  for my $q ( 1 .. $q_max ) {
    my $c_q = $spec->{$q};
    if ( $q == 1 && $nof_seqs ) {
      $c_q += $nof_seqs - $n;
      $n = $nof_seqs;
    }
    next if not $c_q;
    $score += $c_q * $q ** 2;
  }
  $score /= $n ** 2;
  # Rescale X to obtain the score
  $score = $n/($n-1) * ($score - 1/$n);
  return $score;
}


=head2 _naive_assembler

  Title   : _naive_assembler
  Usage   : 
  Function: Reassemble the specified sequences only based on their position in
            the contig. This naive assembly only verifies that the minimum
            overlap length and percentage identity are respected. No actual
            alignment is done
  Returns : arrayref of contigs and singlets
  Args    : Bio::Assembly::Contig
            array reference of sequence IDs to use [optional]
            minimum overlap length (integer)       [optional]
            minimum percentage identity (integer)  [optional]

=cut

sub _naive_assembler {
  my ($self, $contig, $seqlist, $min_overlap, $min_identity) = @_;

  # Use all reads if none was specified:
  if (not defined $seqlist) {
    for my $seq ($contig->each_seq) {
      push @$seqlist, $seq->id;
    }
  }

  # Sanity checks
  if ( ! ref $seqlist || ! ref($seqlist) eq 'ARRAY') {
    $self->throw('Expecting an array reference. Got ['.ref($seqlist)."] \n");
  }
  my $max = scalar @$seqlist;
  $self->throw("Expecting at least 2 sequences as input for _naive_assembler")
    if ($max < 2);

  # Build contig graph
  my %seq_hash = map { $_ => undef } (@$seqlist) if (scalar @$seqlist > 0);
  my ($g, $overlaps) = $self->_contig_graph($contig, \%seq_hash, $min_overlap, $min_identity);

  # Construct sub-contigs
  my @contig_objs;
  my $num = 1;
  if (defined $g) {
    for my $connected_reads ($g->connected_components) { # reads that belong in contigs
      my $sub_id = $contig->id.'_'.$num;
      my $sub_contig = $self->_create_subcontig($contig, $connected_reads, $sub_id);
      push @contig_objs, $sub_contig;
      $num++;
      for my $read_id ( @$connected_reads ) {
        delete $seq_hash{$read_id};
      }
    }
  }

  # Construct sub-singlets
  my @singlet_objs;
  for my $read_id ( keys %seq_hash ) {
    my $read = $contig->get_seq_by_name($read_id);
    my $sub_singlet = Bio::Assembly::Singlet->new(
      -id => $contig->id.'_'.$num,
      -seqref => $self->_obj_copy($read)
    );
    $num++;
    push @singlet_objs, $sub_singlet;
  }

  return [@contig_objs, @singlet_objs];
}


=head2 _create_subcontig

  Title   : _create_subcontig
  Usage   : 
  Function: Create a subcontig from another contig
  Returns : Bio::Assembly::Contig object
  Args    : Bio::Assembly::Contig
            arrayref of the IDs of the reads to includes in the subcontig
            ID to give to the subcontig

=cut

sub _create_subcontig {
  my ($self, $contig, $read_ids, $sub_contig_id) = @_;

  my $sub_contig = Bio::Assembly::Contig->new( -id => $sub_contig_id );

  # Get min and max read coordinates
  my ($min, $max) = (undef, undef);
  for my $read_id ( @$read_ids ) {
    my ($aln_coord) = $contig->get_features_collection->get_features_by_type("_aligned_coord:$read_id");
    my $seq_start = $aln_coord->location->start;
    my $seq_end   = $aln_coord->location->end;
    $min = $seq_start if (not defined $min) || ((defined $min) && ($seq_start < $min));
    $max = $seq_end   if (not defined $max) || ((defined $max) && ($seq_end   > $max));
  }

  # Add reads to subcontig
  for my $read_id (@$read_ids) {
    my $read  = $self->_obj_copy($contig->get_seq_by_name($read_id));
    my $coord = $self->_obj_copy($contig->get_seq_coord($read));
    if ($min > 1) {
      # adjust read coordinates
      $coord->start( $coord->start - $min + 1 );
      $coord->end( $coord->end - $min + 1 );
    }
    $sub_contig->set_seq_coord($coord, $read);
  }

  # Truncate copy of original consensus to new boundaries
  my $cons_seq  = $contig->get_consensus_sequence;
  $sub_contig->set_consensus_sequence( $self->_obj_copy($cons_seq, $min, $max) );
  my $cons_qual = $contig->get_consensus_quality;
  if ($cons_qual) {
     $sub_contig->set_consensus_quality( $self->_obj_copy($cons_qual, $min, $max) );
  }

  return $sub_contig;
}

=head2 _obj_copy

  Title   : _obj_copy
  Usage   : 
  Function: Copy (most of) an object, and optionally truncate it
  Returns : another a Bio::LocatableSeq, Bio::Seq::PrimaryQual, or
              Bio::SeqFeature::Generic object
  Args    : a Bio::LocatableSeq, Bio::Seq::PrimaryQual, or
              Bio::SeqFeature::Generic object
            a start position
            an end position

=cut

sub _obj_copy {
  my ($self, $obj, $start, $end) = @_;
  my $new;
  if ($obj->isa('Bio::Seq::PrimaryQual')) {
    my $qual = [@{$obj->qual}]; # copy of the quality scores
    if (defined $start && defined $end && $start !=1 && $end != scalar(@$qual)) {
      # Truncate the quality scores
      $qual = [ splice @$qual, $start - 1, $end - $start + 1 ];
    }
    $new = Bio::Seq::PrimaryQual->new(
      -qual   => $qual,
      -id     => $obj->id,
    );

  } elsif ($obj->isa('Bio::LocatableSeq')) {
    my $seq = $obj->seq;
    if (defined $start && defined $end && $start != 1 && $end != length($seq)) {
      # Truncate the aligned sequence
      $seq = substr $seq, $start - 1, $end - $start + 1;
    }
    $new = Bio::LocatableSeq->new(
      -seq      => $seq,
      -id       => $obj->id,
      -start    => $obj->start,
      -strand   => $obj->strand,
      -alphabet => $obj->alphabet,
    );
 
  } elsif ($obj->isa('Bio::SeqFeature::Generic')) {
    $new = Bio::SeqFeature::Generic->new(
      -start  => $obj->start,
      -end    => $obj->end
    );
  }

  return $new;
}


=head2 _new_from_assembly

  Title   : _new_from_assembly
  Usage   : 
  Function: Creates a new contig spectrum object based solely on the result of 
            an assembly, contig or singlet
  Returns : Bio::Assembly::Tools::ContigSpectrum object
  Args    : Bio::Assembly::Scaffold, Contig or Singlet object

=cut

sub _new_from_assembly {
  # Create new contig spectrum object based purely on what we can get from the
  # assembly object
  my ($self, $assemblyobj) = @_;
  my $csp = Bio::Assembly::Tools::ContigSpectrum->new();
  # 1: Set id
  $csp->{'_id'} = $assemblyobj->id;
  # 2: Set overlap statistics: nof_overlaps, min_overlap, avg_overlap,
  #  min_identity and avg_identity
  $csp->{'_eff_asm_params'} = $self->{'_eff_asm_params'};
  $csp->{'_min_overlap'}    = $self->{'_min_overlap'};
  $csp->{'_min_identity'}   = $self->{'_min_identity'};
  if ( $csp->{'_eff_asm_params'} > 0 ) {
    ( $csp->{'_avg_overlap'}, $csp->{'_avg_identity'}, $csp->{'_min_overlap'},
      $csp->{'_min_identity'}, $csp->{'_nof_overlaps'} )
      = $csp->_get_assembly_overlap_stats($assemblyobj);
  }
  # 3: Set sequence statistics: nof_seq and avg_seq_len
  ($csp->{'_avg_seq_len'}, $csp->{'_nof_seq'}) = $self->_get_assembly_seq_stats($assemblyobj);
  ### any use in using _naive_assembler here to re-assemble with specific minmum criteria?
  # 4: Set the spectrum: spectrum and max_size
  for my $contigobj ( $self->_get_contig_like($assemblyobj) ) {
    my $size = $contigobj->num_sequences;
    if (defined $csp->{'_spectrum'}{$size}) {
      $csp->{'_spectrum'}{$size}++;
    } else {
      $csp->{'_spectrum'}{$size} = 1;
    }
    $csp->{'_max_size'} = $size if $size > $csp->{'_max_size'};
  }
  # 5: Set list of assembly objects used
  push @{$csp->{'_assembly'}}, $assemblyobj;
  # 6: Set number of repetitions
  $csp->{'_nof_rep'} = 1;
  return $csp;
}


=head2 _new_dissolved_csp

  Title   : _new_dissolved_csp
  Usage   : 
  Function: create a dissolved contig spectrum object
  Returns : dissolved contig spectrum
  Args    : mixed contig spectrum
            header of sequences to keep in this contig spectrum

=cut

sub _new_dissolved_csp {
  my ($self, $mixed_csp, $seq_header) = @_;
  # Sanity checks on the mixed contig spectrum

  # min_overlap and min_identity must be specified if there are some overlaps
  # in the mixed contig
  if ($mixed_csp->{'_nof_overlaps'} > 0) {
    unless ( defined $self->{'_min_overlap'}      || 
             defined $mixed_csp->{'_min_overlap'}  ) {
      $self->throw("min_overlap must be defined in the dissolved contig spectrum".
        " or mixed contig spectrum to dissolve a contig");
    }
    unless ( defined $self->{'_min_identity'}      ||
             defined $mixed_csp->{'_min_identity'}  ) {
      $self->throw("min_identity must be defined in the dissolved contig spectrum".
        " or mixed contig spectrum");
    }
  }

  # there must be at least one assembly in mixed contig spectrum
  if (!defined $mixed_csp->{'_assembly'} ||
      scalar @{$mixed_csp->{'_assembly'}} < 1) {
    $self->throw("The mixed contig spectrum must be based on at least one
    assembly");
  }

  # New dissolved contig spectrum object
  my $dissolved = Bio::Assembly::Tools::ContigSpectrum->new();

  # Take attributes of the parent contig spectrum if they exist, or those of the
  # mixed contig spectrum otherwise
  my ($eff_asm_params, $min_overlap, $min_identity);
  if ($self->{'_eff_asm_params'}) {
    $eff_asm_params = $self->{'_eff_asm_params'};
  } else {
    $eff_asm_params = $mixed_csp->{'_eff_asm_params'};
  }
  if ($self->{'_min_overlap'}) {
    $min_overlap = $self->{'_min_overlap'};
  } else {
    $min_overlap = $mixed_csp->{'_min_overlap'};
  }
  if ($self->{'_min_identity'}) {
    $min_identity = $self->{'_min_identity'};
  } else {
    $min_identity = $mixed_csp->{'_min_identity'};
  }
  ($dissolved->{'_eff_asm_params'}, $dissolved->{'_min_overlap'},
    $dissolved->{'_min_identity'}) = ($eff_asm_params, $min_overlap,
    $min_identity);

  # Dissolve each assembly
  for my $obj (@{$mixed_csp->{'_assembly'}}) {
    for my $contig ( $self->_get_contig_like($obj) ) {

      # Dissolve this assembly/contig/singlet for the given sequences
      my $dissolved_objs = $self->_dissolve_contig( $contig, $seq_header,
        $min_overlap, $min_identity );

      # Add dissolved contigs to contig spectrum
      for my $dissolved_obj (@$dissolved_objs) {
        $dissolved->assembly($dissolved_obj);
        $dissolved->{'_nof_rep'}--;
      }

    }
  }

  # Update nof_rep
  $dissolved->{'_nof_rep'} += $mixed_csp->{'_nof_rep'};

  return $dissolved;
}


=head2 _dissolve_contig

  Title   : _dissolve_contig
  Usage   : 
  Function: dissolve a contig
  Returns : arrayref of contigs and singlets
  Args    : mixed contig spectrum
            header of sequences to keep in this contig spectrum
            minimum overlap
            minimum identity

=cut

sub _dissolve_contig {
  my ($self, $contig, $wanted_origin, $min_overlap, $min_identity) = @_;

  # List of reads
  my @seqs;
  if ($contig->isa('Bio::Assembly::Singlet')) {
    @seqs = $contig->seqref;
  } elsif ($contig->isa('Bio::Assembly::Contig')) {
    @seqs = $contig->each_seq;
  }

  # Get sequences from the desired metagenome
  my @contig_seqs;
  for my $seq (@seqs) {
    my $seq_id = $seq->id;
    # get sequence origin
    next unless $self->_seq_origin($seq_id) eq $wanted_origin;
    # add it to hash
    push @contig_seqs, $seq_id;
  }

  # Update spectrum
  my $size = scalar @contig_seqs;
  my $objs;
  if ($size == 1) {
    # create a singlet and add it to list of objects
    my $id  = $contig_seqs[0]; 
    my $seq = $contig->get_seq_by_name($id);
    push @$objs, Bio::Assembly::Singlet->new(-id => $contig->id, -seqref => $self->_obj_copy($seq) );
  } elsif ($size > 1) {
    # Reassemble good sequences
    my $contig_objs = $self->_naive_assembler( $contig, \@contig_seqs,
      $min_overlap, $min_identity );
    push @$objs, @$contig_objs;
  }

  return $objs;
}


=head2 _new_cross_csp

  Title   : _new_cross_csp
  Usage   : 
  Function: create a cross contig spectrum object
  Returns : cross-contig spectrum
  Args    : mixed contig spectrum

=cut

sub _new_cross_csp {
  my ($self, $mixed_csp) = @_;
  # Sanity check on the mixed contig spectrum
  # There must be at least one assembly
  if (!defined $mixed_csp->{'_assembly'} ||
      scalar @{$mixed_csp->{'_assembly'}} < 1) {
    $self->throw("The mixed contig spectrum must be based on at least one ".
    "assembly.");
  }
  
  # New dissolved contig spectrum object
  my $cross = Bio::Assembly::Tools::ContigSpectrum->new();
  
  # Take attributes from parent or from mixed contig spectrums
  my ($eff_asm_params, $min_overlap, $min_identity);
  if ($self->{'_eff_asm_params'}) {
    $eff_asm_params = $self->{'_eff_asm_params'};
  } else {
    $eff_asm_params = $mixed_csp->{'_eff_asm_params'};
  }
  if ($self->{'_min_overlap'}) {
    $min_overlap = $self->{'_min_overlap'};
  } else {
    $min_overlap = $mixed_csp->{'_min_overlap'};
  }
  if ($self->{'_min_identity'}) {
    $min_identity = $self->{'_min_identity'};
  } else {
    $min_identity = $mixed_csp->{'_min_identity'};
  }
  ($cross->{'_eff_asm_params'},$cross->{'_min_overlap'},$cross->{'_min_identity'})
    = ($eff_asm_params, $min_overlap, $min_identity);

  # Get cross contig spectrum for each assembly
  for my $obj ( @{$mixed_csp->{'_assembly'}} ) {
    for my $contig ( $self->_get_contig_like($obj) ) {

      # Go through contigs and skip the pure ones
      my ($cross_contigs, $nof_cross_singlets) = $self->_cross_contig($contig,
        $min_overlap, $min_identity);

      # Add cross-contig
      for my $cross_contig ( @$cross_contigs ) {
        $cross->assembly($cross_contig);
        $cross->{'_nof_rep'}--;
      }

      # Add cross-singlets
      $cross->{'_spectrum'}->{'1'} += $nof_cross_singlets;

    }
  }

  # Update nof_rep
  $cross->{'_nof_rep'} += $mixed_csp->{'_nof_rep'};

  return $cross;
}


=head2 _cross_contig

  Title   : _cross_contig
  Usage   : 
  Function: calculate cross contigs
  Returns : arrayref of cross-contigs
            number of cross-singlets
  Args    : contig
            minimum overlap
            minimum identity

=cut

sub _cross_contig {
  my ($self, $contig, $min_overlap, $min_identity) = @_;

  my $nof_cross_singlets = 0;
  my @cross_contigs;

  # Weed out pure contigs
  my %all_origins;
  for my $seq ($contig->each_seq) {
    my $seq_id = $seq->id;
    my $seq_origin = $self->_seq_origin($seq_id);
    if (not defined $seq_origin) {
      $self->warn("Sequence $seq_id doesn't have any header. Skipping it...");
      next;
    }
    if ( scalar keys %all_origins > 1 ) {
      # a cross-contig spectrum
      last;
    }
    $all_origins{$seq_origin} = undef;
  }
  if ( scalar keys %all_origins <= 1 ) {
    # a pure contig
    return \@cross_contigs, $nof_cross_singlets;
  }
  %all_origins = ();

  # Break the cross-contigs using the specified stringency
  my $test_contigs = $self->_naive_assembler($contig, undef, $min_overlap, $min_identity);

  # Find cross contigs and singlets
  for my $test_contig ( @$test_contigs ) {

    # Find cross-contigs
    my %origins;
    for my $seq ($test_contig->each_seq) {
      my $seq_id = $seq->id;
      my $seq_origin = $self->_seq_origin($seq_id);
      next if not defined $seq_origin;
      push @{$origins{$seq_origin}}, $seq_id;
    }
    if (scalar keys %origins > 1) {
      # Found a cross-contig
      push @cross_contigs, $test_contig;
    } else {
      next;
    }

    # Find cross-singlets
    for my $origin (keys %origins) {
      my @ori_ids = @{$origins{$origin}};
      if (scalar @ori_ids == 1) {
        $nof_cross_singlets++;
      } elsif (scalar @ori_ids > 1) {
        # Dissolve contig for the given origin

        ### consider using the minimum overlap and identity here again?
        my $ori_contigs = $self->_naive_assembler($test_contig, \@ori_ids, undef, undef);

        for my $ori_contig (@$ori_contigs) {
          $nof_cross_singlets++ if $ori_contig->num_sequences == 1;
        }
      }
    }


  }

  return \@cross_contigs, $nof_cross_singlets;
}


=head2 _seq_origin

  Title   : _seq_origin
  Usage   : 
  Function: determines where a sequence comes from using its header. For example
            the origin of the sequence 'metagenome1|gi|9626988|ref|NC_001508.1|'
            is 'metagenome1'
  Returns : origin
  Args    : sequence ID

=cut

sub _seq_origin {
  # Current sequence origin. Example: sequence with ID
  # 'metagenome1|gi|9626988|ref|NC_001508.1|' has header 'metagenome1'
  my ($self, $seq_id) = @_;
  my $origin;
  if ( $seq_id =~ m/^(.+?)\|/ ) {
    $origin = $1;
  }
  return $origin;
}


=head2 _import_assembly

  Title   : _import_assembly
  Usage   : $csp->_import_assembly($assemblyobj);
  Function: Update a contig spectrum object based on an assembly, contig or
            singlet object
  Returns : 1 for success
  Args    : Bio::Assembly::Scaffold, Contig or Singlet object

=cut

sub _import_assembly {
  my ($self, $assemblyobj) = @_;
  # Sanity check
  if ( ! ref $assemblyobj                                ||
       ( ! $assemblyobj->isa('Bio::Assembly::ScaffoldI') &&
         ! $assemblyobj->isa('Bio::Assembly::Contig')    )) {
    $self->throw("Unable to process non Bio::Assembly::ScaffoldI, Contig or ".
      "Singlet object [".ref($assemblyobj)."]");
  }
  # Create new object from assembly
  my $csp = $self->_new_from_assembly($assemblyobj);

  # Update current contig spectrum object with new one
  $self->add($csp);

  return 1;
}


=head2 _import_spectrum

  Title   : _import_spectrum
  Usage   : $csp->_import_spectrum({ 1 => 90 , 2 => 3 , 4 => 1 })
  Function: update a contig spectrum object based on a contig spectrum
            represented as a hash (key: contig size, value: number of contigs of
            this size)
  Returns : 1 for success
  Args    : contig spectrum as a hash reference

=cut

sub _import_spectrum {
  my ($self, $spectrum) = @_;
  # Sanity check
  if( ! ref $spectrum || ! ref $spectrum eq 'HASH') {
    $self->throw("Spectrum should be a hash reference, but it is [".
      ref($spectrum)."]");
  }

  # Update the spectrum (+ nof_rep, max_size and nof_seq)
  for my $size (keys %$spectrum) {
    # Get the number of contigs of different size
    if (defined $self->{'_spectrum'}{$size}) {
      $self->{'_spectrum'}{$size} += $$spectrum{$size};
    } else {
      $self->{'_spectrum'}{$size} = $$spectrum{$size};
    }
    # Update nof_seq
    $self->{'_nof_seq'} += $size * $$spectrum{$size};
    # Update max_size
    $self->{'_max_size'} = $size if $size > $self->{'_max_size'};
  }

  # If the contig spectrum has only zero 1-contigs, max_size is zero
  $self->{'_max_size'} = 0 if scalar keys %{$self->{'_spectrum'}} == 1 &&
    defined $self->{'_spectrum'}{'1'} && $self->{'_spectrum'}{'1'} == 0;

  # Update nof_rep
  $self->{'_nof_rep'}++;
  return 1;
}

=head2 _import_dissolved_csp

  Title   : _import_dissolved_csp
  Usage   : $csp->_import_dissolved_csp($mixed_csp, $seq_header);
  Function: Update a contig spectrum object by dissolving a mixed contig
            spectrum based on the header of the sequences
  Returns : 1 for success
  Args    : Bio::Assembly::Tools::ContigSpectrum
            sequence header string

=cut

sub _import_dissolved_csp {
  my ($self, $mixed_csp, $seq_header) = @_;
  # Sanity check
  if (not defined $mixed_csp || not defined $seq_header) {
    $self->throw("Expecting a contig spectrum reference and sequence header as".
    " arguments");
  }
  # Create new object from assembly
  my $dissolved_csp = $self->_new_dissolved_csp($mixed_csp, $seq_header);
  # Update current contig spectrum object with new one
  $self->add($dissolved_csp);
  return 1;
}


=head2 _import_cross_csp

  Title   : _import_cross_csp
  Usage   : $csp->_import_cross_csp($mixed_csp);
  Function: Update a contig spectrum object by calculating the cross contig
            spectrum based on a mixed contig spectrum
  Returns : 1 for success
  Args    : Bio::Assembly::Tools::ContigSpectrum

=cut

sub _import_cross_csp {
  my ($self, $mixed_csp) = @_;
  # Sanity check
  if (not defined $mixed_csp) {
    $self->throw("Expecting a contig spectrum reference as argument");
  }

  # Create new object from assembly
  my $cross_csp = $self->_new_cross_csp($mixed_csp);
  my $nof_1_cross_contigs = $cross_csp->spectrum->{1};

  # Update current contig spectrum object with new one
  $self->add($cross_csp);

  # Remove 1-contigs
  $self->{'_nof_seq'} -= $nof_1_cross_contigs;

  return 1;
}

=head2 _get_contig_like

  Title   : _get_contig_like
  Usage   : my @contig_like_objs = $csp->_get_contig_like($assembly_obj);
  Function: Get contigs and singlets from an assembly, contig or singlet
  Returns : array of Bio::Assembly::Contig and Singlet objects
  Args    : a Bio::Assembly::Scaffold, Contig or singlet object

=cut

sub _get_contig_like {
  my ($self, $assembly_obj) = @_;
  my @contig_objs;
  if ($assembly_obj->isa('Bio::Assembly::ScaffoldI')) {
    # all contigs and singlets in the scaffold
    push @contig_objs, ($assembly_obj->all_contigs, $assembly_obj->all_singlets);
  } else {
    # a contig or singlet
    @contig_objs = $assembly_obj;
  }
  return @contig_objs;
}


=head2 _get_assembly_seq_stats

  Title   : _get_assembly_seq_stats
  Usage   : my $seqlength = $csp->_get_assembly_seq_stats($assemblyobj);
  Function: Get sequence statistics from an assembly:
              average sequence length, number of sequences
  Returns : average sequence length (decimal)
            number of sequences (integer)
  Args    : Bio::Assembly::Scaffold, Contig or singlet object
            hash reference with the IDs of the sequences to consider [optional]

=cut

sub _get_assembly_seq_stats {
  my ($self, $assemblyobj, $seq_hash) = @_;

  # Sanity checks
  if ( !defined $assemblyobj ||
       ( !$assemblyobj->isa('Bio::Assembly::ScaffoldI') &&
         !$assemblyobj->isa('Bio::Assembly::Contig')       ) ) {
    $self->throw("Must provide a Bio::Assembly::Scaffold, Contig or Singlet object");
  }
  $self->throw("Expecting a hash reference. Got [".ref($seq_hash)."]")
    if (defined $seq_hash && ! ref($seq_hash) eq 'HASH');

  # Update sequence stats
  my @asm_stats = (0,0);
  # asm_stats = (avg_seq_len, nof_seq)
  for my $contigobj ( $self->_get_contig_like($assemblyobj) ) {
    @asm_stats = $self->_update_seq_stats( @asm_stats,
      $self->_get_contig_seq_stats($contigobj, $seq_hash) );
  }

  return @asm_stats;
}

=head2 _get_contig_seq_stats

  Title   : _get_contig_seq_stats
  Usage   : my $seqlength = $csp->_get_contig_seq_stats($contigobj);
  Function: Get sequence statistics from a contig:
              average sequence length, number of sequences
  Returns : average sequence length (decimal)
            number of sequences (integer)
  Args    : contig object reference
            hash reference with the IDs of the sequences to consider [optional]

=cut

sub _get_contig_seq_stats {
  my ($self, $contigobj, $seq_hash) = @_;
  my @contig_stats = (0, 0);
  # contig_stats = (avg_length, nof_seq)
  for my $seqobj ($contigobj->each_seq) {
    next if defined $seq_hash && !defined $$seq_hash{$seqobj->id};
    my $seq_string;
    if ($contigobj->isa('Bio::Assembly::Singlet')) { # a singlet
      $seq_string = $contigobj->seqref->seq;
    } else { # a contig
      $seq_string = $seqobj->seq;
    }
    # Number of non-gap characters in the sequence
    my $seq_len = $seqobj->_ungapped_len;
    my @seq_stats = ($seq_len);
    @contig_stats = $self->_update_seq_stats(@contig_stats, @seq_stats);
  }
  return @contig_stats;
}


=head2 _update_seq_stats

  Title   : _update_seq_stats
  Usage   : 
  Function: Update the number of sequences and their average length 1
            average identity 1
            minimum length 1
            minimum identity 1
            number of overlaps 1 average sequence length
  Returns : average sequence length
            number of sequences
  Args    : average sequence length 1
            number of sequences 1
            average sequence length 2
            number of sequences 2           

=cut

sub _update_seq_stats {
  my ($self, $p_avg_length, $p_nof_seq, $n_avg_length, $n_nof_seq) = @_;
  # Defaults
  if (not defined $n_nof_seq) {
    $n_nof_seq = 1;
  }
  # Update overlap statistics
  my $avg_length = 0;
  my $nof_seq = $p_nof_seq + $n_nof_seq;
  if ($nof_seq != 0) {
    $avg_length = ($p_avg_length * $p_nof_seq + $n_avg_length * $n_nof_seq) / $nof_seq;
  }
  return $avg_length, $nof_seq;
}


=head2 _get_assembly_overlap_stats

  Title   : _get_assembly_overlap_stats
  Usage   : my ($avglength, $avgidentity, $minlength, $min_identity, $nof_overlaps)
              = $csp->_get_assembly_overlap_stats($assemblyobj);
  Function: Get statistics about pairwise overlaps in contigs of an assembly
  Returns : average overlap length
            average identity percent
            minimum overlap length
            minimum identity percent
            number of overlaps
  Args    : Bio::Assembly::Scaffold, Contig or Singlet object
            hash reference with the IDs of the sequences to consider [optional]

=cut

sub _get_assembly_overlap_stats {
  my ($self, $assembly_obj, $seq_hash) = @_;

  # Sanity check
  if ( !defined $assembly_obj ||
       ( !$assembly_obj->isa('Bio::Assembly::ScaffoldI') &&
         !$assembly_obj->isa('Bio::Assembly::Contig')       ) ) {
    $self->throw("Must provide a Bio::Assembly::ScaffoldI, Contig or Singlet object");
  }
  $self->throw("Expecting a hash reference. Got [".ref($seq_hash)."]")
    if (defined $seq_hash && ! ref($seq_hash) eq 'HASH');

  # Look at all the contigs (no singlets!)
  my @asm_stats = (0, 0, undef, undef, 0);
  # asm_stats = (avg_length, avg_identity, min_length, min_identity, nof_overlaps)
  for my $contig_obj ( $self->_get_contig_like($assembly_obj) ) {
    @asm_stats = $self->_update_overlap_stats(  @asm_stats,
      $self->_get_contig_overlap_stats($contig_obj, $seq_hash) );
  }

  return @asm_stats;
}


=head2 _get_contig_overlap_stats

  Title   : _get_contig_overlap_stats
  Usage   : my ($avglength, $avgidentity, $minlength, $min_identity, $nof_overlaps)
              = $csp->_get_contig_overlap_stats($contigobj);
  Function: Get statistics about pairwise overlaps in a contig or singlet. The
              statistics are obtained using graph theory: each read is a node
              and the edges between 2 reads are weighted by minus the number of
              conserved residues in the alignment between the 2 reads. The
              minimum spanning tree of this graph represents the overlaps that
              form the contig. Overlaps that do not satisfy the minimum overlap
              length and similarity get a malus on their score.
              Note: This function requires the optional BioPerl dependency
              module called 'Graph'
  Returns : average overlap length
            average identity percent
            minimum overlap length
            minimum identity percent
            number of overlaps
  Args    : Bio::Assembly::Contig or Singlet object
            hash reference with the IDs of the sequences to consider [optional]

=cut

sub _get_contig_overlap_stats {
  my ($self, $contig_obj, $seq_hash) = @_;

  # Sanity check
  $self->throw("Must provide a Bio::Assembly::Contig object")
    if (!defined $contig_obj || !$contig_obj->isa("Bio::Assembly::Contig"));
  $self->throw("Expecting a hash reference. Got [".ref($seq_hash)."]")
    if (defined $seq_hash && ! ref($seq_hash) eq 'HASH');   

  my @contig_stats = (0, 0, undef, undef, 0);
  # contig_stats = (avg_length, avg_identity, min_length, min_identity, nof_overlaps)

  # Build contig graph
  ### consider providing the minima to _contig_graph here too?
  my ($g, $overlaps) = $self->_contig_graph($contig_obj, $seq_hash);

  if ( defined $g ) {
    # Graph minimum spanning tree (tree that goes through strongest overlaps)
    $g = $g->MST_Kruskal();

    # Calculate minimum overlap length and identity for this contig
    for my $edge ( $g->edges ) {
      # Retrieve overlap information
      my ($id1, $id2) = @$edge;
      if (not exists $$overlaps{$id1}{$id2}) {
        ($id2, $id1) = @$edge;
      }
      my ($score, $length, $identity) = @{$$overlaps{$id1}{$id2}};
      # Update contig stats
      my @overlap_stats = ($length, $identity);
      @contig_stats = $self->_update_overlap_stats(@contig_stats, @overlap_stats);
    }
  }

  return @contig_stats;
}


=head2 _update_overlap_stats

  Title   : _update_overlap_stats
  Usage   : 
  Function: update the number of overlaps and their minimum and average length
            and identity
  Returns : 
  Args    : average length 1
            average identity 1
            minimum length 1
            minimum identity 1
            number of overlaps 1
            average length 2
            average identity 2
            minimum length 2
            minimum identity 2
            number of overlaps 2

=cut

sub _update_overlap_stats {
  my ($self,
    $p_avg_length, $p_avg_identity, $p_min_length, $p_min_identity, $p_nof_overlaps,
    $n_avg_length, $n_avg_identity, $n_min_length, $n_min_identity, $n_nof_overlaps)
    = @_;

  # Defaults
  if (not defined $n_nof_overlaps) { $n_nof_overlaps = 1 };
  if ((not defined $n_min_length) && ($n_avg_length != 0))   { $n_min_length = $n_avg_length };
  if ((not defined $n_min_identity) && ($n_avg_identity != 0)) { $n_min_identity = $n_avg_identity };

  # Update overlap statistics
  my ($avg_length, $avg_identity, $min_length, $min_identity, $nof_overlaps)
    = (0, 0, undef, undef, 0);
  $nof_overlaps = $p_nof_overlaps + $n_nof_overlaps;
  if ($nof_overlaps > 0) {
    $avg_length = ($p_avg_length * $p_nof_overlaps + $n_avg_length * $n_nof_overlaps) / $nof_overlaps;
    $avg_identity = ($p_avg_identity * $p_nof_overlaps + $n_avg_identity * $n_nof_overlaps) / $nof_overlaps;
  }

  if ( not defined $p_min_length ) {
    $min_length = $n_min_length;
  } elsif ( not defined $n_min_length ) {
    $min_length = $p_min_length;
  } else { # both values are defined
    if ($n_min_length < $p_min_length) {
      $min_length = $n_min_length;
    } else {
      $min_length = $p_min_length;
    }
  }

  if ( not defined $p_min_identity ) {
    $min_identity = $n_min_identity;
  } elsif ( not defined $n_min_identity ) {
    $min_identity = $p_min_identity;
  } else { # both values are defined
    if ($n_min_identity < $p_min_identity) {
      $min_identity = $n_min_identity;
    } else {
      $min_identity = $p_min_identity;
    }
  }

  return $avg_length, $avg_identity, $min_length, $min_identity, $nof_overlaps;
}


=head2 _overlap_alignment

  Title   : _overlap_alignment
  Usage   : 
  Function: Produce an alignment of the overlapping section of two sequences of
            a contig. Minimum overlap length and percentage identity can be
            specified. Return undef if the sequences do not overlap or do not
            meet the minimum overlap criteria.
  Return  : Bio::SimpleAlign object reference
            alignment overlap length
            alignment overlap identity
  Args    : Bio::Assembly::Contig object reference
            Bio::LocatableSeq contig sequence 1
            Bio::LocatableSeq contig sequence 2
            minium overlap length [optional]
            minimum overlap identity percentage[optional]

=cut

sub _overlap_alignment {
  my ($self, $contig, $qseq, $tseq, $min_overlap, $min_identity) = @_;
  # get query and target sequence position
  my $qpos   = $contig->get_seq_coord($qseq);
  my $tpos   = $contig->get_seq_coord($tseq);
  # check that there is an overlap
  my $qend   = $qpos->end;
  my $tstart = $tpos->start;
  return if $qend < $tstart;
  my $qstart = $qpos->start;
  my $tend   = $tpos->end;
  return if $qstart > $tend;
  # get overlap boundaries and check overlap length
  my $left;
  if ($qstart >= $tstart) {
    $left = $qstart
  } else {
    $left = $tstart;
  }
  my $right;
  if ($qend > $tend) {
    $right = $tend;
  } else {
    $right = $qend;
  }
  my $overlap = $right - $left + 1;
  return if defined $min_overlap && $overlap < $min_overlap;
  # slice query and target sequence to overlap boundaries
  my $qleft =
    $contig->change_coord('gapped consensus', "aligned ".$qseq->id, $left);
  my $qstring = substr($qseq->seq, $qleft - 1, $overlap);
  my $tleft =
    $contig->change_coord('gapped consensus', "aligned ".$tseq->id, $left);
  my $tstring = substr($tseq->seq, $tleft - 1, $overlap);
  # remove gaps present in both sequences at the same position
  for (my $pos = 0 ; $pos < $overlap ; $pos++) {
    my $qnt = substr($qstring, $pos, 1);
    if ($qnt eq '-') {
      my $tnt = substr($tstring, $pos, 1);
      if ($tnt eq '-') {
        substr($qstring, $pos, 1, '');
        substr($tstring, $pos, 1, '');
        $pos--;
        $overlap--;
      }
    }
  }
  return if defined $min_overlap && $overlap < $min_overlap;
  # count the number of gaps remaining in each sequence
  my $qgaps = ($qstring =~ tr/-//);
  my $tgaps = ($tstring =~ tr/-//);
  # make an alignment object with the query and target sequences
  my $aln = Bio::SimpleAlign->new;
  my $alseq = Bio::LocatableSeq->new(
        -id       => 1,
        -seq      => $qstring,
        -start    => 1,
        -end      => $overlap - $qgaps,
        -alphabet => 'dna',
  );
  $aln->add_seq($alseq);
  $alseq = Bio::LocatableSeq->new(
        -id       => 2,
        -seq      => $tstring,
        -start    => 1,
        -end      => $overlap - $tgaps,
        -alphabet => 'dna',
  );
  $aln->add_seq($alseq);

  # check overlap percentage identity
  my $identity = $aln->overall_percentage_identity;
  return if defined $min_identity && $identity < $min_identity;

  # all checks passed, return alignment
  return $aln, $overlap, $identity;
}


=head2 _contig_graph

  Title   : _contig_graph
  Usage   : 
  Function: Creates a graph data structure of the contig.The graph is undirected.
            The vertices are the reads of the contig and edges are the overlap
            between the reads. The edges are weighted by the opposite of the
            overlap, so it is negative and the better the overlap, the lower the
            weight.
  Return  : Graph object or undef
            hashref of overlaps (score, length, identity) for each read pair
  Args    : Bio::Assembly::Contig object reference
            hash reference with the IDs of the sequences to consider [optional]
            minimum overlap length (integer)                         [optional]
            minimum percentage identity (integer)                    [optional]

=cut

sub _contig_graph {
  my ($self, $contig_obj, $seq_hash, $min_overlap, $min_identity) = @_;

  # Sanity checks
  if( !ref $contig_obj || ! $contig_obj->isa('Bio::Assembly::Contig') ) {
        $self->throw("Unable to process non Bio::Assembly::Contig ".
        "object [".ref($contig_obj)."]");
  }

  if (not eval { require Graph::Undirected }) {
    $self->throw("Error: the module 'Graph' is needed by the method ".
      "_contig_graph but could not be found\n$@");
  }

  # Skip contigs of 1 sequence (they have no overlap)
  my @seq_objs = $contig_obj->each_seq;
  my $nof_seqs = scalar @seq_objs;

  return if ($nof_seqs <= 1);

  # Calculate alignment between all pairs of reads
  my %overlaps;
  for my $i (0 .. $nof_seqs-1) {
    my $seq_obj = $seq_objs[$i];
    my $seq_id  = $seq_obj->id;

    # Skip this read if not in list of wanted sequences
    next if defined $seq_hash && !exists $$seq_hash{$seq_id};

    # What is the best sequence to align to?
    my ($best_score, $best_length, $best_identity);
    for my $j ($i+1 .. $nof_seqs-1) {

      # Skip this sequence if not in list of wanted sequences
      my $target_obj = $seq_objs[$j];
      my $target_id = $target_obj->id;
      next if defined $seq_hash && !exists $$seq_hash{$target_id};

      # How much overlap with this sequence?
      my ($aln_obj, $length, $identity)
        = $self->_overlap_alignment($contig_obj, $seq_obj, $target_obj, $min_overlap, $min_identity);
      next if ! defined $aln_obj; # there was no sequence overlap or overlap not good enough

      # Score the overlap as the number of conserved residues. In practice, it
      # seems to work better than giving +1 for match and -3 for errors
      # (mismatch or indels)
      my $score = $length * $identity / 100;

      # Apply a malus (square root) for scores that do not satisfy the minimum
      # overlap length similarity. It is necessary for overlaps that get a high
      # score without satisfying both the minimum values.
      if ( ( $min_overlap  && ($length   < $min_overlap ) ) ||
           ( $min_identity && ($identity < $min_identity) ) ) {
          $score = sqrt($score);
      }
      $overlaps{$seq_id}{$target_id} = [$score, $length, $identity];

    }

  }

  # Process overlaps
  my $g; # the Graph object
  if (scalar keys %overlaps >= 1) {
    # At least 1 overlap. Create a weighted undirected graph
    $g = Graph::Undirected->new();
    for my $seq_id (keys %overlaps) {
      for my $target_id (keys %{$overlaps{$seq_id}}) {
        my $score  = @{$overlaps{$seq_id}{$target_id}}[0];
        my $weight = -$score;
        $g->add_weighted_edge($seq_id, $target_id, $weight);
      }
    }

  }

  return $g, \%overlaps;

}


=head2 _draw_graph

  Title   : _draw_graph
  Usage   : 
  Function: Generates a PNG picture of the contig graph. It is mostly for
            debugging purposes.
  Return  : 1 for success
  Args    : a Graph object
            hashref of overlaps (score, length, identity) for each read pair
            name of output file
            overlap info to display: 'score' (default), 'length' or 'identity'

=cut

sub _draw_graph {
  my ($self, $g, $overlaps, $outfile, $edge_type) = @_;

  $self->throw("Error: need to provide a graph as input\n") if not defined $g;

  if (not eval { require GraphViz }) {
    $self->throw("Error: the module 'GraphViz' is needed by the method ".
      "_draw_graph but could not be found\n$@");
  }

  $edge_type ||= 'score';

  my $viz = GraphViz->new( directed => 0 );

  for my $edge ( $g->edges ) {
    # Retrieve overlap information
    my ($id1, $id2) = @$edge;
    if (not exists $$overlaps{$id1}{$id2}) {
      ($id2, $id1) = @$edge;
    }
    my ($score, $length, $identity) = @{$$overlaps{$id1}{$id2}};

    my $edge_val;
    if ($edge_type eq 'score') {
      $edge_val = $score;
    } elsif ($edge_type eq 'length') {
      $edge_val = $length;
    } elsif ($edge_type eq 'identity') {
      $edge_val = $identity;
    } else {
      $self->throw("Error: invalid edge type to display, '$edge_val'");
    }
    $viz->add_edge($id1 => $id2, label => $edge_val);
  }
  open my $fh, '>', $outfile or $self->throw("Error: Could not write file '$outfile': $!");
  print $fh $viz->as_png;
  close $fh;
  return 1;
}


1;

__END__
