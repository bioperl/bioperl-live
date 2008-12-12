#
# BioPerl module for Bio::Assembly::Tools::ContigSpectrum
#
# Copyright by Florent Angly
#
# You may distribute this module under the same terms as Perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Assembly::Tools::ContigSpectrum

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
  print "Average sequence length: ".$from_assembly->avg_seq_length." bp\n";
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

=head2 Implemention

The simplest representation of a contig spectrum is as a hash
representation where the key is the contig size (number of sequences
making up the contig) and the value the number of contigs of this
size.

In fact, it is useful to have more information associated with the
contig spectrum, hence the Bio::Assembly::Tools::ContigSpectrum module
implements an object containing a contig spectrum hash and additional
information. The get/set methods to access them are:

    id              contig spectrum ID
    nof_seq         number of sequences
    nof_rep         number of repetitions (assemblies) used
    max_size        size of (number of sequences in) the largest contig
    nof_overlaps    number of overlaps
    min_overlap     minimum overlap length for building a contig
    min_identity    minimum sequence identity over the overlap length
    avg_overlap     average overlap length
    avg_identity    average overlap identity
    avg_seq_length  average sequence length
    eff_asm_params  effective assembly parameters
    spectrum        hash representation of a contig spectrum

  Operations on the contig spectra:

    to_string       create a string representation of the spectrum
    spectrum        import a hash contig spectrum
    assembly        determine a contig spectrum from an assembly
    dissolve        calculate a dissolved contig spectrum (based on assembly)
    cross           produce a cross contig spectrum (based on assembly)
    add             add a contig spectrum to an existing one
    average         make an average of several contig spectra

When using operations that rely on knowing "where" (from what
metagenomes) a sequence came from (i.e. when creating a dissolved or
cross contig spectrum), make sure that the sequences used for the
assembly have a name header, e.g.  E<gt>metagenome1|seq1,
E<gt>metagenome2|seq1, ...

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
BioPerl modules. Send your comments and suggestions preferably to the
BioPerl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://bio.perl.org/MailList.html     - About the mailing lists

=head2 Reporting Bugs

Report bugs to the BioPerl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

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
use Bio::Align::PairwiseStatistics;

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
    $assembly, $eff_asm_params, $dissolve, $cross) = $self->_rearrange( [qw(ID
    NOF_SEQ NOF_REP MAX_SIZE NOF_OVERLAPS MIN_OVERLAP MIN_IDENTITY AVG_OVERLAP
    AVG_IDENTITY AVG_SEQ_LEN SPECTRUM ASSEMBLY EFF_ASM_PARAMS DISSOLVE CROSS)],
    @args );

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
  $self->{'_spectrum'}       = {1 => 0};  # contig spectrum hash representation
  $self->{'_assembly'}       = []; # list of assembly objects used
  
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
  
  # Finally get stuff that can be gotten in an automated way
  $self->_import_spectrum($spectrum) if defined($spectrum);
  $self->_import_assembly($assembly) if defined($assembly);
  if (defined($dissolve)) {
    my ($mixed_csp, $header) = (@$dissolve[0], @$dissolve[1]);
    $self->_import_dissolved_csp($mixed_csp, $header);
  }
  $self->_import_cross_csp($cross) if defined($cross);
  
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
            spectrum based or derived from an assembly is calulated. The
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
  Usage   : my @asm_list = $csp->assembly();
  Function: Get a reference to the list of assembly object reference used to
            make the contig spectrum object / Update the contig spectrum object
            based on an assembly object.
  Returns : array of Bio::Assembly::Scaffold
  Args    : Bio::Assembly::Scaffold

=cut

sub assembly {
  my ($self, $assembly) = @_;
  if (defined $assembly) {
    $self->_import_assembly($assembly);
  }
  my @asm_list = @{$self->{'_assembly'}} if defined $self->{'_assembly'};
  return \@asm_list;
}

=head2 drop_assembly

  Title   : drop_assembly
  Usage   : $csp->drop_assembly();
  Function: Remove all assembly objects associated with a contig spectrum.
            Assembly objects can be big. This method allows to free some memory
            when assembly information is not needed anymore.
  Returns : 1 for success, 0 for failure
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
            these sequences to the mixed contig spectrum based on the assembly.
            The mixed contig spectrum object must have been created based on one
            (or several) assembly object(s). Additionally, min_overlap and
            min_identity must have been set (either manually using min_overlap
            or automatically by switching on the eff_asm_params option).
  Returns : 1 for success, 0 for failure
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
  Returns : 1 for success, 0 for failure
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
  Returns : 1 for success, 0 for failure
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
    my $tot_num_overlaps = $csp->{'_nof_overlaps'} + $self->{'_nof_overlaps'};
    $self->{'_min_overlap'} = $csp->{'_min_overlap'} if
      defined $csp->{'_min_overlap'} && ( ! defined $self->{'_min_overlap'} ||
      $csp->{'_min_overlap'} < $self->{'_min_overlap'} );
    $self->{'_min_identity'} = $csp->{'_min_identity'} if
      defined $csp->{'_min_identity'} && ( ! defined $self->{'_min_identity'} ||
      $csp->{'_min_identity'} < $self->{'_min_identity'} );
    if ($tot_num_overlaps != 0) {
      $self->{'_avg_overlap'} =
        ($csp->{'_avg_overlap'} * $csp->{'_nof_overlaps'}
        + $self->{'_avg_overlap'} * $self->{'_nof_overlaps'})
        / $tot_num_overlaps;
      $self->{'_avg_identity'} =
        ($csp->{'_avg_identity'} * $csp->{'_nof_overlaps'}
        + $self->{'_avg_identity'} * $self->{'_nof_overlaps'})
        / $tot_num_overlaps;
    }
    $self->{'_nof_overlaps'} = $tot_num_overlaps;
  }
  # Update sequence statistics
  my $tot_nof_seq = $csp->{'_nof_seq'} + $self->{'_nof_seq'};
  if (not $tot_nof_seq == 0) {
    $self->{'_avg_seq_len'} = ($csp->{'_avg_seq_len'} * $csp->{'_nof_seq'} +
      $self->{'_avg_seq_len'} * $self->{'_nof_seq'}) / $tot_nof_seq;
  }
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
  $avg->{'_eff_asm_params'} = 1;
  
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


=head2 _naive_assembler

  Title   : _naive_assembler
  Usage   : 
  Function: Determines the contig spectrum (hash representation) of a subset of
            sequences from a mixed contig spectrum by "reassembling" the
            specified sequences only based on their position in the contig. This
            naive assembly only verifies that the minimum overlap length and
            percentage identity are respected. There is no actual alignment done
  Returns : contig spectrum hash reference
  Args    : Bio::Assembly::Contig
            sequence ID array reference
            minimum overlap length (integer) [optional]
            minimum percentage identity (integer) [optional]

=cut

sub _naive_assembler {
  my ($self, $contig, $seqlist, $min_overlap, $min_identity) = @_;
  # Sanity checks
  if ( ! ref $seqlist || ! ref($seqlist) eq 'ARRAY') {
    $self->throw('Expecting an array reference. Got ['.ref($seqlist)."] \n");
  }
  my $max = scalar @$seqlist;
  $self->throw("Expecting at least 2 sequences as input for _naive_assembler")
    if ($max < 2);
  # Assembly
  my %spectrum = (1 => 0);
  my %overlap_map;
  my %has_overlap;
  # Map what sequences overlap with what sequences
  for (my $i = 0 ; $i < $max-1 ; $i++) {
    # query sequence
    my $qseqid = $$seqlist[$i];
    my $qseq   = $contig->get_seq_by_name($qseqid);
    my $is_singlet = 1;
    for (my $j = $i+1 ; $j < $max ; $j++) {
      # target sequence
      my $tseqid = $$seqlist[$j];
      my $tseq = $contig->get_seq_by_name($tseqid);
      # try to align sequences
      my ($aln, $overlap, $identity)
        = $self->_overlap_alignment($contig, $qseq, $tseq, $min_overlap,
        $min_identity);
      # if there is no valid overlap, go to next sequence
      next if ! defined $aln;
      # the overlap is valid
      $is_singlet = 0;
      push @{$overlap_map{$qseqid}}, $tseqid;
      $has_overlap{$tseqid} = 1;
      $has_overlap{$qseqid} = 1;
    }
    # check if sequence is in previously seen overlap
    if (exists $has_overlap{$qseqid}) {
      $is_singlet = 0;
    }
    if ($is_singlet == 1) {
      $spectrum{1}++;
    }
  }
  # take care of last sequence
  my $last_is_singlet = 1;
  if (exists $has_overlap{$$seqlist[$max-1]}) {
    $last_is_singlet = 0;
  }
  if ($last_is_singlet == 1) {
    $spectrum{1}++;
  }
  # Parse overlap map
  for my $seqid (@$seqlist) {
    # list of sequences that should go in the contig
    next if not exists $overlap_map{$seqid};
    my @overlist = @{$overlap_map{$seqid}};
    for (my $j = 0 ; $j < scalar(@overlist) ; $j++) {
      my $otherseqid = $overlist[$j];
      if (exists $overlap_map{$otherseqid}) {
        push @overlist, @{$overlap_map{$otherseqid}};
        delete $overlap_map{$otherseqid};
      }
    }
    # remove duplicates from list
    @overlist = sort @overlist;
    for (my $j = 0 ; $j < scalar(@overlist)-1 ; $j++) {
      if ( $overlist[$j] eq $overlist[$j+1] ) {
        splice @overlist, $j, 1;
        $j--;
      }
    }
    # update spectrum with size of contig
    my $qsize = scalar(@overlist) + 1;
    if (defined $spectrum{$qsize}) {
      $spectrum{$qsize}++;
    } else {
      $spectrum{$qsize} = 1;
    }
  }
  return \%spectrum;
}


=head2 _new_from_assembly

  Title   : _new_from_assembly
  Usage   : 
  Function: Creates a new contig spectrum object based solely on the result of 
            an assembly
  Returns : Bio::Assembly::Tools::ContigSpectrum
  Args    : Bio::Assembly::Scaffold

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
  if ($csp->{'_eff_asm_params'} > 0) {
     my ($nover, $minl, $avgl, $minid, $avgid)
       = $csp->_get_overlap_stats($assemblyobj);
     $csp->{'_min_overlap'}  = $minl;
     $csp->{'_min_identity'} = $minid;
     $csp->{'_avg_overlap'}  = $avgl;
     $csp->{'_avg_identity'} = $avgid;
     $csp->{'_nof_overlaps'} = $nover;
  }
  # 3: Set sequence statistics: nof_seq and avg_seq_len
  my ($nseq, $avgseql) = $self->_get_seq_stats($assemblyobj);
  $csp->{'_avg_seq_len'} = $avgseql;
  $csp->{'_nof_seq'}     = $nseq;
  # 4: Set the spectrum: spectrum and max_size
  for my $contigobj ($assemblyobj->all_contigs) {
    my $size = $contigobj->no_sequences;
    if (defined $csp->{'_spectrum'}{$size}) {
      $csp->{'_spectrum'}{$size}++;
    } else {
      $csp->{'_spectrum'}{$size} = 1;
    }
    $csp->{'_max_size'} = $size if $size > $csp->{'_max_size'};
  }
  my $nof_singlets = $assemblyobj->get_nof_singlets();
  if (defined $nof_singlets) {
    $csp->{'_spectrum'}{1} += $nof_singlets;
    $csp->{'_max_size'} = 1 if $nof_singlets >= 1 && $csp->{'_max_size'} < 1;
  }
  # 5: Set list of assembly objects used
  push @{$csp->{'_assembly'}}, $assemblyobj;
  # 6: Set number of repetitions
  $csp->{'_nof_rep'} = 1;
  return $csp;
}



=head2 _new_dissolved_csp

  Title   : 
  Usage   : create a dissolved contig spectrum object
  Function: 
  Returns : 
  Args    : 


=cut

sub _new_dissolved_csp {
  my ($self, $mixed_csp, $seq_header) = @_;
  # Sanity checks on the mixed contig spectrum

  # min_overlap and min_identity must be specified if there are some overlaps
  # in the mixed contig
  unless ($mixed_csp->{'_nof_overlaps'} == 0) {
    unless ( defined $self->{'_min_overlap'} || 
      defined $mixed_csp->{'_min_overlap'} ) {
      $self->throw("min_overlap must be defined in the dissolved contig spectrum".
        " or mixed contig spectrum to dissolve a contig");
    }
    unless ( defined $self->{'_min_identity'} ||
      defined $mixed_csp->{'_min_identity'} ) {
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
  
  # take parent attributes if existent or mixed attributes otherwise
  if ($self->{'_eff_asm_params'}) {
    $dissolved->{'_eff_asm_params'} = $self->{'_eff_asm_params'};
  } else {
    $dissolved->{'_eff_asm_params'} = $mixed_csp->{'_eff_asm_params'};
  }
  if ($self->{'_min_overlap'} && $self->{'_min_identity'}) {
    ( $dissolved->{'_min_overlap'}, $dissolved->{'_min_identity'} ) = 
      ( $self->{'_min_overlap'}, $self->{'_min_identity'} );
  } else {
    ( $dissolved->{'_min_overlap'}, $dissolved->{'_min_identity'} ) = 
      ( $mixed_csp->{'_min_overlap'}, $mixed_csp->{'_min_identity'} );
  }
  
  # Dissolve each assembly
  for my $assembly (@{$mixed_csp->{'_assembly'}}) {
    # Dissolve this assembly for the given sequences
    my %asm_spectrum = (1 => 0);
    my %good_seqs;
    # For each contig
    for my $contig ($assembly->all_contigs) {
      # Get good sequences
      my @contig_seqs;
      for my $seq ($contig->each_seq) {
        my $seq_id = $seq->id;
        # get sequence origin
        next unless $seq_id =~ m/^$seq_header\|/;
        # add it to hash
        push @contig_seqs, $seq_id;
        $good_seqs{$seq_id} = 1;
      }
      # Update spectrum
      my $size = scalar @contig_seqs;
      if ($size == 0) {
        next;
      } elsif ($size == 1) {
        $asm_spectrum{1}++;
      } elsif ($size > 1) {
        # Reassemble good sequences
        my $contig_spectrum = $dissolved->_naive_assembler(
          $contig, \@contig_seqs, $dissolved->{'_min_overlap'},
          $dissolved->{'_min_identity'});
        # update spectrum
        for my $qsize (keys %$contig_spectrum) {
          $asm_spectrum{$qsize} += $$contig_spectrum{$qsize};
        }
      } else {
        $self->throw("The size is not valid... how could that happen?");
      }
    }
    # For each singlet
    for my $singlet ($assembly->all_singlets) {
      my $seq_id = $singlet->seqref->id;
      # get sequence origin
      next unless $seq_id =~ m/^$seq_header\|/;
      # add it to hash
      $good_seqs{$seq_id} = 1;
      # update spectrum
      $asm_spectrum{1}++;
    }
    # Update spectrum
    $dissolved->_import_spectrum(\%asm_spectrum);
    # Update nof_rep
    $dissolved->{'_nof_rep'}--;
    $dissolved->{'_nof_rep'} += $mixed_csp->{'_nof_rep'};

    # Get sequence stats
    my ($nseq, $avgseql) = $dissolved->_get_seq_stats($assembly, \%good_seqs);
    $dissolved->{'_avg_seq_len'} = $avgseql;
    $dissolved->{'_nof_seq'}     = $nseq;
  
    # Get eff_asm_param for these sequences
    if ($dissolved->{'_eff_asm_params'} > 0) {
      my ($nover, $minl, $avgl, $minid, $avgid)
        = $dissolved->_get_overlap_stats($assembly, \%good_seqs);
      $dissolved->{'_min_overlap'}  = $minl;
      $dissolved->{'_min_identity'} = $minid;
      $dissolved->{'_avg_overlap'}  = $avgl;
      $dissolved->{'_avg_identity'} = $avgid;
      $dissolved->{'_nof_overlaps'} = $nover;
    }

  }
  return $dissolved;
}


=head2 _new_cross_csp

  Title   : 
  Usage   : 
  Function: create a cross contig spectrum object
  Returns : 
  Args    : 


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
  my %spectrum = (1 => 0);
  
  # Take parent or mixed attributes
  if ($self->{'_eff_asm_params'}) {
    $cross->{'_eff_asm_params'} = $self->{'_eff_asm_params'};
  } else {
    $cross->{'_eff_asm_params'} = $mixed_csp->{'_eff_asm_params'};
  }
  if ($self->{'_min_overlap'} && $self->{'_min_identity'}) {
    ( $cross->{'_min_overlap'}, $cross->{'_min_identity'} ) = 
      ( $self->{'_min_overlap'}, $self->{'_min_identity'} );
  } else {
    ( $cross->{'_min_overlap'}, $cross->{'_min_identity'} ) = 
      ( $mixed_csp->{'_min_overlap'}, $mixed_csp->{'_min_identity'} );
  }
  
  # Get cross contig spectrum for each assembly
  for my $assembly (@{$mixed_csp->{'_assembly'}}) {
    # Go through contigs and skip the pure ones
    my %good_seqs;
    for my $contig ($assembly->all_contigs) {
      # Get origins
      my @seq_origins;
      my @seq_ids;
      for my $seq ($contig->each_seq) {
        # current sequence origin
        my $seq_id = $seq->id;
        $seq_id =~ m/^(.+)\|/;
        my $seq_header = $1;
        $self->warn("Sequence $seq_id does not seem to have a header. Skipping".
          " it...") if not defined $seq_header;
        $seq_header ||= '';
        push @seq_origins, $seq_header;
        push @seq_ids, $seq_id;
      }
      my $qsize = scalar(@seq_ids);
      my @origins = sort { $a cmp $b } @seq_origins;
      my $size = scalar(@origins);
      for (my $i = 1 ; $i < $size ; $i++) {
        if ($origins[$i] eq $origins[$i-1]) {
          splice @origins, $i, 1;
          $i--;
          $size--;
        }
      }
      # Update cross-contig number in spectrum
      if ($size > 1) { # cross-contig detected
        # update good sequences
        for my $seq_id (@seq_ids) {
          $good_seqs{$seq_id} = 1;
        }
        # update number of cross q-contigs in spectrum
        if (defined $spectrum{$qsize}) {
          $spectrum{$qsize} = 1;
        } else {
          $spectrum{$qsize}++;
        }
      }
      # Update number of cross 1-contigs
      if ($size > 1) { # cross-contig detected
        for my $origin (@origins) {
          # sequences to use
          my @ids;
          for (my $i = 0 ; $i < $qsize ; $i++) {
            my $seq_origin = $seq_origins[$i];
            my $seq_id = $seq_ids[$i];
            push @ids, $seq_id if $seq_origin eq $origin;
          }
          if (scalar @ids == 1) {
            $spectrum{1}++;
          } elsif (scalar @ids > 1) {
            my $contig_spectrum = $cross->_naive_assembler(
              $contig, \@ids, $cross->{'_min_overlap'},
              $cross->{'_min_identity'});
            $spectrum{1} += $$contig_spectrum{1};
          } else {
            $self->throw("The size is <= 0. How could such a thing happen?");
          }
        }
      }
    }
    # Get sequence stats
    my ($nseq, $avgseql) = $cross->_get_seq_stats($assembly, \%good_seqs);
    $cross->{'_avg_seq_len'} = $avgseql;
    $cross->{'_nof_seq'}     = $nseq;
    # Get eff_asm_param for these sequences
    if ($cross->{'_eff_asm_params'} > 0) {
      my ($nover, $minl, $avgl, $minid, $avgid)
        = $cross->_get_overlap_stats($assembly, \%good_seqs);
      $cross->{'_min_overlap'}  = $minl;
      $cross->{'_min_identity'} = $minid;
      $cross->{'_avg_overlap'}  = $avgl;
      $cross->{'_avg_identity'} = $avgid;
      $cross->{'_nof_overlaps'} = $nover;
    }
  }
  
  $cross->_import_spectrum(\%spectrum);
  # Update nof_rep
  $cross->{'_nof_rep'}--;
  $cross->{'_nof_rep'} += $mixed_csp->{'_nof_rep'};
  
  return $cross;
}

=head2 _import_assembly

  Title   : _import_assembly
  Usage   : $csp->_import_assembly($assemblyobj);
  Function: Update a contig spectrum object based on an assembly object
  Returns : 1 for success, 0 for error
  Args    : Bio::Assembly::Scaffold assembly object

=cut

sub _import_assembly {
  my ($self, $assemblyobj) = @_;
  # Sanity check
  if( !ref $assemblyobj || ! $assemblyobj->isa('Bio::Assembly::Scaffold') ) {
        $self->throw("Unable to process non Bio::Assembly::Scaffold assembly ".
        "object [".ref($assemblyobj)."]");
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
  Returns : 1 for success, 0 for error
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
  Returns : 1 for success, 0 for error
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
  Returns : 1 for success, 0 for error
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
  
  # Update current contig spectrum object with new one
  $self->add($cross_csp);
  
  return 1;
}


=head2 _get_seq_stats

  Title   : _get_seq_stats
  Usage   : my $seqlength = $csp->_get_seq_stats($assemblyobj);
  Function: Get sequence statistics from an assembly:
              number of sequences, average sequence length
  Returns : number of sequences (integer)
            average sequence length (decimal)
  Args    : assembly object reference
            hash reference with the IDs of the sequences to consider [optional]

=cut

sub _get_seq_stats {
  my ($self, $assemblyobj, $seq_hash) = @_;

  # sanity check
  $self->throw("Must provide a Bio::Assembly::Scaffold object")
    if (!defined $assemblyobj || !$assemblyobj->isa("Bio::Assembly::Scaffold"));
  $self->throw("Expecting a hash reference. Got [".ref($seq_hash)."]")
    if (defined $seq_hash && ! ref($seq_hash) eq 'HASH');

  my $avg_seq_len = 0;
  my $nof_seq = 0;
  for my $contigobj ($assemblyobj->all_contigs) {
    for my $seqobj ($contigobj->each_seq) {
      my $seq_id = $seqobj->id;
      next if defined $seq_hash && !defined $$seq_hash{$seq_id};
      $nof_seq++;
      my $seq_string = $seqobj->seq;
      $seq_string =~ s/-//g;
      $avg_seq_len += length($seq_string);
    }
  }
  for my $singletobj ($assemblyobj->all_singlets) {
    my $seq_id = $singletobj->seqref->id;
    next if defined $seq_hash && !defined $$seq_hash{$seq_id};
    $nof_seq++;
    my $seq_string = $singletobj->seqref->seq;
    $seq_string =~ s/-//g;
    $avg_seq_len += length($seq_string);
  }
  $avg_seq_len /= $nof_seq unless $nof_seq == 0;
  return $nof_seq, $avg_seq_len;
}


=head2 _get_overlap_stats

  Title   : _get_overlap_stats
  Usage   : my ($minlength, $min_identity, $avglength, $avgidentity)
              = $csp->_get_overlap_stats($assemblyobj);
  Function: Get statistics about pairwise overlaps in contigs of an assembly
  Returns : number of overlaps
            minimum overlap length
            average overlap length
            minimum identity percent
            average identity percent
  Args    : assembly object reference
            hash reference with the IDs of the sequences to consider [optional]

=cut

sub _get_overlap_stats {
  my ($self, $assembly_obj, $seq_hash) = @_;

  # sanity check
  $self->throw("Must provide a Bio::Assembly::Scaffold object")
    if (!defined $assembly_obj || !$assembly_obj->isa("Bio::Assembly::Scaffold"));
  $self->throw("Expecting a hash reference. Got [".ref($seq_hash)."]")
    if (defined $seq_hash && ! ref($seq_hash) eq 'HASH');
  
  my $matchdef = $self->{'_eff_asm_params'};
  my ($min_length, $avg_length, $min_identity, $avg_identity, $nof_overlaps)
    = (undef, 0, undef, 0, 0);
  
  # Look at all the contigs (and I really mean no singlets!)
  for my $contig_obj ($assembly_obj->all_contigs) {
    my $nof_seq = 0;

    # Look at best overlap possible with previous sequences in contig
    my @all_seq_objs = $contig_obj->each_seq;
    # sequences should be ordered by starting position
    for (my $i = 0 ; $i < scalar(@all_seq_objs) ; $i++) {
      my $seq_obj    = $all_seq_objs[$i];
      my $seq_id    = $seq_obj->id;
      
      # skip this sequence if not in list of wanted sequences
      next if defined $seq_hash && !defined $$seq_hash{$seq_id};
      $nof_seq++;
      
      # skip the first sequence (no other sequence to compare against)
      next if $nof_seq <= 1;
      
      # what is the best previous sequence to align to?
      my $stats = Bio::Align::PairwiseStatistics->new;
      my $target_obj;
      my $target_id;
      my $best_score;
      my $best_length;
      my $best_identity;
      
      for (my $j = $i-1 ; $j >= 0 ; $j--) {
        my $tmp_target_obj = $all_seq_objs[$j];
        my $tmp_target_id = $tmp_target_obj->id;
        
        # skip this sequence if not in list of wanted sequences
        next if defined $seq_hash && !defined $$seq_hash{$tmp_target_id};
        
        # find overlap with that sequence
        my ($aln_obj, $tmp_length, $tmp_identity)
          = $self->_overlap_alignment($contig_obj, $seq_obj, $tmp_target_obj);
        next if ! defined $aln_obj; # there was no sequence overlap
        my $tmp_score = $stats->score_nuc($aln_obj);
        
        # update score and best sequence for overlap
        if (!defined $best_score || $best_score < $tmp_score) {
          $best_score    = $tmp_score;
          $best_length   = $tmp_length;
          $best_identity = $tmp_identity;
          $target_obj    = $tmp_target_obj;
          $target_id     = $tmp_target_id;
        }
      }
      
      # Update our overlap statistics
      if (defined $best_score) {
        $avg_length += $best_length;
        $avg_identity += $best_identity;
        $min_length = $best_length if ! defined $min_length ||
          $best_length < $min_length;
        $min_identity = $best_identity if ! defined $min_identity ||
          $best_identity < $min_identity;
        $nof_overlaps++;
      }
    }
  }
  
  # averaging
  unless ($nof_overlaps == 0) {
    $avg_length /= $nof_overlaps;
    $avg_identity /= $nof_overlaps;
  }
  
  return $nof_overlaps, $min_length, $avg_length, $min_identity, $avg_identity;
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
            minimum overlap percentage identity [optional]

=cut

sub _overlap_alignment {
  my ($self, $contig, $qseq, $tseq, $min_overlap, $min_identity) = @_;
  # get query sequence position  
  my $qpos   = $contig->get_seq_coord($qseq);
  my $qstart = $qpos->start;
  my $qend   = $qpos->end;
  # get target sequence position
  my $tpos   = $contig->get_seq_coord($tseq);
  my $tstart = $tpos->start;
  my $tend   = $tpos->end;
  # check that there is an overlap
  return if $qstart > $tend || $qend < $tstart;
  # get overlap boundaries and check overlap length
  my $left = $qstart;
  $left = $tstart if $qstart < $tstart;
  my $right = $qend;
  $right = $tend if $qend > $tend;
  my $overlap = $right - $left + 1;
  return if defined $min_overlap && $overlap < $min_overlap;
  # slice query and target sequence to overlap boundaries
  my $qleft = $contig->change_coord('gapped consensus', "aligned ".$qseq->id,
    $left);
  my $qright = $qleft + $overlap - 1;
  my $qstring = $qseq->seq;
  $qstring = substr($qstring, $qleft - 1, $overlap);
  my $tleft = $contig->change_coord('gapped consensus', "aligned ".$tseq->id, 
    $left);
  my $tright = $tleft + $overlap - 1;
  my $tstring = $tseq->seq;
  $tstring = substr($tstring, $tleft - 1, $overlap);
  # remove gaps present in both sequences at the same position
  for (my $pos = 0 ; $pos < $overlap ; $pos++) {
    my $qnt = substr($qstring, $pos, 1);
    my $tnt = substr($tstring, $pos, 1);
    if ($qnt eq '-' && $tnt eq '-') {
      substr($qstring, $pos, 1, '');
      substr($tstring, $pos, 1, '');
      $pos--;
      $overlap--;
    }
  }
  return if defined $min_overlap && $overlap < $min_overlap;
  # make an aligned object
  my $aln = Bio::SimpleAlign->new;
  my $qalseq = Bio::LocatableSeq->new(
        -id       => 1,
        -seq      => $qstring,
        -start    => 1,
        -alphabet => 'dna'
  );
  $aln->add_seq($qalseq);
  my $talseq = Bio::LocatableSeq->new(
        -id       => 2,
        -seq      => $tstring,
        -start    => 1,
        -alphabet => 'dna'
  );
  $aln->add_seq($talseq);
  # check overlap percentage identity
  my $identity = $aln->overall_percentage_identity;
  return if defined $min_identity && $identity < $min_identity;
  # all checks passed, return alignment
  return $aln, $overlap, $identity;
}

1;

__END__
