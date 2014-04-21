# BioPerl module for Bio::Tools::AmpliconSearch
#
# Copyright Florent Angly
#
# You may distribute this module under the same terms as perl itself


package Bio::Tools::AmpliconSearch;

use strict;
use warnings;
use Bio::Tools::IUPAC;
use Bio::SeqFeature::Amplicon;
use Bio::Tools::SeqPattern;
# we require Bio::SeqIO
# and Bio::SeqFeature::Primer

use base qw(Bio::Root::Root);

my $template_str;


=head1 NAME

Bio::Tools::AmpliconSearch - Find amplicons in a template using degenerate PCR primers

=head1 SYNOPSIS

   use Bio::PrimarySeq;
   use Bio::Tools::AmpliconSearch;

   my $template = Bio::PrimarySeq->new(
      -seq => 'aaaaaCCCCaaaaaaaaaaTTTTTTaaaaaCCACaaaaaTTTTTTaaaaaaaaaa',
   );
   my $fwd_primer = Bio::PrimarySeq->new(
      -seq => 'CCNC',
   );
   my $rev_primer = Bio::PrimarySeq->new(
      -seq => 'AAAAA',
   );

   my $search = Bio::Tools::AmpliconSearch->new(
      -template   => $template,
      -fwd_primer => $fwd_primer,
      -rev_primer => $rev_primer,
   );
   
   while (my $amplicon = $search->next_amplicon) {
      print "Found amplicon at position ".$amplicon->start.'..'.$amplicon->end.":\n";
      print $amplicon->seq->seq."\n\n";
   }

   # Now change the template (but you could change the primers instead) and look
   # for amplicons again

   $template = Bio::PrimarySeq->new(
      -seq => 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa',
   );
   $search->template($template);

   while (my $amplicon = $search->next_amplicon) {
      print "Found amplicon at position ".$amplicon->start.'..'.$amplicon->end.":\n";
      print $amplicon->seq->seq."\n\n";
   }

=head1 DESCRIPTION

Perform an in silico PCR reaction, i.e. search for amplicons in a given template
sequence using the specified degenerate primer.

The template sequence is a sequence object, e.g. L<Bio::Seq>, and the primers
can be a sequence or a L<Bio::SeqFeature::Primer> object and contain ambiguous
residues as defined in the IUPAC conventions. The primer sequences are converted
into regular expressions using L<Bio::Tools::IUPAC> and the matching regions of
the template sequence, i.e. the amplicons, are returned as L<Bio::Seq::PrimedSeq>
objects.

AmpliconSearch will look for amplicons on both strands (forward and reverse-
complement) of the specified template sequence. If the reverse primer is not
provided, an amplicon will be returned and span a match of the forward primer to
the end of the template. Similarly, when no forward primer is given, match from
the beginning of the template sequence. When several amplicons overlap, only the
shortest one to more accurately represent the biases of PCR. Future improvements
may include modelling the effects of the number of PCR cycles or temperature on
the PCR products.

=head1 TODO

Future improvements may include:

=over

=item *

Allowing a small number of primer mismatches

=item *

Reporting all amplicons, including overlapping ones

=item *

Putting a limit on the length of amplicons, in accordance with the processivity
of the polymerase used

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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

=head1 AUTHOR

Florent Angly <florent.angly@gmail.com>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=head2 new

 Title    : new
 Usage    : my $search = Bio::Tools::AmpliconSearch->new( );
 Function : Initialize an amplicon search
 Args     : -template       Sequence object for the template sequence. This object
                            will be converted to Bio::Seq if needed in since features
                            (amplicons and primers) will be added to this object.
            -fwd_primer     A sequence object representing the forward primer
            -rev_primer     A sequence object representing the reverse primer
            -primer_file    Read primers from a sequence file. It replaces
                            -fwd_primer and -rev_primer (optional)
            -attach_primers Whether or not to attach primers to Amplicon objects. Default: 0 (off)
 Returns  : A Bio::Tools::AmpliconSearch object

=cut

sub new {
   my ($class, @args) = @_;
   my $self = $class->SUPER::new(@args);
   my ($template, $primer_file, $fwd_primer, $rev_primer, $attach_primers) =
      $self->_rearrange([qw(TEMPLATE PRIMER_FILE FWD_PRIMER REV_PRIMER ATTACH_PRIMERS)],
      @args);

   # Get primers
   if (defined $primer_file) {
      $self->primer_file($primer_file);
   } else {
      $self->fwd_primer($fwd_primer || '');
      $self->rev_primer($rev_primer || '');
   }

   # Get template sequence
   $self->template($template) if defined $template;

   $self->attach_primers($attach_primers) if defined $attach_primers;

   return $self;
}


=head2 template

 Title    : template
 Usage    : my $template = $search->template;
 Function : Get/set the template sequence. Setting a new template resets any
            search in progress.
 Args     : Optional Bio::Seq object
 Returns  : A Bio::Seq object

=cut

sub template {
   my ($self, $template) = @_;
   if (defined $template) {
      if ( not(ref $template) || not $template->isa('Bio::PrimarySeqI') ) {
         # Not a Bio::Seq or Bio::PrimarySeq
         $self->throw("Expected a sequence object as input but got a '".ref($template)."'\n");
      }
      if (not $template->isa('Bio::SeqI')) {
         # Convert sequence object to Bio::Seq Seq so that features can be added
         my $primary_seq = $template;
         $template = Bio::Seq->new();
         $template->primary_seq($primary_seq);
      }
      $self->{template} = $template;
      # Reset search in progress
      $template_str = undef;
   }
   return $self->{template};
}


=head2 fwd_primer

 Title    : fwd_primer
 Usage    : my $primer = $search->fwd_primer;
 Function : Get/set the forward primer. Setting a new forward primer resets any
            search in progress.
 Args     : Optional sequence object or primer object or '' to match beginning
            of sequence.            
 Returns  : A sequence object or primer object or undef

=cut

sub fwd_primer {
   my ($self, $primer) = @_;
   if (defined $primer) {
      $self->_set_primer('fwd', $primer);
   }
   return $self->{fwd_primer};
}


=head2 rev_primer

 Title    : rev_primer
 Usage    : my $primer = $search->rev_primer;
 Function : Get/set the reverse primer. Setting a new reverse primer resets any
            search in progress.
 Args     : Optional sequence object or primer object or '' to match end of
            sequence. 
 Returns  : A sequence object or primer object or undef

=cut

sub rev_primer {
   my ($self, $primer) = @_;
   if (defined $primer) {
      $self->_set_primer('rev', $primer);
   }
   return $self->{rev_primer};
}


sub _set_primer {
   # Save a primer (sequence object) and convert it to regexp. Type is 'fwd' for
   # the forward primer or 'rev' for the reverse primer.
   my ($self, $type, $primer) = @_;
   my $re;
   my $match_rna = 1;
   if ($primer eq '') {
      $re = $type eq 'fwd' ? '^' : '$';
   } else {
      if ( not(ref $primer) || (
           not($primer->isa('Bio::PrimarySeqI')) &&
           not($primer->isa('Bio::SeqFeature::Primer')) ) ) {
         $self->throw('Expected a sequence or primer object as input but got a '.ref($primer)."\n");
      }
      $self->{$type.'_primer'} = $primer;
      my $seq = $primer->isa('Bio::SeqFeature::Primer') ? $primer->seq : $primer;
      $re = Bio::Tools::IUPAC->new(
         -seq => $type eq 'fwd' ? $seq : $seq->revcom,
      )->regexp($match_rna);
   }
   $self->{$type.'_regexp'} = $re;
   # Reset search in progress
   $template_str = undef;
   $self->{regexp} = undef;
   return $self->{$type.'_primer'};
}


=head2 primer_file

 Title    : primer_file
 Usage    : my ($fwd, $rev) = $search->primer_file;
 Function : Get/set a sequence file to read the primer from. The first sequence
            must be the forward primer, and the second is the optional reverse
            primer. After reading the file, the primers are set using fwd_primer()
            and rev_primer() and returned.
 Args     : Sequence file
 Returns  : Array containing forward and reverse primers as sequence objects.

=cut

sub primer_file {
   my ($self, $primer_file) = @_;
   # Read primer file and convert primers into regular expressions to catch
   # amplicons present in the database

   if (not defined $primer_file) {
      $self->throw("Need to provide an input file\n");
   }

   # Mandatory first primer
   require Bio::SeqIO;
   my $in = Bio::SeqIO->new( -file => $primer_file );
   my $fwd_primer = $in->next_seq;
   if (not defined $fwd_primer) {
      $self->throw("The file '$primer_file' contains no primers\n");
   }
   $fwd_primer->alphabet('dna'); # Force the alphabet since degenerate primers can look like protein sequences

   # Optional reverse primers
   my $rev_primer = $in->next_seq;
   if (defined $rev_primer) {
      $rev_primer->alphabet('dna');
   } else {
      $rev_primer = '';
   }
   
   $in->close;

   $self->fwd_primer($fwd_primer);
   $self->rev_primer($rev_primer);

   return ($fwd_primer, $rev_primer);
}


=head2 attach_primers

 Title    : attach_primers
 Usage    : my $attached = $search->attach_primers;
 Function : Get/set whether or not to attach primer objects to the amplicon
            objects.
 Args     : Optional integer (1 for yes, 0 for no)
 Returns  : Integer (1 for yes, 0 for no)

=cut

sub attach_primers {
   my ($self, $attach) = @_;
   if (defined $attach) {
      $self->{attach_primers} = $attach;
      require Bio::SeqFeature::Primer;
   }
   return $self->{attach_primers} || 0;
}


=head2 next_amplicon

 Title    : next_amplicon
 Usage    : my $amplicon = $search->next_amplicon;
 Function : Get the next amplicon
 Args     : None
 Returns  : A Bio::SeqFeature::Amplicon object

=cut

sub next_amplicon {
   my ($self) = @_;

   # Initialize search
   if (not defined $template_str) {
      $self->_init;
   }

   my $re = $self->_regexp;

   my $amplicon;
   if ($template_str  =~ m/$re/g) {
      my ($match, $rev_match) = ($1, $2);
      my $strand = $rev_match ? -1 : 1;
      $match = $match || $rev_match;
      my $end   = pos($template_str);
      my $start = $end - length($match) + 1;
      $amplicon = $self->_attach_amplicon($start, $end, $strand);
   }

   # If no more matches. Make sure calls to next_amplicon() will return undef.
   if (not $amplicon) {
      $template_str = '';
   }

   return $amplicon;
}


sub _init {
   my ($self) = @_;
   # Sanity checks
   if ( not $self->template ) {
      $self->throw('Need to provide a template sequence');
   }
   if ( not($self->fwd_primer) && not($self->rev_primer) ) {
      $self->throw('Need to provide at least a primer');
   }
   # Set the template sequence string
   $template_str = $self->template->seq;
   # Set the regular expression to match amplicons
   $self->_regexp;

   return 1;
}


sub _regexp {
   # Get the regexp to match amplicon. If the regexp is not set, initialize it.
   my ($self, $regexp) = @_;

   if ( not defined $self->{regexp} ) {
      # Build regexp that matches amplicons on both strands and reports shortest
      # amplicon when there are several overlapping amplicons

      my $fwd_regexp = $self->_fwd_regexp;
      my $rev_regexp = $self->_rev_regexp;

      my ($fwd_regexp_rc, $basic_fwd_match, $rev_regexp_rc, $basic_rev_match);
      if ($fwd_regexp eq '^') {
         $fwd_regexp_rc = '';
         $basic_fwd_match = "(?:.*?$rev_regexp)";
      } else {
         $fwd_regexp_rc = Bio::Tools::SeqPattern->new(
            -seq  => $fwd_regexp,
            -type => 'dna',
         )->revcom->str;
         $basic_fwd_match = "(?:$fwd_regexp.*?$rev_regexp)";
      }

      if ($rev_regexp eq '$') {
         $rev_regexp_rc = '';
         $basic_rev_match = "(?:.*?$fwd_regexp_rc)";
      } else {
         $rev_regexp_rc = Bio::Tools::SeqPattern->new(
            -seq  => $rev_regexp,
            -type => 'dna',
         )->revcom->str;
         $basic_rev_match = "(?:$rev_regexp_rc.*?$fwd_regexp_rc)";
      }

      my $fwd_exclude     = "(?!$basic_rev_match".
                            ($fwd_regexp eq '^' ? '' : "|$fwd_regexp").
                            ")";

      my $rev_exclude     = "(?!$basic_fwd_match".
                            ($rev_regexp eq '$' ? '' : "|$rev_regexp_rc").
                            ')';

      $self->{regexp} = qr/
                      ( $fwd_regexp    (?:$fwd_exclude.)*? $rev_regexp    ) |
                      ( $rev_regexp_rc (?:$rev_exclude.)*? $fwd_regexp_rc )
      /xi;
   }

   return $self->{regexp};
}


=head2 annotate_template

 Title    : annotate_template
 Usage    : my $template = $search->annotate_template;
 Function : Search for all amplicons and attach them to the template.
            This is equivalent to running:
               while (my $amplicon = $self->next_amplicon) {
                  # do something
               }
               my $annotated = $self->template;
 Args     : None
 Returns  : A Bio::Seq object with attached Bio::SeqFeature::Amplicons (and
            Bio::SeqFeature::Primers if you set -attach_primers to 1).

=cut

sub annotate_template {
   my ($self) = @_;
   # Search all amplicons and attach them to template
   1 while $self->next_amplicon;
   # Return annotated template
   return $self->template;
}


sub _fwd_regexp {
   my ($self) = @_;
   return $self->{fwd_regexp};
}


sub _rev_regexp {
   my ($self) = @_;
   return $self->{rev_regexp};
}


sub _attach_amplicon {
   # Create an amplicon object and attach it to template
   my ($self, $start, $end, $strand) = @_;

   # Create Bio::SeqFeature::Amplicon feature and attach it to the template
   my $amplicon = Bio::SeqFeature::Amplicon->new(
      -start    => $start,
      -end      => $end,
      -strand   => $strand,
      -template => $self->template,
   );

   # Create Bio::SeqFeature::Primer feature and attach them to the amplicon
   if ($self->attach_primers) {
      for my $type ('fwd', 'rev') {
         my ($pstart, $pend, $pstrand, $primer_seq);

         # Coordinates relative to amplicon
         if ($type eq 'fwd') {
            # Forward primer
            $primer_seq = $self->fwd_primer;
            next if not defined $primer_seq;
            $pstart  = 1;
            $pend    = $primer_seq->length;
            $pstrand = $amplicon->strand;
         } else {
            # Optional reverse primer
            $primer_seq = $self->rev_primer;
            next if not defined $primer_seq;
            $pstart  = $end - $primer_seq->length + 1;
            $pend    = $end;
            $pstrand = -1 * $amplicon->strand;
         }

         # Absolute coordinates needed
         $pstart += $start - 1;
         $pend   += $start - 1;

         my $primer = Bio::SeqFeature::Primer->new(
            -start    => $pstart,
            -end      => $pend,
            -strand   => $pstrand,
            -template => $amplicon,
         );

         # Attach primer to amplicon
         if ($type eq 'fwd') {
            $amplicon->fwd_primer($primer);
         } else {
            $amplicon->rev_primer($primer);
         }

      }
   }

   return $amplicon;
}


1;
