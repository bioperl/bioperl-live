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
# we require Bio::SeqIO
# and Bio::SeqFeature::Primer

use base qw(Bio::Root::Root);

my $template_str;


=head1 NAME

Bio::Tools::AmpliconSearch - Use PCR primers to identify amplicons in a template sequence

=head1 SYNOPSIS

XXX

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

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR

Florent Angly <florent.angly@gmail.com>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=head2 new

 Title    : new
 Usage    : my $search = Bio::Tools::AmpliconSearch->new( );
 Function : 
 Args     : -template       Sequence object for the template sequence. This object
                            will be converted to Bio::Seq if needed since features
                            (amplicons and primers) will be added to this object.
            -fwd_primer     A sequence object representing the forward primer
            -rev_primer     A sequence object representing the reverse primer
            -primer_file    Read primers from a sequence file. It replaces -fwd_primer and -rev_primer (optional)
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
      ($fwd_primer, $rev_primer) = $self->_get_primers_from_file($primer_file);
   }
   $self->_set_fwd_primer($fwd_primer);
   $self->_set_rev_primer($rev_primer);

   # Get template sequence
   $self->_set_template($template) if defined $template;
   if ( $template && not($fwd_primer) && not($rev_primer) ) {
      $self->throw('Need to provide at least a primer');
   }

   $self->_set_attach_primers($attach_primers) if defined $attach_primers;

   return $self;
}


=head2 template

 Title    : template
 Usage    : my $template = $search->template;
 Function : Get the template sequence
 Args     : None
 Returns  : A Bio::Seq object

=cut

sub template {
   my ($self) = @_;
   return $self->{template};
}

sub _set_template {
   my ($self, $template) = @_;
   if ( not(ref $template) || not $template->isa('Bio::PrimarySeqI')) {
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
   $template_str = $self->template->seq;
   $self->_set_strand(1);
   return $self->template;
}


=head2 fwd_primer

 Title    : fwd_primer
 Usage    : my $primer = $search->fwd_primer;
 Function : Get the forward primer.
 Args     : None
 Returns  : A sequence object or primer object or undef

=cut

sub fwd_primer {
   my ($self) = @_;
   return $self->{fwd_primer};
}

sub _set_fwd_primer {
   my ($self, $primer) = @_;
   return $self->_set_primer('fwd', $primer);
}


=head2 rev_primer

 Title    : rev_primer
 Usage    : my $primer = $search->rev_primer;
 Function : Get the reverse primer.
 Args     : None
 Returns  : A sequence object or primer object or undef

=cut

sub rev_primer {
   my ($self) = @_;
   return $self->{rev_primer};
}

sub _set_rev_primer {
   my ($self, $primer) = @_;
   return $self->_set_primer('rev', $primer);
}


sub _set_primer {
   # Save a primer (sequence object) and convert it to regexp. Type is 'fwd' or 'rev'.
   my ($self, $type, $primer) = @_;
   my $re;
   if (defined $primer) {
      if ( not(ref $primer) || (
           not($primer->isa('Bio::PrimarySeqI')) &&
           not($primer->isa('Bio::SeqFeature::Primer')) ) ) {
         $self->throw('Expected a sequence or primer object as input but got a '.ref($primer)."\n");
      }
      $self->{$type.'_primer'} = $primer;
      my $seq = $primer->isa('Bio::SeqFeature::Primer') ? $primer->seq : $primer;
      $re = Bio::Tools::IUPAC->new(
         -seq => $type eq 'fwd' ? $seq : $seq->revcom,
      )->regexp;
   } else {
      $re = $type eq 'fwd' ? qr/^/ : qr/$/;
   }
   $self->{$type.'_regexp'} = $re;
   return $self->{$type.'_primer'};
}


sub _get_primers_from_file {
   my ($self, $primer_file) = @_;
   # Read primer file and convert primers into regular expressions to catch
   # amplicons present in the database

   if (not defined $primer_file) {
      $self->throw("Need to provide an input file\n");
   }

   # Mandatory first primer
   require Bio::SeqIO;
   my $in = Bio::SeqIO->newFh( -file => $primer_file );
   my $fwd_primer = <$in>;
   if (not defined $fwd_primer) {
      $self->throw("The file '$primer_file' contains no primers\n");
   }
   $fwd_primer->alphabet('dna'); # Force the alphabet since degenerate primers can look like protein sequences

   # Optional reverse primers
   my $rev_primer = <$in>;
   if (defined $rev_primer) {
      $rev_primer->alphabet('dna');
   }
   
   #### $in->close;
   #### close $in;
   undef $in;

   return ($fwd_primer, $rev_primer);
}


sub fwd_regexp {
   my ($self) = @_;
   return $self->{fwd_regexp};
}


sub rev_regexp {
   my ($self) = @_;
   return $self->{rev_regexp};
}


=head2 attach_primers

 Title    : attach_primers
 Usage    : my $primers_attached = $search->attach_primers;
 Function : Get whether or not primer objects will be attached to the amplicon
            objects.
 Args     : None
 Returns  : Integer (1 for yes, 0 for no)

=cut

sub attach_primers {
   my ($self) = @_;
   return $self->{attach_primers} || 0;
}

sub _set_attach_primers {
   my ($self, $val) = @_;
   $self->{attach_primers} = $val;
   require Bio::SeqFeature::Primer;
   return $self->attach_primers;
}


sub _cur_strand {
   my ($self) = @_;
   return $self->{cur_strand};
}

sub _set_strand {
   my ($self, $strand) = @_;
   $self->{cur_strand} = $strand;
   return $self->_cur_strand;
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
   my $amplicon;

   my $strand = $self->_cur_strand;
   my $fwd_regexp = $self->fwd_regexp;
   my $rev_regexp = $self->rev_regexp;

   if ($template_str  =~ m/($fwd_regexp.*?$rev_regexp)/g) {
      my $end   = pos($template_str);
      my $start = $end - length($1) + 1;
      # Now trim the left end to obtain the shortest amplicon
      my $ampliconstr = substr $template_str, $start - 1, $end - $start + 1;
      if ($ampliconstr =~ m/$fwd_regexp.*($fwd_regexp)/g) {
         $start += pos($ampliconstr) - length($1);
      }
      $amplicon = $self->_attach_amplicon($start, $end, $strand);
   }

   if ( not $amplicon ) {
      if ( $strand == 1 ) {
         # Exhausted all matches in forward strand. Search in the reverse strand.
         $template_str = $self->template->revcom->seq;
         $self->_set_strand(-1);
         $amplicon = $self->next_amplicon;
      } else {
         # No more matches. Make sure calls to next_amplicon() will return undef.
         $template_str = '';
      }
   }

#  # May want to use Bio::Range intersection() to check for overlap

#  # Get amplicons from forward and reverse strand
#  my $fwd_amplicons = _extract_amplicons_from_strand($seq, $fwd_regexp, $rev_regexp, 1);
#  my $rev_amplicons = _extract_amplicons_from_strand($seq, $fwd_regexp, $rev_regexp, -1);
#
#  # Deal with nested amplicons by removing the longest of the two
#  my $re = qr/(\d+)\.\.(\d+)/;
#  for (my $rev = 0; $rev < scalar @$rev_amplicons; $rev++) {
#    my ($rev_start, $rev_end) = ( $rev_amplicons->[$rev]->{_amplicon} =~ m/$re/ );
#    for (my $fwd = 0; $fwd < scalar @$fwd_amplicons; $fwd++) {
#      my ($fwd_start, $fwd_end) = ( $fwd_amplicons->[$fwd]->{_amplicon} =~ m/$re/ );
#      if ( ($fwd_start < $rev_start) && ($rev_end < $fwd_end) ) {
#        splice @$fwd_amplicons, $fwd, 1; # Remove forward amplicon
#        $fwd--;
#        next;
#      }
#      if ( ($rev_start < $fwd_start) && ($fwd_end < $rev_end) ) {
#        splice @$rev_amplicons, $rev, 1; # Remove reverse amplicon
#        $rev--;
#      }
#    }
#  }
#  
#  my $amplicons = [ @$fwd_amplicons, @$rev_amplicons ];

   return $amplicon;
}


sub _attach_amplicon {
   # Create an amplicon object and attach it to template
   my ($self, $start, $end, $strand) = @_;

   if ($strand == -1) {
      # Calculate coordinates relative to forward strand. For example, given a
      # read starting at 10 and ending at 23 on the reverse complement of a 100
      # bp, give a strand or 77 and end of 90 on forward strand.
      my $length = length $template_str;
      $start = $length - $start + 1;
      $end   = $length - $end + 1;
      ($start, $end) = ($end, $start);
   }

   # Create Bio::SeqFeature::Amplicon feature and attach it to the template
   my $amplicon = Bio::SeqFeature::Amplicon->new(
      -start    => $start,
      -end      => $end,
      -strand   => $strand,
      -template => $self->template,
   );

   ####
   use Data::Dumper;
   print "Amplicon: ".$amplicon->start." .. ".$amplicon->end." -> ".$amplicon->seq->seq."\n";
   print "AMPLICON: ".Dumper($amplicon);
   ####

   # Create Bio::SeqFeature::Primer feature and attach them to the amplicon
   if ($self->attach_primers) {
      for my $type ('fwd', 'rev') {
         my ($pstart, $pend, $pstrand, $primer_seq);
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

         ####
         $pstart += $start - 1;
         $pend   += $start - 1;
         ####

         Bio::SeqFeature::Primer->new(
            -start    => $pstart,
            -end      => $pend,
            -strand   => $pstrand,
            -template => $amplicon,
         );
      }
   }

   ####
   use Data::Dumper;
   print "Amplicon: ".$amplicon->start." .. ".$amplicon->end." -> ".$amplicon->seq->seq."\n";
   print "AMPLICON: ".Dumper($amplicon);
   ####

   return $amplicon;
}


=head2 annotate_template

 Title    : annotate_template
 Usage    : my $template = $search->annotate_template;
 Function : Search for all amplicons and attach them to the template.
            This is equivalent to running:
               XXX
               
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

1;
