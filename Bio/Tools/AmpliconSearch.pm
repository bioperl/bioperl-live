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
complement) of the specified template sequence When two amplicons overlap, an
option allows to discard the longest one to more accurately represent the biases
of PCR. Future improvements may include modelling the effects of the number of
PCR cycles or temperature on the PCR products.

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
 Args     : -template        Sequence object. Note that the template sequence will be converted to Bio::Seq if needed. Features (amplicons and primers will be added to this object).
            -forward_primer
            -reverse_primer
            -primer_file    replaces -forward_primer and -reverse_primer (optional)
            -attach_primers whether or not to attach primers to Amplicon objects. Default: 0 (off)
 Returns  : A Bio::Tools::AmpliconSearch object

=cut

sub new {
   my ($class, @args) = @_;
   my $self = $class->SUPER::new(@args);
   my ($template, $primer_file, $forward_primer, $reverse_primer,
       $attach_primers) = $self->_rearrange([qw(TEMPLATE PRIMER_FILE
       FORWARD_PRIMER REVERSE_PRIMER ATTACH_PRIMERS)], @args);

   # Get primers
   if (defined $primer_file) {
      ($forward_primer, $reverse_primer) = $self->_get_primers_from_file($primer_file);
   }

   $self->_set_forward_primer($forward_primer) if defined $forward_primer; ###
   $self->_set_reverse_primer($reverse_primer);

   # Get template sequence
   $self->_set_template($template) if defined $template;

   $self->_set_attach_primers($attach_primers) if defined $attach_primers;

   return $self;
}


=head2 template

 Title    : template
 Usage    : my $template = $search->template;
 Function : Get the template sequence
 Args     :                   
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
   #if ($self->attach_primer && not $template->isa('Bio::SeqI')) {
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


sub strand {
   my ($self) = @_;
   return $self->{strand};
}

sub _set_strand {
   my ($self, $strand) = @_;
   $self->{strand} = $strand;
   return $self->strand;
}


sub forward_primer {
   my ($self) = @_;
   return $self->{forward_primer};
}

sub _set_forward_primer {
   my ($self, $primer) = @_;
   if (not(ref $primer) || not $primer->isa('Bio::PrimarySeqI') || not $primer->isa('Bio::SeqFeature::Primer') ) { 
      # Not a sequence or a primer object
      $self->throw("Expected a sequence or primer object as input but got a ".ref($primer)."\n");
   }
   $self->{forward_primer} = $primer;
   $self->_set_forward_regexp( Bio::Tools::IUPAC->new( -seq => $primer )->regexp );
   return $self->forward_primer;
}


sub reverse_primer {
   my ($self) = @_;
   return $self->{reverse_primer};
}

sub _set_reverse_primer {
   my ($self, $primer) = @_;
   my $re;

   # Set the reverse primer
   if (defined $primer) {
      if (not(ref $primer) || not $primer->isa('Bio::PrimarySeqI') || not $primer->isa('Bio::SeqFeature::Primer') ) { 
         # Not a sequence or a primer object
         $self->throw("Expected a sequence or primer object as input but got a ".ref($primer)."\n");
      }
      $self->{reverse_primer} = $primer;
      $re = Bio::Tools::IUPAC->new( -seq => $primer->revcom )->regexp;
   }

   # No reverse primer given, match end of string
   else {
      $re = qr/$/;
   }

   $self->_set_reverse_regexp($re);
   return $self->reverse_primer;
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


sub forward_regexp {
   my ($self) = @_;
   return $self->{forward_regexp};
}


sub _set_forward_regexp {
   my ($self, $regexp) = @_;
   $self->{forward_regexp} = $regexp;
   return $self->forward_regexp;
}


sub reverse_regexp {
   my ($self) = @_;
   return $self->{reverse_regexp};
}


sub _set_reverse_regexp {
   my ($self, $regexp) = @_;
   $self->{reverse_regexp} = $regexp;
   return $self->reverse_regexp;
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


=head2 next_amplicon

 Title    : next_amplicon
 Usage    : my $amplicon = $search->next_amplicon;
 Function : Get the next amplicon
 Args     : None
 Returns  : A Bio::SeqFeature::Amplicon object

=cut

sub next_amplicon {
   my ($self, @args) = @_;
   my $amplicon;

   my $strand = $self->strand;
   my $fwd_regexp = $self->forward_regexp ||
      $self->throw("Need to provide at least a primer\n");
   my $rev_regexp = $self->reverse_regexp;

   if ($template_str  =~ m/($fwd_regexp.*?$rev_regexp)/g) {
      my $end   = pos($template_str);
      my $start = $end - length($1) + 1;
      # Now trim the left end to obtain the shortest amplicon
      my $ampliconstr = substr $template_str, $start - 1, $end - $start + 1;
      if ($ampliconstr =~ m/$fwd_regexp.*($fwd_regexp)/g) {
         $start += pos($ampliconstr) - length($1);
      }
      $amplicon = $self->_create_amplicon($start, $end, $strand);
   }

   if ( ($strand == 1) && not($amplicon) ) {
      # No more matches in the forward strand. Search the reverse-complement.
      $template_str = $self->template->revcom->seq;
      $self->_set_strand(-1);
      $amplicon = $self->next_amplicon;
   }

   return $amplicon;
}


sub _create_amplicon {
   # Create an amplicon sequence and register its coordinates
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

   # Create Bio::SeqFeature::Primer feature and attach them to the amplicon
   if ($self->attach_primers) {
      for my $type ('fwd', 'rev') {
         my ($pstart, $pend, $pstrand);

         if ($type eq 'fwd') {
            # Forward primer
            $pstart  = 1;
            $pend    = $self->forward_primer->length;
            $pstrand = $amplicon->strand;
         }
         
         else {
            # Optional reverse primer
            my $primer_seq = $self->reverse_primer;
            next if not defined $primer_seq;
            $pstart  = $end - $primer_seq->length + 1;
            $pend    = $end;
            $pstrand = -1 * $amplicon->strand;
         }

         my $primer = Bio::SeqFeature::Primer->new(
            -start    => $pstart,
            -end      => $pend,
            -strand   => $pstrand,
            -template => $amplicon,
         );
      }
   }

   return $amplicon;
}


1;
