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

When two amplicons overlap, an option allows to discard the longest one to more
accurately represent the biases of PCR. Future improvements may include
modelling the effects of the number of PCR cycles or temperature on the PCR
products.

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
      $self->_set_primers_from_file($primer_file);
   } else {
      $self->_set_forward_primer($forward_primer) if defined $forward_primer;
      $self->_set_reverse_primer($reverse_primer) if defined $reverse_primer;
   }

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
   return $self->template;
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
   if (not(ref $primer) || not $primer->isa('Bio::PrimarySeqI') || not $primer->isa('Bio::SeqFeature::Primer') ) { 
      # Not a sequence or a primer object
      $self->throw("Expected a sequence or primer object as input but got a ".ref($primer)."\n");
   }
   $self->{reverse_primer} = $primer;
   $self->_set_reverse_regexp( Bio::Tools::IUPAC->new( -seq => $primer->revcom )->regexp );
   return $self->reverse_primer;
}


sub _set_primers_from_file {
   my ($self, $primer_file) = @_;
   # Read primer file and convert primers into regular expressions to catch
   # amplicons present in the database

   if (not defined $primer_file) {
      $self->throw("Need to provide an input file\n");
   }

   # Mandatory first primer
   require Bio::SeqIO;
   my $in = Bio::SeqIO->newFh( -file => $primer_file );
   my $primer = <$in>;
   if (not defined $primer) {
      $self->throw("The file '$primer_file' contains no primers\n");
   }
   $primer->alphabet('dna'); # Force the alphabet since degenerate primers can look like protein sequences
   $self->_set_forward_primer($primer);

   # Optional reverse primers
   $primer = undef;
   $primer = <$in>;
   if (defined $primer) {
      $primer->alphabet('dna');
      $self->_set_reverse_primer($primer);
   }
   
   #### $in->close;
   #### close $in;
   undef $in;

   return 1;
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
 Usage    : my $amplicon = $search->
 Function : Get the next amplicon
 Args     :
 Returns  : an amplicon object

=cut

sub next_amplicon {
   my ($self, @args) = @_;

   my $amplicon;

   my $seqstr = $self->template->seq;
   my $fwd_regexp = $self->forward_regexp;
   my $rev_regexp = $self->reverse_regexp;
 
   #### orientation? need to revcom the template if we have exhausted the fwd orientation
   my $strand = 1;

   if ( defined($fwd_regexp) && not(defined $rev_regexp) ) {
      # From forward primer to end of template
      if ($seqstr  =~ m/($fwd_regexp)/g) {
         my $start = pos($seqstr) - length($1) + 1;
         my $end   = length($seqstr);
         $amplicon = $self->_create_amplicon($start, $end, $strand);
      }

   } elsif ( defined($fwd_regexp) && defined($rev_regexp) ) {
      # From forward to reverse primer
#  while ( $seqstr =~ m/($fwd_regexp.*?$rev_regexp)/g ) {
#      if ($seqstr  =~ m/($fwd_regexp.*?$rev_regexp)/g) {
      if ($seqstr  =~ m/\G.*?($fwd_regexp.*?$rev_regexp)/g) {
 
         #####
         print "match: '$1'\n";
         #####

         my $end   = pos($seqstr);
         my $start = $end - length($1) + 1;

         ####
         pos($seqstr) = $start + 2;
         ####

         # Now trim the left end to obtain the shortest amplicon
         my $ampliconstr = substr $seqstr, $start - 1, $end - $start + 1;
         if ($ampliconstr =~ m/$fwd_regexp.*($fwd_regexp)/g) {
            $start += pos($ampliconstr) - length($1);
         }
         $amplicon = $self->_create_amplicon($start, $end, $strand);
      }

   } else {
      $self->throw("Need to provide at least a forward primer\n");
   }

   return $amplicon;
}


sub _create_amplicon {
   # Create an amplicon sequence and register its coordinates
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

#   my $coord;

#   if ($want_a eq 'Bio::PrimarySeq') {
#      require Bio::PrimarySeq;
#      $amplicon = Bio::PrimarySeq->new( -seq => XXX );
#   } elsif ($want_a eq 'Bio::LocatableSeq') {
#      require Bio::LocatableSeq;
#      $amplicon = Bio::PrimarySeq->new( -seq => XXX );
#   } elsif ($want_a eq 'Bio::SeqFeature::Amplicon') {
#      require Bio::SeqFeature::Primer;
#      $amplicon = Bio::PrimarySeq->new( -start => , -end => , -strand => );
#   } else {
#      $self->throw("'$want_a' is not a valid class to use to return amplicons.");
#   }

   #$amplicon = $want_a->new( -seq => );
#   if (not $want_a eq 'Bio::PrimarySeq') {
#      # add coordinates
#      $amplicon->start( XXX );
#      $amplicon->end( XXX );
#      $amplicon->strand( XXX );
#      if ( ($want_a eq 'Bio::SeqFeature::Amplicon') && ($self->attach_primers) ) {
#         # add primers to
#         $amplicon->fwd_primer( XXX );
#         $amplicon->rev_primer( XXX );
#      }
#   }

  
  #my $fwd_primer = Bio::SeqFeature::Primer(
  #    -seq      => $seq_object,
  #    -location =>
  #    -strand   =>
  #);
  #my $rev_primer = Bio::SeqFeature::Primer(
  #    -seq      => $seq_object,
  #    -location =>
  #    -strand   =>
  #);

  # Should be able to give a Bio::SeqFeature::Primer or simply a Bio::Seq (but
  # may need to ensure that reverse primer has reverse orientation??)
  #my $primedseq = Bio::Seq::PrimedSeq->new(-seq => $template, 
  #                                         -left_primer => $fwd_primer,
  #                                         -right_primer => $rev_primer);

  #if ($orientation == -1) {
  #  # Calculate coordinates relative to forward strand. For example, given a
  #  # read starting at 10 and ending at 23 on the reverse complement of a 100 bp
  #  # sequence, return complement(77..90).
  #  $amplicon = $seq->revcom->trunc($start, $end);
  #  my $seq_len = $seq->length;
  #  $start = $seq_len - $start + 1;
  #  $end   = $seq_len - $end + 1;
  #  ($start, $end) = ($end, $start);
  #  $coord = "complement($start..$end)";
  #} else {
  #  $amplicon = $seq->trunc($start, $end);
  #  $coord = "$start..$end";
  #}
  #$amplicon->{_amplicon} = $coord;

  return $amplicon;
}


1;
