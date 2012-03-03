# BioPerl module for Bio::Tools::AmpliconSearch
#
# Copyright Florent Angly
#
# You may distribute this module under the same terms as perl itself


package Bio::Tools::AmpliconSearch;

use strict;
use warnings;
use Bio::Tools::IUPAC;

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
 Args     : -template
            -forward_primer
            -reverse_primer
            -primer_file    replaces -forward_primer and -reverse_primer (optional)
            -want_a         type of object to have next_amplicon() return: PrimarySeq, LocatableSeq, Amplicon
            -attach_primers when returning amplicon objects as Amplicons, whether or not to attach primers to it
 Returns  : 

=cut

sub new {
   my ($class, @args) = @_;
   my $self = $class->SUPER::new(@args);
   my ($template, $primer_file, $forward_primer, $reverse_primer, $want_a,
       $attach_primers) = $self->_rearrange([qw(TEMPLATE PRIMER_FILE
       FORWARD_PRIMER REVERSE_PRIMER WANT_A ATTACH_PRIMERS)], @args);

   if (defined $primer_file) {
      $self->_set_primers($primer_file);
   } else {
      $self->_set_forward_primer($forward_primer) if defined $forward_primer;
      $self->_set_reverse_primer($reverse_primer) if defined $reverse_primer;
   }

   $self->_set_template($template) if defined $template;
   $self->_set_want_a($want_a) if defined $want_a;
   $self->_set_attach_primers($attach_primers) if defined $attach_primers;

   return $self;
}


sub get_template {
   my ($self) = @_;
   return $self->{'template'};
}


sub _set_template {
   my ($self, $val) = @_;
   $self->{'template'} = $val;
   if (not $val->isa('Bio::PrimarySeqI') ) { 
      # not a Bio::Seq or Bio::PrimarySEq
      $self->throw("Expected a sequence object as input but got a ".ref($val)."\n");
   }
   return $self->get_template;
}


sub _set_primers {
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


sub get_forward_primer {
   my ($self) = @_;
   return $self->{'forward_primer'};
}


sub _set_forward_primer {
   my ($self, $val) = @_;
   $self->{'forward_primer'} = $val;
   if (not $val->isa('Bio::PrimarySeqI') || not $val->isa('Bio::SeqFeature::Primer') ) { 
      # Not a sequence or a primer object
      $self->throw("Expected a sequence or primer object as input but got a ".ref($val)."\n");
   }
   $self->_set_forward_regexp( Bio::Tools::IUPAC->new( -seq => $val )->regexp );
   return $self->get_forward_primer;
}


sub get_reverse_primer {
   my ($self) = @_;
   return $self->{'reverse_primer'};
}


sub _set_reverse_primer {
   my ($self, $val) = @_;
   $self->{'reverse_primer'} = $val;
   if (not $val->isa('Bio::PrimarySeqI') || not $val->isa('Bio::SeqFeature::Primer') ) { 
      # Not a sequence or a primer object
      $self->throw("Expected a sequence or primer object as input but got a ".ref($val)."\n");
   }
   $self->_set_reverse_regexp( Bio::Tools::IUPAC->new( -seq => $val->revcom )->regexp );
   return $self->get_reverse_primer;
}


sub get_forward_regexp {
   my ($self) = @_;
   return $self->{'forward_regexp'};
}


sub _set_forward_regexp {
   my ($self, $val) = @_;
   $self->{'forward_regexp'} = $val;
   return $self->get_forward_regexp;
}


sub get_reverse_regexp {
   my ($self) = @_;
   return $self->{'reverse_regexp'};
}


sub _set_reverse_regexp {
   my ($self, $val) = @_;
   $self->{'reverse_regexp'} = $val;
   return $self->get_reverse_regexp;
}


sub get_want_a {
   my ($self) = @_;
   return $self->{'want_a'} || 'Bio::PrimarySeq';
}


sub _set_want_a {
   my ($self, $val) = @_;
   $self->{'want_a'} = $val;
   return $self->get_want_a;
}


sub get_attach_primers {
   my ($self) = @_;
   return $self->{'attach_primers'} || 0;
}


sub _set_attach_primers {
   my ($self, $val) = @_;
   $self->{'attach_primers'} = $val;
   return $self->get_attach_primers;
}


=head2 next_amplicon

 Title    : next_amplicon
 Usage    : my $amplicon = $search->
 Function : Get the next amplicon
 Args     : -want_a :  Bio::PrimarySeq (default), Bio::LocatableSeq or Bio::SeqFeature::Amplicon
            -attach_primers : 0 (no, default) or 1 (yes). Only possible is -want_a => 'Bio::SeqFeature::Amplicon'
 Returns  : an amplicon object

=cut

sub next_amplicon {
   my ($self, @args) = @_;

   my $amplicon;

   my $seq = $self->get_template;
   my $seqstr = $seq->seq;
   my $fwd_regexp = $self->get_forward_regexp;
   my $rev_regexp = $self->get_reverse_regexp;
 
   #### orientation?
   my $orientation = 1;

   if ( defined($fwd_regexp) && not(defined $rev_regexp) ) {
       # From forward primer to end of template
       $seqstr   =~ m/($fwd_regexp)/g;
       my $start = pos($seqstr) - length($1) + 1;
       my $end   = length($seqstr); ### $seq->length;
       $amplicon = $self->_create_amplicon($start, $end, $orientation);

   } elsif ( defined($fwd_regexp) && defined($rev_regexp) ) {
       # From forward to reverse primer
       $seqstr   =~ m/($fwd_regexp.*?$rev_regexp)/g;
       my $end   = pos($seqstr);
       my $start = $end - length($1) + 1;
       # Now trim the left end to obtain the shortest amplicon
       my $ampliconstr = substr $seqstr, $start - 1, $end - $start + 1;
       if ($ampliconstr =~ m/$fwd_regexp.*($fwd_regexp)/g) {
          $start += pos($ampliconstr) - length($1);
       }
       $amplicon = $self->_create_amplicon($start, $end, $orientation);

   } else {
      $self->throw("Need to provide at least a forward primer\n");
   }

   #### return an amplicon? primedseq?
   return 1;
}


sub _create_amplicon {
   # Create an amplicon sequence and register its coordinates
   my ($self, $start, $end, $orientation) = @_;

   my $want_a = $self->want_a;

   my $template = $self->get_template;
   my $amplicon;
   my $coord;

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

  return $amplicon
}


1;
