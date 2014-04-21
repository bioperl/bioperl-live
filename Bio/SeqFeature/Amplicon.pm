#
# BioPerl module for Bio::SeqFeature::Amplicon
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Copyright Florent Angly
#
# You may distribute this module under the same terms as perl itself


=head1 NAME

Bio::SeqFeature::Amplicon - Amplicon feature

=head1 SYNOPSIS

  # Amplicon with explicit sequence
  use Bio::SeqFeature::Amplicon;
  my $amplicon = Bio::SeqFeature::Amplicon->new( 
      -seq        => $seq_object,
      -fwd_primer => $primer_object_1,
      -rev_primer => $primer_object_2,
  );

  # Amplicon with implicit sequence
  use Bio::Seq;
  my $template = Bio::Seq->new( -seq => 'AAAAACCCCCGGGGGTTTTT' );
  $amplicon = Bio::SeqFeature::Amplicon->new(
      -start => 6,
      -end   => 15,
  );
  $template->add_SeqFeature($amplicon);
  print "Amplicon start   : ".$amplicon->start."\n";
  print "Amplicon end     : ".$amplicon->end."\n";
  print "Amplicon sequence: ".$amplicon->seq->seq."\n";
  # Amplicon sequence should be 'CCCCCGGGGG'

=head1 DESCRIPTION

Bio::SeqFeature::Amplicon extends L<Bio::SeqFeature::Subseq> to represent an
amplicon sequence and optional primer sequences.

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
the bugs and their resolution.  Bug reports can be submitted via 
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Florent Angly <florent.angly@gmail.com>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


package Bio::SeqFeature::Amplicon;

use strict;

use base qw(Bio::SeqFeature::SubSeq);

=head2 new

 Title   : new()
 Usage   : my $amplicon = Bio::SeqFeature::Amplicon( -seq => $seq_object );
 Function: Instantiate a new Bio::SeqFeature::Amplicon object
 Args    : -seq        , the sequence object or sequence string of the amplicon (optional)
           -fwd_primer , a Bio::SeqFeature primer object with specified location on amplicon (optional)
           -rev_primer , a Bio::SeqFeature primer object with specified location on amplicon (optional)
 Returns : A Bio::SeqFeature::Amplicon object

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($fwd_primer, $rev_primer) =
        $self->_rearrange([qw(FWD_PRIMER REV_PRIMER)], @args);
    $fwd_primer && $self->fwd_primer($fwd_primer);
    $rev_primer && $self->rev_primer($rev_primer);
    return $self;
}  


sub _primer {
    # Get or set a primer. Type is either 'fwd' or 'rev'.
    my ($self, $type, $primer) = @_;
    if (defined $primer) {
        if ( not(ref $primer) || not $primer->isa('Bio::SeqFeature::Primer') ) {
            $self->throw("Expected a primer object but got a '".ref($primer)."'\n");
        }
        if ( not defined $self->location ) {
            $self->throw("Location of $type primer on amplicon is not known. ".
                "Use start(), end() or location() to set it.");
        }
        $primer->primary_tag($type.'_primer');
        $self->add_SeqFeature($primer);
    }
    return (grep { $_->primary_tag eq $type.'_primer' } $self->get_SeqFeatures)[0];
}    


=head2 fwd_primer

 Title   : fwd_primer
 Usage   : my $primer = $feat->fwd_primer();
 Function: Get or set the forward primer. When setting it, the primary tag
           'fwd_primer' is added to the primer and its start, stop and strand
           attributes are set if needed, assuming that the forward primer is
           at the beginning of the amplicon and the reverse primer at the end.
 Args    : A Bio::SeqFeature::Primer object (optional)
 Returns : A Bio::SeqFeature::Primer object

=cut

sub fwd_primer {
    my ($self, $primer) = @_;
    return $self->_primer('fwd', $primer);
}


=head2 rev_primer

 Title   : rev_primer
 Usage   : my $primer = $feat->rev_primer();
 Function: Get or set the reverse primer. When setting it, the primary tag
          'rev_primer' is added to the primer.
 Args    : A Bio::SeqFeature::Primer object (optional)
 Returns : A Bio::SeqFeature::Primer object

=cut

sub rev_primer {
    my ($self, $primer) = @_;
    return $self->_primer('rev', $primer);
}


1;
