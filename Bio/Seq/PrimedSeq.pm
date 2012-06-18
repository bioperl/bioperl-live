#
# This is the original copyright statement. I have relied on Chad's module
# extensively for this module.
#
# Copyright (c) 1997-2001 bioperl, Chad Matsalla. All Rights Reserved.
#           This module is free software; you can redistribute it and/or
#           modify it under the same terms as Perl itself. 
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code
#
# But I have modified lots of it, so I guess I should add:
#
# Copyright (c) 2003 bioperl, Rob Edwards. All Rights Reserved.
#           This module is free software; you can redistribute it and/or
#           modify it under the same terms as Perl itself. 
#
# Copyright Rob Edwards
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code


=head1 NAME

Bio::Seq::PrimedSeq - A sequence and a pair of primers matching on it

=head1 SYNOPSIS

  use Bio::Seq;
  use Bio::Seq::PrimedSeq;

  my $template = Bio::Seq->new( -seq => 'AGCTTTTCATTCTGACTGCAAC' );
  my $left     = Bio::Seq->new( -seq => 'AGCT' );
  my $right    = Bio::Seq->new( -seq => 'GTTGC' );

  my $primedseq = Bio::Seq::PrimedSeq->new(
          -seq          => $template,  # sequence object
          -left_primer  => $left,      # sequence or primer object
          -right_primer => $right      # sequence or primer object
  );

  # Get the primers as Bio::SeqFeature::Primer objects
  my @primer_objects = $primedseq->get_primer();

  # Sequence object representing the PCR product, i.e. the section of the target
  # sequence contained between the primers (primers included)
  my $amplicon_seq = $primedseq->amplicon();

  print 'Got amplicon sequence: '.$amplicon_seq->seq."\n";
  # Amplicon should be: agctTTTCATTCTGACTgcaac

=head1 DESCRIPTION

This module was created to address the fact that a primer is more than a
SeqFeature and that there is a need to represent the primer-sequence complex and
the attributes that are associated with the complex.

A PrimedSeq object is initialized with a target sequence and two primers. The
location of the primers on the target sequence is calculated so that one can
then get the PCR product, or amplicon sequence. From the PrimedSeq object you
can also retrieve information about melting temperatures and what not on each
of the primers and the amplicon. The L<Bio::Tools::Primer3> uses PrimedSeq
objects extensively.

Note that this module does not simulate PCR. It assumes that the primers will
match in a single location on the target sequence and does not understand
degenerate primers.

The L<Bio::Seq::PrimedSeq> was initially started by Chad Matsalla, and later
improved on by Rob Edwards.

=head1 RECIPES

=over

=item 1.

Run Primer3 to get PrimedSeq objects:

  use Bio::SeqIO;
  use Bio::Tools::Run::Primer3;

  # Read a target sequences from a FASTA file
  my $file = shift || die "Need a file to read";
  my $seqin = Bio::SeqIO->new(-file => $file);
  my $seq = $seqin->next_seq;

  # Use Primer3 to design some primers
  my $primer3 = Bio::Tools::Run::Primer3->new(-seq => $seq);
  my $results = $primer3->run; # default parameters

  # Write all the results in a Genbank file
  my $seqout = Bio::SeqIO->new(-file => ">primed_sequence.gbk", 
                               -format => 'genbank');
  while (my $primedseq = $results->next_primer) {
     $seqout->write_seq( $primedseq->annotated_seq );
  }

=item 2.

Create a genbank file for a sequence and its cognate primers:

  use Bio::SeqIO;
  use Bio::Seq::PrimedSeq;

  # Read a FASTA file that contains the target sequence and its two primers
  my $file = shift || die "$0 <file>";
  my $seqin = Bio::SeqIO->new(-file => $file);
  my ($template, $leftprimer, $rightprimer) = 
        ($seqin->next_seq, $seqin->next_seq, $seqin->next_seq);

  # Set up a PrimedSeq object
  my $primedseq = Bio::Seq::PrimedSeq->new(-seq => $template, 
                                           -left_primer => $leftprimer,
                                           -right_primer => $rightprimer);

  # Write the sequences in an output Genbank file
  my $seqout = Bio::SeqIO->new(-file => ">primed_sequence.gbk",
                               -format => 'genbank');
  $seqout->write_seq($primedseq->annotated_sequence);

=back

=head1 FEEDBACK

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

Rob Edwards, redwards@utmem.edu

Based on a module written by Chad Matsalla, bioinformatics1@dieselwurks.com

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Seq::PrimedSeq;

use strict;
use Bio::SeqFeature::Primer;

use base qw(Bio::Root::Root Bio::SeqFeature::Generic);


=head2 new

 Title   : new()
 Usage   : my $primedseq = Bio::SeqFeature::Primer->new( 
                            -seq => $sequence,
                            -left_primer => $left_primer,
                            -right_primer => $right_primer
           );
 Function: Construct a primed sequence.
 Returns : A Bio::Seq::PrimedSeq object
 Args    :  -seq => a Bio::Seq object (required)
            -left_primer => a Bio::SeqFeature::Primer or sequence object (required)
            -right_primer => a Bio::SeqFeature::Primer or sequence object (required)

           If you pass a sequence object to specify a primer, it will be used to
           construct a Bio::SeqFeature::Primer that you can retrieve with the
           L<get_primer> method.

           Many other parameters can be included including all of the output
           parameters from the primer3 program (see L<Bio::Tools::Primer3>). At
           the moment these parameters will simply be stored and do anything.

=cut

sub new {
    my($class,%args) = @_;
    my $self = $class->SUPER::new(%args);

    # Need an amplicon sequence
    $self->{seq} = delete $args{-seq} || delete $args{-target_sequence} ||
        $self->throw("Need to provide a sequence during PrimedSeq object construction");
    if (! ref($self->{seq}) || ! $self->{seq}->isa('Bio::SeqI') ) {
        $self->throw("The target_sequence must be a Bio::Seq to create this object.");
    }

    # Need a left and right primers
    for my $primer ( 'left', 'right' ) {
        $self->{$primer} = delete $args{'-'.$primer.'_primer'} ||
            $self->throw("Need to provide both primers during PrimedSeq object construction");
        if ( ref $self->{$primer} && $self->{$primer}->isa('Bio::PrimarySeqI') ) {
            # Convert Bio::Seq or Bio::PrimarySeq objects to Bio::SeqFeature::Primer
            $self->{$primer} = Bio::SeqFeature::Primer->new(-seq => $self->{$primer});
        }
        if (not (ref $self->{$primer} && $self->{$primer}->isa("Bio::SeqFeature::Primer"))) {
            $self->throw("Primers must be Bio::SeqFeature::Primer objects but got a ".ref($self->{$primer}));
        }
    }

    # Save remaining arguments
    while (my ($arg, $val) = each %args) {
        $self->{$arg} = $val;
    }
 
    # Now we have the sequences, let's find out where they are
    $self->_place_primers();

    return $self;
}


=head2 get_primer

 Title   : get_primer();
 Usage   :  my @primers = $primedseq->get_primer();
              or
            my $primer = $primedseq->get_primer('-left_primer');
 Function: Get the primers associated with the PrimedSeq object.
 Returns : A Bio::SeqFeature::Primer object
 Args    : For the left primer, use: l, left, left_primer or -left_primer
           For the right primer, use: r, right, right_primer or -right_primer
           For both primers [default], use: b, both, both_primers or -both_primers

=cut

sub get_primer {
    my ($self, $arg) = @_;
    if (! defined $arg ) {
        return ($self->{'left'}, $self->{'right'});
    } elsif ( $arg =~ /^-?l/ ) { 
        # What a cheat, I couldn't be bothered to write all the exact statements!
        # Hah, now you can write 'leprechaun' to get the left primer.
        return $self->{'left'}
    } elsif ( $arg =~ /^-?r/ ) {
        return $self->{'right'};
    } elsif ( $arg =~ /^-?b/ ) {
        return ($self->{'left'}, $self->{'right'});
    }
}


=head2 annotated_sequence

 Title   : annotated_sequence
 Usage   : my $annotated_sequence_object = $primedseq->annotated_sequence();
 Function: Get an annotated sequence object containg the left and right 
           primers
 Returns : An annotated sequence object or 0 if not defined.
 Args    : 
 Note    : Use this method to return a sequence object that you can write
           out (e.g. in GenBank format). See the example above.

=cut

sub annotated_sequence {
    my $self = shift;
    if (not exists $self->{annotated_sequence}) {
        $self->_set_seqfeature;
    }
    return $self->{annotated_sequence} || 0;
}


=head2 amplicon

 Title   : amplicon
 Usage   : my $amplicon = $primedseq->amplicon();
 Function: Retrieve the amplicon as a sequence object. The amplicon sequnce is
           only the section of the target sequence between the primer matches
           (primers included).
 Returns : A Bio::Seq object. To get the sequence use $amplicon->seq
 Args    : None
 Note    : 

=cut

sub amplicon {
    my ($self) = @_;
    my $seqobj = Bio::Seq->new(
        -id  => 'Amplicon_from_'.($self->seq->id || 'unidentified'),
        -seq => $self->{'amplicon_sequence'},
    );
    return $seqobj;
}


=head2 seq

 Title   : seq
 Usage   : my $seqobj = $primedseq->seq();
 Function: Retrieve the target sequence as a sequence object
 Returns : A seq object. To get the sequence use $seqobj->seq
 Args    : None
 Note    : 

=cut

sub seq {
    my $self = shift;
    return $self->{seq};
}


=head2 _place_primers

 Title   : _place_primers
 Usage   : $self->_place_primers();
 Function: An internal method to place the primers on the sequence and 
           set up the ranges of the sequences
 Returns : Nothing
 Args    : None
 Note    : Internal use only

=cut

sub _place_primers {
    my $self = shift;

    # Get the target and primer sequence strings, all in uppercase
    my $target_seq = uc $self->{'seq'}->seq();
    my $left_seq   = uc $self->{'left'}->seq()->seq();
    my $right_seq  = uc $self->{'right'}->seq()->revcom()->seq();

    # Locate primers on target sequence
    my ($before, $middle, $after) = ($target_seq =~ /^(.*)$left_seq(.*)$right_seq(.*)$/);
    if (not defined $before || not defined $middle || not defined $after) {
        if ($target_seq !~ /$left_seq/) {
            $self->throw("Could not place left primer on target");
        }
        if ($target_seq !~ /$right_seq/) {
            $self->throw("Could not place right primer on target")
        }
    }

    # Bio::Seq::PrimedSeq positions
    my $left_location  = length($before).','.length($left_seq);
    my $right_location = (length($target_seq) - length($after) - 1).','.length($right_seq);
    my $amplicon_size  = length($left_seq) + length($middle) + length($right_seq);
 
    # Primer3 positions
    my $primer_left         = $self->{'left'}->{'PRIMER_LEFT'};
    my $primer_right        = $self->{'right'}->{'PRIMER_RIGHT'};
    my $primer_product_size = $self->{'PRIMER_PRODUCT_SIZE'};

    # Compare Bio::Seq::PrimedSeq positions to Primer3 positions
    if (defined $primer_left) {
        if (not $primer_left eq $left_location) {
            $self->warn("Note got |$primer_left| from primer3 and |$left_location| for the left primer.");
        }
    } else {
        $self->{'left'}->{'PRIMER_LEFT'} = $left_location;
    }

    if (defined $primer_right) {
        if (not $primer_right eq $right_location) {
            $self->warn("Note got |$primer_right| from primer3 and |$right_location| for the right primer.");
        }
    } else {
        $self->{'right'}->{'PRIMER_RIGHT'} = $right_location;
    }

    if (defined $primer_product_size) {
        if (not $primer_product_size eq $amplicon_size) {
            $self->warn("Note got |$primer_product_size| from primer3 and |$amplicon_size| for the size.");
        }
    } else {
        $self->{'PRIMER_PRODUCT_SIZE'} = $amplicon_size;
    }

    $self->{'amplicon_sequence'} = lc($left_seq).uc($middle).lc($right_seq);

}


=head2 _set_seqfeature

 Title   : _set_seqfeature
 Usage   : $self->_set_seqfeature()
 Function: An internal method to create Bio::SeqFeature::Generic objects
           for the primed seq
 Returns : Nothing
 Args    : None
 Note    : Internal use only. Should only call this once left and right 
           primers have been placed on the sequence. This will then set 
           them as sequence features so hopefully we can get a nice output 
           with write_seq.

=cut

sub _set_seqfeature {
    my $self = shift;
    unless ($self->{'left'}->{'PRIMER_LEFT'} && 
              $self->{'right'}->{'PRIMER_RIGHT'}) {
        $self->warn('Hmmm. Haven\'t placed primers, but trying to make annotated sequence');
        return 0;
    }
    my ($start, $length) = split /,/, $self->{'left'}->{'PRIMER_LEFT'};
    my $tm = $self->{'left'}->{'PRIMER_LEFT_TM'} || $self->{'left'}->Tm || 0;

    my $seqfeatureL = Bio::SeqFeature::Generic->new(
        -start       => $start + 1,
        -end         => $start + $length,
        -strand      => 1,
        -primary_tag => 'left_primer',
        -source      => 'primer3',
        -tag         => {new => 1, author => 'Bio::Seq::PrimedSeq', Tm => $tm},
    );
 
    ($start, $length) = split /,/, $self->{'right'}->{'PRIMER_RIGHT'};
    $tm = $self->{'right'}->{'PRIMER_RIGHT_TM'} || $self->{'right'}->Tm || 0;
 
    my $seqfeatureR = Bio::SeqFeature::Generic->new(
        -start       => $start - $length + 2,
        -end         => $start + 1,
        -strand      => -1,
        -primary_tag => 'right_primer',
        -source      => 'primer3',
        -tag         => {new => 1, author => 'Bio::Seq::PrimedSeq', Tm => $tm},
    );

    # Now add the sequences to an annotated sequence
    $self->{annotated_sequence} = $self->{seq};
    $self->{annotated_sequence}->add_SeqFeature($seqfeatureL);
    $self->{annotated_sequence}->add_SeqFeature($seqfeatureR);
}

1;
