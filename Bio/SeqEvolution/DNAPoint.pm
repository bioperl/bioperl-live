#
# BioPerl module for Bio::SeqEvolution::DNAPoint
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki at bioperl dot org>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqEvolution::DNAPoint - evolve a sequence by point mutations

=head1 SYNOPSIS

  # $seq is a Bio::PrimarySeqI to mutate
  $evolve = Bio::SeqEvolution::Factory->new (-rate => 5,
                                             -seq => $seq,
                                             -identity => 50
                                             );
  $newseq = $evolve->next_seq;


=head1 DESCRIPTION

Bio::SeqEvolution::DNAPoint implements the simplest evolution model:
nucleotides change by point mutations, only. Transition/transversion
rate of the change, rate(), can be set.

The new sequences are named with the id of the reference sequence
added with a running number. Placing a new sequence into a factory to
be evolved resets that counter. It can also be called directly with
L<reset_sequence_counter>.

The default sequence type returned is Bio::PrimarySeq. This can be
changed to any Bio::PrimarySeqI compliant sequence class.

Internally the probability of the change of one nucleotide is mapped
to scale from 0 to 100. The probability of the transition occupies
range from 0 to some value. The remaining range is divided equally
among the two transversion nucleotides. A random number is then
generated to pick up one change.

Not that the default transition/transversion rate, 1:1, leads to
observed transition/transversion ratio of 1:2 simply because there is
only one transition nucleotide versus two transversion nucleotides.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

  Heikki Lehvaslaiho E<lt>heikki at bioperl dot orgE<gt>

=head1 CONTRIBUTORS

Additional contributor's names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqEvolution::DNAPoint;
use strict;
use Bio::Root::Root;
use Bio::SeqEvolution::Factory;

use Bio::Variation::DNAMutation;
use Bio::Variation::Allele;
use Bio::SimpleAlign;

use base qw(Bio::SeqEvolution::Factory);


sub _initialize {
    my($self, @args) = @_;

    $self->SUPER::_initialize(@args);
    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param; # lowercase keys

    exists $param{'-rate'} && $self->rate($param{'-rate'});

    $self->_init_mutation_engine;
}


sub _init_mutation_engine {
    my $self = shift;

    # arrays of possible changes have transitions as first items
    my %changes;
    $self->{'_changes'}->{'a'} = ['t', 'c', 'g'];
    $self->{'_changes'}->{'t'} = ['a', 'c', 'g'];
    $self->{'_changes'}->{'c'} = ['g', 'a', 't'];
    $self->{'_changes'}->{'g'} = ['c', 'a', 't'];


    # given the desired rate, find out where cut off points need to be
    # when random numbers are generated from 0 to 100
    # we are ignoring identical mutations (e.g. A->A) to speed things up
    my $bin_size = 100/($self->rate + 2);
    $self->{'_transition'} = 100 - (2*$bin_size);
    $self->{'_first_transversion'} = $self->{'_transition'} + $bin_size;

    $self->_init_alignment;
}

sub _init_alignment {
    my $self = shift;

    # put the initial sequence into the alignment object
    $self->{'_align'} = Bio::SimpleAlign->new(-verbose => -1);
    return unless $self->seq;
    $self->{'_ori_locatableseq'} = Bio::LocatableSeq->new(-id => 'ori',
                                                         -seq=> $self->seq->seq);
    $self->{'_mut_locatableseq'} = Bio::LocatableSeq->new(-id => 'mut',
                                                         -seq=> $self->seq->seq);
    $self->{'_align'}->add_seq($self->{'_ori_locatableseq'});
    $self->{'_align'}->add_seq($self->{'_mut_locatableseq'});
}

=head2 seq

 Title   : seq
 Usage   : $obj->seq($newval)
 Function: Set the sequence object for the original sequence
 Returns : The sequence object
 Args    : newvalue (optional)

Setting this will reset mutation and generated mutation counters.

=cut

sub seq{
   my $self = shift;

   if (@_) {
       my $seq = shift;
       $self->throw('Need a valid Bio::PrimarySeqI, not [', ref($seq), ']')
           unless $seq->isa('Bio::PrimarySeqI');
       
       $self->throw('Only nucleotide sequences are supported')
           if $seq->alphabet eq 'protein';
       $self->throw('No ambiquos nucleotides allowed in the input sequence')
           if $seq->seq =~ m/[^acgt]/;

       $self->{'_seq'} = $seq;

       # unify the look of sequence strings and cache the information
       $self->{'_ori_string'} = lc $seq->seq; # lower case
       $self->{'_ori_string'} =~ s/u/t/; # simplyfy our life; modules should deal with the change anyway
       $self->{'_seq_length'} = $seq->length;

       $self->reset_sequence_counter;
   }
   return $self->{'_seq'};
}

=head2 set_mutated_seq

  Title   : seq_mutated_seq
  Usage   : $obj->set_mutated_seq($newval)
  Function: In case of mutating a sequence with multiple evolvers, this
  Returns : set_mutated_seq
  Args    : newvalue (optional)

=cut

sub set_mutated_seq {
    my $self = shift;

    if (@_) {
        my $seq = shift;
        $self->throw('Need a valid Bio::PrimarySeqI, not [', ref($seq), ']')
            unless $seq->isa('Bio::PrimarySeqI');
        $self->throw('Only nucleotide sequences are supported')
            if $seq->alphabet eq 'protein';
        $self->throw('No ambiquos nucleotides allowed in the input sequence')
            if $seq->seq =~ m/[^acgt]/;

        $self->{'_seq_mutated'} = $seq;

        # unify the look of sequence strings and cache the information
        $self->{'_mut_string'} = lc $seq->seq; # lower case
        $self->{'_mut_string'} =~ s/u/t/; # simplyfy our life; modules should deal with the change anyway


        $self->reset_sequence_counter;
    }
    #set returned sequence to be the last mutated string
    $self->{'_seq'}->seq($self->{'_mut_string'});
    return $self->{'_seq'};
}


=head2 rate

  Title   : rate
  Usage   : $obj->rate($newval)
  Function: Set the transition/transversion rate.
  Returns : value of rate
  Args    : newvalue (optional)

Transition/transversion ratio is an observed attribute of an sequence
comparison. We are dealing here with the transition/transversion rate
that we set for our model of sequence evolution.

Note that we are using standard nucleotide alphabet and that there can
there is only one transition versus two possible transversions. Rate 2
is needed to have an observed transition/transversion ratio of 1.

=cut

sub rate{
   my $self = shift;
   if (@_) {
       $self->{'_rate'} = shift @_;
       $self->_init_mutation_engine;
   }
   return $self->{'_rate'} || 1;
}

=head2 next_seq

  Title   : next_seq
  Usage   : $obj->next_seq
  Function: Evolve the reference sequence to desired level
  Returns : A new sequence object mutated from the reference sequence
  Args    : -

=cut

sub next_seq {
    my $self = shift;
    $self->{'_mut_string'} = $self->{'_ori_string'};
    $self->reset_mutation_counter;

    $self->{'_mutations'} = [];

    while (1) {
        # find the location in the string to change
        my $loc = int (rand length($self->{'_mut_string'})) + 1;

        $self->mutate($loc); # for modularity

        # stop evolving if any of the limit has been reached
        last if $self->identity && $self->get_alignment_identity <= $self->identity;
        last if $self->pam &&  100*$self->get_mutation_counter/$self->{'_seq_length'} >= $self->pam;
        last if $self->mutation_count &&  $self->get_mutation_counter >= $self->mutation_count;
    }
    $self->_increase_sequence_counter;

    my $type = $self->seq_type;
    return $type->new(-id => $self->seq->id.  "-". $self->get_sequence_counter,
                      -description => $self->seq->description,
                      -seq => $self->{'_mut_string'}
                     )
}

=head2 mutate

  Title   : mutate
  Usage   : $obj->mutate
  Function: mutate the sequence at the given location according to the model
  Returns : true
  Args    : integer, start location of the mutation, required argument

Called from next_seq().

=cut

sub mutate {
    my $self = shift;
    my $loc = shift;
    $self->throw('the first argument is the location of the mutation') unless $loc;

    # nucleotide to change
    my $oldnuc = substr $self->{'_mut_string'}, $loc-1, 1;
    my $newnuc;


    # find the nucleotide it is changed to
    my $choose = rand(100);     # scale is 0-100
    if ($choose < $self->{'_transition'} ) {
        $newnuc =  $self->{'_changes'}->{$oldnuc}[0];
    } elsif ($choose < $self->{'_first_transversion'} ) {
        $newnuc =  $self->{'_changes'}->{$oldnuc}[1];
    } else {
        $newnuc =  $self->{'_changes'}->{$oldnuc}[2];
    }

    # do the change
    substr $self->{'_mut_string'}, $loc-1, 1 , $newnuc;
    $self->_increase_mutation_counter;
    $self->{'_mut_locatableseq'}->seq($self->{'_mut_string'});

    print STDERR "$loc$oldnuc>$newnuc\n" if $self->verbose > 0;

    push @{$self->{'_mutations'}}, "$loc$oldnuc>$newnuc";
}


1;
