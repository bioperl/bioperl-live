## $Id$

# BioPerl module for Bio::Tools::HMM
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Yee Man Chan <ymc@yahoo.com>
#
# Copyright Yee Man Chan
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::HMM - Perl extension to perform Hidden Markov Model calculations

=head1 SYNOPSIS

  use Bio::Tools::HMM;
  use Bio::SeqIO;
  use Bio::Matrix::Scoring;

  # create a HMM object
  # ACGT are the bases NC mean non-coding and coding
  $hmm = Bio::Tools::HMM->new('-symbols' => "ACGT", '-states' => "NC");

  # initialize some training observation sequences
  $seq1 = Bio::SeqIO->new(-file => $ARGV[0], -format => 'fasta');
  $seq2 = Bio::SeqIO->new(-file => $ARGV[1], -format => 'fasta');
  @seqs = ($seq1, $seq2);

  # train the HMM with the observation sequences
  $hmm->baum_welch_training(\@seqs);

  # get parameters
  $init = $hmm->init_prob; # returns an array reference
  $matrix1 = $hmm->transition_prob; # returns Bio::Matrix::Scoring
  $matrix2 = $hmm->emission_prob; # returns Bio::Matrix::Scoring

  # initialize training hidden state sequences
  $hs1 = "NCNCNNNNNCCCCNNCCCNNNNC";
  $hs2 = "NCNNCNNNNNNCNCNCNNNCNCN";
  @hss = ($hs1, $hs2);

  # train the HMM with both observation sequences and hidden state
  # sequences
  $hmm->statistical_training(\@seqs, \@hss);

  # with the newly calibrated HMM, we can use viterbi algorithm
  # to obtain the hidden state sequence underlying an observation 
  # sequence
  $hss = $hmm->viterbi($seq); # returns a string of hidden states

=head1 DESCRIPTION

Hidden Markov Model (HMM) was first introduced by Baum and his colleagues
in a series of classic papers in the late 1960s and 1970s. It was first
applied to the field of speech recognition with great success in the 1970s.

Explosion in the amount sequencing data in the 1990s opened the field
of Biological Sequence Analysis. Seeing HMM's effectiveness in detecing
signals in biological sequences, Krogh, Mian and Haussler used HMM to find
genes in E. coli DNA in a classical paper in 1994. Since then, there have
been extensive application of HMM to other area of Biology, for example,
multiple sequence alignment, CpG island detection and so on.

=head1 DEPENDENCIES

This package comes with the main bioperl distribution. You also need
to install the lastest bioperl-ext package which contains the XS code
that implements the algorithms. This package won't work if you haven't
compiled the bioperl-ext package.

=head1 TO-DO


=over 3

=item 1.

Allow people to set and get the tolerance level in the EM algorithm.

=item 2.

Allow people to set and get the maximum number of iterations 
to run in the EM algorithm.

=item 3.

A function to calculate the probability of an observation sequence

=item 4.

A function to do posterior decoding, ie to find the probabilty of
seeing a certain observation symbol at position i.

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 
 
Please direct usage questions or support issues to the mailing list:
  
L<bioperl-l@bioperl.org>
  
rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR

        This implementation was written by Yee Man Chan (ymc@yahoo.com).
        Copyright (c) 2005 Yee Man Chan. All rights reserved. This program
        is free software; you can redistribute it and/or modify it under
        the same terms as Perl itself. All the code are written by Yee
        Man Chan without borrowing any code from anywhere.

=cut

package Bio::Tools::HMM;

use strict;
use warnings;

use base qw(Bio::Root::Root);

BEGIN {
    eval {
        require Bio::Ext::HMM;
    };
    if ( $@ ) {
        die("\nThe C-compiled engine for Hidden Markov Model (HMM) has not been installed.\n Please read the install the bioperl-ext package\n\n");
        exit(1);
    }
}

sub new {
   my ($class, @args) = @_;

   my $self = $class->SUPER::new(@args);

   my ($symbols, $states, $init, $a_mat, $e_mat) = $self->_rearrange([qw(SYMBOLS
								STATES
                                                                INIT
                                                                AMAT
                                                                EMAT
							)], @args);

   $self->throw("Observation Symbols are not defined!") unless defined $symbols; 
   $self->throw("Hidden States are not defined!") unless defined $states; 

   if (defined $symbols) {
      if (scalar($symbols)) {
         # check duplicate symbols
         if ($self->symbols($symbols) < 0) {
            $self->throw("Duplicate symbols!\n");
         }
      }
      else {
         $self->throw("We don't support list of symbols in this version.\n");
      }
   }

   if (defined $states) {
      if (scalar($states)) {
         # check duplicate states
         if ($self->states($states) < 0) {
            $self->throw("Duplicate states!\n");
         }
      }
      else {
         $self->throw("We don't support list of states in this version.\n");
      }
   }

   $self->{'hmm'} = Bio::Ext::HMM::HMM->new($symbols, $states);
   return $self;
}

=head2 likelihood

 Title   : likelihood
 Usage   : $prob = $hmm->likelihood($seq)
 Function: Calculate the probability of an observation sequence given an HMM
 Returns : An floating point number that is the logarithm of the probability
           of an observation sequence given an HMM
 Args    : The only argument is a string that is the observation sequence
           you are interested in. Note that the sequence must not contain
           any character that is not in the alphabet of observation symbols.

=cut

sub likelihood {
   my ($self, $seq) = @_;
   my $valid_symbols;

   if( ! defined $seq) {
      $self->warn("Cannot calculate without supply an observation sequence!");
      return;
   }
   my $s = $self->{'symbols'};
   $_ = $seq;
   $valid_symbols = eval "tr/$s//;"; 
   if ($valid_symbols != length($seq)) {
      $self->throw("Observation Sequence contains characters that is not in the
                    alphabet of observation symbols!\n");
   }
   return Bio::Ext::HMM->HMM_likelihood($self->{'hmm'}, $seq);
}

=head2 statistical_training

 Title   : statistical_training
 Usage   : $hmm->statistical_training(\@seqs, \@hss)
 Function: Estimate the parameters of an HMM given an array of observation 
           sequence and an array of the corresponding hidden state 
           sequences
 Returns : Returns nothing. The parameters of the HMM will be set to the 
           estimated values
 Args    : The first argument is a reference to an array of observation 
           sequences. The second argument is a reference to an array of
           hidden state sequences. Note that the lengths of an observation
           sequence and a hidden state sequence must be the same.

=cut

sub statistical_training {
   my ($self, $seqs, $hss) = @_;
   my $valid_symbols;
   my ($seq_cnt, $hs_cnt);
   my $i;

   if( ! defined $seqs or ! defined $hss) {
      $self->warn("Cannot calculate without supply an observation and a hidden state sequence!");
      return;
   }
   $seq_cnt = @{$seqs};
   $hs_cnt = @{$seqs};
   if ($seq_cnt != $hs_cnt) {
      $self->throw("There must be the same number of observation sequences and 
                    hidden state sequences!\n");
   }
   for ($i = 0; $i < $seq_cnt; ++$i) {
      if (length(@{$seqs}[$i]) != length(@{$hss}[$i])) {
         $self->throw("The corresponding observation sequences and hidden state sequences must be of the same length!\n");
      }
   }
   foreach my $seq (@{$seqs}) {
      my $s = $self->{'symbols'};
      $_ = $seq;
      $valid_symbols = eval "tr/$s//;"; 
      if ($valid_symbols != length($seq)) {
         $self->throw("Observation Sequence contains characters that is not in the
                alphabet of observation symbols!\n");
      }
   }
   foreach my $seq (@{$hss}) {
      my $s = $self->{'states'};
      $_ = $seq;
      $valid_symbols = eval "tr/$s//;"; 
      if ($valid_symbols != length($seq)) {
         $self->throw("Hidden State Sequence contains characters that is not in the
                alphabet of hidden states!\n");
      }
   }
   Bio::Ext::HMM->HMM_statistical_training($self->{'hmm'}, $seqs, $hss);
}

=head2 baum_welch_training

 Title   : baum_welch_training
 Usage   : $hmm->baum_welch_training(\@seqs)
 Function: Estimate the parameters of an HMM given an array of observation 
           sequence
 Returns : Returns nothing. The parameters of the HMM will be set to the 
           estimated values
 Args    : The only argument is a reference to an array of observation 
           sequences. 

=cut

sub baum_welch_training {
   my ($self, $seqs) = @_;
   my $valid_symbols;

   if( ! defined $seqs) {
      $self->warn("Cannot calculate without supply an observation sequence!");
      return;
   }
   foreach my $seq (@{$seqs}) {
      my $s = $self->{'symbols'};
      $_ = $seq;
      $valid_symbols = eval "tr/$s//;"; 
      if ($valid_symbols != length($seq)) {
         $self->throw("Observation Sequence contains characters that is not in the
                alphabet of observation symbols!\n");
      }
   }
   Bio::Ext::HMM->HMM_baum_welch_training($self->{'hmm'}, $seqs);
}

=head2 viterbi

 Title   : viterbi
 Usage   : $hss = $hmm->viterbi($seq)
 Function: Find out the hidden state sequence that can maximize the 
           probability of seeing observation sequence $seq.
 Returns : Returns a string that is the hidden state sequence that maximizes
           the probability of seeing $seq.
 Args    : The only argument is an observation sequence.

=cut

sub viterbi {
   my ($self, $seq) = @_;
   my $valid_symbols;

   if( ! defined $seq) {
      $self->warn("Cannot calculate without supply an observation sequence!");
      return;
   }
   my $s = $self->{'symbols'};
   $_ = $seq;
   $valid_symbols = eval "tr/$s//;"; 
   if ($valid_symbols != length($seq)) {
      $self->throw("Observation Sequence contains characters that is not in the
             alphabet of observation symbols!\n");
   }
   return Bio::Ext::HMM->HMM_viterbi($self->{'hmm'}, $seq);
}

=head2 symbols

 Title     : symbols 
 Usage     : $symbols = $hmm->symbols() #get
           : $hmm->symbols($value) #set
 Function  : the set get for the observation symbols
 Example   :
 Returns   : symbols string
 Arguments : new value

=cut

sub symbols {
   my ($self,$val) = @_;
   my %alphabets = ();
   my $c;

   if ( defined $val ) {
# find duplicate
      
      for (my $i = 0; $i < length($val); ++$i) {
         $c = substr($val, $i, 1);
         if (defined $alphabets{$c}) {
            $self->throw("Can't have duplicate symbols!");
         }
         else {
            $alphabets{$c} = 1;
         }
      }
      $self->{'symbols'} = $val;
   }
   return $self->{'symbols'};
}


=head2 states

 Title     : states
 Usage     : $states = $hmm->states() #get
           : $hmm->states($value) #set
 Function  : the set get for the hidden states
 Example   :
 Returns   : states string
 Arguments : new value

=cut

sub states {
   my ($self,$val) = @_;
   my %alphabets = ();
   my $c;

   if ( defined $val ) {
# find duplicate
      
      for (my $i = 0; $i < length($val); ++$i) {
         $c = substr($val, $i, 1);
         if (defined $alphabets{$c}) {
            $self->throw("Can't have duplicate states!");
         }
         else {
            $alphabets{$c} = 1;
         }
      }
      $self->{'states'} = $val;
   }
   return $self->{'states'};
}

=head2 init_prob

 Title     : init_prob
 Usage     : $init = $hmm->init_prob() #get
           : $hmm->transition_prob(\@init) #set
 Function  : the set get for the initial probability array
 Example   :
 Returns   : reference to double array
 Arguments : new value

=cut

sub init_prob {
   my ($self, $init) = @_;
   my $i;
   my @A;

   if (defined $init) {
      if (ref($init)) {
         my $size = @{$init};
         my $sum = 0.0;
         foreach (@{$init}) {
            $sum += $_;
         }
         if ($sum != 1.0) {
            $self->throw("The sum of initial probability array must be 1.0!\n");
         }
         if ($size != length($self->{'states'})) {
            $self->throw("The size of init array $size is different from the number of HMM's hidden states!\n");
         }
         for ($i = 0; $i < $size; ++$i) {
            Bio::Ext::HMM::HMM->set_init_entry($self->{'hmm'}, substr($self->{'states'}, $i, 1), @{$init}[$i]);
         }
      }
      else {
         $self->throw("Initial Probability array must be a reference!\n");
      }
   }
   else {
      for ($i = 0; $i < length($self->{'states'}); ++$i) {
         $A[$i] = Bio::Ext::HMM::HMM->get_init_entry($self->{'hmm'}, substr($self->{'states'}, $i, 1));
      }
      return \@A;
   } 
}

=head2 transition_prob

 Title     : transition_prob
 Usage     : $transition_matrix = $hmm->transition_prob() #get
           : $hmm->transition_prob($matrix) #set
 Function  : the set get for the transition probability mairix
 Example   :
 Returns   : Bio::Matrix::Scoring 
 Arguments : new value

=cut

sub transition_prob {
   my ($self, $matrix) = @_;
   my ($i, $j);
   my @A;

   if (defined $matrix) {
      if ($matrix->isa('Bio::Matrix::Scoring')) {
         my $row = join("", $matrix->row_names);
         my $col = join("", $matrix->column_names);
         if ($row ne $self->{'states'}) {
            $self->throw("Names of the rows ($row) is different from the states of HMM " . $self->{'states'});
         } 
         if ($col ne $self->{'states'}) {
            $self->throw("Names of the columns ($col) is different from the states of HMM " . $self->{'states'});
         }
         for ($i = 0; $i < length($self->{'states'}); ++$i) {
            my $sum = 0.0;
            my $a = substr($self->{'states'}, $i, 1);
            for ($j = 0; $j < length($self->{'states'}); ++$j) {
               my $b = substr($self->{'states'}, $j, 1);
               $sum += $matrix->entry($a, $b);
            }
            if ($sum != 1.0) {
               $self->throw("Sum of probabilities for each from-state must be 1.0; got $sum\n");
            }
         }
         for ($i = 0; $i < length($self->{'states'}); ++$i) {
            my $a = substr($self->{'states'}, $i, 1);
            for ($j = 0; $j < length($self->{'states'}); ++$j) {
               my $b = substr($self->{'states'}, $j, 1);
               Bio::Ext::HMM::HMM->set_a_entry($self->{'hmm'}, $a, $b, $matrix->entry($a, $b));
            }
         }
      }
      else {
         $self->throw("Transition Probability matrix must be of type Bio::Matrix::Scoring.\n");
      }
   }
   else {
      for ($i = 0; $i < length($self->{'states'}); ++$i) {
         for ($j = 0; $j < length($self->{'states'}); ++$j) {
            $A[$i][$j] = Bio::Ext::HMM::HMM->get_a_entry($self->{'hmm'}, substr($self->{'states'}, $i, 1), substr($self->{'states'}, $j, 1));
         }
      }
      my @rows = split(//, $self->{'states'});
      return $matrix = Bio::Matrix::Scoring->new(-values => \@A, -rownames => \@rows, -colnames => \@rows);
   } 
}

=head2 emission_prob

 Title     : emission_prob
 Usage     : $emission_matrix = $hmm->emission_prob() #get
           : $hmm->emission_prob($matrix) #set
 Function  : the set get for the emission probability mairix
 Example   :
 Returns   : Bio::Matrix::Scoring 
 Arguments : new value

=cut

sub emission_prob {
   my ($self, $matrix) = @_;
   my ($i, $j);
   my @A;

   if (defined $matrix) {
      if ($matrix->isa('Bio::Matrix::Scoring')) {
         my $row = join("", $matrix->row_names);
         my $col = join("", $matrix->column_names);
         if ($row ne $self->{'states'}) {
            $self->throw("Names of the rows ($row) is different from the states of HMM " . $self->{'states'});
         } 
         if ($col ne $self->{'symbols'}) {
            $self->throw("Names of the columns ($col) is different from the symbols of HMM " . $self->{'symbols'});
         }
         for ($i = 0; $i < length($self->{'states'}); ++$i) {
            my $sum = 0.0;
            my $a = substr($self->{'states'}, $i, 1);
            for ($j = 0; $j < length($self->{'symbols'}); ++$j) {
               my $b = substr($self->{'symbols'}, $j, 1);
               $sum += $matrix->entry($a, $b);
            }
            if ($sum != 1.0) {
               $self->throw("Sum of probabilities for each state must be 1.0; got $sum\n");
            }
         }
         for ($i = 0; $i < length($self->{'states'}); ++$i) {
            my $a = substr($self->{'states'}, $i, 1);
            for ($j = 0; $j < length($self->{'symbols'}); ++$j) {
               my $b = substr($self->{'symbols'}, $j, 1);
               Bio::Ext::HMM::HMM->set_e_entry($self->{'hmm'}, $a, $b, $matrix->entry($a, $b));
            }
         }
      }
      else {
         $self->throw("Emission Probability matrix must be of type Bio::Matrix::Scoring.\n");
      }
   }
   else {
      for ($i = 0; $i < length($self->{'states'}); ++$i) {
         for ($j = 0; $j < length($self->{'symbols'}); ++$j) {
            $A[$i][$j] = Bio::Ext::HMM::HMM->get_e_entry($self->{'hmm'}, substr($self->{'states'}, $i, 1), substr($self->{'symbols'}, $j, 1));
         }
      }
      my @rows = split(//, $self->{'states'});
      my @cols = split(//, $self->{'symbols'});
      return $matrix = Bio::Matrix::Scoring->new(-values => \@A, -rownames => \@rows, -colnames => \@cols);
   } 
}

1;
