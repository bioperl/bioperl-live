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

Bio::Seq::PrimedSeq - A representation of a sequence and two primers 
flanking a target region

=head1 SYNOPSIS

  # The easiest way to use this is probably either, (i), get the
  # output from Bio::Tools::Run::Primer3, Bio::Tools::Primer3, or 
  # Bio::Tools::PCRSimulation:

      # For example, start with a fasta file
      use Bio::SeqIO;
      use Bio::Tools::Run::Primer3;

      my $file = shift || die "need a file to read";
      my $seqin = Bio::SeqIO->new(-file => $file);
      my $seq = $seqin->next_seq;

      # use primer3 to design some primers
      my $primer3run = Bio::Tools::Run::Primer3->new(-seq => $seq);
      $primer3run -> run; # run it with the default parameters

      # create a file to write the results to
      my $seqout = Bio::SeqIO->new(-file => ">primed_sequence.gbk", 
                                   -format => 'genbank');

      # now just get all the results and write them out.
      while (my $results = $primer3run->next_primer) {
         $seqout->write_seq($results->annotated_seq);
      }

  # Or, (ii), to create a genbank file for a sequence and its cognate
  # primers:

     use Bio::SeqIO;
     use Bio::Seq::PrimedSeq;

     # have a sequence file ($file) with the template, and two primers
     # that match it, in fasta format

     my $file = shift || die "$0 <file>";
     my $seqin = Bio::SeqIO->new(-file => $file);

     # read three sequences
     my ($template, $leftprimer, $rightprimer) =
           ($seqin->next_seq, $seqin->next_seq, $seqin->next_seq);
     # set up the primed sequence object
     my $primedseq = Bio::Seq::PrimedSeq->new(-seq => $template, 
                                              -left_primer => $leftprimer,
                                              -right_primer => $rightprimer);
     # open a file for output
     my $seqout = Bio::SeqIO->new(-file => ">primed_sequence.gbk",
                                  -format => 'genbank');
     # print the sequence out
     $seqout->write_seq($primedseq->annotated_sequence);

  # This should output a genbank file with the two primers labeled.

=head1 DESCRIPTION

This module is a slightly glorified capsule containing a primed sequence. 
It was created to address the fact that a primer is more than a seqfeature 
and there need to be ways to represent the primer-sequence complex and 
the behaviors and attributes that are associated with the complex.

The primers are represented as Bio::SeqFeature::Primer objects, and should
be instantiated first.

A simple way to create a PrimedSeq object is as follows:

  my $primedseq = Bio::Seq::PrimedSeq->new(
          -seq          => $seq,  # Bio::Seq object,
          -left_primer  => $left, # Bio::SeqFeature::Primer object,
          -right_primer => $right # Bio::SeqFeature::Primer object,
  );

From the PrimedSeq object you should be able to retrieve
information about melting temperatures and what not on each of the primers 
and the amplicon.

This is based on the PrimedSeq.pm module started by Chad Matsalla, with 
additions/improvements by Rob Edwards.

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

use vars qw ($AUTOLOAD @RES %OK_FIELD $ID);

use base qw(Bio::Root::Root Bio::SeqFeature::Generic);

BEGIN {
 @RES = qw(); # nothing here yet, not sure what we want!
 foreach my $attr (@RES) {$OK_FIELD{$attr}++}
}

$ID = 'Bio::Tools::Analysis::Nucleotide::PrimedSeq';

sub AUTOLOAD {
 my $self = shift;
 my $attr = $AUTOLOAD;
 $attr =~ s/.*:://;
 $self->throw("Unallowed parameter: $attr !") unless $OK_FIELD{$attr};
 $self->{$attr} = shift if @_;
 return $self->{$attr};
}


=head2 new

 Title   : new()
 Usage   : $primed_sequence = Bio::SeqFeature::Primer->new( 
                                     -seq => $sequence,
                                     -left_primer => $left_primer,
                                     -right_primer => $right_primer);
 Function: A constructor for an object representing a primed sequence 
 Returns : A Bio::Seq::PrimedSeq object
 Args    :  -seq => a Bio::Seq object (required)
            -left_primer => a Bio::SeqFeature::Primer object (required)
            -right_primer => a Bio::SeqFeature::Primer object (required)

           Many other parameters can be included including all of the output
           parameters from the primer3 program. At the moment most of these
           parameters will not do anything.

=cut

sub new {

	# note, I have cleaned up a lot of the script that Chad had written here,
	# and I have removed the part where he removed the - before the tags.
	# Very confusing.

	my($class,%args) = @_;
	my $self = $class->SUPER::new(%args);
   # these are the absolute minimum components required to make
   # a primedseq

	foreach my $key (keys %args) {
		if ($key =~  /^-seq/i) {
			$self->{target_sequence} = $args{$key};
			next;
		} else {
			my $okey;
			($okey = $key) =~ s/^-//;
			if (($okey eq "left_primer" || $okey eq "right_primer") && 
				 ref($args{$key}) && $args{$key}->isa('Bio::SeqI') ) {
				# We have been given a Bio::Seq object, 
				# make it a Bio::SeqFeature::Primer object
				$self->{$okey} = Bio::SeqFeature::Primer->new(-seq => $args{$key});
				push @{$self->{'arguments'}},$okey;
				next;
			}

			$self->{$okey} = $args{$key};
			push @{$self->{'arguments'}},$okey;
		}
	}
	# and now the insurance - make sure that things are ok
	if (!$self->{target_sequence} || !$self->{left_primer} || 
		 !$self->{right_primer} ) {
		$self->throw("You must provide a -seq, -left_primer, and -right_primer to create this object.");
	}
  
	if (! ref($self->{target_sequence}) ||
		 ! $self->{target_sequence}->isa('Bio::SeqI') ) {
		$self->throw("The target_sequence must be a Bio::Seq to create this object.");
 }
	if (! ref($self->{left_primer}) ||
		 ! $self->{left_primer}->isa("Bio::SeqFeature::Primer") || 
		 ! ref($self->{right_primer}) ||
		 ! $self->{right_primer}->isa("Bio::SeqFeature::Primer")) {
		$self->throw("You must provide a left_primer and right_primer, both as Bio::SeqFeature::Primer to create this object.");
	}
 
	# now we have the sequences, lets find out where they are
	$self->_place_seqs();
	return $self;
}


=head2 get_primer

 Title   : get_primer();
 Usage   : $primer = $primedseq->get_primer(l, left, left_primer, 
           -left_primer) to return the left primer or 
	        $primer = $primedseq->get_primer(r, right, right_primer, 
           -right_primer) to return the right primer or
	        $primer = $primedseq->get_primer(b, both, both_primers, 
           -both_primers)
           to return the left primer, right primer array
 Function: A getter for the left primer in thie PrimedSeq object.
 Returns : A Bio::SeqFeature::Primer object
 Args    : Either of (l, left, left_primer, -left_primer) to get left 
           primer.
           Either of (r, right, right_primer, -right_primer) to get 
           right primer
	        Either of (b, both, both_primers, -both_primers) to get 
           both primers. 
           Note that this is plural. [default]

=cut

sub get_primer() {
 my ($self, $arg) = @_;
 if (! defined $arg ) {
  return ($self->{'left_primer'}, $self->{'right_primer'});
 } elsif( $arg =~ /^l/ || $arg =~ /^-l/) { 
  # what a cheat, I couldn't be bothered to write all those or statements!
  # Hah, now you can write leprechaun to get the left primer.
  return $self->{'left_primer'}
 }
 elsif ($arg =~ /^r/ || $arg =~ /^-r/) {return $self->{'right_primer'}}
 elsif ($arg =~ /^b/ || $arg =~ /^-b/) {return ($self->{'left_primer'}, $self->{'right_primer'})}
}



=head2 annotated_sequence

 Title   : annotated_sequence
 Usage   : $annotated_sequence_object = $primedseq->annotated_sequence()
 Function: Get an annotated sequence object containg the left and right 
           primers
 Returns : An annotated sequence object or 0 if not defined.
 Args    : 
 Note    : Use this method to return a sequence object that you can write
           out (e.g. in GenBank format). See the example above.

=cut

sub annotated_sequence {
  my $self = shift;
  if (exists $self->{annotated_sequence}) {return $self->{annotated_sequence}}
  else {return 0}
}

=head2 amplicon

 Title   : amplicon
 Usage   : my $amplicon = $primedseq->amplicon()
 Function: Retrieve the amplicon as a sequence object
 Returns : A seq object. To get the sequence use $amplicon->seq
 Args    : None
 Note    : 

=cut

sub amplicon {
 my ($self,@args) = @_;
 my $id = $self->{'-seq'}->{'id'};
 unless ($id) {$id=""}
 # this just prevents a warning when $self->{'-seq'}->{'id'} is not defined
 $id = "Amplicon from ".$id;
 
 my $seqobj=Bio::Seq->new(-id=>$id, seq=>$self->{'amplicon_sequence'});
 return $seqobj;
}


=head2 seq

 Title   : seq
 Usage   : my $seqobj = $primedseq->seq()
 Function: Retrieve the target sequence as a sequence object
 Returns : A seq object. To get the sequence use $seqobj->seq
 Args    : None
 Note    : 

=cut

sub seq {
 my $self = shift;
 return $self->{target_sequence};
}

=head2 _place_seqs

 Title   : _place_seqs
 Usage   : $self->_place_seqs()
 Function: An internal method to place the primers on the sequence and 
           set up the ranges of the sequences
 Returns : Nothing
 Args    : None
 Note    : Internal use only

=cut

sub _place_seqs {
	my $self = shift;
 
	# we are going to pull out the target sequence, and then the primer sequences
	my $target_sequence = $self->{'target_sequence'}->seq();
 
	# left primer
	my $left_seq = $self->{'left_primer'}->seq()->seq();

	my $rprc = $self->{'right_primer'}->seq()->revcom();
 
	my $right_seq=$rprc->seq();
  
	# now just change the case, because we keep getting screwed on this
	$target_sequence=uc($target_sequence);
	$left_seq=uc($left_seq);
	$right_seq=uc($right_seq);
 
	unless ($target_sequence =~ /(.*)$left_seq(.*)$right_seq(.*)/) {
		unless ($target_sequence =~ /$left_seq/) {$self->throw("Can't place left sequence on target!")}
		unless ($target_sequence =~ /$right_seq/) {$self->throw("Can't place right sequence on target!")}
 }
 
	my ($before, $middle, $after) = ($1, $2, $3); # note didn't use $`, $', and $& because they are bad. Just use length instead.

	# cool now we can figure out lengths and what not.
	# we'll figure out the position and compare it to known positions (e.g. from primer3)
 
	my $left_location = length($before). ",". length($left_seq);
	my $right_location = (length($target_sequence)-length($after)-1).",".length($right_seq);
	my $amplicon_size = length($left_seq)+length($middle)+length($right_seq);
 
	if (exists $self->{'left_primer'}->{'PRIMER_LEFT'}) {
		# this is the left primer from primer3 input
		# just check to make sure it is right
		unless ($self->{'left_primer'}->{'PRIMER_LEFT'} eq $left_location) {
			$self->warn("Note got |".$self->{'left_primer'}->{'PRIMER_LEFT'}."| from primer3 and |$left_location| for the left primer. You should email redwards\@utmem.edu about this.");
		}
	}
	else {
		$self->{'left_primer'}->{'PRIMER_LEFT'}=$left_location;
	}
 
	if (exists $self->{'right_primer'}->{'PRIMER_RIGHT'}) {
		# this is the right primer from primer3 input
		# just check to make sure it is right
		unless ($self->{'right_primer'}->{'PRIMER_RIGHT'} eq $right_location) {
			$self->warn("Note got |".$self->{'right_primer'}->{'PRIMER_RIGHT'}."| from primer3 and |$right_location| for the right primer. You should email redwards\@utmem.edu about this.");
		}
	}
	else {
		$self->{'right_primer'}->{'PRIMER_RIGHT'}=$right_location;
	}
 
	if (exists $self->{'PRIMER_PRODUCT_SIZE'}) {
		# this is the product size from primer3 input
		# just check to make sure it is right
		unless ($self->{'PRIMER_PRODUCT_SIZE'} eq $amplicon_size) {
			$self->warn("Note got |".$self->{'PRIMER_PRODUCT_SIZE'}."| from primer3 and |$amplicon_size| for the size. You should email redwards\@utmem.edu about this.");
		}
	}
	else {
		$self->{'PRIMER_PRODUCT_SIZE'} = $amplicon_size;
	}
 
	$self->{'amplicon_sequence'}= lc($left_seq).uc($middle).lc($right_seq); # I put this in a different case, but I think the seqobj may revert this
	
	$self->_set_seqfeature;
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
	unless ($self->{'left_primer'}->{'PRIMER_LEFT'} && 
			  $self->{'right_primer'}->{'PRIMER_RIGHT'}) {
		$self->warn("hmmm. Haven't placed primers, but trying to make annotated sequence");
		return 0;
	}
	my ($start, $length) = split /,/, $self->{'left_primer'}->{'PRIMER_LEFT'};
	my $tm=$self->{'left_primer'}->{'PRIMER_LEFT_TM'} || $self->{'left_primer'}->Tm || 0;

	my $seqfeatureL=Bio::SeqFeature::Generic->new(
						  -start => $start+1, -end => $start+$length, -strand => 1,
                    -primary_tag => 'left_primer', -source => 'primer3',
                    -tag    => {new => 1, author => 'Bio::Seq::PrimedSeq', Tm => $tm}
															  );
 
	($start, $length) = split /,/, $self->{'right_primer'}->{'PRIMER_RIGHT'};
	$tm=$self->{'right_primer'}->{'PRIMER_RIGHT_TM'} || $self->{'right_primer'}->Tm || 0;
 
	my $seqfeatureR=Bio::SeqFeature::Generic->new(
   -start => $start-$length+2, -end => $start+1, -strand => -1,
   -primary_tag => 'right_primer', -source => 'primer3',
   -tag    => {new => 1, author => 'Bio::Seq::PrimedSeq', Tm => $tm}
															  );

	# now add the sequences to a annotated sequence
	$self->{annotated_sequence} = $self->{target_sequence};
	$self->{annotated_sequence}->add_SeqFeature($seqfeatureL);
	$self->{annotated_sequence}->add_SeqFeature($seqfeatureR);
}

1;
