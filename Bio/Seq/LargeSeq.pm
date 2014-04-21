#
# BioPerl module for Bio::Seq::LargeSeq
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney, Jason Stajich
#
# Copyright Ewan Birney, Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::LargeSeq - SeqI compliant object that stores sequence as
files in /tmp

=head1 SYNOPSIS

  # normal primary seq usage

=head1 DESCRIPTION

This object stores a sequence as a series of files in a temporary
directory. The aim is to allow someone the ability to store very large
sequences (eg, E<gt> 100MBases) in a file system without running out
of memory (eg, on a 64 MB real memory machine!).

Of course, to actually make use of this functionality, the programs
which use this object B<must> not call $primary_seq-E<gt>seq otherwise
the entire sequence will come out into memory and probably paste your
machine. However, calls $primary_seq-E<gt>subseq(10,100) will cause
only 90 characters to be brought into real memory.

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

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::LargeSeq;
use vars qw($AUTOLOAD);
use strict;

# Object preamble

use Bio::Seq::LargePrimarySeq;

use base qw(Bio::Seq Bio::Seq::LargeSeqI);


sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($pseq) = $self->_rearrange([qw(PRIMARYSEQ)], @args);

    if( ! defined $pseq ) {
	$pseq = Bio::Seq::LargePrimarySeq->new(@args);
    }
    $self->primary_seq($pseq);

    return $self;
}


=head2 trunc

 Title   : trunc
 Usage   : $subseq = $myseq->trunc(10,100);
 Function: Provides a truncation of a sequence,

 Example :
 Returns : a fresh Bio::SeqI object
 Args    :

=cut

sub trunc {
    my ($self, $s, $e) = @_;
    return new Bio::Seq::LargeSeq
        ('-display_id' => $self->display_id,
         '-accession_number' => $self->accession_number,
         '-desc' => $self->desc,
         '-alphabet' => $self->alphabet,
         -primaryseq => $self->primary_seq->trunc($s,$e));
}

=head2 Bio::Seq::LargePrimarySeq methods

=cut

=head2 add_sequence_as_string

 Title   : add_sequence_as_string
 Usage   : $seq->add_sequence_as_string("CATGAT");
 Function: Appends additional residues to an existing LargePrimarySeq object.
           This allows one to build up a large sequence without storing
           entire object in memory.
 Returns : Current length of sequence
 Args    : string to append

=cut

sub add_sequence_as_string {
    my ($self,$str) = @_;
    return $self->primary_seq->add_sequence_as_string($str);
}

1;
