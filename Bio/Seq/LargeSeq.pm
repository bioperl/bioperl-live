
#
# BioPerl module for Bio::LargeSeq
#
# Cared for by Jason Stajich
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::LargeSeq - SeqI compliant object that stores sequence as files in /tmp 

=head1 SYNOPSIS

  # normal primary seq usage

=head1 DESCRIPTION

This object stores a sequence as a series of files in a temporary
directory. The aim is to allow someone the ability to store very large
sequences (eg, > 100MBases) in a file system without running out of memory
(eg, on a 64 MB real memory machine!). 

Of course, to actually make use of this functionality, the programs
which use this object B<must> not call $primary_seq->seq otherwise the
entire sequence will come out into memory and probably paste your
machine. However, calls $primary_seq->subseq(10,100) will cause only
90 characters to be brought into real memory.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::LargeSeq;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inherits from Bio::Root::Objecttest 8, 

use Bio::Root::RootI;
use Bio::Seq::LargePrimarySeq;
use Bio::Seq;

@ISA = qw(Bio::Seq);


sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my ($pseq) = $self->_rearrange([qw(PRIMARYSEQ)], @args);

    if( ! defined $pseq ) {
	$pseq = new Bio::Seq::LargePrimarySeq(@args);
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
    return new Bio::Seq::LargeSeq(
				  '-display_id' => $self->display_id,
				  '-accession_number' => $self->accession_number,
				  '-desc' => $self->desc,
				  '-moltype' => $self->moltype,
				  -primaryseq => 
				  $self->primary_seq->trunc($s,$e));

}

1;
