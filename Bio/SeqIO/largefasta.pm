# $Id$
# BioPerl module for Bio::SeqIO::largefasta
#
# Cared for by Jason Stajich
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself
# _history
# 
# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::largefasta - method i/o on very large fasta sequence files

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::SeqIO class.

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from fasta flat
file databases.

This module handles very large sequence files by using the
L<Bio::Seq::LargePrimarySeq> module to store all the sequence data in
a file.  This can be a problem if you have limited disk space on your
computer because this will effectively cause 2 copies of the sequence
file to reside on disk for the life of the
L<Bio::Seq::LargePrimarySeq> object.  The default location for this is
specified by the L<File::Spec>-E<gt>tmpdir routine which is usually /tmp
on UNIX.  If a sequence file is larger than the swap space (capacity
of the /tmp dir) this could cause problems for the machine.  It is
possible to set the directory where the temporary file is located by
adding the following line to your code BEFORE calling next_seq.

    $Bio::Seq::LargePrimarySeq::DEFAULT_TEMP_DIR = 'newdir';

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org            - General discussion
  http://bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHORS - Jason Stajich

Email: jason@chg.mc.duke.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::largefasta;
use vars qw(@ISA $FASTALINELEN);
use strict;
# Object preamble - inherits from Bio::Root::Object

use Bio::SeqIO;
use Bio::Seq::LargePrimarySeq;
use Bio::Seq::LargeSeq;

$FASTALINELEN = 60;
@ISA = qw(Bio::SeqIO);

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    : NONE

=cut

sub next_seq {
    return next_primary_seq( $_[0], 1 );
}

=head2 next_primary_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::PrimarySeq object
 Args    : NONE

=cut

sub next_primary_seq {
  my( $self, $as_next_seq ) = @_;

#  local $/ = "\n";
  my $largepseq = new Bio::Seq::LargePrimarySeq();
  my ($id,$fulldesc,$entry);
  my $count = 0;
  my $seen = 0;
  while( defined ($entry = $self->_readline) ) {      
      if( $seen == 1 && $entry =~ /^\s*>/ ) {
	  $self->_pushback($entry);
	  return $largepseq;
      }
      if ( $entry eq '>' ) { $seen = 1; next; }      
      elsif( $entry =~ /\s*>(.+?)$/ ) {
	  $seen = 1;
	  ($id,$fulldesc) = ($1 =~ /^\s*(\S+)\s*(.*)$/)
	      or $self->warn("Can't parse fasta header");
	  $largepseq->display_id($id);
	  $largepseq->primary_id($id);	  
	  $largepseq->desc($fulldesc);
      } else {
	  $entry =~ s/\s+//g;
	  $largepseq->add_sequence_as_string($entry);
      }
      (++$count % 1000 == 0 && $self->verbose() > 0) && print "line $count\n";
  }
  if( ! $seen ) { return undef; } 
  if( $as_next_seq ) {
      return new Bio::Seq::LargeSeq(-primaryseq => $largepseq );      
  } else {
      return $largepseq;
  }
}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq(@seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq object


=cut

sub write_seq {
   my ($self,@seq) = @_;
   foreach my $seq (@seq) {       
     my $top = $seq->id();
     if ($seq->can('desc') and my $desc = $seq->desc()) {
	 $desc =~ s/\n//g;
	 $top .= " $desc";
     }
     $self->_print (">",$top,"\n");
     my $end = $seq->length();
     my $start = 1;
     while( $start < $end ) {
	 my $stop = $start + $FASTALINELEN - 1;
	 $stop = $end if( $stop > $end );
	 $self->_print($seq->subseq($start,$stop), "\n");
	 $start += $FASTALINELEN;
     }
   }
   return 1;
}

1;
