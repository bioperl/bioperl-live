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

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org            - General discussion
  bioperl-guts-l@bioperl.org       - Technically-oriented discussion
  http://bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHORS - Ewan Birney & Lincoln Stein

Email: birney@ebi.ac.uk
       lstein@cshl.org


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::largefasta;
use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::Object

use Bio::SeqIO;
use Bio::Seq::LargePrimarySeq;

@ISA = qw(Bio::SeqIO);
# override new here to insure we instantiate this class 

sub new {
    my ($class,@args) = @_;    
    my $self = bless {}, $class;
    $self->_initialize(@args);
    return $self;
}

sub _initialize {
  my($self,@args) = @_;
  return unless $self->SUPER::_initialize(@args);
}

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

#  local $/ = "\n>";
  my $largeseq = new Bio::Seq::LargePrimarySeq;(-moltype=>'DNA');
  my ($id,$fulldesc,$entry);
  my $count = 0;
  while( defined ($entry = $self->_readline) ) {
      next if ( $entry eq '>');
      if( $entry =~ /\s*>(.+?)/ ) {
	  ($id,$fulldesc) = ($1 =~ /^\s*(\S+)\s*(.*)/)
	      or $self->warn("Can't parse fasta header");
	  $largeseq->display_id($id);
	  $largeseq->primary_id($id);	  
	  $largeseq->desc($fulldesc);
      } else {
	  $entry =~ s/\s//g;
	  $largeseq->add_sequence_as_string($entry);
      }
      (++$count % 1000 == 0 && $self->verbose() > 0) && print "line $count\n";
  }
  return $largeseq;
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
     my $str = $seq->seq;
     my $top = $seq->id();
     if ($seq->can('desc') and my $desc = $seq->desc()) {
	 $desc =~ s/\n//g;
        $top .= " $desc";
     }
     $str=~ s/(.{1,60})/$1\n/g;
     $self->_print (">",$top,"\n",$str) or return;
   }
   return 1;
}

1;
