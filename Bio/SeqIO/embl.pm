#
# BioPerl module for Bio::SeqIO::fasta
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#          and Lincoln Stein <lstein@cshl.org>
#
# Copyright Ewan Birney & Lincoln Stein
#
# You may distribute this module under the same terms as perl itself
# _history
# October 18, 1999  Largely rewritten by Lincoln Stein

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::embl - EMBL sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::SeqIO class.

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from EMBL flat
file databases.

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

=head1 AUTHORS - Ewan Birney & Lincoln Stein

Email: birney@sanger.ac.uk
       lstein@cshl.org


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::embl;
use vars '@ISA';
use strict;
use Bio::Seq;
# Object preamble - inheriets from Bio::Root::Object

use Bio::SeqIO;

@ISA = 'Bio::SeqIO';

# new() is inherited from Bio::SeqIO

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    :


=cut

sub next_seq{
   my ($self,@args) = @_;
   my ($seq,$fh,$c,$line,$name,$desc,$seqc);
   
   $line = $self->_readline;
   $line =~ /^ID\s+(\S+)/ || $self->throw("EMBL stream with no ID. Not embl in my book");
   $name = $1;

   while( defined ($_ = $self->_readline) ) {
       /^DE\s+(\S.*\S)/ && do { $desc = $1;};
       /^SQ/ && last;
   }

   while( defined ($_ = $self->_readline) ) {
       /^\/\// && last;
       $_ = uc($_);
       s/\W//g;
       $seqc .= $_;
       eof $fh && last;
   }

   $seq = Bio::Seq->new(-seq => $seqc , -id => $name, -desc => $desc);
   
   return $seq;

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
   for my $seq (@seq) {
     my $i;
     my $str = $seq->seq;
     $self->_print("ID   ", $seq->id(), 
		   "\nDE    ", $seq->desc(), 
		   "\nCC   \nCC   Written by Bioperl SeqIO module.\nCC   Only the information in the Sequence object is written in this file\nCC   \nSQ   \n") or return;
     while ($str =~ /(.{1,50})/g) {
       my $line = $1;
       $line =~ s/(.{10})/$1 /g;
       $self->_print('    ',$line,"\n") or return;
     }
     $self->_print("//\n") or return;
   }
   return 1;
}
    
1;
