#
# BioPerl module for Bio::SeqIO::PIR
#
# Cared for by Aaron Mackey <amackey@virginia.edu>
#
# Copyright Aaron Mackey
#
# You may distribute this module under the same terms as perl itself
#
# _history
# October 18, 1999  Largely rewritten by Lincoln Stein

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::pir - PIR sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::SeqIO class.

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from fasta flat
file databases.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

   bioperl-l@bioperl.org             - General discussion
   bioperl-guts-l@bioperl.org        - Automated bug and CVS messages
   http://bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHORS

Aaron Mackey <amackey@virginia.edu>
Lincoln Stein <lstein@cshl.org>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::pir;
use vars qw(@ISA);
use strict;
use Bio::SeqIO;

@ISA = qw(Bio::SeqIO);

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    :

=cut

sub next_seq{
   my ($self,@args) = @_;
   my ($seq,$line,$name,$sfs,$desc);

   return unless $line = $self->_readline;
   $self->throw("PIR stream read attempted without leading '>P1;' [ $line ]")
     unless  $line =~ /^>(?:P|F)1;(\S+)\s*(\|.*)?\s*$/;
   $name = $1;
   $sfs = $2;

   chomp($desc = $self->_readline);
   local $/ = "";
   #my $junk = $self->_readline;  # throw away everything to first empty line
   $seq = $self->_readline;   # everything else is the sequence
   $seq =~ s/\s+//g;
   $seq = Bio::Seq->new(-seq => $seq,
			-id => $name,
			-desc => $desc,
			-names => defined $sfs ? { 'sfnum' => [ split(/\s*\|?\s+/, $sfs) ] } : undef );
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
   my ($self, @seq) = @_;
   for my $seq (@seq) {
     my $str = $seq->seq();
     $str =~ s/(.{10})/$1 /g;
     $str =~ s/(.{66})/$1\n/g;
     return unless $self->_print(">P1;", $seq->id(), 
				 "\n", $seq->desc(), "\n", 
				 "\n",$str, "\n");
   }
   return 1;
}

1;
