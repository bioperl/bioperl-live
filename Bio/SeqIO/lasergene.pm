#-----------------------------------------------------------------------------
# PACKAGE : Bio::SeqIO::Lasergene
# AUTHOR  : Malcolm Cook <mec@stowers-institute.org>
# CREATED : Feb 16 1999
# REVISION: $Id$
#            
# _History_
#
# This code is based on the Bio::SeqIO::raw module with
# the necessary minor tweaks necessary to get it to read (only) Lasergene formatted sequences 

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::Lasergene - Lasergene sequence file input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::SeqIO> class.

=head1 DESCRIPTION

This object can product Bio::Seq::RichSeq objects from Lasergene sequence files.

IT DOES NOT PARSE ANY ATTIBUTE VALUE PAIRS IN THE HEADER OF THE LASERGENE FORMATTED FILE.

IT DOES NOT WRITE THESE FILES EITHER.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://www.bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHORS

  Malcolm Cook  <mec@stowers-institute.org>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::SeqIO::lasergene;
use strict;
use vars qw(@ISA);

use Bio::SeqIO;
use Bio::Seq;

@ISA = qw(Bio::SeqIO);

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    :


=cut
use Bio::Annotation;
use Bio::Annotation::Comment;

sub next_seq{
   my ($self,@args) = @_;
   return undef unless defined($_ = $self->_readline()); # stream is empty! there is no next seq.
   my $seq =  Bio::Seq::RichSeq->new();
   my $comment = "";
   while (! /\^\^/) {
     chomp;
     $comment .= $_ . ";  ";

     #could consider populating BioSeq fields here, but I'm not in that business now!!!
     #my ($att,@val) = split  /\:\w*/  ; #  m/(.*)\:(.*)/;     
     #print STDERR "att:$att;val:@val\n";

     $_ = $self->_readline() or die "unexpected end of file";
   }
   $_ = $self->_readline() or die "unexpected end of file"; # this should be all and only sequence.
   $seq->seq($_);

   my ($annotation,%params,@features) = (new Bio::Annotation::Collection);
   my $commobj = Bio::Annotation::Comment->new();
   $commobj->text($comment);
   $annotation->add_Annotation('comment',$commobj);
   $seq->annotation($annotation);

   return $seq;
}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq($seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq object


=cut


=head2 write_qual

 Title   : write_qual
 Usage   : $stream->write_qual($seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq object


=cut

1;
