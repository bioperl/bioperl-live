# $Id$
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

This object can transform Bio::Seq objects to and from pir flat
file databases.

Note: This does not completely preserve the PIR format - quality 
information about sequence is currently discarded since bioperl 
does not have a mechanism for handling these encodings in sequence 
data.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://www.bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHORS

Aaron Mackey E<lt>amackey@virginia.eduE<gt>
Lincoln Stein E<lt>lstein@cshl.orgE<gt>
Jason Stajich E<lt>jason@chg.mc.duke.eduE<gt>


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::pir;
use vars qw(@ISA);
use strict;

use Bio::SeqIO;
use Bio::PrimarySeq;
use Bio::Seq;

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

 Title   : next_primary_seq
 Usage   : $seq = $stream->next_primary_seq()
 Function: returns the next sequence in the stream as a Bio::PrimarySeq
 Returns : Bio::Seq object
 Args    :


=cut

sub next_primary_seq {
    my ($self,$as_next_seq) = @_;
    local $/ = "\n>";
    return unless my $line = $self->_readline;
    if( $line eq '>' ) {	# handle the very first one having no comment
	return unless $line = $self->_readline;
    }
    my ($top, $desc,$seq) = ( $line =~ /^(.+?)\n(.+?)\n([^>]*)/s )  or
	$self->throw("Cannot parse entry PIR entry [$line]");

    
    my ( $type,$id ) = ( $top =~ /^>?([PF])1;(\S+)\s*$/ ) or
	$self->throw("PIR stream read attempted without leading '>P1;' [ $line ]");

    # P - indicates complete protein
    # F - indicates protein fragment
    # not sure how to stuff these into a Bio object 
    # suitable for writing out.
    $seq =~ s/\*//g;
    $seq =~ s/[\(\)\.\/\=\,]//g;
    $seq =~ s/\s+//g;		# get rid of whitespace

    my ($alphabet,$seqobj) = ('protein',undef);
    # TODO - not processing SFS data
    if ($as_next_seq) {
	# Return a Bio::Seq if asked for
	$seqobj = Bio::Seq->new(-seq        => $seq,
				-primary_id => $id,
				-id         => $type. '1;' . $id,
				-desc       => $desc,
				-alphabet    => $alphabet
				);
    } else {
	$seqobj = Bio::PrimarySeq->new(-seq        => $seq,
				       -primary_id => $id,
				       -id         => $type. '1;' . $id,
				       -desc       => $desc,
				       -alphabet    => $alphabet
				       );
    }
    return $seqobj;
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
     return unless $self->_print(">".$seq->id(), 
				 "\n", $seq->desc(), "\n", 
				 $str, "*\n");
   }
   return 1;
}

1;
