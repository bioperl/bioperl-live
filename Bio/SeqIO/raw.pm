#-----------------------------------------------------------------------------
# PACKAGE : Bio::SeqIO::raw
# AUTHOR  : Ewan Birney <birney@ebi.ac.uk>
# CREATED : Feb 16 1999
# REVISION: $Id$
#            
# Copyright (c) 1997-9 bioperl, Ewan Birney. All Rights Reserved.
#           This module is free software; you can redistribute it and/or 
#           modify it under the same terms as Perl itself.
#
# _History_
#
# Ewan Birney <birney@ebi.ac.uk> developed the SeqIO 
# schema and the first prototype modules.
#
# This code is based on his Bio::SeqIO::Fasta module with
# the necessary minor tweaks necessary to get it to read
# and write raw formatted sequences made by
# chris dagdigian <dag@sonsorol.org>
#
# October 18, 1999  Largely rewritten by Lincoln Stein
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::raw - raw sequence file input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::SeqIO> class.

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from raw flat
file databases.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  bioperl-guts-l@bioperl.org            - Technically-oriented discussion
  http://www.bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHORS

  Ewan Birney  <birney@ebi.ac.uk>
  Lincoln Stein <lstein@cshl.org>

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqIO::raw;
use vars '@ISA';
use strict;
use Bio::SeqIO;
# Object preamble - inheriets from Bio::Root::Object

use FileHandle;

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
   ## When its 1 sequence per line with no formatting at all,
   ## grabbing it should be easy :)

   my $nextline = $self->_readline();
   if( !defined $nextline ){ return undef; }

   my $sequence = uc($nextline);
   $sequence =~ s/\W//g;

   return  Bio::Seq->new(-seq => $sequence);

}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq($seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq object


=cut

sub write_seq {
   my ($self,@seq) = @_;
   foreach (@seq) {
     $self->_print($_->seq, "\n") or return;
   }
   return 1;
}

1;
