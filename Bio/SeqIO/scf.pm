#-----------------------------------------------------------------------------
# PACKAGE : Bio::SeqIO::scf
# AUTHOR  : Aaron Mackey <amackey@virginia.edu>
# CREATED : Feb 16 1999
# REVISION: $Id$
#            
# Copyright (c) 1997-9 bioperl, Aaron Mackeyy. All Rights Reserved.
#           This module is free software; you can redistribute it and/or 
#           modify it under the same terms as Perl itself.
#
# _History_
#
# Ewan Birney <birney@sanger.ac.uk> developed the SeqIO 
# schema and the first prototype modules.
#
# This code is based on his Bio::SeqIO::Fasta module with
# the necessary minor tweaks necessary to get it to read
# and write raw formatted sequences made by
# chris dagdigian <dag@sonsorol.org>
#
# October 18, 1999 Largely rewritten by Lincoln Stein
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::scf - SCF sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::SeqIO> class.

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from SCF trace
files.

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

=head1 AUTHORS

  Aaron Mackey
  Lincoln Stein lstein@cshl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::SeqIO::scf;
use vars qw(@ISA);
use strict;
use Bio::SeqIO;
# Object preamble - inherits from Bio::Root::Object

@ISA = qw(Bio::SeqIO);

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  return unless my $make = $self->SUPER::_initialize(@args);
}

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    :


=cut

sub next_seq{
   my ($self,@args) = @_;
   my ($seq, $seqc, $fh, $buffer, %names);

   $fh = $self->_filehandle();
   unless ($fh) {  # simulate the <> function
     if ( !fileno(ARGV) or eof(ARGV) ) {
       return unless my $ARGV = shift;
       open(ARGV,$ARGV) or
	 $self->throw("Could not open $ARGV for SCF stream reading $!");
     }
     $fh = \*ARGV;
   }

   binmode $fh; # for the Win32/Mac crowds
   return unless read $fh, $buffer, 128;  # no exception; probably end of file

   my ($scf,
       $sample,
       $sample_offset,
       $bases,
       $bases_left_clip,
       $bases_right_clip,
       $bases_offset,
       $comment_size,
       $comments_offset,
       $version,
       $sample_size,
       $code_set,
       @header_spare ) = unpack "a4 NNNNNNNN a4 NN N20", $buffer;

   my $bytes_to_read = $bases_offset - 128;
   $self->throw('Unexpected end of file while reading from SCF file')
     unless $bytes_to_read == read $fh, $buffer, $bytes_to_read; # skip over sample info

   $seqc = '';
   for (0 .. ($bases - 1) ) {
     $self->throw('Unexpected end of file while reading from SCF file')
       unless 12 == read(STDIN, $buffer, 12);
     my ($index,
	 $prob_A,
	 $prob_C,
	 $prob_G,
	 $prob_T,
	 $base,
	 @base_spare) = unpack "N C C C C a C3", $buffer;
     $seqc .= $base;
   }

   {
     local ($/) = "\n";
     while (<$fh>) {
       m/([^\=]+)\=(.*)/;
       $names{$1} = $2;
     }
   }

   $seq = Bio::Seq->new(-seq => $seqc,
			-id => $names{'NAME'},
			-names => \%names,
		       );
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
   $self->throw("Sorry, Bioperl cannot yet write SCF format!");
}

1;
