# $Id$
# BioPerl module for Bio::SeqIO::fasta
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#          and Lincoln Stein <lstein@cshl.org>
#
# Copyright Ewan Birney & Lincoln Stein
#
# You may distribute this module under the same terms as perl itself
# _history
# October 18, 1999  Largely rewritten by Lincoln Stein

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::fasta - fasta sequence input/output stream

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

package Bio::SeqIO::fasta;
use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::Object

use Bio::SeqIO;
use Bio::Seq::SeqFactory;

@ISA = qw(Bio::SeqIO);

sub _initialize {
  my($self,@args) = @_;
  $self->SUPER::_initialize(@args);  
  if( ! defined $self->sequence_factory ) {
      $self->sequence_factory(new Bio::Seq::SeqFactory(-verbose => $self->verbose(), -type => 'Bio::Seq'));      
  }
}

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    : NONE

=cut

sub next_seq {
    my( $self ) = @_;
    my $seq;
    my $alphabet;
    local $/ = "\n>";
    return unless my $entry = $self->_readline;

    if ($entry eq '>')  {	# very first one
	return unless $entry = $self->_readline;
    }
    
    my ($top,$sequence) = $entry =~ /^>?(.+?)\n([^>]*)/s
	or $self->throw("Can't parse fasta entry");
    my ($id,$fulldesc) = $top =~ /^\s*(\S+)\s*(.*)/
	or $self->throw("Can't parse fasta header");
    if ($id eq '') {$id=$fulldesc;} # FIX incase no space between > and name \AE
    $sequence =~ s/\s//g;	# Remove whitespace

    # for empty sequences we need to know the mol.type
    $alphabet = $self->alphabet();
    if(length($sequence) == 0) {
	if(! defined($alphabet)) {
	    # let's default to dna
	    $alphabet = "dna";
	}
    } else {
	# we don't need it really, so disable
	$alphabet = undef;
    }

    # create the seq object
    $seq = $self->{'_seqio_seqfactory'}->create_sequence(-seq        => $sequence,
						    -id         => $id,
						    -primary_id => $id,
						    -desc       => $fulldesc,
						    -alphabet    => $alphabet
						    );
    # if there wasn't one before, set the guessed type
    if( ! defined $alphabet ) {
	$self->alphabet($seq->alphabet());
    }
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
   foreach my $seq (@seq) {
     my $str = $seq->seq;
     my $top = $seq->display_id();
     if ($seq->can('desc') and my $desc = $seq->desc()) {
	 $desc =~ s/\n//g;
        $top .= " $desc";
     }
     if(length($str) > 0) {
	 $str =~ s/(.{1,60})/$1\n/g;
     } else {
	 $str = "\n";
     }
     $self->_print (">",$top,"\n",$str) or return;
   }
   return 1;
}

1;
