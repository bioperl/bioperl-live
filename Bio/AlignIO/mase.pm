# $Id$
#
# BioPerl module for Bio::AlignIO::mase

#	based on the Bio::SeqIO::mase module
#       by Ewan Birney <birney@sanger.ac.uk>
#       and Lincoln Stein  <lstein@cshl.org>
#
#       and the SimpleAlign.pm module of Ewan Birney
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself
# _history
# September 5, 2000
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::mase - mase sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::AlignIO class.

=head1 DESCRIPTION

This object can transform Bio::SimpleAlign objects to and from mase flat
file databases.

=head1 FEEDBACK

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHORS - Peter Schattner

Email: schattner@alum.mit.edu


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::mase;
use vars qw(@ISA);
use strict;

use Bio::AlignIO;

@ISA = qw(Bio::AlignIO);


=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
 Returns : SimpleAlign object
 Args    : NONE

=cut

sub next_aln {
    my $self = shift;
    my $entry;
    my $name;
    my $start;
    my $end;
    my $seq;
    my $add;
    my $count = 0;
    my $seq_residues;

    my $aln =  Bio::SimpleAlign->new();


    while( $entry = $self->_readline) {
        $entry =~ /^;/ && next;
	if(  $entry =~ /^(\S+)\/(\d+)-(\d+)/ ) {
	    $name = $1;
	    $start = $2;
	    $end = $3;
	} else {
	    $entry =~ s/\s//g;
	    $name = $entry;
	    $end = -1;
	}

	$seq = "";

	while( $entry = $self->_readline) {
	    $entry =~ /^;/ && last;
	    $entry =~ s/[^A-Za-z\.\-]//g;
	    $seq .= $entry;
	}
	if( $end == -1) {
	    $start = 1;

	    $seq_residues = $seq;
	    $seq_residues =~ s/\W//g;
	    $end = length($seq_residues);
	}

	$add = new Bio::LocatableSeq('-seq'=>$seq,
			    '-id'=>$name,
			    '-start'=>$start,
			    '-end'=>$end,
			    );


       $aln->add_seq($add);


#  If $end <= 0, we have either reached the end of
#  file in <> or we have encountered some other error
#
   if ($end <= 0) { undef $aln;}

   }

   return $aln;
}



=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in mase format  ###Not yet implemented!###
 Returns : 1 for success and 0 for error
 Args    : Bio::SimpleAlign object


=cut

sub write_aln {
    my ($self,@aln) = @_;

    $self->throw("Sorry: mase-format output, not yet implemented! /n");
}

1;
