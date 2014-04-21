#
# BioPerl module for Bio::AlignIO::mase

#   based on the Bio::SeqIO::mase module
#       by Ewan Birney <birney@ebi.ac.uk>
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

Do not use this module directly.  Use it via the L<Bio::AlignIO> class.

=head1 DESCRIPTION

This object can transform L<Bio::Align::AlignI> objects to and from mase flat
file databases.

=head1 FEEDBACK

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHORS - Peter Schattner

Email: schattner@alum.mit.edu


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::mase;
use strict;


use base qw(Bio::AlignIO);


=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
 Returns : L<Bio::Align::AlignI> object
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

    my $aln =  Bio::SimpleAlign->new(-source => 'mase');


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

	$add = Bio::LocatableSeq->new('-seq'        => $seq,
				      '-display_id' => $name,
				      '-start'      => $start,
				      '-end'        => $end,
				      '-alphabet'   => $self->alphabet,
				      );


       $aln->add_seq($add);


#  If $end <= 0, we have either reached the end of
#  file in <> or we have encountered some other error
#
   if ($end <= 0) { undef $aln;}

   }

   return $aln if $aln->num_sequences;
   return;
}



=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in mase format  ###Not yet implemented!###
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object


=cut

sub write_aln {
    my ($self,@aln) = @_;
    $self->throw_not_implemented();
}

1;
