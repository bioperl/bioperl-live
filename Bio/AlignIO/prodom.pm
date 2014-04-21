#
# BioPerl module for Bio::AlignIO::prodom

#   based on the Bio::SeqIO::prodom module
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

Bio::AlignIO::prodom - prodom sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::AlignIO> class.

=head1 DESCRIPTION

This object can transform L<Bio::Align::AlignI> objects to and from prodom flat
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

package Bio::AlignIO::prodom;
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
    my ($acc, $fake_id, $start, $end, $seq, $add, %names);

    my $aln =  Bio::SimpleAlign->new(-source => 'prodom');

    while( $entry = $self->_readline) {

       if ($entry =~ /^AC\s+(\S+)\s*$/) {         #ps 9/12/00
	   $aln->id( $1 );
       }
       elsif ($entry =~ /^AL\s+(\S+)\|(\S+)\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s*$/){    #ps 9/12/00
	   $acc=$1;
	   $fake_id=$2;  # Accessions have _species appended
	   $start=$3;
	   $end=$4;
	   $seq=$5;

	   $names{'fake_id'} = $fake_id;

	   $add = Bio::LocatableSeq->new('-seq'      => $seq,
					 '-id'       => $acc,
					 '-start'    => $start,
					 '-end'      => $end,
					 '-alphabet' => $self->alphabet,
			       );

	   $aln->add_seq($add);
       }
       elsif ($entry =~ /^CO/) {
	   # the consensus line marks the end of the alignment part of the entry
	   last;
       }
	}
	
   return $aln if $aln->num_sequences;
   return;
}



=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in prodom format  ###Not yet implemented!###
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object


=cut

sub write_aln {
    my ($self,@aln) = @_;
    $self->throw_not_implemented();
}

1;
