# $Id$
#
# BioPerl module for Bio::AlignIO::fasta

#	based on the Bio::SeqIO::fasta module
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

Bio::AlignIO::fasta - fasta sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::AlignIO class.

=head1 DESCRIPTION

This object can transform Bio::SimpleAlign objects to and from fasta flat
file databases.

=head1 FEEDBACK

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHORS - Peter Schattner

Email: schattner@alum.mit.edu


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::fasta;
use vars qw(@ISA);
use strict;

use Bio::AlignIO;

@ISA = qw(Bio::AlignIO);


=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
 Returns : SimpleAlign object - returns 0 on end of file
	    or on error
 Args    : NONE

=cut

sub next_aln {
    my $self = shift;
    my $entry;
    my ($start,$end,$name,$seqname,$seq,$seqchar,$tempname,%align);
    my $aln =  Bio::SimpleAlign->new();

    while(defined ($entry = $self->_readline)) {
	if($entry =~ /^>(\S+)/ ) {
	    $tempname = $1;
	    if( defined $name ) {
		# put away last name and sequence

		if( $name =~ /(\S+)\/(\d+)-(\d+)/ ) {
		    $seqname = $1;
		    $start = $2;
		    $end = $3;
		} else {
		    $seqname=$name;
		    $start = 1;
		    $end = length($seqchar);   #ps 9/6/00
		}
		
		$seq = new Bio::LocatableSeq('-seq'=>$seqchar,
				    '-id'=>$seqname,
				    '-start'=>$start,
				    '-end'=>$end,
				    );

		$aln->add_seq($seq);
	     }
	     $name = $tempname;
	     $seqchar  = "";
	     next;
	}
	$entry =~ s/[^A-Za-z\.\-]//g;
	$seqchar .= $entry;

    }
#
#  Next two lines are to silence warnings that
#  otherwise occur at EOF when using <$fh>

   if (!defined $name) {$name="";}
   if (!defined $seqchar) {$seqchar="";}

#  Put away last name and sequence
    if( $name =~ /(\S+)\/(\d+)-(\d+)/ ) {
	$seqname = $1;
	$start = $2;
	$end = $3;
    } else {
	$seqname=$name;
	$start = 1;
	$end = length($seqchar);   #ps 9/6/00
#	$end = length($align{$name});
    }


#  If $end <= 0, we have either reached the end of
#  file in <> or we have encountered some other error
#
   if ($end <= 0) { undef $aln; return $aln;}


    $seq = new Bio::LocatableSeq('-seq'=>$seqchar,
			'-id'=>$seqname,
			'-start'=>$start,
			'-end'=>$end,
			);

    $aln->add_seq($seq);

    return $aln;

}
	

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in fasta format
 Returns : 1 for success and 0 for error
 Args    : Bio::SimpleAlign object


=cut

sub write_aln {
    my ($self,@aln) = @_;
    my ($seq,$rseq,$name,$count,$length,$seqsub);

  foreach my $aln (@aln) {

      foreach $rseq ( $aln->each_seq() ) {
	$name = $aln->displayname($rseq->get_nse());
	$seq  = $rseq->seq();	
	$self->_print (">$name\n") or return ;	
	$count =0;
	$length = length($seq);
	while( ($count * 60 ) < $length ) {
	    $seqsub = substr($seq,$count*60,60);
	    $self->_print ("$seqsub\n") or return ;
	    $count++;
	}
      }
   }
   return 1;
}

1;
