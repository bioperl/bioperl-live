# $Id$
#
# BioPerl module for Bio::AlignIO::clustalw

#	based on the Bio::SeqIO modules
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

Bio::AlignIO::clustalw - clustalw sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::AlignIO class.

=head1 DESCRIPTION

This object can transform Bio::SimpleAlign objects to and from clustalw flat
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

package Bio::AlignIO::clustalw;
use vars qw(@ISA);
use strict;

use Bio::AlignIO;
use Bio::LocatableSeq;

@ISA = qw(Bio::AlignIO);

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream
 Returns : SimpleAlign object
 Args    : NONE

=cut

sub next_aln {
    my ($self) = @_;

    my $first_line;
    if( defined ($first_line  = $self->_readline ) 
	&& $first_line !~ /CLUSTAL/ ) {	
	$self->warn("trying to parse a file which does not start with a CLUSTAL header");
    }
    my %alignments;
    my $aln =  Bio::SimpleAlign->new();
    while( defined ($_ = $self->_readline) ) {
	next if ( /^\s+$/ );	

	my ($seqname, $aln_line) = ('', '');	
	if( /^\s*(\S+)\s*\/\s*(\d+)-(\d+)\s+(\S+)\s*$/ ) {
	    # clustal 1.4 format
	    ($seqname,$aln_line) = ("$1:$2-$3",$4);
	} elsif( /^(\S+)\s+([A-Z\-]+)\s*$/ ) {
	    ($seqname,$aln_line) = ($1,$2);
	} else { next }
	$alignments{$seqname} .= $aln_line;  
    }
    my ($sname,$start,$end);
    foreach my $name ( keys %alignments ) {
	if( $name =~ /(\S+):(\d+)-(\d+)/ ) {
	    ($sname,$start,$end) = ($1,$2,$3);	    
	} else {
	    ($sname, $start) = ($name,1);
	    my $str  = $alignments{$name};
	    $str =~ s/[^A-Za-z]//g;
	    $end = length($str);
	}
	my $seq = new Bio::LocatableSeq('-seq'   => $alignments{$name},
					 '-id'    => $sname,
					 '-start' => $start,
					 '-end'   => $end);
	$aln->addSeq($seq);
    }
    undef $aln if( !defined $end || $end <= 0);
    return $aln;
}



=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the clustalw-format object (.aln) into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::SimpleAlign object


=cut

sub write_aln {
   my ($self,@aln) = @_;
   my ($count,$length,$seq,@seq,$tempcount);

   foreach my $aln (@aln) {

    $self->_print (sprintf("CLUSTAL W(1.4) multiple sequence alignment\n\n\n")) or return;

    $length = $aln->length_aln();
    $count = $tempcount = 0;
    @seq = $aln->eachSeq();

    while( $count < $length ) {
	foreach $seq ( @seq ) {
#
#  Following lines are to suppress warnings
#  if some sequences in the alignment are much longer than others.
	   my $substring;
           my $seqchars = $seq->seq();
	SWITCH: {
		if (length($seqchars) >= ($count + 50)) {
			$substring = substr($seqchars,$count,50); last SWITCH; }
		if (length($seqchars) >= $count) {
			$substring = substr($seqchars,$count); last SWITCH; }
                $substring = "";
	}

	   $self->_print (sprintf("%-22s %s\n",$aln->get_displayname($seq->get_nse()),$substring)) or return;
	}
	$self->_print (sprintf("\n\n")) or return;
	$count += 50;
    }
   }
   return 1;
}

1;
