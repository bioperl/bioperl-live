# $Id$
#
# BioPerl module for Bio::AlignIO::selex

#	based on the Bio::SeqIO::selex module
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

Bio::AlignIO::selex - selex sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::AlignIO> class.

=head1 DESCRIPTION

This object can transform L<Bio::Align::AlignI> objects to and from selex flat
file databases.

=head1 FEEDBACK

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHORS - Peter Schattner

Email: schattner@alum.mit.edu


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::selex;
use vars qw(@ISA);
use strict;
use Bio::AlignIO;

@ISA = qw(Bio::AlignIO);

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream. Tries to read *all* selex
          It reads all non whitespace characters in the alignment
          area. For selexs with weird gaps (eg ~~~) map them by using
          $al->map_chars('~','-')
 Returns : L<Bio::Align::AlignI> object
 Args    : NONE

=cut

sub next_aln {
    my $self = shift;
    my $entry;
    my ($start,$end,%align,$name,$seqname,$seq,$count,%hash,%c2name, %accession, $no);
    my $aln =  Bio::SimpleAlign->new(-source => 'selex');

    # in selex format, every non-blank line that does not start
    # with '#=' is an alignment segment; the '#=' lines are mark up lines.
    # Of particular interest are the '#=GF <name/st-ed> AC <accession>'
    # lines, which give accession numbers for each segment

    while( $entry = $self->_readline) {
        $entry =~ /^\#=GS\s+(\S+)\s+AC\s+(\S+)/ && do {
	    				$accession{ $1 } = $2;
	    				next;
					};
	$entry !~ /^([^\#]\S+)\s+([A-Za-z\.\-]+)\s*/ && next;
	
	$name = $1;
	$seq = $2;

	if( ! defined $align{$name}  ) {
	    $count++;
	    $c2name{$count} = $name;
	}
	$align{$name} .= $seq;
    }

    # ok... now we can make the sequences

    $count = 0;
    foreach $no ( sort { $a <=> $b } keys %c2name ) {
	$name = $c2name{$no};

	if( $name =~ /(\S+)\/(\d+)-(\d+)/ ) {
	    $seqname = $1;
	    $start = $2;
	    $end = $3;
	} else {
	    $seqname=$name;
	    $start = 1;
	    $end = length($align{$name});
	}
	$seq = new Bio::LocatableSeq('-seq'=>$align{$name},
			    '-id'=>$seqname,
			    '-start'=>$start,
			    '-end'=>$end,
			    '-type'=>'aligned',
				     '-accession_number' => $accession{$name},

			    );

	$aln->add_seq($seq);
	$count++;
    }

#  If $end <= 0, we have either reached the end of
#  file in <> or we have encountered some other error
#
   if ($end <= 0) { undef $aln;}

    return $aln;
}


=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in selex format
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object


=cut

sub write_aln {
    my ($self,@aln) = @_;
    my ($namestr,$seq,$add);
    my ($maxn);
    foreach my $aln (@aln) {
	$maxn = $aln->maxdisplayname_length();
	foreach $seq ( $aln->each_seq() ) {
	    $namestr = $aln->displayname($seq->get_nse());
	    $add = $maxn - length($namestr) + 2;
	    $namestr .= " " x $add;
	    $self->_print (sprintf("%s  %s\n",$namestr,$seq->seq())) or return;
	}
    }
    $self->flush if $self->_flush_on_write && defined $self->_fh;
    return 1;
}

1;
