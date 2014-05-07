#
# BioPerl module for Bio::AlignIO::msf
#   based on the Bio::SeqIO::msf module
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

Bio::AlignIO::msf - msf sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::AlignIO> class.

=head1 DESCRIPTION

This object can transform L<Bio::Align::AlignI> objects to and from msf
flat file databases.

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

package Bio::AlignIO::msf;
use vars qw(%valid_type);
use strict;

use Bio::SeqIO::gcg; # for GCG_checksum()
use Bio::SimpleAlign;

use base qw(Bio::AlignIO);

BEGIN {
	%valid_type = qw( dna N rna N protein P );
}

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream. Tries to read *all* MSF
           It reads all non whitespace characters in the alignment
           area. For MSFs with weird gaps (eg ~~~) map them by using
           $aln->map_chars('~','-')
 Returns : Bio::Align::AlignI object
 Args    : NONE

=cut

sub next_aln {
    my $self = shift;
    my $entry;
    my (%hash,$name,$str,@names,$seqname,$start,$end,$count,$seq);

    my $aln =  Bio::SimpleAlign->new(-source => 'gcg' );

    while( $entry = $self->_readline) {
	$entry =~ m{//} && last; # move to alignment section
	$entry =~ /Name:\s+(\S+)/ && do { $name = $1;
					  $hash{$name} = ""; # blank line
					  push(@names,$name); # we need it ordered!
				      };
	# otherwise - skip
    }

    # alignment section

    while( $entry = $self->_readline) {
	next if ( $entry =~ /^\s+(\d+)/ ) ;
	$entry =~ /^\s*(\S+)\s+(.*)$/ && do {
	    $name = $1;
	    $str = $2;
	    if( ! exists $hash{$name} ) {
		$self->throw("$name exists as an alignment line but not in the header. Not confident of what is going on!");
	    }
	    $str =~ s/\s//g;
	    $str =~ s/~/-/g;
	    $hash{$name} .= $str;
	};
    }

    return if @names < 1;

    # now got this as a name - sequence hash. Let's make some sequences!

    for $name ( @names ) {
	if( $name =~ m{(\S+)/(\d+)-(\d+)} ) {
	    $seqname = $1;
	    $start = $2;
	    $end = $3;
	} else {
	    $seqname = $name;
	    $start = 1;
	    $str = $hash{$name};
	    $str =~ s/[^0-9A-Za-z$Bio::LocatableSeq::OTHER_SYMBOLS]//g;

	    $end = length($str);
	}

	$seq = Bio::LocatableSeq->new('-seq'        => $hash{$name},
				      '-display_id' => $seqname,
				      '-start'      => $start,
				      '-end'        => $end,
				      '-alphabet'   => $self->alphabet,
				      );
	$aln->add_seq($seq);

#  If $end <= 0, we have either reached the end of
#  file in <> or we have encountered some other error

    }
    return $aln if $aln->num_sequences;
    return;
}


=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in MSF format
           Sequence type of the alignment is determined by the first sequence.
 Returns : 1 for success and 0 for error
 Args    : Bio::Align::AlignI object


=cut

sub write_aln {
	my ($self,@aln) = @_;
	my $msftag;
	my $type;
	my $count = 0;
	my $maxname;
	my ($length,$date,$name,$seq,$miss,$pad,%hash,@arr,$tempcount,$index);
	foreach my $aln (@aln) {
		if( ! $aln || ! $aln->isa('Bio::Align::AlignI')  ) {
			$self->warn("Must provide a Bio::Align::AlignI object when calling write_aln");
			next;
		}
		$date = localtime(time);
		$msftag = "MSF";
		$type = $valid_type{$aln->get_seq_by_pos(1)->alphabet};
		$maxname = $aln->maxdisplayname_length();
		$length  = $aln->length();
		$name = $aln->id();
		if( !defined $name ) {
			$name = "Align";
		}

		$self->_print (sprintf("\n%s   MSF: %d  Type: %s  %s  Check: 00 ..\n\n",
			       $name,  $aln->num_sequences, $type, $date));

    my $seqCountFormat = "%".($maxname > 20 ? $maxname + 2: 22)."s%-27d%27d\n";
    my $seqNameFormat  = "%-".($maxname > 20 ? $maxname : 20)."s  ";
        
		foreach $seq ( $aln->each_seq() ) {
			$name = $aln->displayname($seq->get_nse());
			$miss = $maxname - length ($name);
			$miss += 2;
			$pad  = " " x $miss;

			$self->_print (sprintf(" Name: %s%sLen:    %d  Check:  %d  Weight:  1.00\n",$name,$pad,length $seq->seq(), Bio::SeqIO::gcg->GCG_checksum($seq)));

			$hash{$name} = $seq->seq();
			push(@arr,$name);
		}
    	# ok - heavy handed, but there you go.
    	#
    	$self->_print ("\n//\n\n\n");

    	while( $count < $length ) {
			# there is another block to go!
			$self->_print (sprintf($seqCountFormat,' ',$count+1,$count+50));
			foreach $name  ( @arr ) {
				$self->_print (sprintf($seqNameFormat,$name));

				$tempcount = $count;
				$index = 0;
				while( ($tempcount + 10 < $length) && ($index < 5)  ) {

					$self->_print (sprintf("%s ",substr($hash{$name},
																	$tempcount,10)));

					$tempcount += 10;
					$index++;
				}	    	#
				# ok, could be the very last guy ;)
				#
				if( $index < 5) {
					# space to print!
					#
					$self->_print (sprintf("%s ",substr($hash{$name},$tempcount)));
					$tempcount += 10;
				}
				$self->_print ("\n");
			}
			$self->_print ("\n\n");
			$count = $tempcount;
    	}
	}
	$self->flush if $self->_flush_on_write && defined $self->_fh;
	return 1;
}

1;
