# $Id$
#
# BioPerl module for Bio::AlignIO::emboss
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::emboss - Parse EMBOSS alignment output (from applications water and needle)

=head1 SYNOPSIS

    # do not use the object directly
    use Bio::AlignIO;
    # read in an alignment from the EMBOSS program water
    my $in = new Bio::AlignIO(-format => 'emboss',
			      -file   => 'seq.water');
    while( my $aln = $in->next_aln ) {
	# do something with the alignment
    }

=head1 DESCRIPTION

This object handles parsing and writing pairwise sequence alignments
from the EMBOSS suite.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::AlignIO::emboss;
use vars qw(@ISA);
use strict;

use Bio::AlignIO;
use Bio::LocatableSeq;

@ISA = qw(Bio::AlignIO );

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
 Returns : Bio::SimpleAlign object - returns 0 on end of file
	    or on error
 Args    : NONE

=cut

sub next_aln {
    my ($self) = @_;
    my %data;    
    while( defined($_ = $self->_readline) ) {
        next if( /^\s+/ );
	if( /(Local|Global):\s*(\S+)\s+vs\s+(\S+)/ ) {
	    $data{'type'} = $1 eq 'Local' ? 'water' : 'needle';
	    $data{'seq1'} = { 
		'start'=> undef,
		'end'=> undef,		
		'name' => $2,
		'data' => '' };
	    $data{'seq2'} = { 
		'start'=> undef,
		'end'=> undef,		
		'name' => $3,
		'data' => '' };
	    $data{'align'} = '';
	} elsif( /Score:\s+(\S+)/ ) {
	    $data{'score'} = $1;		
	} elsif( /^$data{'seq1'}->{'name'}/ ) {
	    my $count = 0;
	    while( defined ($_) ) {
		if($count == 0 || $count == 2 ) {
		    my ($seq,$start,$align,$end) = split;
		    my $seqname = sprintf("seq%d", ($count == 0) ? '1' : '2'); 
		    $data{$seqname}->{'data'} .= $align;
		    $data{$seqname}->{'start'} ||= $start;
		    $data{$seqname}->{'end'} = $end;
		} else { 
		    s/^\s+//;
		    s/\s+$//;
		    $data{'align'} .= $_;
		}
		last if( $count++ == 2);
		$_ = $self->_readline();
	    }
	}
    }
    my $aln =  Bio::SimpleAlign->new();
    
    foreach my $seqname ( qw(seq1 seq2) ) { 
	return undef unless ( defined $data{$seqname} );
	my $seq = new Bio::LocatableSeq('-seq' => $data{$seqname}->{'data'},
					'-id'  => $data{$seqname}->{'name'},
					'-start'=> $data{$seqname}->{'start'},
					'-end' => $data{$seqname}->{'end'},
					'-type' => 'aligned');
	$aln->add_seq($seq);
    }
    return $aln;
}
    

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in emboss format
 Returns : 1 for success and 0 for error
 Args    : Bio::Align::AlignI object


=cut

sub write_aln {
    my ($self,@aln) = @_;

    $self->throw("Sorry: writing emboss output is not currently available! \n");
}

1;
