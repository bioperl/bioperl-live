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

This object can transform Bio::Align::AlignI objects to and from clustalw flat
file databases.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org               - General discussion
  http://bio.perl.org/MailList.html   - About the mailing lists


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

package Bio::AlignIO::clustalw;
use vars qw(@ISA $LINELENGTH);
use strict;

use Bio::AlignIO;
use Bio::LocatableSeq;
use Bio::SimpleAlign; # to be Bio::Align::Simple

$LINELENGTH = 60;

@ISA = qw(Bio::AlignIO);

=head2 new

 Title   : new
 Usage   : $alignio = new Bio::AlignIO(-format => 'clustalw', 
				       -file => 'filename');
 Function: returns a new Bio::AlignIO object to handle clustalw files
 Returns : Bio::AlignIO::clustalw object
 Args    : -verbose => verbosity setting (-1,0,1,2)
           -file    => name of file to read in or with ">" - writeout
           -fh      => alternative to -file param - provide a filehandle 
                       to read from/write to 
           -format  => type of Alignment Format to process or produce
           -percentages => (clustalw only) display a percentage of identity
                           in each line of the alignment.

           -linelength=> Set the alignment output line length (default 60)
=cut

sub _initialize {
    my ($self, @args) = @_;
    $self->SUPER::_initialize(@args);
    my ($percentages,
	$ll) = $self->_rearrange([qw(PERCENTAGES LINELENGTH)], @args);
    defined $percentages && $self->percentages($percentages);
    $self->line_length($ll || $LINELENGTH);
}

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream
 Returns : Bio::Align::AlignI object
 Args    : NONE

See L<Bio::Align::AlignI> for details

=cut

sub next_aln {
    my ($self) = @_;

    my $first_line;
    if( defined ($first_line  = $self->_readline ) 
	&& $first_line !~ /CLUSTAL/ ) {	
	$self->warn("trying to parse a file which does not start with a CLUSTAL header");
    }
    my %alignments;
    my $aln =  Bio::SimpleAlign->new(-source => 'clustalw');
    my $order = 0;
    my %order;
    $self->{_lastline} = '';
    while( defined ($_ = $self->_readline) ) {
	next if ( /^\s+$/ );	

	my ($seqname, $aln_line) = ('', '');	
	if( /^\s*(\S+)\s*\/\s*(\d+)-(\d+)\s+(\S+)\s*$/ ) {
	    # clustal 1.4 format
	    ($seqname,$aln_line) = ("$1:$2-$3",$4);
	} elsif( /^(\S+)\s+([A-Z\-]+)\s*$/ ) {
	    ($seqname,$aln_line) = ($1,$2);
	} else { $self->{_lastline} = $_; next }
	
	if( !exists $order{$seqname} ) {
	    $order{$seqname} = $order++;
	}

	$alignments{$seqname} .= $aln_line;  
    }
    my ($sname,$start,$end);
    foreach my $name ( sort { $order{$a} <=> $order{$b} } keys %alignments ) {
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
	$aln->add_seq($seq);
    }
    undef $aln if( !defined $end || $end <= 0);
    return $aln;
}

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the clustalw-format object (.aln) into the stream
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object


=cut

sub write_aln {
    my ($self,@aln) = @_;
    my ($count,$length,$seq,@seq,$tempcount,$line_len);
    $line_len = $self->line_length || $LINELENGTH;
    foreach my $aln (@aln) {
	if( ! $aln || ! $aln->isa('Bio::Align::AlignI')  ) { 
	    $self->warn("Must provide a Bio::Align::AlignI object when calling write_aln");
	    next;
	}
	my $matchline = $aln->match_line;
    
	$self->_print (sprintf("CLUSTAL W(1.81) multiple sequence alignment\n\n\n")) or return;

	$length = $aln->length();
	$count = $tempcount = 0;
	@seq = $aln->each_seq();
	my $max = 22;
	foreach $seq ( @seq ) {
	    $max = length ($aln->displayname($seq->get_nse())) 
		if( length ($aln->displayname($seq->get_nse())) > $max );
	}
	while( $count < $length ) {
	    foreach $seq ( @seq ) {
#
#  Following lines are to suppress warnings
#  if some sequences in the alignment are much longer than others.

		my ($substring);
		my $seqchars = $seq->seq();		
	      SWITCH: {
		  if (length($seqchars) >= ($count + $line_len)) {
		      $substring = substr($seqchars,$count,$line_len); 
		      last SWITCH; 
		  } elsif (length($seqchars) >= $count) {
		      $substring = substr($seqchars,$count); 
		      last SWITCH; 
		  }
		  $substring = "";
	      }
		
		$self->_print (sprintf("%-".$max."s %s\n",
				       $aln->displayname($seq->get_nse()),
				       $substring)) or return;
	    }		
	    
	    my $linesubstr = substr($matchline, $count,$line_len);
	    my $percentages = '';
	    if( $self->percentages ) {
		my ($strcpy) = ($linesubstr);
		my $count = ($strcpy =~ tr/\*//);
		$percentages = sprintf("\t%d%%", 100 * ($count / length($linesubstr)));
	    }
	    $self->_print (sprintf("%-".$max."s %s%s\n", '', $linesubstr,
				   $percentages));	    
	    $self->_print (sprintf("\n\n")) or return;
	    $count += $line_len;
	}
    }
    $self->_fh->flush if $self->_flush_on_write && defined $self->_fh;
    return 1;
}

=head2 percentages

 Title   : percentages
 Usage   : $obj->percentages($newval)
 Function: Set the percentages flag - whether or not to show percentages in 
           each output line
 Returns : value of percentages
 Args    : newvalue (optional)


=cut

sub percentages { 
    my ($self,$value) = @_; 
    if( defined $value) {
	$self->{'_percentages'} = $value; 
    } 
    return $self->{'_percentages'}; 
}

=head2 line_length

 Title   : line_length
 Usage   : $obj->line_length($newval)
 Function: Set the alignment output line length
 Returns : value of line_length
 Args    : newvalue (optional)


=cut

sub line_length {
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_line_length'} = $value;
    }
    return $self->{'_line_length'};
}

1;
