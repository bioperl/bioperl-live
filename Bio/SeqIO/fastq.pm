# BioPerl module for Bio::SeqIO::fastq
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Tony Cox <avc@sanger.ac.uk>
#
# Copyright Tony Cox
#
# You may distribute this module under the same terms as perl itself
#
# _history
#
# October 29, 2001  incept data
# June 20, 2009 updates for Illumina variant FASTQ formats for Solexa and later

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::fastq - fastq sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::SeqIO class.

=head1 DESCRIPTION

This object can transform Bio::Seq and Bio::Seq::Quality objects to and from
fastq flat file databases.

Fastq is a file format used frequently at the Sanger Centre to bundle a fasta
sequence and its quality data. A typical fastaq entry takes the from:

  @HCDPQ1D0501
  GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT.....
  +HCDPQ1D0501
  !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65.....

Fastq files have sequence and quality data on a single line and the
quality values are single-byte encoded.

This parser now has preliminary support for Illumina v 1.0 and 1.3 FASTQ, though
it needs extensive testing.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 
 
Please direct usage questions or support issues to the mailing list:
  
L<bioperl-l@bioperl.org>
  
rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHORS - Tony Cox

Email: avc@sanger.ac.uk


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::fastq;
use strict;

use Bio::Seq::SeqFactory;

use base qw(Bio::SeqIO);

sub _initialize {
    my($self,@args) = @_;
    $self->SUPER::_initialize(@args);
    my ($variant, $validate) = $self->_rearrange([qw(VARIANT
                                                   VALIDATE)], @args);
    $variant ||= 'sanger';
    $self->variant($variant);
    $validate   && $self->validate($validate);
    
    if( ! defined $self->sequence_factory ) {
        $self->sequence_factory(Bio::Seq::SeqFactory->new(-verbose => $self->verbose(), -type => 'Bio::Seq::Quality'));      
    }
}


=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq::Quality object
 Args    : NONE
 Status  : Stable
 
=cut

sub next_seq {
    my( $self ) = @_;
    while (defined(my $data = $self->next_dataset)) {
        # Are FASTQ sequences w/o any sequence valid?  Removing for now
        # -cjfields 6.22.09
        my $seq = $self->sequence_factory->create(%$data);
        return $seq;
    }
    return;
}

# pure perl version
sub next_dataset {
    my $self = shift;
    local $/ = "\n";
    my ($data, $ct);
    my $mode = '-seq';
    
    # speed this up by directly accessing the filehandle and in-lining the
    # _readline stuff vs. making the repeated method calls. Tradeoff is speed
    # over repeated code.
    
    # we can probably normalize line endings using PerlIO::eol or
    # Encode::Newlines
    
    my $fh = $self->_fh;
    my $line = $self->{lastline} || <$fh>;
    
    FASTQ:
    while ($line) {
        $line =~ s/\015\012/\012/;
        $line =~ tr/\015/\n/;
        if ($mode eq '-seq' && $line =~ m{^@([^\n]+)$}xmso) {
            $data->{-descriptor} = $1;
            my ($id,$fulldesc);
            if ($data->{-descriptor} =~ /^\s*(\S+)\s*(.*)/) {
                ($id,$fulldesc) = ($1, $2);
            } else {
                $self->throw("Can't parse fastq header");
            }
            $data->{-id} = $id;
            $data->{-desc} = $fulldesc;
            $ct->{-seq} = 0;
        } elsif ($mode eq '-seq' && $line =~ m{^\+([^\n]*)}xmso) {
			my $desc = $1;
            $self->throw("No description line parsed") unless $data->{-descriptor};
            if ($desc && $data->{-descriptor} ne $desc) {
                $self->throw("Quality descriptor [$desc] doesn't match seq description ".$data->{-descriptor} );
            }
            $mode = '-raw_quality';
            $ct->{-raw_quality} = 0;
        } else {
            # this uses a fairly loose check where the number of lines of
            # both qual and seq match before bailing from the loop
            if ($mode eq '-raw_quality' &&
                $ct->{-raw_quality} == $ct->{-seq}) {
                $self->{lastline} = $line;
                last FASTQ
            }
            $ct->{$mode}++;
            chomp $line;
            $data->{$mode} .= $line;
        }
        $line = <$fh>;
        if (!defined $line) {
            delete $self->{lastline};
            last FASTQ;
        }
    }
    
    # simple quality control tests
    if (defined $data) {
        if (length $data->{-seq} != length $data->{-raw_quality}) {
            $self->throw("Quality string\n".$data->{-raw_quality}."\nlength [".length($data->{-raw_quality})."] doesn't match ".
                         "length of sequence\n".$data->{-seq}."\n[".length($data->{-seq})."], $.");
        }
        
		# need to benchmark this vs the prior unpack("C", ...) version
        my @qual = map {$self->{chr2phred}->{$_}}
            unpack("A1" x length($data->{-raw_quality}), $data->{-raw_quality});
            
        # this should be somewhat rarer now, but needs to be present JIC
        if ($self->variant eq 'solexa') {
            # The conversion needs to be to PHRED score, but solexa (aka illumina 1.0)
            # has Solexa qual units, not PHRED qual units.  Convert over...
            # this doesn't account for very low scores yet!
			# NOTE: code is kludged from MAQ (maq.sourceforge.net)
            @qual = map {sprintf("%.0f",(10 * log(1 + 10 ** ($_ / 10.0)) / log(10)))} @qual;
        }
        
        $data->{-qual} = \@qual;
    }
    
    return $data;
}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq(@seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq::Quality
 Note    : This now conforms to SeqIO spec (module output is same format as
           next_seq)
 Status  : Stable

=cut

# This should be creating fastq output only.  Bio::SeqIO::fasta and 
# Bio::SeqIO::qual should be used for that output

sub write_seq {
    my ($self,@seq) = @_;
    foreach my $seq (@seq) {
		unless ($seq->isa("Bio::Seq::Quality")){
			$self->warn("You can't write FASTQ without supplying a Bio::Seq::Quality object! ", ref($seq), "\n");
			next;
		}
		my $str = $seq->seq;
		my @qual = @{$seq->qual};
		my $top = $seq->display_id();
		if ($seq->can('desc') and my $desc = $seq->desc()) {
			$desc =~ s/\n//g;
			$top .= " $desc";
		}
		if(length($str) == 0) {
			$str = "\n";
		}
		my $qual = join('', map {$self->{phred2chr}->{$_}} @qual);
		$self->_print ("\@",$top,"\n",$str,"\n") or return;
		$self->_print ("+",$top,"\n",$qual,"\n") or return;
    }
    return 1;
}

=head2 write_fastq

 Title   : write_fastq
 Usage   : $stream->write_fastq(@seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq::Quality object
 Status  : Deprecated (delegates to write_seq)

=cut

sub write_fastq {
    my ($self,@seq) = @_;
    return $self->write_seq(@seq);
}

=head2 write_fasta

 Title   : write_fasta
 Usage   : $stream->write_fasta(@seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq object
 Note    : This method does not currently delegate to Bio::SeqIO::fasta
           (maybe it should?).  Not sure whether we should keep this as a
		   convenience method.
 Status  : Unstable

=cut

sub write_fasta {
    my ($self,@seq) = @_;
	if (!exists($self->{fasta_proxy})) {
		$self->{fasta_proxy} = Bio::SeqIO->new(-format => 'fasta', -fh => $self->_fh);
	}
	return $self->{fasta_proxy}->write_seq(@seq);
}

=head2 write_qual

 Title   : write_qual
 Usage   : $stream->write_qual(@seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq::Quality object
 Note    : This method does not currently delegate to Bio::SeqIO::qual
           (maybe it should?).  Not sure whether we should keep this as a
		   convenience method.
 Status  : Unstable
 
=cut

sub write_qual {
	my ($self,@seq) = @_;
	if (!exists($self->{qual_proxy})) {
		$self->{qual_proxy} = Bio::SeqIO->new(-format => 'qual', -fh => $self->_fh);
	}
	
	return $self->{qual_proxy}->write_seq(@seq);
}

=head2 variant

 Title   : variant
 Usage   : $format  = $obj->variant();
 Function: Get and set method for the quality sequence variant.  This is
           important for indicating the encoding/decoding to be used for
           quality data.
           
           Current values accepted are:
            'sanger'   (orginal FASTQ)
                ASCII encoding from 33-126, PHRED quality score from 0 to 93
            'solexa'   (aka illumina1.0)
                ASCII encoding from 59-104, SOLEXA quality score from -5 to 40
            'illumina' (aka illumina1.3)
                ASCII encoding from 64-104, PHRED quality score from 0 to 40
            
            (Derived from the MAQ website):
            For 'solexa', scores are converted to PHRED qual scores using:
                $Q = 10 * log(1 + 10 ** (ord($sq) - 64) / 10.0)) / log(10)
            

 Returns : string
 Args    : new value, string

=cut

{
    my %VARIANT = (
        sanger     => {
            'ord_start'  => 33,
            'ord_end'    => 126,
            'qual_start' => 0,
            'qual_end'   => 93
            },
        solexa     => {
            'ord_start' => 59,
            'ord_end'   => 104,
            'qual_start' => -5,
            'qual_end'   => 40
            },
        illumina   => {
            'ord_start' => 64,
            'ord_end'   => 104,
            'qual_start' => 0,
            'qual_end'   => 40
            },
    );
    
sub variant {
    my ($self, $enc) = @_;
    if (defined $enc) {
        $enc = lc $enc;
        $self->throw('Not a valid FASTQ variant format') unless exists $VARIANT{$enc};
        # cache encode/decode values for quicker accession
		my $ct = 0;
		my ($qs, $qe, $os, $oe) = @{ $VARIANT{$enc} }{qw(qual_start qual_end ord_start ord_end)};
		for my $c ($os..$oe) {
			my $char = chr($c);
			my $score = $qs + $ct++;
			if ($enc eq 'solexa') {
				$score = sprintf("%.0f",(10 * log(1 + 10 ** ($score / 10.0)) / log(10)));
			}
			$self->{chr2phred}->{$char} = $score;
			$self->{phred2chr}->{$score} = $char;
		}
        $self->{qualtype} = $enc;
    }
    return $self->{qualtype};
}

}

sub validate {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{_validate_qual} = $val;
    }
    return $self->{_validate_qual};
}

1;
