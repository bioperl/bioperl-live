# BioPerl module for Bio::SeqIO::fastq
#
# Cared for by Tony Cox <avc@sanger.ac.uk>
#
# Copyright Tony Cox
#
# You may distribute this module under the same terms as perl itself
# _history
# October 29, 2001  incept data

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::fastq - fastq sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::SeqIO class.

=head1 DESCRIPTION

This object can transform Bio::Seq and Bio::Seq::Quality
objects to and from fastq flat file databases.

Fastq is a file format used frequently at the Sanger Centre to bundle
a fasta sequence and its quality data. A typical fastaq entry takes
the from:

  @HCDPQ1D0501
  GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT.....
  +HCDPQ1D0501
  !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65.....

Fastq files have sequence and quality data on a single line and the
quality values are single-byte encoded. To retrieve the decimal values
for qualities you need to subtract 33 (or Octal 41) from each byte and
then convert to a '2 digit + 1 space' integer. You can check if 33 is
the right number because the first byte which is always '!'
corresponds to a quality value of 0.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

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

=cut

sub next_seq {
  my( $self ) = @_;
  my $seq;
  my $alphabet;
  local $/ = "\n";
  my $seqdata;
  my @datatype = qw(seqdesc seq qualdesc qual);
  while (@datatype) {
    return unless my $line = $self->_readline; # bail if any data is incomplete
    chomp $line;
    my $type = shift @datatype;
    if ($type eq 'seqdesc' || $type eq 'qualdesc') {
      $self->throw("$line doesn't match fastq descriptor line type") unless
        $line =~ m{^[\+@](.*)$};
        $line = $1;
    }
    $seqdata->{$type} = $line;
  }
  $self->warn("Seq/Qual descriptions don't match; using sequence description\n")
    unless $seqdata->{seqdesc} eq $seqdata->{qualdesc};
  my ($id,$fulldesc) = $seqdata->{seqdesc} =~ /^\s*(\S+)\s*(.*)/
    or $self->throw("Can't parse fastq header");
  if ($id eq '') {$id=$fulldesc;}   # FIX incase no space between \@ and name
  $seqdata->{seq} =~ s/\s//g;             # Remove whitespace
  $seqdata->{qual} =~ s/\s//g;
  
  if(length($seqdata->{seq}) != length($seqdata->{qual})){
    $self->warn("Fastq sequence/quality data length mismatch error\n",
                "Sequence: ",$seqdata->{seqdesc},", seq length: ",length($seqdata->{seq}), " Qual length: ", length($seqdata->{qual}), " \n",
                $seqdata->{seq},"\n",$seqdata->{qual},"\n");
  }

  my @qual = split('', $seqdata->{qual});

  my $qual;
  foreach (@qual) {$qual .=  (unpack("C",$_) - 33) ." "};
  
  # for empty sequences we need to know the mol.type
  $alphabet = $self->alphabet();
  if(length($seqdata->{seq}) == 0) {
      if(! defined($alphabet)) {
	  # let's default to dna
	  $alphabet = "dna";
      }
  } else {
      # we don't need it really, so disable
      $alphabet = undef;
  }

  # create the Quality object
  $seq = $self->sequence_factory->create(
					 -qual         => $qual,
					 -seq          => $seqdata->{seq},
					 -id           => $id,
					 -primary_id   => $id,
					 -desc         => $fulldesc,
					 -alphabet     => $alphabet
					 );
  
  # if there wasn't one before, set the guessed type
  $self->alphabet($seq->alphabet());
  
  return $seq;
}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq(@seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq::Quality or Bio::Seq object


=cut

sub write_seq {
   my ($self,@seq) = @_;
   foreach my $seq (@seq) {
     my $str = $seq->seq;
     my $top = $seq->display_id();
     if ($seq->can('desc') and my $desc = $seq->desc()) {
	 $desc =~ s/\n//g;
        $top .= " $desc";
     }
     if(length($str) > 0) {
	    $str =~ s/(.{1,60})/$1\n/g;
     } else {
	    $str = "\n";
     }
     
     $self->_print (">",$top,"\n",$str) or return;
   }

   $self->flush if $self->_flush_on_write && defined $self->_fh;
   return 1;
}

=head2 write_qual

 Title   : write_qual
 Usage   : $stream->write_qual(@seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq::Quality object


=cut

sub write_qual {
   my ($self,@seq) = @_;
   foreach my $seq (@seq) {
     unless ($seq->isa("Bio::Seq::Quality")){
        $self->warn("You can write FASTQ without supplying a Bio::Seq::Quality object! ", ref($seq), "\n");
        next;
     } 
     my @qual = @{$seq->qual};
     my $top = $seq->display_id();
     if ($seq->can('desc') and my $desc = $seq->desc()) {
	 $desc =~ s/\n//g;
        $top .= " $desc";
     }
     my $qual = "" ;
     if(scalar(@qual) > 0) {
        my $max = 60;
        for (my $q = 0;$q<scalar(@qual);$q++){
            $qual .= $qual[$q] . " ";
            if(length($qual) > $max){
                $qual .= "\n";
                $max += 60;
            }
        }
     } else {
	    $qual = "\n";
     }
     
     $self->_print (">",$top,"\n",$qual,"\n") or return;
   }
   return 1;
}

=head2 write_fastq

 Title   : write_fastq
 Usage   : $stream->write_fastq(@seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq::Quality object


=cut

sub write_fastq {
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
     my $qual = "" ;
     if(scalar(@qual) > 0) {
        for (my $q = 0;$q<scalar(@qual);$q++){
            $qual .= chr($qual[$q] + 33);
        }
     } else {
	    $qual = "\n";
     }
     
     $self->_print ("\@",$top,"\n",$str,"\n") or return;
     $self->_print ("+",$top,"\n",$qual,"\n") or return;
   }
   return 1;
}
1;
