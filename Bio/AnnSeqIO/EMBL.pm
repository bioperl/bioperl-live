

#
# BioPerl module for Bio::SeqIO::EMBL
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::AnnSeqIO::EMBL - EMBL sequence input/output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the AnnSeqIO handler system. Go:

    $stream = Bio::AnnSeqIO->new(-file => $filename, -format => 'EMBL');

    while my $annseq ( $stream->next_annseq() ) {
	# do something with $annseq
    }

=head1 DESCRIPTION

This object can transform Bio::AnnSeq objects to and from EMBL flat
file databases.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::AnnSeqIO::EMBL;
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::AnnSeq;
use Bio::AnnSeqIO::FTHelper;
use Bio::SeqFeature::Generic;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use FileHandle;

@ISA = qw(Bio::Root::Object Exporter);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;
   
  my ($file,$fh) = $self->_rearrange([qw(
					 FILE
					 FH
					 )],
				     @args,
				     );
  if( $file && $fh ) {
      $self->throw("Providing both a file and a filehandle for reading from - oly one please!");
  }

  if( !$file && !$fh ) {
      $self->throw("Neither a file (-file) nor a filehandle (-fh) provided to EMBL opening");
  }

  if( $file ) {
      $fh = new FileHandle;
      $fh->open($file) || $self->throw("Could not open $file for EMBL stream reading $!");
  }
  
  # hash for functions for decoding keys.
  $self->{'_func_ftunit_hash'} = {}; 
  $self->_filehandle($fh);
      

# set stuff in self from @args
 return $make; # success - we hope!
}


=head2 next_annseq

 Title   : next_annseq
 Usage   : $seq = $stream->next_annseq()
 Function: returns the next sequence in the stream
 Returns : Bio::AnnSeq object
 Args    :


=cut

sub next_annseq{
   my ($self,@args) = @_;
   my ($seq,$fh,$c,$line,$name,$desc,$seqc);

   $fh = $self->_filehandle();

   if( eof $fh ) {
       return undef; # no throws - end of file
   }

   $line = <$fh>;
   $line =~ /^ID\s+(\S+)/ || $self->throw("EMBL stream with no ID. Not embl in my book");
   $name = $1;

   while( <$fh> ) {
       /^DE\s+(\S.*\S)/ && do { $desc = $1;};
       # we need to do the upper level annotation

       /^FH/ && last;
   }

   while( <$fh> ) {
       /FT   \w/ && last;
   }

   my $buffer = $_;
   
   my $annseq = Bio::AnnSeq->new();
   print "Buffer starts at $buffer\n";
   FEATURE_TABLE :   
   while( !eof($fh) ) {
       my $ftunit = &_read_FTHelper_EMBL($fh,\$buffer);
       # process ftunit
       print "FT unit key ", $ftunit->key, " and loc ",$ftunit->loc,"\n";
       &_generic_seqfeature($annseq,$ftunit);
       if( $buffer !~ /^FT/ ) {
	   last;
       }
	   
   }
	
   if( $buffer !~ /^SQ/  ) {
       while( <$fh> ) {
	   /^SQ/ && last;
       }
   }

	       
   while( <$fh> ) {
       /^\/\// && last;
       $_ = uc($_);
       s/\W//g;
       $seqc .= $_;
       eof $fh && last;
   }

   $seq = Bio::Seq->new(-seq => $seqc , -id => $name, -desc => $desc);
   $annseq->seq($seq);
   return $annseq;

}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq($seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq object


=cut

sub write_seq {
   my ($self,$seq) = @_;
   my $fh = $self->_filehandle();
   my $i;
   my $str = $seq->seq;

   print $fh "ID   ", $seq->id(), "\nDE    ", $seq->desc(), "\nCC   \nCC   Written by Bioperl SeqIO module.\nCC   Only the information in the Sequence object is written in this file\nCC   \nSQ   \n";
   print $fh "    ";
   for ($i = 10; $i < length($str); $i += 10) {
       print $fh substr($str,$i,10), " ";
       if( $i%50 == 0 ) {
	   print $fh "\n    ";
       }
   }
   print $fh "\n//\n";
   return 1;
}

=head2 _filehandle

 Title   : _filehandle
 Usage   : $obj->_filehandle($newval)
 Function: 
 Example : 
 Returns : value of _filehandle
 Args    : newvalue (optional)


=cut

sub _filehandle{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_filehandle'} = $value;
    }
    return $obj->{'_filehandle'};

}

=head2 _read_FTHelper_EMBL

 Title   : _read_FTHelper_EMBL
 Usage   : &_read_FTHelper_EMBL($fh,$buffer)
 Function: reads the next FT key line
 Example :
 Returns : Bio::AnnSeqIO::FTHelper object 
 Args    : filehandle and reference to a scalar


=cut

sub _read_FTHelper_EMBL{
   my ($fh,$buffer) = @_;
   my ($key,$loc,$out);

   if( $$buffer =~ /^FT\s+(\S+)\s+(\S+)/ ) {
       $key = $1;
       $loc = $2;
   }

   if( $loc =~ /\(/ && $loc !~ /complement\(\d+\.\.\d+\)/ ) {
       # more location to read
       while( <$fh> ) {
	   /^FT\s+(\S+)/ || do { &Bio::Root::Object::throw("Weird location line in EMBL feature table"); };
	   $loc .= $1;
	   $loc =~ /\)/ && last;
       }

   }

   # should read in other fields

   while( <$fh> ) {
       /^FT   \w/ && last;
   }

   $$buffer = $_;
   
   $out = new Bio::AnnSeqIO::FTHelper();
   $out->key($key);
   $out->loc($loc);
   return $out;
}

=head2 _generic_seqfeature

 Title   : _generic_seqfeature
 Usage   : &_generic_seqfeature($annseq,$fthelper)
 Function: processes fthelper into a generic seqfeature
 Returns : nothing (places a new seqfeature into annseq)
 Args    : Bio::AnnSeq,Bio::AnnSeqIO::FTHelper


=cut

sub _generic_seqfeature{
   my ($annseq,$fth) = @_;
   my ($sf);

   $sf = new Bio::SeqFeature::Generic;
   if( $fth->loc =~ /join/ ) {
       my $strand;
       if ( $fth->loc =~ /complement/ ) {
	   $strand = -1;
       } else {
	   $strand = 1;
       }
       $sf->strand($strand);
       $sf->primary_tag($fth->key . "_parent");
       $sf->source_tag('EMBL');
       $sf->has_tag("parent",1);

       # we need to make sub features
       my $loc = $fth->loc;
       while( $loc =~ /(\d+)\.\.(\d+)/g ) {
	   my $start = $1;
	   my $end   = $2;
	   print "Processing $start-$end\n";
	   my $sub = new Bio::SeqFeature::Generic;
	   $sub->primary_tag($fth->key);
	   $sub->start($start);
	   $sub->end($end);
	   $sub->strand($strand);
	   $sub->source_tag('EMBL');
	   $sf->add_sub_SeqFeature($sub,'EXPAND');
       }

   } else {
       $fth->loc =~ /(\d+)\.\.(\d+)/ || do {
	   $annseq->throw("Weird location line [" . $fth->loc . "] in reading EMBL");
	   last;
       };
       $sf->start($1);
       $sf->end($2);
       $sf->source_tag('EMBL');
       $sf->primary_tag($fth->key);
       if( $fth->loc =~ /complement/ ) {
	   $sf->strand(-1);
       } else {
	   $sf->strand(1);
       }
   }

   $annseq->add_SeqFeature($sf);
}

sub DESTROY {
    my $self = shift;
    my $fh;
    $fh = $self->_filehandle();

    if( defined $fh ) {
	$fh->close();
    }

    $self->{'_filehandle'} = '';
}
    


