

#
# BioPerl module for Bio::SeqIO::Fasta
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::Fasta - Fasta sequence input/output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the SeqIO handler system. Go:

    $stream = Bio::SeqIO->new(-file => $filename, -format => 'Fasta');

    while $seq ( <$stream> ) {
	# $seq is a Bio::Seq object
    }


=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from fasta flat
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


package Bio::SeqIO::Fasta;
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::Seq;
# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use IO::File;

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
      $self->throw("Neither a file (-file) nor a filehandle (-fh) provided to Fasta opening");
  }


  if( $file ) {

      $fh = new IO::File;
      $fh->open($file) || $self->throw("Could not open $file for Fasta stream reading $!");
  }

#  print "Setting filehandle to $fh\n";
  $self->_filehandle($fh);

# print "Initializing push buffer\n";
  $self->_pushbuffer('');

# set stuff in self from @args
  return $make; # success - we hope!
}

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    :


=cut

sub next_seq{
   my ($self,@args) = @_;
   my ($seq,$fh,$c,$line,$name,$desc,$seqc);

   $fh = $self->_filehandle();

   if( eof $fh ) {
       return undef; # no throws - end of file
   }

   $line = $self->_popbuffer(); # may be '>' character or undef.
   $line .= <$fh>;

   if( $line !~ /^>\s*(\S+)\s*(.*?)\s*$/ ) {
       $self->throw("Fasta stream read attempted with no '>' as first character[ $line ]");
   }
   $name = $1;
   $desc = $2;

   while( <$fh> ) {
       $_ = uc($_);
       s/\W//g;
       $seqc .= $_;
     GETC:
       eof $fh && last;
       $c = getc $fh;
       if( $c eq '>' ) {
	   $self->_pushbuffer($c);
	   last;
       } elsif( $c eq $/ ) { # don't slurp past newline (record separator).
	   goto GETC;
       } else {
	   $seqc .= uc($c) unless $c =~ m/\W/;
       }
   }

   $seq = Bio::Seq->new(-seq => $seqc , -id => $name, -desc => $desc);

   return $seq;

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

#  for ($i = 60; $i < length($str); $i += 60+1) {
#      # this is not ideal.
#      substr($str,$i,0) = "\n";
#  }
   $str = join("\n", grep { length } split(/(.{60})/, $str)); # how's that? -AJM

   print $fh ">", $seq->id(), " ", $seq->desc(), "\n", $str, "\n";
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

sub _pushbuffer {
    my ($obj, $value) = @_;
    if ( defined $value) {
	if ( exists ($obj->{'_pushbuffer'}) ) {
	    $obj->{'_pushbuffer'} .= $value;
	} else {
	$obj->{'_pushbuffer'} = $value;
	}
    }
    return $obj->{'_pushbuffer'};

}

sub _popbuffer {
    my ($obj) = @_;
    return delete ($obj->{'_pushbuffer'}) if exists ($obj->{'_pushbuffer'});
    return undef;
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
    



