#----------------------------------------------------------------------------
# PACKAGE : Bio::SeqIO::GCG
# AUTHOR  : Ewan Birney <birney@sanger.ac.uk>, 
#         : Chris Dagdigian <dag@sonsorol.org>
# CREATED : Feb 16 1999
# REVISION: $Id$
#            
# Copyright (c) 1997-9 bioperl, Ewan Birney. All Rights Reserved.
#           This module is free software; you can redistribute it and/or 
#           modify it under the same terms as Perl itself.
#
# _History_
#
#
# Ewan Birney <birney@sanger.ac.uk> developed the SeqIO 
# schema and the first prototype modules.
#
# This code is based on his Bio::SeqIO::Fasta module with
# the necessary alterations needed to get it to read
# and write GCG formatted sequences made by
# chris dagdigian <dag@sonsorol.org>
#
#
#-----------------------------------------------------------------------------
## 
##
## POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::GCG - Raw sequence input/output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the SeqIO handler system. Go:

    $stream = Bio::SeqIO->new(-file => $filename, -format => 'GCG');

    while $seq ( <$stream> ) {
	# $seq is a Bio::Seq object
    }


=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from GCG-formatted
sequence files. 

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

=head1 AUTHORS - Ewan Birney & Chris Dagdigian

Email birney@sanger.ac.uk, dag@sonsorol.org

Ewan Birney wrote the all the base SeqIO stuff. This module is based
on his SeqIO::Fasta and has been tweaked by Chris Dagdigian,
to support GCG formatted sequence input/output.

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqIO::GCG;
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::Seq;
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
      $self->throw("Providing both a file and a filehandle for reading from - only one please!");
  }

  if( $file ) {
      $fh = new FileHandle;
      $fh->open($file) || $self->throw("Could not open $file for GCG stream reading $!");
   }

   $self->_filehandle($fh);
      

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
   my ($self,@args)    = @_;
   my($id,$type,$desc) = undef;
   my ($seq,$fh,$line,$sequence,$chksum);


   $fh = $self->_filehandle();

   if( $fh->eof ) {
       return undef; # no throws - end of file
   }

   while( <$fh> ) {

       ## Get the descriptive info (anything before the line with '..')
       unless( /\.\.$/ ) { $desc.= $_; }

       ## Pull ID, Checksum & Type from the line containing '..'
       /\.\.$/ && do     { $line = $_; chomp; 
                           if(/Check\:\s(\d+)\s/) { $chksum = $1; }
                           if(/Type:\s(\w)\s/)    { $type   = $1; }
                           if(/(.*)\s+Length/)    { $id     = $1; }
                           last; 
                         }
   }


   while( <$fh> ) {

       ## This is where we grab the sequence info.

       if( /\.\.$/ ) { 
        $self->throw("Looks like start of another sequence. See documentation. "); 
       }

       next if($_ eq "\n");       ## skip whitespace lines in formatted seq
       s/[^a-zA-Z]//g;            ## remove anything that is not alphabet char
       $_ = uc($_);               ## uppercase sequence
       $sequence .= $_;
       $fh->eof && last;
   }


   ##If we parsed out a checksum, we might as well test it

   if(defined $chksum) { 
       unless(_validate_checksum($sequence,$chksum)) {
	   $self->throw("Checksum failure on parsed sequence.");
       }
   }

   ## Remove whitespace from identifier because the constructor
   ## will throw a warning otherwise...
   if(defined $id) { $id =~ s/\s+//g;}

   ## Turn our parsed "Type: N" or "Type: P" (if found) into the appropriate
   ## keyword that the constructor expects...
   if(defined $type) {
       if($type eq "N") { $type = "dna";      }
       if($type eq "P") { $type = "amino";    }
   }


   $seq = Bio::Seq->new(-seq  => $sequence, 
                        -id   => $id, 
                        -desc => $desc, 
                        -type => $type );
   return $seq;
}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq($seq)
 Function: writes the fromatted $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq object, sequence string


=cut

sub write_seq {
   my ($self,$seq) = @_;
   my $fh          = $self->_filehandle();
   my $str         = $seq->seq;
   my $comment     = $seq->desc; 
   my $id          = $seq->id;
   my($timestamp)  = localtime;
   my($sum,$offset,$len,$type,$i,$j,$cnt,@out);

   $len = length($str);

   $offset=0;
   ## Set the offset if we have any non-standard numbering going on
   if($seq->start < 0)   { $offset = ( 0 + $self->start); }
   if($seq->start >= 1)  { $offset = $self->start;        }  
   if($seq->start == 0)  { $offset = -1;                  }

   $sum=0;

   #Generate the GCG Checksum value
   for($i=0; $i<$len ;$i++) {             
       $cnt++;
       $sum += $cnt * ord(substr($str,$i,1));
       ($cnt == 57) && ($cnt=0);
   }
   $sum %= 10000;

   #Output the sequence header info
   push(@out,"$comment\n");                        
   push(@out," $id Length: $len  $timestamp  $type Check: $sum  ..\n\n");

   #Format the sequence
   $i = $#out + 1;
   for($j = 0 ; $j < $len ; ) {
       if( $j % 50 == 0) {
	   $out[$i] = sprintf("%8d  ",($j+$offset)); #numbering 
       }
       $out[$i] .= sprintf("%s",substr($str,$j,10));
       $j += 10;
       if( $j < $len && $j % 50 != 0 ) {
	   $out[$i] .= " ";
       }elsif($j % 50 == 0 ) {
	   $out[$i++] .= "\n";
       }                           
   }
   local($^W) = 0;
   if($j % 50 != 0 ) {
       $out[$i] .= "\n";
   }
   $out[$i] .= "\n";


   print $fh join("",@out);
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


=head2 _validate_checksum

 Title   : _validate_checksum
 Usage   : n/a - internal method
 Function: if parsed gcg sequence contains a checksum field
         : we compare it to a value computed here on the parsed
         : sequence. A checksum mismatch would indicate some
         : type of parsing failure occured.
         :
 Returns : 1 for success, 0 for failure
 Args    : string containing parsed seq, value of parsed cheksum


=cut

sub _validate_checksum {
    my($seq,$parsed_sum) = @_;
    my($i,$len,$computed_sum,$cnt);

    $len = length($seq);

    #Generate the GCG Checksum value

    for($i=0; $i<$len ;$i++) {             
	$cnt++;
	$computed_sum += $cnt * ord(substr($seq,$i,1));
	($cnt == 57) && ($cnt=0);
    }
    $computed_sum %= 10000;

    ## Compare and decide if success or failure

    if($parsed_sum = $computed_sum) {
	return 1;
    } else { return 0; }


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
    






