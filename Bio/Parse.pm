package Bio::Parse;
require Exporter;
@ISA         = qw(Exporter);
@EXPORT      = qw();
@EXPORT_OK   = qw($VERSION $OK);
$VERSION     = 0.05;

# $Id$ #

############  MANUAL EDITING NEEDED HERE!
############

## Uncomment the $OK stuff (it is mentioned twice to
## shut strict up. Set $READSEQ_PATH to the full path
## to the readseq executable.
 
## $OK = "Y"; 
## $OK = "Y"; 
my $READSEQ_PATH=undef;

############
############


# Copyright (c) 1997-1998
# Chris Dagdigian and others. All Rights Reserved.
# This module is free software; you can redistribute it and/or modify 
# it under the same terms as Perl itself.

require 5.003;
use Carp;
use POSIX;


## POD-formatted documentation

=head1 NAME

Seq::Parse - The Bioperl ReadSeq interface

=head1 SYNOPSIS

Simple perl interface/wrapper to D.G. Gilbert's ReadSeq program.
Used by Seq.pm when internal parsing/formatting code fails.

**NOTE** Not currently used by any of the core bioperl modules.
It can be used as a standalone interface to the readseq package but
manual editing of is required. See the first few lines of the .pm
file for details.

=head1 DESCRIPTION

This package was called upon by Seq.pm when internal attemts to
format or parse a sequence fail. It is currently not used by
any bioperl core module. Basically we decided to deal with
sequence formatting in a different way.

Parse.pm is a simple interface to D.G. Gilbert's ReadSeq program,
it is not meant to be particularly elegant or efficient. The
interface should be abstract enough to allow future versions
to seamlessly access other sequence conversion programs besides
ReadSeq.

At this time the interface methods have not been fully thought
out or implemented. Suggestions are welcome.

If ReadSeq is not on the local system, or this package is not
properly configured, Seq.pm will (hopefully) realize this
and not attempt to use this code.

=head1 USAGE

The ReadSeq executable needs to be installed on your system.

Readseq is freely distributed and is available in
shell archive (.shar) form via FTP from
ftp.bio.indiana.edu (129.79.224.25) in the
molbio/readseq directory.
(URL) ftp://ftp.bio.indiana.edu/molbio/readseq/

=head2 Standalone

 use Parse;

=head2 With Seq.pm

If properly configured, Seq.pm will automatically use this module when
internal methods at parsing or formatting fail. 

The correct path to the readseq executable is configured into
this module during the 'make Makefile.PL' phase of installation.

Manual edits needed in Parse.pm if auto-configuration does not happen:

- Change the value of B<$READSEQ_PATH> so that it defines a path to the
  ReadSeq executable on your system.

- Uncomment the line(s) containing $OK = "Y"


=head2 As a standalone module

Parse.pm should be usable is a standalone module. See the usage instructions.

=head2 Sequence Conversion/Formatting

ReadSeq has trouble with raw sequences so an explicit
convert_from_raw() method has been written.
The following code will return the sequence "GAATTCGTT" as a GCG formatted string.

 $reply  = &Parse::convert_from_raw(-sequence=>'GAATTCGTT',
                                    -fmt=>'gcg'); 


The "fmt" named-parameter field can be set for the following formats:

 IG        (or 'Stanford')
 GenBank   (or 'GB')
 NBRF
 EMBL
 GCG
 Strider
 Fitch
 Fasta
 Zuker
 Phylip3.2 (use 'Phylip3')
 Phylip
 Plain     (or 'Raw')
 PIR       (or 'CODATA')
 MSF
 ASN.1     (use 'ASN1')
 PAUP
 Pretty

The "options" named-parameter field can be used to pass switches directly
to the ReadSeq executable. This option should only be used by people
familiar with operating ReadSeq on the command-line. Use at your own
risk as this has not been fully tested. 

As an example, the ReadSeq switch C<-c> will cause all of the
characters in the formatted sequence to be returned in lowercase.

 $reply  = &Parse::convert_from_raw(-sequence=>"$seq_string",
                                    -options=>'-c', 
                                    -fmt=>'gcg'); 

=head1 Appendix

The following documentation describes the various functions
contained in this package. Some functions
are for internal use and are not meant to be called by
the user; they are preceded by an underscore ("_").
 

=cut

#
##
###
#### END of main POD documentation. Let the code begin..
###
##
#


## Package Global
## Map sequence formats to ReadSeq 
## command line switches
my(%Format);
$Format{ig}        = 1;
$Format{stanford}  = 1;
$Format{genbank}   = 2;
$Format{gb}        = 2;
$Format{nbrf}      = 3;
$Format{embl}      = 4;
$Format{gcg}       = 5;
$Format{strider}   = 6;
$Format{fitch}     = 7;
$Format{fasta}     = 8;
$Format{zuker}     = 9;
$Format{olsen}     = 10;
$Format{phylip2}   = 11;
$Format{phylip}    = 12;
$Format{plain}     = 13;
$Format{raw}       = 13;
$Format{pir}       = 14;
$Format{codata}    = 14;
$Format{msf}       = 15;
$Format{asn1}      = 16;
$Format{paup}      = 17;




#############################
=head2 ## Internal methods ##
=cut

#_______________________________________________________________________

=head2 _rearrange()

 
 Title     : _rearrange
 Usage     : n/a (internal function)
 Function  : Rearranges named parameters to requested order.
 Example   : &_rearrange([SEQUENCE,ID,DESC],@p);
 Returns   : @params - an array of parameters in the requested order.
 Argument  : $order : a reference to an array which describes the desired
                      order of the named parameters.
             @param : an array of parameters, either as a list (in
                      which case the function simply returns the list),
                      or as an associative array (in which case the
                      function sorts the values according to @{$order}
                      and returns that new array.

=cut

sub _rearrange {
  # This function was taken from CGI.pm, written by Dr. Lincoln
  # Stein, and adapted for use in Bio::* by Richard Resnick.
  my($order,@param) = @_;

  # If there are no parameters, we simply wish to return
  # an empty array which is the size of the @{$order} array.
  return ('') x $#{$order} unless @param;

  # If we've got parameters, we need to check to see whether
  # they are named or simply listed. If they are listed, we
  # can just return them.
  return @param unless (defined($param[0]) && $param[0]=~/^-/);

  # Now We've got to do some work on the named parameters.
  # The next few lines strip out the '-' characters which
  # preceed the keys, and capitalizes them.
  my $i;
  for ($i=0;$i<@param;$i+=2) {
        $param[$i]=~s/^\-//;
        $param[$i]=~tr/a-z/A-Z/;
  }
  
  # Now we'll convert the @params variable into an associative array.
  my(%param) = @param;

  my(@return_array);
  
  # What we intend to do is loop through the @{$order} variable,
  # and for each value, we use that as a key into our associative
  # array, pushing the value at that key onto our return array.
  my($key);

  foreach $key (@{$order}) {
        my($value) = $param{$key};
        delete $param{$key};
        push(@return_array,$value);
  }

  my(@restkeys) = keys %param;
  if (scalar(@restkeys) > 0) {
       carp("@restkeys not processed in _rearrange(), did you use a
       nonexisting parameter name ? ");
  }

  return (@return_array);
}




#_______________________________________________________________________

=head2 _write_tmp_file()

 
 Title     : _write_tmp_file
 Usage     : n/a (internal function)
 Function  : Writes a temporary file to disk. Uses
           : the POSIX tmpnam() call to get path &
           : filename. Should be more portable than
           : just writing to /tmp. 
           :
 Example   : &_write_tmp_file("$formatted_sequence");
 Returns   : string containing the temp file path 
 Argument  : string that is to be written to disk
 

=cut

sub _write_tmp_file {
my($contents) = @_;
my($fh,$tmp_path);

     #Get temporary filename to write sequence
     $tmp_path = POSIX::tmpnam();

     #Write sequence to file
     open(FILE_OUT,">$tmp_path") || carp "unable to open filehandle: $!\n";
       print FILE_OUT $contents;
     close(FILE_OUT);
 
$tmp_path;   #Return the POSIX-derived temp pathname 
}


 
#################
## Main Methods

#_______________________________________________________________________

=head2 version()

 
 Title     : version
 Usage     : &Parse::version;
 Function  : Prints current package version 
 Example   : &Parse::version;
 Returns   : none
 Argument  : none
           :

=cut

sub version {
  print "Bio::Parse Version is ", $Bio::Parse::VERSION, ".\n";
  return 1;
}



#_______________________________________________________________________

=head2 convert_from_raw()

 
 Title     : convert_from_raw()
 Usage     : print &Parse::convert_from_raw(-sequence=>$raw_seq,
           :                                -fmt=>'asn1');
           :
           : $reply  = &Parse::convert_from_raw(-sequence=>'GAATTCGTT',
           :                                    -options=>'-c',
           :                                    -fmt=>'gcg'); 
           :
 Function  : ReadSeq does not function well when called upon 
           : to read or convert "raw" or unformatted sequence 
           : strings or files. This code will take a given 
           : raw sequence and manipulate it into FASTA
           : format before invoking ReadSeq.
           :
           : The following named paramters may be used as
           : arguments:
           :
           :  -sequence=>     Sequence string.
           :  -fmt=>          Format sequence will be converted to. 
           :  -options=>      String containing command-line
           :                  switches for ReadSeq. Passed
           :                  directly.
           :
 Example   : see usage
 Returns   : Formatted sequence string 
 Argument  : named parameters, see function
           :

=cut

sub convert_from_raw {
my(@params) = @_;
my($form,$filename,$reply) = "";
my($seq,$location,$fmt,$opts,$temp_seqpath) = "";
my(@response) = "";

 ## Get named parameters
 if(defined(@params)) {
   ($seq,$fmt,$opts) = 
             &_rearrange(['SEQUENCE',
                         'FMT',
                         'OPTIONS',],
                        @params);
  }


 ##Write user submitted seq to a temp file
 if(defined($seq)) {
    $temp_seqpath  =  &_write_tmp_file(">\n$seq\n");
    $filename = $temp_seqpath;
 }

 unless(defined($opts)) { $opts="";}

 ##Error check user submitted format choice
 $fmt =~ tr/A-Z/a-z/;
  unless($Format{$fmt}) { carp "Error: unknown or improper format selected.\n"; }
  else {
        $form = $Format{$fmt};
       }

       open(READSEQ,"$READSEQ_PATH -pipe -format=$form $opts $filename |") || 
         carp "Unable to sucessfuly open pipe to ReadSeq";
         @response = <READSEQ>;
       close(READSEQ);

    ## Delete any temporary files...
    if(defined($temp_seqpath)) { unlink($temp_seqpath); }

$reply = join('',@response);
}

#_______________________________________________________________________

=head2 convert()

 
 Title     : convert
           :
 Usage     : print &Parse::convert(-sequence=>$raw_seq,
           :                       -fmt=>'asn1');
           :
           : $reply  = &Parse::convert(-sequence=>'GAATTCGTT',
           :                           -options=>'-c',
           :                           -fmt=>'gcg'); 
           :
           : $reply  = &Parse::convert(-location=>'/tmp/a.seq',
           :                           -fmt=>'raw'); 
           :
 Note      : ReadSeq does not function well when called upon 
           : to read or convert "raw" or unformatted sequence 
           : strings or files. User beware.
           : 
 Function  : Will read/parse a given sequence string *OR* a given
           : sequence file.
           :
           : If a sequence string AND a sequence file path are
           : both passed in, the file path will be used with no
           : complaint.
           :
           : The following named paramters may be used as
           : arguments:
           : 
           :  -sequence=>     Sequence string.
           :  -location=>     Sequence file path.
           :  -fmt=>          Format sequence will be converted to. 
           :  -options=>      String containing command-line
           :                  switches for ReadSeq. Passed
           :                  directly.
           :
 Example   : see usage
 Returns   : Formatted sequence string 
 Argument  : named parameters, see function
           :

=cut

sub convert {
my(@params) = @_;
my($form,$filename,$reply) = "";
my($seq,$location,$fmt,$opts,$temp_seqpath) = "";
my(@response) = "";

 ## Get named parameters
 if(defined(@params)) {
     ($seq,$location,$fmt,$opts) = 
             &_rearrange(['SEQUENCE',
                         'LOCATION',
                         'FMT',
                         'OPTIONS',],
                        @params);
 }

 ##Write user submitted seq to a temp file
 if(defined($seq)) {
    $temp_seqpath  =  &_write_tmp_file($seq);
    $filename = $temp_seqpath;
 }

 ##Deal with a user submitted sequence file path
 ## (takes precedence over submitted seq string)
 if(defined($location)) {
     unless(-r $location) { carp "Error: Sequence file $location is not readable\n"; }
     $filename = $location;
 }

 unless(defined($opts)) { $opts="";}

 ##Error check user submitted format choice
 $fmt =~ tr/A-Z/a-z/;
  unless($Format{$fmt}) { carp "Error: unknown or improper format selected.\n"; }
  else {
        $form = $Format{$fmt};
       }

       open(READSEQ,"$READSEQ_PATH -pipe -format=$form $opts $filename |") || 
         carp "Unable to sucessfuly open pipe to ReadSeq";
         @response = <READSEQ>;
       close(READSEQ);

    ## Delete any temporary files...
    if(defined($temp_seqpath)) { unlink($temp_seqpath); }

$reply = join('',@response);
}



#-----------------------------------

=head1 ACKNOWLEDGEMENTS


=head1 SEE ALSO

 Core bioperl modules

=head1 REFERENCES

Bioperl Project
http://bio.perl.org

=head1 COPYRIGHT

Copyright (c) 1997-1998 Chris Dagdigian, Georg Fuellen,
Steven E. Brenner and others. All Rights Reserved.
This module is free software; you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut



