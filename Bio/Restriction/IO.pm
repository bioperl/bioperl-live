# $Id$
#
# BioPerl module for Bio::Restriction::IO
#
# Cared for by Rob Edwards <redwards@utmem.edu>
#
# Copyright Rob Edwards
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Restriction::IO - Handler for sequence variation IO Formats

=head1 SYNOPSIS

    use Bio::Restriction::IO;

    $in  = Bio::Restriction::IO->new(-file => "inputfilename" ,
                                     -format => 'withrefm');
    $out = Bio::Restriction::IO->new(-file => ">outputfilename" ,
                                     -format => 'bairoch');
    my $res = $in->read; # a Bio::Restriction::EnzymeCollection
    $out->write($res);

  # or

  #    use Bio::Restriction::IO;
  #
  #    #input file format can be read from the file extension (dat|xml)
  #    $in  = Bio::Restriction::IO->newFh(-file => "inputfilename");
  #    $out = Bio::Restriction::IO->newFh('-format' => 'xml');
  #
  #    # World's shortest flat<->xml format converter:
  #    print $out $_ while <$in>;

=head1 DESCRIPTION

Bio::Restriction::IO is a handler module for the formats in the
Restriction IO set (eg, Bio::Restriction::IO::XXX). It is the
officially sanctioned way of getting at the format objects, which most
people should use.

The structure, conventions and most of the code is inherited from
L<Bio::SeqIO> module. The main difference is that instead of using
methods next_seq and write_seq, you drop '_seq' from the method names.

Also, instead of dealing only with individual Bio::Restriction::Enzyme
objects, write() will go through all enzymes from a
Bio::Restriction::EnzymeCollection objects and read() will slurp in all
enzymes into a Collection.

For more details, see documentation in Bio::SeqIO.

=head1 TO DO

This inherits from SeqIOfor convinience sake. Get rid of it by copying
relevant methods in. Bio::SeqIO has too many sequence IO specific
tweeks in it.

=head1 SEE ALSO

L<Bio::SeqIO>, 
L<Bio::Restriction::Enzyme>, 
L<Bio::Restriction::EnzymeCollection>

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists Your participation is much appreciated.

  bioperl-l@bioperl.org                     - General discussion
  http://bio.perl.org/MailList.html         - About the mailing lists

=head2 Reporting Bugs

report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR

Rob Edwards, redwards@utmem.edu

=head1 CONTRIBUTORS

Heikki Lehvaslaiho, heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Restriction::IO;

use strict;
use vars qw(@ISA %FORMAT);
use Bio::SeqIO;
@ISA = 'Bio::SeqIO';

%FORMAT = (
            'itype2'    => 'itype2',
            '8'         => 'itype2',
            'withrefm'  => 'withrefm',
            '31'        => 'withrefm',
            'base'      => 'base',
            '0'         => 'base',
	    'bairoch'   => 'bairoch',
	    '19'        => 'bairoch',
	    'macvector' => 'bairoch',
	    'vectorNTI' => 'bairoch'
);

=head2 new

 Title   : new
 Usage   : $stream = Bio::Restriction::IO->new(-file => $filename,
                                               -format => 'Format')
 Function: Returns a new seqstream
 Returns : A Bio::Restriction::IO::Handler initialised with
           the appropriate format
 Args    : -file => $filename
           -format => format
           -fh => filehandle to attach to

=cut

sub new {
   my ($class, %param) = @_;
   my ($format);

   @param{ map { lc $_ } keys %param } = values %param;  # lowercase keys

   $format = $FORMAT{$param{'-format'}} if defined $param{'-format'};
   $format ||= $class->_guess_format( $param{-file} || $ARGV[0] )
             || 'base';
   $format = "\L$format"; # normalize capitalization to lower case

   return undef unless $class->_load_format_module($format);
   return "Bio::Restriction::IO::$format"->new(%param);
}


sub _load_format_module {
  my ($class, $format) = @_;
  my $module = "Bio::Restriction::IO::" . $format;
  my $ok;
  eval {
      $ok = $class->_load_module($module);
  };
  if ( $@ ) {
    print STDERR <<END;
$class: $format cannot be found
Exception $@
For more information about the IO system please see the IO docs.
This includes ways of checking for formats at compile time, not run time
END
  ;
  }
  return $ok;
}

=head2 read

 Title   : read
 Usage   : $renzs = $stream->read
 Function: reads all the restrction enzymes from the stream
 Returns : a Bio::Restriction::EnzymeCollection object
 Args    :

=cut

sub read {
   my ($self, $seq) = @_;
   $self->throw("Not implemented");
}

sub next {
   my ($self, $seq) = @_;
   $self->throw("Not implemented");
}

sub next_seq {
   my ($self, $seq) = @_;
   $self->throw("Not implemented");
}

=head2 write

 Title   : write
 Usage   : $stream->write($seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Restriction::EnzymeCOllection object

=cut

sub write {
    my ($self, $seq) = @_;
    $self->throw("Sorry, you cannot write to a generic ".
                 "Bio::Restricion::IO object.");
}

sub write_seq {
   my ($self, $seq) = @_;
   $self->warn("These are not sequence objects. ".
               "Use method 'write' instead of 'write_seq'.");
   $self->write($seq);
}

=head2 _guess_format

 Title   : _guess_format
 Usage   : $obj->_guess_format($filename)
 Function:
 Example :
 Returns : guessed format of filename (lower case)
 Args    :

=cut

sub _guess_format {
   my $class = shift;
   return  unless $_ = shift;
   return 'flat'     if /\.dat$/i;
   return 'xml'     if /\.xml$/i;
}


1;
