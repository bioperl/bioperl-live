#
# BioPerl module for Bio::Variation::IO
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Heikki Lehvaslaiho <heikki-at-bioperl-dot-org>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Variation::IO - Handler for sequence variation IO Formats

=head1 SYNOPSIS

    use Bio::Variation::IO;

    $in  = Bio::Variation::IO->new(-file => "inputfilename" , 
                                   -format => 'flat');
    $out = Bio::Variation::IO->new(-file => ">outputfilename" ,
                                   -format => 'xml');

    while ( my $seq = $in->next() ) {
	   $out->write($seq);
    }

  # or

    use Bio::Variation::IO;

    #input file format can be read from the file extension (dat|xml)
    $in  = Bio::Variation::IO->newFh(-file => "inputfilename");
    $out = Bio::Variation::IO->newFh(-format => 'xml');

    # World's shortest flat<->xml format converter:
    print $out $_ while <$in>;

=head1 DESCRIPTION

Bio::Variation::IO is a handler module for the formats in the 
Variation IO set (eg, Bio::Variation::IO::flat). It is the officially 
sanctioned way of getting at the format objects, which most people 
should use.

The structure, conventions and most of the code is inherited from
L<Bio::SeqIO> module. The main difference is that instead of using
methods next_seq and write_seq, you drop '_seq' from the method names.

The idea is that you request a stream object for a particular format.
All the stream objects have a notion of an internal file that is read
from or written to. A particular SeqIO object instance is configured
for either input or output. A specific example of a stream object is
the Bio::Variation::IO::flat object.

Each stream object has functions

   $stream->next();

and

   $stream->write($seqDiff);

also

   $stream->type() # returns 'INPUT' or 'OUTPUT'

As an added bonus, you can recover a filehandle that is tied to the
SeqIO object, allowing you to use the standard E<lt>E<gt> and print 
operations to read and write sequence objects:

    use Bio::Variation::IO;

    $stream = Bio::Variation::IO->newFh(-format => 'flat'); 
    # read from standard input

    while ( $seq = <$stream> ) {
	   # do something with $seq
    }

and

    print $stream $seq; # when stream is in output mode

This makes the simplest ever reformatter

    #!/usr/local/bin/perl

    $format1 = shift;
    $format2 = shift;

    use Bio::Variation::IO;

    $in  = Bio::Variation::IO->newFh(-format => $format1 );
    $out = Bio::Variation::IO->newFh(-format => $format2 );

    print $out $_ while <$in>;


=head1 CONSTRUCTORS

=head2 Bio::Variation::IO-E<gt>new()

   $seqIO = Bio::Variation::IO->new(-file => 'filename',   -format=>$format);
   $seqIO = Bio::Variation::IO->new(-fh   => \*FILEHANDLE, -format=>$format);
   $seqIO = Bio::Variation::IO->new(-format => $format);

The new() class method constructs a new Bio::Variation::IO object.  The
returned object can be used to retrieve or print BioSeq objects. new()
accepts the following parameters:

=over 4

=item -file

A file path to be opened for reading or writing.  The usual Perl
conventions apply:

   'file'       # open file for reading
   '>file'      # open file for writing
   '>>file'     # open file for appending
   '+<file'     # open file read/write
   'command |'  # open a pipe from the command
   '| command'  # open a pipe to the command

=item -fh

You may provide new() with a previously-opened filehandle.  For
example, to read from STDIN:

   $seqIO = Bio::Variation::IO->new(-fh => \*STDIN);

Note that you must pass filehandles as references to globs.

If neither a filehandle nor a filename is specified, then the module
will read from the @ARGV array or STDIN, using the familiar E<lt>E<gt>
semantics.

=item -format

Specify the format of the file.  Supported formats include:

   flat        pseudo EMBL format
   xml         seqvar xml format

If no format is specified and a filename is given, then the module
will attempt to deduce it from the filename.  If this is unsuccessful,
Fasta format is assumed.

The format name is case insensitive.  'FLAT', 'Flat' and 'flat' are
all supported.

=back

=head2 Bio::Variation::IO-E<gt>newFh()

   $fh = Bio::Variation::IO->newFh(-fh   => \*FILEHANDLE, -format=>$format);
   $fh = Bio::Variation::IO->newFh(-format => $format);
   # etc.

   #e.g.
   $out = Bio::Variation::IO->newFh( '-FORMAT' => 'flat');
   print $out $seqDiff;

This constructor behaves like new(), but returns a tied filehandle
rather than a Bio::Variation::IO object.  You can read sequences from this
object using the familiar E<lt>E<gt> operator, and write to it using print().
The usual array and $_ semantics work.  For example, you can read all
sequence objects into an array like this:

  @mutations = <$fh>;

Other operations, such as read(), sysread(), write(), close(), and printf() 
are not supported.

=head1 OBJECT METHODS

See below for more detailed summaries.  The main methods are:

=head2 $sequence = $seqIO-E<gt>next()

Fetch the next sequence from the stream.

=head2 $seqIO-E<gt>write($sequence [,$another_sequence,...])

Write the specified sequence(s) to the stream.

=head2 TIEHANDLE(), READLINE(), PRINT()

These provide the tie interface.  See L<perltie> for more details.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the 
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

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

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Variation::IO;

use strict;


use base qw(Bio::SeqIO Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : $stream = Bio::Variation::IO->new(-file => $filename, -format => 'Format')
 Function: Returns a new seqstream
 Returns : A Bio::Variation::IO::Handler initialised with the appropriate format
 Args    : -file => $filename
           -format => format
           -fh => filehandle to attach to

=cut


sub new {
   my ($class, %param) = @_;
   my ($format);

   @param{ map { lc $_ } keys %param } = values %param;  # lowercase keys
   $format = $param{'-format'}
             || $class->_guess_format( $param{-file} || $ARGV[0] )
             || 'flat';
   $format = "\L$format"; # normalize capitalization to lower case

   return unless $class->_load_format_module($format);
   return "Bio::Variation::IO::$format"->new(%param);
}


=head2 format

 Title   : format
 Usage   : $format = $stream->format()
 Function: Get the variation format
 Returns : variation format
 Args    : none

=cut

# format() method inherited from Bio::Root::IO


sub _load_format_module {
  my ($class, $format) = @_;
  my $module = "Bio::Variation::IO::" . $format;
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

=head2 next

 Title   : next
 Usage   : $seqDiff = $stream->next
 Function: reads the next $seqDiff object from the stream
 Returns : a Bio::Variation::SeqDiff object
 Args    :

=cut

sub next {
   my ($self, $seq) = @_;
   $self->throw("Sorry, you cannot read from a generic Bio::Variation::IO object.");
}

sub next_seq {
   my ($self, $seq) = @_;
   $self->throw("These are not sequence objects. Use method 'next' instead of 'next_seq'.");
   $self->next($seq);
}

=head2 write

 Title   : write
 Usage   : $stream->write($seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Variation::SeqDiff object

=cut

sub write {
    my ($self, $seq) = @_;
    $self->throw("Sorry, you cannot write to a generic Bio::Variation::IO object.");
}

sub write_seq {
   my ($self, $seq) = @_;
   $self->warn("These are not sequence objects. Use method 'write' instead of 'write_seq'.");
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
   return unless $_ = shift;
   return 'flat'     if /\.dat$/i;
   return 'xml'     if /\.xml$/i;
}


1;
