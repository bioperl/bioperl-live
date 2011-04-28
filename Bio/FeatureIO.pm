#
# BioPerl module for Bio::FeatureIO
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Allen Day <allenday@ucla.edu>
#
# Copyright Allen Day
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::FeatureIO - Handler for FeatureIO

=head1 SYNOPSIS

  use Bio::FeatureIO;

  #read from a file
  $in  = Bio::FeatureIO->new(-file => "my.gff" , -format => 'GFF');

  #read from a filehandle
  $in  = Bio::FeatureIO->new(-fh => \*GFF , -format => 'GFF');

  #read features already attached to a sequence
  my $feat = Bio::FeatureIO->new(-seq => $seq , -format => 'features');

  #read new features for existing sequence
  my $seq = Bio::FeatureIO->new(-seq => $seq , -format => 'Das');

  #write out features
  $out = Bio::FeatureIO->new(-file    => ">outputfilename" ,
                             -format  => 'GFF' ,
                             -version => 3);

  while ( my $feature = $in->next_feature() ) {
    $out->write_feature($feature);
  }

=head1 DESCRIPTION

NOTE:  This is a highly experimental I/O iterator subsystem for genomic sequence
features. It is not complete, and is now undergoing a significant refactoring in
a separate branch of BioPerl to address specific issues regarding the current
implementation. Use at your own risk.

Bio::FeatureIO is a handler module for the formats in the FeatureIO set (eg,
Bio::FeatureIO::GFF).

The Bio::FeatureIO system can be thought of like biological file handles.
They are attached to filehandles with smart formatting rules (eg,
GFF format, or BED format) and 
can either read or write feature objects (Bio::SeqFeature objects, or
more correctly, Bio::FeatureHolderI implementing objects, of which 
Bio::SeqFeature is one such object). If you want to know what to 
do with a Bio::SeqFeatureI object, read L<Bio::SeqFeatureI>.

The idea is that you request a stream object for a particular format.
All the stream objects have a notion of an internal file that is read
from or written to. A particular FeatureIO object instance is configured
for either input or output. A specific example of a stream object is
the Bio::FeatureIO::gff object.

Each stream object has functions:

  $stream->next_feature();
  $stream->write_feature($feature);

=head1 SUPPORTED FORMATS

 name                         module
 -----------------------------------
 BED                          bed.pm
 GFF                          gff.pm
 GTF                          gtf.pm
 InterPro (IPRScan 4.0)  interpro.pm
 PTT (NCBI protein table)     ptt.pm


=head1 CONSTRUCTORS

=head2 Bio::FeatureIO-E<gt>new()

   $featureIO = Bio::FeatureIO->new(-file => 'filename',   -format=>$format);
   $featureIO = Bio::FeatureIO->new(-fh   => \*FILEHANDLE, -format=>$format);
   $featureIO = Bio::FeatureIO->new(-seq  => $seq,         -format=>$format);

The new() class method constructs a new Bio::FeatureIO object.  The
returned object can be used to retrieve or print Seq objects. new()
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

   $featio = Bio::FeatureIO->new(-fh => \*STDIN);

Note that you must pass filehandles as references to globs.

If neither a filehandle nor a filename is specified, then the module
will read from the @ARGV array or STDIN, using the familiar E<lt>E<gt>
semantics.

A string filehandle is handy if you want to modify the output in the
memory, before printing it out. The following program reads in EMBL
formatted entries from a file and prints them out in fasta format with
some HTML tags:

  use Bio::FeatureIO;
  use IO::String;
  my $in  = Bio::FeatureIO->new('-file' => "my.gff" ,
  			    '-format' => 'EMBL');
  while ( my $f = $in->next_feature() ) {
      # the output handle is reset for every file
      my $stringio = IO::String->new($string);
      my $out = Bio::FeatureIO->new('-fh' => $stringio,
  			        '-format' => 'gtf');
      # output goes into $string
      $out->write_feature($f);
      # modify $string
      $string =~ s|(>)(\w+)|$1<font color="Red">$2</font>|g;
      # print into STDOUT
      print $string;
  }

=item -format

Specify the format of the file.  See above for list of supported formats

=item -flush

By default, all files (or filehandles) opened for writing sequences
will be flushed after each write_seq() (making the file immediately
usable).  If you don't need this facility and would like to marginally
improve the efficiency of writing multiple sequences to the same file
(or filehandle), pass the -flush option '0' or any other value that
evaluates as defined but false:

  my $f1 = Bio::FeatureIO->new -file   => "<a.f1",
                              -format => "f1";
  my $f2 = Bio::FeatureIO->new -file   => ">a.f2",
                              -format => "f2",
                              -flush  => 0; # go as fast as we can!

  while($feature = $f1->next_feature) { $f2->write_feature($feature) }

=back

=head2 Bio::FeatureIO-E<gt>newFh()

   $fh = Bio::FeatureIO->newFh(-fh   => \*FILEHANDLE, -format=>$format);
   $fh = Bio::FeatureIO->newFh(-format => $format);
   # etc.

This constructor behaves like new(), but returns a tied filehandle
rather than a Bio::FeatureIO object.  You can read sequences from this
object using the familiar E<lt>E<gt> operator, and write to it using
print().  The usual array and $_ semantics work.  For example, you can
read all sequence objects into an array like this:

  @features = <$fh>;

Other operations, such as read(), sysread(), write(), close(), and printf()
are not supported.

=head1 OBJECT METHODS

See below for more detailed summaries.  The main methods are:

=head2 $feature = $featureIO-E<gt>next_feature()

Fetch the next feature from the stream.

=head2 $featureIO-E<gt>write_feature($feature [,$another_feature,...])

Write the specified feature(s) to the stream.

=head2 TIEHANDLE(), READLINE(), PRINT()

These provide the tie interface.  See L<perltie> for more details.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
to one of the Bioperl mailing lists.

Your participation is much appreciated.

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

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Allen Day

Email allenday@ucla.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#' Let the code begin...

package Bio::FeatureIO;

use strict;

use Symbol;

use base qw(Bio::Root::Root Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : $stream = Bio::FeatureIO->new(-file => $filename, -format => 'Format')
 Function: Returns a new feature stream
 Returns : A Bio::FeatureIO stream initialised with the appropriate format
 Args    : Named parameters:
             -file => $filename
             -fh => filehandle to attach to
             -format => format

=cut

my $entry = 0;

sub new {
  my ($caller,@args) = @_;
  my $class = ref($caller) || $caller;

  # or do we want to call SUPER on an object if $caller is an
  # object?
  if( $class =~ /Bio::FeatureIO::(\S+)/ ) {

    my ($self) = $class->SUPER::new(@args);	
    $self->_initialize(@args);
    return $self;

  } else {

	my %param = @args;

	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	my $format = $param{'-format'} ||
      $class->_guess_format( $param{-file} || $ARGV[0] );
	
	if( ! $format ) {
      if ($param{-file}) {
        $format = $class->_guess_format($param{-file});
      } elsif ($param{-fh}) {
        $format = $class->_guess_format(undef);
      }
	}
	$format = "\L$format";	# normalize capitalization to lower case
	return unless( $class->_load_format_module($format) );
	return "Bio::FeatureIO::$format"->new(@args);

  }
}

=head2 newFh

 Title   : newFh
 Usage   : $fh = Bio::FeatureIO->newFh(-file=>$filename,-format=>'Format')
 Function: does a new() followed by an fh()
 Example : $fh = Bio::FeatureIO->newFh(-file=>$filename,-format=>'Format')
           $feature = <$fh>;   # read a feature object
           print $fh $feature; # write a feature object
 Returns : filehandle tied to the Bio::FeatureIO::Fh class
 Args    :

See L<Bio::FeatureIO::Fh>

=cut

sub newFh {
  my $class = shift;
  return unless my $self = $class->new(@_);
  return $self->fh;
}

=head2 fh

 Title   : fh
 Usage   : $obj->fh
 Function:
 Example : $fh = $obj->fh;      # make a tied filehandle
           $feature = <$fh>;   # read a feature object
           print $fh $feature; # write a feature object
 Returns : filehandle tied to Bio::FeatureIO class
 Args    : none

=cut


sub fh {
  my $self = shift;
  my $class = ref($self) || $self;
  my $s = Symbol::gensym;
  tie $$s,$class,$self;
  return $s;
}

# _initialize is chained for all FeatureIO classes

sub _initialize {
    my($self, %arg) = @_;

    # flush is initialized by the Root::IO init

    # initialize the IO part
    $self->seq($arg{-seq});
    $self->_initialize_io(%arg);
}

=head2 next_feature

 Title   : next_feature
 Usage   : $feature = stream->next_feature
 Function: Reads the next feature object from the stream and returns it.

           Certain driver modules may encounter entries in the stream
           that are either misformatted or that use syntax not yet
           understood by the driver. If such an incident is
           recoverable, e.g., by dismissing a feature of a feature
           table or some other non-mandatory part of an entry, the
           driver will issue a warning. In the case of a
           non-recoverable situation an exception will be thrown.  Do
           not assume that you can resume parsing the same stream
           after catching the exception. Note that you can always turn
           recoverable errors into exceptions by calling
           $stream->verbose(2).

 Returns : a Bio::SeqFeatureI feature object
 Args    : none

See L<Bio::Root::RootI>, L<Bio::SeqFeatureI>

=cut

sub next_feature {
   my ($self, $seq) = @_;
   $self->throw("Sorry, you cannot read from a generic Bio::FeatureIO object.");
}

=head2 write_feature

 Title   : write_feature
 Usage   : $stream->write_feature($feature)
 Function: writes the $feature object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::SeqFeature object

=cut

sub write_feature {
    my ($self, $seq) = @_;
    if(ref($self) eq __PACKAGE__){
      $self->throw("Sorry, you cannot write to a generic Bio::FeatureIO object.");
    } else {
      $self->throw_not_implemented();
    }
}

=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL FeatureIO stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns :
 Args    :

=cut

sub _load_format_module {
    my ($self, $format) = @_;
    my $class = ref($self) || $self;
    my $module = $class."::$format";#"Bio::Feature::" . $format;
    my $ok;

    eval {
	$ok = $self->_load_module($module);
    };
    if ( $@ ) {
    print STDERR <<END;
$self: $format cannot be found
Exception $@
For more information about the FeatureIO system please see the FeatureIO docs.
This includes ways of checking for formats at compile time, not run time
END
  ;
  }
  return $ok;
}

=head2 seq

 Title   : seq
 Usage   : $obj->seq() OR $obj->seq($newSeq)
 Example :
 Returns : Bio::SeqI object
 Args    : newSeq (optional)

=cut

sub seq {
  my $self = shift;
  my $val = shift;

  $self->{'seq'} = $val if defined($val);
  return $self->{'seq'};
}

=head2 _filehandle

 Title   : _filehandle
 Usage   : $obj->_filehandle($newval)
 Function: This method is deprecated. Call _fh() instead.
 Example :
 Returns : value of _filehandle
 Args    : newvalue (optional)


=cut

sub _filehandle {
    my ($self,@args) = @_;
    return $self->_fh(@args);
}

=head2 _guess_format

 Title   : _guess_format
 Usage   : $obj->_guess_format($filename)
 Function: guess format based on file suffix
 Example :
 Returns : guessed format of filename (lower case)
 Args    :
 Notes   : See "SUPPORTED FORMATS"

=cut

sub _guess_format {
   my $class = shift;
   return unless $_ = shift;
   return 'gff'     if /\.gff3?$/i;
   return 'gff'     if /\.gtf$/i;
   return 'bed'     if /\.bed$/i;
   return 'ptt'     if /\.ptt$/i;

   return 'gff'; #the default
}

sub DESTROY {
    my $self = shift;
    $self->close();
}

sub TIEHANDLE {
    my ($class,$val) = @_;
    return bless {'featio' => $val}, $class;
}

sub READLINE {
  my $self = shift;
  return $self->{'featio'}->next_feature() unless wantarray;
  my (@list, $obj);
  push @list, $obj while $obj = $self->{'featio'}->next_feature();
  return @list;
}

sub PRINT {
  my $self = shift;
  $self->{'featio'}->write_feature(@_);
}

1;

