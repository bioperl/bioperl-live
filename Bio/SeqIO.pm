
# $Id$
#
# BioPerl module for Bio::SeqIO
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#       and Lincoln Stein  <lstein@cshl.org>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself
#
# _history
# October 18, 1999  Largely rewritten by Lincoln Stein

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO - Handler for SeqIO Formats

=head1 SYNOPSIS

    use Bio::SeqIO;

    $in  = Bio::SeqIO->new(-file => "inputfilename" , '-format' => 'Fasta');
    $out = Bio::SeqIO->new(-file => ">outputfilename" , '-format' => 'EMBL');
    # note: we quote -format to keep older Perls from complaining.

    while ( my $seq = $in->next_seq() ) {
	$out->write_seq($seq);
    }

Now, to actually get at the sequence object, use the standard Bio::Seq
methods (look at L<Bio::Seq> if you don't know what they are)

    use Bio::SeqIO;

    $in  = Bio::SeqIO->new(-file => "inputfilename" , '-format' => 'genbank');

    while ( my $seq = $in->next_seq() ) {
       print "Sequence ",$seq->id," first 10 bases ",$seq->subseq(1,10),"\n";
    }


The SeqIO system does have a filehandle binding. Most people find this
a little confusing, but it does mean you write the world's smallest
reformatter

    use Bio::SeqIO;

    $in  = Bio::SeqIO->newFh(-file => "inputfilename" , '-format' => 'Fasta');
    $out = Bio::SeqIO->newFh('-format' => 'EMBL');

    # World's shortest Fasta<->EMBL format converter:
    print $out $_ while <$in>;


=head1 DESCRIPTION

Bio::SeqIO is a handler module for the formats in the SeqIO set (eg,
Bio::SeqIO::fasta). It is the officially sanctioned way of getting at
the format objects, which most people should use.

The Bio::SeqIO system can be thought of like biological file handles.
They are attached to filehandles with smart formatting rules (eg,
genbank format, or EMBL format, or binary trace file format) and 
can either read or write sequence objects (Bio::Seq objects, or
more correctly, Bio::SeqI implementing objects, of which Bio::Seq is
one such object). If you want to know what to do with a Bio::Seq
object, read L<Bio::Seq>.

The idea is that you request a stream object for a particular format.
All the stream objects have a notion of an internal file that is read
from or written to. A particular SeqIO object instance is configured
for either input or output. A specific example of a stream object is
the Bio::SeqIO::fasta object.

Each stream object has functions

   $stream->next_seq();

and

   $stream->write_seq($seq);

also

   $stream->type() # returns 'INPUT' or 'OUTPUT'

As an added bonus, you can recover a filehandle that is tied to the
SeqIO object, allowing you to use the standard E<lt>E<gt> and print operations
to read and write sequence objects:

    use Bio::SeqIO;

    $stream = Bio::SeqIO->newFh(-format => 'Fasta'); # read from standard input

    while ( $seq = <$stream> ) {
	# do something with $seq
    }

and

    print $stream $seq; # when stream is in output mode

This makes the simplest ever reformatter

    #!/usr/local/bin/perl

    $format1 = shift;
    $format2 = shift || die "Usage: reformat format1 format2 < input > output";

    use Bio::SeqIO;

    $in  = Bio::SeqIO->newFh(-format => $format1 );
    $out = Bio::SeqIO->newFh(-format => $format2 );
    #note: you might want to quote -format to keep older perl's from complaining.

    print $out $_ while <$in>;


=head1 CONSTRUCTORS

=head2 Bio::SeqIO-E<gt>new()

   $seqIO = Bio::SeqIO->new(-file => 'filename',   -format=>$format);
   $seqIO = Bio::SeqIO->new(-fh   => \*FILEHANDLE, -format=>$format);
   $seqIO = Bio::SeqIO->new(-format => $format);

The new() class method constructs a new Bio::SeqIO object.  The
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

   $seqIO = Bio::SeqIO->new(-fh => \*STDIN);

Note that you must pass filehandles as references to globs.

If neither a filehandle nor a filename is specified, then the module
will read from the @ARGV array or STDIN, using the familiar E<lt>E<gt>
semantics.

A string filehandle is handy if you want to modify the output in the
memory, before printing it out. The following program reads in EMBL
formatted entries from a file and prints them out in fasta format with
some HTML tags:

  use Bio::SeqIO;
  use IO::String;
  my $in  = Bio::SeqIO->new('-file' => "emblfile" , 
  			    '-format' => 'EMBL');
  while ( my $seq = $in->next_seq() ) {
      # the output handle is reset for every file
      my $stringio = IO::String->new($string);
      my $out = Bio::SeqIO->new('-fh' => $stringio,
  			        '-format' => 'fasta');
      # output goes into $string
      $out->write_seq($seq);
      # modify $string
      $string =~ s|(>)(\w+)|$1<font color="Red">$2</font>|g;
      # print into STDOUT
      print $string;
  }

=item -format

Specify the format of the file.  Supported formats include:

   Fasta       FASTA format
   EMBL        EMBL format
   GenBank     GenBank format
   swiss       Swissprot format
   PIR         Protein Information Resource format
   GCG         GCG format
   raw         Raw format (one sequence per line, no ID)
   ace         ACeDB sequence format
   game        GAME XML format
   phd         phred output
   qual        Quality values (get a sequence of quality scores)
   Fastq       Fastq format
   SCF         SCF tracefile format
   ABI         ABI tracefile format
   ALF         ALF tracefile format
   CTF         CTF tracefile format
   ZTR         ZTR tracefile format
   PLN         Staden plain tracefile format
   EXP         Staden tagged experiment tracefile format

If no format is specified and a filename is given then the module
will attempt to deduce the format from the filename suffix.  If this
is unsuccessful then Fasta format is assumed.

The format name is case insensitive.  'FASTA', 'Fasta' and 'fasta' are
all valid suffixes.

Currently, the tracefile formats (except for SCF) require installation
of the external Staden "io_lib" package, as well as the
Bio::SeqIO::staden::read package available from the bioperl-ext
repository.

=item -flush

By default, all files (or filehandles) opened for writing sequences
will be flushed after each write_seq() (making the file immediately
usable).  If you don't need this facility and would like to marginally
improve the efficiency of writing multiple sequences to the same file
(or filehandle), pass the -flush option '0' or any other value that
evaluates as defined but false:

  my $gb = new Bio::SeqIO -file   => "<gball.gbk",
                          -format => "gb";
  my $fa = new Bio::SeqIO -file   => ">gball.fa",
                          -format => "fasta",
                          -flush  => 0; # go as fast as we can!
  while($seq = $gb->next_seq) { $fa->write_seq($seq) }


=back

=head2 Bio::SeqIO-E<gt>newFh()

   $fh = Bio::SeqIO->newFh(-fh   => \*FILEHANDLE, -format=>$format);
   $fh = Bio::SeqIO->newFh(-format => $format);
   # etc.

This constructor behaves like new(), but returns a tied filehandle
rather than a Bio::SeqIO object.  You can read sequences from this
object using the familiar E<lt>E<gt> operator, and write to it using
print().  The usual array and $_ semantics work.  For example, you can
read all sequence objects into an array like this:

  @sequences = <$fh>;

Other operations, such as read(), sysread(), write(), close(), and printf() 
are not supported.

=head1 OBJECT METHODS

See below for more detailed summaries.  The main methods are:

=head2 $sequence = $seqIO-E<gt>next_seq()

Fetch the next sequence from the stream.

=head2 $seqIO-E<gt>write_seq($sequence [,$another_sequence,...])

Write the specified sequence(s) to the stream.

=head2 TIEHANDLE(), READLINE(), PRINT()

These provide the tie interface.  See L<perltie> for more details.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
to one of the Bioperl mailing lists.

Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/MailList.shtml      - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney, Lincoln Stein

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#' Let the code begin...

package Bio::SeqIO;

use strict;
use vars qw(@ISA);

use Bio::Root::Root;
use Bio::Root::IO;
use Bio::Factory::SequenceStreamI;
use Symbol();

@ISA = qw(Bio::Root::Root Bio::Root::IO Bio::Factory::SequenceStreamI);

sub BEGIN {
    eval { require Bio::SeqIO::staden::read; };
}

=head2 new

 Title   : new
 Usage   : $stream = Bio::SeqIO->new(-file => $filename, -format => 'Format')
 Function: Returns a new seqstream
 Returns : A Bio::SeqIO stream initialised with the appropriate format
 Args    : -file => $filename
           -format => format
           -fh => filehandle to attach to

See L<Bio::SeqIO::Handler>

=cut

my $entry = 0;

sub new {
    my ($caller,@args) = @_;
    my $class = ref($caller) || $caller;
    
    # or do we want to call SUPER on an object if $caller is an
    # object?
    if( $class =~ /Bio::SeqIO::(\S+)/ ) {
	my ($self) = $class->SUPER::new(@args);	
	$self->_initialize(@args);
	return $self;
    } else { 

	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	my $format = $param{'-format'} || 
	    $class->_guess_format( $param{-file} || $ARGV[0] ) ||
		'fasta';
	$format = "\L$format";	# normalize capitalization to lower case

	# normalize capitalization
	return undef unless( &_load_format_module($format) );
	return "Bio::SeqIO::$format"->new(@args);
    }
}

=head2 newFh

 Title   : newFh
 Usage   : $fh = Bio::SeqIO->newFh(-file=>$filename,-format=>'Format')
 Function: does a new() followed by an fh()
 Example : $fh = Bio::SeqIO->newFh(-file=>$filename,-format=>'Format')
           $sequence = <$fh>;   # read a sequence object
           print $fh $sequence; # write a sequence object
 Returns : filehandle tied to the Bio::SeqIO::Fh class
 Args    :

See L<Bio::SeqIO::Fh>

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
           $sequence = <$fh>;   # read a sequence object
           print $fh $sequence; # write a sequence object
 Returns : filehandle tied to Bio::SeqIO class
 Args    : none

=cut


sub fh {
  my $self = shift;
  my $class = ref($self) || $self;
  my $s = Symbol::gensym;
  tie $$s,$class,$self;
  return $s;
}

# _initialize is chained for all SeqIO classes

sub _initialize {
    my($self, @args) = @_;

    # flush is initialized by the Root::IO init

    my ($seqfact) = $self->_rearrange([qw(SEQFACTORY )], @args);

    $seqfact && $self->sequence_factory($seqfact);

    # initialize the IO part
    $self->_initialize_io(@args);
}

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = stream->next_seq
 Function: Reads the next sequence object from the stream and returns it.

           Certain driver modules may encounter entries in the stream that
           are either misformatted or that use syntax not yet understood
           by the driver. If such an incident is recoverable, e.g., by
           dismissing a feature of a feature table or some other non-mandatory
           part of an entry, the driver will issue a warning. In the case
           of a non-recoverable situation an exception will be thrown.
           Do not assume that you can resume parsing the same stream after
           catching the exception. Note that you can always turn recoverable
           errors into exceptions by calling $stream->verbose(2).
 Returns : a Bio::Seq sequence object
 Args    : none

See L<Bio::Root::RootI>, L<Bio::Factory::SeqStreamI>, L<Bio::Seq>

=cut

sub next_seq {
   my ($self, $seq) = @_;
   $self->throw("Sorry, you cannot read from a generic Bio::SeqIO object.");
}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq($seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq object

=cut

sub write_seq {
    my ($self, $seq) = @_;
    $self->throw("Sorry, you cannot write to a generic Bio::SeqIO object.");
}


=head2 alphabet

 Title   : alphabet
 Usage   : $self->alphabet($newval)
 Function: Set/get the molecule type for the Seq objects to be created.
 Example : $seqio->alphabet('protein')
 Returns : value of alphabet: 'dna', 'rna', or 'protein'
 Args    : newvalue (optional)
 Throws  : Exception if the argument is not one of 'dna', 'rna', or 'protein'

=cut

sub alphabet {
   my ($self, $value) = @_;

   if ( defined $value) {
       # instead of hard-coding the allowed values once more, we check by
       # creating a dummy sequence object
       eval {
	   require Bio::PrimarySeq;
	   my $seq = Bio::PrimarySeq->new('-alphabet' => $value);
       };
       if($@) {
	   $self->throw("Invalid alphabet: $value\n. See Bio::PrimarySeq for allowed values.");
       }
       $self->{'alphabet'} = "\L$value";
   }
   return $self->{'alphabet'};
}

=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL SeqIO stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns :
 Args    :

=cut

sub _load_format_module {
  my ($format) = @_;
  my ($module, $load, $m);

  # untaint operation for safe web-based running
  if ($format =~ /^([\w:]+)$/) {
    $format = $1;
  } else {
    print STDERR <<END;
_load_format_module: $format is an illegal perl package name
For more information about the SeqIO system please see the SeqIO docs.
This includes ways of checking for formats at compile time, not run time
END
  ;
    return;
  }

  $module = "_<Bio/SeqIO/$format.pm";
  $load = "Bio/SeqIO/$format.pm";

  return 1 if $main::{$module};
  eval {
    require $load;
  };
  if ( $@ ) {
    print STDERR <<END;
$load: $format cannot be found
Exception $@
For more information about the SeqIO system please see the SeqIO docs.
This includes ways of checking for formats at compile time, not run time
END
  ;
    return;
  }
  return 1;
}

=head2 _concatenate_lines

 Title   : _concatenate_lines
 Usage   : $s = _concatenate_lines($line, $continuation_line)
 Function: Private. Concatenates two strings assuming that the second stems
           from a continuation line of the first. Adds a space between both
           unless the first ends with a dash.

           Takes care of either arg being empty.
 Example :
 Returns : A string.
 Args    :

=cut

sub _concatenate_lines {
    my ($self, $s1, $s2) = @_;

    $s1 .= " " if($s1 && ($s1 !~ /-$/) && $s2);
    return ($s1 ? $s1 : "") . ($s2 ? $s2 : "");
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
 Notes   : formats that _filehandle() will guess include fasta,
           genbank, scf, pir, embl, raw, gcg, ace, bsml, swissprot,
           fastq and phd/phred

=cut

sub _guess_format {
   my $class = shift;
   return unless $_ = shift;
   return 'fasta'   if /\.(fasta|fast|seq|fa|fsa|nt|aa)$/i;
   return 'genbank' if /\.(gb|gbank|genbank)$/i;
   return 'scf'     if /\.scf$/i;
   return 'scf'     if /\.scf$/i;
   return 'abi'     if /\.abi$/i;
   return 'alf'     if /\.alf$/i;
   return 'ctf'     if /\.ctf$/i;
   return 'ztr'     if /\.ztr$/i;
   return 'pln'     if /\.pln$/i;
   return 'exp'     if /\.exp$/i;
   return 'pir'     if /\.pir$/i;
   return 'embl'    if /\.(embl|ebl|emb|dat)$/i;
   return 'raw'     if /\.(txt)$/i;
   return 'gcg'     if /\.gcg$/i;
   return 'ace'     if /\.ace$/i;
   return 'bsml'    if /\.(bsm|bsml)$/i;
   return 'swiss'   if /\.(swiss|sp)$/i;
   return 'phd'     if /\.(phd|phred)$/i;
   return 'fastq'   if /\.fastq$/i;
}

sub DESTROY {
    my $self = shift;

    $self->close();
}

sub TIEHANDLE {
    my ($class,$val) = @_;
    return bless {'seqio' => $val}, $class;
}

sub READLINE {
  my $self = shift;
  return $self->{'seqio'}->next_seq() unless wantarray;
  my (@list, $obj);
  push @list, $obj while $obj = $self->{'seqio'}->next_seq();
  return @list;
}

sub PRINT {
  my $self = shift;
  $self->{'seqio'}->write_seq(@_);
}

=head2 sequence_factory

 Title   : sequence_factory
 Usage   : $seqio->sequence_factory($seqfactory)
 Function: Get/Set the Bio::Factory::SequenceFactoryI
 Returns : Bio::Factory::SequenceFactoryI
 Args    : [optional] Bio::Factory::SequenceFactoryI


=cut

sub sequence_factory{
   my ($self,$obj) = @_;   
   if( defined $obj ) {
       if( ! ref($obj) || ! $obj->isa('Bio::Factory::SequenceFactoryI') ) {
	   $self->throw("Must provide a valid Bio::Factory::SequenceFactoryI object to ".ref($self)." sequence_factory()");
       }
       $self->{'_seqio_seqfactory'} = $obj;
   }
   $self->{'_seqio_seqfactory'};
}

1;

