
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
    # note: we quote -format to keep older perl's from complaining.

    while ( $seq = $in->next_seq() ) {
	$out->write_seq($seq);
    }

or

    use Bio::SeqIO;

    $in  = Bio::SeqIO->newFh(-file => "inputfilename" , '-format' => 'Fasta');
    $out = Bio::SeqIO->newFh('-format' => 'EMBL');

    # World's shortest Fasta<->EMBL format converter:
    print $output $_ while <$in>;

=head1 DESCRIPTION

Bio::SeqIO is a handler module for the formats in the SeqIO set (eg,
Bio::SeqIO::fasta). It is the officially sanctioned way of getting at
the format objects, which most people should use.

The SeqIO system replaces the old parse_XXX functions in the Seq
object.

The idea is that you request a stream object for a particular format.
All the stream objects have a notion of an internal file that is read
from or written to (the same object handles both input and output).  A
physical example of a stream object is the Bio::SeqIO::Fasta object.

Each stream object has functions

   $stream->next_seq();

and

   $stream->write_seq($seq);

also

   $stream->type() # returns 'INPUT' or 'OUTPUT'

As an added bonus, you can recover a filehandle that is tied to the
SeqIO object, allowing you to use the standard <> and print operations
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


Notice that the reformatter will only convert information that is held in
the Seq object, which at the  moment is only the sequence and the id. More
information will be converted through the expanded or larger object which
the bioperl developers are talking about.

It is not good for reformatting genbank to embl therefore, but was never
designed for this task anyway.

=head1 CONSTRUCTORS

=head2 Bio::SeqIO->new()

   $seqIO = Bio::SeqIO->new(-file => 'filename',   -format=>$format);
   $seqIO = Bio::SeqIO->new(-fh   => \*FILEHANDLE, -format=>$format);
   $seqIO = Bio::SeqIO->new(-format => $format);

The new() class method constructs a new Bio::SeqIO object.  The
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

   $seqIO = Bio::SeqIO->new(-fh => \*STDIN);

Note that you must pass filehandles as references to globs.

If neither a filehandle nor a filename is specified, then the module
will read from the @ARGV array or STDIN, using the familiar <>
semantics.

=item -format

Specify the format of the file.  Supported formats include:

   Fasta       FASTA format
   EMBL        EMBL format
   SCF         SCF tracefile format
   PIR         Protein Information Resource format
   GCG         GCG format
   raw         Raw format (one sequence per line, no ID)

If no format is specified and a filename is given, then the module
will attempt to deduce it from the filename.  If this is unsuccessful,
Fasta format is assumed.

The format name is case insensitive.  'FASTA', 'Fasta' and 'fasta' are
all supported.

=back

=head2 Bio::SeqIO->newFh()

   $fh = Bio::SeqIO->newFh(-fh   => \*FILEHANDLE, -format=>$format);
   $fh = Bio::SeqIO->newFh(-format => $format);
   # etc.

This constructor behaves like new(), but returns a tied filehandle
rather than a Bio::SeqIO object.  You can read sequences from this
object using the familiar <> operator, and write to it using print().
The usual array and $_ semantics work.  For example, you can read all
sequence objects into an array like this:

  @sequences = <$fh>;

Other operations, such as read(), sysread(), write() and printf() are
not supported.

=head1 OBJECT METHODS

See below for more detailed summaries.  The main methods are:

=head2 $sequence = $seqIO->next_seq()

Fetch the next sequence from the stream.

=head2 $seqIO->write_seq($sequence [,$another_sequence,...])

Write the specified sequence(s) to the stream.

=head2 TIEHANDLE(), READLINE(), PRINT()

These provide the tie interface.  See L<perltie> for more details.

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

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO;

use strict;
use vars '@ISA';

use Bio::Root::Object;
use Bio::Seq;
use Symbol();

@ISA = 'Bio::Root::Object';

=head2 new

 Title   : new
 Usage   : $stream = Bio::SeqIO->new(-file => $filename, -format => 'Format')
 Function: Returns a new seqstream
 Returns : A Bio::SeqIO::Handler initialsed with the appropiate format
 Args    : -file => $filename 
           -format => format
           -fh => filehandle to attach to

=cut

my $entry = 0;

sub new {
   my ($class,%param) = @_;
   my ($format);
   my ($handler,$stream);

   @param{ map { lc $_ } keys %param } = values %param;  # lowercase keys
   $format = $param{'-format'} 
             || $class->_guess_format( $param{-file} || $ARGV[0] )
             || 'Fasta';
   $format = "\L$format"; # normalize capitalization to lower case

   if( &_load_format_module($format) == 0 ) { # normalize capitalization
       return undef;
   }

   $stream = "Bio::SeqIO::$format"->_new(%param);
   return $stream;
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
 Returns : filehandle tied to the Bio::SeqIO::Fh class
 Args    : 

=cut


sub fh {
  my $self = shift;
  my $class = ref($self) || $self;
  my $s = Symbol::gensym;
  tie $$s,$class,$self;
  return $s;
}

sub _new {
  my $self = shift;
  $self->SUPER::new(@_);
}

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;  # exception handling code

  my ($file,$fh) = $self->_rearrange([qw(
					 FILE
					 FH
					)],
				     @args,
				     );
  if( defined $file and defined $fh ) {
      $self->throw("Providing both a file and a filehandle for reading from - only one please!");
  }

  if( defined $file and $file ne '' ) {
    $fh = Symbol::gensym();
    open ($fh,$file) || $self->throw("Could not open $file for Fasta stream reading $!");
  }
  $self->_filehandle($fh) if defined $fh;
  return $make; # success - we hope!
}

=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL SeqIO stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns : 
 Args    :

=cut
  
sub _load_format_module{
  my ($format) = @_;
  my ($module,$load,$m);
  
  $module = "_<Bio/SeqIO/$format.pm";
  $load = "Bio/SeqIO/$format.pm";
  
  return 1 if $main::{$module};
  eval {
    require $load;
  };
  if( $@ ) {
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

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = stream->next_seq
 Function: reads the next $seq object from the stream
 Returns : the Bio::Seq sequence object
 Args    : 

=cut

sub next_seq {
   my ($self,$seq) = @_;
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
   my ($self,$seq) = @_;
   $self->throw("Sorry, you cannot write to a generic Bio::SeqIO object.");
}

=head2 _print

 Title   : _print
 Usage   : $obj->_print(@lines)
 Function: 
 Example : 
 Returns : writes output

=cut

sub _print {
  my $self = shift;
  my $fh = $self->_filehandle || \*STDOUT;
  print $fh @_;
}

=head2 _readline

 Title   : _readline
 Usage   : $obj->_readline($newval)
 Function: 
 Example : 
 Returns : reads a line of input

=cut

sub _readline {
  my $self = shift;
  my $fh = $self->_filehandle;
  my $line = defined($fh) ? <$fh> : <>;
  return $line;
}

=head2 _filehandle

 Title   : _filehandle
 Usage   : $obj->_filehandle($newval)
 Function: 
 Example : 
 Returns : value of _filehandle
 Args    : newvalue (optional)


=cut

sub _filehandle {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_filehandle'} = $value;
    }
    return $obj->{'_filehandle'};
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
   return 'fasta'   if /\.(fasta|fast|seq|fa|fsa|nt|aa)$/i;
   return 'genbank' if /\.(gb|gbank|genbank)$/i;
   return 'scf'     if /\.scf$/i;
   return 'pir'     if /\.pir$/i;
   return 'embl'    if /\.(embl|ebl|emb)$/i;
   return 'raw'     if /\.(txt)$/i;
   return 'gcg'     if /\.gcg$/i;
}

sub TIEHANDLE {
  my $class = shift;
  return bless {seqio => shift},$class;
}

sub READLINE {
  my $self = shift;
  return $self->{seqio}->next_seq() unless wantarray;
  my (@list,$obj);
  push @list,$obj  while $obj = $self->{seqio}->next_seq();
  return @list;
}

sub PRINT {
  my $self = shift;
  $self->{seqio}->write_seq(@_);
}

1;

