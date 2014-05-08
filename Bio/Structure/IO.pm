#
# BioPerl module for Bio::Structure::IO
#
# Copyright 2001, 2002 Kris Boulez
#
# You may distribute this module under the same terms as perl itself
#
# _history
# October 18, 1999  Largely rewritten by Lincoln Stein
# November 16, 2001 Copied Bio::SeqIO to Bio::Structure::IO and modified
# 			where needed. Factoring out common methods
# 			(to Bio::Root::IO) might be a good idea.

# POD documentation - main docs before the code

=head1 NAME

Bio::Structure::IO - Handler for Structure Formats

=head1 SYNOPSIS

    use Bio::Structure::IO;

    $in  = Bio::Structure::IO->new(-file => "inputfilename",
                                   -format => 'pdb');

    while ( my $struc = $in->next_structure() ) {
       print "Structure ", $struc->id, " number of models: ",
             scalar $struc->model,"\n";
    }

=head1 DESCRIPTION

Bio::Structure::IO is a handler module for the formats in the
Structure::IO set (e.g. L<Bio::Structure::IO::pdb>). It is the officially
sanctioned way of getting at the format objects, which most people
should use.

The Bio::Structure::IO system can be thought of like biological file
handles.  They are attached to filehandles with smart formatting rules
(e.g. PDB format) and can either read or write structure objects
(Bio::Structure objects, or more correctly, Bio::Structure::StructureI
implementing objects, of which Bio::Structure is one such object). If
you want to know what to do with a Bio::Structure object, read
L<Bio::Structure>.

The idea is that you request a stream object for a particular format.
All the stream objects have a notion of an internal file that is read
from or written to. A particular Structure::IO object instance is
configured for either input or output. A specific example of a stream
object is the Bio::Structure::IO::pdb object.

Each stream object has functions

   $stream->next_structure();

and

   $stream->write_structure($struc);

also

   $stream->type() # returns 'INPUT' or 'OUTPUT'

As an added bonus, you can recover a filehandle that is tied to the
Structure::IOIO object, allowing you to use the standard E<lt>E<gt>
and print operations to read and write structure::IOuence objects:

    use Bio::Structure::IO;

    $stream = Bio::Structure::IO->newFh(-format => 'pdb'); # read from standard input

    while ( $structure = <$stream> ) {
   	# do something with $structure
    }

and

    print $stream $structure; # when stream is in output mode


=head1 CONSTRUCTORS

=head2 Bio::Structure::IO-E<gt>new()

   $stream = Bio::Structure::IO->new(-file => 'filename',   -format=>$format);
   $stream = Bio::Structure::IO->new(-fh   => \*FILEHANDLE, -format=>$format);
   $stream = Bio::Structure::IO->new(-format => $format);

The new() class method constructs a new Bio::Structure::IO object. The
returned object can be used to retrieve or print Bio::Structure
objects.  new() accepts the following parameters:

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

   $strucIO = Bio::Structure::IO->new(-fh => \*STDIN);

Note that you must pass filehandles as references to globs.

If neither a filehandle nor a filename is specified, then the module
will read from the @ARGV array or STDIN, using the familiar E<lt>E<gt>
semantics.

=item -format

Specify the format of the file.  Supported formats include:

   pdb         Protein Data Bank format

If no format is specified and a filename is given, then the module
will attempt to deduce it from the filename.  If this is unsuccessful,
PDB format is assumed.

The format name is case insensitive.  'PDB', 'Pdb' and 'pdb' are
all supported.

=back

=head2 Bio::Structure::IO-E<gt>newFh()

   $fh = Bio::Structure::IO->newFh(-fh   => \*FILEHANDLE, -format=>$format);
   $fh = Bio::Structure::IO->newFh(-format => $format);
   # etc.

This constructor behaves like new(), but returns a tied filehandle
rather than a Bio::Structure::IO object.  You can read structures from this
object using the familiar E<lt>E<gt> operator, and write to it using
print().  The usual array and $_ semantics work.  For example, you can
read all structure objects into an array like this:

  @structures = <$fh>;

Other operations, such as read(), sysread(), write(), close(), and printf()
are not supported.

=head1 OBJECT METHODS

See below for more detailed summaries.  The main methods are:

=head2 $structure = $structIO-E<gt>next_structure()

Fetch the next structure from the stream.

=head2 $structIO-E<gt>write_structure($struc [,$another_struc,...])

Write the specified structure(s) to the stream.

=head2 TIEHANDLE(), READLINE(), PRINT()

These provide the tie interface.  See L<perltie> for more details.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.
Bug reports can be submitted via the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHORS - Ewan Birney, Lincoln Stein, Kris Boulez

Email birney@ebi.ac.uk, lstein@cshl.org, kris.boulez@algonomics.com


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Structure::IO;

use strict;

use Bio::PrimarySeq;
use Symbol;

use base qw(Bio::Root::Root Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : $stream = Bio::Structure::IO->new(-file => $filename, -format => 'Format')
 Function: Returns a new structIOstream
 Returns : A Bio::Structure::IO handler initialised with the appropriate format
 Args    : -file => $filename
           -format => format
           -fh => filehandle to attach to

=cut

my $entry = 0;

sub new {
    my ($caller,@args) = @_;
    my $class = ref($caller) || $caller;

    # or do we want to call SUPER on an object if $caller is an
    # object?
    if( $class =~ /Bio::Structure::IO::(\S+)/ ) {
	my ($self) = $class->SUPER::new(@args);
	$self->_initialize(@args);
	return $self;
    } else {

	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	my $format = $param{'-format'} ||
	    $class->_guess_format( $param{-file} || $ARGV[0] ) ||
		'pdb';
	$format = "\L$format";	# normalize capitalization to lower case

	# normalize capitalization
	return unless( &_load_format_module($format) );
	return "Bio::Structure::IO::$format"->new(@args);
    }
}

=head2 newFh

 Title   : newFh
 Usage   : $fh = Bio::Structure::IO->newFh(-file=>$filename,-format=>'Format')
 Function: does a new() followed by an fh()
 Example : $fh = Bio::Structure::IO->newFh(-file=>$filename,-format=>'Format')
           $structure = <$fh>;   # read a structure object
           print $fh $structure; # write a structure object
 Returns : filehandle tied to the Bio::Structure::IO::Fh class
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
           $structure = <$fh>;   # read a structure object
           print $fh $structure; # write a structure object
 Returns : filehandle tied to the Bio::Structure::IO::Fh class
 Args    :

=cut


sub fh {
  my $self = shift;
  my $class = ref($self) || $self;
  my $s = Symbol::gensym;
  tie $$s,$class,$self;
  return $s;
}


=head2 format

 Title   : format
 Usage   : $format = $obj->format()
 Function: Get the structure format
 Returns : structure format
 Args    : none

=cut

# format() method inherited from Bio::Root::IO


# _initialize is chained for all SeqIO classes

sub _initialize {
    my($self, @args) = @_;

    # not really necessary unless we put more in RootI
    $self->SUPER::_initialize(@args);

    # initialize the IO part
    $self->_initialize_io(@args);
}

=head2 next_structure

 Title   : next_structure
 Usage   : $structure = stream->next_structure
 Function: Reads the next structure object from the stream and returns a
           Bio::Structure::Entry object.

           Certain driver modules may encounter entries in the stream that
           are either misformatted or that use syntax not yet understood
           by the driver. If such an incident is recoverable, e.g., by
           dismissing a feature of a feature table or some other non-mandatory
           part of an entry, the driver will issue a warning. In the case
           of a non-recoverable situation an exception will be thrown.
           Do not assume that you can resume parsing the same stream after
           catching the exception. Note that you can always turn recoverable
           errors into exceptions by calling $stream->verbose(2) (see
           Bio::RootI POD page).
 Returns : a Bio::Structure::Entry object
 Args    : none

=cut

sub next_structure {
   my ($self, $struc) = @_;
   $self->throw("Sorry, you cannot read from a generic Bio::Structure::IO object.");
}

# Do we want people to read out the sequence directly from a $structIO stream
#
##=head2 next_primary_seq
##
## Title   : next_primary_seq
## Usage   : $seq = $stream->next_primary_seq
## Function: Provides a primaryseq type of sequence object
## Returns : A Bio::PrimarySeqI object
## Args    : none
##
##
##=cut
##
##sub next_primary_seq {
##   my ($self) = @_;
##
##   # in this case, we default to next_seq. This is because
##   # Bio::Seq's are Bio::PrimarySeqI objects. However we
##   # expect certain sub classes to override this method to provide
##   # less parsing heavy methods to retrieving the objects
##
##   return $self->next_seq();
##}

=head2 write_structure

 Title   : write_structure
 Usage   : $stream->write_structure($structure)
 Function: writes the $structure object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Structure object

=cut

sub write_seq {
    my ($self, $struc) = @_;
    $self->throw("Sorry, you cannot write to a generic Bio::Structure::IO object.");
}


# De we need this here
#
##=head2 alphabet
##
## Title   : alphabet
## Usage   : $self->alphabet($newval)
## Function: Set/get the molecule type for the Seq objects to be created.
## Example : $seqio->alphabet('protein')
## Returns : value of alphabet: 'dna', 'rna', or 'protein'
## Args    : newvalue (optional)
## Throws  : Exception if the argument is not one of 'dna', 'rna', or 'protein'
##
##=cut
##
##sub alphabet {
##   my ($self, $value) = @_;
##
##   if ( defined $value) {
##       # instead of hard-coding the allowed values once more, we check by
##       # creating a dummy sequence object
##       eval {
##	   my $seq = Bio::PrimarySeq->new('-alphabet' => $value);
##       };
##       if($@) {
##	   $self->throw("Invalid alphabet: $value\n. See Bio::PrimarySeq for allowed values.");
##       }
##       $self->{'alphabet'} = "\L$value";
##   }
##   return $self->{'alphabet'};
##}

=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL Structure::IO stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns :
 Args    :

=cut

sub _load_format_module {
  my ($format) = @_;
  my ($module, $load, $m);

  $module = "_<Bio/Structure/IO/$format.pm";
  $load = "Bio/Structure/IO/$format.pm";

  return 1 if $main::{$module};
  eval {
    require $load;
  };
  if ( $@ ) {
    print STDERR <<END;
$load: $format cannot be found
Exception $@
For more information about the Structure::IO system please see the
Bio::Structure::IO docs.  This includes ways of checking for formats at
compile time, not run time
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
   return 'embl'    if /\.(embl|ebl|emb|dat)$/i;
   return 'raw'     if /\.(txt)$/i;
   return 'gcg'     if /\.gcg$/i;
   return 'ace'     if /\.ace$/i;
   return 'bsml'    if /\.(bsm|bsml)$/i;
   return 'pdb'     if /\.(ent|pdb)$/i;
}

sub DESTROY {
    my $self = shift;

    $self->close();
}

sub TIEHANDLE {
    my ($class,$val) = @_;
    return bless {'structio' => $val}, $class;
}

sub READLINE {
  my $self = shift;
  return $self->{'structio'}->next_seq() || undef unless wantarray;
  my (@list, $obj);
  push @list, $obj while $obj = $self->{'structio'}->next_seq();
  return @list;
}

sub PRINT {
  my $self = shift;
  $self->{'structio'}->write_seq(@_);
}

1;

