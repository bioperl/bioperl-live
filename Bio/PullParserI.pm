#
# BioPerl module for Bio::PullParserI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PullParserI - A base module for fast 'pull' parsing

=head1 SYNOPSIS

    # do not use this class, it is intended for parser module
    # writers only

=head1 DESCRIPTION

If you are writing a module to parse some new format, you may wish to use
a 'pull' approach whereby you only do work (reading file data, parsing it,
turning the parsed data in an object) when absolutely necessary.

PullParserI provides a system for doing exactly that. As a PullParser you
need a chunk. A chunk is just a Bio::Root::IO that contains all the raw data
you would want to parse. You can use the chunk() method to create a chunk from
a filename, existing filehandle or even a string. If you make a chunk from a
large file, but actually only want your chunk to be some portion of the whole
file, supply start and end amounts in bytes to chunk() at the same time.
The methods _chunk_seek() and _chunk_tell() provide seeks and tells that are
relative to the start and end of your chunk, not the whole file.

The other thing you will need to decide when making a chunk is how to handle
piped input. A PullParser typically needs seekable data to parse, so if your
data is piped in and unseekable, you must decide between creating a temp file
or reading the input into memory, which will be done before the chunk becomes
usable and you can begin any parsing. Alternatively you can choose to force
a sequential read, in which case you can make use of _dependencies() to define
the linear order of methods that would result in the file being read
sequentially. The return value of _sequential() is also useful here, if you
would need to cache some data or otherwise behave differently during a
sequential read.

The main method in the system is get_field(). This method relies on the
existance of a private hash reference accessible to it with the method
_fields(). That hash ref should have as keys all the sorts of data you will want
to parse (eg. 'score'), and prior to parsing the values would be undefined. A
user of your module can then call either $module-E<gt>get_field('score') or
$module-E<gt>score and get_field will either return the answer from
$self-E<gt>_fields-E<gt>{score} if it is defined, or call a method _discover_score()
first if not. So for the system to work you need to define a _discover_*()
method for every field in the fields hash, and ensure that the method stores an
answer in the fields hash.

How you implement your _discover_* methods is up to you, though you should never
call a _discover_* method directly yourself; always use get_field(), since
get_field() will deal with calling dependent methods for you if a forced
sequenctial read is in progress due to piped input. You will almost certainly
want to make use of the various chunk-related methods of this class (that are
denoted private by the leading '_'; this means you can use them as the author of
a parser class, but users of your parser should not). 

Primary amongst them is _*_chunk_by_end() to which you provide text that
represents the end of your desired chunk and it does a readline with your
argument as $/. The chunk knows about its line-endings, so if you want your
end definition to include a new line, just always use "\n" and PullParserI will
do any necessary conversion for you.

If your input data is hierarchical (eg. report-E<gt>many results-E<gt>many hits-E<gt>many
hsps), and you want an object at the leaf of the hierarchy to have access to
information that is shared amongst all of them (is parsed in the root), you
don't have to copy the data to each leaf object; simply by defining parent(),
when you call get_field() and the requested field isn't in your leaf's fields
hash, the leaf's parent will be asked for the field instead, and so on till
root.

See Bio::SearchIO::hmmer_pull for an example of implementing a parser using
PullParserI.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 CONTRIBUTORS

Inspired by a posting by Aaron J. Mackey

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::PullParserI;

use vars qw($AUTOLOAD $FORCE_TEMP_FILE);
use strict;

use Bio::Root::IO;

use base qw(Bio::Root::RootI);

BEGIN {
    # chunk() needs perl 5.8 feature for modes other than temp_file, so will
    # workaround by forcing temp_file mode in <5.8. Could also rewrite using
    # IO::String, but don't want to.
    if ($] < 5.008) {
        $FORCE_TEMP_FILE = 1;
    }
}

=head2 _fields

 Title   : _fields
 Usage   : $obj->_fields( { field1 => undef } );
           my $fields_ref = $obj->_fields;
 Function: Get/set the hash reference containing all the fields for this parser
 Returns : hash ref
 Args    : none to get, OR hash ref to set

=cut

sub _fields {
    my $self = shift;
    if (@_) {
        $self->{_fields} = shift;
    }
    unless (defined $self->{_fields}) {
        $self->{_fields} = { };
    }
    return $self->{_fields};
}

=head2 has_field

 Title   : has_field
 Usage   : if ($obj->has_field('field_name') {...}
 Function: Ask if a particular object has a given field (doesn't ask ancestors)
 Returns : boolean
 Args    : string (the field name to test)

=cut

sub has_field {
    my ($self, $desired) = @_;
    $desired || return;
    return exists $self->_fields->{$desired};
}

=head2 get_field

 Title   : get_field
 Usage   : my $field_value = $obj->get_field('field_name');
 Function: Get the value of a given field. If this $obj doesn't have the field,
           it's parent() will be asked, and so on until there are no more
           parents.
 Returns : scalar, warns if a value for the field couldn't be found and returns
           undef.
 Args    : string (the field to get)

=cut

sub get_field {
    my $self = shift;
    my $desired = shift || return keys %{$self->_fields};
    if (exists $self->_fields->{$desired}) {
        unless (defined $self->_fields->{$desired}) {
            my $method = '_discover_'.$desired;
            
            my $dependency = $self->_dependencies($desired);
            if ($dependency && ! defined $self->_fields->{$dependency}) {
                $self->get_field($dependency);
            }
            
            # it might exist now
            $self->$method unless defined $self->_fields->{$desired};
        }
        return $self->_fields->{$desired};
    }
    
    # is it a field of our parent? (checks all ancestors)
    if (my $parent = $self->parent) {
        return $parent->get_field($desired);
    }
    
    $desired =~ s/_discover_//;
    $self->warn("This report does not hold information about '$desired'");
    return;
}

=head2 parent

 Title   : parent
 Usage   : $obj->parent($parent_obj);
           my $parent_obj = $obj->parent;
 Function: Get/set the parent object of this one.
 Returns : Bio::PullParserI
 Args    : none to get, OR Bio::PullParserI to set

=cut

sub parent {
    my $self = shift;
    if (@_) { $self->{parent} = shift }
    return $self->{parent} || return;
}

=head2 chunk

 Title   : chunk
 Usage   : $obj->chunk($filename);
           my $chunk = $obj->chunk;
 Function: Get/set the chunk of this parser.
 Returns : Bio:Root::IO
 Args    : none to get, OR
           First argument of a GLOB reference, filename string, string data to
           treat as the chunk, or Bio::Root::IO.
           Optionally, also provide:
           -start => int : the byte position within the thing described by the
                           first argument to consider as the start of this
                           chunk (default 0)
           -end   => int : the byte position to consider as the end (default
                           true end)
           -piped_behaviour => 'memory'|'temp_file'|'sequential_read'

           The last option comes into effect when the first argument is
           something that cannot be seeked (eg. piped input filehandle).
            'memory'          means read all the piped input into a string
                              first, then set the chunk to that string.
            'temp_file'       means read all the piped input and output it to
                              a temp file, then set the chunk to that temp file.
            'sequential_read' means that the piped input should be read
                              sequentially and your parsing code must cope with
                              not being able to seek.
           'memory' is the fastest but uses the most memory. 'temp_file' and
           'sequential_read' can be slow, with 'temp_file' being the most memory
           efficient but requiring disc space. The default is 'sequential_read'.
           Note that in versions of perl earlier than 5.8 only temp_file works
           and will be used regardless of what value is supplied here.

=cut

sub chunk {
    my $self = shift;
    
    if (@_) {
        my $thing = shift || $self->throw("Trying to set chunk() to an undefined value");
        if (ref($thing) eq 'GLOB') {
            $self->{_chunk} = Bio::Root::IO->new(-fh => $thing);
        }
        elsif (ref(\$thing) eq 'SCALAR') {
            if ($thing !~ /\n/ && -e $thing) {
                $self->{_chunk} = Bio::Root::IO->new(-file => $thing);
            }
            else {
                unless ($FORCE_TEMP_FILE) {
                    # treat a string as a filehandle
                    open my $fake_fh, "+<", \$thing or $self->throw("Could not open file '$thing': $!"); # requires perl 5.8
                    $self->{_chunk} = Bio::Root::IO->new(-fh => $fake_fh);
                }
                else {
                    my ($handle) = $self->{_chunk}->tempfile();
                    print $handle $thing;
                    $self->{_chunk} = Bio::Root::IO->new(-fh => $handle);
                }
            }
        }
        elsif ($thing->isa('Bio::Root::IO')) {
            $self->{_chunk} = $thing;
        }
        else {
            $self->throw("Unknown input into chunk()");
        }
        
        my ($piped_behaviour, $start, $end);
        if (@_) {
            ($piped_behaviour, $start, $end) =
                $self->_rearrange([qw(PIPED_BEHAVIOUR START END)], @_);
        }
        $piped_behaviour ||= 'sequential_read';
        $FORCE_TEMP_FILE && ($piped_behaviour = 'temp_file');
        $start ||= 0;
        $self->_chunk_true_start($start);
        $self->_chunk_true_end($end);
        
        # determine if the chunk is seekable
        my $fh = $self->{_chunk}->_fh;
        seek($fh, 0, 0);
        my $first_line = <$fh>;
        seek($fh, 0, 0);
        my $seekable = tell($fh) == 0;
        unless ($seekable) {
            if ($piped_behaviour eq 'memory') {
                my $string = $first_line;
                while (<$fh>) {
                    $string .= $_;
                }
                $self->chunk($string);
            }
            elsif ($piped_behaviour eq 'temp_file') {
                my ($handle) = $self->{_chunk}->tempfile();
                print $handle $first_line;
                while (<$fh>) {
                    print $handle $_;
                }
                seek($handle, 0, 0);
                $self->chunk($handle);
            }
            elsif ($piped_behaviour eq 'sequential_read') {
                $self->{_chunk}->_pushback($first_line);
                $self->_sequential(1);
            }
            else {
                $self->throw("Unknown piped behaviour type '$piped_behaviour'");
            }
        }
        
        # determine our line ending
        if ($first_line =~ /\r\n/) {
            $self->_line_ending("\r\n");
        }
        elsif ($first_line =~ /\r/) {
            $self->_line_ending("\r");
        }
        else {
            $self->_line_ending("\n");
        }
    }
    
    return $self->{_chunk} || return;
}

=head2 _sequential

 Title   : _sequential
 Usage   : if ($obj->_sequential) {...}
 Function: Ask if we have to do operations such that the input is read
           sequentially.
 Returns : boolean
 Args    : none to get, OR boolean to set (typically, you should never set this
           yourself)

=cut

sub _sequential {
    my $self = shift;
    if (@_) {
        $self->{_sequential} = shift;
    }
    return $self->{_sequential} || 0;
}

=head2 _dependencies

 Title   : _dependencies
 Usage   : $obj->_dependencies( { field1 => field2 } );
           my $dependancy = $obj->_dependencies('field_name');
 Function: Set the fields that are dependent on each other, or get the field
           than another is dependent upon.
 Returns : string (a field name)
 Args    : string (a field name) to get, OR hash ref to initially set, with
           field names as keys and values, key field being dependent upon value
           field.

=cut

sub _dependencies {
    my ($self, $thing) = @_;
    $thing || return;
    if (ref($thing) eq 'HASH') {
        $self->{_dependencies} = $thing;
    }
    else {
        return $self->{_dependencies}->{$thing};
    }
}

=head2 _chunk_true_start

 Title   : _chunk_true_start
 Usage   : my $true_start = $obj->_chunk_true_start;
 Function: Get/set the true start position of the chunk within the filehandle
           it is part of.
 Returns : int
 Args    : none to get, OR int to set (typically, you won't set this yourself)

=cut

sub _chunk_true_start {
    my $self = shift;
    if (@_) {
        $self->{_chunk_start} = shift;
    }
    return $self->{_chunk_start} || 0;
}

=head2 _chunk_true_end

 Title   : _chunk_true_end
 Usage   : my $true_end = $obj->_chunk_true_end;
 Function: Get/set for the true end position of the chunk within the filehandle
           it is part of.
 Returns : int
 Args    : none to get, OR int to set (typically, you won't set this yourself)

=cut

sub _chunk_true_end {
    my $self = shift;
    if (@_) {
        $self->{_chunk_end} = shift;
    }
    return $self->{_chunk_end};
}

=head2 _line_ending

 Title   : _line_ending
 Usage   : my $line_ending = $obj->_line_ending;
 Function: Get/set for the line ending for the chunk.
 Returns : string
 Args    : none to get, OR string to set (typically, you won't set this
           yourself)

=cut

sub _line_ending {
    my $self = shift;
    if (@_) {
        $self->{_chunk_line_ending} = shift;
    }
    return $self->{_chunk_line_ending};
}

=head2 _chunk_seek

 Title   : _chunk_seek
 Usage   : $obj->_chunk_seek($pos);
 Function: seek() the chunk to the provided position in bytes, relative to the
           defined start of the chunk within its filehandle.

           In _sequential() mode, this function does nothing.

 Returns : n/a
 Args    : int

=cut

sub _chunk_seek {
    my ($self, $pos) = @_;
    $self->throw("Undefined position passed") unless defined $pos;
    return if $self->_sequential;
    
    my $fh = $self->chunk->_fh;
    
    # seek to the defined start
    seek($fh, $self->_chunk_true_start, 0);
    
    # now seek to desired position relative to defined start
    seek($fh, $pos, 1);
}

=head2 _chunk_tell

 Title   : _chunk_seek
 Usage   : my $pos = $obj->_chunk_tell;
 Function: Get the current tell() position within the chunk, relative to the
           defined start of the chunk within its filehandle.

           In _sequential() mode, this function does nothing.

 Returns : int
 Args    : none

=cut

sub _chunk_tell {
    my $self = shift;
    return if $self->_sequential;
    
    my $fh = $self->chunk->_fh;
    return tell($fh) - $self->_chunk_true_start;
}

=head2 _get_chunk_by_nol

 Title   : _chunk_seek
 Usage   : my $string = $obj->_get_chunk_by_nol;
 Function: Get a chunk of chunk() from the current position onward for the given
           number of lines.
 Returns : string
 Args    : int (number of lines you want)

=cut

sub _get_chunk_by_nol {
    my ($self, $nol) = @_;
    $nol > 0 || $self->throw("Can't request a chunk of fewer than 1 lines");
    
    # hope that $/ is \n
    
    my ($line, $count);
    while (defined($_ = $self->chunk->_readline)) {
        $line .= $_;
        $count++;
        last if $count == $nol;
    }
    
    my $current = $self->_chunk_tell;
    my $end = ($current || 0) + $self->_chunk_true_start;
    if (! $current || ($self->_chunk_true_end ? $end <= $self->_chunk_true_end : 1)) {
        return $line;
    }
    return;
}

=head2 _get_chunk_by_end

 Title   : _get_chunk_by_end
 Usage   : my $string = $obj->_get_chunk_by_end;
 Function: Get a chunk of chunk() from the current position onward till the end
           of the line, as defined by the supplied argument.
 Returns : string
 Args    : string (line ending - if you want the line ending to include a new
           line, always use \n)

=cut

sub _get_chunk_by_end {
    my ($self, $chunk_ending) = @_;
    
    my $start = $self->_chunk_tell;
    
    my $line_ending = $self->_line_ending;
    $chunk_ending =~ s/\n/$line_ending/g;
    local $/ = $chunk_ending || '';
    my $line = $self->chunk->_readline;
    
    my $current = $self->_chunk_tell;
    my $end = ($current || 0) + $self->_chunk_true_start;
    if (! $current || ($self->_chunk_true_end ? $end <= $self->_chunk_true_end : 1)) {
        return $line;
    }
    
    $self->_chunk_seek($start);
    return;
}

=head2 _find_chunk_by_end

 Title   : _find_chunk_by_end
 Usage   : my $string = $obj->_find_chunk_by_end;
 Function: Get the start and end of what would be a chunk of chunk() from the
           current position onward till the end of the line, as defined by the
           supplied argument.

           In _sequential() mode, this function does nothing.

 Returns : _chunk_tell values for start and end in 2 element list
 Args    : string (line ending - if you want the line ending to include a new
           line, always use \n)

=cut

sub _find_chunk_by_end {
    my ($self, $chunk_ending) = @_;
    return if $self->_sequential;
    
    my $line_ending = $self->_line_ending;
    $chunk_ending =~ s/\n/$line_ending/g;
    local $/ = $chunk_ending || '';
    
    my $start = $self->_chunk_tell;
    $self->chunk->_readline;
    my $end = $self->_chunk_tell;
    
    my $comp_end = $end + $self->_chunk_true_start;
    if ($self->_chunk_true_end ? $comp_end <= $self->_chunk_true_end : 1) {
        return ($start, $end);
    }
    
    $self->_chunk_seek($start);
    return;
}

=head2 AUTOLOAD

 Title   : AUTOLOAD
 Usage   : n/a
 Function: Assumes that any unknown method called should be treated as
           get_field($method_name).
 Returns : n/a
 Args    : n/a

=cut

sub AUTOLOAD {
    my $self = shift;
    ref($self) || return;
    
	my $name = $AUTOLOAD;
	$name =~ s/.*://; # strip fully-qualified portion
    
    # is it one of our fields?
    return $self->get_field($name);
}

1;
