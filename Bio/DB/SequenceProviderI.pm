package Bio::DB::SequenceProviderI;

# $Id $
# A provider of Bio::PrimarySeqI objects.

=head1 NAME

Bio::DB::SequenceProviderI - A provider of L<Bio::PrimarySeqI> objects.

=head1 SYNOPSIS

  # First we need a provider object, such as Bio::DB::SimpleSequenceProvider.
  use Bio::DB::SimpleSequenceProvider;
  # Create it and add the (given) Bio::PrimarySeqI object $roa1_human_seq.
  my $provider = new Bio::DB::SimpleSequenceProvider( $roa1_human_seq );
  # Now we can use the provider to retrieve it by its unique_id.
  my $seq = $db->get_Seq_by_id( 'ROA1_HUMAN' );
  # We could also have retrieved it by its accession.

=head1 DESCRIPTION

This is a pure interface class - in other words, all this does is define
methods which other (concrete) classes will actually implement.

The Bio::DB::SequenceProviderI class defines what methods a generic
database class should have.  At the moment it is just the ability to
make L<Bio::PrimarySeqI> objects from a name (id) or an accession number.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR

Ewan Birney originally wrote this class, as Bio::DB::RandomAccessI.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=head1 CONTRIBUTORS

Paul Edlefsen E<lt>paul@systemsbiology.orgE<gt>.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

use vars qw(@ISA);
use strict;

use Bio::Root::RootI;

## TODO: For now we're leaving Bio::DB::RandomAccessI and
## Bio::DB::SeqI.  Later we should replace them with this (IMO).  For
## now we'll subclass them.
#@ISA = qw( Bio::Root::RootI );
use Bio::DB::SeqI;
# TODO: Also integrate updating methods (from Bio::DB::UpdateableSeqI)
@ISA = qw( Bio::Root::RootI Bio::DB::SeqI );

## Note: This is to support Bio::DB::SeqI and might be removeable
## later on.  A grep turns up no usage references (there are a few
## implementations though, just no calls to this method).
=head2 get_all_primary_ids

 Title   : get_all_primary_ids
 Usage   : my @primary_ids = $provider->get_all_primary_ids()
 Function: Returns an array of all the primary_ids of the 
           sequence objects in the database. These
           maybe ids (display style) or accession numbers
           or something else completely different - they
           *are not* meaningful outside of this database
           implementation.
 Returns : an array of strings
 Args    : none
 Status  : Public

=cut

sub get_all_primary_ids {
  shift->throw_not_implemented();
}

=head2 sequences

 Title   : sequences
 Usage   : my @seqs = $provider->sequences( @names );
           OR
           my @seqs = $provider->sequences( %args );
 Function: Retrieves a list of L<Bio::PrimarySeqI> objects.
 Returns : a list of L<Bio::PrimarySeqI> objects
           OR
           (when the -iterator option is true)
             a L<Bio::Factory::SequenceStreamI> object
           OR
           (when the -callback option is true) true iff the callbacks
             completed.
 Args    : see below
 Status  : Public
 Throws  : "$some_text does not exist" ($some_text might be anything)
            if a version is given and a matching sequence exists,
            but not of that version.

Two types of argument lists are accepted.  In the positional argument
form, the arguments are treated as a list of sequence names (as if they
were given as -names => \@_).  In the named parameter form, the
arguments are a series of -name=E<gt>value pairs.  Note that the table
below is not exhaustive; implementations must support these but may
support other arguments as well (and are responsible for documenting the
difference).

  Argument       Description
  --------       ------------

  -unique_id     A (string) unique_id.  See also -namespace and -id.
  -unique_ids    An array reference to multiple unique_id values.

  -primary_id    A (string) primary_id.  See also -namespace and -id.
  -primary_ids   An array reference to multiple primary_id values.

  -display_id    A (string) display_id.  See also -namespace and -id.
  -display_ids   An array reference to multiple display_id values.

  -id            A (string) unique_id or primary_id or display_id.
                 See also -namespace.
  -ids           An array reference to multiple id values.

  -accession     A (string) accession.  See also -namespace and -name.
  -accessions    An array reference to multiple accession values.

  -name          A (string) accession or id.  See also -namespace and -id.
  -names         An array reference to multiple accession/id values.

  -namespace     A (string) namespace qualifier to help resolve the names.
  -class         same as -namespace

  -iterator      Return a L<Bio::Factory::SequenceStreamI>
  -stream        same as -iterator

  -callback      A callback to invoke on each sequence

All plural arguments are interchangeable with their singular counterparts.

The -unique_ids argument is a reference to a list of strings.  Every
returned sequence must have its unique_id value in this list if the
list is provided*.  Iff a sequence has no unique_id value then its
primary_id will be used instead.  A -unique_id argument is treated as
a single-element list of unique_ids.

The -ids argument is a reference to a list of strings.  Every returned
sequence must have its unique_id, primary_id, or display_id value in
this list if the list is provided*.  An -id argument is treated as a
single-element list of ids.

The -accessions argument is a reference to a list of strings.  Every
returned sequence must have its accession value in this list if the
list is provided*.  An -accession argument is treated as a
single-element list of names.  If the accession value contains a dot
('.') then the suffix after the final dot will be interpreted as a
version number.  If the accession is available but not in the given
version then an exception (ending in ' does not exist') will be
thrown.  An empty version value is acceptable and means 'latest version',
which is also the default, so -accession => 'foo' and -accession =>
'foo.' mean the same thing.

NOTE: If your accession value contains a dot (unrelated to version
number) then you should postpend the value with a dot to disambiguate:
'accession.that.contains.dots' will be interpreted as the accession
'accession.that.contains' with version 'dots' unless you postpend it
with a dot, like so: -accession => 'accession.that.contains.dots.'

The -names argument is a reference to a list of strings.  Every
returned sequence must have its accession, unique_id, primary_id,
or display_id value in this list if the list is provided*.  A -name
argument is treated as a single-element list of names.

If a -namespace is provided then names and ids (both queries and
targets) will be prepended with "$namespace:" as a bonus.  So if you
do sequences( -names => [ 'foo', 'bar' ], -namespace => 'ns' ) then any
sequence with the accession, unique_id, primary_id, or display_id
'foo', 'ns:foo', 'bar', or 'ns:bar' will be returned.

If -iterator is true, then the method returns an object of type
L<Bio::Factory::SequenceStreamI>.  Each call to next_seq() on this
object returns a L<Bio::PrimarySeqI> object from this provider.

If -callback is passed a code reference, the code reference will be
invoked on each sequence returned.  The code will be passed two
arguments consisting of the current sequence and this SequenceProviderI
object, and must return a true value.  If the code returns a false
value, sequence retrieval will be aborted.

-callback and -iterator are mutually exclusive options.  If -iterator
is defined, then -callback is ignored.

Footnote (*): -unique_ids, -primary_ids, -display_ids, -ids,
-accessions, and -names will be ORed together if they appear together,
so if you do $provider->sequences( -unique_id => 'foo', -accession =>
'bar' ) then you will get all sequences with unique_id 'foo' and
(also) all sequences with accession 'bar'.

NOTE: the following methods all build on top of sequences(), and do not
need to be explicitly implemented.

    get_Seq_by_id()
    get_Seq_by_primary_id()
    get_Seq_by_accession()
    get_Seq_by_version()
    get_Seq_stream()
    get_PrimarySeq_stream()

=cut

sub sequences {
  shift->throw_not_implemented();
}

=head2 get_Seq_by_id

 Title   : get_Seq_by_id
 Usage   : my $seq = $provider->get_Seq_by_id( $id, %args );
           OR
           my @seqs = $provider->get_Seq_by_id( $id, %args );
 Function: Gets a L<Bio::PrimarySeqI> object or list of objects by id.
 Returns : a L<Bio::PrimarySeqI> object (in scalar context)
           OR
           a list of L<Bio::PrimarySeqI> objects (in list context)
           OR
           undef if there are none with that id (in either context)
 Args    : the unique_id or primary_id or display_id (string) of a sequence,
           or a reference to a list thereof.
 Status  : Public
 Throws  : "more than one sequence corresponds" if there are multiple sequences
            to return and this method is called in a scalar context.
           "$some_text does not exist" ($some_text might be anything)
            if a version is given and a matching sequence exists,
            but not of that version.

This method is identical to sequences() except that it takes the -id
(or -ids) value as its first argument.  Also, if called from a scalar
context it attempts to return a single sequence (and throws an
exception if it is unable to do so).

NOTE: This is defined in the interface in terms of sequences().  You do not
have to implement it.

=cut

sub get_Seq_by_id {
  my $self = shift;
  my $id = shift;
  if( $_[ 0 ] =~ /^-/ ) {
    if( wantarray ) {
      return $self->sequences( @_, -id => $id );
    } else {
      my @result_array = $self->sequences( @_, -id => $id );
      if( scalar( @result_array ) > 1 ) {
        $self->throw( "more than one sequence corresponds" );
      } elsif( scalar( @result_array ) == 1 ) {
        return shift @result_array;
      } else {
        return undef;
      }
    }
  } else {
    if( wantarray ) {
      return $self->sequences( -names => \@_, -id => $id );
    } else {
      my @result_array = $self->sequences( -names => \@_, -id => $id );
      if( scalar( @result_array ) > 1 ) {
        $self->throw( "more than one sequence corresponds" );
      } elsif( scalar( @result_array ) == 1 ) {
        return shift @result_array;
      } else {
        return undef;
      }
    }
  }
} # get_Seq_by_id(..)

=head2 get_Seq_by_primary_id

 Title   : get_Seq_by_primary_id
 Usage   : my $seq = $provider->get_Seq_by_primary_id( $id, %args );
           OR
           my @seqs = $provider->get_Seq_by_primary_id( $id, %args );
 Function: Gets a L<Bio::PrimarySeqI> object or list of objects by primary_id.
 Returns : a L<Bio::PrimarySeqI> object (in scalar context)
           OR
           a list of L<Bio::PrimarySeqI> objects (in list context)
           OR
           undef if there are none with that primary id (in either context)
 Args    : the primary_id (string) of a sequence, or a reference to a list
           thereof.
 Status  : Public
 Throws  : "more than one sequence corresponds" if there are multiple sequences
            to return and this method is called in a scalar context.
           "$some_text does not exist" ($some_text might be anything)
            if a version is given and a matching sequence exists,
            but not of that version.

This method is identical to sequences() except that it takes the -primary_id
(or -primary_ids) value as its first argument.  Also, if called from a scalar
context it attempts to return a single sequence (and throws an
exception if it is unable to do so).

NOTE: This is defined in the interface in terms of sequences().  You do not
have to implement it.

=cut

sub get_Seq_by_primary_id {
  my $self = shift;
  my $id = shift;
  if( $_[ 0 ] =~ /^-/ ) {
    if( wantarray ) {
      return $self->sequences( @_, -primary_id => $id );
    } else {
      my @result_array = $self->sequences( @_, -primary_id => $id );
      if( scalar( @result_array ) > 1 ) {
        $self->throw( "more than one sequence corresponds" );
      } elsif( scalar( @result_array ) == 1 ) {
        return shift @result_array;
      } else {
        return undef;
      }
    }
  } else {
    if( wantarray ) {
      return $self->sequences( -names => \@_, -id => $id );
    } else {
      my @result_array = $self->sequences( -names => \@_, -id => $id );
      if( scalar( @result_array ) > 1 ) {
        $self->throw( "more than one sequence corresponds" );
      } elsif( scalar( @result_array ) == 1 ) {
        return shift @result_array;
      } else {
        return undef;
      }
    }
  }
} # get_Seq_by_primary_id(..)

=head2 get_Seq_by_acc

 Title   : get_Seq_by_acc
 Usage   : my $seq = $provider->get_Seq_by_acc( $accession_string );
           OR
           my $seq = $provider->get_Seq_by_acc( $namespace, $accession [, $version] );
           OR
           my @seqs = $provider->get_Seq_by_acc( $accession_string );
           OR
           my @seqs = $provider->get_Seq_by_acc( $namespace, $accession [, $version] );
 Function: Gets a L<Bio::PrimarySeqI> object or list of objects by accession.
 Returns : a L<Bio::PrimarySeqI> object (in scalar context)
           OR
           a list of L<Bio::PrimarySeqI> objects (in list context)
           OR
           undef if there are none with that id (in either context)
 Args    : accession (string), or a two- or three-
               element list consisting of namespace, accession, and version.
 Status  : Public
 Throws  : "more than one sequence corresponds" if there are multiple sequences
            to return and this method is called in a scalar context.
           "$some_text does not exist" ($some_text might be anything)
            if a version is given and a matching sequence exists,
            but not of that version.

This method is identical to sequences() except that it takes either
the -accession (or -accessions) value as its sole argument or it takes
3 arguments: the values for -namespace, -accessions, and a special
version value that will get postpended to the accession (after a dot;
see the sequences() docs for more about this).  Also, if called from a
scalar context it attempts to return a single sequence (and throws an
exception if it is unable to do so).

This method differs from get_Seq_by_version only in the following
respect: if the $version argument is undef, a dot will be postpended
anyway, just in case (with the assumption that the $accession value
does not contain version information).

NOTE: This is defined in the interface in terms of sequences().  You do not
have to implement it.

=cut

sub get_Seq_by_acc {
  my $self = shift;
  my $accession = shift;

  my ( $namespace, $version );
  if( @_ ) {
    $namespace = $accession;
    ( $accession, $version ) = @_;
  }

  unless( $accession ) {
    @self->throw( "Bio::DB::SequenceProviderI->get_Seq_by_acc(..) requires either 1 argument (an accession value or a reference to a non-empty list of accession values) or 2 or 3 arguments (namespace, accession(s), and optionally a version).  Somehow it didn't get an $accession." );
  }

  if( ref $accession ) {
    map { $_ .= ".$version" } @$accession;
  } else {
    $accession .= ".$version";
  }
  my @result_array =
    $self->sequences(
      -accession => $accession,
      -namespace => $namespace
    );
  if( wantarray ) {
    return @result_array;
  } else {
    if( scalar( @result_array ) > 1 ) {
      $self->throw( "more than one sequence corresponds" );
    } elsif( scalar( @result_array ) == 1 ) {
      return shift @result_array;
    } else {
      return undef;
    }
  }
} # get_Seq_by_acc(..)

=head2 get_Seq_by_version

 Title   : get_Seq_by_version
 Usage   : my $seq = $provider->get_Seq_by_version( $accession_string );
           OR
           my $seq = $provider->get_Seq_by_version( $namespace, $accession [, $version] );
           OR
           my @seqs = $provider->get_Seq_by_version( $accession_string );
           OR
           my @seqs = $provider->get_Seq_by_version( $namespace, $accession [, $version] );
 Function: Gets a L<Bio::PrimarySeqI> object or list of objects by accession.
 Returns : a L<Bio::PrimarySeqI> object (in scalar context)
           OR
           a list of L<Bio::PrimarySeqI> objects (in list context)
           OR
           undef if there are none with that id (in either context)
 Args    : versioned accession (string), or a two- or three-
               element list consisting of namespace, accession, and version.
 Status  : Public
 Throws  : "more than one sequence corresponds" if there are multiple sequences
            to return and this method is called in a scalar context.
           "$some_text does not exist" ($some_text might be anything)
            if a version is given and a matching sequence exists,
            but not of that version.

This method is identical to sequences() except that it takes either
the -accession (or -accessions) value as its sole argument or it takes
3 arguments: the values for -namespace, -accessions, and a special
version value that will get postpended to the accession (after a dot;
see the sequences() docs for more about this).  Also, if called from a
scalar context it attempts to return a single sequence (and throws an
exception if it is unable to do so).

This method differs from get_Seq_by_acc(..) only in the following
respect: if the $version argument is undef, no dot will be postpended,
with the assumption that the given $accession already contains a
version.

NOTE: This is defined in the interface in terms of sequences().  You do not
have to implement it.

=cut

sub get_Seq_by_version {
  my $self = shift;
  my $accession = shift;

  my ( $namespace, $version );
  if( @_ ) {
    $namespace = $accession;
    ( $accession, $version ) = @_;
  }
  unless( $accession ) {
    @self->throw( "Bio::DB::SequenceProviderI->get_Seq_by_version(..) requires either 1 argument (an accession value or a reference to a non-empty list of accession values) or 2 or 3 arguments (namespace, accession(s), and optionally a version).  Somehow it didn't get an $accession." );
  }

  # NOTE: This if.. is the only diff b/n this method and get_Seq_by_acc(..).
  #       In the other method the contents are always executed (there's no if).
  if( $version ) {
    if( ref $accession ) {
      map { $_ .= ".$version" } @$accession;
    } else {
      $accession .= ".$version";
    }
  }
  my @result_array =
    $self->sequences(
      -accession => $accession,
      -namespace => $namespace
    );
  if( wantarray ) {
    return @result_array;
  } else {
    if( scalar( @result_array ) > 1 ) {
      $self->throw( "more than one sequence corresponds" );
    } elsif( scalar( @result_array ) == 1 ) {
      return shift @result_array;
    } else {
      return undef;
    }
  }
} # get_Seq_by_version(..)

=head2 get_Seq_stream

 Title   : get_Seq_stream
 Usage   : my $iterator = $provider->get_Seq_stream( %args )
 Function: get an iterator over the sequences provided by this provider
 Returns : a L<Bio::Factory::SequenceStreamI>
 Args    : same as sequences()
 Status  : Public

This method is identical to sequences() except that it always generates
an iterator.

NOTE: This is defined in the interface in terms of sequences().  You do not
have to implement it.

=cut

sub get_Seq_stream {
  my $self = shift;
  if( $_[ 0 ] =~ /^-/ ) {
    return $self->sequences( @_, -iterator => 1 );
  } else {
    return $self->sequences( -names => \@_, -iterator => 1 );
  }
} # get_Seq_stream()

=head2 get_PrimarySeq_stream [deprecated]

 Title   : get_PrimarySeq_stream [deprecated]
 Usage   : my $iterator = $provider->get_PrimarySeq_stream( %args )
 Function: get an iterator over the sequences provided by this provider
 Returns : a L<Bio::Factory::SequenceStreamI>
 Args    : same as sequences()
 Status  : Public

This method is identical to sequences() except that it always generates
an iterator.

NOTE: This is a (glob ref) alias to get_Seq_stream().  It is
deprecated.  get_Seq_stream() should be used instead.

=cut

  *get_PrimarySeq_stream = \&get_Seq_stream;

1;

__END__

