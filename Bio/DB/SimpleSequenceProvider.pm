package Bio::DB::SimpleSequenceProvider;

## TODO: Finish this.  This is a hacky start.  See the other TODOs.

# $Id$
## A simple implementation of the Bio::DB::SequenceProviderI
## interface.

=head1 NAME

Bio::DB::SimpleSequenceProvider - An in-memory implementation of the
SequenceProviderI interface.

=head1 SYNOPSIS

=head1 DESCRIPTION

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

Paul Edlefsen E<lt>paul@systemsbiology.orgE<gt>.

Copyright (c) 2003 Institute for Systems Biology

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=head1 CONTRIBUTORS

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

use strict;
use vars qw( $VERSION @ISA );
use overload 
  '""' => 'toString',
  cmp   => '_cmp';

use Bio::DB::GFF::Util::Rearrange; # for 'rearrange'
use Bio::RelRange qw( &absSeqId );

use Bio::Root::Root;
use Bio::DB::SequenceProviderI;

$VERSION = '0.01';
@ISA = qw( Bio::Root::Root Bio::DB::SequenceProviderI );

=head2 new

 Title   : new
 Usage   : my $obj =
             new Bio::DB::SimpleSequenceProvider( @sequences );
 Function: Builds a new Bio::DB::SimpleSequenceProvider object 
 Returns : a new Bio::DB::SimpleSequenceProvider
 Args    : L<Bio::PrimarySeqI> objects to store herein

=cut

sub new {
  my( $class, @args ) = @_;

  my $self = $class->SUPER::new( @args );
  $self->_initialize_simple_sequence_provider( @args );
  return $self;
} # new(..)

sub _initialize_simple_sequence_provider {
  my $self = shift;
  my @args = @_;

  return if( $self->{ '_simple_sequence_provider_initialized' } );
  $self->{ '_unique_id_to_sequence_table' }  = {};
  $self->{ '_primary_id_to_sequence_table' } = {};
  $self->{ '_display_id_to_sequence_table' } = {};
  $self->{ '_accession_to_sequence_table' }  = {};

  my ( $sequences );
  if( scalar( @args ) && ( $args[ 0 ] =~ /^-/ ) ) {
    ( $sequences ) =
      rearrange(
        [
         [ qw( SEQUENCES SEQUENCE SEQS SEQ ) ],
        ],
        @args
      );
  } else {
    # sequences
    $sequences = \@args;
  }
  ## Fix up sequences.
  if( $sequences ) {
    unless( ref $sequences eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $sequences = [ $sequences ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$sequences ) ||
       ( ( scalar( @$sequences ) == 1 ) && !( $sequences->[ 0 ] ) )
      ) {
      undef $sequences;
    }
  }

  if( $sequences ) {
    foreach my $seq ( @$sequences ) {
      ## TODO: REMOVE
      use Data::Dumper;
      $self->throw( "\$seq is $seq, a " . ref( $seq ) );
      next unless( ref( $seq ) && $seq->isa( "Bio::PrimarySeqI" ) );
      unless( $self->_insert_sequence( $seq ) ) {
        $self->throw( "duplicate sequence: $seq" );
      }
    }
  }

  $self->{ '_simple_sequence_provider_initialized' }++;
  return $self;
} # _initialize_simple_sequence_provider(..)

=head2 get_all_primary_ids

 Title   : get_all_primary_ids
 Usage   : my @primary_ids = $provider->get_all_primary_ids()
 Function: Returns an array of all the primary_ids of the sequence
           objects in this data store. These maybe ids (display style)
           or accession numbers or something else completely different
           - they may be anything so long as each sequence has a
           different one.  Note that although some sequences may have
           undefined primary_ids (bad!), the returned list will not
           include undef.
 Returns : an array of strings
 Args    : none
 Status  : Public

=cut

sub get_all_primary_ids {
  my $self = shift;
  return grep { $_ ne 'undef' }
         keys %{ $self->{ '_primary_id_to_sequence_table' } };
} # get_all_primary_ids(..)

=head2 unique_ids

 Title     : unique_ids
 Usage     : my @unique_ids = $sequenceprovider->unique_ids();
 Function  : Return a list of the unique_ids of all sequences provided
             by this SequenceProvider.  Note that although some
             sequences may have undefined unique_ids, the returned
             list will not include undef.
 Returns   : an array of strings.
 Args      : none
 Status    : Public

=cut

sub unique_ids {
  my $self = shift;
  return grep { $_ ne 'undef' }
         keys %{ $self->{ '_unique_id_to_sequence_table' } };
} # unique_ids(..)

=head2 add_sequences

 Title   : add_sequences
 Usage   : my @added = $collection->add_sequences( @sequence_list );
 Function: Adds the given sequences to this provider.
 Returns : The sequences added (or their count, in scalar context).
 Args    : An array of L<Bio::PrimarySeqI>s
 Status  : Public

=cut

sub add_sequences {
  my $self = shift;
  my @sequence_list = @_;

  my $count = 0;
  my @added;
  foreach my $sequence ( @sequence_list ) {
    if( $sequence = $self->_insert_sequence( $sequence ) ) {
      if( wantarray ) {
        push( @added, $sequence );
      } else {
        $count++;
      }
    }
  }

  return ( wantarray ? @added : $count );
} # add_sequences(..)

=head2 remove_sequences

 Title   : remove_sequences
 Usage   : my @removed = $collection->remove_sequences( @sequence_list )
 Function: Removes the requested sequences from this provider.
 Returns : The removed sequences (or their count, in scalar context)
 Args    : An array of L<Bio::PrimarySeqI>s or their ids (see below)
 Status  : Public

  If any argument is a string, it will be taken as either the
  unique_id or the primary_id or the accession of a sequence to be
  removed.  The return list will contain the actual removed sequence
  if there was a match.

=cut

sub remove_sequences {
  my $self = shift;
  my @sequence_list = @_;

  my $count = 0;
  my @removed;
  foreach my $sequence ( @sequence_list ) {
    # Special case, to allow the list to include unique_ids, etc.
    ## TODO: Right now this won't work correctly; if the $sequence
    ## value is for the primary_id then the _remove_sequence(..)
    ## method will remove it from that table but not from the
    ## unique_id table.
    if( ref \$sequence eq 'STRING' ) {
      $sequence =
        Bio::PrimarySeq->new(
          '-unique_id' => $sequence,
          '-primary_id' => $sequence,
          '-accession_number' => $sequence
        );
    }
    if( $sequence = $self->_remove_sequence( $sequence ) ) {
      if( wantarray ) {
        push( @removed, $sequence );
      } else {
        $count++;
      }
    }
  }
  return ( wantarray ? @removed : $count );
} # remove_sequences()

=head2 sequences

 Title   : sequences
 Usage   : my @seqs = $provider->sequences( @names );
           OR
           my @seqs = $provider->sequences( %args );
 Function: Retrieves a list of L<Bio::PrimarySeqI> objects.
 Returns : a list of L<Bio::PrimarySeqI> objects
           OR
           (when the -iterator option is true)
             a L<Bio::Seq::IteratorI> object
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

  -iterator      Return a L<Bio::Seq::IteratorI>
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
L<Bio::Seq::IteratorI>.  Each call to next_seq() on this
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
  my $self = shift;

  my ( $unique_ids, $primary_ids, $display_ids, $ids, $accessions, $names, $namespace, $iterator, $callback );
  if( scalar( @_ ) && $_[ 0 ] =~ /^-/ ) {
    ( $unique_ids, $primary_ids, $display_ids, $ids, $accessions, $names, $namespace, $iterator, $callback ) =
      rearrange(
        [ [ qw( UNIQUE_IDS UNIQUEIDS UNIQUE_ID UNIQUEID ) ],
          [ qw( PRIMARY_IDS PRIMARYIDS PRIMARY_ID PRIMARYID ) ],
          [ qw( DISPLAY_IDS DISPLAYIDS DISPLAY_ID DISPLAYID ) ],
          [ qw( IDS ID ) ],
          [ qw( ACCESSIONS ACCESSION ACCESSION_NUMBERS ACCESSION_NUMBER
                ACCESSIONNUMBERS ACCESSIONNUMBER ) ],
          [ qw( NAME NAMES ) ],
          [ qw( NAMESPACE NAME_SPACE CLASS ) ],
          [ qw( ITERATOR STREAM ) ],
          [ qw( CALLBACK CALL_BACK ) ]
        ],
        @_
      );
  } else {
    ## Names.
    $names = \@_;
  }

  ## Fix up unique_ids.
  if( $unique_ids ) {
    unless( ref $unique_ids eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $unique_ids = [ $unique_ids ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$unique_ids ) ||
       ( ( scalar( @$unique_ids ) == 1 ) && !( $unique_ids->[ 0 ] ) )
      ) {
      undef $unique_ids;
    }
  }

  ## Fix up primary_ids.
  if( $primary_ids ) {
    unless( ref $primary_ids eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $primary_ids = [ $primary_ids ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$primary_ids ) ||
       ( ( scalar( @$primary_ids ) == 1 ) && !( $primary_ids->[ 0 ] ) )
      ) {
      undef $primary_ids;
    }
  }

  ## Fix up display_ids.
  if( $display_ids ) {
    unless( ref $display_ids eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $display_ids = [ $display_ids ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$display_ids ) ||
       ( ( scalar( @$display_ids ) == 1 ) && !( $display_ids->[ 0 ] ) )
      ) {
      undef $display_ids;
    }
  }

  ## Fix up ids.
  if( $ids ) {
    unless( ref $ids eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $ids = [ $ids ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$ids ) ||
       ( ( scalar( @$ids ) == 1 ) && !( $ids->[ 0 ] ) )
      ) {
      undef $ids;
    }
  }

  ## Fix up accessions.
  if( $accessions ) {
    unless( ref $accessions eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $accessions = [ $accessions ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$accessions ) ||
       ( ( scalar( @$accessions ) == 1 ) && !( $accessions->[ 0 ] ) )
      ) {
      undef $accessions;
    }
  }

  ## Fix up names.
  if( $names ) {
    unless( ref $names eq 'ARRAY' ) {
      ## The incoming value might be a simple scalar instead of a list.
      $names = [ $names ];
    }
    ## Just in case it's an array ref to an empty string:
    if(
       !scalar( @$names ) ||
       ( ( scalar( @$names ) == 1 ) && !( $names->[ 0 ] ) )
      ) {
      undef $names;
    }
  }

  # names are just accessions OR ids
  if( $names ) {
    $accessions = [] unless defined( $accessions );
    push( @{ $accessions }, @$names );
    $ids = [] unless defined( $ids );
    push( @{ $ids }, @$names );
  }

  # ids are just unique_ids, primary_ids, or display_ids
  if( $ids ) {
    $unique_ids = [] unless defined( $unique_ids );
    push( @{ $unique_ids }, @$ids );
    $primary_ids = [] unless defined( $primary_ids );
    push( @{ $primary_ids }, @$ids );
    $display_ids = [] unless defined( $display_ids );
    push( @{ $display_ids }, @$ids );
  }

  ## -iterator and -callback are mutually exclusive.
  if( $iterator && $callback ) {
    $self->throw( "The -iterator and -callback options are mutually exclusive, and yet both have been specified.  Perhaps you could apply your callback method to each element returned by the iterator?" );
  }

  ## This may be horribly inefficient, but what can ya do?  That's why
  ## this is called SimpleSequenceProvider.  Make your own if you want
  ## something better.
  my @sequences_to_return;

  ## TODO: Right now we're not supporting namespaces or versions.  We should.

  # By unique_ids
  if( $unique_ids ) {
    foreach my $unique_id ( @$unique_ids ) {
      next unless $self->{ '_unique_id_to_sequence_table' }->{ $unique_id };
      if( defined $callback ) {
        unless( $callback->(
                  $self->{ '_unique_id_to_sequence_table' }->{ $unique_id },
                  $self
                ) ) {
          return 0; # return false when the callback fails.
        }
      } else {
        push( @sequences_to_return,
              $self->{ '_unique_id_to_sequence_table' }->{ $unique_id } );
      }
    }
  }

  # By primary_ids
  if( $primary_ids ) {
    foreach my $primary_id ( @$primary_ids ) {
      next unless $self->{ '_primary_id_to_sequence_table' }->{ $primary_id };
      if( defined $callback ) {
        unless( $callback->(
                  $self->{ '_primary_id_to_sequence_table' }->{ $primary_id },
                  $self
                ) ) {
          return 0; # return false when the callback fails.
        }
      } else {
        push( @sequences_to_return,
              $self->{ '_primary_id_to_sequence_table' }->{ $primary_id } );
      }
    }
  }

  # By display_ids
  if( $display_ids ) {
    foreach my $display_id ( @$display_ids ) {
      next unless $self->{ '_display_id_to_sequence_table' }->{ $display_id };
      foreach my $sequence (
        @{ $self->{ '_display_id_to_sequence_table' }->{ $display_id } }
      ) {
        if( defined $callback ) {
          unless( $callback->( $sequence, $self ) ) {
            return 0; # return false when the callback fails.
          }
        } else {
          push( @sequences_to_return, $sequence );
        }
      }
    }
  }

  # By accessions
  if( $accessions ) {
    foreach my $accession ( @$accessions ) {
      next unless $self->{ '_accession_to_sequence_table' }->{ $accession };
      foreach my $sequence (
        @{ $self->{ '_accession_to_sequence_table' }->{ $accession } }
      ) {
        if( defined $callback ) {
          unless( $callback->( $sequence, $self ) ) {
            return 0; # return false when the callback fails.
          }
        } else {
          push( @sequences_to_return, $sequence );
        }
      }
    }
  }

  # If they've not given anything specific then they want everything
  unless( $unique_ids || $primary_ids || $display_ids || $accessions ) {
    if( defined $callback ) {
      foreach my $sequence (
        map { ( $_ eq 'undef' ) ?
              @{ $self->{ '_unique_id_to_sequence_table' }->{ 'undef' } } :
              $self->{ '_unique_id_to_sequence_table' }->{ $_ }
            }
            keys %{ $self->{ '_unique_id_to_sequence_table' } }
      ) {
        unless( $callback->( $sequence, $self ) ) {
          return 0; # return false when the callback fails.
        }
      }
    } else {
      push( @sequences_to_return,
            map { ( $_ eq 'undef' ) ?
                  @{ $self->{ '_unique_id_to_sequence_table' }->{ 'undef' } } :
                  $self->{ '_unique_id_to_sequence_table' }->{ $_ }
                }
                keys %{ $self->{ '_unique_id_to_sequence_table' } }
          );
    }
  } # End if no ids or accessions have been given, then return all sequences.

  if( $iterator ) {
    return Bio::Seq::SimpleIterator->new( @sequences_to_return );
  } elsif( defined $callback ) {
    return 1;
  } else {
    return @sequences_to_return;
  }
} # sequences(..)

=head2 sequence_count

 Title   : sequence_count
 Usage   : $collection->sequence_count()
 Function: Return the number of L<Bio::PrimarySeqI> objects that would
           be returned by a call to sequences() with no arguments.
 Returns : integer representing the number of sequence objects
 Args    : None

=cut

sub sequence_count {
  my $self = shift;

  my $count = 0;
  map { ( $_ eq 'undef' ) ?
        ( $count +=
          @{ $self->{ '_unique_id_to_sequence_table' }->{ 'undef' } } ) :
        ( $count +=
          ( defined( $self->{ '_unique_id_to_sequence_table' }->{ $_ } ) ?
            1 :
            0
          )
        )
      }
      keys %{ $self->{ '_unique_id_to_sequence_table' } };
  return $count;
} # sequence_count()

#### TODO: _insert_sequence, etc., should probably take some more measures to ensure consistency of the tables.  For one thing, tests should precede modification.  For another, one of the unique tables (unique_id and primary_id are unique within this provider) may be sufficient to identify a sequence even when the other tables are unable to; this fact could be utilized to avoid throwing exceptions when updating a sequence that has, eg., changed its unique_id but not its primary_id.

=head2 _insert_sequence

 Title   : _insert_sequence
 Usage   : $self->_insert_sequence( $sequence )
 Function: Inserts the sequence in the store.
 Returns : The sequence that was inserted, or undef if it was already
           in the store.
 Args    : A L<Bio::PrimarySeqI> object
 Status  : Protected

=cut

sub _insert_sequence {
  my $self = shift;
  my ( $sequence ) = @_;

  unless( defined( $sequence ) ) {
    $self->throw( "\$simple_sequence_provider->_insert_sequence( undef ): \$sequence is undef!" );
  }
  unless( ref( $sequence ) && $sequence->isa( 'Bio::PrimarySeqI' ) ) {
    $self->throw( "\$simple_sequence_provider->_insert_sequence( $sequence ): \$sequence is not a Bio::PrimarySeqI!" );
  }

  # By unique_id
  my $unique_id = $sequence->unique_id();
  unless( defined $unique_id ) {
    $unique_id = 'undef';
  }
  # Unique ids are, like, unique n stuff.
  if( $unique_id eq 'undef' ) {
    # Let's wait till we're sure it's new before we push it onto the
    # 'undef' list in the unique_id_to_sequence_table.
  } elsif( $self->{ '_unique_id_to_sequence_table' }->{ $unique_id } ) {
    return undef;
  } else {
    $self->{ '_unique_id_to_sequence_table' }->{ $unique_id } = $sequence;
  }
  
  # By primary_id
  my $primary_id = $sequence->primary_id();
  unless( defined $primary_id ) {
    $primary_id = 'undef';
  }
  # Primary ids are unique within one provider..
  if( $primary_id eq 'undef' ) {
    # Let's wait till we're sure it's new before we push it onto the
    # 'undef' list in the unique_id_to_sequence_table.
  } elsif( $self->{ '_primary_id_to_sequence_table' }->{ $primary_id } ) {
    return undef;
  } else {
    $self->{ '_primary_id_to_sequence_table' }->{ $primary_id } = $sequence;
  }
  
  # By display_id
  my $accession = $sequence->accession_number();
  my $display_id = $sequence->display_id();
  unless( defined $display_id ) {
    $display_id = 'undef';
  }
  if( defined $self->{ '_display_id_to_sequence_table' }->{ $display_id } ) {
    foreach my $sequence_with_the_same_display_id
            ( @{ $self->{ '_display_id_to_sequence_table' }->{ $display_id } } )
    {
      if( defined $unique_id ) {
        unless( $sequence_with_the_same_display_id->unique_id() eq $unique_id ) {
          next;
        }
      }
      if( defined $primary_id ) {
        unless( $sequence_with_the_same_display_id->primary_id() eq $primary_id ) {
          next;
        }
      }
      if( defined $accession ) {
        unless( $sequence_with_the_same_display_id->accession_number() eq $accession ) {
          next;
        }
      }
      # Okay, then as far as we can tell they're the same..
      return undef;
    }
  }
  push( @{ $self->{ '_display_id_to_sequence_table' }->{ $display_id } },
        $sequence );

  # By accession
  unless( defined $accession ) {
    $accession = 'undef';
  }
  if( defined $self->{ '_accession_to_sequence_table' }->{ $accession } ) {
    foreach my $sequence_with_the_same_accession
            ( @{ $self->{ '_accession_to_sequence_table' }->{ $accession } } )
    {
      if( $sequence_with_the_same_accession eq $sequence ) {
        return undef;
      }
    }
  }
  push( @{ $self->{ '_accession_to_sequence_table' }->{ $accession } },
        $sequence );
  if( $unique_id eq 'undef' ) {
    push( @{ $self->{ '_unique_id_to_sequence_table' }->{ 'undef' } },
          $sequence );
  }
  if( $primary_id eq 'undef' ) {
    push( @{ $self->{ '_primary_id_to_sequence_table' }->{ 'undef' } },
          $sequence );
  }
  return $sequence;
} # _insert_sequence(..)

=head2 _update_sequence

 Title   : _update_sequence
 Usage   : $self->_update_sequence( $sequence )
 Function: Updates the sequence in the store.
 Returns : The sequence that was updated or undef if the sequence is not
           in the store.
 Args    : A L<Bio::PrimarySeqI> object
 Status  : Protected

=cut

sub _update_sequence {
  my $self = shift;
  my ( $sequence ) = @_;

  unless( defined( $sequence ) ) {
    $self->throw( "\$simple_sequence_provider->_update_sequence( undef ): \$sequence is undef!" );
  }
  unless( ref( $sequence ) && $sequence->isa( 'Bio::PrimarySeqI' ) ) {
    $self->throw( "\$simple_sequence_provider->_update_sequence( $sequence ): \$sequence is not a Bio::PrimarySeqI!" );
  }

  my $found_it;

  # By unique_id
  my $unique_id = $sequence->unique_id();
  unless( defined $unique_id ) {
    $unique_id = 'undef';
  }
  # Unique ids are, like, unique n stuff.
  if( $unique_id eq 'undef' ) {
    my $num_seqs = 0;
    if(
       defined( $self->{ '_unique_id_to_sequence_table' }->{ $unique_id } )
      ) {
      $num_seqs =
        scalar( @{ $self->{ '_unique_id_to_sequence_table' }->{ $unique_id } } );
    }
    $found_it = 0;
    for( my $i = 0; $i < $num_seqs; $i++ ) {
      unless( $self->{ '_unique_id_to_sequence_table' }->{ $unique_id }->[ $i ]
              eq $sequence ) {
        next;
      }
      $found_it = 1;
      $self->{ '_unique_id_to_sequence_table' }->{ $unique_id }->[ $i ] =
        $sequence;
      last;
    }
    unless( $found_it ) {
      return undef;
    }
  } elsif( $self->{ '_unique_id_to_sequence_table' }->{ $unique_id } ) {
    $self->{ '_unique_id_to_sequence_table' }->{ $unique_id } = $sequence;
  } else {
    return undef;
  }
  
  # By primary_id
  my $primary_id = $sequence->primary_id();
  unless( defined $primary_id ) {
    $primary_id = 'undef';
  }
  # Primary ids are unique within one provider..
  if( $primary_id eq 'undef' ) {
    my $num_seqs = 0;
    if(
       defined( $self->{ '_primary_id_to_sequence_table' }->{ $primary_id } )
      ) {
      $num_seqs =
        scalar( @{ $self->{ '_primary_id_to_sequence_table' }->{ $primary_id } } );
    }
    $found_it = 0;
    for( my $i = 0; $i < $num_seqs; $i++ ) {
      unless( $self->{ '_primary_id_to_sequence_table' }->{ $primary_id }->[ $i ]
              eq $sequence ) {
        next;
      }
      $found_it = 1;
      $self->{ '_primary_id_to_sequence_table' }->{ $primary_id }->[ $i ] =
        $sequence;
      last;
    }
    unless( $found_it ) {
      return undef;
    }
  } elsif( $self->{ '_primary_id_to_sequence_table' }->{ $primary_id } ) {
    $self->{ '_primary_id_to_sequence_table' }->{ $primary_id } = $sequence;
  } else {
    return undef;
  }
  
  # By display_id
  my $display_id = $sequence->display_id();
  unless( defined $display_id ) {
    $display_id = 'undef';
  }
  my $num_seqs = 0;
  if(
     defined( $self->{ '_display_id_to_sequence_table' }->{ $display_id } )
    ) {
    $num_seqs =
      scalar( @{ $self->{ '_display_id_to_sequence_table' }->{ $display_id } } );
  }
  $found_it = 0;
  for( my $i = 0; $i < $num_seqs; $i++ ) {
    unless( $self->{ '_display_id_to_sequence_table' }->{ $display_id }->[ $i ]
            eq $sequence ) {
      next;
    }
    $found_it = 1;
    $self->{ '_display_id_to_sequence_table' }->{ $display_id }->[ $i ] =
      $sequence;
    last;
  }
  unless( $found_it ) {
    return undef;
  }

  # By accession
  my $accession = $sequence->accession_number();
  unless( defined $accession ) {
    $accession = 'undef';
  }
  my $num_seqs = 0;
  if(
     defined( $self->{ '_accession_to_sequence_table' }->{ $accession } )
    ) {
    $num_seqs =
      scalar( @{ $self->{ '_accession_to_sequence_table' }->{ $accession } } );
  }
  $found_it = 0;
  for( my $i = 0; $i < $num_seqs; $i++ ) {
    unless( $self->{ '_accession_to_sequence_table' }->{ $accession }->[ $i ]
            eq $sequence ) {
      next;
    }
    $found_it = 1;
    $self->{ '_accession_to_sequence_table' }->{ $accession }->[ $i ] =
      $sequence;
    last;
  }
  unless( $found_it ) {
    return undef;
  }

  return $sequence;
} # _update_sequence(..)

=head2 _insert_or_update_sequence

 Title   : _insert_or_update_sequence
 Usage   : $self->_insert_or_update_sequence( $sequence )
 Function: Inserts or updates the sequence in the store.
 Returns : The sequence that was added or updated.
 Args    : A L<Bio::PrimarySeqI> object
 Status  : Protected

=cut

sub _insert_or_update_sequence {
  my $self = shift;
  my ( $sequence ) = @_;

  unless( defined( $sequence ) ) {
    $self->throw( "\$simple_sequence_provider->_update_sequence( undef ): \$sequence is undef!" );
  }
  unless( ref( $sequence ) && $sequence->isa( 'Bio::PrimarySeqI' ) ) {
    $self->throw( "\$simple_sequence_provider->_update_sequence( $sequence ): \$sequence is not a Bio::PrimarySeqI!" );
  }

  my $found_it;

  # By unique_id
  my $unique_id = $sequence->unique_id();
  unless( defined $unique_id ) {
    $unique_id = 'undef';
  }
  # Unique ids are, like, unique n stuff.
  if( $unique_id eq 'undef' ) {
    my $num_seqs = 0;
    if(
       defined( $self->{ '_unique_id_to_sequence_table' }->{ $unique_id } )
      ) {
      $num_seqs =
        scalar( @{ $self->{ '_unique_id_to_sequence_table' }->{ $unique_id } } );
    }
    $found_it = 0;
    for( my $i = 0; $i < $num_seqs; $i++ ) {
      unless( $self->{ '_unique_id_to_sequence_table' }->{ $unique_id }->[ $i ]
              eq $sequence ) {
        next;
      }
      $found_it = 1;
      $self->{ '_unique_id_to_sequence_table' }->{ $unique_id }->[ $i ] =
        $sequence;
      last;
    }
    unless( $found_it ) {
      # Add it.
      push( @{ $self->{ '_unique_id_to_sequence_table' }->{ $unique_id } },
            $sequence );
    }
  } else {
    $self->{ '_unique_id_to_sequence_table' }->{ $unique_id } = $sequence;
  }
  
  # By primary_id
  my $primary_id = $sequence->primary_id();
  unless( defined $primary_id ) {
    $primary_id = 'undef';
  }
  # Primary ids are unique within one provider..
  if( $primary_id eq 'undef' ) {
    my $num_seqs = 0;
    if(
       defined( $self->{ '_primary_id_to_sequence_table' }->{ $primary_id } )
      ) {
      $num_seqs =
        scalar( @{ $self->{ '_primary_id_to_sequence_table' }->{ $primary_id } } );
    }
    $found_it = 0;
    for( my $i = 0; $i < $num_seqs; $i++ ) {
      unless( $self->{ '_primary_id_to_sequence_table' }->{ $primary_id }->[ $i ]
              eq $sequence ) {
        next;
      }
      $found_it = 1;
      $self->{ '_primary_id_to_sequence_table' }->{ $primary_id }->[ $i ] =
        $sequence;
      last;
    }
    unless( $found_it ) {
      # Add it.
      push( @{ $self->{ '_primary_id_to_sequence_table' }->{ $primary_id } },
            $sequence );
    }
  } else {
    $self->{ '_primary_id_to_sequence_table' }->{ $primary_id } = $sequence;
  }
  
  # By display_id
  my $display_id = $sequence->display_id();
  unless( defined $display_id ) {
    $display_id = 'undef';
  }
  my $num_seqs = 0;
  if(
     defined( $self->{ '_display_id_to_sequence_table' }->{ $display_id } )
    ) {
    $num_seqs =
      scalar(
        @{ $self->{ '_display_id_to_sequence_table' }->{ $display_id } }
      );
  }
  $found_it = 0;
  for( my $i = 0; $i < $num_seqs; $i++ ) {
    unless( $self->{ '_display_id_to_sequence_table' }->{ $display_id }->[ $i ]
            eq $sequence ) {
      next;
    }
    $found_it = 1;
    $self->{ '_display_id_to_sequence_table' }->{ $display_id }->[ $i ] =
      $sequence;
    last;
  }
  unless( $found_it ) {
    # Add it.
    push( @{ $self->{ '_display_id_to_sequence_table' }->{ $display_id } },
          $sequence );
  }

  # By accession
  my $accession = $sequence->accession_number();
  unless( defined $accession ) {
    $accession = 'undef';
  }
  my $num_seqs = 0;
  if(
     defined( $self->{ '_accession_to_sequence_table' }->{ $accession } )
    ) {
    $num_seqs =
      scalar( @{ $self->{ '_accession_to_sequence_table' }->{ $accession } } );
  }
  $found_it = 0;
  for( my $i = 0; $i < $num_seqs; $i++ ) {
    unless( $self->{ '_accession_to_sequence_table' }->{ $accession }->[ $i ]
            eq $sequence ) {
      next;
    }
    $found_it = 1;
    $self->{ '_accession_to_sequence_table' }->{ $accession }->[ $i ] =
      $sequence;
    last;
  }
  unless( $found_it ) {
    # Add it.
    push( @{ $self->{ '_accession_to_sequence_table' }->{ $accession } },
          $sequence );
  }

  return $sequence;
} # _insert_or_update_sequence(..)

=head2 _remove_sequence

 Title   : _remove_sequence
 Usage   : $self->_remove_sequence( $sequence )
 Function: Removes the sequence from the store.
 Returns : The sequence that was removed or undef if the sequence was not
           in the store.
 Args    : A L<Bio::PrimarySeqI> object
 Status  : Protected

=cut

sub _remove_sequence {
  my $self = shift;
  my ( $sequence ) = @_;

  unless( defined( $sequence ) ) {
    $self->throw( "\$simple_sequence_provider->_remove_sequence( undef ): \$sequence is undef!" );
  }
  unless( ref( $sequence ) && $sequence->isa( 'Bio::PrimarySeqI' ) ) {
    $self->throw( "\$simple_sequence_provider->_remove_sequence( $sequence ): \$sequence is not a Bio::PrimarySeqI!" );
  }

  my $found_it;

  # By unique_id
  my $unique_id = $sequence->unique_id();
  unless( defined $unique_id ) {
    $unique_id = 'undef';
  }
  # Unique ids are, like, unique n stuff.
  if( $unique_id eq 'undef' ) {
    my $num_seqs = 0;
    if(
       defined( $self->{ '_unique_id_to_sequence_table' }->{ $unique_id } )
      ) {
      $num_seqs =
        scalar( @{ $self->{ '_unique_id_to_sequence_table' }->{ $unique_id } } );
    }
    $found_it = 0;
    for( my $i = 0; $i < $num_seqs; $i++ ) {
      unless( $self->{ '_unique_id_to_sequence_table' }->{ $unique_id }->[ $i ]
              eq $sequence ) {
        next;
      }
      $sequence =
        $self->{ '_unique_id_to_sequence_table' }->{ $unique_id }->[ $i ];
      $found_it = 1;
      splice( @{ $self->{ '_unique_id_to_sequence_table' }->{ $unique_id } },
              $i,
              1 );
      last;
    }
    unless( $found_it ) {
      return undef;
    }
  } elsif( $self->{ '_unique_id_to_sequence_table' }->{ $unique_id } ) {
    delete $self->{ '_unique_id_to_sequence_table' }->{ $unique_id };
  } else {
    return undef;
  }
  
  # By primary_id
  my $primary_id = $sequence->primary_id();
  unless( defined $primary_id ) {
    $primary_id = 'undef';
  }
  # Primary ids are unique within one provider..
  if( $primary_id eq 'undef' ) {
    my $num_seqs = 0;
    if(
       defined( $self->{ '_primary_id_to_sequence_table' }->{ $primary_id } )
      ) {
      $num_seqs =
        scalar( @{ $self->{ '_primary_id_to_sequence_table' }->{ $primary_id } } );
    }
    $found_it = 0;
    for( my $i = 0; $i < $num_seqs; $i++ ) {
      unless( $self->{ '_primary_id_to_sequence_table' }->{ $primary_id }->[ $i ]
              eq $sequence ) {
        next;
      }
      $sequence =
        $self->{ '_primary_id_to_sequence_table' }->{ $unique_id }->[ $i ];
      $found_it = 1;
      splice( @{ $self->{ '_primary_id_to_sequence_table' }->{ $primary_id } },
              $i,
              1 );
      last;
    }
    unless( $found_it ) {
      return undef;
    }
  } elsif( $self->{ '_primary_id_to_sequence_table' }->{ $primary_id } ) {
    delete $self->{ '_primary_id_to_sequence_table' }->{ $primary_id };
  } else {
    return undef;
  }

  # By display_id
  my $display_id = $sequence->display_id();
  unless( defined $display_id ) {
    $display_id = 'undef';
  }
  my $num_seqs = 0;
  if(
     defined( $self->{ '_display_id_to_sequence_table' }->{ $display_id } )
    ) {
    $num_seqs =
      scalar( @{ $self->{ '_display_id_to_sequence_table' }->{ $display_id } } );
  }
  $found_it = 0;
  for( my $i = 0; $i < $num_seqs; $i++ ) {
    unless( $self->{ '_display_id_to_sequence_table' }->{ $display_id }->[ $i ]
            eq $sequence ) {
      next;
    }
    $found_it = 1;
    $sequence =
      $self->{ '_display_id_to_sequence_table' }->{ $unique_id }->[ $i ];
    splice( @{ $self->{ '_display_id_to_sequence_table' }->{ $display_id } },
            $i,
            1 );
    last;
  }
  unless( $found_it ) {
    return undef;
  }

  # By accession
  my $accession = $sequence->accession_number();
  unless( defined $accession ) {
    $accession = 'undef';
  }
  my $num_seqs = 0;
  if(
     defined( $self->{ '_accession_to_sequence_table' }->{ $accession } )
    ) {
    $num_seqs =
      scalar( @{ $self->{ '_accession_to_sequence_table' }->{ $accession } } );
  }
  $found_it = 0;
  for( my $i = 0; $i < $num_seqs; $i++ ) {
    unless( $self->{ '_accession_to_sequence_table' }->{ $accession }->[ $i ]
            eq $sequence ) {
      next;
    }
    $found_it = 1;
    $sequence =
      $self->{ '_accession_to_sequence_table' }->{ $unique_id }->[ $i ];
    splice( @{ $self->{ '_accession_to_sequence_table' }->{ $accession } },
            $i,
            1 );
    last;
  }
  unless( $found_it ) {
    return undef;
  }

  return $sequence;
} # _remove_sequence(..)

## TODO: ERE I AM.  Above is what I've done so far.  need to do remove_sequence and insert_or_update_sequence.  _*_seq(..) are from UpdateableSeqI, and they are almost like ours except that they throw exceptions instead of returning undef.  We'll need to add an equals (override eq) to PrimarySeqI, etc, perhaps something like this:
#  if( defined $unique_id ) {
#    unless( $other_sequence->unique_id() eq $unique_id ) {
#      return 0;
#    }
#  }
#  if( defined $primary_id ) {
#    unless( $other_sequence->primary_id() eq $primary_id ) {
#      return 0;
#    }
#  }
#  if( defined $display_id ) {
#    unless( $other_sequence->display_id() eq $display_id ) {
#      return 0;
#    }
#  }
#  if( defined $accession ) {
#    unless( $other_sequence->accession_number() eq $accession ) {
#      return 0;
#    }
#  }
#  # Okay, then as far as we can tell they're the same..
#  return 1;

=head2 _add_seq

 Title   : _add_seq
 Usage   : _add_seq($seq)
 Function: Adds a new sequence
 Example : 
 Returns : will throw an exception if
           sequences accession number already exists
 Args    : a new seq object - should have an accession number

=cut

# Implements UpdateableSeqI
sub _add_seq {
  my $self = shift;
  my ( $sequence ) = @_;

  unless( $self->_add_sequence( $sequence ) ) {
    $self->throw( "Unable to add a sequence that has already been added or that is identical ( by eq ) to another that has been added." );
  }
}

=head2 _remove_seq

 Title   : _remove_seq
 Usage   : _remove_seq($seq)
 Function: Removes an existing sequence
 Example : 
 Returns : will throw an exception if
           sequence does not exists for the primary_id
 Args    : a seq object that was retrieved from Bio::DB::UpdateableSeqI

=cut

# Implements UpdateableSeqI
sub _remove_seq {
  my $self = shift;
  my ( $sequence ) = @_;

  unless( $self->_remove_sequence( $sequence ) ) {
    $self->throw( "Unable to remove a sequence that has not been added." );
  }
} # _remove_seq(..)

=head2 _update_seq

 Title   : _update_seq
 Usage   : _update_seq($seq)
 Function: Updates a sequence
 Example : 
 Returns : will throw an exception if
           sequence is out of sync from expected val.
 Args    : a seq object that was retrieved from Bio::DB::UpdateableSeqI

=cut

# Implements UpdateableSeqI
sub _update_seq {
  my $self = shift;
  my ( $sequence ) = @_;

  unless( $self->_update_sequence( $sequence ) ) {
    $self->throw( "Unable to update a sequence that has not been added, or that is no longer identical ( by eq ) to itself as-added." );
  }
} # _update_seq(..)

=head2 toString

 Title   : toString
 Usage   : $str_val = $collection->toString()
 Function: returns "A SimpleSequenceProvider.";
 Returns : a String
 Args    : None
 Status  : Public

  This method is a hack.

=cut

sub toString {
  my $self = shift;
  ## TODO: Something cool.
  return "A SimpleSequenceProvider.";
} # toString()

## method for overload for comparing two SimpleSequenceProvider
## objects.  Uses toString().
sub _cmp {
  my $self = shift;
  my ( $b, $reversed ) = @_;
  my $a = $self->toString();
  ( $a, $b ) = ( $b, $a ) if $reversed;
  return ( $a cmp $b );
}

1;

__END__
