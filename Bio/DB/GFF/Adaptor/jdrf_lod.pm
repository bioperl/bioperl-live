# $Id$

=head1 NAME

Bio::DB::GFF::Adaptor::jdrf_lod

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::DB::GFF::Adaptor::jdrf_lod;

use strict;

use DBI;
use Bio::DB::GFF;
use Bio::DB::GFF::Util::Rearrange; # for rearrange()
use Bio::DB::GFF::Util::Binning;
use DBI;
#use Bio::DB::GFF::Adaptor::dbi::caching_handle;
use vars qw($VERSION @ISA);

@ISA =  qw(Bio::DB::GFF);
$VERSION = '0.01';

use Bio::Range;

##############################################################################


=head2 new

 Title   : new
 Usage   : $db = Bio::DB::GFF->new(@args)
 Function: create a new adaptor
 Returns : a Bio::DB::GFF object
 Args    : see below
 Status  : Public

This is the constructor for the adaptor.  It is called automatically
by Bio::DB::GFF-E<gt>new.  In addition to arguments that are common among
all adaptors, the following class-specific arguments are recgonized:

  Argument       Description
  --------       -----------

  -dsn           the DBI data source, e.g. 'dbi:mysql:ens0040'

  -user          username for authentication

  -pass          the password for authentication

=cut

# Create a new Bio::DB::GFF::Adaptor::jdrf_lod object
sub new {
  my $class = shift;
  my ($features_db,$username,$auth,$other) = rearrange([
							[qw(FEATUREDB DB DSN)],
							[qw(USERNAME USER)],
							[qw(PASSWORD PASS)],
						       ],@_);

  $features_db  || $class->throw("new(): Provide a data source or DBI database");

  if (!ref($features_db)) {
    my $dsn = $features_db;
    my @args;
    push @args,$username if defined $username;
    push @args,$auth     if defined $auth;
    $features_db = #Bio::DB::GFF::Adaptor::dbi::caching_handle->new($dsn,@args)
      DBI->connect( $dsn, @args )
        || $class->throw("new(): Failed to connect to $dsn: "
                         . #Bio::DB::GFF::Adaptor::dbi::caching_handle->errstr);
                         DBI->errstr);
  } else {
    $features_db->isa('DBI::db') 
      || $class->throw("new(): $features_db is not a DBI handle");
  }

  # fill in object
  return bless {
		features_db => $features_db
	       },$class;
}

sub debug {
  my $self = shift;
  $self->features_db->debug(@_);
  $self->SUPER::debug(@_);
}

=head2 features_db

 Title   : features_db
 Usage   : $dbh = $db->features_db
 Function: get database handle
 Returns : a DBI handle
 Args    : none
 Status  : Public

=cut

sub features_db { shift->{features_db} }
sub dbh         { shift->{features_db} }

=head2 get_dna

 Title   : get_dna
 Usage   : $string = $db->get_dna($name,$start,$stop,$class)
 Function: get DNA string
 Returns : ''
 Args    : name, class, start and stop of desired segment
 Status  : Public

returns ''.

=cut

sub get_dna {
  return '';
}


=head2 get_abscoords

 Title   : get_abscoords
 Usage   : ($refseq,$refclass,$start,$stop,$strand) = $db->get_abscoords($name,$class)
 Function: get absolute coordinates for landmark
 Returns : an array ref -- see below
 Args    : name and class of desired landmark
 Status  : Public

This method performs the low-level resolution of a landmark into a
reference sequence and position.

The result is an array ref, each element of which is a five-element
list containing reference sequence name, class, start, stop and strand.

=cut

sub get_abscoords {
  my $self = shift;
  my ( $name, $class, $refseq )  = @_;
  my ( $chr, $start, $end ) =
    ( $name =~ /^(.+)\:(\d+)\,(\d+)$/ );
  my $strand = '+';
  if( $start > $end ) {
    ( $start, $end ) = ( $end, $start );
    $strand = '-';
  }

  ## TODO: REMOVE
  #warn "jdrf_lod:get_abscoords( $chr:$start-$end ), name is $name;";

  ## TODO: Dehackify
  if( ( defined( $chr ) && ( $chr !~ /^\d/ ) ) ||
      ( defined( $name ) && ( $name !~ /^\d/ ) ) ) {
    return;
  }
  ## Paul's note: I figured out how to do this by trial and error.  I
  ## can't say I understand *why* I had to do it like I did, but what
  ## can ya do?
  return [ [ $chr|$name, undef, 1, ( $end - $start + 1 ), $strand ] ];
}


=head2 get_features

 Title   : get_features
 Usage   : $db->get_features($search,$options,$callback)
 Function: get list of features for a region
 Returns : count of number of features retrieved
 Args    : see below
 Status  : protected

The first argument is a hash reference containing search criteria for
retrieving features.  It contains the following keys:

   rangetype One of "overlaps", "contains" or "contained_in".  Indicates
              the type of range query requested.

   refseq    ID of the landmark that establishes the absolute 
              coordinate system.

   refclass  Class of this landmark.  Can be ignored by implementations
              that don't recognize such distinctions.

   start     Start of the range, inclusive.

   stop      Stop of the range, inclusive.

   types     Array reference containing the list of annotation types
              to fetch from the database.  Each annotation type is an
              array reference consisting of [source,method].

The second argument is a hash reference containing certain options
that affect the way information is retrieved:

   sort_by_group
             A flag.  If true, means that the returned features should be
             sorted by the group that they're in.

   sparse    A flag.  If true, means that the expected density of the 
             features is such that it will be more efficient to search
             by type rather than by range.  If it is taking a long
             time to fetch features, give this a try.

   binsize   A true value will create a set of artificial features whose
             start and stop positions indicate bins of the given size, and
             whose scores are the number of features in the bin.  The
             class of the feature will be set to "bin", and its name to
             "method:source".  This is a handy way of generating histograms
             of feature density.

The third argument, the $callback, is a code reference to which
retrieved features are passed.  It is described in more detail below.

This routine is responsible for getting arrays of GFF data out of the
database and passing them to the callback subroutine.  The callback
does the work of constructing a Bio::DB::GFF::Feature object out of
that data.  The callback expects a list of 13 fields:

  $refseq      The reference sequence
  $start       feature start
  $stop        feature stop
  $source      feature source
  $method      feature method
  $score       feature score
  $strand      feature strand
  $phase       feature phase
  $groupclass  group class (may be undef)
  $groupname   group ID (may be undef)
  $tstart      target start for similarity hits (may be undef)
  $tstop       target stop for similarity hits (may be undef)
  $feature_id  A unique feature ID (may be undef)

These fields are in the same order as the raw GFF file, with the
exception that the group column has been parsed into group class and
group name fields.

The feature ID, if provided, is a unique identifier of the feature
line.  The module does not depend on this ID in any way, but it is
available via Bio::DB::GFF-E<gt>id() if wanted.  In the dbi::mysql and
dbi::mysqlopt adaptor, the ID is a unique row ID.  In the acedb
adaptor it is not used.

=cut

# Given sequence name, range, and optional filter, retrieve list of
# all features.  Passes features through callback.
sub get_features {
  my $self = shift;
  my ($search,$options,$callback) = @_;
  $callback || $self->throw('must provide a callback argument');

  ## TODO: We don't yet use refclass at all!  Should we?
  my $rangetype = $search->{ 'rangetype' };
  my $chr       = $search->{ 'refseq' };
  my $start     = $search->{ 'start' };
  my $end       = $search->{ 'stop' };
  my $types     = $search->{ 'types' };

  ## TODO: REMOVE
  #warn "Hey.  got chr $chr.  types is [ ".$types->[ 0 ]." ]";

  if( $chr =~ /\:/ ) {
    my ( $ref_chr, $ref_start, $ref_end ) =
      ( $chr =~ /^[^\d]*(\d+)\:(\d+)\,(\d+)$/ );
    $chr = $ref_chr;
    ## TODO: Adjust the given start..
    if( $start ) {
      $self->throw( "Sorry, but we didn't yet implement the ability in jdrf_lod to query on subranges." );
    }
    $start = $ref_start;
    $end = $ref_end;
  } elsif( $chr =~ /[^\d]/ ) {
    my ( $just_the_digits ) = ( $chr =~ /^[^\d]+(\d+)$/ );
    unless( $just_the_digits ) {
      warn "Unable to parse out the digits from this chromosome id: '$chr'; expecting digits at the end of the string, preceded by non-digits, such as 'chr4'";
      return;
    } else {
      $chr = $just_the_digits;
    }
  }
  my $query_range =
    Bio::Range->new( '-seq_id' => $chr, '-start' => $start, '-end' => $end );

  ## TODO: Use all types.  For now, use just the first.
  if( !$types || !@$types || !@{ $types->[ 0 ] } ) {
    warn "No types!";
    warn $self->stack_trace_dump();
    return 0;
  }
  ## Note that for now we only really have one type, so we'll return
  ## after returning the lod features, even if multiple of the
  ## arguments have 'lod' in them.
  for( my $type_i; $type_i < @$types; $type_i++ ) {
    my ( $method, $source ) = @{ $types->[ $type_i ] };
    ## TODO: REMOVE
    #warn "Got type $method:$source";
    next unless ( $method eq 'lod' );

    ## TODO: REMOVE
    #warn "Hey.  Gonna do chr $chr, start $start, end $end.";

    my $dbh = $self->dbh();
    # First get the locus id.
    ## TODO: Also consider the rangetype...
    my ( $loc_id ) = $dbh->selectrow_array("
        SELECT	locus_id
        FROM		locus
        WHERE		chromosome = '$chr'
                AND (	(locus_start BETWEEN '$start' AND '$end')
                        OR ( locus_end BETWEEN '$start' AND '$end') 
                        OR ( ( locus_start <= '$start' ) AND ( locus_end >= '$end') ) )
                AND 	record_status_id <> 'D'
    ");
    unless( $loc_id ) {
      ## Ain't no lod data if we're not in a locus of some sort.
      return undef;
    }
    # Now get the strat ids
    my @str_ids = $self->_get_str_ids_for_segment( $loc_id );
    unless( @str_ids ) {
      ## Ain't no lod data if there ain't no str_id.
      return undef;
    }
    ## TODO: Deal with the fact that there will often be multiple str_ids.
    my $str_id = shift @str_ids; 
    ## TODO: REMOVE when we figure out what to do if there's multiple str_ids.
    if( @str_ids ) {
      warn "Oop.  There's multiple str_ids in this region: $chr:$start-$end.  Using only the first ($str_id), and ignoring these: ( ".join( ', ', @str_ids )." )";
    }

    ## TODO: REMOVE
    #warn "Hello again.  str_id is $str_id, loc_id is $loc_id.";
    
    my $sth =
      $dbh->prepare("
           SELECT marker_name,
               lod,
               chromosome_start,
               chromosome_end,
               fuzzy
           FROM marker_warehouse
           WHERE stratification_id = '$str_id'
      ");
    $sth->execute() ||
      warn "The statement failed!  errstr is ".$dbh->errstr;
    
    my %markers;
    my @last_fuzzy;
    my $last_end;
    my $marker_order = 1;
    while( my ( $mrk_name, $mrk_lod, $chrom_start, $chrom_end, $fuzzy ) =
           $sth->fetchrow_array() ) {

      ## TODO: Use $fuzzy.
      if( $chrom_start ) { # If there's a definite start position
        $markers{ $mrk_name }{ 'lod' }   = $mrk_lod;
        $markers{ $mrk_name }{ 'score' } = int( $mrk_lod * 100 );
        $markers{ $mrk_name }{ 'start' } = $chrom_start;
        $markers{ $mrk_name }{ 'end' }   = $chrom_end;
        $markers{ $mrk_name }{ 'fuzzy' } = 0;
        $markers{ $mrk_name }{ 'order' } = $marker_order;
        foreach my $fuzzy_marker ( @last_fuzzy ) {
          $markers{ $fuzzy_marker }{ 'end' } = $chrom_start;
        }
    
        undef @last_fuzzy;
        $last_end	= $chrom_end;
      } else {
        $markers{ $mrk_name }{ 'lod' }   = $mrk_lod;
        $markers{ $mrk_name }{ 'score' } = int( $mrk_lod * 100 );
        $markers{ $mrk_name }{ 'start' } = $last_end;
        $markers{ $mrk_name }{ 'fuzzy' } = 1;
        $markers{ $mrk_name }{ 'order' } = $marker_order;
        push( @last_fuzzy, $mrk_name );
      }
      ++$marker_order;
    }
    
    foreach my $fuzzy_marker ( @last_fuzzy ) {
      $markers{ $fuzzy_marker }{ 'end' } = $end;
    }

    ## We return one feature on either side of the range as well.
    my $handy_range =
      Bio::Range->new( '-seq_id' => $chr, '-start' => 1, '-end' => 1 );
    my $previous_mrk_name;
    my $overlaps;
    my $count;
    foreach my $mrk_name ( sort { $markers{ $a }{ 'start' } <=> $markers{ $b }{ 'start' } } keys %markers ) {
      $handy_range->start( $markers{ $mrk_name }{ 'start' } );
      $handy_range->end( $markers{ $mrk_name }{ 'end' } );
      ## TODO: REMOVE
      #warn "That marker is at ".$handy_range->start().", ".$handy_range->end().".";
      unless( $overlaps = $query_range->overlaps( $handy_range ) ) {
        ## We still need to return the one that's to the right of our range.
        if( $markers{ $mrk_name }{ 'start' } > $start ) {
          # If this is the first one we're returning..
          if( !$count && defined( $previous_mrk_name ) ) {
            ## Return the previous one too.
            $self->_do_callback( $callback, $previous_mrk_name, $chr, $str_id, $method, \%markers );
            $count++;
          }
          $self->_do_callback( $callback, $mrk_name, $chr, $str_id, $method, \%markers );
          $count++;
          # But only the first one to the right; now we're done.
          last;
        }
   
        $previous_mrk_name = $mrk_name;
        next;
      }
      # If this is the first one we're returning..
      if( !$count && defined( $previous_mrk_name ) ) {
        ## Return the previous one too.
        $self->_do_callback( $callback, $previous_mrk_name, $chr, $str_id, $method, \%markers );
        $count++;
      }
      $self->_do_callback( $callback, $mrk_name, $chr, $str_id, $method, \%markers );
      $previous_mrk_name = $mrk_name;
      $count++;
    }

    # We can return now, because we're only capable of returning
    # features of type 'lod', and we just did that...
    return $count;
  } # End for each of the given types..
  return 0;
} # get_features(..)

sub _do_callback {
  my $self = shift;
  my ( $callback, $mrk_name, $chr, $str_id, $method, $markers ) = @_;
  ## TODO: REMOVE
  #warn "This next one has a score of ".$markers->{ $mrk_name }{ 'score' }.".";
  $callback->( # 13 arguments, plus a hash of gsf tag values (ala Bio::SeqFeature::Generic)
    $chr,                             # $refseq      The reference sequence
    $markers->{ $mrk_name }{ 'start' }, # $start       feature start
    $markers->{ $mrk_name }{ 'end' },   # $stop        feature stop
    $str_id,                          # $source      feature source
    $method,                          # $method      feature method
    $markers->{ $mrk_name }{ 'score' }, # $score       feature score
    0,                                # $strand      feature strand
    1,                                # $phase       feature phase
    undef,                            # $groupclass  group class (may be undef)
    undef,                            # $groupname   group ID (may be undef)
    undef,                            # $tstart      target start for similarity hits (may be undef)
    undef,                            # $tstop       target stop for similarity hits (may be undef)
    $mrk_name,                        # $feature_id  A unique feature ID (may be undef)
    'lod:$str_id',                    # $group_id    The group id..
    ## Tags and their values
    'fuzzy' => $markers->{ $mrk_name }{ 'fuzzy' },
    'link'  => "http://jdrfdev.systemsbiology.net/cgi-bin/jdrf_publication.cgi?str_id=$str_id",
    'description'  => ( $mrk_name.': '.$markers->{ $mrk_name }{ 'score' } )
  );
} # _do_callback(..)

=head2 classes

 Title   : classes
 Usage   : $db->classes
 Function: return list of landmark classes in database
 Returns : a list of classes
 Args    : none
 Status  : public

Returns the empty list.

=cut

sub classes {
  return ();
}

=head2 _feature_by_name

 Title   : _feature_by_name
 Usage   : $db->get_features_by_name($name,$class,$callback)
 Function: get a list of features by name and class
 Returns : count of number of features retrieved
 Args    : name of feature, class of feature, and a callback
 Status  : protected

Not implemented.

=cut

sub _feature_by_name {
  shift->throw_not_implemented( @_ );
}

=head2 _feature_by_id

 Title   : _feature_by_id
 Usage   : $db->_feature_by_id($ids,$type,$callback)
 Function: get a list of features by ID
 Returns : count of number of features retrieved
 Args    : arrayref containing list of IDs to fetch and a callback
 Status  : protected

Not implemented

=cut

sub _feature_by_id {
  shift->throw_not_implemented( @_ );
}

sub _feature_by_attribute {
  shift->throw_not_implemented( @_ );
}

=head2 get_types

 Title   : get_types
 Usage   : $db->get_types($refseq,$refclass,$start,$stop,$count)
 Function: get list of types
 Returns : a list of Bio::DB::GFF::Typename objects
 Args    : see below
 Status  : Public

This method is responsible for fetching the list of feature type names
from the database.  The query may be limited to a particular range, in
which case the range is indicated by a landmark sequence name and
class and its subrange, if any.  These arguments may be undef if it is
desired to retrieve all feature types in the database (which may be a
slow operation in some implementations).

If the $count flag is false, the method returns a simple list of
vBio::DB::GFF::Typename objects.  If $count is true, the method returns
a list of $name=E<gt>$count pairs, where $count indicates the number of
times this feature occurs in the range.

Not Implemented (yet)
=cut

sub get_types {
  my $self = shift;
  my ( $chr, $class, $start, $stop, $count ) = @_;
  ## TODO: Dehackify
  if( $chr !~ /^\d/ ) {
    return;
  }
  unless( exists( $self->{ '__valid_ranges_hack' } ) ) {
    foreach my $valid_range qw(
           11:65695316,70981556
           14:83890000,91290000
           2:202980000,213140000
           10:109010000,116960000
           1:228150000,229170000
                               ) {
      my ( $vr_chr, $vr_start, $vr_end ) =
        ( $valid_range =~ /^(\d+):(\d+),(\d+)$/ );
      push( @{ $self->{ '__valid_ranges_hack' } },
            Bio::Range->new( '-seq_id' => $vr_chr, '-start' => $vr_start, '-end' => $vr_end ) );
    }
  }
  my $query_range =
    Bio::Range->new( '-seq_id' => $chr, '-start' => $start, '-end' => $stop );
  foreach my $valid_range ( @{ $self->{ '__valid_ranges_hack' } } ) {
    if( $valid_range->overlaps( $query_range ) ) {
      ## TODO: REMOVE
      #warn "Yes, a lod region";
      return ( 'lod' => 1 );
    }
  }
  ## TODO: REMOVE
  #warn "No, not a lod region";
  return undef;
} # get_types(..)

##### GET STRATIFICATION IDS FOR A SEGMENT
sub _get_str_ids_for_segment {
  my $self = shift;
  my ( $loc_id ) = @_;

  my $dbh = $self->dbh();
  my $sth = $dbh->prepare("
          SELECT	DISTINCT str.stratification_id
          FROM		stratification str,
                                  stratification_marker_linkage strmrk,
                                  marker mrk
          WHERE		strmrk.stratification_id = str.stratification_id
                  AND	strmrk.marker_id = mrk.marker_id
                  AND	mrk.locus_id = '$loc_id'
                  AND	mrk.record_status_id <> 'D'
                  AND	strmrk.record_status_id <> 'D'
                  AND	str.record_status_id <> 'D'
  ");
  $sth->execute() ||
      warn "The statement failed!  errstr is ".$dbh->errstr;
  
  my @str_ids;
  my $str_id;
  while( ( $str_id ) = $sth->fetchrow_array() ) {
    push( @str_ids, $str_id );
  }
  return @str_ids;
} # _get_str_ids_for_segment(..)

=head2 meta

 Title   : meta
 Usage   : $value = $db->meta($name [,$newval])
 Function: get or set a meta variable
 Returns : a string
 Args    : meta variable name and optionally value
 Status  : public

Get or set a named metavariable for the database.  Metavariables can
be used for database-specific settings.  This method calls two
class-specific methods which must be implemented:

  make_meta_get_query()   Returns a sql fragment which given a meta
                          parameter name, returns its value.  One bind
                          variable.
  make_meta_set_query()   Returns a sql fragment which takes two bind
                          arguments, the parameter name and its value


Don't make changes unless you know what you're doing!  It will affect the
persistent database.

Not Implemented.

=cut

sub meta {
  shift->throw_not_implemented( @_ );
}

sub dna_chunk_size {
  shift->throw_not_implemented( @_ );
}

=head2 get_features_iterator

 Title   : get_features_iterator
 Usage   : $iterator = $db->get_features_iterator($search,$options,$callback)
 Function: create an iterator on a features() query
 Returns : A Bio::DB::GFF::Adaptor::dbi::iterator object
 Args    : see get_features()
 Status  : public

This method is similar to get_features(), except that it returns an
iterator across the query.  See
L<Bio::DB::GFF::Adaptor::dbi::iterator>.

Not Implemented.

=cut

sub get_features_iterator {
  shift->throw_not_implemented( @_ );
}

########################## loading and initialization  #####################

=head2 do_initialize

 Title   : do_initialize
 Usage   : $success = $db->do_initialize($drop_all)
 Function: initialize the database
 Returns : a boolean indicating the success of the operation
 Args    : a boolean indicating whether to delete existing data
 Status  : protected

This method will load the schema into the database.  If $drop_all is
true, then any existing data in the tables known to the schema will be
deleted.

Implemented to just silently return 1.

=cut

# Create the schema from scratch.
# You will need create privileges for this.
sub do_initialize {
  1;
}

=head2 DESTROY

 Title   : DESTROY
 Usage   : $db->DESTROY
 Function: disconnect database at destruct time
 Returns : void
 Args    : none
 Status  : protected

This is the destructor for the class.

=cut

sub DESTROY {
  my $self = shift;
  $self->features_db->disconnect if defined $self->features_db;
}

1;

__END__

=head1 BUGS

=head1 SEE ALSO

L<Bio::DB::GFF>, L<bioperl>

=head1 AUTHOR

Michelle Whiting E<lt>mwhiting@systemsbiology.orgE<gt>.

Copyright (c) 2003 Institute for Systems Biology.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

