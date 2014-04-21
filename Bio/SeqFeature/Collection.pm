#
# BioPerl module for Bio::SeqFeature::Collection
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Collection - A container class for SeqFeatures
suitable for performing operations such as finding features within a
range, that match a certain feature type, etc.

=head1 SYNOPSIS

  use Bio::SeqFeature::Collection;
  use Bio::Location::Simple;
  use Bio::Tools::GFF;
  use Bio::Root::IO;
  use File::Spec;
  # let's first input some features
  my $gffio = Bio::Tools::GFF->new(-file => File::Spec->catfile
  				 ("t","data","myco_sites.gff"),
  				 -gff_version => 2);
  my @features = ();
  # loop over the input stream
  while(my $feature = $gffio->next_feature()) {
      # do something with feature
      push @features, $feature;
  }
  $gffio->close();
  # build the Collection object
  my $col = Bio::SeqFeature::Collection->new();
  # add these features to the object
  my $totaladded = $col->add_features(\@features);

  my @subset = $col->features_in_range(-start => 1,
  				     -end => 25000,
  				     -strand => 1,
  				     -contain => 0);
  # subset should have 18 entries for this dataset
  print "size is ", scalar @subset, "\n";
  @subset = $col->features_in_range(-range => Bio::Location::Simple->new
  				  (-start => 70000,
  				   -end => 150000,
  				   -strand => -1),
  				  -contain => 1,
  				  -strandmatch => 'strong');

  # subset should have 22 entries for this dataset
  print "size is ", scalar @subset, "\n";
  print "total number of features in collection is ",
         $col->feature_count(),"\n";

=head1 DESCRIPTION

This object will efficiently allow one for query subsets of ranges
within a large collection of sequence features (in fact the objects
just have to be Bio::RangeI compliant).  This is done by the creation
of bins which are stored in order in a B-Tree data structure as
provided by the DB_File interface to the Berkeley DB.

This is based on work done by Lincoln for storage in a mysql instance
- this is intended to be an embedded in-memory implementation for
easily querying for subsets of a large range set.

Collections can be made persistent by keeping the indexfile and
passing in the -keep flag like this:

  my $collection = Bio::SeqFeature::Collection->new(-keep => 1,
                                                   -file => 'col.idx');
  $collaction->add_features(\@features);
  undef $collection;

  # To reuse this collection, next time you initialize a Collection object
  # specify the filename and the index will be reused.
  $collection = Bio::SeqFeature::Collection->new(-keep => 1,
                                                -file => 'col.idx');



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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 CONTRIBUTORS

Using code and strategy developed by Lincoln Stein (lstein@cshl.org)
in Bio::DB::GFF implementation.  Credit also to Lincoln for suggesting
using Storable to serialize features rather than my previous implementation
which kept the features in memory.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Collection;
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::DB::GFF::Util::Binning;
use DB_File;
use Bio::Location::Simple;
use Bio::SeqFeature::Generic;
use Storable qw(freeze thaw);

use base qw(Bio::Root::Root Bio::SeqFeature::CollectionI);


# This may need to get re-optimized for BDB usage as these
# numbers were derived empirically by Lincoln on a mysql srv
# running on his laptop

# this is the largest that any reference sequence can be (100 megabases)
use constant MAX_BIN    => 100_000_000;

# this is the smallest bin (1 K)
use constant MIN_BIN    => 1_000;

=head2 new

 Title   : new
 Usage   : my $obj = Bio::SeqFeature::Collection->new();
 Function: Builds a new Bio::SeqFeature::Collection object
 Returns : Bio::SeqFeature::Collection
 Args    :

           -minbin        minimum value to use for binning
                          (default is 100,000,000)
           -maxbin        maximum value to use for binning
                          (default is 1,000)
           -file          filename to store/read the
                          BTREE from rather than an in-memory structure
                          (default is false and in-memory).
           -keep          boolean, will not remove index file on
                          object destruction.
           -features      Array ref of features to add initially

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($maxbin,$minbin, $file, $keep,
      $features) = $self->_rearrange([qw(MAXBIN MINBIN FILE KEEP
					 FEATURES)],@args);

  defined $maxbin && $self->max_bin($maxbin);
  defined $minbin && $self->min_bin($minbin);

  defined $features &&  $self->add_features($features);
  $DB_BTREE->{'flags'} = R_DUP ;
  $DB_BTREE->{'compare'} = \&_compare;
  $self->{'_btreehash'} = {};
  if( $file ) {
      $self->debug("using file $file");
      $self->indexfile($file);
  }
  $self->keep($keep);
  $self->{'_btree'} = tie %{$self->{'_btreehash'}},
  'DB_File', $self->indexfile, O_RDWR|O_CREAT, 0640, $DB_BTREE;
  $self->{'_btree'} || $self->throw("Unable to tie DB_File handle");
  return $self;
}


=head2 add_features

 Title   : add_features
 Usage   : $collection->add_features(\@features);
 Function:
 Returns : number of features added
 Args    : arrayref of Bio::SeqFeatureI objects to index


=cut

sub add_features{
   my ($self,$feats) = @_;
   if( ref($feats) !~ /ARRAY/i ) {
       $self->warn("Must provide a valid Array reference to add_features");
       return 0;
   }
   my $count = 0;
   foreach my $f ( @$feats ) {
       if( ! $f || ! ref($f) || ! $f->isa('Bio::RangeI') ) {
	   $self->warn("Must provide valid Bio::RangeI objects to add_features, skipping object '$f'\n");
	   next;
       }
       my $bin = bin($f->start,$f->end,$self->min_bin);
       my $serialized = &feature_freeze($f);
       $self->{'_btree'}->put($bin,$serialized);
       if( $f->isa('Bio::SeqFeature::Generic') ) {
	   $self->debug( "$bin for ". $f->location->to_FTstring(). " matches ".$#{$self->{'_features'}}. "\n");
       }
       $count++;
   }
   return $count;
}


=head2 features_in_range

 Title   : features_in_range
 Usage   : my @features = $collection->features_in_range($range)
 Function: Retrieves a list of features which were contained or overlap the
           the requested range (see Args for way to specify overlap or
				only those containe)d
 Returns : List of Bio::SeqFeatureI objects
 Args    : -range => Bio::RangeI object defining range to search,
           OR
           -start  => start,
           -end    => end,
           -strand  => strand

           -contain => boolean - true if feature must be completely
                       contained with range
                       OR false if should include features that simply overlap
                       the range. Default: true.
           -strandmatch =>  'strong',  ranges must have the same strand
                            'weak',    ranges must have the same
                                           strand or no strand
                            'ignore', ignore strand information
                           Default. 'ignore'.

=cut

sub features_in_range{
   my $self = shift;
   my (@args) = @_;
   my ($range, $contain, $strandmatch,$start,$end,$strand);
   if( @args == 1 ) {
       $range = shift @args;
   } else {
       ($start,$end,$strand,$range,
	$contain,$strandmatch) = $self->_rearrange([qw(START END
						       STRAND
						       RANGE CONTAIN
						       STRANDMATCH)],
						   @args);
       $contain = 1 unless defined $contain;
   }
   $strand = 1 unless defined $strand;
   if( $strand !~ /^([\-\+])$/ &&
       $strand !~ /^[\-\+]?1$/ ) {
       $self->warn("must provide a valid numeric or +/- for strand");
       return ();
   }
   if( defined $1 ) { $strand .= 1; }

   if( !defined $start && !defined $end ) {
       if( ! defined $range || !ref($range) || ! $range->isa("Bio::RangeI") )
       {
	   $self->warn("Must defined a valid Range for the method feature_in_range");
	   return ();
       }
       ($start,$end,$strand) = ($range->start,$range->end,$range->strand);
   }
   my $r = Bio::Location::Simple->new(-start => $start,
				     -end   => $end,
				     -strand => $strand);

   my @features;
   my $maxbin = $self->max_bin;
   my $minbin = $self->min_bin;
   my $tier = $maxbin;
   my ($k,$v,@bins) = ("",undef);
   while ($tier >= $minbin) {
	my ($tier_start,$tier_stop) = (bin_bot($tier,$start),
				       bin_top($tier,$end));
       if( $tier_start == $tier_stop ) {
	   my @vals = $self->{'_btree'}->get_dup($tier_start);
	   if( scalar @vals > 0 ) {
	       push @bins, map { thaw($_) } @vals;
	   }
       } else {	
	   $k = $tier_start;
	   my @vals;
	   for( my $rc = $self->{'_btree'}->seq($k,$v,R_CURSOR);
	        $rc == 0;
	        $rc = $self->{'_btree'}->seq($k,$v, R_NEXT) ) {
	       last if( $k > $tier_stop || $k < $tier_start);
	       push @bins, thaw($v);
	   }
       }
       $tier /= 10;
   }
   my %seen = ();
   foreach my $t ( map { ref($_) } @bins) {
       next if $seen{$t}++;
       eval "require $t";

       if( $@ ) {
	   $self->warn("Trying to thaw a stored feature $t which does not appear in your Perl library. $@");
	   next;
       }
   }
   $strandmatch = 'ignore' unless defined $strandmatch;
   return ( $contain ) ? grep { $r->contains($_,$strandmatch) } @bins :
       grep { $r->overlaps($_,$strandmatch)} @bins;
}

=head2 remove_features

 Title   : remove_features
 Usage   : $collection->remove_features(\@array)
 Function: Removes the requested sequence features (based on features
	   which have the same location)
 Returns : Number of features removed
 Args    : Arrayref of Bio::RangeI objects


=cut

sub remove_features{
   my ($self,$feats) = @_;
   if( ref($feats) !~ /ARRAY/i ) {
       $self->warn("Must provide a valid Array reference to remove_features");
       return 0;
   }
   my $countprocessed = 0;

   foreach my $f ( @$feats ) {
       next if ! ref($f) || ! $f->isa('Bio::RangeI');
       my $bin = bin($f->start,$f->end,$self->min_bin);
       my @vals = $self->{'_btree'}->get_dup($bin);
       my $vcount = scalar @vals;

       foreach my $v ( @vals )  {
	   # Once we have uniquely identifiable field
	   # I think it will work better.
	   if( $v eq &feature_freeze($f) ) {
	       $self->{'_btree'}->del_dup($bin,$v);
	       $vcount--;
	       $countprocessed++;
	   }
       }
       if( $vcount == 0 ) {
	   $self->{'_btree'}->del($bin);
       }
   }
   $countprocessed;

}

=head2 get_all_features

 Title   : get_all_features
 Usage   : my @f = $col->get_all_features()
 Function: Return all the features stored in this collection (Could be large)
 Returns : Array of Bio::RangeI objects
 Args    : None


=cut

sub get_all_features{
   my ($self) = @_;
   my @features;
   my ($key,$value);
   for (my $status = $self->{'_btree'}->seq($key, $value, R_FIRST) ;
	$status == 0 ;
	$status = $self->{'_btree'}->seq($key, $value, R_NEXT) )
   {   next unless defined $value;
       push @features, &thaw($value);
   }
   if( scalar @features !=  $self->feature_count() ) {
       $self->warn("feature count does not match actual count\n");
   }
   return @features;
}


=head2 min_bin

 Title   : min_bin
 Usage   : my $minbin= $self->min_bin;
 Function: Get/Set the minimum value to use for binning
 Returns : integer
 Args    : [optional] minimum bin value


=cut

sub min_bin {
  my ($self,$min) = @_;
  if( defined $min ) {
      $self->{'_min_bin'} = $min;
  }
  return $self->{'_min_bin'}  || MIN_BIN;
}

=head2 max_bin

 Title   : max_bin
 Usage   : my $maxbin= $self->max_bin;
 Function: Get/Set the maximum value to use for binning
 Returns : integer
 Args    : [optional] maximum bin value


=cut

sub max_bin {
  my ($self,$max) = @_;
  if( defined $max ) {
      $self->{'_max_bin'} = $max;
  }
  return $self->{'max_bin'} || MAX_BIN;
}

=head2 feature_count

 Title   : feature_count
 Usage   : my $c = $col->feature_count()
 Function: Retrieve the total number of features in the collection
 Returns : integer
 Args    : none


=cut

sub feature_count {
    my $self = shift;
    my $count = 0;
    for ( keys %{$self->{'_btreehash'}} ) {
	my $v = $self->{'_btreehash'}->{$_};
	next unless defined  $v;
	$count++;
    }
    $count;
}

=head2 indexfile

 Title   : indexfile
 Usage   : $obj->indexfile($newval)
 Function: Get/set the filename where index is kept
 Returns : value of indexfile (a filename string)
 Args    : on set, new value (a filename string )


=cut

sub indexfile{
    my $self = shift;

    return $self->{'indexfile'} = shift if @_;
    return $self->{'indexfile'};
}

=head2 keep

 Title   : keep
 Usage   : $obj->keep($newval)
 Function: Get/set boolean flag to keep the indexfile after
           exiting program
 Example :
 Returns : value of keep (boolean)
 Args    : on set, new value (boolean)


=cut

sub keep{
    my $self = shift;

    return $self->{'keep'} = shift if @_;
    return $self->{'keep'};
}

sub _compare{
    if( defined $_[0] && ! defined $_[1]) {
	return -1;
    } elsif ( defined $_[1] && ! defined $_[0]) {
	return 1;
    }
    $_[0] <=> $_[1];
}

sub feature_freeze {
    my $obj = shift;
    _remove_cleanup_methods($obj);
    return freeze($obj);
}

sub _remove_cleanup_methods {
    my $obj = shift;
    
    # we have to remove any cleanup methods here for Storable
    for my $funcref ( $obj->_cleanup_methods ) {
        $obj->_unregister_for_cleanup($funcref);
    }
    
    # ... and the same for any contained features; hopefully any implementations
    # adhere to implementing Bio::SeqFeatureI::sub_SeqFeature
    
    for my $contained ($obj->sub_SeqFeature) {
        _remove_cleanup_methods($contained);
    }
    
    1;
}

sub feature_thaw {
    return thaw(shift);
}

sub DESTROY {
    my $self = shift;
    $self->{'_btree'} = undef;
    untie(%{$self->{'_btreehash'}});
    if( ! $self->keep && $self->indexfile ) {
        my $f = $self->indexfile;
        $self->debug( "unlinking ".$f. "\n");
        close($f);
        unlink($f);
    }
    $self->SUPER::DESTROY();    
}

1;
