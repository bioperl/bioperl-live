# $Id$
#
# BioPerl module for Bio::SeqFeature::Collection
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
  # let's first input some features
  my $gffio = Bio::Tools::GFF->new(-file => Bio::Root::IO->catfile
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
  my $col = new Bio::SeqFeature::Collection();
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

=head1 DESCRIPTION

This object will efficiently allow one for query subsets of ranges
within a large collection of sequence features (in fact the objects
just have to be Bio::RangeI compliant).  This is done by the creation
of bins which are stored in order in a B-Tree data structure as
provided by the DB_File interface to the Berkeley DB.

This is based on work done by Lincoln for storage in a mysql instance
- this is intended to be an embedded in-memory implementation for
easily quering for subsets of a large range set.  All features are
held in memory even if the -usefile flag is provided.

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
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 CONTRIBUTORS

Using code and strategy developed by Lincoln Stein (lstein@cshl.org)
in Bio::DB::GFF implementation.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Collection;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Root::IO;
use Bio::DB::GFF::Util::Binning;
use DB_File;
use Bio::Location::Simple;

@ISA = qw(Bio::Root::Root );


# This may need to get re-optimized for BDB usage as these
# numbers were derived empirically by Lincoln on a mysql srv
# running on his laptop

# this is the largest that any reference sequence can be (100 megabases)
use constant MAX_BIN    => 100_000_000;

# this is the smallest bin (1 K)
use constant MIN_BIN    => 1_000;

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::SeqFeature::Collection();
 Function: Builds a new Bio::SeqFeature::Collection object 
 Returns : Bio::SeqFeature::Collection
 Args    :

           -minbin        minimum value to use for binning
                          (default is 100,000,000)
           -maxbin        maximum value to use for binning
                          (default is 1,000)
           -usefile       boolean to use a file to store
                          BTREE rather than an in-memory structure 
                          (default is false or in-memory).

           -features      Array ref of features to add initially

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($maxbin,$minbin,$usefile,$features) = $self->_rearrange([qw(MAXBIN MINBIN
							 USEFILE
							 FEATURES)],@args);

  defined $maxbin && $self->max_bin($maxbin);
  defined $minbin && $self->min_bin($minbin);

  defined $features &&  $self->add_features($features);

  my $tmpname = undef;
  if( $usefile ) { 
      $self->{'_io'} = new Bio::Root::IO;
      $tmpname = $self->{'_io'}->tempfile();
      unlink($tmpname);
      $self->debug("tmpfile is $tmpname");
  }
  $DB_BTREE->{'flags'} = R_DUP ;
  $DB_BTREE->{'compare'} = \&_compare;
  $self->{'_btreehash'} = {};
  $self->{'_btree'} = tie %{$self->{'_btreehash'}}, 
      'DB_File', $tmpname, O_RDWR|O_CREAT, 0640, $DB_BTREE;

#  possibly storing/retrieving as floats for speed improvement?
#  $self->{'_btree'}->filter_store_key  ( sub { $_ = pack ("f", $_) } );
#  $self->{'_btree'}->filter_fetch_key  ( sub { $_ = unpack("f", $_) } );

  $self->{'_features'} = [];
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
       if( ! $f->isa('Bio::RangeI') ) {
	   $self->warn("Must provide valid Bio::RangeI objects to add_features, skipping object '$f'\n");
	   next;
       }
       my $bin = bin($f->start,$f->end,$self->min_bin);       
       $self->debug( "$bin for ". $f->location->to_FTstring(). "\n");
       push @{$self->{'_features'}}, $f;
       $self->{'_btreehash'}->{$bin} = $#{$self->{'_features'}};
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
   $strand |= 1;
   if( !defined $start && !defined $end ) {
       if( ! defined $range || !ref($range) || ! $range->isa("Bio::RangeI") ) 
       { 
	   $self->warn("Must defined a valid Range for the method feature_in_range");
	   return ();
       }
       ($start,$end,$strand) = ($range->start,$range->end,$range->strand);
   }
   my $r = new Bio::Location::Simple(-start => $start,
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
	   if( @vals ) {
	       push @bins, map { $self->{'_features'}->[$_] } @vals;
	   } 
       } else {	   
	   $k = $tier_start;
	   $self->{'_btree'}->seq($k,$v, R_CURSOR);	   
	   my @vals;
	   while( $self->{'_btree'}->seq($k,$v, R_NEXT) == 0 ) {
	       last if( $k > $tier_stop || $k < $tier_start);
	       push @vals, $v;
	   }
	   push @bins, map { $self->{'_features'}->[$_] } @vals;
       }
       $tier /= 10;
   }   
   
   $strandmatch = 'ignore' unless defined $strandmatch;
   return ( $contain ) ? grep { $r->contains($_,$strandmatch) } @bins : 
       grep { $r->overlaps($_,$strandmatch)} @bins;
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
sub _compare{ $_[0] <=> $_[1]}

sub _comparepack { unpack("f", $_[0]) <=> unpack("f", $_[1]) ;}

sub DESTROY { 
    my $self = shift;
    $self->SUPER::DESTROY();
    if( defined $self->{'_io'} )  {
	$self->{'_io'}->_io_cleanup();    
	$self->{'_io'} = undef;
    }
}

1;
