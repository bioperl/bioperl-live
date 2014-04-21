#-------------------------------------------------------------------------------
#
# BioPerl module Bio::Restriction::EnzymeCollection
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Rob Edwards <redwards@utmem.edu>
#
# You may distribute this module under the same terms as perl itself
#-------------------------------------------------------------------------------

## POD Documentation:

=head1 NAME

Bio::Restriction::EnzymeCollection - Set of restriction endonucleases

=head1 SYNOPSIS

  use Bio::Restriction::EnzymeCollection;

  # Create a collection with the default enzymes.
  my $default_collection = Bio::Restriction::EnzymeCollection->new();

  # Or create a collection from a REBASE 'withrefm' file obtained from
  # ftp://ftp.neb.com/pub/rebase/. (See Bio::Restriction::IO for more
  # information.)
  my $rebase = Bio::Restriction::IO->new(
      -file   => 'withrefm.610',
      -format => 'withrefm' );
  my $rebase_collection = $rebase->read();

  # Or create an empty collection and set the enzymes later. See
  # 'CUSTOM COLLECTIONS' below for more information.
  my $empty_collection =
    Bio::Restriction::EnzymeCollection->new( -empty => 1 );

  # Get an array of Bio::Restriction::Enzyme objects from the collection.
  my @enzymes = $default_collection->each_enzyme();

  # Get a Bio::Restriction::Enzyme object for a particular enzyme by name.
  my $enz = $default_collection->get_enzyme( 'EcoRI' );

  # Get a Bio::Restriction::EnzymeCollection object containing the enzymes
  # that have the equivalent of 6-bp recognition sequences.
  my $six_cutters = $default_collection->cutters( 6 );

  # Get a Bio::Restriction::EnzymeCollection object containing the enzymes
  # that are rare cutters.
  my $rare_cutters = $default_collection->cutters( -start => 6, -end => 8 );

  # Get a Bio::Restriction::EnzymeCollection object that contains enzymes
  # that generate blunt ends:
  my $blunt_cutters = $default_collection->blunt_enzymes();

  # See 'CUSTOM COLLECTIONS' below for an example of creating a
  # Bio::Restriction::EnzymeCollection object with a specified subset of
  # enzymes using methods provided by the Bio::RestrictionEnzyme class.

=head1 DESCRIPTION

Bio::Restriction::EnzymeCollection represents a collection of
restriction enzymes.

If you create a new collection directly rather than from a REBASE
file using L<Bio::Restriction::IO>, it will be populated by a
default set of enzymes with site and cut information
only.

Use L<Bio::Restriction::Analysis> to figure out which enzymes are
available and where they cut your sequence.

=head1 CUSTOM COLLECTIONS

Note that the underlying L<Bio::Restriction::Enzyme> objects have a rich
variety of methods that allow more complicated selections than the methods
that are defined by Bio::Restriction::EnzymeCollection.

For example, the way to create a custom collection of Type II enzymes
is as follows:

  my $complete_collection =
      Bio::Restriction::EnzymeCollection->new();
  my $type_ii_collection  =
      Bio::Restriction::EnzymeCollection->new( -empty => 1 );
  $type_ii_collection->enzymes(
      grep { $_->type() eq 'II' } $complete_collection->each_enzyme() );

=head1 SEE ALSO

L<Bio::Restriction::IO> - read in enzymes from REBASE files

L<Bio::Restriction::Analysis> - figure out what enzymes cut a sequence

L<Bio::Restriction::Enzyme> - define a single restriction enzyme

=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

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
the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Rob Edwards, redwards@utmem.edu

=head1 CONTRIBUTORS

Heikki Lehvaslaiho, heikki-at-bioperl-dot-org

=head1 COPYRIGHT

Copyright (c) 2003 Rob Edwards.

Some of this work is Copyright (c) 1997-2002 Steve A. Chervitz. All
Rights Reserved.

This module is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=head1 APPENDIX

Methods beginning with a leading underscore are considered private and
are intended for internal use by this module. They are not considered
part of the public interface and are described here for documentation
purposes only.

=cut


package Bio::Restriction::EnzymeCollection;
use strict;

use Bio::Restriction::Enzyme;
use Bio::Restriction::IO;

use Data::Dumper;

use base qw(Bio::Root::Root);

=head2 new

 Title     : new
 Function  : Initializes the Restriction::EnzymeCollection object
 Returns   : The Restriction::EnzymeCollection object
 Arguments : optional named parameter -empty

Set parameter -empty to true if you do NOT want the collection be
populated by the default set of prototype type II enzymes.

Alternatively, pass an array of enzymes to -enzymes parameter.

=cut

sub new {
    my($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($empty, $enzymes) =
            $self->_rearrange([qw(
                                  EMPTY
                                  ENZYMES
                                 )], @args);

    $self->{'_all_enzymes'} = [];
    $self->{'_enzymes'} = {};

    return $self if $empty;


    if ($enzymes) {
	# as advertised in pod/maj
	$self->throw( "Arg to -enzymes must be an arrayref to Bio::Restriction::Enzyme objects") unless ref($enzymes) eq 'ARRAY';
	$self->enzymes(@$enzymes);
	return $self;
    }
    else {
	# the default set of enzymes
	my $in  = Bio::Restriction::IO->new(-verbose => $self->verbose);
	return $in->read;
    }
}

=head2 Manipulate the enzymes within the collection

=cut

=head2 enzymes

 Title     : enzyme
 Function  : add/get method for enzymes and enzyme collections
 Returns   : object itself
 Arguments : array of Bio::Restriction::Enzyme and
             Bio::Restriction::EnzymeCollection objects

=cut

sub enzymes {
    my ($self, @enzs)=@_;
    foreach my $e (@enzs) {
        if ( ref $e eq '') {
            print "|$e|\n";
        }
        elsif ($e->isa('Bio::Restriction::EnzymeI')) {
            push(@{$self->{'_all_enzymes'}},$e);
            $self->{'_enzymes'}->{$e->name} = $e;
        }
        elsif ($e->isa('Bio::Restriction::EnzymeCollection')) {
           $self->enzymes($e->each_enzyme);
        } else {
            my $r = 1;
            $self->warn("EnzymeCollection can not deal with ".
                        ref($e)." objects");
        }
    }
    return $self;
}

#
# method to remove duplicates?
#

=head2 each_enzyme

 Title     : each_enzyme
 Function  : get an array of enzymes
 Returns   : array of Bio::Restriction::Enzyme objects
 Arguments : -

=cut

sub each_enzyme {
    my $self = shift;
    return @{$self->{'_all_enzymes'}};
}

=head2 get_enzyme

 Title     : get_enzyme
 Function  : Gets a Bio::Restriction::Enzyme object for the enzyme name
 Returns   : A Bio::Restriction::Enzyme object or undef
 Arguments : An enzyme name that is in the collection

=cut

sub get_enzyme {
    my ($self, $name)=@_;
    return $self->{'_enzymes'}->{$name};
}


=head2 available_list

 Title     : available_list
 Function  : Gets a list of all the enzymes that we know about
 Returns   : A reference to an array with all the enzyme names
             that we have defined or 0 if none are defined
 Arguments : Nothing
 Comments  : Note, I maintain this for backwards compatibility,
             but I don't like the name as it is very ambiguous

=cut

sub available_list {
    my ($self, $size)=@_;
    my @keys = sort keys %{$self->{'_enzymes'}};
    return @keys;
}

=head2 longest_cutter

 Title     : longest_cutter
 Function  : Gets the enzyme with the longest recognition site
 Returns   : A Bio::Restriction::Enzyme object
 Arguments : Nothing
 Comments  : Note, this is used by Bio::Restriction::Analysis
             to figure out what to do with circular sequences

=cut

sub longest_cutter {
    my ($self)=@_;
    my $longest=0; my $longest_enz='.';
    foreach my $enz ($self->each_enzyme) {
     my $len=$enz->recognition_length;
     if ($len > $longest) {$longest=$len; $longest_enz=$enz}
    }
    return $longest_enz;
}

=head2 Filter enzymes

=cut

=head2 blunt_enzymes

  Title     : blunt_enzymes
  Function  : Gets a list of all the enzymes that are blunt cutters
  Returns   : A reference to an array with all the enzyme names that
              are blunt cutters or 0 if none are defined
  Arguments : Nothing
  Comments  : 

This is an example of the kind of filtering better done by the scripts
using the rich collection of methods in Bio::Restriction::Enzyme.

=cut

sub blunt_enzymes {
    my $self=shift;
    my $bs = Bio::Restriction::EnzymeCollection->new(-empty => 1);
    return $bs->enzymes(  grep { $_->overhang eq 'blunt' }  $self->each_enzyme );
}


=head2 cutters

  Title     : cutters
  Function  : Gets a list of all the enzymes that recognize a
              certain size, e.g. 6-cutters
  Usage     : $cutters = $collection->cutters(6);
  Returns   : A reference to an array with all the enzyme names
              that are x cutters or 0 if none are defined
  Arguments : A positive number for the size of cutters to return
              OR
              A range: (-start => 6, -end => 8,
                        -inclusive => 1, -exclusive = 0 )

The default for a range is 'inclusive'


=cut

sub cutters {
    my ($self) = shift;

    return unless @_; # no argument

    if (scalar @_ == 1 ) {
        my $size = shift;
        my @sizes;
        (ref $size eq 'ARRAY') ? push @sizes, @{$size} : push @sizes, $size;
        my $bs = Bio::Restriction::EnzymeCollection->new(-empty => 1);
        for my $size (@sizes) {
            $self->throw("Need a positive number [$size]")
                unless $size =~ /[+]?[\d\.]+/;
            foreach my $e ($self->each_enzyme) {
                ##print $e->name, ": ", $e->cutter, "\n"  if $e->cutter == $size;
                $bs->enzymes($e) if $e->cutter == $size;
            }
        }
        return $bs;

    } else { # named arguments

        my ($start, $end, $inclusive, $exclusive ) =
            $self->_rearrange([qw(
                                  START
                                  END
                                  INCLUSIVE
                                  EXCLUSIVE
                                 )], @_);

        $self->throw("Start needs a positive number [$start]")
            unless $start =~ /[+]?[\d\.]+/;
        $self->throw("End needs a positive number [$end]")
            unless $end =~ /[+]?[\d\.]+/;

        my $limits;
        $inclusive = 1 if $inclusive or not $exclusive;
        $inclusive = 0 if $exclusive;

        my $bs = Bio::Restriction::EnzymeCollection->new(-empty => 1);
        if ($inclusive) {
            foreach my $e ($self->each_enzyme) {
                $bs->enzymes($e) if $e->cutter >= $start and $e->cutter <= $end;
            }
        } else {
            foreach my $e ($self->each_enzyme) {
                $bs->enzymes($e) if $e->cutter > $start and $e->cutter < $end;
            }
        }
        return $bs;
    }
}


1;
