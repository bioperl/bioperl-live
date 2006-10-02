# $Id$
#-------------------------------------------------------------------------------
#
# BioPerl module Bio::Restriction::EnzymeCollection
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

  # create a set with a default enzymes
  my $collection = Bio::Restriction::EnzymeCollection;


  print "No of enzymes: ", scalar $collection->each_enzyme, "\n";

  # find something about a particular enzyme
  my $enz=$collection->get_enzyme('EcoRI');

  # and we know about some special types of enzymes.
  my $blunt_collection = $collection->blunt_enzymes; # enzymes without an overhang

  # see 'CUSTOM COLLECTIONS' below  for more options

  # the most common selection criteria is how many times the enzymes
  # cuts. This can be estimated using the length and specificity of
  # the recognition site

  # enzymes that have euivalent of  6bp recognition sequence
  my $six_cutters=$collection->cutters(6);

  # all rare cutters
  my $rare_cutters=$collection->cutters(-start=>6, -end=>8);


=head1 DESCRIPTION

Bio::Restriction::EnzymeCollection represents a collection of
restriction enzymes.

If you create a new collection directly rather than from a REBASE
format file using Bio::RestrictionIO, it will be populated by a
default set of protype, typeII enzymes with site and cut information
only.

Use Bio::Restriction::Analysis to figure out which enzymes are
available and where they cut your sequence.


=head1 CUSTOM COLLECTIONS

Note, that the underlying Enzyme objects are much more rich and allow
more complicated selections than the predefinend methods. The way to
create a custom subset is as follows:

  my $initial_collection;
  my $new_collection = Bio::Restriction::EnzymeCollection(-empty => 1);
  foreach $enzyme ($initial_collection) {
      # this selects only type II enzymes
      $new_collection($enzyme) if $enzyme->type eq 'II';
  }

=head1 COMMENTS

I am trying to make this backwards compatible with
Bio::Tools::RestrictionEnzyme. Undoubtedly some things will break, but
we can fix things as we progress.....!

I have added another comments section at the end of this POD that
discusses a couple of areas I know are broken (at the moment)

=head1 SEE ALSO

L<Bio::Restriction::IO> - read in enzymes from REBASE files
(s

L<Bio::Restriction::Analysis> - figure out what enzymes cut a sequence
(start here)

L<Bio::Restriction::Enzyme> - defining a single restriction enzyme


=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

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

use vars qw ();
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

    my ($empty) =
            $self->_rearrange([qw(
                                  EMPTY
                                 )], @args);

    $self->{'_all_enzymes'} = [];
    $self->{'_enzymes'} = {};

    return $self if $empty;

    # the default set of enzymes
    my $in  = Bio::Restriction::IO->new(-verbose => $self->verbose);
    return $in->read;

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
    my $bs = new Bio::Restriction::EnzymeCollection(-empty => 1);
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
        $self->throw("Need a positive number [$size]")
            unless $size =~ /[+]?[\d\.]+/;

        my $bs = new Bio::Restriction::EnzymeCollection(-empty => 1);

        foreach my $e ($self->each_enzyme) {
            ##print $e->name, ": ", $e->cutter, "\n"  if $e->cutter == $size;
            $bs->enzymes($e) if $e->cutter == $size;
        }
        return $bs;
        #return $bs->enzymes(  grep { ($_->cutter == $size) }  $self->each_enzyme );

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

        my $bs = new Bio::Restriction::EnzymeCollection(-empty => 1);
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
