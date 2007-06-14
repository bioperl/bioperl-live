# $Id$
# BioPerl module for Bio::Restriction::IO::itype2
#
# Cared for by Rob Edwards <redwards@utmem.edu>
#
# Copyright Rob Edwards
#
# You may distribute this module under the same terms as perl itself
#

# POD documentation - main docs before the code

=head1 NAME

Bio::Restriction::IO::itype2 - itype2 enzyme set

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::Restriction::IO class.

=head1 DESCRIPTION

This is tab delimited, entry per line format which is fast to process.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR

Rob Edwards, redwards@utmem.edu

=head1 CONTRIBUTORS

Heikki Lehvaslaiho, heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Restriction::IO::itype2;

use strict;

use Bio::Restriction::Enzyme;
use Bio::Restriction::EnzymeCollection;

use Data::Dumper;

use base qw(Bio::Restriction::IO::base);



sub new {
    my($class, @args) = @_;
    my $self = bless {}, $class;

    $self->_initialize(@args);
    return $self;
}

sub _initialize {
    my($self,@args) = @_;

    my ($verbose) =
            $self->_rearrange([qw(
                                  VERBOSE
                                 )], @args);
    $verbose || 0;
    $self->verbose($verbose);

    $self->_companies;
    return unless $self->SUPER::_initialize(@args);
}

=head2 read

 Title   : read
 Usage   : $renzs = $stream->read
 Function: reads all the restrction enzymes from the stream
 Returns : a Bio::Restriction::Restriction object
 Args    : none

Internally creates a hash of enzyme information which is passed on to
_create_enzyme method. See
L<Bio::Restriction::IO::base::_create_enzyme>.

=cut

sub read {
    my $self = shift;

    my $renzs = Bio::Restriction::EnzymeCollection->new(-empty => 1);

    # read until start of data
    while (defined( my $line = $self->_readline()) ) {
        next if $line =~ /^[ R]/;
        $self->_pushback($line);
        last;
    }

    # enzyme name [tab] prototype [tab] recognition sequence with
    # cleavage site [tab] methylation site and type [tab] commercial
    # source [tab] references

    while (defined(my $line = $self->_readline()) ) {
        $self->debug($line);
        chomp $line;

        my ($name, $prototype, $site, $meth, $vendor, $refs) = split /\t/, $line;
        # we need mininum name and site
        unless ($site) {
            $self->warn("Can not parse line with name [$name]") if $self->verbose > 0;
            next;
        }
        next unless $name;

        # four cut enzymes are not in this format
        my $precut;
        if ($site =~ m/^\((\d+\/\d+)\)[ATGCN]+/) {
            $precut=$1;
            $site =~ s/\($precut\)//;
        }
        # -------------- cut ---------------

        my @sequences;
        if ($site =~ /\,/) {
            @sequences = split /\,/, $site;
            $site=shift @sequences;
        }

        my ($cut, $comp_cut);
        ($site, $cut, $comp_cut) = $self->_cuts_from_site($site);


        my $re = Bio::Restriction::Enzyme->new(-name=>$name,
                                              -site => $site
                                             );
        $renzs->enzymes($re);

        if ($cut) {
            $re->cut($self->_coordinate_shift_to_cut(length($site), $cut));
            $re->complementary_cut($self->_coordinate_shift_to_cut(length($site), $comp_cut));
        }


        #
        # prototype 
        #

        $prototype ? $re->prototype_name($prototype) : $re->is_prototype;


        #
        # methylation
        #
        my @meths;
        if ($meth) {
            # this can be either X(Y) or X(Y),X2(Y2)
            # where X is the base and y is the type of methylation
            if ( $meth =~ /(\S+)\((\d+)\),(\S+)\((\d+)\)/ ) { # two msites per site
                #my ($p1, $m1, $p2, $m2) = ($1, $2, $3, $4);
                $re->methylation_sites($self->_meth($re,$1, $2),
                                       $self->_meth($re,$3,$4));
            }
            elsif ($meth =~ /(\S+)\((\d+)\)/ ) { # one msite per site or more sites
                #print Dumper $meth;
                $re->methylation_sites( $self->_meth($re,$1,$2) );
                @meths = split /, /, $meth;
                $meth=shift @meths;
            } else {
                $self->warn("Unknown methylation format [$meth]") if $self->verbose >0;
            }
        }

        #
        # vendors
        #
        if ($vendor) {
            $re->vendors( split / */, $vendor);
        }

        #
        # references
        #
        $re->references(map {split /\n+/} $refs) if $refs;

        #
        # create special types of Enzymes
        #
        $self->_make_multisites($renzs, $re, \@sequences, \@meths) if @sequences;
        $self->_make_multicuts($renzs, $re, $precut) if $precut;

    }

    return $renzs;
}



=head2 write

 Title   : write
 Usage   : $stream->write($renzs)
 Function: writes restriction enzymes into the stream
 Returns : 1 for success and 0 for error
 Args    : a Bio::Restriction::Enzyme
           or a Bio::Restriction::EnzymeCollection object

=cut

sub write {
    my ($self,@h) = @_;
    $self->throw_not_implemented;
}

1;
