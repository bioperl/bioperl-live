# $Id$
# BioPerl module for Bio::Restriction::IO::withrefm
#
# Cared for by Rob Edwards <redwards@utmem.edu>
#
# Copyright Rob Edwards
#
# You may distribute this module under the same terms as perl itself
#

# POD documentation - main docs before the code

=head1 NAME

Bio::Restriction::IO::withrefm - withrefm enzyme set

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::Restriction::IO class.

=head1 DESCRIPTION

This is the most complete format of the REBASE files, and basically
includes all the data on each of the restriction enzymes.


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

package Bio::Restriction::IO::withrefm;

use vars qw(%WITH_REFM_FIELD);
use strict;

#use Bio::Restriction::IO;
use Bio::Restriction::Enzyme;
use Bio::Restriction::EnzymeCollection;

use Data::Dumper;

use base qw(Bio::Restriction::IO::base);

=head2 read

 Title   : read
 Usage   : $renzs = $stream->read
 Function: reads all the restrction enzymes from the stream
 Returns : a Bio::Restriction::Restriction object
 Args    : none

=cut

sub read {
    my $self = shift;

    my $renzs = Bio::Restriction::EnzymeCollection->new(-empty => 1);

    local $/ = '<1>';
    while (defined(my $entry=$self->_readline()) ) {

        # not an entry.
        next unless $entry =~ /<2>/;

        #$self->debug("|$entry|\n");

        #
        # Minimal information
        #
        my ($name) = $entry =~ /^(\S+)/;
        my ($site) = $entry =~ /\<3\>([^\n]+)/;


        if ( ! defined $site || $site eq '' or $site eq '?') {
            $self->warn("$name: no site. Skipping") if $self->verbose > 1;
            next;
        }

        my $precut;
        if ($site =~ m/^\((\w+\/\w+)\)\w+\((\w+\/\w+)\)/) {
            $precut = $1;
            $site =~ s/\($precut\)//;
        }

        # there are a couple of sequences that have multiple
        # recognition sites eg M.PhiBssHII: ACGCGT,CCGCGG,RGCGCY,RCCGGY,GCGCGC


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
        # prototype / isoschizomers
        #

        my ($isoschizomers) = $entry =~ /<2>([^\n]+)/;

        if ($isoschizomers) {
            # bug 2179
            # here's the trick; if there are no enzymes here, the enzyme in <1>
            # is the prototype (see withref format for this).  However, one
            # can't unequivicably assign prototype based on the presence of
            # enzymes or which one is first without building a logic kit on
            # determining how these are assigned.
            
            # we could add in a hook to check an outside prototype file here
            # at some point; for now we'll just warn if is_prototype() is called
            my @isos = split /\,/, $isoschizomers;
            $re->isoschizomers(@isos);
            #$re->is_prototype(0);
        } else {
            $re->is_prototype(1);
        }

        #
        # methylation
        #

        my ($meth) = $entry =~ /<4>([^\n]+)/;
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
        # microbe
        #
        my ($microbe) = $entry =~ /<5>([^\n]+)/;
        $re->microbe($microbe) if $microbe;

        #
        # source
        #
        my ($source) = $entry =~ /<6>([^\n]+)/;
        $re->source($source) if $source;

        #
        # vendors
        #
        my ($vendors) = $entry =~ /<7>([^\n]+)/;
        $re->vendors(split / */, $vendors) if $vendors;

        #
        # references
        #
        my ($refs) = $entry =~ /<8>(.+)/s;
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
