# BioPerl module for Bio::Restriction::IO::withrefm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Rob Edwards, redwards@utmem.edu

=head1 CONTRIBUTORS

Heikki Lehvaslaiho, heikki-at-bioperl-dot-org
Mark A. Jensen, maj-at-fortinbras-dot-us

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

        # there are a couple of sequences that have multiple
        # recognition sites eg M.PhiBssHII: ACGCGT,CCGCGG,RGCGCY,RCCGGY,GCGCGC
	# TaqII : GACCGA(11/9),CACCCA(11/9)

        my @sequences;
        if ($site =~ /\,/) {
            @sequences = split (/\,/, $site);
            $site=shift @sequences;
        }

	# this regexp now parses all possible components
	# $1 : (s/t) or undef 
	# $2 : [site]
	# $3 : (m/n) or undef /maj

	no warnings; # avoid faulty 'uninitialized value' warnings
                     # occurring against the variables set by 
	             # regexp matching (unless anyone has other ideas...)

	my ($precut, $recog, $postcut) = ( $site =~ m/^(?:\((-?\w+\/-?\w+)\))?([\w^]+)(?:\((-?\w+\/-?\w+)\))?/ );


        #
        # prototype / isoschizomers
        #

        my ($isoschizomers) = $entry =~ /<2>([^\n]+)/;
	my @isos = split(/\,/,$isoschizomers);
	my $is_prototype = (@isos ? 1 : 0);

        #
        # microbe
        #
        my ($microbe) = $entry =~ /<5>([^\n]+)/;

        #
        # source
        #
        my ($source) = $entry =~ /<6>([^\n]+)/;

        #
        # vendors
        #
        my ($vendors) = $entry =~ /<7>([^\n]+)/;
	my @vendors = split(/ */, $vendors);


        #
        # references
        #
        my ($refs) = $entry =~ /<8>(.+)<1>/s;
	my @refs = map {split /\n+/} $refs;

	use warnings; 
	       
	# when enz is constructed, site() will contain original characters,
	# but recog() will contain a regexp if required.../maj
        my $re = Bio::Restriction::Enzyme->new(
	    -name          => $name,
	    -site          => $recog,
	    -recog         => $recog,
	    -precut        => $precut,
	    -postcut       => $postcut,
	    -is_prototype  => $is_prototype,
	    -isoschizomers => [@isos],
	    -source        => $source,
	    -vendors       => [@vendors],
	    -references    => [@refs],
	    -xln_sub       => \&_xln_sub
	    );

        #
        # methylation: easier to set here during parsing/maj
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
                $re->methylation_sites( $self->_meth($re,$1,$2) );
                @meths = split (/\, /, $meth);
                $meth=shift @meths;
            } else {
                $self->warn("Unknown methylation format [$meth]") if $self->verbose >0;
            }
        }

	# the _make_multicuts function now takes place in the 
	# Enzyme constructor / maj

        #
        # create special types of Enzymes
	# (because of object cloning in _make_multisites, this happens
	#  after everything else is set /maj)
        # (with the removal of the collection from the arglist, this
	# call (or its code) could now be placed in the constructor,
	# which is safer (since this has to happen last), 
	# but it requires the methylation info, which 
	# is more natural to get out here in the parsing./maj

        $self->_make_multisites($re, \@sequences, \@meths, \&_xln_sub) if @sequences;

        $renzs->enzymes($re);


    }

    return $renzs;
}

=head2 _xln_sub

 Title   : _xln_sub
 Function: Translates withrefm coords to Bio::Restriction coords
 Args    : Bio::Restriction::Enzyme object, scalar integer (cut posn)
 Note    : Used internally; pass as a coderef to the B:R::Enzyme 
           constructor

=cut

sub _xln_sub { 
    my ($z,$c) = @_; 
    return ($c < 0 ? $c : length($z->string)+$c);
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
