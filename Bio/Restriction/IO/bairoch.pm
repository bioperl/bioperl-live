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

Bio::Restriction::IO::bairoch - bairoch enzyme set

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::Restriction::IO class.

=head1 DESCRIPTION

This is the most complete format of the REBASE files, and basically
includes all the data on each of the restriction enzymes.

This parser is for the Bairoch format (aka MacVector, Vector NTI, PC/Gene 
(Bairoch) format), REBASE format #19

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists Your participation is much appreciated.

  bioperl-l@bioperl.org                    - General discussion
  http://bio.perl.org/MailList.html        - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR

Rob Edwards, redwards@utmem.edu

=head1 CONTRIBUTORS

Heikki Lehvaslaiho, heikki@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Restriction::IO::bairoch;

use vars qw(@ISA %WITH_REFM_FIELD);
use strict;

#use Bio::Restriction::IO;
use Bio::Restriction::IO::base;
use Bio::Restriction::Enzyme;
use Bio::Restriction::Enzyme::MultiCut;
use Bio::Restriction::Enzyme::MultiSite;
use Bio::Restriction::EnzymeCollection;

use Data::Dumper;

@ISA = qw(Bio::Restriction::IO::base);


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

    return unless $self->SUPER::_initialize(@args);
}

=head2 read

 Title   : read
 Usage   : $renzs = $stream->read
 Function: reads all the restrction enzymes from the stream
 Returns : a Bio::Restriction::Restriction object
 Args    : none

=cut

sub read {
    my $self = shift;

    my $renzs = new Bio::Restriction::EnzymeCollection(-empty => 1);

    local $/ = '//';
    while (defined(my $entry=$self->_readline()) ) {
        $self->debug("|$entry|\n");

        #
        # Minimal information
        #
        my ($name) = $entry =~ /ID\s+(\S+)/;
        my ($site) = $entry =~ /RS\s+([^\n]+)/;
        next unless ($name && $site);
       
        # the standard sequence format for these guys is:
	# GATC, 2;
	# or, for enzymes that cut more than once
	# GATC, 2; GTAC, 2; 

        # there are a couple of sequences that have multiple
        # recognition sites. 


        my @sequences;
        if ($site =~ /\;/) {
            @sequences = split /\;/, $site;
            $site=shift @sequences;
        }
        
	my ($seq, $cut)=split /,\s+/, $site;
	print STDERR "SITE: |$site| GAVE: |$seq| and |$cut|\n";
	if ($seq eq '?') {
	   $self->warn("$name: no site. Skipping") if $self->verbose > 1;
	   next;
	}
	
        # this is mainly an error check to make sure that I am adding what I think I am!	
	if ($seq !~ /[NGATC]/i) {
	  $self->throw("Sequence $name has weird sequence: |$seq|");
	}
	my $re;
	if ($cut eq "?") {
          $re = new Bio::Restriction::Enzyme(-name=>$name, -seq => $seq);
	}
	else {
           if ($cut !~ /^-?\d+$/) {
	     $self->throw("Cut site from $name is weird: |$cut|\n");
           }
	
           $re = new Bio::Restriction::Enzyme(-name=>$name,
                                              -cut => $cut,
                                              -seq => $seq
                                              );
	}
        $renzs->enzymes($re);

        #
        # create special types of Enzymes
        #
        $self->_make_multisites($renzs, $re, \@sequences, \@meths) if @sequences;


        #
        # prototype / isoschizomers
        #
        my ($prototype) = $entry =~ /PT\s+([^\n]+)/;

        if ($prototype) {
            #$re->isoschizomers(split /\,/, $isoschizomers);
            #NOTE: Need to add a method so that we can add isoschosomers to enzymes that may not exist!
	    $re->is_prototype(0);
        } else {
            $re->is_prototype(1);
        }

        #
        # methylation
        #

        my ($meth) = $entry =~ /MS\s+([^\n]+)/;
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
        my ($microbe) = $entry =~ /OS\s+([^\n]+)/;
        $re->microbe($microbe) if $microbe;

        #
        # source
        #
        #my ($source) = $entry =~ /<6>([^\n]+)/;
        #$re->source($source) if $source;

        #
        # vendors
        #
        my ($vendors) = $entry =~ /CR\s+([^\n]+)/;
        $re->vendors(split /,\s*/, $vendors) if $vendors;

        #
        # references
        #
        #my ($refs) = $entry =~ /<8>(.+)/s;
        #$re->references(map {split /\n+/} $refs) if $refs;



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
