# $Id$
#
# BioPerl module for Bio::Factory::FTLocationFactory
#
# Cared for by Hilmar Lapp <hlapp at gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

#
# (c) Hilmar Lapp, hlapp at gnf.org, 2002.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2002.
#
# You may distribute this module under the same terms as perl itself.
# Refer to the Perl Artistic License (see the license accompanying this
# software package, or see http://www.perl.com/language/misc/Artistic.html)
# for the terms under which you may use, modify, and redistribute this module.
# 
# THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#

# POD documentation - main docs before the code

=head1 NAME

Bio::Factory::FTLocationFactory - DESCRIPTION of Object

=head1 SYNOPSIS

    # parse a string into a location object
    $loc = Bio::Factory::FTLocationFactory->from_string("join(100..200, 400..500);
=head1 DESCRIPTION

Implementation of string-encoded location parsing for the Genbank feature table
encoding of locations.

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

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Factory::FTLocationFactory;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Factory::LocationFactoryI;
use Bio::Location::Simple;
use Bio::Location::Split;
use Bio::Location::Fuzzy;


@ISA = qw(Bio::Root::Root Bio::Factory::LocationFactoryI);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Factory::FTLocationFactory();
 Function: Builds a new Bio::Factory::FTLocationFactory object 
 Returns : an instance of Bio::Factory::FTLocationFactory
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  return $self;
}

=head2 from_string

 Title   : from_string
 Usage   : $loc = $locfactory->from_string("100..200");
 Function: Parses the given string and returns a Bio::LocationI implementing
           object representing the location encoded by the string.

           This implementation parses the Genbank feature table
           encoding of locations.
 Example :
 Returns : A Bio::LocationI implementing object.
 Args    : A string.


=cut

sub from_string{
    # the third parameter is purely optional and indicates a recursive
    # call if set
    my ($self,$locstr,$is_rec) = @_;
    my $loc;

    # there is no place in FT-formatted location strings where whitespace 
    # carries meaning, so strip it off entirely upfront
    $locstr =~ s/\s+//g if ! $is_rec;

    # does it contain an operator?
    if($locstr =~ /^([A-Za-z]+)\((.*)\)$/) {
	# yes:
	my $op = $1;
	my $oparg = $2;
	if($op eq "complement") {
	    # parse the argument recursively, then set the strand to -1
	    $loc = $self->from_string($oparg, 1);
	    $loc->strand(-1);
	} elsif(($op eq "join") || ($op eq "order")) {
	    # This is a split location. Split into components and parse each
	    # one recursively, then gather into a SplitLocationI instance.
	    #
	    # Note: The following code will /not/ work with nested
	    # joins (you want to have grammar-based parsing for that).
	    $loc = Bio::Location::Split->new(-splittype => $op);
	    foreach my $substr (split(/,/, $oparg)) {
		$loc->add_sub_Location($self->from_string($substr, 1));
	    }
	} else {
	    $self->throw("operator \"$op\" unrecognized by parser");
	}
    } else {
	# no operator, parse away
	$loc = $self->_parse_location($locstr);
    }
    return $loc;
}

=head2 _parse_location

 Title   : _parse_location
 Usage   : $loc = $locfactory->_parse_location( $loc_string)

 Function: Parses the given location string and returns a location object 
           with start() and end() and strand() set appropriately.
           Note that this method is private.
 Returns : A Bio::LocationI implementing object or undef on failure
 Args    : location string

=cut

sub _parse_location {
    my ($self, $locstr) = @_;
    my ($loc, $seqid);

    $self->debug( "Location parse, processing $locstr\n");

    # 'remote' location?
    if($locstr =~ /^(\S+):(.*)$/) {
	# yes; memorize remote ID and strip from location string
	$seqid = $1;
	$locstr = $2;
    }

    # split into start and end
    my ($start, $end) = split(/\.\./, $locstr);
    # remove enclosing parentheses if any
    $start =~ s/\((.*)\)/$1/;
    $end =~ s/\((.*)\)/$1/;

    # Is this a simple (exact) or a fuzzy location? Simples have exact start
    # and end, or is between two adjacent bases. Everything else is fuzzy.
    my $loctype = ".."; # exact with start and end as default
    my $locclass = "Bio::Location::Simple";
    if(! defined($end)) {
	if($locstr =~ /(\d+)([\.\^])(\d+)/) {
	    $start = $1;
	    $end = $3;
	    $loctype = $2;
	    $locclass = "Bio::Location::Fuzzy"
		unless (($end-1) == $start) && ($loctype eq "^");
	} else {
	    $end = $start;
	}
    }
    if ( ($start =~ /[\>\<\?\.\^]/) || ($end   =~ /[\>\<\?\.\^]/) ) {
	$locclass = 'Bio::Location::Fuzzy';
    } 

    # instantiate location and initialize
    $loc = $locclass->new(-start => $start, -end  => $end, -strand => 1,
			  -location_type => $loctype);
    # set remote ID if remote location
    if($seqid) {
	$loc->seq_id($seqid);
	$loc->is_remote(1);
    }

    # done (hopefully)
    return $loc;
    
}

1;
