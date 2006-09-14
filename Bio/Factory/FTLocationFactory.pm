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

Bio::Factory::FTLocationFactory - A FeatureTable Location Parser

=head1 SYNOPSIS

    # parse a string into a location object
    $loc = Bio::Factory::FTLocationFactory->from_string("join(100..200, 
                                                         400..500");

=head1 DESCRIPTION

Implementation of string-encoded location parsing for the Genbank feature
table encoding of locations.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net

=head1 CONTRIBUTORS

Jason Stajich, jason-at-bioperl-dot-org

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
    if($locstr =~ /^(\w+)\((.*)\)$/) {  
        # yes:
        my $op = lc($1);    
        my $oparg = $2;
        if($op eq "complement") {
            # parse the argument recursively, then set the strand to -1
            $loc = $self->from_string($oparg, 1);
            $loc->strand(-1);
        } elsif($op eq "join" || $op eq "order" || $op eq "bond" ) {
            # This is a split location. Split into components and parse each
            # one recursively, then gather into a SplitLocationI instance.
            #
            # Note: The following code will /not/ work with nested
            # joins (you want to have grammar-based parsing for that).
            $loc = Bio::Location::Split->new(-verbose   => $self->verbose,
                                                        -splittype => $op);

            # have to do this to capture nested joins, something like this
            # join(11..21,join(100..300,complement(150..230)))
            # This fixes bug #1674
            my $re;
            $re = qr{
                        \(                 # <<--- paired parens required
                        (?:
                        (?> [^()]+ )    # Non-parens without backtracking
                        |
                        (??{ $re })     # Group with matching parens
                      )*
                        \)                 # ---->> paired parens required
                    }x;
            my $oparg_orig = $oparg;
            my @sections;
            # lets capture and remove all the sections which
            # are groups        
            while( $oparg =~ s/(join|complement|bond|order)$re//ig ) {
                # oh man this is SUUUCCCH a hack
                # I don't know what else to do though
                # s// seems to be dropping the whole 
                # warn("rematch is $&   $` $'\n");
                # the code use to just be this line
                push @sections, $&;
                # but I recognized join(...,complement(join(..)))
                # was failling
                my $before = $`;
                if( $oparg ne $before . $') { #'
                    $oparg = $before . $'; # '
                }
            }
            push @sections, split(/,/,$oparg) if length($oparg);
            # because we don't necessarily process the string in-order
            # as we are pulling the data from the string out for
            # groups first, then pulling out data, comma delimited
            # I am re-sorting the sections based on their position
            # in the original string, using the index function to figure
            # out their position in the string
            # --jason
            # resort based on input order, schwartzian style!
            @sections = map { shift @$_ } sort { $a->[1] <=> $b->[1] }
                        map { [$_, index($oparg_orig, $_)] } @sections;
            # end of fix for bug #1674
            foreach my $s (@sections) {
                next unless length($s);     
                $loc->add_sub_Location($self->from_string($s, 1));
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
    
    #$self->debug( "Location parse, processing $locstr\n");
    # 'remote' location?
    if($locstr =~ /^(\S+):(.*)$/) {
        # yes; memorize remote ID and strip from location string
        $seqid = $1;
        $locstr = $2;
    }
    
    # split into start and end
    my ($start, $end) = split(/\.\./, $locstr);
    # remove enclosing parentheses if any; note that because of parentheses
    # possibly surrounding the entire location the parentheses around start
    # and/or may be asymmetrical
    $start =~ s/^\(+//;
    $start =~ s/\)+$//;
    $end   =~ s/^\(+// if $end;
    $end   =~ s/\)+$// if $end;

    # Is this a simple (exact) or a fuzzy location? Simples have exact start
    # and end, or is between two adjacent bases. Everything else is fuzzy.
    my $loctype = ".."; # exact with start and end as default

    $loctype = '?' if ( ($locstr =~ /\?/) &&
                              ($locstr !~ /\?\d+/) );
    $loctype = '?' if ( ($locstr =~ /\?/) &&
                              ($locstr !~ /\?\d+/) );

    my $locclass = "Bio::Location::Simple";
    if(! defined($end)) {
        if($locstr =~ /(\d+)([\.\^])(\d+)/) {
            $start = $1;
            $end = $3;
            $loctype = $2;
            $locclass = "Bio::Location::Fuzzy"
              unless (abs($end-$start) <= 1) && ($loctype eq "^");
        } else {
            $end = $start;
        }
    }
    # start_num and end_num are for the numeric only versions of 
    # start and end so they can be compared
    # in a few lines
    my ($start_num, $end_num) = ($start,$end);
    if ( ($start =~ /[\>\<\?\.\^]/) || ($end   =~ /[\>\<\?\.\^]/) ) {
        $locclass = 'Bio::Location::Fuzzy';
        if($start =~ /(\d+)/) {
            ($start_num) = $1;
        } else { 
            $start_num = 0
        }
        if ($end =~ /(\d+)/) {
            ($end_num)   = $1;
        } else { $end_num = 0 }
    } 
    my $strand = 1;

    if( $start_num > $end_num && $loctype ne '?') {
        ($start,$end,$strand) = ($end,$start,-1);
    }
    # instantiate location and initialize
    $loc = $locclass->new(-verbose => $self->verbose,
                                 -start   => $start, 
                                 -end     => $end, 
                                 -strand  => $strand, 
                                 -location_type => $loctype);
    # set remote ID if remote location
    if($seqid) {
        $loc->is_remote(1);
        $loc->seq_id($seqid);
    }

    # done (hopefully)
    return $loc;
}

1;
