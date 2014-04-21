#
# BioPerl module for Bio::Factory::FTLocationFactory
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net

=head1 CONTRIBUTORS

Jason Stajich, jason-at-bioperl-dot-org
Chris Fields, cjfields-at-uiuc-dot-edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Factory::FTLocationFactory;
use vars qw($LOCREG);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Location::Simple;
use Bio::Location::Split;
use Bio::Location::Fuzzy;


use base qw(Bio::Root::Root Bio::Factory::LocationFactoryI);

BEGIN {
    # the below is an optimized regex obj. from J. Freidl's Mastering Reg Exp.
    $LOCREG = qr{
                (?>
                [^()]+
                |
                \(
                (??{$LOCREG})
                \)
                )*
                }x;     
}

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Factory::FTLocationFactory->new();
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

sub from_string {
    my ($self,$locstr,$op) = @_;
    my $loc;
    
    #$self->debug("$locstr\n");
    
    # $op for operator (error handling)
    
    # run on first pass only
    # Note : These location types are now deprecated in GenBank (Oct. 2006)
    if (!defined($op)) {
        # convert all (X.Y) to [X.Y]
        $locstr =~ s{\((\d+\.\d+)\)}{\[$1\]}g;
        # convert ABC123:(X..Y) to ABC123:[X..Y]
        # we should never see the above
        $locstr =~ s{:\((\d+\.{2}\d+)\)}{:\[$1\]}g;
    }
    
    if ($locstr =~ m{(.*?)\(($LOCREG)\)(.*)}o) { # any matching parentheses?
        my ($beg, $mid, $end) = ($1, $2, $3);
        my (@sublocs) = (split(q(,),$beg), $mid, split(q(,),$end));
        
        my @loc_objs;
        my $loc_obj;
        my @gl_subloc_strands;
        
        SUBLOCS:
        while (@sublocs) {
            my $subloc = shift @sublocs;
            next if !$subloc;
            my $oparg = ($subloc eq 'join'   || $subloc eq 'bond' ||
                         $subloc eq 'order'  || $subloc eq 'complement') ? $subloc : undef;
            # has operator, requires further work (recurse)
            if ($oparg) {
                my $sub = shift @sublocs;
                # simple split operators (no recursive calls needed)
                if (($oparg eq 'join' || $oparg eq 'order' || $oparg eq 'bond' )
                     && $sub !~ m{(?:join|order|bond)}) {
                    my @splitlocs = split(q(,), $sub);
                    $loc_obj = Bio::Location::Split->new(-verbose => 1,
                                                         -splittype => $oparg);
                    # Store strand values for later consistency check
                    my @subloc_strands;
                    my @s_objs;
                    foreach my $splitloc (@splitlocs) {
                        next unless $splitloc;
                        my $sobj;
                        if ($splitloc =~ m{\(($LOCREG)\)}) {
                            my $comploc = $1;
                            $sobj = $self->_parse_location($comploc);
                            $sobj->strand(-1);
                            push @subloc_strands,    -1;
                            push @gl_subloc_strands, -1;
                        } else {
                            $sobj = $self->_parse_location($splitloc);
                            push @subloc_strands,    1;
                            push @gl_subloc_strands, 1;
                        }
                        push @s_objs, $sobj;
                    }

                    # Sublocations strand values consistency check to set
                    # Guide Strand and sublocations adding order
                    if (scalar @s_objs > 0) {
                        my $identical    = 0;
                        my $gl_identical = 0;

                        my $first_value = $subloc_strands[0];
                        foreach my $strand (@subloc_strands) {
                            $identical++ if ($strand == $first_value);
                        }

                        my $first_gl_value = $gl_subloc_strands[0];
                        foreach my $gl_strand (@gl_subloc_strands) {
                            $gl_identical++ if ($gl_strand == $first_gl_value);
                        }

                        if ($identical == scalar @subloc_strands) {
                            # Set guide_strand if all sublocations have the same strand
                            $loc_obj->guide_strand($first_value);

                            # Reverse sublocation order for negative strand locations in cases like this:
                            # join(1..11,join(complement(40..50),complement(60..70)))
                            # But not this:
                            # join(complement(10..20),complement(30..40))
                            if (    $gl_identical != scalar @gl_subloc_strands
                                and $first_value  == -1
                                ) {
                                @s_objs = reverse @s_objs;
                            }
                        }
                        else {
                            # Mixed strand values
                            $loc_obj->guide_strand(undef);
                        }

                        # Add sublocations
                        foreach my $s_obj (@s_objs) {
                            $loc_obj->add_sub_Location($s_obj);
                        }
                    }
                } else {
                    $loc_obj = $self->from_string($sub, $oparg);
                    # reinsure the operator is set correctly for this level
                    # unless it is complement
                    $loc_obj->splittype($oparg) unless $oparg eq 'complement';
                }
            }
            # no operator, simple or fuzzy 
            else {
                $loc_obj = $self->from_string($subloc,1);
            }
            if ($op && $op eq 'complement') {
                $loc_obj->strand(-1);
                push @gl_subloc_strands, -1;
            }
            else {
                push @gl_subloc_strands, 1;
            }

            push @loc_objs, $loc_obj;
        }
        my $ct = @loc_objs;
        if ($op && !($op eq 'join' || $op eq 'order' || $op eq 'bond')
                && $ct > 1 ) {
            $self->throw("Bad operator $op: had multiple locations ".
                         scalar(@loc_objs).", should be SplitLocationI");
        }
        if ($ct > 1) {
            $loc = Bio::Location::Split->new();
            $loc->add_sub_Location(shift @loc_objs) while (@loc_objs);
            return $loc;
        } else {
            $loc = shift @loc_objs;
            return $loc;
        }
    } else { # simple location(s)
        $loc = $self->_parse_location($locstr);
        $loc->strand(-1) if ($op && $op eq 'complement');
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
    if($locstr =~ m{^(\S+):(.*)$}o) {
        # yes; memorize remote ID and strip from location string
        $seqid = $1;
        $locstr = $2;
    }
    
    # split into start and end
    my ($start, $end) = split(/\.\./, $locstr);
    # remove enclosing parentheses if any; note that because of parentheses
    # possibly surrounding the entire location the parentheses around start
    # and/or may be asymmetrical
    # Note: these are from X.Y fuzzy locations, which are deprecated!
    $start =~ s/(?:^\[+|\]+$)//g if $start;
    $end   =~ s/(?:^\[+|\]+$)//g if $end;

    # Is this a simple (exact) or a fuzzy location? Simples have exact start
    # and end, or is between two adjacent bases. Everything else is fuzzy.
    my $loctype = ".."; # exact with start and end as default

    $loctype = '?' if ( ($locstr =~ /\?/) && ($locstr !~ /\?\d+/) );

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
