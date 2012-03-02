#
# BioPerl module for Bio::SeqFeature::Primer
#
# This is the original copyright statement. I have relied on Chad's module
# extensively for this module.
#
# Copyright (c) 1997-2001 bioperl, Chad Matsalla. All Rights Reserved.
#           This module is free software; you can redistribute it and/or
#           modify it under the same terms as Perl itself. 
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code
#
# But I have modified lots of it, so I guess I should add:
#
# Copyright (c) 2003 bioperl, Rob Edwards. All Rights Reserved.
#           This module is free software; you can redistribute it and/or
#           modify it under the same terms as Perl itself. 
#
# Copyright Rob Edwards
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Primer - Primer Generic SeqFeature

=head1 SYNOPSIS

 # Set up a single primer that can be used in a PCR reaction

 use Bio::SeqFeature::Primer;

 # Initiate a primer with raw sequence
 my $primer = Bio::SeqFeature::Primer->new( -seq => 'CTTTTCATTCTGACTGCAACG' );

 # Get the primary tag for the primer. This is always 'Primer'.
 my $tag = $primer->primary_tag;

 # Get or set the start and end of the primer match on the template
 $primer->start(2);
 my $start = $primer->start;
 $primer->end(19);
 my $end = $primer->end;

 # Get or set the strand of the primer. Strand should be 1, 0, or -1
 $primer->strand(-1);
 my $strand = $primer->strand;

 # Get or set the ID of the primer
 $primer->display_id('test_id');
 my $id = $primer->display_id;

 # Calculate the Tm (melting temperature) for the primer
 my $tm = $primer->Tm;

 print "These are the details of the primer:\n\tTag:\t\t$tag\n\tLocation\t$location\n\tStart:\t\t$start\n";
 print "\tEnd:\t\t$end\n\tStrand:\t\t$strand\n\tID:\t\t$id\n\tTm:\t\t$tm\n";

=head1 DESCRIPTION

This module handle PCR primer sequences. The L<Bio::SeqFeature::Primer> object
is a L<Bio::SeqFeature::Generic> object that can additionally contain a primer
sequence and its coordinates on a template sequence. The primary_tag() for this
object is 'Primer'. A method is provided to calculate the melting temperature Tm
of the primer. L<Bio::SeqFeature::Primer> objects are useful to build
L<Bio::Seq::PrimedSeq> amplicon objects such as the ones returned by
L<Bio::Tools::Primer3>.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR 

Rob Edwards, redwards@utmem.edu

The original concept and much of the code was written by
Chad Matsalla, bioinformatics1@dieselwurks.com

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Primer;

use strict;
use Bio::PrimarySeq;
use Bio::Tools::SeqStats;

use base qw(Bio::SeqFeature::Generic);


=head2 new()

 Title   : new()
 Usage   : $primer = Bio::SeqFeature::Primer(-id => 'primerX', -seq => $seq_object);
 Function: Instantiate a new Bio::SeqFeature::Primer object
 Returns : A Bio::SeqFeature::Primer object
 Args    : You must pass either a sequence object (preferable) or a sequence string.
           You can also pass arbitray arguments, e.g. -TARGET => '5,3' which will
           be stored in $primer->{'-TARGET'}

=cut

sub new {
    my ($class, %args) = @_;
    my $self = $class->SUPER::new(%args);

    # Set the display ID
    my $id = delete $args{'-id'} || 'SeqFeature Primer object';
    ####$self->SUPER::display_name($id);

    # Set the primer sequence
    my $seq = delete $args{'-seq'} || delete $args{'-sequence'} ||
        $self->throw("Need to provide a sequence during Primer object construction\n");
    if (not ref $seq) {
         $seq = Bio::PrimarySeq->new( -seq => $seq, -id => $id );
    } else {
        if (not $seq->isa('Bio::PrimarySeqI')) {
            $self->throw("Expected a sequence object but got a [".ref($seq)."]\n");
        }
    }
    $self->{seq} = $seq;
    
    # Save arbitrary parameters like:
    #   TARGET=513,26
    #   PRIMER_FIRST_BASE_INDEX=1
    #   PRIMER_LEFT=484,20
    while (my ($arg, $val) = each %args) {
        $self->{$arg} = $val;
    }

    $self->primary_tag('Primer');

    return $self;
}


=head2 seq()

 Title   : seq()
 Usage   : $seq = $primer->seq();
 Function: Get the sequence object of this Primer.
 Returns : A Bio::PrimarySeq object
 Args    : None.

=cut

sub seq {
    my $self = shift;
    return $self->{seq};
}


#=head2 primary_tag()

# Title   : primary_tag()
# Usage   : $tag = $primer->primary_tag();
# Function: Get the primary tag associated with this Primer.
# Returns : A string, always 'Primer'.
# Args    : None.

#=cut

##sub primary_tag {
##    return 'Primer';
##}


#=head2 source_tag()

# Title   : source_tag()
# Usage   : $tag = $feature->source_tag();
# Function: Get or set the source tag associated with this Primer.
# Returns : A string.
# Args    : If an argument is provided, the source of this SeqFeature
#           is set to that argument.

#=cut

##sub source_tag {
##    my ($self, $insource) = @_;
##    if ($insource) {
##        $self->{source} = $insource;
##    }
##    return $self->{source};
##}


#=head2 location()

# Title   : location()
# Usage   : $tag = $primer->location();
# Function: Get or set the location of the primer on the template sequence  
# Returns : If the location is set, return that, if not, return 0. 
#           Note: You can set the location to whatever you want, but at the
#           moment I am using the primer3 notation of location. In this form, the
#           location of both primers is given as their start from the 5' end, and
#           a length. In this case, the left primer is given from the leftmost
#           end, but the right primer is given from the rightmost end. You could
#           also specify a location as a L<Bio::Location::Simple> object or other.
#           You can also use start() and end() to get the leftmost and rightmost
#           base of each primer sequence.
# Args    : If supplied will set a location

#=cut

##sub location {
##    my ($self, $location) = @_;
##    if ($location) {
##        $self->{location} = $location;
##    }
##    return $self->{location} || 0;
##}
# Use Bio::SeqFeature::Generic's location() method, but just have an extra conversion
# layer for compatibility with previous usages of location() in Bio::SeqFeature::Primer.

sub location {
    my ($self, $location) = @_;
    if ($location) {
        unless (ref($location) and $location->isa('Bio::LocationI')) {
            # Convert location string for backward compatibility
            $self->warn('Passing a string to location() is deprecated. Try to pass a Bio::Location::Simple object or use start() and end() instead.');
            my ($start, $end) = ($location =~ m/^(\d+).*?(\d*)$/);
            if (not $start) {
                $self->throw("Wrong location entered: $location");
            }
            $end = $start if not $end;
            $location = Bio::Location::Simple->new( -start => $start, -end => $end );
        }
        $self->SUPER::location($location);
    }
    return $self->SUPER::location;
}


#=head2 start()

# Title   : start()
# Usage   : $start_position = $primer->start($new_position);
# Function: Get or set the start position of this Primer on the template.
#           This is the leftmost base, regardless of whether it is a left or
#           right primer, i.e. start() < end().
# Returns : The start position of this primer or 0 if not set.
# Args    : If supplied will set a start position.

#=cut

##sub start {
##    my ($self, $start) = @_;
##    if ($start) {
##        $self->{start_position} = $start;
##    }
##    return $self->{start_position} || 0;
##}


#=head2 end()

# Title   : end()
# Usage   : $end_position = $primer->end($new_position);
# Function: Get or set the end position of this Primer on the template.
#           This is the rightmost base, regardless of whether it is a left or
#           right primer, i.e. start() < end().
# Returns : The end position of this primer or 0 if not set.
# Args    : If supplied will set an end position.

#=cut

##sub end {
##    my ($self, $end) = @_;
##    if ($end) {
##        $self->{end_position} = $end;
##    }
##    return $self->{end_position} || 0;
##}


#=head2 strand()

# Title   : strand()
# Usage   : $strand = $primer->strand();
# Function: Get or set the strand. It can be 1, 0 (not set), or -1.
# Returns : The strand that the primer binds to.
# Args    : If an argument is supplied will set the strand, otherwise will return it.

#=cut

##sub strand {
##    my ($self, $strand) = @_;
##    if ($strand) {
##        unless ($strand == -1 || $strand == 0 || $strand == 1) {
##            $self->throw("Strand must be either 1, 0, or -1, not $strand");
##        }
##        $self->{strand} = $strand;
##    }
##    return $self->{strand} || 0;
##}


#=head2 display_id()

# Title   : display_id()
# Usage   : $id = $primer->display_id($new_id);
# Function: Get or set the display ID for this Primer.
# Returns : A scalar.
# Args    : If an argument is provided, the display_id of this primer is
#           set to that value.

#=cut

##sub display_id {
##    my ($self, $newid) = @_;
##    if ($newid) {
##        $self->seq()->display_id($newid);
##    }
##    return $self->seq()->display_id();
##}


=head2 Tm()

 Title   : Tm()
 Usage   : $tm = $primer->Tm(-salt => 0.05, -oligo => 0.0000001);
 Function: Calculate the Tm (melting temperature) of the primer
 Returns : A scalar containing the Tm.
 Args    : -salt  : set the Na+ concentration on which to base the calculation
                    (default=0.05 molar).
         : -oligo : set the oligo concentration on which to base the
                    calculation (default=0.00000025 molar).
 Notes   : Calculation of Tm as per Allawi et. al Biochemistry 1997
           36:10581-10594. Also see documentation at
           http://www.idtdna.com/Scitools/Scitools.aspx as they use this
           formula and have a couple nice help pages. These Tm values will be
           about are about 0.5-3 degrees off from those of the idtdna web tool.
           I don't know why.

           This was suggested by Barry Moore (thanks!). See the discussion on
           the bioperl-l with the subject "Bio::SeqFeature::Primer Calculating
           the PrimerTM"

=cut

sub Tm {
    my ($self, %args) = @_;
    my $salt_conc = 0.05; # salt concentration (molar units)
    my $oligo_conc = 0.00000025; # oligo concentration (molar units)
    if ($args{'-salt'}) {
        # Accept object defined salt concentration
        $salt_conc = $args{'-salt'};
    } 
    if ($args{'-oligo'}) {
        # Accept object defined oligo concentration
        $oligo_conc = $args{'-oligo'};
    }
    my $seqobj = $self->seq();
    my $length = $seqobj->length();
    my $sequence = uc $seqobj->seq();
    my @dinucleotides;
    my $enthalpy;
    my $entropy;
    # Break sequence string into an array of all possible dinucleotides
    while ($sequence =~ /(.)(?=(.))/g) {
        push @dinucleotides, $1.$2;
    }
    # Build a hash with the thermodynamic values
    my %thermo_values = ('AA' => {'enthalpy' => -7.9,
                                  'entropy'  => -22.2},
                         'AC' => {'enthalpy' => -8.4,
                                  'entropy'  => -22.4},
                         'AG' => {'enthalpy' => -7.8,
                                  'entropy'  => -21},
                         'AT' => {'enthalpy' => -7.2,
                                  'entropy'  => -20.4},
                         'CA' => {'enthalpy' => -8.5,
                                  'entropy'  => -22.7},
                         'CC' => {'enthalpy' => -8,
                                  'entropy'  => -19.9},
                         'CG' => {'enthalpy' => -10.6,
                                  'entropy'  => -27.2},
                         'CT' => {'enthalpy' => -7.8,
                                  'entropy'  => -21},
                         'GA' => {'enthalpy' => -8.2,
                                  'entropy'  => -22.2},
                         'GC' => {'enthalpy' => -9.8,
                                  'entropy'  => -24.4},
                         'GG' => {'enthalpy' => -8,
                                  'entropy'  => -19.9},
                         'GT' => {'enthalpy' => -8.4,
                                  'entropy'  => -22.4},
                         'TA' => {'enthalpy' => -7.2,
                                  'entropy'  => -21.3},
                         'TC' => {'enthalpy' => -8.2,
                                  'entropy'  => -22.2},
                         'TG' => {'enthalpy' => -8.5,
                                  'entropy'  => -22.7},
                         'TT' => {'enthalpy' => -7.9,
                                  'entropy'  => -22.2},
                         'A' =>  {'enthalpy' => 2.3,
                                  'entropy'  => 4.1},
                         'C' =>  {'enthalpy' => 0.1,
                                  'entropy'  => -2.8},
                         'G' =>  {'enthalpy' => 0.1,
                                  'entropy'  => -2.8},
                         'T' =>  {'enthalpy' => 2.3,
                                  'entropy'  => 4.1}
                        );
    # Loop through dinucleotides and calculate cumulative enthalpy and entropy values
    for (@dinucleotides) {
        $enthalpy += $thermo_values{$_}{enthalpy};
        $entropy += $thermo_values{$_}{entropy};
    }
    # Account for initiation parameters
    $enthalpy += $thermo_values{substr($sequence,  0, 1)}{enthalpy};
    $entropy  += $thermo_values{substr($sequence,  0, 1)}{entropy};
    $enthalpy += $thermo_values{substr($sequence, -1, 1)}{enthalpy};
    $entropy  += $thermo_values{substr($sequence, -1, 1)}{entropy};
    # Symmetry correction
    $entropy -= 1.4;
    my $r = 1.987; # molar gas constant
    my $tm = $enthalpy * 1000 / ($entropy + ($r * log($oligo_conc))) - 273.15 + (12* (log($salt_conc)/log(10)));

    return $tm;
 }

=head2 Tm_estimate

 Title   : Tm_estimate
 Usage   : $tm = $primer->Tm_estimate(-salt => 0.05);
 Function: Estimate the Tm (melting temperature) of the primer
 Returns : A scalar containing the Tm.
 Args    : -salt set the Na+ concentration on which to base the calculation.
 Notes   : This is only an estimate of the Tm that is kept in for comparative
           reasons. You should probably use Tm instead!

           This Tm calculations are taken from the Primer3 docs: They are
           based on Bolton and McCarthy, PNAS 84:1390 (1962) 
           as presented in Sambrook, Fritsch and Maniatis,
           Molecular Cloning, p 11.46 (1989, CSHL Press).

           Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) - 600/length

           where [Na+] is the molar sodium concentration, %GC is the
           %G+C of the sequence, and length is the length of the sequence.

           However.... I can never get this calculation to give me the same result
           as primer3 does. Don't ask why, I never figured it out. But I did 
           want to include a Tm calculation here because I use these modules for 
           other things besides reading primer3 output.

           The primer3 calculation is saved as 'PRIMER_LEFT_TM' or 'PRIMER_RIGHT_TM'
           and this calculation is saved as $primer->Tm so you can get both and
           average them!

=cut

sub Tm_estimate {

    # This should probably be put into seqstats as it is more generic, but what the heck.

    my ($self, %args) = @_;
    my $salt = 0.2;
    if ($args{'-salt'}) {
        $salt = $args{'-salt'}
    };
    my $seqobj = $self->seq();
    my $length = $seqobj->length();
    my $seqdata = Bio::Tools::SeqStats->count_monomers($seqobj);
    my $gc=$$seqdata{'G'} + $$seqdata{'C'};
    my $percent_gc = ($gc/$length)*100;

    my $tm = 81.5+(16.6*(log($salt)/log(10)))+(0.41*$percent_gc) - (600/$length);

    # and now error check compared to primer3
    # note that this NEVER gives me the same values, so I am ignoring it
    # you can get these out separately anyway

    #if ($self->{'PRIMER_LEFT_TM'}) {
    # unless ($self->{'PRIMER_LEFT_TM'} == $tm) {
    #  $self->warn("Calculated $tm for Left primer but received ".$self->{'PRIMER_LEFT_TM'}." from primer3\n");
    # }
    #}
    #elsif ($self->{'PRIMER_RIGHT_TM'}) {
    # unless ($self->{'PRIMER_RIGHT_TM'} == $tm) {
    #   $self->warn("Calculated $tm for Right primer but received ".$self->{'PRIMER_RIGHT_TM'}." from primer3\n");
    # }
    #}

    return $tm; 
} 

=head2 primary_tag, source_tag, location, start, end, strand...

The documentation of L<Bio::SeqFeature::Generic> describes all the methods that
L<Bio::SeqFeature::Primer> object inherit.

=cut

1;
