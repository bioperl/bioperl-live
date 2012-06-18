=pod

=head1 NAME

Bio::FeatureIO::vecscreen_simple - read/write features from NCBI vecscreen -f 3
output

=head1 SYNOPSIS

    # read features 
    my $fin = Bio::FeatureIO->new(-file=>'vecscreen.out',
                                  -format=>'vecscreen_simple');
    my @vec_regions;
    while (my $f = $fin->next_feature) {
      push @vec_regions, $f;
    }
    
    # write features NOT IMPLEMENTED

=head1 DESCRIPTION

vecscreen is a system for quickly identifying segments of a nucleic
acid sequence that may be of vector origin.  NCBI developed vecscreen
to minimize the incidence and impact of vector contamination in public
sequence databases.  GenBank Annotation Staff use vecscreen to verify
that sequences submitted for inclusion in the database are free from
contaminating vector sequence. Any sequence can be screened for vector
contamination using vecscreen.

This module provides parsing for vecscreen '-f 3' output, described in
the vecscreen documentation as 'Text list, no alignments'

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

 https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Robert Buels

Email rmb32 AT cornell.edu

=head1 CONTRIBUTORS

Based on ptt.pm by Torsten Seeman

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::FeatureIO::vecscreen_simple;

use strict;
use base qw(Bio::FeatureIO);
use Bio::SeqFeature::Generic;

=head2 _initialize

Title   : _initialize
Function: Reading? parses the header of the input
          Writing? 

=cut

sub _initialize {
 my($self,%arg) = @_;

 $self->SUPER::_initialize(%arg);

 if ($self->mode eq 'r') {
   $self->{parse_state}->{seqname}   = '';
   $self->{parse_state}->{matchtype} = '';
 }
 else {
   $self->throw('vecscreen_simple feature writing not implemented');
 }
}

=head2 next_feature

  Title   : next_feature
  Usage   : $io->next_feature()
  Function: read the next feature from the vecscreen output file
  Args    : none
  Returns : Bio::SeqFeatureI object

=cut

sub next_feature {
 my $self = shift;
 return unless $self->mode eq 'r'; # returns if can't read next_feature when we're in write mode

 while ( my $line = $self->_readline() ) {
   chomp $line;
   if ( $line =~ /^>Vector (\S+)/ ) {
     $self->{parse_state}{seqname} = $1;
   } elsif ( $line =~ /^\s*WARNING/ ) {
     $self->warn("$self->{parse_state}{seqname}: vecscreen says: $line\n");
   } elsif ( $line =~ /\S/ ) {

     $self->{parse_state}{seqname}
	or $self->throw("Unexpected line in vecscreen output '$line'");

     #if it's not a vector line, it should be either a match type or
     #a coordinates line
     my $lcline = lc $line;

     if ( $line =~ /^(\d+)\s+(\d+)\s*$/ ) {
	my ($s,$e) = ($1,$2);

	my $matchtype = $self->{parse_state}{matchtype};
	$matchtype =~ s/\s/_/g; #replace whitespace with underscores for the primary tag
	return Bio::SeqFeature::Generic->new( -start    => $s,
					      -end      => $e,
					      -primary  => $matchtype,
					      -seq_id   => $self->{parse_state}{seqname},
					    );
     } elsif ( $lcline eq 'no hits found' ) {
	$self->{parse_state}{seqname} = '';
     } elsif ( grep $lcline eq $_, 'strong match', 'moderate match', 'weak match', 'suspect origin') {
	$self->{parse_state}{matchtype} = $lcline;
     } else {
	$self->throw("Parse error.  Expected a match type or coordinate line but got '$line'");
     }
   } else {
     #blank line, ignore it and reset parser

     $self->{parse_state}{seqname} = ''; #< a line with whitespace
     #indicates a boundary
     #between output for
     #different sequences
     $self->{parse_state}{matchtype} = '';
   }
 }

 return;
}

=head2 write_feature (NOT IMPLEMENTED)

  Title   : write_feature
  Usage   : $io->write_feature($feature)
  Function: write a Bio::SeqFeatureI object in vecscreen -f 3 format
  Example :
  Args    : Bio::SeqFeatureI object
  Returns :

=cut

sub write_feature {
 shift->throw_not_implemented;
}


###
1;#do not remove
###
