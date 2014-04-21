# $Id: Match.pm,v 1.2 2007/06/14 18:01:52 nathan Exp $
#
# BioPerl module for Bio::Tools::Match
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Match - Parses output from Transfac's match(TM)

=head1 SYNOPSIS

  use strict;

  use Bio::Tools::Match;

  my $parser = Bio::Tools::Match->new(-file => "match.out");
  
  while (my $feat = $parser->next_result) {
    my $start = $feat->start;
    my $end = $feat->end;
    my $core_score = $feat->score;
    my $matrix_score = ($feat->annotation->get_Annotations('matrix_score'))[0]->value;
    my $matrix_id = ($feat->annotation->get_Annotations('matrix_id'))[0]->value;
  }

=head1 DESCRIPTION

This module is used to parse the output from Transfac's match(TM) program. It
doesn't support the histogram output of match.

Each result is a Bio::SeqFeature::Annotated representing a single matrix match.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list. Your participation is much appreciated.

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

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::Match;
use strict;

use Bio::SeqFeature::Generic;
use Bio::Annotation::SimpleValue;

use base qw(Bio::Root::Root Bio::Root::IO);


=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Match->new();
 Function: Builds a new Bio::Tools::Match object
 Returns : Bio::Tools::Match
 Args    : -file (or -fh) should contain the contents of a standard match output

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    $self->_initialize_io(@args);
    
    return $self;
}

=head2 next_result

 Title   : next_result
 Usage   : $result = $obj->next_result();
 Function: Returns the next result available from the input, or undef if there
           are no more results.
 Returns : Bio::SeqFeature::Annotated object. Features are annotated with tags
           for 'matrix_score', 'matrix_id' and a 'predicted' tag.
 Args    : none

=cut

sub next_result {
    my ($self) = @_;
    
    my $line = $self->_readline || return;
    
    if (! $self->{found_seq_id} && $line =~ /^Inspecting sequence ID\s+(.+)/) {
        $self->{found_seq_id} = $1;
    }
    
    while ($line !~ /^\s\S+\s+\|\s+\d+/) {
        $line = $self->_readline || return;
    }
    
    
    # The first column gives the TRANSFAC(r) identifier of the matching matrix,
    # then comes the position and the strand where the respective match has been
    # found. The core similarity score is given in column three, the matrix
    # similarity score in column four. The last column gives the matching
    # sequence.
    #
    #
    #Search for sites by WeightMatrix library: /home/sendu/files/programs/transfac/cgi-bin/data/matrix.dat
    #Sequence file: sequence.fa
    #Site selection profile: mxprf Profile generated from /home/sendu/files/programs/transfac/cgi-bin/data/matrix.dat with default values.
    #
    #
    #Inspecting sequence ID   Homo_sapiens
    #
    # V$MYOD_01              |        5 (+) |  0.751 |  0.784 | ttaGAGGTggcg
    # V$MYOD_01              |        5 (-) |  0.778 |  0.580 | ttagAGGTGgcg
    # V$MYOD_01              |       30 (+) |  0.751 |  0.581 | gctCAGGCaccc
    #[...]
    # V$RORA_Q4              |    53610 (+) |  0.775 |  0.668 | tgtgggGGCCA
    # V$RORA_Q4              |    53639 (+) |  0.775 |  0.636 | gtcgggGGACA
    #
    # Total sequences length=53654
    #
    # Total number of found sites=1735559
    #
    # Frequency of sites per nucleotide=32.347243
    
    my ($matrix_id, $start, $strand, $core_score, $matrix_score, $seq) = $line =~ /^\s(\S+)\s+\|\s+(\d+)\s+\(([+-])\)\s+\|\s+(\S+)\s+\|\s+(\S+)\s+\|\s+(\S+)/;
    my $feat = Bio::SeqFeature::Generic->new(
        -seq_id => $self->{found_seq_id},
        -start  => $start, 
        -end    => $start + length($seq) - 1,
        -strand => 1,
        -score  => $core_score,
        -source => 'transfac_match');
    
    my $sv = Bio::Annotation::SimpleValue->new(-tagname => 'predicted', -value => 1);
    $feat->annotation->add_Annotation($sv);
    $sv = Bio::Annotation::SimpleValue->new(-tagname => 'matrix_score', -value => $matrix_score);
    $feat->annotation->add_Annotation($sv);
    $sv = Bio::Annotation::SimpleValue->new(-tagname => 'matrix_id', -value => $matrix_id);
    $feat->annotation->add_Annotation($sv);
    
    return $feat;
}

1;
