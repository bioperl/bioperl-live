#
# BioPerl module for Bio::Search::HSP::PSLHSP
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::HSP::PSLHSP - A HSP for PSL output

=head1 SYNOPSIS

  # get a PSLHSP somehow (SearchIO::psl)

=head1 DESCRIPTION

This is a HSP for PSL output so we can handle seq_inds differently.

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

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::HSP::PSLHSP;
use strict;

# Object preamble - inherits from Bio::Root::Root


use base qw(Bio::Search::HSP::GenericHSP);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::HSP::PSLHSP->new();
 Function: Builds a new Bio::Search::HSP::PSLHSP object 
 Returns : an instance of Bio::Search::HSP::PSLHSP
 Args    : -gapblocks => arrayref of gap locations which are [start,length]
                         of gaps


=cut

sub new { 
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($qgaplocs,
	$hgaplocs,
	$mismatches) = $self->_rearrange([qw(QUERY_GAPBLOCKS
					     HIT_GAPBLOCKS
					     MISMATCHES)],
				       @args);
    $self->gap_blocks('query',$qgaplocs) if defined $qgaplocs;
    $self->gap_blocks('hit',  $hgaplocs) if defined $hgaplocs;
    $self->mismatches($mismatches) if defined $mismatches;
    return $self;
}

=head2 gap_blocks

 Title   : gap_blocks
 Usage   : $obj->gap_blocks($seqtype,$blocks)
 Function: Get/Set the gap blocks
 Returns : value of gap_blocks (a scalar)
 Args    : sequence type - 'query' or 'hit'
           blocks - arrayref of block start,length


=cut

sub gap_blocks {
    my ($self,$seqtype,$blocks) = @_;
    if( ! defined $seqtype ) { $seqtype = 'query' }
    $seqtype = lc($seqtype);
    $seqtype = 'hit' if $seqtype eq 'sbjct';
    if( $seqtype !~ /query|hit/i ) { 
	$self->warn("Expect either 'query' or 'hit' as argument 1 for gap_blocks");
    }

    unless( defined $blocks ) {
	return $self->{'_gap_blocks'}->{$seqtype};
    } else { 
	return $self->{'_gap_blocks'}->{$seqtype} = $blocks;
    }
}

=head2 mismatches

 Title   : mismatches
 Usage   : $obj->mismatches($newval)
 Function: Get/Set the number of mismatches
 Returns : value of mismatches (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub mismatches{
    my $self = shift;
    return $self->{'mismatches'} = shift if @_;
    return $self->{'mismatches'};
}

1;
