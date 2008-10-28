# $Id$
#
# BioPerl module for Bio::Search::HSP::FastaHSP
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::HSP::FastaHSP - HSP object for FASTA specific data

=head1 SYNOPSIS

  # get a FastaHSP from a SearchIO stream
  my $in = Bio::SearchIO->new(-format => 'fasta', -file => 'filename.fasta');

  while( my $r = $in->next_result) {
      while( my $hit = $r->next_result ) {
           while( my $hsp = $hit->next_hsp ) {
              print "smith-waterman score (if available): ", 
                    $hsp->sw_score(),"\n";
           }
      }
  }

=head1 DESCRIPTION

Describe the object here

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

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Search::HSP::FastaHSP;
use strict;


use base qw(Bio::Search::HSP::GenericHSP);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::HSP::FastaHSP->new();
 Function: Builds a new Bio::Search::HSP::FastaHSP object 
 Returns : Bio::Search::HSP::FastaHSP
 Args    : -swscore => smith-waterman score

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  
  my ($swscore) = $self->_rearrange([qw(SWSCORE)], @args);

  defined $swscore && $self->sw_score($swscore);

  return $self;
}


=head2 sw_score

 Title   : sw_score
 Usage   : $obj->sw_score($newval)
 Function: Get/Set Smith-Waterman score
 Returns : value of sw_score
 Args    : newvalue (optional)


=cut

sub sw_score{
    my ($self,$value) = @_;
    if( defined $value || ! defined $self->{'_sw_score'} ) {
	$value = 0 unless defined $value; # default value
	$self->{'_sw_score'} = $value;
    }
    return $self->{'_sw_score'};
}


sub get_aln {
    my ($self) = @_;
    require Bio::LocatableSeq;
    require Bio::SimpleAlign;
    my $aln = Bio::SimpleAlign->new();
    my $hs = $self->hit_string();
    my $qs = $self->query_string();

    # fasta reports some extra 'regional' sequence information
    # we need to clear out first
    # this seemed a bit insane to me at first, but it appears to 
    # work --jason
    
    # we infer the end of the regional sequence where the first
    # non space is in the homology string
    # then we use the HSP->length to tell us how far to read
    # to cut off the end of the sequence
        
    my ($start) = 0;
    if( $self->homology_string() =~ /^(\s+)/ ) {
	$start = CORE::length($1);
    }
    $self->debug("hs seq is '$hs'\n");
    $self->debug("qs seq is '$qs'\n");

    $hs = substr($hs, $start,$self->length('total'));
    $qs = substr($qs, $start,$self->length('total'));
    foreach my $seq ( $qs,$hs)  {
	foreach my $f ( '\\', '/', ' ') {
	    my $index =  index($seq,$f);
	    while( $index >=0 && length($seq) > 0 ) {
		substr($hs,$index,1) = '';
		substr($qs,$index,1) = '';
		$self->debug( "$f, $index+1, for ".length($seq). " ($seq)\n");
		$index = index($seq,$f,$index+1);
	    }
	}
    }

    my $seqonly = $qs;
    $seqonly =~ s/\s+//g;
    my ($q_nm,$s_nm) = ($self->query->seq_id(),
			$self->hit->seq_id());
    unless( defined $q_nm && CORE::length ($q_nm) ) {
	$q_nm = 'query';
    }
    unless( defined $s_nm && CORE::length ($s_nm) ) {
	$s_nm = 'hit';
    }
    my $query = Bio::LocatableSeq->new('-seq'   => $seqonly,
				      '-id'    => $q_nm,
				      '-start' => $self->query->start,
				      '-end'   => $self->query->end,
                      '-mapping' => [1, $self->{_query_mapping}]                      
				      );
    $seqonly = $hs;
    $seqonly =~ s/\s+//g;
    my $hit =  Bio::LocatableSeq->new('-seq'    => $seqonly,
				      '-id'    => $s_nm,
				      '-start' => $self->hit->start,
				      '-end'   => $self->hit->end,
                      '-mapping' => [1, $self->{_hit_mapping}]
				      );
    $aln->add_seq($query);
    $aln->add_seq($hit);
    return $aln;
}


1;
