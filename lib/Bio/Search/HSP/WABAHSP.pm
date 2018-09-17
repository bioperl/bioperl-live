#
# BioPerl module for Bio::Search::HSP::WABAHSP
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Search::HSP::WABAHSP - HSP object suitable for describing WABA alignments

=head1 SYNOPSIS

# use this object as you would a GenericHSP
# a few other methods have been added including state

=head1 DESCRIPTION

This object implements a few of the extra methods such as
hmmstate_string which returns the HMM state representation for the
WABA alignment.  We also must implement a method to calculate
homology_string since it is not returned by the algorithm in the
machine readable format.

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

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::Search::HSP::WABAHSP;
use strict;
use Bio::Root::RootI;

use base qw(Bio::Search::HSP::GenericHSP);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::HSP::WABAHSP->new();
 Function: Builds a new Bio::Search::HSP::WABAHSP object 
 Returns : Bio::Search::HSP::WABAHSP
 Args    : -hmmstate_seq => the string representing the state output from WABA

=cut

sub new {
  my($class,@args) = @_;

  # gotta do some preprocessing before we send the arguments to the superclass
  my ($len,$qs,$hs) = Bio::Root::RootI->_rearrange([qw(HSP_LENGTH
						       QUERY_SEQ 
						       HIT_SEQ)],@args);  
  if( $len != length($qs) ) {
    Bio::Root::RootI->warn("HSP_LENGTH must equal length of query_seq string, using value from QUERY_SEQ\n");
      $len = length($qs);
  }
  my( $homol_seq,$gapct,$identical) = ('',0,0);
  
  for(my $i=0;$i<$len;$i++) {
      my $q = substr($qs,$i,1);
      my $h = substr($hs,$i,1);
      if( $q eq '-' || $h eq '-' ) {
	  $homol_seq .= ' ';
	  $gapct ++;
      } elsif( $q eq $h ) { 
	  $homol_seq .= '|';
	  $identical++;
      } else { 
	  $homol_seq .= ' ';
      }
  }
  my $self = $class->SUPER::new('-conserved' => $identical,
				'-identical' => $identical,
				'-gaps'      => $gapct,
				'-homology_seq' => $homol_seq,
				@args);
    
  my ($hmmst) = $self->_rearrange([qw(HMMSTATE_SEQ)],@args);
  defined $hmmst && $self->hmmstate_string($hmmst);
  
  $self->add_tag_value('Target' , join(" ","Sequence:".$self->hit->seq_id, 
				       $self->hit->start, $self->hit->end));

  return $self;
}

=head2 hmmstate_string

 Title   : hmmstate_string
 Usage   : my $hmmseq = $wabahsp->hmmstate_string();
 Function: Get/Set the WABA HMM stateseq
 Returns : string
 Args    : [optional] string


=cut

sub hmmstate_string{
   my ($self,$val) = @_;
   if( defined $val ) { 
       $self->{'_hmmstate_string'} = $val;
   }
   return $self->{'_hmmstate_string'};
}

=head2 homology_string

 Title   : homolgy_string
 Usage   : my $homology_str = $hsp->homology_string();
 Function: Homology string must be calculated for a WABA HSP so we can do
           so here and cache the result so it is only done once
 Returns : string
 Args    : none


=cut

sub homology_string{
   my ($self) = @_;
   return '';
}


1;
