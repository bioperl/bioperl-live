# BioPerl module for Bio::Tools::Promoterwise
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Shawn Hoon <shawnh@fugu-sg.org>
#
# Copyright Shawn Hoon
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Promoterwise - parser for Promoterwise tab format output

=head1 SYNOPSIS

  use Bio::Tools::Promoterwise;

  my $pw = Bio::Tools::Promoterwise->new(-file=>"out",
                                         -query1_seq=>$seq1,
                                         -query2_seq=>$seq2);
  while (my $fp = $pw->next_result){
    print "Hit Length: ".$fp->feature1->length."\n";
    print "Hit Start: ".$fp->feature1->start."\n";
    print "Hit End: ".$fp->feature1->end."\n";
    print "Hsps: \n";
    my @first_hsp = $fp->feature1->sub_SeqFeature;
    my @second_hsp = $fp->feature2->sub_SeqFeature;
    foreach my $i (0..$#first_hsp){
      print $first_hsp[$i]->start. " ".$first_hsp[$i]->end." ".
            $second_hsp[$i]->start. " ".$second_hsp[$i]->end."\n";
    }
  }

=head1 DESCRIPTION

Promoteriwise is an alignment algorithm that relaxes the constraint
that local alignments have to be co-linear. Otherwise it provides a
similar model to DBA, which is designed for promoter sequence
alignments.  Promoterwise is written by Ewan Birney.  It is part of
the wise2 package available at
L<ftp://ftp.ebi.ac.uk/pub/software/unix/wise2/>

This module is the parser for the Promoterwise output in tab format.

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

=head1 AUTHOR - Shawn Hoon

Email shawnh@fugu-sg.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Promoterwise;
use strict;

use Bio::SeqFeature::FeaturePair;
use Bio::SeqFeature::Generic;

use base qw(Bio::Root::Root Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Promoterwise->new();
 Function: Builds a new Bio::Tools::Promoterwise object
 Returns : L<Bio::Tools::Promoterwise>
 Args    : -fh/-file => $val, # for initing input, see Bio::Root::IO


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->_initialize_io(@args);
  my ($query1,$query2) = $self->_rearrange([qw(QUERY1_SEQ QUERY2_SEQ)],@args);
  $self->query1_seq($query1) if ($query1);
  $self->query2_seq($query2) if ($query2);

  return $self;
}

=head2 next_result

 Title   : next_result
 Usage   : my $r = $rpt_masker->next_result
 Function: Get the next result set from parser data
 Returns : an  L<Bio::SeqFeature::FeaturePair>
 Args    : none


=cut

sub next_result {
  my ($self) = @_;
  $self->_parse unless $self->_parsed;
  return $self->_next_result;
}

sub _parse{
   my ($self) = @_;
   my (%hash,@fp);
   while (defined($_ = $self->_readline()) ) {
       chomp;
       my @array = split;
       push @{$hash{$array[-1]}}, \@array;
   }
   foreach my $key(keys %hash){
    my $sf1 = Bio::SeqFeature::Generic->new(-primary=>"conserved_element",
                                            -source_tag=>"promoterwise");
    $sf1->attach_seq($self->query1_seq) if $self->query1_seq;
    my $sf2 = Bio::SeqFeature::Generic->new(-primary=>"conserved_element",
                                            -source_tag=>"promoterwise");
    $sf2->attach_seq($self->query2_seq) if $self->query2_seq;
    foreach my $info(@{$hash{$key}}){
	
      my ($score,$id1,$start_1,$end_1, $strand_1,$s1_len,
	  $id2,$start_2,$end_2,$strand_2,$s2_len, $group);
      if( @{$info} == 12 ) {
	  ($score,$id1,$start_1,$end_1, $strand_1,$s1_len,
	   $id2,$start_2,$end_2,$strand_2,$s2_len, $group) = @{$info};
      } elsif( @{$info} == 10 ) {
	  ($score,$id1,$start_1,$end_1, $strand_1,
	   $id2,$start_2,$end_2,$s2_len, $group) = @{$info};
      } else {
	  $self->throw("unknown promoterwise output, ", scalar @{$info},
		       " columns, expected 10 or 12\n");
      }
      if(!$sf1->strand && !$sf2->strand){
        $sf1->strand($strand_1);
        $sf2->strand($strand_2);
        $sf1->seq_id($id1);
        $sf2->seq_id($id2);
        $sf1->score($score);
        $sf2->score($score);
      }

      my $sub1 = Bio::SeqFeature::Generic->new(-start=>$start_1,
                                              -seq_id=>$id1,
                                              -end  =>$end_1,
                                              -strand=>$strand_1,
                                              -primary=>"conserved_element",
                                              -source_tag=>"promoterwise",
                                              -score=>$score);
      $sub1->attach_seq($self->query1_seq) if $self->query1_seq;

      my $sub2 = Bio::SeqFeature::Generic->new(-start=>$start_2,
                                              -seq_id=>$id2,
                                              -end  =>$end_2,
                                              -strand=>$strand_2,
                                              -primary=>"conserved_element",
                                              -source_tag=>"promoterwise",
                                              -score=>$score);
      $sub2->attach_seq($self->query2_seq) if $self->query2_seq;
      $sf1->add_SeqFeature($sub1,'EXPAND');
      $sf2->add_SeqFeature($sub2,'EXPAND');
    }

    my $fp = Bio::SeqFeature::FeaturePair->new(-feature1=>$sf1,
                                               -feature2=>$sf2);
    push @fp, $fp;
  }
    $self->_feature_pairs(\@fp);
    $self->_parsed(1);
    return;
}

sub _feature_pairs {
  my ($self,$fp) = @_;
  if($fp){
    $self->{'_feature_pairs'} = $fp;
  }
  return  $self->{'_feature_pairs'};
}

sub _next_result {
  my ($self) = @_;
  return unless (exists($self->{'_feature_pairs'}) && @{$self->{'_feature_pairs'}});
  return shift(@{$self->{'_feature_pairs'}});
}
sub _parsed {
  my ($self,$flag) = @_;
  if($flag){
    $self->{'_flag'} = 1;
  }
  return $self->{'_flag'};
}

sub query1_seq {
  my ($self,$val) = @_;
  if($val){
    $self->{'query1_seq'} = $val;
  }
  return $self->{'query1_seq'};
}
sub query2_seq {
  my ($self,$val) = @_;
  if($val){
    $self->{'query2_seq'} = $val;
  }
  return $self->{'query2_seq'};
}
1;
