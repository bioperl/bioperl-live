# BioPerl module for Bio::Tools::FootPrinter
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

Bio::Tools::FootPrinter - write sequence features in FootPrinter format

=head1 SYNOPSIS

    use Bio::Tools::FootPrinter;

    my $tool = Bio::Tools::FootPrinter->new(-file=>"footprinter.out");

    while (my $result = $tool->next_feature){
      foreach my $feat($result->sub_SeqFeature){
        print $result->seq_id."\t".$feat->start."\t".$feat->end."\t".$feat->seq->seq."\n";
      }
    }

=head1 DESCRIPTION

This module writes sequence features in FootPrinter format. 
See L<http://bio.cs.washington.edu/software.html> for more details.

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


package Bio::Tools::FootPrinter;
use strict;

use Bio::SeqFeature::Generic;
use Bio::PrimarySeq;

use base qw(Bio::Root::Root Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::FootPrinter->new();
 Function: Builds a new Bio::Tools::FootPrinter object 
 Returns : Bio::Tools::FootPrinter
 Args    : -fh/-file => $val, # for initing input, see Bio::Root::IO

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->_initialize_io(@args);

  return $self;
}

=head2 next_feature

 Title   : next_feature
 Usage   : my $r = $footprint->next_feature
 Function: Get the next feature from parser data
 Returns : L<Bio::SeqFeature::Generic> 
 Args    : none

=cut

sub next_feature{
   my ($self) = @_;
   $self->_parse_predictions() unless $self->_predictions_parsed();
   return shift @{$self->{'_feature'}};

}

=head2 _add_feature

 Title   : _add_feature
 Usage   : $footprint->_add_feature($feat)
 Function: Add feature to array
 Returns : none
 Args    : none

=cut

sub _add_feature {
    my ($self,$feat) = @_;
    if($feat){
        push @{$self->{'_feature'}},$feat;
    }
}

=head2 _parse_predictions

 Title   : _parse_predictions
 Usage   : my $r = $footprint->_parse_predictions
 Function: do the parsing 
 Returns : none 
 Args    : none

=cut

sub _parse_predictions {
  my ($self) = @_;
  $/="";
  my ($seq,$second,$third,$name);
  while ($_ = $self->_readline) {
    chomp;
    my @array = split("\n",$_);
    if ($#array == 5) {
      # get rid of header
      shift(@array); shift(@array);
    }
    if($#array == 3){
        if($name){
            $name=~s/>//;
            my $feat = $self->_parse($name,$seq,$second,$third);
            $self->_add_feature($feat);
        }
        $name    = shift @array;
        $seq     = $array[0];
        $second  = $array[1];
        $third   = $array[2];
        next;
    }
    $seq        .= $array[0];
    $third      .= $array[2];
  }
  
  $seq || return;
  
  $name=~s/>//;
  my $feat = $self->_parse($name,$seq,$second,$third);
  $self->_add_feature($feat);

  $self->_predictions_parsed(1);
}

=head2 _predictions_parsed

 Title   : _predictions_parsed
 Usage   : $footprint->_predictions_parsed(1)
 Function: Get/Set for whether predictions parsed
 Returns : 1/0
 Args    : none

=cut

sub _predictions_parsed {
    my ($self,$val) = @_;
    if($val){
        $self->{'_predictions_parsed'} = $val;
    }
    return $self->{'_predictions_parsed'};
}


=head2 _parse

 Title   : _parse
 Usage   : $footprint->_parse($name,$seq,$pattern)
 Function: do the actual parsing
 Returns : L<Bio::SeqFeature::Generic>
 Args    : none

=cut

sub _parse {
    my ($self,$name,$seq,$score,$pattern) = @_;
    my @char  = split('',$pattern);
    my @score = split('',$score);

    my ($prev,$word,@words,@word_scores,$word_score);

    my $i = 0;
    for my $c ( @char ) {
        if( ! $word) {
            $word .= $c;
            $prev = $c;
	    defined $score[$i] && 
		($score[$i] =~ /\d/) && ($word_score += $score[$i]);
        } elsif ($c eq $prev){
	    $word .=$c;
	    $prev  = $c;
	    defined $score[$i] && 
		($score[$i] =~ /\d/) && ($word_score += $score[$i]);
        } else {
            # remove words with only \s
            $word=~s/\s+//g;
            if ($word ne ''){
		push @words, $word;
		push @word_scores, ($word_score/length($word));
            }
            $word =$c;
	    $prev = $c;
	    $word_score = 0;
	    defined $score[$i] &&
		($score[$i] =~ /\d/) && ($word_score += $score[$i]);
        }
	$i++;
    }
    $word =~s/\s+//g;
    if( length($word) ){
	push @words, $word;
    }
    my $last;
    my $feat = Bio::SeqFeature::Generic->new(-seq_id=>$name);
    my $offset = $i = 0;
    my $count = 1;
    for my $w (@words){
        if(length($w) ) { 
	    my $index = index($pattern,$w,$offset);
	    $offset = $index + length($w);
	    my $subfeat = Bio::SeqFeature::Generic->new 
		( -seq_id  =>"$name-motif".$count++,
		  -start   => $index+1, 
		  -end     => $index+length($w),
		  -source  =>"FootPrinter",
		  -score   => $word_scores[$i]
		  );
	    # ugh - not sure the sub_SeqFeature situation will
	    # be around forever- things should probably be
	    # grouped by a 'group' tag instead ala GFF3 
	    # perhaps when Lincoln's API changes are 
	    # made to SeqFeatures this will get changed
	    $feat->add_sub_SeqFeature($subfeat,'EXPAND');
        }
	$i++;
    }
    my $priseq = Bio::PrimarySeq->new(-id=>$name,-seq=>$seq);
    $feat->attach_seq($priseq);
    return $feat;

}

1;
