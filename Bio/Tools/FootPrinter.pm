# BioPerl module for Bio::Tools::FootPrinter
#
# Cared for by Shawn Hoon <shawnh@fugu-sg.org>
#
# Copyright Shawn Hoon
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::FootPrinter - DESCRIPTION of Object

=head1 SYNOPSIS

    use Bio::Tools::FootPrinter;

    my $tool = Bio::Tools::FootPrinter->new(-file=>"footprinter.out");

    while (my $result = $tool->next_feature){
      foreach my $feat($result->sub_SeqFeature){
        print $result->seq_id."\t".$feat->start."\t".$feat->end."\t".$feat->seq->seq."\n";
      }
    }

=head1 DESCRIPTION

A parser for FootPrinter  output 

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Shawn Hoon 

Email shawnh@fugu-sg.org 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::FootPrinter;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::SeqFeature::Generic;
use Bio::PrimarySeq;
use Bio::Root::IO;

@ISA = qw(Bio::Root::Root Bio::Root::IO );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::FootPrinter();
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
  my ($seq,$third,$name);
  while ($_ = $self->_readline) {
    chomp;
    my @array = split("\n",$_);
    if($#array == 3){
        if($name){
            $name=~s/>//;
            my $feat = $self->_parse($name,$seq,$third);
            $self->_add_feature($feat);
        }
        $name = shift @array;
        $seq=$array[0];
        $third=$array[2];
        next;
    }
    $seq.=$array[0];
    $third.=$array[2];
  }
  $name=~s/>//;
  my $feat = $self->_parse($name,$seq,$third);
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
  my ($self,$name,$seq,$pattern) = @_;
    my @char = split('',$pattern);
    my $prev;
   my $word;
    my @words;
    foreach my $c(@char){
        if(!$word){
            $word .= $c;
            $prev = $c;
            next;
        }
        if ($c eq $prev){
          $word.=$c;
          $prev = $c;
        }
        else {
            #remove words with only \s
            $word=~s/\s+//g;
            if ($word ne ''){
              push @words, $word;
            }
            $word=$c;
           $prev = $c;
           
        }
    }
    $word=~s/\s+//g;
    if($word ne ''){
      push @words, $word;
    }
    my $last;
    my $feat = new Bio::SeqFeature::Generic(-seq_id=>$name);
    my $offset=0;
    foreach my $w(@words){
        if($w !~ /^$/){
          my $index = index($pattern,$w,$offset);
          $offset = $index + length($w);
          my $subfeat = new Bio::SeqFeature::Generic ( -seq_id=>$name,
                                                    -start => $index+1, 
                                                    -end =>$index+length($w),
                                                    -source=>"FootPrinter");
          $feat->add_sub_SeqFeature($subfeat,'EXPAND');
        }
    }
    my $priseq = Bio::PrimarySeq->new(-id=>$name,-seq=>$seq);
    $feat->attach_seq($priseq);
    return $feat;
    
}

1;
