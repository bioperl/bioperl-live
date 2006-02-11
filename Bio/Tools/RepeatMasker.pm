# $Id$
#
# BioPerl module for Bio::Tools::RepeatMasker
#
# Cared for by Shawn Hoon <shawnh@fugu-sg.org>
#
# Copyright Shawn Hoon
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::RepeatMasker

=head1 SYNOPSIS

    use Bio::Tools::RepeatMasker;
    my $parser = new Bio::Tools::RepeatMasker(-file => 'seq.fa.out');
    while( my $result = $parser->next_result ) {
      # get some value
    }

=head1 DESCRIPTION

A parser for RepeatMasker output

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
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Shawn Hoon

Email shawnh@fugu-sg.org

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::RepeatMasker;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::SeqFeature::FeaturePair;
use Bio::Root::IO;

@ISA = qw(Bio::Root::Root Bio::Root::IO );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::RepeatMasker();
 Function: Builds a new Bio::Tools::RepeatMasker object
 Returns : Bio::Tools::RepeatMasker
 Args    : -fh/-file => $val, for initing input, see Bio::Root::IO

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->_initialize_io(@args);

  return $self;
}

=head2 next_result

 Title   : next_result
 Usage   : my $r = $rpt_masker->next_result
 Function: Get the next result set from parser data
 Returns : Bio::SeqFeature::FeaturePair
 Args    : none

=cut

sub next_result{
   my ($self) = @_;
   while ($_=$self->_readline()) {
        if (/no repetitive sequences detected/) {
           $self->warn( "RepeatMasker didn't find any repetitive sequences\n");
           return ;
        }
        if (/\d+/) { #ignore introductory lines
          my @element = split;
          # ignore features with negatives
          next if ($element[11-13] =~ /-/);
          my (%feat1, %feat2);
          my ($score, $query_name, $query_start, $query_end, $strand,
          $repeat_name, $repeat_class ) = (split)[0, 4, 5, 6, 8, 9, 10];

          my ($hit_start,$hit_end);
          if ($strand eq '+') {
            ($hit_start, $hit_end) = (split)[11, 12];
            $strand = 1;
          }
          elsif ($strand eq 'C') {
            ($hit_start, $hit_end) = (split)[12, 13];
            $strand = -1;
          }
          my $rf = Bio::SeqFeature::Generic->new
	      (-seq_id      => $query_name,
	       -score       => $score,
	       -start       => $query_start,
	       -end         => $query_end,
	       -strand      => $strand,
	       -source_tag  => 'RepeatMasker',
	       -primary_tag => $repeat_class);

          my $rf2 = Bio::SeqFeature::Generic->new
	      (-seq_id         => $repeat_name,
	       -score          => $score,
	       -start          => $hit_start,
	       -end            => $hit_end,
	       -strand         => $strand,
	       -source_tag     => "RepeatMasker",
	       -primary_tag    => $repeat_class);
	  
          my $fp = Bio::SeqFeature::FeaturePair->new(-feature1=>$rf2,
                                                     -feature2=>$rf);
          return $fp;
        }
    }
}

1;
