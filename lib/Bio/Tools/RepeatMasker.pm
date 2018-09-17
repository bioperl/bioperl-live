#
# BioPerl module for Bio::Tools::RepeatMasker
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

Bio::Tools::RepeatMasker - a parser for RepeatMasker output

=head1 SYNOPSIS

    use Bio::Tools::RepeatMasker;
    my $parser = Bio::Tools::RepeatMasker->new(-file => 'seq.fa.out');
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

=head1 AUTHOR - Shawn Hoon

Email shawnh@fugu-sg.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::RepeatMasker;
use strict;

use Bio::SeqFeature::FeaturePair;

use base qw(Bio::Root::Root Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::RepeatMasker->new();
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
           Feature1 is the Query coordinates and Feature2 is the Hit
 Args    : none

=cut

sub next_result {
    my ($self) = @_;
    local $_;
    while (defined($_=$self->_readline()) ) {
	if (/no repetitive sequences detected/) {
	    $self->warn( "RepeatMasker didn't find any repetitive sequences\n");
	    return ;
	}
	#ignore introductory lines
	if (/\d+/) {
	    my @element = split;
	    # ignore features with negatives
	    next if ($element[11-13] =~ /-/);
	    my (%feat1, %feat2);
	    my @line = split;
	    my ($score, $query_name, $query_start, $query_end, $strand,
		$repeat_name, $repeat_class ) = @line[0, 4, 5, 6, 8, 9, 10];

	    my ($hit_start,$hit_end);

	    if ($strand eq '+') {
		($hit_start, $hit_end) = @line[11, 12];
		$strand = 1;
	    } elsif ($strand eq 'C') {
		($hit_end, $hit_start) = @line[12, 13];
		$strand = -1;
	    }
	    my $rf = Bio::SeqFeature::Generic->new
		(-seq_id      => $query_name,
		 -score       => $score,
		 -start       => $query_start,
		 -end         => $query_end,
		 -strand      => $strand,
		 -source_tag  => 'RepeatMasker',
		 -primary_tag => $repeat_class,
		 -tag => { 'Target'=> [$repeat_name, $hit_start, $hit_end]},
		);

	    my $rf2 = Bio::SeqFeature::Generic->new
		(-seq_id         => $repeat_name,
		 -score          => $score,
		 -start          => $hit_start,
		 -end            => $hit_end,
		 -strand         => $strand,
		 -source_tag     => "RepeatMasker",
		 -primary_tag    => $repeat_class,
		 -tag => { 'Target'=> [$query_name,$query_start,$query_end] },
		);

	    my $fp = Bio::SeqFeature::FeaturePair->new(-feature1 => $rf,
						       -feature2 => $rf2);
	    return $fp;
	}
    }
}

1;
