#
# BioPerl module for Bio::Tools::Seg
#
# Copyright Balamurugan Kumarasamy
# Totally re-written, added docs and tests -- Torsten Seemann, Sep 2006
#
# Copyright 
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Seg - parse C<seg> output 

=head1 SYNOPSIS

  use Bio::Tools::Seg;
  my $parser = Bio::Tools::Seg->(-file => 'seg.fasta');
  while ( my $f = $parser->next_result ) {
    if ($f->score < 1.5) {
      print $f->location->to_FTstring, " is low complexity\n";
    }
  }

=head1 DESCRIPTION

C<seg> identifies low-complexity regions on a protein sequence.
It is usually part of the C<WU-BLAST> and C<InterProScan> packages.

The L<Bio::Tools::Seg> module will only parse the "fasta" output 
modes of C<seg>, i.e. C<seg -l> (low complexity regions only), 
C<seg -h> (high complexity regions only), or C<seg -a> (both low 
and high). 

It creates a L<Bio::SeqFeature::Generic> for each FASTA-like entry 
found in the input file. It is up to the user to appropriately filter 
these using the feature's score.

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

=head1 AUTHOR - Torsten Seemann

Email - torsten.seemann AT infotech.monash.edu.au

=head1 CONTRIBUTOR - Bala

Email - savikalpa@fugu-sg.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tools::Seg;
use strict;

use Bio::SeqFeature::Generic;
use base qw(Bio::Root::Root Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Seg->new();
 Function: Builds a new Bio::Tools::Seg object
 Returns : Bio::Tools::Seg
 Args    : -fh/-file => $val, # for initing input, see Bio::Root::IO

=cut


sub new {
	my($class,@args) = @_;
	my $self = $class->SUPER::new(@args);
	$self->_initialize_io(@args);
	return $self;
}

=head2 next_result

 Title   : next_result
 Usage   : my $feat = $seg->next_result
 Function: Get the next result set from parser data
 Returns : Bio::SeqFeature::Generic
 Args    : none

=cut

sub next_result {
	my ($self) = @_;

	# For example in this line 
	# test_prot(214-226) complexity=2.26 (12/2.20/2.50)
	# $1 is test_prot  
	# $2 is 214 
	# $3 is 226 
	# $4 is 2.26

	while (my $line = $self->_readline) {
		if ($line =~ /^\>\s*?(\S+)?\s*?\((\d+)\-(\d+)\)\s*complexity=(\S+)/) {
			return Bio::SeqFeature::Generic->new(
				-seq_id     => $1,
				-start      => $2,
				-end        => $3,
				-score      => $4,
				-source_tag => 'Seg',
				-primary    => 'low_complexity'
			);
		}
	}
}


1;
