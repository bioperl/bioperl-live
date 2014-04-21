#
#
# BioPerl module for Bio::Tools::SiRNA::Ruleset::tuschl
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Donald Jackson, donald.jackson@bms.com
#
# Copyright Bristol-Myers Squibb
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::SiRNA::Ruleset::tuschl - Perl object implementing the
tuschl group's rules for designing small inhibitory RNAs

=head1 SYNOPSIS

Do not use this module directly.  Instead, use Bio::Tools::SiRNA and
specify the tuschl ruleset:

  use Bio::Tools::SiRNA;

  my $sirna_designer = Bio::Tools::SiRNA->new( -target => $bio_seq,
                                               -rules  => 'tuschl'
    );
  my @pairs = $sirna_designer->design;

  foreach $pair (@pairs) {
      my $sense_oligo_sequence = $pair->sense->seq;
      my $antisense_oligo_sequence = $pair->antisense->seq;

      # print out results
      print join ("\t", $pair->start, $pair->end, $pair->rank,
                  $sense_oligo_sequence, $antisense_oligo_sequence), "\n";
  }

=head1 DESCRIPTION

This package implements the rules for designing siRNA reagents
developed by Tuschl and colleagues (see
http://www.rockefeller.edu/labheads/tuschl/sirna.html). It looks for
oligos that match the following patterns in the target sequence:

  1. AA(N19)TT (rank 1)
  2. AA(N21) (rank 2)
  3. NA(N21) (rank 3)

The package also supports selection of siRNA seqences that can be
transcribed by pol3:

    A[A,G]N17[C,T]

=head1 SEE ALSO

L<Bio::Tools::SiRNA>, L<Bio::SeqFeature::SiRNA::Pair>,
L<Bio::SeqFeature::SiRNA::Oligo>.

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

=head1 AUTHOR

Donald Jackson (donald.jackson@bms.com)

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _


=cut

package Bio::Tools::SiRNA::Ruleset::tuschl;

use strict;
use warnings;

use base qw(Bio::Tools::SiRNA);

our %PATTERNS = ( 1 	=> '(AA.{19}TT)',
		  2 	=> '(AA.{19}[ACG][ACG])',
		  3 	=> '([CGT]A.{21})',
		  Pol3	=> '(.A[AG].{17}[CT]..)'
		  );

our $DEFAULT_CUTOFF = 2;

=head2 new

  Title	: new
  Usage	: Do not call directly - use Bio::Tools::SiRNA->new instead.
  Returns : Bio::Tools::SiRNA::Ruleset::saigo object
  Args	: none

=cut

sub new {
    my ($proto, %args) = @_;
    my $class = ref($proto) || $proto;
    
    $args{'RULES'} = 'tuschl';

    return $class->SUPER::new(%args);
}

sub _regex {
    my ($self, $rank) = @_;
    return $PATTERNS{$rank};
}

sub cutoff {
    my ($self, $cutoff) = @_;
    if ($cutoff) {
	$self->{'cutoff'} = $cutoff;
    }
    elsif (!$self->{'cutoff'}) {
	$self->{'cutoff'} = $DEFAULT_CUTOFF;
    }
    return $self->{'cutoff'};
}


sub _get_oligos {
    #use regular expressions to pull out oligos
    my ($self) = @_;

    my @ranks;
    if ($self->cutoff eq 'pol3') {
	@ranks = ('pol3');
    }
    else {
	@ranks = (1 .. $self->cutoff);
    }
    
    foreach my $rank (@ranks) {
	my $regex = $self->_regex($rank);
	#my @exclude;


# 	my ($targregion) = grep { $_->primary_tag eq 'Target' } $self->target->top_SeqFeatures;
# 	my $seq = $targregion->seq->seq;
# 	# but this way I loose start info
# 	my $targstart = $targregion->start;
	my ($seq, $targstart) = $self->_get_targetregion();

	while ( $seq =~ /(.*?)$regex/gi ) {
	    my $target = $2;

	    # check for too many Gs (or Cs on the other strand)
	    next if ( $target =~ /G{ $self->gstring,}/io );
	    next if ( $target =~ /C{ $self->gstring,}/io );
# 	skip Ns (for filtering)
	    next if ( $target =~ /N/i);

	    my $start = length($1) + $targstart;
	    my $stop = $start + length($target) -1;

	    my @gc = ( $target =~ /G|C/gi);
	    my $fxGC = sprintf("%2.2f", (scalar(@gc) / length($target)));
	    next if ($fxGC < $self->min_gc);
	    next if ($fxGC > $self->max_gc);
	    
	    $self->add_oligos($target, $start, $rank);
	}
    }
}

	 
sub _get_sense {
    my ($self, $target) = @_;
    # trim off 1st 2 nt to get overhang
    $target =~ s/^..//;
    # convert T's to U's (transcribe)
    $target =~ s/T/U/gi;
    # force last 2 nt to be T's
    $target =~ s/..$/TT/;

    return $target;
}

sub _get_anti {
    my ($self, $target) = @_;
    my @target = split(//, $target);
    my ($nt,@antitarget);

    while ($nt = pop @target) {
	push(@antitarget, $self->_comp($nt));
    }
    my $anti = join('', @antitarget);
    # trim off 1st 2 nt to get overhang
    $anti =~ s/^..//;
    # convert T's to U's
    $anti =~ s/T/U/gi;
    # convert last 2 NT's to T
    $anti =~ s/..$/TT/;

    return $anti;
}


1;
