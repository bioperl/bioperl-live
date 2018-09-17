# BioPerl module for Bio::Tools::SiRNA::Ruleset::saigo
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

Bio::Tools::SiRNA::Ruleset::saigo - Perl object implementing the Saigo
group's rules for designing small inhibitory RNAs

=head1 SYNOPSIS

Do not use this module directly.  Instead, use Bio::Tools::SiRNA and
specify the saigo ruleset:

  use Bio::Tools::SiRNA;

  my $sirna_designer = Bio::Tools::SiRNA->new( -target => $bio_seq,
                                               -rules  => 'saigo'
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
published by Ui-Tei et al (2004).  The rules are:

=over 5

=item 1.

The first base in the sense strand of the duplex must be a G or C

=item 2.

The first base in the antisense strand of the duplex must be an A or U

=item 3.

The first 7 nucleotides in the antisense strand of the duplex must be
A or U

=item 4.

There cannot be more than 9 consecutive G or C nucleotides

=item 5.

The first 12 nucleotides in the sense strand of the duplex should have
33-66% GC

=back

The module inherits from Bio::Tools::SiRNA.  See the documentation for
that module for information on how to specify the target and recover
the SiRNA duplex information.

=head2 EXPORT

None.

=head1 SEE ALSO

L<Bio::Tools::SiRNA>, 
L<Bio::SeqFeature::SiRNA::Pair>,
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

package Bio::Tools::SiRNA::Ruleset::saigo;

use strict;
use warnings;

use base qw(Bio::Tools::SiRNA);

=head2 new

  Title	: new
  Usage  : Do not call directly - use Bio::Tools::SiRNA->new instead.
  Returns : Bio::Tools::SiRNA::Ruleset::saigo object
  Args	: none

=cut

sub new {
    my ($proto, %args) = @_;
    my $class = ref($proto) || $proto;
    
    $args{'RULES'} = 'saigo';
    
    return $class->SUPER::new(%args);
 }

sub _get_oligos {
    my ($self) = @_;

    my ($targseq, $targstart) = $self->_get_targetregion;

    foreach my $i (0 .. (length($targseq) - 23)) {
	my $testseq = substr($targseq, $i, 23);		
	$self->add_oligos($testseq, $targstart + $i + 1) if ($self->_oligo_ok($testseq));
    }
}


sub _get_sense {
    my ($self, $target) = @_;
    #trim off 1st 2 nt to get overhang
    $target =~ s/^..//;
    #convert T's to U's (transcribe)
    $target =~ s/T/U/gi;

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
    #trim off 1st 2 nt to get overhang
    $anti =~ s/^..//;
    #convert T's to U's
    $anti =~ s/T/U/gi;

    return $anti;
}

sub _oligo_ok {
    my ($self, $testseq) = @_;

    $self->debug("Testing $testseq...\n");

    my @testseq = split(//, $testseq);
    # is 5p end of sense strand a G/C?
    unless ($testseq[2] =~ /[GC]/i) {
	$self->debug("No G/C at sense 5' end\n");
	return 0;
    }
    # is 5p end of antisense strand an A/T?
    unless ($testseq[20] =~ /[AT]/i) {
	$self->debug("No A/T at antisense 5' end\n");
	return 0;
    }

    # are 4 of the last 7 bases in the duplex A/T?
    my $atcount_3p = grep { /[AT]/i } @testseq[14 .. 20];
    unless ($atcount_3p >= 4) {
	$self->debug("Found $atcount_3p A/T in last 7 bases of duplex\n");
	return 0;
    }
    # what is gc fraction in rest of duplex? Target: 33 to 66 pct gc (4-8 of 12)
    my $gccount_5p = grep { /[GC]/i } @testseq[2 .. 13];
    if ($gccount_5p < 4) {
	$self->debug("Found only $gccount_5p GCs in 5p end of duplex\n");
	return 0;
    }
    if ($gccount_5p > 8) {
	$self->debug("Found only $gccount_5p GCs in 5p end of duplex\n");
	return 0;
    }
    
    # no more than 9 consecutive GC
    if ($testseq =~ /[GC]{9,}?/i) {
	$self->debug("Found more than 9 consecutive GCs\n");
	return 0;
    }

    $self->debug("Oligo passed \n");
    return 1;
}

1;
