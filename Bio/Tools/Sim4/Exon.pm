
#
# BioPerl module for Bio::Tools::Sim4::Exon
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
# and Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Ewan Birney, Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Sim4::Exon - A single exon determined by an alignment

=head1 SYNOPSIS

  # See Bio::Tools::Sim4::Results for a description of the context.

  # an instance of this class is-a Bio::SeqFeature::SimilarityPair

  # coordinates of the exon (recommended way):
  print "exon from ", $exon->start(),
  	" to ", $exon->end(), "\n";

  # the same (feature1() inherited from Bio::SeqFeature::FeaturePair)
  print "exon from ", $exon->feature1()->start(),
  	" to ", $exon->feature1()->end(), "\n";
  # also the same (query() inherited from Bio::SeqFeature::SimilarityPair):
  print "exon from ", $exon->query()->start(),
  	" to ", $exon->query()->end(), "\n";

  # coordinates on the matching EST (recommended way):
  print "matches on EST from ", $exon->est_hit()->start(),
  	" to ", $exon->est_hit()->end(), "\n";

  # the same (feature2() inherited from Bio::SeqFeature::FeaturePair)
  print "matches on EST from ", $exon->feature2()->start(),
  	" to ", $exon->feature2()->end(), "\n";
  # also the same (subject() inherited from Bio::SeqFeature::SimilarityPair):
  print "exon from ", $exon->subject()->start(),
  	" to ", $exon->subject()->end(), "\n";

=head1 DESCRIPTION

This class inherits from Bio::SeqFeature::SimilarityPair and represents an
exon on a genomic sequence determined by similarity, that is, by aligning an
EST sequence (using Sim4 in this case). Consequently, the notion of query and
subject is always from the perspective of the genomic sequence: query refers
to the genomic seq, subject to the aligned EST hit. Because of this,
$exon-E<gt>start(), $exon-E<gt>end() etc will always return what you expect. 

To get the coordinates on the matching EST, refer to the properties of the
feature returned by L<est_hit>().

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Ewan Birney, Hilmar Lapp

Email birney@sanger.ac.uk
Hilmar Lapp E<lt>hlapp@gmx.netE<gt> or E<lt>hilmar.lapp@pharma.novartis.comE<gt>.

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Sim4::Exon;
use vars qw(@ISA);
use strict;

use Bio::SeqFeature::FeaturePair;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::SimilarityPair;

@ISA = qw(Bio::SeqFeature::SimilarityPair);

sub new {
    my ($class,@args) = @_;
    my %param = @args;
    my $self = $class->SUPER::new(@args);

    my ($prim, $source) = $self->_rearrange([qw(PRIMARY SOURCE)], @args);

    $self->primary_tag('exon') unless $prim;
    $self->source_tag('Sim4') unless $source;
    $self->strand(0) unless defined($self->strand());
    $self->query();
    return $self; 
}

=head2 percentage_id

 Title   : percentage_id
 Usage   : $obj->percentage_id($newval)
 Function: This is a synonym for 100 * $obj->est_hit()->frac_identical().
 Returns : value of percentage_id
 Args    : newvalue (optional)


=cut

sub percentage_id {
    my ($self, @args) = @_;
    my $frac;
    my $val;
    my $delegated = 0;
    
    if(@args) {
	$frac = $args[0];
	$frac /= 100.0 if defined($frac);
    }
    if($self->query()->can('frac_identical')) {
	if(defined($frac)) {
	    $self->query()->frac_identical($frac);
	}
	$val = 100.0 * $self->query()->frac_identical();
	$delegated = 1;
    }
    if($self->est_hit()->can('frac_identical')) {
	if(defined($frac)) {
	    $self->est_hit()->frac_identical($frac);
	}
	# this intentiously overwrites previous $val
	$val = 100.0 * $self->est_hit()->frac_identical();
	$delegated = 1;
    }
    if(! $delegated) {
	if(@args) {
	    $val = shift(@args);
	    $self->{'percentage_id'} = $val;
	} else {
	    $val = $self->{'percentage_id'};
	}
    }
    return $val;
}

=head2 est_hit

 Title   : est_hit
 Usage   : $est_feature = $obj->est_hit();
 Function: Returns the EST hit pointing to (i.e., aligned to by Sim4) this
           exon (i.e., genomic region). At present, merely a synonym for
           $obj->feature2().
 Returns : An Bio::SeqFeatureI implementing object.
 Args    : 


=cut

sub est_hit {
    my $self = shift;
    return $self->feature2(@_);
}

1;
