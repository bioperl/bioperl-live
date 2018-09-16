#
# BioPerl module for Bio::Tools::Spidey::Exon
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ryan Golhar <golharam@umdnj.edu>
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Spidey::Exon - A single exon determined by an alignment

=head1 SYNOPSIS

  # See Bio::Tools::Spidey::Results for a description of the context.

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
EST sequence (using Spidey in this case). Consequently, the notion of query and
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
the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Ryan Golhar

Email golharam@umdnj.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Spidey::Exon;
use strict;


use base qw(Bio::SeqFeature::SimilarityPair);

sub new {
    my ($class,@args) = @_;
    my %param = @args;
    my $self = $class->SUPER::new(@args);

    my ($prim, $prim_tag, $source, $source_tag) = 
	$self->_rearrange([qw(PRIMARY
			      PRIMARY_TAG 
			      SOURCE
			      SOURCE_TAG)], 
			  @args);

    $self->primary_tag('exon') unless $prim || $prim_tag;
    $self->source_tag('Spidey') unless $source || $source_tag;
    $self->strand(0) unless defined($self->strand());
    $self->query();
    return $self; 
}

=head2 percentage_id

 Title   : percentage_id
 Usage   : $obj->percentage_id
 Function: This is the percent id as reported by Spidey
 Returns : value of percentage_id
 Args    : 


=cut

sub percentage_id {
	my ($self, @args) = @_;
	my $val;
    
	if(@args) {
	    $val = shift(@args);
	    $self->{'percentage_id'} = $val;
	} else {
	    $val = $self->{'percentage_id'};
	}
	return $val;
}

=head2 est_hit

 Title   : est_hit
 Usage   : $est_feature = $obj->est_hit();
 Function: Returns the EST hit pointing to (i.e., aligned to by Spidey) this
           exon (i.e., genomic region). At present, merely a synonym for
           $obj->feature2().
 Returns : An Bio::SeqFeatureI implementing object.
 Args    : 


=cut

sub est_hit {
    my $self = shift;
    return $self->feature2(@_);
}

=head2 mismatches

 Title   : mismatches
 Usage   : $obj->mismatches;
 Function: Returns the mismatches of the cDNA to (i.e., aligned to by Spidey) this
           exon (i.e., genomic region). 
 Returns : value of mismatches.
 Args    : 


=cut

sub mismatches {
	my ($self, @args) = @_;
	my $val;
    
	if(@args) {
	    $val = shift(@args);
	    $self->{'mismatches'} = $val;
	} else {
	    $val = $self->{'mismatches'};
	}
	return $val;
}

=head2 gaps

 Title   : gaps
 Usage   : $obj->gaps;
 Function: Returns the gaps of the cDNA to (i.e., aligned to by Spidey) this
           exon (i.e., genomic region). 
 Returns : value of gaps.
 Args    : 


=cut

sub gaps {
	my ($self, @args) = @_;
	my $val;
    
	if(@args) {
	    $val = shift(@args);
	    $self->{'gaps'} = $val;
	} else {
	    $val = $self->{'gaps'};
	}
	return $val;
}

=head2 donor

 Title   : donor
 Usage   : $obj->donor;
 Function: Returns 0 if a splice donor site does not exist, or 
           1 if a splice donor site exists
 Returns : value of existence of donor splice site (0 or 1)
 Args    :


=cut

sub donor {
	my ($self, @args) = @_;
	my $val;

	if (@args) {
		$val = shift @args;
		$self->{'donor'} = $val;
	} else {
		$val = $self->{'donor'};
	}
	return $val;
}

=head2 acceptor

 Title   : acceptor
 Usage   : $obj->acceptor;
 Function: Returns 0 if a splice acceptor site does not exist, or 
           1 if a splice acceptor site exists
 Returns : value of existence of acceptor splice site (0 or 1)
 Args    :


=cut

sub acceptor {
	my ($self, @args) = @_;
	my $val;

	if (@args) {
		$val = shift @args;
		$self->{'acceptor'} = $val;
	} else {
		$val = $self->{'acceptor'};
	}
	return $val;
}


1;
