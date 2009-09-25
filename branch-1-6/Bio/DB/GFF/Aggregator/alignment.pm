=head1 NAME

Bio::DB::GFF::Aggregator::alignment -- Alignment aggregator

=head1 SYNOPSIS

  use Bio::DB::GFF;

  # Open the sequence database
  my $db      = Bio::DB::GFF->new( -adaptor => 'dbi:mysql',
                                   -dsn     => 'dbi:mysql:elegans42',
				   -aggregator => ['alignment'],
				 );

 -----------------------------
 Aggregator method: alignment
 Main method:       (none)
 Sub methods:       nucleotide_match,EST_match,cDNA_match,expressed_sequence_match,
                    translated_nucleotide_match,protein_match,HSP
 -----------------------------

=head1 DESCRIPTION

Bio::DB::GFF::Aggregator::alignment is one of the default aggregators,
and was written to be compatible with the C elegans GFF files.  It
aggregates raw "similarity" features into composite features of type
"alignment".  A better name for this class might be
"gapped_alignment."

This aggregator does not insist that there be a single top-level
feature that spans one end of the alignment to the other.  As a
result, it can produce truncated alignments if the entire alignment is
not contained within the segment of interest.

=cut

package Bio::DB::GFF::Aggregator::alignment;

use strict;


use base qw(Bio::DB::GFF::Aggregator);

=head2 aggregate

 Title   : aggregate
 Usage   : $features = $a->aggregate($features,$factory)
 Function: aggregate a feature list into composite features
 Returns : an array reference containing modified features
 Args    : see L<Bio::DB::GFF::Aggregator>
 Status  : Public

Because of the large number of similarity features, the aggregate()
method is overridden in order to perform some optimizations.

=cut

# we look for features of type Sequence and add them to a pseudotype transcript
sub aggregate {
  my $self = shift;
  my $features = shift;
  my $factory  = shift;

  my $matchsub = $self->match_sub($factory) or return;
  my $passthru = $self->passthru_sub($factory);
  my $method   = $self->get_method;

  my (%alignments,%targets,@result);

  warn "running alignment aggregator" if $factory->debug;
  for my $feature (@$features) {

    if ($matchsub->($feature)) {

      my $group  = $feature->{group};
      my $source = $feature->source;
      unless (exists $alignments{$group,$source}) {
	my $type = Bio::DB::GFF::Typename->new($method,$source);
	
	my $f = $feature->clone;
	# this is a violation of OO encapsulation, but need to do it this way
	# to achieve desired performance
	@{$f}{qw(type score phase)} = ($type,undef,undef);

	$alignments{$group,$source} = $f or next;
      }

      my $main = $alignments{$group,$source};
      $main->add_subfeature($feature);
      push @result,$feature if $passthru && $passthru->($feature);
    } else {
      push @result,$feature;
    }
  }

  warn "running aligner adjuster" if $factory->debug;
  for my $alignment (values %alignments) {
    $alignment->adjust_bounds;
    $alignment->compound(1);
    push @result,$alignment;
  }
  warn "aligner done" if $factory->debug;
  @$features = @result;
}

=head2 method

 Title   : method
 Usage   : $aggregator->method
 Function: return the method for the composite object
 Returns : the string "alignment"
 Args    : none
 Status  : Public

=cut

sub method { 'alignment' }

=head2 part_names

 Title   : part_names
 Usage   : $aggregator->part_names
 Function: return the methods for the sub-parts
 Returns : the full list of aggregated methods
 Args    : none
 Status  : Public

=cut

sub part_names {
  my $self = shift;
 return qw(nucleotide_match EST_match cDNA_match
	   expressed_sequence_match
	   translated_nucleotide_match
	   protein_match HSP);
}

1;
