=head1 NAME

Bio::DB::GFF::Aggregator::clone -- Clone aggregator

=head1 SYNOPSIS

  use Bio::DB::GFF;

  # Open the sequence database
  my $db      = Bio::DB::GFF->new( -adaptor => 'dbi:mysql',
                                   -dsn     => 'dbi:mysql:elegans42',
				   -aggregator => ['transcript','clone'],
				 );

 ----------------------------------------------------------------------------
 Aggregator method: clone
 Main method:       -none-
 Sub methods:       Clone_left_end Clone_right_end region:Genomic_canonical
 ----------------------------------------------------------------------------

=head1 DESCRIPTION

Bio::DB::GFF::Aggregator::clone is one of the default aggregators, and
was written to be compatible with the C elegans GFF files.  It
aggregates raw "Clone_left_end", "Clone_right_end", and
"region:Genomic_canonical" features into composite features of type
"clone".

=cut

package Bio::DB::GFF::Aggregator::clone;

use strict;


use base qw(Bio::DB::GFF::Aggregator);

=head2 aggregate

 Title   : aggregate
 Usage   : $features = $a->aggregate($features,$factory)
 Function: aggregate a feature list into composite features
 Returns : an array reference containing modified features
 Args    : see L<Bio::DB::GFF::Aggregator>
 Status  : Public

The WormBase GFF model is unusual in that clones aren't identified as
a single feature with start and stop positions, but as two features, a
"left end" and a "right end".  One or both of these features may be
absent.  In order to accommodate this, the aggregator will return undef
for the start and/or stop if one or both of the ends are missing.

=cut

#'

# we look for features of type Sequence and add them to a pseudotype transcript
sub aggregate {
  my $self = shift;
  my $features = shift;
  my $factory  = shift;

  my $matchsub    = $self->match_sub($factory) or return;
  my $passthru    = $self->passthru_sub($factory);
  my $method      = $self->get_method;

  my (%clones,%types,@result);
  for my $feature (@$features) {

    if ($feature->group && $matchsub->($feature)) {

      if ($feature->method =~ /^region|Sequence$/ && $feature->source eq 'Genomic_canonical') {
	$clones{$feature->group}{canonical} = $feature;
      } elsif ($feature->method eq 'Clone_left_end') {
	$clones{$feature->group}{left} = $feature;
      } elsif ($feature->method eq 'Clone_right_end') {
	$clones{$feature->group}{right} = $feature;
      }
      push @result,$feature if $passthru && $passthru->($feature);
    } else {
      push @result,$feature;
    }
  }

  for my $clone (keys %clones) {
    my $canonical = $clones{$clone}{canonical} or next;

    # the genomic_canonical doesn't tell us where the clone starts and stops
    # so don't assume it
    my $duplicate = $canonical->clone;   # make a duplicate of the feature
    # munge the method and source fields
    my $source    = $duplicate->source;
    my $type = $types{$method,$source} ||= Bio::DB::GFF::Typename->new($method,$source);
    $duplicate->type($type);

    my ($start,$stop) = $duplicate->strand > 0 ? ('start','stop') : ('stop','start');
    @{$duplicate}{$start,$stop} =(undef,undef);

    $duplicate->{$start} = $clones{$clone}{left}{$start}  if exists $clones{$clone}{left};
    $duplicate->{$stop}  = $clones{$clone}{right}{$stop}  if exists $clones{$clone}{right};
    $duplicate->method($self->method);
    push @result,$duplicate;
  }

  @$features = @result;
}

=head2 method

 Title   : method
 Usage   : $aggregator->method
 Function: return the method for the composite object
 Returns : the string "clone"
 Args    : none
 Status  : Public

=cut

sub method { 'clone' }

=head2 part_names

 Title   : part_names
 Usage   : $aggregator->part_names
 Function: return the methods for the sub-parts
 Returns : the list ("Clone_left_end", "Clone_right_end", "region:Genomic_canonical")
 Args    : none
 Status  : Public

=cut

sub part_names {
  my $self = shift;
  return qw(Clone_left_end Clone_right_end region:Genomic_canonical Sequence:Genomic_canonical);
}

1;

__END__

=head1 BUGS

None reported.


=head1 SEE ALSO

L<Bio::DB::GFF>, L<Bio::DB::GFF::Aggregator>

=head1 AUTHOR

Lincoln Stein E<lt>lstein@cshl.orgE<gt>.

Copyright (c) 2001 Cold Spring Harbor Laboratory.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

