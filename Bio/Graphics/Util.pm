package Bio::Graphics::Util;

# $Id$
# Non object-oriented utilities used here-and-there in Bio::Graphics modules

=head1 NAME

Bio::Graphics::Util - non-object-oriented utilities used in Bio::Graphics modules

=cut

use strict;
require Exporter;
use base qw(Exporter);
use vars '@EXPORT','@EXPORT_OK';
@EXPORT = 'frame_and_offset';
use Bio::Root::Version;

=over 4

=item ($frame,$offset) = frame_and_offset($pos,$strand,$phase)

Calculate the reading frame for a given genomic position, strand and
phase.  The offset is the offset from $pos to the first nucleotide
of the reading frame.

In a scalar context, returns the frame only.

=back

=cut

sub frame_and_offset {
  my ($pos,$strand,$phase) = @_;
  $strand ||= +1;
  $phase  ||= 0;
  my $codon_start =  $strand >= 0
                   ? $pos + $phase
	           : $pos - $phase;  # probably wrong
  my $frame  = ($codon_start-1) % 3;
#  my $frame = $strand >= 0
#    ? ($pos - $phase - 1) % 3
#    : (1 - $pos - $phase) % 3;
  my $offset = $strand >= 0 ? $phase : -$phase;
  return wantarray ? ($frame,$offset) : $frame;
}


1;
