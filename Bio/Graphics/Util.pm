package Bio::Graphics::Util;

# $Id$
# Non object-oriented utilities used here-and-there in Bio::Graphics modules

use strict;
require Exporter;
use vars '@ISA','@EXPORT','@EXPORT_OK';
@ISA = 'Exporter';
@EXPORT = 'frame_and_offset';

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
  my $frame = $strand >= 0 
    ? ($pos - $phase - 1) % 3
    : (1 - $pos - $phase) % 3;
  my $offset = -$phase % 3;
  $offset   *= -1 if $strand < 0;
  return wantarray ? ($frame,$offset) : $frame;
}


1;
