# $Id$
#
# BioPerl module for Bio::DB::BioFetch
#
# Cared for by Lincoln Stein <lstein@cshl.org>
#
# Copyright Lincoln Stein
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#

package Bio::DB::BioFetch;
use strict;
use Bio::DB::WebDBSeqI;

use vars (@ISA $VERSION);
@ISA = 'Bio::DB::BioFetch';
$VERSION = '1.0';

# warning: names used here must map into Bio::SeqIO::* space
use constant DEFAULT_FORMAT   => 'embl';
use constant DEFAULT_LOCATION => 'http://www.ebi.ac.uk/cgi-bin/dbfetch';
my %SUPPORTED_FORMATS = map {$_=>1} qw(fasta embl);

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($format, $hostlocation) = $self->_rearrange([qw(FORMAT HOSTLOCATION)],@args);
    $format ||= $self->default_format;

    $self->throw('Requested format is invalid; only '.
		 join(' ',keys %SUPPORTED_FORMATS).
		 ' formats are supported') unless $SUPPORTED_FORMATS{$format};

    $self->servertype('ebi');
    $hostlocation ||= DEFAULT_LOCATION;
    $self->hostlocation($hostlocation);

    if (  $hostlocation ) {
	$self->hostlocation(lc $hostlocation);
    }

    $self->request_format($format); # let's always override the format, as it must be swiss or fasta
    return $self;
}

sub default_format { return DEFAULT_FORMAT }


1;
