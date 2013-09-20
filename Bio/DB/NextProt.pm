package Bio::DB::NextProt;

use strict;
use warnings;


sub new {
	my ($class, @args) = @_;
	#my $self = $class->SUPER::new(@args);
	my $self = {};
	$self->{_fileInput}	= undef,
	bless($self, $class);
	return $self;
}


1;
