package Bio::DB::NextProt;

use strict;
use warnings;
use JSON;
use REST::Client;

sub new {
	my ($class, @args) = @_;
	#my $self = $class->SUPER::new(@args);
	my $self = {};
	$self->{_client}	= REST::Client->new({host=> "http://www.nextprot.org", timeout => 10,});
	$self->{_query}		= undef;
	$self->{_filter}	= undef;
	$self->{_format}	= "json";
	bless($self, $class);
	return $self;
}

sub search_protein() {
	my $self  = shift;
    my %param = @_;

	my $path = "/rest/protein/list";
	my $string = '';

	if (defined $param{'-query'} && defined $param{'-filter'}) {
		
		$self->{_client}->GET($path."?query=".$param{'-query'}."&filter=".$param{'-filter'}."&format=".$self->{_format});

	} elsif (defined $param{'-query'}) {
		
		$self->{_client}->GET($path."?query=".$param{'-query'}."&format=".$self->{_format});

	} elsif (defined $param{'-filter'}) {
		
		$self->{_client}->GET($path."?filter=".$param{'-filter'}."&format=".$self->{_format});
	}

	&reset_params();

	return $self->{_client}->responseContent();

}

sub search_cv() {
	my $self  = shift;
	my %param = @_;
	
	my $path   = "/rest/cv/list";
	my $string = '';

    if (defined $param{'-query'} && defined $param{'-filter'}) {

        $self->{_client}->GET($path."?query=".$param{'-query'}."&filter=".$param{'-filter'}."&format=".$self->{_format});

    } elsif (defined $param{'-query'}) {

        $self->{_client}->GET($path."?query=".$param{'-query'}."&format=".$self->{_format});

    } elsif (defined $param{'-filter'}) {

        $self->{_client}->GET($path."?filter=".$param{'-filter'}."&format=".$self->{_format});
    }

	&reset_params();

	return $self->{_client}->responseContent();

}

sub get_protein_info() {
	my $self  = shift;
	my %param = @_;

	my $path   = "/rest/entry/";
	my $string = '';

}

sub reset_params() {
	my $self = shift;

	$self->{_query}  = undef;
	$self->{_filter} = undef;
}


1;





























