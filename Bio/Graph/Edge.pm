#!/bin/perl -w
use strict;
package Bio::Graph::Edge;
use Bio::Root::Root;
use Bio::IdentifiableI;
use vars qw(@ISA);
@ISA = qw(Bio::Root::Root Bio::IdentifiableI);

##array based, not hash based ##..., therefore does not use 
#Bio::Root::Root->new().

sub new {

	my ($caller, @args) = @_;
	my $class  = ref ($caller) || $caller;
	my $self    = [];
	bless ($self, $class);

	my ($weight, $id, $nodes) = $self->_rearrange([qw(
													 WEIGHT 
													 ID
													 NODES)],
												 @args);
	$self->[0] = $nodes->[0];
	$self->[1] = $nodes->[1];
	$self->[2] = defined($weight)?$weight:undef; 
	$self->[3] = defined($id)?$id:undef; 
	return $self;

}

sub weight {
	my $self = shift;
	if (@_) {$self->[2] = shift;}
	return defined($self->[2])?$self->[2]:undef;
}

sub object_id {
	my $self            = shift;
	if (@_) {
		my $v  = shift;
		if (ref ($v)) {
			$self->throw ("Edge ID must be a text value, not a [".
							ref($v). "].");
			} 
		$self->[3] = shift;
	}
	return defined($self->[3])?$self->[3]:undef;
}

sub nodes {
	my ($self, @args) = @_;
	if (@args >= 2 ) {
		$self->[0] =  $args[0];
		$self->[1] =  $args[1];
		}
	return ($self->[0], $self->[1]);
	
}
