# This template file is in the Public Domain.
# You may do anything you want with this file.
#
# $Id$
#

package Bio::SeqIO::game::idHandler;
use vars qw{ $AUTOLOAD };
use strict;

sub new {
    my ($caller,@args) = @_;
    my $class = ref($caller) || $caller;
    my $self = bless( {
	'ids'          =>  [],
    }, $class);
    return $self;
}

# Basic PerlSAX
sub start_document            {
    my ($self, $document) = @_;
}

sub end_document              {
    my ($self, $document) = @_;
    return $self->{'ids'};
}
sub start_element             {
     my ($self, $element) = @_;

     if ($element->{'Name'} eq 'bx-seq:seq') {
       if ($element->{'Attributes'}->{'bx-seq:id'}) {
	 push @{$self->{'ids'}}, $element->{'Attributes'}->{'bx-seq:id'};
       } else {
	 if ($self->can('warn')) {
	   $self->warn('WARNING: Attribute bx-seq:id is required on bx-seq:seq. Sequence will not be parsed.');
	 } else {
	   warn('WARNING: Attribute bx-seq:id is required on bx-seq:seq. Sequence will not be parsed.');
	 }
       }
     }
    return 0;
  }

sub end_element               {
    my ($self, $element) = @_;

}
sub characters   {
    my ($self, $text) = @_;
  }


# Others
sub AUTOLOAD {
    my $self = shift;

    my $method = $AUTOLOAD;
    $method =~ s/.*:://;
    return if $method eq 'DESTROY';

    print "UNRECOGNIZED $method\n";
}

1;

__END__
