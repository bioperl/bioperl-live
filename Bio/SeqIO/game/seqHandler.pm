# This template file is in the Public Domain.
# You may do anything you want with this file.
#
# $Id$
#

package Bio::SeqIO::game::seqHandler;
use vars qw{ $AUTOLOAD @ISA };
use XML::Handler::Subs;

@ISA = qw(XML::Handler::Subs);

sub new {
    my ($class,$seq) = @_;
    my $self = bless ( {
	seq  => $seq,
    }, $class);
    return $self;
}

# Basic PerlSAX
sub start_document            {
    my ($self, $document) = @_;
    
    $self->{Names} = [];

}
sub end_document              {
    my ($self, $document) = @_;

    return new Bio::PrimarySeq( -seq => $self->{residues},
				-moltype => $self->{moltype},
				-id => $self->{seq},
				-accession => $self->{accession},
				-desc => $self->{desc},
				-length => $self->{length}
			      );

    delete $self->{Names};
  }

sub start_element             {
     my ($self, $element) = @_;

     push @{$self->{Names}}, $element->{Name};
     $string = '';

     if ($element->{Name} eq 'bx-seq:seq') {
       if ($element->{Attributes}->{'bx-seq:id'} == $self->{seq}) {
	 $self->{in_current_seq} = 'true';
	 $self->{moltype} = $element->{Attributes}->{'bx-seq:type'};
	 $self->{length} = $element->{Attributes}->{'bx-seq:length'};
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

    if ($self->{in_current_seq} eq 'true') {
      
      if ($self->in_element('bx-seq:residues')) {
	while ($string =~ s/\s+//) {};
	$self->{residues} = $string;
      }
      
      
      if ($self->in_element('bx-seq:name')) {
	$string =~ s/^\s+//g;
	$string =~ s/\s+$//;
	$string =~ s/\n//g;
	$self->{name} = $string;
      }
      
      
      if ($self->in_element('bx-link:id')  && $self->within_element('bx-link:dbxref')) {
      	$string =~ s/^\s+//g;
	$string =~ s/\s+$//;
	$string =~ s/\n//g;
	$self->{accession} = $string;
      }

      if ($self->in_element('bx-seq:description')) {
	$self->{desc} = $string;
      }
      
      if ($self->in_element('bx-seq:seq')) {
	$self->{in_current_seq} = 'false';
      }
    }
    
    pop @{$self->{Names}};
    
  }

sub characters   {
    my ($self, $text) = @_;

    $string .= $text->{Data};
  }

sub in_element {
    my ($self, $name) = @_;

    return ($self->{Names}[-1] eq $name);
}

sub within_element {
    my ($self, $name) = @_;

    my $count = 0;
    foreach my $el_name (@{$self->{Names}}) {
	$count ++ if ($el_name eq $name);
    }

    return $count;
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
