#!/usr/bin/perl

# This template file is in the Public Domain.
# You may do anything you want with this file.
#
# $Id$
#

package Bio::SeqIO::game::featureHandler;
use Bio::SeqFeature::Generic;
use vars qw{ $AUTOLOAD @ISA };
use strict;

@ISA = qw(XML::Handler::Subs);

sub new {
    my ($class,$seq,$length,$type) = @_;
    my $self = bless {
	seq      => $seq,
	type     => $type,
	length   => $length,
	string   => '',
	feat     => {},
	feats    => [],
	comp_id  => 1,
    }, $class;
    return $self;
}

# Basic PerlSAX
sub start_document            {
    my ($self, $document) = @_;
    
    $self->{Names} = [];
    $self->{Nodes} = [];
    $self->{feats} = [];

}
sub end_document              {
    my ($self, $document) = @_;

    delete $self->{Names};

    return $self->{feats};
  }

sub start_element             {
     my ($self, $element) = @_;

     push @{$self->{Names}}, $element->{Name};
     $self->{string} = '';
         
     if ($self->in_element('bx-feature:seq_relationship')) {
	 if (defined $element->{Attributes}->{'bx-feature:seq'} && 
	     defined $self->{seq} &&
	     $element->{Attributes}->{'bx-feature:seq'} eq $self->{seq}) {
	     $self->{in_current_seq} = 'true';
	 } 
     }
     

     if ($self->in_element('bx-computation:computation')) {
       $self->{feat} = {};
       if (defined $element->{Attributes}->{'bx-computation:id'}) {
	 $self->{feat}->{computation_id} = $element->{Attributes}->{'bx-computation:id'};
       }  else {
	 $self->{feat}->{computation_id} = $self->{comp_id};
	 $self->{comp_id}++;
       }
     }
     
     if ($self->in_element('bx-feature:feature')) {
       if (defined $element->{Attributes}->{'bx-feature:id'}) {
	 $self->{feat}->{id} = $element->{Attributes}->{'bx-feature:id'};
       }
     }
     
     if ($self->in_element('bx-annotation:annotation')) {
       $self->{feat} = {};
       $self->{feat}->{annotation_id} = $element->{Attributes}->{'bx-annotation:id'};
     }


    return 0;
  }

sub end_element               {
  my ($self, $element) = @_;
  
  if ($self->in_element('bx-computation:program')) {
    $self->{string} =~ s/^\s+//g;
    $self->{string} =~ s/\s+$//;
    $self->{string} =~ s/\n//g;
    $self->{feat}->{source_tag} = $self->{string};
    }
  
  if ($self->in_element('bx-annotation:author')) {
    $self->{string} =~ s/^\s+//g;
    $self->{string} =~ s/\s+$//;
    $self->{string} =~ s/\n//g;
    $self->{feat}->{source_tag} = "Annotated by $self->{string}.";
  }
  
  if ($self->in_element('bx-feature:type')) {
    $self->{string} =~ s/^\s+//g;
    $self->{string} =~ s/\s+$//;
    $self->{string} =~ s/\n//g;
    $self->{feat}->{primary_tag} = $self->{string};
  }
  
  if ($self->in_element('bx-feature:start')) {
    $self->{string} =~ s/^\s+//g;
    $self->{string} =~ s/\s+$//;
    $self->{string} =~ s/\n//g;
    $self->{feat}->{start} = $self->{string};
  }
  
  if ($self->in_element('bx-feature:end')) {
    $self->{string} =~ s/^\s+//g;
    $self->{string} =~ s/\s+$//;
    $self->{string} =~ s/\n//g;
    $self->{feat}->{end} = $self->{string};
  }
  
  if ($self->in_element('bx-computation:score')) {
    $self->{string} =~ s/^\s+//g;
    $self->{string} =~ s/\s+$//;
    $self->{string} =~ s/\n//g;
    $self->{feat}->{score} = $self->{string};
  }
  
  if ($self->in_element('bx-feature:seq_relationship')) {
    if ($self->{feat}->{start} > $self->{feat}->{end}) {
      my $new_start = $self->{feat}->{end};
      $self->{feat}->{end} = $self->{feat}->{start};
      $self->{feat}->{start} = $new_start;
      $self->{feat}->{strand} = -1;
    } else {
      $self->{feat}->{strand} = 1;
    }
    my $new_feat = new Bio::SeqFeature::Generic(
						-start=>$self->{feat}->{start},
						-end=>$self->{feat}->{end},
						-source=>$self->{feat}->{source_tag},
						-primary=>$self->{feat}->{primary_tag},
						-score=>$self->{feat}->{score}
					       );
    if (defined $self->{feat}->{computation_id}) {
      $new_feat->add_tag_value('computation_id', $self->{feat}->{computation_id} );
    } elsif (defined $self->{feat}->{annotation_id}) {
      $new_feat->add_tag_value('annotation_id', $self->{feat}->{annotation_id} );
    }
    if (defined $self->{feat}->{id}) {
      $new_feat->add_tag_value('id', $self->{feat}->{id} );
    }

    push @{$self->{feats}}, $new_feat;
    $self->{feat} = { 
		     seqid => $self->{feat}->{curr_seqid},
		     primary_tag => $self->{feat}->{primary_tag},
		     source_tag => $self->{feat}->{source_tag},
		     computation_id => $self->{feat}->{computation_id},
		     annotation_id => $self->{feat}->{annotation_id}
		    }
  }
  
  
  pop @{$self->{Names}};
  pop @{$self->{Nodes}};
  
}

sub characters   {
    my ($self, $text) = @_;
    
    $self->{string} .= $text->{Data};
  }

sub in_element {
    my ($self, $name) = @_;
    
    return (defined $self->{Names}[-1] && 
	    $self->{Names}[-1] eq $name);
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
