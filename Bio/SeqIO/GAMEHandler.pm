# This template file is in the Public Domain.
# You may do anything you want with this file.
#
# $Id$
#
=head1 NAME

Bio::SeqIO::Bioxml  Parses bioxml 0.3 and higher into Bio::Seq objects. - 
Sorry, no namespaces yet!

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::SeqIO class.

=head1 DESCRIPTION

This object can transform Bio::Seq objects from (and eventually to) bioxml game
versions 0.3 and higher.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioxml-dev@bioxml.org        - Technical discussion - Moderate volume
  bioxml-announce@bioxml.org   - General Announcements - Pretty dead
  http://www.bioxml.org/MailingLists/         - About the mailing lists

=head1 AUTHORS - Brad Marshall & Ewan Birney & Lincoln Stein

Email: bradmars@yahoo.com
       birney@sanger.ac.uk
       lstein@cshl.org


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut



package Bio::SeqIO::GAMEHandler;
use vars qw{ $AUTOLOAD };

sub new {
    my $class = shift;
    my $self = {};
    my $seqs;
    my @feats = [];
    my $string;
    my $feat = {};
    return bless ($self, $class);
}

# Basic PerlSAX
sub start_document            {
    my ($self, $document) = @_;

    $self->{Names} = [];
    $self->{Nodes} = [];
    $self->{Seqs} = {};
    $self->{feats} = [];
}
sub end_document              {
    my ($self, $document) = @_;

    delete $self->{Names};
    delete $self->{Nodes};

    return();
}
sub start_element             {
     my ($self, $element) = @_;

    $string = '';
    push @{$self->{Names}}, $element->{Name};
    push @{$self->{Nodes}}, $element;

    if ($self->in_element('seq')) {
      if ($self->{Seqs}->{$self->{Nodes}[-1]->{Attributes}->{id}}->{type} eq  'dna') {
	  $feat->{strand} = 1;
	} elsif ($self->{Seqs}->{$self->{Nodes}[-1]->{Attributes}->{id}}->{type} eq  'aa') {
	  $feat->{strand} = 0;
	  $self->{Seqs}->{$self->{Nodes}[-1]->{Attributes}->{id}}->{type} = 'protein';
	} elsif ($self->{Seqs}->{$self->{Nodes}[-1]->{Attributes}->{id}}->{type} eq  'rna') {
	  $feat->{strand} = 0;
	} else {
	  if ($self->can('warn')) {
	    $self->warn('Sequence type is either missing or not one of dna, rna, or aa.  bioperl will guess.');
	  } else {
	    warn('Sequence type is either missing or not one of dna, rna, or aa.  bioperl will guess.');
	  }

	    $self->{Seqs}->{$self->{Nodes}[-1]->{Attributes}->{id}}->{type} = '';
	}

      if ($self->{Nodes}[-1]->{Attributes} -> {length}) {
	$self->{Seqs}->{$self->{Nodes}[-1]->{Attributes}->{id}} -> {length} = $self->{Nodes}[-1]->{Attributes}->{length};
      }
    }

    if ($self->in_element('seq_relationship')) {
      if ($self->{Nodes}->[-1]->{Attributes}->{seq}) {
	$feat->{seqid} = $self->{Nodes}->[-1]->{Attributes}->{seq};
      } 
    }

     if ($self->in_element('computation')) {
       $feat->{curr_seqid} = $self->{Nodes}[-1]->{Attributes}->{seq};
       $feat->{seqid} = $feat->{curr_seqid}
     }
     
     if ($self->in_element('annotation')) {
       $feat->{curr_seqid} = $self->{Nodes}[-1]->{Attributes}->{seq};
       $feat->{seqid} = $feat->{curr_seqid}
     }
     


    return 0;
  }

sub end_element               {
    my ($self, $element) = @_;

    if ($self->in_element('residues')) {
      while ($string =~ s/\s+//) {};
       $self->{Seqs}->{$self->{Nodes}[-2]->{Attributes}->{id}} -> {seq} = $string;
   }
    
    if ($self->in_element('program')) {
	$string =~ s/^\s+//g;
	$string =~ s/\s+$//;
	$string =~ s/\n//g;
	$feat->{source_tag} = $string;
    }

    if ($self->in_element('author')) {
      if ($self->within_element('annotation')){
	$string =~ s/^\s+//g;
	$string =~ s/\s+$//;
	$string =~ s/\n//g;
	$feat->{source_tag} = "Annotated by $string .";
      }
    }
    
    if ($self->in_element('type')) {
      $string =~ s/^\s+//g;
      $string =~ s/\s+$//;
      $string =~ s/\n//g;
      $feat->{primary_tag} = $string;
    }

    if ($self->in_element('start')) {
      $string =~ s/^\s+//g;
      $string =~ s/\s+$//;
      $string =~ s/\n//g;
      $feat->{start} = $string;
    }

    if ($self->in_element('end')) {
      $string =~ s/^\s+//g;
      $string =~ s/\s+$//;
      $string =~ s/\n//g;
      $feat->{end} = $string;
    }

    if ($self->in_element('score')) {
      $string =~ s/^\s+//g;
      $string =~ s/\s+$//;
      $string =~ s/\n//g;
      $feat->{score} = $string;
    }

    if ($self->in_element('seq_relationship')) {
      if ($feat->{start} > $feat->{end}) {
	$feat->{start} = $self->{Seqs}->{$feat->{seqid}}->{length} - $feat->{start};
	$feat->{end} = $self->{Seqs}->{$feat->{seqid}}->{length} - $feat->{end};
	if ($self->{Seqs}->{$feat->{seqid}}->{type} eq  'dna') {
	  $feat->{strand} = -1;
	} elsif ($self->{Seqs}->{$feat->{seqid}}->{type} eq  'protein'){
	  $self->warn("adding amino acid $feat->{seqid} with start less than end");
	  $feat->{strand} = 0;
	} elsif ($self->{Seqs}->{$feat->{seqid}}->{type} eq  'rna') {
	  warn("adding rna $feat->{seqid} with start less than end");
	  $feat->{strand} = 0;
	} else {
	  $feat->{strand} = 0;
	}
      } else {
	$feat->{start}++;
	$feat->{end}++;
	if ($self->{Seqs}->{$feat->{seqid}}->{type} eq  'dna') {
	  $feat->{strand} = 1;
	} elsif ($self->{Seqs}->{$feat->{seqid}}->{type} eq  'protein') {
	  $feat->{strand} = 0;
	} elsif ($self->{Seqs}->{$feat->{seqid}}->{type} eq  'rna') {
	  $feat->{strand} = 0;
	} else {
	  $feat->{strand} = 0;
	}
      }
      push @{$self->{Seqs}->{$feat->{seqid}}->{features}}, $feat;
      $feat = { seqid => $feat->{curr_seqid},
		primary_tag => $feat->{primary_tag},
		source_tag => $feat->{source_tag}
	      }
    }

    pop @{$self->{Names}};
    pop @{$self->{Nodes}};

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

sub getSeqs {
  my $self=shift;
  return $self->{Seqs};
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
