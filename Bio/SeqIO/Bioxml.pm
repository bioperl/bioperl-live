#
# BioPerl module for Bio::SeqIO::fasta
#
# Cared for by Brad Marshall <bradmars@yahoo.com>
#         
#
# Copyright Ewan Birney & Lincoln Stein & Brad Marshall
#
# You may distribute this module under the same terms as perl itself
# _history
# June 25, 2000     written by Brad Marshall
#
# POD documentation - main docs before the code

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

# Let the code begin...

package Bio::SeqIO::bioxml;
use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::Object

use Bio::SeqIO;
use Bio::SeqIO::GAMEHandler;
use XML::Parser::PerlSAX;
use Bio::SeqFeature::Generic;
use Bio::PrimarySeq;
use Bio::Seq;

@ISA = qw(Bio::SeqIO);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  my $xmlfile ="";

  $self->{counter} = 0;
  
  open (FILE, @args[1]);
  while (<FILE>){
    $xmlfile .= $_;
  }

  eval {
    my $handler = Bio::SeqIO::GAMEHandler->new();
    my $options = {Handler=>$handler};
    my $parser = XML::Parser::PerlSAX->new($options);
    $parser->parse($xmlfile);
    $self->{seqs} = $handler->getSeqs();
    foreach my $key (keys %{$self->{seqs}}) {
      push @{$self->{seq_keys}}, $key;
    }
  };
  if ($@) {
    $self->throw("There was an error parsing the xml document $args[1].  It may not be well-formed");
    exit(0);
  }
  
  
  return unless my $make = $self->SUPER::_initialize(@args);
}

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    : NONE

=cut

sub next_seq {
  my $self = shift;
  
  return 0 if (scalar @{$self->{seq_keys}} < $self->{counter});
  
  my $annseq= Bio::Seq->new();
  
  
  my $pseq = 
    Bio::PrimarySeq->new(-seq => $self->{seqs}->{$self->{seq_keys}[$self->{counter}]}->{seq},
			 #-accession => $self->{seq_keys}[$self->{counter}],
			 -id => $self->{seq_keys}[$self->{counter}],
			 -moltype => $self->{seqs}->{$self->{seq_keys}[$self->{counter}]}->{type}
			);
  
  foreach my $feature (@{$self->{seqs}->{$self->{seq_keys}[$self->{counter}]}->{features}}) {
    #   if ($feature->{start} <= $feature->{end}) {
      
    my $feat = 
      new Bio::SeqFeature::Generic ( -start => $feature->{start}, 
				     -end => $feature->{end},
				     -strand => $feature->{strand},
				     -primary => $feature->{primary_tag},
				     -source => $feature->{source_tag},
				     -score  => $feature->{score}
				   );
    
    
    $annseq->add_SeqFeature($feat);
  }  
  
  $annseq->primary_seq($pseq);
  ++$self->{counter};
  return $annseq;
}


=head2 write_seq

 Title   : write_seq
 Usage   : Not Yet Implemented! $stream->write_seq(@seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq object


=cut

sub write_seq {
   my ($self,@seq) = @_;
   $self->throw("Sorry, write seq is not yet implemented.");
   return 1;
}

1;
