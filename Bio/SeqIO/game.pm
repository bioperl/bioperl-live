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

Bio::SeqIO::Bioxml  Parses game 0.1 and higher into and out of 
 Bio::Seq objects.  

=head1 SYNOPSIS

To use this module you need XML::Parser, XML::Parser::PerlSAX
and XML::Writer.

Do not use this module directly.  Use it via the Bio::SeqIO class.

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from bioxml seq,
computation, feature and annotation dtds,versions 0.1 and higher.
These can be found at http://www.bioxml.org/dtds/current.  It does
this using the idHandler, seqHandler and featureHandler modules you
should have gotten with this one.

The idea is that any bioxml features can be turned into bioperl
annotations.  When Annotations and computations are parsed in, they
gain additional info in the bioperl SeqFeature tag attribute.  These
can be used to reconstitute a computation or annotation by the bioxml
with the bx-handler module when write_seq is called.

If you use this to write SeqFeatures that were not generated from
computations or annotations, it will output a list of bioxml features.
Some data may be lost in this step, since bioxml features just have a
span, type and description - nothing about the anlysis performed.

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

package Bio::SeqIO::game;
use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::Object

use Bio::SeqIO;
use Bio::SeqIO::game::seqHandler;
use Bio::SeqIO::game::featureHandler;
use Bio::SeqIO::game::idHandler;
use XML::Parser::PerlSAX;
use Bio::SeqFeature::Generic;
use XML::Writer;

use Bio::Seq;

@ISA = qw(Bio::SeqIO);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

#----------
sub new {
#----------
    my($class, @param) = @_;
    my $self = {};
    bless $self, ref($class) || $class;
    $self->_initialize(@param);
    $self;
}


sub _initialize {
  my($self,@args) = @_;
  my $xmlfile ="";
  
  $self->{counter} = 0;
  $self->{id_counter} = 1;  
  
  ($self->{file} ) = $self->_rearrange( [ qw(FILE) ], @args);
  $self->throw("did not specify a file to read, Filehandle suport is not implemented currently") if( !defined $self->{file});
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

  unless (defined @{$self->{seqs}}) {

      eval {
	my $handler = Bio::SeqIO::game::idHandler->new();
	my $options = {Handler=>$handler};
	my $parser = XML::Parser::PerlSAX->new($options);
	$self->{seqs} = $parser->parse(Source => { SystemId => $self->{file} });
      };
      if ($@) {
	$self->warn("There was an error parsing the xml document.  It may not be well-formed or empty.\n$@");
	return 0;
      }
  }
  my $seq = shift(@{$self->{seqs}});
  if ($seq) {
    my $handler = Bio::SeqIO::game::seqHandler->new($seq);
    my $options = {Handler=>$handler};
    my $parser = XML::Parser::PerlSAX->new($options);
    my $pseq = $parser->parse(Source => { SystemId => $self->{file} });

    my $handler = Bio::SeqIO::game::featureHandler->new($pseq->length(), $pseq->moltype());
    my $options = {Handler=>$handler};

    my $parser = XML::Parser::PerlSAX->new($options);
    my $features = $parser->parse(Source => { SystemId => $self->{file} });
    my $seq = Bio::Seq->new();
    foreach my $feature (@{$features}) {
      $seq->add_SeqFeature($feature);
    }
    $seq->primary_seq($pseq);
    return $seq;
  }
}

=head2 next_primary_seq

 Title   : next_primary_seq
 Usage   : $seq = $stream->next_primary_seq()
 Function: returns the next primary sequence (ie no seq_features) in the stream
 Returns : Bio::Primary_Seq object
 Args    : NONE

=cut

sub next_primary_seq {
  my $self=shift;

  unless (defined @{$self->{seqs}}) {
    
    eval {
      my $handler = Bio::SeqIO::game::idHandler->new();
      my $options = {Handler=>$handler};
      my $parser = XML::Parser::PerlSAX->new($options);
      $self->{seqs} = $parser->parse(Source => { SystemId => $self->{file} });
    };
    if ($@) {
      $self->warn("There was an error parsing the xml document.  It may not be well-formed or empty.");
      return 0;
    }
  }

  my $seq = shift(@{$self->{seqs}});
  if ($seq) {
    my $handler = Bio::SeqIO::game::seqHandler->new($seq);
    my $options = {Handler=>$handler};
    my $parser  = XML::Parser::PerlSAX->new($options);
    my $pseq = $parser->parse(Source => { SystemId => $self->{file} });
    return $pseq;
  }
}


=head2 write_seq

 Title   : write_seq
 Usage   : Not Yet Implemented! $stream->write_seq(@seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq object


=cut

sub write_seq {
  my ($self,@seqs) = @_;

  my $bxfeat = "http://www.bioxml.org/dtds/feature/v0_1";
  my $bxann = "http://www.bioxml.org/dtds/annotation/v0_1";
  my $bxcomp = "http://www.bioxml.org/dtds/computation/v0_1";
  my $bxgame = "http://www.bioxml.org/dtds/game/v0_1";
  my $bxlink = "http://www.bioxml.org/dtds/link/v0_1";
  my $bxseq = "http://www.bioxml.org/dtds/seq/v0_1";
  
  my $writer = new XML::Writer(FILENAME => $self->{file},
  		               NAMESPACES => 1,
			       DATA_MODE => 1,
			       DATA_INDENT => 4,
			       PREFIX_MAP => {$bxfeat => 'bx-feature',
					      $bxann => 'bx-annotation',
					      $bxcomp => 'bx-computation',
					      $bxgame => 'bx-game',
					      $bxlink => 'bx-link',
					      $bxseq => 'bx-seq'
					     });
  
  
  $writer->xmlDecl("UTF-8");
  $writer->doctype("bx-game:game", 'game', "http://www.bioxml.org/dtds/current/game.dtd");
  
  $writer ->startTag([$bxgame, 'game']);
  foreach my $seq (@seqs) {
     $writer ->startTag([$bxseq, 'seq'], 
      		       [$bxseq, 'id'] => $seq->display_id,
		       [$bxseq, 'length'] => $seq->length,
		       [$bxseq, 'type'] => $seq->moltype);
    if ($seq->length > 0) {
      $writer ->startTag([$bxseq, 'residues']);
      $writer->characters($seq->seq);
      $writer ->endTag([$bxseq, 'residues']);
    }
    $writer->endTag([$bxseq, 'seq']);
    
    my @feats = $seq->all_SeqFeatures;
    
     my $features;
    
    foreach my $feature (@feats) {
      if ($feature->has_tag('annotation_id')) {
	my @ann_id = $feature->each_tag_value('annotation_id');
	push (@{$features->{annotations}->{@ann_id[0]}}, $feature);
      } elsif ($feature->has_tag('computation_id')) {
	my @comp_id = $feature->each_tag_value('computation_id');
	push (@{$features->{computations}->{@comp_id[0]}}, $feature);
      } else {
	push (@{$features->{everybody_else}}, $feature);
      }
    }
    
    foreach my $key (keys %{$features->{annotations}}) {
      $writer->startTag([$bxann, 'annotation'],
			[$bxann, 'id']=>$key
		       );
      $writer->startTag([$bxann, 'seq_link']);
      $writer->startTag([$bxlink, 'link']);
      $writer->emptyTag([$bxlink, 'ref_link'],
		       [$bxlink, 'ref'] => $seq->display_id());
      $writer->endTag([$bxlink, 'link']);
      $writer->endTag([$bxann, 'seq_link']);
      $writer->startTag([$bxann, 'feature_set']);
      foreach my $feature(@{$features->{annotations}->{$key}}) {
	$self->__draw_feature($writer, $feature, $seq);
      }
      $writer->endTag([$bxann, 'feature_set']);
      $writer->endTag([$bxann, 'annotation']);
    }
    
    foreach my $key (keys %{$features->{computations}}) {
      $writer->startTag([$bxcomp, 'computation'],
			[$bxcomp, 'id']=>$key
		       );
      $writer->startTag([$bxcomp, 'seq_link']);
      $writer->startTag([$bxlink, 'link']);
      $writer->emptyTag([$bxlink, 'ref_link'],
		       [$bxlink, 'ref'] => $seq->display_id());
      $writer->endTag([$bxlink, 'link']);
      $writer->endTag([$bxcomp, 'seq_link']);
      $writer->startTag([$bxcomp, 'feature_set']);
      foreach my $feature(@{$features->{computations}->{$key}}) {
	$self->__draw_feature($writer, $feature, $seq);
      }
      $writer->endTag([$bxcomp, 'feature_set']);
      $writer->endTag([$bxcomp, 'computation']);
    }

    foreach my $feature(@{$features->{everybody_else}}) {
      $self->__draw_feature($writer, $feature, $seq);
    }
  }
   $writer->endTag([$bxgame, 'game']);
}

sub __draw_feature {
  my ($self, $writer, $feature, $seq, $parent) = @_;
  my ($subfeature, $subset, @subfeatures, $score, $score_val, $score_no);
  my $bxfeat = "http://www.bioxml.org/dtds/feature/v0_1";

  if (!$feature->has_tag('id')) {
    #$feature->remove_tag('id');
    $feature->add_tag_value('id', $self->{id_counter});
    $self->{id_counter}++;
  }

  my @id = $feature->each_tag_value('id');
  if ($parent == "") {
	  $writer->startTag([$bxfeat, 'feature'],
			    [$bxfeat, 'id'] => $id[0]
			   );
  } else {
	  $writer->startTag([$bxfeat, 'feature'],
			    [$bxfeat, 'id'] => $id[0],
			    [$bxfeat, 'parent'] => $parent
			   );
  }

  $writer->startTag([$bxfeat, 'type']);
  $writer->characters($feature->primary_tag());
  $writer->endTag([$bxfeat, 'type']);

  foreach $score ($feature->all_tags()) {
    next if ($score eq 'id');
    $writer->startTag([$bxfeat, 'score'],
                      [$bxfeat, 'type'] => $score 
                     );
    $score_no = 0;
    foreach $score_val ($feature->each_tag_value($score)) {
       $writer->characters(' ') if ($score_no > 0);
       $writer->characters($score_val);
       $score_no++;
    }
    $writer->endTag([$bxfeat, 'score']);
  }

  $writer->startTag([$bxfeat, 'seq_relationship'],
		    [$bxfeat, 'seq'] => $seq->display_id,
		    [$bxfeat, 'type'] => 'query'
		   );
  
  $writer->startTag([$bxfeat, 'span']);
  $writer->startTag([$bxfeat, 'start']);
  $writer->characters($feature->start());
  $writer->endTag([$bxfeat, 'start']);
  $writer->startTag([$bxfeat, 'end']);
  $writer->characters($feature->end());
  $writer->endTag([$bxfeat, 'end']);
  $writer->endTag([$bxfeat, 'span']);
  $writer->endTag([$bxfeat, 'seq_relationship']);
  $writer->endTag([$bxfeat, 'feature']);

  #proces subseqfeature's, exons, introns, promotors, whatever...

  foreach $subset (('exons', 'promoters','poly_A_sites','utrs','introns','sub_SeqFeature')) {
    #determine if it exists
    if (@subfeatures = eval ( '$feature->' . $subset . '()' )) {
      foreach $subfeature (@subfeatures) {
         $self->__draw_feature ($writer, $subfeature, $seq, $id[0]);
      }        
    }
  }
}

1;

