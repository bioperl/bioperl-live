# BioPerl module for Bio::SeqIO::interpro
#
# Author: Jared Fox (jaredfox@ucla.edu)
# Contributor: Allen Day (allenday@ucla.edu)
# Copyright Jared Fox and Allen Day
#
# You may distribute this module under the same terms as perl itself
#
# Bio::SeqIO::interpro will parse interpro scan XML (version 1.2) and create 
# Bio::SeqFeature::Generic objects based on the contents of the XML document. 
# Bio::SeqIO::Interpro will also attach the annotation given in the XML file to the
# Bio::SeqFeature::Generic objects that it creates.

package Bio::SeqIO::interpro;
use vars qw(@ISA);
use strict;

use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use XML::DOM;
use XML::DOM::XPath;
use Bio::Seq::SeqFactory;
use Bio::Annotation::Collection;
use Bio::Annotation::DBLink;
use Data::Dumper;

@ISA = qw(Bio::SeqIO);

my $idcounter = {};  # Used to generate unique id values
my $nvtoken = ": ";  # The token used if a name/value pair has to be stuffed
                     # into a single line

=head1 METHODS

=cut

=head2 next_seq

 Title   : next_seq
 Usage   : my $bioSeqObj = $stream->next_seq
 Function: Retrieves the next sequence from a SeqIO::interpro stream.
 Returns : A reference to a Bio::Seq::RichSeq object
 Args    : 

=cut

sub next_seq {
  my $self = shift;
  my ($desc);
  my $bioSeq = $self->sequence_factory->create(-verbose =>$self->verbose());

  my $zinc = "(\"zincins\")";
  my $wing = "\"Winged helix\"";
  my $finger = "\"zinc finger\"";


  my $xml_fragment = undef;
  while(my $line = $self->_readline()){

    my $where = index($line, $zinc);
    my $wherefinger = index($line, $finger);
    my $finishedline = $line;
    my $wingwhere = index($line, $wing);

    #the interpro XML is not fully formed, so we need to convert the extra double quotes
    #and ampersands into the appropriate XML chracter codes
    if($where > 0){
      my @linearray = split /$zinc/, $line;
      $finishedline = join "&quot;zincins&quot;", $linearray[0], $linearray[2];
    }
    if(index($line, "&") > 0){
      my @linearray = split /&/, $line;
      $finishedline = join "&amp;", $linearray[0], $linearray[1];
    }
    if($wingwhere > 0){
      my @linearray = split /$wing/, $line;
      $finishedline = join "&quot;Winged helix&quot;", $linearray[0], $linearray[1];
    }

    $xml_fragment .= $finishedline;
    last if $finishedline =~ m!</protein>!;
  }

  return undef unless $xml_fragment =~ /<protein/;

  $self->parse_xml($xml_fragment);

  my $dom = $self->dom;

  my ($protein_node) = $dom->findnodes('/protein');
  my @interproNodes = $protein_node->findnodes('/protein/interpro');
  for(my $interpn=0; $interpn<scalar(@interproNodes); $interpn++){
    my $ipnlevel = join "", "/protein/interpro[", $interpn + 1, "]";
    my @matchNodes = $protein_node->findnodes($ipnlevel);
    for(my $match=0; $match<scalar(@matchNodes); $match++){
      my $matlevel = join "", "/protein/interpro[", $interpn+1, "]/match[", $match+1, "]/location";
      my @locNodes = $protein_node->findnodes($matlevel);

#      $self->warn(join '*', map { $_->getAttribute('score') } @locnodes);

        my @seqFeatures = map { Bio::SeqFeature::Generic->new(
                               -start => $_->getAttribute('start'), 
                               -end => $_->getAttribute('end'), 
                               -score => $_->getAttribute('score'), 
                               -source_tag => 'IPRscan',
                               -primary_tag => 'region',
  #                            -source_tag => $interproNodes[$interpn]->getAttribute('id'), 
                               -display_name => $interproNodes[$interpn]->getAttribute('name'),
                               -seq_id => $protein_node->getAttribute('id'),
                                ),
        } @locNodes;
        foreach my $seqFeature (@seqFeatures){
          #my $annotationCollection = Bio::Annotation::Collection->new;
          my $annotation1 = Bio::Annotation::DBLink->new;
          $annotation1->database($matchNodes[$match]->getAttribute('dbname'));
          $annotation1->primary_id($matchNodes[$match]->getAttribute('id'));
          $annotation1->comment($matchNodes[$match]->getAttribute('name'));
          $seqFeature->annotation->add_Annotation('dblink',$annotation1);

          my $annotation2 = Bio::Annotation::DBLink->new;
          $annotation2->database('INTERPRO');
          $annotation2->primary_id($interproNodes[$interpn]->getAttribute('id'));
          $annotation2->comment($interproNodes[$interpn]->getAttribute('name'));
          $seqFeature->annotation->add_Annotation('dblink',$annotation2)
        }
        $bioSeq->add_SeqFeature(@seqFeatures);
     }
  }
  my $accession = $protein_node->getAttribute('id');
  my $displayname = $protein_node->getAttribute('name');
  $bioSeq->accession($accession);
  $bioSeq->display_name($displayname);
  return $bioSeq;
}

sub _initialize {
  my($self,@args) = @_;

  $self->SUPER::_initialize(@args);
  # hash for functions for decoding keys.
  $self->{'_func_ftunit_hash'} = {}; 

  my %param = @args;  # From SeqIO.pm
  @param{ map { lc $_ } keys %param } = values %param; # lowercase keys


  my $line = undef;
  #fast forward to first <protein/> record.
  while($line = $self->_readline()){
    if($line =~ /<protein/){
      $self->_pushback($line);
      last;
    }
  }

  $self->xml_parser( XML::DOM::Parser->new() );

  $self->sequence_factory( new Bio::Seq::SeqFactory
                           ( -verbose => $self->verbose(),
                             -type => 'Bio::Seq::RichSeq'))
    if( ! defined $self->sequence_factory );
}

sub sequence_factory {
  my $self = shift;
  my $val = shift;

  $self->{'sequence_factory'} = $val if defined($val);
  return $self->{'sequence_factory'};
}

sub xml_parser {
  my $self = shift;
  my $val = shift;

  $self->{'xml_parser'} = $val if defined($val);
  return $self->{'xml_parser'};
}

sub parse_xml {
  my ($self,$xml) = @_;
  $self->dom( $self->xml_parser->parse($xml) );
  return 1;
}

sub dom {
  my $self = shift;
  my $val = shift;

  $self->{'dom'} = $val if defined($val);
  return $self->{'dom'};
}


1;
