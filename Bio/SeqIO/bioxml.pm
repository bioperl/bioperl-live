#
# BioPerl module for Bio::SeqIO::bioxml
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#          and Lincoln Stein <lstein@cshl.org>
#
# Copyright Ewan Birney & Lincoln Stein
#
# You may distribute this module under the same terms as perl itself
# _history
# October 18, 1999  Largely rewritten by Lincoln Stein

# POD documentation - main docs before the code
# $Id$ 

=head1 NAME

Bio::SeqIO::bioxml - fasta sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::SeqIO class.

=head1 REQUIREMENTS

To use this code you need the CPAN modules XML::Parser and XML::DOM
installed locally.

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from xml sequence
files that comply to the bioxml seq tag.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to the bioxml developer's mailing list.
Your participation is much appreciated.

  bioxml-dev@bioxml.org                 the bioxml list.
  bioperl-l@bioperl.org                - General bioperl discussion
  bioperl-guts-l@bioperl.org           - Technically-oriented bioperl discussion
  http://www.bioxml.org/MailList.html  - About the bioxml mailing list
  http://www.bioperl.org/MailList.html - About the bioperl mailing lists

=head2 Reporting Bugs

For now send bug reports to the bioxml list above.

=head1 AUTHORS - Brad Marshall, Ewan Birney, Lincoln Stein

Email: bradmars@yahoo.com
       birney@ebi.ac.uk
       lstein@cshl.org


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut
#'

# Let the code begin...

package Bio::SeqIO::bioxml;
use vars qw(@ISA);
use strict;
# Object preamble - inherits from Bio::Root::Object

use Bio::SeqIO;
use Bio::Root::Object;
use XML::DOM;
use XML::Handler::BuildDOM;

@ISA = qw(Bio::SeqIO Bio::Root::Object);

sub new {
    my ($class,@args) = @_;    
    my $self = bless {}, $class;
    $self->_initialize(@args);
    return $self;
}

sub _initialize {
    my($self,@args) = @_;
    $self->SUPER::_initialize(@args);

    if (@args[2] ne '') {
	my $xmlfile = "";
	while (my $next_line = $self->_readline) {
	    $xmlfile= $xmlfile.$next_line;
	}
	my ($doc, $seq);
	eval {
	    my $parser = new XML::DOM::Parser();
	    $doc = $parser->parse($xmlfile);
	    $seq = $doc->getElementsByTagName("seq");
	};
	if ($@) {
	    $self->throw("There was an error parsing the xml document $args[1].  It may not be well-formed");
	    return 0;
	}
    }
    $self->{'_seqctr'} = 0;
}

sub _residues {
    my ($self,$sequence) = @_;
    if (defined $sequence->item($seqnum)->getElementsByTagName("residues")->item(0)){
	my $sequence=$seq->item($seqnum)->getElementsByTagName("residues")->item(0)->getFirstChild()->getNodeValue();
	$sequence =~ s/ //g;
	$sequence =~ s/\n//g;
	my $type = $seq->item($seqnum)->getElementsByTagName("residues")->item(0)->getAttribute("type");
	return ($sequence, $type);
    } 
}

sub _id {
    my $seq = shift;
    my $name = $seq->item($seqnum)->getAttribute("name");
    if (defined $name){
        return $name;
    } 
}

sub _accession {
    my $seq = shift;
    my $accession = $seq->item($seqnum)->getElementsByTagName("unique_id")->item(0);
    if (defined  $accession){
        my $accession=$accession->getFirstChild()->getNodeValue();
        $accession =~ s/^\s+//;
        $accession =~ s/\s+$//;
        return $accession;
    } 
}
    
=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    :

=cut

sub next_seq{
    my ($self ) = @_;
    if ($self->{'_seqctr'} < $seq->getLength()) {
	my $id = &_id($seq);
	my ($sequence, $type)=&_residues($seq);
	my $accession = &_accession($seq);

	$seqnum++;


	return Bio::Seq->new(-seq => $sequence,
			     -accession_number => $accession,
			     -display_id=>$id,
			     -moltype=>$type
			     );
    }
}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq(@seq)
 Function: writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq object


=cut

sub write_seq {
   my ($self,@seq) = @_;
   my $id=1;
   foreach $seq (@seq) {
     my $doc = new XML::DOM::Document();
     my $node = $doc->createElement("seq");
     $doc->appendChild($node);
     my $curr_node = $doc->getFirstChild();
     my $attr = $doc->createAttribute("id");
     $curr_node->setAttributeNode($attr);
     $curr_node->setAttribute("id", $id);
     if ($seq->display_id()) {
       my $attr = $doc->createAttribute("name");
       $curr_node->setAttributeNode($attr);
       $curr_node->setAttribute("name", $seq->display_id());
     }
     $id++;
     if ($seq->accession_number() ne "unknown") {
       my $node = $doc->createElement("dbxref");
       $curr_node->appendChild($node);
       my $curr_node= $curr_node->getFirstChild();
       my $node = $doc->createElement("database");
       $curr_node->appendChild($node);
       my $curr_node= $curr_node->getFirstChild();
       my $node = $doc->createTextNode("unknown");
       $curr_node->appendChild($node);
       my $curr_node= $curr_node->getParentNode();
       my $node = $doc->createElement("accession");
       $curr_node->appendChild($node);
       my $curr_node= $curr_node->getLastChild();
       my $node = $doc->createTextNode($seq->accession_number());
       $curr_node->appendChild($node);
       my $curr_node = $curr_node->getParentNode();
       my $curr_node = $curr_node->getParentNode();
     }
     if ($seq->seq()){
       my $node = $doc->createElement("residues");
       $curr_node->appendChild($node);
       my $curr_node = $curr_node->getLastChild();
       my $node = $doc->createTextNode($seq->seq);
       $curr_node->appendChild($node);
     }
     print $doc->toString;
   }
   return 1;
}

1;
