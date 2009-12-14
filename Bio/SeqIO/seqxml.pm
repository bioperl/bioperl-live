# $Id$
# BioPerl module for Bio::SeqIO::seqxml
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Dave Messina <dmessina@cpan.org>
#
# Copyright Dave Messina
#
# You may distribute this module under the same terms as perl itself
# _history
# December 2009  initial version

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::seqxml - SeqXML sequence input/output stream

=head1 SYNOPSIS

  # Do not use this module directly.  Use it via the Bio::SeqIO class.

  use Bio::SeqIO;
  my $seqio = Bio::SeqIO->new(-format => 'seqxml',
                              -file   => 'my_seqs.xml');

  my $seq_object = $seqio->next_seq;

  print join("\t", 
             $seq_object->display_id,
             $seq_object->description,
             $seq_object->seq,           
            ), "\n";

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from SeqXML format.
For more information on the SeqXML standard, visit L<http://www.seqxml.org>.

This module is based in part (particularly the XML-parsing part) on 
Bio::TreeIO::phyloxml by Mira Han.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHORS - Dave Messina

Email: I<dmessina@cpan.org>

=head1 CONTRIBUTORS


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::seqxml;

use strict;

use Bio::Seq;
use Bio::Seq::SeqFactory;
use Bio::Species;
use XML::LibXML;
use XML::LibXML::Reader;

use base qw(Bio::SeqIO);

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : L<Bio::Seq::RichSeq> object, or nothing if no more available
 Args    : none

=cut

sub next_seq {
    my ($self) = @_;
    my $reader = $self->{'_reader'};
    my $entry;

    while ( $reader->read ) {

        # we're done if we hit </entry>
        if ( $reader->nodeType == XML_READER_TYPE_END_ELEMENT ) {
            if ( $reader->name eq 'entry' ) {
                $entry = $self->end_element_entry();
                last;
            }
        }
        $self->processXMLnode;
    }

    return $entry;
}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq(@seq)
 Function: Writes the $seq object into the stream
 Returns : 1 for success and 0 for error
 Args    : Array of 1 or more L<Bio::PrimarySeqI> objects

=cut

sub write_seq {
    my ( $self, $seq ) = @_;
    $self->throw("SeqXML writing not implemented yet.");
}

=head2 _initialize

 Title   : _initialize
 Usage   : $self->_initialize(@args) 
 Function: constructor (for internal use only)
 Returns : none
 Args    : none
 Throws  : Exception if XML::LibXML::Reader is not initialized

=cut

sub _initialize {
    my ( $self, @args ) = @_;

    $self->SUPER::_initialize(@args);
    if ( !defined $self->sequence_factory ) {
        $self->sequence_factory(
            Bio::Seq::SeqFactory->new(
                -verbose => $self->verbose(),
                -type    => 'Bio::Seq::RichSeq',
            )
        );
    }

    # reading in SeqXML
    if ( $self->mode eq 'r' ) {
        if ( $self->_fh ) {
            $self->{'_reader'} = XML::LibXML::Reader->new(
                IO        => $self->_fh,
                no_blanks => 1
            );
        }
        if ( !$self->{'_reader'} ) {
            $self->throw("XML::LibXML::Reader not initialized");
        }

        # data structures used during parsing
        ## holds version and source data
        $self->{'_seqxml_metadata'} = {};
        ## holds data temporarily during parsing
        $self->{'_current_entry_data'} = {};

        $self->_initialize_seqxml_node_methods();
    }

    # writing out SeqXML
    elsif ( $self->mode eq 'w' ) {

        # print default lines
        $self->_print( '<?xml version="1.0" encoding="UTF-8"?>', "\n" );
        $self->_print(
'<seqxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.seqxml.org http://www.seqxml.org/1.00/seqxml.xsd" xmlns="http://www.seqxml.org">',
            "\n"
        );
    }

}

=head2 _initialize_seqxml_node_methods

 Title   : _initialize_seqxml_node_methods
 Usage   : $self->_initialize_xml_node_methods
 Function: sets up code ref mapping of each seqXML node type
           to a method for processing that node type 
 Returns : none
 Args    : none

=cut

sub _initialize_seqxml_node_methods {
    my ($self) = @_;

    my %start_elements = (
        'seqXML'        => \&element_seqXML,
        'entry'         => \&element_entry,
        'species'       => \&element_species,
        'description'   => \&element_description,
        'rnaSeq'        => \&element_rnaSeq,
        'dnaSeq'        => \&element_dnaSeq,
        'aaSeq'         => \&element_aaSeq,
        'alternativeID' => \&element_alternativeID,
        'property'      => \&element_property,
    );
    $self->{'_start_elements'} = \%start_elements;

    my %end_elements = (
        'seqXML'        => \&end_element_seqXML,
        'entry'         => \&end_element_entry,
        'species'       => \&end_element_default,
        'description'   => \&end_element_default,
        'rnaSeq'        => \&end_element_rnaSeq,
        'dnaSeq'        => \&end_element_dnaSeq,
        'aaSeq'         => \&end_element_aaSeq,
        'alternativeID' => \&end_element_alternativeID,
        'property'      => \&end_element_default,
    );
    $self->{'_end_elements'} = \%end_elements;

}

=head1 Methods for parsing the XML document

=cut

=head2 processXMLNode

 Title   : processXMLNode
 Usage   : $seqio->processXMLNode
 Function: reads the XML node and processes according to the node type
 Returns : none
 Args    : none

=cut

sub processXMLnode {
    my ($self)   = @_;
    my $reader   = $self->{'_reader'};
    my $nodetype = $reader->nodeType;

    if ( $nodetype == XML_READER_TYPE_ELEMENT ) {
        $self->{'_current_element_name'} = $reader->name;

        if ( exists $self->{'_start_elements'}->{ $reader->name } ) {
            my $method = $self->{'_start_elements'}->{ $reader->name };
            $self->$method();
        }
        else {
            my $name = $reader->name;
            $self->warn("unexpected start element encountered: $name");
        }

        # if ( $reader->isEmptyElement ) {
        #     # element is complete
        #     # set nodetype so it can jump and
        #     # do procedures for XML_READER_TYPE_END_ELEMENT
        #     $nodetype = XML_READER_TYPE_END_ELEMENT;
        # }

    }
    elsif ( $nodetype == XML_READER_TYPE_TEXT ) {

        # store key-value pair of element name and the corresponding text
        my $name = $self->{'_current_element_name'};
        $self->{'_current_entry_data'}->{$name} = $reader->value;

    }
    elsif ( $nodetype == XML_READER_TYPE_END_ELEMENT ) {
        if ( exists $self->{'_end_elements'}->{ $reader->name } ) {
            my $method = $self->{'_end_elements'}->{ $reader->name };
            $self->$method();
        }
        else {
            my $name = $reader->name;
            $self->warn("unexpected end element encountered: $name");
        }
        $self->{'_current_element_name'} = {};    # empty current element name
    }
    else {
        $self->throw(
            "unexpected node type " . $nodetype,
            " encountered (name: ",
            $reader->name, ")\n"
        );
    }

    if ( $self->debug ) {
        printf "%d %d %s %d\n",
          (
            $reader->depth, $reader->nodeType,
            $reader->name,  $reader->isEmptyElement
          );
    }
}

=head2 processAttribute

 Title   : processAttribute
 Usage   : $seqio->processAttribute(\%hash_for_attribute);
 Function: reads the attributes of the current element into a hash
 Returns : none
 Args    : hash reference where the attributes will be stored.

=cut

sub processAttribute {
    my ( $self, $data ) = @_;
    my $reader = $self->{'_reader'};

    # several ways of reading attributes:
    # read all attributes:
    if ( $reader->moveToFirstAttribute ) {
        do {
            $data->{ $reader->name() } = $reader->value;
        } while ( $reader->moveToNextAttribute );
        $reader->moveToElement;
    }
}

=head2 element_seqXML

 Title   : element_seqXML
 Usage   : $self->element_seqXML
 Function: processes the opening <seqXML> node
 Returns : none
 Args    : none

=cut

sub element_seqXML {
    my ($self) = @_;
    my $reader = $self->{'_reader'};

    # reset for every new <seqXML></seqXML> block

    if ( $reader->hasAttributes() ) {
        $self->processAttribute( $self->{'_seqxml_metadata'} );
    }
    else {
        $self->throw("no SeqXML metadata!");
    }
}

=head2 element_entry

 Title   : element_entry
 Usage   : $self->element_entry
 Function: processes a sequence <entry> node
 Returns : none
 Args    : none
 Throws  : Exception if sequence ID is not present in <entry> element

=cut

sub element_entry {
    my ($self) = @_;
    my $reader = $self->{'_reader'};

    if ( $reader->hasAttributes() ) {
        $self->processAttribute( $self->{'_current_entry_data'} );
    }
    else {
        $self->throw("no sequence ID!");
    }
}

=head2 element_species

 Title   : element_entry
 Usage   : $self->element_entry
 Function: processes a <species> node, creating a Bio::Species object
 Returns : none
 Args    : none
 Throws  : Exception if <species> tag exists but is empty,
           or if the attributes 'name' or 'ncbiTaxID' are undefined

=cut

sub element_species {
    my ($self) = @_;
    my $reader = $self->{'_reader'};
    my $data   = $self->{'_current_entry_data'};

    my $species_data = {};
    my $species_obj;

    if ( $reader->hasAttributes() ) {
        $self->processAttribute($species_data);
    }
    else {
        $self->throw("no species information!");
    }

    if (   defined $species_data->{'name'}
        && defined $species_data->{'ncbiTaxID'} )
    {
        $species_obj =
          Bio::Species->new( -ncbi_taxid => $species_data->{'ncbiTaxID'}, );
        $species_obj->species( $species_data->{'name'} );
        $data->{'species'} = $species_obj;
    }
    else {
        $self->throw("<species> attributes name and ncbiTaxID are undefined");
    }

}

=head2 element_description

 Title   : element_description
 Usage   : $self->element_description
 Function: processes a sequence <description> node
 Returns : none
 Args    : none

=cut

sub element_description {
    my ($self) = @_;
    my $reader = $self->{'_reader'};
}

=head2 element_rnaSeq

 Title   : element_rnaSeq
 Usage   : $self->element_rnaSeq
 Function: processes a sequence <rnaSeq> node
 Returns : none
 Args    : none

=cut

sub element_rnaSeq {
    my ($self) = @_;
    my $reader = $self->{'_reader'};

    my $data = $self->{'_current_entry_data'};
    $data->{'alphabet'} = 'rna';
    $data->{'sequence'} = $data->{'rnaSeq'};

}

=head2 element_dnaSeq

 Title   : element_dnaSeq
 Usage   : $self->element_dnaSeq
 Function: processes a sequence <dnaSeq> node
 Returns : none
 Args    : none

=cut

sub element_dnaSeq {
    my ($self) = @_;
    my $reader = $self->{'_reader'};

    my $data = $self->{'_current_entry_data'};
    $data->{'alphabet'} = 'dna';
    $data->{'sequence'} = $data->{'dnaSeq'};

}

=head2 element_aaSeq

 Title   : element_aaSeq
 Usage   : $self->element_aaSeq
 Function: processes a sequence <aaSeq> node
 Returns : none
 Args    : none

=cut

sub element_aaSeq {
    my ($self) = @_;
    my $reader = $self->{'_reader'};

    my $data = $self->{'_current_entry_data'};
    $data->{'alphabet'} = 'protein';
    $data->{'sequence'} = $data->{'aaSeq'};

}

=head2 element_alternativeID

 Title   : element_alternativeID
 Usage   : $self->element_alternativeID
 Function: processes a sequence <alternativeID> node
 Returns : none
 Args    : none

=cut

sub element_alternativeID {
    my ($self) = @_;
    my $reader = $self->{'_reader'};

}

=head2 element_property

 Title   : element_property
 Usage   : $self->element_property
 Function: processes a sequence <property> node
 Returns : none
 Args    : none

=cut

sub element_property {
    my ($self) = @_;
    my $reader = $self->{'_reader'};

}

=head2 end_element_seqXML

 Title   : end_element_seqXML
 Usage   : $self->end_element_seqXML
 Function: processes the closing </seqXML> node
 Returns : none
 Args    : none

=cut

sub end_element_seqXML {
    my ($self) = @_;
    my $reader = $self->{'_reader'};

}

=head2 end_element_dnaSeq

 Title   : end_element_dnaSeq
 Usage   : $self->end_element_dnaSeq
 Function: processes a sequence <dnaSeq> node
 Returns : none
 Args    : none

=cut

sub end_element_dnaSeq {
    my ($self) = @_;
    my $reader = $self->{'_reader'};

    my $data = $self->{'_current_entry_data'};
    $data->{'alphabet'} = 'dna';
    $data->{'sequence'} = $data->{'dnaSeq'};

}

=head2 end_element_rnaSeq

 Title   : end_element_rnaSeq
 Usage   : $self->end_element_rnaSeq
 Function: processes a sequence <rnaSeq> node
 Returns : none
 Args    : none

=cut

sub end_element_rnaSeq {
    my ($self) = @_;
    my $reader = $self->{'_reader'};

    my $data = $self->{'_current_entry_data'};
    $data->{'alphabet'} = 'rna';
    $data->{'sequence'} = $data->{'rnaSeq'};
}

=head2 end_element_aaSeq

 Title   : end_element_aaSeq
 Usage   : $self->end_element_aaSeq
 Function: processes a sequence <aaSeq> node
 Returns : none
 Args    : none

=cut

sub end_element_aaSeq {
    my ($self) = @_;
    my $reader = $self->{'_reader'};

    my $data = $self->{'_current_entry_data'};
    $data->{'alphabet'} = 'protein';
    $data->{'sequence'} = $data->{'aaSeq'};

}

=head2 end_element_entry

 Title   : end_element_entry
 Usage   : $self->end_element_entry
 Function: processes the closing </entry> node, creating the Seq object
 Returns : a Bio::Seq::RichSeq object
 Args    : none
 Throws  : Exception if sequence, sequence ID, or alphabet are missing

=cut

sub end_element_entry {
    my ($self) = @_;
    my $reader = $self->{'_reader'};

    my $data = $self->{'_current_entry_data'};

    # make sure we've got at least a seq, an ID, and an alphabet
    unless ( $data->{'sequence'} ) {
        $self->throw("this entry lacks a sequence");
    }
    unless ( $data->{'id'} ) {
        $self->throw("this entry lacks an id");
    }
    unless ( $data->{'alphabet'} ) {
        $self->throw("this entry lacks an alphabet");
    }

    # create new sequnce object with minimum necessary parameters
    my $seq_obj = $self->sequence_factory->create(
        -seq        => $data->{'sequence'},
        -alphabet   => $data->{'alphabet'},
        -id         => $data->{'id'},
        -primary_id => $data->{'id'},
    );

    # add additional parameters if available
    if ( $data->{'description'} ) {
        $seq_obj->desc( $data->{'description'} );
    }
    if ( $data->{'species'} ) {
        $seq_obj->species( $data->{'species'} );
    }

    # empty the temporary data store
    $self->{'_current_entry_data'} = {};

    return $seq_obj;
}

=head2 end_element_default

 Title   : end_element_default
 Usage   : $self->end_element_default
 Function: processes all other nodes
 Returns : none
 Args    : none

=cut

sub end_element_default {
    my ($self) = @_;
    my $reader = $self->{'_reader'};

}

1;
