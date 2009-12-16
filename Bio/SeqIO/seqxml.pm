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
  
  # read a SeqXML file
  my $seqio = Bio::SeqIO->new(-format => 'seqxml',
                              -file   => 'my_seqs.xml');

  while (my $seq_object = $seqio->next_seq) {
      print join("\t", 
                 $seq_object->display_id,
                 $seq_object->description,
                 $seq_object->seq,           
                ), "\n";
  }
  
  # write a SeqXML file
  #
  # Note that you can (optionally) specify the source
  # (usually a database) and source version.
  my $seqwriter = Bio::SeqIO->new(-format        => 'seqxml',
                                  -file          => ">outfile.xml",
                                  -source        => 'Ensembl',
                                  -sourceVersion => '56');
  $seqwriter->write_seq($seq_object);

=head1 DESCRIPTION

This object can transform Bio::Seq objects to and from SeqXML format.
For more information on the SeqXML standard, visit L<http://www.seqxml.org>.

In short, SeqXML is a lightweight sequence format that takes advantage
of the validation capabilities of XML while not overburdening the user
with a strict and complicated schema.

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
use Bio::Annotation::DBLink;
use XML::LibXML;
use XML::LibXML::Reader;
use XML::Writer;

use base qw(Bio::SeqIO);

=head2 _initialize

 Title   : _initialize
 Usage   : $self->_initialize(@args) 
 Function: constructor (for internal use only).
 
           Besides the usual SeqIO arguments (-file, -fh, etc.),
           Bio::SeqIO::seqxml accepts three arguments which are used
           when writing out a seqxml file. They are all optional.
 Returns : none
 Args    : -source         => source string (usually a database name)
           -sourceVersion  => source version. The version number of the source
           -seqXMLversion  => the version of seqXML that will be used
 Throws  : Exception if XML::LibXML::Reader or XML::Writer
           is not initialized

=cut

sub _initialize {
    my ( $self, @args ) = @_;

    $self->SUPER::_initialize(@args);
    if ( !defined $self->sequence_factory ) {
        $self->sequence_factory(
            Bio::Seq::SeqFactory->new(
                -verbose => $self->verbose(),
                -type    => 'Bio::Seq',
            )
        );
    }

    # holds version and source data
    $self->{'_seqxml_metadata'} = {};

    # load any passed parameters
    my %params = @args;
    if ($params{'-sourceVersion'}) {
        $self->sourceVersion($params{'-sourceVersion'});
    }
    if ($params{'-source'}) {
        $self->source($params{'-source'});
    }
    if ($params{'-seqXMLversion'}) {
        $self->seqXMLversion($params{'-seqXMLversion'});
    }
    # reading in SeqXML
    if ( $self->mode eq 'r' ) {
        if ( $self->_fh ) {
            $self->{'_reader'} = XML::LibXML::Reader->new(
                IO        => $self->_fh,
                no_blanks => 1,
            );
        }
        if ( !$self->{'_reader'} ) {
            $self->throw("XML::LibXML::Reader not initialized");
        }

        # holds data temporarily during parsing
        $self->{'_current_entry_data'} = {};

        $self->_initialize_seqxml_node_methods();
    }

    # writing out SeqXML
    elsif ( $self->mode eq 'w' ) {
        if ( $self->_fh ) {
            $self->{'_writer'} = XML::Writer->new(
                OUTPUT      => $self->_fh,
                DATA_MODE   => 1,
                DATA_INDENT => 1,
            );
            if ( !$self->{'_writer'} ) {
                $self->throw("XML::Writer not initialized");
            }

            # write SeqXML header
            $self->{'_writer'}->xmlDecl("UTF-8");
            $self->{'_writer'}->startTag(
                'seqXML',
                'source'        => $self->source,
                'sourceVersion' => $self->sourceVersion,
                'seqXMLversion' => $self->seqXMLversion,
            );
        }
    }

}

=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : L<Bio::Seq> object, or nothing if no more available
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
    my ( $self, @seqs ) = @_;
    my $writer = $self->{'_writer'};

    foreach my $seqobj (@seqs) {
        $self->throw("Trying to write with no seq!") unless defined $seqobj;

        if ( !ref $seqobj || !$seqobj->isa('Bio::SeqI') ) {
            $self->warn(
" $seqobj is not a SeqI compliant module. Attempting to dump, but may fail!"
            );
        }

        # opening tag and ID
        my $id = $seqobj->display_id;
        if ($id) {
            $writer->startTag( 'entry', 'id' => $id );
        }
        else {
            $self->throw(" $seqobj has no ID!");
        }

        # species and NCBI taxID
        if ( $seqobj->species ) {
            my $name  = $seqobj->species->species;
            my $taxid = $seqobj->species->ncbi_taxid;
            if ( $name && ( $taxid =~ /[0-9]+/ ) ) {
                $writer->emptyTag(
                    'species',
                    'name'      => $seqobj->species->species,
                    'ncbiTaxID' => $seqobj->species->ncbi_taxid
                );
            }
            else {
                $self->throw("$seqobj has malformed species data");
            }
        }

        # description
        if ( $seqobj->desc ) {
            $writer->dataElement( 'description', $seqobj->desc );
        }

        # sequence
        if ( $seqobj->seq ) {
            my $alphabet = $seqobj->alphabet;
            my %seqtype  = (
                'rna'     => 'rnaSeq',
                'dna'     => 'dnaSeq',
                'protein' => 'aaSeq'
            );
            unless ( exists( $seqtype{$alphabet} ) ) {
                $self->throw("invalid sequence alphabet $alphabet!");
            }
            $writer->dataElement( $seqtype{$alphabet}, $seqobj->seq );
        }

        # alternative IDs
        my @dblinks = $seqobj->get_Annotations('dblink');
        foreach my $dblink (@dblinks) {
            unless ( $dblink->database && $dblink->primary_id ) {
                $self->throw("dblink $dblink is malformed");
            }
            $writer->emptyTag(
                'alternativeID',
                'source' => $dblink->database,
                'id'     => $dblink->primary_id,
            );
        }

        # properties
        my @annotations = $seqobj->get_Annotations();
        foreach my $annot_obj (@annotations) {
            next if ( $annot_obj->tagname eq 'dblink' );
            unless ( $annot_obj->tagname ) {
                $self->throw("property $annot_obj is missing a tagname");
            }
            if ( $annot_obj->value ) {
                $writer->emptyTag(
                    'property',
                    'name'  => $annot_obj->tagname,
                    'value' => $annot_obj->value,
                );
            }
            else {
                $writer->emptyTag(
                    'property',
                    'name' => $annot_obj->tagname,
                );
            }

        }

        # closing tag
        $writer->endTag('entry');

        # make sure it gets written to the file
        $self->flush if $self->_flush_on_write && defined $self->_fh;
        return 1;
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
        'seqXML'        => \&end_element_default,
        'entry'         => \&end_element_entry,
        'species'       => \&end_element_default,
        'description'   => \&end_element_default,
        'rnaSeq'        => \&end_element_rnaSeq,
        'dnaSeq'        => \&end_element_dnaSeq,
        'aaSeq'         => \&end_element_aaSeq,
        'alternativeID' => \&end_element_default,
        'property'      => \&end_element_default,
    );
    $self->{'_end_elements'} = \%end_elements;

}

=head2 source

 Title   : source
 Usage   : $self->source
 Function: gets/sets the data source in the <seqXML> header
 Returns : the data source string
 Args    : To set the source, call with a source string as the argument.

=cut

sub source {
    my ( $self, $value ) = @_;
    my $metadata = $self->{'_seqxml_metadata'};

    # set if a value is supplied
    if ($value) {
        $metadata->{'source'} = $value;
    }

    return $metadata->{'source'};
}

=head2 sourceVersion

 Title   : sourceVersion
 Usage   : $self->sourceVersion
 Function: gets/sets the data source version in the <seqXML> header
 Returns : the data source version string
 Args    : To set the source version, call with a source version string
           as the argument.

=cut

sub sourceVersion {
    my ( $self, $value ) = @_;
    my $metadata = $self->{'_seqxml_metadata'};

    # set if a value is supplied
    if ($value) {
        $metadata->{'sourceVersion'} = $value;
    }

    return $metadata->{'sourceVersion'};
}

=head2 seqXMLversion

 Title   : seqXMLversion
 Usage   : $self->seqXMLversion
 Function: gets/sets the seqXML version in the <seqXML> header
 Returns : the seqXML version string.
 Args    : To set the seqXML version, call with a seqXML version string
           as the argument.

=cut

sub seqXMLversion {
    my ( $self, $value ) = @_;
    my $metadata = $self->{'_seqxml_metadata'};

    # set if a value is supplied
    if ($value) {
        $metadata->{'seqXMLversion'} = $value;
    }

    return $metadata->{'seqXMLversion'};
}

=head1 Methods for parsing the XML document

=cut

=head2 processXMLNode

 Title   : processXMLNode
 Usage   : $seqio->processXMLNode
 Function: reads the XML node and processes according to the node type
 Returns : none
 Args    : none
 Throws  : Exception on unexpected XML node type, warnings on unexpected
           XML element names.

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

    # reset for every new <seqXML> block
    $self->{'_seqxml_metadata'} = {};

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
 Function: processes a sequence <description> node;
           a no-op -- description text is read by
           processXMLnode
 Returns : none
 Args    : none

=cut

sub element_description {
    my ($self) = @_;
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
 Function: processes a sequence <alternativeID> node,
           creating a Bio::Annotation::DBLink object
 Returns : none
 Args    : none

=cut

sub element_alternativeID {
    my ($self) = @_;
    my $reader = $self->{'_reader'};
    my $data   = $self->{'_current_entry_data'};

    my $alternativeID = {};
    my $annotation_obj;

    if ( $reader->hasAttributes() ) {
        $self->processAttribute($alternativeID);
    }
    else {
        $self->throw("no alternative ID data!");
    }

    if (   defined $alternativeID->{'source'}
        && defined $alternativeID->{'id'} )
    {
        $annotation_obj = Bio::Annotation::DBLink->new(
            -primary_id => $alternativeID->{'id'},
            -database   => $alternativeID->{'source'},
            -tagname    => 'dblink',
        );
        push @{ $data->{'alternativeIDs'} }, $annotation_obj;
    }
    else {
        $self->throw("malformed alternative ID data!");
    }
}

=head2 element_property

 Title   : element_property
 Usage   : $self->element_property
 Function: processes a sequence <property> node, creating a
           Bio::Annotation::SimpleValue object
 Returns : none
 Args    : none

=cut

sub element_property {
    my ($self) = @_;
    my $reader = $self->{'_reader'};
    my $data   = $self->{'_current_entry_data'};

    my $property = {};
    my $annotation_obj;

    if ( $reader->hasAttributes() ) {
        $self->processAttribute($property);
    }
    else {
        $self->throw("no property data!");
    }

    if ( defined $property->{'name'} ) {
        $annotation_obj =
          Bio::Annotation::SimpleValue->new( -tagname => $property->{'name'} );

        if ( defined $property->{'value'} ) {
            $annotation_obj->value( $property->{'value'} );
        }

        push @{ $data->{'properties'} }, $annotation_obj;
    }
    else {
        $self->throw("malformated property!");
    }
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
 Returns : a Bio::Seq object
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
    if ( $data->{'alternativeIDs'} ) {
        foreach my $annotation_obj ( @{ $data->{'alternativeIDs'} } ) {
            $seq_obj->add_Annotation($annotation_obj);
        }
    }
    if ( $data->{'properties'} ) {
        foreach my $annotation_obj ( @{ $data->{'properties'} } ) {
            $seq_obj->add_Annotation($annotation_obj);
        }
    }

    # empty the temporary data store
    $self->{'_current_entry_data'} = {};

    return $seq_obj;
}

=head2 end_element_default

 Title   : end_element_default
 Usage   : $self->end_element_default
 Function: processes all other closing tags;
           a no-op.
 Returns : none
 Args    : none

=cut

sub end_element_default {
    my ($self) = @_;
}

=head2 DESTROY

 Title   : DESTROY
 Usage   : called automatically by Perl just before object
           goes out of scope
 Function: writes closing </seqXML> tag and performs a write flush
 Returns : none
 Args    : none

=cut

sub DESTROY {
    my $self = shift;
    if ( $self->mode eq 'w' ) {
        $self->{'_writer'}->endTag("seqXML");
        $self->{'_writer'}->end();
        $self->flush if $self->_flush_on_write && defined $self->_fh;
    }
    $self->SUPER::DESTROY;
}

1;
