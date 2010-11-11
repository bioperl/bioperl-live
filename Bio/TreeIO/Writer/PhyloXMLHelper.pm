# $Id: phyloxml.pm,v 0.01 2007-03-27 12:43:27 heikki Exp $
#
# BioPerl module for Bio::TreeIO::Writer::PhyloXMLHelper
#
# Cared for by Chris Fields cjfields at bioperl dot org
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::TreeIO::Writer::PhyloXMLHelper - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists
  
=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.open-bio.org/
  
=head1 AUTHOR - Chris Fields

  Email cjfields at bioperl dot org

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::TreeIO::Writer::PhyloXMLHelper;

use strict;
use warnings;
use XML::LibXML;
use Scalar::Util qw(blessed reftype);

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : my $obj = new
 Function: Builds a new Bio::TreeIO::Writer::PhyloXMLHelper object 
 Returns : an instance of Bio::TreeIO::Writer::PhyloXMLHelper
 Args    : 
 
=cut

# TODO: get rid of this code stink!

no strict 'refs';

our %TABLE = map {
    $_ => \&{"${_}_el"}
} qw(name branch_length confidence width color taxonomy sequence events
binary_characters distribution date reference property clade symbol accession
name location mol_seq uri annotation domain_architecture id code scientific_name
authority common_name synonym rank uri phyloxml phylogeny document generic);

use strict;

# Just an external generic bag of methods for generating XML node elements,
# ayttributes, etc for phyloXML. Hopefully some usable parts that could be used
# elsewhere

# TODO: maybe wrap this up in an event-based writer

sub create_node {
    my ($self, $data, $parent) = @_;
    if ($parent) {
        my $class = blessed($parent);
        if (!defined($class) || !$parent->isa('XML::LibXML::Node')) {
            $self->throw('The parent must be a XML::LibXML::Node, got '.ref($parent));
        } 
    }
    $self->throw("Must supply a defined data 'Name' key") unless defined($data->{Name});
    if (exists $TABLE{$data->{Name}}) {
        return $TABLE{$data->{Name}}->($data, $parent);
    } else {
        $self->debug("Unknown event: ".$data->{Name}."\n");
    }
}

####################### OUTPUT HANDLERS #######################

# These are split up into very specific methods for each piece of data; we'll
# optimize from there (simplify down)

# basically these just create the XML::LibXML::Note/Attribute as needed

# All methods get a string and an optional XML::LibXML-compliant parent node,
# trying to be very non-BioPerl-reliant

# TODO: determine if $str should be a simple string or hold more data (hashref)

sub document_el {
    my ($data) = @_; # no parent
    XML::LibXML::Document->new( '1.0', 'UTF-8' );
}

# catch-all
sub generic_el {
    my ($data, $parent) = @_;
    my ($tag, $val, $attributes, $type) =
        @{$data->{Data}}{qw(Tag Value Attributes Class)};
    
    # not using class for now

    my $node = XML::LibXML::Element->new($tag);
    
    # should the attributes be assigned to the parent or the new node? Latter...
    if ($attributes) {
        while (my ($name, $att_value) = each(%$attributes)) {
            $node->setAttribute( $name, $att_value);
        }
    }
    
    $node->appendTextNode($val) if $val;
    if ($parent) {
        $parent->addChild($node);
    }
    $node;
    #my $node = $parent ? $parent->addNewChild( undef, $str) : XML::LibXML::Element->new($str);
}

sub phyloxml_el {
    my ($data, $parent) = @_;
    my $root_el = $parent->createElement('phyloxml');
    $root_el->setNamespace("http://www.w3.org/2001/XMLSchema-instance", 'xsi', 0);
    $root_el->setNamespace("http://www.phyloxml.org");
    $parent->setDocumentElement($root_el);
    $root_el;
}

sub phylogeny_el {
    my ($data, $parent) = @_;
    Bio::Root::Root->throw("Parent must be a XML::LibXML::Node") if !$parent
        || !$parent->isa('XML::LibXML::Node') ;
    $parent->addNewChild(undef, $data->{Name});
}

sub branch_length_el {}
sub confidence_el {}
sub width_el {}
sub color_el {}
sub taxonomy_el {}
sub sequence_el {}
sub events_el {}
sub binary_characters_el {}
sub distribution_el {}
sub date_el {}
sub reference_el {}
sub property_el {}

# For <sequence>, the order is:sub _symbol_el {}
sub accession_el {}
sub name_el {}
sub location_el {}
sub mol_seq_el {}
sub uri_el {}
sub annotation_el {}
sub domain_architecture_el {}

   
# For <taxonomy>, the order is:szub _id_el {}
sub code_el {}
sub scientific_name_el {}
sub authority_el {}
sub common_name_el {}
sub synonym_el {}
sub rank_el {}

#sub _uri_el {}



1;
