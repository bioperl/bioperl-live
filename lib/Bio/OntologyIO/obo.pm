#
# BioPerl module for Bio::OntologyIO::obo
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Sohel Merchant, s-merchant at northwestern.edu
#
# Copyright Sohel Merchant
#
# You may distribute this module under the same terms as perl itself


=head1 NAME

Bio::OntologyIO::obo - parser for OBO flat-file format

=head1 SYNOPSIS

  use Bio::OntologyIO;

  # do not use directly -- use via Bio::OntologyIO
  my $parser = Bio::OntologyIO->new
        ( -format => "obo",
          -file   =>  "gene_ontology.obo");

  while(my $ont = $parser->next_ontology()) {
  print "read ontology ",$ont->name()," with ",
               scalar($ont->get_root_terms), " root terms, and ",
               scalar($ont->get_all_terms),  " total terms, and ",
               scalar($ont->get_leaf_terms), " leaf terms\n";
  }

=head1 DESCRIPTION

Parser for OBO flat-file format. 'obo' example:

 format-version: 1.2
 ontology: so/dev/externalDerived
 property_value: owl:versionInfo "$Revision: 80 $" xsd:string
 default-namespace: SO

 [Term]
 id: SO_0000343
 name: match
 def: "A region of sequence, aligned to another sequence." []

 [Term]
 id: SO_0000039
 name: match_part
 def: "A part of a match." []
 is_a: SO_0000343

Specification: L<http://www.geneontology.org/GO.format.obo-1_2.shtml>.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

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

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Sohel Merchant

Email: s-merchant@northwestern.edu

Address:

  Northwestern University
  Center for Genetic Medicine (CGM), dictyBase
  Suite 1206,
  676 St. Clair st
  Chicago IL 60611

=head2 CONTRIBUTOR

 Hilmar Lapp, hlapp at gmx.net
 Chris Mungall, cjm at fruitfly.org
 Brian Osborne, briano@bioteam.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package  Bio::OntologyIO::obo;

use strict;

use Bio::Root::IO;
use Bio::Ontology::OBOEngine;
use Bio::Ontology::Ontology;
use Bio::Ontology::OntologyStore;
use Bio::Ontology::TermFactory;
use Bio::Annotation::Collection;
use Text::Balanced qw(extract_quotelike extract_bracketed);

use constant TRUE  => 1;
use constant FALSE => 0;

use base qw(Bio::OntologyIO);

=head2 new

 Title   : new
 Usage   : $parser = Bio::OntologyIO->new(
                             -format => "obo",
                             -file => "gene_ontology.obo");
 Function: Creates a new dagflat parser.
 Returns : A new dagflat parser object, implementing Bio::OntologyIO.
 Args    : -file      => a single ontology flat file holding the
                         terms, descriptions and relationships
           -ontology_name => the name of the ontology; if not specified the
                          parser will assign the name of the ontology as the
                          default-namespace header value from the OBO file.
           -engine     => the Bio::Ontology::OntologyEngineI object
                          to be reused (will be created otherwise); note
                          that every Bio::Ontology::OntologyI will
                          qualify as well since that one inherits from the
                          former.

See L<Bio::OntologyIO>.

=cut

# Let OntologyIO::new() do the instantiation, and override
# _initialize for all initialization work
sub _initialize {
    my ( $self, %arg ) = @_;

    my ( $file, $name, $eng ) = $self->_rearrange(
        [ qw( FILE ONTOLOGY_NAME ENGINE) ], %arg
    );

    $self->SUPER::_initialize(%arg);
    delete $self->{'_ontologies'};

    # Ontology engine (and possibly name if it's an OntologyI)
    $eng = Bio::Ontology::OBOEngine->new() unless $eng;
    if ( $eng->isa("Bio::Ontology::OntologyI") ) {
        $self->ontology_name( $eng->name() );
        $eng = $eng->engine() if $eng->can('engine');
    }
    $self->_ont_engine($eng);

    $self->ontology_name($name) if $name;
}

=head2 ontology_name

 Title   : ontology_name
 Usage   : $obj->ontology_name($newval)
 Function: Get/set the name of the ontology parsed by this module.
 Example :
 Returns : value of ontology_name (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub ontology_name {
    my $self = shift;

    return $self->{'ontology_name'} = shift if @_;
    return $self->{'ontology_name'};
}

=head2 parse

 Title   : parse()
 Usage   : $parser->parse();
 Function: Parses the files set with "new" or with methods
           defs_file and _flat_files.

           Normally you should not need to call this method as it will
           be called automatically upon the first call to
           next_ontology().

 Returns : Bio::Ontology::OntologyEngineI
 Args    :

=cut

sub parse {
    my $self = shift;

    # Setup the default term factory if not done by anyone yet
    $self->term_factory(Bio::Ontology::TermFactory->new( -type => "Bio::Ontology::OBOterm" ) )
      unless $self->term_factory();

    # Parse the file header
    my $annotations_collection = $self->_header();

    # Create the default ontology object itself
    my $ont = Bio::Ontology::Ontology->new(
        -name   => $self->ontology_name(),
        -engine => $self->_ont_engine()
    );

    # Assign the file headers
    $ont->annotation($annotations_collection);

    # Set up the ontology of the relationship types
    for (
        $self->_part_of_relationship(),
        $self->_is_a_relationship(),
        $self->_related_to_relationship(),
        $self->_regulates_relationship(),
        $self->_positively_regulates_relationship(),
        $self->_negatively_regulates_relationship(),
      )
    {
        $_->ontology($ont);
    }

    $self->_add_ontology($ont);

    # Adding new terms
    while ( my $term = $self->_next_term() ) {

        # Check if the terms has a valid ID and NAME otherwise ignore the term
        if ( !$term->identifier() || !$term->name() ) {
            $self->throw( "OBO File Format Error on line "
                  . $self->{'_current_line_no'}
                  . "\nThe term does not have a id/name tag. This term will be ignored."
            );
            next;
        }

        my $new_ontology_flag = 1;
        my $ontologies_array_ref = $self->{'_ontologies'};

        for my $ontology ( @$ontologies_array_ref ) {
            my ($oname, $t_ns) = ( $ontology->name, $term->namespace );
            next unless ( defined($oname) && defined($t_ns) );
            if ( $oname eq $t_ns ) {
                # No need to create new ontology
                $new_ontology_flag = 0;
                $ont = $ontology;
            }
        }

        if ( $new_ontology_flag && $term->namespace ) {
            my $new_ont = Bio::Ontology::Ontology->new(
                -name   => $term->namespace,
                -engine => $self->_ont_engine
            );
            $new_ont->annotation($annotations_collection);
            $self->_add_ontology($new_ont);
            $ont = $new_ont;
        }

        $self->_add_term( $term, $ont );

        # Adding the IS_A relationship
        for my $parent_term ( @{$self->{'_isa_parents'}} ) {
            # Check if parent exists, if not then add the term to the graph.
            if ( ! $self->_has_term($parent_term) ) {
                $self->_add_term( $parent_term, $ont ); # !
            }

            $self->_add_relationship( $parent_term, $term,
                $self->_is_a_relationship(), $ont );
        }

        # Adding the other relationships like part_of, related_to, develops_from
        for my $relationship ( keys %{$self->{'_relationships'}} ) {
            my $reltype;
            # Check if relationship exists, if not add it
            if ( $self->_ont_engine->get_relationship_type($relationship) ) {
                $reltype = $self->_ont_engine->get_relationship_type($relationship);
            }
            else {
                $self->_ont_engine->add_relationship_type( $relationship, $ont );
                $reltype = $self->_ont_engine->get_relationship_type($relationship);
            }

            # Check if the id already exists in the graph
            for my $id ( @{$self->{'_relationships'}->{$relationship}} ) {
                my $parent_term = $self->_create_term_object();
                $parent_term->identifier($id);
                $parent_term->ontology($ont);

                if ( ! $self->_has_term($parent_term) ) {
                    $self->_add_term( $parent_term, $ont );
                }

                $self->_add_relationship( $parent_term, $term, $reltype, $ont );
            }
        }
    }

    return $self->_ont_engine();
}

=head2 next_ontology

 Title   : next_ontology
 Usage   :
 Function: Get the next available ontology from the parser. This is the
           method prescribed by Bio::OntologyIO.
 Example :
 Returns : An object implementing Bio::Ontology::OntologyI, and nothing if
           there is no more ontology in the input.
 Args    :

=cut

sub next_ontology {
    my $self = shift;

    # Parse if not done already
    $self->parse() unless exists( $self->{'_ontologies'} );

    # Return next available ontology
    if ( exists( $self->{'_ontologies'} ) ) {
        my $ont = shift( @{ $self->{'_ontologies'} } );
        if ($ont) {
            my $store = Bio::Ontology::OntologyStore->new();
            $store->register_ontology($ont);

            return $ont;
        }
    }
    return;
}

=head2 close

 Title   : close
 Usage   :
 Function: Closes this ontology stream and associated file handles.

           Clients should call this method especially when they write
           ontologies.

           We need to override this here in order to close the file
           handle for the term definitions file.

 Example :
 Returns : none
 Args    : none

=cut

sub close {
    my $self = shift;

    # first call the inherited implementation
    $self->SUPER::close();
}

# INTERNAL METHODS

sub _add_ontology {
    my $self = shift;
    $self->{'_ontologies'} = [] unless exists( $self->{'_ontologies'} );
    for my $ont (@_) {
        $self->throw( ref($ont) . " does not implement Bio::Ontology::OntologyI" )
            unless ref($ont) && $ont->isa("Bio::Ontology::OntologyI");
        # The ontology name may have been auto-discovered while parsing
        # the file
        $ont->name( $self->ontology_name ) unless $ont->name();
        push( @{ $self->{'_ontologies'} }, $ont );
    }
}

# This simply delegates. See Ontology::OBOEngine::add_term.
sub _add_term {
    my ( $self, $term, $ont ) = @_;
    $term->ontology($ont) if $ont && ( !$term->ontology );
    $self->_ont_engine()->add_term($term);
}

# This simply delegates. See OBOEngine
sub _part_of_relationship {
    my $self = shift;

    return $self->_ont_engine()->part_of_relationship(@_);
}

# This simply delegates. See OBOEngine
sub _is_a_relationship {
    my $self = shift;

    return $self->_ont_engine()->is_a_relationship(@_);
}

# This simply delegates. See OBOEngine
sub _related_to_relationship {
    my $self = shift;

    return $self->_ont_engine()->related_to_relationship(@_);
}

# This simply delegates. See OBOEngine
sub _regulates_relationship {
    my $self = shift;

    return $self->_ont_engine()->regulates_relationship(@_);
}

# This simply delegates. See OBOEngine
sub _positively_regulates_relationship {
    my $self = shift;

    return $self->_ont_engine()->positively_regulates_relationship(@_);
}

# This simply delegates. See OBOEngine
sub _negatively_regulates_relationship {
    my $self = shift;

    return $self->_ont_engine()->negatively_regulates_relationship(@_);
}

# This simply delegates. See OBOEngine
sub _add_relationship {
    my ( $self, $parent, $child, $type, $ont ) = @_;
    # Note the triple terminology (subject,predicate,object) corresponds to
    # (child,type,parent)
    $self->_ont_engine()->add_relationship( $child, $type, $parent, $ont );
}

# This simply delegates. See OBOEngine
sub _has_term {
    my $self = shift;

    return $self->_ont_engine()->has_term(@_);
}

# Holds the OBO engine to be parsed into
sub _ont_engine {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{"_ont_engine"} = $value;
    }

    $self->{"_ont_engine"};
}

# Removes the escape characters from the file
sub _filter_line {
    my ( $self, $line ) = @_;

    chomp($line);
    $line =~ tr [\200-\377] [\000-\177];
                               # see 'man perlop', section on tr/
                               # weird ascii characters should be excluded
    $line =~ tr/\0-\10//d;     # remove weird characters; ascii 0-8
                               # preserve \11 (9 - tab) and \12 (10-linefeed)
    $line =~ tr/\13\14//d;     # remove weird characters; 11,12
                               # preserve \15 (13 - carriage return)
    $line =~ tr/\16-\37//d;    # remove 14-31 (all rest before space)
    $line =~ tr/\177//d;       # remove DEL character

    $line =~ s/^\!.*//;
    $line =~ s/[^\\]\!.*//;
    $line =~ s/[^\\]\#.*//;
    $line =~ s/^\s+//;
    $line =~ s/\s+$//;

    return $line;
}

# Parses the header
sub _header {
    my $self = shift;
    my $annotation_collection = Bio::Annotation::Collection->new();
    my ( $tag, $value );
    my $line_counter = 0;
    $self->{'_current_line_no'} = 0;
    my $format_version_header_flag = 0;
    my $default_namespace_header_flag = 0;

    while ( my $line = $self->_readline() ) {
        ++$line_counter;
        my $line = $self->_filter_line($line);

        if ( !$line ) {
            if ( !$format_version_header_flag ) {
                $self->throw("Format Error - Cannot find tag format-version." .
                             "This is required in header" );
            }

            $self->{'_current_line_no'} = $line_counter;
            return $annotation_collection;
        }

        # Check if there is a header
        if ( $line =~ /\[\w*\]/ ) {
                $self->throw("Format Error - Cannot find tag format-version." .
                             "This is required in header." );
          }

        # If the line is not null, check it contains at least one colon
        $self->_check_colon( $line, $line_counter );

        # These are the allowed headers. Any other headers will be ignored
        if ( $line =~ /^(\[|format-version:
                           |data-version:
                           |typeref:
                           |version:
                           |date:
                           |saved-by:
                           |auto-generated-by:
                           |default-namespace:
                           |remark:
                           |subsetdef:
                           |import:
                           |synonymtypedef:
                           |idspace:
                           |default-relationship-id-prefix:
                           |id-mapping:
                           )/x
            )
        {
            if ( $line =~ /^([\w\-]+)\:\s*(.*)/ ) {
                ( $tag, $value ) = ( $1, $2 );
            }

            if ( $tag =~ /format-version/) {
                $format_version_header_flag = 1;
            }elsif( $tag =~ /default-namespace/ ) {
                $default_namespace_header_flag = 1;
            }

            my $header = Bio::Annotation::SimpleValue->new( -value => $value );
            $annotation_collection->add_Annotation( $tag, $header );

            # Assign the Ontology name as the value of the default-namespace header
            if ( $tag =~ /default-namespace/i ) {
                $self->ontology_name($value);
            }
        }
    }
}

# Parses each stanza of the file
sub _next_term {
    my $self = shift;
    my $term;
    my $skip_stanza_flag = 1;
    my $line_counter     = $self->{'_current_line_no'};

    while ( my $line = $self->_readline() ) {
        ++$line_counter;
        my $line = $self->_filter_line($line);

        if ( !$line && $term ) {
            $self->{'_current_line_no'} = $line_counter;
            return $term;
        }

        if ( ( $line =~ /^\[(\w+)\]\s*(.*)/ ) ) { # New stanza
            if ( uc($1) eq "TERM" ) {
                $term             = $self->_create_term_object;
                $skip_stanza_flag = 0;

                # Reset the relationships after each stanza
                $self->{'_relationships'} = {};
                $self->{'_isa_parents'}   = undef;
            }
            elsif ( uc($1) eq "TYPEDEF" ) {
                $skip_stanza_flag = 1;
                # Check if this typedef is already defined by the relationship
            }
            else {
                $skip_stanza_flag = 1;
                $self->warn("OBO File Format Warning on line $line_counter $line\n"
                          . "Unrecognized stanza type found. Skipping this stanza." );
            }
            next;
        }

        # If the line is not null, check it contains at least one colon
        $self->_check_colon( $line, $line_counter );

        # If there is any tag value other than the list below move to the next tag
        next if (( $line !~ /^(\[|id:
                                 |is_anonymous:
                                 |name:
                                 |namespace:
                                 |alt_id:
                                 |def:
                                 |comment:
                                 |subset:
                                 |synonym:
                                 |xref:
                                 |is_a:
                                 |intersection_of:
                                 |union_of:
                                 |disjoint_from:
                                 |relationship:
                                 |is_obsolete:
                                 |replaced_by:
                                 |consider:
                                 |created_by:
                                 |creation_date:
                               )/x
                ) || $skip_stanza_flag );

        # Tag/value pair
        if ( $line =~ /^([\w\-]+)\:\s*(.*)/ ) {
            my ( $tag, $val ) = ( $1, $2 );

            # If no value for the tag throw a warning
            if ( !$val ) {
                $self->warn("OBO File Format Warning on line $line_counter $line\n" .
                            "Tag has no value."
                );
            }

            my $qh;
            ( $val, $qh ) = $self->_extract_quals($val);
            my $val2 = $val;
            $val2 =~ s/\\,/,/g;
            $tag = uc($tag);
            if ( $tag eq "ID" ) {
                $term->identifier($val);
                if ( $self->_has_term($term) ) {
                    $term = $self->_ont_engine()->get_terms($val);
                }
            }
            elsif ( $tag eq "NAME" ) {
                $term->name($val);
            }
            elsif ( $tag eq "XREF_ANALOG" ) {
                if ( !$term->has_dbxref($val) ) {
                    $term->add_dbxref(-dbxrefs => $self->_to_annotation( [$val] ) );
                }
            }
            elsif ( $tag eq "XREF_UNKNOWN" ) {
                $term->add_dbxref(-dbxrefs => $self->_to_annotation( [$val] ) );
            }
            elsif ( $tag eq "NAMESPACE" ) {
                $term->namespace($val);
            }
            elsif ( $tag eq "DEF" ) {
                my ( $defstr, $parts ) = $self->_extract_qstr($val);
                $term->definition($defstr);
                my $ann = $self->_to_annotation($parts);
                $term->add_dbxref( -dbxrefs => $ann );
            }
            elsif ( $tag eq "SYNONYM" ) {
                $term->add_synonym($val);
            }
            elsif ( $tag eq "ALT_ID" ) {
                $term->add_secondary_id($val);
            }
            elsif ( $tag =~ /XREF/i ) {
                $term->add_secondary_id($val);
            }
            elsif ( $tag eq "IS_OBSOLETE" ) {
                if ( $val eq 'true' ) {
                    $val = 1;
                }
                elsif ( $val eq 'false' ) {
                    $val = 0;
                }
                $term->is_obsolete($val);
            }
            elsif ( $tag eq "COMMENT" ) {
                $term->comment($val);
            }
            elsif ( $tag eq "RELATIONSHIP" ) {
                $self->_handle_relationship_tag($val);
            }
            elsif ( $tag eq "IS_A" ) {
                $val =~ s/ //g;
                my $parent_term = $self->_create_term_object();
                $parent_term->identifier($val);
                push @{ $self->{'_isa_parents'} }, $parent_term;
            }
        }
    }

    $term;
}


# Creates a Bio::Ontology::OBOterm object
sub _create_term_object {
    my ($self) = @_;
    my $term = $self->term_factory->create_object();
    $term;
}

sub _extract_quals {
    my ( $self, $str ) = @_;

    my %q = ();
    if ( $str =~ /(.*)\s+(\{.*\})\s*$/ ) {
        my $return_str = $1;
        my $extr       = $2;
        if ($extr) {
            my @qparts = $self->_split_on_comma($extr);
            foreach (@qparts) {
                if (/(\w+)=\"(.*)\"/) {
                    $q{$1} = $2;
                }
                elsif (/(\w+)=\'(.*)\'/) {
                    $q{$1} = $2;
                }
                else {
                    warn("$_ in $str");
                }
            }
        }
        return ( $return_str, \%q );
    }
    else {
        return ( $str, {} );
    }
}

sub _extract_qstr {
    my ( $self, $str ) = @_;

    my ( $extr, $rem, $prefix ) = extract_quotelike($str);
    my $txt = $extr;
    $txt =~ s/^\"//;
    $txt =~ s/\"$//;
    if ($prefix) {
        warn("illegal prefix: $prefix in: $str");
    }

    my @extra = ();

    # e.g. synonym: "foo" EXACT [...]
    if ( $rem =~ /(\w+)\s+(\[.*)/ ) {
        $rem = $2;
        push( @extra, split( ' ', $1 ) );
    }

    my @parts = ();
    while ( ( $extr, $rem, $prefix ) = extract_bracketed( $rem, '[]' ) ) {
        last unless $extr;
        $extr =~ s/^\[//;
        $extr =~ s/\]$//;
        push( @parts, $extr ) if $extr;
    }
    @parts =
      map { $self->_split_on_comma($_) } @parts;

    $txt =~ s/\\//g;

    ( $txt, \@parts, \@extra );
}

sub _split_on_comma {
    my ( $self, $str ) = @_;
    my @parts = ();
    while ( $str =~ /(.*[^\\],\s*)(.*)/ ) {
        $str = $1;
        my $part = $2;
        unshift( @parts, $part );
        $str =~ s/,\s*$//;
    }
    unshift( @parts, $str );

    return map { s/\\//g; $_ } @parts;
}

# This method checks for an existing colon in a line
sub _check_colon {
    my ( $self, $line, $line_no ) = @_;
    if ( $line && !( $line =~ /:/ ) ) {
        $self->throw("OBO File Format Error on line $line_no $line\n" .
                     "Cannot find key-terminating colon"
        );
    }
}

# This method handles relationship tags
sub _handle_relationship_tag {
    my ( $self, $val ) = @_;
    my @parts        = split( / /, $val );
    my $relationship = uc($parts[0]);
    my $id           = $parts[1] =~ /\^(w+)\s+\!/ ? $1 : $parts[1];
    push @{$self->{_relationships}->{$relationship}}, $id;
}

# Convert simple strings to Bio::Annotation::DBLinks
sub _to_annotation {
    my ($self , $links) = @_;
    return unless $links;
    my @dbxrefs;
    for my $string (@{$links}) {
        my ($db, $id) = split(':',$string,2);
        push @dbxrefs, Bio::Annotation::DBLink->new(-database => $db, -primary_id => $id);
    }

    \@dbxrefs;
}

1;
