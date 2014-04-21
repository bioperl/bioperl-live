#
# BioPerl module for Bio::OntologyIO
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp at gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

#
# (c) Hilmar Lapp, hlapp at gmx.net, 2003.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2003.
#
# You may distribute this module under the same terms as perl itself.
# Refer to the Perl Artistic License (see the license accompanying this
# software package, or see http://www.perl.com/language/misc/Artistic.html)
# for the terms under which you may use, modify, and redistribute this module.
#
# THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#

# POD documentation - main docs before the code

=head1 NAME

Bio::OntologyIO - Parser factory for Ontology formats

=head1 SYNOPSIS

    use Bio::OntologyIO;

    my $parser = Bio::OntologyIO->new(-format => "go",
                                      -file=> $file);

    while(my $ont = $parser->next_ontology()) {
         print "read ontology ",$ont->name()," with ",
               scalar($ont->get_root_terms)," root terms, and ",
               scalar($ont->get_leaf_terms)," leaf terms\n";
    }

=head1 DESCRIPTION

This is the parser factory for different ontology sources and
formats. Conceptually, it is very similar to L<Bio::SeqIO>, but the
difference is that the chunk of data returned as an object is an
entire ontology.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::OntologyIO;
use strict;

# Object preamble - inherits from Bio::Root::Root


use base qw(Bio::Root::Root Bio::Root::IO);

#
# Maps from format name to driver suitable for the format.
#
my %format_driver_map = (
                         "go"          => "goflat",
                         "so"          => "soflat",
                         "interpro"    => "InterProParser",
                         "interprosax" => "Handlers::InterPro_BioSQL_Handler",
                         "evoc"        => "simplehierarchy",
                         "obo"         => "obo"
                         );

=head2 new

 Title   : new
 Usage   : my $parser = Bio::OntologyIO->new(-format => 'go', @args);
 Function: Returns a stream of ontologies opened on the specified input
           for the specified format.
 Returns : An ontology parser (an instance of Bio::OntologyIO) initialized
           for the specified format.
 Args    : Named parameters. Common parameters are

              -format    - the format of the input; the following are
                           presently supported:
                  goflat: DAG-Edit Gene Ontology flat files
                  go    : synonymous to goflat
                  soflat: DAG-Edit Sequence Ontology flat files
                  so    : synonymous to soflat
                  simplehierarchy: text format with one term per line
                          and indentation giving the hierarchy
                  evoc  : synonymous to simplehierarchy
                  interpro: InterPro XML
                  interprosax: InterPro XML - this is actually not a
                          Bio::OntologyIO compliant parser; instead it
                          persists terms as they are encountered.
                          L<Bio::OntologyIO::Handlers::InterPro_BioSQL_Handler>
                  obo   : OBO format style from Gene Ontology Consortium
              -file      - the file holding the data
              -fh        - the stream providing the data (-file and -fh are
                          mutually exclusive)
              -ontology_name - the name of the ontology
              -engine    - the L<Bio::Ontology::OntologyEngineI> object
                          to be reused (will be created otherwise); note
                          that every L<Bio::Ontology::OntologyI> will
                          qualify as well since that one inherits from the
                          former.
              -term_factory - the ontology term factory to use. Provide a
                          value only if you know what you are doing.

           DAG-Edit flat file parsers will usually also accept the
           following parameters.

              -defs_file - the name of the file holding the term
                          definitions
              -files     - an array ref holding the file names (for GO,
                          there will usually be 3 files: component.ontology,
                          function.ontology, process.ontology)

           Other parameters are specific to the parsers.

=cut

sub new {
    my ($caller,@args) = @_;
    my $class = ref($caller) || $caller;
    # or do we want to call SUPER on an object if $caller is an
    # object?
    if( $class =~ /Bio::OntologyIO::(\S+)/ ) {
        my ($self) = $class->SUPER::new(@args);
        $self->_initialize(@args);
        return $self;
    } else {
        my %param = @args;
        @param{ map { lc $_ } keys %param } = values %param; # lowercase keys
        my $format = $class->_map_format($param{'-format'});

        # normalize capitalization
        return unless( $class->_load_format_module($format) );
        return "Bio::OntologyIO::$format"->new(@args);
    }

}


=head2 format

 Title   : format
 Usage   : $format = $parser->format()
 Function: Get the ontology format
 Returns : ontology format
 Args    : none

=cut

# format() method inherited from Bio::Root::IO


sub _initialize {
    my($self, @args) = @_;

    # initialize factories etc
    my ($eng,$fact,$ontname) =
        $self->_rearrange([qw(TERM_FACTORY)
                           ], @args);
    # term object factory
    $self->term_factory($fact) if $fact;

    # initialize the Bio::Root::IO part
    $self->_initialize_io(@args);
}

=head2 next_ontology

 Title   : next_ontology
 Usage   : $ont = $stream->next_ontology()
 Function: Reads the next ontology object from the stream and returns it.
 Returns : a L<Bio::Ontology::OntologyI> compliant object, or undef at the
           end of the stream
 Args    : none


=cut

sub next_ontology {
    shift->throw_not_implemented();
}

=head2 term_factory

 Title   : term_factory
 Usage   : $obj->term_factory($newval)
 Function: Get/set the ontology term factory to use.

           As a user of this module it is not necessary to call this
           method as there will be default. In order to change the
           default, the easiest way is to instantiate
           L<Bio::Ontology::TermFactory> with the proper -type
           argument. Most if not all parsers will actually use this
           very implementation, so even easier than the aforementioned
           way is to simply call
           $ontio->term_factory->type("Bio::Ontology::MyTerm").

 Example :
 Returns : value of term_factory (a Bio::Factory::ObjectFactoryI object)
 Args    : on set, new value (a Bio::Factory::ObjectFactoryI object, optional)


=cut

sub term_factory{
    my $self = shift;

    return $self->{'term_factory'} = shift if @_;
    return $self->{'term_factory'};
}

=head1 Private Methods

  Some of these are actually 'protected' in OO speak, which means you
  may or will want to utilize them in a derived ontology parser, but
  you should not call them from outside.

=cut

=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL OntologyIO stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns :
 Args    :

=cut

sub _load_format_module {
    my ($self, $format) = @_;
    my $module = "Bio::OntologyIO::" . $format;
    my $ok;

    eval {
        $ok = $self->_load_module($module);
    };
    if ( $@ ) {
        print STDERR <<END;
$self: $format cannot be found
Exception $@
For more information about the OntologyIO system please see the docs.
This includes ways of checking for formats at compile time, not run time
END
    }
    return $ok;
}

sub DESTROY {
    my $self = shift;

    $self->close();
}

sub _map_format {
    my $self = shift;
    my $format = shift;
    my $mod;

    if($format) {
        $mod = $format_driver_map{lc($format)};
        $mod = lc($format) unless $mod;
    } else {
        $self->throw("unable to guess ontology format, specify -format");
    }
    return $mod;
}

sub unescape {
  my( $self, $ref ) = @_;
  $ref =~ s/&lt\\;/\</g;
  $ref =~ s/&gt\\;/\>/g;
  $ref =~ s/&pct\\;/\%/g;
  $ref =~ s/\\n/\n/g;
  $ref =~ s/\\t/\t/g;
  return $ref;
}

1;
