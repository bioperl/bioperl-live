# $Id$
#
# BioPerl module for Bio::OntologyIO
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

    my $parser = Bio::OntologyIO->new(-format => "go", ...);

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

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::OntologyIO;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Root::IO;

@ISA = qw(Bio::Root::Root Bio::Root::IO);

#
# Maps from format name to driver suitable for the format.
#
my %format_driver_map = (
			 "go"       => "simpleGOparser",
			 "interpro" => "InterProParser",
			 );

=head2 new

 Title   : new
 Usage   : my $parser = Bio::OntologyIO->new(-format => 'go', @args);
 Function: Returns a stream of ontologies opened on the specified input
           for the specified format.
 Returns : An ontology parser (an instance of Bio::OntologyIO) initialized
           for the specified format.
 Args    : Named parameters. Common parameters are

              -format  the format of the input; supported right now are
                       'go' and 'interpro'
              -file    the file holding the data
              -fh      the stream providing the data (-file and -fh are
                       mutually exclusive)
           
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
	return undef unless( $class->_load_format_module($format) );
	return "Bio::OntologyIO::$format"->new(@args);
    }
}

sub _initialize {
    my($self, @args) = @_;
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

1;
