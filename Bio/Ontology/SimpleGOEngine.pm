# $Id$
#
# BioPerl module for Bio::Ontology::SimpleGOEngine
#
# Cared for by Christian M. Zmasek <czmasek@gnf.org> or <cmzmasek@yahoo.com>
#
# (c) Christian M. Zmasek, czmasek@gnf.org, 2002.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2002.
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
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Ontology::SimpleGOEngine - a Ontology Engine for GO implementing OntologyEngineI

=head1 SYNOPSIS

  use Bio::Ontology::SimpleGOEngine;

  my $parser = Bio::Ontology::SimpleGOEngine->new
        ( -defs_file => "/home/czmasek/GO/GO.defs",
          -files     => ["/home/czmasek/GO/component.ontology",
                         "/home/czmasek/GO/function.ontology",
                         "/home/czmasek/GO/process.ontology"] );

  my $engine = $parser->parse();

  my $IS_A       = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
  my $PART_OF    = Bio::Ontology::RelationshipType->get_instance( "PART_OF" );
  my $RELATED_TO = Bio::Ontology::RelationshipType->get_instance( "RELATED_TO" );

=head1 DESCRIPTION

This class is deprecated and instead Bio::Ontology::OBOEngine should be used.

Needs Graph.pm from CPAN.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via
the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR

Christian M. Zmasek

Email: czmasek@gnf.org  or  cmzmasek@yahoo.com

WWW:   http://www.genetics.wustl.edu/eddy/people/zmasek/

Address:

  Genomics Institute of the Novartis Research Foundation
  10675 John Jay Hopkins Drive
  San Diego, CA 92121

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...
package Bio::Ontology::SimpleGOEngine;

use strict;

use base qw(Bio::Ontology::OBOEngine);


# Internal methods
# ----------------

## Overiding this method from OBOEngine
# Checks the correct format of a GOBO-formatted id
# Gets the id out of a term or id string
sub _get_id {
    my ( $self, $term ) = @_;
    my $id = $term;

    if(ref($term)) {
        # use TermI standard API
        $self->throw("Object doesn't implement Bio::Ontology::TermI. ".
                     "Bummer.")
            unless $term->isa("Bio::Ontology::TermI");
        $id = $term->identifier();
        # if there is no ID, we need to fake one from ontology name and name
        # in order to achieve uniqueness
        if(!$id) {
            $id = $term->ontology->name() if $term->ontology();
            $id = $id ? $id.'|' : '';
            $id .= $term->name();
        }
    }
    # don't fuss if it looks remotely standard, and we trust GO terms
    return $id
#        if $term->isa("Bio::Ontology::GOterm")||($id =~ /^[A-Z_]{1,8}:\d{1,}$/);
        if $term->isa("Bio::Ontology::GOterm")||($id =~ /^\w+:\w+$/);
    # prefix with something if only numbers
    if($id =~ /^\d+$/) {
        $self->warn(ref($self).": identifier [$id] is only numbers - ".
                    "prefixing with 'GO:'");
        return "GO:" . $id;
    }
    # we shouldn't have gotten here if it's at least a remotely decent ID
    $self->throw(ref($self).": non-standard identifier '$id'\n")
        unless $id =~ /\|/;
    return $id;
} # _get_id


1;
