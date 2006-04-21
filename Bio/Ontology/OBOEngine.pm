
=head1 NAME

OBOEngine - An Ontology Engine for OBO style flat file format from Gene Ontology Consortium which inherits from Bio::Ontology::SimpleGOEngine

=head1 SYNOPSIS

  use Bio::Ontology::OBOEngine;

  my $parser = Bio::Ontology::OBOEngine->new
        ( -file     => "gene_ontology.obo" );

  my $engine = $parser->parse();

=head1 DESCRIPTION

Needs Graph.pm from CPAN.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                         - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via
the web:

  http://bugzilla.bioperl.org/

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
 Chris Mungall,   cjm at fruitfly.org


=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Ontology::OBOEngine;

use strict;
use vars qw( @ISA );
use Bio::Ontology::SimpleGOEngine;

@ISA = qw( Bio::Ontology::SimpleGOEngine );

=head2 new

 Title   : new
 Usage   : $engine = Bio::Ontology::OBOEngine->new()
 Function: Creates a new OBOEngine
 Returns : A new OBOEngine object
 Args    :

=cut

sub new {
    my ( $class, @args ) = @_;

    my $self = $class->SUPER::new(@args);

    $self->init();

    return $self;
}    # new

# Internal methods
# ----------------

## Overirding this method from SimpleGoEngine
# Checks the correct format of a GOBO-formatted id
# Gets the id out of a term or id string
sub _get_id {
    my ( $self, $term ) = @_;
    my $id = $term;

    if ( ref($term) ) {

        # use TermI standard API
        $self->throw(
            "Object doesn't implement Bio::Ontology::TermI. " . "Bummer." )
          unless $term->isa("Bio::Ontology::TermI");
        $id = $term->identifier();

        # if there is no ID, we need to fake one from ontology name and name
        # in order to achieve uniqueness
        if ( !$id ) {
            $id = $term->ontology->name() if $term->ontology();
            $id = $id ? $id . '|' : '';
            $id .= $term->name();
        }
    }

    return $id

#        if $term->isa("Bio::Ontology::GOterm")||($id =~ /^[A-Z_]{1,8}:\d{1,}$/);
      if $term->isa("Bio::Ontology::OBOterm") || ( $id =~ /^\w+:\w+$/ );

    # prefix with something if only numbers
    #     if($id =~ /^\d+$/) {
    #         $self->warn(ref($self).": identifier [$id] is only numbers - ".
    #                     "prefixing with 'GO:'");
    #         return "GO:" . $id;
    #     }
    # we shouldn't have gotten here if it's at least a remotely decent ID
    $self->throw( ref($self) . ": non-standard identifier '$id'\n" )
      unless $id =~ /\|/;
    return $id;
}    # _get_id

1;
