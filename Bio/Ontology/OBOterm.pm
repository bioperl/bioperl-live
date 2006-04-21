
=head1 NAME

OBOterm - representation of OBO terms

=head1 SYNOPSIS

  $term = Bio::Ontology::OBOterm->new
    ( -identifier       => "GO:0005623",
      -name        => "Cell",
      -definition  => "The basic structural and functional unit ...",
      -is_obsolete => 0,
      -comment     => "" );

  $term->add_reference( @refs );
  $term->add_secondary_id( @ids );
  $term->add_synonym( @synonym );

  # etc.

=head1 DESCRIPTION

This is class for OBO terms

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org             - General discussion
  http://bio.perl.org/MailList.html - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the web:

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

=head1 APPENDIX

The rest of the documentation details each of the object
methods.

=cut

# Let the code begin...

package Bio::Ontology::OBOterm;
use vars qw( @ISA );
use strict;
use Bio::Ontology::Term;

use constant TRUE  => 1;
use constant FALSE => 0;

@ISA = qw( Bio::Ontology::Term );

=head2 new

 Title   : new
 Usage   :   $term = Bio::Ontology::OBOterm->new
     ( -identifier       => "GO:0005623",
      -name        => "Cell",
      -definition  => "The basic structural and functional unit ...",
      -is_obsolete => 0,
      -comment     => "" );

 Function: Creates a new Bio::Ontology::OBOterm.
 Returns : A new Bio::Ontology::OBOterm object.
 Args    : -identifier    => the id of this OBO term [GO:nnnnnnn]
                             integer of seven digits)
           -name          => the name of this OBO term [scalar]
           -definition    => the definition of this OBO term [scalar]
           -ontology      => the ontology for this term (a
                             Bio::Ontology::OntologyI compliant object)
           -version       => version information [scalar]
           -is_obsolete   => the obsoleteness of this OBO term [0 or 1]
           -comment       => a comment [scalar]

=cut

sub new {

    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    return $self;
}    # new

=head2 has_dblink

  Title   : has_dblink
  Usage   : $term->has_dblink($dblink);
  Function: Checks if a DBXref is already existing in the OBOterm object
  Return  : TRUE/FALSE
  Args    : [arg1] A DBxref identifier

=cut

sub has_dblink {
    my ( $self, $value ) = @_;
    return unless defined $value;
    my $context = "_default";
    $self->throw("'all' is a reserved word for context.") if $context eq 'all';
    $context ||= '_default';
    if ( ( $self->{_dblinks}->{$context} ) && grep { $_ eq $value }
        @{ $self->{_dblinks}->{$context} } )
    {
        return TRUE;
    }
    else {
        return FALSE;
    }
}

1;
