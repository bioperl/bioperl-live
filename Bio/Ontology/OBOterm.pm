
=head1 NAME

Bio::Ontology::OBOterm - representation of OBO terms

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

This is data holder class for OBO terms. It is currently a dummy class since we anticipate that the
OBO term will become more richer with more features being added to OBO flat-file format.

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
the bugs and their resolution.  Bug reports can be submitted via the web:

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

=head1 APPENDIX

The rest of the documentation details each of the object
methods.

=cut

# Let the code begin...

package Bio::Ontology::OBOterm;
use strict;

use constant TRUE  => 1;
use constant FALSE => 0;

use base qw(Bio::Ontology::Term);

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



1;