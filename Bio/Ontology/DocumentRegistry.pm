# $Id$
#
# BioPerl module for Bio::Ontology::DocumentRegistry
#
# Cared for by Allen Day <allenday@ucla.edu>
#
# Copyright Allen Day
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Ontology::DocumentRegistry - Keep track of where to find ontologies.
Allows lookups by name.

=head1 SYNOPSIS

my $registry = Bio::Ontology::DocumentRegistry->get_instance();
my($ont,$def,$fmt) = $registry->documents('Sequence Ontology');

my $io = Bio::OntologyIO->new(-url => $ont , -defs_url => $def , -format => $fmt);
my $so = $io->next_ontology();
#...

=head1 DESCRIPTION

Don't use this directly, use Bio::Ontology::OntologyStore instead.  Bio::Ontology::OntologyStore
uses Bio::Ontology::DocumentRegistry to load and cache ontologies as object graphs, you can just
ask it for what you want by name.  See L<Bio::Ontology::OntologyStore> for details.

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

=head1 AUTHOR - Allen Day

Email allenday@ucla.edu

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Ontology::DocumentRegistry;
use strict;
use base qw(Bio::Root::Root);
use Data::Dumper;

my $instance;

BEGIN {
$instance = {
             'Sequence Ontology' => {
                ontology => "http://cvs.sourceforge.net/viewcvs.py/*checkout*/song/ontology/so.ontology?rev=HEAD",
                definitions => "http://cvs.sourceforge.net/viewcvs.py/*checkout*/song/ontology/so.definition?rev=HEAD",
                format => 'soflat',
                                    },
            };

bless $instance, __PACKAGE__;
}


sub new {
  return shift->get_instance(@_);
}

=head2 get_instance()

=over

=item Usage

  my $singleton = Bio::Ontology::DocumentRegistry->get_instance();

=item Function

constructor

=item Returns

the Bio::Ontology::DocumentRegistry singleton.

=item Arguments

None

=back

=cut

sub get_instance {
  my ($self,%arg) = @_;


}


sub get_instance {
  return $instance;
}

=head2 documents()

=over

=item Usage

  my($ontology_url, $definitions_url, $format) = $obj->documents('Sequence Ontology');

=item Function

  maps an ontology name to a list of (local or) remote URIs where the files can be located.

=item Returns

  a 3-item list:

    (1) URI for the ontology file
    (2) URI for the ontology definitions file
    (3) format of the files (dagedit, obo, etc)

=item Arguments

  name of an ontology, e.g. 'Sequence Ontology', or 'Cellular Component (Gene Ontology)'

=back

=cut


sub documents {
  my($self,$name) = @_;

  if(defined($self->{$name})){
    return ($self->{$name}{ontology} , $self->{$name}{definitions}, $self->{$name}{format});
  } else {
    return ();
  }
}
