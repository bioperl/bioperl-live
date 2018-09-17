#
# BioPerl module for Bio::Ontology::DocumentRegistry
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
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

  my $io = Bio::OntologyIO->new(-url => $ont,
                                -defs_url => $def,
                                -format => $fmt);
  my $so = $io->next_ontology();
  #...

=head1 DESCRIPTION

Do not use this directly, use Bio::Ontology::OntologyStore instead.
Bio::Ontology::OntologyStore uses Bio::Ontology::DocumentRegistry to
load and cache ontologies as object graphs, you can just ask it for
what you want by name.  See L<Bio::Ontology::OntologyStore> for
details.

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

=head1 AUTHOR - Allen Day

Email allenday@ucla.edu

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
	    ontology => "http://song.cvs.sourceforge.net/*checkout*/song/ontology/so.ontology?rev=HEAD",
        definitions => "http://song.cvs.sourceforge.net/*checkout*/song/ontology/so.definition?rev=HEAD",
        format => 'soflat',
                                    },
   'Sequence Ontology OBO' => {
	    ontology => "http://song.cvs.sourceforge.net/*checkout*/song/ontology/so.obo?rev=HEAD",
        definitions => "http://song.cvs.sourceforge.net/*checkout*/song/ontology/so.definition?rev=HEAD",
        format => 'obo',
                                    },
   
   #### TODO Server http://umn.dl.sourceforge.net/ does not respond, are there
   #### alternative sources? 
   'Sequence Ontology Feature Annotation' => {
        ontology => 'http://umn.dl.sourceforge.net/sourceforge/song/sofa.ontology',
        definitions =>'http://umn.dl.sourceforge.net/sourceforge/song/sofa.definition',
        format => 'soflat',
                                    },
    'Gene Ontology' => {
         ontology => [
							 'http://www.geneontology.org/ontology/function.ontology',
							 'http://www.geneontology.org/ontology/process.ontology',
							 'http://www.geneontology.org/ontology/component.ontology'
							],
			definitions => 'http://www.geneontology.org/ontology/GO.defs',
         format => 'soflat',
							  },
            };

#aliases
$instance->{Gene_Ontology} = $instance->{'Gene Ontology'};

bless $instance, __PACKAGE__;
}


sub new {
  return shift->get_instance(@_);
}

=head2 get_instance

 Title   : get_instance
 Usage   : my $singleton = Bio::Ontology::DocumentRegistry->get_instance();
 Function: constructor
 Returns : The Bio::Ontology::DocumentRegistry singleton.
 Args    : None
 Usage

=cut

sub get_instance {
  return $instance;
}

=head2 documents

 Title   : documents
 Usage   : my($ontology_url, $definitions_url, $format) = $obj->documents('Sequence Ontology');
 Function: Maps an ontology name to a list of (local or) remote URIs where the
           files can be located.
 Returns : A 3-item list:
           (1) URI for the ontology file
           (2) URI for the ontology definitions file
           (3) format of the files (dagedit, obo, etc)
 Args    : Name of an ontology, e.g. 'Sequence Ontology', or 'Cellular Component 
           (Gene Ontology)'

=cut


sub documents {
  my($self,$name) = @_;

  if(defined($self->{$name})){
    return ($self->{$name}{ontology} , $self->{$name}{definitions}, $self->{$name}{format});
  } else {
    return ();
  }
}

1;
