# $Id$
#
# BioPerl module for Bio::DB::Taxonomy
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::Taxonomy - Access to a taxonomy database

=head1 SYNOPSIS

  use Bio::DB::Taxonomy;
  my $db = new Bio::DB::Taxonomy(-source => 'entrez');
  # use NCBI Entrez over HTTP
  my $taxonid = $db->get_taxonid('Homo sapiens');

=head1 DESCRIPTION

This is a front end module for access to a taxonomy database.

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

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

Describe contact details here

=head1 CONTRIBUTORS

Sendu Bala: bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::Taxonomy;
use vars qw(@ISA $DefaultSource);
use strict;

use Bio::Root::HTTPget;
$DefaultSource = 'entrez';

@ISA = qw(Bio::Root::HTTPget);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::DB::Taxonomy(-source => 'entrez');
 Function: Builds a new Bio::DB::Taxonomy object.
 Returns : an instance of Bio::DB::Taxonomy
 Args    : -source => which database source 'entrez' or 'flatfile'

=cut

sub new {
  my($class,@args) = @_;

  if( $class =~ /Bio::DB::Taxonomy::(\S+)/ ) {
      my ($self) = $class->SUPER::new(@args);
      $self->_initialize(@args);
      return $self;
  } else { 
      my %param = @args;
      @param{ map { lc $_ } keys %param } = values %param; # lowercase keys
      my $source = $param{'-source'} || $DefaultSource;

      $source = "\L$source";	# normalize capitalization to lower case

      # normalize capitalization
      return unless( $class->_load_tax_module($source) );
      return "Bio::DB::Taxonomy::$source"->new(@args);
  }
}

# empty for now
sub _initialize { }

=head2 get_Taxonomy_Node

 Title   : get_Taxonomy_Node
 Usage   : my $species = $db->get_Taxonomy_Node(-taxonid => $taxonid)
 Function: Get a Bio::Taxonomy::Taxon object for a taxonid.
 Returns : Bio::Taxonomy::Taxon object
 Args    : -taxonid => taxonomy id (to query by taxonid)
            OR
           -name   => string (to query by a taxonomy name: common name, 
                              species, genus, etc)

=cut

sub get_Taxonomy_Node {
    shift->throw_not_implemented();
}

=head2 get_taxonids

 Title   : get_taxonids
 Usage   : my @taxonids = $db->get_taxonids('Homo sapiens');
 Function: Searches for a taxonid (typically ncbi_taxon_id) based on a query
           string. Note that multiple taxonids can match to the same supplied
           name.
 Returns : array of integer ids in list context, one of these in scalar context
 Args    : string representing taxanomic (node) name

=cut

sub get_taxonids {
    shift->throw_not_implemented();
}

*get_taxonid = \&get_taxonids;
*get_taxaid = \&get_taxonids;

=head2 _load_tax_module

 Title   : _load_tax_module
 Usage   : *INTERNAL Bio::DB::Taxonomy stuff*
 Function: Loads up (like use) a module at run time on demand

=cut

sub _load_tax_module {
    my ($self, $source) = @_;
    my $module = "Bio::DB::Taxonomy::" . $source;
    my $ok;

    eval { $ok = $self->_load_module($module) };
    if ( $@ ) {
	print STDERR $@;
	print STDERR <<END;
$self: $source cannot be found
Exception $@
For more information about the Bio::DB::Taxonomy system please see
the Bio::DB::Taxonomy docs.  This includes ways of checking for 
formats at compile time, not run time.
END
  ;
    }
    return $ok;
}

1;
