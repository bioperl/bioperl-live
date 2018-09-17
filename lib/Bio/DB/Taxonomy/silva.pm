#
# BioPerl module for Bio::DB::Taxonomy::silva
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Copyright Florent Angly
#
# You may distribute this module under the same terms as perl itself


=head1 NAME

Bio::DB::Taxonomy::silva - Use the Silva taxonomy

=head1 SYNOPSIS

  use Bio::DB::Taxonomy;

  my $db = Bio::DB::Taxonomy->new(
     -source   => 'silva',
     -taxofile => 'SSURef_108_tax_silva_trunc.fasta',
  );

=head1 DESCRIPTION

This is an implementation of Bio::DB::Taxonomy which stores and accesses the
Silva taxonomy. Internally, Bio::DB::Taxonomy::silva keeps the taxonomy
into memory by using Bio::DB::Taxonomy::list. As a consequence, note that the
IDs assigned to the taxonomy nodes, e.g. sv72, are arbitrary, contrary to the
pre-defined IDs that NCBI assigns to taxons. Note also that no rank names or
common names are assigned to the taxa of Bio::DB::Taxonomy::silva.

The latest Silva taxonomy (2011) contains about 126,000 taxa and occupies 
about 124 MB of memory once parsed into a Bio::DB::Taxonomy::silva object.
Obviously, it can take a little while to load.

The taxonomy file SSURef_108_tax_silva_trunc.fasta that this module uses is
available from L<http://www.arb-silva.de/no_cache/download/archive/release_108/Exports/>.

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

=head1 AUTHOR - Florent Angly

florent.angly@gmail.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


package Bio::DB::Taxonomy::silva;

use strict;
use Bio::SeqIO;

use base qw(Bio::DB::Taxonomy Bio::DB::Taxonomy::list);

$Bio::DB::Taxonomy::list::prefix = 'sv';


=head2 new

 Title   : new
 Usage   : my $obj = Bio::DB::Taxonomy::silva->new();
 Function: Builds a new Bio::DB::Taxonomy::silva object 
 Returns : an instance of Bio::DB::Taxonomy::silva
 Args    : -taxofile  => name of the FASTA file containing the taxonomic information,
                         typically 'SSURef_108_tax_silva_trunc.fasta' (mandatory)

=cut

sub new {
   # Override Bio::DB::Taxonomy
   my($class, @args) = @_;
   my $self = $class->SUPER::new(@args);
   my ($taxofile) = $self->_rearrange([qw(TAXOFILE)], @args);
  
   if ( $taxofile ) {
      $self = $self->_build_taxonomy($taxofile);
   }

   return $self;
}


sub _build_taxonomy {
   my ($self, $taxofile) = @_;

   my $taxonomy = Bio::DB::Taxonomy::list->new();
   my %taxas;
   my $desc_re = qr/^>\S+?(?:\s+(.+))?$/;

   # One could open the file using Bio::SeqIO::fasta, but it is slower and we
   # only need the sequence descriptions

   open my $in, '<', $taxofile or $self->throw("Could not read file '$taxofile': $!");

   # Populate taxonomy with taxonomy obtained from sequence description
   while (my $line = <$in>) {

      next if $line !~ $desc_re;
      my $taxo_string = $1;
      next if not $taxo_string;

      # Example of taxonomy string:
      # 1/ Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus faecium DO
      # 2/ Eukaryota;Metazoa;Chordata;Craniata;Vertebrata;Euteleostomi;Mammalia;Eutheria;Euarchontoglires;Glires;
      #       Rodentia;Sciurognathi;Muroidea;Muridae;Murinae;Rattus;;Rattus norvegicus (Norway rat)
      
      # Skip already seen taxas
      next if exists $taxas{$taxo_string};
      $taxas{$taxo_string} = undef;

      # Strip the common name (could save it if Bio::DB::Taxonomy::list supported it)
      $taxo_string =~ s/ \(.*\)$//;

      # Save lineage
      # Unfortunately, we cannot easily add ranks since they vary from 2 to 23 for every entry
      my @names = split /;/, $taxo_string;
      $taxonomy->add_lineage(
         -names => \@names,
      );

   }

   close $in;

   return $taxonomy;
}


1;
