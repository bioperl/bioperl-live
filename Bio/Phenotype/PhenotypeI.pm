#
# BioPerl module for Bio::Phenotype::PhenotypeI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Christian M. Zmasek <czmasek-at-burnham.org> or <cmzmasek@yahoo.com>
#
# (c) Christian M. Zmasek, czmasek-at-burnham.org, 2002.
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

Bio::Phenotype::PhenotypeI - An interface for classes modeling phenotypes

=head1 SYNOPSIS

  #get Bio::Phenotype::PhenotypeI somehow

  print $phenotype->name(), "\n";
  print $phenotype->description(), "\n";

  my @keywords = ( "achondroplasia", "dwarfism" );
  $phenotype->add_keywords( @keywords ); 
  foreach my $keyword ( $phenotype->each_keyword() ) {
       print $keyword, "\n";
  }
  $phenotype->remove_keywords();


  foreach my $gene_symbol ( $phenotype->each_gene_symbol() ) {
       print $gene_symbol, "\n";
  }

  foreach my $corr ( $phenotype->each_Correlate() ) {
       # Do something with $corr
  }

  foreach my $var ( $phenotype->each_Variant() ) {
       # Do something with $var (mutation)
  }

  foreach my $measure ( $phenotype->each_Measure() ) {
       # Do something with $measure
  }


=head1 DESCRIPTION

This superclass defines common methods for classes modelling phenotypes.
Bio::Phenotype::OMIM::OMIMentry is an example of an instantiable phenotype
class (the design of this interface was partially guided by the need
to model OMIM entries).
Please note. This interface provides methods to associate mutations
(methods "each_Variant", ...) and genotypes (methods "each_Genotype", ...) 
with phenotypes. Yet, these aspects might need some future enhancements,
especially since there is no "genotype" class yet.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the 
Bioperl mailing lists  Your participation is much appreciated.

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

report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Christian M. Zmasek

Email: czmasek-at-burnham.org  or  cmzmasek@yahoo.com

WWW:   http://monochrome-effect.net/

Address: 

  Genomics Institute of the Novartis Research Foundation
  10675 John Jay Hopkins Drive
  San Diego, CA 92121

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Phenotype::PhenotypeI;
use base qw(Bio::Root::RootI);



=head2 name

 Title   : name
 Usage   : $obj->name( "r1" );
           or
           print $obj->name();
 Function: Set/get for the name or id of this phenotype.
 Returns : A name or id [scalar].
 Args    : A name or id [scalar] (optional).

=cut

sub name {
    my ( $self ) = @_;

    $self->throw_not_implemented();

} # name




=head2 description

 Title   : description
 Usage   : $obj->description( "This is ..." );
           or
           print $obj->description();
 Function: Set/get for the description of this phenotype.
 Returns : A description [scalar].
 Args    : A description [scalar] (optional).

=cut

sub description {
     my ( $self ) = @_;

    $self->throw_not_implemented();

} # description




=head2 species

 Title   : species
 Usage   : $obj->species( $species );
           or
           $species = $obj->species();
 Function: Set/get for the species of this phenotype.
 Returns : A species [Bio::Species].
 Args    : A species [Bio::Species] (optional).

=cut

sub species {
    my ( $self ) = @_;

    $self->throw_not_implemented();

} # species




=head2 comment

 Title   : comment
 Usage   : $obj->comment( "putative" );
           or
           print $obj->comment();
 Function: Set/get for a comment about this phenotype.
 Returns : A comment [scalar].
 Args    : A comment [scalar] (optional).

=cut

sub comment {
    my ( $self ) = @_;

    $self->throw_not_implemented();

} # comment




=head2 each_gene_symbol

 Title   : each_gene_symbol()
 Usage   : @gs = $obj->each_gene_symbol();                 
 Function: Returns a list of gene symbols [scalars, most likely Strings]
           associated with this phenotype.
 Returns : A list of scalars.
 Args    :

=cut

sub each_gene_symbol {
    my ( $self ) = @_;

    $self->throw_not_implemented();

} # each_gene_symbol


=head2 add_gene_symbols

 Title   : add_gene_symbols
 Usage   : $obj->add_gene_symbols( @gs );
           or
           $obj->add_gene_symbols( $gs );                  
 Function: Pushes one or more gene symbols [scalars, most likely Strings]
           into the list of gene symbols.
 Returns : 
 Args    : scalar(s).

=cut

sub add_gene_symbols {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # add_gene_symbols


=head2 remove_gene_symbols

 Usage   : $obj->remove_gene_symbols();
 Function: Deletes (and returns) the list of gene symbols [scalars,
           most likely Strings] associated with this phenotype.
 Returns : A list of scalars.
 Args    :

=cut

sub remove_gene_symbols {
    my ( $self ) = @_;

    $self->throw_not_implemented();

} # remove_gene_symbols




=head2 each_Variant

 Title   : each_Variant()
 Usage   : @vs = $obj->each_Variant();                 
 Function: Returns a list of Bio::Variation::VariantI implementing objects
           associated with this phenotype.
           This is for representing the actual mutation(s) causing this 
           phenotype.
           {* The "variants" data member and its methods will/might need to be
           changed/improved in one way or another, CZ 09/06/02 *}
 Returns : A list of Bio::Variation::VariantI implementing objects.
 Args    :

=cut

sub each_Variant {
    my ( $self ) = @_;

    $self->throw_not_implemented();

} # each_Variant


=head2 add_Variants

 Usage   : $obj->add_Variants( @vs );
           or
           $obj->add_Variants( $v );                  
 Function: Pushes one or more Bio::Variation::VariantI implementing objects
           into the list of Variants.
 Returns : 
 Args    : Bio::Variation::VariantI implementing object(s).

=cut

sub add_Variants {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # add_Variants


=head2 remove_Variants

 Title   : remove_Variants
 Usage   : $obj->remove_Variants();
 Function: Deletes (and returns) the list of Bio::Variation::VariantI implementing
           objects associated with this phenotype.
 Returns : A list of Bio::Variation::VariantI implementing objects.
 Args    :

=cut

sub remove_Variants {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # remove_Variants




=head2 each_Reference

 Title   : each_Reference()
 Usage   : @refs = $obj->each_Reference();                 
 Function: Returns a list of Bio::Annotation::Reference objects
           associated with this phenotype.
 Returns : A list of Bio::Annotation::Reference objects.
 Args    :

=cut

sub each_Reference {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # each_Reference


=head2 add_References 

 Title   : add_References
 Usage   : $obj->add_References( @refs );
           or
           $obj->add_References( $ref );                  
 Function: Pushes one or more Bio::Annotation::Reference objects
           into the list of References.
 Returns : 
 Args    : Bio::Annotation::Reference object(s).

=cut

sub add_References {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # add_References


=head2 remove_References

 Title   : remove_References()
 Usage   : $obj->remove_References();
 Function: Deletes (and returns) the list of Bio::Annotation::Reference objects
           associated with this phenotype.
 Returns : A list of Bio::Annotation::Reference objects.
 Args    :

=cut

sub remove_References {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # remove_References




=head2 each_CytoPosition

 Title   : each_CytoPosition()
 Usage   : @cps = $obj->each_CytoPosition();                 
 Function: Returns a list of Bio::Map::CytoPosition objects
           associated with this phenotype.
 Returns : A list of Bio::Map::CytoPosition objects.
 Args    :

=cut

sub each_CytoPosition {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # each_CytoPosition


=head2 add_CytoPositions

 Title   : add_CytoPositions
 Usage   : $obj->add_CytoPositions( @cps );
           or
           $obj->add_CytoPositions( $cp );                  
 Function: Pushes one or more Bio::Map::CytoPosition objects
           into the list of CytoPositions.
 Returns : 
 Args    : Bio::Map::CytoPosition object(s).

=cut

sub add_CytoPositions {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # add_CytoPositions


=head2 remove_CytoPositions

 Title   : remove_CytoPositions
 Usage   : $obj->remove_CytoPositions();
 Function: Deletes (and returns) the list o fBio::Map::CytoPosition objects
           associated with this phenotype.
 Returns : A list of Bio::Map::CytoPosition objects.
 Args    :

=cut

sub remove_CytoPositions {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # remove_CytoPositions




=head2 each_Correlate

 Title   : each_Correlate()
 Usage   : @corrs = $obj->each_Correlate();                 
 Function: Returns a list of Bio::Phenotype::Correlate objects
           associated with this phenotype.
           (Correlates are correlating phenotypes in different species;
           inspired by mouse correlates of human phenotypes in the OMIM
           database.)
 Returns : A list of Bio::Phenotype::Correlate objects.
 Args    :

=cut

sub each_Correlate {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # each_Correlate




=head2 add_Correlates

 Title   : add_Correlates
 Usage   : $obj->add_Correlates( @corrs );
           or
           $obj->add_Correlates( $corr );                  
 Function: Pushes one or more Bio::Phenotype::Correlate objects
           into the list of Correlates.
 Returns : 
 Args    : Bio::Phenotype::Correlate object(s).

=cut

sub add_Correlates {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # add_Correlates


=head2 remove_Correlates

 Title   : remove_Correlates
 Usage   : $obj->remove_Correlates();
 Function: Deletes (and returns) the list of Bio::Phenotype::Correlate objects
           associated with this phenotype.
 Returns : A list of Bio::Phenotype::Correlate objects.
 Args    :

=cut

sub remove_Correlates {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # remove_Correlates




=head2 each_Measure

 Title   : each_Measure()
 Usage   : @ms = $obj->each_Measure();                 
 Function: Returns a list of Bio::Phenotype::Measure objects
           associated with this phenotype.
           (Measure is for biochemically defined phenotypes
           or any other types of measures.)
 Returns : A list of Bio::Phenotype::Measure objects.
 Args    :

=cut

sub each_Measure {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # each_Measure


=head2 add_Measures

 Title   : add_Measures
 Usage   : $obj->add_Measures( @ms );
           or
           $obj->add_Measures( $m );                  
 Function: Pushes one or more Bio::Phenotype::Measure objects
           into the list of Measures.
 Returns : 
 Args    : Bio::Phenotype::Measure object(s).

=cut

sub add_Measures {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # add_Measures


=head2 remove_Measures

 Title   : remove_Measures
 Usage   : $obj->remove_Measures();
 Function: Deletes (and returns) the list of Bio::Phenotype::Measure objects
           associated with this phenotype.
 Returns : A list of Bio::Phenotype::Measure objects.
 Args    :

=cut

sub remove_Measures {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # remove_Measures




=head2 each_keyword

 Title   : each_keyword()
 Usage   : @kws = $obj->each_keyword();                 
 Function: Returns a list of key words [scalars, most likely Strings]
           associated with this phenotype.
 Returns : A list of scalars.
 Args    :

=cut

sub each_keyword {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # each_keyword


=head2 add_keywords

 Title   : add_keywords
 Usage   : $obj->add_keywords( @kws );
           or
           $obj->add_keywords( $kw );                  
 Function: Pushes one or more keywords [scalars, most likely Strings]
           into the list of key words.
 Returns : 
 Args    : scalar(s).

=cut

sub add_keywords {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # add_keywords


=head2 remove_keywords

 Title   : remove_keywords
 Usage   : $obj->remove_keywords();
 Function: Deletes (and returns) the list of key words [scalars,
           most likely Strings] associated with this phenotype.
 Returns : A list of scalars.
 Args    :

=cut

sub remove_keywords {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # remove_keywords




=head2 each_DBLink

 Title   : each_DBLink()
 Usage   : @dbls = $obj->each_DBLink();                 
 Function: Returns a list of Bio::Annotation::DBLink objects
           associated with this phenotype.
 Returns : A list of Bio::Annotation::DBLink objects.
 Args    :

=cut

sub each_DBLink {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
}


=head2 add_DBLinks

 Title   : add_DBLinks
 Usage   : $obj->add_DBLinks( @dbls );
           or
           $obj->add_DBLinks( $dbl );                  
 Function: Pushes one or more Bio::Annotation::DBLink objects
           into the list of DBLinks.
 Returns : 
 Args    : Bio::Annotation::DBLink object(s).

=cut

sub add_DBLinks {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # add_DBLinks


=head2 remove_DBLinks

 Title   : remove_DBLinks
 Usage   : $obj->remove_DBLinks();
 Function: Deletes (and returns) the list of Bio::Annotation::DBLink objects
           associated with this phenotype.
 Returns : A list of Bio::Annotation::DBLink objects.
 Args    :

=cut

sub remove_DBLinks {
    my ( $self ) = @_;

    $self->throw_not_implemented();

} # remove_DBLinks




=head2 each_Genotype

 Title   : each_Reference()
 Usage   : @gts = $obj->each_Reference();                 
 Function: Returns a list of "Genotype" objects
           associated with this phenotype.
           {* the "genotypes" data member and its methods certainly will/needs to be
           changed/improved in one way or another since there is
           no "Genotype" class yet, CZ 09/06/02 *}
 Returns : A list of "Genotype" objects.
 Args    :

=cut

sub each_Genotype {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # each_Genotype


=head2 add_Genotypes

 Title   : add_Genotypes
 Usage   : $obj->add_Genotypes( @gts );
           or
           $obj->add_Genotypes( $gt );                  
 Function: Pushes one or more "Genotypes"
           into the list of "Genotypes".
 Returns : 
 Args    : "Genotypes(s)".

=cut

sub add_Genotypes {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # add_Genotypes


=head2 remove_Genotypes

 Title   : remove_Genotypes
 Usage   : $obj->remove_Genotypes();
 Function: Deletes (and returns) the list of "Genotype" objects
           associated with this phenotype.
 Returns : A list of "Genotype" objects.
 Args    :

=cut

sub remove_Genotypes {
    my ( $self ) = @_;

    $self->throw_not_implemented();
    
} # remove_Genotypes





1;
