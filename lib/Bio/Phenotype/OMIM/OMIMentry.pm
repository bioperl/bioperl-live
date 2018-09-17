#
# BioPerl module for Bio::Phenotype::OMIM::OMIMentry
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

Bio::Phenotype::OMIM::OMIMentry - represents OMIM (Online Mendelian
Inheritance in Man) database entries

=head1 SYNOPSIS

  $obj = Bio::Phenotype::OMIM::OMIMentry->new( -mim_number          => 200000,
                                               -description         => "This is ...",
                                               -more_than_two_genes => 1 );

=head1 DESCRIPTION

Inherits from Bio::Phenotype::PhenotypeI.
Bio::Phenotype::OMIM::OMIMparser parses the flat file representation
of OMIM (i.e. files "omim.txt" and "genemap") returning OMIMentry objects. 

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


package Bio::Phenotype::OMIM::OMIMentry;
use strict;

use Bio::Phenotype::OMIM::MiniMIMentry;
use Bio::Phenotype::OMIM::OMIMentryAllelicVariant;

use constant TRUE              => 1;
use constant FALSE             => 0;
use constant DEFAULT_MIM_NUMER => 0;

use base qw(Bio::Phenotype::Phenotype);




=head2 new

 Title   : new
 Usage   : $obj = Bio::Phenotype::OMIM::OMIMentry->new( -mim_number          => 200000,
                                                        -description         => "This is ...",
                                                        -more_than_two_genes => 1 );                      
 Function: Creates a new OMIMentry object.
 Returns : A new OMIMentry object.
 Args    : -mim_number                     => the MIM number
           -title                          => the title or name
           -alternative_titles_and_symbols => the "alternative titles and symbols"    
           -more_than_two_genes            => can phenotype can be caused by mutation in any of two or more genes?       
           -is_separate                    => is this phenotype separate from those represented by other entries  
           -description                    => the description of this phenotype
           -mapping_method                 => the mapping method      
           -gene_status                    => the gene status of this       
           -comment                        => a comment        
           -species                        => ref to the the species (human)
           -created                        => created by whom/when       
           -edited                         => edited by whom/when    
           -contributors                   => contributed by whom/when 
           -additional_references          => "see also"     
           -clinical_symptoms              => the clinical symptoms
           -minimim                        => the Mini MIM associated with this OMIM antry

=cut

sub new {

    my( $class,@args ) = @_;
    
    my $self = $class->SUPER::new( @args );
   
    my ( $mim_number,
         $title,
         $alternative_titles_and_symbols,     
         $more_than_two_genes,       
         $is_separate,    
         $description,
         $mapping_method,     
         $gene_status,       
         $comment,        
         $species,
         $created,       
         $edited,    
         $contributors,
         $additional_references,     
         $clinical_symptoms, 
         $miniMIM )
    = $self->_rearrange( [ qw( MIM_NUMBER
                               TITLE
                               ALTERNATIVE_TITLES_AND_SYMBOLS
                               MORE_THAN_TWO_GENES
                               IS_SEPARATE
                               DESCRIPTION
                               MAPPING_METHOD
                               GENE_STATUS
                               COMMENT
                               SPECIES
                               CREATED
                               EDITED
                               CONTRIBUTORS
                               ADDITIONAL_REFERENCES
                               CLINICAL_SYMPTOMS
                               MINIMIM ) ], @args );
   
    $self->init(); 
    
    $mim_number                     && $self->MIM_number( $mim_number );
    $title                          && $self->title( $title );
    $alternative_titles_and_symbols && $self->alternative_titles_and_symbols( $alternative_titles_and_symbols );     
    $more_than_two_genes            && $self->more_than_two_genes( $more_than_two_genes );      
    $is_separate                    && $self->is_separate( $is_separate );   
    $description                    && $self->description( $description );
    $mapping_method                 && $self->mapping_method( $mapping_method );     
    $gene_status                    && $self->gene_status( $gene_status );       
    $comment                        && $self->comment( $comment );        
    $species                        && $self->species( $species );
    $created                        && $self->created( $created );       
    $edited                         && $self->edited( $edited );    
    $contributors                   && $self->contributors( $contributors );
    $additional_references          && $self->additional_references( $additional_references );     
    $clinical_symptoms              && $self->clinical_symptoms_raw( $clinical_symptoms );
    $miniMIM                        && $self->miniMIM( $miniMIM );
                                                    
    return $self;
    
} # new



=head2 init

 Title   : init()
 Usage   : $obj->init();   
 Function: Initializes this OMIMentry to all "" and empty lists.
 Returns : 
 Args    :

=cut

sub init {

    my( $self ) = @_;

    $self->MIM_number( DEFAULT_MIM_NUMER );
    $self->title( "" );
    $self->alternative_titles_and_symbols( "" );
    $self->more_than_two_genes( FALSE );
    $self->is_separate( FALSE );
    $self->description( "" );
    $self->mapping_method( "" );
    $self->gene_status( "" );
    $self->comment( "" );
    my $species = Bio::Species->new();
    $species->classification( qw( sapiens Homo ) );
    $self->species( $species );
    $self->created( "" );
    $self->edited( "" );
    $self->contributors( "" );
    $self->additional_references( "" );
    $self->clinical_symptoms( {} );
    $self->remove_Correlates();
    $self->remove_References();
    $self->remove_AllelicVariants();
    $self->remove_CytoPositions();
    $self->remove_gene_symbols();
    $self->remove_Genotypes();
    $self->remove_DBLinks();
    $self->remove_keywords();
    $self->remove_Variants();
    $self->remove_Measures();
    $self->miniMIM( Bio::Phenotype::OMIM::MiniMIMentry->new() );
  
} # init



sub to_string {

    my( $self ) = @_;

    my $s = "";

    $s .= "-- MIM number:\n";
    $s .= $self->MIM_number()."\n\n";
    $s .= "-- Title:\n";
    $s .= $self->title()."\n\n";
    $s .= "-- Alternative Titles and Symbols:\n";
    $s .= $self->alternative_titles_and_symbols()."\n\n";
    $s .= "-- Can be caused by Mutation in any of two or more Genes:\n";
    $s .= $self->more_than_two_genes()."\n\n";
    $s .= "-- Phenotype is separate:\n";
    $s .= $self->is_separate()."\n\n"; 
    $s .= "-- Description:\n";
    $s .= $self->description()."\n\n";
    $s .= "-- Species:\n";
    $s .= $self->species()->binomial()."\n\n";
    $s .= "-- Clinical Symptoms:\n";
    $s .= $self->clinical_symptoms()."\n\n";
    $s .= "-- Allelic Variants:\n";
    $s .= $self->_array_to_string( $self->each_AllelicVariant() )."\n";
    $s .= "-- Cyto Positions:\n";
    $s .= $self->_array_to_string( $self->each_CytoPosition() )."\n";
    $s .= "-- Gene Symbols:\n";
    $s .= $self->_array_to_string( $self->each_gene_symbol() )."\n";
    $s .= "-- Correlates:\n";
    $s .= $self->_array_to_string( $self->each_Correlate() )."\n";
    $s .= "-- References:\n";
    $s .= $self->_array_to_string( $self->each_Reference() )."\n";
    $s .= "-- Additional References:\n";
    $s .= $self->additional_references()."\n\n";
    $s .= "-- Mapping Method:\n";
    $s .= $self->mapping_method()."\n\n";
    $s .= "-- Gene status:\n";
    $s .= $self->gene_status()."\n\n";
    $s .= "-- Created:\n";
    $s .= $self->created()."\n\n";
    $s .= "-- Contributors:\n";
    $s .= $self->contributors()."\n\n";
    $s .= "-- Edited:\n";
    $s .= $self->edited()."\n\n";
    $s .= "-- Comment:\n";
    $s .= $self->comment()."\n\n";
    $s .= "-- MiniMIM:\n";
    $s .= $self->miniMIM()->to_string()."\n\n";
    return $s;
    

} # to_string



=head2 MIM_number

 Title   : MIM_number
 Usage   : $omim->MIM_number( "100050" );
           or
           print $omim->MIM_number();
 Function: Set/get for the MIM number of this OMIM entry.
 Returns : The MIM number [an integer larger than 100000].
 Args    : The MIM number [an integer larger than 100000] (optional).

=cut

sub MIM_number {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        if ( $value =~ /\D/
        || ( $value < 100000 && $value != DEFAULT_MIM_NUMER ) ) {
            $self->throw( "Found [$value]" 
            . " where [integer larger than 100000] expected" );
        }
        $self->{ "_MIM_number" } = $value;
    }

    return $self->{ "_MIM_number" };

} # MIM_number




=head2 title

 Title   : title
 Usage   : $omim->title( "AARSKOG SYNDROME" );
           or
           print $omim->title();
 Function: Set/get for the title or name of this OMIM entry.
           This method is an alias to the method "name" of
           Bio::Phenotype::PhenotypeI.
 Returns : The title [scalar].
 Args    : The title [scalar] (optional).

=cut

sub title {
    my $self = shift;
    
    $self->name(@_);
    
} # title




=head2 alternative_titles_and_symbols

 Title   : alternative_titles_and_symbols
 Usage   : $omim->alternative_titles_and_symbols( "AORTIC ANEURYSM, ABDOMINAL" );
           or
           print $omim->alternative_titles_and_symbols();
 Function: Set/get for the "alternative titles and symbols" of this OMIM entry.
           Currently, everything after the first line of title (TI) field is
           considered "alternative titles and symbols".
 Returns : "alternative titles and symbols" [scalar].
 Args    : "alternative titles and symbols" [scalar] (optional).

=cut

sub alternative_titles_and_symbols {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_alternative_titles_and_symbols" } = $value;
    }

    return $self->{ "_alternative_titles_and_symbols" };

} # alternative_titles_and_symbols




=head2 more_than_two_genes

 Title   : more_than_two_genes
 Usage   : $omim->more_than_two_genes( 1 );
           or
           print $omim->more_than_two_genes();
 Function: This is true if this phenotype can be caused
           by mutation in any of two or more genes.
           In OMIM, this is indicated by a number symbol (#)
           before an entry number (e.g. #114480 -- BREAST CANCER).
 Returns : [1 or 0].
 Args    : [1 or 0] (optional).

=cut

sub more_than_two_genes {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->_is_true_or_false( $value );
        $self->{ "_more_than_two_genes" } = $value;
    }

    return $self->{ "_more_than_two_genes" };

} # more_than_two_genes




=head2 is_separate

 Title   : is_separate
 Usage   : $omim->is_separate( 1 );
           or
           print $omim->is_separate();
 Function: This is true if the phenotype determined by the gene at
           the given locus is separate from those represented by
           other entries where "is_separate" is true and if the mode
           of inheritance of the phenotype has been proved
           (in the judgment of the authors and editors).
           In OMIM, this is indicated by a asterisk  (*)
           before an entry number (e.g. *113705 BREAST CANCER,
           TYPE 1; BRCA1).
 Returns : [1 or 0].
 Args    : [1 or 0] (optional).

=cut

sub is_separate {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->_is_true_or_false( $value );
        $self->{ "_is_separate" } = $value;
    }

    return $self->{ "_is_separate" };

} # is_separate




=head2 mapping_method

 Title   : mapping_method
 Usage   : $omim->mapping_method( "PCR of somatic cell hybrid DNA" );
           or
           print $omim->mapping_method();
 Function: Set/get for the mapping method of this OMIM entry.
 Returns : The mapping method [scalar].
 Args    : The mapping method [scalar] (optional).

=cut

sub mapping_method {
    my $self = shift;
    return $self->{ "_mapping_method" } = shift if(@_);
    return $self->{ "_mapping_method" };
} # mapping_method

=head2 gene_status

 Title   : gene_status
 Usage   : $omim->gene_status( "C" );
           or
           print $omim->gene_status();
 Function: Set/get for the gene status of this OMIM entry.
           The certainty with which assignment of loci to chromosomes or the linkage
           between two loci has been established has been graded into the following
           classes:
           <L>C = confirmed - observed in at least two laboratories or in several families.
           <L>P = provisional - based on evidence from one laboratory or one family.
           <L>I = inconsistent - results of different laboratories disagree.
           <L>L = limbo - evidence not as strong as that provisional, but included for
           heuristic reasons. (Same as `tentative'.)

 Returns :  [C, P, I, or L].
 Args    :  [C, P, I, or L] (optional).

=cut

sub gene_status {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        #unless ( $value eq "C"
        #      || $value eq "P"
        #      || $value eq "I"
        #      || $value eq "L"
        #      || $value eq "A"  # !?
        #      || $value eq "H"  # !?
        #      || $value eq "U"  # !?
        #      || $value eq "" ) {
        #    $self->throw( "Found [$value]" 
        #    . " where [C, P, I, or L] expected" );
        #}
        unless ( $value eq "C"
              || $value eq "P"
              || $value eq "I"
              || $value eq "L"
              || $value eq "" ) {
            $value = "";
        }
        
        $self->{ "_gene_status" } = $value;
    }

    return $self->{ "_gene_status" };

} # gene_status


=head2 clinical_symptoms

 Title   : clinical_symptoms
 Usage   : $omim->clinical_symptoms({});
 Function: Set/get for the clinical symptoms of this OMIM entry.
 Returns : [hash reference].
 Args    : [hash reference]. Suggested not to assign alone. Parser will do.

=cut

sub clinical_symptoms {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        unless(ref($value) eq 'HASH'){
            $self->throw('a hash referenced needed');
        }
        $self->{ "_clinical_symptoms" } = $value;
    }

    return $self->{ "_clinical_symptoms" };

} # clinical_symptoms

=head2 clinical_symptoms_raw

  Title     : clinical_symptoms_raw
  Usage     : $omim->clinical_symptoms( "Patients with ..." );
              print $omim->clinical_symptoms();
  Functions : Get/set for text information of clinical symptoms
  Returns   : The clinical symptoms [scalar].
  Args      : The clinical symptoms [scalar] (optional).

=cut 

sub clinical_symptoms_raw {
    my $self = shift;
    return $self->{_clinical_symptoms_raw} = shift if @_;
    return $self->{_clinical_symptoms_raw};
}

=head2 add_clinical_symptoms

  Title     : add_clinical_symptoms
  Usage     : $entry->add_clinical_symptoms('Ears', 'Floppy ears', 'Lop-ears');
  Function  : add one or more symptoms on one part of body.
  Returns   : [none]
  Args      : ($part, @symptoms)
              $part, the text name of part/organism of human
              @symptoms, an array of text description

=cut

sub add_clinical_symptoms {
    my ($self, $part, @symptoms) = @_;
    unless(defined $part){
        $self->throw('a part/organism must be assigned');
    }
    $self->{_clinical_symptoms} = {} unless $self->{_clinical_symptoms};
    $self->{_clinical_symptoms}->{$part} = [] 
        unless $self->{_clinical_symptoms}->{$part};
    push @{$self->{_clinical_symptoms}->{$part}}, @symptoms;
}

=head2 query_clinical_symptoms

  Title     : get_clinical_symptoms
  Usage     : @symptoms = $self->query_clinical_symptoms('Ears');
  Function  : get all symptoms specific to one part/organism.
  Returns   : an array of text
  Args      : $organ

=cut

sub query_clinical_symptoms {
    my ($self, $organ)=@_;
    my $symptoms=$self->{_clinical_symptoms}->{$organ};
    @$symptoms;
}

sub get_clinical_symptom_organs {
    my ($self)=@_;
    keys %{$self->{_clinical_symptoms}};
}

=head2 created

 Title   : created
 Usage   : $omim->created( "Victor A. McKusick: 6/4/1986" );
           or
           print $omim->created();
 Function: Set/get for the created field of the OMIM database.
 Returns : Name(s) and date(s) [scalar - free form].
 Args    : Name(s) and date(s) [scalar - free form] (optional).

=cut

sub created {
    my $self = shift;
    return $self->{ "_created" } = shift if(@_);
    return $self->{ "_created" };

} # created




=head2 contributors

 Title   : contributors
 Usage   : $omim->contributors( "Kelly A. Przylepa - revised: 03/18/2002" );
           or
           print $omim->contributors();
 Function: Set/get for the contributors field of the OMIM database.
 Returns : Name(s) and date(s) [scalar - free form].
 Args    : Name(s) and date(s) [scalar - free form] (optional).

=cut

sub contributors {
    my  $self = shift;
    $self->{ "_contributors" } = shift if(@_);
    return $self->{ "_contributors" };

} # contributors




=head2 edited

 Title   : edited
 Usage   : $omim->edited( "alopez: 06/03/1997" );
           or
           print $omim->edited();
 Function: Set/get for the edited field of the OMIM database.
 Returns : Name(s) and date(s) [scalar - free form].
 Args    : Name(s) and date(s) [scalar - free form] (optional).

=cut

sub edited {
    my $self = shift;
    return $self->{ "_edited" } = shift if(@_);
    return $self->{ "_edited" };

} # edited




=head2 additional_references

 Title   : additional_references
 Usage   : $omim->additional_references( "Miller er al." );
           or
           print $omim->additional_references();
 Function: Set/get for the additional references of this OMIM antry
           (see also).
 Returns : additional reference [scalar].
 Args    : additional reference [scalar] (optional).

=cut

sub additional_references {
    my $self = shift;
    return $self->{ "_additional_references" } = shift if(@_);
    return $self->{ "_additional_references" };

} # additional_references

=head2 miniMIM

 Title   : miniMIM
 Usage   : $omim->miniMIM( $MM );
           or
           $MM = $omim->miniMIM();
 Function: Set/get for the Mini MIM associated with this OMIM antry
           (see also).
 Returns : [Bio::Phenotype::OMIM::MiniMIMentry].
 Args    : [Bio::Phenotype::OMIM::MiniMIMentry] (optional).

=cut

sub miniMIM {

    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->_check_ref_type( $value, "Bio::Phenotype::OMIM::MiniMIMentry" );
        $self->{ "_mini_mim" } = $value;
    }
    
    return $self->{ "_mini_mim" };
}

=head2 each_AllelicVariant

 Title   : each_AllelicVariant()
 Usage   : @avs = $obj->each_AllelicVariant();                 
 Function: Returns a list of Bio::Phenotype::OMIM::OMIMentryAllelicVariant objects
           associated with this OMIM entry.
 Returns : A list of Bio::Phenotype::OMIM::OMIMentryAllelicVariant objects.
 Args    :

=cut

sub each_AllelicVariant {
    my ( $self ) = @_;
    
    return @{$self->{"_allelic_variants"}} if exists($self->{"_allelic_variants"});
    return ();    
} # each_AllelicVariant


=head2 add_AllelicVariants

 Title   : add_AllelicVariants
 Usage   : $obj->add_AllelicVariants( @avs );
           or
           $obj->add_AllelicVariants( $av );                  
 Function: Pushes one or more OMIMentryAllelicVariant
           into the list of OMIMentryAllelicVariants.
 Returns : 
 Args    : Bio::Phenotype::OMIM::OMIMentryAllelicVariant object(s).

=cut

sub add_AllelicVariants {
    my ( $self, @values ) = @_;
    
    return unless( @values );

    foreach my $value ( @values ) {  
        $self->_check_ref_type( $value, "Bio::Phenotype::OMIM::OMIMentryAllelicVariant" );
    }
        
    push( @{ $self->{ "_allelic_variants" } }, @values );
    
} # add_AllelicVariants


=head2 remove_AllelicVariants

 Title   : remove_AllelicVariants
 Usage   : $obj->remove_AllelicVariants();
 Function: Deletes (and returns) the list of OMIMentryAllelicVariant objects
           associated with this OMIM entry.
 Returns : A list of OMIMentryAllelicVariant objects.
 Args    :

=cut

sub remove_AllelicVariants {
    my ( $self ) = @_;
     
    my @a = $self->each_AllelicVariant();
    $self->{ "_allelic_variants" } = [];
    return @a;

} # remove_AllelicVariants


# Title   : _array_to_string         
# Function:
# Returns : 
# Args    : 
sub _array_to_string {
    my( $self, @value ) = @_;

    my $s = "";
    
    for ( my $i = 0; $i < scalar( @value ); ++$i ) {
        if ( ! ref( $value[ $i ] ) ) {
            $s .= "#" . $i . "\n-- Value:\n" . $value[ $i ] . "\n";
        }
        elsif ( $value[ $i ]->isa( "Bio::Phenotype::OMIM::OMIMentryAllelicVariant" ) 
        ||      $value[ $i ]->isa( "Bio::Phenotype::Correlate" ) ) {
            $s .= "#" . $i . "\n" . ( $value[ $i ] )->to_string() . "\n";
        }
        elsif ( $value[ $i ]->isa( "Bio::Annotation::Reference" ) ) {
            $s .= "#".$i."\n-- Authors:\n".( $value[ $i ] )->authors()."\n";
            $s .= "-- Title:\n".( $value[ $i ] )->title()."\n";
            $s .= "-- Location:\n".( $value[ $i ] )->location()."\n";
        }
        elsif ( $value[ $i ]->isa( "Bio::Map::CytoPosition" ) ) {
            $s .= "#" . $i . "\n-- Value:\n" . ( $value[ $i ] )->value() . "\n";
        }
    }
    
    return $s;
    
} # _array_to_string


# Title   :_is_true_or_false              
# Function: Checks whether the argument is 1 or 0.
# Returns : 
# Args    : The value to be checked.
sub _is_true_or_false {
    my ( $self, $value ) = @_;
    unless ( $value !~ /\D/ && ( $value == TRUE || $value == FALSE ) ) {
        $self->throw( "Found [" . $value
        . "] where " . TRUE . " or " . FALSE . " expected" );
    }
} # _is_true_or_false


1;
