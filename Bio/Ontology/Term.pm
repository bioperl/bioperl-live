# $Id$
#
# BioPerl module for Bio::Ontology::Term
#
# Cared for by Christian M. Zmasek <czmasek@gnf.org> or <zmasek@yahoo.com>
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

Term - interface for ontology terms

=head1 SYNOPSIS

#get Bio::Ontology::TermI somehow.
  
  print $term->identifier(), "\n";
  print $term->name(), "\n";
  print $term->definition(), "\n";
  print $term->is_obsolete(), "\n";
  print $term->comment(), "\n";
  
  foreach my $synonym ( $term->each_synonym() ) {
      print $synonym, "\n";
  }
    

=head1 DESCRIPTION

This is "dumb" interface for ontology terms providing basic methods
(it provides no functionality related to graphs).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org             - General discussion
  http://bio.perl.org/MailList.html - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR

Christian M. Zmasek

Email: czmasek@gnf.org  or  zmasek@yahoo.com

WWW:   http://www.genetics.wustl.edu/eddy/people/zmasek/

Address: 

  Genomics Institute of the Novartis Research Foundation
  10675 John Jay Hopkins Drive
  San Diego, CA 92121

=head1 APPENDIX

The rest of the documentation details each of the object
methods.

=cut


# Let the code begin...

package Bio::Ontology::Term;
use vars qw( @ISA );
use strict;
use Bio::Root::Object;
use Bio::Ontology::TermI;

use constant TRUE    => 1;
use constant FALSE   => 0;

@ISA = qw( Bio::Root::Root Bio::Ontology::TermI );



=head2 new

 Title   : new
 Usage   : $term = Bio::Ontology::Term->new( -identifier  => "16847",
                                             -name        => "1-aminocyclopropane-1-carboxylate synthase",
                                             -definition  => "Catalysis of ...",
                                             -is_obsolete => 0,
                                             -comment     => "" );                   
 Function: Creates a new Bio::Ontology::Term.
 Returns : A new Bio::Ontology::Term object.
 Args    : -identifier            => the identifier of this term [scalar]
           -name                  => the name of this term [scalar]
           -definition            => the definition of this term [scalar]  
           -category              => a relationship between this Term and another Term [TermI or scalar]
           -version               => version information [scalar]
           -is_obsolete           => the obsoleteness of this term [0 or 1]   
           -comment               => a comment [scalar]

=cut

sub new {

    my( $class,@args ) = @_;
    
    my $self = $class->SUPER::new( @args );
   
    my ( $identifier,
         $name,
         $definition,
         $category,
         $version,     
         $is_obsolete,       
         $comment )
	= $self->_rearrange( [ qw( IDENTIFIER
				   NAME
				   DEFINITION
				   CATEGORY 
				   VERSION    
				   IS_OBSOLETE      
				   COMMENT ) ], @args );
   
    $self->init(); 
    
    $identifier            && $self->identifier( $identifier );
    $name                  && $self->name( $name );
    $definition            && $self->definition( $definition );
    $category              && $self->category( $category );   
    $version               && $self->version( $version );   
    $is_obsolete           && $self->is_obsolete( $is_obsolete );      
    $comment               && $self->comment( $comment  ); 
  
                                                    
    return $self;
    
} # new



sub init {

    my( $self ) = @_;

    $self->identifier( "" );
    $self->name( "" );
    $self->definition( "" );
    $self->version( "" );
    $self->is_obsolete( FALSE );
    $self->comment( "" );
    $self->remove_synonyms();
  
} # init



=head2 identifier

 Title   : identifier
 Usage   : $term->identifier( "0003947" );
           or
           print $term->identifier();
 Function: Set/get for the identifier of this Term.
 Returns : The identifier [scalar].
 Args    : The identifier [scalar] (optional).

=cut

sub identifier {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_identifier" } = $value ? $value : undef; # no empty string
    }

    return $self->{ "_identifier" };

} # identifier




=head2 name

 Title   : name
 Usage   : $term->name( "N-acetylgalactosaminyltransferase" );
           or
           print $term->name();
 Function: Set/get for the name of this Term.
 Returns : The name [scalar].
 Args    : The name [scalar] (optional).

=cut

sub name {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_name" } = $value ? $value : undef; # no empty string
    }

    return $self->{ "_name" };

} # name





=head2 definition

 Title   : definition
 Usage   : $term->definition( "Catalysis of ..." );
           or
           print $term->definition();
 Function: Set/get for the definition of this Term.
 Returns : The definition [scalar].
 Args    : The definition [scalar] (optional).

=cut

sub definition {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_definition" } = $value ? $value : undef; # no empty string
    }

    return $self->{ "_definition" };

} # definition



=head2 category

 Title   : category
 Usage   : $term->category( $top );
           or 
           $top = $term->category();
 Function: Set/get for a relationship between this Term and
           another Term (e.g. the top level of the ontology).
 Returns : The category of this Term [TermI].
 Args    : The category of this Term [TermI or scalar -- which
           becomes the name of the catagory term] (optional).

=cut

sub category {
     my ( $self, $value ) = @_;
    
    if ( defined $value ) {
	if(! $value) {
            # no empty string here
	    $self->{"_category"} = undef;
	}
	elsif ( ! ref( $value ) ) {
            my $term = $self->new();
            $term->name( $value );
            $self->{ "_category" } = $term; 
        }
        elsif ( $value->isa( "Bio::Ontology::TermI" ) ) {
            $self->{ "_category" } = $value; 
        } 
        else {
            $self->throw( "Found [". ref( $value ) 
            . "] where [Bio::Ontology::TermI] or [scalar] expected" );
        }
    }
    
    return $self->{ "_category" };
    
} # category



=head2 version

 Title   : version
 Usage   : $term->version( "1.00" );
           or 
           print $term->version();
 Function: Set/get for version information.
 Returns : The version [scalar].
 Args    : The version [scalar] (optional).

=cut

sub version {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_version" } = $value ? $value : undef; # no empty string
    }

    return $self->{ "_version" };
    
} # version



=head2 is_obsolete

 Title   : is_obsolete
 Usage   : $term->is_obsolete( 1 );
           or
           if ( $term->is_obsolete() )
 Function: Set/get for the obsoleteness of this Term.
 Returns : the obsoleteness [0 or 1].
 Args    : the obsoleteness [0 or 1] (optional).

=cut

sub is_obsolete {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->_is_true_or_false( $value );
        $self->{ "_is_obsolete" } = $value;
    }

    return $self->{ "_is_obsolete" };

} # is_obsolete





=head2 comment

 Title   : comment
 Usage   : $term->comment( "Consider the term ..." );
           or 
           print $term->comment();
 Function: Set/get for an arbitrary comment about this Term.
 Returns : A comment.
 Args    : A comment (optional).

=cut

sub comment {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_comment" } = $value ? $value : undef; # no empty string
    }
   
    return $self->{ "_comment" };
    
} # comment




=head2 each_synonym

 Title   : each_synonym()
 Usage   : @aliases = $term->each_synonym();                 
 Function: Returns a list of aliases of this Term.
 Returns : A list of aliases [array of [scalar]].
 Args    :

=cut

sub each_synonym {
    my ( $self ) = @_;
    
    if ( $self->{ "_synonyms" } ) {
        return @{ $self->{ "_synonyms" } };
    }
    else {
        return my @a = (); 
    }
    
} # each_synonym


=head2 add_synonyms

 Title   : add_synonyms
 Usage   : $term->add_synonyms( @asynonyms );
           or
           $term->add_synonyms( $synonym );                  
 Function: Pushes one or more synonyms into the list of synonyms.
 Returns : 
 Args    : One synonym [scalar] or a list of synonyms [array of [scalar]].

=cut

sub add_synonyms {
    my ( $self, @values ) = @_;
    
    return unless( @values );
        
    push( @{ $self->{ "_synonyms" } }, @values );
    
} # add_synonyms


=head2 remove_synonyms

 Title   : remove_synonyms()
 Usage   : $term->remove_synonyms();
 Function: Deletes (and returns) the synonyms of this Term.
 Returns : A list of synonyms [array of [scalar]].
 Args    :

=cut

sub remove_synonyms {
    my ( $self ) = @_;
     
    my @a = $self->each_synonym();
    $self->{ "_synonyms" } = [];
    return @a;

} # remove_synonyms




# Title   :_is_true_or_false
# Function: Checks whether the argument is TRUE or FALSE.
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
