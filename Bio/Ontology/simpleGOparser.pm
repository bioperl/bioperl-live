# $Id Exp $
#
# BioPerl module for Bio::Ontology::simpleGOparser
#
# Cared for by Christian M. Zmasek <czmasek@gnf.org> or <cmzmasek@hotmail.com>
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

simpleGOparser - a simple GO parser returning a simpleGOengine

=head1 SYNOPSIS

use Bio::Ontology::simpleGOparser;



my $parser = Bio::Ontology::simpleGOparser->new( -go_defs_file_name    => "/home/czmasek/GO/GO.defs",
                                                 -components_file_name => "/home/czmasek/GO/component.ontology",
                                                 -functions_file_name  => "/home/czmasek/GO/function.ontology",
                                                 -processes_file_name  => "/home/czmasek/GO/process.ontology" );


my $engine = $parser->parse();

my $IS_A    = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
my $PART_OF = Bio::Ontology::RelationshipType->get_instance( "PART_OF" );

                       

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

report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR

Christian M. Zmasek

Email: czmasek@gnf.org  or  cmzmasek@hotmail.com

WWW:   http://www.genetics.wustl.edu/eddy/people/zmasek/

Address: 

  Genomics Institute of the Novartis Research Foundation
  10675 John Jay Hopkins Drive
  San Diego, CA 92121

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package  Bio::Ontology::simpleGOparser;

use Graph;

use vars qw( @ISA );
use strict;

use Bio::Root::IO;
use Bio::Root::Root;
use Bio::Ontology::GOterm;
use Bio::Ontology::simpleGOengine;

use constant TRUE         => 1;
use constant FALSE        => 0;


@ISA = qw( Bio::Root::Root );





=head2 new

 Title   : new
 Usage   : $parser = Bio::Ontology::simpleGOparser->new( -go_defs_file_name    => "/path/to/GO.defs",
                                                         -components_file_name => "/path/to/component.ontology"
                                                         -functions_file_name  => "/path/to/function.ontology"
                                                         -processes_file_name  => "/path/to/process.ontology" );                      
 Function: Creates a new simpleGOparser.
 Returns : A new simpleGOparser object.
 Args    : -go_defs_file_name    => the GO defs-file name
           -components_file_name => the component.ontology-file name
           -functions_file_name  => the function.ontology-file name
           -processes_file_name  => the process.ontology-file name

=cut

sub new {
    my( $class, @args ) = @_;
    
    my $self = $class->SUPER::new( @args );

    my ( $go_defs_file_name,
         $components,
         $functions,
         $processes, ) 
    = $self->_rearrange( [ qw( GO_DEFS_FILE_NAME
                               COMPONENTS_FILE_NAME
                               FUNCTIONS_FILE_NAME
                               PROCESSES_FILE_NAME ) ], @args );

    $self->init(); 
    
    $go_defs_file_name && $self->go_defs_file_name( $go_defs_file_name );
    $components        && $self->components_file_name( $components );
    $functions         && $self->functions_file_name( $functions );
    $processes         && $self->processes_file_name( $processes );
    
    
                         
    return $self;
} # new



=head2 init

 Title   : init()
 Usage   : $parser->init();   
 Function: Initializes this object.
 Returns : 
 Args    :

=cut

sub init {
    my ( $self ) = @_;
    
    $self->{ "_go_defs_file_name" }    = undef;
    $self->{ "_components_file_name" } = undef;
    $self->{ "_functions_file_name" }  = undef;
    $self->{ "_processes_file_name" }  = undef;
    $self->_done( FALSE );
    $self->_not_first_record( FALSE );
    $self->_term( "" );
    
    $self->_go_engine( Bio::Ontology::simpleGOengine->new() );

} # init



=head2 parse

 Title   : parse()
 Usage   : $parser->parse();   
 Function: Parses the files set wirh "new" or with methods
           go_defs_file_name, components_file_name, functions_file_name,
           processes_file_name.
 Returns : [Bio::Ontology::simpleGOengine]
 Args    :

=cut

sub parse {
    my ( $self ) = @_;
    my $x = 0;
    
    while( my $goterm = $self->_next_term() ) {
        $self->_add_term( $goterm );
    }

    if ( $self->components_file_name() ) {
        $self->_parse_relationships_file( $self->components_file_name(), "components ontology" );
    }
    if ( $self->functions_file_name() ) {
        $self->_parse_relationships_file( $self->functions_file_name(), "functions ontology" );
    }
    if ( $self->processes_file_name() ) {
        $self->_parse_relationships_file( $self->processes_file_name(), "processes ontology" );
    }
    
    return $self->_go_engine();
    
} # parse



=head2 go_defs_file_name

 Title   : go_defs_file_name
 Usage   : $parser->go_defs_file_name( "GO.defs" );
 Function: Set/get for the GO defs-file_name.
 Returns : The GO defs-file_name [string].
 Args    : The GO defs-file_name [string] (optional).

=cut

sub go_defs_file_name {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_go_defs_file_name" } = $value;
        if ( $value =~ /\W/ ) {
            $self->_go_defs_file( new Bio::Root::IO->new( -file => $value ) );
        } 
    }
    
    return $self->{ "_go_defs_file_name" };
} # go_defs_file_name



=head2 components_file_name

 Title   : components_file_name
 Usage   : $parser-> components_file_name( "function.ontology" );
 Function: Set/get for the function ontology file name.
 Returns : The function ontology file name [string].
 Args    : The function ontology file name [string] (optional).

=cut

sub components_file_name {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_components_file_name" } = $value;
    }
    
    return $self->{ "_components_file_name" };
} # components_file_name



=head2 functions_file_name

 Title   : functions_file_name
 Usage   : $parser->functions_file_name( "function.ontology" );
 Function: Set/get for functions file name.
 Returns : The functions file name [string].
 Args    : The functions file name [string] (optional).

=cut

sub functions_file_name {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_functions_file_name" } = $value;
    }
    
    return $self->{ "_functions_file_name" };
} # functions_file_name



=head2 processes_file_name

 Title   : processes_file_name
 Usage   : $parser->processes_file_name( "GO.defs" );
 Function: Set/get for the processes file name.
 Returns : The processes file name [string].
 Args    : The processes file name [string] (optional).

=cut

sub processes_file_name {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_processes_file_name" } = $value;
    }
    
    return $self->{ "_processes_file_name" };
} # processes_file_name










# INTERNAL METHODS
# ----------------


# This simply delegates. See simpleGOengine.
sub _add_term {
    my ( $self, $term ) = @_;

    $self->_go_engine()->add_term( $term );


} # _add_term 



# This simply delegates. See simpleGOengine
sub _part_of_relationship {
    my ( $self, $term ) = @_;

    return $self->_go_engine()->part_of_relationship();


} # _part_of_relationship 



# This simply delegates. See simpleGOengine
sub _is_a_relationship {
    my ( $self, $term ) = @_;

    return $self->_go_engine()->is_a_relationship();


} # _is_a_relationship 



# This simply delegates. See simpleGOengine
sub _add_relationship {
    my ( $self, $parent, $child, $type ) = @_;

   
    $self->_go_engine()->add_relationship( $parent, $child, $type );


} # _add_term 


# This simply delegates. See simpleGOengine
sub _has_term {
    my ( $self, $term ) = @_;

    
    return $self->_go_engine()->has_term( $term );


} # _add_term



# This parses the relationships files
# Arguments: filename, category-name (ie: function ontology)
sub _parse_relationships_file {
    my ( $self, $file_name, $category ) = @_;
    
    my $file = new Bio::Root::IO->new( -file => $file_name );
    my $cat = Bio::Ontology::Term->new( -name => $category );
    
    my @stack       = ();
    my $prev_spaces = -1;
    my $prev_term   = "";
    
    
    while( my $line = $file->_readline() ) {
        
        if ( $line =~ /^!/ ) {
            next;
        }
        
        my $current_term   = $self->_get_first_goid( $line );
        my @isa_parents    = $self->_get_isa_goids( $line );
        my @partof_parents = $self->_get_partof_goids( $line );
        my @syns           = $self->_get_synonyms( $line );
        my @sec_go_ids     = $self->_get_secondary_goids( $line );
        my @cross_refs     = $self->_get_db_cross_refs( $line );
        
        
        if ( ! $self->_has_term( $current_term ) ) {
            my $goterm = $self->_create_GOentry( $self->_get_name( $line, $current_term ), $current_term );
            $self->_add_term( $goterm );
        }
        
        my $current_term_object = $self->_go_engine()->get_term( $current_term );
        
        $current_term_object->add_dblinks( @cross_refs );
        $current_term_object->add_secondary_GO_ids( @sec_go_ids );
        $current_term_object->add_synonyms( @syns );
        unless ( $line =~ /^\$/ ) {
            $current_term_object->category( $cat );
        }
        foreach my $parent ( @isa_parents ) {
            if ( ! $self->_has_term( $parent ) ) {
                my $goterm = $self->_create_GOentry( $self->_get_name( $line, $parent ), $parent );
                $self->_add_term( $goterm );
            }
            
            $self->_add_relationship( $parent,
                                      $current_term,
                                      $self->_is_a_relationship() );
             
        }
        foreach my $parent ( @partof_parents ) {
            if ( ! $self->_has_term( $parent ) ) {
                my $goterm = $self->_create_GOentry( $self->_get_name( $line, $parent ), $parent );
                $self->_add_term( $goterm );
            }
           
            $self->_add_relationship( $parent,
                                      $current_term,
                                      $self->_part_of_relationship() );
        }
        
        my $current_spaces = $self->_count_spaces( $line );
        
        if ( $current_spaces != $prev_spaces  ) {
       
            if ( $current_spaces == $prev_spaces + 1 ) {
                push( @stack, $prev_term ); 
            }
            elsif ( $current_spaces < $prev_spaces ) {
                my $n = $prev_spaces -  $current_spaces;
                for ( my $i = 0; $i < $n; ++$i ) {
                    pop( @stack );
                }
            }
            else {
                die( "format error" );
            } 
        }
        
        my $parent = $stack[ @stack - 1 ];
        
        
        if ( $line =~ /^\$/ ) {
        }
        elsif ( $line =~ /^\s*</ ) {
            $self->_add_relationship( $parent,
                                      $current_term,
                                      $self->_part_of_relationship() );
        }
        elsif ( $line =~ /^\s*%/ ) {
            $self->_add_relationship( $parent,
                                      $current_term,
                                      $self->_is_a_relationship() );
        }
        else {
            die( "format error" );
        }
        
        
        $prev_spaces = $current_spaces;
        
        $prev_term = $current_term;
        
    } 
} # _parse_relationships_file



# Parses the 1st GO id number out of line.
sub _get_first_goid {
    my ( $self, $line ) = @_;
    
    if ( $line =~ /;\s*(GO:\d{7})/ ) {
        return $1;
    }
    else {
        die( "format error" );
    }
    
} # _get_first_goid



# Parses the name out of line.
sub _get_name {
    my ( $self, $line, $goid ) = @_;
    
    if ( $line =~ /([^;^<^%^,]+);\s*$goid/ ) {
        my $name = $1;
        $name =~ s/\s+$//;
        $name =~ s/^\s+//;
        return $name;
    }
    else {
        return undef;
    }
} # _get_name   


# Parses the synonyms out of line.
sub _get_synonyms {
    my ( $self, $line ) = @_;
   
    my @synonyms = ();
   
    while ( $line =~ /synonym\s*:\s*([^;^<^%]+)/g ) {
        my $syn = $1;
        $syn =~ s/\s+$//;
        $syn =~ s/^\s+//;
        push( @synonyms, $syn );
    }
    return @synonyms;
    
} # _get_synonyms



# Parses the db cross refs out of line.
sub _get_db_cross_refs {
    my ( $self, $line ) = @_;
   
    my @refs = ();
   
    while ( $line =~ /;([^;^<^%^:]+:[^;^<^%^:]+)/g ) {
        my $ref = $1;
        if ( $ref =~ /synonym/ || $ref =~ /GO:\d{7}/ ) {
            next;
        }
        $ref =~ s/\s+$//;
        $ref =~ s/^\s+//;
        push( @refs, $ref );
    }
    return @refs;
    
} 


# Parses the secondary go ids out of a line
sub _get_secondary_goids {
    my ( $self, $line ) = @_;
    my @secs = ();
   
    while ( $line =~ /,\s*(GO:\d{7})/g ) {
        my $sec = $1;
        push( @secs, $sec );
    }
    return @secs;
    
} # _get_secondary_goids 



# Parses the is a ids out of a line
sub _get_isa_goids {
    my ( $self, $line ) = @_;
    
    my @ids = ();
    
    $line =~ s/GO:\d{7}//;
    
    while ( $line =~ /%[^<^,]*?(GO:\d{7})/g ) {
        push( @ids, $1 );
    }
    return @ids; 
} # _get_isa_goids



# Parses the part of ids out of a line
sub _get_partof_goids {
    my ( $self, $line ) = @_;
    
    my @ids = ();
    
    $line =~ s/GO:\d{7}//;
    
    while ( $line =~ /<[^%^,]*?(GO:\d{7})/g ) {
        push( @ids, $1 );
    }
    return @ids; 
    
    
} # _get_partof_goids




# Counts the spaces at the beginning of a line in the relationships files
sub _count_spaces {
    my ( $self, $line ) = @_;
     
    if ( $line =~ /^([ ]+)/ ) {
         return length( $1 );
    }
    else {
         return 0;
    }
} # _count_spaces




# "next" method for parsing the defintions file
sub _next_term {
    my ( $self ) = @_;

    if ( $self->_done() == TRUE ) {
        return undef;
    }
    
    my $line      = "";
    my $goid      = "";
    my $next_term = "";
    my $def       = "";
    my $comment   = "";
    my @def_refs  = ();
    
    while( $line = ( $self->_go_defs_file )->_readline() ) {
    
        if ( $line !~ /\S/ 
        ||   $line =~ /^\s*!/ ) {
            next;
        }
        
        elsif ( $line =~ /^\s*term:\s*(.+)/ ) {
            $next_term = $1;
            if ( $self->_not_first_record() == TRUE ) {
                my $entry = $self->_create_GOentry( $self->_term(), $goid, $def, $comment, \@def_refs );
                $self->_term( $next_term );
                return $entry;
            }
            else {
                $self->_term( $next_term );
                $self->_not_first_record( TRUE );
            }
        }
        elsif ( $line =~ /^\s*goid:\s*(.+)/ ) {
            $goid = $1;
        }
        elsif ( $line =~ /^\s*definition:\s*(.+)/ ) {
            $def = $1;   
        }
        elsif ( $line =~ /^\s*definition_reference:\s*(.+)/ ) {
            push( @def_refs, $1 );  
        }
        elsif ( $line =~ /^\s*comment:\s*(.+)/ ) {
            $comment = $1;  
        }
    }
    $self->_done( TRUE );
    return $self->_create_GOentry( $self->_term(), $goid, $def, $comment, \@def_refs );
} # _next_term





# Holds the GO engine to be parsed into
sub _go_engine {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_go_engine" } = $value;
    }
    
    return $self->{ "_go_engine" };
} # _go_enginee




# Used to create GO terms.
# Arguments: name, id
sub _create_GOentry {
    my ( $self, $name, $goid ) = @_;

    my $term = Bio::Ontology::GOterm->new();

    $term->GO_id( $goid );
    $term->name( $name );
    
    return $term;

} # _create_GOentry



# Holds whether first record or not
sub _not_first_record {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        unless ( $value == FALSE || $value == TRUE ) {
            $self->throw( "Argument to method \"_not_first_record\" must be either ".TRUE." or ".FALSE );
        }
        $self->{ "_not_first_record" } = $value;
    }
    
    return $self->{ "_not_first_record" };
} # _not_first_record



# Holds whether done or not
sub _done {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        unless ( $value == FALSE || $value == TRUE ) {
            $self->throw( "Found [$value] where [" . TRUE
            ." or " . FALSE . "] expected" );
        }
        $self->{ "_done" } = $value;
    }
    
    return $self->{ "_done" };
} # _done



# Holds the GO definitions file
sub _go_defs_file {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        unless ( $value->isa( "Bio::Root::IO" ) ) {
            $self->throw( "Argument to method \"_go_defs_File\" is not a valid \"Bio::Root::IO\"" );
        }
        $self->{ "_go_defs_file" } = $value;
     
    }
    
    return $self->{ "_go_defs_file" };
} # _go_defs_file



# Holds a term.  
sub _term {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_term" } = $value;
    }
    
    return $self->{ "_term" };
} # _term  
  
  
