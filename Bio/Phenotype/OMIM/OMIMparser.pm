# $Id$
#
# BioPerl module for Bio::Phenotype::OMIM::OMIMparser
#
# Cared for by Christian M. Zmasek <czmasek@gnf.org>
#
# Copyright Christian M. Zmasek
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

OMIMparser - parser for the OMIM database

=head1 SYNOPSIS

  use Bio::Phenotype::OMIM::OMIMparser;

  $omim_parser = Bio::Phenotype::OMIM::OMIMparser->new( -genemap  => "/path/to/genemap",
                                                        -omimtext => "/path/to/omim.txt" );
                                                        
  while ( my $omim_entry = $omim_parser->next() ) {
    # This prints everything.
    print( $omim_entry->to_string() );
    print "\n\n";
    
    # This gets individual data (some of them object-arrays)
    # (and illustrates the relevant methods of OMIMentry).
    my $number     = $omim_entry->MIM_number();
    my $title      = $omim_entry->title();
    my $alt        = $omim_entry->alternative_titles_and_symbols();
    my $mtt        = $omim_entry->more_than_two_genes();
    my $is_separat = $omim_entry->is_separate();
    my $desc       = $omim_entry->description();
    my $mm         = $omim_entry->mapping_method();
    my $gs         = $omim_entry->gene_status();
    my $cr         = $omim_entry->created();
    my $cont       = $omim_entry->contributors();
    my $sa         = $omim_entry->additional_references();
    my $cs         = $omim_entry->clinical_symptoms();
    my $mini_mim   = $omim_entry->miniMIM();
    my @corrs      = $omim_entry->each_Correlate();
    my @refs       = $omim_entry->each_Reference();
    my @avs        = $omim_entry->each_AllelicVariant();
    my @cps        = $omim_entry->each_CytoPosition();
    my @gss        = $omim_entry->each_gene_symbol();
    # do something ...
  }                                                      

=head1 DESCRIPTION

This parser returns Bio::Phenotype::OMIM::OMIMentry objects
(which inherit from Bio::Phenotype::PhenotypeI).
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

Email: czmasek@gnf.org  or  zmasek@yahoo.com

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


package Bio::Phenotype::OMIM::OMIMparser;

use vars qw( @ISA );
use strict;

use Bio::Root::IO;
use Bio::Root::Root;
use Bio::Species;
use Bio::Annotation::Reference;
use Bio::Map::CytoPosition;
use Bio::Phenotype::OMIM::OMIMentry;
use Bio::Phenotype::OMIM::OMIMentryAllelicVariant;
use Bio::Phenotype::Correlate;

@ISA = qw( Bio::Root::Root );


use constant DEFAULT_STATE               => 0;
use constant MIM_NUMBER_STATE            => 1;
use constant TITLE_STATE                 => 2;
use constant TEXT_STATE                  => 3;
use constant MINI_MIM_TEXT_STATE         => 4;
use constant ALLELIC_VARIANT_STATE       => 5;
use constant SEE_ALSO_STATE              => 6;
use constant REF_STATE                   => 7;
use constant SYMPT_STATE                 => 8;
use constant CONTRIBUTORS_STATE          => 9;
use constant CREATED_BY_STATE            => 10;
use constant EDITED_BY_STATE             => 11;
use constant MINI_MIM_EDITED_BY_STATE    => 12;
use constant MINI_MIM_CREATED_BY_STATE   => 13;
use constant MINI_MIM_CONTRIBUTORS_STATE => 14;
use constant TRUE                        => 1;
use constant FALSE                       => 0;



=head2 new

 Title   : new
 Usage   : $omim_parser = Bio::Phenotype::OMIM::OMIMparser->new( -genemap  => "/path/to/genemap",
                                                                 -omimtext => "/path/to/omim.txt" );                      
 Function: Creates a new OMIMparser.
 Returns : A new OMIMparser object.
 Args    : -genemap  => the genemap file name (optional)
           -omimtext => the omim text file name

=cut

sub new {
    my( $class, @args ) = @_;
    
    my $self = $class->SUPER::new( @args );

    my ( $genemap_file_name, $omimtxt_file_name ) 
    = $self->_rearrange( [ qw( GENEMAP OMIMTEXT ) ], @args );

    $self->init(); 
    
    $genemap_file_name && $self->genemap_file_name( $genemap_file_name );
    
    $omimtxt_file_name && $self->omimtxt_file_name( $omimtxt_file_name);
                         
    return $self;
}




=head2 next

 Title   : next()
 Usage   : while ( my $omim_entry = $omim_parser->next() ) {
               # do something with $omim_entry
           }    
 Function: Returns an Bio::Phenotype::OMIM::OMIMentry or
           undef once the end of the omim text file is reached.
 Returns : A Bio::Phenotype::OMIM::OMIMentry.
 Args    :

=cut

sub next  {
  
    my ( $self ) = @_;
    
    unless( defined( $self->_omimtxtFile() ) ) {
        $self->_noOmimTextFile();
    }
    
    if ( $self->_done() == TRUE ) {
        return undef;
    }

    my $fieldtag          = "";
    my $contents          = "";
    my $line              = "";
    my $state             = DEFAULT_STATE;
    my $saw_mini_min_flag = FALSE;
    my %record            = ();
    
    while( $line = ( $self->_omimtxtFile )->_readline() ) {
        if ( $line =~ /^\s*\*RECORD\*/ ) {
            if ( $self->_notFirstRecord() == TRUE ) {
                $record{ $state } = $contents;
                my $omim_entry = $self->_createOMIMentry( \%record );
                return $omim_entry;
            }
            else {
                $self->_notFirstRecord( TRUE );
            }
            
        }
        elsif ( $line =~ /^\s*\*FIELD\*\s*(\S+)/ ) {
            $fieldtag = $1;
            if ( $state != DEFAULT_STATE ) {
                if ( exists( $record{ $state } ) ) {
                    chomp( $record{ $state } );
                    $record{ $state } = $record{ $state }.$contents;
                }
                else {
                    $record{ $state } = $contents;
                }
            }
            $contents = "";
            
            if ( $fieldtag eq "NO" ) {
                $state = MIM_NUMBER_STATE;
                $saw_mini_min_flag = FALSE;   
            }
            elsif ( $fieldtag eq "TI" ) {
                $state = TITLE_STATE;
                $saw_mini_min_flag = FALSE;   
            }
            elsif ( $fieldtag eq "TX" ) {
                $state = TEXT_STATE;
                $saw_mini_min_flag = FALSE;   
            }
            elsif ( $fieldtag eq "MN" ) {
                $state = MINI_MIM_TEXT_STATE;
                $saw_mini_min_flag = TRUE;           
            }
            elsif ( $fieldtag eq "AV" ) {
                $state = ALLELIC_VARIANT_STATE;
                $saw_mini_min_flag = FALSE;     
            }
            elsif ( $fieldtag eq "SA" ) { 
                $state = SEE_ALSO_STATE;
                $saw_mini_min_flag = FALSE;   
            }
            elsif ( $fieldtag eq "RF" ) {
                $state = REF_STATE;
                $saw_mini_min_flag = FALSE;   
            }
            elsif ( $fieldtag eq "CS" ) {
                $state = SYMPT_STATE;
                $saw_mini_min_flag = FALSE;   
            }
            elsif ( $fieldtag eq "CN" ) {
                if ( $saw_mini_min_flag == TRUE ) {
                    $state = MINI_MIM_CONTRIBUTORS_STATE;
                }
                else {
                    $state = CONTRIBUTORS_STATE;
                }     
            }
            elsif ( $fieldtag eq "CD" ) {
                if ( $saw_mini_min_flag == TRUE ) {
                    $state = MINI_MIM_CREATED_BY_STATE;
                }
                else {
                    $state = CREATED_BY_STATE;
                }     
            }
            elsif ( $fieldtag eq "ED" ) {
                if ( $saw_mini_min_flag == TRUE ) {
                    $state = MINI_MIM_EDITED_BY_STATE;
                }
                else {
                    $state = EDITED_BY_STATE;
                }     
            }
            else {
                print "new tag: $fieldtag\n";
                die;
            }

        }
        else {
            $contents .= $line;
        }
    }

    $self->_omimtxtFile()->close();
    $self->_done( TRUE );

    unless( %record ) {
        $self->_notOmim();
    }

    $record{ $state } = $contents;
    my $omim_entry = $self->_createOMIMentry( \%record );
    
    return $omim_entry;

} # next




=head2 init

 Title   : init()
 Usage   : $omim_parser->init();   
 Function: Initializes this OMIMparser to all "".
 Returns : 
 Args    :

=cut

sub init {
    my ( $self ) = @_;
    
    $self->genemap_file_name( "" );
    $self->omimtxt_file_name( "" );
    $self->_genemapHash( {} );
    $self->_omimtxtFile( undef );
    $self->_notFirstRecord( FALSE );
    $self->_done( FALSE );

} # init




=head2 genemap_file_name

 Title   : genemap_file_name
 Usage   : $omimparser->genemap_file_name( "genemap" );
 Function: Set/get for the genemap file name.
 Returns : The genemap file name [string].
 Args    : The genemap file name [string] (optional).

=cut

sub genemap_file_name {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_genemap_file_name" } = $value;
        if ( $value =~ /\W/ ) {
            _genemapHash( $self->_read_genemap( $value ) );
        }
    }
    
    return $self->{ "_genemap_file_name" };
} # genemap_file_name




=head2 omimtxt_file_name

 Title   : omimtxt_file_name
 Usage   : $omimparser->omimtxt_file_name( "omim.txt" );
 Function: Set/get for the omim text file name.
 Returns : The the omim text file name [string].
 Args    : The the omim text file name [string] (optional).

=cut

sub omimtxt_file_name {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{ "_omimtxt_file_name" } = $value;
        if ( $value =~ /\W/ ) {
            $self->_omimtxtFile( new Bio::Root::IO->new( -file => $value ) );
        } 
    }
    
    return $self->{ "_omimtxt_file_name" };
} # omimtxt_file_name





sub _createOMIMentry {
    my ( $self, $record_ref ) = @_;
    
    my $omim_entry = Bio::Phenotype::OMIM::OMIMentry->new();
    my $mini_mim   = Bio::Phenotype::OMIM::MiniMIMentry->new();
    
    while ( ( my $key, my $val ) = each( %$record_ref ) ) {
        
        if ( $key == MIM_NUMBER_STATE ) {
            $val =~ s/\s+//g;
            $val =~ s/\D//g;
           
            $omim_entry->MIM_number( $val );
            
            my $gm = $self->_genemapHash();
            if ( exists( $$gm{ $val } ) ) {
                $self->_parse_genemap( $omim_entry, $val );
            }
            
        }
        elsif ( $key == TITLE_STATE ) {
            my ( $title, $alt_titles ) = $self->_parse_title( $val );
            $omim_entry->title( $title );
            $omim_entry->alternative_titles_and_symbols( $alt_titles );
            if ( $title =~ /^\*/ ) {
                 $omim_entry->is_separate( TRUE );
            }
            elsif ( $title =~ /^#/ ) {
                 $omim_entry->more_than_two_genes( TRUE );
            } 
        }
        elsif ( $key == TEXT_STATE ) {
            $omim_entry->description( $val );
        }
        elsif ( $key == ALLELIC_VARIANT_STATE ) {
            my @allelic_variants =  $self->_parse_allelic_variants( $val );
            $omim_entry->add_AllelicVariants( @allelic_variants );
        }
        elsif ( $key == SEE_ALSO_STATE ) {
            $omim_entry->additional_references( $val );
        }
        elsif ( $key == REF_STATE ) {
            my @refs =  $self->_parse_references( $val );
            $omim_entry->add_References( @refs );
        }
        elsif ( $key == SYMPT_STATE ) {
            $omim_entry->clinical_symptoms( $val );
        }
        elsif ( $key == CONTRIBUTORS_STATE ) {
            $omim_entry->contributors( $val );
        }
        elsif ( $key == CREATED_BY_STATE ) {
            $omim_entry->created( $val );
        }
        elsif ( $key == EDITED_BY_STATE ) {
            $omim_entry->edited( $val );
        }
        elsif ( $key == MINI_MIM_TEXT_STATE ) {
            $mini_mim->description( $val );
        }
        elsif ( $key == MINI_MIM_CONTRIBUTORS_STATE ) {
            $mini_mim->contributors( $val );
        }
        elsif ( $key == MINI_MIM_CREATED_BY_STATE ) {
            $mini_mim->created( $val );
        }
        elsif ( $key == MINI_MIM_EDITED_BY_STATE ) {
            $mini_mim->edited( $val );
        }
    
    }
    
    my $man = Bio::Species->new();
    $man->classification( qw( sapiens Homo ) );
    $man->common_name( "man" );
    $omim_entry->species( $man );
    $omim_entry->miniMIM( $mini_mim );
    
    return $omim_entry;

} # _createOMIMentry




sub _parse_genemap {
     my ( $self, $omim_entry, $val ) = @_;
     
     my $genemap_line = ${ $self->_genemapHash() }{ $val };
     my @a = split( /\|/, $genemap_line );

     my $locations = $a[ 4 ];
     if ( defined ( $locations ) ) {
          $locations =~ s/\s+//g;
          my @ls = split( /[,;]/, $locations );
          my @cps;
          foreach my $l ( @ls ) {
               my $cp = Bio::Map::CytoPosition->new( -value => $l );
               push( @cps, $cp ); 
          }
          $omim_entry->add_CytoPositions( @cps );
     }

     my $gene_symbols = $a[ 5 ];
     if ( defined ( $gene_symbols ) ) {
          $gene_symbols =~ s/\s+//g;
          my @gss = split( /[,;]/, $gene_symbols );
          $omim_entry->add_gene_symbols( @gss );
     }

     my $mouse_correlates = $a[ 16 ];
     if ( defined ( $mouse_correlates ) ) {
          $mouse_correlates =~ s/\s+//g;
          my @mcs = split( /[,;]/, $mouse_correlates );
          my @cs;
          foreach my $mc ( @mcs ) {
               my $mouse = Bio::Species->new();
               $mouse->classification( qw( musculus Mus ) );
               $mouse->common_name( "mouse" );
               my $c = Bio::Phenotype::Correlate->new();
               $c->name( $mc );
               $c->species( $mouse );
               $c->type( "OMIM mouse correlate" );

               push( @cs, $c ); 
          }
          $omim_entry->add_Correlates( @cs );
     }

     $omim_entry->gene_status( $a[ 6 ] ) if defined $a[ 6 ];
     $omim_entry->mapping_method( $a[ 10 ] ) if defined $a[ 10 ];
     $omim_entry->comment( $a[ 11 ] ) if defined $a[ 11 ];

} # _parse_genemap




sub _parse_allelic_variants {
    my ( $self, $text ) = @_;
    
    my @allelic_variants;
    my $number      = "";
    my $title       = "";
    my $symbol      = "";
    my $mutation    = "";
    my $description = "";
     
    my @lines = split( /\n/, $text );
    
    foreach my $line ( @lines ) {
        if ( $line !~ /\W/ ) {
            next;
        }
        elsif ( $line =~ /^(\s*\.\d+)/ ) {
            if ( $number ne "" ) {
                my $allelic_variant = $self->_create_allelic_variant( $number, $title, 
                                                    $symbol, $mutation, $description );
                
                push( @allelic_variants, $allelic_variant );
            }
            $number      = $1;
            $title       = "";
            $symbol      = "";
            $mutation    = "";
            $description = "";
        }
        elsif ( $title eq "" ) {
            $title = $line;
        }
        elsif ( $symbol eq "" ) {
            if ( $line =~ /\s*(.+)\s*,\s*(.+)/ ) {
                $symbol   = $1;
                $mutation = $2;
            }
            else {
                $symbol = $line;
            }
        }
        else {
            $description .= $line;
        }
    }
    
    my $allelic_variant = $self->_create_allelic_variant( $number, $title, 
                                        $symbol, $mutation, $description );
    
    push( @allelic_variants, $allelic_variant );
    
    return @allelic_variants;
    
} # _parse_allelic_variants




sub _create_allelic_variant {
    my ( $self, $number, $title, $symbol, $mutation, $description ) = @_;
    
    $mutation =~ /(\d+)([a-z]+)(\d+)/i;
    
    my $aa_ori   = $1;
    my $aa_mut   = $3;
    my $position = $2;
    
    my $allelic_variant = Bio::Phenotype::OMIM::OMIMentryAllelicVariant->new();
    $allelic_variant->number( $number );
    $allelic_variant->aa_ori( $aa_ori );
    $allelic_variant->aa_mut( $aa_mut );
    $allelic_variant->position( $position );
    $allelic_variant->title( $title );
    $allelic_variant->symbol( $symbol );
    $allelic_variant->description( $description );
     
    return $allelic_variant; 
    
} # _create_allelic_variant




sub _parse_title {
    my ( $self, $text ) = @_;
     
    $text =~ /^(.+)\n/;
    my $title = $1;
    $text =~ s/^(.+)\n//;
    $text =~ s/;+/;/g;
    #$text =~ s/\n+/ /g;
    $text =~ s/^;//;
    $text =~ s/;$//;
    
    return ( $title, $text ); 
} # _parse_title




sub _parse_references {
    my ( $self, $text ) = @_;
    
    $text =~ s/\A\s+//;
    $text =~ s/\s+\z//;
    $text =~ s/\A\d+\.\s*//;
    
    my @references;
    
    my @texts = split( /\s*\n\s*\n\s*\d+\.\s*/, $text );
    
    my $authors   = "";
    my $title     = "";
    my $location  = "";
    
    foreach my $t ( @texts ) {
    
        $authors   = "";
        $title     = "";
        $location  = "";
        
        $t =~ s/\s+/ /g;
        
        if ( $t =~ /(.+?)\s*:\s*(.+?[.?!])\s+(.+?)\s+(\S+?)\s*:\s*(\w?\d+.*)\s*,\s*(\d+)/ ) {
            $authors    = $1;
            $title      = $2;
            my $journal = $3;
            my $volume  = $4;
            my $fromto  = $5;
            my $year    = $6;
            my $from    = "",
            my $to      = "";
            if ( $fromto =~ /(\d+)-+(\d+)/ ) {
                $from = $1;
                $to   = "-".$2;
            }
            elsif ( $fromto =~ /\A(\w+)/ ) {
                $from = $1;
            }
            $location = $journal." ".$volume." ".$from.$to." (".$year.")";
        }
       
            
        elsif ( $t =~ /(.+?)\s*:\s*(.+?[.?!])\s*(.+?)\z/ ) {
            $authors   = $1;
            $title     = $2;
            $location  = $3;
        }
        else {
            $title = $t;  
        }
        
        $authors =~ s/;//g;
        
        $title =~ s/\.\z//;
         
        my $ref = Bio::Annotation::Reference->new( -title    => $title,
                                                   -location => $location,
                                                   -authors  => $authors );
         
        push( @references, $ref );
       
    }
    return @references;
    
} # _parse_references




sub _genemapHash {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        unless ( ref( $value ) eq "HASH" ) {
            $self->throw( "Argument to method \"_genemapHash\" is not a reference to an Hash" );
        }
        $self->{ "_genemap_hash" } = $value;
     
    }
    
    return $self->{ "_genemap_hash" };
} # _genemapHash




sub _notFirstRecord {

    my ( $self, $value ) = @_;

    if ( defined $value ) {
        unless ( $value == FALSE || $value == TRUE ) {
            $self->throw( "Argument to method \"_notFirstRecord\" must be either ".TRUE." or ".FALSE );
        }
        $self->{ "_not_first_record" } = $value;
    }
    
    return $self->{ "_not_first_record" };
} # _notFirstRecord




sub _done {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        unless ( $value == FALSE || $value == TRUE ) {
            $self->throw( "Argument to method \"_done\" must be either ".TRUE." or ".FALSE );
        }
        $self->{ "_done" } = $value;
    }
    
    return $self->{ "_done" };
} # _done




sub _omimtxtFile {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        unless ( $value->isa( "Bio::Root::IO" ) ) {
            $self->throw( "Argument to method \"_omimtxtFile\" is not a valid \"Bio::Root::IO\"" );
        }
        $self->{ "_omimtxt_file" } = $value;
     
    }
    
    return $self->{ "_omimtxt_file" };
} # _omimtxtFile




sub _read_genemap {
    my ( $self, $genemap_file_name ) = @_;
    
    my $line         = "";
    my %genemap_hash = ();
    my $genemap_file = new Bio::Root::IO->new( -file => $genemap_file_name );
    my @a            = ();
    my %gm           = ();
    
    while( $line = $genemap_file->_readline() ) {
        @a = split( /\|/, $line );
        unless( scalar( @a ) == 18 ) {
            $self->throw( "Gene map file \"".$self->genemapFileName()."\" is not in the expected format" );
        }
        $gm{ $a[ 9 ] } = $line;
    }
    $genemap_file->close();
    $self->_genemapHash( \%gm );
  
} #_read_genemap 




sub _noOmimTextFile {
    my ( $self ) = @_;

    my $msg =  "Need to indicate a OMIM text file to read from with\n";
    $msg .= "either \"OMIMparser->new( -omimtext => \"path/to/omim.txt\" );\"\n";
    $msg .= "or \"\$omim_parser->omimtxtFileName( \"path/to/omim.txt\" );\"";
    $self->throw( $msg );
} # _noOmimTextFile




sub _notOmim {
    my ( $self ) = @_;

    my $msg =  "File \"".$self->omimtxtFileName()."\" appears not to be a OMIM text file";
    $self->throw( $msg );
} # _notOmim

1;
