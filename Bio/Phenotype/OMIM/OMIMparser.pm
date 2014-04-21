#
# BioPerl module for Bio::Phenotype::OMIM::OMIMparser
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

Bio::Phenotype::OMIM::OMIMparser - parser for the OMIM database

=head1 SYNOPSIS

  use Bio::Phenotype::OMIM::OMIMparser;

  # The OMIM database is available as textfile at:
  # ftp://ncbi.nlm.nih.gov/repository/OMIM/omim.txt.Z
  # The genemap is available as textfile at:
  # ftp://ncbi.nlm.nih.gov/repository/OMIM/genemap

  $omim_parser = Bio::Phenotype::OMIM::OMIMparser->new( -genemap  => "/path/to/genemap",
                                                        -omimtext => "/path/to/omim.txt" );

  while ( my $omim_entry = $omim_parser->next_phenotype() ) {
    # This prints everything.
    print( $omim_entry->to_string() );
    print "\n\n";

    # This gets individual data (some of them object-arrays)
    # (and illustrates the relevant methods of OMIMentry).
    my $numb  = $omim_entry->MIM_number();                     # *FIELD* NO
    my $title = $omim_entry->title();                          # *FIELD* TI - first line
    my $alt   = $omim_entry->alternative_titles_and_symbols(); # *FIELD* TI - additional lines
    my $mtt   = $omim_entry->more_than_two_genes();            # "#" before title
    my $sep   = $omim_entry->is_separate();                    # "*" before title
    my $desc  = $omim_entry->description();                    # *FIELD* TX
    my $mm    = $omim_entry->mapping_method();                 # from genemap
    my $gs    = $omim_entry->gene_status();                    # from genemap
    my $cr    = $omim_entry->created();                        # *FIELD* CD
    my $cont  = $omim_entry->contributors();                   # *FIELD* CN
    my $ed    = $omim_entry->edited();                         # *FIELD* ED
    my $sa    = $omim_entry->additional_references();          # *FIELD* SA
    my $cs    = $omim_entry->clinical_symptoms_raw();              # *FIELD* CS
    my $comm  = $omim_entry->comment();                        # from genemap

    my $mini_mim   = $omim_entry->miniMIM();                   # *FIELD* MN
      # A Bio::Phenotype::OMIM::MiniMIMentry object.
      # class Bio::Phenotype::OMIM::MiniMIMentry
      # provides the following:
      # - description()
      # - created()
      # - contributors()
      # - edited() 
      #
    # Prints the contents of the MINI MIM entry (most OMIM entries do
    # not have MINI MIM entries, though).
    print $mini_mim->description()."\n";
    print $mini_mim->created()."\n";
    print $mini_mim->contributors()."\n";
    print $mini_mim->edited()."\n";

    my @corrs      = $omim_entry->each_Correlate();            # from genemap
      # Array of Bio::Phenotype::Correlate objects.
      # class Bio::Phenotype::Correlate
      # provides the following:
      # - name()
      # - description() (not used)
      # - species() (always mouse)
      # - type() ("OMIM mouse correlate")
      # - comment() 

    my @refs       = $omim_entry->each_Reference();            # *FIELD* RF
      # Array of Bio::Annotation::Reference objects.


    my @avs        = $omim_entry->each_AllelicVariant();       # *FIELD* AV
      # Array of Bio::Phenotype::OMIM::OMIMentryAllelicVariant objects.
      # class Bio::Phenotype::OMIM::OMIMentryAllelicVariant
      # provides the following:
      # - number (e.g ".0001" )
      # - title (e.g "ALCOHOL INTOLERANCE" )
      # - symbol (e.g "ALDH2*2" )
      # - description (e.g "The ALDH2*2-encoded protein has a change ..." )
      # - aa_ori  (used if information in the form "LYS123ARG" is found)
      # - aa_mut (used if information in the form "LYS123ARG" is found)
      # - position (used if information in the form "LYS123ARG" is found)
      # - additional_mutations (used for e.g. "1-BP DEL, 911T")

    my @cps        = $omim_entry->each_CytoPosition();         # from genemap
      # Array of Bio::Map::CytoPosition objects.

    my @gss        = $omim_entry->each_gene_symbol();          # from genemap
      # Array of strings.

    # do something ...
  }

=head1 DESCRIPTION

This parser returns Bio::Phenotype::OMIM::OMIMentry objects
(which inherit from Bio::Phenotype::PhenotypeI).
It parses the OMIM database available as 
ftp://ncbi.nlm.nih.gov/repository/OMIM/omim.txt.Z 
together with (optionally) the gene map file at
ftp://ncbi.nlm.nih.gov/repository/OMIM/genemap.


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


package Bio::Phenotype::OMIM::OMIMparser;

use strict;

use Bio::Root::IO;
use Bio::Species;
use Bio::Annotation::Reference;
use Bio::Map::CytoPosition;
use Bio::Phenotype::OMIM::OMIMentry;
use Bio::Phenotype::OMIM::OMIMentryAllelicVariant;
use Bio::Phenotype::Correlate;

use base qw(Bio::Root::Root);


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




=head2 next_phenotype

 Title   : next_phenotype()
 Usage   : while ( my $omim_entry = $omim_parser->next_phenotype() ) {
               # do something with $omim_entry
           }    
 Function: Returns an Bio::Phenotype::OMIM::OMIMentry or
           undef once the end of the omim text file is reached.
 Returns : A Bio::Phenotype::OMIM::OMIMentry.
 Args    :

=cut

sub next_phenotype {
    my ( $self ) = @_;
    
    unless( defined( $self->_OMIM_text_file() ) ) {
        $self->_no_OMIM_text_file_provided_error();
    }
    
    if ( $self->_done() == TRUE ) {
        return;
    }

    my $fieldtag          = "";
    my $contents          = "";
    my $line              = "";
    my $state             = DEFAULT_STATE;
    my $saw_mini_min_flag = FALSE;
    my %record            = ();
    
    while( $line = ( $self->_OMIM_text_file )->_readline() ) {
        if ( $line =~ /^\s*\*RECORD\*/ ) {
            if ( $self->_is_not_first_record() == TRUE ) {
                $self->_add_to_hash( $state, $contents,\%record );
                my $omim_entry = $self->_createOMIMentry( \%record );
                return $omim_entry;
            }
            else {
                $self->_is_not_first_record( TRUE );
            }
            
        }
        elsif ( $line =~ /^\s*\*FIELD\*\s*(\S+)/ ) {
            $fieldtag = $1;
            if ( $state != DEFAULT_STATE ) {
                $self->_add_to_hash( $state, $contents,\%record );
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
                print "Warning: Unknown tag: $fieldtag\n";
            }

        }
        else {
            $contents .= $line;
        }
    }

    $self->_OMIM_text_file()->close();
    $self->_done( TRUE );

    unless( %record ) {
        $self->_not_a_OMIM_text_file_error();
    }

    $self->_add_to_hash( $state, $contents,\%record );
    
    my $omim_entry = $self->_createOMIMentry( \%record );
    
    return $omim_entry;

} # next_phenotype




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
    $self->_genemap_hash( {} );
    $self->_OMIM_text_file( undef );
    $self->_is_not_first_record( FALSE );
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
        $self->_genemap_hash( $self->_read_genemap( $value ) );
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
            $self->_OMIM_text_file( Bio::Root::IO->new->new( -file => $value ) );
        } 
    }
    
    return $self->{ "_omimtxt_file_name" };
} # omimtxt_file_name





sub _createOMIMentry {
    my ( $self, $record_ref ) = @_;
    
    my $omim_entry = Bio::Phenotype::OMIM::OMIMentry->new();
    my $mini_mim   = Bio::Phenotype::OMIM::MiniMIMentry->new();
    
    while ( ( my $key, my $val ) = each( %$record_ref ) ) {
        
        $val =~ s/^\s+//;
        $val =~ s/\s+$//;
        
        if ( $key == MIM_NUMBER_STATE ) {
            $val =~ s/\s+//g;
            $val =~ s/\D//g;
           
            $omim_entry->MIM_number( $val );
            
            my $gm = $self->_genemap_hash();
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
            $val = undef if($val =~ /DESCRIPTION1\nDESCRIPTION2/);
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
            $val = '' if($val eq 'clinical symptoms');
            $omim_entry->clinical_symptoms_raw( $val );
        }
        elsif ( $key == CONTRIBUTORS_STATE ) {
            $val = undef if($val =~ /cn1\ncn2\ncn3/);
            $omim_entry->contributors( $val );
        }
        elsif ( $key == CREATED_BY_STATE ) {
            $val = undef if($val =~ /cd1\ncd2\ncd3/);
            $omim_entry->created( $val );
        }
        elsif ( $key == EDITED_BY_STATE ) {
            $val = undef if($val =~ /ed1\ned2\ned3/);
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

    # parse the symptoms text into a hash-based structure.
    $self->_finer_parse_symptoms($omim_entry);
    
    return $omim_entry;

} # _createOMIMentry


sub _finer_parse_symptoms {
    my ($self, $omim_entry) = @_;
    my $text = $omim_entry->clinical_symptoms_raw;
    if( $text ) { 
	my $part;
	for my $line (split /\n/, $text){
		if ($line =~ /^([\w\s,]+)\:\s*$/) {
		$part = $1;
	    } elsif( $line =~ /^\s+$/ ) {
	    } elsif($line =~ /^(\s+)([^;]+)\;?\s*$/){
		my $symptom = $2;
		if( ! $part ) { 
		    # $self->warn("$text\nline='$line'\n");
		    next;
		}
		$omim_entry->add_clinical_symptoms($part, $symptom);
	    }
	}
    }
    $omim_entry->clinical_symptoms_raw('');
}

sub _parse_genemap {
     my ( $self, $omim_entry, $val ) = @_;
     
     my $genemap_line = ${ $self->_genemap_hash() }{ $val };
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
    my $number          = "";
    my $title           = "";
    my $symbol_mut_line = "";
    my $prev_line       = "";
    my $description     = "";
    my $saw_empty_line  = FALSE;
     
    my @lines = split( /\n/, $text );
    
    foreach my $line ( @lines ) {
        if ( $line !~ /\w/ ) {
             $saw_empty_line = TRUE;
        }
        elsif ( $line =~ /^\s*(\.\d+)/ ) {
            my $current_number = $1;
            if ( $number ne "" ) {
                my $allelic_variant = $self->_create_allelic_variant( $number, $title, 
                                                    $symbol_mut_line, $description );
                
                push( @allelic_variants, $allelic_variant );
            }
            $number          = $current_number;
            $title           = "";
            $prev_line       = "";
            $symbol_mut_line = "";
            $description     = "";
            $saw_empty_line  = FALSE;
        }
        elsif ( $title eq "" ) {
            $title = $line;
        }
        elsif ( $saw_empty_line == FALSE ) {
            $prev_line = $line;
        }
        elsif ( $saw_empty_line == TRUE ) {
            if ( $prev_line ne "" ) {
                $symbol_mut_line = $prev_line;
                $prev_line       = "";
            }
            if ( $description ne "" ) {
                $description .= "\n" . $line;
            }
            else {
                $description = $line;
            }
        }
    }
    
    my $allelic_variant = $self->_create_allelic_variant( $number, $title, 
                                        $symbol_mut_line, $description );
    
    push( @allelic_variants, $allelic_variant );
    
    return @allelic_variants;
    
} # _parse_allelic_variants




sub _create_allelic_variant {
    my ( $self, $number, $title, $symbol_mut_line, $description ) = @_;
    
    my $symbol   = "";
    my $mutation = "";
    my $aa_ori   = "";
    my $aa_mut   = "";
    my $position = "";
   
    if ( $symbol_mut_line =~ /\s*(.+?)\s*,\s*([a-z]{3})(\d+)([a-z]{3})/i ) {
         $symbol   = $1;
         $aa_ori   = $2;
         $aa_mut   = $4;
         $position = $3;
    }
    elsif ( $symbol_mut_line =~ /\s*(.+?)\s*,\s*(.+)/ ) {
         $symbol   = $1;
         $mutation = $2;
    }
    else {
         $symbol = $symbol_mut_line;
    }
    
    if ( ! defined( $description ) ) { $self->throw("undef desc"); }
    if ( ! defined( $mutation ) )   { $self->throw("undef mutation"); }
  
    
    my $allelic_variant = Bio::Phenotype::OMIM::OMIMentryAllelicVariant->new();
    $allelic_variant->number( $number );
    $allelic_variant->aa_ori( $aa_ori );
    $allelic_variant->aa_mut( $aa_mut );
    $allelic_variant->position( $position );
    $allelic_variant->title( $title );
    $allelic_variant->symbol( $symbol );
    $allelic_variant->description( $description );
    $allelic_variant->additional_mutations( $mutation );
     
    return $allelic_variant; 
    
} # _create_allelic_variant




sub _parse_title {
    my ( $self, $text ) = @_;
    my $title = "";
    if ( $text =~ /^(.+)\n/ ) {
        $title = $1;
        $text  =~ s/^.+\n//;
    }
    else {
        $title = $text;
        $text  = "";
    
    }
    
    return ( $title, $text ); 
} # _parse_title




sub _parse_references {
    my ( $self, $text ) = @_;
    
    $text =~ s/\A\s+//;
    $text =~ s/\s+\z//;
    $text =~ s/\A\d+\.\s*//;
    
    my @references;
    
    my @texts = split( /\s*\n\s*\n\s*\d+\.\s*/, $text );
    
    foreach my $t ( @texts ) {
    
        my $authors   = "";
        my $title     = "";
        my $location  = "";
        
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
         
        my $ref = Bio::Annotation::Reference->new( -title    => $title,
                                                   -location => $location,
                                                   -authors  => $authors );
        push( @references, $ref );
       
    }
    return @references;
    
} # _parse_references




sub _genemap_hash {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        unless ( ref( $value ) eq "HASH" ) {
            $self->throw( "Argument to method \"_genemap_hash\" is not a reference to an Hash" );
        }
        $self->{ "_genemap_hash" } = $value;
     
    }
    
    return $self->{ "_genemap_hash" };
} # _genemap_hash




sub _is_not_first_record {

    my ( $self, $value ) = @_;

    if ( defined $value ) {
        unless ( $value == FALSE || $value == TRUE ) {
            $self->throw( "Found [$value] where [" . TRUE
            ." or " . FALSE . "] expected" );
        }
        $self->{ "_not_first_record" } = $value;
    }
    
    return $self->{ "_not_first_record" };
} # _is_not_first_record




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




sub _OMIM_text_file {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        unless ( $value->isa( "Bio::Root::IO" ) ) {
            $self->throw( "[$value] is not a valid \"Bio::Root::IO\"" );
        }
        $self->{ "_omimtxt_file" } = $value;
     
    }
    
    return $self->{ "_omimtxt_file" };
} # _OMIM_text_file




sub _read_genemap {
    my ( $self, $genemap_file_name ) = @_;
    
    my $line         = "";
    my %genemap_hash = ();
    my $genemap_file = Bio::Root::IO->new( -file => $genemap_file_name );
    my @a            = ();
    my %gm           = ();
    
    while( $line = $genemap_file->_readline() ) {
        @a = split( /\|/, $line );
        unless( scalar( @a ) == 18 ) {
            $self->throw( "Gene map file \"".$self->genemap_file_name()
            . "\" is not in the expected format."
            . " Make sure there is a linebreak after the final line." );
        }
        $gm{ $a[ 9 ] } = $line;
    }
    $genemap_file->close();
    $self->_genemap_hash( \%gm );
  
} #_read_genemap 




sub _no_OMIM_text_file_provided_error {
    my ( $self ) = @_;

    my $msg =  "Need to indicate a OMIM text file to read from with\n";
    $msg .= "either \"OMIMparser->new( -omimtext => \"path/to/omim.txt\" );\"\n";
    $msg .= "or \"\$omim_parser->omimtxt_file_name( \"path/to/omim.txt\" );\"";
    $self->throw( $msg );
} # _no_OMIM_text_file_provided_error




sub _not_a_OMIM_text_file_error {
    my ( $self ) = @_;

    my $msg =  "File \"".$self->omimtxt_file_name() . 
    "\" appears not to be a OMIM text file";
    $self->throw( $msg );
} # _not_a_OMIM_text_file_error




sub _add_to_hash {
    my ( $self, $state, $contents, $record_ref ) = @_;
    
    if ( exists( $record_ref->{ $state } ) ) {
        chomp( $record_ref->{ $state } );
        $record_ref->{ $state } = $record_ref->{ $state } . $contents;
    }
    else {
        $record_ref->{ $state } = $contents;
    }
} # _add_to_hash



1;
