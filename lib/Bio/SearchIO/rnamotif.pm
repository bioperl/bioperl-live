#
# BioPerl module for Bio::SearchIO::rnamotif
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chris Fields <cjfields-at-uiuc-dot-edu>
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::rnamotif - SearchIO-based RNAMotif parser

=head1 SYNOPSIS

  # do not call this module directly. Use Bio::SearchIO.

=head1 DESCRIPTION

This is a highly experimental SearchIO-based parser for output from the rnamotif
program (one of the programs in the RNAMotif suite). It currently parses only
raw rnamotif output for RNAMotif versions 3.0 and above; older versions may work
but will not be supported. rmfmt output will not be supported at this time.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Chris Fields

Email cjfields-at-uiuc-dot-edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SearchIO::rnamotif;
use strict;

use base qw(Bio::SearchIO);

my %MODEMAP = (
	    'Result'             => 'result',
	    'Hit'                => 'hit',
	    'Hsp'                => 'hsp'
	    );

my %MAPPING = ( 
        # commented out tags have not been assigned
        
        'Hsp_score'        => 'HSP-score',
        'Hsp_custom-data'  => 'HSP-custom_score',
        
        # rnamotif has no evalue
        
        # descriptor has no start, end; same as hit start, end
        'Hsp_query-from'  => 'HSP-query_start',
        'Hsp_query-to'    => 'HSP-query_end',
        'Hsp_hit-from'    => 'HSP-hit_start', 
        'Hsp_hit-to'      => 'HSP-hit_end',
        
        # descriptor has no start, end
        
        'Hsp_hseq'        => 'HSP-hit_seq',
        'Hsp_align-len'   => 'HSP-hsp_length',
        
        # build this from scratch, simple WUSS-format
        'Hsp_structure'   => 'HSP-meta',
        'Hsp_stranded'    => 'HSP-stranded',        
        
        # not supported for RNAMotif

        'Hit_id'        => 'HIT-name',
        'Hit_accession' => 'HIT-accession',
        'Hit_gi'        => 'HIT-ncbi_gi',
        'Hit_def'       => 'HIT-description',
        'Hit_score'     => 'HIT-score', # best HSP score
 
        'RNAMotif_program'  => 'RESULT-algorithm_name', # get/set 
        'RNAMotif_version'  => 'RESULT-algorithm_version', # get/set 
        'RNAMotif_query-def'=> 'RESULT-query_name', # get/set
        # No length (query is a descriptor)
        'RNAMotif_query-acc'=> 'RESULT-query_accession', # get/set 
        'RNAMotif_querydesc'=> 'RESULT-query_description', # get/set
        'RNAMotif_db'       => 'RESULT-database_name',  # get/set 
	     );

# use structure_delimiters to set custom delimiters

my @VALID_SYMBOLS = qw(5-prime 3-prime single-strand unknown);
my %STRUCTURE_SYMBOLS = (
                   '5-prime'        => '<',
                   '3-prime'        => '>',
                   'single-strand'  => '.',
                   'unknown'        => '?'
                    # may add more for quartets, triplets
                  );

my $DEFAULT_VERSION = '3.0.3';

=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO->new();
 Function: Builds a new Bio::SearchIO::rnamotif object 
 Returns : Bio::SearchIO::rnamotif parser
 Args    : -fh/-file     => RNAMotif filename
           -format       => 'rnamotif'
           -model        => query model (or descriptor, in this case)
           -database     => database name (default undef)
           -query_acc    => query accession (default undef)
           -hsp_minscore => minimum HSP score cutoff
           -hsp_maxscore => maximum HSP score cutoff
           -symbols      => hash ref of structure symbols to use
                            (default symbols in %STRUCTURE_SYMBOLS hash)

=cut

sub _initialize {
    my ( $self, @args ) = @_;
    $self->SUPER::_initialize(@args);
    my ($version, $model, $database, $maxcutoff, $mincutoff, $seqdistance,
        $accession, $symbols) =
       $self->_rearrange([qw(VERSION MODEL DATABASE HSP_MAXSCORE 
                          HSP_MINSCORE SEQ_DISTANCE QUERY_ACC SYMBOLS)],@args);
    my $handler = $self->_eventHandler;
    $handler->register_factory(
        'result',
        Bio::Factory::ObjectFactory->new(
            -type      => 'Bio::Search::Result::GenericResult',
            -interface => 'Bio::Search::Result::ResultI',
            -verbose => $self->verbose
        )
    );

    $handler->register_factory(
        'hit',
        Bio::Factory::ObjectFactory->new(
            -type      => 'Bio::Search::Hit::ModelHit',
            -interface => 'Bio::Search::Hit::HitI',
            -verbose => $self->verbose
        )
    );

    $handler->register_factory(
        'hsp',
        Bio::Factory::ObjectFactory->new(
            -type      => 'Bio::Search::HSP::ModelHSP',
            -interface => 'Bio::Search::HSP::HSPI',
            -verbose => $self->verbose
        )
    );
    $model      && $self->model($model);
    $database   && $self->database($database);
    $accession  && $self->query_accession($accession);
    $version ||= $DEFAULT_VERSION;
    $self->algorithm_version($version);
    $self->throw("Cannot define both a minimal and maximal cutoff")
           if (defined($mincutoff) && defined($maxcutoff));
    defined($mincutoff)   && $self->hsp_minscore($mincutoff);
    defined($maxcutoff)   && $self->hsp_maxscore($maxcutoff);
    $symbols  ||= \%STRUCTURE_SYMBOLS;
    $self->structure_symbols($symbols);
}

=head2 next_result

 Title   : next_result
 Usage   : my $hit = $searchio->next_result;
 Function: Returns the next Result from a search
 Returns : Bio::Search::Result::ResultI object
 Args    : none

=cut

sub next_result {
    my ($self) = @_;
    my $seentop = 0;
    local $/ = "\n";
    local $_;
    my ($rm, $d, $descriptor, $file, $oktobuild);
    my ($hitid, $hitdesc, $hspid, $lastid, $lastscore);
    my $sprintf;
    
    # user-determined Result data
    my ($accession, $db, $model) =
       ($self->query_accession, $self->database, $self->model);
    # HSP building options
    my $hsp_min = $self->hsp_minscore;
    my $hsp_max = $self->hsp_maxscore;
    my $version = $self->algorithm_version;
    my $laststart;
    
    my $verbose = $self->verbose;    # cache for speed?
    $self->start_document();
    PARSER:
    while ( defined( my $line = $self->_readline ) ) {
        # start of report
        next if $line =~ m{^\s+$};
        if (index($line,'#RM') == 0) {
            if (index($line,'#RM scored') == 0 ) {
                if ($seentop) {
                    $self->_pushback($line);
                    last PARSER;
                }
                $self->start_element({'Name' => 'Result'});
                $self->element_hash({
                    'RNAMotif_program'      => 'rnamotif',
                    'RNAMotif_version'      => $version,
                    'RNAMotif_query-acc'    => $accession,
                    'RNAMotif_db'           => $db
                    });                
                $seentop = 1;
                #$self->debug("Start result\n");
            } elsif (index($line,'#RM descr') == 0) {
                ($rm, $d, $descriptor) = split ' ', $line, 3;
                # toss $rm, $d; keep $descr
                chomp $descriptor;
                $self->{'_descriptor'} = $descriptor;
                $self->element(
                               {'Name' => 'RNAMotif_querydesc',
                                'Data' => $descriptor}
                              );
            } elsif(index($line,'#RM dfile') == 0) {
                ($rm, $d, $file) = split ' ', $line, 3;
                # toss $rm, $d; keep $file
                chomp $file;
                $self->element(
                               {'Name' => 'RNAMotif_query-def',
                                'Data' => $file}
                              );
            } else {
                $self->debug("Unrecognized line: $line");
            }
        } elsif ($line =~ s{^>}{}) {
            chomp $line;
            ($hitid, $hitdesc) = split ' ',$line,2;
            
            if ($self->within_element('hit') && ($hitid ne $lastid)) {
                $self->element(
                       {'Name' => 'Hit_score',
                        'Data' => $lastscore}
                      ) if $lastscore;
                $self->end_element({'Name' => 'Hit'});
                $self->start_element({'Name' => 'Hit'});
            } elsif (!$self->within_element('hit')) {
                $self->start_element({'Name' => 'Hit'});
            }
            my ($gi, $acc, $ver) = $self->_get_seq_identifiers($hitid);
            
            $self->element_hash({
                'Hit_id'        => $hitid,
                'Hit_gi'        => $gi,
                'Hit_accession' => $ver ? "$acc.$ver" :
                                    $acc ? $acc : $hitid,
                'Hit_def'       => $hitdesc}
              );
            $lastid = $hitid;
        } elsif ($line =~ m{^(\S+)\s+(.+?)\s+(\d+)\s+(\d+)\s+(\d+)\s(.*)$}xm) {
            chomp $line;
            my $hspid = $1;
            my ($score, $strand, $start, $length , $seq) = ($2, $3, $4, $5, $6);
            $score *= 1;  # implicitly cast any odd '0.000' to float
            # sanity check ids
            unless ($hitid eq $hspid) {
                $self->throw("IDs do not match!");
            }
            # check score for possible sprintf data, mark as such, cache result
            if (!defined($sprintf)) {
                if ($score =~ m{[^0-9.-]+}gxms) {
                    if (defined $hsp_min || defined $hsp_max ) {
                        $self->warn("HSP data likely contains custom score; ".
                                    "ignoring min/maxscore");
                    }
                    $sprintf = $oktobuild = 1;
                } else {
                    $sprintf = 0;
                }
            }
            
            if (!$sprintf) {
                if (($hsp_min && $score <= $hsp_min) 
                          || ($hsp_max && ($score >= $hsp_max)) ) {
                    # do not build HSP
                    $oktobuild = 0;
                } else {
                    $oktobuild = 1;
                    
                    # store best hit score based on the hsp min/maxscore only
                    if (defined $hsp_min && $score > $hsp_min) {
                        $lastscore = $score if !$lastscore || $score > $lastscore;
                    } elsif (defined $hsp_max && $score < $hsp_max) {
                        $lastscore = $score if !$lastscore || $score < $lastscore;
                    } 
                }
            }
            
            # build HSP
            if ($oktobuild) {
                my $end;
                # calculate start/end
                if( $strand==0 ) {
                    $end = $start + $length -1;
                } else {
                    $end = $start - $length + 1;
                }
                
                my ($rna, $meta) = $self->_motif2meta($seq, $descriptor);
                
                $self->start_element({'Name' => 'Hsp'});
                my $rnalen = $rna =~ tr{ATGCatgc}{ATGCatgc};
                $self->element_hash({
                        'Hsp_stranded'      => 'HIT', 
                        'Hsp_hseq'          => $rna,
                        'Hsp_query-from'    => 1,
                        'Hsp_query-to'      =>length($rna),
                        'Hsp_hit-from'      => $start,
                        'Hsp_hit-to'        => $end,
                        'Hsp_structure'     => $meta,
                        'Hsp_align-len'     => length($rna),
                        'Hsp_score'         => $sprintf ? undef : $score,
                        'Hsp_custom-data'   => $sprintf ? $score : undef,
                        });
                $self->end_element({'Name' => 'Hsp'});
                $oktobuild = 0 if (!$sprintf);
            }
        }
    }
    if ($self->within_element('hit')) {
        $self->element(
               {'Name' => 'Hit_score',
                'Data' => $lastscore}
              ) if $lastscore;
        $self->end_element( { 'Name' => 'Hit' } );
    }
    if ($seentop) {
        $self->end_element( { 'Name' => 'Result' } );
    }
    return $self->end_document();
}

=head2 start_element

 Title   : start_element
 Usage   : $eventgenerator->start_element
 Function: Handles a start element event
 Returns : none
 Args    : hashref with at least 2 keys 'Data' and 'Name'


=cut

sub start_element {
    my ( $self, $data ) = @_;

    # we currently don't care about attributes
    my $nm   = $data->{'Name'};
    my $type = $MODEMAP{$nm};
    if ($type) {
        if ( $self->_eventHandler->will_handle($type) ) {
            my $func = sprintf( "start_%s", lc $type );
            $self->_eventHandler->$func( $data->{'Attributes'} );
        }
        unshift @{ $self->{'_elements'} }, $type;
    }
    if ( defined $type
        && $type eq 'result' )
    {
        $self->{'_values'} = {};
        $self->{'_result'} = undef;
    }
}

=head2 end_element

 Title   : start_element
 Usage   : $eventgenerator->end_element
 Function: Handles an end element event
 Returns : none
 Args    : hashref with at least 2 keys, 'Data' and 'Name'


=cut

sub end_element {
    my ( $self, $data ) = @_;
    my $nm   = $data->{'Name'};
    my $type = $MODEMAP{$nm};
    my $rc;

    if ($type) {
        if ( $self->_eventHandler->will_handle($type) ) {
            my $func = sprintf( "end_%s", lc $type );
            $rc = $self->_eventHandler->$func( $self->{'_reporttype'},
                $self->{'_values'} );
        }
        my $lastelem = shift @{ $self->{'_elements'} };
    }
    elsif ( $MAPPING{$nm} ) {
        if ( ref( $MAPPING{$nm} ) =~ /hash/i ) {
            my $key = ( keys %{ $MAPPING{$nm} } )[0];
            $self->{'_values'}->{$key}->{ $MAPPING{$nm}->{$key} } =
              $self->{'_last_data'};
        }
        else {
            $self->{'_values'}->{ $MAPPING{$nm} } = $self->{'_last_data'};
        }
    }
    else {
        $self->debug("unknown nm $nm, ignoring\n");
    }
    $self->{'_last_data'} = '';    # remove read data if we are at
                                   # end of an element
    $self->{'_result'} = $rc if ( defined $type && $type eq 'result' );
    return $rc;
}

=head2 element

 Title   : element
 Usage   : $eventhandler->element({'Name' => $name, 'Data' => $str});
 Function: Convenience method that calls start_element, characters, end_element
 Returns : none
 Args    : Hash ref with the keys 'Name' and 'Data'


=cut

sub element {
    my ( $self, $data ) = @_;
    # simple data calls (%MAPPING) do not need start_element
    $self->characters($data);
    $self->end_element($data);
}


=head2 element_hash

 Title   : element
 Usage   : $eventhandler->element_hash({'Hsp_hit-from' => $start,
                                        'Hsp_hit-to'   => $end,
                                        'Hsp_score'    => $lastscore});
 Function: Convenience method that takes multiple simple data elements and
           maps to appropriate parameters
 Returns : none
 Args    : Hash ref with the mapped key (in %MAPPING) and value

=cut

sub element_hash {
    my ($self, $data) = @_;
    $self->throw("Must provide data hash ref") if !$data || !ref($data);
    for my $nm (sort keys %{$data}) {
        next if $data->{$nm} && $data->{$nm} =~ m{^\s*$}o;
        if ( $MAPPING{$nm} ) {
            if ( ref( $MAPPING{$nm} ) =~ /hash/i ) {
                my $key = ( keys %{ $MAPPING{$nm} } )[0];
                $self->{'_values'}->{$key}->{ $MAPPING{$nm}->{$key} } =
                  $data->{$nm};
            }
            else {
                $self->{'_values'}->{ $MAPPING{$nm} } = $data->{$nm};
            }
        }
    }
}

=head2 characters

 Title   : characters
 Usage   : $eventgenerator->characters($str)
 Function: Send a character events
 Returns : none
 Args    : string


=cut

sub characters {
    my ( $self, $data ) = @_;
    return unless ( defined $data->{'Data'} && $data->{'Data'} !~ /^\s+$/o );
    $self->{'_last_data'} = $data->{'Data'};
}

=head2 within_element

 Title   : within_element
 Usage   : if( $eventgenerator->within_element($element) ) {}
 Function: Test if we are within a particular element
           This is different than 'in' because within can be tested
           for a whole block.
 Returns : boolean
 Args    : string element name 

=cut

sub within_element {
    my ( $self, $name ) = @_;
    return 0
      if ( !defined $name
        || !defined $self->{'_elements'}
        || scalar @{ $self->{'_elements'} } == 0 );
    foreach ( @{ $self->{'_elements'} } ) {
        return 1 if ( $_ eq $name );
    }
    return 0;
}

=head2 in_element

 Title   : in_element
 Usage   : if( $eventgenerator->in_element($element) ) {}
 Function: Test if we are in a particular element
           This is different than 'within' because 'in' only 
           tests its immediate parent.
 Returns : boolean
 Args    : string element name 

=cut

sub in_element {
    my ( $self, $name ) = @_;
    return 0 if !defined $self->{'_elements'}->[0];
    return ( $self->{'_elements'}->[0] eq $name );
}

=head2 start_document

 Title   : start_document
 Usage   : $eventgenerator->start_document
 Function: Handle a start document event
 Returns : none
 Args    : none

=cut

sub start_document {
    my ($self) = @_;
    $self->{'_lasttype'} = '';
    $self->{'_values'}   = {};
    $self->{'_result'}   = undef;
    $self->{'_elements'} = [];
}

=head2 end_document

 Title   : end_document
 Usage   : $eventgenerator->end_document
 Function: Handles an end document event
 Returns : Bio::Search::Result::ResultI object
 Args    : none

=cut

sub end_document {
    my ($self) = @_;
    return $self->{'_result'};
}

=head2 result_count

 Title   : result_count
 Usage   : my $count = $searchio->result_count
 Function: Returns the number of results we have processed
 Returns : integer
 Args    : none

=cut

sub result_count {
    my $self = shift;
    return $self->{'_result_count'};
}

=head2 descriptor

 Title   : descriptor
 Usage   : my $descr = $parser->descriptor();
 Function: Get/Set descriptor name.  Some versions of RNAMotif do not add the
           descriptor name to the output.  This also overrides any name found
           while parsing.
 Returns : String (name of model)
 Args    : [optional] String (name of model)

=cut

sub descriptor {
    my $self = shift;
    return $self->{'_descriptor'} = shift if @_;
    return $self->{'_descriptor'};
}

=head2 model

 Title   : model
 Usage   : my $model = $parser->model();
 Function: Get/Set model; Infernal currently does not output
           the model name (Rfam ID)
 Returns : String (name of model)
 Args    : [optional] String (name of model)
 Note    : this is a synonym for descriptor()

=cut

sub model { shift->descriptor(@_) }

=head2 database

 Title   : database
 Usage   : my $database = $parser->database();
 Function: Get/Set database; Infernal currently does not output
           the database name
 Returns : String (database name)
 Args    : [optional] String (database name)

=cut

sub database {
    my $self = shift;
    return $self->{'_database'} = shift if @_;
    return $self->{'_database'};
}

=head2 query_accession

 Title   : query_accession
 Usage   : my $acc = $parser->query_accession();
 Function: Get/Set query (model) accession; RNAMotif currently does not output
           the accession number
 Returns : String (accession)
 Args    : [optional] String (accession)

=cut

sub query_accession {
    my $self = shift;
    return $self->{'_query_accession'} = shift if @_;
    return $self->{'_query_accession'};
}

=head2 algorithm_version

 Title   : algorithm_version
 Usage   : my $ver = $parser->algorithm_version();
 Function: Get/Set algorithm version (not defined in RNAMotif output)
 Returns : String (accession)
 Args    : [optional] String (accession)

=cut

sub algorithm_version {
    my $self = shift;
    return $self->{'_algorithm'} = shift if @_;
    return $self->{'_algorithm'};
}

=head2 hsp_minscore

 Title   : hsp_minscore
 Usage   : my $cutoff = $parser->hsp_minscore();
 Function: Get/Set min score cutoff (for generating Hits/HSPs).
 Returns : score (number)
 Args    : [optional] score (number)
 Note    : Cannot be set along with hsp_maxscore()

=cut

sub hsp_minscore {
    my ($self, $score) = shift;
    $self->throw('Minscore not set to a number') if
        ($score && $score !~ m{[0-9.]+});
    return $self->{'_hsp_minscore'} = shift if @_;
    return $self->{'_hsp_minscore'};
}

=head2 hsp_maxscore

 Title   : hsp_maxscore
 Usage   : my $cutoff = $parser->hsp_maxscore();
 Function: Get/Set max score cutoff (for generating Hits/HSPs).
 Returns : score (number)
 Args    : [optional] score (number)
 Note    : Cannot be set along with hsp_minscore()

=cut

sub hsp_maxscore {
    my ($self, $score) = shift;
    $self->throw('Maxscore not set to a number') if
        ($score && $score !~ m{[0-9.]+});
    return $self->{'_hsp_maxscore'} = shift if @_;
    return $self->{'_hsp_maxscore'};
}

=head2 structure_symbols

 Title   : structure_symbols
 Usage   : my $hashref = $parser->structure_symbols();
 Function: Get/Set RNA structure symbols
 Returns : Hash ref of delimiters (5' stem, 3' stem, single-strand, etc)
         : default = < (5-prime)
                     > (3-prime)
                     . (single-strand)
                     ? (unknown) 
 Args    : Hash ref of substitute delimiters, using above keys.

=cut

sub structure_symbols {
    my ($self, $delim) = @_;
    if ($delim) {
        if (ref($delim) =~ m{HASH}) {
            my %data = %{ $delim };
            for my $d (@VALID_SYMBOLS) {
                if ( exists $data{$d} ) {
                    $self->{'_delimiter'}->{$d} = $data{$d};
                }
            }
        } else {
            $self->throw("Args to helix_delimiters() should be in a hash reference");
        }
    }
    return $self->{'_delimiter'};
}

#Private methods

=head2 _motif2meta

 Title   : _motif2meta
 Usage   : my ($rna, $meta) = $parser->_motif2meta($str, $descr);
 Function: Creates meta string from sequence and descriptor
 Returns : array of sequence, meta strings
 Args    : Array of string data and descriptor data

 Note: This is currently a quick and simple way of making simple
 RNA structures (stem-loops, helices, etc) from RNAMotif descriptor
 data in the output file.  It does not currently work with pseudoknots,
 triplets, G-quartets, or other more complex RNA structural motifs.

=cut

sub _motif2meta {
    my ($self, $str, $descriptor) = @_;
    my ($rna, $meta);
    my @desc_el = split ' ',$descriptor;
    my @seq_el = split ' ',$str;
    my $symbol = $self->structure_symbols();
    if ($#desc_el != $#seq_el) {
        $self->throw("Descriptor elements and seq elements do not match");
    }
    while (@desc_el) {
        my $struct;
        my ($seq, $motif) = (shift @seq_el, shift @desc_el);
        $struct = (index($motif,'h5') == 0) ? $symbol->{'5-prime'} :
                  (index($motif,'h3') == 0) ? $symbol->{'3-prime'} :
                  (index($motif,'ss') == 0) ? $symbol->{'single-strand'}  :
                  (index($motif,'ctx')== 0) ? $symbol->{'single-strand'}  :
                  $symbol->{'unknown'};
        $meta .= $struct x (length($seq));
        $rna .= $seq;
    }
    return ($rna, $meta);
}

1;
