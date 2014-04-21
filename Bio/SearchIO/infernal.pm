#
# BioPerl module for Bio::SearchIO::infernal
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

Bio::SearchIO::infernal - SearchIO-based Infernal parser

=head1 SYNOPSIS

  my $parser = Bio::SearchIO->new(-format => 'infernal',
                                  -file => 'purine.inf');
  while( my $result = $parser->next_result ) {
        # general result info, such as model used, Infernal version
        while( my $hit = $result->next_hit ) {
            while( my $hsp = $hit->next_hsp ) {
                # ...
            }
        }
  }

=head1 DESCRIPTION

This is a SearchIO-based parser for Infernal output from the cmsearch program.
It currently parses cmsearch output for Infernal versions 0.7-1.0; older
versions may work but will not be supported.

As the first stable version has been released (and output has stabilized) it is
highly recommended that users upgrade to using the latest Infernal release.
Support for the older pre-v.1 developer releases will be dropped for future core
1.6 releases. 

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

=head1 CONTRIBUTORS

  Jeffrey Barrick, Michigan State University

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SearchIO::infernal;
use strict;

use Data::Dumper;
use base qw(Bio::SearchIO);

our %MODEMAP = (
	    'Result'             => 'result',
	    'Hit'                => 'hit',
	    'Hsp'                => 'hsp'
	    );

our %MAPPING = ( 
        'Hsp_bit-score'   => 'HSP-bits',
        'Hsp_score'       => 'HSP-score',
        'Hsp_evalue'      => 'HSP-evalue', # evalues only in v0.81, optional
        'Hsp_pvalue'      => 'HSP-pvalue', # pvalues only in v0.81, optional
        'Hsp_query-from'  => 'HSP-query_start',
        'Hsp_query-to'    => 'HSP-query_end',
        'Hsp_hit-from'    => 'HSP-hit_start', 
        'Hsp_hit-to'      => 'HSP-hit_end', 
        'Hsp_gaps'        => 'HSP-hsp_gaps', 
        'Hsp_hitgaps'     => 'HSP-hit_gaps',
        'Hsp_querygaps'   => 'HSP-query_gaps',
        'Hsp_qseq'        => 'HSP-query_seq',
        'Hsp_hseq'        => 'HSP-hit_seq',
        'Hsp_midline'     => 'HSP-homology_seq',
        'Hsp_structure'   => 'HSP-meta',
        'Hsp_align-len'   => 'HSP-hsp_length',
        'Hsp_stranded'    => 'HSP-stranded',
        
        'Hit_id'        => 'HIT-name',
        'Hit_len'       => 'HIT-length',
        'Hit_gi'        => 'HIT-ncbi_gi',
        'Hit_accession' => 'HIT-accession',
        'Hit_def'       => 'HIT-description',
        'Hit_signif'    => 'HIT-significance', # evalues only in v0.81, optional
        'Hit_p'         => 'HIT-p',            # pvalues in 1.0, optional
        'Hit_score'     => 'HIT-score', # best HSP bit score
        'Hit_bits'      => 'HIT-bits', # best HSP bit score
 
        'Infernal_program'  => 'RESULT-algorithm_name', # get/set 
        'Infernal_version'  => 'RESULT-algorithm_version', # get/set 
        'Infernal_query-def'=> 'RESULT-query_name', # get/set 
        'Infernal_query-len'=> 'RESULT-query_length', 
        'Infernal_query-acc'=> 'RESULT-query_accession', # get/set 
        'Infernal_querydesc'=> 'RESULT-query_description', # get/set
        'Infernal_db'       => 'RESULT-database_name',  # get/set 
        'Infernal_db-len'   => 'RESULT-database_entries', # none yet
        'Infernal_db-let'   => 'RESULT-database_letters', # none yet
	     );

my $MINSCORE = 0;
my $DEFAULT_ALGORITHM = 'cmsearch';
my $DEFAULT_VERSION = '1.0';

my @VALID_SYMBOLS = qw(5-prime 3-prime single-strand unknown gap);
my %STRUCTURE_SYMBOLS = (
                   '5-prime'        => '<',
                   '3-prime'        => '>',
                   'single-strand'  => ':',
                   'unknown'        => '?',
                   'gap'            => '.'
                   );

=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::infernal->new();
 Function: Builds a new Bio::SearchIO::infernal object 
 Returns : Bio::SearchIO::infernal
 Args    : -fh/-file      => cmsearch (infernal) filename
           -format        => 'infernal'
           -model         => query model (Rfam ID) (default undef)
           -database      => database name (default undef)
           -query_acc     => query accession, eg. Rfam accession RF####
           -query_desc    => query description, eg. Rfam description
           -hsp_minscore  => minimum HSP score cutoff
           -convert_meta  => boolean, set to convert meta string to simple WUSS format
           -symbols       => hash ref of structure symbols to use
                             (default symbols in %STRUCTURE_SYMBOLS hash)

=cut

sub _initialize {
    my ( $self, @args ) = @_;
    $self->SUPER::_initialize(@args);
    my ($model, $database, $convert, $symbols, $cutoff,
        $desc, $accession, $algorithm, $version) =
        $self->_rearrange([qw(MODEL
                          DATABASE
                          CONVERT_META
                          SYMBOLS
                          HSP_MINSCORE
                          QUERY_DESC
                          QUERY_ACC
                          ALGORITHM
                          VERSION)],@args);
    my $handler = $self->_eventHandler;
    $handler->register_factory(
        'result',
        Bio::Factory::ObjectFactory->new(
            -type      => 'Bio::Search::Result::GenericResult',
            -interface => 'Bio::Search::Result::ResultI',
            -verbose   => $self->verbose
        )
    );

    $handler->register_factory(
        'hit',
        Bio::Factory::ObjectFactory->new(
            -type      => 'Bio::Search::Hit::ModelHit',
            -interface => 'Bio::Search::Hit::HitI',
            -verbose   => $self->verbose
        )
    );

    $handler->register_factory(
        'hsp',
        Bio::Factory::ObjectFactory->new(
            -type      => 'Bio::Search::HSP::ModelHSP',
            -interface => 'Bio::Search::HSP::HSPI',
            -verbose   => $self->verbose
        )
    );
	
    defined $model     && $self->model($model);
    defined $database  && $self->database($database);
    defined $accession && $self->query_accession($accession);
    defined $convert   && $self->convert_meta($convert);
    defined $desc      && $self->query_description($desc);
    
    $version ||= $DEFAULT_VERSION;
    $self->version($version);
    $symbols ||= \%STRUCTURE_SYMBOLS;
    $self->structure_symbols($symbols);
    $cutoff ||= $MINSCORE;
    $self->hsp_minscore($cutoff);
    $algorithm ||= $DEFAULT_ALGORITHM;
    $self->algorithm($algorithm);
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
    unless (exists $self->{'_handlerset'}) {
        my $line;
        while ($line = $self->_readline) {
            # advance to first line
            next if $line =~ m{^\s*$};
            # newer output starts with model name
            if ($line =~ m{^\#\s+cmsearch\s}) {
                $self->{'_handlerset'} = 'latest';
			} elsif ($line =~ m{^CM\s\d+:}) {
                $self->{'_handlerset'} = 'pre-1.0';
            } else {
                $self->{'_handlerset'} ='old';
            }
            last;
        }
        $self->_pushback($line);
		#if ($self->{'_handlerset'} ne '1.0') {
		#	$self->deprecated(
		#	-message => "Parsing of Infernal pre-1.0 release is deprecated;\n".
		#		"upgrading to Infernal 1.0 or above is highly recommended",
		#	-version => 1.007);
		#}
    }
    return ($self->{'_handlerset'} eq 'latest')  ? $self->_parse_latest :
		   ($self->{'_handlerset'} eq 'pre-1.0') ? $self->_parse_pre :
			$self->_parse_old;
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

=head2 model

 Title   : model
 Usage   : my $model = $parser->model();
 Function: Get/Set model; Infernal currently does not output
           the model name (Rfam ID)
 Returns : String (name of model)
 Args    : [optional] String (name of model)

=cut

sub model {
    my $self = shift;
    return $self->{'_model'} = shift if @_;
    return $self->{'_model'};
}

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

=head2 algorithm

 Title   : algorithm
 Usage   : my $algorithm = $parser->algorithm();
 Function: Get/Set algorithm; current versions of Infernal do not output
           the algorithm name
 Returns : String (algorithm name)
 Args    : [optional] String (algorithm name)

=cut

sub algorithm {
    my $self = shift;
    return $self->{'_algorithm'} = shift if @_;
    return $self->{'_algorithm'};
}

=head2 query_accession

 Title   : query_accession
 Usage   : my $acc = $parser->query_accession();
 Function: Get/Set query (model) accession; Infernal currently does not output
           the accession number (Rfam accession #)
 Returns : String (accession)
 Args    : [optional] String (accession)

=cut

sub query_accession {
    my $self = shift;
    return $self->{'_query_accession'} = shift if @_;
    return $self->{'_query_accession'};
}

=head2 query_description

 Title   : query_description
 Usage   : my $acc = $parser->query_description();
 Function: Get/Set query (model) description; Infernal currently does not output
           the Rfam description
 Returns : String (description)
 Args    : [optional] String (description)

=cut

sub query_description {
    my $self = shift;
    return $self->{'_query_description'} = shift if @_;
    return $self->{'_query_description'};
}

=head2 hsp_minscore

 Title   : hsp_minscore
 Usage   : my $cutoff = $parser->hsp_minscore();
 Function: Get/Set min bit score cutoff (for generating Hits/HSPs)
 Returns : score (number)
 Args    : [optional] score (number)

=cut

sub hsp_minscore {
    my $self = shift;
    return $self->{'_hsp_minscore'} = shift if @_;
    return $self->{'_hsp_minscore'};
}

=head2 convert_meta

 Title   : convert_meta
 Usage   : $parser->convert_meta(1);
 Function: Get/Set boolean flag for converting Infernal WUSS format
           to a simple bracketed format (simple WUSS by default) 
 Returns : boolean flag (TRUE or FALSE)
 Args    : [optional] boolean (eval's to TRUE or FALSE)

=cut

sub convert_meta {
    my $self = shift;
    return $self->{'_convert_meta'} = shift if @_;
    return $self->{'_convert_meta'};
}

=head2 version

 Title   : version
 Usage   : $parser->version();
 Function: Set the Infernal cmsearch version
 Returns : version
 Args    : [optional] version

=cut

sub version {
    my $self = shift;
    return $self->{'_version'} = shift if @_;
    return $self->{'_version'};
}

=head2 structure_symbols

 Title   : structure_symbols
 Usage   : my $hashref = $parser->structure_symbols();
 Function: Get/Set RNA structure symbols
 Returns : Hash ref of delimiters (5' stem, 3' stem, single-strand, etc)
         : default = < (5-prime)
                     > (3-prime)
                     : (single-strand)
                     ? (unknown)
                     . (gap)
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

=head2 simple_meta

 Title   : simple_meta
 Usage   : my $string = $parser->simple_meta($str);
 Function: converts more complex WUSS meta format into simple bracket format
           using symbols defined in structure_symbols()
 Returns : converted string
 Args    : [required] string to convert
 Note    : This is a very simple conversion method to get simple bracketed
           format from Infernal data.  If the convert_meta() flag is set,
           this is the method used to convert the strings.

=cut

sub simple_meta {
    my ($self, $str) = @_;
    $self->throw("No string arg sent!") if !$str;
    my $structs = $self->structure_symbols();
    my ($ls, $rs, $ss, $unk, $gap) = ($structs->{'5-prime'}, $structs->{'3-prime'},
                                $structs->{'single-strand'}, $structs->{'unknown'},
                                $structs->{'gap'});
    $str =~ s{[\(\<\[\{]}{$ls}g;
    $str =~ s{[\)\>\]\}]}{$rs}g;
    $str =~ s{[:,_-]}{$ss}g;
    $str =~ s{\.}{$gap}g;
    # unknown not handled yet
    return $str;
}

## private methods

# this is a hack which guesses the format and sets the handler for parsing in
# an instance; it'll be taken out when infernal 1.0 is released

sub _parse_latest {
    my ($self) = @_;
    my $seentop = 0;
    local $/ = "\n";
    my ($accession, $description) = ($self->query_accession, $self->query_description);
    my ($maxscore, $mineval, $minpval);
    $self->start_document();
    my ($lasthit, $lastscore, $lasteval, $lastpval, $laststart, $lastend);
    PARSER:
    while (my $line = $self->_readline) {
        next if $line =~ m{^\s+$};
        # stats aren't parsed yet...
		if ($line =~ m{^\#\s+cmsearch}xms) {
			$seentop = 1;
			$self->start_element({'Name' => 'Result'});
			$self->element_hash({
					'Infernal_program'   => 'CMSEARCH'
				});
		}
		elsif ($line =~ m{^\#\sINFERNAL\s+(\d+\.\d+)}xms) {
			$self->element_hash({
					'Infernal_version'   => $1,
				});
		}
		elsif ($line =~ m{^\#\scommand:.*?\s(\S+)$}xms) {
			$self->element_hash({
					'Infernal_db'   => $1,
				});
		}
		elsif ($line =~ m{^\#\s+dbsize\(Mb\):\s+(\d+\.\d+)}xms) {
			# store absolute DB length
			$self->element_hash({
					'Infernal_db-let'   => $1 * 1e6
				});		
		}		
		elsif ($line =~ m{^CM(?:\s(\d+))?:\s*(\S+)}xms) {
			# not sure, but it's possible single reports may contain multiple
			# models; if so, they should be rolled over into a new ResultI
			#print STDERR "ACC: $accession\nDESC: $description\n";
			$self->element_hash({
			        'Infernal_query-def' => $2, # present in output now
			        'Infernal_query-acc' => $accession,
			        'Infernal_querydesc' => $description
			    });
        }
		elsif ($line =~ m{^>\s*(\S+)} ){
            #$self->debug("Start Hit: Found hit:$1\n");
            if ($self->in_element('hit')) {
                $self->element_hash({'Hit_score' => $maxscore,
                                     'Hit_bits'  => $maxscore});
                ($maxscore, $minpval, $mineval) = undef;
                $self->end_element({'Name' => 'Hit'});
            }
            $lasthit = $1;
        }
        elsif ($line =~ m{
            ^\sQuery\s=\s\d+\s-\s\d+,\s   # Query start/end
            Target\s=\s(\d+)\s-\s(\d+)    # Target start/end
            }xmso) {
            # Query (model) start/end always the same, determined from
            # the HSP length
            ($laststart, $lastend) = ($1, $2);
            #$self->debug("Found hit coords:$laststart - $lastend\n");
        } elsif ($line =~ m{
            ^\sScore\s=\s([\d\.]+),\s   # Score = Bitscore (for now)
            (?:E\s=\s([\d\.e-]+),\s     # E-val optional
             P\s=\s([\d\.e-]+),\s)?     # P-val optional
            GC\s=                       # GC not captured
            }xmso
                 ) {
            ($lastscore, $lasteval, $lastpval) = ($1, $2, $3);
            #$self->debug(sprintf("Found hit data:Score:%s,Eval:%s,Pval:%s\n",$lastscore, $lasteval||'', $lastpval||''));
            $maxscore ||= $lastscore;
            if ($lasteval && $lastpval) {
                $mineval ||= $lasteval;
                $minpval ||= $lastpval;
                $mineval = ($mineval > $lasteval)  ? $lasteval :
                        $mineval;
                $minpval = ($minpval > $lastpval)  ? $lastpval :
                        $minpval;
            }
            $maxscore = ($maxscore < $lastscore)  ? $lastscore :
                        $maxscore;
            if (!$self->within_element('hit')) {
                my ($gi, $acc, $ver) = $self->_get_seq_identifiers($lasthit);
                $self->start_element({'Name' => 'Hit'});
                $self->element_hash({
                    'Hit_id'           => $lasthit,
                    'Hit_accession'    => $ver ? "$acc.$ver" :
                                           $acc ? $acc : $lasthit,
                    'Hit_gi'           => $gi
                    });
            }
            if (!$self->in_element('hsp')) {
                $self->start_element({'Name' => 'Hsp'});
            }
            
            # hsp is similar to older output
        } elsif ($line =~ m{^(\s+)[<>\{\}\(\)\[\]:_,-\.]+}xms) { # start of HSP
            $self->_pushback($line); # set up for loop
            #$self->debug("Start HSP\n");
            # what is length of the gap to the structure data?
            my $offset = length($1);
            my ($ct, $strln) = 0;
            my $hsp;
            HSP:
            my %hsp_key = ('0' => 'meta',
               '1' => 'query',
               '2' => 'midline',
               '3' => 'hit');
            HSP:
            while (defined ($line = $self->_readline)) {
                chomp $line;
          		next if (!$line); # toss empty lines
                # next if $line =~ m{^\s*$}; # toss empty lines
                # it is possible to have homology lines consisting
                # entirely of spaces if the subject has a large
                # insertion where nothing matches the model
                
                # exit loop if at end of file or upon next hit/HSP
                if ($line =~ m{^\s{0,2}\S+}) {
                    $self->_pushback($line);
                    last HSP;
                }
                # iterate to keep track of each line (4 lines per hsp block)
                my $iterator = $ct % 4;
                # strlen set only with structure lines (proper length)
                $strln = length($line) if $iterator == 0;
                # only grab the data needed (hit start and stop in hit line above)
												
                my $data = substr($line, $offset, $strln-$offset);
                $hsp->{ $hsp_key{$iterator} } .= $data;
                
                $ct++;
            }
            
            # query start, end are from the actual query length (entire hit is
            # mapped to CM data, so all CM data is represented)
            # works for now...
            if ($self->in_element('hsp')) {				
				# In some cases with HSPs unaligned residues are present in
				# the hit or query (Ex: '*[ 8]*' is 8 unaligned residues).
				# This info needs to be passed on unmodifed to the HSP class
				# and handled there as it is subjectively changed based on
				# use.
                my $strlen = 0;
				
				# catch any insertions and add them into the actual length
				while ($hsp->{'query'} =~ m{\*\[\s*(\d+)\s*\]\*}g) {
					$strlen += $1;
				}
				# add on the actual residues
				$strlen += $hsp->{'query'} =~ tr{A-Za-z}{A-Za-z};
				
                my $metastr = ($self->convert_meta) ? ($self->simple_meta($hsp->{'meta'})) :
                            $hsp->{'meta'};
                $self->element_hash(
                               {'Hsp_stranded'  => 'HIT',
                                'Hsp_qseq'      => $hsp->{'query'},
                                'Hsp_hseq'      => $hsp->{'hit'},
                                'Hsp_midline'   => $hsp->{'midline'},
                                'Hsp_structure' => $metastr,
                                'Hsp_query-from' => 1,
                                'Infernal_query-len' => $strlen,
                                'Hsp_query-to'   => $strlen,
                                'Hsp_hit-from'  => $laststart,
                                'Hsp_hit-to'    => $lastend,
                                'Hsp_score'     => $lastscore,
                                'Hsp_bit-score' => $lastscore,
                            });
                $self->element_hash(
                               {'Hsp_evalue'    => $lasteval,
                                'Hsp_pvalue'    => $lastpval,
                            }) if ($lasteval && $lastpval);
                $self->end_element({'Name' => 'Hsp'});
            }
        # result now ends with // and 'Fin'
        } elsif ($line =~ m{^//}xms )  {
            if ($self->within_element('result') && $seentop) {
                if ($self->in_element('hit')) {
                    $self->element_hash({'Hit_score'    => $maxscore,
                                         'Hit_bits'     => $maxscore});
                    # don't know where to put minpval yet
                    $self->element_hash({'Hit_signif'   => $mineval}) if $mineval;
                    $self->element_hash({'Hit_p'        => $minpval}) if $minpval;
                    $self->end_element({'Name' => 'Hit'});
                }
                last PARSER;
            }
        }
    }
    $self->within_element('hit') && $self->end_element( { 'Name' => 'Hit' } );
    $self->end_element( { 'Name' => 'Result' } ) if $seentop;
    return $self->end_document();
}

# cmsearch 0.81 (pre-1.0)
sub _parse_pre {
    my ($self) = @_;
    my $seentop = 0;
    local $/ = "\n";
    my ($accession, $db, $algorithm, $description, $version) =
       ($self->query_accession, $self->database, $self->algorithm,
        $self->query_description, '0.81');
    my ($maxscore, $mineval, $minpval);
    $self->start_document();
    my ($lasthit, $lastscore, $lasteval, $lastpval, $laststart, $lastend);
    PARSER:
    while (my $line = $self->_readline) {
        next if $line =~ m{^\s+$};
        # stats aren't parsed yet...
        if ($line =~ m{CM\s\d+:\s*(\S+)}xms) {
            #$self->debug("Start Result: Found model:$1\n");
            if (!$self->within_element('result')) {
                $seentop = 1;
                $self->start_element({'Name' => 'Result'});
                $self->element_hash({
                        'Infernal_program'   => $algorithm,
                        'Infernal_query-def' => $1, # present in output now
                        'Infernal_query-acc' => $accession,
                        'Infernal_querydesc' => $description,
                        'Infernal_db'        => $db
                    });
            }
        } elsif ($line =~ m{^>\s*(\S+)} ){
            #$self->debug("Start Hit: Found hit:$1\n");
            if ($self->in_element('hit')) {
                $self->element_hash({'Hit_score' => $maxscore,
                                     'Hit_bits'  => $maxscore});
                ($maxscore, $minpval, $mineval) = undef;
                $self->end_element({'Name' => 'Hit'});
            }
            $lasthit = $1;
        }
        elsif ($line =~ m{
            ^\sQuery\s=\s\d+\s-\s\d+,\s   # Query start/end
            Target\s=\s(\d+)\s-\s(\d+)    # Target start/end
            }xmso) {
            # Query (model) start/end always the same, determined from
            # the HSP length
            ($laststart, $lastend) = ($1, $2);
            #$self->debug("Found hit coords:$laststart - $lastend\n");
        } elsif ($line =~ m{
            ^\sScore\s=\s([\d\.]+),\s   # Score = Bitscore (for now)
            (?:E\s=\s([\d\.e-]+),\s     # E-val optional
             P\s=\s([\d\.e-]+),\s)?     # P-val optional
            GC\s=                       # GC not captured
            }xmso
                 ) {
            ($lastscore, $lasteval, $lastpval) = ($1, $2, $3);
            #$self->debug(sprintf("Found hit data:Score:%s,Eval:%s,Pval:%s\n",$lastscore, $lasteval||'', $lastpval||''));
            $maxscore ||= $lastscore;
            if ($lasteval && $lastpval) {
                $mineval ||= $lasteval;
                $minpval ||= $lastpval;
                $mineval = ($mineval > $lasteval)  ? $lasteval :
                        $mineval;
                $minpval = ($minpval > $lastpval)  ? $lastpval :
                        $minpval;
            }
            $maxscore = ($maxscore < $lastscore)  ? $lastscore :
                        $maxscore;
            if (!$self->within_element('hit')) {
                my ($gi, $acc, $ver) = $self->_get_seq_identifiers($lasthit);
                $self->start_element({'Name' => 'Hit'});
                $self->element_hash({
                    'Hit_id'           => $lasthit,
                    'Hit_accession'    => $ver ? "$acc.$ver" :
                                           $acc ? $acc : $lasthit,
                    'Hit_gi'           => $gi
                    });
            }
            if (!$self->in_element('hsp')) {
                $self->start_element({'Name' => 'Hsp'});
            }
            
            # hsp is similar to older output
        } elsif ($line =~ m{^(\s+)[<>\{\}\(\)\[\]:_,-\.]+}xms) { # start of HSP
            $self->_pushback($line); # set up for loop
            #$self->debug("Start HSP\n");
            # what is length of the gap to the structure data?
            my $offset = length($1);
            my ($ct, $strln) = 0;
            my $hsp;
            HSP:
            my %hsp_key = ('0' => 'meta',
               '1' => 'query',
               '2' => 'midline',
               '3' => 'hit');
            HSP:
            while (defined ($line = $self->_readline)) {
                chomp $line;
          		next if (!$line); # toss empty lines
                # next if $line =~ m{^\s*$}; # toss empty lines
                # it is possible to have homology lines consisting
                # entirely of spaces if the subject has a large
                # insertion where nothing matches the model
                
                # exit loop if at end of file or upon next hit/HSP
                if ($line =~ m{^\s{0,2}\S+}) {
                    $self->_pushback($line);
                    last HSP;
                }
                # iterate to keep track of each line (4 lines per hsp block)
                my $iterator = $ct%4;
                # strlen set only with structure lines (proper length)
                $strln = length($line) if $iterator == 0;
                # only grab the data needed (hit start and stop in hit line above)
												
                my $data = substr($line, $offset, $strln-$offset);
                $hsp->{ $hsp_key{$iterator} } .= $data;
                
                $ct++;
            }
            
            # query start, end are from the actual query length (entire hit is
            # mapped to CM data, so all CM data is represented)
            # works for now...
            if ($self->in_element('hsp')) {
                my $strlen = $hsp->{'query'} =~ tr{A-Za-z}{A-Za-z};
                
                my $metastr;
                $metastr = ($self->convert_meta) ? ($self->simple_meta($hsp->{'meta'})) :
                            ($hsp->{'meta'});
                $self->element_hash(
                               {'Hsp_stranded'  => 'HIT',
                                'Hsp_qseq'      => $hsp->{'query'},
                                'Hsp_hseq'      => $hsp->{'hit'},
                                'Hsp_midline'   => $hsp->{'midline'},
                                'Hsp_structure' => $metastr,
                                'Hsp_query-from' => 1,
                                'Infernal_query-len' => $strlen,
                                'Hsp_query-to'   => $strlen,
                                'Hsp_hit-from'  => $laststart,
                                'Hsp_hit-to'    => $lastend,
                                'Hsp_score'     => $lastscore,
                                'Hsp_bit-score' => $lastscore,
                            });
                $self->element_hash(
                               {'Hsp_evalue'    => $lasteval,
                                'Hsp_pvalue'    => $lastpval,
                            }) if ($lasteval && $lastpval);
                $self->end_element({'Name' => 'Hsp'});
            }
        # result now ends with // and 'Fin'
        } elsif ($line =~ m{^//}xms )  {
            if ($self->within_element('result') && $seentop) {
                $self->element(
                            {'Name' => 'Infernal_version',
                             'Data' => $version}
                            );
                if ($self->in_element('hit')) {
                    $self->element_hash({'Hit_score'    => $maxscore,
                                         'Hit_bits'     => $maxscore});
                    # don't know where to put minpval yet
                    $self->element_hash({'Hit_signif'   => $mineval}) if $mineval; 
                    $self->end_element({'Name' => 'Hit'});
                }
                last PARSER;
            }
        }
    }
    $self->within_element('hit') && $self->end_element( { 'Name' => 'Hit' } );
    $self->end_element( { 'Name' => 'Result' } ) if $seentop;
    return $self->end_document();
}

# cmsearch 0.72 and below; will likely be dropped when Infernal 1.0 is released
sub _parse_old {
    my ($self) = @_;
    my $seentop = 0;
    local $/ = "\n";
    my ($accession, $db, $algorithm, $model, $description, $version) =
       ($self->query_accession, $self->database, $self->algorithm,
        $self->model, $self->query_description, $self->version);
    my $maxscore;
    my $cutoff = $self->hsp_minscore;
    $self->start_document();
    local ($_);
    my $line;
    my ($lasthit, $lastscore, $laststart, $lastend);
    my $hitline;
    PARSER:
    while ( defined( $line = $self->_readline ) ) {
        next if $line =~ m{^\s+$};
        # bypass this for now...
        next if $line =~ m{^HMM\shit};
        # pre-0.81
        if ($line =~ m{^sequence:\s+(\S+)} ){
            if (!$self->within_element('result')) {
                $seentop = 1;
                $self->start_element({'Name' => 'Result'});
                $self->element_hash({
                        'Infernal_program'   => $algorithm,
                        'Infernal_query-def' => $model,
                        'Infernal_query-acc' => $accession,
                        'Infernal_querydesc' => $description,
                        'Infernal_db'        => $db
                    });
            }
            if ($self->in_element('hit')) {
                $self->element_hash({'Hit_score' => $maxscore,
                                     'Hit_bits'  => $maxscore});
                $maxscore = undef;
                $self->end_element({'Name' => 'Hit'});
            }            
            $lasthit = $1;
        } elsif ($line =~ m{^hit\s+\d+\s+:\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+bits}xms) {
            ($laststart, $lastend, $lastscore) = ($1, $2, $3);
            $maxscore = $lastscore unless $maxscore;
            if ($lastscore > $cutoff) {
                if (!$self->within_element('hit')) {
                    my ($gi, $acc, $ver) = $self->_get_seq_identifiers($lasthit);
                    $self->start_element({'Name' => 'Hit'});
                    $self->element_hash({
                        'Hit_id'           => $lasthit,
                        'Hit_accession'    => $ver ? "$acc.$ver" :
                                               $acc ? $acc : $lasthit,
                        'Hit_gi'           => $gi
                        });
                }
                # necessary as infernal 0.71 has repeated hit line
                if (!$self->in_element('hsp')) {
                    $self->start_element({'Name' => 'Hsp'});
                }
                $maxscore = ($maxscore < $lastscore)  ? $lastscore :
                            $maxscore;
            }
        } elsif ($line =~ m{^(\s+)[<>\{\}\(\)\[\]:_,-\.]+}xms) { # start of HSP
            $self->_pushback($line); # set up for loop
            # what is length of the gap to the structure data?
            my $offset = length($1);
            my ($ct, $strln) = 0;
            my $hsp;
            HSP:
            my %hsp_key = ('0' => 'meta',
               '1' => 'query',
               '2' => 'midline',
               '3' => 'hit');
            HSP:
            while ($line = $self->_readline) {
                next if $line =~ m{^\s*$}; # toss empty lines
                chomp $line;
                # exit loop if at end of file or upon next hit/HSP
                if (!defined($line) || $line =~ m{^\S+}) {
                    $self->_pushback($line);
                    last HSP;
                }
                # iterate to keep track of each line (4 lines per hsp block)
                my $iterator = $ct%4;
                # strlen set only with structure lines (proper length)
                $strln = length($line) if $iterator == 0;
                # only grab the data needed (hit start and stop in hit line above)

                my $data = substr($line, $offset, $strln-$offset);
                $hsp->{ $hsp_key{$iterator} } .= $data;
                $ct++;
            }
            # query start, end are from the actual query length (entire hit is
            # mapped to CM data, so all CM data is represented)
            # works for now...
            if ($self->in_element('hsp')) {
                my $strlen = $hsp->{'query'} =~ tr{A-Za-z}{A-Za-z};
                
		my $metastr;
                # Ugh...these should be passed in a hash
                $metastr = ($self->convert_meta) ? ($self->simple_meta($hsp->{'meta'})) :
                            ($hsp->{'meta'});
                $self->element_hash(
                               {'Hsp_stranded'  => 'HIT',
                                'Hsp_qseq'      => $hsp->{'query'},
                                'Hsp_hseq'      => $hsp->{'hit'},
                                'Hsp_midline'   => $hsp->{'midline'},
                                'Hsp_structure' => $metastr,
                                'Hsp_query-from' => 1,
                                'Infernal_query-len' => $strlen,
                                'Hsp_query-to'   => $strlen,
                                'Hsp_hit-from'  => $laststart,
                                'Hsp_hit-to'    => $lastend,
                                'Hsp_score'     => $lastscore,
                                'Hsp_bit-score' => $lastscore
                            });
                $self->end_element({'Name' => 'Hsp'});
            }
        } elsif ($line =~ m{^memory}xms || $line =~ m{^CYK\smemory}xms )  {
            if ($self->within_element('result') && $seentop) {
                $self->element(
                            {'Name' => 'Infernal_version',
                             'Data' => $version}
                            );
                if ($self->in_element('hit')) {
                    $self->element_hash({'Hit_score'    => $maxscore,
                                         'Hit_bits'     => $maxscore});
                    $self->end_element({'Name' => 'Hit'});
                }
                last PARSER;
            }
        }
    }
    $self->within_element('hit') && $self->end_element( { 'Name' => 'Hit' } );
    $self->end_element( { 'Name' => 'Result' } ) if $seentop;
    return $self->end_document();
}

1;
