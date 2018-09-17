#
# BioPerl module for Bio::SearchIO::erpin
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

Bio::SearchIO::erpin - SearchIO-based ERPIN parser

=head1 SYNOPSIS

  # do not call this module directly. Use Bio::SearchIO.

=head1 DESCRIPTION

This is an experimental SearchIO-based parser for output from
the erpin program.  It currently parses erpin output for ERPIN
versions 4.2.5 and above; older versions may work but will not be supported.

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

package Bio::SearchIO::erpin;
use strict;

use Data::Dumper;
use base qw(Bio::SearchIO);

my %MODEMAP = (
	    'Result'             => 'result',
	    'Hit'                => 'hit',
	    'Hsp'                => 'hsp'
	    );

my %MAPPING = ( 
        'Hsp_bit-score'   => 'HSP-bits',
        'Hsp_score'       => 'HSP-score',
        'Hsp_evalue'      => 'HSP-evalue', # no evalues yet
        'Hsp_query-from'  => 'HSP-query_start',
        'Hsp_query-to'    => 'HSP-query_end',
        'Hsp_hit-from'    => 'HSP-hit_start', #
        'Hsp_hit-to'      => 'HSP-hit_end', #
        'Hsp_gaps'        => 'HSP-hsp_gaps', 
        'Hsp_hitgaps'     => 'HSP-hit_gaps',
        'Hsp_querygaps'   => 'HSP-query_gaps',
        'Hsp_qseq'        => 'HSP-query_seq',
        'Hsp_hseq'        => 'HSP-hit_seq',
        'Hsp_midline'     => 'HSP-homology_seq',
        'Hsp_structure'   => 'HSP-meta',
        'Hsp_align-len'   => 'HSP-hsp_length',
        'Hsp_stranded'    => 'HSP-stranded',
        
        # not supported yet
        'Hsp_positive'    => 'HSP-conserved',
        'Hsp_identity'    => 'HSP-identical',

        'Hit_id'        => 'HIT-name',
        'Hit_len'       => 'HIT-length',
        'Hit_gi'        => 'HIT-ncbi_gi',
        'Hit_accession' => 'HIT-accession',
        'Hit_def'       => 'HIT-description',
        'Hit_signif'    => 'HIT-significance', # none yet
        'Hit_score'     => 'HIT-score', # best HSP bit score
        'Hit_bits'      => 'HIT-bits', # best HSP bit score
 
        'ERPIN_program'  => 'RESULT-algorithm_name', # get/set 
        'ERPIN_version'  => 'RESULT-algorithm_version', # get/set 
        'ERPIN_query-def'=> 'RESULT-query_name', # get/set 
        'ERPIN_query-len'=> 'RESULT-query_length', 
        'ERPIN_query-acc'=> 'RESULT-query_accession', # get/set 
        'ERPIN_querydesc'=> 'RESULT-query_description', # get/set
        'ERPIN_db'       => 'RESULT-database_name',  # get/set 
        'ERPIN_db-len'   => 'RESULT-database_entries', # none yet
        'ERPIN_db-let'   => 'RESULT-database_letters', # none yet
        
        'Parameters_cutoff'      => { 'RESULT-parameters' => 'cutoff' },
        'Parameters_expect'      => { 'RESULT-parameters' => 'expect' },
        
        'Parameters_include'     => { 'RESULT-parameters' => 'include' },
        'Parameters_sc-match'    => { 'RESULT-parameters' => 'match' },
        'Parameters_sc-mismatch' => { 'RESULT-parameters' => 'mismatch' },
        'Parameters_gap-open'    => { 'RESULT-parameters' => 'gapopen' },
        'Parameters_gap-extend'  => { 'RESULT-parameters' => 'gapext' },
        'Parameters_filter'      => { 'RESULT-parameters' => 'filter' },
        'Parameters_allowgaps'   => { 'RESULT-parameters' => 'allowgaps' },
        'Parameters_full_dbpath' => { 'RESULT-parameters' => 'full_dbpath' },        
        'Statistics_db-let'      => { 'RESULT-statistics' => 'dbletters' },
	     );

my $MINSCORE = 0;
my $DEFAULT_VERSION = '4.2.5';
my $DEFAULT_ALGORITHM = 'erpin';

=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::infernal->new();
 Function: Builds a new Bio::SearchIO::infernal object 
 Returns : Bio::SearchIO::infernal
 Args    : -fh/-file     => cmsearch (infernal) filename
           -format       => 'erpin'
           -algorithm    => algorithm (default 'Infernal')
           -query_acc    => query accession, eg. Rfam accession (default undef)
           -hsp_minscore => minimum HSP score cutoff
           -version      => ERPIN version (not reported in output)

=cut

sub _initialize {
    my ( $self, @args ) = @_;
    $self->SUPER::_initialize(@args);
    my ($cutoff, $accession, $version) =
       $self->_rearrange([qw(HSP_MINSCORE QUERY_ACC VERSION)],@args);
    my $handler = $self->_eventHandler;
    $handler->register_factory(
        'result',
        Bio::Factory::ObjectFactory->new(
            -type      => 'Bio::Search::Result::GenericResult',
            -interface => 'Bio::Search::Result::ResultI',
            -verbose   => $self->verbose()
        )
    );

    $handler->register_factory(
        'hit',
        Bio::Factory::ObjectFactory->new(
            -type      => 'Bio::Search::Hit::ModelHit',
            -interface => 'Bio::Search::Hit::HitI',
            -verbose   => $self->verbose()
        )
    );

    $handler->register_factory(
        'hsp',
        Bio::Factory::ObjectFactory->new(
            -type      => 'Bio::Search::HSP::ModelHSP',
            -interface => 'Bio::Search::HSP::HSPI',
            -verbose   => $self->verbose()
        )
    );
    $accession && $self->query_accession($accession);
    $cutoff ||= $MINSCORE;
    $self->hsp_minscore($cutoff);
    $version ||= $DEFAULT_VERSION;
    $self->algorithm_version($version);
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
    my $accession = $self->query_accession;
    my $minscore = $self->hsp_minscore;
    my $version = $self->algorithm_version;
    my $verbose = $self->verbose;    # cache for speed?
    $self->start_document();
    my ($lasthit, $lastscore, $lastlen, $lasteval);
    #my $hitline;
    PARSER:
    while ( defined( my $line = $self->_readline ) ) {
        next if $line =~ m{^\s*$};
        if ($line =~ m{^Training\sset:\s+"(.*)"}xmso) {
            if ($seentop) {
                $self->_pushback($line);
                last PARSER;
            }
            $self->start_element({'Name' => 'Result'});
            $self->element_hash( {
                'ERPIN_query-def'   => $1,
                'ERPIN_program'     =>'erpin',
                'ERPIN_version'     => $version,
                'ERPIN_query-acc'   => $accession,
                });
            $seentop = 1;
            # parse rest of header here
            HEADER:
            while (defined ($line = $self->_readline) ) {
                next if $line =~ m{^\s*$};
                if (index($line, '>') == 0 ||
                    index($line, '-------- at level 1 --------') == 0) {
                    $self->_pushback($line);
                    last HEADER;
                }
                if ($line =~ m{^\s+(\d+\ssequences\sof\slength\s\d+)}xmso) {
                    $self->element(
                        {'Name' => 'ERPIN_querydesc',
                         'Data' => $1}
                       );
                } elsif ($line =~ m{^Cutoff:\s+(\S+)}xmso) {
                    $self->element(
                        {'Name' => 'Parameters_cutoff',
                         'Data' => $1}
                                  );
                } elsif ($line =~ m{^Database:\s+"(.*)"}xmso) {
                    $self->element(
                        {'Name' => 'ERPIN_db',
                         'Data' => $1}
                       );
                } elsif ($line =~ m{^\s+(\d+)\snucleotides\sto\sbe\sprocessed\sin\s(\d+)\ssequences}xmso) {
                    $self->element_hash(
                        {'ERPIN_db-len' => $2,
                         'ERPIN_db-let' => $1}
                       );
                } elsif ($line =~ m{^E-value\sat\scutoff\s\S+\sfor\s\S+\sdouble\sstrand\sdata:\s+(\S+)}xmso) {
                    $self->element(
                                   {'Name' => 'Parameters_expect',
                                    'Data' => $1}
                                  );
                } elsif ($line =~ m{^\s+(ATGC\sratios:\s+(?:\S+\s+\S+\s+\S+\s+\S+))}) {
                    $self->element(
                                   {'Name' => 'Statistics_db-let',
                                    'Data' => $1}
                                  );
                }
            }
        } elsif ($line =~ m{^>(\S+)\s+(.*)}xmso ) {
            my ($id, $desc) = ($1, $2);
            chomp $desc;
            # desc line is repeated for each strand, so must check
            # prior to starting a new hit
            if (!$lasthit || $id ne $lasthit) {
                if ($self->within_element('hit') ) {
                    $self->element_hash({
                        'Hit_signif' => $lasteval,
                        'Hit_score'  => $lastscore,
                        'Hit_bits'   => $lastscore
                        });
                    $self->end_element({'Name' => 'Hit'});
                }                
                $self->start_element({'Name' => 'Hit'});
                my ($gi, $acc, $ver) = $self->_get_seq_identifiers($id);
            
                $self->element_hash({
                    'Hit_id'        => $id,
                    'Hit_gi'        => $gi,
                    'Hit_accession' => $ver ? "$acc.$ver" :
                                        $acc ? $acc : $id,
                    'Hit_def'       => $desc
                    });
            }
            $lasthit = $id;
        } elsif ( (index($line, 'FW') == 0) || (index($line, 'RC') == 0)) {
            my ($str, $hn, $pos, $score, $eval) = split ' ', $line;
            if ($minscore < $score) {
                $self->start_element({'Name' => 'Hsp'});
                
                my ($start, $end) = split m{\.\.}, $pos, 2;
                ($start, $end) = ($end, $start) if ($str eq 'RC');
                $line = $self->_readline;
                chomp $line;
                $self->element_hash({
                    'Hsp_stranded'     => 'HIT',
                    'Hsp_hit-from'     => $start,
                    'Hsp_hit-to'       => $end,
                    'Hsp_score'        => $score,
                    'Hsp_bit-score'    => $score,
                    'Hsp_evalue'       => $eval,
                    'Hsp_query-from'   => 1,
                    'Hsp_query-to'     => length($line),
                    'Hsp_align-len'    => length($line),
                    'Hsp_hseq'         =>$line
                    });
                $self->end_element({'Name' => 'Hsp'});
                $lastscore = $score if (!$lastscore || $lastscore < $score);
                $lasteval = $eval if (!$lasteval || $lasteval > $eval);
            }
        } else {
            #$self->debug("Dropped data: $line");
        }
    }
    if ($seentop) {
        if ($self->within_element('hit')) {
            $self->element_hash({
                    'Hit_signif'    => $lasteval,
                    'Hit_score'     => $lastscore,
                    'Hit_bits'      => $lastscore
                    });
            $self->end_element({'Name' => 'Hit'}); 
        }
        $self->end_element({'Name' => 'Result'});
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

1;
