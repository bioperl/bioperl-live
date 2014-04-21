#
# BioPerl module for Bio::SearchIO::XML::BlastHandler
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich, Chris Fields
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::XML::BlastHandler - XML Handler for NCBI Blast XML parsing.

=head1 SYNOPSIS

  # This is not to be used directly.

=head1 DESCRIPTION

This is the XML handler for BLAST XML parsing. Currently it passes elements off
to the event handler, which is ultimately responsible for Bio::Search object
generation.

This was recently split off from the original code for Bio::SearchIO::blastxml
primarily for maintenance purposes.

=head1 DEPENDENCIES

In addition to parts of the Bio:: hierarchy, this module uses:

 XML::SAX::Base

which comes with the XML::SAX distribution.

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

=head1 AUTHOR - Jason Stajich, Chris Fields

Email jason-at-bioperl.org
Email cjfields-at-uiuc dot edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::SearchIO::XML::BlastHandler;
use base qw(Bio::Root::Root XML::SAX::Base);

my %MODEMAP = (
                'Iteration'   => 'result',
                'Hit'         => 'hit',
                'Hsp'         => 'hsp'
);

# major post 2.2.12 BLAST XML changes
# 1) moved XML Handler to it's own class
# 2) reconfigure blastxml to deal with old and new BLAST XML output

my %MAPPING = (
                # Result-specific fields
                'BlastOutput_program'   => 'RESULT-algorithm_name',
                'BlastOutput_version'   => 'RESULT-algorithm_version',
                'BlastOutput_db'        => 'RESULT-database_name',
                'BlastOutput_reference' => 'RESULT-program_reference',
                'BlastOutput_query-def' => 'RESULT-query_description',
                'BlastOutput_query-len' => 'RESULT-query_length',
                'BlastOutput_query-ID'  => 'runid',                
                'Parameters_matrix'     => { 'RESULT-parameters' => 'matrix'},
                'Parameters_expect'     => { 'RESULT-parameters' => 'expect'},
                'Parameters_include'    => { 'RESULT-parameters' => 'include'},
                'Parameters_sc-match'   => { 'RESULT-parameters' => 'match'},
                'Parameters_sc-mismatch' => { 'RESULT-parameters' => 'mismatch'},
                'Parameters_gap-open'   => { 'RESULT-parameters' => 'gapopen'},
                'Parameters_gap-extend' => { 'RESULT-parameters' => 'gapext'},
                'Parameters_filter'     => {'RESULT-parameters' => 'filter'},
                'Statistics_db-num'     => 'RESULT-database_entries',
                'Statistics_db-len'     => 'RESULT-database_letters',
                'Statistics_hsp-len'    => { 'RESULT-statistics' => 'hsplength'},
                'Statistics_eff-space'  => { 'RESULT-statistics' => 'effectivespace'},
                'Statistics_kappa'      => { 'RESULT-statistics' => 'kappa' },
                'Statistics_lambda'     => { 'RESULT-statistics' => 'lambda' },
                'Statistics_entropy'    => { 'RESULT-statistics' => 'entropy'},
                
                # HSP specific fields
                'Hsp_bit-score'  => 'HSP-bits',
                'Hsp_score'      => 'HSP-score',
                'Hsp_evalue'     => 'HSP-evalue',
                'Hsp_query-from' => 'HSP-query_start',
                'Hsp_query-to'   => 'HSP-query_end',
                'Hsp_hit-from'   => 'HSP-hit_start',
                'Hsp_hit-to'     => 'HSP-hit_end',
                'Hsp_positive'   => 'HSP-conserved',
                'Hsp_identity'   => 'HSP-identical',
                'Hsp_gaps'       => 'HSP-gaps',
                'Hsp_hitgaps'    => 'HSP-hit_gaps',
                'Hsp_querygaps'  => 'HSP-query_gaps',
                'Hsp_qseq'       => 'HSP-query_seq',
                'Hsp_hseq'       => 'HSP-hit_seq',
                'Hsp_midline'    => 'HSP-homology_seq',
                'Hsp_align-len'  => 'HSP-hsp_length',
                'Hsp_query-frame'=> 'HSP-query_frame',
                'Hsp_hit-frame'  => 'HSP-hit_frame',

                # Hit specific fields
                'Hit_id'               => 'HIT-name',
                'Hit_len'              => 'HIT-length',
                'Hit_accession'        => 'HIT-accession',
                'Hit_def'              => 'HIT-description',
                'Hit_num'              => 'HIT-order',
                'Iteration_iter-num'   => 'HIT-iteration',
                'Iteration_stat'       => 'HIT-iteration_statistic',
                
                # if these tags are present, they will overwrite the
                # above with more current data (i.e. multiquery hits)
                'Iteration_query-def'   => 'RESULT-query_description',
                'Iteration_query-len'   => 'RESULT-query_length',       
                'Iteration_query-ID'    => 'runid',
               );

# these XML tags are ignored for now
my %IGNOREDTAGS = (
                'Hsp_num'              => 1,#'HSP-order',
                'Hsp_pattern-from'     => 1,#'patternend',
                'Hsp_pattern-to'       => 1,#'patternstart',
                'Hsp_density'          => 1,#'hspdensity',
                'Iteration_message'    => 1,
                'Hit_hsps'             => 1,
                'BlastOutput_param'    => 1,
                'Iteration_hits'       => 1,
                'Statistics'           => 1,
                'Parameters'           => 1,
                'BlastOutput'          => 1,
                'BlastOutput_iterations' => 1,     
                   );

=head2 SAX methods

=cut

=head2 start_document

 Title   : start_document
 Usage   : $parser->start_document;
 Function: SAX method to indicate starting to parse a new document
 Returns : none
 Args    : none

=cut

sub start_document{
    my ($self) = @_;
    $self->{'_lasttype'} = '';
    $self->{'_values'} = {};
    $self->{'_result'}= [];
}

=head2 end_document

 Title   : end_document
 Usage   : $parser->end_document;
 Function: SAX method to indicate finishing parsing a new document
 Returns : Bio::Search::Result::ResultI object
 Args    : none

=cut

sub end_document{
   my ($self,@args) = @_;
   
   # reset data carried throughout parse
   $self->{'_resultdata'} = undef;
   
   # pass back ref to results queue; caller must reset handler results queue
   return $self->{'_result'};
}

=head2 start_element

 Title   : start_element
 Usage   : $parser->start_element($data)
 Function: SAX method to indicate starting a new element
 Returns : none
 Args    : hash ref for data

=cut

sub start_element{
    my ($self,$data) = @_;
    # we currently don't care about attributes
    my $nm = $data->{'Name'};

    if( my $type = $MODEMAP{$nm} ) {
        if( $self->eventHandler->will_handle($type) ) {
            my $func = sprintf("start_%s",lc $type);
            $self->eventHandler->$func($data->{'Attributes'});
        }                                                    
    }
}

=head2 end_element

 Title   : end_element
 Usage   : $parser->end_element($data)
 Function: Signals finishing an element
 Returns : Bio::Search object dpending on what type of element
 Args    : hash ref for data

=cut

sub end_element{
    my ($self,$data) = @_;

    my $nm = $data->{'Name'};
    my $rc;
    if($nm eq 'BlastOutput_program' &&
       $self->{'_last_data'} =~ /(t?blast[npx])/i ) {
        $self->{'_type'} = uc $1; 
    }
    if ($nm eq 'Iteration') {
        map {
            $self->{'_values'}->{$_} = $self->{'_resultdata'}->{$_};
            } keys %{ $self->{'_resultdata'} };
    }
    if( my $type = $MODEMAP{$nm} ) {
        if( $self->eventHandler->will_handle($type) ) {
            my $func = sprintf("end_%s",lc $type);
            $rc = $self->eventHandler->$func($self->{'_type'},
                                              $self->{'_values'});
        }
    }
    elsif( exists $MAPPING{$nm} ) { 
        if ( ref($MAPPING{$nm}) =~ /hash/i ) {
            my $key = (keys %{$MAPPING{$nm}})[0];
            $self->{'_values'}->{$key}->{$MAPPING{$nm}->{$key}} = $self->{'_last_data'};
        } else {
            $self->{'_values'}->{$MAPPING{$nm}} = $self->{'_last_data'};
        }
    }
    elsif( exists $IGNOREDTAGS{$nm} ){
        # ignores these elements for now
    }
    else {      
        $self->debug("ignoring unrecognized element type $nm\n");
    }
    $self->{'_last_data'} = ''; # remove read data if we are at 
                                # end of an element
                                
    # add to ResultI array
    $self->{'_result'} = $rc if( $nm eq 'Iteration' );
    # reset values for each Result round
    if ($nm eq 'Iteration') {
        $self->{'_values'} = {};
    }
}

=head2 characters

 Title   : characters
 Usage   : $parser->characters($data)
 Function: Signals new characters to be processed
 Returns : characters read
 Args    : hash ref with the key 'Data'


=cut

sub characters{
   my ($self,$data) = @_;
   return unless ( defined $data->{'Data'} && $data->{'Data'} !~ /^\s+$/ );
   $self->{'_last_data'} .= $data->{'Data'};
}

sub eventHandler {
    my $self = shift;
    return $self->{'_handler'} = shift if @_;
    return $self->{'_handler'};
}

1;
