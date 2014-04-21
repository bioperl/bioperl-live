#
# BioPerl module for Bio::SearchIO::fasta
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::fasta - A SearchIO parser for FASTA results

=head1 SYNOPSIS

  # Do not use this object directly, use it through the SearchIO system
   use Bio::SearchIO;
   my $searchio = Bio::SearchIO->new(-format => 'fasta',
                    -file   => 'report.FASTA');
   while( my $result = $searchio->next_result ) {
    # ... do what you would normally doi with Bio::SearchIO.
   }

=head1 DESCRIPTION

This object contains the event based parsing code for FASTA format
reports.  It creates L<Bio::Search::HSP::FastaHSP> objects instead of
L<Bio::Search::HSP::GenericHSP> for the HSP objects. 

This module will parse -m 9 -d 0 output as well as default m 1 output
from FASTA as well as SSEARCH.

Also see the SearchIO HOWTO:
L<http://bioperl.open-bio.org/wiki/HOWTO:SearchIO>.

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

=head1 AUTHOR - Jason Stajich, Aaron Mackey, William Pearson

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SearchIO::fasta;
use vars qw(%MODEMAP %MAPPING $IDLENGTH);
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::Factory::ObjectFactory;

BEGIN {

    # Set IDLENGTH to a new value if you have
    # compile FASTA with a different ID length
    # (actually newest FASTA allows the setting of this
    #  via -C parameter, default is 6)
    $IDLENGTH = 6;

    # mapping of NCBI Blast terms to Bioperl hash keys
    %MODEMAP = (
        'FastaOutput' => 'result',
        'Hit'         => 'hit',
        'Hsp'         => 'hsp'
    );

    # This should really be done more intelligently, like with
    # XSLT

    %MAPPING = (
        'Hsp_bit-score'   => 'HSP-bits',
        'Hsp_score'       => 'HSP-score',
        'Hsp_sw-score'    => 'HSP-swscore',
        'Hsp_evalue'      => 'HSP-evalue',
        'Hsp_evalue2'     => 'HSP-evalue2',
        'Hsp_query-from'  => 'HSP-query_start',
        'Hsp_query-to'    => 'HSP-query_end',
        'Hsp_hit-from'    => 'HSP-hit_start',
        'Hsp_hit-to'      => 'HSP-hit_end',
        'Hsp_positive'    => 'HSP-conserved',
        'Hsp_identity'    => 'HSP-identical',
        'Hsp_gaps'        => 'HSP-hsp_gaps',
        'Hsp_hitgaps'     => 'HSP-hit_gaps',
        'Hsp_querygaps'   => 'HSP-query_gaps',
        'Hsp_qseq'        => 'HSP-query_seq',
        'Hsp_hseq'        => 'HSP-hit_seq',
        'Hsp_midline'     => 'HSP-homology_seq',
        'Hsp_align-len'   => 'HSP-hsp_length',
        'Hsp_query-frame' => 'HSP-query_frame',
        'Hsp_hit-frame'   => 'HSP-hit_frame',

        'Hit_id'        => 'HIT-name',
        'Hit_len'       => 'HIT-length',
        'Hit_accession' => 'HIT-accession',
        'Hit_def'       => 'HIT-description',
        'Hit_signif'    => 'HIT-significance',
        'Hit_score'     => 'HIT-score',

        'FastaOutput_program'   => 'RESULT-algorithm_name',
        'FastaOutput_version'   => 'RESULT-algorithm_version',
        'FastaOutput_query-def' => 'RESULT-query_name',
        'FastaOutput_querydesc' => 'RESULT-query_description',
        'FastaOutput_query-len' => 'RESULT-query_length',
        'FastaOutput_db'        => 'RESULT-database_name',
        'FastaOutput_db-len'    => 'RESULT-database_entries',
        'FastaOutput_db-let'    => 'RESULT-database_letters',

        'Parameters_matrix'      => { 'RESULT-parameters' => 'matrix' },
        'Parameters_expect'      => { 'RESULT-parameters' => 'expect' },
        'Parameters_include'     => { 'RESULT-parameters' => 'include' },
        'Parameters_sc-match'    => { 'RESULT-parameters' => 'match' },
        'Parameters_sc-mismatch' => { 'RESULT-parameters' => 'mismatch' },
        'Parameters_gap-open'    => { 'RESULT-parameters' => 'gapopen' },
        'Parameters_gap-ext'     => { 'RESULT-parameters' => 'gapext' },
        'Parameters_word-size'   => { 'RESULT-parameters' => 'wordsize' },
        'Parameters_ktup'        => { 'RESULT-parameters' => 'ktup' },
        'Parameters_filter'      => { 'RESULT-parameters' => 'filter' },
        'Statistics_db-num'      => { 'RESULT-statistics' => 'dbentries' },
        'Statistics_db-len'      => { 'RESULT-statistics' => 'dbletters' },
        'Statistics_hsp-len'     => { 'RESULT-statistics' => 'hsplength' },
        'Statistics_eff-space'   => { 'RESULT-statistics' => 'effectivespace' },
        'Statistics_kappa'       => { 'RESULT-statistics' => 'kappa' },
        'Statistics_lambda'      => { 'RESULT-statistics' => 'lambda' },
        'Statistics_entropy'     => { 'RESULT-statistics' => 'entropy' },
    );
}

use base qw(Bio::SearchIO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::fasta->new();
 Function: Builds a new Bio::SearchIO::fasta object 
 Returns : Bio::SearchIO::fasta
 Args    : -idlength - set ID length to something other 
                       than the default (6), this is only
                       necessary if you have compiled FASTA
                       with a new default id length to display
                       in the HSP alignment blocks

=cut

sub _initialize {
    my ( $self, @args ) = @_;
    $self->SUPER::_initialize(@args);
    return unless @args;
    my ($idlength) = $self->_rearrange( [qw(IDLENGTH)], @args );
    $self->idlength( $idlength || $IDLENGTH );
    $self->_eventHandler->register_factory(
        'hsp',
        Bio::Factory::ObjectFactory->new(
            -type      => 'Bio::Search::HSP::FastaHSP',
            -interface => 'Bio::Search::HSP::HSPI'
        )
    );
    return 1;
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
    local $/ = "\n";
    local $_;

    my $data    = '';
    my $seentop = 0;
    my $current_hsp;
    my $m9HSP = 0;
    $self->start_document();
    my @hit_signifs;
    while ( defined( $_ = $self->_readline ) ) {
        next if ( !$self->in_element('hsp')
            && /^\s+$/ );    # skip empty lines
        if (
               m/(\S+)\s+searches\s+a\s+(protein\s+or\s+DNA\s+)?sequence/oxi
            || /(\S+)\s+compares\s+a/
	    || /(\S+)\s+performs\s+a/
	    || /(\S+)\s+produces\s/
	    || /(\S+)\s+finds\s+/	# for lalign, but does not work because no "The best scores are:"
            || (   m/^\#\s+/	# has a command log line
                && ( $_ = $self->_readline )
                && /(\S+)\s+searches\s+a\s+(protein\s+or\s+DNA\s+)?sequence/oxi
		   || /(\S+)\s+compares\s+a/
		   || /(\S+)\s+performs\s+a/
		   || /(\S+)\s+produces\s/
		   || /(\S+)\s+finds\s+/	# for lalign, but does not work because no "The best scores are:"
	    )
          )
        {
            if ($seentop) {
                $self->_pushback($_);
                $self->end_element( { 'Name' => 'FastaOutput' } );
                return $self->end_document();
            }
            $self->{'_reporttype'} = $1;
            $self->start_element( { 'Name' => 'FastaOutput' } );
            $self->{'_result_count'}++;
            $seentop = 1;
            #$self->debug( "reporttype is " . $self->{'_reporttype'} . "\n" );
            $self->element(
                {
                    'Name' => 'FastaOutput_program',
                    'Data' => $self->{'_reporttype'}
                }
            );
            my $version;
            # version 35 version string on same line
            if (/version/) {
                ($version) = (/version\s+(\S+)/);
            }
            # earlier versions, it's on the next line
            else {
                $_ = $self->_readline();
                ($version) = (/version\s+(\S+)/);
            }
            $version = '' unless defined $version;
            $self->{'_version'} = $version;
            $self->element(
                {
                    'Name' => 'FastaOutput_version',
                    'Data' => $version
                }
            );

            my ( $last, $leadin, $type, $querylen, $querytype, $querydef );

            while ( defined( $_ = $self->_readline() ) ) {
                if (
                    /^ (
                       (?:\s+>) |             # fa33 lead-in
                       (?:\s*\d+\s*>>>)       # fa34 mlib lead-in
                      )
                      (.*)
                   /x
                  )
                {
                    ( $leadin, $querydef ) = ( $1, $2 );
                    if ( $leadin =~ m/>>>/ ) {
                        if ( $querydef =~
                            /^(.*?)\s+(?:\-\s+)?(\d+)\s+(aa|nt).*$/o )
                        {
                            ( $querydef, $querylen, $querytype ) =
                              ( $1, $2, $3 );
                            last;
                        }
                    }
                    else {
                        if ( $last =~ /(\S+)[:,]\s*(\d+)\s+(aa|nt)/ ) {
                            ( $querylen, $querytype ) = ( $2, $3 );
                            $querydef ||= $1;
                            last;
                        }
                    }
                }
                elsif (m/^\s*vs\s+\S+/o) {
                    if ( $last =~ /(\S+)[,:]\s+(\d+)\s+(aa|nt)/o ) {
                        ( $querydef, $querylen, $querytype ) = ( $1, $2, $3 );
                        last;
                    }
                }
                $last = $_;
            }
            if (   $self->{'_reporttype'}
                && $self->{'_reporttype'} eq 'FASTA' )
            {
                if ( $querytype eq 'nt' ) {
                    $self->{'_reporttype'} = 'FASTN';
                }
                elsif ( $querytype eq 'aa' ) {
                    $self->{'_reporttype'} = 'FASTP';
                }
            }
            my ( $name, $descr ) = $querydef =~ m/^(\S+)\s*(.*?)\s*$/o;
            $self->element(
                {
                    'Name' => 'FastaOutput_query-def',
                    'Data' => $name
                }
            );
            $self->element(
                {
                    'Name' => 'FastaOutput_querydesc',
                    'Data' => $descr
                }
            );
            if ($querylen) {
                $self->element(
                    {
                        'Name' => 'FastaOutput_query-len',
                        'Data' => $querylen
                    }
                );
            }
            else {
                $self->warn("unable to find and set query length");
            }
            if (
                   $last =~ /^\s*vs\s+(\S+)/
                || ( $last =~ /^searching\s+(\S+)\s+library/ )
                || ( $last =~ /^Library:\s+(\S+)\s+/ )
                || (
                    defined $_
                    && (   /^\s*vs\s+(\S+)/
                        || /^Library:\s+(\S+)\s+/ )
                )
                || ( defined( $_ = $self->_readline() )
                    && ( /^\s*vs\s+(\S+)/ || /^Library:\s+(\S+)/ ) )
              )
            {
                $self->element(
                    {
                        'Name' => 'FastaOutput_db',
                        'Data' => $1
                    }
                );
            }
            elsif (m/^\s+opt(?:\s+E\(\))?$/o) {

           # histogram ... read over it more rapidly than the larger outer loop:
                while ( defined( $_ = $self->_readline ) ) {
                    last if m/^>\d+/;
                }
            }
        }
        elsif (/(\d+)\s+residues\s+in\s+(\d+)\s+(?:library\s+)?sequences/) {
            $self->element(
                {
                    'Name' => 'FastaOutput_db-let',
                    'Data' => $1
                }
            );
            $self->element(
                {
                    'Name' => 'FastaOutput_db-len',
                    'Data' => $2
                }
            );
            $self->element(
                {
                    'Name' => 'Statistics_db-len',
                    'Data' => $1
                }
            );
            $self->element(
                {
                    'Name' => 'Statistics_db-num',
                    'Data' => $2
                }
            );
        }
        elsif (/Lambda=\s*(\S+)/) {
            $self->element(
                {
                    'Name' => 'Statistics_lambda',
                    'Data' => $1
                }
            );
        }
        elsif (/K=\s*(\S+)/) {
            $self->element(
                {
                    'Name' => 'Statistics_kappa',
                    'Data' => $1
                }
            );
        }
        elsif (/^\s*(Smith-Waterman)/) {

            $self->{'_reporttype'} = $1;

	    m/\[\s*(\S+)\s+matrix \([^\)]+\)(xS)?\],/;

            $self->element(
                {
                    'Name' => 'Parameters_matrix',
                    'Data' => $1
                }
            );
            $self->element(
                {
                    'Name' => 'Parameters_filter',
                    'Data' => defined $2 ? 1 : 0,
                }
            );
	    if (/\s*gap\-penalty:\s*(\-?\d+)\/(\-?\d+)/) {
		$self->element(
		    {
			'Name' => 'Parameters_gap-open',
			'Data' => $1,
		    }
		    );

		$self->element(
		    {
			'Name' => 'Parameters_gap-ext',
			'Data' => $2,
		    }
		    );
	    }
	    elsif (/\s*open\/ext:\s*(\-?\d+)\/(\-?\d+)/) {
		$self->element(
		    {
			'Name' => 'Parameters_gap-open',
			'Data' => $1,
		    }
		    );

		$self->element(
		    {
			'Name' => 'Parameters_gap-ext',
			'Data' => $2,
		    }
		    );
	    }

            $self->element(
                {
                    'Name' => 'FastaOutput_program',
                    'Data' => $self->{'_reporttype'}
                }
            );
        }
        elsif (/The best( related| unrelated)? scores are:/) {
            my $rel    = $1;
            my @labels = split;
            @labels = map {
                if ( $_ =~ m/^E\((\d+)\)$/o )
                {
                    $self->element(
                        { 'Name' => 'Statistics_eff-space', 'Data' => $1 } );
                    "evalue";
                }
                else {
                    # canonicalize changed column headers
                    if ($_ eq "gapl") {
                        $_ = "lgaps";
                    } elsif ($_ eq "gapq") {
                        $_ = "qgaps";
                    } elsif ($_ eq "E2()") {
                        $_ = "evalue2";
                    }
                    
                    $_;
                }
            } @labels[ $rel ? 5 : 4 .. $#labels ];

            while ( defined( $_ = $self->_readline() )
                && !/^\s+$/ )
            {
                my @line = split;

                if ( $line[-1] =~ m/\=/o && $labels[-1] ne 'aln_code' ) {
                    # unlabelled alignment hit;
                    push @labels, "aln_code";
                }

                if ($line[0] eq "+-") {
		    $m9HSP = 1;
                    # parse HSP, add to last parsed Hit
                    my %hspData;
                    
                    my @line = split;

                    @hspData{@labels} = splice( @line, @line - @labels );
                    $hspData{lframe} = $hit_signifs[-1]->{lframe};
                    
                    push @{$hit_signifs[-1]->{HSPs}}, \%hspData;

                    next;
                }
		elsif ($line[0] eq '>>><<<') {
		    last;
		}

                my (%data, %hspData);
                @data{@labels} = @hspData{@labels} = splice( @line, @line - @labels );
                if ( $line[-1] =~ m/\[([1-6rf])\]/o ) {
                    my $fr = $1;
                    $hspData{lframe} = $data{lframe} = (
                        $fr =~ /\d/o
                        ? ( $fr <= 3 ? "+$fr" : "-@{[$fr-3]}" )
                        : ( $fr eq 'f' ? '+1' : '-1' )
                    );
                    pop @line;
                }
                else {
                    $hspData{lframe} = $data{lframe} = '0';
                }

                if ( $line[-1] =~ m/^\(?(\d+)\)$/ ) {
                    $data{hit_len} = $1;
                    pop @line;
                    if ( $line[-1] =~ m/^\($/ ) {
                        pop @line;
                    }
                }
                else {
                    $data{hit_len} = 0;
                }

                # rebuild the first part of the line, preserving spaces:
                ($_) = m/^(\S+(?:\s+\S+){$#line})/;

                my ( $id, $desc ) = split( /\s+/, $_, 2 );
                my @pieces = split( /\|/, $id );
                my $acc = pop @pieces;
                $acc =~ s/\.\d+$//;

                @data{qw(id desc acc)} = ( $id, $desc, $acc );
                $data{HSPs} = [ \%hspData ];
                
                push @hit_signifs, \%data;
            }
        }
        elsif (
/^\s*([T]?FAST[XYAF]).+,\s*(\S+)\s*matrix[^\]]+?(xS)?\]\s*ktup:\s*(\d+)/
          )
        {

            $self->element(
                {
                    'Name' => 'Parameters_matrix',
                    'Data' => $2
                }
            );
            $self->element(
                {
                    'Name' => 'Parameters_filter',
                    'Data' => defined $3 ? 1 : 0,
                }
            );
            $self->element(
                {
                    'Name' => 'Parameters_ktup',
                    'Data' => $4
                }
            );
            $self->{'_reporttype'} = $1
              if ( $self->{'_reporttype'} !~ /FAST[PN]/i );

#
# get gap-pen line for FASTA33, which is not on the matrix line
#
# FASTA (3.36 June 2000) function [optimized, BL50 matrix (15:-5)] ktup: 2
#  join: 36, opt: 24, gap-pen: -12/ -2, width:  16
#
	    $_ = $self->_readline();
	    if (/(?:gap\-pen|open\/ext):\s+([\-\+]?\d+)\s*\/\s*([\-\+]?\d+)/) {
		$self->element(
		    {
			'Name' => 'Parameters_gap-open',
			'Data' => $1
		    }
		    );
		$self->element(
		    {
			'Name' => 'Parameters_gap-ext',
			'Data' => $2
		    }
		    );
	    }

            $self->element(
                {
                    'Name' => 'FastaOutput_program',
                    'Data' => $self->{'_reporttype'}
                }
            );

        }
        elsif (/^Algorithm:\s+(\S+)\s+.*\s*\(([^)]+)\)\s+(\S+)/) {
            $self->{'_reporttype'} = $1
              if ( $self->{'_reporttype'} !~ /FAST[PN]/i );
        }
        elsif ( /^Parameters:/ ) {    # FASTA 35.04/FASTA 36
	    m/Parameters:\s+(\S+)\s+matrix\s+\([^\)]+\)(xS)?,?\s/;
            $self->element(
                {
                    'Name' => 'Parameters_matrix',
                    'Data' => $1
                }
            );
            $self->element(
                {
                    'Name' => 'Parameters_filter',
                    'Data' => defined $2 ? $2 : 0,
                }
            );
	    if (/ktup:\s(\d+)/) {
		$self->element(
		    {
			'Name' => 'Parameters_ktup',
			'Data' => $1
		    }
		    );
		if (/ktup:\s\d+$/) {
		    $_ = $self->_readline();
		}

	    }
	    if (/(?:gap\-pen|open\/ext):\s+([\-\+]?\d+)\s*\/\s*([\-\+]?\d+)/) {
		$self->element(
		    {
			'Name' => 'Parameters_gap-open',
			'Data' => $1
		    }
		    );
		$self->element(
		    {
			'Name' => 'Parameters_gap-ext',
			'Data' => $2
		    }
		    );
	    }
            $self->element(
                {
                    'Name' => 'FastaOutput_program',
                    'Data' => $self->{'_reporttype'}
                }
            );
        }
        elsif (
	    /^\s+ktup:\s*(\d+),/
          )
        {
            $self->element(
                {
                    'Name' => 'Parameters_ktup',
                    'Data' => $1
                }
            );
        }
        elsif (/^(>--)$/ || /^>>(?!>)(.+?)\s+(?:\((\d+)\s*(aa|nt)\))?$/) {
            if ( $self->in_element('hsp') ) {
                $self->end_element( { 'Name' => 'Hsp' } );
            }
            
            my $firstHSP = 0;
            if ($1 ne ">--") {
                $firstHSP = 1;

                my ($hit_id, $len, $alphabet) = ($1, $2, $3);
                if (!$len || !$alphabet) {
                    WRAPPED:
                    while (defined($_ = $self->_readline)) {
                        if (/(.*?)\s+\((\d+)\s*(aa|nt)\)/) {
                            ($len, $alphabet) = ($2, $3);
                            $hit_id .= $1 ? " ".$1 : '';
                            last WRAPPED;
                        }
                        if (/^>>(?!>)/) { # too far, throw
                            $self->throw("Couldn't find length, bailing");
                        }
                    }
                }
                if ( $self->in_element('hit') ) {
                    $self->end_element( { 'Name' => 'Hit' } );
                    shift @hit_signifs if @hit_signifs;
                }
    
                $self->start_element( { 'Name' => 'Hit' } );
                $self->element(
                    {
                        'Name' => 'Hit_len',
                        'Data' => $len
                    }
                );
                my ( $id, $desc ) = split( /\s+/, $hit_id, 2 );
                $self->element(
                    {
                        'Name' => 'Hit_id',
                        'Data' => $id
                    }
                );
    
                #$self->debug("Hit ID is $id\n");
                my @pieces = split( /\|/, $id );
                my $acc = pop @pieces;
                $acc =~ s/\.\d+$//;
                $self->element(
                    {
                        'Name' => 'Hit_accession',
                        'Data' => $acc
                    }
                );
                $self->element(
                    {
                        'Name' => 'Hit_def',
                        'Data' => $desc
                    }
                );
            }
	    else {
#		push @{$hit_signifs[0]->{HSPs}}, $current_hsp;
	    }


            $_ = $self->_readline();
            my ( $score, $bits, $e, $e2 ) = /Z-score: \s* (\S+) \s*
                               (?: bits: \s* (\S+) \s+ )?
                               (?: E|expect ) \s* \((?:\d+)?\) :? \s*(\S+)
                               (?: \s* E2 \s* \(\) :? \s*(\S+) )?
                               /ox;
            $bits = $score unless defined $bits;

            my ($v);

	    if ($firstHSP && !$m9HSP) {
	      $v = shift @{$hit_signifs[0]->{HSPs}}
		if (@hit_signifs && @{$hit_signifs[0]->{HSPs}});
	      $current_hsp = $v;
	    }
	    else {
	      $v = $current_hsp;
	    }

            if ( defined $v ) {
                @{$v}{qw(evalue evalue2 bits z-sc)} = ( $e, $e2, $bits, $score );
            }

            if ($firstHSP) {
                $self->element(
                    {
                        'Name' => 'Hit_signif',
                        'Data' => $v ? $v->{evalue} : $e
                    }
                );
                $self->element(
                    {
                        'Name' => 'Hit_score',
                        'Data' => $v ? $v->{bits} : $bits
                    }
                );
            }
            
            $self->start_element( { 'Name' => 'Hsp' } );

            $self->element(
                {
                    'Name' => 'Hsp_score',
                    'Data' => $v ? $v->{'z-sc'} : $score
                }
            );
            $self->element(
                {
                    'Name' => 'Hsp_evalue',
                    'Data' => $v ? $v->{evalue} : $e
                }
            );
            $self->element(
                {
                    'Name' => 'Hsp_evalue2',
                    'Data' => $v && exists($v->{evalue2}) ? $v->{evalue2} : $e2
                }
            ) if (($v && exists($v->{evalue2})) || defined $e2);

            $self->element(
                {
                    'Name' => 'Hsp_bit-score',
                    'Data' => $v ? $v->{bits} : $bits
                }
            );
            $_ = $self->_readline();

            if (s/global\/.* score:\s*(\d+)\;?//) {
                $self->element(
                    {
                        'Name' => 'Hsp_sw-score',
                        'Data' => $1
                    }
                );
            }
            if (s/Smith-Waterman score:\s*(\d+)\;?//) {
                $self->element(
                    {
                        'Name' => 'Hsp_sw-score',
                        'Data' => $1
                    }
                );
            }
            if (
                / (\d*\.?\d+)\% \s* identity
                 (?:\s* \(\s*(\S+)\% \s* (?:ungapped|similar) \) )?
                 \s* in \s* (\d+) \s+ (?:aa|nt) \s+ overlap \s*
                 \( (\d+) \- (\d+) : (\d+) \- (\d+) \)
               /x
              )
            {
                my ( $identper, $gapper, $len, $querystart, $queryend,
                    $hitstart, $hitend )
                  = ( $1, $2, $3, $4, $5, $6, $7 );
                my $ident = sprintf( "%.0f", ( $identper / 100 ) * $len );
                my $positive = sprintf( "%.0f", ( $gapper / 100 ) * $len );

                $self->element(
                    {
                        'Name' => 'Hsp_identity',
                        'Data' => $ident
                    }
                );
                $self->element(
                    {
                        'Name' => 'Hsp_positive',
                        'Data' => $positive
                    }
                );
                $self->element(
                    {
                        'Name' => 'Hsp_align-len',
                        'Data' => $len
                    }
                );

                $self->element(
                    {
                        'Name' => 'Hsp_query-from',
                        'Data' => $querystart
                    }
                );
                $self->element(
                    {
                        'Name' => 'Hsp_query-to',
                        'Data' => $queryend
                    }
                );
                $self->element(
                    {
                        'Name' => 'Hsp_hit-from',
                        'Data' => $hitstart
                    }
                );
                $self->element(
                    {
                        'Name' => 'Hsp_hit-to',
                        'Data' => $hitend
                    }
                );

            }

            if ($v) {
                $self->element(
                    { 'Name' => 'Hsp_querygaps', 'Data' => $v->{qgaps} } )
                  if exists $v->{qgaps};
                $self->element(
                    { 'Name' => 'Hsp_hitgaps', 'Data' => $v->{lgaps} } )
                  if exists $v->{lgaps};

                if ( $self->{'_reporttype'} =~ m/^FAST[NXY]$/o ) {
                    if ( 8 == scalar grep { exists $v->{$_} }
                        qw(an0 ax0 pn0 px0 an1 ax1 pn1 px1) )
                    {
                        if ( $v->{ax0} < $v->{an0} ) {
                            $self->element(
                                {
                                    'Name' => 'Hsp_query-frame',
                                    'Data' =>
                                      "-@{[(($v->{px0} - $v->{ax0}) % 3) + 1]}"
                                }
                            );
                        }
                        else {
                            $self->element(
                                {
                                    'Name' => 'Hsp_query-frame',
                                    'Data' =>
                                      "+@{[(($v->{an0} - $v->{pn0}) % 3) + 1]}"
                                }
                            );
                        }
                        if ( $v->{ax1} < $v->{an1} ) {
                            $self->element(
                                {
                                    'Name' => 'Hsp_hit-frame',
                                    'Data' =>
                                      "-@{[(($v->{px1} - $v->{ax1}) % 3) + 1]}"
                                }
                            );
                        }
                        else {
                            $self->element(
                                {
                                    'Name' => 'Hsp_hit-frame',
                                    'Data' =>
                                      "+@{[(($v->{an1} - $v->{pn1}) % 3) + 1]}"
                                }
                            );
                        }
                    }
                    else {
                        $self->element(
                            {
                                'Name' => 'Hsp_query-frame',
                                'Data' => $v->{lframe}
                            }
                        );
                        $self->element(
                            { 'Name' => 'Hsp_hit-frame', 'Data' => 0 } );
                    }
                }
                else {
                    $self->element(
                        { 'Name' => 'Hsp_query-frame', 'Data' => 0 } );
                    $self->element(
                        { 'Name' => 'Hsp_hit-frame', 'Data' => $v->{lframe} } );
                }

            }
            else {
                $self->warn("unable to parse FASTA score line: $_");
            }
        }
        elsif (/\d+\s*residues\s*in\s*\d+\s*query\s*sequences/) {
            if ( $self->in_element('hsp') ) {
                $self->end_element( { 'Name' => 'Hsp' } );
            }
            if ( $self->in_element('hit') ) {
                $self->end_element( { 'Name' => 'Hit' } );
                shift @hit_signifs if @hit_signifs;

            }

           #       $_ = $self->_readline();
           #       my ( $liblen,$libsize) = /(\d+)\s+residues\s*in(\d+)\s*library/;
           # fast forward to the end of the file as there is
           # nothing else left to do with this file and want to be sure and
           # reset it
            while ( defined( $_ = $self->_readline() ) ) {
                last if (/^Function used was/);
                if (
                    /(\S+)\s+searches\s+a\s+(protein\s+or\s+DNA\s+)?
           sequence/oxi || /(\S+)\s+compares\s+a/oi
                  )
                {
                    $self->_pushback($_);
                }
            }

            $self->_processHits(@hit_signifs) if @hit_signifs;
            
            $self->end_element( { 'Name' => 'FastaOutput' } );
            return $self->end_document();
        }
        elsif (/^\s*\d+\s*>>>/) {
            if ( $self->within_element('FastaOutput') ) {
                if ( $self->in_element('hsp') ) {
                    $self->end_element( { 'Name' => 'Hsp' } );
                }
                if ( $self->in_element('hit') ) {
                    $self->end_element( { 'Name' => 'Hit' } );
                    shift @hit_signifs if @hit_signifs;

                }

                $self->_processHits(@hit_signifs) if (@hit_signifs);
                
                $self->end_element( { 'Name' => 'FastaOutput' } );
                $self->_pushback($_);
                return $self->end_document();
            }
            else {
                $self->start_element( { 'Name' => 'FastaOutput' } );
                $self->{'_result_count'}++;
                $seentop = 1;
                $self->element(
                    {
                        'Name' => 'FastaOutput_program',
                        'Data' => $self->{'_reporttype'}
                    }
                );
                $self->element(
                    {
                        'Name' => 'FastaOutput_version',
                        'Data' => $self->{'_version'}
                    }
                );

                my ( $type, $querylen, $querytype, $querydef );

                if (/^\s*\d+\s*>>>(.*)/) {
                    $querydef = $1;
                    if ( $querydef =~ /^(.*?)\s+(?:\-\s+)?(\d+)\s+(aa|nt).*$/o )
                    {
                        ( $querydef, $querylen, $querytype ) = ( $1, $2, $3 );
                    }
                }

                if (   $self->{'_reporttype'}
                    && $self->{'_reporttype'} eq 'FASTA' )
                {
                    if ( $querytype eq 'nt' ) {
                        $self->{'_reporttype'} = 'FASTN';
                    }
                    elsif ( $querytype eq 'aa' ) {
                        $self->{'_reporttype'} = 'FASTP';
                    }
                }
                my ( $name, $descr ) =
                  ( $querydef =~ m/^(\S+)(?:\s+(.*))?\s*$/o );
                $self->element(
                    {
                        'Name' => 'FastaOutput_query-def',
                        'Data' => $name
                    }
                );
                $self->element(
                    {
                        'Name' => 'FastaOutput_querydesc',
                        'Data' => $descr
                    }
                );
                if ($querylen) {
                    $self->element(
                        {
                            'Name' => 'FastaOutput_query-len',
                            'Data' => $querylen
                        }
                    );
                }
                else {
                    $self->warn("unable to find and set query length");
                }
                if ( defined( $_ = $self->_readline() )
                    && ( /^\s*vs\s+(\S+)/ || /^Library:\s+(\S+)/ ) )
                {
                    $self->element(
                        {
                            'Name' => 'FastaOutput_db',
                            'Data' => $1
                        }
                    );
                }

            }
        }
        elsif ( $self->in_element('hsp') ) {
            my @data  = ( [], [], [] );
            my $count = 0;
            my $len   = $self->idlength + 1;
            my ($seq1_id);
            while ( defined($_) ) {
                chomp;
                #$self->debug("$count $_\n");
                if (/residues in \d+\s+query\s+sequences/o) {
                    $self->_pushback($_);
                    last;
                }
                elsif (/^>>>\*\*\*/o) {
                    $self->end_element( { Name => "Hsp" } );
                    last;
                }
                elsif (/^>>/o) {
                    $self->_pushback($_);
                    last;
                }
                elsif (/^\s*\d+\s*>>>/o) {
                    $self->_pushback($_);
                    last;
                }
                if ( $count == 0 ) {
                    if (/^(\S+)\s+/) {
                        $self->_pushback($_);
                        $count = 2;
                    }
                    elsif ( /^\s+\d+/ || /^\s+$/ ) {

                        # do nothing, this is really a 0 line
                    }
                    elsif ( length($_) == 0 ) {
                        $count = -1;
                    }
                    else {
                        $self->_pushback($_);
                        $count = 0;
                    }
                }
                elsif ( $count == 1 || $count == 3 ) {
                    if (/^(\S+)\s+/) {
                        $len = CORE::length($1) if $len < CORE::length($1);
                        s/\s+$//;    # trim trailing spaces,we don't want them
                        push @{ $data[ $count - 1 ] }, substr( $_, $len );
                    }
                    elsif (/^\s+(\d+)/) {
                        $count = -1;
                        $self->_pushback($_);
                    }
                    elsif ( /^\s+$/ || length($_) == 0 ) {
                        $count = 5;

                        # going to skip these
                    }
                    elsif ( /\s+\S+fasta3\d\s+/) {
                        # this is something that looks like a path but contains
                        # the fasta3x executable string, such as:
                        # /usr/local/fasta3/bin/fasta35 -n -U -Q -H -A -E 2.0 -C 19 -m 0 -m 9i test.fa ../other_mirs.fa -O test.fasta35
                        last;
                    }
                    else {
                        $self->throw(
                            "Unrecognized alignment line ($count) '$_'");
                    }
                }
                elsif ( $count == 2 ) {
                    if (/^\s+\d+\s+/) {
                        $self->warn("$_\n") if $self->verbose > 0;

                        # we are on a Subject part of the alignment
                        # but we THOUGHT we were on the Query
                        # move that last line to the proper place
                        push @{ $data[2] }, pop @{ $data[0] };
                        $count = 4;
                    }
                    else {

                        # toss the first IDLENGTH characters of the line
                        if ( length($_) >= $len ) {
                            push @{ $data[ $count - 1 ] }, substr( $_, $len );
                        }
                    }
                }
                last if ( $count++ >= 5 );
                $_ = $self->_readline();
            }
            if ( @{ $data[0] } || @{ $data[2] } ) {
                $self->characters(
                    {
                        'Name' => 'Hsp_qseq',
                        'Data' => join( '', @{ $data[0] } )
                    }
                );
                $self->characters(
                    {
                        'Name' => 'Hsp_midline',
                        'Data' => join( '', @{ $data[1] } )
                    }
                );
                $self->characters(
                    {
                        'Name' => 'Hsp_hseq',
                        'Data' => join( '', @{ $data[2] } )
                    }
                );
            }
        }
        else {
            if ( !$seentop ) {
                $self->debug($_);
                #$self->warn("unrecognized FASTA Family report file!");
                #return;
            }
        }
    }
    if ( $self->in_element('result') ) {
        if ( $self->in_element('hsp') ) {
            $self->end_element( { 'Name' => 'Hsp' } );
        }
        if ( $self->in_element('hit') ) {
            $self->end_element( { 'Name' => 'Hit' } );
            shift @hit_signifs if @hit_signifs;
        }
        $self->end_element( { 'Name' => 'FastaOutput' } );
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
    my $nm = $data->{'Name'};
    if ( my $type = $MODEMAP{$nm} ) {
        $self->_mode($type);
        if ( my $handler = $self->_will_handle($type) ) {
            my $func = sprintf( "start_%s", lc $type );
            $handler->$func( $data->{'Attributes'} );
        }
        unshift @{ $self->{'_elements'} }, $type;
    }
    if ( $nm eq 'FastaOutput' ) {
        $self->{'_values'} = {};
        $self->{'_result'} = undef;
        $self->{'_mode'}   = '';
    }

}

=head2 end_element

 Title   : start_element
 Usage   : $eventgenerator->end_element
 Function: Handles an end element event
 Returns : none
 Args    : hashref with at least 2 keys 'Data' and 'Name'


=cut

sub end_element {
    my ( $self, $data ) = @_;
    my $nm = $data->{'Name'};
    my $rc;

    # Hsp are sort of weird, in that they end when another
    # object begins so have to detect this in end_element for now
    if ( $nm eq 'Hsp' ) {
        foreach (qw(Hsp_qseq Hsp_midline Hsp_hseq)) {
            $self->element(
                {
                    'Name' => $_,
                    'Data' => $self->{'_last_hspdata'}->{$_}
                }
            );
        }
        $self->{'_last_hspdata'} = {};
    }

    if ( my $type = $MODEMAP{$nm} ) {
        if ( my $handler = $self->_will_handle($type) ) {
            my $func = sprintf( "end_%s", lc $type );
            $rc = $handler->$func( $self->{'_reporttype'}, $self->{'_values'} );
        }
        shift @{ $self->{'_elements'} };

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
        $self->warn("unknown nm $nm, ignoring\n");
    }
    $self->{'_last_data'} = '';    # remove read data if we are at
                                   # end of an element
    $self->{'_result'} = $rc if ( $nm eq 'FastaOutput' );
    return $rc;

}

=head2 element

 Title   : element
 Usage   : $eventhandler->element({'Name' => $name, 'Data' => $str});
 Function: Convience method that calls start_element, characters, end_element
 Returns : none
 Args    : Hash ref with the keys 'Name' and 'Data'


=cut

sub element {
    my ( $self, $data ) = @_;
    $self->start_element($data);
    $self->characters($data);
    $self->end_element($data);
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

    return unless ( defined $data->{'Data'} );
    if ( $data->{'Data'} =~ /^\s+$/ ) {
        return unless $data->{'Name'} =~ /Hsp\_(midline|qseq|hseq)/;
    }

    if (   $self->in_element('hsp')
        && $data->{'Name'} =~ /Hsp\_(qseq|hseq|midline)/ )
    {

        $self->{'_last_hspdata'}->{ $data->{'Name'} } .= $data->{'Data'};
    }

    $self->{'_last_data'} = $data->{'Data'};
}

=head2 _mode

 Title   : _mode
 Usage   : $obj->_mode($newval)
 Function: 
 Example : 
 Returns : value of _mode
 Args    : newvalue (optional)


=cut

sub _mode {
    my ( $self, $value ) = @_;
    if ( defined $value ) {
        $self->{'_mode'} = $value;
    }
    return $self->{'_mode'};
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
      if (!defined $name && !defined $self->{'_elements'}
        || scalar @{ $self->{'_elements'} } == 0 );
    foreach ( @{ $self->{'_elements'} } ) {
        if ( $_ eq $name || $_ eq $MODEMAP{$name} ) {
            return 1;
        }
    }
    return 0;
}

=head2 in_element

 Title   : in_element
 Usage   : if( $eventgenerator->in_element($element) ) {}
 Function: Test if we are in a particular element
           This is different than 'in' because within can be tested
           for a whole block.
 Returns : boolean
 Args    : string element name 


=cut

sub in_element {
    my ( $self, $name ) = @_;
    return 0 if !defined $self->{'_elements'}->[0];
    return (
        $self->{'_elements'}->[0] eq $name
          || ( exists $MODEMAP{$name}
            && $self->{'_elements'}->[0] eq $MODEMAP{$name} )
    );
}

=head2 start_document

 Title   : start_document
 Usage   : $eventgenerator->start_document
 Function: Handles a start document event
 Returns : none
 Args    : none


=cut

sub start_document {
    my ($self) = @_;
    $self->{'_lasttype'} = '';
    $self->{'_values'}   = {};
    $self->{'_result'}   = undef;
    $self->{'_mode'}     = '';
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
    my ( $self, @args ) = @_;
    return $self->{'_result'};
}

=head2 idlength

 Title   : idlength
 Usage   : $obj->idlength($newval)
 Function: Internal storage of the length of the ID desc
           in the HSP alignment blocks.  Defaults to
           $IDLENGTH class variable value
 Returns : value of idlength
 Args    : newvalue (optional)


=cut

sub idlength {
    my ( $self, $value ) = @_;
    if ( defined $value ) {
        $self->{'_idlength'} = $value;
    }
    return $self->{'_idlength'} || $IDLENGTH;
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

sub attach_EventHandler {
    my ( $self, $handler ) = @_;

    $self->SUPER::attach_EventHandler($handler);

    # Optimization: caching the EventHandler since it is used a lot
    # during the parse.

    $self->{'_handler_cache'} = $handler;
    return;
}

=head2 _will_handle

 Title   : _will_handle
 Usage   : Private method. For internal use only.
              if( $self->_will_handle($type) ) { ... }
 Function: Provides an optimized way to check whether or not an element of a 
           given type is to be handled.
 Returns : Reference to EventHandler object if the element type is to be handled.
           undef if the element type is not to be handled.
 Args    : string containing type of element.

Optimizations:

=over 2

=item 1

Using the cached pointer to the EventHandler to minimize repeated
lookups.

=item 2

Caching the will_handle status for each type that is encountered so
that it only need be checked by calling
handler-E<gt>will_handle($type) once.

=back

This does not lead to a major savings by itself (only 5-10%).  In
combination with other optimizations, or for large parse jobs, the
savings good be significant.

To test against the unoptimized version, remove the parentheses from
around the third term in the ternary " ? : " operator and add two
calls to $self-E<gt>_eventHandler().

=cut

sub _will_handle {
    my ( $self, $type ) = @_;
    my $handler = $self->{'_handler_cache'};
    my $will_handle =
      defined( $self->{'_will_handle_cache'}->{$type} )
      ? $self->{'_will_handle_cache'}->{$type}
      : ( $self->{'_will_handle_cache'}->{$type} =
          $handler->will_handle($type) );

    return $will_handle ? $handler : undef;
}

=head2 _processHits

 Title   : _processHits
 Usage   : Private method. For internal use only.
 Function: Process/report any hits/hsps we saw in the top table, not in alignments.
 Returns : nothing.
 Args    : array of hits to process.

=cut

sub _processHits {

  my ($self, @hit_signifs) = @_;

  # process remaining best hits
  for my $hit (@hit_signifs) {

    # Hsp_score Hsp_evalue Hsp_bit-score
    # Hsp_sw-score Hsp_gaps Hsp_identity Hsp_positive
    # Hsp_align-len Hsp_query-from Hsp_query-to
    # Hsp_hit-from Hsp_hit-to Hsp_qseq Hsp_midline

    $self->start_element( { 'Name' => 'Hit' } );
    $self->element(
		   {
		    'Name' => 'Hit_len',
		    'Data' => $hit->{hit_len}
		   }
		  ) if exists $hit->{hit_len};
    $self->element(
		   {
		    'Name' => 'Hit_id',
		    'Data' => $hit->{id}
		   }
		  ) if exists $hit->{id};
    $self->element(
		   {
		    'Name' => 'Hit_accession',
		    'Data' => $hit->{acc}
		   }
		  ) if exists $hit->{acc};
    $self->element(
		   {
		    'Name' => 'Hit_def',
		    'Data' => $hit->{desc}
		   }
		  ) if exists $hit->{desc};
    $self->element(
		   {
		    'Name' => 'Hit_signif',
		    'Data' => $hit->{evalue}
		   }
		  ) if exists $hit->{evalue};
    $self->element(
		   {
		    'Name' => 'Hit_score',
		    'Data' => $hit->{bits}
		   }
		  ) if exists $hit->{bits};

    for my $hsp (@{$hit->{HSPs}}) {

      $self->start_element( { 'Name' => 'Hsp' } );
      $self->element({'Name' => 'Hsp_score', 'Data' => $hsp->{'z-sc'}})
	if exists $hsp->{'z-sc'};
      $self->element({'Name' => 'Hsp_evalue', 'Data' => $hsp->{evalue} } )
	if exists $hsp->{evalue};
      $self->element({'Name' => 'Hsp_evalue2', 'Data' => $hsp->{evalue2} } )
	if exists $hsp->{evalue2};

      $self->element({'Name' => 'Hsp_bit-score', 'Data' => $hsp->{bits} } )
	if exists $hsp->{bits};
      $self->element({'Name' => 'Hsp_sw-score', 'Data' => $hsp->{'n-w'} } )
	if exists $hsp->{'n-w'};
      $self->element({'Name' => 'Hsp_sw-score', 'Data' => $hsp->{sw} } )
	if exists $hsp->{sw};
      $self->element({'Name' => 'Hsp_gaps', 'Data' => $hsp->{'%_gid'} } )
	if exists $hsp->{'%_gid'};
      $self->element({
		      'Name' => 'Hsp_identity',
		      'Data' =>
		      sprintf( "%.0f", $hsp->{'%_id'} * $hsp->{alen} )
		     }) if ( exists $hsp->{'%_id'} && exists $hsp->{alen} );

      if ( exists $hsp->{'%_gid'} ) {
	$self->element(
		       {
			'Name' => 'Hsp_positive',
			'Data' =>
			sprintf( "%.0f", $hsp->{'%_gid'} * $hsp->{alen} )
		       }
		      ) if exists $hsp->{'%_gid'} && exists $hsp->{alen};
      } else {
	$self->element(
		       {
			'Name' => 'Hsp_positive',
			'Data' =>
			sprintf( "%.0f", $hsp->{'%_id'} * $hsp->{alen} )
		       }
		      ) if ( exists $hsp->{'%_id'} && exists $hsp->{alen} );
      }

      $self->element(
		     {
		      'Name' => 'Hsp_align-len', 'Data' => $hsp->{alen} } )
	if exists $hsp->{alen};
      $self->element(
		     {
		      'Name' => 'Hsp_query-from', 'Data' => $hsp->{an0} } )
	if exists $hsp->{an0};
      $self->element(
		     {
		      'Name' => 'Hsp_query-to', 'Data' => $hsp->{ax0} } )
	if exists $hsp->{ax0};
      $self->element(
		     {
		      'Name' => 'Hsp_hit-from', 'Data' => $hsp->{an1} } )
	if exists $hsp->{an1};
      $self->element(
		     {
		      'Name' => 'Hsp_hit-to', 'Data' => $hsp->{ax1} } )
	if exists $hsp->{ax1};

      $self->element(
		     {
		      'Name' => 'Hsp_querygaps', 'Data' => $hsp->{qgaps} } )
	if exists $hsp->{qgaps};
      $self->element(
		     {
		      'Name' => 'Hsp_hitgaps', 'Data' => $hsp->{lgaps} } )
	if exists $hsp->{lgaps};

      if ( $self->{'_reporttype'} =~ m/^FAST[NXY]$/o ) {
	if ( 8 == scalar grep { exists $hsp->{$_} }
	     qw(an0 ax0 pn0 px0 an1 ax1 pn1 px1) ) {
	  if ( $hsp->{ax0} < $hsp->{an0} ) {
	    $self->element(
			   {
			    'Name' => 'Hsp_query-frame',
			    'Data' =>
			    "-@{[(($hsp->{px0} - $hsp->{ax0}) % 3) + 1]}"
			   }
			  );
	  } else {
	    $self->element(
			   {
			    'Name' => 'Hsp_query-frame',
			    'Data' =>
			    "+@{[(($hsp->{an0} - $hsp->{pn0}) % 3) + 1]}"
			   }
			  );
	  }
	  if ( $hsp->{ax1} < $hsp->{an1} ) {
	    $self->element(
			   {
			    'Name' => 'Hsp_hit-frame',
			    'Data' =>
			    "-@{[(($hsp->{px1} - $hsp->{ax1}) % 3) + 1]}"
			   }
			  );
	  } else {
	    $self->element(
			   {
			    'Name' => 'Hsp_hit-frame',
			    'Data' =>
			    "+@{[(($hsp->{an1} - $hsp->{pn1}) % 3) + 1]}"
			   }
		);
	  }
	} else {
	  $self->element(
			 {
			  'Name' => 'Hsp_query-frame',
			  'Data' => $hsp->{lframe}
			 }
			);
	  $self->element(
			 {
			  'Name' => 'Hsp_hit-frame', 'Data' => 0 } );
	}
      } else {
	$self->element(
		       {
			'Name' => 'Hsp_query-frame', 'Data' => 0 } );
	$self->element(
		       {
			'Name' => 'Hsp_hit-frame',
			'Data' => $hsp->{lframe}
		       }
		      );
      }

      $self->end_element( { 'Name' => 'Hsp' } );
    }

    $self->end_element( { 'Name' => 'Hit' } );

  }
}



1;

