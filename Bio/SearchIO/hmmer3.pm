#
# BioPerl module for Bio::SearchIO::hmmer3
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Thomas Sharpton <thomas.sharpton@gmail.com>
#
# Copyright Thomas Sharpton
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::hmmer3

=head1 SYNOPSIS

use Bio::SearchIO;

my $searchio = Bio::SearchIO->new(
    -format  => 'hmmer',
    -version => 3,
    -file    => 'hmmsearch.out'
);

my $result = $searchio->next_result;
my $hit = $result->next_hit;
print $hit->name, $hit->description, $hit->significance, 
      $hit->score, "\n";

my $hsp = $hit->next_hsp;
print $hsp->start('hit'), $hsp->end('hit'), $hsp->start('query'), 
      $hsp->end('query'), "\n";

=head1 DESCRIPTION

Code to parse output from hmmsearch, hmmscan, and nhmmer, compatible with
both version 2 and version 3 of the HMMER package from L<http://hmmer.org>.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support

Please direct usage questions or support issues to the mailing list:

L<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and
reponsive experts will be able look at the problem and quickly
address it. Please include a thorough description of the problem
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Thomas Sharpton

Email thomas.sharpton@gmail.com

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

briano at bioteam.net

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SearchIO::hmmer3;

use strict;
use Data::Dumper;
use Bio::Factory::ObjectFactory;
use Bio::Tools::IUPAC;
use vars qw(%MAPPING %MODEMAP);
use base qw(Bio::SearchIO::hmmer);

BEGIN {

    # mapping of HMMER items to Bioperl hash keys
    %MODEMAP = (
        'HMMER_Output' => 'result',
        'Hit'          => 'hit',
        'Hsp'          => 'hsp'
    );

    %MAPPING = (
        'Hsp_bit-score'    => 'HSP-bits',
        'Hsp_score'        => 'HSP-score',
        'Hsp_evalue'       => 'HSP-evalue',
        'Hsp_query-from'   => 'HSP-query_start',
        'Hsp_query-to'     => 'HSP-query_end',
        'Hsp_query-strand' => 'HSP-query_strand',
        'Hsp_hit-from'     => 'HSP-hit_start',
        'Hsp_hit-to'       => 'HSP-hit_end',
        'Hsp_hit-strand'   => 'HSP-hit_strand',
        'Hsp_positive'     => 'HSP-conserved',
        'Hsp_identity'     => 'HSP-identical',
        'Hsp_gaps'         => 'HSP-hsp_gaps',
        'Hsp_hitgaps'      => 'HSP-hit_gaps',
        'Hsp_querygaps'    => 'HSP-query_gaps',
        'Hsp_qseq'         => 'HSP-query_seq',
        'Hsp_csline'       => 'HSP-cs_seq',
        'Hsp_hseq'         => 'HSP-hit_seq',
        'Hsp_midline'      => 'HSP-homology_seq',
        'Hsp_pline'        => 'HSP-pp_seq',
        'Hsp_align-len'    => 'HSP-hsp_length',
        'Hsp_query-frame'  => 'HSP-query_frame',
        'Hsp_hit-frame'    => 'HSP-hit_frame',

        'Hit_id'        => 'HIT-name',
        'Hit_len'       => 'HIT-length',
        'Hit_accession' => 'HIT-accession',
        'Hit_desc'      => 'HIT-description',
        'Hit_signif'    => 'HIT-significance',
        'Hit_score'     => 'HIT-score',

        'HMMER_program'   => 'RESULT-algorithm_name',
        'HMMER_version'   => 'RESULT-algorithm_version',
        'HMMER_query-def' => 'RESULT-query_name',
        'HMMER_query-len' => 'RESULT-query_length',
        'HMMER_query-acc' => 'RESULT-query_accession',
        'HMMER_querydesc' => 'RESULT-query_description',
        'HMMER_hmm'       => 'RESULT-hmm_name',
        'HMMER_seqfile'   => 'RESULT-sequence_file',
        'HMMER_db'        => 'RESULT-database_name',
    );
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
    my ( $last, @hit_list, @hsp_list, %hspinfo, %hitinfo, %domaincounter );
    local $/ = "\n";
    local $_;

    my @ambiguous_nt = keys %Bio::Tools::IUPAC::IUB;
    my $ambiguous_nt = join '', @ambiguous_nt;

    my $verbose = $self->verbose;    # cache for speed? JES's idea in hmmer.pm
    $self->start_document();
    local ($_);

    # This is here to ensure that next_result doesn't produce infinite loop
    if ( !defined( $_ = $self->_readline ) ) {
        return undef;
    }
    else {
        $self->_pushback($_);
    }

    my $hit_counter = 0; # helper variable for non-unique hit IDs

    # Regex goes here for HMMER3
    # Start with hmmsearch processing
    while ( defined( $_ = $self->_readline ) ) {
        my $lineorig = $_;
        chomp;

        # Grab the program name
        if ( $_ =~ m/^\#\s(\S+)\s\:\:\s/ ) {
            my $prog = $1;

            # TO DO: customize the above regex to adapt to other
            # program types (hmmscan, etc)
            $self->start_element( { 'Name' => 'HMMER_Output' } );
            $self->{'_result_count'}++;   #Might need to move to another block
            $self->element(
                {   'Name' => 'HMMER_program',
                    'Data' => uc($prog)
                }
            );
        }

        # Get the HMMER package version and release date
        elsif ( $_ =~ m/^\#\sHMMER\s+(\S+)\s+\((.+)\)/ ) {
            my $version     = $1;
            my $versiondate = $2;
            $self->{'_hmmidline'} = $_;
            $self->element(
                {   'Name' => 'HMMER_version',
                    'Data' => $version
                }
            );
        }

        # Get the query info
        elsif ( $_ =~ /^\#\squery (?:\w+ )?file\:\s+(\S+)/ ) {
            if (   $self->{'_reporttype'} eq 'HMMSEARCH'
                || $self->{'_reporttype'} eq 'NHMMER' )
            {
                $self->{'_hmmfileline'} = $lineorig;
                $self->element(
                    {   'Name' => 'HMMER_hmm',
                        'Data' => $1
                    }
                );
            }
            elsif ( $self->{'_reporttype'} eq 'HMMSCAN' ) {
                $self->{'_hmmseqline'} = $lineorig;
                $self->element(
                    {   'Name' => 'HMMER_seqfile',
                        'Data' => $1
                    }
                );
            }
        }

        # If this is a report without alignments
        elsif ( $_ =~ m/^\#\sshow\salignments\sin\soutput/ ) {
            $self->{'_alnreport'} = 0;
        }

        # Get the database info
        elsif ( $_ =~ m/^\#\starget\s\S+\sdatabase\:\s+(\S+)/ ) {

            if (   $self->{'_reporttype'} eq 'HMMSEARCH'
                || $self->{'_reporttype'} eq 'NHMMER' )
            {
                $self->{'_hmmseqline'} = $lineorig;
                $self->element(
                    {   'Name' => 'HMMER_seqfile',
                        'Data' => $1
                    }
                );
            }
            elsif ( $self->{'_reporttype'} eq 'HMMSCAN' ) {
                $self->{'_hmmfileline'} = $lineorig;
                $self->element(
                    {   'Name' => 'HMMER_hmm',
                        'Data' => $1
                    }
                );
            }
        }

        # Get query data
        elsif ( $_ =~ s/^Query:\s+// ) {
            # For  multi-query reports
            if (    (   not exists $self->{_values}->{"RESULT-algorithm_name"}
                     or not exists $self->{_values}->{"RESULT-algorithm_version"}
                     )
                and exists $self->{_hmmidline}
                ) {
                my ($version, $versiondate) = $self->{_hmmidline} =~ m/^\#\sHMMER\s+(\S+)\s+\((.+)\)/;
                $self->element(
                    {   'Name' => 'HMMER_program',
                        'Data' => $self->{_reporttype}
                    }
                );
                $self->element(
                    {   'Name' => 'HMMER_version',
                        'Data' => $version
                    }
                );
            }
            if (    (   not exists $self->{_values}->{"RESULT-hmm_name"}
                     or not exists $self->{_values}->{"RESULT-sequence_file"}
                     )
                and (   exists $self->{_hmmfileline}
                     or exists $self->{_hmmseqline}
                     )
                ) {
                if (   $self->{'_reporttype'} eq 'HMMSEARCH'
                    or $self->{'_reporttype'} eq 'NHMMER'
                    ) {
                    my ($qry_file)    = $self->{_hmmfileline} =~ m/^\#\squery (?:\w+ )?file\:\s+(\S+)/;
                    my ($target_file) = $self->{_hmmseqline}  =~ m/^\#\starget\s\S+\sdatabase\:\s+(\S+)/;
                    $self->element(
                        {   'Name' => 'HMMER_hmm',
                            'Data' => $qry_file
                        }
                    );
                    $self->element(
                        {   'Name' => 'HMMER_seqfile',
                            'Data' => $target_file
                        }
                    );
                }
                elsif ( $self->{'_reporttype'} eq 'HMMSCAN' ) {
                    my ($qry_file)    = $self->{_hmmseqline}  =~ m/^\#\squery \w+ file\:\s+(\S+)/;
                    my ($target_file) = $self->{_hmmfileline} =~ m/^\#\starget\s\S+\sdatabase\:\s+(\S+)/;
                    $self->element(
                        {   'Name' => 'HMMER_seqfile',
                            'Data' => $qry_file
                        }
                    );
                    $self->element(
                        {   'Name' => 'HMMER_hmm',
                            'Data' => $target_file
                        }
                    );
                }
            }

            unless (s/\s+\[[L|M]\=(\d+)\]$//) {
                warn "Error parsing length for query, offending line $_\n";
                exit(0);
            }
            my $querylen = $1;
            $self->element(
                {   'Name' => 'HMMER_query-len',
                    'Data' => $querylen
                }
            );
            $self->element(
                {   'Name' => 'HMMER_query-def',
                    'Data' => $_
                }
            );
        }

        # Get Accession data
        elsif ( $_ =~ s/^Accession:\s+// ) {
            s/\s+$//;
            $self->element(
                {   'Name' => 'HMMER_query-acc',
                    'Data' => $_
                }
            );
        }

        # Get description data
        elsif ( $_ =~ s/^Description:\s+// ) {
            s/\s+$//;
            $self->element(
                {   'Name' => 'HMMER_querydesc',
                    'Data' => $_
                }
            );
        }

        # hmmsearch, nhmmer, and hmmscan-specific formatting here
        elsif (
            defined $self->{'_reporttype'}
            && (   $self->{'_reporttype'} eq 'HMMSEARCH'
                || $self->{'_reporttype'} eq 'HMMSCAN'
                || $self->{'_reporttype'} eq 'NHMMER' )
            )
        {
            # Complete sequence table data above inclusion threshold,
            # hmmsearch or hmmscan
            if ( $_ =~ m/Scores for complete sequence/ ) {
                while ( defined( $_ = $self->_readline ) ) {
                    if (   $_ =~ m/inclusion threshold/
                        || m/Domain( and alignment)? annotation for each/
                        || m/\[No hits detected/
                        || m!^//! )
                    {
                        $self->_pushback($_);
                        last;
                    }
                    elsif (   $_ =~ m/^\s+E-value\s+score/
                           || $_ =~ m/\-\-\-/
                           || $_ =~ m/^$/
                        )
                    {
                        next;
                    }

                    # Grab table data
                    $hit_counter++;
                    my ($eval_full,  $score_full, $bias_full, $eval_best,
                        $score_best, $bias_best,  $exp,       $n,
                        $hitid,      $desc,       @hitline
                    );
                    @hitline    = split( " ", $_ );
                    $eval_full  = shift @hitline;
                    $score_full = shift @hitline;
                    $bias_full  = shift @hitline;
                    $eval_best  = shift @hitline;
                    $score_best = shift @hitline;
                    $bias_best  = shift @hitline;
                    $exp        = shift @hitline;
                    $n          = shift @hitline;
                    $hitid      = shift @hitline;
                    $desc       = join " ", @hitline;

                    $desc = '' if ( !defined($desc) );

                    push @hit_list,
                        [ $hitid, $desc, $eval_full, $score_full ];
                    $hitinfo{"$hitid.$hit_counter"} = $#hit_list;
                }
            }

            # nhmmer
            elsif ( /Scores for complete hits/ ) {
                while ( defined( $_ = $self->_readline ) ) {
                    if (   /inclusion threshold/
                        || /Annotation for each hit/
                        || /\[No hits detected/
                        || m!^//! )
                    {
                        $self->_pushback($_);
                        last;
                    }
                    elsif (   $_ =~ m/^\s+E-value\s+score/
                           || $_ =~ m/\-\-\-/
                           || $_ =~ m/^$/
                        )
                    {
                        next;
                    }

                    # Grab table data
                    $hit_counter++;
                    my ($eval,  $score, $bias, $hitid,
                        $start, $end,   $desc, @hitline
                    );
                    @hitline = split( " ", $_ );
                    $eval    = shift @hitline;
                    $score   = shift @hitline;
                    $bias    = shift @hitline;
                    $hitid   = shift @hitline;
                    $start   = shift @hitline;
                    $end     = shift @hitline;
                    $desc    = join ' ', @hitline;

                    $desc = '' if ( !defined($desc) );

                    push @hit_list, [ $hitid, $desc, $eval, $score ];
                    $hitinfo{"$hitid.$hit_counter"} = $#hit_list;
                }
            }

            # Complete sequence table data below inclusion threshold
            elsif ( /inclusion threshold/ ) {
                while ( defined( $_ = $self->_readline ) ) {
                    if (   /Domain( and alignment)? annotation for each/
                        || /Internal pipeline statistics summary/ 
                        || /Annotation for each hit\s+\(and alignments\)/
                        )
                    {
                        $self->_pushback($_);
                        last;
                    }
                    elsif (   $_ =~ m/inclusion threshold/
                           || $_ =~ m/^$/
                        )
                    {
                        next;
                    }

                    # Grab table data
                    $hit_counter++;
                    my ($eval_full,  $score_full, $bias_full, $eval_best,
                        $score_best, $bias_best,  $exp,       $n,
                        $hitid,      $desc,       @hitline
                    );
                    @hitline    = split( " ", $_ );
                    $eval_full  = shift @hitline;
                    $score_full = shift @hitline;
                    $bias_full  = shift @hitline;
                    $eval_best  = shift @hitline;
                    $score_best = shift @hitline;
                    $bias_best  = shift @hitline;
                    $exp        = shift @hitline;
                    $n          = shift @hitline;
                    $hitid      = shift @hitline;
                    $desc       = join " ", @hitline;

                    $desc = '' if ( !defined($desc) );

                    push @hit_list,
                        [ $hitid, $desc, $eval_full, $score_full ];
                    $hitinfo{"$hitid.$hit_counter"} = $#hit_list;
                }
            }

            # Domain annotation for each sequence table data,
            # for hmmscan, hmmsearch & nhmmer
            elsif (   /Domain( and alignment)? annotation for each/
                   or /Annotation for each hit\s+\(and alignments\)/
                   ) {
                @hsp_list = ();    # Here for multi-query reports
                my $name;
                my $annot_counter = 0;

                while ( defined( $_ = $self->_readline ) ) {
                    if (   /\[No targets detected/
                        || /Internal pipeline statistics/ )
                    {
                        $self->_pushback($_);
                        last;
                    }

                    if ( $_ =~ m/^\>\>\s(\S*)\s+(.*)/ ) {
                        $name    = $1;
                        my $desc = $2;
                        $annot_counter++;
                        $domaincounter{"$name.$annot_counter"} = 0;

                        # The Hit Description from the Scores table can be truncated if
                        # its too long, so use the '>>' line description when its longer
                        if (length $hit_list[
                                             $hitinfo{"$name.$annot_counter"}
                                             ]
                                             [1] < length $desc
                            ) {
                            $hit_list[ $hitinfo{"$name.$annot_counter"} ][1] = $desc;
                        }

                        while ( defined( $_ = $self->_readline ) ) {
                            if (   $_ =~ m/Internal pipeline statistics/
                                || $_ =~ m/Alignments for each domain/
                                || $_ =~ m/^\s+Alignment:/
                                || $_ =~ m/^\>\>/ )
                            {
                                $self->_pushback($_);
                                last;
                            }
                            elsif (   $_ =~ m/^\s+score\s+bias/
                                   || $_ =~ m/^\s+\#\s+score/
                                   || $_ =~ m/^\s+------\s+/
                                   || $_ =~ m/^\s\-\-\-\s+/
                                   || $_ =~ m/^$/
                                )
                            {
                                next;
                            }

                            # Grab hsp data from table, push into @hsp;
                            if ($self->{'_reporttype'} =~ m/(?:HMMSCAN|HMMSEARCH|NHMMER)/) {
                                my ( $domain_num, $score,    $bias,
                                     $ceval,      $ieval,
                                     $hmm_start,  $hmm_stop, $hmm_cov,
                                     $seq_start,  $seq_stop, $seq_cov,
                                     $env_start,  $env_stop, $env_cov,
                                     $hitlength,  $acc );
                                my @vals;

                                if ( # HMMSCAN & HMMSEARCH
                                    ( $domain_num, $score,    $bias,
                                      $ceval,      $ieval,
                                      $hmm_start,  $hmm_stop, $hmm_cov,
                                      $seq_start,  $seq_stop, $seq_cov,
                                      $env_start,  $env_stop, $env_cov,
                                      $acc ) = (
                                            m|^\s+(\d+)\s\!*\?*\s+     # domain number
                                              (\S+)\s+(\S+)\s+         # score, bias
                                              (\S+)\s+(\S+)\s+         # c-eval, i-eval
                                              (\d+)\s+(\d+)\s+(\S+)\s+ # hmm start, stop, coverage
                                              (\d+)\s+(\d+)\s+(\S+)\s+ # seq start, stop, coverage
                                              (\d+)\s+(\d+)\s+(\S+)\s+ # env start, stop, coverage
                                              (\S+)                    # posterior probability accuracy
                                               \s*$|ox
                                        )
                                    ) {
                                    # Values assigned when IF succeeded

                                    # Try to get the Hit length from the alignment information
                                    $hitlength = 0;
                                    if ($self->{'_reporttype'} eq 'HMMSEARCH') {
                                        # For Hmmsearch, if seq coverage ends in ']' it means that the alignment
                                        # runs until the end. In that case add the END coordinate to @hitinfo
                                        # to use it as Hit Length
                                        if ( $seq_cov =~ m/\]$/ ) {
                                            $hitlength = $seq_stop;
                                        }
                                    }
                                    elsif ($self->{'_reporttype'} eq 'HMMSCAN') {
                                        # For Hmmscan, if hmm coverage ends in ']' it means that the alignment
                                        # runs until the end. In that case add the END coordinate to @hitinfo
                                        # to use it as Hit Length
                                        if ( $hmm_cov =~ m/\]$/ ) {
                                            $hitlength = $hmm_stop;
                                        }
                                    }
                                }
                                elsif ( # NHMMER
                                       ( $score,     $bias,     $ceval,
                                         $hmm_start, $hmm_stop, $hmm_cov,
                                         $seq_start, $seq_stop, $seq_cov,
                                         $env_start, $env_stop, $env_cov,
                                         $hitlength, $acc ) = (
                                            m|^\s+[!?]\s+
                                              (\S+)\s+(\S+)\s+(\S+)\s+ # score, bias, evalue
                                              (\d+)\s+(\d+)\s+(\S+)\s+ # hmm start, stop, coverage
                                              (\d+)\s+(\d+)\s+(\S+)\s+ # seq start, stop, coverage
                                              (\d+)\s+(\d+)\s+(\S+)\s+ # env start, stop, coverage
                                              (\d+)\s+(\S+)            # target length, pp accuracy
                                               .*$|ox
                                        )
                                    ) {
                                    # Values assigned when IF succeeded
                                }
                                else {
                                    print "Missed this line: $_\n";
                                    next;
                                }

                                my $info = $hit_list[ $hitinfo{"$name.$annot_counter"} ];
                                if ( !defined $info ) {
                                    $self->warn(
                                        "Incomplete information: can't find HSP $name in list of hits\n"
                                    );
                                    next;
                                }

                                $domaincounter{"$name.$annot_counter"}++;
                                my $hsp_key
                                    = $name . "_" . $domaincounter{"$name.$annot_counter"};

                                # Keep it simple for now. let's customize later
                                @vals = (
                                    $hmm_start, $hmm_stop,
                                    $seq_start, $seq_stop,
                                    $score,     $ceval,
                                    $hitlength, '',
                                    '',         '',
                                    '',         ''
                                );
                                push @hsp_list, [ $name, @vals ];
                                $hspinfo{"$hsp_key.$annot_counter"} = $#hsp_list;
                            }
                        }
                    }
                    elsif ( /Alignment(?:s for each domain)?:/ ) {
                        #line counter
                        my $count = 0;

                        # There's an optional block, so we sometimes need to
                        # count to 3, and sometimes to 4.
                        my $max_count = 3;
                        my $lastdomain;
                        my $hsp;
                        my ( $csline, $hline, $midline, $qline, $pline );

                        # To avoid deleting whitespaces from the homology line,
                        # keep track of the position and length of the alignment
                        # in each individual hline/qline, to take them as reference
                        # and use them in the homology line
                        my $align_offset = 0;
                        my $align_length = 0;

                        while ( defined( $_ = $self->_readline ) ) {
                            if (   $_ =~ m/^\>\>/
                                || $_ =~ m/Internal pipeline statistics/ )
                            {
                                $self->_pushback($_);
                                last;
                            }
                            elsif ($_ =~ m/^$/ )
                            {
                                next;
                            }

                            if (   $_ =~ /\s\s\=\=\sdomain\s(\d+)\s+/
                                or $_ =~ /\s\sscore:\s\S+\s+/
                                ) {
                                my $domainnum = $1 || 1;
                                $count = 0;
                                my $key = $name . "_" . $domainnum;
                                $hsp        = $hsp_list[ $hspinfo{"$key.$annot_counter"} ];
                                $csline     = $$hsp[-5];
                                $hline      = $$hsp[-4];
                                $midline    = $$hsp[-3];
                                $qline      = $$hsp[-2];
                                $pline      = $$hsp[-1];
                                $lastdomain = $name;
                            }
                            # model data track, some reports don't have
                            elsif ( $_ =~ m/\s+\S+\s(?:CS|RF)$/ ) {
                                my @data = split( " ", $_ );
                                $csline .= $data[-2];
                                $max_count++;
                                $count++;
                                next;
                            }
                            # Query line and Hit line swaps positions
                            # depending of the program
                            elsif (    $count == $max_count - 3
                                   or  $count == $max_count - 1
                                   ) {
                                my @data = split( " ", $_ );

                                my $line_offset = 0;
                                while ($_ =~ m/$data[-2]/g) {
                                    $line_offset = pos;
                                }
                                if ($line_offset != 0) {
                                    $align_length = length $data[-2];
                                    $align_offset = $line_offset - $align_length;
                                }

                                if ($self->{'_reporttype'} eq 'HMMSCAN') {
                                    # hit sequence
                                    $hline .= $data[-2] if ($count == $max_count - 3);
                                    # query sequence
                                    $qline .= $data[-2] if ($count == $max_count - 1);
                                }
                                else { # hmmsearch & nhmmer
                                    # hit sequence
                                    $hline .= $data[-2] if ($count == $max_count - 1);
                                    # query sequence
                                    $qline .= $data[-2] if ($count == $max_count - 3);
                                }

                                $count++;
                                next;
                            }
                            # conservation track
                            # storage isn't quite right - need to remove
                            # leading/lagging whitespace while preserving
                            # gap data (latter isn't done, former is)
                            elsif ( $count == $max_count - 2 ) {
                                $midline .= substr $_, $align_offset, $align_length;
                                $count++;
                                next;
                            }
                            # posterior probability track
                            elsif ( $count == $max_count ) {
                                my @data   = split(" ", $_);
                                $pline    .= $data[-2];
                                $count     = 0;
                                $max_count = 3;
                                $$hsp[-5]  = $csline;
                                $$hsp[-4]  = $hline;
                                $$hsp[-3]  = $midline;
                                $$hsp[-2]  = $qline;
                                $$hsp[-1]  = $pline;
                                next;
                            }
                            else {
                                print "Missed this line: $_\n";
                            }
                        }
                    }
                }
            }

            # End of report
            elsif ( m/Internal pipeline statistics/ || m!^//! ) {
                # If within hit, hsp close;
                if ( $self->within_element('hit') ) {
                    if ( $self->within_element('hsp') ) {
                        $self->end_element( { 'Name' => 'Hsp' } );
                    }
                    $self->end_element( { 'Name' => 'Hit' } );
                }

                # Grab summary statistics of run
                while ( defined( $_ = $self->_readline ) ) {
                    last if ( $_ =~ m/^\/\/$/ );
                }

                # Do a lot of processing of hits and hsps here
                my $index = 0;
                while ( my $hit = shift @hit_list ) {
                    $index++;
                    my $hit_name    = shift @$hit;
                    my $hit_desc    = shift @$hit;
                    my $hit_signif  = shift @$hit;
                    my $hit_score   = shift @$hit;
                    my $num_domains = $domaincounter{"$hit_name.$index"} || 0;

                    $self->start_element( { 'Name' => 'Hit' } );
                    $self->element(
                        {   'Name' => 'Hit_id',
                            'Data' => $hit_name
                        }
                    );
                    $self->element(
                        {   'Name' => 'Hit_desc',
                            'Data' => $hit_desc
                        }
                    );
                    $self->element(
                        {   'Name' => 'Hit_signif',
                            'Data' => $hit_signif
                        }
                    );
                    $self->element(
                        {   'Name' => 'Hit_score',
                            'Data' => $hit_score
                        }
                    );

                    for my $i ( 1 .. $num_domains ) {
                        my $key = $hit_name . "_" . $i;
                        my $hsp = $hsp_list[ $hspinfo{"$key.$index"} ];
                        if ( defined $hsp ) {
                            my $hsp_name = shift @$hsp;
                            $self->start_element( { 'Name' => 'Hsp' } );
                            # Since HMMER doesn't print some data explicitly,
                            # calculate it from the homology line (midline)
                            if ($$hsp[-3] ne '') {
                                my $length    = length $$hsp[-3];
                                my $identical = ($$hsp[-3] =~ tr/a-zA-Z//);
                                my $positive  = ($$hsp[-3] =~ tr/+//) + $identical;
                                $self->element(
                                    {
                                        'Name' => 'Hsp_align-len',
                                        'Data' => $length
                                    }
                                );
                                $self->element(
                                    {   'Name' => 'Hsp_identity',
                                        'Data' => $identical
                                    }
                                );
                                $self->element(
                                    {   'Name' => 'Hsp_positive',
                                        'Data' => $positive
                                    }
                                );
                            }
                            else {
                                $self->element(
                                    {   'Name' => 'Hsp_identity',
                                        'Data' => 0
                                    }
                                );
                                $self->element(
                                    {   'Name' => 'Hsp_positive',
                                        'Data' => 0
                                    }
                                );
                            }
                            if ( $self->{'_reporttype'} eq 'HMMSCAN' ) {
                                $self->element(
                                    {   'Name' => 'Hsp_hit-from',
                                        'Data' => shift @$hsp
                                    }
                                );
                                $self->element(
                                    {   'Name' => 'Hsp_hit-to',
                                        'Data' => shift @$hsp
                                    }
                                );
                                $self->element(
                                    {   'Name' => 'Hsp_query-from',
                                        'Data' => shift @$hsp
                                    }
                                );
                                $self->element(
                                    {   'Name' => 'Hsp_query-to',
                                        'Data' => shift @$hsp
                                    }
                                );
                            }
                            elsif (   $self->{'_reporttype'} eq 'HMMSEARCH'
                                   or $self->{'_reporttype'} eq 'NHMMER'
                                   ) {
                                $self->element(
                                    {   'Name' => 'Hsp_query-from',
                                        'Data' => shift @$hsp
                                    }
                                );
                                $self->element(
                                    {   'Name' => 'Hsp_query-to',
                                        'Data' => shift @$hsp
                                    }
                                );
                                $self->element(
                                    {   'Name' => 'Hsp_hit-from',
                                        'Data' => shift @$hsp
                                    }
                                );
                                $self->element(
                                    {   'Name' => 'Hsp_hit-to',
                                        'Data' => shift @$hsp
                                    }
                                );
                            }
                            $self->element(
                                {   'Name' => 'Hsp_score',
                                    'Data' => shift @$hsp
                                }
                            );
                            $self->element(
                                {   'Name' => 'Hsp_evalue',
                                    'Data' => shift @$hsp
                                }
                            );
                            my $hitlength = shift @$hsp;
                            if ( $hitlength != 0 ) {
                                $self->element(
                                    {   'Name' => 'Hit_len',
                                        'Data' => $hitlength
                                    }
                                );
                            }
                            $self->element(
                                {   'Name' => 'Hsp_csline',
                                    'Data' => shift @$hsp
                                }
                            );
                            $self->element(
                                {   'Name' => 'Hsp_hseq',
                                    'Data' => shift @$hsp
                                }
                            );
                            $self->element(
                                {   'Name' => 'Hsp_midline',
                                    'Data' => shift @$hsp
                                }
                            );
                            $self->element(
                                {   'Name' => 'Hsp_qseq',
                                    'Data' => shift @$hsp
                                }
                            );
                            $self->element(
                                {   'Name' => 'Hsp_pline',
                                    'Data' => shift @$hsp
                                }
                            );

                            # Only nhmmer output has strand information
                            if ( $self->{'_reporttype'} eq 'NHMMER' ) {
                                my $hstart = $self->get_from_element('HSP-hit_start');
                                my $hend   = $self->get_from_element('HSP-hit_end');
                                my $hstrand = ( $hstart < $hend ) ? 1 : -1;

                                my $qstart = $self->get_from_element('HSP-query_start');
                                my $qend   = $self->get_from_element('HSP-query_end');
                                my $qstrand = ( $qstart < $qend ) ? 1 : -1;

                                $self->element(
                                    {   'Name' => 'Hsp_query-strand',
                                        'Data' => $qstrand
                                    }
                                );
                                $self->element(
                                    {   'Name' => 'Hsp_hit-strand',
                                        'Data' => $hstrand
                                    }
                                );
                            }

                            $self->end_element( { 'Name' => 'Hsp' } );
                        }
                    }
                    $self->end_element( { 'Name' => 'Hit' } );
                }
                @hit_list = ();
                %hitinfo  = ();
                last;
            }
        }
        else {
            print "Missed line: $_\n";
            $self->debug($_);
        }
        $last = $_;
    }
    $self->end_element( { 'Name' => 'HMMER_Output' } );
    my $result = $self->end_document();
    return $result;
}

=head2 start_element

 Title   : start_element
 Usage   : $eventgenerator->start_element
 Function: Handles a start event
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

 Title   : end_element
 Usage   : $eventgeneartor->end_element
 Function: Handles and end element event
 Returns : none
 Args    : hashref with at least 2 keys 'Data' and 'Name'

=cut

sub end_element {

    my ( $self, $data ) = @_;
    my $nm   = $data->{'Name'};
    my $type = $MODEMAP{$nm};
    my $rc;

    if ( $nm eq 'HMMER_program' ) {
        if ( $self->{'_last_data'} =~ /(N?HMM\S+)/i ) {
            $self->{'_reporttype'} = uc $1;
        }
    }

    # Hsp are sort of weird, in that they end when another
    # object begins so have to detect this in end_element for now
    if ( $nm eq 'Hsp' ) {
        foreach (qw(Hsp_csline Hsp_qseq Hsp_midline Hsp_hseq Hsp_pline)) {
            my $data = $self->{'_last_hspdata'}->{$_};
            if ( $data && $_ eq 'Hsp_hseq' ) {

                # replace hmm '.' gap symbol by '-'
                $data =~ s/\./-/g;
            }
            $self->element(
                {   'Name' => $_,
                    'Data' => $data
                }
            );
        }
        $self->{'_last_hspdata'} = {};
    }
    if ($type) {
        if ( $self->_eventHandler->will_handle($type) ) {
            my $func = sprintf( "end_%s", lc $type );
            $rc = $self->_eventHandler->$func( $self->{'_reporttype'},
                $self->{'_values'} );
        }
        my $lastelem = shift @{ $self->{'_elements'} };

        # Flush corresponding values from the {_values} buffer
        my $name = uc $type;
        foreach my $key (keys %{ $self->{_values} }) {
            delete $self->{_values}->{$key} if ($key =~ m/^$name-/);
        }
    }
    elsif ( $MAPPING{$nm} ) {
        if ( ref( $MAPPING{$nm} ) =~ /hash/i ) {
            my $key = ( keys %{ $MAPPING{$nm} } )[0];
            $self->{'_values'}->{$key}->{ $MAPPING{$nm}->{$key} }
                = $self->{'_last_data'};
        }
        else {
            $self->{'_values'}->{ $MAPPING{$nm} } = $self->{'_last_data'};

            # print "lastdata is " . $self->{'_last_data'} . "\n";
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
    $self->start_element($data);
    $self->characters($data);
    $self->end_element($data);
}

=head2 get_from_element

 Title   : get_from_element
 Usage   : $self->get_from_element('HSP-hit_start');
 Function: Convenience method to retrieve data from '_values' hash
 Returns : string
 Args    : key

=cut

sub get_from_element {
    my ($self,$key) = @_;
    my $values = $self->{_values};
    $values->{$key};
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

    if (   $self->in_element('hsp')
        && $data->{'Name'} =~ /Hsp\_(?:qseq|hseq|csline|pline|midline)/o
        && defined $data->{'Data'} )
    {
        $self->{'_last_hspdata'}->{ $data->{'Name'} } .= $data->{'Data'};
    }
    return unless ( defined $data->{'Data'} && $data->{'Data'} !~ /^\s+$/o );

    $self->{'_last_data'} = $data->{'Data'};
}

=head2 within_element

 Title   : within_element
 Usage   : if( $eventgenerator->within_element( $element ) ) {}
 Function: Test if we are within a particular element
           This is different than 'in' because within can be tested for
           a whole block
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
 Usage   : if( $eventgenerator->in_element( $element ) ) {}
 Function: Test if we are in a particular element
           This is different than 'within' because 'in' only
           tests its immediate parent
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
 Function: Returns the number of results processed
 Returns : interger
 Args    : none

=cut

sub result_count {
    my $self = shift;
    return $self->{'_result_count'};
}

1;
