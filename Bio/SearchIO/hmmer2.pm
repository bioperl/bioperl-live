#
# BioPerl module for Bio::SearchIO::hmmer2
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::hmmer2 - A parser for HMMER output (hmmpfam, hmmsearch)

=head1 SYNOPSIS

    # do not use this class directly it is available through Bio::SearchIO
    use Bio::SearchIO;
    my $in = Bio::SearchIO->new(-format => 'hmmer2',
                               -file   => 't/data/L77119.hmmer');
    while( my $result = $in->next_result ) {
        # this is a Bio::Search::Result::HMMERResult object
        print $result->query_name(), " for HMM ", $result->hmm_name(), "\n";
        while( my $hit = $result->next_hit ) {
            print $hit->name(), "\n";
            while( my $hsp = $hit->next_hsp ) {
                print "length is ", $hsp->length(), "\n";
            }
        }
    }

=head1 DESCRIPTION

This object implements a parser for HMMER output.

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

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SearchIO::hmmer2;

use strict;

use Bio::Factory::ObjectFactory;

use vars qw(%MAPPING %MODEMAP
);

use base qw(Bio::SearchIO::hmmer);

BEGIN {

    # mapping of HMMER items to Bioperl hash keys
    %MODEMAP = (
        'HMMER_Output' => 'result',
        'Hit'          => 'hit',
        'Hsp'          => 'hsp'
    );

    %MAPPING = (
        'Hsp_bit-score'   => 'HSP-bits',
        'Hsp_score'       => 'HSP-score',
        'Hsp_evalue'      => 'HSP-evalue',
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
        'Hsp_csline'      => 'HSP-cs_seq',
        'Hsp_hseq'        => 'HSP-hit_seq',
        'Hsp_midline'     => 'HSP-homology_seq',
        'Hsp_align-len'   => 'HSP-hsp_length',
        'Hsp_query-frame' => 'HSP-query_frame',
        'Hsp_hit-frame'   => 'HSP-hit_frame',

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
    my $seentop = 0;
    my $reporttype;
    my ( $last, @hitinfo, @hspinfo, %hspinfo, %hitinfo );
    local $/ = "\n";
    local $_;

    my $verbose = $self->verbose;    # cache for speed?
    $self->start_document();
    local ($_);
    while ( defined( $_ = $self->_readline ) ) {
        my $lineorig = $_;

        chomp;
        if (/^HMMER\s+(\S+)\s+\((.+)\)/o) {
            my ( $prog, $version ) = split;
            if ($seentop) {
                $self->_pushback($_);
                $self->end_element( { 'Name' => 'HMMER_Output' } );
                return $self->end_document();
            }
            $self->{'_hmmidline'} = $_;
            $self->start_element( { 'Name' => 'HMMER_Output' } );
            $self->{'_result_count'}++;
            $seentop = 1;
            if ( defined $last ) {
                ($reporttype) = split( /\s+/, $last );
                $reporttype = uc($reporttype) if defined $reporttype;
                $self->element(
                    {
                        'Name' => 'HMMER_program',
                        'Data' => $reporttype
                    }
                );
            }
            $self->element(
                {
                    'Name' => 'HMMER_version',
                    'Data' => $version
                }
            );
        }
        elsif (s/^HMM file:\s+//o) {
            $self->{'_hmmfileline'} = $lineorig;
            $self->element(
                {
                    'Name' => 'HMMER_hmm',
                    'Data' => $_
                }
            );
        }
        elsif (s/^Sequence\s+(file|database):\s+//o) {
            $self->{'_hmmseqline'} = $lineorig;
            if ( $1 eq 'database' ) {
                $self->element(
                    {
                        'Name' => 'HMMER_db',
                        'Data' => $_
                    }
                );
            }
            $self->element(
                {
                    'Name' => 'HMMER_seqfile',
                    'Data' => $_
                }
            );
        }
        elsif (s/^Query(?:\s+(?:sequence|HMM))?(?:\s+\d+)?:\s+//o) {
            if ( !$seentop ) {

                # we're in a multi-query report
                $self->_pushback($lineorig);
                $self->_pushback( $self->{'_hmmseqline'} );
                $self->_pushback( $self->{'_hmmfileline'} );
                $self->_pushback( $self->{'_hmmidline'} );
                next;
            }
            s/\s+$//;
            $self->element(
                {
                    'Name' => 'HMMER_query-def',
                    'Data' => $_
                }
            );
        }
        elsif (s/^Accession:\s+//o) {
            s/\s+$//;
            $self->element(
                {
                    'Name' => 'HMMER_query-acc',
                    'Data' => $_
                }
            );
        }
        elsif (s/^Description:\s+//o) {
            s/\s+$//;
            $self->element(
                {
                    'Name' => 'HMMER_querydesc',
                    'Data' => $_
                }
            );
        }
        elsif (   defined $self->{'_reporttype'}
              && (   $self->{'_reporttype'} eq 'HMMSEARCH'
                  || $self->{'_reporttype'} eq 'HMMPFAM' )
              ) {
            # PROCESS RESULTS HERE
            if (/^Scores for (?:complete sequences|sequence family)/o) {
                while ( defined( $_ = $self->_readline ) ) {
                    last if (/^\s+$/);
                    next if (   /^Model\s+Description/o
                             || /^Sequence\s+Description/o
                             || /^\-\-\-/o );

                    chomp;
                    my @line = split;
                    my ( $name, $domaintotal, $evalue, $score ) =
                      ( shift @line, pop @line, pop @line, pop @line );
                    my $desc = join( ' ', @line );
                    push @hitinfo, [ $name, $desc, $score, $evalue, $domaintotal ];
                    $hitinfo{$name} = $#hitinfo;
                }
            }
            elsif (/^Parsed for domains:/o) {
                @hspinfo = ();

                while ( defined( $_ = $self->_readline ) ) {
                    last if (/^\s+$/);
                    if (m!^//!) {
                        $self->_pushback($_);
                        last;
                    }
                    next if ( /^(?:Model|Sequence)\s+Domain/ || /^\-\-\-/ );

                    chomp;
                    if (
                        my ( $name, $domainct, $domaintotal,
                             $seq_start, $seq_end, $seq_cov,
                             $hmm_start, $hmm_end, $hmm_cov,
                             $score, $evalue ) = (
                                m!^(\S+)\s+          # domain name
                                   (\d+)/(\d+)\s+    # domain num out of num
                                   (\d+)\s+(\d+)\s+  # seq start, end
                                   (\S+)\s+          # seq coverage
                                   (\d+)\s+(\d+)\s+  # hmm start, end
                                   (\S+)\s+          # hmm coverage
                                   (\S+)\s+          # score
                                   (\S+)             # evalue
                                    \s*$!ox
                            )
                        ) {
                        my $hindex = $hitinfo{$name};
                        if ( !defined $hindex ) {
                            push @hitinfo,
                              [ $name, '', $score, $evalue, $domaintotal ];
                            $hitinfo{$name} = $#hitinfo;
                            $hindex = $#hitinfo;
                        }

                        my $info = $hitinfo[$hindex];
                        if ( !defined $info ) {
                            if ($self->{'_reporttype'} eq 'HMMSEARCH') {
                                $self->warn(
                                    "Incomplete Sequence information, can't find $name hitinfo says $hitinfo{$name}"
                                );
                            }
                            elsif ($self->{'_reporttype'} eq 'HMMPFAM') {
                                $self->warn(
                                    "Incomplete Domain information, can't find $name hitinfo says $hitinfo{$name}"
                                );
                            }
                            next;
                        }

                        # Try to get HMM and Sequence lengths from the alignment information
                        if ($self->{'_reporttype'} eq 'HMMSEARCH') {
                            # For Hmmsearch, if seq coverage ends in ']' it means that the alignment
                            # runs until the end. In that case add the END coordinate to @hitinfo
                            # to use it as Hit Length
                            if (    $seq_cov =~ m/\]$/
                                and scalar @{ $hitinfo[$hindex] } == 5
                                ) {
                                push @{ $hitinfo[$hindex] }, $seq_end ;
                            }
                            # For Hmmsearch, if hmm coverage ends in ']', it means that the alignment
                            # runs until the end. In that case use the END coordinate as Query Length
                            if (    $hmm_cov =~ m/\]$/
                                and not exists $self->{_values}->{'RESULT-query_length'}
                                ) {
                                $self->element(
                                    {   'Name' => 'HMMER_query-len',
                                        'Data' => $hmm_end
                                    }
                                );
                            }
                        }
                        elsif ($self->{'_reporttype'} eq 'HMMPFAM') {
                            # For Hmmpfam, if hmm coverage ends in ']' it means that the alignment
                            # runs until the end. In that case add the END coordinate to @hitinfo
                            # to use it as Hit Length
                            if (    $hmm_cov =~ m/\]$/
                                and scalar @{ $hitinfo[$hindex] } == 5
                                ) {
                                push @{ $hitinfo[$hindex] }, $hmm_end ;
                            }
                            # For Hmmpfam, if seq coverage ends in ']', it means that the alignment
                            # runs until the end. In that case use the END coordinate as Query Length
                            if (    $seq_cov =~ m/\]$/
                                and not exists $self->{_values}->{'RESULT-query_length'}
                                ) {
                                $self->element(
                                    {   'Name' => 'HMMER_query-len',
                                        'Data' => $seq_end
                                    }
                                );
                            }
                        }

                        my @vals = ($seq_start, $seq_end,
                                    $hmm_start, $hmm_end,
                                    $score,     $evalue);
                        push @hspinfo, [ $name, @vals ];
                    }
                }
            }
            elsif (/^Alignments of top/o) {
                my ( $prelength, $count, $width );
                $count = 0;
                my %domaincounter;
                my $second_tier = 0;
                my $csline      = '';

                while ( defined( $_ = $self->_readline ) ) {
                    next if ( /^Align/o );

                    if (   m/^Histogram/o
                        || m!^//!o
                        || m/^Query(?:\s+(?:sequence|HMM))?(?:\s+\d+)?:/o
                        ) {
                        if ( $self->in_element('hsp') ) {
                            $self->end_element( { 'Name' => 'Hsp' } );
                        }
                        if ( $self->within_element('hit') ) {
                            $self->end_element( { 'Name' => 'Hit' } );
                        }
                        $self->_pushback($_);
                        last;
                    }

                    chomp;
                    if (
                        my ( $name, $domainct, $domaintotal,
                             $from, $to ) = (
                                m/^\s*(.+):
                                   \s+ domain \s+ (\d+) \s+ of \s+ (\d+) ,
                                   \s+ from   \s+ (\d+) \s+ to \s+ (\d+)/x
                            )
                        ) {
                        $domaincounter{$name}++;
                        if ( $self->within_element('hit') ) {
                            if ( $self->within_element('hsp') ) {
                                $self->end_element( { 'Name' => 'Hsp' } );
                            }
                            $self->end_element( { 'Name' => 'Hit' } );
                        }

                        my $info = [ @{ $hitinfo[ $hitinfo{$name} ] } ];
                        if (   !defined $info
                            || $info->[0] ne $name
                            ) {
                            $self->warn(
                                  "Somehow the Model table order does not match the order in the domains (got "
                                .  $info->[0]
                                . ", expected $name). We're back loading this from the alignment information instead"
                            );
                            $info = [
                                $name, '',
                                /score \s+ ([^,\s]+), \s+E\s+=\s+ (\S+)/ox,
                                $domaintotal
                            ];
                            push @hitinfo, $info;
                            $hitinfo{$name} = $#hitinfo;
                        }

                        $self->start_element( { 'Name' => 'Hit' } );
                        $self->element(
                            {
                                'Name' => 'Hit_id',
                                'Data' => shift @{$info}
                            }
                        );
                        $self->element(
                            {
                                'Name' => 'Hit_desc',
                                'Data' => shift @{$info}
                            }
                        );
                        $self->element(
                            {
                                'Name' => 'Hit_score',
                                'Data' => shift @{$info}
                            }
                        );
                        $self->element(
                            {
                                'Name' => 'Hit_signif',
                                'Data' => shift @{$info}
                            }
                        );
                        my $dom_total = shift @{$info};
                        if (my $hit_end = shift @{$info}) {
                            $self->element(
                                {
                                    'Name' => 'Hit_len',
                                    'Data' => $hit_end
                                }
                            );
                        }

                        $self->start_element( { 'Name' => 'Hsp' } );
                        my $HSPinfo = shift @hspinfo;
                        my $id      = shift @$HSPinfo;

                        if ( $id ne $name ) {
                            $self->throw(
                                  "Somehow the domain list details do not match "
                                . "the table (got $id, expected $name)"
                            );
                        }

                        if ($self->{'_reporttype'} eq 'HMMSEARCH') {
                            $self->element(
                                {
                                    'Name' => 'Hsp_hit-from',
                                    'Data' => shift @$HSPinfo
                                }
                            );
                            $self->element(
                                {
                                    'Name' => 'Hsp_hit-to',
                                    'Data' => shift @$HSPinfo
                                }
                            );
                            $self->element(
                                {
                                    'Name' => 'Hsp_query-from',
                                    'Data' => shift @$HSPinfo
                                }
                            );
                            $self->element(
                                {
                                    'Name' => 'Hsp_query-to',
                                    'Data' => shift @$HSPinfo
                                }
                            );
                        }
                        elsif ($self->{'_reporttype'} eq 'HMMPFAM') {
                            $self->element(
                                {
                                    'Name' => 'Hsp_query-from',
                                    'Data' => shift @$HSPinfo
                                }
                            );
                            $self->element(
                                {
                                    'Name' => 'Hsp_query-to',
                                    'Data' => shift @$HSPinfo
                                }
                            );
                            $self->element(
                                {
                                    'Name' => 'Hsp_hit-from',
                                    'Data' => shift @$HSPinfo
                                }
                            );
                            $self->element(
                                {
                                    'Name' => 'Hsp_hit-to',
                                    'Data' => shift @$HSPinfo
                                }
                            );
                        }
                        $self->element(
                            {
                                'Name' => 'Hsp_score',
                                'Data' => shift @$HSPinfo
                            }
                        );
                        $self->element(
                            {
                                'Name' => 'Hsp_evalue',
                                'Data' => shift @$HSPinfo
                            }
                        );

                        if ( $domaincounter{$name} == $domaintotal ) {
                            $hitinfo[ $hitinfo{$name} ] = undef;
                        }
                    }
                    else {

                        # Might want to change this so that it
                        # accumulates all the of the alignment lines into
                        # three array slots and then tests for the
                        # end of the line
                        if ($_ =~ m/^\s+(?:CS|RF)\s+/o && $count == 0) {
                            # Buffer the CS line now and process it later at
                            # midline point, where $prelength and width will be known
                            $csline = $_;
                            next;
                        }
                        elsif (/^(\s+ \*->) (\S+)/ox) {
                            # start of domain
                            $prelength = CORE::length($1);
                            $width     = 0;

                            # deal with fact that start and stop is on same line
                            my $data = $2;
                            if ($data =~ s/<-?\*?\s*$//)
                            {
                                $width = CORE::length($data);
                            }
 
                            if ($self->{'_reporttype'} eq 'HMMSEARCH') {
                                $self->element(
                                    {
                                        'Name' => 'Hsp_qseq',
                                        'Data' => $data
                                    }
                                );
                            }
                            elsif ($self->{'_reporttype'} eq 'HMMPFAM') {
                                $self->element(
                                    {
                                        'Name' => 'Hsp_hseq',
                                        'Data' => $data
                                    }
                                );
                            }
                            $count       = 0;
                            $second_tier = 0;
                        }
                        elsif (/^(\s+) (\S+) <-?\*? \s*$/ox) {
                            # end of domain
                            $prelength -= 3 unless ( $second_tier++ );
                            if ($self->{'_reporttype'} eq 'HMMSEARCH') {
                                $self->element(
                                    {
                                        'Name' => 'Hsp_qseq',
                                        'Data' => $2
                                    }
                                );
                            }
                            elsif ($self->{'_reporttype'} eq 'HMMPFAM') {
                                $self->element(
                                    {
                                        'Name' => 'Hsp_hseq',
                                        'Data' => $2
                                    }
                                );
                            }
                            $width = CORE::length($2);
                            $count = 0;
                        }
                        elsif ( ( $count != 1 && /^\s+$/o )
                               || CORE::length($_) == 0
                               || /^\s+\-?\*\s*$/
                               || /^\s+\S+\s+\-\s+\-\s*$/ )
                        {
                            next;
                        }
                        elsif ( $count == 0 ) {
                            $prelength -= 3 unless ( $second_tier++ );
                            unless ( defined $prelength ) {

                                # $self->warn("prelength not set");
                                next;
                            }
                            if ($self->{'_reporttype'} eq 'HMMSEARCH') {
                                $self->element(
                                    {
                                        'Name' => 'Hsp_qseq',
                                        'Data' => substr( $_, $prelength )
                                    }
                                );
                            }
                            elsif ($self->{'_reporttype'} eq 'HMMPFAM') {
                                $self->element(
                                    {
                                        'Name' => 'Hsp_hseq',
                                        'Data' => substr( $_, $prelength )
                                    }
                                );
                            }
                        }
                        elsif ( $count == 1 ) {
                            if ( !defined $prelength ) {
                                $self->warn("prelength not set");
                            }
                            if ($width) {
                                $self->element(
                                    {
                                        'Name' => 'Hsp_midline',
                                        'Data' => substr( $_, $prelength, $width )
                                    }
                                );
                                if ($csline ne '') {
                                    $self->element(
                                        {
                                            'Name' => 'Hsp_csline',
                                            'Data' => substr( $csline, $prelength, $width )

                                        }
                                    );
                                    $csline = '';
                                }
                            }
                            else {
                                $self->element(
                                    {
                                        'Name' => 'Hsp_midline',
                                        'Data' => substr( $_, $prelength )
                                    }
                                );
                                if ($csline ne '') {
                                    $self->element(
                                        {
                                            'Name' => 'Hsp_csline',
                                            'Data' => substr( $csline, $prelength )
                                        }
                                    );
                                    $csline = '';
                                }
                            }
                        }
                        elsif ( $count == 2 ) {
                            if (   /^\s+(\S+)\s+(\d+|\-)\s+(\S*)\s+(\d+|\-)/o) {
                                if ($self->{'_reporttype'} eq 'HMMSEARCH') {
                                    $self->element(
                                        {
                                            'Name' => 'Hsp_hseq',
                                            'Data' => $3
                                        }
                                    );
                                }
                                elsif ($self->{'_reporttype'} eq 'HMMPFAM') {
                                    $self->element(
                                        {
                                            'Name' => 'Hsp_qseq',
                                            'Data' => $3
                                        }
                                    );
                                }
                            }
                            else {
                                $self->warn("unrecognized line ($count): $_\n");
                            }
                        }
                        $count = 0 if $count++ >= 2;
                    }
                }
            }
            elsif ( /^Histogram/o || m!^//!o ) {
                my %domaincounter;

                while ( my $HSPinfo = shift @hspinfo ) {
                    my $id   = shift @$HSPinfo;
                    $domaincounter{$id}++;

                    my $info = [ @{ $hitinfo[ $hitinfo{$id} ] } ];
                    next unless defined $info;

                    $self->start_element( { 'Name' => 'Hit' } );
                    $self->element(
                        {
                            'Name' => 'Hit_id',
                            'Data' => shift @{$info}
                        }
                    );
                    $self->element(
                        {
                            'Name' => 'Hit_desc',
                            'Data' => shift @{$info}
                        }
                    );
                    $self->element(
                        {
                            'Name' => 'Hit_score',
                            'Data' => shift @{$info}
                        }
                    );
                    $self->element(
                        {
                            'Name' => 'Hit_signif',
                            'Data' => shift @{$info}
                        }
                    );
                    my $domaintotal = shift @{$info};
                    if (my $hit_end = shift @{$info}) {
                        $self->element(
                            {
                                'Name' => 'Hit_len',
                                'Data' => $hit_end
                            }
                        );
                    }

                    # Histogram is exclusive of Hmmsearch, not found in Hmmpfam,
                    # so just use Hmmsearch start/end order (first hit, then query)
                    $self->start_element( { 'Name' => 'Hsp' } );
                    $self->element(
                        {
                            'Name' => 'Hsp_hit-from',
                            'Data' => shift @$HSPinfo
                        }
                    );
                    $self->element(
                        {
                            'Name' => 'Hsp_hit-to',
                            'Data' => shift @$HSPinfo
                        }
                    );
                    $self->element(
                        {
                            'Name' => 'Hsp_query-from',
                            'Data' => shift @$HSPinfo
                        }
                    );
                    $self->element(
                        {
                            'Name' => 'Hsp_query-to',
                            'Data' => shift @$HSPinfo
                        }
                    );
                    $self->element(
                        {
                            'Name' => 'Hsp_score',
                            'Data' => shift @$HSPinfo
                        }
                    );
                    $self->element(
                        {
                            'Name' => 'Hsp_evalue',
                            'Data' => shift @$HSPinfo
                        }
                    );
                    $self->end_element( { 'Name' => 'Hsp' } );
                    $self->end_element( { 'Name' => 'Hit' } );

                    if ( $domaincounter{$id} == $domaintotal ) {
                        $hitinfo[ $hitinfo{$id} ] = undef;
                    }
                }
                @hitinfo = ();
                %hitinfo = ();
                last;
            }
            # uncomment to see missed lines with verbose on
            #else {
            #    $self->debug($_);
            #}
        }
        $last = $_;
    }
    $self->end_element( { 'Name' => 'HMMER_Output' } ) unless !$seentop;
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
 Args    : hashref with at least 2 keys 'Data' and 'Name'


=cut

sub end_element {
    my ( $self, $data ) = @_;
    my $nm   = $data->{'Name'};
    my $type = $MODEMAP{$nm};
    my $rc;

    if ( $nm eq 'HMMER_program' ) {
        if ( $self->{'_last_data'} =~ /(HMM\S+)/i ) {
            $self->{'_reporttype'} = uc $1;
        }
    }

    # Hsp are sort of weird, in that they end when another
    # object begins so have to detect this in end_element for now
    if ( $nm eq 'Hsp' ) {
        foreach (qw(Hsp_csline Hsp_qseq Hsp_midline Hsp_hseq)) {
            my $data = $self->{'_last_hspdata'}->{$_};
            if ($data && $_ eq 'Hsp_hseq') {
                # replace hmm '.' gap symbol by '-'
                $data =~ s/\./-/g;
            }
            $self->element(
                {
                    'Name' => $_,
                    'Data' => $data
                }
            );
            # Since HMMER doesn't print some data explicitly,
            # calculate it from the homology line (midline)
            if ($_ eq 'Hsp_midline') {
                if ($data) {
                    my $length    = length $data;
                    my $identical = ($data =~ tr/a-zA-Z//);
                    my $positive  = ($data =~ tr/+//) + $identical;
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
            }
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

    if (   $self->in_element('hsp')
        && $data->{'Name'} =~ /Hsp\_(?:qseq|hseq|csline|midline)/o
        && defined $data->{'Data'} )
    {
        $self->{'_last_hspdata'}->{ $data->{'Name'} } .= $data->{'Data'};
    }
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
           tests its immediete parent.
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

1;
