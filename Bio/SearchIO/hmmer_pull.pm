#
# BioPerl module for Bio::SearchIO::hmmer_pull
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::hmmer_pull - A parser for HMMER output

=head1 SYNOPSIS

    # do not use this class directly it is available through Bio::SearchIO
    use Bio::SearchIO;
    my $in = Bio::SearchIO->new(-format => 'hmmer_pull',
                               -file   => 't/data/hmmpfam.bigout');
    while (my $result = $in->next_result) {
        # this is a Bio::Search::Result::HmmpfamResult object
        print $result->query_name(), " for HMM ", $result->hmm_name(), "\n";
        while (my $hit = $result->next_hit) {
            print $hit->name(), "\n";
            while (my $hsp = $hit->next_hsp) {
                print "length is ", $hsp->length(), "\n";
            }
        }
    }

=head1 DESCRIPTION

This object implements a pull-parser for HMMER output. It is fast since it
only does work on request (hence 'pull').

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

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SearchIO::hmmer_pull;

use strict;


use base qw(Bio::SearchIO Bio::PullParserI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::hmmer_pull->new();
 Function: Builds a new Bio::SearchIO::hmmer_pull object 
 Returns : Bio::SearchIO::hmmer_pull
 Args    : -fh/-file => HMMER output filename
           -format   => 'hmmer_pull'
           -evalue   => float or scientific notation number to be used
                        as an evalue cutoff for hits
           -score    => integer or scientific notation number to be used
                        as a score value cutoff for hits
           -hsps     => integer minimum number of hsps (domains) a hit must have
           -piped_behaviour => 'temp_file'|'memory'|'sequential_read'

           -piped_behaviour defines what the parser should do if the input is
            an unseekable filehandle (eg. piped input), see
            Bio::PullParserI::chunk for details. Default is 'sequential_read'.

=cut

sub _initialize {
    my ($self, @args) = @_;
    
    # don't do normal SearchIO initialization
    
    my ($writer, $file, $fh, $piped_behaviour, $evalue, $score, $hsps) =
                            $self->_rearrange([qw(WRITER
                                                  FILE FH
                                                  PIPED_BEHAVIOUR
                                                  EVALUE
                                                  SCORE
                                                  HSPS)], @args);
    $self->writer($writer) if $writer;
    
    $self->_fields( { ( header => undef,
                        algorithm => undef,
                        algorithm_version => undef,
                        algorithm_reference => '',
                        hmm_file => undef,
                        hmm_name => undef,
                        sequence_file => undef,
                        sequence_database => undef,
                        database_name => undef,
                        database_letters => undef,
                        database_entries => undef,
                        next_result => undef,
                        evalue_cutoff => '[unset]',
                        score_cutoff => '[unset]',
                        hsps_cutoff => '[unset]' ) } );
    
    $self->_fields->{evalue_cutoff} = $evalue if $evalue;
    $self->_fields->{score_cutoff} = $score if $score;
    $self->_fields->{hsps_cutoff} = $hsps if $hsps;
    
    $self->_dependencies( { ( algorithm => 'header',
                              algorithm_version => 'header',
                              hmm_file => 'header',
                              hmm_name => 'header',
                              sequence_file => 'header',
                              sequence_database => 'header' ) } );
    
    $self->chunk($file || $fh || $self->throw("-file or -fh must be supplied"),
                 -piped_behaviour => $piped_behaviour || 'sequential_read');
}

sub _discover_header {
    my $self = shift;
    $self->_chunk_seek(0);
    my $header = $self->_get_chunk_by_nol(8);
    $self->{_after_header} = $self->_chunk_tell;
    
    my ($algo) = $header =~ /^(hmm\S+) - search/m;
    $self->_fields->{algorithm} = uc $algo;
    
    ($self->_fields->{algorithm_version}) = $header =~ /^HMMER\s+?(\S+)/m;
    
    ($self->_fields->{hmm_file}) = $header =~ /^HMM file:\s.+?(\S+)$/m;
    $self->_fields->{hmm_name} = $self->_fields->{hmm_file};
    
    ($self->_fields->{sequence_file}) = $header =~ /^Sequence (?:file|database):\s.+?(\S+)$/m;
    $self->_fields->{sequence_database} = $self->_fields->{sequence_file};
    
    $self->_fields->{header} = 1;
}

sub _discover_database_name {
    my $self = shift;
    my $type = $self->get_field('algorithm');
    
    if ($type eq 'HMMPFAM') {
        $self->_fields->{database_name} = $self->get_field('hmm_file');
    }
    elsif ($type eq 'HMMSEARCH') {
        $self->_fields->{database_name} = $self->get_field('sequence_file');
    }
}

sub _discover_next_result {
    my $self = shift;
    my $type = $self->get_field('algorithm'); # also sets _after_header if not set
    
    if ($type eq 'HMMPFAM') {
        use Bio::Search::Result::HmmpfamResult;
        
        unless ($self->_sequential) {
            $self->_chunk_seek($self->{_end_of_previous_result} || $self->{_after_header});
            
            my ($start, $end) = $self->_find_chunk_by_end("//\n");
            return if $start == $end;
            $self->_fields->{next_result} = Bio::Search::Result::HmmpfamResult->new(-chunk => [($self->chunk, $start, $end)],
                                                                                   -parent => $self);
            
            $self->{_end_of_previous_result} = $end;
        }
        else {
            # deliberatly don't cache these, which means rewind won't work;
            # if we cached we may as well have used 'memory' option to
            # -piped_behaviour
            my $chunk = $self->_get_chunk_by_end("//\n");
            $chunk || return;
            $self->_fields->{next_result} = Bio::Search::Result::HmmpfamResult->new(-chunk => [$chunk],
                                                                                   -parent => $self);
        }
    }
    elsif ($type eq 'HMMSEARCH') {
        $self->throw("Can't handle hmmsearch yet\n");
    }
    else {
        $self->throw("Unknown report type");
    }
}

=head2 next_result

 Title   : next_result
 Usage   : my $hit = $searchio->next_result;
 Function: Returns the next Result from a search
 Returns : Bio::Search::Result::ResultI object
 Args    : none

=cut

sub next_result {
    my $self = shift;
    my $result = $self->get_field('next_result') || return;
    
    undef $self->_fields->{next_result};
    
    $self->{'_result_count'}++;
    return $result;
}

=head2 result_count

 Title   : result_count
 Usage   : my $count = $searchio->result_count
 Function: Returns the number of results we have processed.
 Returns : integer
 Args    : none

=cut

sub result_count {
    my $self = shift;
    return $self->{'_result_count'};
}

=head2 rewind

 Title   : rewind
 Usage   : $searchio->rewind;
 Function: Allow one to reset the Result iterator to the beginning, so that
           next_result() will subsequently return the first result and so on.

           NB: result objects are not cached, so you will get new result objects
           each time you rewind. Also, note that result_count() counts the
           number of times you have called next_result(), so will not be able
           tell you how many results there were in the file if you use rewind().

 Returns : n/a
 Args    : none

=cut

sub rewind {
	my $self = shift;
    if ($self->_sequential) {
        $self->warn("rewind has no effect on piped input when you have chosen 'sequential_read' mode");
    }
	delete $self->{_end_of_previous_result};
}

1;
