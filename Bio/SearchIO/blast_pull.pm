#
# BioPerl module for Bio::SearchIO::blast_pull
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

Bio::SearchIO::blast_pull - A parser for BLAST output

=head1 SYNOPSIS

    # do not use this class directly it is available through Bio::SearchIO
    use Bio::SearchIO;
    my $in = Bio::SearchIO->new(-format => 'blast_pull',
                               -file   => 't/data/new_blastn.txt');
    while (my $result = $in->next_result) {
        # this is a Bio::Search::Result::BlastPullResult object
        print "Results for ", $result->query_name(), "\n";
        while (my $hit = $result->next_hit) {
            print $hit->name(), "\n";
            while (my $hsp = $hit->next_hsp) {
                print "length is ", $hsp->length(), "\n";
            }
        }
    }

=head1 DESCRIPTION

This object implements a pull-parser for BLAST output. It is fast since it
only does work on request (hence 'pull').

Currently only NCBI BLASTN and BLASTP are supported.

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

package Bio::SearchIO::blast_pull;

use strict;
use Bio::Search::Result::BlastPullResult;

use base qw(Bio::SearchIO Bio::PullParserI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::blast_pull->new();
 Function: Builds a new Bio::SearchIO::blast_pull object 
 Returns : Bio::SearchIO::blast_pull
 Args    : -fh/-file => BLAST output filename
           -format   => 'blast_pull'
           -evalue   => float or scientific notation number to be used
                        as an evalue cutoff for hits
           -score    => integer or scientific notation number to be used
                        as a score value cutoff for hits
           -piped_behaviour => 'temp_file'|'memory'|'sequential_read'

           -piped_behaviour defines what the parser should do if the input is
            an unseekable filehandle (eg. piped input), see
            Bio::PullParserI::chunk for details. Default is 'memory'.

=cut

sub _initialize {
    my ($self, @args) = @_;
    
    # don't do normal SearchIO initialization
    
    my ($writer, $file, $fh, $piped_behaviour, $evalue, $score) =
                            $self->_rearrange([qw(WRITER
                                                  FILE FH
                                                  PIPED_BEHAVIOUR
                                                  EVALUE
                                                  SCORE)], @args);
    $self->writer($writer) if $writer;
    
    $self->_fields( { ( header => undef,
                        algorithm => undef,
                        algorithm_version => undef,
                        algorithm_reference => '',
                        database_name => undef,
                        database_letters => undef,
                        database_entries => undef,
                        next_result => undef,
                        evalue_cutoff => '[unset]',
                        score_cutoff => '[unset]' ) } );
    
    $self->_fields->{evalue_cutoff} = $evalue if $evalue;
    $self->_fields->{score_cutoff} = $score if $score;
    
    $self->_dependencies( { ( algorithm => 'header',
                              algorithm_version => 'header',
                              database_name => 'header',
                              database_letters => 'header',
                              database_entries => 'header' ) } );
    
    $self->chunk($file || $fh || $self->throw("-file or -fh must be supplied"),
                 -piped_behaviour => $piped_behaviour || 'memory');
}

sub _discover_header {
    my $self = shift;
    $self->_chunk_seek(0);
    my $header = $self->_get_chunk_by_end("\nQuery=");
    $self->{_after_header} = $self->_chunk_tell;
    
    #*** won't catch all types? only support blastn/p now anyway
    $header =~ /^(\S+) (\S+\s+\S+)/;
    $self->_fields->{algorithm} = $1;
    $self->_fields->{algorithm_version} = $2;
    
    my ($database) = $header =~ /^Database: (.+)/sm;
    
    unless ($database) {
        # earlier versions put query before database?
        my $header2 = $self->_get_chunk_by_end(".done\n");
        ($database) = $header2 =~ /^Database: (.+)/sm;
    }
    
    $database =~ s/\s+(\d\S+) sequences; (\d\S+) total letters.*//s;
    my $entries = $1;
    my $letters = $2;
    $database =~ s/\n//g;
    $entries =~ s/,//g;
    $letters =~ s/,//g;
    $self->_fields->{database_name} = $database;
    $self->_fields->{database_entries} = $entries;
    $self->_fields->{database_letters} = $letters;
    
    $self->_fields->{header} = 1;
}

sub _discover_next_result {
    my $self = shift;
    return if $self->{_after_results};
    my $type = $self->get_field('algorithm'); # also sets _after_header if not set
    
    if ($type eq 'BLASTN' || $type eq 'BLASTP') {
        unless ($self->_sequential) {
            $self->_chunk_seek($self->{_end_of_previous_result} || $self->{_after_header});
            
            my ($start, $end) = $self->_find_chunk_by_end("\nQuery=");
            return if ($start == $end);
            
            unless ($end) {
                $start = $self->{_end_of_previous_result} || $self->{_after_header};
                $end = undef;
            }
            
            $self->_fields->{next_result} = Bio::Search::Result::BlastPullResult->new(-chunk => [($self->chunk, $start, $end)],
                                                                                     -parent => $self);
            
            $self->{_end_of_previous_result} = $end;
        }
        else {
            #*** doesn't work for the last result, needs fixing - try getting the database end chunk on failure?...
            $self->throw("sequential mode not yet implemented");
            my $chunk = $self->_get_chunk_by_end("\nQuery=");
            $chunk || return;
            $self->_fields->{next_result} = Bio::Search::Result::BlastPullResult->new(-chunk => [$chunk],
                                                                                   -parent => $self);
        }
    }
    else {
        $self->throw("Can only handle NCBI BLASTN and BLASTP right now");
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
	delete $self->{_end_of_previous_result};
}

1;
