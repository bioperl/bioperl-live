#
# BioPerl module for Bio::Search::Result::BlastPullResult
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

Bio::Search::Result::BlastPullResult - A parser and result object for BLASTN
                                     results

=head1 SYNOPSIS

    # generally we use Bio::SearchIO to build these objects
    use Bio::SearchIO;
    my $in = Bio::SearchIO->new(-format => 'blast_pull',
							   -file   => 'result.blast');

    while (my $result = $in->next_result) {
		print $result->query_name, " ", $result->algorithm, " ", $result->num_hits(), " hits\n";
    }

=head1 DESCRIPTION

This object implements a parser for NCBI BLASTN result output.

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

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Search::Result::BlastPullResult;

use strict;

use Bio::Search::Hit::BlastPullHit;

use base qw(Bio::Root::Root Bio::Search::Result::PullResultI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::Result::hmmpfam->new();
 Function: Builds a new Bio::SearchIO::Result::hmmpfam object 
 Returns : Bio::SearchIO::Result::hmmpfam
 Args    : -chunk  => [Bio::Root::IO, $start, $end] (required if no -parent)
           -parent => Bio::PullParserI object (required if no -chunk)
           -parameters => hash ref of search parameters (key => value), optional
           -statistics => hash ref of search statistics (key => value), optional

		   where the array ref provided to -chunk contains an IO object
           for a filehandle to something representing the raw data of the
           result, and $start and $end define the tell() position within the
           filehandle that the result data starts and ends (optional; defaults
           to start and end of the entire thing described by the filehandle)

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
	
	$self->_setup(@args);
	
	foreach my $field (qw( header hit_table hsp_table alignments next_model
	models query_length stats_and_params)) {
		$self->_fields->{$field} = undef;
	}
	
	$self->_dependencies( { ( query_name => 'header',
                              query_accession => 'header',
                              query_description => 'header',
			    query_length => 'header',
                              hit_table => 'header',
			    num_hits => 'hit_table',
			    no_hits_found => 'hit_table' ) } );
    
    return $self;
}

#
# PullParserI discovery methods so we can answer all ResultI questions
#

sub _discover_header {
	my $self = shift;
	$self->_chunk_seek(0);
	my $header = $self->_get_chunk_by_end("Value\n");
	if (!$header) {
	    $header = $self->_get_chunk_by_end("***** No hits found ******\n");
	    $self->{_no_hits_found} = 1;
	}
	$self->throw("Invalid header returned") unless $header;
	$self->{_after_header} = $self->_chunk_tell;
	
	($self->_fields->{query_name}) = $header =~ /^\s*(\S+)/;
	$self->_fields->{query_accession} = '';
	$self->_fields->{query_description} = '';
	
	if ($header =~ /^Length=(\d+)/m) {
		$self->_fields->{query_length} = $1;
	}
	elsif ($header =~ /\((\d+) letters\)/) {
		# older form?
		$self->_fields->{query_length} = $1;
		
		if ($header =~ /^\s*\(\d+ letters/) {
			# there wasn't a query sequence name
			$self->_fields->{query_name} = '';
		}
	}
	
	$self->_fields->{header} = 1;
}

sub _discover_hit_table {
	my $self = shift;
	$self->_chunk_seek($self->{_after_header});
	
	my $table = $self->_get_chunk_by_end("\n>");
	
	if (!$table && !$self->{_no_hits_found}) {
		# no alignments, presumably; hit table comprises the remainder of this
		# result
		while (my $line = $self->_get_chunk_by_nol(1)) {
			$table .= $line;
		}
	}
    
	$table ||= '';
	
	$self->{_after_hit_table} = $self->_chunk_tell;
	
	my $evalue_cutoff = $self->get_field('evalue_cutoff');
	undef $evalue_cutoff if $evalue_cutoff eq '[unset]';
	my $score_cutoff = $self->get_field('score_cutoff');
	undef $score_cutoff if $score_cutoff eq '[unset]';
	
	my @table;
	my $no_hit = 1;
	while ($table =~ /^(\S+)\s+(\S.*?)?\s+(\S+)\s+([\de]\S*)\s*\n/gm) {
		$no_hit = 0;
		my ($name, $desc, $score, $evalue) = ($1, $2, $3, $4);
		$desc ||= '';
        if ($evalue =~ /^e/) {
            $evalue = '1'.$evalue;
        }
		next if ($evalue_cutoff && $evalue > $evalue_cutoff);
		next if ($score_cutoff && $score < $score_cutoff);
		push(@table, [$name, $desc, $score, $evalue]);
	}
	$self->_fields->{hit_table} = \@table;
	$self->{_next_hit_index} = @table > 0 ? 0 : -1;
	
	$self->_fields->{no_hits_found} = $no_hit;
	$self->_fields->{num_hits} = @table;
}

sub _discover_next_hit {
	my $self = shift;
	my $hit_table = $self->get_field('hit_table');
	return if $self->{_next_hit_index} == -1;
	my $hit_row = ${$hit_table}[$self->{_next_hit_index}];
	
	$self->_chunk_seek($self->{_end_of_previous_hit} || $self->{_after_hit_table});
	
	my ($start, $end) = $self->_find_chunk_by_end("\n>");
	unless ($end) {
		$start = $self->{_end_of_previous_hit} || $self->{_after_hit_table};
		$end = $self->_chunk_true_end;
	}
	else {
		$end += $self->_chunk_true_start;
	}
	$start += $self->_chunk_true_start;
	
	$self->{_end_of_previous_hit} = $end - $self->_chunk_true_start;
	
	#*** needs to inherit piped_behaviour, and we need to deal with _sequential
	#    ourselves
	$self->_fields->{next_hit} = Bio::Search::Hit::BlastPullHit->new(-parent => $self,
								-chunk => [$self->chunk, $start, $end],
								-hit_data => $hit_row);
	
	$self->{_next_hit_index}++;
	
	if ($self->{_next_hit_index} > $#{$hit_table}) {
		$self->{_next_hit_index} = -1;
	}
}

sub _discover_stats_and_params {
	my $self = shift;
	$self->_chunk_seek(0);
	my ($start, $end) = $self->_find_chunk_by_end("\n  Database: ");
	$self->_chunk_seek($end);
	
	my $gapped = 0;
	while ($self->_get_chunk_by_nol(1)) {
		if (/Number of letters in database:\s+(\S+)/) {
			my $stat = $1;
			$stat =~ s/,//g;
			$self->add_statistic('dbletters', $stat);
		}
		elsif (/Number of sequences in database:\s+(\S+)/) {
			my $stat = $1;
			$stat =~ s/,//g;
			$self->add_statistic('dbentries', $stat);
		}
		elsif (/^Lambda/) {
			my $line = $self->_get_chunk_by_nol(1);
			$line =~ /\s+(\S+)\s+(\S+)\s+(\S+)/;
			$self->add_statistic($gapped ? 'lambda_gapped' : 'lambda', $1);
			$self->add_statistic($gapped ? 'kappa_gapped' : 'kappa', $2);
			$self->add_statistic($gapped ? 'entropy_gapped' : 'entropy', $3);
			$gapped = 1;
		}
		elsif (/^Matrix: (.+)\n/) {
			$self->add_parameter('matrix', $1);
		}
		elsif (/^Gap Penalties: Existence: (\d+), Extension: (\d+)/) {
			$self->add_parameter('gapopen', $1);
			$self->add_parameter('gapext', $2);
		}
		elsif (/^Number of Hits to DB: (\d+)/) {
			$self->add_statistic('Hits_to_DB', $1);
		}
		elsif (/^Number of extensions: (\d+)/) {
			$self->add_statistic('num_extensions', $1);
		}
		elsif (/^Number of successful extensions: (\d+)/) {
			$self->add_statistic('num_successful_extensions', $1);
		}
		elsif (/^Number of sequences better than (\S+):/) {
			$self->add_parameter('expect', $1);
		}
		elsif (/^[Ll]ength of query: (\d+)/) {
			$self->add_statistic('querylength', $1);
		}
        elsif (/^[Ee]ffective HSP length: (\d+)/) {
			$self->add_statistic('effective_hsplength', $1);
		}
		elsif (/^[Ee]ffective length of database: (\d+)/) {
			$self->add_statistic('effectivedblength', $1);
		}
		elsif (/^[Ee]ffective search space: (\d+)/) {
			$self->add_statistic('effectivespace', $1);
		}
		elsif (/^[Ee]ffective search space used: (\d+)/) {
			$self->add_statistic('effectivespaceused', $1);
		}
		elsif (/^([TAXS]\d?): (\d+)(?: \((\S+))?/) {
			$self->add_statistic($1, $2);
			$self->add_statistic($1.'_bits', $3) if $3;
		}
	}
	
	$self->_fields->{stats_and_params} = 1;
}

=head2 next_hit

 Title   : next_hit
 Usage   : while( $hit = $result->next_hit()) { ... }
 Function: Returns the next available Hit object, representing potential
           matches between the query and various entities from the database.
 Returns : a Bio::Search::Hit::HitI object or undef if there are no more.
 Args    : none

=cut

sub next_hit {
	my $self = shift;
    my $hit = $self->get_field('next_hit');
	undef $self->_fields->{next_hit};
	return $hit;
}

=head2 hits

 Title   : hits
 Usage   : my @hits = $result->hits
 Function: Returns the HitI objects contained within this Result
 Returns : Array of Bio::Search::Hit::HitI objects
 Args    : none

See Also: L<Bio::Search::Hit::HitI>

=cut

sub hits {
	my $self = shift;
	my $old = $self->{_next_hit_index} || 0;
	$self->rewind;
	my @hits;
	while (defined(my $hit = $self->next_hit)) {
		push(@hits, $hit);
	}
	$self->{_next_hit_index} = @hits > 0 ? $old : -1;
	return @hits;
}

=head2 sort_hits

 Title		: sort_hits
 Usage		: $result->sort_hits('<score')
 Function	: Sorts the hits so that they come out in the desired order when
              hits() or next_hit() is called.
 Returns	: n/a
 Args		: A coderef for the sort function. See the documentation on the Perl
              sort() function for guidelines on writing sort functions.
			  By default the sort order is ascending significance value (ie.
			  most significant hits first).
			  *** example

=cut

sub sort_hits {
    my ($self, $code_ref) = @_;
	$self->throw("Not implemented yet");
}

=head2 rewind

 Title   : rewind
 Usage   : $result->rewind;
 Function: Allow one to reset the Hit iterator to the beginning, so that
           next_hit() will subsequently return the first hit and so on.
 Returns : n/a
 Args    : none

=cut

sub rewind {
	my $self = shift;
	unless ($self->_fields->{hit_table}) {
		$self->get_field('hit_table');
	}
	$self->{_next_hit_index} = @{$self->_fields->{hit_table}} > 0 ? 0 : -1;
}

=head2 get_statistic

 Title   : get_statistic
 Usage   : my $gap_ext = $result->get_statistic('kappa')
 Function: Returns the value for a specific statistic available 
           from this result
 Returns : string
 Args    : name of statistic (string)

=cut

sub get_statistic {
	my $self = shift;
	$self->get_field('stats_and_params');
	return $self->SUPER::get_statistic(@_);
}

=head2 get_parameter

 Title   : get_parameter
 Usage   : my $gap_ext = $result->get_parameter('gapext')
 Function: Returns the value for a specific parameter used
           when running this result
 Returns : string
 Args    : name of parameter (string)

=cut

sub get_parameter {
	my $self = shift;
	$self->get_field('stats_and_params');
	return $self->SUPER::get_parameter(@_);
}

1;
