#
# BioPerl module for Bio::Search::Result::HmmpfamResult
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

Bio::Search::Result::HmmpfamResult - A parser and result object for hmmpfam
                                     results

=head1 SYNOPSIS

    # generally we use Bio::SearchIO to build these objects
    use Bio::SearchIO;
    my $in = Bio::SearchIO->new(-format => 'hmmer_pull',
							   -file   => 'result.hmmer');

    while (my $result = $in->next_result) {
		print $result->query_name, " ", $result->algorithm, " ", $result->num_hits(), " hits\n";
    }

=head1 DESCRIPTION

This object implements a parser for hmmpfam result output, a program in the HMMER
package.

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

package Bio::Search::Result::HmmpfamResult;

use strict;

use Bio::Search::Hit::HmmpfamHit;

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
	
	foreach my $field (qw( header hit_table hsp_table alignments next_model models query_length )) {
		$self->_fields->{$field} = undef;
	}
	
	$self->_dependencies( { ( query_name => 'header',
                              query_accession => 'header',
                              query_description => 'header',
                              hit_table => 'header',
							  num_hits => 'hit_table',
							  no_hits_found => 'hit_table',
                              hsp_table => 'hit_table',
                              next_alignment => 'hsp_table' ) } );
    
    return $self;
}

#
# PullParserI discovery methods so we can answer all ResultI questions
#

sub _discover_header {
	my $self = shift;
	$self->_chunk_seek(0);
	my $header = $self->_get_chunk_by_end("all domains):\n");
	$self->{_after_header} = $self->_chunk_tell;
	
	$header || $self->throw("Could not find hmmer header, is file really hmmer format?");
	
	($self->_fields->{query_name}) = $header =~ /^Query(?:\s+sequence)?:\s+(\S+)/m;
	($self->_fields->{query_accession}) = $header =~ /^Accession:\s+(\S+)/m;
	($self->_fields->{query_description}) = $header =~ /^Description:\s+(\S.+)/m;
	$self->_fields->{query_accession} ||= '';
	$self->_fields->{query_description} ||= '';
	
	$self->_fields->{header} = 1; # stop this method being called again
}

sub _discover_hit_table {
	my $self = shift;
	
	$self->_chunk_seek($self->{_after_header});
	my $table = $self->_get_chunk_by_end("for domains:\n");
	$self->{_after_hit_table} = $self->_chunk_tell;
	
	my $evalue_cutoff = $self->get_field('evalue_cutoff');
	undef $evalue_cutoff if $evalue_cutoff eq '[unset]';
	my $score_cutoff = $self->get_field('score_cutoff');
	undef $score_cutoff if $score_cutoff eq '[unset]';
	my $hsps_cutoff = $self->get_field('hsps_cutoff');
	undef $hsps_cutoff if $hsps_cutoff eq '[unset]';
	
	my @table;
	my $no_hit = 1;
	while ($table =~ /^(\S+)\s+(\S.+?)?\s+(\S+)\s+(\S+)\s+(\d+)\n/gm) {
		$no_hit = 0;
		my $evalue = abs($4); # consistency for tests under Windows
		next if ($evalue_cutoff && $evalue > $evalue_cutoff);
		next if ($score_cutoff && $3 < $score_cutoff);
		next if ($hsps_cutoff && $5 < $hsps_cutoff);
		push(@table, [$1, $2, $3, $evalue, $5]);
	}
	$self->_fields->{hit_table} = \@table;
	$self->{_next_hit_index} = @table > 0 ? 0 : -1;
	
	$self->_fields->{no_hits_found} = $no_hit;
	$self->_fields->{num_hits} = @table;
}

sub _discover_hsp_table {
	my $self = shift;
	
	$self->_chunk_seek($self->{_after_hit_table});
	my $table = $self->_get_chunk_by_end("top-scoring domains:\n");
	$table ||= $self->_get_chunk_by_end("//\n"); # A0 reports
	$self->{_after_hsp_table} = $self->_chunk_tell;
	
	my %table;
	# can't save this regex work for the hsp object because the hit object needs
	# its length, so may as well just do all the work here
	while ($table =~ /^(\S+)\s+(\d+)\/\d+\s+(\d+)\s+(\d+)\s+\S\S\s+(\d+)\s+(\d+)\s+\S(\S)\s+(\S+)\s+(\S+)/gm) {
		# rank query_start query_end hit_start hit_end score evalue
		my $evalue = abs($9); # consistency for tests under Windows
		push(@{$table{$1}->{hsp_data}}, [$2, $3, $4, $5, $6, $8, $evalue]);
		if ($7 eq ']') {
			$table{$1}->{hit_length} = $6;
		}
	}
	$self->_fields->{hsp_table} = \%table;
}

sub _discover_alignments {
	my $self = shift;
	$self->_fields->{alignments} = { };
}

sub _next_alignment {
	my $self = shift;;
	return if $self->{_no_more_alignments};
	
	my $aligns = $self->_fields->{alignments};
	
	unless (defined $self->{_after_previous_alignment}) {
		$self->_chunk_seek($self->{_after_hsp_table});
		my $chunk = $self->_get_chunk_by_end(": domain");
		unless ($chunk) {
			$self->{_no_more_alignments} = 1;
			return;
		}
		
		$self->{_after_previous_alignment} = $self->_chunk_tell;
		$self->{_next_alignment_start_text} = $chunk;
		return $self->_next_alignment;
	}
	
	$self->_chunk_seek($self->{_after_previous_alignment});
	my $chunk = $self->_get_chunk_by_end(": domain");
	unless ($chunk) {
		$self->_chunk_seek($self->{_after_previous_alignment});
		$chunk = $self->_get_chunk_by_end("//");
		
		unless ($chunk) {
			$self->{_no_more_alignments} = 1;
			return;
		}
	}
	
	$self->{_after_previous_alignment} = $self->_chunk_tell;
	
	if (defined $self->{_next_alignment_start_text}) {
		$chunk = $self->{_next_alignment_start_text}.$chunk;
	}
	
	$chunk =~ s/(\S+: domain)$//;
	$self->{_next_alignment_start_text} = $1;
	
	my ($name, $domain) = $chunk =~ /^(\S+): domain (\d+)/;
	$aligns->{$name.'~~~~'.$domain} = $chunk;
	return 1;
}

sub _discover_next_hit {
	my $self = shift;
	my @hit_table = @{$self->get_field('hit_table')};
	return if $self->{_next_hit_index} == -1;
	
	#[name description score significance num_hsps rank]
	my @hit_data = (@{$hit_table[$self->{_next_hit_index}++]}, $self->{_next_hit_index});
	
	$self->_fields->{next_hit} = Bio::Search::Hit::HmmpfamHit->new(-parent => $self,
																  -hit_data => \@hit_data);
	
	if ($self->{_next_hit_index} > $#hit_table) {
		$self->{_next_hit_index} = -1;
	}
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

=head2 next_model

 Title   : next_model
 Usage   : my $domain = $result->next_model
 Function: Returns the next domain - this is an alias for next_hit()
 Returns : L<Bio::Search::Hit::HitI> object
 Args    : none

=cut

*next_model = \&next_hit;

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

=head2 models

 Title   : models
 Usage   : my @domains = $result->models;
 Function: Returns the list of HMM models seen - this is an alias for hits()
 Returns : Array of L<Bio::Search::Hit::HitI> objects
 Args    : none

=cut

*models = \&hits;

=head2 sort_hits

 Title		: sort_hits
 Usage		: $result->sort_hits('<score')
 Function	: Sorts the hits so that they come out in the desired order when
              hits() or next_hit() is called.
 Returns	: n/a
 Args		: A coderef for the sort function. See the documentation on the Perl
              sort() function for guidelines on writing sort functions.
			  You will be sorting array references, not HitI objects. The
			  references contain name as element 0, description as element 1,
			  score as element 2, significance as element 3 and number of hsps
			  as element 4.
			  By default the sort order is ascending significance value (ie.
			  most significant hits first).
 Note		: To access the special variables $a and $b used by the Perl sort()
              function the user function must access
			  Bio::Search::Result::HmmpfamResult namespace. 
              For example, use : 
              $result->sort_hits(
				sub{$Bio::Search::Result::HmmpfamResult::a->[2]
				                         <=> 
					$Bio::Search::Result::HmmpfamResult::b->[2]});
              NOT $result->sort_hits($a->[2] <=> $b->[2]);

=cut

sub sort_hits {
    my ($self, $code_ref) = @_;
	$code_ref ||= sub { $a->[3] <=> $b->[3] };
	
	# avoid creating hit objects just to sort, hence force user to sort on
	# the array references in hit table
	my $table_ref = $self->get_field('hit_table');
	@{$table_ref} > 1 || return;
	
	my @sorted = sort $code_ref @{$table_ref};
	@sorted == @{$table_ref} || $self->throw("Your sort routine failed to give back all hits!");
	$self->_fields->{hit_table} = \@sorted;
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

1;
