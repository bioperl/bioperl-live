#
# BioPerl module for Bio::Search::Hit::HmmpfamHit
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

Bio::Search::Hit::HmmpfamHit - A parser and hit object for hmmpfam hits

=head1 SYNOPSIS

    # generally we use Bio::SearchIO to build these objects
    use Bio::SearchIO;
    my $in = Bio::SearchIO->new(-format => 'hmmer_pull',
							   -file   => 'result.hmmer');

    while (my $result = $in->next_result) {
		while (my $hit = $result->next_hit) {
			print $hit->name, "\n";
			print $hit->score, "\n";
			print $hit->significance, "\n";

			while (my $hsp = $hit->next_hsp) {
				# process HSPI objects
			}
		}
    }

=head1 DESCRIPTION

This object implements a parser for hmmpfam hit output, a program in the HMMER
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

package Bio::Search::Hit::HmmpfamHit;

use strict;

use Bio::Search::HSP::HmmpfamHSP;

use base qw(Bio::Root::Root Bio::Search::Hit::PullHitI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::Hit::HmmpfamHit->new();
 Function: Builds a new Bio::Search::Hit::HmmpfamHit object.
 Returns : Bio::Search::Hit::HmmpfamHit
 Args    : -chunk    => [Bio::Root::IO, $start, $end] (required if no -parent)
           -parent   => Bio::PullParserI object (required if no -chunk)
           -hit_data => array ref with [name description score significance
		                                num_hsps rank]

           where the array ref provided to -chunk contains an IO object
           for a filehandle to something representing the raw data of the
           hit, and $start and $end define the tell() position within the
           filehandle that the hit data starts and ends (optional; defaults
           to start and end of the entire thing described by the filehandle)

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
	
	$self->_setup(@args);
	
	my $fields = $self->_fields;
	foreach my $field (qw( next_domain domains hsp_data )) {
		$fields->{$field} = undef;
	}
	
	my $hit_data = $self->_raw_hit_data;
	if ($hit_data && ref($hit_data) eq 'ARRAY') {
		foreach my $field (qw(name description score significance num_hsps rank)) {
			$fields->{$field} = shift(@{$hit_data});
		}
	}
	$fields->{hit_start} = 1;
	
	delete $self->_fields->{accession};
	
	$self->_dependencies( { ( length => 'hsp_data' ) } );
	
    return $self;
}

#
# PullParserI discovery methods so we can answer all HitI questions
#

sub _discover_description {
	# this should be set when this object is created, but if it was undef as is
	# possible, this _discover method will be called: just return and keep the
	# return value undef
	return;
}

sub _discover_hsp_data {
	my $self = shift;
	my $hsp_table = $self->get_field('hsp_table');
	my $hsp_data = $hsp_table->{$self->get_field('name')} || undef;
	if ($hsp_data) {
		if (defined $hsp_data->{hit_length}) {
			$self->_fields->{length} = $hsp_data->{hit_length};
		}
		
		# rank query_start query_end hit_start hit_end score evalue
		$self->_fields->{hsp_data} = $hsp_data->{hsp_data};
	}
}

sub _discover_query_start {
	my $self = shift;
	my $hsp_data = $self->get_field('hsp_data') || return;
	
	my ($this_hsp) = sort { $a->[1] <=> $b->[1] } @{$hsp_data};
	$self->_fields->{query_start} = $this_hsp->[1];
}

sub _discover_query_end {
	my $self = shift;
	my $hsp_data = $self->get_field('hsp_data') || return;
	
	my ($this_hsp) = sort { $b->[2] <=> $a->[2] } @{$hsp_data};
	$self->_fields->{query_end} = $this_hsp->[2];
}

sub _discover_hit_start {
	my $self = shift;
	my $hsp_data = $self->get_field('hsp_data') || return;
	
	my ($this_hsp) = sort { $a->[3] <=> $b->[3] } @{$hsp_data};
	$self->_fields->{hit_start} = $this_hsp->[3];
}

sub _discover_hit_end {
	my $self = shift;
	my $hsp_data = $self->get_field('hsp_data') || return;
	
	my ($this_hsp) = sort { $b->[4] <=> $a->[4] } @{$hsp_data};
	$self->_fields->{hit_end} = $this_hsp->[4];
}

sub _discover_next_hsp {
	my $self = shift;
	my $hsp_data = $self->get_field('hsp_data') || return;
	unless (defined $self->{_next_hsp_index}) {
		$self->{_next_hsp_index} = 0;
	}
	return if $self->{_next_hsp_index} == -1;
	
	$self->_fields->{next_hsp} = Bio::Search::HSP::HmmpfamHSP->new(-parent => $self,
																  -hsp_data => $hsp_data->[$self->{_next_hsp_index}++]);
	
	if ($self->{_next_hsp_index} > $#{$hsp_data}) {
		$self->{_next_hsp_index} = -1;
	}
}

=head2 next_hsp

 Title    : next_hsp
 Usage    : while( $hsp = $obj->next_hsp()) { ... }
 Function : Returns the next available High Scoring Pair
 Example  : 
 Returns  : L<Bio::Search::HSP::HSPI> object or null if finished
 Args     : none

=cut

sub next_hsp {
    my $self = shift;
    my $hsp = $self->get_field('next_hsp');
	undef $self->_fields->{next_hsp};
	return $hsp;
}

=head2 next_domain

 Title   : next_domain 
 Usage   : my $domain = $hit->next_domain();
 Function: An alias for L<next_hsp()>, this will return the next HSP
 Returns : L<Bio::Search::HSP::HSPI> object
 Args    : none

=cut

*next_domain = \&next_hsp;

=head2 hsps

 Usage     : $hit_object->hsps();
 Purpose   : Get a list containing all HSP objects.
 Example   : @hsps = $hit_object->hsps();
 Returns   : list of L<Bio::Search::HSP::BlastHSP> objects.
 Argument  : none

=cut

sub hsps {
    my $self = shift;
	my $old = $self->{_next_hsp_index} || 0;
	$self->rewind;
	my @hsps;
	while (defined(my $hsp = $self->next_hsp)) {
		push(@hsps, $hsp);
	}
	$self->{_next_hsp_index} =  @hsps > 0 ? $old : -1;
	return @hsps;
}

=head2 domains

 Title   : domains
 Usage   : my @domains = $hit->domains();
 Function: An alias for L<hsps()>, this will return the full list of hsps
 Returns : array of L<Bio::Search::HSP::HSPI> objects
 Args    : none

=cut

*domains = \&hsps;

=head2 hsp

 Usage     : $hit_object->hsp( [string] );
 Purpose   : Get a single HSPI object for the present HitI object.
 Example   : $hspObj  = $hit_object->hsp;  # same as 'best'
           : $hspObj  = $hit_object->hsp('best');
           : $hspObj  = $hit_object->hsp('worst');
 Returns   : Object reference for a L<Bio::Search::HSP::HSPI> object.
 Argument  : String (or no argument).
           :   No argument (default) = highest scoring HSP (same as 'best').
           :   'best'  = highest scoring HSP.
           :   'worst' = lowest scoring HSP.
 Throws    : Exception if an unrecognized argument is used.

See Also   : L<hsps()|hsps>, L<num_hsps>()

=cut

sub hsp {
    my ($self, $type) = @_;
	$type ||= 'best';
	my $hsp_data = $self->get_field('hsp_data') || return;
	
	my $sort;
	if ($type eq 'best') {
		$sort = sub { $a->[6] <=> $b->[6] };
	}
	elsif ($type eq 'worst') {
		$sort = sub { $b->[6] <=> $a->[6] };
	}
	else {
		$self->throw("Unknown arg '$type' given to hsp()");
	}
	
	my ($this_hsp) = sort $sort @{$hsp_data};
	return Bio::Search::HSP::HmmpfamHSP->new(-parent => $self, -hsp_data => $this_hsp);
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
	my $hsp_data = $self->get_field('hsp_data') || return;
	$self->{_next_hsp_index} = @{$hsp_data} > 0 ? 0 : -1;
}

# have p() a synonym of significance()
sub p {
	return shift->significance;
}

=head2 strand

 Usage     : $sbjct->strand( [seq_type] );
 Purpose   : Gets the strand(s) for the query, sbjct, or both sequences.
           : For hmmpfam, the answers are always 1 (forward strand).
 Example   : $qstrand = $sbjct->strand('query');
           : $sstrand = $sbjct->strand('hit');
           : ($qstrand, $sstrand) = $sbjct->strand();
 Returns   : scalar context: integer '1'
           : array context without args: list of two strings (1, 1)
           : Array context can be "induced" by providing an argument of 'list'
		   : or 'array'.
 Argument  : In scalar context: seq_type = 'query' or 'hit' or 'sbjct' (default
           : = 'query') ('sbjct' is synonymous with 'hit')

=cut

sub strand {
    my ($self, $type) = @_;
	$type ||= (wantarray ? 'list' : 'query');
    $type = lc($type);
	if ($type eq 'list' || $type eq 'array') {
		return (1, 1);
	}
	return 1;
}

=head2 frac_aligned_query

 Usage     : $hit_object->frac_aligned_query();
 Purpose   : Get the fraction of the query sequence which has been aligned
           : across all HSPs (not including intervals between non-overlapping
           : HSPs).
 Example   : $frac_alnq = $hit_object->frac_aligned_query();
 Returns   : undef (the length of query sequences is unknown in Hmmpfam reports)
 Argument  : none

=cut

# noop
sub frac_aligned_query { }

1;
