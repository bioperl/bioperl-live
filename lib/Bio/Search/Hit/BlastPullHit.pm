#
# BioPerl module for Bio::Search::Hit::BlastPullHit
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

Bio::Search::Hit::BlastPullHit - A parser and hit object for BLASTN hits

=head1 SYNOPSIS

    # generally we use Bio::SearchIO to build these objects
    use Bio::SearchIO;
    my $in = Bio::SearchIO->new(-format => 'blast_pull',
							   -file   => 'result.blast');

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

This object implements a parser for BLASTN hit output.

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

package Bio::Search::Hit::BlastPullHit;

use strict;

use Bio::Search::HSP::BlastPullHSP;

use base qw(Bio::Root::Root Bio::Search::Hit::PullHitI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::Hit::BlastNHit->new();
 Function: Builds a new Bio::Search::Hit::BlastNHit object.
 Returns : Bio::Search::Hit::BlastNHit
 Args    : -chunk    => [Bio::Root::IO, $start, $end] (required if no -parent)
           -parent   => Bio::PullParserI object (required if no -chunk)
           -hit_data => array ref with [name description score significance]

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
	foreach my $field (qw( header start_end )) {
		$fields->{$field} = undef;
	}
	
    my $hit_data = $self->_raw_hit_data;
	if ($hit_data && ref($hit_data) eq 'ARRAY') {
		foreach my $field (qw(name description score significance)) {
			$fields->{$field} = shift(@{$hit_data});
		}
	}
    
	$self->_dependencies( { ( name => 'header',
							  length => 'header',
							  description => 'header',
							  accession => 'header',
                              next_hsp => 'header',
                              query_start => 'start_end',
                              query_end => 'start_end',
                              hit_start => 'start_end',
                              hit_end => 'start_end' ) } );
    
    return $self;
}

#
# PullParserI discovery methods so we can answer all HitI questions
#

sub _discover_header {
	my $self = shift;
	$self->_chunk_seek(0);
	my $header = $self->_get_chunk_by_end("\n Score = ");
    
    unless ($header) {
        # no alignment or other data; all information was in the hit table of
        # the result
        $self->_calculate_accession_from_name;
        
        $self->_fields->{header} = 1;
        return;
    }
    
	$self->{_after_header} = $self->_chunk_tell;
	
    ($self->_fields->{name}, $self->_fields->{description}, $self->_fields->{length}) = $header =~ /^(\S+)\s+(\S.+?)?\s+Length\s*=\s*(\d+)/sm;
    if ($self->_fields->{description}) {
        $self->_fields->{description} =~ s/\n//g;
    }
    else {
        $self->_fields->{description} = '';
    }
	
    $self->_calculate_accession_from_name;
	
	$self->_fields->{header} = 1;
}

sub _calculate_accession_from_name {
    my $self = shift;
    my $name = $self->get_field('name');
    if ($name =~ /.+?\|.+?\|.+?\|(\w+)/) {
        $self->_fields->{accession} = $1;
    }
    elsif ($self->_fields->{name} =~ /.+?\|(\w+)?\./) {
        # old form?
        $self->_fields->{accession} = $1;
    }
    else {
        $self->_fields->{accession} = $name;
    }
}

sub _discover_start_end {
    my $self = shift;
    
    my ($q_start, $q_end, $h_start, $h_end);
    foreach my $hsp ($self->hsps) {
        my ($this_q_start, $this_h_start) = $hsp->start;
        my ($this_q_end, $this_h_end) = $hsp->end;
        
        if (! defined $q_start || $this_q_start < $q_start) {
            $q_start = $this_q_start;
        }
        if (! defined $h_start || $this_h_start < $h_start) {
            $h_start = $this_h_start;
        }
        
        if (! defined $q_end || $this_q_end > $q_end) {
            $q_end = $this_q_end;
        }
        if (! defined $h_end || $this_h_end > $h_end) {
            $h_end = $this_h_end;
        }
    }
    
    $self->_fields->{query_start} = $q_start;
    $self->_fields->{query_end} = $q_end;
    $self->_fields->{hit_start} = $h_start;
    $self->_fields->{hit_end} = $h_end;
}

sub _discover_next_hsp {
	my $self = shift;
    my $pos = $self->{_end_of_previous_hsp} || $self->{_after_header};
    return unless $pos;
    $self->_chunk_seek($pos);
    
    my ($start, $end) = $self->_find_chunk_by_end("\n Score = ");
    if ((defined $end && ($end + $self->_chunk_true_start) > $self->_chunk_true_end) || ! $end) {
		$start = $self->{_end_of_previous_hsp} || $self->{_after_header};
		$end = $self->_chunk_true_end;
	}
	else {
		$end += $self->_chunk_true_start;
	}
	$start += $self->_chunk_true_start;
    
    return if $start >= $self->_chunk_true_end;
    
    $self->{_end_of_previous_hsp} = $end - $self->_chunk_true_start;
    
    #*** needs to inherit piped_behaviour, and we need to deal with _sequential
	#    ourselves
	$self->_fields->{next_hsp} = Bio::Search::HSP::BlastPullHSP->new(-parent => $self,
                                                                    -chunk => [$self->chunk, $start, $end]);
}

sub _discover_num_hsps {
    my $self = shift;
    $self->_fields->{num_hsps} = $self->hsps;
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

=head2 hsps

 Usage     : $hit_object->hsps();
 Purpose   : Get a list containing all HSP objects.
 Example   : @hsps = $hit_object->hsps();
 Returns   : list of L<Bio::Search::HSP::BlastHSP> objects.
 Argument  : none

=cut

sub hsps {
    my $self = shift;
	my $old = $self->{_end_of_previous_hsp};
	$self->rewind;
	my @hsps;
	while (defined(my $hsp = $self->next_hsp)) {
		push(@hsps, $hsp);
	}
	$self->{_end_of_previous_hsp} = $old;
	return @hsps;
}

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
	$self->throw_not_implemented;
}

=head2 rewind

 Title   : rewind
 Usage   : $result->rewind;
 Function: Allow one to reset the HSP iterator to the beginning, so that
           next_hsp() will subsequently return the first hsp and so on.
 Returns : n/a
 Args    : none

=cut

sub rewind {
	my $self = shift;
	delete $self->{_end_of_previous_hsp};
}

# have p() a synonym of significance()
sub p {
	return shift->significance;
}

1;
