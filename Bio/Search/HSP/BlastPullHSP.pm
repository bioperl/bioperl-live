#
# BioPerl module for Bio::Search::HSP::BlastPullHSP
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

Bio::Search::HSP::BlastPullHSP - A parser and HSP object for BlastN hsps

=head1 SYNOPSIS

    # generally we use Bio::SearchIO to build these objects
    use Bio::SearchIO;
    my $in = Bio::SearchIO->new(-format => 'hmmer_pull',
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

This object implements a parser for BlastN hsp output.

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

package Bio::Search::HSP::BlastPullHSP;

use strict;
use base qw(Bio::Search::HSP::PullHSPI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::HSP::BlastNHSP->new();
 Function: Builds a new Bio::Search::HSP::BlastNHSP object.
 Returns : Bio::Search::HSP::BlastNHSP
 Args    : -chunk  => [Bio::Root::IO, $start, $end] (required if no -parent)
           -parent => Bio::PullParserI object (required if no -chunk)

           where the array ref provided to -chunk contains an IO object
           for a filehandle to something representing the raw data of the
           hsp, and $start and $end define the tell() position within the
           filehandle that the hsp data starts and ends (optional; defaults
           to start and end of the entire thing described by the filehandle)

=cut

sub new {
    my ($class, @args) = @_;
	my $self = $class->SUPER::new(@args);
	
	$self->_setup(@args);
	
	my $fields = $self->_fields;
	foreach my $field (qw( header alignment query_strand hit_strand )) {
		$fields->{$field} = undef;
	}
	
	$self->_dependencies( { ( score => 'header',
                              bits => 'header',
                              evalue => 'header',
                              total_gaps => 'header',
                              query_strand => 'header',
                              hit_strand => 'header',
                              alignment => 'header',
                              query_string => 'alignment',
                              hit_string => 'alignment',
                              homology_string => 'alignment',
                              query_start => 'alignment',
                              query_end => 'alignment',
                              hit_start => 'alignment',
                              hit_end => 'alignment',
                              hit_identical_inds => 'seq_inds',
							  hit_conserved_inds => 'seq_inds',
							  hit_nomatch_inds => 'seq_inds',
                              hit_gap_inds => 'seq_inds',
                              query_identical_inds => 'seq_inds',
							  query_conserved_inds => 'seq_inds',
							  query_nomatch_inds => 'seq_inds',
							  query_gap_inds => 'seq_inds' ) } );
	
    return $self;
}

#
# PullParserI discovery methods so we can answer all HitI questions
#

sub _discover_header {
    my $self = shift;
    $self->_chunk_seek(0);
	my $header = $self->_get_chunk_by_end("\nQuery");
	$self->{_after_header} = $self->_chunk_tell;
	
	($self->_fields->{bits}, $self->_fields->{score}, $self->_fields->{evalue},
     $self->_fields->{total_gaps}, $self->_fields->{query_strand}, $self->_fields->{hit_strand})
     = $header =~ /^\s*(\S+) bits \((\d+)\),\s+Expect = (\S+)(?:\s+.+Gaps = (\d+))?(?:.+Strand\s*=\s*(\w+)\s*\/\s*(\w+))?/sm;
    
    if ($self->_fields->{query_strand}) {
        # protein blasts don't have strand
        for my $strand_type ('query_strand', 'hit_strand') {
            $self->_fields->{$strand_type} = $self->_fields->{$strand_type} eq 'Plus' ? 1 : -1;
        }
    }
    else {
        $self->_fields->{query_strand} = 0;
        $self->_fields->{hit_strand} = 0;
    }
	
	if ($self->_fields->{evalue} =~ /^e/) {
		$self->_fields->{evalue} = '1'.$self->_fields->{evalue};
	}
    
    # query_gaps isn't always given
    $self->_fields->{total_gaps} = '[unset]' unless $self->_fields->{total_gaps};
    
	$self->_fields->{header} = 1;
}

sub _discover_alignment {
    my $self = shift;
    $self->_chunk_seek($self->{_after_header});
    
    # work out various basic fields for the hsp
    # (quicker to do this all at once instead of each method working on
    # alignment itself)
    my ($query_string, $hit_string, $homology_string, $q_start, $h_start, $q_end, $h_end);
    while (my $strip = $self->_get_chunk_by_end("\nQuery") || $self->_get_chunk_by_nol(4)) {
        $strip =~ /\s+(\d+)\s+(\S+)\s+(\d+)\s+(\S.+)\nSbjct:?\s+(\d+)\s+(\S+)\s+(\d+)/gm || last;
        my $q1 = $1;
        $query_string .= $2;
        my $q2 = $3;
        my $hom = $4;
        my $h1 = $5;
        $hit_string .= $6;
        my $h2 = $7;
        
        $hom = ' 'x(length($6) - length($hom)).$hom;
        $homology_string .= $hom;
        
        for my $q ($q1, $q2) {
            if (! defined $q_start || $q < $q_start) {
                $q_start = $q;
            }
            if (! defined $q_end || $q > $q_end) {
                $q_end = $q;
            }
        }
        for my $h ($h1, $h2) {
            if (! defined $h_start || $h < $h_start) {
                $h_start = $h;
            }
            if (! defined $h_end || $h > $h_end) {
                $h_end = $h;
            }
        }
    }
    
    $self->_fields->{query_string} = $query_string;
    $self->_fields->{hit_string} = $hit_string;
    $self->_fields->{homology_string} = $homology_string;
    
    $self->_fields->{query_start} = $q_start;
    $self->_fields->{query_end} = $q_end;
    $self->_fields->{hit_start} = $h_start;
    $self->_fields->{hit_end} = $h_end;
    
    ($self->{_query_gaps}) = $query_string =~ tr/-//;
    ($self->{_hit_gaps}) = $hit_string =~ tr/-//;
    ($self->{_total_gaps}) = $self->{_query_gaps} + $self->{_hit_gaps};
    
    $self->_fields->{alignment} = 1; # stop this method being called again
}

# seq_inds related methods, all just need seq_inds field to have been gotten
sub _discover_seq_inds {
    my $self = shift;
    my ($seqString, $qseq, $sseq) = ( $self->get_field('homology_string'),
                                      $self->get_field('query_string'),
                                      $self->get_field('hit_string') );
    
    # (code largely lifted from GenericHSP)
    
    # Using hashes to avoid saving duplicate residue numbers.
    my %identicalList_query = ();
    my %identicalList_sbjct = ();
    my %conservedList_query = ();
    my %conservedList_sbjct = ();
    my @gapList_query = ();
    my @gapList_sbjct = ();
    my %nomatchList_query = ();
    my %nomatchList_sbjct = ();
    
    my $resCount_query = $self->get_field('query_end');
    my $resCount_sbjct = $self->get_field('hit_end');
    
    my ($mchar, $schar, $qchar);
    while ($mchar = chop($seqString) ) {
        ($qchar, $schar) = (chop($qseq), chop($sseq));
        
        if ($mchar eq '+' || $mchar eq '.' || $mchar eq ':') { 
            $conservedList_query{ $resCount_query } = 1; 
            $conservedList_sbjct{ $resCount_sbjct } = 1;
        }
        elsif ($mchar eq ' ') { 
            $nomatchList_query{ $resCount_query } = 1;
            $nomatchList_sbjct{ $resCount_sbjct } = 1;
        }
        else { 
            $identicalList_query{ $resCount_query } = 1; 
            $identicalList_sbjct{ $resCount_sbjct } = 1;
        }
        
        if ($qchar eq '-') {
            push(@gapList_query, $resCount_query);
        }
        else { 	    
            $resCount_query -= 1;
        }
        if ($schar eq '-') {
            push(@gapList_sbjct, $resCount_sbjct);
        }
        else { 	    
            $resCount_sbjct -= 1;
        }
    }
    
    my $fields = $self->_fields;
    $fields->{hit_identical_inds} = [ sort { $a <=> $b } keys %identicalList_sbjct ];
    $fields->{hit_conserved_inds} = [ sort { $a <=> $b } keys %conservedList_sbjct ];
    $fields->{hit_nomatch_inds} = [ sort { $a <=> $b } keys %nomatchList_sbjct ];
    $fields->{hit_gap_inds} = [ reverse @gapList_sbjct ];
    $fields->{query_identical_inds} = [ sort { $a <=> $b } keys %identicalList_query ];
    $fields->{query_conserved_inds} = [ sort { $a <=> $b } keys %conservedList_query ];
    $fields->{query_nomatch_inds} = [ sort { $a <=> $b } keys %nomatchList_query ];
    $fields->{query_gap_inds} = [ reverse @gapList_query ];
    
    $fields->{seq_inds} = 1;
}

=head2 query

 Title   : query
 Usage   : my $query = $hsp->query
 Function: Returns a SeqFeature representing the query in the HSP
 Returns : L<Bio::SeqFeature::Similarity>
 Args    : none

=cut

sub query {
    my $self = shift;
    unless ($self->{_created_query}) {
        $self->SUPER::query( new  Bio::SeqFeature::Similarity
                  ('-primary'  => $self->primary_tag,
                   '-start'    => $self->get_field('query_start'),
                   '-end'      => $self->get_field('query_end'),
                   '-expect'   => $self->get_field('evalue'),
                   '-score'    => $self->get_field('score'),
                   '-strand'   => $self->get_field('query_strand'),
                   '-seq_id'   => $self->get_field('query_name'),
                   '-seqlength'=> $self->get_field('query_length'),
                   '-source'   => $self->get_field('algorithm'),
                   '-seqdesc'  => $self->get_field('query_description'),
                   '-frame'    => 0 # not known?
                   ) );
		$self->{_created_query} = 1;
    }
    return $self->SUPER::query(@_);
}

=head2 hit

 Title   : hit
 Usage   : my $hit = $hsp->hit
 Function: Returns a SeqFeature representing the hit in the HSP
 Returns : L<Bio::SeqFeature::Similarity>
 Args    : [optional] new value to set

=cut

sub hit {
    my $self = shift;
    unless ($self->{_created_hit}) {
        $self->SUPER::hit( new  Bio::SeqFeature::Similarity
                  ('-primary'  => $self->primary_tag,
                   '-start'    => $self->get_field('hit_start'),
                   '-end'      => $self->get_field('hit_end'),
                   '-expect'   => $self->get_field('evalue'),
                   '-score'    => $self->get_field('score'),
                   '-strand'   => $self->get_field('hit_strand'),
                   '-seq_id'   => $self->get_field('name'),
                   '-seqlength'=> $self->get_field('length'),
                   '-source'   => $self->get_field('algorithm'),
                   '-seqdesc'  => $self->get_field('description'),
                   '-frame'    => 0 # not known?
                   ) );
		$self->{_created_hit} = 1;
    }
    return $self->SUPER::hit(@_);
}

=head2 gaps

 Title    : gaps
 Usage    : my $gaps = $hsp->gaps( ['query'|'hit'|'total'] );
 Function : Get the number of gap characters in the query, hit, or total alignment.
 Returns  : Integer, number of gap characters or 0 if none
 Args     : 'query' = num conserved / length of query seq (without gaps)
            'hit'   = num conserved / length of hit seq (without gaps)
            'total' = num conserved / length of alignment (with gaps)
            default = 'total' 

=cut

sub gaps {
    my ($self, $type) = @_;
    
    $type = lc $type if defined $type;
    $type = 'total' if (! defined $type || $type eq 'hsp' || $type !~ /query|hit|subject|sbjct|total/); 
    $type = 'hit' if $type =~ /sbjct|subject/;
    
    if ($type eq 'total') {
        my $answer = $self->get_field('total_gaps');
        return $answer unless $answer eq '[unset]';
    }
    
    $self->get_field('alignment'); # make sure gaps have been calculated
    
    return $self->{'_'.$type.'_gaps'};
}

=head2 strand

 Title   : strand
 Usage   : $hsp->strand('query')
 Function: Retrieves the strand for the HSP component requested
 Returns : +1 or -1 (0 if unknown)
 Args    : 'hit' or 'subject' or 'sbjct' to retrieve the strand of the subject
           'query' to retrieve the query strand (default)
           'list' or 'array' to retreive both query and hit together

=cut

sub strand {
    my $self = shift;
    my $val = shift;
    $val = 'query' unless defined $val;
    $val =~ s/^\s+//;

    if ($val =~ /^q/i) {
        return $self->get_field('query_strand');
    }
    elsif ($val =~ /^hi|^s/i) {
        return $self->get_field('hit_strand');
    }
    elsif ($val =~ /^list|array/i) {
        return ($self->get_field('query_strand'), $self->get_field('hit_strand'));
    }
    else { 
        $self->warn("unrecognized component '$val' requested\n");
    }
    return 0;
}

=head2 start

 Title   : start
 Usage   : $hsp->start('query')
 Function: Retrieves the start for the HSP component requested
 Returns : integer, or list of two integers (query start and subject start) in
           list context
 Args    : 'hit' or 'subject' or 'sbjct' to retrieve the start of the subject
           'query' to retrieve the query start (default)

=cut

sub start {
    my $self = shift;
    my $val = shift;
    $val = (wantarray ? 'list' : 'query') unless defined $val;
    $val =~ s/^\s+//;
    
    if ($val =~ /^q/i) { 
        return $self->get_field('query_start');
    }
    elsif ($val =~ /^(hi|s)/i) {
        return $self->get_field('hit_start');
    }
    elsif ($val =~ /^list|array/i) {
        return ($self->get_field('query_start'), $self->get_field('hit_start') );
    }
    else { 
        $self->warn("unrecognized component '$val' requested\n");
    }
    return 0;
}

=head2 end

 Title   : end
 Usage   : $hsp->end('query')
 Function: Retrieves the end for the HSP component requested
 Returns : integer, or list of two integers (query end and subject end) in
           list context
 Args    : 'hit' or 'subject' or 'sbjct' to retrieve the end of the subject
           'query' to retrieve the query end (default)

=cut

sub end {
    my $self = shift;
    my $val = shift;
    $val = (wantarray ? 'list' : 'query') unless defined $val;
    $val =~ s/^\s+//;
    
    if ($val =~ /^q/i) { 
        return $self->get_field('query_end');
    }
    elsif ($val =~ /^(hi|s)/i) {
        return $self->get_field('hit_end');
    }
    elsif ($val =~ /^list|array/i) {
        return ($self->get_field('query_end'), $self->get_field('hit_end'));
    }
    else {
        $self->warn("unrecognized end component '$val' requested\n");
    }
    return 0;
}

=head2 pvalue

 Title   : pvalue
 Usage   : my $pvalue = $hsp->pvalue();
 Function: Returns the P-value for this HSP
 Returns : undef (Hmmpfam reports do not have p-values)
 Args    : none

=cut

sub pvalue { }

1;
