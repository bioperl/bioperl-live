#
# BioPerl module for Bio::Search::HSP::HmmpfamHSP
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

Bio::Search::HSP::HmmpfamHSP - A parser and HSP object for hmmpfam hsps

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

This object implements a parser for hmmpfam hsp output, a program in the HMMER
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

package Bio::Search::HSP::HmmpfamHSP;

use strict;
use base qw(Bio::Search::HSP::PullHSPI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Search::HSP::HmmpfamHSP->new();
 Function: Builds a new Bio::Search::HSP::HmmpfamHSP object.
 Returns : Bio::Search::HSP::HmmpfamHSP
 Args    : -chunk  => [Bio::Root::IO, $start, $end] (required if no -parent)
           -parent => Bio::PullParserI object (required if no -chunk)
           -hsp_data => array ref with [rank query_start query_end hit_start
										hit_end score evalue]

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
	foreach my $field (qw( alignment )) {
		$fields->{$field} = undef;
	}
	
	my $hsp_data = $self->_raw_hsp_data;
	if ($hsp_data && ref($hsp_data) eq 'ARRAY') {
		my @hsp_data = @{$hsp_data}; # don't alter the reference
		foreach my $field (qw(rank query_start query_end hit_start hit_end score evalue)) {
			$fields->{$field} = shift(@hsp_data);
		}
	}
	
	$self->_dependencies( { ( query_string => 'alignment',
                              hit_string => 'alignment',
                              homology_string => 'alignment',
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

sub _discover_alignment {
    my $self = shift;
    my $alignments_hash = $self->get_field('alignments');
	
    my $identifier = $self->get_field('name').'~~~~'.$self->get_field('rank');
	
    while (! defined $alignments_hash->{$identifier}) {
		last unless $self->parent->parent->_next_alignment;
    }
    my $alignment = $alignments_hash->{$identifier};
	
    if ($alignment) {
		# work out query, hit and homology strings, and some stats
        # (quicker to do this all at once instead of each method working on
        # $alignment string itself)
        
        my ($query_string, $hit_string, $homology_string);
        while ($alignment =~ /\s+(\S+)\n\s+(\S.+)\n\s+\S+\s+\d+\s+(\S+)\s+\d/gm) {
            my $hi = $1;
            my $ho = $2;
            $query_string .= $3;
            
            $hi =~ s/\*\-\>//;
            $ho = ' 'x(length($hi) - length($ho)).$ho;
			$hi =~ s/\<\-\*//;
            
            $hit_string .= $hi;
            $homology_string .= $ho;
        }
        
        $self->_fields->{query_string} = $query_string;
        $self->_fields->{hit_string} = $hit_string;
        $homology_string =~ s/   $//;
        $self->_fields->{homology_string} = $homology_string;
        
        ($self->{_query_gaps}) = $query_string =~ tr/-//;
        ($self->{_hit_gaps}) = $hit_string =~ tr/.//;
        ($self->{_total_gaps}) = $self->{_query_gaps} + $self->{_hit_gaps};
    }
    
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
        if ($schar eq '.') {
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
                   '-strand'   => 1,
                   '-seq_id'   => $self->get_field('query_name'),
                   #'-seqlength'=> $self->get_field('query_length'),  (not known)
                   '-source'   => $self->get_field('algorithm'),
                   '-seqdesc'  => $self->get_field('query_description')
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
		# the full length isn't always known (given in the report), but don't
		# warn about the missing info all the time
		my $verbose = $self->parent->parent->parent->verbose;
		$self->parent->parent->parent->verbose(-1);
		my $seq_length = $self->get_field('length');
		$self->parent->parent->parent->verbose($verbose);
		
        $self->SUPER::hit( new  Bio::SeqFeature::Similarity
                  ('-primary'  => $self->primary_tag,
                   '-start'    => $self->get_field('hit_start'),
                   '-end'      => $self->get_field('hit_end'),
                   '-expect'   => $self->get_field('evalue'),
                   '-score'    => $self->get_field('score'),
                   '-strand'   => 1,
                   '-seq_id'   => $self->get_field('name'),
                   $seq_length ? ('-seqlength' => $seq_length) : (),
                   '-source'   => $self->get_field('algorithm'),
                   '-seqdesc'  => $self->get_field('description')
                   ) );
		$self->{_created_hit} = 1;
    }
    return $self->SUPER::hit(@_);
}

=head2 gaps

 Title    : gaps
 Usage    : my $gaps = $hsp->gaps( ['query'|'hit'|'total'] );
 Function : Get the number of gaps in the query, hit, or total alignment.
 Returns  : Integer, number of gaps or 0 if none
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
    
    $self->get_field('alignment'); # make sure gaps have been calculated
    
    return $self->{'_'.$type.'_gaps'};
}

=head2 pvalue

 Title   : pvalue
 Usage   : my $pvalue = $hsp->pvalue();
 Function: Returns the P-value for this HSP
 Returns : undef (Hmmpfam reports do not have p-values)
 Args    : none

=cut

# noop
sub pvalue {  }

1;
