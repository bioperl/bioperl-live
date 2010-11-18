#
# BioPerl module for Bio::SearchIO::psl
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SearchIO::psl - A parser for PSL output (UCSC)

=head1 SYNOPSIS

  use Bio::SearchIO;
  my $parser = Bio::SearchIO->new(-file   => 'file.psl',
                                 -format => 'psl');
  while( my $result = $parser->next_result ) {
  }

=head1 DESCRIPTION

This is a SearchIO driver for PSL format.
PSL format is documented here:
http://genome.ucsc.edu/goldenPath/help/customTrack.html#PSL

By default it assumes PSL output came from BLAT you can override that
by specifying -program_name =E<gt> 'BLASTZ' when initializing the
SearchIO object.

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

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SearchIO::psl;
use vars qw(%MAPPING %MODEMAP $DEFAULT_WRITER_CLASS $DefaultProgramName);

use strict;
use Bio::Search::HSP::HSPFactory;
use Bio::Search::Hit::HitFactory;
use Bio::Search::Result::ResultFactory;

$DefaultProgramName   = 'BLAT';
$DEFAULT_WRITER_CLASS = 'Bio::SearchIO::Writer::HitTableWriter';

# mapping of terms to Bioperl hash keys
%MODEMAP = (
    'PSLOutput' => 'result',
    'Result'    => 'result',
    'Hit'       => 'hit',
    'Hsp'       => 'hsp'
);

%MAPPING = (
    'Hsp_bit-score'   => 'HSP-bits',
    'Hsp_score'       => 'HSP-score',
    'Hsp_evalue'      => 'HSP-evalue',
    'Hsp_query-from'  => 'HSP-query_start',
    'Hsp_query-to'    => 'HSP-query_end',
    'Hsp_query-strand'=> 'HSP-query_strand',
    'Hsp_hit-from'    => 'HSP-hit_start',
    'Hsp_hit-to'      => 'HSP-hit_end',
    'Hsp_hit-strand'  => 'HSP-hit_strand',
    'Hsp_positive'    => 'HSP-conserved',
    'Hsp_identity'    => 'HSP-identical',
    'Hsp_mismatches'  => 'HSP-mismatches',
    'Hsp_qgapblocks'  => 'HSP-query_gapblocks',
    'Hsp_hgapblocks'  => 'HSP-hit_gapblocks',
    'Hsp_gaps'        => 'HSP-hsp_gaps',
    'Hsp_hitgaps'     => 'HSP-hit_gaps',
    'Hsp_querygaps'   => 'HSP-query_gaps',
    'Hsp_align-len'   => 'HSP-hsp_length',
    'Hsp_query-frame' => 'HSP-query_frame',
    'Hsp_hit-frame'   => 'HSP-hit_frame',

    'Hit_id'        => 'HIT-name',
    'Hit_len'       => 'HIT-length',
    'Hit_accession' => 'HIT-accession',
    'Hit_def'       => 'HIT-description',
    'Hit_signif'    => 'HIT-significance',
    'Hit_score'     => 'HIT-score',
    'Hit_bits'      => 'HIT-bits',

    'PSLOutput_program'   => 'RESULT-algorithm_name',
    'PSLOutput_version'   => 'RESULT-algorithm_version',
    'PSLOutput_query-def' => 'RESULT-query_name',
    'PSLOutput_query-len' => 'RESULT-query_length',
    'PSLOutput_query-acc' => 'RESULT-query_accession',
    'PSLOutput_querydesc' => 'RESULT-query_description',
    'PSLOutput_db'        => 'RESULT-database_name',
    'PSLOutput_db-len'    => 'RESULT-database_entries',
    'PSLOutput_db-let'    => 'RESULT-database_letters',
);

use base qw(Bio::SearchIO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::SearchIO::psl->new();
 Function: Builds a new Bio::SearchIO::psl object 
 Returns : an instance of Bio::SearchIO::psl
 Args    :


=cut

sub _initialize {
    my ( $self, @args ) = @_;
    $self->SUPER::_initialize(@args);
    my ($pname, $distinct, $qt, $ht) = $self->_rearrange( [qw(PROGRAM_NAME
                                                 DISTINCT_LOCATIONS
                                                 QUERY_TYPE
                                                 HIT_TYPE)], @args );
    $self->program_name( $pname || $DefaultProgramName );
    $distinct && $self->distinct_locations(1);
    $qt ||= 'dna';
    $ht ||= 'dna';
    $self->query_type($qt);
    $self->hit_type($ht);
    
    $self->_eventHandler->register_factory(
        'result',
        Bio::Search::Result::ResultFactory->new(
            -type => 'Bio::Search::Result::GenericResult'
        )
    );

    $self->_eventHandler->register_factory(
        'hit',
        Bio::Search::Hit::HitFactory->new(
            -type => 'Bio::Search::Hit::GenericHit'
        )
    );
    $self->_eventHandler->register_factory(
        'hsp',
        Bio::Search::HSP::HSPFactory->new(
            -type => 'Bio::Search::HSP::PSLHSP'
        )
    );
}

=head2 next_result

 Title   : next_result
 Usage   : my $result = $parser->next_result
 Function: Parse the next result from the data stream
 Returns : L<Bio::Search::Result::ResultI> or undef if no more results
 Args    : none

=cut

sub next_result {
    my ($self) = @_;
    unless (exists $self->{'_handlerset'}) {
        my $line;
        while ($line = $self->_readline) {
            # advance to first line
            next if $line =~ m{^\s*$};
            # newer output starts with model name
            if ($line =~ m{^\#\s+cmsearch\s}) {
                $self->{'_handlerset'} = 'latest';
			} elsif ($line =~ m{^CM\s\d+:}) {
                $self->{'_handlerset'} = 'pre-1.0';
            } else {
                $self->{'_handlerset'} ='old';
            }
            last;
        }
        $self->_pushback($line);
    }
    return $self->{'distinct_locations'} ? $self->_parse_distinct :
		   $self->_parse_simple;
}

sub _parse_simple {
    my ($self) = @_;
    my ( $lastquery, $lasthit );
    local $/ = "\n";
    local $_;

    # skip over any header lines
    while( defined($_ = $self->_readline) and ! /^\d+\s+\d+\s+/ ) {}
    $self->_pushback($_);

    # now start the main parsing loop
    while ( defined( $_ = $self->_readline ) ) {
        my (
            $matches,      $mismatches,    $rep_matches,  $n_count,
            $q_num_insert, $q_base_insert, $t_num_insert, $t_base_insert,
            $strand,       $q_name,        $q_length,     $q_start,
            $q_end,        $t_name,        $t_length,     $t_start,
            $t_end,        $block_count,   $block_sizes,  $q_starts,
            $t_starts
        ) = split(/\s/,$_);
        
        # TODO: what do we do with strand of '++' or similar?

        $q_length > 0 or $self->throw("parse error, invalid query length '$q_length'");
        my $score = sprintf( "%.2f",  100 * ( $matches + $mismatches + $rep_matches ) / $q_length );

        # this is overall percent identity...
        my $match_total  = $matches + $mismatches + $rep_matches;
        $match_total > 0
            or $self->throw("parse error, matches + mismatches + rep_matches must be > 0!");
        my $percent_id = sprintf("%.2f", 100 * ( $matches + $rep_matches ) / $match_total );

        # Remember Jim's code is 0 based
        if ( defined $lastquery
            && $lastquery ne $q_name )
        {
            $self->end_element( { 'Name' => 'Hit' } );
            $self->end_element( { 'Name' => 'PSLOutput' } );
            $self->_pushback($_);
            return $self->end_document;
        }
        elsif ( !defined $lastquery ) {
            $self->{'_result_count'}++;
            $self->start_element( { 'Name' => 'PSLOutput' } );
            $self->element_hash({
                    'PSLOutput_program' => $self->program_name,
                    'PSLOutput_query-def' => $q_name,
                    'PSLOutput_query-len' => $q_length
                }
            );
            $self->start_element( { 'Name' => 'Hit' } );
            $self->element_hash({
                    'Hit_id' => $t_name,
                    'Hit_len' => $t_length,
                    'Hit_score' => $score
                }
            );
        }
        elsif ( $lasthit ne $t_name ) {
            $self->end_element( { 'Name' => 'Hit' } );
            $self->start_element( { 'Name' => 'Hit' } );
            $self->element_hash({
                    'Hit_id' => $t_name,
                    'Hit_len' => $t_length,
                    'Hit_score' => $score
                }
            );
        }

        my $identical = $matches + $rep_matches;
        $self->start_element( { 'Name' => 'Hsp' } );
        $self->element_hash({
                'Hsp_score'  => $score,
                'Hsp_identity' => $identical,
                'Hsp_positive' => $identical,
                'Hsp_mismatches' => $mismatches,
                'Hsp_gaps' => $q_base_insert + $t_base_insert,
                'Hsp_querygaps' => $t_base_insert,
                'Hsp_hitgaps' => $q_base_insert,
                'Hsp_query-from' => ( $strand eq '+' ) ? $q_start + 1 : $q_end,
                'Hsp_query-to' => ( $strand eq '+' ) ? $q_end :$q_start + 1,
            }
        );
        my $hsplen =
          $q_base_insert +
          $t_base_insert +
          abs( $t_end - $t_start ) +
          abs( $q_end - $q_start );
        $self->element_hash({
                'Hsp_hit-from' => $t_start + 1,
                'Hsp_hit-to'   => $t_end,
                'Hsp_align-len' => $hsplen
            }
        );

        # cleanup trailing commas in some output
        $block_sizes =~ s/\,$//;
        $q_starts    =~ s/\,$//;
        $t_starts    =~ s/\,$//;
        my @blocksizes = split( /,/, $block_sizes );    # block sizes
        my @qstarts = split( /,/, $q_starts ); # starting position of each block
                                               # in query
        my @tstarts = split( /,/, $t_starts ); # starting position of each block
                                               # in target
        my ( @qgapblocks, @hgapblocks );

        for ( my $i = 0 ; $i < $block_count ; $i++ ) {
            if ( $strand eq '+' ) {
                push @qgapblocks, [ $qstarts[$i] + 1, $blocksizes[$i] ];
            }
            else {
                push @qgapblocks, [ $q_length - $qstarts[$i], $blocksizes[$i] ];
            }
            push @hgapblocks, [ $tstarts[$i] + 1, $blocksizes[$i] ];
        }
        $self->element_hash({
                'Hsp_qgapblocks' => \@qgapblocks,
                'Hsp_hgapblocks' => \@hgapblocks
            }
        );
        $self->end_element( { 'Name' => 'Hsp' } );
        $lastquery = $q_name;
        $lasthit   = $t_name;
    }
    if ( defined $lasthit || defined $lastquery ) {
        $self->end_element( { 'Name' => 'Hit' } );
        $self->end_element( { 'Name' => 'Result' } );
        return $self->end_document;
    }
}

# will be merging this and the above into a more coherent group, but
# keeping separate for now.  Noticably, this wasn't handling more
# complex output (BLATX, TBLATN, TBLATX, etc)

sub _parse_distinct {
    my ($self) = @_;
    my ( $lastquery );
    local $/ = "\n";
    local $_;
    
    # if these aren't passed in, we can only assume both of these are
    # dna (the default)
    my ($qt, $ht) = ($self->query_type, $self->hit_type);
    
    # skip over any header lines
    while( defined($_ = $self->_readline) and ! /^\d+\s+\d+\s+/ ) {}
    $self->_pushback($_);

    # now start the main parsing loop
    while ( defined( $_ = $self->_readline ) ) {
        my (
            $matches,      $mismatches,    $rep_matches,  $n_count,
            $q_num_insert, $q_base_insert, $t_num_insert, $t_base_insert,
            $strand,       $q_name,        $q_length,     $q_start,
            $q_end,        $t_name,        $t_length,     $t_start,
            $t_end,        $block_count,   $block_sizes,  $q_starts,
            $t_starts
        ) = split(/\s/,$_);

        $q_length > 0 or $self->throw("parse error, invalid query length '$q_length'");
        my $score = sprintf( "%.2f",  100 * ( $matches + $mismatches + $rep_matches ) / $q_length );

        # this is overall percent identity...
        my $match_total  = $matches + $mismatches + $rep_matches;
        $match_total > 0
            or $self->throw("parse error, matches + mismatches + rep_matches must be > 0!");
        my $percent_id = sprintf("%.2f", 100 * ( $matches + $rep_matches ) / $match_total );

        # Remember Jim's code is 0 based
        if ( defined $lastquery && $lastquery ne $q_name )
        {
            $self->end_element( { 'Name' => 'PSLOutput' } );
            $self->_pushback($_);
            return $self->end_document;
        }
        else {
            if (!$self->within_element('result')) {
                $self->{'_result_count'}++;
                $self->start_element( { 'Name' => 'PSLOutput' } );
                $self->element_hash({
                        'PSLOutput_program' => $self->program_name,
                        'PSLOutput_query-def' => $q_name,
                        'PSLOutput_query-len' => $q_length
                    }
                );
            }
            $self->start_element( { 'Name' => 'Hit' } );
            $self->element_hash({
                    'Hit_id' => $t_name,
                    'Hit_len' => $t_length,
                    'Hit_score' => $score
                }
            );

            #my $identical = $matches + $rep_matches;
            
            my @blocksizes = split( /,/, $block_sizes );    # block sizes
            my @qstarts = split( /,/, $q_starts ); # starting position of each block
                                                   # in query
            my @tstarts = split( /,/, $t_starts ); # starting position of each block
                                                   # in target
            my ($qstrand, $hstrand) = length($strand) > 1 ?
                split('',$strand,2) : ($strand, '+');
                
            for my $i ( 0..$#blocksizes ) {
                
                # TODO: start/end coordinates may need to be adjusted based on
                # the type of BLAT run; I don't think there is a way to tell
                # this via PSL output, so this will have to be specified as
                # an object attribute or based off that.
                
                $self->start_element( { 'Name' => 'Hsp' } );
                                
                my ($qstart, $qend) = ( $qstrand eq '+' ) ?
                    ($qstarts[$i] + 1, $qstarts[$i] + $blocksizes[$i]) :
                    ($q_length - $qstarts[$i] - $blocksizes[$i] + 1, $q_length - $qstarts[$i]);
                my ($hstart, $hend) = ($tstarts[$i] + 1, $tstarts[$i] + $blocksizes[$i]);
                
                $self->debug(sprintf("Q: %d-%d:%s\nH: %d-%d:%s\n",
                                     $qstart, $qend, $qstrand,
                                     $hstart, $hend, $hstrand, $strand)) if $qstrand eq '-';
                
                # no score for HSPs in this case
                
                $self->element_hash({
                        'Hsp_query-from' => $qstart,
                        'Hsp_query-to'   => $qend,
                        'Hsp_query-strand'   => $qstrand,
                        'Hsp_hit-from'   => $hstart,
                        'Hsp_hit-to'     => $hend,
                        'Hsp_hit-strand'   => $hstrand,
                    }
                );
                #$self->element(
                #    {
                #        'Name' => 'Hsp_score',
                #        'Data' => $score
                #    }
                #);
                #$self->element(
                #    {
                #        'Name' => 'Hsp_identity',
                #        'Data' => $identical
                #    }
                #);
                #$self->element(
                #    {
                #        'Name' => 'Hsp_positive',
                #        'Data' => $identical
                #    }
                #);
                #$self->element(
                #    {
                #        'Name' => 'Hsp_mismatches',
                #        'Data' => $mismatches
                #    }
                #);
                #$self->element(
                #    {
                #        'Name' => 'Hsp_gaps',
                #        'Data' => $q_base_insert + $t_base_insert
                #    }
                #);
                #
                ## query gaps are the number of target inserts and vice-versa
                #$self->element(
                #    {
                #        'Name' => 'Hsp_querygaps',
                #        'Data' => $t_base_insert
                #    }
                #);
                #$self->element(
                #    {
                #        'Name' => 'Hsp_hitgaps',
                #        'Data' => $q_base_insert
                #    }
                #);
                #if ( $strand eq '+' ) {
                #    $self->element(
                #        {
                #            'Name' => 'Hsp_query-from',
                #            'Data' => $q_start + 1
                #        }
                #    );
                #    $self->element(
                #        {
                #            'Name' => 'Hsp_query-to',
                #            'Data' => $q_end
                #        }
                #    );
                #}
                #else {
                #    $self->element(
                #        {
                #            'Name' => 'Hsp_query-to',
                #            'Data' => $q_start + 1
                #        }
                #    );
                #    $self->element(
                #        {
                #            'Name' => 'Hsp_query-from',
                #            'Data' => $q_end
                #        }
                #    );
                #}
                #my $hsplen =
                #  $q_base_insert +
                #  $t_base_insert +
                #  abs( $t_end - $t_start ) +
                #  abs( $q_end - $q_start );
                #$self->element(
                #    {
                #        'Name' => 'Hsp_hit-from',
                #        'Data' => $t_start + 1
                #    }
                #);
                #$self->element(
                #    {
                #        'Name' => 'Hsp_hit-to',
                #        'Data' => $t_end
                #    }
                #);
                #$self->element(
                #    {
                #        'Name' => 'Hsp_align-len',
                #        'Data' => $hsplen
                #    }
                #);
                
                $self->end_element( { 'Name' => 'Hsp' } );
            }
            $lastquery = $q_name;
            $self->end_element( { 'Name' => 'Hit' } );
        }
    }
    if ( defined $lastquery ) {
        $self->end_element( { 'Name' => 'Result' } );
        return $self->end_document;
    }

}

=head2 element_hash

 Title   : element
 Usage   : $eventhandler->element_hash({'Hsp_hit-from' => $start,
                                        'Hsp_hit-to'   => $end,
                                        'Hsp_score'    => $lastscore});
 Function: Convenience method that takes multiple simple data elements and
           maps to appropriate parameters
 Returns : none
 Args    : Hash ref with the mapped key (in %MAPPING) and value

=cut

sub element_hash {
    my ($self, $data) = @_;
    $self->throw("Must provide data hash ref") if !$data || !ref($data);
    for my $nm (sort keys %{$data}) {
        next if $data->{$nm} && $data->{$nm} =~ m{^\s*$}o;
        if ( $MAPPING{$nm} ) {
            if ( ref( $MAPPING{$nm} ) =~ /hash/i ) {
                my $key = ( keys %{ $MAPPING{$nm} } )[0];
                $self->{'_values'}->{$key}->{ $MAPPING{$nm}->{$key} } =
                  $data->{$nm};
            }
            else {
                $self->{'_values'}->{ $MAPPING{$nm} } = $data->{$nm};
            }
        }
    }
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
    my $nm = $data->{'Name'};
    if ( my $type = $MODEMAP{$nm} ) {
        $self->_mode($type);
        if ( $self->_eventHandler->will_handle($type) ) {
            my $func = 'start_'.lc $type;
            $self->_eventHandler->$func( $data->{'Attributes'} );
        }
        unshift @{ $self->{'_elements'} }, $type;
    }
    if ( $nm eq 'PSLOutput' ) {
        $self->{'_values'} = {};
        $self->{'_result'} = undef;
        $self->{'_mode'}   = '';
    }
}

=head2 end_element

 Title   : end_element
 Usage   : $eventgenerator->end_element
 Function: Handles an end element event
 Returns : return value from the associated end_$type event handler
 Args    : hashref with at least 2 keys 'Data' and 'Name'


=cut

sub end_element {
    my ( $self, $data ) = @_;
    my $nm = $data->{'Name'};
    my $rc;

    # Hsp are sort of weird, in that they end when another
    # object begins so have to detect this in end_element for now

    if ( my $type = $MODEMAP{$nm} ) {
        if ( $self->_eventHandler->will_handle($type) ) {
            my $func = 'end_'.lc $type;
            $rc = $self->_eventHandler->$func( $self->{'_reporttype'},
                $self->{'_values'} );
        }
        shift @{ $self->{'_elements'} };

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
        $self->warn(
            __PACKAGE__ . "::end_element: unknown nm '$nm', ignoring\n" );
    }
    $self->{'_last_data'} = '';    # remove read data if we are at
                                   # end of an element
    $self->{'_result'}    = $rc
      if ( defined $nm
        && defined $MODEMAP{$nm}
        && $MODEMAP{$nm} eq 'result' );
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

    return unless ( defined $data->{'Data'} );
    if ( $data->{'Data'} =~ /^\s+$/ ) {
        return unless $data->{'Name'} =~ /Hsp\_(midline|qseq|hseq)/;
    }

    if (   $self->in_element('hsp')
        && $data->{'Name'} =~ /Hsp\_(qseq|hseq|midline)/ )
    {

        $self->{'_last_hspdata'}->{ $data->{'Name'} } .= $data->{'Data'};
    }

    $self->{'_last_data'} = $data->{'Data'};
}

=head2 _mode

 Title   : _mode
 Usage   : $obj->_mode($newval)
 Function: 
 Example : 
 Returns : value of _mode
 Args    : newvalue (optional)


=cut

sub _mode {
    my ( $self, $value ) = @_;
    if ( defined $value ) {
        $self->{'_mode'} = $value;
    }
    return $self->{'_mode'};
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
    return 0 if (!defined $name || !defined($self->{'_elements'})
        || scalar @{ $self->{'_elements'} } == 0 );
    foreach ( @{ $self->{'_elements'} } ) {
        if ( $_ eq $name ) {
            return 1;
        }
    }
    return 0;
}

=head2 in_element

 Title   : in_element
 Usage   : if( $eventgenerator->in_element($element) ) {}
 Function: Test if we are in a particular element
           This is different than 'in' because within can be tested
           for a whole block.
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
 Function: Handles a start document event
 Returns : none
 Args    : none


=cut

sub start_document {
    my ($self) = @_;
    $self->{'_lasttype'} = '';
    $self->{'_values'}   = {};
    $self->{'_result'}   = undef;
    $self->{'_mode'}     = '';
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
    my ( $self, @args ) = @_;
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

sub report_count { shift->result_count }

=head2 program_name

 Title   : program_name
 Usage   : $obj->program_name($newval)
 Function: Get/Set the program name
 Returns : value of program_name (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub program_name {
    my $self = shift;

    $self->{'program_name'} = shift if @_;
    return $self->{'program_name'} || $DefaultProgramName;
}

=head2 distinct_locations

The original PSL parser parses very simple start/end information, but 
doesn't disambiguate blocks into HSPs (though the data is stored)
Therefore, every time the query name changes, a new Result object is
generated; if the target ID stays the same, these are clustered together
as one Hit object with two distinct HSP regions.

This has the advantage of clustering Results/Hits/HSPs similar to BLAST
reports, but the important disadvantage is that the blocks cannot be further
disambiguated into subparts below that of HSPs. Since BLAT can 

=cut

sub distinct_locations {
    my $self = shift;
    return $self->{distinct_locations} = shift if @_;
    return $self->{distinct_locations};
}

sub hit_type {
    my $self = shift;
    return $self->{query_type} = shift if @_;
    return $self->{query_type};
}

sub query_type {
    my $self = shift;
    return $self->{hit_type} = shift if @_;
    return $self->{hit_type};
}

1;
