
=head1 NAME

Bio::Tools::TandemRepeatsFinder - a parser for Tandem Repeats Finder output

=head1 SYNOPSIS

    use Bio::Tools::TandemRepeatsFinder;

    # create parser
    my $parser = Bio::Tools::Bio::Tools::TandemRepeatsFinder->new(-file => 'tandem_repeats.out');

    # loop through results
    while( my $feature = $parser->next_result ) {

        # print the source sequence id, start, end, percent matches, and the consensus sequence
        my ($percent_matches)    = $feat->get_tag_values('percent_matches');
        my ($consensus_sequence) = $feat->get_tag_values('consensus_sequence');
        print $feat->seq_id()."\t".$feat->start()."\t".$feat->end()."\t$percent_matches\t$consensus_sequence\n"; 

    }

=head1 DESCRIPTION

A parser for Tandem Repeats Finder output.  
Written and tested for version 4.00

Location, seq_id, and score are stored in Bio::SeqFeature::Generic feature.
All other data is stored in tags.  The availabale tags are

        period_size
        copy_number
        consensus_size
        percent_matches
        percent_indels
        percent_a
        percent_c
        percent_g
        percent_t
        entropy
        consensus_sequence
        repeat_sequence
        run_parameters
        sequence_description

The run_parameters are stored in a hashref with the following key:

        match_weight
        mismatch_weight
        indel_weight
        match_prob
        indel_prob
        min_score
        max_period_size

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Eric Just

Email e-just@northwestern.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::Tools::TandemRepeatsFinder;
use strict;
use constant DEBUG => 0;
use Bio::SeqFeature::Generic;

use base qw(Bio::Root::Root Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::TandemRepeatsFinder->new();
 Function: Builds a new Bio::Tools::TandemRepeatsFinder object
 Returns : Bio::Tools::TandemRepeatsFinder
 Args    : -fh/-file => $val, for initing input, see Bio::Root::IO

=cut

sub new {
    my ( $class, @args ) = @_;

    my $self = $class->SUPER::new(@args);
    $self->_initialize_io(@args);

    return $self;
}

=head2 version

 Title   : version
 Usage   : $self->version( $version )
 Function: get/set the version of Tandem Repeats finder that was used in analysis
 Returns : value of version of 
 Args    : new value (optional)

=cut

sub version {
    my ( $self, $value ) = @_;
    if ( defined $value ) {
        $self->{'version'} = $value;
    }
    return $self->{'version'};
}

=head2 _current_seq_id

 Title   : _current_seq_id
 Usage   : $self->_current_seq_id( $current_seq_id )
 Function: get/set the _current_seq_id
 Returns : value of _current_seq_id
 Args    : new value (optional)

=cut

sub _current_seq_id {
    my ( $self, $value ) = @_;
    if ( defined $value ) {
        $self->{'_current_seq_id'} = $value;
    }
    return $self->{'_current_seq_id'};
}

=head2 _current_seq_description

 Title   : _current_seq_description
 Usage   : $self->_current_seq_description( $current_seq_id )
 Function: get/set the _current_seq_description
 Returns : value of _current_seq_description
 Args    : new value (optional)

=cut

sub _current_seq_description {
    my ( $self, $value ) = @_;
    if ( defined $value ) {
        $self->{'_current_seq_description'} = $value;
    }
    return $self->{'_current_seq_description'};
}

=head2 _current_parameters

 Title   : _current_parameters
 Usage   : $self->_current_parameters( $parameters_hashref )
 Function: get/set the _current_parameters
 Returns : hashref representing current parameters parsed from results file
         : keys are
               match_weight
               mismatch_weight
               indel_weight
               match_prob
               indel_prob
               min_score
               max_period_size
 Args    : parameters hashref (optional)

=cut

sub _current_parameters {
    my ( $self, $value ) = @_;
    if ( defined $value ) {
        $self->{'_current_parameters'} = $value;
    }
    return $self->{'_current_parameters'};
}

=head2 next_result

 Title   : next_result
 Usage   : my $r = $trf->next_result()
 Function: Get the next result set from parser data
 Returns : Bio::SeqFeature::Generic
 Args    : none

=cut

sub next_result {
    my ($self) = @_;
    while ( defined( $_ = $self->_readline() ) ) {

        # Parse Version line
        if (/^Version (.+)/) {
            my $version = $1;
            $self->warn("parsed version: $version\n") if DEBUG;
            $self->warn( qq{ Bio::Tools::TandemRepeatsFinder was written and tested for Tandem Repeats Masker Version 4.00 output
You appear to be using Verion $version.  Use at your own risk.}) if ($version != 4);
            $self->version($version);
        }

        # Parse Sequence identifier
        # i.e. Sequence: DDB0215018 |Masked Chromosomal Sequence| Chr 2f
        elsif ( /^Sequence: ([^\s]+)\s(.+)?/ ) {
            my $seq_id          = $1;
            my $seq_description = $2;
            $self->warn("parsed sequence_id: $seq_id\n") if DEBUG;
            $self->_current_seq_id($seq_id);
            $self->_current_seq_description($seq_description);
        }

        # Parse Parameters
        # i.e. Parameters: 2 7 7 80 10 50 12
        elsif (/^Parameters: (.+)/) {
            my $params = $1;
            $self->warn("parsed parameters: $params\n") if DEBUG;

            my @param_array = split /\s/, $params;

            my $param_hash = {
                match_weight    => $param_array[0],
                mismatch_weight => $param_array[1],
                indel_weight    => $param_array[2],
                match_prob      => $param_array[3],
                indel_prob      => $param_array[4],
                min_score       => $param_array[5],
                max_period_size => $param_array[6]
            };
            $self->_current_parameters($param_hash);
        }

        # Parse Data
        # i.e. 13936 13960 12 2.1 12 100 0 50 16 8 52 24 1.70 T TTTTTTTTTT
        elsif (/^\d+\s\d+\s\d+/) {

            # call internal method to create Bio::SeqFeature::Generic
            # to represent tandem repeat
            return $self->_create_feature($_);
        }

        elsif (DEBUG) {
            $self->warn( "UNPARSED LINE:\n" . $_ );
        }
    }
    return;
}

=head2 _create_feature

 Title   : _create_feature
 Usage   : internal method used by 'next_feature'
 Function: Takes a line from the results file and creates a bioperl object
 Returns : Bio::SeqFeature::Generic
 Args    : none

=cut

sub _create_feature {
    my ( $self, $line ) = @_;

    # split the line and store into named variables
    my @element = split /\s/, $line;
    my (
        $start,          $end,                $period_size,
        $copy_number,    $consensus_size,     $percent_matches,
        $percent_indels, $score,              $percent_a,
        $percent_c,      $percent_g,          $percent_t,
        $entropy,        $consensus_sequence, $repeat_sequence
    ) = @element;

    # create tag hash from data in line
    my $tags = {
        period_size          => $period_size,
        copy_number          => $copy_number,
        consensus_size       => $consensus_size,
        percent_matches      => $percent_matches,
        percent_indels       => $percent_indels,
        percent_a            => $percent_a,
        percent_c            => $percent_c,
        percent_g            => $percent_g,
        percent_t            => $percent_t,
        entropy              => $entropy,
        consensus_sequence   => $consensus_sequence,
        repeat_sequence      => $repeat_sequence,
        run_parameters       => $self->_current_parameters(),
        sequence_description => $self->_current_seq_description()
    };

    # create feature from start/end etc
    my $feat = Bio::SeqFeature::Generic->new(
        -seq_id      => $self->_current_seq_id(),
        -score       => $score,
        -start       => $start,
        -end         => $end,
        -source_tag  => 'Tandem Repeats Finder',
        -primary_tag => 'tandem repeat',
        -tag         => $tags
    );

    return $feat;

}

1;

