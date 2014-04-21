#
# BioPerl module for Bio::Tools::Blat
#
# Written by Balamurugan Kumarasamy
#
# You may distribute this module under the same terms as perl itself
#

=head1 NAME

Bio::Tools::Blat - parser for Blat program

=head1 SYNOPSIS

  use Bio::Tools::Blat;
  my $blat_parser = Bio::Tools::Blat->new(-fh =>$filehandle );
  while( my $blat_feat = $blat_parser->next_result ) {
        push @blat_feat, $blat_feat;
  }

=head1 DESCRIPTION

 Parser for Blat program

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
 of the bugs and their resolution. Bug reports can be submitted the
 web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Balamurugan Kumarasamy

 Email: bala@tll.org.sg

=head1 APPENDIX

 The rest of the documentation details each of the object methods.
 Internal methods are usually preceded with a _

=cut

package Bio::Tools::Blat;
use strict;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::FeaturePair;
use Bio::SeqFeature::Generic;
use base qw(Bio::Root::Root Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Blat->new(-fh=>$filehandle);
 Function: Builds a new Bio::Tools::Blat object
 Returns : Bio::Tools::Blat
 Args    : -filename
           -fh (filehandle)

=cut

sub new {
      my($class,@args) = @_;

      my $self = $class->SUPER::new(@args);
      $self->_initialize_io(@args);

      return $self;
}

=head2 next_result

 Title   : next_result
 Usage   : my $feat = $blat_parser->next_result
 Function: Get the next result set from parser data
 Returns : L<Bio::SeqFeature::Generic>
 Args    : none

=cut

sub next_result {
	my ($self) = @_;
	my $filehandle;
	my $line;
	my $id;

	while ($_=$self->_readline()){
		# first split on spaces:
		$line = $_;
		chomp $line;

		my ($matches, $mismatches, $rep_matches, $n_count, $q_num_insert,
			 $q_base_insert, $t_num_insert, $t_base_insert, $strand, $q_name,
			 $q_length, $q_start, $q_end, $t_name, $t_length, 
			 $t_start, $t_end, $block_count, $block_sizes, $q_starts,
			 $t_starts
			) = split;

		my $superfeature = Bio::SeqFeature::Generic->new();

		# ignore any preceeding text
		next unless ( $matches =~/^\d+$/ );

		# create as many features as blocks there are in each output line
		my (%feat1, %feat2);
		$feat1{name} = $t_name;
		$feat2{name} = $q_name;

		$strand = $1 if ($strand =~/([+-])[+-]/);

		$feat2{strand} = 1;
		$feat1{strand} = $strand;

		my $percent_id = sprintf "%.2f",
		(100 * ($matches + $rep_matches)/( $matches + $mismatches + $rep_matches));

		unless ( $q_length ){
			$self->warn("length of query is zero, something is wrong!");
			next;
		}

		my $score   = sprintf "%.2f",
		(100 * ( $matches + $mismatches + $rep_matches ) / $q_length);

		# size of each block of alignment (inclusive)
		my @block_sizes     = split ",",$block_sizes;

		# start position of each block (you must add 1 as psl output 
		# is off by one in the start coordinate)
		my @q_start_positions = split ",",$q_starts;
		my @t_start_positions = split ",",$t_starts;

		$superfeature->seq_id($q_name);
		$superfeature->score( $score );
		$superfeature->add_tag_value('percent_id',$percent_id);

		# each line of output represents one possible entire aligment 
		# of the query (feat1) and the target(feat2)

		for (my $i=0; $i<$block_count; $i++ ){

			my ($query_start,$query_end);

			if ( $strand eq '+' ){
				$query_start = $q_start_positions[$i] + 1;
				$query_end   = $query_start + $block_sizes[$i] - 1;
			}else{
				$query_end   = $q_length  - $q_start_positions[$i];
				$query_start = $query_end - $block_sizes[$i] + 1;
			}

			#$feat2 {start} = $q_start_positions[$i] + 1;
			#$feat2 {end}   = $feat2{start} + $block_sizes[$i] - 1;
			$feat2 {start} = $query_start;
			$feat2 {end}   = $query_end;
			if ( $query_end <  $query_start ){
				$self->warn("dodgy feature coordinates: end = $query_end, start = $query_start. Reversing...");
				$feat2 {end}   = $query_start;
				$feat2 {start} = $query_end;
			}

			$feat1 {start} = $t_start_positions[$i] + 1;
			$feat1 {end}   = $feat1{start} + $block_sizes[$i] - 1;

			# we put all the features with the same score and percent_id
			$feat2 {score}   = $score;
			$feat1 {score}   = $feat2 {score};
			$feat2 {percent} = $percent_id;
			$feat1 {percent} = $feat2 {percent};

			# other stuff:
			$feat1 {db}         = undef;
			$feat1 {db_version} = undef;
			$feat1 {program}    = 'blat';
			$feat1 {p_version}  = '1';
			$feat1 {source}     = 'blat';
			$feat1 {primary}    = 'similarity';
			$feat2 {source}     = 'blat';
			$feat2 {primary}    = 'similarity';

			my $feature_pair = $self->create_feature(\%feat1, \%feat2);
			$superfeature->add_sub_SeqFeature( $feature_pair,'EXPAND');
		}
		return $superfeature;
	}
}

=head2 create_feature

 Title   : create_feature
 Usage   : my $feat=$blat_parser->create_feature($feature,$seqname)
 Function: creates a SeqFeature Generic object
 Returns : L<Bio::SeqFeature::Generic>
 Args    :


=cut

sub create_feature {
    my ($self, $feat1,$feat2) = @_;
    my $feature1= Bio::SeqFeature::Generic->new(
							  -seq_id     =>$feat1->{name},
							  -start      =>$feat1->{start},
                       -end        =>$feat1->{end},
                       -strand     =>$feat1->{strand},
                       -score      =>$feat1->{score},
                       -source     =>$feat1->{source},
                       -primary    =>$feat1->{primary} );

    my $feature2= Bio::SeqFeature::Generic->new(
                       -seq_id     =>$feat2->{name},
							  -start      =>$feat2->{start},
                       -end        =>$feat2->{end},
                       -strand     =>$feat2->{strand},
                       -score      =>$feat2->{score},
                       -source     =>$feat2->{source},
                       -primary    =>$feat2->{primary} );

    my $featurepair = Bio::SeqFeature::FeaturePair->new;
    $featurepair->feature1 ($feature1);
    $featurepair->feature2 ($feature2);

	 $featurepair->add_tag_value('evalue',$feat2->{p});
	 $featurepair->add_tag_value('percent_id',$feat2->{percent});
	 $featurepair->add_tag_value("hid",$feat2->{primary});
    return  $featurepair;
}

1;
