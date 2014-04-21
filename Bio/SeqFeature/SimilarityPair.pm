#
# BioPerl module for Bio::SeqFeature::SimilarityPair
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::SimilarityPair - Sequence feature based on the similarity
                  of two sequences.

=head1 SYNOPSIS

  $sim_pair = Bio::SeqFeature::SimilarityPair->from_searchResult($blastHit);

  $sim = $sim_pair->query(); # a Bio::SeqFeature::Similarity object - the query
  $sim = $sim_pair->hit();   # dto - the hit.

  # some properties for the similarity pair
  $expect = $sim_pair->significance();
  $score = $sim_pair->score();
  $bitscore = $sim_pair->bits();

  # this will not write the description for the sequence (only its name)
  print $sim_pair->query()->gff_string(), "\n";

=head1 DESCRIPTION

Lightweight similarity search result as a pair of Similarity
features. This class inherits off Bio::SeqFeature::FeaturePair and
therefore implements Bio::SeqFeatureI, whereas the two features of the
pair are descendants of Bio::SeqFeature::Generic, with better support
for representing similarity search results in a cleaner way.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Hilmar Lapp

Email hlapp@gmx.net or hilmar.lapp@pharma.novartis.com

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::SimilarityPair;
use strict;

use Bio::SeqFeature::Similarity;
use Bio::Factory::ObjectFactory;

use base qw(Bio::SeqFeature::FeaturePair);

=head2 new

 Title   : new
 Usage   : my $similarityPair = Bio::SeqFeature::SimilarityPair->new
                                 (-hit   => $hit,
                                  -query => $query,
                                  -source => 'blastp');
 Function: Initializes a new SimilarityPair object
 Returns : Bio::SeqFeature::SimilarityPair
 Args    : -query => The query in a Feature pair 
           -hit   => (formerly '-subject') the subject/hit in a Feature pair


=cut

sub new {
    my($class,@args) = @_;

    if(! grep { lc($_) eq "-feature_factory"; } @args) {
	# if no overriding factory is provided, provide our preferred one
	my $fact = Bio::Factory::ObjectFactory->new(
                                    -type => "Bio::SeqFeature::Similarity",
				    -interface => "Bio::SeqFeatureI");
	push(@args, '-feature_factory', $fact);
    }
    my $self = $class->SUPER::new(@args);

    my ($primary, $hit, $query, $fea1, $source,$sbjct) =
        $self->_rearrange([qw(PRIMARY
                              HIT
                              QUERY
                              FEATURE1
                              SOURCE
                              SUBJECT
                              )],@args);
    
    if( $sbjct ) { 
        # undeprecated by Jason before 1.1 release 
        # $self->deprecated("use of -subject deprecated: SimilarityPair now uses 'hit'");
        if(! $hit) { $hit = $sbjct } 
        else { 
            $self->warn("-hit and -subject were specified, using -hit and ignoring -subject");
        }
    }

    # set the query and subject feature if provided
    $self->query( $query) if $query && ! $fea1;
    $hit && $self->hit($hit);

    # the following refer to feature1, which is guaranteed to exist
    if( defined $primary || ! defined $self->primary_tag) { 
        $primary = 'similarity' unless defined $primary;
        $self->primary_tag($primary);
    } 

    $source && $self->source_tag($source);

    return $self;
}

#
# Everything else is just inherited from SeqFeature::FeaturePair.
#

=head2 query

 Title   : query
 Usage   : $query_feature = $obj->query();
           $obj->query($query_feature);
 Function: The query object for this similarity pair
 Returns : Bio::SeqFeature::Similarity
 Args    : [optional] Bio::SeqFeature::Similarity

See L<Bio::SeqFeature::Similarity>, L<Bio::SeqFeature::FeaturePair>

=cut

sub query {
    return shift->feature1(@_);
}




=head2 subject

 Title   : subject
 Usage   : $sbjct_feature = $obj->subject();
           $obj->subject($sbjct_feature);
 Function: Get/Set Subject for a SimilarityPair 
 Returns : Bio::SeqFeature::Similarity
 Args    : [optional] Bio::SeqFeature::Similarity
 Notes   : Deprecated.  Use the method 'hit' instead

=cut

sub subject { 
    my $self = shift;
#    $self->deprecated("Method subject deprecated: use hit() instead");
    $self->hit(@_); 
}

=head2 hit

 Title   : hit
 Usage   : $sbjct_feature = $obj->hit();
           $obj->hit($sbjct_feature);
 Function: Get/Set Hit for a SimilarityPair 
 Returns : Bio::SeqFeature::Similarity
 Args    : [optional] Bio::SeqFeature::Similarity


=cut

sub hit {
    return shift->feature2(@_);
}

=head2 source_tag

 Title   : source_tag
 Usage   : $source = $obj->source_tag(); # i.e., program
           $obj->source_tag($evalue);
 Function: Gets the source tag (program name typically) for a feature 
 Returns : string
 Args    : [optional] string


=cut

sub source_tag {
    my ($self, @args) = @_;

    if(@args) {
        $self->hit()->source_tag(@args);
    }
    return $self->query()->source_tag(@args);
}

=head2 significance

 Title   : significance
 Usage   : $evalue = $obj->significance();
           $obj->significance($evalue);
 Function: 
 Returns : 
 Args    : 


=cut

sub significance {
    my ($self, @args) = @_;

    if(@args) {
        $self->hit()->significance(@args);
    }
    return $self->query()->significance(@args);
}

=head2 score

 Title   : score
 Usage   : $score = $obj->score();
           $obj->score($value);
 Function: 
 Returns : 
 Args    : 


=cut

sub score {
    my ($self, @args) = @_;

    if(@args) {
        $self->hit()->score(@args);
    }
    # Note: You might think it's only getting set on the hit object.
    # Actually, it's getting set on both hit and query.

    return $self->query()->score(@args);
}

=head2 bits

 Title   : bits
 Usage   : $bits = $obj->bits();
           $obj->bits($value);
 Function: 
 Returns : 
 Args    : 


=cut

sub bits {
    my ($self, @args) = @_;

    if(@args) {
        $self->hit()->bits(@args);
    }
    return $self->query()->bits(@args);
}

#################################################################
# aliases for backwards compatibility or convenience            #
#################################################################

*sbjct = \&subject;

1;
