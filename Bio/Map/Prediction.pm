# $Id: Prediction.pm,v 1.10 2006/09/28 14:09:40 sendu Exp $
#
# BioPerl module for Bio::Map::Prediction
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

Bio::Map::Prediction - An object representing the predictions of something
that can have multiple locations in several maps.

=head1 SYNOPSIS

  use Bio::Map::Prediction;
  use Bio::Map::Position;

  # normally you would get predictions from a run wrapper like
  # Bio::Tools::Run::Meme, but here we create some manually:
  my $pred1 = Bio::Map::Prediction->new(-source => 'meme');
  Bio::Map::Position->new(-element => $prediction1,
						  -map => Bio::Map::GeneMap->get(-gene => 'gene1',
						                                 -species => 'species1'),
						  -start => 950,
						  -end => 960);
  Bio::Map::Position->new(-element => $prediction1,
                          -map => Bio::Map::GeneMap->get(-gene => 'gene1',
						                                 -species => 'species2'),
                          -start => 1950,
                          -end => 1960);
  Bio::Map::Position->new(-element => $prediction1,
                          -map => Bio::Map::GeneMap->get(-gene => 'gene2',
						                                 -species => 'species1'),
                          -start => 955,
                          -end => 965);
  Bio::Map::Position->new(-element => $prediction1,
                          -map => Bio::Map::GeneMap->get(-gene => 'gene2',
	      				                                 -species => 'species2'),
                          -start => 1955,
                          -end => 1965);

  my $pred2 = Bio::Map::Prediction->new(-source => 'gerp');
  Bio::Map::Position->new(-element => $prediction2,
                          -map => Bio::Map::GeneMap->get(-gene => 'gene1',
						                                 -species => 'species1'),
                          -start => 950,
                          -end => 960);
  # etc.

  # find the places where predictions agree
  use Bio::Map::GeneRelative;
  my $rel = Bio::Map::GeneRelative->new(-gene => 0);
  my $di = Bio::Map::Mappable->disconnected_intersections([$pred1, $pred2],
                                                -min_mappables_percent => 100,
                                                -min_map_percent => 100,
                                                -relative => $rel);
  my @positions = $di->get_positions;

=head1 DESCRIPTION

For example, used to model transcription factor binding site predictions, which
can have multiple locations in several maps.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing list.  Your participation is much appreciated.

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

package Bio::Map::Prediction;
use strict;

use base qw(Bio::Map::Mappable);

=head2 new

 Title   : new
 Usage   : my $prediction = Bio::Map::Prediction->new();
 Function: Builds a new Bio::Map::Prediction object
 Returns : Bio::Map::Prediction
 Args    : -name   => string : name of the mappable element
           -id     => string : id of the mappable element
           -source => string : name of the prediction program

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    my ($source) = $self->_rearrange([qw(SOURCE)], @args);
    $self->source($source) if $source;
    
    return $self;
}

=head2 source

 Title   : name
 Usage   : $mappable->name($new_name);
	       my $name = $mappable->name();
 Function: Get/Set the name for this Mappable
 Returns : A scalar representing the current name of this Mappable
 Args    : none to get
           string to set

=cut

sub source {
    my $self = shift;
    if (@_) { $self->{_source} = shift }
    return $self->{_source} || '';
}

1;
