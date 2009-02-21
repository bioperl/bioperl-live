# $Id$
# BioPerl module for Bio::Expression::FeatureGroup
#
# Copyright Allen Day <allenday@ucla.edu>, Stanley Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::FeatureGroup - a set of DNA/RNA features.  ISA
Bio::Expression::FeatureI

=head1 SYNOPSIS

#

=head1 DESCRIPTION

A set of DNA/RNA features.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 
 
Please direct usage questions or support issues to the mailing list:
  
L<bioperl-l@bioperl.org>
  
rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR

Allen Day E<lt>allenday@ucla.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::Expression::FeatureGroup;

use strict;

use base qw(Bio::Root::Root Bio::Expression::FeatureI);
use vars qw($DEBUG);

=head2 new

 Title   : new
 Usage   : $featuregroup = Bio::Expression::FeatureGroup->new(%args);
 Function: create a new featuregroup object
 Returns : a Bio::Expression::FeatureGroup object
 Args    : an optional hash of parameters to be used in initialization:
           -id    --  the featuregroup ID
           -type  --  the featuregroup type

=cut

sub new {
  my($class,@args) = @_;
  my $self = bless {}, $class;
  $self->_initialize(@args);
  return $self;
}

=head2 _initialize

 Title   : _initialize
 Usage   : $featuregroup->_initialize(@args);
 Function: initialize the featuregroup object
 Returns : nothing
 Args    : @args

=cut

sub _initialize{
  my ($self,@args) = @_;
  my %param = @args;

  $self->type($param{-type});
  $self->id($param{-id}    );

  $self->SUPER::_initialize(@args);
  $DEBUG = 1 if( ! defined $DEBUG && $self->verbose > 0);
}

=head2 type

 Title   : type
 Usage   : $featuregroup->type($optional_arg);
 Function: get/set the type of the featuregroup
 Comments: this is probably going to be a string like
           "quality control", "mismatch blah blah", etc.
 Returns : the featuregroup type
 Args    : a new value for the featuregroup type

=cut

sub type {
  my $self = shift;
  $self->{type} = shift if @_;
  return $self->{type};
}

=head2 id

 Title   : id
 Usage   : $featuregroup->id($optional_arg);
 Function: get/set the id of the featuregroup
 Returns : the featuregroup id
 Args    : a new value for the featuregroup id

=cut

sub id {
  my $self = shift;
  $self->{id} = shift if @_;
  return $self->{id};
}


=head2 standard_deviation

 Title   : standard_deviation
 Usage   : $featuregroup->standard_deviation($optional_arg);
 Function: get/set the standard deviation of the featuregroup value
 Returns : the featuregroup standard deviation
 Args    : a new value for the featuregroup standard deviation
 Notes   : this method does no calculation, it merely holds a value

=cut

sub standard_deviation {
  my $self = shift;
  $self->{standard_deviation} = shift if @_;
  return $self->{standard_deviation};
}

=head2 quantitation

 Title   : quantitation
 Usage   : $featuregroup->quantitation($optional_arg);
 Function: get/set the quantitation of the featuregroup
 Returns : the featuregroup's quantitated value
 Args    : a new value for the featuregroup's quantitated value
 Notes   : this method does no calculation, it merely holds a value

=cut

sub quantitation {
  my $self = shift;
  $self->{quantitation} = shift if @_;
  return $self->{quantitation};
}

=head2 quantitation_units

 Title   : quantitation_units
 Usage   : $featuregroup->quantitation_units($optional_arg);
 Function: get/set the quantitation units of the featuregroup
 Returns : the featuregroup's quantitated value units
 Args    : a new value for the featuregroup's quantitated value units

=cut

sub quantitation_units {
  my $self = shift;
  $self->{quantitation_units} = shift if @_;
  return $self->{quantitation_units};
}

=head2 presence

 Title   : presence
 Usage   : $featuregroup->presence($optional_arg);
 Function: get/set the presence call of the featuregroup
 Returns : the featuregroup's presence call
 Args    : a new value for the featuregroup's presence call

=cut

sub presence {
  my $self = shift;
  $self->{presence} = shift if @_;
  return $self->{presence};
}

=head2 add_feature

 Title   : add_feature
 Usage   : $feature_copy = $featuregroup->add_feature($feature);
 Function: add a feature to the featuregroup
 Returns : see this_feature()
 Args    : a Bio::Expression::FeatureI compliant object

=cut

sub add_feature {
  my($self,@args) = @_;
  foreach my $feature (@args){
	$self->throw('Features must be Bio::Expression::FeatureI compliant') unless $feature->isa('Bio::Expression::FeatureI');
    push @{$self->{features}}, $feature;
  }

  return $self->{features} ? $self->{features}->[-1] : undef;
}

=head2 this_feature

 Title   : this_feature
 Usage   : $feature = $featuregroup->this_feature
 Function: access the last feature added to the featuregroup
 Returns : the last feature added to the featuregroup
 Args    : none

=cut

sub this_feature {
  my $self = shift;
  return $self->{features} ? $self->{features}->[-1] : undef;
}

=head2 each_feature

 Title   : each_feature
 Usage   : @features = $featuregroup->each_feature
 Function: returns a list of Bio::Expression::FeatureI compliant
           objects
 Returns : a list of objects
 Args    : none

=cut

sub each_feature {
  my $self = shift;
  return @{$self->{features}} if defined($self->{features});
  return ();
}

=head2 each_feature_quantitation

 Title   : each_feature_quantitation
 Usage   : @featurequantitions = $featuregroup->each_feature_quantitation;
 Function: returns an list of quantitations of the features in the featuregroup
 Returns : a list of numeric values
 Args    : none

=cut

sub each_feature_quantitation {
  my $self = shift;
  my @values = ();
  push @values, $_->value foreach $self->each_feature;
  return @values;
}

=head2 is_qc

 Title   : is_qc
 Usage   : $is_quality_control = $featuregroup->is_qc
 Function: get/set whether or not the featuregroup is used for quality control purposes
 Returns : a boolean (equivalent)
 Args    : a new value

=cut

sub is_qc {
  my $self = shift;
  $self->{is_qc} = shift if defined @_;
  return $self->{is_qc};
}

1;
