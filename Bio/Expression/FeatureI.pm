# $Id$
# BioPerl module for Bio::Expression::FeatureI
#
# Copyright Allen Day <allenday@ucla.edu>, Stan Nelson <snelson@ucla.edu>
# Human Genetics, UCLA Medical School, University of California, Los Angeles

# POD documentation - main docs before the code

=head1 NAME

Bio::Expression::FeatureI - an interface class for DNA/RNA features

=head1 SYNOPSIS

Do not use this module directly

=head1 DESCRIPTION

This provides a standard bioperl interface class for representing
DNA and RNA features.  It cannot be instantiated directly, but serves
as an abstract base class for implementors.

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
package Bio::Expression::FeatureI;

use strict;
use Bio::Root::RootI;

use base qw(Bio::Root::RootI Bio::PrimarySeqI);
use vars qw($DEBUG);

=head2 quantitation()

  Title   : value
  Usage   : $val = $ftr->quantitation()
  Function: get/set the feature's quantitation
  Returns : A numeric value
  Args    : a new numeric value (optional)

=cut

sub quantitation {
  shift->throw_not_implemented();
}

=head2 quantitation_units()

  Title   : quantitation_units
  Usage   : $units = $ftr->quantitation_units()
  Function: get/set the units of the feature's quantitation
  Returns : A string or undef
  Args    : a new string (optional)

=cut

sub quantitation_units {
  shift->throw_not_implemented();
}

=head2 standard_deviation()

  Title   : standard_deviation
  Usage   : $std_dev = $ftr->standard_deviation()
  Function: get/set the feature's standard deviation of quantitation()
  Returns : A numeric value
  Args    : a new numeric value (optional)
  Comments: no calculation is done here

=cut

sub standard_deviation {
  shift->throw_not_implemented();
}

=head2 sample_count()

  Title   : sample_count
  Usage   : $sample_count = $ftr->sample_count()
  Function: get/set the number of samples used to calculate
            quantitation()
  Returns : An integer
  Args    : a new integer (optional)

=cut

sub sample_count {
  shift->throw_not_implemented();
}

1;
