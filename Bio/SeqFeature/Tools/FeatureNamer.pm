#
# bioperl module for Bio::SeqFeature::Tools::FeatureNamer
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chris Mungall <cjm@fruitfly.org>
#
# Copyright Chris Mungall
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Tools::FeatureNamer - generates unique persistent names for features

=head1 SYNOPSIS

  use Bio::SeqIO;
  use Bio::SeqFeature::Tools::FeatureNamer;

  # first fetch a genbank SeqI object
  $seqio =
    Bio::SeqIO->new(-file=>'AE003644.gbk',
                    -format=>'GenBank');
  $seq = $seqio->next_seq();

  $namer = Bio::SeqFeature::Tools::FeatureNamer->new;
  my @features = $seq->get_SeqFeatures;
  foreach my $feature (@features) {
    $namer->name_feature($feature) unless $feature->display_name;
  }  

=head1 DESCRIPTION

This is a helper class for providing names for SeqFeatures

The L<Bio::SeqFeatureI> class provides a display_name
method. Typically the display_name is not set when parsing formats
such as genbank - instead properties such as B<label>, B<product> or
B<gene> are set in a somewhat inconsistent manner.

In addition, when generating subfeatures (for example, exons that are
subfeatures of a transcript feature), it is often desirable to name
these subfeatures before either exporting to another format or
reporting to the user.

This module is intended to help given uniform display_names to
features and their subfeatures.

=head1 TODO

Currently the naming policy is hardcoded. It may be desirable to allow
plugging in variations on naming policies; this could be done either
by subclassing, anonymous subroutines (closures) or
parameterization. Contact the author if you feel you have need for a
different naming policy


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

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

report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Chris Mungall

Email:  cjm AT fruitfly DOT org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::SeqFeature::Tools::FeatureNamer;
use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : $unflattener = Bio::SeqFeature::Tools::FeatureNamer->new();
 Function: constructor
 Example : 
 Returns : a new Bio::SeqFeature::Tools::FeatureNamer
 Args    : see below


=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

#    my($typemap) =
#	$self->_rearrange([qw(TYPEMAP
#			     )],
#                          @args);#

#    $typemap  && $self->typemap($typemap);
    return $self; # success - we hope!
}

=head2 name_feature

 Title   : name_feature
 Usage   : $namer->name_feature($sf);
 Function: sets display_name
 Example :
 Returns : 
 Args    : L<Bio::SeqFeatureI>

This method calls generate_feature_name() and uses the returned value
to set the display_name of the feature

=cut

sub name_feature {
    my ($self, $sf) = @_;
    my $name = $self->generate_feature_name($sf);
    $sf->display_name($name);
}

=head2 name_contained_features

 Title   : name_contained_features
 Usage   : $namer->name_contained_features($sf);
 Function: sets display_name for all features contained by sf
 Example :
 Returns : 
 Args    : L<Bio::SeqFeatureI>

iterates through all subfeatures of a certain feature (using
get_all_SeqFeatures) and names each subfeatures, based on the
generated name for the holder feature

A subfeature is named by concatenating the generated name of the
container feature with the type and a number.

For example, if the containing feature is a gene with display name
B<dpp>, subfeatures will be named dpp-mRNA-1 dpp-mRNA2 dpp-exon1
dpp-exon2 etc

=cut

sub name_contained_features{
   my ($self,$sf) = @_;
   my $cname = $self->generate_feature_name($sf);
   my @subsfs = $sf->get_all_SeqFeatures;
   my %num_by_type = ();
   foreach my $ssf (@subsfs) {
       my $type = $ssf->primary_tag;
       my $num = $num_by_type{$type} || 0;
       $num++;
       $num_by_type{$type} = $num;
       $ssf->display_name("$cname-$type-$num");
   }
   return;
}

=head2 generate_feature_name

 Title   : generate_feature_name
 Usage   : $name = $namer->generate_feature_name($sf);
 Function: derives a sensible human readable name for a $sf
 Example :
 Returns : str
 Args    : L<Bio::SeqFeatureI>

returns a generated name (but does not actually set display_name).

If display_name is already set, the method will return this

Otherwise, the name will depend on the property:

=over

=item label

=item product

=item gene

=item locus_tag

=back

(in order of priority)

=cut

sub generate_feature_name {
    my ($self, $sf) = @_;

    my $name = $sf->display_name;
    if (!$name) {
	if ($sf->has_tag("label")) {
	    ($name) = $sf->get_tag_values("label");
	}
	elsif ($sf->has_tag("product")) {
	    ($name) = $sf->get_tag_values("product");
	}
	elsif ($sf->primary_tag eq 'gene' &&
	       $sf->has_tag("gene")) {
	    ($name) = $sf->get_tag_values("gene");
	}
	elsif ($sf->primary_tag eq 'gene' &&
	       $sf->has_tag("locus_tag")) {
	    ($name) = $sf->get_tag_values("locus_tag");
	}
	else {
	    $name =  $sf->display_name;
	}
    }
    return $name;
}

1;
