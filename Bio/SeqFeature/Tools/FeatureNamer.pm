# $Id$
#
# bioperl module for Bio::SeqFeature::Tools::FeatureNamer
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


=head1 DESCRIPTION


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                         - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Chris Mungall

Email:  cjm@fruitfly.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::SeqFeature::Tools::FeatureNamer;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root
use Bio::Root::Root;

@ISA = qw(Bio::Root::Root);

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
 Args    :


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
 Args    :


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
 Args    : SeqFeatureI


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
