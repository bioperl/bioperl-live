# $Id$
#
# BioPerl module for Bio::FeatureHolderI
#
# Cared for by Hilmar Lapp <hlapp at gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::FeatureHolderI - the base interface an object with features must implement

=head1 SYNOPSIS

    use Bio::SeqIO;
    # get a feature-holding object somehow: for example, Bio::SeqI objects
    # have features
    my $seqio = Bio::SeqIO->new(-fh => \*STDIN, -format => 'genbank);
    while (my $seq = $seqio->next_seq()) {
        # $seq is-a Bio::FeatureHolderI, hence:
        my @feas = $seq->get_SeqFeatures();
        # each element is-a Bio::SeqFeatureI
        foreach my $fea (@feas) {
            # do something with the feature objects
        }
    }

=head1 DESCRIPTION

This is the base interface that all feature-holding objects must
implement.

Popular feature-holders are for instance L<Bio::Seq> objects. Since
L<Bio::SeqFeatureI> defines a sub_SeqFeature() method, most
Bio::SeqFeatureI implementations like L<Bio::SeqFeature::Generic> will
implement the feature holder interface as well.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Hilmar Lapp

Email hlapp at gmx.net

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...
#'


package Bio::FeatureHolderI;
use vars qw(@ISA);
use strict;
use Carp;
use Bio::Root::RootI;

@ISA = qw( Bio::Root::RootI );

=head2 get_SeqFeatures

 Title   : get_SeqFeatures
 Usage   :
 Function: Get the feature objects held by this feature holder.
 Example :
 Returns : an array of Bio::SeqFeatureI implementing objects
 Args    : none

At some day we may want to expand this method to allow for a feature
filter to be passed in.

=cut

sub get_SeqFeatures{
    shift->throw_not_implemented();
}

=head2 feature_count

 Title   : feature_count
 Usage   : $obj->feature_count()
 Function: Return the number of SeqFeatures attached to a feature holder.

           This is before flattening a possible sub-feature tree.

           We provide a default implementation here that just counts
           the number of objects returned by get_SeqFeatures().
           Implementors may want to override this with a more
           efficient implementation.

 Returns : integer representing the number of SeqFeatures
 Args    : None

At some day we may want to expand this method to allow for a feature
filter to be passed in.

Our default implementation allows for any number of additional
arguments and will pass them on to get_SeqFeatures(). I.e., in order to
support filter arguments, just support them in get_SeqFeatures().

=cut

sub feature_count {
    return scalar(shift->get_SeqFeatures(@_));
}

=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   :
 Function: Get the flattened tree of feature objects held by this
           feature holder. The difference to get_SeqFeatures is that
           the entire tree of sub-features will be flattened out.

           We provide a default implementation here, so implementors
           don''t necessarily need to implement this method.

 Example :
 Returns : an array of Bio::SeqFeatureI implementing objects
 Args    : none

At some day we may want to expand this method to allow for a feature
filter to be passed in.

Our default implementation allows for any number of additional
arguments and will pass them on to any invocation of
get_SeqFeatures(), wherever a component of the tree implements
FeatureHolderI. I.e., in order to support filter arguments, just
support them in get_SeqFeatures().

=cut

sub get_all_SeqFeatures{
    my $self = shift;
    my @flatarr;

    foreach my $feat ( $self->get_SeqFeatures(@_) ){
	push(@flatarr,$feat);
	&_add_flattened_SeqFeatures(\@flatarr,$feat,@_);
    }
    return @flatarr;
}

sub _add_flattened_SeqFeatures {
    my ($arrayref,$feat,@args) = @_;
    my @subs = ();

    if($feat->isa("Bio::FeatureHolderI")) {
	@subs = $feat->get_SeqFeatures(@args);
    } elsif($feat->isa("Bio::SeqFeatureI")) {
	@subs = $feat->sub_SeqFeature();
    } else {
	confess ref($feat)." is neither a FeatureHolderI nor a SeqFeatureI. ".
	    "Don't know how to flatten.";
    }
    foreach my $sub (@subs) {
	push(@$arrayref,$sub);
	&_add_flattened_SeqFeatures($arrayref,$sub);
    }

}

=head2 set_ParentIDs_from_hierarchy()

 Title   : set_ParentIDs_from_hierarchy()
 Usage   : $seq->set_ParentIDs_from_hierarchy()
 Function: populates tags Parent and ID via holder hierarchy
 Example :
 Returns : 
 Args    :

This is mainly for GFF3 export

GFF3 uses the tags ID and Parent to represent the feature containment
hierarchy; it does NOT use the feature holder tree

This method sets Parent (and ID for any parents not set) based on
feature holder/containement hierarchy, ready for GFF3 output

=cut

# method author: cjm@fruitfly.org
sub set_ParentIDs_from_hierarchy(){
   my ($self) = @_;

   # we will traverse the tree of contained seqfeatures
   # (a seqfeature is itself a holder)

   # start with the top-level features
   my @sfs = $self->get_SeqFeatures;

   # clear existing parent tags
   # (we assume this is the desired behaviour)
   my @all_sfs = $self->get_all_SeqFeatures;
   foreach (@all_sfs) {
       if ($_->has_tag('Parent')) {
           $_->remove_tag('Parent');
       }
   }
   

   # iterate until entire tree traversed
   while (@sfs) {
       my $sf = shift @sfs;
       my @subsfs = $sf->get_SeqFeatures;

       # see if the ID tag 
       my $id = $sf->primary_id;
       if (!$id) {
           # the skolem function feature(seq,start,end,type)
           # is presumed to uniquely identify this feature, and
           # to also be persistent
           $id = $sf->generate_unique_persistent_id;
       }
       foreach my $subsf (@subsfs) {
           $subsf->add_tag_value('Parent', $id);
       }
       
       # push children on to end of stack (breadth first search)
       push(@sfs, @subsfs);
   }
   return;
}

=head2 create_hierarchy_from_ParentIDs

 Title   : create_hierarchy_from_ParentIDs
 Usage   :
 Function:
 Example :
 Returns : list of top SeqFeatures
 Args    :

inverse of set_ParentIDs_from_hierarchy



=cut

sub create_hierarchy_from_ParentIDs{
   my ($self,@args) = @_;

   my @sfs = $self->get_all_SeqFeatures;
   my %sf_by_ID = ();
   foreach (@sfs) {
       my $id = $_->primary_id;
       next unless $id;
       if ($sf_by_ID{$id}) {
           $self->throw("DUPLICATE ID: $id");
       }
       $sf_by_ID{$id} = $_;
       $_->remove_SeqFeatures; # clear existing hierarchy (assume this is desired)
   }
   if (!%sf_by_ID) {
       # warn??
       # this is actually expected behaviour for some kinds of data;
       # eg lists of STSs - no containment hierarchy
       return;
   }

   my @topsfs = 
     grep {
         my @parents = $_->get_tagset_values('Parent');
         foreach my $parent (@parents) {
             $parent->add_SeqFeature($_);
         }
         !@parents;
     } @sfs;
   $self->remove_SeqFeatures;
   $self->add_SeqFeature($_) foreach @topsfs;
   return @topsfs;
}



1;
