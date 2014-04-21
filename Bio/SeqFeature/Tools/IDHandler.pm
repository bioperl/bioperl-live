#
# bioperl module for Bio::SeqFeature::Tools::IDHandler
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

Bio::SeqFeature::Tools::IDHandler - maps $seq_feature-E<gt>primary_tag

=head1 SYNOPSIS

  use Bio::SeqIO;
  use Bio::SeqFeature::Tools::IDHandler;


=head1 DESCRIPTION

Class to map $seq_feature-E<gt>primary_tag


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

Email:  cjm@fruitfly.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::SeqFeature::Tools::IDHandler;
use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : $unflattener = Bio::SeqFeature::Tools::IDHandler->new();
 Function: constructor
 Example : 
 Returns : a new Bio::SeqFeature::Tools::IDHandler
 Args    : see below


=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my($generate_id_sub) =
	$self->_rearrange([qw(GENERATE_ID_SUB
			     )],
                          @args);

    return $self; # success - we hope!
}

=head2 set_ParentIDs_from_hierarchy()

 Title   : set_ParentIDs_from_hierarchy()
 Usage   : $idhandler->set_ParentIDs_from_hierarchy($fholder)
 Function: populates tags Parent and ID via holder hierarchy
 Example :
 Returns : 
 Args    : Bio::featureHolderI (either a SeqFeature or a Seq)

This is mainly for GFF3 export

GFF3 uses the tags ID and Parent to represent the feature containment
hierarchy; it does NOT use the feature holder tree

This method sets Parent (and ID for any parents not set) based on
feature holder/containement hierarchy, ready for GFF3 output

=cut

# method author: cjm@fruitfly.org
sub set_ParentIDs_from_hierarchy(){
   my $self = shift;
   my ($featholder) = @_;

   # we will traverse the tree of contained seqfeatures
   # (a seqfeature is itself a holder)

   # start with the top-level features
   my @sfs = $featholder->get_SeqFeatures;

   # clear existing parent tags
   # (we assume this is the desired behaviour)
   my @all_sfs = $featholder->get_all_SeqFeatures;
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
 Usage   : $idhandler->set_ParentIDs_from_hierarchy($fholder)
 Function: inverse of set_ParentIDs_from_hierarchy
 Example :
 Returns : list of top SeqFeatures
 Args    :


=cut

sub create_hierarchy_from_ParentIDs{
   my ($self,$featholder,@args) = @_;

   my @sfs = $featholder->get_all_SeqFeatures;
   my %sf_by_ID = ();
   foreach (@sfs) {
       my $id = $_->primary_id;
       next unless $id;
       if ($sf_by_ID{$id}) {
           $featholder->throw("DUPLICATE ID: $id");
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
             $sf_by_ID{$parent}->add_SeqFeature($_)
		 if exists $sf_by_ID{$parent};
         }
         !@parents;
     } @sfs;
   $featholder->remove_SeqFeatures;
   $featholder->add_SeqFeature($_) foreach @topsfs;
   return @topsfs;
}


=head2 generate_unique_persistent_id

 Title   : generate_unique_persistent_id
 Usage   :
 Function: generates a unique and persistent identifier for this
 Example :
 Returns : value of primary_id (a scalar)
 Args    :

Will generate an ID, B<and> set primary_id() (see above)

The ID is a string generated from 

  seq_id
  primary_tag
  start
  end

There are three underlying assumptions: that all the above accessors
are set; that seq_id is a persistent and unique identifier for the
sequence containing this feature; and that 

  (seq_id, primary_tag, start, end) 

is a "unique constraint" over features

The ID is persistent, so long as none of these values change - if they
do, it is considered a separate entity

=cut

# method author: cjm@fruitfly.org
sub generate_unique_persistent_id{
   my ($self,$sf,@args) = @_;

   my $id;
   if (!$sf->isa("Bio::SeqFeatureI")) {
       $sf->throw("not a Bio::SeqFeatureI");
   }
   my $seq_id = $sf->seq_id || $sf->throw("seq_id must be set: ".$sf->display_name);
   #my $seq_id = $sf->seq_id || 'unknown_seq';
   if ($sf->has_tag('transcript_id')) {
       ($id) = $sf->get_tag_values('transcript_id');
   }
   elsif ($sf->has_tag('protein_id')) {
       ($id) = $sf->get_tag_values('protein_id');
   }
   else {
       my $source = $sf->source_tag || $sf->throw("source tag must be set: ".$sf->display_name);
       #my $source = $sf->source_tag || 'unknown_source';
       my $start = $sf->start || $sf->throw("start must be set or is zero: ".$sf->display_name);
       my $end = $sf->end || $sf->throw("end must be set");
       my $type = $sf->primary_tag || $sf->throw("primary_tag/type must be set: ".$sf->display_name);

       $id = "$source:$type:$seq_id:$start:$end";
   }
   $sf->primary_id($id);
   return $id;
}

1;
