# $Id$
#
# bioperl module for Bio::SeqFeature::Tools::Unflattener
#
# Cared for by Chris Mungall <cjm@fruitfly.org>
#
# Copyright Chris Mungall
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Tools::Unflattener - nests a flat list of genbank-sourced features

=head1 SYNOPSIS

  # standard / generic use - unflatten a genbank record
  use Bio::SeqIO;
  use Bio::SeqFeature::Tools::Unflattener;

  # generate an Unflattener object
  $unflattener = Bio::SeqFeature::Tools::Unflattener->new;

  # first fetch a genbank SeqI object
  $seqio =
    Bio::SeqIO->new(-file=>'AE003644.gbk',
                    -format=>'GenBank');
  while ($seq = $seqio->next_seq()) {

    # get top level unflattended SeqFeatureI objects
    $unflattener->unflatten_seq(-seq=>$seq,
                                -use_magic=>1);
    @top_sfs = $seq->get_SeqFeatures;
  }


=head1 DESCRIPTION

Most GenBank entries for annotated genomic DNA contain a B<flat> list
of features. These features can be parsed into an equivalent flat list
of L<Bio::SeqFeatureI> objects using the standard L<Bio::SeqIO>
classes. However, it is often desirable to B<unflatten> this list into
something resembling actual B<gene models>, whereby genes, mRNAs and CDSs
are B<nested> according to the nature of the gene model.

The BioPerl object model allows us to store these kind of associations
between SeqFeatures in B<containment hierarchies> -- any SeqFeatureI
object can contain nested SeqFeatureI objects. The
Bio::SeqFeature::Tools::Unflattener object facilitates construction of
these hierarchies from the underlying GenBank flat-feature-list
representation.

For example, if you were to look at a typical GenBank DNA entry, say,
B<AE003644>, you would see a flat list of features:

  source
  gene
  mRNA CG4491-RA
  CDS CG4491-PA
  gene
  tRNA tRNA-Pro
  gene
  mRNA CG32954-RA
  mRNA CG32954-RC
  mRNA CG32954-RB
  CDS CG32954-PA
  CDS CG32954-PB
  CDS CG32954-PC

(shown as [type, product-name] pairs)

We would like to convert the above list into the B<containment
hierarchy>, shown below:

  source
  gene
    mRNA CG4491-RA
      CDS CG4491-PA
  gene
    tRNA tRNA-Pro
  gene
    mRNA CG32954-RA
      CDS CG32954-PA
    mRNA CG32954-RC
      CDS CG32954-PC
    mRNA CG32954-RB
      CDS CG32954-PB

We do this using a call on a Bio::SeqFeature::Tools::Unflattener object

  @sfs = $unflattener->unflatten_seq(-seq=>$seq);

This would return a list of the B<top level> (i.e. containing)
SeqFeatureI objects - in this case, genes. Other top level features
are possible; for instance, the B<source> feature which is always
present, and other features such as B<variation> or B<misc_feature>
types.

The containment hierarchy can be accessed using the get_SeqFeature()
call on any feature object - see L<Bio::SeqFeature::FeatureHolderI>.
The following code will traverse the containment hierarchy for a
feature:

  sub traverse {
    $sf = shift;   #  $sf isa Bio::SeqfeatureI
    
    # ...do something with $sf!

    # depth first traversal of containment tree
    @contained_sfs = $sf->get_SeqFeatures;
    traverse($_) foreach @contained_sfs;
  }

Once you have built the hierarchy, you can do neat stuff like turn the
features into 'rich' feature objects (eg
L<Bio::SeqFeature::Gene::GeneStructure>) or convert to a suitable
format such as GFF3 or chadoxml (after mapping to the Sequence
Ontology); this step is not described here.

=head1 USING MAGIC

Due to the quixotic nature of how features are stored in
GenBank/EMBL/DDBJ, there is no guarantee that the default behaviour of
this module will produce perfect results. Sometimes it is hard or
impossible to build a correct containment hierarchy if the information
provided is simply too lossy, as is often the case. If you care deeply
about your data, you should always manually inspect the resulting
containment hierarchy; you may have to customise the algorithm for
building the hierarchy, or even manually tweak the resulting
hierarchy. This is explained in more detail further on in the document.

However, if you are satisfied with the default behaviour, then you do
not need to read any further. Just make sure you set the parameter
B<use_magic> - this will invoke incantations which will magically
produce good results no matter what the idiosyncracies of the
particular GenBank record in question.

For example

  $unflattener->unflatten_seq(-seq=>$seq,
                              -use_magic=>1);

The success of this depends on the phase of the moon at the time the
entry was submitted to GenBank. Note that the magical recipe is being
constantly improved, so the results of invoking magic may vary
depending on the bioperl release.

If you are skeptical of magic, or you wish to exact fine grained
control over how the entry is unflattened, or you simply wish to
understand more about how this crazy stuff works, then read on!

=head1 ALGORITHM

This is the default algorithm; you should be able to override any part
of it to customise.

The core of the algorithm is in two parts

=over

=item Partitioning into groups

=item Resolving the containment hierarchy

=back

There are other optional steps after the completion of these two
steps; we now describe in more detail what is going on.

=head2 Partitioning into groups

First of all the flat feature list is partitioned into B<group>s.

The default way of doing this is to use the B<gene> attribute; if we
look at two features from GenBank accession AE003644:

     gene            20111..23268
                     /gene="noc"
                     /locus_tag="CG4491"
                     /note="last curated on Thu Dec 13 16:51:32 PST 2001"
                     /map="35B2-35B2"
                     /db_xref="FLYBASE:FBgn0005771"
     mRNA            join(20111..20584,20887..23268)
                     /gene="noc"
                     /locus_tag="CG4491"
                     /product="CG4491-RA"
                     /db_xref="FLYBASE:FBgn0005771"

Both these features share the same /gene tag which is "noc", so they
correspond to the same gene model (the CDS feature is not shown, but
this also has a tag-value /gene="noc").

Not all groups need to correspond to gene models, but this is the most
common use case; later on we shall describe how to customise this.

Sometimes other tags have to be used; for instance, if you look at the
entire record for AE003644 you will see you actually need the use the
/locus_tag attribute. This attribute is actually not present in most
records!

You can override this like this:

  $collection->unflatten_seq(-seq=>$seq, group_tag=>'locus_tag');

Alternatively, if you B<-use_magic>, the object will try and make a
guess as to what the correct group_tag should be.

=head2 Resolving the containment mapping

After the grouping is done, we end up with a list of groups which
probably contain features of type 'gene', 'mRNA', 'CDS'.

Each group is itself flat; we need to add an extra level of
organisation. Usually this is because different spliceforms
(represented by the 'mRNA' feature) can give rise to different
translations (represented by the 'CDS' feature). We want to correctly
associate mRNAs to CDSs.

We want to go from a group like this:

  [ gene mRNA mRNA mRNA CDS CDS CDS ]

to a containment hierarchy like this:

  gene
    mRNA
      CDS
    mRNA
      CDS
    mRNA
      CDS

In which each CDS corresponds to its containing mRNA.

How can we do this? The bad news is that there is no guaranteed way of
doing this correctly for all of GenBank. Occasionally the submission
will have been done in such a way as to reconstruct the containment
hierarchy. However, this is not consistent across databank entries, so
no generic solution can be provided witin bioperl. This module does
provide the framework within which you can customise a solution for
the particular dataset you are interested in - see later.

The good news is that there is an inference we can do that should
produce pretty good results most of the time. It uses splice
coordinate data - this is the default behaviour of this module.

=head2 Using splice site coordinates to infer containment

If an mRNA is to be the container for a CDS, then the splice site
coordinates of the CDS must fit inside the splice site coordinates of
the mRNA.

Ambiguities can still arise, but the results produced should still be
reasonable and consistent at the sequence level. For example

  mRNA    XXX---XX--XXXXXX--XXXX         join(1..3,7..8,11..16,19..23)
  mRNA    XXX-------XXXXXX--XXXX         join(1..3,11..16,19..23)
  CDS                 XXXX--XX           join(13..16,19..20)
  CDS                 XXXX--XX           join(13..16,19..20)

[obviously the positions have been scaled down]

We cannot unambiguously match mRNA with CDS based on splice sites,
since both CDS share the splice site locations 16^17 and
18^19. However, the consequences of making a wrong match are probably
not that severe. Any annotation data attached to the first CDS is
probably identical to the seconds CDS, other than identifiers.

The default behaviour of this module is to make an arbitrary call
where it is ambiguous (the mapping will always be bijective).

[NOTE: not tested on EMBL data, which may not be bijective; ie two
mRNAs can share the same CDS??]

Of course, if you are dealing with an organism with no alternate
splicing, you have nothing to worry about here! There is no ambiguity
possible, so you will always get a tree that looks like this:

  gene
    mRNA
      CDS

=head1 ADVANCED

=head2 Customising the grouping of features

The default behaviour is suited mostly to building models of protein
coding genes and noncoding genes from genbank genomic DNA submissions.

You can change the tag used to partition the feature by passing in a
different group_tag argument - see the unflatten_seq() method

Other behaviour may be desirable. For example, even though SNPs
(features of type 'variation' in GenBank) are not actually part of the
gene model, it may be desirable to group SNPs that overlap or are
nearby gene models.

It should certainly be possible to extend this module to do
this. However, I have yet to do this!

In the meantime, you could write your own grouping subroutine, and
feed the results into unflatten_groups() [see the method documentation
below]


=head2 Customising the resolution of the containment hierarchy

Once the flat list of features has been partitioned into groups, the
method unflatten_group() is called on each group to build a tree.

The algorithm for doing this is described above; ambiguities are
resolved by using splice coordinates. As discussed, this can be
ambiguous.

Some submissions may contain information in tags/attributes that hint
as to the mapping that needs to be made between the features.

For example, with the Drosophila Melanogaster release 3 submission, we
see that CDS features in alternately spliced mRNAs have a form like
this:

     CDS             join(145588..145686,145752..146156,146227..146493)
                     /locus_tag="CG32954"
                     /note="CG32954 gene product from transcript CG32954-RA"
                                                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^
                     /codon_start=1
                     /product="CG32954-PA"
                     /protein_id="AAF53403.1"
                     /db_xref="GI:7298167"
                     /db_xref="FLYBASE:FBgn0052954"
                     /translation="MSFTLTNKNVIFVAGLGGIGLDTSKELLKRDLKNLVILDRIENP..."

Here the /note tag provides the clue we need to link CDS to mRNA
(highlighted with ^^^^). We just need to find the mRNA with the tag

  /product="CG32954-RA"

I have no idea how consistent this practice is across submissions; it
is consistent for the fruitfly genome submission.

We can customise the behaviour of unflatten_group() by providing our
own resolver method. This obviously requires a bit of extra
programming, but there is no way to get around this.

Here is an example of how to pass in your own resolver; this example
basically checks the parent (container) /product tag to see if it
matches the required string in the child (contained) /note tag.

       $unflattener->unflatten_seq(-seq=>$seq,
                                 -group_tag=>'locus_tag',
                                 -resolver_method=>sub {
                                     my $self = shift;
                                     my ($sf, @candidate_container_sfs) = @_;
                                     if ($sf->has_tag('note')) {
                                         my @notes = $sf->get_tag_values('note');
                                         my @trnames = map {/from transcript\s+(.*)/;
                                                            $1} @notes;
                                         @trnames = grep {$_} @trnames;
                                         my $trname;
                                         if (@trnames == 0) {
                                             $self->throw("UNRESOLVABLE");
                                         }
                                         elsif (@trnames == 1) {
                                             $trname = $trnames[0];
                                         }
                                         else {
                                             $self->throw("AMBIGUOUS: @trnames");
                                         }
                                         my @container_sfs =
                                           grep {
                                               my ($product) =
                                                 $_->has_tag('product') ?
                                                   $_->get_tag_values('product') :
                                                     ('');
                                               $product eq $trname;
                                           } @candidate_container_sfs;
                                         if (@container_sfs == 0) {
                                             $self->throw("UNRESOLVABLE");
                                         }
                                         elsif (@container_sfs == 1) {
                                             # we got it!
                                             return $container_sfs[0];
                                         }
                                         else {
                                             $self->throw("AMBIGUOUS");
                                         }
                                     }
                                 });

the resolver method is only called when there is more than one spliceform.

=head2 Parsing mRNA records

Some of the entries in sequence databanks are for mRNA sequences as
well as genomic DNA. We may want to build models from these too.

NOT YET DONE - IN PROGRESS!!!

Open question - what would these look like?

=head1 SEE ALSO

Feature table description

  http://www.ebi.ac.uk/embl/Documentation/FT_definitions/feature_table.htm


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

package Bio::SeqFeature::Tools::Unflattener;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root
use Bio::Root::Root;

@ISA = qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : $unflattener = Bio::SeqFeature::Tools::Unflattener->new();
 Function: constructor
 Example : 
 Returns : a new Bio::SeqFeature::Tools::Unflattener
 Args    : see below


=cut


sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my($seq, $group_tag) =
	$self->_rearrange([qw(SEQ
                              GROUP_TAG
			     )],
                          @args);

    $seq  && $self->seq($seq);
    $group_tag  && $self->group_tag($group_tag);
    return $self; # success - we hope!
}



=head2 seq

 Title   : seq
 Usage   : $obj->seq($newval)
 Function: 
 Example : 
 Returns : value of seq (a Bio::SeqI)
 Args    : on set, new value (a Bio::SeqI, optional)

The Bio::SeqI object should hold a flat list of Bio::SeqFeatureI
objects; this is the list that will be unflattened.

=cut

sub seq{
    my $self = shift;

    return $self->{'seq'} = shift if @_;
    return $self->{'seq'};
}

=head2 group_tag

 Title   : group_tag
 Usage   : $obj->group_tag($newval)
 Function: 
 Example : 
 Returns : value of group_tag (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

This is the tag that will be used to collect elements from the flat
feature list into groups; for instance, if we look at two typical
GenBank features:

     gene            20111..23268
                     /gene="noc"
                     /locus_tag="CG4491"
                     /note="last curated on Thu Dec 13 16:51:32 PST 2001"
                     /map="35B2-35B2"
                     /db_xref="FLYBASE:FBgn0005771"
     mRNA            join(20111..20584,20887..23268)
                     /gene="noc"
                     /locus_tag="CG4491"
                     /product="CG4491-RA"
                     /db_xref="FLYBASE:FBgn0005771"

We can see that these comprise the same gene model because they share
the same /gene attribute; we want to collect these together in groups.

Setting group_tag is optional. The default is to use 'gene'. In the
example above, we could also use /locus_tag

=cut

sub group_tag{
    my $self = shift;

    return $self->{'group_tag'} = shift if @_;
    return $self->{'group_tag'};
}

=head2 containment_hierarchy

 Title   : containment_hierarchy
 Usage   : $obj->containment_hierarchy({mRNA=>'gene', CDS=>'mRNA')
 Function: 
 Example : 
 Returns : value of containment_hierarchy (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

A hash representing the containment structure that the seq_feature
nesting should conform to; each key represents the contained (child)
type; each value represents the container (parent) type.

=cut

sub containment_hierarchy{
    my $self = shift;

    return $self->{'containment_hierarchy'} = shift if @_;
    return $self->{'containment_hierarchy'} || $self->_default_containment_hierarchy;
}

sub _default_containment_hierarchy{
    return {
            mRNA => 'gene',
            tRNA => 'gene',
            CDS => 'mRNA',
	    exon => 'mRNA',
	    intron => 'mRNA',
           };
}

=head2 structure_type

 Title   : structure_type
 Usage   : $obj->structure_type($newval)
 Function: 
 Example : 
 Returns : value of structure_type (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=over

=item Type 0 (DEFAULT)

typically contains

  source
  gene
  mRNA
  CDS

with this structure type, we want the seq_features to be nested like this

  gene
    mRNA
    CDS
      exon

exons and introns are implicit from the mRNA 'join' location

to get exons from the mRNAs, you will need this call (see below)

  $unflattener->feature_from_splitloc(-seq=>$seq);

=item Type 1

typically contains

  source
  gene
  CDS
  exon
  intron

there are no mRNA features

with this structure type, we want the seq_features to be nested like this

  gene
    CDS
      exon
      intron

exon and intron may or may not be present; they may be implicit from
the mRNA 'join' location



=back

=cut

sub structure_type{
    my $self = shift;

    return $self->{'structure_type'} = shift if @_;
    return $self->{'structure_type'};
}


=head2 get_container_type

 Title   : get_container_type
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_container_type{
   my ($self,$type) = @_;
   my @roots = $self->_get_containment_hierarchy_roots;
   if (grep {$_ eq $type} @roots) {
       # it is a root - no parents/containers
       return;
   }
   my $ch = $self->containment_hierarchy;
   my $ctype = $ch->{$type};
   if (!$ctype) {
       # asterix acts as a wild card
       $ctype = $ch->{'*'};
   }
   return $ctype;
}

sub _get_containment_hierarchy_roots {
    my $self = shift;
    my $ch = $self->containment_hierarchy;
    my @parents = values %$ch;
    # find parents that do not have parents themselves
    return grep {!$ch->{$_}} @parents;
}

=head2 unflatten_seq

 Title   : unflatten_seq
 Usage   : @sfs = $unflattener->unflatten_seq($seq);
 Function: turns a flat list of features into a list of holder features
 Example :
 Returns : list of Bio::SeqFeatureI objects
 Args    : see below

partitions a list of features then arranges them in a nested tree; see
above for full explanation.

note - the original Bio::SeqI object passed in will be modified

Arguments

  -seq   :          a Bio::SeqI object; must contain Bio::SeqFeatureI objects

  -resolver_method: a CODE reference
                    see the documentation above for an example of
                    a subroutine that can be used to resolve hierarchies
                    within groups.

                    this is optional - if nothing is supplied, a default
                    subroutine will be used (see below)

  -group_tag:       a string
                    [ see the group_tag() method ]
                    this overrides the default group_tag which is 'gene'

=cut

sub unflatten_seq{
   my ($self,@args) = @_;

    my($seq, $resolver_method, $group_tag, $containment_hierarchy, 
       $structure_type, $use_magic) =
	$self->_rearrange([qw(SEQ
                              RESOLVER_METHOD
                              GROUP_TAG
                              CONTAINMENT_HIERARCHY
			      STRUCTURE_TYPE
			      USE_MAGIC
			     )],
                          @args);

   # seq we want to unflatten
   $seq = $seq || $self->seq;


   # prevent bad argument combinations
   if ($containment_hierarchy &&
       defined($structure_type)) {
       $self->throw("You cannot set both -containment_hierarchy and -structure_type\n".
		    "(the former is implied by the latter)");
   }

   # remember the current value of containment_hierarchy, to reset later
   my $old_containment_hierarchy = $self->containment_hierarchy;
   $self->containment_hierarchy($containment_hierarchy) if defined $containment_hierarchy;

   # remember old structure_type
   my $old_structure_type = $self->structure_type;
   $self->structure_type($structure_type) if defined $structure_type;

   # if we are sourcing our data from genbank, all the
   # features should be flat (eq no sub_SeqFeatures)
   my @flat_seq_features = $seq->get_SeqFeatures;
   my @all_seq_features = $seq->get_all_SeqFeatures;

   # sanity checks
   if (@all_seq_features > @flat_seq_features) {
       $self->throw("It looks as if this sequence has already been unflattened");
   }
   if (@all_seq_features < @flat_seq_features) {
       $self->throw("ASSERTION ERROR: something is seriously wrong with your features");
   }

   # tag for ungrouping; usually /gene or /locus_tag
   #     for example:        /gene="foo"
   $group_tag = $group_tag || $self->group_tag;
   if ($use_magic) {
       # use magic to guess the group tag
       my @sfs_with_locus_tag =
	 grep {$_->has_tag("locus_tag")} @flat_seq_features;
       if (@sfs_with_locus_tag) {
	   if ($group_tag && $group_tag ne 'locus_tag') {
	       $self->throw("You have explicitly set group_tag to be '$group_tag'\n".
			    "However, I detect that some features use /locus_tag\n".
			    "I believe that this is the correct group_tag to use\n".
			    "You can resolve this by either NOT setting -group_tag\n".
			    "OR you can unset -use_magic to regain control");
	   }
	   # I know what's best for you.
	   # You want to be using /locus_tag=foo here
	   #
	   # see GenBank entry AE003677 (version 3) for an example
	   $group_tag = 'locus_tag';
       }
   }
   if (!$group_tag) {
       $group_tag = 'gene';
   }

   # ------------------------------
   # GROUP FEATURES using $group_tag
   #     collect features into unstructured groups
   # ------------------------------

   # -------------
   # we want to generate a list of groups;
   # each group is a list of SeqFeatures; this
   # group probably (but not necessarily)
   # corresponds to a gene model.
   #
   # this array will look something like this:
   # ([$f1], [$f2, $f3, $f4], ...., [$f97, $f98, $f99])
   #
   # there are also 'singleton' groups, with one member.
   # for instance, the 'source' feature is in a singleton group;
   # the same with others such as 'misc_feature'
   my @groups = ();
   # -------------

   # --------------------
   # we hope that the genbank record allows us to group by some grouping
   # tag.
   # for instance, most of the time a gene model can be grouped using
   # the gene tag - that is where you see
   #                    /gene="foo"
   # in a genbank record
   # --------------------
   
   # keep an index of groups by their
   # grouping tag
   my %group_by_tag = ();
   

   # iterate through all features, putting them into groups
   foreach my $sf (@flat_seq_features) {
       if (!$sf->has_tag($group_tag)) {
	   # SINGLETON
           # this is an ungroupable feature;
           # add it to a group of its own
           push(@groups, [$sf]);
       }
       else {
	   # NON-SINGLETON
           my @group_tagvals = $sf->get_tag_values($group_tag);
           if (@group_tagvals > 1) {
	       # sanity check:
               # currently something can only belong to one group
               $self->throw(">1 value for /$group_tag: @group_tagvals\n".
                            "At this time this module is not equipped to handle this");
           }
	   # get value of group tag
           my $gtv = shift @group_tagvals;
           $gtv || $self->throw("Empty /$group_tag vals not allowed!");

           # is this a new group?
           my $group = $group_by_tag{$gtv};
           if ($group) {
               # this group has been encountered before - add current
               # sf to the end of the group
               push(@$group, $sf);
           }
           else {
               # new group; add to index and create new group
               $group = [$sf];  # currently one member; probably more to come
               $group_by_tag{$gtv} = $group;
               push(@groups, $group);
           }
       }
   }
   
   foreach my $group (@groups) {
       $self->_remove_duplicates_from_group($group);
   }

   # LOGGING
   if ($self->verbose) {
       print "GROUPS:\n";
       foreach my $group (@groups) {
	   $self->_write_group($group, $group_tag);
       }
   }
   # -

   # --- MAGIC ---
   my $need_to_infer_exons = 0;
   my $need_to_infer_mRNAs = 0;
   if ($use_magic) {
       if (defined($structure_type)) {
	   $self->throw("Can't combine use_magic AND setting structure_type");
       }
       foreach my $group (@groups) {
	   my @introns = grep {$_->primary_tag eq 'exon'} @$group;
	   my @exons = grep {$_->primary_tag eq 'exon'} @$group;
	   my @mrnas = grep {$_->primary_tag eq 'mRNA'} @$group;
	   my @cdss = grep {$_->primary_tag eq 'CDS'} @$group;
	   
	   if (@cdss) {
	       # a pc gene model should contain at the least a CDS
	       if (!@mrnas) {
		   # looks like structure_type == 1
		   $structure_type = 1;
		   $need_to_infer_mRNAs = 1;
	       }
	       # we always infer exons in magic mode
	       $need_to_infer_exons = 1;
	   }
	   else {
	       # this doesn't seem to be any kind of protein coding gene model
	   }
       }
   }
   # --- END OF MAGIC ---

   # >>>>>>>>>                   <<<<<<<<<<<<<
   # --------- FINISHED GROUPING -------------
   # >>>>>>>>>                   <<<<<<<<<<<<<

   # set the containment hierarchy if desired
   # see docs for structure_type() method
   if ($structure_type) {
       if ($structure_type == 1) {
	   $self->containment_hierarchy(
					{CDS => 'gene',
					 exon => 'CDS',
					 intron => 'CDS',
					}
				       );
       }
       else {
	   $self->throw("structure_type == $structure_type is currently unknown");
       }
   }


   # we have done the first part of the unflattening.
   # we now have a list of groups; each group is a list of seqfeatures.
   # the actual group itself is flat; we may want to unflatten this further;
   # for instance, a gene model can contain multiple mRNAs and CDSs. We may want
   # to link the correct mRNA to the correct CDS via the bioperl sub_SeqFeature tree.
   #
   # what we would end up with would be
   #  gene1
   #    mRNA-a
   #      CDS-a
   #    mRNA-b
   #      CDS-b
   my @top_sfs = $self->unflatten_groups(-groups=>\@groups,
                                         -resolver_method=>$resolver_method);
   
   # restore settings
   $self->containment_hierarchy($old_containment_hierarchy);

   # restore settings
   $self->structure_type($old_structure_type);

   # modify the original Seq object - the top seqfeatures are now
   # the top features from each group
   $seq->remove_SeqFeatures;
   $seq->add_SeqFeature(@top_sfs);

   # extra stuff
   if ($need_to_infer_mRNAs) {
       $self->infer_mRNA_from_CDS(-seq=>$seq);
   }
   if ($need_to_infer_exons) {
       # get rid of existing exons and introns - they
       # are redundant
       $self->remove_exons_and_introns(-seq=>$seq);
       
       # infer exons
       $self->feature_from_splitloc(-seq=>$seq);
   }

   # return new top level features; this can also 
   # be retrieved via
   #   $seq->get_SeqFeatures();
   return @top_sfs;
}

sub _remove_duplicates_from_group {
    my $self = shift;
    my $group = shift;

    # ::: WEIRD BOUNDARY CASE CODE :::
    # for some reason, there are some gb records with two gene
    # features for one gene; for example, see ATF14F8.gbk
    # in the t/data directory
    #
    # in this case, we get rid of one of the genes

    my @genes = grep {$_->primary_tag eq 'gene'} @$group;
    if (@genes > 1) {
	# OK, if we look at ATF14F8.gbk we see that some genes
	# just exist as a single location, some exist as a multisplit location;
	#
	# eg

	#     gene            16790..26395
	#                     /gene="F14F8_60"
	#     ...
	#     gene            complement(join(16790..19855,20136..20912,21378..21497,
	#                     21654..21876,22204..22400,22527..23158,23335..23448,
	#                     23538..23938,24175..24536,24604..24715,24889..24984,
	#                     25114..25171,25257..25329,25544..25589,25900..26018,
	#                     26300..26395))
	#                     /gene="F14F8_60"

	# the former is the 'standard' way of representing the gene in genbank;
	# the latter is redundant with the CDS entry. So we shall get rid of
	# the latter with the following filter

	@genes =
	  grep {
	      my $loc = $_->location;
	      if ($loc->isa("Bio::Location::SplitLocationI")) {
		  my @locs = $loc->each_Location;		  
		  if (@locs > 1) {
		      0;
		  }
		  else {
		      1;
		  }
	      }
	      else {
		  1;
	      }
	  } @genes;

	if (@genes > 1) {
	    # OK, that didn't work. Our only resort is to just pick one at random
	    @genes = ($genes[0]);
	}
	if (@genes) {
	    @genes == 1 || $self->throw("ASSERTION ERROR");
	    @$group =
	      ($genes[0], grep {$_->primary_tag ne 'gene'} @$group);
	}
    }
    # its a dirty job but someone's gotta do it
    return;
}


=head2 unflatten_groups

 Title   : unflatten_groups
 Usage   :
 Function: nests groups of features in containment hierarchies
 Example :
 Returns : list of Bio::SeqFeatureI objects that are holders
 Args    : see below

Arguments

  -groups:          list of list references; inner list is of Bio::SeqFeatureI objects
                    e.g.  ( [$sf1], [$sf2, $sf3, $sf4], [$sf5, ...], ...)

  -resolver_method: a CODE reference
                    see the documentation above for an example of
                    a subroutine that can be used to resolve hierarchies
                    within groups.

                    this is optional - a default subroutine will be used


NOTE: You should not need to call this method, unless you want fine
grained control over how the unflattening has been performed

=cut

sub unflatten_groups{
   my ($self,@args) = @_;
   my($groups, $resolver_method) =
     $self->_rearrange([qw(GROUPS
                           RESOLVER_METHOD
                          )],
                          @args);

   # this is just a simple wrapper for unflatten_group()
   return 
     map {
         $self->unflatten_group(-group=>$_, 
                                -resolver_method=>$resolver_method)
     } @$groups;
}

=head2 unflatten_group

 Title   : unflatten_group
 Usage   :
 Function: nests a group of features
 Example :
 Returns : Bio::SeqFeatureI objects that holds other features
 Args    : see below

Arguments

  -group:           reference to list of Bio::SeqFeatureI objects

  -resolver_method: a CODE reference
                    see the documentation above for an example of
                    a subroutine that can be used to resolve hierarchies
                    within groups

                    this is optional - a default subroutine will be used


You should not need to call this method, unless you want fine grained
control over how the unflattening has been performed

=cut

sub unflatten_group{
   my ($self,@args) = @_;

   my($group, $resolver_method) =
     $self->_rearrange([qw(GROUP
                           RESOLVER_METHOD
                          )],
                          @args);

   my @sfs = @$group;
   my $containment_hierarchy = $self->containment_hierarchy;

   # use default function for resolving, unless one passed in
   $resolver_method = $resolver_method || \&_resolve_container_for_sf;

   # find all the features for which there is no
   # containing feature type (eg genes)
   my @top_sfs =
     grep { 
         !$self->get_container_type($_->primary_tag);
     } @sfs;

   if (@top_sfs > 1) {
       $self->_write_group($group, $self->group_tag);
       $self->throw("multiple sfs @top_sfs in group");
   }


   my %sfs_by_type = ();
   foreach my $sf (@sfs) {
       push(@{$sfs_by_type{$sf->primary_tag}}, $sf);
   }

   my %container = ();
   my %has_been_chosen = ();

   # ALGORITHM: build containment graph
   #
   # find all possible containers for each SF;
   # for instance, for a CDS, the possible containers are all
   # the mRNAs in the same group. For a mRNA, the possible
   # containers are any SFs of type 'gene' (should only be 1).
   #
   # contention is resolved by checking coordinates of splice sites
   foreach my $sf (@sfs) {
       my $type = $sf->primary_tag;
       my $container_type = 
         $self->get_container_type($type);
       if ($container_type) {
           my @possible_container_sfs =
             @{$sfs_by_type{$container_type} || []};
           @possible_container_sfs =
             grep {!$has_been_chosen{$_}} @possible_container_sfs; 
           if (@possible_container_sfs == 0) {
           }
           elsif (@possible_container_sfs == 1) {
#               push(@{$sf_parents{$sf}}, $possible_container_sfs[0]);
               $container{$sf} = $possible_container_sfs[0];
           }
           else {
               # contention for container; use resolver
               # (default is to check splice sites)
               my @container_sfs =
                 $resolver_method->($self, $sf, @possible_container_sfs);
               if (@container_sfs == 0) {
                   # none found
                   printf "SF:%s\n\n", $sf->gff_string;
                   printf "POSS:%s\n",
                     join("\n\n", map{$_->gff_string} @possible_container_sfs);
                   $self->throw("No container for SF");
               }
               elsif (@container_sfs == 1) {
                   # perfect - exactly one possible container
                   $container{$sf} = $container_sfs[0];
               }
               else {
                   # ambiguous - loss has occurred in the GenBank
                   # representation
                   my $c = shift @container_sfs;
                   $has_been_chosen{$c} = 1;
                   $container{$sf} = $c;                   
               }
           }
       }
   }

   my @top = ();
   foreach my $sf (@sfs) {
       my $container_sf = $container{$sf};
       if ($container_sf) {
           $container_sf->add_SeqFeature($sf);
       }
       else {
           push(@top, $sf);
       }
   }
   return @top;

}

# ----------------------------------------------
# writes a group to stdout
#
# mostly for logging/debugging
# ----------------------------------------------
sub _write_group {
    my $self = shift;
    my $group = shift;
    my $group_tag = shift;

    my $f = $group->[0];
    my $label = '';
    if ($f->has_tag($group_tag)) {
	($label) = $f->get_tag_values($group_tag);
    }
    printf("  GROUP [%s]:%s\n",
	   $label,
	   join(' ',
		map { $_->primary_tag } @$group));

}

# -----------------------------------------------
#
# returns all possible containers for an SF based
# on splice site coordinates; splice site coords
# must be contained
# -----------------------------------------------
sub _resolve_container_for_sf{
   my ($self, $sf, @possible_container_sfs) = @_;

   my @coords = $self->_get_splice_coords_for_sf($sf);
   my $splice_uniq_str = "@coords";
   
   my @sfs =
     grep {
         my @coords = $self->_get_splice_coords_for_sf($_);
         index("@coords", $splice_uniq_str) > -1;
     } @possible_container_sfs;
   return @sfs;
}

sub _get_splice_coords_for_sf {
    my $self = shift;
    my $sf = shift;

   my @locs = $sf->location;
   if ($sf->location->isa("Bio::Location::SplitLocationI")) {
       @locs = $sf->location->each_Location;
   }

   # get an ordered list of (start, end) positions
   my @coords =
     map {
         $_->strand > 0 ? ($_->start, $_->end) : ($_->end, $_->start)
     } @locs;

   # remove first and last leaving only splice sites
   pop @coords;
   shift @coords;
    return @coords;
}

use Bio::Location::Simple;
use Bio::SeqFeature::Generic;

=head2 feature_from_splitloc

 Title   : feature_from_splitloc
 Usage   : $unflattener->feature_from_splitloc(-feature=>$sf);
 Function:
 Example :
 Returns : 
 Args    : see below

At this time all this method does is generate exons for mRNA ot tRNA features

Arguments:

  -feature:    a Bio::SeqFeatureI object (that conforms to Bio::FeatureHolderI)
  -seq:        a Bio::SeqI object that contains Bio::SeqFeatureI objects


=cut

sub feature_from_splitloc{
   my ($self,@args) = @_;

   my($sf, $seq) =
     $self->_rearrange([qw(FEATURE
                           SEQ
                          )],
                          @args);
   my @sfs = ($sf);
   if ($seq) {
       $seq->isa("Bio::SeqI") || $self->throw("$seq NOT A SeqI");
       @sfs = $seq->get_all_SeqFeatures;
   }
   foreach my $sf (@sfs) {

       $sf->isa("Bio::SeqFeatureI") || $self->throw("$sf NOT A SeqFeatureI");
       $sf->isa("Bio::FeatureHolderI") || $self->throw("$sf NOT A FeatureHolderI");

       my $type = $sf->primary_tag;
       next unless $type eq 'mRNA' or $type eq 'tRNA';

       my @locs = $sf->location;
       if ($sf->location->isa("Bio::Location::SplitLocationI")) {
           @locs = $sf->location->each_Location;
       }
       
       my @subsfs =
         map {
             my $subsf = Bio::SeqFeature::Generic->new(-location=>$_,
                                                       -primary_tag=>'exon');
             $subsf;
         } @locs;
       
       $sf->location(Bio::Location::Simple->new());
       $sf->add_SeqFeature($_, 'EXPAND') foreach @subsfs;
       if (!$sf->location->strand) {
	   # correct weird bioperl bug
	   $sf->location->strand($subsfs[0]->location->strand);
       }

       
   }
   return;
}

=head2 infer_mRNA_from_CDS

 Title   : infer_mRNA_from_CDS
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

given a "type 1" containment hierarchy

  gene
    CDS
      exon

this will infer the uniform "type 0" containment hierarchy

  gene
    mRNA
      CDS
      exon

all the children of the CDS will be moved to the mRNA

=cut

sub infer_mRNA_from_CDS{
   my ($self,@args) = @_;

   my($sf, $seq) =
     $self->_rearrange([qw(FEATURE
                           SEQ
                          )],
                          @args);
   my @sfs = ($sf);
   if ($seq) {
       $seq->isa("Bio::SeqI") || $self->throw("$seq NOT A SeqI");
       @sfs = $seq->get_all_SeqFeatures;
   }
   foreach my $sf (@sfs) {

       $sf->isa("Bio::SeqFeatureI") || $self->throw("$sf NOT A SeqFeatureI");
       $sf->isa("Bio::FeatureHolderI") || $self->throw("$sf NOT A FeatureHolderI");

       if ($sf->primary_tag eq 'mRNA') {
	   $self->throw("you cannot infer mRNAs if there are already mRNAs present");
       }

       my @cdsl = grep {$_->primary_tag eq 'CDS' } $sf->get_SeqFeatures;
       if (@cdsl) {
	   my @children = grep {$_->primary_tag ne 'CDS'} $sf->get_SeqFeatures;
	   my @mrnas = ();
	   foreach my $cds (@cdsl) {
	       
	       my $loc = Bio::Location::Split->new;
	       foreach my $cdsexonloc ($cds->location->each_Location) {
		   my $subloc =
		     Bio::Location::Simple->new(-start=>$cdsexonloc->start,
						-end=>$cdsexonloc->end,
						-strand=>$cdsexonloc->strand);
		   $loc->add_sub_Location($subloc);
	       }
	       # share the same location
	       my $mrna =
		 Bio::SeqFeature::Generic->new(-location=>$loc,
					       -primary_tag=>'mRNA');
	       
	       $mrna->add_SeqFeature($cds);
	       # mRNA steals children of CDS
	       foreach my $subsf ($cds->get_SeqFeatures) {
		   $mrna->add_SeqFeature($subsf);
	       }
	       $cds->remove_SeqFeatures;
	       push(@mrnas, $mrna);
	   }
	   # change gene/CDS to gene/mRNA
	   $sf->remove_SeqFeatures;
	   $sf->add_SeqFeature($_) foreach (@mrnas, @children);
       }
   }
   return;
   

}

sub remove_exons_and_introns {
   my ($self,@args) = @_;

   my($sf, $seq) =
     $self->_rearrange([qw(FEATURE
                           SEQ
                          )],
                          @args);
   my @sfs = ($sf);
   if ($seq) {
       $seq->isa("Bio::SeqI") || $self->throw("$seq NOT A SeqI");
       @sfs = $seq->get_all_SeqFeatures;
   }
   foreach my $sf (@sfs) {

       $sf->isa("Bio::SeqFeatureI") || $self->throw("$sf NOT A SeqFeatureI");
       $sf->isa("Bio::FeatureHolderI") || $self->throw("$sf NOT A FeatureHolderI");

       my %bad = (exon=>1, intron=>1);
       my @togo = grep {$bad{$_->primary_tag}} $sf->get_SeqFeatures;

       if (@togo) {
	   my @tostay = grep {!$bad{$_->primary_tag}} $sf->get_SeqFeatures;

	   $sf->remove_SeqFeatures;
	   $sf->add_SeqFeature($_) foreach (@tostay);
       }
   }
   return;
   

}



1;
