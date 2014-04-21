#
# bioperl module for Bio::SeqFeature::Tools::Unflattener
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

Bio::SeqFeature::Tools::Unflattener - turns flat list of genbank-sourced features into a nested SeqFeatureI hierarchy

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
  my $out =
    Bio::SeqIO->new(-format=>'asciitree');
  while ($seq = $seqio->next_seq()) {

    # get top level unflattended SeqFeatureI objects
    $unflattener->unflatten_seq(-seq=>$seq,
                                -use_magic=>1);
    $out->write_seq($seq);

    @top_sfs = $seq->get_SeqFeatures;
    foreach my $sf (@top_sfs) {
	# do something with top-level features (eg genes)
    }
  }


=head1 DESCRIPTION

Most GenBank entries for annotated genomic DNA contain a B<flat> list
of features. These features can be parsed into an equivalent flat list
of L<Bio::SeqFeatureI> objects using the standard L<Bio::SeqIO>
classes. However, it is often desirable to B<unflatten> this list into
something resembling actual B<gene models>, in which genes, mRNAs and CDSs
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

  gene CG4491
  mRNA CG4491-RA
  CDS CG4491-PA

  gene tRNA-Pro
  tRNA tRNA-Pro

  gene CG32954
  mRNA CG32954-RA
  mRNA CG32954-RC
  mRNA CG32954-RB
  CDS CG32954-PA
  CDS CG32954-PB
  CDS CG32954-PC

These features have sequence locations, but it is not immediately
clear how to write code such that each mRNA is linked to the
appropriate CDS (other than relying on IDs which is very bad)

We would like to convert the above list into the B<containment
hierarchy>, shown below:

  source
  gene
    mRNA CG4491-RA
      CDS CG4491-PA
      exon
      exon
  gene
    tRNA tRNA-Pro
      exon
  gene
    mRNA CG32954-RA
      CDS CG32954-PA
      exon
      exon
    mRNA CG32954-RC
      CDS CG32954-PC
      exon
      exon
    mRNA CG32954-RB
      CDS CG32954-PB
      exon
      exon

Where each feature is nested underneath its container. Note that exons
have been automatically inferred (even for tRNA genes).

We do this using a call on a L<Bio::SeqFeature::Tools::Unflattener>
object

  @sfs = $unflattener->unflatten_seq(-seq=>$seq);

This would return a list of the B<top level> (i.e. container)
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

=head1 PROBLEMATIC DATA AND INCONSISTENCIES

Occasionally the Unflattener will have problems with certain
records. For example, the record may contain inconsistent data - maybe
there is an B<exon> entry that has no corresponding B<mRNA> location. 

The default behaviour is to throw an exception reporting the problem,
if the problem is relatively serious - for example, inconsistent data.

You can exert more fine grained control over this - perhaps you want
the Unflattener to do the best it can, and report any problems. This
can be done - refer to the methods.

  error_threshold()

  get_problems()

  report_problems()

  ignore_problems()

=head1 ALGORITHM

This is the default algorithm; you should be able to override any part
of it to customise.

The core of the algorithm is in two parts

=over

=item Partitioning the flat feature list into groups

=item Resolving the feature containment hierarchy for each group

=back

There are other optional steps after the completion of these two
steps, such as B<inferring exons>; we now describe in more detail what
is going on.

=head2 Partitioning into groups

First of all the flat feature list is partitioned into B<group>s.

The default way of doing this is to use the B<gene> attribute; if we
look at two features from GenBank accession AE003644.3:

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
common use case; later on we shall describe how to customise the
grouping.

Sometimes other tags have to be used; for instance, if you look at the
entire record for AE003644.3 you will see you actually need the use the
/locus_tag attribute. This attribute is actually B<not present> in
most records!

You can override this:

  $collection->unflatten_seq(-seq=>$seq, -group_tag=>'locus_tag');

Alternatively, if you B<-use_magic>, the object will try and make a
guess as to what the correct group_tag should be.

At the end of this step, we should have a list of groups - there is no
structure within a group; the group just serves to partition the flat
features. For the example data above, we would have the following groups.

  [ source ]
  [ gene mRNA CDS ]
  [ gene mRNA CDS ]
  [ gene mRNA CDS ]
  [ gene mRNA mRNA mRNA CDS CDS CDS ]

=head3 Multicopy Genes

Multicopy genes are usually rRNAs or tRNAs that are duplicated across
the genome. Because they are functionally equivalent, and usually have
the same sequence, they usually have the same group_tag (ie gene
symbol); they often have a /note tag giving copy number. This means
they will end up in the same group. This is undesirable, because they
are spatially disconnected.

There is another step, which involves splitting spatially disconnected
groups into distinct groups

this would turn this

 [gene-rrn3 rRNA-rrn3 gene-rrn3 rRNA-rrn3]

into this

 [gene-rrn3 rRNA-rrn3] [gene-rrn3 rRNA-rrn3]

based on the coordinates

=head3 What next?

The next step is to add some structure to each group, by making
B<containment hierarchies>, trees that represent how the features
interrelate

=head2 Resolving the containment hierarchy

After the grouping is done, we end up with a list of groups which
probably contain features of type 'gene', 'mRNA', 'CDS' and so on.

Singleton groups (eg the 'source' feature) are ignored at this stage.

Each group is itself flat; we need to add an extra level of
organisation. Usually this is because different spliceforms
(represented by the 'mRNA' feature) can give rise to different
protein products (indicated by the 'CDS' feature). We want to correctly
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

In which each CDS is nested underneath the correct corresponding mRNA.

For entries that contain no alternate splicing, this is simple; we
know that the group

  [ gene mRNA CDS ]

Must resolve to the tree

  gene
    mRNA
      CDS

How can we do this in entries with alternate splicing? The bad
news is that there is no guaranteed way of doing this correctly for
any GenBank entry. Occasionally the submission will have been done in
such a way as to reconstruct the containment hierarchy. However, this
is not consistent across databank entries, so no generic solution can
be provided by this object. This module does provide the framework
within which you can customise a solution for the particular dataset
you are interested in - see later.

The good news is that there is an inference we can do that should
produce pretty good results the vast majority of the time. It uses
splice coordinate data - this is the default behaviour of this module,
and is described in detail below.

=head2 Using splice site coordinates to infer containment

If an mRNA is to be the container for a CDS, then the splice site
coordinates (or intron coordinates, depending on how you look at it)
of the CDS must fit inside the splice site coordinates of the mRNA.

Ambiguities can still arise, but the results produced should still be
reasonable and consistent at the sequence level. Look at this fake
example:

  mRNA    XXX---XX--XXXXXX--XXXX         join(1..3,7..8,11..16,19..23)
  mRNA    XXX-------XXXXXX--XXXX         join(1..3,11..16,19..23)
  CDS                 XXXX--XX           join(13..16,19..20)
  CDS                 XXXX--XX           join(13..16,19..20)

[obviously the positions have been scaled down]

We cannot unambiguously match mRNA with CDS based on splice sites,
since both CDS share the splice site locations 16^17 and
18^19. However, the consequences of making a wrong match are probably
not very severe. Any annotation data attached to the first CDS is
probably identical to the seconds CDS, other than identifiers.

The default behaviour of this module is to make an arbitrary call
where it is ambiguous (the mapping will always be bijective; i.e. one
mRNA -E<gt> one CDS).

[TODO: NOTE: not tested on EMBL data, which may not be bijective; ie two
mRNAs can share the same CDS??]

This completes the building of the containment hierarchy; other
optional step follow

=head1 POST-GROUPING STEPS

=head2 Inferring exons from mRNAs

This step always occurs if B<-use_magic> is invoked.

In a typical GenBank entry, the exons are B<implicit>. That is they
can be inferred from the mRNA location.

For example:

     mRNA            join(20111..20584,20887..23268)

This tells us that this particular transcript has two exons. In
bioperl, the mRNA feature will have a 'split location'.

If we call

  $unflattener->feature_from_splitloc(-seq=>$seq);

This will generate the necessary exon features, and nest them under
the appropriate mRNAs. Note that the mRNAs will no longer have split
locations - they will have simple locations spanning the extent of the
exons. This is intentional, to avoid redundancy.

Occasionally a GenBank entry will have both implicit exons (from the
mRNA location) B<and> explicit exon features.

In this case, exons will still be transferred. Tag-value data from the
explicit exon will be transfered to the implicit exon. If exons are
shared between mRNAs these will be represented by different
objects. Any inconsistencies between implicit and explicit will be
reported.

=head3 tRNAs and other noncoding RNAs

exons will also be generated from these features

=head2 Inferring mRNAs from CDS

Some GenBank entries represent gene models using features of type
gene, mRNA and CDS; some entries just use gene and CDS.

If we only have gene and CDS, then the containment hierarchies will
look like this:

  gene
    CDS

If we want the containment hierarchies to be uniform, like this

  gene
    mRNA
      CDS

Then we must create an mRNA feature. This will have identical
coordinates to the CDS. The assumption is that there is either no
untranslated region, or it is unknown.

To do this, we can call

   $unflattener->infer_mRNA_from_CDS(-seq=>$seq);

This is taken care of automatically, if B<-use_magic> is invoked.

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
this. However, I have yet to code this part!!! If anyone would find
this useful let me know.

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

Ideally we would like a way of combining a mRNA record with the
corresponding SeFeature entry from the appropriate genomic DNA
record. This could be problemmatic in some cases - for example, the
mRNA sequences may not match 100% (due to differences in strain,
assembly problems, sequencing problems, etc). What then...?

=head1 SEE ALSO

Feature table description

  http://www.ebi.ac.uk/embl/Documentation/FT_definitions/feature_table.html

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                         - General discussion
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

package Bio::SeqFeature::Tools::Unflattener;
use strict;

# Object preamble - inherits from Bio::Root::Root
use Bio::Location::Simple;
use Bio::SeqFeature::Generic;
use Bio::Range;


use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : $unflattener = Bio::SeqFeature::Tools::Unflattener->new();
           $unflattener->unflatten_seq(-seq=>$seq);
 Function: constructor
 Example : 
 Returns : a new Bio::SeqFeature::Tools::Unflattener
 Args    : see below

Arguments

  -seq       : A L<Bio::SeqI> object (optional)
               the sequence to unflatten; this can also be passed in
               when we call unflatten_seq()

  -group_tag : a string representing the /tag used to partition flat features
               (see discussion above)

=cut


sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my($seq, $group_tag, $trust_grouptag) =
	$self->_rearrange([qw(SEQ
                              GROUP_TAG
                              TRUST_GROUPTAG
			     )],
                          @args);

    $seq  && $self->seq($seq);
    $group_tag  && $self->group_tag($group_tag);
    # $self->{'trust_grouptag'}= $trust_grouptag if($trust_grouptag); #dgg suggestion
    return $self; # success - we hope!
}

sub DESTROY {
    my $self = shift;
    return if $self->{_reported_problems};
    return if $self->{_ignore_problems};
    my @probs = $self->get_problems;
    if (!$self->{_problems_reported} &&
	scalar(@probs)) {
	$self->warn(
	    "WARNING: There are UNREPORTED PROBLEMS.\n".
	    "You may wish to use the method report_problems(), \n",
	    "or ignore_problems() on the Unflattener object\n");
    }
    return;
}

=head2 seq

 Title   : seq
 Usage   : $unflattener->seq($newval)
 Function: 
 Example : 
 Returns : value of seq (a Bio::SeqI)
 Args    : on set, new value (a Bio::SeqI, optional)

The Bio::SeqI object should hold a flat list of Bio::SeqFeatureI
objects; this is the list that will be unflattened.

The sequence object can also be set when we call unflatten_seq()

=cut

sub seq{
    my $self = shift;

    return $self->{'seq'} = shift if @_;
    return $self->{'seq'};
}

=head2 group_tag

 Title   : group_tag
 Usage   : $unflattener->group_tag($newval)
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

=head2 partonomy

 Title   : partonomy
 Usage   : $unflattener->partonomy({mRNA=>'gene', CDS=>'mRNA')
 Function: 
 Example : 
 Returns : value of partonomy (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

A hash representing the containment structure that the seq_feature
nesting should conform to; each key represents the contained (child)
type; each value represents the container (parent) type.

=cut

sub partonomy{
    my $self = shift;

    return $self->{'partonomy'} = shift if @_;
    if (!$self->{'partonomy'}) {
	$self->{'partonomy'} = $self->_default_partonomy;
    }
    return $self->{'partonomy'};
}

sub _default_partonomy{
    return {
            mRNA => 'gene',
            tRNA => 'gene',
            rRNA => 'gene',
            scRNA => 'gene',
            snRNA => 'gene',
            snoRNA => 'gene',
            misc_RNA => 'gene',
            CDS => 'mRNA',
	    exon => 'mRNA',
	    intron => 'mRNA',

            pseudoexon => 'pseudogene',
            pseudointron => 'pseudogene',
            pseudotranscript => 'pseudogene',
           };
}

=head2 structure_type

 Title   : structure_type
 Usage   : $unflattener->structure_type($newval)
 Function: 
 Example : 
 Returns : value of structure_type (a scalar)
 Args    : on set, new value (an int or undef, optional)

GenBank entries conform to different flavours, or B<structure
types>. Some have mRNAs, some do not.

Right now there are only two base structure types defined. If you set
the structure type, then appropriate unflattening action will be
taken.  The presence or absence of explicit exons does not affect the
structure type.

If you invoke B<-use_magic> then this will be set automatically, based
on the content of the record.

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
  exon [optional]
  intron [optional]

there are no mRNA features

with this structure type, we want the seq_features to be nested like this

  gene
    CDS
      exon
      intron

exon and intron may or may not be present; they may be implicit from
the CDS 'join' location

=back

=cut

sub structure_type{
    my $self = shift;

    return $self->{'structure_type'} = shift if @_;
    return $self->{'structure_type'};
}

=head2 get_problems

 Title   : get_problems
 Usage   : @probs = get_problems()
 Function: Get the list of problem(s) for this object.
 Example :
 Returns : An array of [severity, description] pairs
 Args    :

In the course of unflattening a record, problems may occur. Some of
these problems are non-fatal, and can be ignored.

Problems are represented as arrayrefs containing a pair [severity,
description]

severity is a number, the higher, the more severe the problem

the description is a text string

=cut

sub get_problems{
    my $self = shift;

    return @{$self->{'_problems'}} if exists($self->{'_problems'});
    return ();
}

=head2 clear_problems

 Title   : clear_problems
 Usage   :
 Function: resets the problem list to empty
 Example :
 Returns : 
 Args    :


=cut

sub clear_problems{
   my ($self,@args) = @_;
   $self->{'_problems'} = [];
   return;
}


# PRIVATE
# see get_problems
sub add_problem{
    my $self = shift;

    $self->{'_problems'} = [] unless exists($self->{'_problems'});
    if ($self->verbose > 0) {
        warn( "PROBLEM: $_\n") foreach @_;
    }
    push(@{$self->{'_problems'}}, @_);
}

# PRIVATE
# see get_problems
sub problem {
    my $self = shift;
    my ($severity, $desc, @sfs) = @_;
    if (@sfs) {
	foreach my $sf (@sfs) {
	    $desc .=
	      sprintf("\nSF [$sf]: ". $sf->location->to_FTstring . "; %s\n",
		      join('; ',
                           $sf->primary_tag,
			   map {
			       $sf->has_tag($_) ?
				 $sf->get_tag_values($_) : ()
			     } qw(locus_tag gene product label)));
	}
    }
    my $thresh = $self->error_threshold;
    if ($severity > $thresh) {
	$self->{_problems_reported} = 1;
	$self->throw("PROBLEM, SEVERITY==$severity\n$desc");
    }
    $self->add_problem([$severity, $desc]);
    return;
}

=head2 report_problems

 Title   : report_problems
 Usage   : $unflattener->report_problems(\*STDERR);
 Function:
 Example :
 Returns : 
 Args    : FileHandle (defaults to STDERR)


=cut

sub report_problems{
   my ($self, $fh) = @_;

   if (!$fh) {
       $fh = \*STDERR;
   }
   foreach my $problem ($self->get_problems) {
       my ($sev, $desc) = @$problem;
       printf $fh "PROBLEM, SEVERITY==$sev\n$desc\n";
   }
   $self->{_problems_reported} = 1;
   return;
}

=head2 ignore_problems

 Title   : ignore_problems
 Usage   : $obj->ignore_problems();
 Function:
 Example :
 Returns : 
 Args    :

Unflattener is very particular about problems it finds along the
way. If you have set the error_threshold such that less severe
problems do not cause exceptions, Unflattener still expects you to
report_problems() at the end, so that the user of the module is aware
of any inconsistencies or problems with the data. In fact, a warning
will be produced if there are unreported problems. To silence, this
warning, call the ignore_problems() method before the Unflattener
object is destroyed.

=cut

sub ignore_problems{
   my ($self) = @_;
   $self->{_ignore_problems} = 1;
   return;
}


=head2 error_threshold

 Title   : error_threshold
 Usage   : $obj->error_threshold($severity)
 Function: 
 Example : 
 Returns : value of error_threshold (a scalar)
 Args    : on set, new value (an integer)

Sets the threshold above which errors cause this module to throw an
exception. The default is 0; all problems with a severity E<gt> 0 will
cause an exception.

If you raise the threshold to 1, then the unflattening process will be
more lax; problems of severity==1 are generally non-fatal, but may
indicate that the results should be inspected, for example, to make
sure there is no data loss.

=cut

sub error_threshold{
    my $self = shift;

    return $self->{'error_threshold'} = shift if @_;
    return $self->{'error_threshold'} || 0;
}



# PRIVATE
#
# given a type (eg mRNA), will return the container type (eg gene)
sub get_container_type{
   my ($self,$type) = @_;
   my @roots = $self->_get_partonomy_roots;
   if (grep {$_ eq $type} @roots) {
       # it is a root - no parents/containers
       return;
   }
   my $ch = $self->partonomy;
   my $ctype = $ch->{$type};
   if (!$ctype) {
       # asterix acts as a wild card
       $ctype = $ch->{'*'};
   }
   return $ctype;
}

# get root node of partonomy hierarchy (usually gene)
sub _get_partonomy_roots {
    my $self = shift;
    my $ch = $self->partonomy;
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

note - the Bio::SeqI object passed in will be modified

Arguments

  -seq   :          a Bio::SeqI object; must contain Bio::SeqFeatureI objects
                    (this is optional if seq has already been set)

  -use_magic:       if TRUE (ie non-zero) then magic will be invoked;
                    see discussion above.

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

    my($seq, $resolver_method, $group_tag, $partonomy, 
       $structure_type, $resolver_tag, $use_magic, $noinfer) =
	$self->_rearrange([qw(SEQ
                              RESOLVER_METHOD
                              GROUP_TAG
                              PARTONOMY
			      STRUCTURE_TYPE
			      RESOLVER_TAG
			      USE_MAGIC
			      NOINFER
			     )],
                          @args);

   # seq we want to unflatten
   $seq = $seq || $self->seq;
   if (!$self->seq) {
       $self->seq($seq);
   }


   # prevent bad argument combinations
   if ($partonomy &&
       defined($structure_type)) {
       $self->throw("You cannot set both -partonomy and -structure_type\n".
		    "(the former is implied by the latter)");
   }

   # remember the current value of partonomy, to reset later
   my $old_partonomy = $self->partonomy;
   $self->partonomy($partonomy) if defined $partonomy;

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
       my @sfs_with_gene_tag =
	 grep {$_->has_tag("gene")} @flat_seq_features;
       my @sfs_with_product_tag =
	 grep {$_->has_tag("product")} @flat_seq_features;
	 
#        if ($group_tag && $self->{'trust_grouptag'}) { # dgg suggestion
# 
#         }
#        elsif
       if (@sfs_with_locus_tag) {
        # dgg note: would like to -use_magic with -group_tag = 'gene' for ensembl genomes
        # where ensembl gene FT have both /locus_tag and /gene, but mRNA, CDS have /gene only
	   if ($group_tag && $group_tag ne 'locus_tag') {
	       $self->throw("You have explicitly set group_tag to be '$group_tag'\n".
			    "However, I detect that some features use /locus_tag\n".
			    "I believe that this is the correct group_tag to use\n".
			    "You can resolve this by either NOT setting -group_tag\n".
			    "OR you can unset -use_magic to regain control");
	   }

	   # use /locus_tag instead of /gene tag for grouping
	   # see GenBank entry AE003677 (version 3) for an example
	   $group_tag = 'locus_tag';
           if ($self->verbose > 0) {
               warn "Set group tag to: $group_tag\n";
           }
       }

       # on rare occasions, records will have no /gene or /locus_tag
       # but it WILL have /product tags. These serve the same purpose
       # for grouping. For an example, see AY763288 (also in t/data)
       if (@sfs_with_locus_tag==0 &&
           @sfs_with_gene_tag==0 &&
           @sfs_with_product_tag>0 &&
           !$group_tag) {
	   $group_tag = 'product';
           if ($self->verbose > 0) {
               warn "Set group tag to: $group_tag\n";
           }
           
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
               $self->problem(2,
			      ">1 value for /$group_tag: @group_tagvals\n".
			      "At this time this module is not equipped to handle this adequately", $sf);
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
   
   # as well as having the same group_tag, a group should be spatially
   # connected. if not, then the group should be split into subgroups.
   # this turns out to be necessary in the case of multicopy genes.
   # the standard way to represent these is as spatially disconnected
   # gene models (usually a 'gene' feature and some kind of RNA feature)
   # with the same group tag; the code below will split these into 
   # seperate groups, one per copy.
   @groups = map { $self->_split_group_if_disconnected($_) } @groups;

   # remove any duplicates; most of the time the method below has
   # no effect. there are some unusual genbank records for which
   # duplicate removal is necessary. see the comments in the
   # _remove_duplicates_from_group() method if you want to know
   # the ugly details
   foreach my $group (@groups) {
       $self->_remove_duplicates_from_group($group);
   }

   # -

   # PSEUDOGENES, PSEUDOEXONS AND PSEUDOINTRONS
   # these are indicated with the /pseudo tag
   # these are mapped to a different type; they should NOT
   # be treated as normal genes
   foreach my $sf (@all_seq_features) {
       if ($sf->has_tag('pseudo')) {
           my $type = $sf->primary_tag;
           # SO type is typically the same as the normal
           # type but preceeded by "pseudo"
           if ($type eq 'misc_RNA' || $type eq 'mRNA') { 
            # dgg: see TypeMapper; both pseudo mRNA,misc_RNA should be pseudogenic_transcript
               $sf->primary_tag("pseudotranscript");
           }
           else {
               $sf->primary_tag("pseudo$type");
           }
       }
   }
   # now some of the post-processing that follows which applies to
   # genes will NOT be applied to pseudogenes; this is deliberate
   # for example, gene models are normalised to be gene-transcript-exon
   # for pseudogenes we leave them as pseudogene-pseudoexon

   # --- MAGIC ---
   my $need_to_infer_exons = 0;
   my $need_to_infer_mRNAs = 0;
   my @removed_exons = ();
   if ($use_magic) {
       if (defined($structure_type)) {
	   $self->throw("Can't combine use_magic AND setting structure_type");
       }
       my $n_introns =
	 scalar(grep {$_->primary_tag eq 'exon'} @flat_seq_features);
       my $n_exons =
	 scalar(grep {$_->primary_tag eq 'exon'} @flat_seq_features);
       my $n_mrnas =
	 scalar(grep {$_->primary_tag eq 'mRNA'} @flat_seq_features);
       my $n_mrnas_attached_to_gene =
	 scalar(grep {$_->primary_tag eq 'mRNA' &&
			$_->has_tag($group_tag)} @flat_seq_features);
       my $n_cdss =
	 scalar(grep {$_->primary_tag eq 'CDS'} @flat_seq_features);
       my $n_rnas =
	 scalar(grep {$_->primary_tag =~ /RNA/} @flat_seq_features);  
       # Are there any CDS features in the record?
       if ($n_cdss > 0) {
           # YES
           
	   # - a pc gene model should contain at the least a CDS

           # Are there any mRNA features in the record?
	   if ($n_mrnas == 0) {
               # NO mRNAs:
	       # looks like structure_type == 1
	       $structure_type = 1;
	       $need_to_infer_mRNAs = 1;
	   }
	   elsif ($n_mrnas_attached_to_gene == 0) {
               # $n_mrnas > 0
               # $n_mrnas_attached_to_gene = 0
               #
               # The entries _do_ contain mRNA features,
               # but none of them are part of a group/gene, i.e. they
               # are 'floating'

	       # this is an annoying weird file that has some floating
	       # mRNA features; 
	       # eg ftp.ncbi.nih.gov/genomes/Schizosaccharomyces_pombe/
               
               if ($self->verbose) {
                   my @floating_mrnas =
                     grep {$_->primary_tag eq 'mRNA' &&
                             !$_->has_tag($group_tag)} @flat_seq_features;
                   printf STDERR "Unattached mRNAs:\n";
                   foreach my $mrna (@floating_mrnas) {
                       $self->_write_sf_detail($mrna);
                   }
                   printf STDERR "Don't know how to deal with these; filter at source?\n";
               }

	       foreach (@flat_seq_features) {
		   if ($_->primary_tag eq 'mRNA') {
		       # what should we do??
		       
		       # I think for pombe we just have to filter
		       # out bogus mRNAs prior to starting
		   }
	       }

	       # looks like structure_type == 2
	       $structure_type = 2;
	       $need_to_infer_mRNAs = 1;
	   }
	   else {
	   }

	   # we always infer exons in magic mode
	   $need_to_infer_exons = 1;
       }
       else {
	   # this doesn't seem to be any kind of protein coding gene model
	   if ( $n_rnas > 0 ) {
	       $need_to_infer_exons = 1;
	   }
       }

       $need_to_infer_exons = 0 if $noinfer; #NML

       if ($need_to_infer_exons) {
	   # remove exons and introns from group -
	   # we will infer exons later, and we
	   # can always infer introns from exons
	   foreach my $group (@groups) {
	       @$group = 
		 grep {
		     my $type = $_->primary_tag();
		     if ($type eq 'exon') {
			 # keep track of all removed exons,
			 # so we can do a sanity check later
			 push(@removed_exons, $_);
		     }
		     $type ne 'exon' && $type ne 'intron'
		 } @$group;
	   }
	   # get rid of any groups that have zero members
	   @groups = grep {scalar(@$_)} @groups;
       }
   }
   # --- END OF MAGIC ---
   
   # LOGICAL ASSERTION
   if (grep {!scalar(@$_)} @groups) {
       $self->throw("ASSERTION ERROR: empty group");
   }

   # LOGGING
   if ($self->verbose > 0) {
       printf STDERR "GROUPS:\n";
       foreach my $group (@groups) {
	   $self->_write_group($group, $group_tag);
       }
   }
   # -

   # --------- FINISHED GROUPING -------------


   # TYPE CONTAINMENT HIERARCHY (aka partonomy)
   # set the containment hierarchy if desired
   # see docs for structure_type() method
   if ($structure_type) {
       if ($structure_type == 1) {
	   $self->partonomy(
                            {CDS => 'gene',
                             exon => 'CDS',
                             intron => 'CDS',
                            }
                           );
       }
       else {
	   $self->throw("structure_type $structure_type is currently unknown");
       }
   }

   # see if we have an obvious resolver_tag
   if ($use_magic) {
       foreach my $sf (@all_seq_features) {
	   if ($sf->has_tag('derived_from')) {
	       $resolver_tag = 'derived_from';
	   }
       }
   }

   if ($use_magic) {
       # point all feature types without a container type to the root type.
       #
       # for example, if we have an unanticipated feature_type, say
       # 'aberration', this should by default point to the parent 'gene'
       foreach my $group (@groups) {
	   my @sfs = @$group;
	   if (@sfs > 1) {
	       foreach my $sf (@sfs) {
		   my $type = $sf->primary_tag;
		   next if $type eq 'gene';
		   my $container_type = $self->get_container_type($type);
		   if (!$container_type) {
		       $self->partonomy->{$type} = 'gene';
		   }
	       }
	   }
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
                                         -resolver_method=>$resolver_method,
					 -resolver_tag=>$resolver_tag);
   
   # restore settings
   $self->partonomy($old_partonomy);

   # restore settings
   $self->structure_type($old_structure_type);

   # modify the original Seq object - the top seqfeatures are now
   # the top features from each group
   $seq->remove_SeqFeatures;
   $seq->add_SeqFeature($_) foreach @top_sfs;

   # --------- FINISHED UNFLATTENING -------------

   # lets see if there are any post-unflattening tasks we need to do

   

   # INFERRING mRNAs
   if ($need_to_infer_mRNAs) {
       if ($self->verbose > 0) {
	   printf STDERR "** INFERRING mRNA from CDS\n";
       }
       $self->infer_mRNA_from_CDS(-seq=>$seq, -noinfer=>$noinfer);
   }

   # INFERRING exons
   if ($need_to_infer_exons) {

       # infer exons, one group/gene at a time
       foreach my $sf (@top_sfs) {
	   my @sub_sfs = ($sf, $sf->get_all_SeqFeatures);
	   $self->feature_from_splitloc(-features=>\@sub_sfs);
       }

       # some exons are stated explicitly; ie there is an "exon" feature
       # most exons are inferred; ie there is a "mRNA" feature with
       # split locations
       #
       # if there were exons explicitly stated in the entry, we need to
       # do two things:
       #
       # make sure these exons are consistent with the inferred exons
       #  (you never know)
       #
       # transfer annotation (tag-vals) from the explicit exon to the
       # new inferred exon
       if (@removed_exons) {
	   my @allfeats = $seq->get_all_SeqFeatures;

	   # find all the inferred exons that are children of mRNA
	   my @mrnas =  grep {$_->primary_tag eq 'mRNA'} @allfeats;
	   my @exons =  
	     grep {$_->primary_tag eq 'exon'}
	       map {$_->get_SeqFeatures} @mrnas;

	   my %exon_h = (); 	   # index of exons by location;

	   # there CAN be >1 exon at a location; we can represent these redundantly
	   # (ie as a tree, not a graph)
	   push(@{$exon_h{$self->_locstr($_)}}, $_) foreach @exons;
	   my @problems = ();      # list of problems;
	                           # each problem is a 
	                           # [$severity, $description] pair
	   my $problem = '';
	   my ($n_exons, $n_removed_exons) =
	     (scalar(keys %exon_h), scalar(@removed_exons));
	   foreach my $removed_exon (@removed_exons) {
	       my $locstr = $self->_locstr($removed_exon);
	       my $inferred_exons = $exon_h{$locstr};
	       delete $exon_h{$locstr};
	       if ($inferred_exons) {
		   my %exons_done = ();
		   foreach my $exon (@$inferred_exons) {

		       # make sure we don't move stuff twice
		       next if $exons_done{$exon};
		       $exons_done{$exon} = 1;

		       # we need to tranfer any tag-values from the explicit
		       # exon to the implicit exon
		       foreach my $tag ($removed_exon->get_all_tags) {
			   my @vals = $removed_exon->get_tag_values($tag);
			   if (!$exon->can("add_tag_value")) {
			       # I'm puzzled as to what should be done here;
			       # SeqFeatureIs are not necessarily mutable,
			       # but we know that in practice the implementing
			       # class is mutable
			       $self->throw("The SeqFeature object does not ".
					    "implement add_tag_value()");
			   }
			   $exon->add_tag_value($tag, @vals);
		       }
		   }
	       } 
               else {
                   # no exons inferred at $locstr
		   push(@problems,
			[1, 
			 "there is a conflict with exons; there was an explicitly ".
			 "stated exon with location $locstr, yet I cannot generate ".
			 "this exon from the supplied mRNA locations\n"]);
	       }
	   }
	   # do we have any inferred exons left over, that were not
	   # covered in the explicit exons?
	   if (keys %exon_h) {
	       # TODO - we ignore this problem for now
	       push(@problems,
		    [1,
		     sprintf("There are some inferred exons that are not in the ".
			     "explicit exon list; they are the exons at locations:\n".
			     join("\n", keys %exon_h)."\n")]);
	   }

	   # report any problems
	   if (@problems) {
	       my $thresh = $self->error_threshold;
	       my @bad_problems = grep {$_->[0] > $thresh} @problems;
	       if (@bad_problems) {
		   printf STDERR "PROBLEM:\n";
		   $self->_write_hier(\@top_sfs);
		   # TODO - allow more fine grained control over this
		   $self->{_problems_reported} = 1;
		   $self->throw(join("\n",
				     map {"@$_"} @bad_problems));
	       }
	       $self->problem(@$_) foreach @problems;
	   }
       }
   }    
   # --- end of inferring exons --

   # return new top level features; this can also 
   # be retrieved via
   #   $seq->get_SeqFeatures();
#   return @top_sfs;
   return $seq->get_SeqFeatures;
}

# _split_group_if_disconnected([@sfs])
#
# as well as having the same group_tag, a group should be spatially
# connected. if not, then the group should be split into subgroups.
# this turns out to be necessary in the case of multicopy genes.
# the standard way to represent these is as spatially disconnected
# gene models (usually a 'gene' feature and some kind of RNA feature)
# with the same group tag; the code below will split these into 
# seperate groups, one per copy.

sub _split_group_if_disconnected {
    my $self = shift;
    my $group = shift;
    my @sfs = @$group;
    my @ranges =
      Bio::Range->disconnected_ranges(@sfs);
    my @groups;
    if (@ranges == 0) {
	$self->throw("ASSERTION ERROR");
    }
    elsif (@ranges == 1) {
	# no need to split the group
	@groups = ($group);
    }
    else {
	# @ranges > 1
	# split the group into disconnected ranges
	if ($self->verbose > 0) {
	    printf STDERR "GROUP PRE-SPLIT:\n";
	    $self->_write_group($group, $self->group_tag);
	}
	@groups =
	  map {
	      my $range = $_;
	      [grep {
		  $_->intersection($range);
	      } @sfs]
	  } @ranges;
	if ($self->verbose > 0) {
	    printf STDERR "SPLIT GROUPS:\n";
	    $self->_write_group($_, $self->group_tag) foreach @groups;	    
	}
    }
    return @groups;
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

	if ($self->verbose > 0) {
	    printf STDERR "REMOVING DUPLICATES:\n";
	}

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
 Function: iterates over groups, calling unflatten_group() [see below]
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
grained control over how the unflattening process.

=cut

sub unflatten_groups{
   my ($self,@args) = @_;
   my($groups, $resolver_method, $resolver_tag) =
     $self->_rearrange([qw(GROUPS
                           RESOLVER_METHOD
			   RESOLVER_TAG
                          )],
                          @args);

   # this is just a simple wrapper for unflatten_group()
   return 
     map {
         $self->unflatten_group(-group=>$_,
                                -resolver_method=>$resolver_method,
				-resolver_tag=>$resolver_tag)
     } @$groups;
}

=head2 unflatten_group

 Title   : unflatten_group
 Usage   :
 Function: nests a group of features into a feature containment hierarchy
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


NOTE: You should not need to call this method, unless you want fine
grained control over how the unflattening process.

=cut

sub unflatten_group{
   my ($self,@args) = @_;

   my($group, $resolver_method, $resolver_tag) =
     $self->_rearrange([qw(GROUP
                           RESOLVER_METHOD
			   RESOLVER_TAG
                          )],
                          @args);

   if ($self->verbose > 0) {
       printf STDERR "UNFLATTENING GROUP:\n";
       $self->_write_group($group, $self->group_tag);
   }

   my @sfs = @$group;

   # we can safely ignore singletons (e.g. [source])
   return $sfs[0] if @sfs == 1;

   my $partonomy = $self->partonomy;

   # $resolver_method is a reference to a SUB that will resolve
   # ambiguous parent/child containment; for example, determining
   # which mRNAs go with which CDSs
   $resolver_method = $resolver_method || \&_resolve_container_for_sf;

   # TAG BASED RESOLVING OF HIERARCHIES
   #
   # if the user specifies $resolver_tag, then we use this tag
   # to pair up ambiguous parents and children;
   #
   # for example, the CDS feature may have a resolver tag of /derives_from
   # which is a 'foreign key' into the /label tag of the mRNA feature
   #
   # this kind of tag-based resolution is possible for a certain subset
   # of genbank records
   #
   # if no resolver tag is specified, we revert to the normal
   # resolver_method
   if ($resolver_tag) {
       my $backup_resolver_method = $resolver_method;
       # closure: $resolver_tag is remembered by this sub
       my $sub = 
	 sub {
	     my ($self, $sf, @possible_container_sfs) = @_;
	     my @container_sfs = ();
	     if ($sf->has_tag($resolver_tag)) {
		 my ($resolver_tagval) = $sf->get_tag_values($resolver_tag);
		 # if a feature has a resolver_tag (e.g. /derives_from)
		 # this specifies the /product, /symbol or /label for the
		 # parent feature
		 @container_sfs = 
		   grep {
		       my $match = 0;
		       $self->_write_sf($_) if $self->verbose > 0;
		       foreach my $tag (qw(product symbol label)) {
			   if ($_->has_tag($tag)) {
			       my @vals =
				 $_->get_tag_values($tag);
			       if (grep {$_ eq $resolver_tagval} @vals) {
				   $match = 1;
				   last;
			       }
			   }   
		       }
		       $match;
		   } @possible_container_sfs;
	     } 
	     else {
		 return $backup_resolver_method->($sf, @possible_container_sfs);
	     }
	     return map {$_=>0} @container_sfs;
	 };
       $resolver_method = $sub;
   }
   else {
       # CONDITION: $resolver_tag is NOT set
       $self->throw("assertion error") if $resolver_tag;
   }
   # we have now set $resolver_method to a subroutine for
   # disambiguatimng parent/child relationships. we will
   # now build the whole containment hierarchy for this group


   # FIND TOP/ROOT SEQFEATURES
   #
   # find all the features for which there is no
   # containing feature type (eg genes)
   my @top_sfs =
     grep { 
         !$self->get_container_type($_->primary_tag);
     } @sfs;

   # CONDITION: there must be at most one root
   if (@top_sfs > 1) {
       $self->_write_group($group, $self->group_tag);
       printf STDERR "TOP SFS:\n";
       $self->_write_sf($_) foreach @top_sfs;
       $self->throw("multiple top-sfs in group");
   }
   my $top_sf = $top_sfs[0];

   # CREATE INDEX OF SEQFEATURES BY TYPE
   my %sfs_by_type = ();
   foreach my $sf (@sfs) {
       push(@{$sfs_by_type{$sf->primary_tag}}, $sf);
   }

   # containment index; keyed by child; lookup parent
   # note: this index uses the stringified object reference of
   # the object as a surrogate lookup key

   my %container = ();   # child -> parent

   # ALGORITHM: build containment graph
   #
   # find all possible containers for each SF;
   # for instance, for a CDS, the possible containers are all
   # the mRNAs in the same group. For a mRNA, the possible
   # containers are any SFs of type 'gene' (should only be 1).
   # (these container-type mappings can be overridden)
   #
   # contention is resolved by checking coordinates of splice sites
   # (this is the default, but can be overridden)
   #
   # most of the time, there is no problem identifying a unique
   # parent for every child; this can be ambiguous when constructing
   # CDS to mRNA relationships with lots of alternate splicing
   #
   # a hash of child->parent relationships is constructed (%container)
   # any mappings that need further resolution (eg CDS to mRNA) are
   # placed in %unresolved

   # %unresolved index
   # (keyed by stringified object reference of child seqfeature)
   my %unresolved = ();    # child -> [parent,score] to be resolved
                           
   # index of seqfeatures by their stringified object reference;
   # this is essentially a way of 'reviving' an object from its stringified
   # reference
   # (see NOTE ON USING OBJECTS AS KEYS IN HASHES, below)
   my %idxsf = map {$_=>$_} @sfs;

   foreach my $sf (@sfs) {
       my $type = $sf->primary_tag;

       # container type (e.g. the container type for CDS is usually mRNA)
       my $container_type = 
         $self->get_container_type($type);
       if ($container_type) {

           my @possible_container_sfs =
             @{$sfs_by_type{$container_type} || []};
           # we now have a list of possible containers
           # (eg for a CDS in an alternately spliced gene, this
           #  would be a list of all the mRNAs for this gene)

	   if (!@possible_container_sfs) {
	       # root of hierarchy
	   }
	   else {
	       if (@possible_container_sfs == 1) {
                   # this is the easy situation, whereby the containment
                   # hierarchy is unambiguous. this will probably be the
                   # case if the genbank record has no alternate splicing
                   # within it

		   # ONE OPTION ONLY - resolved!
		   $container{$sf} = $possible_container_sfs[0];

	       }
	       else {
		   # MULTIPLE CONTAINER CHOICES
		   $self->throw("ASSERTION ERROR") unless @possible_container_sfs > 1;

                   # push this onto the %unresolved graph, and deal with it
                   # later

                   # for now we hardcode things such that the only type 
                   # with ambiguous parents is a CDS; if this is violated,
                   # it has a weak problem class of '1' so the API user
                   # can easily set things to ignore these
		   if ($sf->primary_tag ne 'CDS') {
		       $self->problem(1,
				      "multiple container choice for non-CDS; ".
				      "CDS to mRNA should be the only ".
				      "relationships requiring resolving",
				      $sf);
		   }

                   # previously we set the SUB $resolver_method
                   $self->throw("ASSERTION ERROR")
                     unless $resolver_method;

                   # $resolver_method will assign scores to
                   # parent/child combinations; later on we
                   # will use these scores to find the optimal
                   # parent/child pairings

                   # the default $resolver_method uses splice sites to
                   # score possible parent/child matches

		   my %container_sfh =
		     $resolver_method->($self, $sf, @possible_container_sfs);
                   if (!%container_sfh) {
                       $self->problem(2,
                                      "no containers possible for SeqFeature of ".
                                      "type: $type; this SF is being placed at ".
                                      "root level",
                                      $sf);
                       # RESOLVED! (sort of - placed at root/gene level)
                       $container{$sf} = $top_sf;

                       # this sort of thing happens if the record is
                       # badly messed up and there is absolutely no indication
                       # of where to put the CDS. Perhaps we should just
                       # place it with a random mRNA?
                   }
		   foreach my $jsf (keys %container_sfh) {

                       # add [score, parent] pairs to the %unresolved
                       # lookup table/graph
		       push(@{$unresolved{$sf}}, 
			    [$idxsf{$jsf}, $container_sfh{$jsf} || 0]);
		   }
	       }
	   }
       }
       else {
           # CONDITION:
           # not container type for $sf->primary_tag
           
           # CONDITION:
	   # $sf must be a root/top node (eg gene)
       }
   }

   if (0) {

       # CODE CURRENTLY DISABLED

       # we require a 1:1 mapping between mRNAs and CDSs;
       # create artificial duplicates if we can't do this...
       if (%unresolved) {
           my %childh = map {$_=>1} keys %unresolved;
           my %parenth = map {$_->[0]=>1} map {@$_} values %unresolved;
           if ($self->verbose > 0) {
               printf STDERR "MATCHING %d CHILDREN TO %d PARENTS\n",
                 scalar(keys %childh), scalar(keys %parenth);
           }
           # 99.99% of the time in genbank genomic record of structure type 0, we
           # see one CDS for every mRNA; one exception is the S Pombe
           # genome, which is all CDS, bar a few spurious mRNAs; we have to
           # filter out the spurious mRNAs in this case
           #
           # another strange case is in the mouse genome, NT_078847.1
           # for Pcdh13 you will notice there is 4 mRNAs and 5 CDSs.
           # most unusual! 
           # I'm at a loss for a really clever thing to do here. I think the
           # best thing is to create duplicate features to preserve the 1:1 mapping
           #       my $suffix_id = 1;
           #       while (keys %childh > keys %parenth) {
           #           
           #       }
       }
   }

   # DEBUGGING CODE
   if ($self->verbose > 0 && scalar(keys %unresolved)) {
       printf STDERR "UNRESOLVED PAIRS:\n";
       foreach my $childsf (keys %unresolved) {
	   my @poss = @{$unresolved{$childsf}};
	   foreach my $p (@poss) {
	       my $parentsf = $p->[0];
	       $childsf = $idxsf{$childsf};
               my @clabels = ($childsf->get_tagset_values(qw(protein_id label product)), "?");
               my @plabels = ($parentsf->get_tagset_values(qw(transcript_id label product)), "?");
	       printf STDERR
                      ("  PAIR: $clabels[0] => $plabels[0]  (of %d)\n", 
                       scalar(@poss));
	   }
       }
   } # -- end of verbose

   # Now we have to fully resolve the containment hierarchy; remember,
   # the graph %container has the fully resolved child->parent links;
   #
   # the graph %unresolved is keyed by children missing parents; we
   # need to put all these orphans in the %container graph
   #
   # we do this using the scores in %unresolved, with the
   # find_best_matches() algorithm
   my $unresolved_problem_reported = 0;
   if (%unresolved) {
       my $new_pairs =
	 $self->find_best_matches(\%unresolved, []);
       if (!$new_pairs) {
           my ($g) = $sfs[0]->get_tagset_values($self->group_tag || 'gene');
	   $self->problem(2,
			  "Could not resolve hierarchy for $g");
           $new_pairs = [];
           $unresolved_problem_reported = 1;
       }
       foreach my $pair (@$new_pairs) {
	   if ($self->verbose > 0) {
	       printf STDERR "  resolved pair @$pair\n";
	   }
	   $container{$pair->[0]} = $pair->[1];
           delete $unresolved{$pair->[0]};
       }
   }

   # CONDITION: containment hierarchy resolved
   if (%unresolved) {
       $self->throw("UNRESOLVED: %unresolved")
         unless $unresolved_problem_reported;
   }

   # make nested SeqFeature hierarchy from @containment_pairs
   # ie put child SeqFeatures into parent SeqFeatures
   my @top = ();
   foreach my $sf (@sfs) {
       my $container_sf = $container{$sf};
       if ($container_sf) {
           # make $sf nested inside $container_sf

           # first check if the container spatially contains the containee
           if ($container_sf->contains($sf)) {
               # add containee
	       $container_sf->add_SeqFeature($sf);
           }
           else {
               # weird case - the container does NOT spatially
               # contain the containee;
               # we expand and throw a warning
               #
               # for an example of this see ZFP91-CNTF dicistronic gene
               # in NCBI chrom 11 build 34.3
	       $self->problem(1,
			      "Container feature does not spatially contain ".
                              "subfeature. Perhaps this is a dicistronic gene? ".
                              "I am expanding the parent feature",
			      $container_sf,
			      $sf);
	       $container_sf->add_SeqFeature($sf, 'EXPAND');
           }
       }
       else {
           push(@top, $sf);
       }
   }
   return @top;
} # -- end of unflatten_group

# -------
# A NOTE ON USING OBJECTS AS KEYS IN HASHES (stringified objects)
#
# Often we with to use seqfeatures as keys in a hashtable; because seqfeatures
# in bioperl have no unique ID, we use a surrogate ID in the form of the
# stringified object references - this is just what you get if you say
#
#  print "$sf\n";
#
# this is guaranteed to be unique (within a particular perl execution)
#
# often we want to 'revive' the objects used as keys in a hash - once the
# objects are used as keys, remember it is the *strings* used as keys and
# not the object itself, so the object needs to be revived using another
# hashtable that looks like this
#
#    %sfidx = map { $_ => $_ } @sfs
#
# -------


# recursively finds the best set of pairings from a matrix of possible pairings
#
# tries to make sure nothing is unpaired
#
# given a matrix of POSSIBLE matches
#  (matrix expressed as hash/lookup; keyed by child object; val = [parent, score]
#
# 
sub find_best_matches {
    my $self = shift;
    my $matrix = shift;
    my $pairs = shift;        # [child,parent] pairs already selected

    my $verbose = $self->verbose;
    #################################print "I";
    if ($verbose > 0) {
	printf STDERR "find_best_matches: (/%d)\n", scalar(@$pairs);
    }

    my %selected_children = map {($_->[0]=>1)} @$pairs;
    my %selected_parents = map {($_->[1]=>1)} @$pairs;
    
    # make a copy of the matrix with the portions still to be
    # resolved
    my %unresolved_parents = ();
    my %unresolved =
      map {
          if ($verbose > 0) {
              printf STDERR "  $_ : %s\n", join("; ", map {"[@$_]"} @{$matrix->{$_}});
          }
	  if ($selected_children{$_}) {
	      ();
	  }
	  else {
	      my @parents =
		grep {
		    !$selected_parents{$_->[0]}
		} @{$matrix->{$_}};
              $unresolved_parents{$_} = 1 foreach @parents;
              # new parents
	      ($_ => [@parents]);
	  }
      } keys %$matrix;
    
    my @I = keys %unresolved;

    return $pairs if !scalar(keys %unresolved_parents);
    # NECESSARY CONDITION:
    # all possible parents have a child match

    return $pairs if !scalar(@I);
    # NECESSARY CONDITION:
    # all possible children have a parent match

    # give those with fewest choices highest priority
    @I = sort {
	# n possible parents
	scalar(@{$unresolved{$a}}) 
	  <=>
	    scalar(@{$unresolved{$b}}) ;
    } @I;
    
    my $csf = shift @I;

    my @J = @{$unresolved{$csf}};  # array of [parent, score]

    # sort by score, highest first
    @J =
      sort {
	  $b->[1] <=> $a->[1]
      } @J;

    # select pair(s) from remaining matrix of possible pairs
    # by iterating through possible parents

    my $successful_pairs;
    foreach my $j (@J) {
	my ($psf, $score) = @$j;
	# would selecting $csf, $psf as a pair
	# remove all choices from another?
	my $bad = 0;
	foreach my $sf (@I) {
	    if (!grep {$_->[0] ne $psf} @{$unresolved{$sf}}) {
		# $psf was the only parent choice for $sf
		$bad = 1;
		last;
	    }
	}
	if (!$bad) {
	    my $pair = [$csf, $psf];
	    my $new_pairs = [@$pairs, $pair];
	    my $set = $self->find_best_matches($matrix, $new_pairs);
	    if ($set) {
		$successful_pairs = $set;
		last;
	    }
	}
    }
    # success
    return $successful_pairs if $successful_pairs;
    # fail
    return 0;
}

# ----------------------------------------------
# writes a group to stdout
#
# mostly for logging/debugging
# ----------------------------------------------
sub _write_group {
    my $self = shift;
    my $group = shift;
    my $group_tag = shift || 'gene';

    my $f = $group->[0];
    my $label = '?';
    if ($f->has_tag($group_tag)) {
	($label) = $f->get_tag_values($group_tag);
    }
    if( $self->verbose > 0 ) { 
	printf STDERR ("  GROUP [%s]:%s\n",
	       $label,
	       join(' ',
		    map { $_->primary_tag } @$group));
    }

}

sub _write_sf {
    my $self = shift;
    my $sf = shift;
    printf STDERR "TYPE:%s\n", $sf->primary_tag;
    return;
}

sub _write_sf_detail {
    my $self = shift;
    my $sf = shift;
    printf STDERR "TYPE:%s\n", $sf->primary_tag;
    my @locs = $sf->location->each_Location;
    printf STDERR "  %s,%s [%s]\n", $_->start, $_->end, $_->strand foreach @locs;
    return;
}

sub _write_hier {
    my $self = shift;
    my @sfs = @{shift || []};
    my $indent = shift || 0;
    if( $self->verbose > 0 ) {
	foreach my $sf (@sfs) {
	    my $label = '?';
	    if ($sf->has_tag('product')) {
		($label) = $sf->get_tag_values('product');
	    }
	    printf STDERR "%s%s $label\n", '  ' x $indent, $sf->primary_tag;
	    my @sub_sfs = $sf->sub_SeqFeature;
	    $self->_write_hier(\@sub_sfs, $indent+1);
	}
    }
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
   my $start = $sf->start;
   my $end = $sf->end;
   my $splice_uniq_str = "@coords";
   
   my @sf_score_pairs = ();
   # a CDS is contained by a mRNA if the locations of the splice
   # coordinates are identical
   foreach (@possible_container_sfs) {
       my @container_coords = $self->_get_splice_coords_for_sf($_);
       my $inside = 
	 !$splice_uniq_str || 
	   index("@container_coords", $splice_uniq_str) > -1;
       if ($inside) {
           # the container cannot be smaller than the thing contained
           if ($_->start > $start || $_->end < $end) {
               $inside = 0;
           }
       }


       # SPECIAL CASE FOR /ribosomal_slippage
       # See: http://www.ncbi.nlm.nih.gov/collab/FT/
       if (!$inside && $sf->has_tag('ribosomal_slippage')) {
	   if ($self->verbose > 0) {
	       printf STDERR "    Checking for ribosomal_slippage\n";
	   }

           # TODO: rewrite this to match introns;
           #  each slippage will be a "fake" small CDS exon
	   my @transcript_splice_sites = @container_coords;
	   my @cds_splice_sites = @coords;
           ##printf STDERR "xxTR SSs: @transcript_splice_sites :: %s\n", $_->get_tag_values('product');
           ##printf STDERR "xxCD SSs: @cds_splice_sites :: %s\n\n", $sf->get_tag_values('product');

	   # find the the first splice site within the CDS
	   while (scalar(@transcript_splice_sites) &&
		  $transcript_splice_sites[0] < $cds_splice_sites[0]) {
	       shift @transcript_splice_sites;
	   }

           ##print STDERR "TR SSs: @transcript_splice_sites\n";
           ##print STDERR "CD SSs: @cds_splice_sites\n\n";

	   if (!(scalar(@transcript_splice_sites)) ||
                 $transcript_splice_sites[0] == $cds_splice_sites[0]) {

               # we will now try and align all splice remaining sites in the transcript and CDS;
               # any splice site that can't be aligned is assumed to be a ribosomal slippage

	       my @slips = ();
	       my $in_exon = 1;
	       $inside = 1;   # innocent until proven guilty..
	       while (@cds_splice_sites) {
		   if (!@transcript_splice_sites) {

                       # ribosomal slippage is after the last transcript splice site
                       # Example: (NC_00007, isoform 3 of PEG10)
                       #     mRNA            join(85682..85903,92646..99007)
                       #     mRNA            join(85682..85903,92646..99007)
                       #     CDS             join(85899..85903,92646..93825,93825..94994)

                       # OR: None of the splice sites align;
                       #  may be a single CDS exon with one slippage inside it.
                       # Example: (NC_00007, isoform 4 of PEG10)
                       #     mRNA            join(85637..85892,92646..99007)
                       #     CDS             join(92767..93825,93825..94994)
                       
                       # Yes, this code is repeated below...
                       my $p1 = shift @cds_splice_sites;
                       my $p2 = shift @cds_splice_sites;
                       if ($self->verbose > 0) {
                           printf STDERR "    Found the ribosomal_slippage: $p1..$p2\n";
                       }
                       push(@slips, ($p2-$p1)-1);
		   }
		   elsif ($cds_splice_sites[0] == $transcript_splice_sites[0]) {
                       # splice sites align: this is not the slippage
		       shift @cds_splice_sites;
		       shift @transcript_splice_sites;
                       ##print STDERR "MATCH\n";
		   }
		   else {
		       # mismatch
		       if ($cds_splice_sites[0] < $transcript_splice_sites[0]) {
			   # potential slippage
			   #             v
			   # ---TTTTTTTTTT----
			   # ---CCCC--CCCC----
			   #       ^

			   my $p1 = shift @cds_splice_sites;
			   my $p2 = shift @cds_splice_sites;
			   if ($self->verbose > 0) {
			       printf STDERR "    Found the ribosomal_slippage: $p1..$p2\n";
			   }
			   push(@slips, ($p2-$p1)-1);
		       }
		       else {
			   # not a potential ribosomal slippage
			   $inside = 0; # guilty!
                           ##print STDERR "FAIL\n";
			   last;
		       }
		   }
	       }
	       if ($inside) {
		   # TODO: this is currently completely arbitrary. How many ribosomal slippages do we allow?
		   # perhaps we need some mini-statistical model here....?
		   if (@slips > 1) {
		       $inside = 0;
		   }
		   # TODO: this is currently completely arbitrary. What is the maximum size of a ribosomal slippage?
		   # perhaps we need some mini-statistical model here....?
		   if (grep {$_ > 2} @slips) {
		       $inside = 0;
		   }
	       }
	   }
	   else {
	       # not a ribosomal_slippage, sorry
	   }
       }
       if ($self->verbose > 0) {
	   printf STDERR "    Checking containment:[$inside] (@container_coords) IN ($splice_uniq_str)\n";
       }
       if ($inside) {
	   # SCORE: matching (ss-scoords+2)/(n-container-ss-coords+2)
	   my $score =
	     (scalar(@coords)+2)/(scalar(@container_coords)+2);
	   push(@sf_score_pairs,
		$_=>$score);
       }
   }
   return @sf_score_pairs;
}

sub _get_splice_coords_for_sf {
    my $self = shift;
    my $sf = shift;

   my @locs = $sf->location;
   if ($sf->location->isa("Bio::Location::SplitLocationI")) {
       @locs = $sf->location->each_Location;
   }

   # get an ordered list of (start, end) positions

#   my @coords =
#     map {
#         $_->strand > 0 ? ($_->start, $_->end) : ($_->end, $_->start)
#     } @locs;

    my @coords = map {($_->start, $_->end)} @locs;

   # remove first and last leaving only splice sites
   pop @coords;
   shift @coords;
    return @coords;
}

=head2 feature_from_splitloc

 Title   : feature_from_splitloc
 Usage   : $unflattener->feature_from_splitloc(-features=>$sfs);
 Function:
 Example :
 Returns : 
 Args    : see below

At this time all this method does is generate exons for mRNA or other RNA features

Arguments:

  -feature:    a Bio::SeqFeatureI object (that conforms to Bio::FeatureHolderI)
  -seq:        a Bio::SeqI object that contains Bio::SeqFeatureI objects
  -features:   an arrayref of Bio::SeqFeatureI object


=cut

sub feature_from_splitloc{
   my ($self,@args) = @_;

   my($sf, $seq, $sfs) =
     $self->_rearrange([qw(FEATURE
                           SEQ
			   FEATURES
                          )],
                          @args);
   my @sfs = (@{$sfs || []});
   push(@sfs, $sf) if $sf;
   if ($seq) {
       $seq->isa("Bio::SeqI") || $self->throw("$seq NOT A SeqI");
       @sfs = $seq->get_all_SeqFeatures;
   }
   my @exons = grep {$_->primary_tag eq 'exon'} @sfs;
   if (@exons) {
       $self->problem(2,
		      "There are already exons, so I will not infer exons");
   }

   # index of features by type+location
   my %loc_h = ();

   # infer for every feature
   foreach my $sf (@sfs) {

       $sf->isa("Bio::SeqFeatureI") || $self->throw("$sf NOT A SeqFeatureI");
       $sf->isa("Bio::FeatureHolderI") || $self->throw("$sf NOT A FeatureHolderI");

       my $type = $sf->primary_tag;
       next unless $type eq 'mRNA' or $type =~ /RNA/;

       # an mRNA from genbank will have a discontinuous location,
       # with each sub-location being equivalent to an exon
       my @locs = $sf->location;

       if ($sf->location->isa("Bio::Location::SplitLocationI")) {
           @locs = $sf->location->each_Location;
       }

       if (!@locs) {
           use Data::Dumper;
           print Dumper $sf;
	   $self->throw("ASSERTION ERROR: sf has no location objects");
       }

       # make exons from locations
       my @subsfs =
         map {
             my $subsf = Bio::SeqFeature::Generic->new(-location=>$_,
                                                       -primary_tag=>'exon');
             ## Provide seq_id to new feature:
             $subsf->seq_id($sf->seq_id) if $sf->seq_id;
             $subsf->source_tag($sf->source_tag) if $sf->source_tag;
             ## Transfer /locus_tag and /gene tag values to inferred
             ## features.  TODO: Perhaps? this should not be done
             ## indiscriminantly but rather by virtue of the setting
             ## of group_tag.
             foreach my $tag (grep /gene|locus_tag/, $sf->get_all_tags) {
                 my @vals = $sf->get_tag_values($tag);
                 $subsf->add_tag_value($tag, @vals);
             }

	     my $locstr = 'exon::'.$self->_locstr($subsf);

	     # re-use feature if type and location the same
	     if ($loc_h{$locstr}) {
		 $subsf = $loc_h{$locstr};
	     }
	     else {
		 $loc_h{$locstr} = $subsf;
	     }
             $subsf;
         } @locs;
       
       # PARANOID CHECK
       $self->_check_order_is_consistent($sf->location->strand,@subsfs);
       #----

       $sf->location(Bio::Location::Simple->new());

       # we allow the exons to define the boundaries of the transcript
       $sf->add_SeqFeature($_, 'EXPAND') foreach @subsfs;


       if (!$sf->location->strand) {
	   # correct weird bioperl bug in previous versions;
	   # strand was not being set correctly
	   $sf->location->strand($subsfs[0]->location->strand);
       }

       
   }
   return;
}

#sub merge_features_with_same_loc {
#   my ($self,@args) = @_;

#   my($sfs, $seq) =
#     $self->_rearrange([qw(FEATURES
#                           SEQ
#                          )],
#                          @args);
#   my @sfs = (@$sfs);
#   if ($seq) {
#       $seq->isa("Bio::SeqI") || $self->throw("$seq NOT A SeqI");
#       @sfs = $seq->get_all_SeqFeatures;
#   }

   
#   my %loc_h = ();
#   foreach my $sf (@sfs) {
#       my $type = $sf->primary_tag;
#       my $locstr = $self->_locstr($sf);
##       $loc_h{$type.$locstr}
#       push(@{$exon_h{$self->_locstr($_)}}, $_) foreach @exons;
#   }
#}

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

a "type 2" containment hierarchy is mixed type "0" and "1" (for
example, see ftp.ncbi.nih.gov/genomes/Schizosaccharomyces_pombe/)

=cut

sub infer_mRNA_from_CDS{
   my ($self,@args) = @_;

   my($sf, $seq, $noinfer) =
     $self->_rearrange([qw(FEATURE
                           SEQ
			   NOINFER
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
       if ($self->verbose > 0) {
           printf STDERR "    Checking $sf %s\n", $sf->primary_tag;
       }
       
       if ($sf->primary_tag eq 'mRNA') {
	   $self->problem(2,
			  "Inferring mRNAs when there are already mRNAs present");
       }

       my @cdsl = grep {$_->primary_tag eq 'CDS' } $sf->get_SeqFeatures;
       if (@cdsl) {
	   my @children = grep {$_->primary_tag ne 'CDS'} $sf->get_SeqFeatures;
	   my @mrnas = ();


	   foreach my $cds (@cdsl) {
	       
               if ($self->verbose > 0) {
                   print "    Inferring mRNA from CDS $cds\n";
               }
               $self->_check_order_is_consistent($cds->location->strand,$cds->location->each_Location);
               
	       my $loc = Bio::Location::Split->new;
	       foreach my $cdsexonloc ($cds->location->each_Location) {
		   my $subloc =
		     Bio::Location::Simple->new(-start=>$cdsexonloc->start,
						-end=>$cdsexonloc->end,
						-strand=>$cdsexonloc->strand);
		   $loc->add_sub_Location($subloc);
	       }
		if ($noinfer) {
		    push(@mrnas, $cds);
		}
		else {
#		    share the same location
		    my $mrna =
			Bio::SeqFeature::Generic->new(-location=>$loc,
				-primary_tag=>'mRNA');

##		    Provide seq_id to new feature:
		    $mrna->seq_id($cds->seq_id) if $cds->seq_id;
		    $mrna->source_tag($cds->source_tag) if $cds->source_tag;

		    $self->_check_order_is_consistent($mrna->location->strand,$mrna->location->each_Location);

#		    make the mRNA hold the CDS; no EXPAND option,
#		    the CDS cannot be wider than the mRNA
		    $mrna->add_SeqFeature($cds);

#		    mRNA steals children of CDS
		    foreach my $subsf ($cds->get_SeqFeatures) {
			$mrna->add_SeqFeature($subsf);
		    }
		    $cds->remove_SeqFeatures;
		    push(@mrnas, $mrna);
		}
	   }
#	   change gene/CDS to gene/mRNA
	   $sf->remove_SeqFeatures;
	   $sf->add_SeqFeature($_) foreach (@mrnas, @children);
       }
   }
   return;
   

}

=head2 remove_types

 Title   : remove_types
 Usage   : $unf->remove_types(-seq=>$seq, -types=>["mRNA"]);
 Function:
 Example :
 Returns : 
 Args    :

removes features of a set type

useful for pre-filtering a genbank record; eg to get rid of STSs

also, there is no way to unflatten
ftp.ncbi.nih.gov/genomes/Schizosaccharomyces_pombe/ UNLESS the bogus
mRNAs in these records are removed (or changed to a different type) -
they just confuse things too much

=cut

sub remove_types{
   my ($self,@args) = @_;

   my($seq, $types) =
     $self->_rearrange([qw(
                           SEQ
			   TYPES
                          )],
                          @args);
   $seq->isa("Bio::SeqI") || $self->throw("$seq NOT A SeqI");
   my @sfs = $seq->get_all_SeqFeatures;
   my %rh = map {$_=>1} @$types;
   @sfs = grep {!$rh{$_->primary_tag}} @sfs;
   $seq->remove_SeqFeatures;
   $seq->add_SeqFeature($_) foreach @sfs;
   return;
}


# _check_order_is_consistent($strand,$ranges) RETURNS BOOL
#
# note: the value of this test is moot - there are many valid,
# if unusual cases where it would flag an anomaly. for example
# transpliced genes such as mod(mdg4) in dmel on AE003744, and
# the following spliced gene on NC_001284:
#
#     mRNA            complement(join(20571..20717,21692..22086,190740..190761,
#                     140724..141939,142769..142998))
#                     /gene="nad5"
#                     /note="trans-splicing, RNA editing"
#                     /db_xref="GeneID:814567"
#
# note how the exons are not in order
#  this will flag a level-3 warning, the user of this module
#  can ignore this and deal appropriately with the resulting
#  unordered exons
sub _check_order_is_consistent {
    my $self = shift;

    my $parent_strand = shift; # this does nothing..?
    my @ranges = @_;
    return unless @ranges;
    my $rangestr =
      join(" ",map{sprintf("[%s,%s]",$_->start,$_->end)} @ranges);
    my $strand = $ranges[0]->strand;
    for (my $i=1; $i<@ranges;$i++) {
	if ($ranges[$i]->strand != $strand) {
            $self->problem(1,"inconsistent strands. Trans-spliced gene? Range: $rangestr");
	    return 1; 
            # mixed ranges - autopass
            # some mRNAs have exons on both strands; for
            # example, the dmel mod(mdg4) gene which is
            # trans-spliced (in actual fact two mRNAs)
	}
    }
    my $pass = 1;
    for (my $i=1; $i<@ranges;$i++) {
	my $rangeP = $ranges[$i-1];
	my $range = $ranges[$i];
	    if ($rangeP->start > $range->end) {
                if ($self->seq->is_circular) {
                    # see for example NC_006578.gbk
                    # we make exceptions for circular genomes here.
                    # see Re: [Gmod-ajax] flatfile-to-json.pl error with GFF
                    # 2010-07-26
                }
                else {
                    # failed - but still get one more chance..
                    $pass = 0;
                    $self->problem(2,"Ranges not in correct order. Strange ensembl genbank entry? Range: $rangestr");
                    last;
                }
	    }
    }
    
    if (!$pass) {
        # sometimes (eg ensembl flavour genbank files)
        # exons on reverse strand listed in reverse order
        # eg join(complement(R1),...,complement(Rn))
        # where R1 > R2
        for (my $i=1; $i<@ranges;$i++) {
            my $rangeP = $ranges[$i-1];
            my $range = $ranges[$i];
	    if ($rangeP->end < $range->start) {
                $self->problem(3,"inconsistent order. Range: $rangestr");
                return 0;
	    }
        }
    }
    return 1; # pass
}

# PRIVATE METHOD: _locstr($sf)
#
# returns a location string for a feature; just the outer boundaries
sub _locstr {
    my $self = shift;
    my $sf = shift;
    return
      sprintf("%d..%d", $sf->start, $sf->end);
}

sub iterate_containment_tree {
    my $self = shift;
    my $feature_holder = shift;
    my $sub = shift;
    $sub->($feature_holder);
    my @sfs = $feature_holder->get_SeqFeatures;
    $self->iterate_containment_tree($_) foreach @sfs;
}

sub find_best_pairs {
    my $matrix = shift;
    my $size = shift;
    my $i = shift || 0;

    for (my $j=0; $j < $size; $j++) {
	my $score = $matrix->[$i][$j];
	if (!defined($score)) {
	    next;
	}
	
    }
    
}

1;
