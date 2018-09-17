#
# bioperl module for Bio::SeqFeature::Tools::TypeMapper
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

Bio::SeqFeature::Tools::TypeMapper - maps $seq_feature-E<gt>primary_tag

=head1 SYNOPSIS

  use Bio::SeqIO;
  use Bio::SeqFeature::Tools::TypeMapper;

  # first fetch a genbank SeqI object
  $seqio =
    Bio::SeqIO->new(-file=>'AE003644.gbk',
                    -format=>'GenBank');
  $seq = $seqio->next_seq();

  $tm = Bio::SeqFeature::Tools::TypeMapper->new;

  # map all the types in the sequence
  $tm->map_types(-seq=>$seq,
		 {CDS=>'ORF',
		  variation=>sub {
		      my $f = shift;
		      $f->length > 1 ?
			'variation' : 'SNP'
		  },
		 });

   # alternatively, use the hardcoded SO mapping
   $tm->map_types_to_SO(-seq=>$seq);

=head1 DESCRIPTION

This class implements an object for mapping between types; for
example, the types in a genbank feature table, and the types specified
in the Sequence Ontology.

You can specify your own mapping, either as a simple hash index, or by
providing your own subroutines.

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

package Bio::SeqFeature::Tools::TypeMapper;
use strict;

# Object preamble - inherits from Bio::Root::Root

use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : $unflattener = Bio::SeqFeature::Tools::TypeMapper->new();
 Function: constructor
 Example : 
 Returns : a new Bio::SeqFeature::Tools::TypeMapper
 Args    : see below


=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my($typemap) =
	$self->_rearrange([qw(TYPEMAP
			     )],
                          @args);

    $typemap  && $self->typemap($typemap);
    return $self; # success - we hope!
}

=head2 typemap

 Title   : typemap
 Usage   : $obj->typemap($newval)
 Function: 
 Example : 
 Returns : value of typemap (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub typemap{
    my $self = shift;

    return $self->{'typemap'} = shift if @_;
    return $self->{'typemap'};
}

=head2 map_types

 Title   : map_types
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

 dgg: added -undefined => "region" option to produce all valid SO mappings.

=cut

sub map_types{
   my ($self,@args) = @_;

   my($sf, $seq, $type_map, $undefmap) =
     $self->_rearrange([qw(FEATURE
                           SEQ
			   TYPE_MAP
			   UNDEFINED
                          )],
                          @args);
   if (!$sf && !$seq) {
       $self->throw("you need to pass in either -feature or -seq");
   }

   my @sfs = ($sf);
   if ($seq) {
       $seq->isa("Bio::SeqI") || $self->throw("$seq NOT A SeqI");
       @sfs = $seq->get_all_SeqFeatures;
   }
   $type_map = $type_map || $self->typemap; # dgg: was type_map;
   foreach my $sf (@sfs) {

       $sf->isa("Bio::SeqFeatureI") || $self->throw("$sf NOT A SeqFeatureI");
       $sf->isa("Bio::FeatureHolderI") || $self->throw("$sf NOT A FeatureHolderI");

       my $type = $sf->primary_tag;
       my $mtype = $type_map->{$type};
       if ($mtype) {
	   if (ref($mtype)) {
	       if (ref($mtype) eq 'CODE') {
		   $mtype = $mtype->($sf);
	       }
	       else {
		   $self->throw('type_map values must be scalar or CODE ref. You said: '.$mtype.' for type: '.$type);
	       }
	   }
	   elsif ($undefmap && $mtype eq 'undefined') { # dgg
	      $mtype= $undefmap;
           }
	   $sf->primary_tag($mtype);
       }
   }
   return;
}

=head2 map_types_to_SO

 Title   : map_types_to_SO
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

hardcodes the genbank to SO mapping

Based on revision 1.22 of SO

Please see the actual code for the mappings

Taken from

L<http://sequenceontology.org/resources/mapping/FT_SO.txt>

dgg: separated out FT_SO_map for caller changes. Update with:

  open(FTSO,"curl -s http://sequenceontology.org/resources/mapping/FT_SO.txt|");
  while(<FTSO>){
    chomp; ($ft,$so,$sid,$ftdef,$sodef)= split"\t";
    print "     '$ft' => '$so',\n" if($ft && $so && $ftdef);             
  }

=cut

sub FT_SO_map  {
  # $self= shift;
  # note : some of the ft_so mappings are commented out and overriden...
    return {
	"-" => ["located_sequence_feature", "so:0000110"],
	"-10_signal" => ["minus_10_signal", "so:0000175"],
	"-35_signal" => ["minus_35_signal", "so:0000176"],
	"3'utr" => ["three_prime_utr", "so:0000205"],
	"3'clip" => ["three_prime_clip", "so:0000557"],
	"5'utr" => ["five_prime_utr", "so:0000204"],
	"5'clip" => ["five_prime_clip", "so:0000555"],
	"caat_signal" => ["caat_signal", "so:0000172"],
	"cds" => ["cds", "so:0000316"],
	"c_region" => ["undefined", ""],
	"d-loop" => ["d_loop", "so:0000297"],
	"d_segment" => ["d_gene", "so:0000458"],
	"gc_signal" => ["gc_rich_region", "so:0000173"],
	"j_segment" => ["undefined", ""],
	"ltr" => ["long_terminal_repeat", "so:0000286"],
	"n_region" => ["undefined", ""],
	"rbs" => ["ribosome_entry_site", "so:0000139"],
	"sts" => ["sts", "so:0000331"],
	"s_region" => ["undefined", ""],
	"tata_signal" => ["tata_box", "so:0000174"],
	"v_region" => ["undefined", ""],
	"v_segment" => ["undefined", ""],
	"attenuator" => ["attenuator", "so:0000140"],
	"conflict" => ["undefined", ""],
	"enhancer" => ["enhancer", "so:0000165"],
	"exon" => ["exon", "so:0000147"],
	"gap" => ["gap", "so:0000730"],
	"gene" => ["gene", "so:0000704"],
	"idna" => ["idna", "so:0000723"],
	"intron" => ["intron", "so:0000188"],
	"mRNA" => ["mRNA", "so:0000234"],
	"mat_peptide" => ["mature_protein_region", "so:0000419"],
	"mature_peptide" => ["mature_protein_region", "so:0000419"],
	#"misc_RNA" => ["transcript", "so:0000673"],
	"misc_binding" => ["binding_site", "so:0000409"],
	"misc_difference" => ["sequence_difference", "so:0000413"],
	"misc_feature" => ["region", undef],
	"misc_recomb" => ["recombination_feature", "so:0000298"],
	"misc_signal" => ["regulatory_region", "so:0005836"],
	"misc_structure" => ["sequence_secondary_structure", "so:0000002"],
	"modified_base" => ["modified_base_site", "so:0000305"],
	"old_sequence" => ["undefined", ""],
	"operon" => ["operon", "so:0000178"],
	"oriT" => ["origin_of_transfer", "so:0000724"],
	"polya_signal" => ["polyA_signal_sequence", "so:0000551"],
	"polya_site" => ["polyA_site", "so:0000553"],
	"precursor_RNA" => ["primary_transcript", "so:0000185"],
	"prim_transcript" => ["primary_transcript", "so:0000185"],
	"primer_bind" => ["primer_binding_site", "so:0005850"],
	"promoter" => ["promoter", "so:0000167"],
	"protein_bind" => ["protein_binding_site", "so:0000410"],
	"rRNA" => ["rRNA", "so:0000252"],
	"repeat_region" => ["repeat_region", "so:0000657"],
	"repeat_unit" => ["repeat_unit", "so:0000726"],
	"satellite" => ["satellite_dna", "so:0000005"],
	"scRNA" => ["scRNA", "so:0000013"],
	"sig_peptide" => ["signal_peptide", "so:0000418"],
	"snRNA" => ["snRNA", "so:0000274"],
	"snoRNA" => ["snoRNA", "so:0000275"],
	#"source" => ["databank_entry", "so:2000061"],
	"stem_loop" => ["stem_loop", "so:0000313"],
	"tRNA" => ["tRNA", "so:0000253"],
	"terminator" => ["terminator", "so:0000141"],
	"transit_peptide" => ["transit_peptide", "so:0000725"],
	"unsure" => "undefined",
	"variation" => ["sequence_variant", "so:0000109"],

	# manually added 
	## has parent = pseudogene ; dgg
	"pseudomRNA" => ["pseudogenic_transcript", "so:0000516"],
	## from unflattener misc_rna ; dgg
	"pseudotranscript" => ["pseudogenic_transcript", "so:0000516"],
	"pseudoexon" => ["pseudogenic_exon", "so:0000507"],
	"pseudoCDS" => ["pseudogenic_exon", "so:0000507"],
	"pseudomisc_feature" => ["pseudogenic_region", "so:0000462"],
	"pseudointron" => ["pseudogenic_region", "so:0000462"],


	## "undefined" => "region",

	# this is the most generic form for rnas;
	# we always represent the processed form of
	# the transcript
	misc_RNA => ['mature_transcript',"so:0000233"],

	# not sure about this one...
	source=>['contig', "SO:0000149"],

	rep_origin=>['origin_of_replication',"SO:0000296"],

	Protein=>['polypeptide',"SO:0000104"],
    };
#  return {
     #"FT term" => "SO term",
     #"-" => "located_sequence_feature",
     #"-10_signal" => "minus_10_signal",
     #"-35_signal" => "minus_35_signal",
     #"3'UTR" => "three_prime_UTR",
     #"3'clip" => "three_prime_clip",
     #"5'UTR" => "five_prime_UTR",
     #"5'clip" => "five_prime_clip",
     #"CAAT_signal" => "CAAT_signal",
     #"CDS" => "CDS",
     #"C_region" => "undefined",
     #"D-loop" => "D_loop",
     #"D_segment" => "D_gene",
     #"GC_signal" => "GC_rich_region",
     #"J_segment" => "undefined",
     #"LTR" => "long_terminal_repeat",
     #"N_region" => "undefined",
     #"RBS" => "ribosome_entry_site",
     #"STS" => "STS",
     #"S_region" => "undefined",
     #"TATA_signal" => "TATA_box",
     #"V_region" => "undefined",
     #"V_segment" => "undefined",
     #"attenuator" => "attenuator",
     #"conflict" => "undefined",
     #"enhancer" => "enhancer",
     #"exon" => "exon",
     #"gap" => "gap",
     #"gene" => "gene",
     #"iDNA" => "iDNA",
     #"intron" => "intron",
     #"mRNA" => "mRNA",
     #"mat_peptide" => "mature_protein_region",
     #"mature_peptide" => "mature_protein_region",
##                     "misc_RNA" => "transcript",
     #"misc_binding" => "binding_site",
     #"misc_difference" => "sequence_difference",
     #"misc_feature" => "region",
     #"misc_recomb" => "recombination_feature",
     #"misc_signal" => "regulatory_region",
     #"misc_structure" => "sequence_secondary_structure",
     #"modified_base" => "modified_base_site",
     #"old_sequence" => "undefined",
     #"operon" => "operon",
     #"oriT" => "origin_of_transfer",
     #"polyA_signal" => "polyA_signal_sequence",
     #"polyA_site" => "polyA_site",
     #"precursor_RNA" => "primary_transcript",
     #"prim_transcript" => "primary_transcript",
     #"primer_bind" => "primer_binding_site",
     #"promoter" => "promoter",
     #"protein_bind" => "protein_binding_site",
     #"rRNA" => "rRNA",
     #"repeat_region" => "repeat_region",
     #"repeat_unit" => "repeat_unit",
     #"satellite" => "satellite_DNA",
     #"scRNA" => "scRNA",
     #"sig_peptide" => "signal_peptide",
     #"snRNA" => "snRNA",
     #"snoRNA" => "snoRNA",
##                     "source" => "databank_entry",
     #"stem_loop" => "stem_loop",
     #"tRNA" => "tRNA",
     #"terminator" => "terminator",
     #"transit_peptide" => "transit_peptide",
     #"unsure" => "undefined",
     #"variation" => "sequence_variant",

      #"pseudomRNA" => "pseudogenic_transcript", ## has parent = pseudogene ; dgg
      #"pseudotranscript" => "pseudogenic_transcript", ## from Unflattener misc_RNA ; dgg
      #"pseudoexon" => "pseudogenic_exon",
      #"pseudoCDS"  => "pseudogenic_exon",
      #"pseudomisc_feature" => "pseudogenic_region",
      #"pseudointron" => "pseudogenic_region",
      
      ### "undefined" => "region",

      ## this is the most generic form for RNAs;
      ## we always represent the processed form of
      ## the transcript
      #misc_RNA=>'processed_transcript',
      
      ## not sure about this one...
      #source=>'contig',
      
      #rep_origin=>'origin_of_replication',
      
      #Protein=>'protein',
      #};
}

sub map_types_to_SO{
   my ($self,@args) = @_;

   push(@args, (-type_map=> $self->FT_SO_map() ) );
   return $self->map_types(@args);
}

=head2 get_relationship_type_by_parent_child

 Title   : get_relationship_type_by_parent_child
 Usage   : $type = $tm->get_relationship_type_by_parent_child($parent_sf, $child_sf);
 Usage   : $type = $tm->get_relationship_type_by_parent_child('mRNA', 'protein');
 Function: given two features where the parent contains the child,
           will determine what the relationship between them in
 Example :
 Returns : 
 Args    : parent SeqFeature, child SeqFeature OR
           parent type string, child type string OR

bioperl Seq::FeatureHolderI hierarchies are equivalent to unlabeled
graphs (where parent nodes are the containers, and child nodes are the
features being contained). For example, a feature of type mRNA can
contain features of type exon.

Some external representations (eg chadoxml or chaosxml) require that
the edges in the feature relationship graph are labeled. For example,
the type between mRNA and exon would be B<part_of>. Although it
stretches the bioperl notion of containment, we could have a CDS
contained by an mRNA (for example, the
L<Bio::SeqFeature::Tools::Unflattener> module takes genbank records
and makes these kind of links. The relationship here would be
B<produced_by>

In chado speak, the child is the B<subject> feature and the parent is
the B<object> feature

=cut

sub get_relationship_type_by_parent_child {
   my ($self,$parent,$child) = @_;
   $parent = ref($parent) ? $parent->primary_tag : $parent;
   $child = ref($child) ? $child->primary_tag : $child;

   my $type = 'part_of'; # default

   # TODO - do this with metadata, or infer via SO itself

   if (lc($child) eq 'protein') {
       $type = 'derives_from';
   }
   if (lc($child) eq 'polypeptide') {
       $type = 'derives_from';
   }
   return $type;
}


1;
