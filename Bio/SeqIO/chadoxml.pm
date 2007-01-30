# $Id$
#
# BioPerl module for Bio::SeqIO::chadoxml
#
# Peili Zhang   <peili@morgan.harvard.edu>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::chadoxml - chadoxml sequence output stream

=head1 SYNOPSIS

It is probably best not to use this object directly, but
rather go through the SeqIO handler system:

    $writer = Bio::SeqIO->new(-file => ">chado.xml",
                              -format => 'chadoxml');

    # assume you already have a Sequence object
    $writer->write_seq($seq_obj);

=head1 DESCRIPTION

This object can transform Bio::Seq objects to chadoxml flat
file databases (for chadoxml DTD, see
http://gmod.cvs.sourceforge.net/gmod/schema/chado/dat/chado.dtd).

This is currently a write-only module.

    $seqio = Bio::SeqIO->new(-file => '>outfile.xml',
                             -format => 'chadoxml');

    # we have a Bio::Seq object $seq which is a gene located on
    # chromosome arm 'X', to be written out to chadoxml
    # before converting to chadoxml, $seq object B<must> be transformed
    # so that all the coordinates in $seq are against the source
    # feature to be passed into Bio::SeqIO::chadoxml->write_seq()
    # -- chromosome arm X in the example below.

    $seqio->write_seq(-seq=>$seq,
                      -genus   => 'Homo',
                      -species => 'sapiens',
                      -seq_so_type=>'gene',
                      -src_feature=>'X',
                      -src_feat_type=>'chromosome_arm',
		                -nounflatten=>1,
                      -is_analysis=>'true',
                      -data_source=>'GenBank');

The chadoxml output of Bio::SeqIO::chadoxml-E<gt>write_seq() method can be
passed to the loader utility in XORT package
(http://gmod.cvs.sourceforge.net/gmod/schema/XMLTools/XORT/)
to be loaded into chado.

This object is currently implemented to work with sequence and
annotation data from whole genome projects deposited in GenBank. It
may not be able to handle all different types of data from all
different sources.

In converting a Bio::Seq object into chadoxml, a top-level feature is
created to represent the object and all sequence features inside the
Bio::Seq object are treated as subfeatures of the top-level
feature. The Bio::SeqIO::chadoxml object calls
Bio::SeqFeature::Tools::Unflattener to unflatten the flat feature list
contained in the subject Bio::Seq object, to build gene model
containment hierarchy conforming to chado central dogma model: gene
--E<gt> mRNA --E<gt> exons and protein.

Destination of data in the subject Bio::Seq object $seq is as following:

	*$seq->display_id:  name of the top-level feature;

	*$seq->accession_number: if defined, uniquename and
				 feature_dbxref of the top-level
				 feature if not defined,
				 $seq->display_id is used as the
				 uniquename of the top-level feature;

	*$seq->molecule: transformed to SO type, used as the feature
			type of the top-level feature if -seq_so_type
			argument is supplied, use the supplied SO type
			as the feature type of the top-level feature;

	*$seq->species: organism of the top-level feature;

	*$seq->seq: residues of the top-level feature;

	*$seq->is_circular, $seq->division: feature_cvterm;

	*$seq->keywords, $seq->desc, comments: featureprop;

	*references: pub and feature_pub;
		medline/pubmed ids: pub_dbxref;
		comments: pubprop;

	*feature "source" span: featureloc for top-level feature;

	*feature "source" db_xref: feature_dbxref for top-level feature;

	*feature "source" other tags: featureprop for top-level feature;

	*subfeature 'symbol' or 'label' tag: feature uniquename, if
                     none of these is present, the chadoxml object
                     generates feature uniquenames as:
                     <gene>-<feature_type>-<span>
                     (e.g. foo-mRNA--1000..3000);

	*gene model: feature_relationship built based on the
                     containment hierarchy;

	*feature span: featureloc;

	*feature accession numbers: feature_dbxref;

	*feature tags (except db_xref, symbol and gene): featureprop;

Things to watch out for:

	*chado schema change: this version works with the chado
                               version tagged chado_1_01 in GMOD CVS.

	*feature uniquenames: especially important if using XORT
                              loader to do incremental load into
                              chado. may need pre-processing of the
                              source data to put the correct
                              uniquenames in place.

	*pub uniquenames: chadoxml->write_seq() has the FlyBase policy
                          on pub uniquenames hard-coded, it assigns
                          pub uniquenames in the following way: for
                          journals and books, use ISBN number; for
                          published papers, use MEDLINE ID; for
                          everything else, use FlyBase unique
                          identifier FBrf#. need to modify the code to
                          implement your policy. look for the comments
                          in the code.

	*for pubs possibly existing in chado but with no knowledge of
         its uniquename:put "op" as "match", then need to run the
                        output chadoxml through a special filter that
                        talks to chado database and tries to find the
                        pub by matching with the provided information
                        instead of looking up by the unique key. after
                        matching, the filter also resets the "match"
                        operation to either "force" (default), or
                        "lookup", or "insert", or "update". the
                        "match" operation is for a special FlyBase use
                        case. please modify to work according to your
                        rules.

	*chado initialization for loading:

		cv & cvterm: in the output chadoxml, all cv's and
                             cvterm's are lookup only. Therefore,
                             before using XORT loader to load the
                             output into chado, chado must be
                             pre-loaded with all necessary CVs and
                             CVterms, including "SO" , "property
                             type", "relationship type", "pub type",
                             "pubprop type", "pub relationship type",
                             "sequence topology", "GenBank feature
                             qualifier", "GenBank division". A pub by
                             the uniquename 'nullpub' of type 'null
                             pub' needs to be inserted.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.
Bug reports can be submitted via the web:

  http://bugzilla.bioperl.org

=head1 AUTHOR - Peili Zhang

Email peili@morgan.harvard.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::chadoxml;
use strict;
use English;

use XML::Writer;
use IO::File;
use IO::Handle;
use Bio::Seq;
use Bio::Seq::RichSeq;
use Bio::SeqIO::FTHelper;
use Bio::Species;
use Bio::Seq::SeqFactory;
use Bio::Factory::SequenceStreamI;
use Bio::SeqFeature::Generic;
use Bio::Annotation::Collection;
use Bio::Annotation::Comment;
use Bio::Annotation::Reference;
use Bio::Annotation::DBLink;
use Bio::SeqFeature::Tools::Unflattener;

#global variables
undef(my %finaldatahash); #data from Bio::Seq object stored in a hash
undef(my %datahash); #data from Bio::Seq object stored in a hash

my $chadotables = 'feature featureprop feature_relationship featureloc feature_cvterm cvterm cv feature_pub pub pub_dbxref pub_author author pub_relationship pubprop feature_dbxref dbxref db';

my %fkey = (
	"cvterm.cv_id"			=> "cv",
	"dbxref.db_id"			=> "db",
	"feature.type_id" 		=> "cvterm",
	"feature.organism_id" 		=> "organism",
	"feature.dbxref_id" 		=> "dbxref",
	"featureprop.type_id" 		=> "cvterm",
	"feature_pub.pub_id" 		=> "pub",
	"feature_cvterm.cvterm_id"	=> "cvterm",
	"feature_cvterm.pub_id"		=> "pub",
	"feature_dbxref.dbxref_id"	=> "dbxref",
	"feature_relationship.object_id"	=> "feature",
	"feature_relationship.subject_id"	=> "feature",
	"feature_relationship.type_id" 	=> "cvterm",
	"featureloc.srcfeature_id"	=> "feature",
	"pub.type_id"			=> "cvterm",
	"pub_dbxref.dbxref_id"		=> "dbxref",
	"pub_author.author_id"		=> "author",
	"pub_relationship.obj_pub_id"	=> "pub",
	"pub_relationship.subj_pub_id"	=> "pub",
	"pub_relationship.type_id"	=> "cvterm",
	"pubprop.type_id"		=> "cvterm",
);

my %feattype_args2so = (
	"aberr"				=> "aberration_junction",
#	"conflict"			=> "sequence_difference",
#	"polyA_signal"			=> "polyA_signal_sequence",
	"variation"			=> "sequence_variant",
	"mutation1"			=> "point_mutation",		#for single-base mutation
	"mutation2"			=> "sequence_variant",		#for multi-base mutation
	"rescue"			=> "rescue_fragment",
#	"rfrag"				=> "restriction_fragment",
	"protein_bind"			=> "protein_binding_site",
	"misc_feature"			=> "region",
#	"prim_transcript"		=> "primary_transcript",
	"CDS"				=> "protein",
	"reg_element"			=> "regulatory_region",
	"seq_variant"			=> "sequence_variant",
	"mat_peptide"			=> "mature_peptide",
	"sig_peptide"			=> "signal_peptide",
);

undef(my %organism);

use base qw(Bio::SeqIO);

sub _initialize {

    my($self,@args) = @_;

    $self->SUPER::_initialize(@args);
    unless( defined $self->sequence_factory ) {
        $self->sequence_factory(new Bio::Seq::SeqFactory
                                (-verbose => $self->verbose(),
                                 -type => 'Bio::Seq::RichSeq'));
    }

}

=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq(-seq=>$seq, -seq_so_type=>$seqSOtype,
			      -src_feature=>$srcfeature,
			      -src_feat_type=>$srcfeattype,
			      -nounflatten=>0 or 1,
			      -is_analysis=>'true' or 'false',
			      -data_source=>$datasource)
 Function: writes the $seq object (must be seq) into chadoxml.
	   Current implementation:
	   1. for non-mRNA records,
	   a top-level feature of type $seq->alphabet is
	   generated for the whole GenBank record, features listed
           are unflattened for DNA records to build gene model
	   feature graph, and for the other types of records all
	   features in $seq are treated as subfeatures of the top-level
	   feature.
	   2. for mRNA records,
	   if a 'gene' feature is present, it B<must> have a /symbol
	   or /label tag to contain the uniquename of the gene. a top-
	   level feature of type 'gene' is generated. the mRNA is written
	   as a subfeature of the top-level gene feature, and the other
	   sequence features listed in $seq are treated as subfeatures
	   of the mRNA feature.
 Returns : 1 for success and 0 for error


 Args     : A Bio::Seq object $seq, optional $seqSOtype, $srcfeature,
	         $srcfeattype, $nounflatten, $is_analysis and $data_source.
           when $srcfeature (a string, the uniquename of the source
           feature) is given, the location and strand information of
           the top-level feature against the source feature will be
           derived from the sequence feature called 'source' of the
           $seq object, a featureloc record is generated for the top
           -level feature on $srcfeature. when $srcfeature is given,
           $srcfeattype must also be present. All feature coordinates
           in $seq should be against $srcfeature.  $seqSOtype is the
           optional SO term to use as the type of the top-level feature.
           For example, a GenBank data file for a Drosophila melanogaster
           genome scaffold has the molecule type of "DNA", when
           converting to chadoxml, a $seqSOtype argument of
           "golden_path_region" can be supplied to save the scaffold
           as a feature of type "golden_path_region" in chadoxml, instead
           of "DNA".  a feature with primary tag of 'source' must be
           present in the sequence feature list of $seq, to decribe the
           whole sequence record.


=cut

sub write_seq {
	my $usage = <<EOUSAGE;
Bio::SeqIO::chadoxml->write_seq()
Usage   : \$stream->write_seq(-seq=>\$seq,
			      -seq_so_type=>\$SOtype,
			      -src_feature=>\$srcfeature,
			      -src_feat_type=>\$srcfeattype,
			      -nounflatten=>0 or 1,
                              -is_analysis=>'true' or 'false',
                              -data_source=>\$datasource)
Args    : \$seq		: a Bio::Seq object
	  \$SOtype	: the SO term to use as the feature type of
	                  the \$seq record, optional
	  \$srcfeature	: unique name of the source feature, a string
			  containing at least one alphabetical letter
			  (a-z, A-Z), optional
	  \$srcfeattype	: feature type of \$srcfeature. one of SO terms.
			  optional
	  when \$srcfeature is given, \$srcfeattype becomes mandatory,
	  \$datasource	: source of the sequence annotation data,
			  e.g. 'GenBank' or 'GFF'.
EOUSAGE

	my ($self,@args) = @_;

	my ($seq, $seq_so_type, $srcfeature, $srcfeattype, $nounflatten, $isanalysis, $datasource, $genus, $species) =
	   $self->_rearrange([qw(SEQ
				 SEQ_SO_TYPE
				 SRC_FEATURE
				 SRC_FEAT_TYPE
				 NOUNFLATTEN
				 IS_ANALYSIS
				 DATA_SOURCE
                                 GENUS
                                 SPECIES
				 )],
			      @args);
	#print "$seq_so_type, $srcfeature, $srcfeattype\n";

	if( !defined $seq ) {
	    $self->throw("Attempting to write with no seq!");
	}

	if( ! ref $seq || ! $seq->isa('Bio::Seq::RichSeqI') ) {
	    $self->warn(" $seq is not a RichSeqI compliant module. Attempting to dump, but may fail!");
	}

	#$srcfeature, when provided, should contain at least one alphabetical letter
	if (defined $srcfeature)
	{
	    if ($srcfeature =~ /[a-zA-Z]/)
	    {
		chomp($srcfeature);
	    } else {
		$self->throw( $usage );
	    }

	    #check for mandatory $srcfeattype
	    if (! defined $srcfeattype)
	    {
		$self->throw( $usage );
		#$srcfeattype must be a string of non-whitespace characters
	    } else {
		if ($srcfeattype =~ /\S+/) {
		    chomp($srcfeattype);
		} else {
		    $self->throw( $usage );
		}
	    }
	}

	# variables local to write_seq()
        my $div = undef;
	my $hkey = undef;
	undef(my @top_featureprops);
	my $name = $seq->display_id if $seq->can('display_id');
        $name = $seq->display_name  if $seq->can('display_name');
	undef(my @feature_cvterms);
	undef(my %sthash);
	undef(my %dvhash);
	undef(my %h1);
	undef(my %h2);
	my $temp = undef;
	my $ann = undef;
	undef(my @references);
	undef(my @feature_pubs);
	my $ref = undef;
	my $location = undef;
	my $fbrf = undef;
	my $journal = undef;
	my $issue = undef;
	my $volume = undef;
	my $volumeissue = undef;
	my $pages = undef;
	my $year = undef;
	my $pubtype = undef;
#	my $miniref= undef;
	my $uniquename = undef;
	my $refhash = undef;
	my $feat = undef;
	my $tag = undef;
	my $tag_cv = undef;
	my $ftype = undef;
	my $subfeatcnt = undef;
	undef(my @top_featrels);
	undef (my %srcfhash);

	local($^W) = 0; # supressing warnings about uninitialized fields.

	if ($seq->can('accession_number') && defined $seq->accession_number && $seq->accession_number ne 'unknown') {
		$uniquename = $seq->accession_number;
	} elsif ($seq->can('accession') && defined $seq->accession && $seq->accession ne 'unknown') {
		$uniquename = $seq->accession;
	} else {
		$uniquename = $name;
	}
        my $len = $seq->length();
	if ($len == 0) {
		$len = undef;
	}

	undef(my $gb_type);
	if (!$seq->can('molecule') || ! defined ($gb_type = $seq->molecule()) ) {
		$gb_type = $seq->can('alphabet') ? $seq->alphabet : 'DNA';
	}
	$gb_type = 'DNA' if $ftype eq 'dna';
	$gb_type = 'RNA' if $ftype eq 'rna';

	if (defined $seq_so_type) {
		$ftype = $seq_so_type;
	}
	else {
		$ftype = $gb_type;
	}

	my %ftype_hash = ( "name" => $ftype, "cv_id" => {"name" => 'SO'});

        if ($species) {
            %organism = ("genus"=>$genus, "species" => $species);
        }
        else {
	    my $spec = $seq->species();
	    if (!defined $spec) {
		$self->throw("$seq does not know what organism it is from, which is required by chado. cannot proceed!\n");
	    } else {
		%organism = ("genus"=>$spec->genus(), "species" => $spec->species());
	    }
        }

        my $residues = $seq->seq || '';

	#set is_analysis flag for gene model features
	undef(my $isanal);
	if ($ftype eq 'gene' || $ftype eq 'mRNA' || $ftype eq 'exon' || $ftype eq 'protein') {
		$isanal = $isanalysis;
		$isanal = 'false' if !defined $isanal;
	}

	%datahash = (
		"name"		=> $name,
		"uniquename"	=> $uniquename,
		"seqlen"	=> $len,
		"residues"	=> $residues,
		"type_id"	=> \%ftype_hash,
		"organism_id"	=> \%organism,
		"is_analysis"	=> $isanal,
		);

	#if $srcfeature is not given, use the Bio::Seq object itself as the srcfeature for featureloc's
	if (!defined $srcfeature) {
		$srcfeature = $uniquename;
		$srcfeattype = $ftype;
	}

	#default data source is 'GenBank'
	if (!defined $datasource) {
		$datasource = 'GenBank';
	}

	if ($datasource =~ /GenBank/i) {
		#sequence topology as feature_cvterm
		if ($seq->can('is_circular') && $seq->is_circular) {
			%sthash = (
				"cvterm_id"	=> {'name' => 'circular',
						    'cv_id' => {
							'name' => 'sequence topology',
						    },
						},
				   "pub_id"	=> {'uniquename' => 'nullpub',
						    'type_id' => {
							'name' => 'null pub',
							'cv_id' => {
							    'name'=> 'pub type',
							},
						    },
						},
				);
		} else {
			%sthash = (
				"cvterm_id"	=> { 'name' => 'linear',
						     'cv_id' => {
							 'name' => 'sequence topology',
						     }
						 },
				"pub_id"	=> {'uniquename' => 'nullpub',
						    'type_id' => {
							'name' => 'null pub',
							'cv_id' => {
							    'name'=> 'pub type',
							},
						    },
						},
				   );
		}
		push(@feature_cvterms, \%sthash);

		#division as feature_cvterm
       		if ($seq->can('division') && defined $seq->division()) {
       		 	$div = $seq->division();
			%dvhash = (
				"cvterm_id"	=> {'name' => $div,
						    'cv_id' => {
							'name' => 'GenBank division'}},
				"pub_id"	=> {'uniquename' => 'nullpub',
						    'type_id' => {
							'name' => 'null pub',
							'cv_id' => {
							    'name'=> 'pub type'},
						    	}},
				);
			push(@feature_cvterms, \%dvhash);
		}

		$datahash{'feature_cvterm'} = \@feature_cvterms;
	}

	#featureprop's
	#DEFINITION
	if ($seq->can('desc') && defined $seq->desc()) {
		$temp = $seq->desc();

		my %prophash = (
			"type_id" 	=> {'name' => 'description',
					    'cv_id' => {
						'name' => 'property type'}},
			"value"		=> $temp,
			);

		push(@top_featureprops, \%prophash);
        }

	#KEYWORDS
	if ($seq->can('keywords')) {
	    $temp = $seq->keywords();

	    if (defined $temp && $temp ne '.' && $temp ne '') {
		my %prophash = (
				"type_id"       => {'name' => 'keywords',
						    'cv_id' =>
						    {'name' => 'property type'}
						},
				"value"          => $temp,
                        	);

		push(@top_featureprops, \%prophash);
	    }
        }

	#COMMENT
	if ($seq->can('annotation')) {
		$ann = $seq->annotation();
		foreach my $comment ($ann->get_Annotations('comment')) {
			$temp = $comment->as_text();
			#print "fcomment: $temp\n";
			my %prophash = (
				"type_id"	=> {'name' => 'comment',
						    'cv_id' =>
						    {'name' => 'property type'}},
				"value"		=> $temp,
				);

			push(@top_featureprops, \%prophash);
		}
	}

	#accession and version as feature_dbxref
	my @top_dbxrefs = ();
	if ($seq->can('accession_number') && defined $seq->accession_number && $seq->accession_number ne 'unknown') {
	    my $db = $self->_guess_acc_db($seq, $seq->accession_number);
	    my %acchash = (
			   "db_id"	=> {'name' => $db},
			   "accession"	=> $seq->accession_number,
			   "version"	=> $seq->seq_version,
			   );
	    my %fdbx = ('dbxref_id' => \%acchash);
	    push(@top_dbxrefs, \%fdbx);
	}

	if( $seq->isa('Bio::Seq::RichSeqI') && defined $seq->get_secondary_accessions() ) {
		my @secacc = $seq->get_secondary_accessions();
		my $acc;
		foreach $acc (@secacc) {
			my %acchash = (
				"db_id"         => {'name' => 'GB'},
				"accession"	=> $acc,
				);
			my %fdbx = ('dbxref_id' => \%acchash);
			push(@top_dbxrefs, \%fdbx);
		}
	}

	#GI number
	if( $seq->isa('Bio::Seq::RichSeqI') && defined ($seq->pid)) {
		my $id = $seq->pid;
		#print "reftype: ", ref($id), "\n";

		#if (ref($id) eq 'HASH') {
		my %acchash = (
			"db_id"		=> {'name' => 'GI'},
			"accession"	=> $id,
			);
		my %fdbx = ('dbxref_id' => \%acchash);
		push (@top_dbxrefs, \%fdbx);
	}

	#REFERENCES as feature_pub
	if (defined $ann) {
	    #get the references
	    @references = $ann->get_Annotations('reference');
	    foreach $ref (@references) {
		undef(my %pubhash);
		$refhash = $ref->hash_tree();
		$location = $ref->location || $refhash->{'location'};
		#print "location: $location\n";

		#get FBrf#, special for FlyBase SEAN loading
		if (index($location, ' ==') >= 0) {
		    $location =~ /\s==/;
				#print "match: $MATCH\n";
				#print "prematch: $PREMATCH\n";
				#print "postmatch: $POSTMATCH\n";
		    $fbrf = $PREMATCH;
		    $location = $POSTMATCH;
		    $location =~ s/^\s//;
		}

		#print "location: $location\n";
		#unpublished reference
		if ($location =~ /Unpublished/) {
		    $pubtype = 'unpublished';
		    %pubhash = (
				"title"		=> $ref->title || $refhash->{'title'},
				#"miniref"	=> substr($location, 0, 255),
				#"uniquename"	=> $fbrf,
				"type_id"	=> {'name' => $pubtype, 'cv_id' => {'name' =>'pub type'}}
				);
		}
		#submitted
		elsif ($location =~ /Submitted/) {
		    $pubtype = 'submitted';

		    %pubhash = (
				"title"		=> $ref->title || $refhash->{'title'},
				#"miniref"	=> substr($location, 0, 255),
				#"uniquename"	=> $fbrf,
				"type_id"	=> {'name' => $pubtype, 'cv_id' => {'name' =>'pub type'}}
				);

		    undef(my $pyear);
		    $pyear = $self->_getSubmitYear($location);
		    if (defined $pyear) {
			$pubhash{'pyear'} = $pyear;
		    }
		}

		#published journal paper
		elsif ($location =~ /\D+\s\d+\s\((\d+|\d+-\d+)\),\s(\d+-\d+|\d+--\d+)\s\(\d\d\d\d\)$/) {
		    $pubtype = 'paper';

				#parse location to get journal, volume, issue, pages & year
		    $location =~ /\(\d\d\d\d\)$/;

		    $year = $MATCH;
		    my $stuff = $PREMATCH;
		    $year =~ s/\(//; #remove the leading parenthesis
		    $year =~ s/\)//; #remove the trailing parenthesis

		    $stuff =~ /,\s(\d+-\d+|\d+--\d+)\s$/;

		    $pages = $MATCH;
		    $stuff = $PREMATCH;
		    $pages =~ s/^, //; #remove the leading comma and space
		    $pages =~ s/ $//; #remove the last space

		    $stuff =~ /\s\d+\s\((\d+|\d+-\d+)\)$/;

		    $volumeissue = $MATCH;
		    $journal = $PREMATCH;
		    $volumeissue =~ s/^ //; #remove the leading space
		    $volumeissue =~ /\((\d+|\d+-\d+)\)$/;
		    $issue = $MATCH;
		    $volume = $PREMATCH;
		    $issue =~ s/^\(//; #remove the leading parentheses
		    $issue =~ s/\)$//; #remove the last parentheses
		    $volume =~ s/^\s//;	#remove the leading space
		    $volume =~ s/\s$//;	#remove the last space

		    %pubhash = (
				"title"		=> $ref->title || $refhash->{'title'},
				"volume"	=> $volume,
				"issue"		=> $issue,
				"pyear"		=> $year,
				"pages"		=> $pages,
				#"miniref"	=> substr($location, 0, 255),
				#"miniref"	=> ' ',
				#"uniquename"	=> $fbrf,
				"type_id"	=> {'name' => $pubtype, 'cv_id' => {'name' =>'pub type'}},
				"pub_relationship" => {
				    'obj_pub_id' => {
					'uniquename' => $journal,
					'title' => $journal,
					#'miniref' => substr($journal, 0, 255),
					'type_id' =>{'name' => 'journal',
						     'cv_id' =>
						     {'name' => 'pub type'
						      },
						 },
						     #'pubprop' =>{'value'=> $journal,
						     #	     'type_id'=>{'name' => 'abbreviation', 'cv_id' => {'name' => 'pubprop type'}},
						     #	    },
						 },
					   'type_id' => {
					       'name' => 'published_in',
					       'cv_id' => {
						   'name' => 'pub relationship type'},
					   },
				},
				);
		}

		#other references
		else {
		    $pubtype = 'other';
		    %pubhash = (
				"title"		=> $ref->title || $refhash->{'title'},
				#"miniref"	=> $fbrf,
				"type_id"	=> {
				    'name' => $pubtype,
				    'cv_id' => {'name' =>'pub type'}
				}
				);
		}

		#pub_author
		my $autref = $self->_getRefAuthors($ref);
		if (defined $autref) {
		    $pubhash{'pub_author'} = $autref;
		}
		# if no author and is type 'submitted' and has submitter address, use the first 100 characters of submitter address as the author lastname.
		else {
		    if ($pubtype eq 'submitted') {
			my $autref = $self->_getSubmitAddr($ref);
			if (defined $autref) {
			    $pubhash{'pub_author'} = $autref;
			}
		    }
		}

		#$ref->comment as pubprop
		#print "ref comment: ", $ref->comment, "\n";
		#print "ref comment: ", $refhash->{'comment'}, "\n";
		if (defined $ref->comment || defined $refhash->{'comment'}) {
		    my $comnt = $ref->comment || $refhash->{'comment'};
				#print "remark: ", $comnt, "\n";
		    $pubhash{'pubprop'} = {
			"type_id"       => {'name' => 'comment', 'cv_id' => {'name' => 'pubprop type'}},
			"value"		=> $comnt,
		    };
		}

		#pub_dbxref
		undef(my @pub_dbxrefs);
		if (defined $fbrf) {
		    push(@pub_dbxrefs, {dbxref_id => {accession => $fbrf, db_id => {'name' => 'FlyBase'}}});
		}
		if (defined ($temp = $ref->medline)) {
		    push(@pub_dbxrefs, {dbxref_id => {accession => $temp, db_id => {'name' => 'MEDLINE'}}});
				#use medline # as the pub's uniquename
		    $pubhash{'uniquename'} = $temp;
		}
		if (defined ($temp = $ref->pubmed)) {
		    push(@pub_dbxrefs, {dbxref_id => {accession => $temp, db_id => {'name' => 'PUBMED'}}});
		}
		$pubhash{'pub_dbxref'} = \@pub_dbxrefs;

		#if the pub uniquename is not defined or blank, put its FBrf# as its uniquename
		#this is unique to FlyBase
		#USERS OF THIS MODULE: PLEASE MODIFY HERE TO IMPLEMENT YOUR POLICY
		# ON PUB UNIQUENAME!!!
		if (!defined $pubhash{'uniquename'} || $pubhash{'uniquename'} eq '') {
		    if (defined $fbrf) {
			$pubhash{'uniquename'} = $fbrf;
		    }
				#else {
				#	$pubhash{'uniquename'} = $self->_CreatePubUname($ref);
				#}
		}

		#add to collection of references
		#if the pub covers the entire sequence of the top-level feature, add it to feature_pubs
		if (($ref->start == 1 && $ref->end == $len) || (!defined $ref->start && !defined $ref->end)) {
		    push(@feature_pubs, {"pub_id" => \%pubhash});
		}
		#the pub is about a sub-sequence of the top-level feature
		#create a feature for the sub-sequence and add pub as its feature_pub
		#featureloc of this sub-sequence is against the top-level feature, in interbase coordinates.
		else {
		    my %parf = (
				'uniquename'	=> $uniquename . ':' . $ref->start . "\.\." . $ref->end,
				'organism_id'	=>\%organism,
				'type_id'	=>{'name' =>'region', 'cv_id' => {'name' => 'SO'}},
				);
		    my %parfsrcf = (
				    'uniquename'	=> $uniquename,
				    'organism_id'	=>\%organism,
				    );
		    my %parfloc = (
				   'srcfeature_id'	=> \%parfsrcf,
				   'fmin'		=> $ref->start - 1,
				   'fmax'		=> $ref->end,
				   );
		    $parf{'featureloc'} = \%parfloc;
		    $parf{'feature_pub'} = {'pub_id' => \%pubhash};
		    my %ffr = (
			       'subject_id'	=> \%parf,
			       'type_id'		=> { 'name' => 'partof', 'cv_id' => { 'name' => 'relationship type'}},
			       );
		    push(@top_featrels, \%ffr);
		}
	    }
	    $datahash{'feature_pub'} = \@feature_pubs;
	}

	##construct srcfeature hash for use in featureloc
	if (defined $srcfeature) {
		%srcfhash = ('uniquename' 	=> $srcfeature,
				'organism_id'   => \%organism,
				'type_id' 	=> {'name' => $srcfeattype, 'cv_id' => {'name' => 'SO'}},
			);

	#	my %fr = (
	#		"object_id"	=> \%srcfhash,
	#		"type_id"	=> { 'name' => 'partof', 'cv_id' => { 'name' => 'relationship type'}},
	#		);

	#	push (@top_featrels, \%fr);
	}

	#unflatten the seq features in $seq if $seq is a gene or a DNA sequence
	if (($gb_type eq 'gene' || $gb_type eq 'DNA') &&
	    !$nounflatten) {
		my $u = Bio::SeqFeature::Tools::Unflattener->new;
		$u->unflatten_seq(-seq=>$seq, -use_magic=>1);
	}

	my @top_sfs = $seq->get_SeqFeatures;
	#print $#top_sfs, "\n";

	#SUBFEATURES

	if ($datasource =~ /GenBank/i) {
		$tag_cv = 'GenBank feature qualifier';
	} elsif ($datasource =~ /GFF/i) {
		$tag_cv = 'GFF tag';
	} else {
		$tag_cv = 'property type';
	}

	my $si = 0;
	foreach $feat (@top_sfs) {
		#$feat = $top_sfs[$si];
		#print "si: $si\n";
		my $prim_tag = $feat->primary_tag;
		#print $prim_tag, "\n";

		# get all qualifiers of the 'source' feature, load these as top_featureprops of the top level feature
		if ($prim_tag eq 'source') {
			foreach $tag ($feat->all_tags()) {
				#db_xref
				if ($tag eq 'db_xref')   {
					my @t1 = $feat->each_tag_value($tag);
					foreach $temp (@t1) {
					   $temp =~ /:/;
					   my $db = $PREMATCH;
					   my $xref = $POSTMATCH;
					   my %acchash = (
						"db_id"		=> {'name' => $db},
						"accession"	=> $xref,
						);
					   my %fdbx = ('dbxref_id' => \%acchash);
					   push (@top_dbxrefs, \%fdbx);
					}
				#other tags as featureprops
				} elsif ($tag ne 'gene') {
					my %prophash = undef;
					%prophash = (
                        			"type_id"       => {'name' => $tag, 'cv_id' => {'name' => $tag_cv}},
						"value"		=> join(' ',$feat->each_tag_value($tag)),
						);
					push(@top_featureprops, \%prophash);
				}
			}

			#featureloc for the top-level feature
			my $fmin = undef;
			my $fmax = undef;
			my $strand = undef;
			my %fl = undef;

			$fmin = $feat->start - 1;
			$fmax = $feat->end;
			$strand = $feat->strand;

			%fl = (
				"srcfeature_id"	=> \%srcfhash,
				"fmin"		=> $fmin,
				"fmax"		=> $fmax,
				"strand"	=> $strand,
				);

			$datahash{'featureloc'} = \%fl;

			#delete 'source' feature from @top_sfs
			splice(@top_sfs, $si, 1);
		}
		$si ++;
	#close loop over top_sfs
	}

	#the top-level features other than 'source'
	foreach $feat (@top_sfs) {
		#print $feat->primary_tag, "\n";

		my $r = $self->_subfeat2featrelhash($name, $ftype, $feat, \%srcfhash, $tag_cv, $isanalysis);

		if (!($ftype eq 'mRNA' && $feat->primary_tag eq 'gene')) {
			my %fr = %$r;
			push(@top_featrels, \%fr);
		} else {
			%finaldatahash = %$r;
		}
	}

	if (@top_dbxrefs) {
		$datahash{'feature_dbxref'} = \@top_dbxrefs;
	}

	if (@top_featureprops) {
		$datahash{'featureprop'} = \@top_featureprops;
	}

	if (@top_featrels) {
		$datahash{'feature_relationship'} = \@top_featrels;
	}

	if ($ftype eq 'mRNA' && %finaldatahash) {
		$finaldatahash{'feature_relationship'} = {
						'subject_id'	=> \%datahash,
						'type_id'	=> { 'name' => 'partof', 'cv_id' => { 'name' => 'relationship type'}},
							 };
	} else {
		%finaldatahash = %datahash;
	}

	my $mainTag = 'feature';
	$self->_hash2xml(undef, $mainTag, \%finaldatahash);

	return 1;
}

sub _hash2xml {
    my $self = shift;
    my $isMatch = undef;
    $isMatch = shift;
    my $ult = shift;
    my $ref = shift;
    my %mh = %$ref;
    undef(my $writer);
    $writer = shift if (@_);
    my $key;
    my $v;
    my $sh;
    my $xx;
    my $yy;
    my $nt;
    my $ntref;
    my $output;
    my $root = shift if (@_);
    #print "ult: $ult\n";
    if (!defined $writer) {
	$root = 1;
	$writer = new XML::Writer(OUTPUT => $self->_fh,
				  DATA_MODE => 1,
				  DATA_INDENT => 3);

	#print header
	$writer->xmlDecl("UTF-8");
	$writer->comment("created by Peili Zhang, Flybase, Harvard University");

	#start chadoxml
	$writer->startTag('chado');
    }
    my $temp;
    my %subh = undef;

    #start opeing tag
    #if pub record of type 'journal', form the 'ref' attribute for special pub lookup
    #requires that the journal name itself is also stored as a pubprop record for the journal with value equal
    #to the journal name and type of 'abbreviation'.
    if ($ult eq 'pub' && $mh{'type_id'}->{'name'} eq 'journal') {
	$writer->startTag($ult, 'ref' => $mh{'title'} . ':journal:abbreviation');
    }

    #special pub match if pub uniquename not known
    elsif ($ult eq 'pub' && !defined $mh{'uniquename'}) {
	$writer->startTag($ult, 'op' => 'match');
	#set the match flag, all the sub tags should also have "op"="match"
	$isMatch = 1;
    }

    #if cvterm or cv, lookup only
    elsif (($ult eq 'cvterm') || ($ult eq 'cv')) {
	$writer->startTag($ult, 'op' => 'lookup');
    }

    #if nested tables of match table, match too
    elsif ($isMatch) {
	$writer->startTag($ult, 'op' => 'match');
    }

    else {
	$writer->startTag($ult);
    }

    #first loop to produce xml for all the table columns
    foreach $key (keys %mh)
    {
	#print "key: $key\n";
	$xx = ' ' . $key;
	$yy = $key . ' ';
	if (index($chadotables, $xx) < 0 && index($chadotables, $yy) < 0)
	{
	    if ($isMatch) {
		$writer->startTag($key, 'op' => 'match');
	    } else {
		$writer->startTag($key);
	    }

	    my $x = $ult . '.' . $key;
	    #the column is a foreign key
	    if (defined $fkey{$x})
	    {
		$nt = $fkey{$x};
		$sh = $mh{$key};
		$self->_hash2xml($isMatch, $nt, $sh, $writer, 0);
	    } else
	    {
		#print "$key: $mh{$key}\n";
		$writer->characters($mh{$key});
	    }
	    $writer->endTag($key);
	}
    }

    #second loop to produce xml for all the nested tables
    foreach $key (keys %mh)
    {
	#print "key: $key\n";
	$xx = ' ' . $key;
	$yy = $key . ' ';
	#a nested table
	if (index($chadotables, $xx) > 0 || index($chadotables, $yy) > 0)
	{
	    #$writer->startTag($key);
	    $ntref = $mh{$key};
	    #print "$key: ", ref($ntref), "\n";
	    if (ref($ntref) =~ 'HASH') {
		$self->_hash2xml($isMatch, $key, $ntref, $writer, 0);
	    } elsif (ref($ntref) =~ 'ARRAY') {
		#print "array dim: ", $#$ntref, "\n";
		foreach $ref (@$ntref) {
				#print "\n";
		    $self->_hash2xml($isMatch, $key, $ref, $writer, 0);
		}
	    }
	    #$writer->endTag($key);
	}
    }

    #end tag
    $writer->endTag($ult);

    if ($root == 1) {
	$writer->endTag('chado');
    }
}

sub _guess_acc_db {
	my $self = shift;
	my $seq = shift;
	my $acc = shift;
	#print "acc: $acc\n";

	if ($acc =~ /^NM_\d{6}/ || $acc =~ /^NP_\d{6}/ || $acc =~ /^NT_\d{6}/ || $acc =~ /^NC_\d{6}/) {
		return "RefSeq";
	} elsif ($acc =~ /^XM_\d{6}/ || $acc =~ /^XP_\d{6}/ || $acc =~ /^XR_\d{6}/) {
		return "RefSeq";
	} elsif ($acc =~ /^[a-zA-Z]{1,2}\d{5,6}/) {
		return "GB";
	} elsif ($seq->molecule() eq 'protein' && $acc =~ /^[a-zA-z]\d{5}/) {
		return "PIR";
	} elsif ($seq->molecule() eq 'protein' && $acc =~ /^\d{6,7}[a-zA-Z]/) {
		return "PRF";
	} elsif ($acc =~ /\d+/ && $acc !~ /[a-zA-Z]/) {
		return "LocusID";
	} elsif ($acc =~ /^CG\d+/ || $acc =~ /^FB[a-z][a-z]\d+/) {
		return "FlyBase";
	} else {
		return "unknown";
	}
}

sub _subfeat2featrelhash {
	my $self = shift;
	my $genename = shift;
	my $seqtype = shift;
	my $feat = shift;
	my $r = shift;
	my %srcf = %$r;		#srcfeature hash for featureloc.srcfeature_id
	my $tag_cv = shift;
	my $isanalysis = shift;

	my $prim_tag = $feat->primary_tag;

	my $sfunique = undef;		#subfeature uniquename
	my $sfname = undef;		#subfeature name
	my $sftype = undef;		#subfeature type

	if ($feat->has_tag('symbol')) {
		($sfunique) = $feat->each_tag_value("symbol");
	} elsif ($feat->has_tag('label')) {
		($sfunique) = $feat->each_tag_value("label");
	} else {
		#$self->throw("$prim_tag at " . $feat->start . "\.\." . $feat->end . " does not have symbol or label! To convert into chadoxml, a seq feature must have a /symbol or /label tag holding its unique name.");
		#generate feature unique name as <genename>-<feature-type>-<span>
		$sfunique = $self->_genFeatUniqueName($genename, $feat);
	}

	if ($feat->has_tag('Name')) {
		($sfname) = $feat->each_tag_value("Name");
	}

	#feature type translation
	if (defined $feattype_args2so{$prim_tag}) {
		$sftype = $feattype_args2so{$prim_tag};
	} else {
		$sftype = $prim_tag;
	}

	if ($prim_tag eq 'mutation') {
		if ($feat->start == $feat->end) {
			$sftype = $feattype_args2so{'mutation1'};
		} else {
			$sftype = $feattype_args2so{'mutation2'};
		}
	}

	#set is_analysis flag for gene model features
	undef(my $isanal);
	if ($sftype eq 'gene' || $sftype eq 'mRNA' || $sftype eq 'exon' || $sftype eq 'protein') {
		$isanal = $isanalysis;
	}

	my %sfhash = (
		"name"			=> $sfname,
		"uniquename"		=> $sfunique,
		"organism_id"		=> \%organism,
		"type_id"		=> { 'name' => $sftype, 'cv_id' => { 'name' => 'SO'}},
		"is_analysis"           => $isanal,
		);

	#make a copy of %sfhash for passing to this method when recursively called
	#my %srcfeat = (
        #        "name"                  => $sfname,
        #        "uniquename"            => $sfunique,
        #        "organism_id"           => \%organism,
        #        "type_id"               => { 'name' => $sftype, 'cv_id' => { 'name' => 'SO'}},
        #        );

	#featureloc for subfeatures
	undef(my $sfmin);
	undef(my $sfmax);
	undef(my $is_sfmin_partial);
	undef(my $is_sfmax_partial);
	undef(my $sfstrand);
	$sfmin = $feat->start - 1;
	$sfmax = $feat->end;
	$sfstrand = $feat->strand();

	#if the gene feature in an mRNA record, cannot use its coordinates, omit featureloc
	if ($seqtype eq 'mRNA' && $sftype eq 'gene') {
	} else {
		if ($feat->location->isa('Bio::Location::FuzzyLocationI')) {
			if ($feat->location->start_pos_type() ne 'EXACT') {
				$is_sfmin_partial = 'true';
			}
			if ($feat->location->end_pos_type() ne 'EXACT') {
				$is_sfmax_partial = 'true';
			}
		}

		my %sfl = (
			"srcfeature_id"	=> \%srcf,
			"fmin"		=> $sfmin,
			"is_fmin_partial" => $is_sfmin_partial || '',
			"fmax"		=> $sfmax,
			"is_fmax_partial" => $is_sfmax_partial || '',
			"strand"	=> $sfstrand,
			);

		$sfhash{'featureloc'} = \%sfl;
	}


	#subfeature tags
	undef(my @sfdbxrefs);		#subfeature dbxrefs
	undef(my @sub_featureprops);	#subfeature props
	foreach my $tag ($feat->all_tags()) {
		#feature_dbxref for features
		if ($tag eq 'db_xref')   {
			my @t1 = $feat->each_tag_value($tag);
			#print "# of dbxref: @t1\n";
			for my $temp (@t1) {
			   $temp =~ /:/;
			   my $db = $PREMATCH;
			   my $xref = $POSTMATCH;
			   #print "db: $db; xref: $xref\n";
			   my %acchash = (
				"db_id"		=> {'name' => $db},
				"accession"	=> $xref,
				);
			   my %sfdbx = ('dbxref_id' => \%acchash);
			   push (@sfdbxrefs, \%sfdbx);
			}
		#featureprop for features, excluding GFF Name & Parent tags
		} elsif ($tag ne 'gene' && $tag ne 'symbol' && $tag ne 'Name' && $tag ne 'Parent') {
			foreach my $val ($feat->each_tag_value($tag)) {
				my %prophash = undef;
				%prophash = (
                      			"type_id"       => {'name' => $tag, 'cv_id' => {'name' => $tag_cv}},
					"value"		=> $val,
				);
				push(@sub_featureprops, \%prophash);
			}
		}
	}
	if (@sub_featureprops) {
		$sfhash{'featureprop'} = \@sub_featureprops;
	}
	if (@sfdbxrefs) {
		$sfhash{'feature_dbxref'} = \@sfdbxrefs;
	}

	undef(my @ssfeatrel);
	if ($feat->has_tag('locus_tag')) {
		($genename)= $feat->each_tag_value('locus_tag');
	} elsif ($feat->has_tag('gene')) {
		($genename)= $feat->each_tag_value('gene');
	}

	foreach my $sf ($feat->get_SeqFeatures()) {
		#print $sf->primary_tag, "\n";
		my $rref = $self->_subfeat2featrelhash($genename, $sftype, $sf, \%srcf, $tag_cv, $isanalysis);
		if (defined $rref) {
			push(@ssfeatrel, $rref);
		}
	}

	if (@ssfeatrel) {
		$sfhash{'feature_relationship'} = \@ssfeatrel;
	}

	#subj-obj relationship type
	undef(my $reltypename);
	if ($sftype eq 'protein') {
		$reltypename = 'producedby';
	} else {
		$reltypename = 'partof';
	}

	my %fr = (
		"subject_id"	=> \%sfhash,
		"type_id"		=> { 'name' => $reltypename, 'cv_id' => { 'name' => 'relationship type'}},
		);

	if ($seqtype eq 'mRNA' && $sftype eq 'gene') {
		return \%sfhash;
	} else {
		return \%fr;
	}

}

#generate uniquename for feature as: <genename>-<feature-type>-<span> (foo-mRNA-10..1000)
sub _genFeatUniqueName {
	my $self = shift;
	my $genename = shift;
	my $feat = shift;
	undef(my $uniquename);
	my $ftype = $feat->primary_tag;
	my $start = $feat->start;
	my $end = $feat->end;

	if ($feat->has_tag('locus_tag')) {
		($genename) = $feat->each_tag_value("locus_tag");
	} elsif ($feat->has_tag('gene')) {
		($genename) = $feat->each_tag_value("gene");
	}

	$uniquename = $genename . '-' . $ftype . '-' . $start . "\.\." . $end;

	return $uniquename;
}

#create uniquename for pubs with no medline id and no FBrf#
#use "<authors>, <year>, <type>" as the uniquename (same as miniref)
#<authors> is <sole-author-surname>    if one author,
#  or <first-author-surname> and <second-author-surname>   if two,
#  or <first-author-surname> et al.   if more
#sub _CreatePubUname {
#	my $self = shift;
#	my $pub = shift;
#	undef(my $pubuname);
#
#	return $pubuname;
#}

#get authors of a reference
#returns ref to the array of author hashes
sub _getRefAuthors {
	my $self = shift;
	my $ref = shift;

	my $temp = $ref->authors;
	undef(my @authors);
	undef(my @aut);

	#there are authors
	if ($temp ne '.') {
		if (index($temp, ' and ') > 0) {
			$temp =~ / and /;
			my $lastauthor = $POSTMATCH;
			@authors = split(/\, /, $PREMATCH);
			push (@authors, $lastauthor);
		} else {
			@authors = split(/\, /, $temp);
		}

		my $a;
		my $i = 0;
		foreach $a (@authors) {
			$i ++;
			#parse the author lastname and givennames
			undef(my $last);
			undef(my $given);
			if (index($a, ',') > 0) {	#genbank format, last,f.m.
				($last, $given) = split(/\,/, $a);
			} elsif (index($a, ' ') > 0) {	#embl format, last f.m.
				($last, $given) = split(/ /, $a);
			}
			my %au = (
				'surname'	=> $last,
				'givennames'	=> $given,
				);
			push(@aut, {author_id => \%au, arank => $i});
		}

		return \@aut;
	}

	#no authors, Bio::SeqIO::genbank doesn't pick up 'CONSRTM' line.
	else {
		return;
	}

}

#extract submission year from the citation of the submitted reference
#genbank format for the submitted citation: JOURNAL   Submitted (DD-MON-YYYY) submitter address
sub _getSubmitYear {
    my $self = shift;
    my $citation = shift;

    if ($citation !~ /Submitted/) {
	$self->warn("not citation for a submitted reference. cannot extract submission year.");
	return;
    } else {
	$citation =~ /Submitted \(\d\d-[a-zA-Z]{3}-\d{4}\)/;
	my $a = $MATCH;
	$a =~ /\d{4}/;
	my $year = $MATCH;

	return $year;
    }
}

sub _getSubmitAddr {
    my $self = shift;
    my $ref = shift;
    undef(my %author);

    my $citation = $ref->location;
    if ($citation !~ /Submitted/) {
	$self->warn("not citation for a submitted reference. cannot extract submission year.");
	return;
    } else {
	$citation =~ /Submitted \(\d\d-[a-zA-Z]{3}-\d{4}\)/;
	my $a = $POSTMATCH;
	if (defined $a) {
	    $a =~ s/^\s//;
	    %author = (
		       'author_id'	=> {'surname'	=> substr($a, 0, 100)},
		       );
	    return \%author;
	} else {
	    return;
	}
    }
}

1;
