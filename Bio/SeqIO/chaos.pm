# $Id$
#
# BioPerl module for Bio::SeqIO::chaos
#
# Chris Mungall <cjm@fruitfly.org>
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::chaos - chaos sequence input/output stream

=head1 SYNOPSIS

    #In general you will not want to use this module directly;
    #use the chaosxml format via SeqIO

    $outstream = Bio::SeqIO->new(-file => $filename, -format => 'chaosxml');

    while ( my $seq = $instream->next_seq() ) {
       $outstream->write_seq($seq);
    }

=head1 DESCRIPTION

This is the guts of L<Bio::SeqIO::chaosxml> - please refer to the
documentation for this module

B<CURRENTLY WRITE ONLY>

ChaosXML is an XML mapping of the chado relational database; for more
information, see http://www.fruitfly.org/chaos-xml

chaos can be represented in various syntaxes - XML, S-Expressions or
indented text. You should see the relevant SeqIO file. You will
probably want to use L<Bio::SeqIO::chaosxml>, which is a wrapper to
this module.

=head2 USING STAG OBJECTS

B<non-standard bioperl stuff you dont necessarily need to know follows>

This module (in write mode) is an B<event producer> - it generates XML
events via the L<Data::Stag> module. If you only care about the final
end-product xml, use L<Bio::SeqIO::chaosxml>

You can treat the resulting chaos-xml stream as stag XML objects;
  
    $outstream = Bio::SeqIO->new(-file => $filename, -format => 'chaos');

    while ( my $seq = $instream->next_seq() ) {
       $outstream->write_seq($seq);
    }
    my $chaos = $outstream->handler->stag;
    # stag provides get/set methods for xml elements
    # (these are chaos objects, not bioperl objects)
    my @features = $chaos->get_feature;
    my @feature_relationships = $chaos->get_feature_relationships;
    # stag objects can be queried with functional-programming
    # style queries
    my @features_in_range =
      $chaos->where('feature',
                    sub {
                         my $featureloc = shift->get_featureloc;
                         $featureloc->strand == 1 &&
                         $featureloc->nbeg > 10000 &&
                         $featureloc->nend < 20000;
                    });
    foreach my $feature (@features_in_range) {
      my $featureloc = $feature->get_featureloc;
      printf "%s [%d->%d on %s]\n", 
        $feature->sget_name,
        $featureloc->sget_nbeg,
        $featureloc->sget_end,
        $featureloc->sget_srcfeature_id;
    }

=head1 MODULES REQUIRED

L<Data::Stag>

Downloadable from CPAN; see also http://stag.sourceforge.net

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://www.bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Chris Mungall

Email cjm@fruitfly.org

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::chaos;
use vars qw(@ISA);
use strict;

use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::Species;
use Bio::Seq::SeqFactory;
use Bio::Annotation::Collection;
use Bio::Annotation::Comment;
use Bio::Annotation::Reference;
use Bio::Annotation::DBLink;
use Bio::SeqFeature::Tools::TypeMapper;
use Bio::SeqFeature::Tools::FeatureNamer;
use Bio::SeqFeature::Tools::IDHandler;
use Data::Stag qw(:all);

@ISA = qw(Bio::SeqIO);

our $TM = 'Bio::SeqFeature::Tools::TypeMapper'; 
our $FNAMER = 'Bio::SeqFeature::Tools::FeatureNamer'; 
our $IDH = 'Bio::SeqFeature::Tools::IDHandler'; 

sub _initialize {
    my($self,@args) = @_;
    
    $self->SUPER::_initialize(@args); 
    if( ! defined $self->sequence_factory ) {
	$self->sequence_factory(new Bio::Seq::SeqFactory
				(-verbose => $self->verbose(), 
				 -type => 'Bio::Seq::RichSeq'));
    }
    my $wclass = $self->default_handler_class;
    $self->handler($wclass);
    if ($self->_fh) {
        $self->handler->fh($self->_fh);
    }
    $self->{_end_of_data} = 0;
    $self->_type_by_id_h({});
    my $t = time;
    my $ppt = localtime $t;
    $self->handler->S("chaos");
    $self->handler->ev(chaos_metadata=>[
                                        [chaos_version=>1],
                                        [chaos_flavour=>'bioperl'],
                                        [feature_unique_key=>'feature_id'],
                                        [equiv_chado_release=>'chado_1_01'],
                                        [export_unixtime=>$t],
                                        [export_localtime=>$ppt],
                                        [export_host=>$ENV{HOST}],
                                        [export_user=>$ENV{USER}],
                                        [export_perl5lib=>$ENV{PERL5LIB}],
                                        [export_program=>$0],
                                        [export_module=>'Bio::SeqIO::chaos'],
                                       ]);

    return;
}

sub DESTROY {
    my $self = shift;
    $self->end_of_data();
    $self->SUPER::DESTROY();
}

sub end_of_data {
    my $self = shift;
    return if $self->{_end_of_data};
    $self->{_end_of_data} = 1;
    $self->handler->E("chaos");
}

sub default_handler_class {
    return Data::Stag->makehandler;
} 

=head2 context_namespace

 Title   : context_namespace
 Usage   : $obj->context_namespace($newval)
 Function: 
 Example : 
 Returns : value of context_namespace (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub context_namespace{
    my $self = shift;

    return $self->{'context_namespace'} = shift if @_;
    return $self->{'context_namespace'};
}


=head2 next_seq

 Title   : next_seq
 Usage   : $seq = $stream->next_seq()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq object
 Args    :

=cut

sub next_seq {
    my ($self,@args) = @_;
    my $seq = $self->sequence_factory->create
	(
         #         '-verbose' =>$self->verbose(), 
         #	 %params,
         #	 -seq => $seqc,
         #	 -annotation => $annotation,
         #	 -features => \@features
        );
    return $seq;
}

sub handler {
    my $self = shift;
    $self->{_handler} = shift if @_;
    return $self->{_handler};
}


=head2 write_seq

 Title   : write_seq
 Usage   : $stream->write_seq($seq)
 Function: writes the $seq object (must be seq) to the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::Seq


=cut

sub write_seq {
    my ($self,$seq) = @_;
    
    if( !defined $seq ) {
	$self->throw("Attempting to write with no seq!");
    }
    
    if( ! ref $seq || ! $seq->isa('Bio::SeqI') ) {
	$self->warn(" $seq is not a SeqI compliant module. Attempting to dump, but may fail!");
    }

    # get a handler - must inherit from Data::Stag::BaseHandler;
    my $w = $self->handler;

    # start of data
    ###    $w->S("chaos_block");

    my $seq_chaos_feature_id;
    if ($seq->accession_number) {
        $seq_chaos_feature_id = $seq->accession_number . '.' . ($seq->can('seq_version') ? $seq->seq_version : $seq->version);
    }
    else {
        $seq_chaos_feature_id = $self->get_chaos_feature_id($seq);
    }

#    if ($seq->accession_number eq 'unknown') {
#        $seq_chaos_feature_id = $self->get_chaos_feature_id('contig', $seq);
#    }

    # data structure representing the core sequence for this record
    my $seqnode =
      Data::Stag->new(feature=>[
                                [feature_id=>$seq_chaos_feature_id],
                                [dbxrefstr=>'SEQDB:'.$seq->accession_number],
                                [name=>$seq->display_name],
				[uniquename=>$seq->display_name .'/'. $seq_chaos_feature_id],
                                [residues=>$seq->seq],
                               ]);
    
    # soft properties
    my %prop = ();

    my ($div, $mol);
    my $len = $seq->length();

    if ( $seq->can('division') ) {
	$div=$seq->division;
    } 
    if( !defined $div || ! $div ) { $div = 'UNK'; }

    if( !$seq->can('molecule') || ! defined ($mol = $seq->molecule()) ) {
	$mol = $seq->alphabet || 'DNA';
    }
#    my $seq_type = 'contig';
#    if ($mol eq 'AA') {
#	$seq_type = 'protein';
#    }
#    if ($mol eq 'RNA') {
#	$seq_type = 'cDNA';
#    }
    $seqnode->set_type('databank_entry');
    
    my $circular = 'linear  ';
    $circular = 'circular' if $seq->is_circular;

    # cheeky hack - access symbol table
    no strict 'refs';
    map {
        $prop{$_} = 
          $ {*$_};
    } qw(mol div circular);
    use strict 'refs';

    map {
        $prop{$_} = $seq->$_() if $seq->can($_);
    } qw(desc keywords);

    local($^W) = 0;   # supressing warnings about uninitialized fields.
    
    my $OS;
    # Organism lines
    if (my $spec = $seq->species) {
        my ($species, $genus, @class) = $spec->classification();
        $OS = "$genus $species";

        if (my $ssp = $spec->sub_species) {
            $OS .= " $ssp";
        }
        if( $spec->common_name ) {
	    $OS .= " (".$spec->common_name.")";
        }
    }
    
    # Reference lines
    my $count = 1;
    foreach my $ref ( $seq->annotation->get_Annotations('reference') ) {
        # TODO
    }
    # Comment lines
    
    foreach my $comment ( $seq->annotation->get_Annotations('comment') ) {
        $seqnode->add_featureprop([[type=>'comment'],[value=>$comment->text]]);
    }
    if ($OS) {
	$seqnode->set_organismstr($OS);
        $self->organismstr($OS);
    }

    my @sfs = $seq->get_SeqFeatures;

    # genbank usually includes a 'source' feature - we just
    # migrate the data from this to the actual source feature
    my @sources = grep {$_->primary_tag eq 'source'} @sfs;
    @sfs = grep {$_->primary_tag ne 'source'} @sfs;
    $self->throw(">1 source types") if @sources > 1;
    my $source = shift @sources;
    if ($source) {

	my $tempw = Data::Stag->makehandler;
	$self->write_sf($source, $seq_chaos_feature_id, $tempw);
	my $snode = $tempw->stag;
	$seqnode->add($_->name, $_->data)
	  foreach ($snode->get_featureprop,
		   $snode->get_feature_dbxref);
	
    }
    

    # throw the writer an event
    $w->ev(@$seqnode);

    $seqnode = undef;      # free memory

    # make events for all the features within the record
    foreach my $sf ( @sfs ) {
        $FNAMER->name_feature($sf);
        $FNAMER->name_contained_features($sf);
        $self->write_sf($sf, $seq_chaos_feature_id);
    }

    # data end
    ### $w->E("chaos_block");
    return 1;
}


sub organismstr{
    my $self = shift;

    return $self->{'organismstr'} = shift if @_;
    return $self->{'organismstr'};
}

# maps ID to type
sub _type_by_id_h {
    my $self = shift;
    $self->{_type_by_id_h} = shift if @_;
    return $self->{_type_by_id_h};
}



# ----
# writes a seq feature
# ----

sub write_sf {
    my $self = shift;
    my $sf = shift;
    my $seq_chaos_feature_id = shift;
    my $w = shift || $self->handler;

    my %props =
      map {
          lc($_)=>[$sf->each_tag_value($_)]
      } $sf->all_tags;

    my $loc = $sf->location;
    my $name = $FNAMER->generate_feature_name($sf);
    my $type = $sf->primary_tag;

    # The CDS (eg in a genbank feature) implicitly represents
    # the protein
    $type =~ s/CDS/protein/;

    my @subsfs = $sf->sub_SeqFeature;
    my @locnodes = ();
    my $sid = $loc->is_remote ? $loc->seq_id : $seq_chaos_feature_id;

    my $CREATE_SPLIT_SFS = 0;

    if($CREATE_SPLIT_SFS &&
       $loc->isa("Bio::Location::SplitLocationI") ) {
        # turn splitlocs into subfeatures
        my $n = 1;
        push(@subsfs,
             map {
                 my $ssf =
                   Bio::SeqFeature::Generic->new(

                                                 -start=>$_->start,
                                                 -end=>$_->end,
                                                 -strand=>$_->strand,
                                                 -primary=>$self->subpartof($type),
                                              );
                 if ($_->is_remote) {
                     $ssf->location->is_remote(1);
                     $ssf->location->seq_id($_->seq_id);
                 }
                 $ssf;
               } $loc->each_Location);
    }
    elsif( $loc->isa("Bio::Location::RemoteLocationI") ) {
        # turn splitlocs into subfeatures
        my $n = 1;
        push(@subsfs,
             map {
                 Bio::SeqFeature::Generic->new(
#                                               -name=>$name.'.'.$n++,
                                               -start=>$_->start,
                                               -end=>$_->end,
                                               -strand=>$_->strand,
                                               -primary=>$self->subpartof($type),
                                              )
               } $loc->each_Location);
    }
    else {
        my ($beg, $end, $strand) = $self->bp2ib($loc);
	if (!$strand) {
	    print "($beg, $end, $strand) - no strand\n";
	    use Data::Dumper;
	    print Dumper $sf;
	    die Dumper $loc;
	}
        @locnodes = (
                     [featureloc=>[
                                   [nbeg=>$beg],
                                   [nend=>$end],
                                   [strand=>$strand],
                                   [srcfeature_id=>$sid],
                                   [locgroup=>0],
                                   [rank=>0],
                                  ]
                     ]
                    );
    }
    my $feature_id = $self->get_chaos_feature_id($sf);
    
    # do something with genbank stuff
    my $pid = $props{'protein_id'};
    my $tn = $props{'translation'};
    my @xrefs = @{$props{'db_xref'} || []};
    if ($pid) {
	push(@xrefs, "protein:$pid->[0]");
    }
    
    my $org = $props{organism};
    if (!$org && $self->organismstr) {
        $org = [$self->organismstr];
    }
    my $uname = $name ? $name.'/'.$feature_id : $feature_id;
    if ($org && $name) {
        my $os = $org->[0];
        $os =~ s/\s+/_/g;
        $os =~ s/_*\(.*//;
        $uname = "$os:$name";
    }
    $self->_type_by_id_h->{$feature_id} = $type;
    my $fnode =
      [feature=>[
                 [feature_id=>$feature_id],
                 $name ? ([name=>$name]) : (),
                 [uniquename=>$uname],
                 [type=>$type],
		 $tn ? ([residues=>$tn->[0]], 
			[seqlen=>length($tn->[0])],
			#####[md5checksum=>md5checksum($tn->[0])],
		       ) :(),
		 $org ? ([organismstr=>$org->[0]]) : (),
                 @locnodes,
		 (map {
		     [feature_dbxref=>[
				       [dbxrefstr=>$_]
				      ]
		     ]
		 } @xrefs),
                 (map {
                     my $k = $_;
		     my $rank=0;
                     map { [featureprop=>[[type=>$k],[value=>$_],[rank=>$rank++]]] } @{$props{$k}}
                 } keys %props),
                ]];
    $w->ev(@$fnode);

    my $rank = 0;
    if (@subsfs) {
	# strand is always determined by FIRST feature listed
	# (see genbank entry for trans-spliced mod(mdg4) AE003734)
	my $strand = $subsfs[0];

	# almost all the time, all features are on same strand
	my @sfs_on_main_strand = grep {$_->strand == $strand} @subsfs;
	my @sfs_on_other_strand = grep {$_->strand != $strand} @subsfs;
	
	sort_by_strand($strand, \@sfs_on_main_strand);
	sort_by_strand(0-$strand, \@sfs_on_other_strand);
	@subsfs = (@sfs_on_main_strand, @sfs_on_other_strand);

	foreach my $ssf (@subsfs) {
	    my $ssfid = $self->write_sf($ssf, $sid);
	    #my $rtype = 'part_of';
            my $rtype = 
              $TM->get_relationship_type_by_parent_child($sf,$ssf);
	    if ($ssf->primary_tag eq 'CDS') {
		$rtype = 'produced_by';
	    }
	    $w->ev(feature_relationship=>[
					  [subject_id=>$ssfid],
					  [object_id=>$feature_id],
					  [type=>$rtype],
					  [rank=>$rank++],
					 ]
		  );
	}
    }
    else {
        # parents not stored as 
        my @parent_ids = @{$props{parent} || []};
        foreach my $parent_id (@parent_ids) {
            my $ptype =
              $self->_type_by_id_h->{$parent_id} || 'unknown';
            my $rtype = 
              $TM->get_relationship_type_by_parent_child($ptype,$type);
	    $w->ev(feature_relationship=>[
					  [subject_id=>$feature_id],
					  [object_id=>$parent_id],
					  [type=>$rtype],
					  [rank=>$rank++],
					 ]
		  );
        }
    }
    return $feature_id;
}

sub sort_by_strand {
    my $strand = shift || 1;
    my $sfs = shift;
    @$sfs = sort { ($a->start <=> $b->start) * $strand } @$sfs;
    return;
}

# private;
# an ID for this session that should be
# unique... hmm
sub session_id {
    my $self = shift;
    $self->{_session_id} = shift if @_;
    if (!$self->{_session_id}) {
        $self->{_session_id} = $$.time;
    }
    return $self->{_session_id};
}


sub get_chaos_feature_id {
    my $self = shift;
    my $ob = shift;

    my $id = $ob->can("primary_id") ? $ob->primary_id : undef;
    if ($id) {
        $id = $self->context_namespace ? $self->context_namespace . ":" . $id : $id;
    }
    else {
        if ($ob->isa("Bio::SeqFeatureI")) {
            $id = $IDH->generate_unique_persistent_id($ob);
        }
        else {
            $self->throw("Cannot generate a unique persistent ID for a Seq without either primary_id or accession");
        }
    }
    return $id;
}

# interbase and directional semantics
sub bp2ib {
    my $self = shift;
    my $loc = shift;
    my ($s, $e, $str) = 
      ref($loc) eq "ARRAY" ? (@$loc) : ($loc->start, $loc->end, $loc->strand);
    $s--;
    if ($str < 0) {
        ($s, $e) = ($e, $s);
    }
    return ($s, $e, $str || 1);
}

sub subpartof {
    my $self = shift;
    my $type = 'partof_'.shift;
    $type =~ s/partof_CDS/CDS_exon/;
    $type =~ s/partof_protein/CDS_exon/;
    $type =~ s/partof_\w*RNA/exon/;
    return $type;
}

1;
