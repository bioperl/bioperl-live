# $Id$
#
# $GNF: projects/gi/symgene/src/perl/seqproc/Bio/SeqIO/locuslink.pm,v 1.6 2002/11/07 11:14:34 hlapp Exp $
#
# BioPerl module for Bio::SeqIO::locuslink
#
# Cared for by Keith Ching <kching at gnf.org>
#
# Copyright Keith Ching
#
# You may distribute this module under the same terms as perl itself

#
# (c) Keith Ching, kching at gnf.org, 2002.
# (c) GNF, Genomics Institute of the Novartis Research Foundation, 2002.
#
# You may distribute this module under the same terms as perl itself.
# Refer to the Perl Artistic License (see the license accompanying this
# software package, or see http://www.perl.com/language/misc/Artistic.html)
# for the terms under which you may use, modify, and redistribute this module.
# 
# THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::locuslink - DESCRIPTION of Object

=head1 SYNOPSIS

   # don't instantiate directly - instead do
   my $seqio = Bio::SeqIO->new(-format => "locuslink", -file => \STDIN);

=head1 DESCRIPTION

This module parses LocusLink into Bio::SeqI objects with rich
annotation, but no sequence.

The input file has to be in the LL_tmpl format - the tabular format
will not work.

The way the current implementation populates the object is rather a
draft work than a finished work of art. Note that at this stage the
locuslink entries cannot be round-tripped, because the parser loses
certain information. For instance, most of the alternative transcript
descriptions are not retained. The parser also misses any element
that deals with visual representation (e.g., 'button') except for the
URLs. Almost all of the pieces of the annotation are kept in the
L<Bio::Annotation::Collection> object.

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

=head1 AUTHOR - Keith Ching

Email kching at gnf.org

Describe contact details here

=head1 CONTRIBUTORS

Hilmar Lapp, hlapp at gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::SeqIO::locuslink;

use strict;
use vars qw(@ISA);

use Bio::SeqIO;
use Bio::Seq::SeqFactory;
use Bio::Species;
use Bio::Annotation::DBLink;
#use Bio::Annotation::Reference;
use Bio::Annotation::Comment;
use Bio::Annotation::SimpleValue;
use Bio::Annotation::OntologyTerm;

@ISA = qw(Bio::SeqIO);

# list of all the field names in locuslink
my @locuslink_keys = qw(
		       ACCNUM
		       ALIAS_PROT
		       ALIAS_SYMBOL
		       ASSEMBLY
		       BUTTON
		       CDD
		       CHR
		       COMP
		       CONTIG
		       CURRENT_LOCUSID
		       DB_DESCR
		       DB_LINK
		       ECNUM
		       EVID
		       EXTANNOT
		       GO
		       GRIF
		       LINK
		       LOCUSID
		       LOCUS_CONFIRMED
		       LOCUS_TYPE
		       MAP
		       MAPLINK
		       NC
		       NG
		       NM
		       NP
		       NR
		       OFFICIAL_GENE_NAME
		       OFFICIAL_SYMBOL
		       OMIM
		       ORGANISM
		       PHENOTYPE
		       PHENOTYPE_ID
		       PMID
		       PREFERRED_GENE_NAME
		       PREFERRED_PRODUCT
		       PREFERRED_SYMBOL
		       PRODUCT
		       PROT
		       RELL
		       STATUS
		       STS
		       SUMFUNC
		       SUMMARY
		       TRANSVAR
		       TYPE
		       UNIGENE
		       XG
		       XM
		       XP
		       XR
		       );

# list of fields to make simple annotations from
# fields not listed here or as a key in feature hash are ignored (lost).
my %anntype_map = (
		   SimpleValue => [qw(
				      ALIAS_PROT
				      ALIAS_SYMBOL
				      CDD
				      CHR
				      CURRENT_LOCUSID
				      ECNUM
				      EXTANNOT
				      MAP
				      NC
				      NR
				      OFFICIAL_GENE_NAME
				      OFFICIAL_SYMBOL
				      PHENOTYPE
				      PREFERRED_GENE_NAME
				      PREFERRED_PRODUCT
				      PREFERRED_SYMBOL
				      PRODUCT
				      RELL
				      SUMFUNC
				      )
				   ],
		   Comment     => [qw(
				      SUMMARY
				      )
				   ],
		   );


# certain fields are not named the same as the symgene database list
my %dbname_map = (
		  pfam    => 'Pfam',
		  smart   => 'SMART',
		  NM      => 'RefSeq',
		  NP      => 'RefSeq',
		  XP      => 'RefSeq',
		  XM      => 'RefSeq',
		  NG      => 'RefSeq',
		  XG      => 'RefSeq',
		  XR      => 'RefSeq',
		  PROT    => 'GenBank',
		  ACCNUM  => 'GenBank',
		  CONTIG  => 'GenBank',
		  # certain fields are not named the same as the symgene
		  # database list: rename the fields the symgene database name
		  # key = field name in locuslink
		  # value = database name in sym
		  #GO      => 'GO',
		  OMIM    => 'MIM',
		  GRIF    => 'GRIF',
		  STS     => 'STS',
		  UNIGENE => 'UniGene',
		  );

# certain CDD entries use the wrong prefix for the accession number
# cddprefix will replace the key w/ the value for these entries
my %cddprefix = (
		 pfam     => 'PF',
		 smart    => 'SM',
		 );

# alternate mappings if one field does not exist
my %alternate_map = (
		  OFFICIAL_GENE_NAME => 'PREFERRED_GENE_NAME',
		  OFFICIAL_SYMBOL    => 'PREFERRED_SYMBOL',
		    );

# for these field names, we only care about the first value X in value X|Y|Z
my @ll_firstelements = qw(
			  NM
			  NP
			  NG
			  XG
			  XM
			  XP
			  XR
			  PROT
			  STS
			  ACCNUM
			  CONTIG
			  GRIF
			  );

# set the default search pattern for all the field names
my %feature_pat_map = map { ($_ , "^$_: (.+)\n"); } @locuslink_keys;

sub _initialize {
  my($self,@args) = @_;

  $self->SUPER::_initialize(@args);

  # overwrite the search pattern w/ the first value pattern
  foreach my $key(@ll_firstelements){
      $feature_pat_map{$key}="^$key: ([^|]+)";
  }

  # special search pattern for cdd entries
  foreach my $key(keys %cddprefix) {
      $feature_pat_map{$key}='^CDD: .+\|'.$key.'(\d+)';
  }

  # special patterns for specific fieldsq
  $feature_pat_map{MAP}      = '^MAP: (.+?)\|';
  $feature_pat_map{MAPHTML}  = '^MAP: .+\|(<.+>)\|';
  $feature_pat_map{GO}       = '^GO: .+\|.+\|\w+\|(GO:\d+)\|';
  $feature_pat_map{GO_DESC}  = '^GO: .+\|(.+)\|\w+\|GO:\d+\|';
  $feature_pat_map{GO_CAT}   = '^GO: (.+)\|.+\|\w+\|GO:\d+\|';
  $feature_pat_map{EXTANNOT} = '^EXTANNOT: (.+)\|(.+)\|\w+\|.+\|\d+';

  # set the sequence factory of none has been set already
  if(! $self->sequence_factory()) {
      $self->sequence_factory(Bio::Seq::SeqFactory->new(
					      -type => 'Bio::Seq::RichSeq'));
  }
}


#########################
#
sub search_pattern{
#
#########################
        my $entry=$_[0];		#text to search
        my $searchconfirm =$_[1];	#to make sure you got the right thing
	my $searchpattern=$_[2];
	my $searchtype=$_[3];
        my @query;
        @query=$entry=~/$searchpattern/gm;
        if (!(@query) && ($searchconfirm ne "FALSE")){
                print "No $searchtype found\n$entry\n";
                return;
        }elsif(!(@query)){
		return;
	}
        if ($searchconfirm ne "FALSE"){
		foreach (@query){
		  if (!($_=~/$searchconfirm/)){
		    print "error\n$entry\n$searchtype parse $_ does not match $searchconfirm\n";
                    exit;
		  }
		}#endforeach
        }#endsearchconfirm
        return(@query);
}#endsub
############
#
sub read_species{
#
############
	my ($spline)=@_;
	my $species;
	my $genus;
	($genus,$species)=$spline=~/([^ ]+) ([^ ]+)/;
	my $make = Bio::Species->new();
	$make->classification( ($species,$genus) );
	return $make;
}
################
#
sub read_dblink{
#
################
	my ($ann,$db,$ref)=@_;
	my @results=$ref ? @$ref : ();
	foreach my $id(@results){
	  if($id){
	    $ann->add_Annotation('dblink',
				 Bio::Annotation::DBLink->new(
							  -database =>$db ,
							  -primary_id =>$id));
	  }
	}
	return($ann);
}

################
#
sub read_reference{
#
################
        my ($ann,$db,$results)=@_;

	if($results){	
	    chomp($results);
	    my @ids=split(/,/,$results);
	    $ann = read_dblink($ann,$db,\@ids) if @ids;
	}
	return $ann; 
}#endsub

################
#
sub add_annotation{
#
################
    my ($ac,$type,$text,$anntype)=@_;
    my @args;

    $anntype = 'SimpleValue' unless $anntype;
    SWITCH : {
	$anntype eq 'SimpleValue' && do {
	    push(@args, -value => $text, -tagname => $type);
	    last SWITCH;
	};
	$anntype eq 'Comment'     && do {
	    push(@args, -text  => $text, -tagname => 'comment');
	    last SWITCH;
	};
    }
    $ac->add_Annotation("Bio::Annotation::$anntype"->new(@args));
    return($ac);
}#endsub

################
#
sub add_annotation_ref{
#
################
        my ($ann,$type,$textref)=@_;
	my @text=$textref ? @$textref : ();
	
	foreach my $text(@text){
		$ann->add_Annotation($type,Bio::Annotation::SimpleValue->new(-value => $text));
        }
        return($ann);
}#endsub

################
#
sub make_unique{
#
##############
    my ($ann,$key) = @_;
    
    my %seen = ();
    foreach my $dbl ($ann->remove_Annotations($key)) {
	if(! $seen{$dbl->as_text()}) {
	    $seen{$dbl->as_text()} = 1;
	    $ann->add_Annotation($dbl);
	} 
    }
    return $ann;
}

################
#
sub next_seq{
#
##############
	my ($self, @args)=@_;
	my (%record,@results,$search,$ref,$cddref);
	my ($PRESENT,@keep);
	# my $seq=Bio::Seq->new();
	my $seq = $self->sequence_factory->create(-verbose =>$self->verbose());
	my $ann = $seq->annotation();


	# LOCUSLINK entries begin w/ >>
	local $/=">>";

	# slurp in a whole entry and return if no more entries
	return unless my $entry = $self->_readline;

	# if its the first entry you have to slurp it in again
	if ($entry eq '>>'){ #first entry
	  return unless $entry = $self->_readline;
	}

	if (!($entry=~/LOCUSID/)){
	    $self->throw("No LOCUSID in first line of record. ".
			 "Not LocusLink in my book.");
	}


	# loop through list of features and get field values
	# place into record hash as array refs
	foreach my $key (keys %feature_pat_map){
		$search=$feature_pat_map{$key};
		@results=search_pattern($entry,'FALSE',$search,$search);
		if ($results[0]){
			push @{$record{$key}},@results;
		}else{
			$record{$key}=0;
		}
	}#endfor

#for the purposes of symgene, we don't care about old loci
	if($record{CURRENT_LOCUSID}){
	  #print "current locus\n";
	  return($self->next_seq(@args));
	}


	# special processing for CDD entries like pfam and smart
	foreach my $key(keys %cddprefix){
	  #print "check CDD $key\n";
	  @keep=();
	  $ref=$record{$key};
	  foreach my $list($ref ? @$ref : ()){
	    # replace AC with correct AC number
	    push(@keep,$cddprefix{$key}.$list);	    
	  }
	  # replace CDD ref with correctly prefixed AC number
	  delete $record{$key};
	  push @{$record{$key}},@keep;
       	}
	# modify CDD references	@=();
	$cddref=$record{CDD};
	@keep=();
	foreach ($cddref ? @$cddref : ()){
	  $PRESENT='FALSE';
	  foreach my $key(keys %cddprefix){
	    if ($_=~/$key/){
	      $PRESENT = 'TRUE';
	      last;
	    }
	  }
	  if($PRESENT eq 'FALSE'){
	    push(@keep,$_);
	  }
	}
	delete $record{CDD};
	push @{$record{CDD}},@keep;

	foreach my $field(keys %dbname_map){
	    $ann=read_dblink($ann,$dbname_map{$field},$record{$field});
	}
	
	# add GO link as an OntologyTerm annotation
	if($record{GO}) {
	    for(my $j = 0; $j < @{$record{GO}}; $j++) {
		my $goann = Bio::Annotation::OntologyTerm->new(
					   -identifier => $record{GO}->[$j],
					   -name => $record{GO_DESC}->[$j],
					   -category => $record{GO_CAT}->[$j]);
		$ann->add_Annotation($goann);
	    }
	}

	$ann=add_annotation_ref($ann,'URL',$record{LINK});
	$ann=add_annotation_ref($ann,'URL',$record{DB_LINK});
	$ann=read_reference($ann,'PUBMED',first_element($record{PMID}));
	my @assembly=split(/,/,first_element($record{ASSEMBLY}));
        $ann=read_dblink($ann,'GenBank',\@assembly);

	# presently we can't store types of dblinks - hence make unique
	make_unique($ann,'dblink');

	# everything else gets a simple tag or comment value annotation
	foreach my $anntype (keys %anntype_map) {
	    foreach my $key (@{$anntype_map{$anntype}}){
		if($record{$key}){
		    foreach (@{$record{$key}}){
			#print "$key\t\t$_\n";
			$ann=add_annotation($ann,$key,$_,$anntype);
		    }
		}
	    }
	}

	$seq->accession_number(first_element($record{LOCUSID}));
	# replace fields w/ alternate if original does not exist
	foreach my $fieldval (keys %alternate_map){
	  if($record{$fieldval}){
	    $record{$fieldval}=first_element($record{$fieldval});
	  }elsif($record{$alternate_map{$fieldval}}){
	    $record{$fieldval}=first_element($record{$alternate_map{$fieldval}});
	  }else{
	    #print "Can't find $fieldval\n$entry\n";
	    $record{$fieldval}='None';
	  }
	}
        $seq->desc($record{OFFICIAL_GENE_NAME});
	$seq->display_id($record{OFFICIAL_SYMBOL});
	$seq->species(read_species(first_element($record{ORGANISM})));

	# dump out object contents
	# show_obj([$seq]);

	# set the sequence annotation
	$seq->annotation($ann);
	return($seq);
}
sub first_element{
	my $ref=$_[0];

	return($ref ? @$ref[0] : 0);
}#endsub
################
#
sub show_obj{
#
################
        my ($seqlistref)=@_;
        my @list=@$seqlistref;
        my $out = Bio::SeqIO->new('-fh' => \*STDOUT, -format => 'genbank' );
	my ($ann,@values,$val);

        foreach my $seq(@list){
                $out->write_seq($seq);
                $ann=$seq->annotation;
       		foreach my $key ( $ann->get_all_annotation_keys() ) {
        		@values = $ann->get_Annotations($key);
	           	foreach my $value ( @values ) {
	              		# value is an Bio::AnnotationI, and defines a "as_text" method
				$val=$value->as_text;
	             		print "Annotation ",$key,"\t\t",$val,"\n";
	           	}
		}
        }
}#endsub

1;
