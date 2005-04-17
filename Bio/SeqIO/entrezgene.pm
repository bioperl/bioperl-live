# $Id$
# BioPerl module for Bio::SeqIO::entrezgene
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::entrezgene - Entrez Gene ASN1 parser

=head1 SYNOPSIS

   # don't instantiate directly - instead do
   my $seqio = Bio::SeqIO->new(-format => 'entrezgene',
                               -file => $file);
   my $gene = $seqio->next_seq;

=head1 DESCRIPTION

This is EntrezGene ASN bioperl parser. It is built on top of 
GI::Parser::Entrezgene, a low level ASN parser built by Mingyi Liu 
(sourceforge.net/projetcs/egparser). The easiest way to use it is 
shown above.

You will get most of the EntrezGene annotation such as gene symbol, 
gene name and description, accession numbers associated 
with the gene, etc. Almost all of these are given as Annotation objects.
A comprehensive list of those objects will be available here later.

If you need all the data do:

   my $seqio = Bio::SeqIO->new(-format => 'entrezgene',
                               -file => $file,
                               -debug => 'on');
   my ($gene,$genestructure,$uncaptured) = $seqio->next_seq;

The $genestructure is a Bio::Cluster::SequenceFamily object. It 
contains all refseqs and the genomic contigs that are associated with 
the paricular gene. You can also modify the output $seq to allow back 
compatibility with old LocusLink parser:

   my $seqio = Bio::SeqIO->new(-format => 'entrezgene',
                               -file => $file,
                               -locuslink => 'convert');

The -debug and -locuslink options slow down the parser.

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
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Stefan Kirov

Email skirov at utk.edu

Describe contact details here

=head1 CONTRIBUTORS

Hilmar Lapp, hlapp at gmx.net

=head1 APPENDIX

This parser is based on GI::Parser::EntrezGene module

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::SeqIO::entrezgene;

use strict;
use vars qw(@ISA);
use Bio::ASN1::EntrezGene;
use Bio::Seq;
use Bio::Species;
use Bio::Annotation::SimpleValue;
use Bio::Annotation::DBLink;
use Bio::Annotation::Comment;
use Bio::SeqFeature::Generic;
use Bio::Annotation::Reference;
use Bio::SeqFeature::Gene::Exon;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::GeneStructure;
use Bio::Cluster::SequenceFamily;
#use Bio::Ontology::Ontology; Relationships.... later
use Bio::Ontology::Term;
use Bio::Annotation::OntologyTerm;
use vars qw(@ISA);
@ISA = qw(Bio::SeqIO);

%main::eg_to_ll =('Official Full Name'=>'OFFICIAL_GENE_NAME',
						'chromosome'=>'CHR',
						'cyto'=>'MAP', 
						'Official Symbol'=>'OFFICIAL_SYMBOL');
@main::egonly = keys %main::eg_to_ll;
# We define $xval and some other variables so we don't have 
# to pass them as arguments
my ($seq,$ann,$xval,%seqcollection,$buf);

sub _initialize {
	my($self,@args) = @_;
	$self->SUPER::_initialize(@args);
	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	$self->{_debug}=$param{-debug};
	$self->{_locuslink}=$param{-locuslink};
	$self->{_parser}=Bio::ASN1::EntrezGene->new(file=>$param{-file});
	#Instantiate the low level parser here (it is -file in Bioperl
   #-should tell M.)
	#$self->{_parser}->next_seq; #First empty record- bug in Bio::ASN::Parser
}


sub next_seq {
    my $self=shift;
    my $value = $self->{_parser}->next_seq(-trimopt=>1); 
	 # $value contains data structure for the
	 # record being parsed. 2 indicates the recommended
	 # trimming mode of the data structure
	 #I use 1 as I prefer not to descend into size 0 arrays
	 return undef unless ($value);
    my $debug=$self->{_debug};
    $ann = Bio::Annotation::Collection->new();
    my @alluncaptured;
    # parse the entry
    my @keys=keys %{$value};
    $xval=$value->[0];
    #Basic data
	 #$xval->{summary}=~s/\n//g; 
    $seq = Bio::Seq->new(
                        -display_id  => $xval->{gene}{locus},
                        -accession_number =>$xval->{'track-info'}{geneid},
                        -desc=>$xval->{summary}
                   );
#DEBUG:REMOVE
#if ($seq->accession_number == 1107) {
#	print "debug record\n";
#}
#DEBUG:REMOVE
    #Source data here
    _add_to_ann($xval->{status},'Entrez Gene Status'); 
    my $lineage=$xval->{source}{org}{orgname}{lineage};
    my $sp=$xval->{source}{org}{orgname}{name}{binomial}{species};
    $lineage=~s/[\s\n]//g;
    my ($comp,@lineage);
    while ($lineage) {
        ($comp,$lineage)=split(/;/,$lineage,2);
        unshift @lineage,$comp;
    }
    unshift @lineage,$sp;
     my $specie=new Bio::Species(-classification=>[@lineage],
                                -ncbi_taxid=>$xval->{source}{org}{db}{tag}{id});
    $specie->common_name($xval->{source}{org}{common});
    if (exists($xval->{source}->{subtype}) && ($xval->{source}->{subtype})) {
        if (ref($xval->{source}->{subtype}) eq 'ARRAY') {
            foreach my $subtype (@{$xval->{source}->{subtype}}) {
                _add_to_ann($subtype->{name},$subtype->{subtype});
            }
        }
        else {
            _add_to_ann($xval->{source}->{subtype}->{name},$xval->{source}->{subtype}->{subtype}); 
        }
    }
    #Synonyms
    if (ref($xval->{gene}->{syn}) eq 'ARRAY') {
        foreach my $symsyn (@{$xval->{gene}->{syn}}) {
        _add_to_ann($symsyn,'ALIAS_SYMBOL');
        }
    }
    else {
        _add_to_ann($xval->{gene}->{syn},'ALIAS_SYMBOL');
    }
    
    
    #COMMENTS (STS not dealt with yet)
    if (ref($xval->{comments}) eq 'ARRAY') {
        for my $i (0..$#{$xval->{comments}}) {
            push @alluncaptured,_process_all_comments($xval->{comments}->[$i]);
           }
    }
    else {
            push @alluncaptured,_process_all_comments($xval->{comments});
    }
       #Gene
       if (exists($xval->{gene}->{db})) {
       if (ref($xval->{gene}->{db}) eq 'ARRAY') {
        foreach my $genedb (@{$xval->{gene}->{db}}) {
            _add_to_ann($genedb->{tag}->{id},$genedb->{db});
        }
        }
        else {
            _add_to_ann($xval->{gene}->{db}->{tag}->{id},$xval->{gene}->{db}->{db});
        }
        delete $xval->{gene}->{db} unless ($debug eq 'off');
        }
       #LOCATION To do: uncaptured stuff
       if (exists($xval->{location})) {
        if (ref($xval->{location}) eq 'ARRAY') {
            foreach my $loc (@{$xval->{location}}) {
                _add_to_ann($loc->{'display-str'},$loc->{method}->{'map-type'});
            }
        }
        else {
            _add_to_ann($xval->{location}->{'display-str'},$xval->{location}->{method}->{'map-type'});
        }
        delete $xval->{location} unless ($debug eq 'off');
       }
       #LOCUS
       if (ref($xval->{locus}) eq 'ARRAY') {
       foreach my $locus (@{$xval->{locus}}) {
        push @alluncaptured,_process_locus($locus);
        }
       }
        else {
            push @alluncaptured,_process_locus($xval->{locus});
        }
        #Homology
        my ($uncapt,$hom,$anchor)=_process_src($xval->{homology}->{source});
        foreach my $homann (@$hom) {
            $ann->add_Annotation('dblink',$homann);
        }
        push @alluncaptured,$uncapt;
        #Index terms
        if (exists($xval->{'xtra-index-terms'})) {
        if (ref($xval->{'xtra-index-terms'}) eq 'ARRAY') {
          foreach my $term (@{$xval->{'xtra-index-terms'}}) {
           _add_to_ann($term,'Index terms');
           }
        }
        else {
          _add_to_ann($xval->{'xtra-index-terms'},'Index terms');
        }
        }
        #PROPERTIES
        my @prop;
        if (exists($xval->{properties})) {
        if (ref($xval->{properties}) eq 'ARRAY') {
          foreach my $property (@{$xval->{properties}}) {
            push @alluncaptured,_process_prop($property);
           }
        }
        else {
          push @alluncaptured,_process_prop($xval->{properties});
        }
        }
        $seq->annotation($ann) unless ($self->{_locuslink} eq 'convert');
        $seq->species($specie);
        my @seqs;
        foreach my $key (keys %seqcollection) {
          push @seqs,@{$seqcollection{$key}};
        }
        my $cluster = Bio::Cluster::SequenceFamily->new(-family_id=>$seq->accession_number,
                                                 -description=>"Entrez Gene " . $seq->accession_number,
                                               -members=>\@seqs);#Our EntrezGene object
        #clean
    unless ($debug eq 'off') {
        delete $xval->{homology}->{source};
        delete($xval->{summary});
        delete($xval->{'track-info'});
        delete($xval->{gene}{locus});
        delete($xval->{source}{org}{orgname}{lineage});
        delete $xval->{source}{org}{orgname}{name}{binomial}{species};
        delete $xval->{gene}{syn};
        delete $xval->{source}->{subtype};
        delete $xval->{comments};
        delete $xval->{properties};
        delete $xval->{'xtra-index-terms'};
        $xval->{status};
    }
    push @alluncaptured,$xval;
        undef %seqcollection;
    undef $xval;
    #print 'x';
    &_backcomp_ll if ($self->{_locuslink} eq 'convert');
    $ann->DESTROY;
    return wantarray ? ($seq,$cluster,\@alluncaptured):$seq;#Hilmar's suggestion
  }

sub _process_refseq {
my $products=shift;
my (@uncaptured,@products);
my $nann=Bio::Annotation::Collection->new();
if (ref($products) eq 'ARRAY') { @products=@{$products}; }
else {push @products,$products ;}
foreach my $product (@products) {
    if (($product->{seqs}->{whole}->{gi})||($product->{accession})){#Minimal data required
   my $nseq = Bio::Seq->new(
                        -accession_number => $product->{seqs}->{whole}->{gi},
                        -display_id=>$product->{accession},
                        -desc=> $product->{heading}
                   );
                   if ($product->{source}) {
                    my ($uncapt,$allann)=_process_src($product->{source});
                    delete $product->{source};
                    push @uncaptured,$uncapt;
                    foreach my $annotation (@{$allann}) {
                        $nann->add_Annotation('dblink',$annotation);
                    }
                    }
    delete  $product->{seqs}->{whole}->{gi};
    delete $product->{accession};
    delete $product->{source};
    delete $product->{heading};
    my ($allan,$allfeat,$uncapt)=_process_comments($product->{comment});
    push @uncaptured,$uncapt;
    foreach my $feat (@$allfeat) {
        $nseq->add_SeqFeature($feat);
    }
    foreach my $annotation (@{$allan}) {
        $nann->add_Annotation('dblink',$annotation);
    }
    $nseq->annotation($nann);
   push @{$seqcollection{seq}},$nseq;
    }
    if ($product->{products}) {
       my @uncapt=_process_refseq($product->{products});
       push @uncaptured,@uncapt;
    }
}
return @uncaptured;
}

sub _process_other_seqs {
}
sub _process_links {
 my $links=shift;
 my (@annot,@uncapt);
 if (ref($links) eq 'ARRAY') {
    foreach my $link (@$links) {
        my ($uncapt,$annot)=_process_src($link->{source});
        push @uncapt,$uncapt;
        foreach my $annotation (@$annot) {
          $ann->add_Annotation('dblink',$annotation);
        }
    }
 }
 else { my ($uncapt,$annot)=_process_src($links->{source});         
        push @uncapt,$uncapt;
        foreach my $annotation (@$annot) {
          $ann->add_Annotation('dblink',$annotation);
        }
    }
return @uncapt;
}

sub _add_to_ann {#Highest level only
my ($val,$tag)=@_;
  #  $val=~s/\n//g;#Low level EG parser leaves this so we take care of them here
    unless ($tag) {
     warn "No tagname for value $val, tag $tag ",$seq->id,"\n";
     return;
    }
        my $simann=new Bio::Annotation::SimpleValue(-value=>$val,-tagname=>$tag);
        $ann->add_Annotation($simann);
}

sub _process_comments {
 my $comment=shift;
 return undef unless ($comment);#Should be more careful when calling _process_comment:To do
    my (@uncaptured,@comments);
    my (@nfeat,@anncol);
    if (ref($comment) eq 'ARRAY') { @comments=@{$comment}; }
    else {push @comments,$comment ;}
    for my $i (0..$#comments) {#Each comments is a
        my ($desc,$nfeat,$add,@ann,@comm);
        my $comm=$comments[$i];
        my $nann=Bio::Annotation::Collection->new();
        my $heading=$comm->{heading} if ($comm->{heading});
        unless (exists($comm->{comment})) {undef $comm->{comment}; $add=1;}#Trick in case we miss something
        while ((exists($comm->{comment})&&$comm->{comment})) {
            if ($comm->{source}) {
               my ($uncapt,$allann,$anchor) = _process_src($comm->{source});
           if ($allann) {
               @ann=@{$allann};
            delete $comm->{source};
            push @uncaptured,$uncapt;
                    foreach my $annotation (@{$allann}) {
                        $nann->add_Annotation('dblink',$annotation);
                         if ($annotation->{_anchor}) {$desc.=$annotation->{_anchor}.' ';}
                         $annotation->optional_id($heading);
                    }
        }
            }
            $comm=$comm->{comment};
            if (ref($comm) eq 'ARRAY') {
              @comm=@{$comm};
            }
            else {
                push @comm,$comm;
            }
            foreach my $ccomm (@comm) {
            next unless ($ccomm);
            my @sfann;
            if (exists($ccomm->{source})) {
                my ($uncapt,$allann,$anchor) = _process_src($ccomm->{source});
               if ($allann) {
                   @sfann=@{$allann};
                delete $ccomm->{source};
                push @uncaptured,$uncapt;
            }
            }
            if ((exists($ccomm->{text}))&&($ccomm->{text}=~/Location/i)){
            my ($l1,$rest)=split(/-/,$ccomm->{text});
            $l1=~s/\D//g;
            $rest=~s/^\s//;
            my ($l2,@rest)=split(/\s/,$rest);
            my (%tags,$tag);
            unless ($l1) {
                next;
            }
            $nfeat=Bio::SeqFeature::Generic->new(-start=>$l1, -end=>$l2, -strand=>$tags{strand}, -source=>$ccomm->{type},
                                -seq_id=>$desc, primary=>$heading);
            foreach my $rest (@rest) {
                $rest=~s/[:,]$//;
                if ($rest=~/\D/) {$tag.=lc($rest).' '; next;}
                next unless (($tag)&&($rest));
                 $nfeat->add_tag_value($tag,$rest);
                 undef $tag;
            }
            foreach my $sfann (@sfann) {
                $nann->add_Annotation('dblink',$sfann);
            }
            $nfeat->annotation($nann);
            push @nfeat,$nfeat;
            delete $ccomm->{text};
            delete $ccomm->{type};
            }
        elsif (exists($ccomm->{label})) {
            _add_to_ann($ccomm->{text},$ccomm->{label});
            delete $ccomm->{text};
            delete $ccomm->{label};
            push @uncaptured,$ccomm;
        }
        elsif (exists($ccomm->{text})) {
            _add_to_ann($ccomm->{text},'description');
            delete $ccomm->{text};
            push @uncaptured,$ccomm;
        }
        }#Bit clumsy but that's what we get from the low level parser
            unless ($nfeat){ push @anncol,@ann;}
        push @uncaptured,$comm;
    }
    }
    return \@anncol,\@nfeat,\@uncaptured;
}

sub _process_src {
    my $src=shift;
    return undef unless (exists($src->{src}->{tag}));
    my @ann;
    my $db=$src->{src}->{db};
    delete $src->{src}->{db};
    my $anchor=$src->{anchor};
    delete $src->{anchor};
    my $url;
    if ($src->{url}) {
            $url=$src->{url};
            $url=~s/\n//g;
            delete $src->{url};
        }
        if ($src->{src}->{tag}->{str}) {
            my @sq=split(/[,;]/,$src->{src}->{tag}->{str});
            delete $src->{src}->{tag};
            foreach my $id (@sq) {
                $id=~s/\n//g;
                undef $anchor if ($anchor eq 'id');
                my $simann=new Bio::Annotation::DBLink(-database => $db,
                                        -primary_id => $id,
                    );
                $simann->url($url) if ($url);#DBLink should have URL!
                push @ann, $simann;
            }
        }
        else {
            my $id=$src->{src}->{tag}->{id};
            delete $src->{src}->{tag};
            undef $anchor if ($anchor eq 'id');
            $id=~s/\n//g;
            my $simann=new Bio::Annotation::DBLink(-database => $db,
                                        -primary_id => $id
                    );
            $simann->{_anchor}=$anchor if ($anchor);
            $simann->{_url}=$url if ($url);#DBLink should have URL!
            push @ann, $simann;
        }
        return $src, \@ann,$anchor;
}

sub _add_references {
my $refs=shift;
if (ref($refs) eq 'ARRAY') {
    foreach my $ref(@$refs) {
        my $refan=new Bio::Annotation::Reference(-database => 'Pubmed',
                                        -primary_id => $ref);
        $ann->add_Annotation('Reference',$refan);
    }
}
else {
    my $refan=new Bio::Annotation::Reference(-database => 'Pubmed',
                                        -primary_id => $refs);
        $ann->add_Annotation('Reference',$refan);
}
}

#Should we do this at all if no seq coord are present?
sub _process_locus {
my $locus=shift;
my @uncapt;
my $gseq=new Bio::Seq(-display_id=>$locus->{accession},-version=>$locus->{version},
            -accession_number=>$locus->{seqs}->{'int'}->{id}->{gi});
delete $locus->{accession};
delete $locus->{version};
delete $locus->{'int'}->{id}->{gi};
my ($start,$end,$strand);
if (exists($locus->{seqs}->{'int'}->{from})) {
 $start=$locus->{seqs}->{'int'}->{from};
 delete $locus->{seqs}->{'int'}->{from};
 #unless ($start) {print $locus->{seqs}->{'int'}->{from},"\n",$locus,"\n";}
 $end=$locus->{seqs}->{'int'}->{to};
 delete $locus->{seqs}->{'int'}->{to};
 delete $locus->{seqs}->{'int'}->{strand};
 $strand=$locus->{seqs}->{'int'}->{strand} eq 'minus'?-1:1;
    my $nfeat=Bio::SeqFeature::Generic->new(-start=>$start, -end=>$end, -strand=>$strand);
    $gseq->add_SeqFeature($nfeat);
}
my @products;
if (ref($locus->{products}) eq 'ARRAY') {
    @products=@{$locus->{products}};
}
else {
    push @products,$locus->{products};
}
delete $locus->{products};
my $gstruct=new Bio::SeqFeature::Gene::GeneStructure;
foreach my $product (@products) {
    my ($tr,$uncapt)=_process_products_coordinates($product,$start,$end,$strand);
    $gstruct->add_transcript($tr);
    undef $tr->{parent}; #Because of a cycleG
    push @uncapt,$uncapt;
}
$gseq->add_SeqFeature($gstruct);
push @{$seqcollection{genestructure}},$gseq;
return @uncapt;
}

=head1 _process_products_coordinates
To do:
=cut


sub _process_products_coordinates {
my $coord=shift;
my $start=shift||0;#In case it is not known: should there be an entry at all?
my $end=shift||1;
my $strand=shift||1;
my (@coords,@uncapt);
my $transcript=new Bio::SeqFeature::Gene::Transcript(-primary=>$coord->{accession},
                                          -start=>$start,-end=>$end,-strand=>$strand, -desc=>$coord->{type});

if ((exists($coord->{'genomic-coords'}->{mix}->{'int'}))||(exists($coord->{'genomic-coords'}->{'packed-int'}))) {
@coords=exists($coord->{'genomic-coords'}->{mix}->{'int'})?@{$coord->{'genomic-coords'}->{mix}->{'int'}}:
                                    @{$coord->{'genomic-coords'}->{'packed-int'}};
foreach my $exon (@coords) {
    next unless (exists($exon->{from}));
    my $exonobj=new Bio::SeqFeature::Gene::Exon(-start=>$exon->{from},-end=>$exon->{to},-strand=>$strand);
    $transcript->add_exon($exonobj);
    delete $exon->{from};
    delete $exon->{to};
    delete $exon->{strand};
    push @uncapt,$exon;
}
}
my ($prot,$uncapt);
if (exists($coord->{products})) {
    my ($prot,$uncapt)=_process_products_coordinates($coord->{products},$start,$end,$strand);
    $transcript->add_SeqFeature($prot);
    push @uncapt,$uncapt;
}
return $transcript,\@uncapt;
}

=head1 _process_prop
To do: process GO
=cut
sub _process_prop {
    my $prop=shift;
    my @uncapt;
    if (exists($prop->{properties})) {#Iterate
        if (ref($prop->{properties}) eq 'ARRAY') {
            foreach my $propn (@{$prop->{properties}}) {
               push @uncapt,_process_prop($propn);
            }
        }
        else {
            push @uncapt,_process_prop($prop->{properties});
        }
    }
    unless ((exists($prop->{heading})) && ($prop->{heading} eq 'GeneOntology')) {
        _add_to_ann($prop->{text},$prop->{label}) if (exists($prop->{text})); 
        delete $prop->{text};
        delete $prop->{label};
        push @uncapt,$prop;
        return \@uncapt;
    }
    #Will do GO later
    if (exists($prop->{comment})) {
    push @uncapt,_process_go($prop->{comment});
    }
}


sub _process_all_comments {
my $product=shift;
my @alluncaptured;
my $heading=$product->{heading} if (exists($product->{heading}));
           if ($heading) {
               delete $product->{heading};
               CLASS: {
                   if ($heading =~ 'RefSeq Status') {#IN case NCBI changes slightly the spacing:-)
                    _add_to_ann($product->{label},'RefSeq status');  last CLASS;
                   }
                   if ($heading =~ 'NCBI Reference Sequences') {#IN case NCBI changes slightly the spacing:-)
                    my @uncaptured=_process_refseq($product);
            push @alluncaptured,@uncaptured; last CLASS;
                   }
                   if ($heading =~ 'Related Sequences') {#IN case NCBI changes slightly the spacing:-)
                    my @uncaptured=_process_refseq($product);
                    push @alluncaptured,@uncaptured;  last CLASS;
                   }
                    if ($heading =~ 'Sequence Tagges Sites') {#IN case NCBI changes slightly the spacing:-)
                    my @uncaptured=_process_links($product);
                     push @alluncaptured,@uncaptured;
                     last CLASS;
                   }
                   if ($heading =~ 'Additional Links') {#IN case NCBI changes slightly the spacing:-)
                    push @alluncaptured,_process_links($product->{comment});
                     last CLASS;
                   }
                   if ($heading =~ 'LocusTagLink') {#IN case NCBI changes slightly the spacing:-)
                     _add_to_ann($product->{source}->{src}->{tag}->{id},$product->{source}->{src}->{db}); 
                    last CLASS;
                   }
                   if ($heading =~ 'Sequence Tagged Sites') {#IN case NCBI changes slightly the spacing:-)
                     push @alluncaptured,_process_STS($product); 
                     delete $product->{comment};
                    last CLASS;
                   }
               }
    }
	if (exists($product->{type})&&($product->{type} eq 'generif')) {
		push @alluncaptured,_process_grif($product);
		return @alluncaptured;#Maybe still process the comments?
	}
	if (exists($product->{refs})) {
                _add_references($product->{refs}->{pmid});
                delete $product->{refs}->{pmid}; push @alluncaptured,$product;
            }
	if (exists($product->{comment})) {
                my ($allan,$allfeat,$uncapt)=_process_comments($product->{comment});
                foreach my $cann (@$allan) {$ann->add_Annotation('dblink',$cann);}
                delete $product->{refs}->{comment}; push @alluncaptured,$uncapt;
            }

return @alluncaptured;
}

sub _process_STS {
my $product=shift;
}

sub _process_go {
    my $comm=shift;
    my @comm;
    push @comm,( ref($comm) eq 'ARRAY')? @{$comm}:$comm;
    foreach my $comp (@comm) {
        if (ref($comp->{comment}) eq 'ARRAY') {
            foreach my $go (@{$comp->{comment}}) {
                my $term=_get_go_term($go);
                my $annterm = new Bio::Annotation::OntologyTerm (-tagname => 'Gene Ontology');
                $annterm->term($term);
                $ann->add_Annotation('OntologyTerm',$annterm);
            }
        }
        else {
            my $term=_get_go_term($comp->{comment});
            my $annterm = new Bio::Annotation::OntologyTerm (-tagname => 'Gene Ontology');
            $annterm->term($term);
            $ann->add_Annotation('OntologyTerm',$annterm);
        }
    }
}

sub _process_grif {
my $grif=shift;
if (ref($grif->{comment}) eq 'ARRAY') {#Insane isn't it?
	my @uncapt;
	foreach my $product (@{$grif->{comment}}) {
		next unless (exists($product->{text})); 
		my $uproduct=_process_grif($product);
	#	$ann->add_Annotation($type,$grifobj);
		push @uncapt,$uproduct;
	}
	return \@uncapt;
}
if (exists($grif->{comment}->{comment})) {
	$grif=$grif->{comment};
}
my $ref= (ref($grif->{refs}) eq 'ARRAY') ? shift @{$grif->{refs}}:$grif->{refs};
my $refergene='';
my ($obj,$type);
if ($ref->{pmid}) {
	$refergene=$grif->{comment}->{source}->{src}->{tag}->{id} if (exists($grif->{comment}->{source})); #unfortunatrely we cannot put yet everything in
	my $grifobj=new  Bio::Annotation::Comment(-text=>$grif->{text});
	$obj = new Bio::Annotation::DBLink(-database => 'generif',
                                        -primary_id => $ref->{pmid}, #The pubmed id (at least the first one) which is a base for the conclusion
                                        -version=>$grif->{version},
                                        -optional_id=>$refergene
                    ); 
	$obj->comment($grifobj);
    $type='dblink';
}
else {
	$obj=new  Bio::Annotation::SimpleValue($grif->{text},'generif');
    $type='generif';
}
delete $grif->{text};
delete $grif->{version};
delete $grif->{type};
delete $grif->{refs};
$ann->add_Annotation($type,$obj);
return $grif;
}
sub _get_go_term {
my $go=shift;
    my $refan=new Bio::Annotation::Reference(-database => 'Pubmed', #We expect one ref per GO
        -primary_id => $go->{refs}->{pmid});
    my $term = Bio::Ontology::Term->new( 
        -identifier  => $go->{source}->{src}->{tag}->{id},
        -name        => $go->{source}->{anchor},
        -definition  => $go->{source}->{anchor},
        -comment     => $go->{source}->{'post-text'},
        -references  => [$refan],
        -version     =>$go->{version});
return $term;
}

sub _backcomp_ll {
my $self=shift;
my $newann=Bio::Annotation::Collection->new();
        #$newann->{_annotation}->{ALIAS_SYMBOL}=$ann->{_annotation}->{ALIAS_SYMBOL};
       # $newann->{_annotation}->{CHR}=$ann->{_annotation}->{chromosome};
       # $newann->{_annotation}->{MAP}=$ann->{_annotation}->{cyto};
       foreach my $tagmap (keys %{$ann->{_typemap}->{_type}}) {
	next if (grep(/$tagmap/,@main::egonly));
        $newann->{_annotation}->{$tagmap}=$ann->{_annotation}->{$tagmap};
	}
        #$newann->{_annotation}->{Reference}=$ann->{_annotation}->{Reference};
        #$newann->{_annotation}->{generif}=$ann->{_annotation}->{generif};
        #$newann->{_annotation}->{comment}=$ann->{_annotation}->{comment};
       # $newann->{_annotation}->{OFFICIAL_GENE_NAME}=$ann->{_annotation}->{'Official Full Name'};
        $newann->{_typemap}->{_type}=$ann->{_typemap}->{_type};
        foreach my $ftype (keys %main::eg_to_ll) {
		my $newkey=$main::eg_to_ll{$ftype};
		$newann->{_annotation}->{$newkey}=$ann->{_annotation}->{$ftype};
		$newann->{_typemap}->{_type}->{$newkey}='Bio::Annotation::SimpleValue';
		delete $newann->{_typemap}->{_type}->{$ftype};
		$newann->{_annotation}->{$newkey}->[0]->{tagname}=$newkey;
        }
	foreach my $dblink (@{$newann->{_annotation}->{dblink}}) {
            next unless ($dblink->{_url});
            my $simann=new Bio::Annotation::SimpleValue(-value=>$dblink->{_url},-tagname=>'URL');
            $newann->add_Annotation($simann);
        }

#        my $simann=new Bio::Annotation::SimpleValue(-value=>$seq->desc,-tagname=>'comment');
#        $newann->add_Annotation($simann);
    $seq->annotation($newann);
return 1;
}

1;
