#!/usr/bin/perl
=head1 NAME

biblio-eutils-example.pl

=head1 SYNOPSIS

Script that uses Bio::Biblio, accessing 'eutils' at PubMed.

As of Bioperl version 1.4 there are 3 bibliographic repositories,
stipulated by the -access argument: soap, eutils, and biofetch.
The default is 'soap'. Not all of these repositories support all
the Biblio methods nor are the contents of these repositories
necessarily the same. Choose wisely!

=head2 PubMed Queries

The syntax of the queries is the same as at PubMed, see
http://www.ncbi.nlm.nih.gov/entrez/query/Pmc/pmchelp.html#SearchFieldDescriptionsandTags
for more information on how to construct queries.

=head2 Parsing Results

Bio::Biblio will give you XML when querying eutils so you have
choose a method to parse XML. A fairly simple approach uses
XML::Twig, shown here. This example shows how query by title and
how to retrieve the titles of the abstracts found.

=cut

use strict;
use Bio::Biblio;
use XML::Twig;

# one-liner to get the number of abstracts found
my $num = new Bio::Biblio(-access => "eutils")->find("Osborne","authors")->
  get_count;

my $biblio = Bio::Biblio->new(-access => "eutils");

my $result = $biblio->find("brain [TI] AND MDM2 [TI]");

my $pmids = $result->get_all_ids;

my $parser = XML::Twig->new(twig_roots => {"ArticleTitle" => \&print_title} );

for my $pmid (@$pmids) {
	my $xml = $biblio->get_by_id($pmid);
	eval {
		$parser->parse($xml);
	};
	if ($@) {
		warn "Problem parsing PubMed $pmid XML: $!\n";
	}
}

sub print_title {
	my ($twig, $elt) = @_;
	print $elt->text,"\n";
	$twig->purge;
}

=head1 PubMed XML Example

<?xml version="1.0"?>
<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, 1st November 2004//EN" "http://www.ncbi.nlm.nih.gov/entrez/query/DTD/pubmed_041101.dtd">
<PubmedArticleSet>
<PubmedArticle>
    <MedlineCitation Owner="NLM" Status="MEDLINE">
        <PMID>15815077</PMID>
        <DateCreated>
            <Year>2005</Year>
            <Month>04</Month>
            <Day>07</Day>
        </DateCreated>
        <DateCompleted>
            <Year>2005</Year>
            <Month>08</Month>
            <Day>29</Day>
        </DateCompleted>
        <Article PubModel="Print">
            <Journal>
                <ISSN>0231-5882</ISSN>
                <JournalIssue>
                    <Volume>23</Volume>
                    <Issue>4</Issue>
                    <PubDate>
                        <Year>2004</Year>
                        <Month>Dec</Month>
                    </PubDate>
                </JournalIssue>
            </Journal>
            <ArticleTitle>Rabbit liver microsomal system: study of interaction with two model N-nitrosamines and their metabolism.</ArticleTitle>
            <Pagination>
                <MedlinePgn>423-33</MedlinePgn>
            </Pagination>
            <Abstract>
                <AbstractText>Rabbit liver microsomes of control (non-treated) or animals induced either by ethanol (EtOH) or phenobarbital (PB) were incubated with N-nitrosodimethylamine (NDMA) or N-nitrosomethylaniline (NMA). Difference spectroscopy showed that NMA is bound to the substrate-binding site of cytochrome P-450 (CYP) isoforms as heme ligand in control and EtOH pre-treated microsomes. On the other hand, PB-induced microsomes exhibit with NMA substrate type of spectra. NDMA does not provide any type of binding spectra with used microsomal systems. Oxidative bio-activation of N-nitrosamines by the microsomal CYP isoforms was measured as formaldehyde formation. Analysis of reaction kinetics in control microsomes revealed, for both substrates, two values of Michaelis-Menten constant (K(m)) for, K(m) values of 0.03 and 0.13 mmol/l for NDMA, and 0.30 and 0.82 mmol/l for NMA. Induction of animals with EtOH resulted in a decrease in the K(m) value for both substrates. In contrast, PB treatment caused an elevation of K(m) value for NDMA. Based on these data, we conclude that EtOH-inducible microsomal CYP isoforms (mainly CYP2E1) are responsible for binding and N-demethylation metabolism of both studied N-nitrosamines in rabbit liver microsomal system. The role of the other CYP isoforms involved in the metabolism of mentioned N-nitrosamines is discussed.</AbstractText>
            </Abstract>
            <Affiliation>Department of Biochemistry, Faculty of Science, Charles University, Hlavova 2030, 128 43 Prague 2, Czech Republic. mis@natur.cuni.cz</Affiliation>
            <AuthorList CompleteYN="Y">
                <Author ValidYN="Y">
                    <LastName>Sulc</LastName>
                    <ForeName>B</ForeName>
                    <Initials>B</Initials>
                </Author>
                <Author ValidYN="Y">
                    <LastName>Kubícková</LastName>
                    <ForeName>B</ForeName>
                    <Initials>B</Initials>
                </Author>
                <Author ValidYN="Y">
                    <LastName>Máslová</LastName>
                    <ForeName>F</ForeName>
                    <Initials>B</Initials>
                </Author>
                <Author ValidYN="Y">
                    <LastName>Hodek</LastName>
                    <ForeName>C</ForeName>
                    <Initials>B</Initials>
                </Author>
            </AuthorList>
            <Language>eng</Language>
            <PublicationTypeList>
                <PublicationType>Journal Article</PublicationType>
            </PublicationTypeList>
        </Article>
        <MedlineJournalInfo>
            <Country>Slovakia</Country>
            <MedlineTA>Gen Physiol Biophys</MedlineTA>
            <NlmUniqueID>8400604</NlmUniqueID>
        </MedlineJournalInfo>
        <ChemicalList>
            <Chemical>
                <RegistryNumber>0</RegistryNumber>
                <NameOfSubstance>N-nitrosodimethylamine</NameOfSubstance>
            </Chemical>
            <Chemical>
                <RegistryNumber>0</RegistryNumber>
                <NameOfSubstance>Nitrosamines</NameOfSubstance>
            </Chemical>
            <Chemical>
                <RegistryNumber>50-06-6</RegistryNumber>
                <NameOfSubstance>Phenobarbital</NameOfSubstance>
            </Chemical>
            <Chemical>
                <RegistryNumber>614-00-6</RegistryNumber>
                <NameOfSubstance>N-methyl-N-nitrosoaniline</NameOfSubstance>
            </Chemical>
            <Chemical>
                <RegistryNumber>64-17-5</RegistryNumber>
                <NameOfSubstance>Ethanol</NameOfSubstance>
            </Chemical>
            <Chemical>
                <RegistryNumber>9035-51-2</RegistryNumber>
                <NameOfSubstance>Cytochrome P-450 Enzyme System</NameOfSubstance>
            </Chemical>
        </ChemicalList>
        <CitationSubset>I</CitationSubset>
        <MeshHeadingList>
            <MeshHeading>
                <DescriptorName MajorTopicYN="N">Animals</DescriptorName>
            </MeshHeading>
            <MeshHeading>
                <DescriptorName MajorTopicYN="N">Cytochrome P-450 Enzyme System</DescriptorName>
                <QualifierName MajorTopicYN="Y">metabolism</QualifierName>
            </MeshHeading>
            <MeshHeading>
                <DescriptorName MajorTopicYN="N">Ethanol</DescriptorName>
                <QualifierName MajorTopicYN="Y">administration &amp; dosage</QualifierName>
            </MeshHeading>
            <MeshHeading>
                <DescriptorName MajorTopicYN="N">Liver</DescriptorName>
                <QualifierName MajorTopicYN="N">drug effects</QualifierName>
                <QualifierName MajorTopicYN="Y">metabolism</QualifierName>
            </MeshHeading>
            <MeshHeading>
                <DescriptorName MajorTopicYN="N">Male</DescriptorName>
            </MeshHeading>
            <MeshHeading>
                <DescriptorName MajorTopicYN="N">Microsomes, Liver</DescriptorName>
                <QualifierName MajorTopicYN="N">drug effects</QualifierName>
                <QualifierName MajorTopicYN="Y">metabolism</QualifierName>
            </MeshHeading>
            <MeshHeading>
                <DescriptorName MajorTopicYN="N">Nitrosamines</DescriptorName>
                <QualifierName MajorTopicYN="Y">metabolism</QualifierName>
            </MeshHeading>
            <MeshHeading>
                <DescriptorName MajorTopicYN="N">Phenobarbital</DescriptorName>
                <QualifierName MajorTopicYN="Y">administration &amp; dosage</QualifierName>
            </MeshHeading>
            <MeshHeading>
                <DescriptorName MajorTopicYN="N">Rabbits</DescriptorName>
            </MeshHeading>
            <MeshHeading>
                <DescriptorName MajorTopicYN="N">Research Support, Non-U.S. Govt</DescriptorName>
            </MeshHeading>
        </MeshHeadingList>
    </MedlineCitation>
    <PubmedData>
        <History>
            <PubMedPubDate PubStatus="pubmed">
                <Year>2005</Year>
                <Month>4</Month>
                <Day>9</Day>
                <Hour>9</Hour>
                <Minute>0</Minute>
            </PubMedPubDate>
            <PubMedPubDate PubStatus="medline">
                <Year>2005</Year>
                <Month>8</Month>
                <Day>30</Day>
                <Hour>9</Hour>
                <Minute>0</Minute>
            </PubMedPubDate>
        </History>
        <PublicationStatus>ppublish</PublicationStatus>
        <ArticleIdList>
            <ArticleId IdType="pubmed">15815077</ArticleId>
        </ArticleIdList>
    </PubmedData>
  </PubmedArticle>
 </PubmedArticleSet>

=cut
