# BioPerl module for Bio::DB::NextProt
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Copyright Felipe da Veiga Leprevost
#
# You may distribute this module under the same terms as perl itself.

=head1 NAME

Bio::DB::NextProt - Object interface to NextProt REST API.

=head1 SYNOPSIS

	my $np = Bio::DB::NextProt->new();

	my @result_1 = $np->search_cv(-query => "kinase#");

	my @result_2 = $np->get_isoform_info(-id => "NX_O00142-2");

	my @result_3 = $np->get_protein_cv_info(-id => "PTM-0205", -format => "html");

=head1 DESCRIPTION

The module allows the dynamic retrieval of information from the NextProt Database
through its API service. All the information below was extracted from the API webpage.
For the moment the results obtained from the API are in pure HTML, XML or JSON, so
you will have to parse them yourself. 

=head2 Search functionalities

=head3 Search Protein

Search proteins matching the query or search proteins for which the filter is true. 
Available filter values are: structure, disease, expression, mutagenesis or proteomics.
Note: only one filter parameter at a time is possible for the moment. 

	@result = $np->search_protein(-query => "kinase");
	@result = $np->search_protein(-query => "kinase", -filter => "disease");

=head3 Control Vocabulary Terms

Search control vocabulary terms matching the query or search control vocabulary terms in the category specified by the filter. 
Available filter values are: enzyme, go, mesh, cell, domain, family, tissue, metal, pathway, disease, keyword, ptm, subcell.
Note: only one category at a time is possible. 

	@result = $np->search_cv(-query => "colon");
	@result = $np->search_cv(-query => "colon", -filter => "keyword");

=head3 Format: 

Output format maybe in JSON (default), HTML or XML.	

	@result = $np->search_protein(-query => "kinase", -filter => "disease", -format => "html");


=head2 Find information by protein entry

=head3 Protein ID

ID is neXtProt identifier.
Retrieve gene name as well as main identifier and isoform sequences

	@result = $np->get_protein_info(-query => "NX_P13051");

=head3 Post-translational modifications

For each isoform of the specified entry, retrieve all post-translational modifications. 

	@result = $np->get_protein_info(-query => "NX_P13051", -retrieve => "ptm");

=head3 Variant

For each isoform of the specified entry, retrieve all variants. 

	@result = $np->get_protein_info(-query => "NX_P13051", -retrieve => "variant");

=head3 Localisation

For each isoform of the specified entry, retrieve all subcellular location. 

	@result = $np->get_protein_info(-query => "NX_P13051", -retrieve => "localisation");

=head3 Expression

Retrieve all expression information by tissue for the specified entry. 

	@result = $np->get_protein_info(-query => "NX_P13051", -retrieve => "expression");

=head3 Format: 

Output format maybe in JSON (default), HTML or XML.

	@result = $np->get_protein_info(-query => "NX_P13051", -retrieve => "expression", -format => "html");

=head2 Find information by isoform

=head3 Protein ID

ID is neXtProt identifier.
Retrieve gene name as well as main identifier and isoform sequences

	@result = $np->get_isoform_info(-query => "NX_O00142-2");

=head3 Post-translational modifications

For each isoform of the specified entry, retrieve all post-translational modifications. 

	@result = $np->get_isoform_info(-query => "NX_P01116-2", -retrieve => "ptm");

=head3 Variant

For each isoform of the specified entry, retrieve all variants. 

	@result = $np->get_isoform_info(-query => "NX_P01116-2", -retrieve => "variant");

=head3 Localisation

For each isoform of the specified entry, retrieve all subcellular location. 

	@result = $np->get_isoform_info(-query => "NX_P01116-2", -retrieve => "localisation");

=head2 Find information by controlled vocabulary term

=head3 Protein ID

ID is neXtProt identifier.
Retrieve the accession, the name and the category of the CV term.

	@result = $np->get_protein_cv_info(-query => "PTM-0205");

=head3 Protein List

List all the proteins associated with the term in neXtProt.

	@result = $np->get_protein_cv_info(-query => "PTM-0205", -retrieve => "proteins");

=head3 Format: 

Output format maybe in JSON (default), HTML or XML.

	@result = $np->get_protein_cv_info(-query => "PTM-0205", -retrieve => "proteins", -format => "html");

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

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

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live

=head1 AUTHOR - Felipe da Veiga Leprevost

Email leprevost@cpan.org

=cut

package Bio::DB::NextProt;

use strict;
use warnings;
use REST::Client;
sub new {
	my ($class, @args) = @_;
	#my $self = $class->SUPER::new(@args);
	my $self = {};
	$self->{_client}	= REST::Client->new({host=> "http://www.nextprot.org", timeout => 10,});
	$self->{_query}		= undef;
	$self->{_filter}	= undef;
	$self->{_format}	= "json";
	bless($self, $class);
	return $self;
}

sub search_protein() {
	
	my $self  = shift;
    my %param = @_;

	my $path = "/rest/protein/list";

    $self->{_format} = $param{'-format'} if defined $param{'-format'};

	if (defined $param{'-query'} && defined $param{'-filter'}) {
		
		$self->{_client}->GET($path."?query=".$param{'-query'}."&filter=".$param{'-filter'}."&format=".$self->{_format});

	} elsif (defined $param{'-query'}) {
		
		$self->{_client}->GET($path."?query=".$param{'-query'}."&format=".$self->{_format});

	} elsif (defined $param{'-filter'}) {
		
		$self->{_client}->GET($path."?filter=".$param{'-filter'}."&format=".$self->{_format});
	}

	&reset_params();

	return $self->{_client}->responseContent();

}

sub search_cv() {
	my $self  = shift;
	my %param = @_;
	
	my $path   = "/rest/cv/list";

    $self->{_format} = $param{'-format'} if defined $param{'-format'};

    if (defined $param{'-query'} && defined $param{'-filter'}) {

        $self->{_client}->GET($path."?query=".$param{'-query'}."&filter=".$param{'-filter'}."&format=".$self->{_format});

    } elsif (defined $param{'-query'}) {

        $self->{_client}->GET($path."?query=".$param{'-query'}."&format=".$self->{_format});

    } elsif (defined $param{'-filter'}) {

        $self->{_client}->GET($path."?filter=".$param{'-filter'}."&format=".$self->{_format});
    }

	&reset_params();

	return $self->{_client}->responseContent();

}

sub get_protein_info() {
	my $self  = shift;
	my %param = @_;

	my $path   = "/rest/entry/";

    $self->{_format} = $param{'-format'} if defined $param{'-format'};

	if (defined $param{'-query'} && $param{'-retrieve'}) {

		$self->{_client}->GET($path.$param{'-query'}."/".$param{'-retrieve'}."?format=".$self->{_format});

	} elsif (defined $param{'-query'}) {

		$self->{_client}->GET($path.$param{'-query'}."?format=".$self->{_format});
	}

	&reset_params();

	return $self->{_client}->responseContent();

}

sub get_isoform_info() {
	my $self  = shift;
	my %param = @_;

	my $path = "/rest/isoform/";

    $self->{_format} = $param{'-format'} if defined $param{'-format'};

    if (defined $param{'-query'} && $param{'-retrieve'}) {

        $self->{_client}->GET($path.$param{'-query'}."/".$param{'-retrieve'}."?format=".$self->{_format});

	} elsif (defined $param{'-query'}) {

	    $self->{_client}->GET($path.$param{'-query'}."?format=".$self->{_format});
	}

	&reset_params();

	return $self->{_client}->responseContent();

}

sub get_protein_cv_info() {
	my $self  = shift;
	my %param = @_;

	my $path = "/rest/cv/";

	$self->{_format} = $param{'-format'} if defined $param{'-format'};

	if (defined $param{'-query'} && $param{'-retrieve'}) {
		
		$self->{_client}->GET($path.$param{'-query'}."/".$param{'-retrieve'}."?format=".$self->{_format});

    } elsif (defined $param{'-query'}) {

        $self->{_client}->GET($path.$param{'-query'}."?format=".$self->{_format});
    }

	&reset_params();
    
	return $self->{_client}->responseContent();

}


sub reset_params() {
	my $self = shift;

	$self->{_query}  = undef;
	$self->{_filter} = undef;
}


1;





























