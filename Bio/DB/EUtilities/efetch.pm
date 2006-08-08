# $Id$
#
# BioPerl module for Bio::DB::EUtilities::efetch
#
# Cared for by Chris Fields
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
# 
# Part of the EUtilities BioPerl package

=head1 NAME

Bio::DB::EUtilities::efetch - retrieval of records from a list of IDs or the
user's environment.

=head1 SYNOPSIS

  my $efetch = Bio::DB::EUtilities->new(
                                       -verbose => 1,
                                       -cookie   => $esearch->next_cookie,
                                       -retmax   => $retmax,
                                       -rettype  => 'fasta'
                                        );
    
  print $efetch->get_response->content;

=head1 DESCRIPTION

L<EFetch|Bio::DB::EUtilities::efetch> retrieves data records from a list of
ID's.  This can be accomplished directly (using C<id>) or indirectly
(by using a L<Cookie|Bio::DB::EUtilities::Cookie>.

=head2 NCBI Efetch Parameters

The following are a general list of parameters that can be used to take
advantage of Efetch.  Up-to-date help for Efetch is available at this URL
(the information below is a summary of the options found there):

  http://eutils.ncbi.nlm.nih.gov/entrez/query/static/efetch_help.html

=over 3

=item C<db>

One or more database available through EUtilities.  EFetch currently only
supports database retrieval from the following databases:

B<pubmed>, B<pmc> (PubMed Central), B<journals>, B<omim>, B<nucleotide>,
B<protein>, B<genome>, B<gene>, B<snp> (dbSBP), B<popset>, and B<taxonomy>.

Also supported are B<sequences> (nucleotide, protein, popset and genome), and
the three subsets of nucleotide: B<nuccore>, B<nucest>, B<nucgss>

=item C<id>

a list of primary ID's

Below are a list of IDs which can be used with EFetch:

For sequence databases:

B<NCBI sequence number> (GI), B<accession>, B<accession.version>, B<fasta>,
B<GeneID>, B<genome ID>, B<seqid>

All other databases:

B<PMID> (pubmed), B<MIM number> (omim), B<GI number> (nucleotide, protein),
B<Genome ID> (genome), B<Popset ID> (popset), B<SNP cluster ID> (snp),
B<UniSTS ID> (unists), B<UniGene cluster ID> (unigene), <MMDB-ID> (structure),
B<PSSM-ID> (cdd), B<3D SDI> (domains), B<TAXID> (taxonomy), B<GEO ID> (geo)

=item C<mindate>, C<maxdate>

limits results by dates (C<yyyy/mm/dd> format, or by year)

=item C<rettype>

Output type based on the database.  Not all return types are compatible with
all return modes (-C<retmode>).  For more information, see the specific
literature or sequence database links at URL above.

Literature databases have the below return types:

B<uilist> (all databases),
B<abstract>, B<citation>, B<medline> (not omim),
B<full> (journals and omim)

Literature databases have the below return types:

B<native> (full record, all databases),
B<fasta>, B<seqid>, B<acc> (nucleotide or protein),
B<gb>, B<gbc>, B<gbwithparts> (nucleotide only),
B<est> (dbEST only),
B<gss> (dbGSS only),
B<gp>, B<gpc> (protein only),
B<chr>, B<flt>, B<rsr>, B<brief>, B<docset> (dbSNP only)

=item C<retmode>

EFetch is set, by default, to return a specific format for each Entrez database;
this is set in the %DATABASE hash in L<Bio::DB::EUtilities>.  To override this
format, you can set -C<retmode>.  The normal return modes are text, HTML, XML,
and ASN1.  Error checking for the set return mode is currently not
implemented.

=item C<report>

Used for the output format for Taxonomy; set to B<uilist>, B<brief>, B<docsum>,
B<xml>

=item C<strand> - I<sequence only>

The strand of DNA to show: 1=plus, 2=minus

=item C<seq_start>, C<seq_stop> - I<sequence only>

the start and end coordinates of the sequence to display

=item C<complexity> - I<sequence only>

The GI is often part of a biological blob containing other GIs

    * 0 - get the whole blob
    * 1 - get the bioseq for gi of interest (default in Entrez)
    * 2 - get the minimal bioseq-set containing the gi of interest
    * 3 - get the minimal nuc-prot containing the gi of interest
    * 4 - get the minimal pub-set containing the gi of interest

=back

=head2 Additional (Bioperl-related) Parameters

These are Bioperl-related settings and are not used as CGI parameters when

=over 3

=item C<eutil>

The relevant EUtility to be used (efetch).  

=item C<cookie>

Uses a L<Cookie|Bio::DB::EUtilities::Cookie>-based search (see below)

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the
evolution of this and other Bioperl modules. Send
your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation
is much appreciated.

  bioperl-l@lists.open-bio.org               - General discussion
  http://www.bioperl.org/wiki/Mailing_lists  - About the mailing lists
  
=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to
help us keep track the bugs and their resolution.
Bug reports can be submitted via the web.

  http://bugzilla.open-bio.org/

=head1 AUTHOR 

Email cjfields at uiuc dot edu

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut

# Let the code begin...

package Bio::DB::EUtilities::efetch;
use strict;
use warnings;
use Bio::DB::EUtilities;
use URI::Escape qw(uri_unescape);

use vars qw(@ISA $EUTIL);

@ISA = qw(Bio::DB::EUtilities);

BEGIN {
    $EUTIL = 'efetch';
}

sub _initialize {
    my ($self, @args ) = @_;
    $self->SUPER::_initialize(@args);
    my ($retmode, $reldate, $mindate, $maxdate, $datetype, $rettype, $retstart, 
        $retmax, $report, $seq_start, $seq_stop, $strand, $complexity) = 
      $self->_rearrange([qw(RETMODE RELDATE MINDATE MAXDATE DATETYPE RETTYPE
        RETSTART RETMAX REPORT SEQ_START SEQ_STOP STRAND COMPLEXITY)], @args);    
    # set by default
    $self->_eutil($EUTIL);
    $datetype ||= 'mdat';
    $self->datetype($datetype) if $datetype;
    defined($retstart)       && $self->retstart($retstart);
    $retmode        && $self->retmode($retmode);
    $retmax         && $self->retmax($retmax);
    $rettype        && $self->rettype($rettype);
    $seq_start      && $self->seq_start($seq_start);
    $seq_stop       && $self->seq_stop($seq_stop);
    $strand         && $self->strand($strand);
    defined($complexity) && $self->complexity($complexity);
    $report         && $self->report($report);    
}

=head2 parse_response

 Title   : parse_response
 Usage   : $db->_parse_response($content)
 Function: parse out response for cookie
 Returns : empty
 Args    : none
 Throws  : 'unparseable output exception'

=cut

# this is NOOP b/c efetch returns raw data to be processed or saved;
# HTTP errors caught in get_response 

sub parse_response {
}

1;
__END__