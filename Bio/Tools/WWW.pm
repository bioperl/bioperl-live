#-----------------------------------------------------------------------------
# PACKAGE : Bio::Tools::WWW.pm
# PURPOSE : To encapsulate commonly used URLs for web key websites in bioinformatics.
# AUTHOR  : Steve A. Chervitz
# CREATED : 27 Aug 1996 
# REVISION: $Id$
#
# For documentation, run this module through pod2html 
# (preferably from Perl v5.004 or better).
#
# MODIFIED: 
#  0.014, sac --- Mon Aug 31 19:41:44 1998
#      * Updated and added a few URLs.
#      * Added method strip_html().
#      * Documentation changes.
#
#-----------------------------------------------------------------------------

package	 Bio::Tools::WWW;
use strict;  

use Bio::Root::Object   ();
use Bio::Root::Global  qw($AUTHORITY);
use Exporter      ();

use vars qw(@ISA @EXPORT_OK %EXPORT_TAGS $ID $VERSION $BioWWW $Revision);

@ISA         = qw( Bio::Root::Object Exporter);
@EXPORT_OK   = qw($BioWWW);
%EXPORT_TAGS = ( obj => [qw($BioWWW)],
		 std => [qw($BioWWW)]);

$ID = 'Bio::Tools::WWW';
$VERSION = 0.014;
$Revision = '$Id$'; #'

## Static object.
$BioWWW = {};
bless $BioWWW, $ID;
$BioWWW->{'_name'} = "Static $ID object";


## POD Documentation:

=head1 NAME

Bio::Tools::WWW.pm - Bioperl manager for web resources related to biology.

=head1 SYNOPSIS

=head2 Object Creation

    use Bio::Tools qw(:obj);

    $pdb = $BioWWW->home_url('pdb');

There is no need to create a new Bio::Tools::WWW.pm object when the
C<:obj> tag is used. This tag will import the static $BioWWW object
created by Bio::Tools::WWW.pm into your name space. This saves you
from having to call C<new Bio::Tools::WWW>.

You are free to not use the :obj tag and create the object as you
like, but a Bio::Tools::WWW object is not configurable; any given
script only needs a single copy.

=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

You also need to define URLs for the following variables in this package:

  $Not_found_url : Generic page to show in place of a 404 error.
  $Tmp_url       : Web-accessible site that is Used for scripts that 
                   need to generate temporary, web-accessible files.
                   The files need not necessarily be HTML files, but 
                   being on the same disk as the server will permit 
                   faster IO from server scripts.

=head1 DESCRIPTION

Bio::Tools::WWW is primarily a URL broker for a select set 
of sites related to bioinformatics/genome analysis. It 
definitely represents a biased, unexhaustive set.
It might be more accurate to call this module 
"Bio::Tools::URL.pm". But this module does handle some non-URL
things and it may do more of this in the future. Having one
module to cover all biologically relevant web utilities
makes it more convenient, especially at this early stage
of development. 

Maintaining accurate URLs over time can be challenging as 
new web sites spring up and old sites are re-organized. Because
of this fact, the URLs in this module are not guaranteed to be
correct or exhaustive and will require periodic updating.

=head2 URL Management

By keeping URL management within Bio::Tools::WWW.pm, other generic 
modules can easily access a variety of different web sites without 
having to know about a potential multitude of specific modules 
specialized for one database or another. A specific example
of this is in B<Bio::Tools::Blast.pm> where the function blast_to_html()
needs access to different URLs in order to add database links
to the Blast report. An alternative approach would be to have
multiple blast_to_html() functions defined within modules
specialized for Blast analyses of different datasets. This, however,
may create maintenance headaches when updating the different
versions of the function. 

=head2 Complex Websites

Websites with complex datasets may require special treatment
within this module. As an example,
URLs for the Saccharomyces Genome Database are clustered
separately in this module, due to (1) the different ways to
access information at this database and (2) the familiarity 
of the developer with this database. The Bio::SGD::WWW.pm inherits from
Bio::Tools::WWW.pm to permit access to the URLs provided by Bio::Tools::WWW.pm
and to SGD-specific HTML and images. 

The organization of Bio::Tools::WWW.pm is expected to evolve as 
websites get born, die, and mutate their APIs.

=head1 SEE ALSO

 http://bio.perl.org/Projects/modules.html  - Online module documentation
 http://bio.perl.org/                       - Bioperl Project Homepage

=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other Bioperl modules.
Send your comments and suggestions preferably to one of the Bioperl mailing lists.
Your participation is much appreciated.

    vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
    vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
    http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the bugs and 
their resolution. Bug reports can be submitted via email or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR

Steve A. Chervitz, sac@genome.stanford.edu

=head1 VERSION

Bio::Tools::WWW.pm, 0.014

=head1 COPYRIGHT

Copyright (c) 1996-98 Steve A. Chervitz. All Rights Reserved.
This module is free software; you can redistribute it and/or 
modify it under the same terms as Perl itself.


=cut


#
##
###
#### END of main POD documentation.
###
##
#


############################  DATA ##################################

### Database homepage links.
my %Home_url =
    (
     'bioperl'         =>'http://bio.perl.org/',
     'bioperl-stanford'=>'http://genome-www.stanford.edu/perlOOP/bioperl/',
     'bioperl-schema'  =>'http://bio.perl.org/Projects/Schema/',
     'biomoo'          =>'http://bioinformatics.weizmann.ac.il/BioMOO/',
     'blast_ncbi'      =>'http://www.ncbi.nlm.nih.gov/BLAST/',
     'blast_wu'        =>'http://blast.wustl.edu/',
     'bsm'             =>'http://www.biochem.ucl.ac.uk/bsm/',
     'clustal'         =>'http://www.csc.fi/molbio/progs/clustalw/clustalw.html',
     'ebi'             =>'http://www.ebi.ac.uk/',
     'emotif'          =>'http://motif.Stanford.EDU/emotif',
     'entrez'          =>'http://www3.ncbi.nlm.nih.gov/Entrez/',
     'expasy'          =>'http://www.expasy.ch/',
     'gdb'             =>'http://www.gdb.org/',  # R.I.P. (Jan 1998); site still functional
     'mips'            =>'http://speedy.mips.biochem.mpg.de/',
     'mmdb'            =>'http://www.ncbi.nlm.nih.gov/Structure/',
     'modbase'         =>'http://guitar.rockefeller.edu/',
     'ncbi'            =>'http://www.ncbi.nlm.nih.gov/',
     'pedant'          =>'http://pedant.mips.biochem.mpg.de',
     'phylip'          =>'http://evolution.genetics.washington.edu/phylip.html',
     'pir'             =>'http://www-nbrf.georgetown.edu/pir/',
     'pfam'            =>'http://pfam.wustl.edu/',
     'pfam_uk'         =>'http://www.sanger.ac.uk/Software/Pfam/',
     'pfam_us'         =>'http://pfam.wustl.edu/',
     'pdb'             =>'http://www.pdb.bnl.gov/',
     'presage'         =>'http://presage.stanford.edu/',
     'geneQuiz'        =>'http://www.sander.ebi.ac.uk/genequiz/genomes/sc/',
     'molMov'          =>'http://bioinfo.mbb.yale.edu/MolMovDB/',
#     'protMot'         =>'http://bioinfo.mbb.yale.edu/ProtMotDB/', # old, use molMov instead
     'pubmed'          =>'http://www.ncbi.nlm.nih.gov/PubMed/',
     'sacch3d'         =>'http://genome-www.stanford.edu/Sacch3D/',
     'sgd'             =>'http://genome-www.stanford.edu/Saccharomyces/',
#     'scop'            =>'http://www.pdb.bnl.gov/scop/',
     'scop'            =>'http://scop.stanford.edu/scop/',
     'swissProt'       =>'http://www.expasy.ch/sprot/sprot-top.html',
     'webmol'          =>'http://genome-www.stanford.edu/structure/webmol/',
     'ypd'             =>'http://quest7.proteome.com/YPDhome.html',
     );

### Database access CGI stems. (For some DBs the home URL can be used as the CGI stem)
my %Stem_url = 
    ( 
      'emotif'      =>'http://dna.Stanford.EDU/cgi-bin/emotif/',
      'entrez'      =>'http://www3.ncbi.nlm.nih.gov/htbin-post/Entrez/query?',
      'pdb'         =>'http://www.pdb.bnl.gov/pdb-bin/',
      'pfam_uk'     =>'http://www.sanger.ac.uk/cgi-bin/Pfam/',
      'pfam_us'     =>'http://pfam.wustl.edu/cgi-bin/',
      'pir'         =>'http://www-nbrf.georgetown.edu/cgi-bin/nbrfget?',
      );


### Database access stems/links.
my %Search_url = 
    ( #'3db'       =>'http://pdb.pdb.bnl.gov/cgi-bin/pdbids?3DB_ID=',   # Former stem
      '3db'          =>$Stem_url{'pdb'}.'opdbshort?oPDBid=',  # New stem (aug 1997)
      'embl'         =>$Home_url{'ebi'}.'htbin/emblfetch?',
      'expasy'       =>$Home_url{'expasy'}.'cgi-bin/',  # program name and query string must be supplied.
      'cath'         =>$Home_url{'bsm'}.'cath/CATHSrch.pl?type=PDB&query=',
      'cog_seq'      =>$Home_url{'ncbi'}.'cgi-bin/COG/nph-cognitor?seq=', # add sequence
      # To cog_orf, append ORF name ('YAL005c'). Case-sensitive! YAL005C won't work!
      'cog_orf'      =>$Home_url{'ncbi'}.'cgi-bin/COG/cogeseq?', 
      'ec1'          =>$Home_url{'gdb'}.'bin/bio/wais_q-bio?object_class_key=30&jhu_id=',
      'ec2'          =>$Home_url{'bsm'}.'enzymes/',
      'ec3'          =>$Home_url{'expasy'}.'cgi-bin/get-enzyme-entry?',
      'emotif_id'    =>$Stem_url{'emotif'}.'nph-identify?sequence=',
      'entrez'       =>$Stem_url{'entrez'}."db=p_r?db=1&choseninfo=ORF_NAME%20[Gene%20Name]\@1\@1&form=4&field=Gene%20Name&mode=0&retrievestring=ORF_NAME%20[Gene%20Name]",
      'gb_n'         =>$Stem_url{'entrez'}."db=n&form=6&dopt=g&uid=",
      'gb_p'         =>$Stem_url{'entrez'}."db=p&form=6&dopt=g&uid=",
      'gb_struct'    =>$Stem_url{'entrez'}."db=t&form=6&dopt=s&uid=",
      'pdb'          =>$Stem_url{'pdb'}.'send-text?filename=',
      'medline'      =>$Stem_url{'entrez'}.'form=6&db=m&Dopt=r&uid=',
      'mmdb'         =>$Stem_url{'entrez'}.'db=t&form=6&Dopt=s&uid=',
      'modbase_orf'  =>$Home_url{'modbase'}.'gm-cgi-bin/orf_page.cgi?pg1=0.5&pg2=1.0&orf=',
      # To the modbase_model, append yeast ORF name &pdb=<4-LETTER_CODE>&chain=<UPCASE LETTER, IF ANY>
      'modbase_model' =>$Home_url{'modbase'}.'gm-cgi-bin/model_page.cgi?pg1=0.5&pg2=1.0&orf=',
      'molMov'       =>$Home_url{'molMov'}.'search.cgi?pdb=',
      'pdb'          =>$Stem_url{'pdb'}.'opdbshort?oPDBid=',  # same as 3db
      'pdb_coord'    =>$Stem_url{'pdb'}.'send-pdb?filename=', # retrieves full coordinate file
      'pfam'         =>$Home_url{'pfam'}.'cgi-bin/nph-hmm_search?evalue=1.0&protseq=',  # default: seq search, US
      'pfam_sp_uk'   =>$Stem_url{'pfam_uk'}.'swisspfamget.pl?name=',
      'pfam_seq_uk'  =>$Stem_url{'pfam_uk'}.'nph-search.cgi?evalue=1.0&type=normal&protseq=',
      'pfam_sp_us'   =>$Stem_url{'pfam_us'}.'getswisspfam?key=',
      'pfam_seq_us'  =>$Stem_url{'pfam_us'}.'nph-hmm_search?evalue=1.0&protseq=',
      'pfam_form'    =>$Home_url{'pfam'}.'cgi-bin/hmm_page.cgi', # interactive search form
      'pir_id'       =>$Stem_url{'pir'}.'fmt=c&xref=0&id=',
      'pir_acc'      =>$Stem_url{'pir'}.'fmt=c&xref=1&id=',
      'pir_uid'      =>$Stem_url{'pir'}.'uid=',
      'pdbSum'       =>$Home_url{'bsm'}.'cath/GetPDBSUMCODE.pl?code=',
#      'protMot'      =>$Home_url{'protMot'}.'search.cgi?pdb=', # old, use molMov instead
      'presage_sp'   =>$Home_url{'presage'}.'search.cgi?spac=',
      'swpr'         =>$Home_url{'expasy'}.'cgi-bin/get-sprot-entry?',
      'swModel'      =>$Home_url{'expasy'}.'cgi-bin/sprot-swmodel-sub?',
      'swprSearch'   =>$Home_url{'expasy'}.'cgi-bin/sprot-search-ful?',
      
      ###  SCOP tlev options can be appended to the stem after adding a PDB ID.
      ###  tlev options are: 'dm'(domain), 'sf'(superfamily), 'fa'(family), 'cf'(common fold), 'cl'(class)
      ###  E.g., search.cgi?pdb=1ARD;tlev=dm

      'scop'         =>$Home_url{'scop'}.'search.cgi?pdb=',  ### better to use scop_pdb.
      'scop_pdb'     =>$Home_url{'scop'}.'search.cgi?pdb=',
      'scop_data'    =>$Home_url{'scop'}.'data/scop.',  ### Deprecated: frequent changes.

      ## Search URLs for SGD/Sacch3D are contained %SGD_url and %S3d_url (below).

      # For wormpep, the query string MUST end with "&keyword=" (after appending a sequence ID)
      'wormpep'        =>'http://www.sanger.ac.uk/cgi-bin/wormpep_fetch.pl?entry=', 
      'wormace'        =>'http://webace.sanger.ac.uk/cgi-bin/webace?db=wormace&class=Sequence&text=yes&object=',

      ### YPD: You must use a valid gene name or ORF name (IFF there is no gene name).
      ###      For this reason it is most convenient to use SGD's Protein_Info link
      ###      which can accept either and will provide a proper link to YPD.
      'ypd'          =>'http://quest7.proteome.com/YPD/',  
      );



### CGI stems for SGD and Sacch3D.
my %SGD_stem_url =
    ('stanford'      =>'http://genome-www.stanford.edu/',
     'sgd'           =>'http://genome-www.stanford.edu/cgi-bin/SGD/',  
     'sgd2'          =>'http://genome-www2.stanford.edu/cgi-bin/SGD/', 
     's3d'           =>'http://genome-www.stanford.edu/cgi-bin/SGD/Sacch3D/',  
     's3d2'          =>'http://genome-www2.stanford.edu/cgi-bin/SGD/Sacch3D/',  
     's3d3'          =>'http://genome-www3.stanford.edu/cgi-bin/SGD/Sacch3D/',  
     'sacchdb'       =>'http://genome-www.stanford.edu/cgi-bin/dbrun/SacchDB?',  
     );

### SGD stems and links.
my %SGD_url = 
    ('home'         =>$Home_url{'sgd'},
     'help'         =>$Home_url{'sgd'}.'help/',
     'mammal'       =>$Home_url{'sgd'}.'mammal/',  
     'worm'         =>$Home_url{'sgd'}.'worm/',  
     'gene'         =>$SGD_stem_url{'sacchdb'}.'find+Locus+',
     'locus'        =>$SGD_stem_url{'sacchdb'}.'find+Locus+',
     'orf'          =>$SGD_stem_url{'sacchdb'}.'find+Locus+',
     'mipsorf'      =>$SGD_stem_url{'sgd'}."mips-orfs?",
     'gene_info'    =>$SGD_stem_url{'sacchdb'}.'find+Gene_Info+',
     'prot_info'    =>$SGD_stem_url{'sacchdb'}.'find+Protein_Info+',
     'seq'          =>$SGD_stem_url{'sgd'}.'seqDisplay?seq=',
     'gi'           =>$SGD_stem_url{'sacchdb'}.'find+Sequence+Database+=+GenPept+AND+NEXT+=+',
     'chr'          =>$SGD_stem_url{'sgd2'}.'seqTools?chr=',
     'chr_old'      =>$SGD_stem_url{'sgd'}.'dnaredir?chr=',
     'seq_an'       =>$SGD_stem_url{'sgd2'}.'seqTools?seqname=',
     'seq_an_old'   =>$SGD_stem_url{'sgd'}.'dnaredir?seqname=',
     'map_chr'      =>$SGD_stem_url{'sgd'}.'ORFMAP/ORFmap?chr=',
     'map_orf'      =>$SGD_stem_url{'sgd'}.'ORFMAP/ORFmap?seq=',
#     'chr'          =>$SGD_stem_url{'sgd2'}.'seqform?chr=',
#     'seg'          =>$SGD_stem_url{'sgd2'}.'seqform?seg=',
#     'fea'          =>$SGD_stem_url{'sgd2'}.'featureform?seg=',
     'feature'      =>$SGD_stem_url{'sgd2'}.'featureform?chr=', # complete with "5&beg=100&end=400"
     'search'       =>$SGD_stem_url{'sgd'}.'search?',
     'images'       =>$SGD_stem_url{'stanford'}.'images/',
     'suggest'      =>$SGD_stem_url{'stanford'}.'forms/sgd-suggestion.html',
     'tmp'          =>$SGD_stem_url{'stanford'}.'tmp/',
     );


### Sacch3D stems and links.
my %S3d_url =
    ('home'          =>$Home_url{'sacch3d'},
     'search'        =>$Home_url{'sacch3d'}.'search.html',
     'help'          =>$Home_url{'sacch3d'}.'help/',
     'new'           =>$Home_url{'sacch3d'}.'new/',
     'chrm'          =>$Home_url{'sacch3d'}.'data/chr',  
     'domains'       =>$Home_url{'sacch3d'}.'domains/',  
     'genequiz'      =>$Home_url{'sacch3d'}.'genequiz/',  
     'analysis'      =>$Home_url{'sacch3d'}.'analysis/',  
     'scop'          =>$SGD_stem_url{'s3d3'}.'getscop?data=',  
     'scop_fold'     =>$SGD_stem_url{'s3d3'}.'getscop?type=fold&data=',  
     'scop_class'    =>$SGD_stem_url{'s3d3'}.'getscop?type=class&data=',  
     'scop_gene'     =>$SGD_stem_url{'s3d3'}.'getscop?type=gene&data=',  
     'gene'          =>$SGD_stem_url{'s3d'}.'get?class=gene&item=',
     'orf'           =>$SGD_stem_url{'s3d'}.'get?class=orf&item=',
     'text'          =>$SGD_stem_url{'s3d'}.'get?class=text&item=',
     'pdb'           =>$SGD_stem_url{'s3d'}.'get?class=pdb&item=',
     'pdb_coord'     =>$SGD_stem_url{'s3d'}.'pdbcoord.pl?id=',
     'dsc'           =>$SGD_stem_url{'s3d'}.'dsc.pl?gene=',
     'emotif'        =>$SGD_stem_url{'s3d'}.'seq_search.pl?db=emotif&gene=',
     'pfam'          =>$SGD_stem_url{'s3d'}.'seq_search.pl?db=pfam&gene=',
     'pfam_uk'       =>$SGD_stem_url{'s3d'}.'seq_search.pl?db=pfam&loc=uk&gene=',
     'pfam_us'       =>$SGD_stem_url{'s3d'}.'seq_search.pl?db=pfam&loc=us&gene=',
     'blast_pdb'     =>$SGD_stem_url{'s3d'}.'getblast?db=pdb&name=',
     'blast_nr'      =>$SGD_stem_url{'s3d'}.'getblast?db=nr&name=',
     'blast_est'     =>$SGD_stem_url{'s3d'}.'getblast?db=est&name=',
     'blast_mammal'  =>$SGD_stem_url{'s3d'}.'getblast?db=mammal&name=',
     'blast_human'   =>$SGD_stem_url{'s3d'}.'getblast?db=human&name=',
     'blast_worm'    =>$SGD_stem_url{'s3d'}.'getblast?db=worm&name=',
     'blast_yeast'   =>$SGD_stem_url{'s3d'}.'getblast?db=yeast&name=',
     'blast_worm_yeast'=>$SGD_stem_url{'s3d'}.'getblast?db=worm&query=worm&name=',
     'patmatch'      =>$SGD_stem_url{'s3d2'}.'grepmatch?',  ## deprecated
     'grepmatch'     =>$SGD_stem_url{'s3d2'}.'grepmatch?',
     'pdb_neighbors' =>$SGD_stem_url{'s3d'}.'pdb_neighbors?id=CHAIN&gene=ORF_NAME',
     );


### 3D viewer stems.
my %Viewer_url = 
#    ('java'     =>$SGD_stem_url{'sgd'}.'Sacch3D/pdbViewer.pl?pdbCode=PDB&orf=',
    (
     'java'     =>$SGD_stem_url{'sgd'}.'Sacch3D/pdbViewer.pl?pdbCode=',  # Default java viewer
     'webmol'   =>$SGD_stem_url{'sgd'}.'Sacch3D/pdbViewer.pl?pdbCode=', 
     'codebase' =>$SGD_stem_url{'stanford'}.'structure/webmol/lib',
     'rasmol'   =>$Stem_url{'pdb'}.'send-ras?filename=',
     'chime'    =>$Stem_url{'pdb'}.'ccpeek?id=',
     'cn3d'     =>$Stem_url{'entrez'}.'db=t&form=6&Dopt=i&Complexity=Cn3D+Subset&uid=',
     'kinemage' =>'http://prosci.org/Kinemage',
     );


### Stock HTML
# The error reporting HTML strings represent some experiments in human psychology: 
# how do you induce users to report errors that you should know about yet not
# get flooded with trivial problems caused by novices?
my %Html = 
    ('authority'  =>qq|<A HREF="mailto:$AUTHORITY"><b>$AUTHORITY</b></A>|,
     'trouble'    => <<"QQ_TROUBLE_QQ",
<p>If this problem persists, <A HREF="mailto:$AUTHORITY"><b>please notify us.</b></A>
Include a copy of this error page with your message. Thanks.<p>
QQ_TROUBLE_QQ
     'notify'     => <<"QQ_NOTIFY_QQ",
<A HREF="mailto:$AUTHORITY"><b>Please notify us.</b></A>
Include a copy of this error page with your message. Thanks.<p>
QQ_NOTIFY_QQ
     'ourFault'   => <<"QQ_FAULT_QQ",
<p><b>This is our fault!</b> There is apparently a problem with our software
that we may not know about. <A HREF="mailto:$AUTHORITY"><b>Please notify us!</b></A>
Include a copy of this error page with your message. Thanks.<p>
QQ_FAULT_QQ
     'techDiff'   => <<"QQ_TECH_QQ",
<p><big>We are experiencing technical difficulties now.<br>
We will have the problem fixed soon. Sorry for any inconvenience.</big><p>
QQ_TECH_QQ

     );


### Miscellaneous URLs. Configure as desired for your site.
my $Not_found_url = 'http://genome-www.stanford.edu/Sacch3D/notfound.html';
my $Tmp_url       = 'http://genome-www.stanford.edu/tmp/';



=head1 APPENDIX

Methods beginning with a leading underscore are considered private
and are intended for internal use by this module. They are
B<not> considered part of the public interface and are described here
for documentation purposes only.

=cut

#########################################################################
##                          ACCESSOR METHODS                            
#########################################################################


=head2 home_url

 Usage     : $BioWWW->home_url(<string>)
 Purpose   : To obtain the homepage URL for a biological database or resource.
 Returns   : String containing the URL (including "http://")
 Argument  : String
           : Currently acceptable arguments are:
           :    bioperl  bioperl-schema  biomoo  bsm  ebi  emotif  entrez 
           :    expasy  mips  mmdb  ncbi  pir  pfam  pdb  geneQuiz  
           :    molMov  pubmed  sacch3d  sgd  scop  swissProt  webmol  ypd
 Throws    : Warns if argument cannot be resolved to a URL.
 Comments  : The URLs listed here do not represent a complete list.
           : Expect this to evolve and grow with time.

See Also   : L<search_url>()

=cut

#-------------
sub home_url { 
#-------------
    my($self,$arg) = @_; 
    $arg eq 'all' and return %Home_url;
    (exists $Home_url{$arg}) ? $Home_url{$arg} 
                             : ($self->warn("Can't resolve argument to URL: $arg"), 
				$Not_found_url);
}



=head2 search_url

 Usage     : $BioWWW->search_url(<string>)
 Purpose   : To provide a URL stem for a search engine at a biological database 
           : or resource.
 Returns   : String containing the URL (including "http://")
 Argument  : String
           : Currently acceptable arguments are:
           :   3db  embl  cath  ec1  ec2  ec3  emotif_id  entrez  gb1  gb2  
           :   gb3  gb4  gb5  pdb  medline  mmdb  pdb  pdb_coord  pfam  pir_acc  
           :   pdbSum  molMov  swpr  swModel  swprSearch  scop  scop_pdb  scop_data 
           :   ypd
 Throws    : Warns if argument cannot be resolved to a URL.
 Comments  : Unlike the homepage URLs, this method does not return a complete
           : URL but a stem which must be further modified, typically by
           : appending data to it, before it can be used. The data appended
           : depends on the specific URL; typically, it is a database ID or
           : other unique identifier.
           : The requirements for each URL will be described here eventually.
           : 
           : The URLs listed here do not represent a complete list.
           : Expect this to evolve and grow with time.
           :
           : Given this complexity, it may be useful to provide special methods
           : for these different URLs. This would however result in an 
           : explosion of methods that might make this module less 
           : maintainable and harder to use.

See Also   : L<home_url>()

=cut

#--------------
sub search_url { 
#--------------
    my($self,$arg) = @_; 
    $arg eq 'all' and return %Search_url;
    (exists $Search_url{$arg}) ? $Search_url{$arg} 
                             : ($self->warn("Can't resolve argument to URL: $arg"), 
				$Not_found_url);
}



=head2 stem_url

 Usage     : $BioWWW->stem_url(<string>)
 Purpose   : To obtain the minimal stem URL for searching a biological database or resource.
 Returns   : String containing the URL (including "http://")
 Argument  : String
           : Currently acceptable arguments are:
           :    emotif  entrez  pdb
 Throws    : Warns if argument cannot be resolved to a URL.
 Comments  : The URLs stems returned by this method are much more minimal than
           : this provided by search_url(). Use of these stems requires knowledge
           : of the CGI scripts which they invoke.

See Also   : L<search_url>()

=cut

#--------------
sub stem_url {
#--------------
    my($self,$arg) = @_; 
    $arg eq 'all' and return %Stem_url;
    (exists $Stem_url{$arg}) ? $Stem_url{$arg}
                             : ($self->warn("Can't resolve argument to URL: $arg"), 
				$Not_found_url);
}

	      

=head2 viewer_url

 Usage     : $BioWWW->viewer_url(<string>)
 Purpose   : To obtain the stem URL for a 3D viewer (RasMol, WebMol, Cn3D)
 Returns   : String containing the URL (including "http://")
 Argument  : String
           : Currently acceptable arguments are:
           :    rasmol webmol cn3d java  (java is an alias for webmol)
 Throws    : Warns if argument cannot be resolved to a URL.
 Comments  : The 4-letter Brookhaven PDB identifier must be appended to the
           : URL provided by this method.
           : The URLs listed here do not represent a complete list.
           : Expect this to evolve and grow with time.

=cut

#---------------
sub viewer_url { 
#---------------
    my($self,$arg) = @_; 
    $arg eq 'all' and return %Viewer_url;
    (exists $Viewer_url{$arg}) ? $Viewer_url{$arg} 
                             : ($self->warn("Can't resolve argument to URL: $arg"), 
				$Not_found_url);
}



=head2 not_found_url

 Usage     : $BioWWW->not_found_url()
 Purpose   : To obtain the URL for a web page to be shown in place of a 404 error.
 Returns   : String containing the URL (including "http://")
 Argument  : n/a
 Throws    : n/a
 Comments  : This URL should be customized as desired.

=cut

#-----------------
sub not_found_url {  my $self = shift; $Not_found_url; }
#-----------------


=head2 tmp_url

 Usage     : $BioWWW->tmp_url()
 Purpose   : To obtain the URL for a temporary, web-accessible directory.
 Returns   : String containing the URL (including "http://")
 Argument  : n/a
 Throws    : n/a
 Comments  : This URL should be customized  as desired.

=cut

#-----------
sub tmp_url {  my $self = shift; $Tmp_url; }
#-----------



=head2 search_link

 Usage     : $BioWWW->search_link(<site>, <value>, <text>)
 Purpose   : Wrapper for search_url() that returns the URL within an HTML anchor.
 Returns   : String containing the HTML anchor ( qq|<A HREF="http://..."</A>|)
 Argument  : <site>  = string to be used as argument for search_url()
           : <value> = string to be appended to the search URL stem.
           : <text>  = string to be shown as the link text (default = <value>).
 Throws    : n/a
 Status    : Experimental

See Also   : L<search_url>()

=cut

#---------------
sub search_link { 
#---------------
    my($self,$arg,$value,$text) = @_; 
    my $url = $self->search_url($arg);
    $text ||= $value;
    qq|<A HREF="$url$value">$text</A>|;
}



=head2 viewer_link

 Usage     : $BioWWW->viewer_link(<site>, <value>, <text>)
 Purpose   : Wrapper for viewer_url() that returns the complete URL within an HTML anchor.
 Returns   : String containing the HTML anchor ( qq|<A HREF="http://..."</A>|)
 Argument  : <site>  = string to be used as argument for viewer_url()
           : <value> = string to be appended to the viewer URL stem.
           : <text>  = string to be shown as the link text (default = <value>).
 Throws    : n/a
 Status    : Experimental

See Also   : L<viewer_url>()

=cut

#----------------
sub viewer_link { 
#----------------
    my($self,$arg,$value,$text) = @_; 
    my $url = $self->viewer_url($arg);
    $text ||= $value;
    qq|<A HREF="$url$value">$text</A>|;
}



=head2 html

 Usage     : $BioWWW->html(<string>)
 Purpose   : To obtain HTML-formatted text for frequently needed web-page messages.
 Returns   : String containing the HTML anchor ( qq|<A HREF="http://..."</A>|)
 Argument  : String.
           : Currently acceptable arguments are:
           :   authority  (mailto: link for webmaster; shows e-mail address as link)
           :   notify     (wraps mailto:authority link with text for link "please notify us")
           :   ourFault   ("this problem is our fault. If it persists <notify-link>")
           :   trouble    (same as ourFault but doesn't blame us for the problem)
           :   techDiff   ("we are experiencing technical difficulties. Please stand by.")
 Throws    : n/a
 Comments  : The authority (webmaster) is imported from the Bio::Root::Global.pm
           : module. The value for $AUTHORITY should be set there, or
           : customize this module so that it doesn't use Bio::Root::Global.pm.

=cut

#----------
sub html { 
#----------
    my($self,$arg) = @_; 
    $arg eq 'all' and return %Html;
    (exists $Html{$arg}) ? $Html{$arg} : "<pre>(missing HTML for \"$arg\")</pre>";
}


###
### Below are accessors specialized for the Saccharomyces Genome Database
### It is possible that they will be moved to Bio::SGD::WWW.pm in the future.
### 


=head2 sgd_url

 Usage     : $BioWWW->sgd_url(<string>)
 Purpose   : To obtain the webpage URL or search stem for SGD.
 Returns   : String containing the URL (including "http://")
 Argument  : String
           : Currently acceptable arguments (TODO).
 Throws    : Warns if argument cannot be resolved to a URL.
 Comments  : This accessor is specialized for the Saccharomyces Genome Database.
           : It is possible that it will be moved to SGD::WWW.pm in the future.

See Also   : L<search_url>()

=cut

#------------
sub sgd_url { 
#------------
    my($self,$arg) = @_; 
    $arg eq 'all' and return %SGD_url;
    (exists $SGD_url{$arg}) ? $SGD_url{$arg} 
                             : ($self->warn("Can't resolve argument to URL: $arg"), 
				$Not_found_url);
}



=head2 s3d_url

 Usage     : $BioWWW->s3d_url(<string>)
 Purpose   : To obtain the webpage URL or search stem for Sacch3D.
 Returns   : String containing the URL (including "http://")
 Argument  : String
           : Currently acceptable arguments (TODO).
 Throws    : Warns if argument cannot be resolved to a URL.
 Comments  : This accessor is specialized for the Saccharomyces Genome Database.
           : It is possible that it will be moved to SGD::WWW.pm in the future.

See Also   : L<search_url>()

=cut

#-----------
sub s3d_url { 
#-----------
    my($self,$arg) = @_; 
    $arg eq 'all' and return %S3d_url;
    (exists $S3d_url{$arg}) ? $S3d_url{$arg} 
                             : ($self->warn("Can't resolve argument to URL: $arg"), 
				$Not_found_url);
}



=head2 sgd_stem_url

 Usage     : $BioWWW->sgd_stem_url(<string>)
 Purpose   : To obtain the minimal stem URL for a SGD/Sacch3D CGI script.
 Returns   : String containing the URL (including "http://")
 Argument  : String
           : Currently acceptable arguments (TODO).
 Throws    : Warns if argument cannot be resolved to a URL.
 Comments  : This accessor is specialized for the Saccharomyces Genome Database.
           : It is possible that it will be moved to SGD::WWW.pm in the future.

See Also   : L<search_url>()

=cut

#-----------------
sub sgd_stem_url { 
#-----------------
    my($self,$arg) = @_; 
    $arg eq 'all' and return %SGD_stem_url;
    (exists $SGD_stem_url{$arg}) ? $SGD_stem_url{$arg} 
                             : ($self->warn("Can't resolve argument to URL: $arg"), 
				$Not_found_url);
}



=head2 s3d_link

 Usage     : $BioWWW->s3d_link(<site>, <value>, <text>)
 Purpose   : Wrapper for s3d_url() that returns the complete URL within an HTML anchor.
 Returns   : String containing the URL (including "http://")
 Argument  : <site>  = string to be used as argument for s3d_url()
           : <value> = string to be appended to the s3d URL stem.
           : <text>  = string to be shown as the link text (default = <value>).
 Throws    : n/a
 Status    : Experimental
 Comments  : This accessor is specialized for the Saccharomyces Genome Database.
           : It is possible that it will be moved to SGD::WWW.pm in the future.

See Also   : L<s3d_url>(), L<sgd_link>()

=cut

#--------------
sub s3d_link { 
#--------------
    my($self,$arg,$value,$text) = @_; 
    my $url = $self->s3d_url($arg);
    $text ||= $value;
    qq|<A HREF="$url$value">$text</A>|;
}
	      


=head2 sgd_link

 Usage     : $BioWWW->sgd_link(<site>, <value>, <text>)
 Purpose   : Wrapper for sgd_url() that returns the complete URL within an HTML anchor.
 Returns   : String containing the URL (including "http://")
 Argument  : <site>  = string to be used as argument for sgd_url()
           : <value> = string to be appended to the sgd URL stem.
           : <text>  = string to be shown as the link text (default = <value>).
 Throws    : n/a
 Status    : Experimental
 Comments  : This accessor is specialized for the Saccharomyces Genome Database.
           : It is possible that it will be moved to SGD::WWW.pm in the future.

See Also   : L<sgd_url>(), L<s3d_link>()

=cut

#--------------
sub sgd_link { 
#--------------
    my($self,$arg,$value,$text) = @_; 
    my $url = $self->sgd_url($arg);
    $text ||= $value;
    qq|<A HREF="$url$value">$text</A>|;
}


#########################################################################
##                        INSTANCE METHODS                              
#########################################################################

## Note that similar functions to those presented below are also availble 
## via L. Stein's CGI.pm. These are more experimental versions.

=head2 start_html

 Usage     : $BioWWW->start_html()
 Purpose   : Prints the "Content-type: text/html\n\n<HTML>\n" header.
 Returns   : n/a; This method prints the Content-type string shown above.
 Argument  : n/a
 Throws    : n/a
 Status    : Experimental
 Comments  : This method prevents redundant invocations thus avoiding th
           : accidental printing of the "content-type..." on the page.
           : If using L. Stein's CGI.pm, this is similar to $query->header()
           : (Does CGI.pm prevent redundant invocation?)

=cut

#---------------'
sub start_html { 
#---------------
    my $self=shift; 
    if(!$self->{'_started_html'}) {
	print "Content-type: text/html\n\n<HTML>\n";
	$self->{'_started_html'} = 1;
    }
}


=head2 redirect

 Usage     : $BioWWW->redirect(<string>)
 Purpose   : Prints the header needed to redirect a web browser to a supplied URL. 
 Returns   : n/a; Prints the redirection header.
 Argument  : String containing the URL to be redirected to.
 Throws    : n/a
 Status    : Experimental

=cut

#-------------
sub redirect {
#-------------
    my($self,$url) = @_;

    print "Location: $url\n";
    print "Content-type: text/html\n\n"; 
}



=head2 pre

 Usage     : $BioWWW->pre("text to be pre-formatted");
 Purpose   : To produce HTML for text that is not to be formated by the brower.
 Returns   : String containing the "<pre>" formatted html.
 Argument  : n/a
 Throws    : n/a
 Status    : Experimental

=cut

#--------
sub pre { 
#--------
    my $self = shift; 
    "<PRE>\n".shift()."\n</PRE>";
}


#----------------
sub html_footer {
#----------------
    my( $self, @param ) = @_;  

    my( $linkTo, $linkText, $modified, $mail, $mailText, $top) = 
	$self->_rearrange([qw(LINKTO LINKTEXT MODIFIED MAIL MAILTEXT TOP)], @param);

    $modified = (scalar $modified) 
	? qq|<center><small><b>Last modified: $modified </b></small></center>| 
        : '';

    $linkTo ||= '';

#    $top = (defined $top) ? qq|<a href="top">Top</a><br>| : '';
    $top = qq|<a href="#top">Top</a>|;  ## Utilizing the HTML bug/feature wherein 
                                          ## a bogus name anchor defaults to the 
                                          ## top of the page. 

    return <<"HTML";	
<p>
<hr size=3 noshade width=95%>
$top | <a href="$linkTo"> $linkText</a><br>
$modified
<small><i><a href="mailto:$mail">$mailText</a></i></small>
</body></html>

HTML
}


=head2 strip_html

 Usage     : $boolean = &strip_html( string_ref, [fast] );
 Purpose   : Removes HTML formatting from a supplied string.
 Returns   : Boolean: true if string was stripped, false if not.
 Argument  : string_ref = reference to a string containing the whole 
           :              web page to be stripped.
           : fast = a non-zero value. Optional. If set, a faster 
           :        but perhaps less thorough procedure is used for
           :        stripping. Default = not fast.
 Throws    : Exception if the argument is not a scalar reference.
 Comments  : Based on code originally written by Alex Dong Li
           : (ali@genet.sickkids.on.ca).
           : This is a more generic version of the function that appears 
           : in Bio::Tools::Blast::HTML.pm
           : This version does not perform any Blast-specific stripping.
           :
           : This employs a simple method for removing tags that
           : will fail under following conditions:
           :  1) if quoted > appears in a tag  (does this ever happen?)
           :  2) if a tag is split over multiple lines and this method is
           :     used to process one line at a time.
           :
           : Without fast mode, large HTML files can take exceedingly long times to
           : strip (e.g., 1Meg file with many tags can take 10 minutes versus 5 seconds
           : in fast mode. Try the swissprot yeast table). If you know the HTML to be
           : well-behaved (i.e., tags are not split across mutiple lines), use fast
           : mode for large, dense files.

=cut

#---------------
sub strip_html {
#---------------
    my ($self, $string_ref, $fast) = @_;

    ref $string_ref eq 'SCALAR' or 
	$self->throw("Can't strip HTML: ".
		     "Argument is should be a SCALAR reference not a ${\ref $string_ref}");

    my $str = $$string_ref;
    my $stripped = 0;

    if($fast) {
	# MULTI-STRING-MODE: Much faster than single-string mode
	# but will miss tags that span multiple lines.
	# This is fine if you know the HTML to be "well-behaved".
	
	my @lines = split("\n", $str);
	foreach (@lines) {
	    s/<[^>]+>|&nbsp//gi and $stripped = 1;
	}
	
	# This regexp likely won't work properly in this mode.
	foreach (@lines) {
	    s/(\A|\n)>\s+/\n\n>/gi and $stripped = 1;
	}
	$$string_ref = join ("\n", @lines);

    } else {

	# SINGLE-STRING-MODE: Can be very slow for long strings with many substitutions.
	
	# Removing all "<>" tags. 
	$str =~ s/<[^>]+>|&nbsp//sgi and $stripped = 1;
	
	# Re-uniting any lone '>' characters. Not really necessary for functional HTML
	$str =~ s/(\A|\n)>\s+/\n\n>/sgi and $stripped = 1;
	
	$$string_ref = $str;
    }
    $stripped;
}


1;
__END__

########################################################################
##                            END OF CLASS                             
########################################################################
 
=head1 FOR DEVELOPERS ONLY

=head2 Data Members

An instance of Bio::Tools::WWW.pm is a blessed reference to a hash containing
all or some of the following fields:

 FIELD           VALUE
 --------------------------------------------------------------
 _started_html   Defined the on the initial invocation of start_html()
                 to avoid duplicate printing out the "Content-type..." header.


=cut

1;


