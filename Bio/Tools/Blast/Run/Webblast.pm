#--------------------------------------------------------------------
# PACKAGE : Bio::Tools::Blast::Run::Webblast.pm
# PURPOSE : To run a Blast analysis on a remote server 
#           and save the results locally.
# AUTHOR  : Steve Chervitz (sac@bioperl.org) 
#               - Webblast.pm: modularized version of webblast.
#           Alex Dong Li (ali@genet.sickkids.on.ca) 
#               - original webblast script.
#           Ross N. Crowhurst (RCrowhurst@hort.cri.nz)
#               - modified Webblast.pm to use LWP to give proxy server 
#                 support.
# CREATED : 4 May 1998
# STATUS  : Alpha
# REVISION: $Id$
#
# For documentation, run this module through pod2html
# (preferably from Perl v5.004 or better).
#
# MODIFICATION HISTORY: See bottom of file.
#
# Copyright (c) 1998 Alex Dong Li/Steve Chervitz/Ross N. Crowhurst. 
#           All Rights Reserved.
#           This module is free software; you can redistribute it and/or 
#           modify it under the same terms as Perl itself.
#------------------------------------------------------------------------

package Bio::Tools::Blast::Run::Webblast;
use strict;

# rnc: now uses HTTP and LWP.
# sac: Make sure HTTP and LWP are available.
#      This check is also performed during bioperl installation.
#      Warning instead of die-ing since exported variables can
#      still be examined. Dieing only if blast_remote() is attempted.

BEGIN {
  use vars qw($Loaded_LWP $Loaded_IOScalar);
  $Loaded_LWP = 1;
  unless( eval "require HTTP::Request::Common" and
	  eval "require LWP::UserAgent") {
    warn "\a\n".'='x50, "\n".
      "WARNING: COULDN'T LOAD THE LWP MODULE.\n\n".
      "   LWP (libwww-perl) is now required to run remote Blasts.\n".
      "   Download it from CPAN: http://www.perl.com/CPAN/.".
	"\n".'='x50, "\n\n";
    $Loaded_LWP = 0;
  }

  HTTP::Request::Common->import(qw(POST)) if $Loaded_LWP;

  $Loaded_IOScalar = 1;
  unless( eval "require IO::Scalar") {
    warn "\a\n".'='x50, "\n".
      "WARNING: COULDN'T LOAD THE IO::Scalar MODULE.\n\n".
      "   IO::Scalar is now required to run remote Blasts.\n".
      "   This module is included in the IO-stringy collection\n".
      "   from CPAN: http://www.perl.com/CPAN/.".
	"\n".'='x50, "\n\n";
    $Loaded_IOScalar = 0;
  }
}

use Bio::SeqIO;
use Bio::Root::Global    qw(:devel);
use Bio::Root::Utilities qw(:obj);
use Bio::Tools::Blast::HTML qw(&strip_html);
use IO::Scalar;
use Carp;

use Exporter;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
@ISA        = qw(Exporter);
@EXPORT     = qw();
@EXPORT_OK  = qw(&blast_remote @Blast_dbp_remote @Blast_dbn_remote 
		 @Blast_matrix_remote);
%EXPORT_TAGS = ( std => [qw(&blast_remote @Blast_dbp_remote @Blast_dbn_remote 
			    @Blast_matrix_remote)] );

use strict;
use vars qw( $ID $VERSION $revision);

$ID       = 'Bio::Tools::Blast::Run::Webblast';
$VERSION  = 1.24; 

# SAC: grouped database names.
# May want to put these sorts of things in a
# Blast::Resources.pm type of module.
use vars qw(@Blast_dbp_remote @Blast_dbn_remote @Blast_matrix_remote);

# Peptide datasets
@Blast_dbp_remote = qw(nr month swissprot yeast ecoli kabat pdb alu );

# Nucleotide datasets
@Blast_dbn_remote = qw(nr month swissprot yeast ecoli kabat pdb alu 
		       dbest dbsts htgs vector mito epd gss
		       est_mouse est_human est_others );

# Not necessarily an exhaustive list.
# rnc:	added PAM30, PAM70, BLOSUM45, BLOSUM80 from NCBI BLAST2 
@Blast_matrix_remote  = qw(BLOSUM45 BLOSUM62 BLOSUM80 PAM30 PAM40 PAM70 PAM120 PAM250 IDENTITY);

# SAC: Consolidated the various default values.
# These are the defaults specified by the NCBI server.
my %_default = (
		'expect'                => 10,  
		'expectPSI'             => 10,	# rnc: PSI expect default is 10 not 0.01,  
		'descriptions'          => 500,	# rnc: DESCRIPTION default is 500 not 100,
		'alignments'            => 500,	# rnc: ALIGNMENTS default is 500 not 100,
		'queryGeneticCode'      => 1,
		'advancedOptionG'       => 11,  # gap creation
		'advancedOptionE'       => 1,   # gap extension
		'advancedOptionQ'       => -3,  # mismatch penalty,  blastn only
		'advancedOptionR'       => 1,   # match reward,      blastn only
		'advancedOptionW'       => 11,  # word size (11 blastn, 3 others)
		'resultReturnProtocol'  => 'www',   # or 'email'
		'minimumLength'         => 30,
		'minimumLengthProtein'  => 10,
		'maximumLength'         => 100000, 
		'maximumLengthProtein'  => 50000, 
# Blast1 no longer supported at NCBI
#		'blastServerURL1'       => 'www.ncbi.nlm.nih.gov/',
		'blastServerURL2'       => 'www.ncbi.nlm.nih.gov/blast/blast.cgi',
		'blastServerURLpsi'     => 'www.ncbi.nlm.nih.gov/blast/psiblast.cgi',
		'blastServerURLphi'     => 'www.ncbi.nlm.nih.gov/blast/phiblast.cgi',
		# rnc:	added WashU-Blast2 server address, ..WashU - first remote WashU Blast Server to use
		# rnc:	added proxy server address option
		'blastServerURLWashU'	=> 'www2.ebi.ac.uk/cgi-bin/newblast2.pl',
		'entryServerURLroot'    => 'www.ncbi.nlm.nih.gov',
		'proxyServerURL'        => '',
		'imageMapURL'           => 'www.ncbi.nlm.nih.gov/cgi-bin/imagemap/BLAST/',

                # rnc: added gi_list , list_org, alignment_view, cutoff, input_type,
		#      strand, histogram,
                'gi_list'               => '(None)',  
                'list_org'              => '',
                'alignment_view'        => '0',
                'cutoff'                => 'default',
                'input_type'            => 'Sequence in FASTA format',
                'strand'                => 'Both',
                'histogram'             => 'Both',

                # sequenceFormat no longer an issue since we now work with seq objects.
		#'sequenceFormat'        => 'GCG', 
                # notAsk is no longer an option since we are trying to automate this procedure
                # at a higher level.
		#'notAsk'                => 'on', 
		);

my $program                 ='blastp';
my $version                 = 2;  # default version of Blast to run
# rnc:  no longer need 'GAPPED_ALIGNMENT=is_set' just 'is_set' #my $gappedAlignmentOption   ='GAPPED_ALIGNMENT=is_set';
my $gappedAlignmentOption   ='is_set';
my $gappedAlignmentFlag     ='on';
my $database                ='nr';
my $extension               = $program;
my $expect                  = $_default{'expect'};      
my $expectPSI               = $_default{'expectPSI'};      
my $filter                  ='default';             # protein-protein searches only
my $matrix                  ='BLOSUM62';            # scoring matrix 
my $NCBI_giOption           ='';                    # default is off
my $NCBI_giFlag             ='off';                 # default is off
my $descriptions            = $_default{'descriptions'};
my $alignments              = $_default{'alignments'};

my $gi_list                 = $_default{'gi_list'}; # default is '(None)'   
my $list_org                = $_default{'list_org'}; #default is ''
my $alignment_view          = $_default{'alignment_view'};
    
# rnc:  no longer need 'OVERVIEW=is_set' just 'is_set'
#       changed default of $graphicalOverviewOption to 'is_set' #my $graphicalOverviewOption = 'OVERVIEW=is_set'; # or '';     # not available for email
#       NOT SURE HOW setting graphicalOverviewOption to on will affect subsequent parsing but its a useful
#       feature for the actual blast report if you load that report into netscape
#       Maybe BioPerl could have an option added to indicate users desire to send HTML version of
#       blast reports to netscape when there is only going to be one or 2 reports and they are
#       running X - this can be controlled elsewhere however but would increase user-friendliness
my $graphicalOverviewOption = ''; # or 'is_set';     # not available for email so is set to '' below when email option on

my $cutoff              = $_default{'cutoff'};
my $input_type          = $_default{'input_type'};
my $strand              = $_default{'strand'};
my $histogram           = $_default{'histogram'};

my $graphicalOverviewFlag   ='off'; # or 'on';
my $queryGeneticCode        = $_default{'queryGeneticCode'};

## All advanced options are left undefined unless specified by user.
my $advancedOptionG         = ''; 
my $advancedOptionE         = ''; 
my $advancedOptionQ         = ''; 
my $advancedOptionR         = ''; 
my $advancedOptionW         = ''; 

my $blastServerURL          = '';
my $resultReturnProtocol    = $_default{'resultReturnProtocol'};
my $resultTypeOption        = '';       # null for plain text or html=is_set
my $resultTypeInHTMLflag    ='on';
#$imageMapURL               =  $_default{'imageMapURL'};
#$imageURL=$imageURLdefault ='www.ncbi.nlm.nih.gov/BLAST/';
my $minimumLength           = $_default{'minimumLength'};
my $maximumLength           = $_default{'maximumLength'};
# this is the root of server for retrieving record, might be different from blastServer
my $entryServerURLroot      = $_default{'entryServerURLroot'};

# rnc: PROXY SERVER 
#       enter valid proxy server name:port eg my $proxyServerURL = "http://proxy.akl.ihug.co.nz:8080";
#       or leave as default (empty string).
# rnc:  could use $pserver_name and $pserver_port then $proxyServerURL = $pserver_name . ":" . $pserver_port
#       if better to separate
# rnc:  my first ISP proxyserver
#my $proxyServerURL = "http://proxy.akl.ihug.co.nz:8080"; # if no proxyserver required then use my $proxyServerURL = "";  
# rnc:  my work proxyserver
#my $proxyServerURL = "http://proxy.marc.cri.nz:8080"; # if no proxyserver required then use my $proxyServerURL = "";  

my $proxyServerURL          = $_default{'proxyserver'};

my $emailAddress            = '';
my $emailOption             = '';
my $outputFileNamePost      = "posting_parameters";
my $currentDate             = '';
my @queue                   = ();  # list of sequence objects

# SAC: new private, class vars
my @_outFiles               = ();
my $_out_dir                = '';    # where to save the Blast output files.
my $_blastObj               = undef;
my $_errmsg                 = '';


# SAC: defunct options and variables:
#my $notAsk             =  $_default{'notAsk'};
#my $sequenceFormat     = $_default{'sequenceFormat'};
#my $queue='';          #A string of tab delimited clone names.
#my $currentTime = '';
#my $username = getpwuid($<);            # get userid
#my @username   = split(//, $username);  # convert the first case into upper
#$username[0] =~ tr/a-z/A-Z/;
#$username    = join("", @username);



## POD Documentation

=head1 NAME

Bio::Tools::Blast::Run::Webblast.pm - Bioperl module for running Blast analyses using a HTTP interface.

=head1 SYNOPSIS

    # Run a Blast
    use Bio::Tools::Blast::Run::Webblast qw(&blast_remote);

    @out_file_names = &blast_remote($object, %named_parameters);

L<blast_remote> is the only exported method of this module
and it returns a list of local file names containing the Blast
reports. C<$object> is a reference to a B<Bio::Root::Object.pm> object or
subclass. See L<blast_remote>() for a description of available parameters.

    # Obtain a list of available databases

    use Bio::Tools::Blast::Run::Webblast qw(@Blast_dbp_remote
					    @Blast_dbn_remote);

    @amino_dbs      = @Blast_dbp_remote;
    @nucleotide_dbs = @Blast_dbn_remote;


=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

=head1 DESCRIPTION

B<Bio::Tools::Blast::Run::Webblast.pm> contains methods and data necessary for
running Blast sequence analyses using a remote server and saving the results locally.

B<Bio::Tools::Blast::run()> provides an interface for Webblast.pm,
so, ideally, you shouldn't use Webblast.pm directly, but via Blast.pm. 

B<FEATURES:>

=over 2

=item * Supports NCBI Blast1, Blast2, and PSI-Blast servers as well as WashU-Blast servers.

=item * Can operate through a proxy server enabling operation from behind a firewall.

=item * Can save reports with and without HTML formatting.

=item * Uses LWP.

=back

In principle, this module can be customized to use different servers
that provide a Blast interface like the NCBI or WashU style servers. 
Such servers could be remote or local. This hasn't been well-tested however. 


=head1 DEPENDENCIES

Bio::Tools::Blast::Run::Webblast.pm is used by B<Bio::Tools::Blast.pm>.
The development of this is thus linked with the Blast.pm module.

=head1 SEE ALSO

 Bio::Tools::Blast.pm                    - Blast object.
 Bio::Tools::Blast::Run::LocalBlast.pm   - Utility module for running Blasts locally.
 Bio::Tools::Blast::HTML.pm              - Blast HTML-formating utility class.
 Bio::Seq.pm                             - Biosequence object  
 Bio::Root::Object.pm                    - Bioperl base object class.


 http://bio.perl.org/Projects/modules.html  - Online module documentation
 http://bio.perl.org/Projects/Blast/        - Bioperl Blast Project     
 http://bio.perl.org/                       - Bioperl Project Homepage

=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bio.perl.org/MailList.html  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via email
or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR

=over 0

=item Steve Chervitz E<lt>sac@bioperl.orgE<gt>

    - Webblast.pm modularized version of webblast script.

=item Alex Dong Li E<lt>ali@genet.sickkids.on.caE<gt>

    - original webblast script.

=item Ross N. Crowhurst E<lt>RCrowhurst@hort.cri.nzE<gt>

    - modified Webblast.pm to use LWP to give proxy server support.

=back 

=head1 VERSION

Bio::Tools::Blast::Run::Webblast.pm, 1.24

=head1 COPYRIGHT

Copyright (c) 1998, 1999 Steve Chervitz, Alex Dong Li,  Ross N. Crowhurst. 
All Rights Reserved.This module is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

=cut

#
##
###
#### END of main POD documentation.
###
##
#'

######################  BEGIN FUNCTIONS  ########################

=head1 APPENDIX

Methods beginning with a leading underscore are considered private
and are intended for internal use by this module. They are
B<not> considered part of the public interface and are described here
for documentation purposes only.


=head2 blast_remote

 Usage     : @files = blast_remote( $blast_object,  %namedParameters);
           : This method is exported.
 Purpose   : Run a remote Blast analysis on one or more sequences.
           : NOTE: The name of this method is potentially misleading
           :       since the a local server could be specified.
  	   :       Probably should be called blast_http.
 Returns   : Array containing a list of filenames of the Blast reports.
 Argument  : First argument should be a Bio::Tools::Blast.pm object reference.
           : This object is primarily used for error reporting
           : Remaining arguments are named parameters: 
           : (PARAMETER TAGS CAN BE UPPER OR LOWER CASE).
           :
           :   -ALIGN      => integer, number of alignments (B, 100)
           :   -ALIGN_VIEW => alignment view option (see below)
           :   -CUTOFF     => Blast score cutoff (60-110 or 'default')
           :   -DATABASE   => name of database (see below)
           :   -DESCR      => integer, number of on-line descriptions (V, 100)
           :   -EXPECT     => expect value cutoff
           :   -EXPECT_PSI => expect value for inclusion in PSI-BLAST iteration 1 
           :   -FILTER     => sequence complexity filter ('default' or 'none')
           :   -GAP        => 'on' or 'off'
           :   -GAP_CREATE => gap creation penalty (G, 5) 
           :   -GAP_EXTEND => gap extension penalty (E, 2)
           :   -GEN_CODE   => integer for special genetic code (see below) blastx only
           :   -GRAPH      => 'on' or 'off' (graphical overview not yet supported)
           :   -HISTOGRAM  => 'on' or 'off' or 'both'
           :   -HTML       => 'on' or 'off' or 'both'
           :   -INPUT_TYPE => 'Sequence in FASTA format' or 'Accession or GI'
           :   -MATRIX     => substitution scoring matrix (blast1 only for NCBI server)
           :   -NCBI_GI    => 'on' or 'off'
           :   -MATCH      => match reward (r, 1)       (blastn only)
           :   -MAX_LEN    => max query sequence length to blast
           :   -MIN_LEN    => min query sequence length to blast
           :   -MISMATCH   => mismatch penalty (q, -3)  (blastn only)
           :   -ORGANISM   => organism name to limit Blast2 search.
           :   -ORGANISM_CUSTOM  => custom organism or taxon name.
           :   -OUT_DIR    => output directory to store blast result files
           :   -PROG       => name of blast program (blastp, blastx, etc.)
           :   -SEQS       => ref to an array of Bio::Seq.pm objects. 
           :   -SERVER     => blast server to use (default is NCBI Blast2)
           :   -STRAND     => Default = 'Both' (not used by NCBI servers)
           :   -VERSION    => blast version (1, 2, PSI, WashU)
           :   -WORD       => word size (W, 11 for blastn, 3 for all others)
#rnc:   LIST_ORG
#       valid list_org entries for blast2 are a string of 50 chars max, default is empty string
           :
 Throws    : Exception if:
           :   - Cannot obtain parameters by calling _rearrange() on the
           :     first argument, which should be a Bio::Tools::Blast.pm object ref.
           :   - No sequences are provided.
           :   - Sequence type is incompatible with Blast program type.
           :   - Database name is not one of the valid names.
           :   - Supplied e-mail address looks invalid.
 Comments  :
  -------------------------------------------------------------
  Available programs: blastn, blastx, dbest, blastp, tblastn, tblastx
  Program versions: 1, 2, PSI, WashU (or WU)

  -------------------------------------------------------------
  Available databases:
        nr, month, swissprot, dbest, dbsts, 
        est_mouse, est_human, est_others, pdb, vector, kabat,
        mito, alu, epd, yeast, ecoli, gss, htgs.

    These are exported by this module in the @Blast_dbp_remote
    and @Blast_dbn_remote arrays.
  -------------------------------------------------------------
  Available Genetic Codes are (blastx only): 

        (1) Standard                    (2) Vertebrate Mitochondrial
        (3) Yeast Mitochondrial         (4) Mold Mitochondrial; ... 
        (5) Invertebrate Mitochondrial  (6) Ciliate Nuclear; ...
        (9) Echinoderm Mitochondrial    (10) Euplotid Nuclear
        (11) Bacterial                  (12) Alternative Yeast Nuclear
        (13) Ascidian Mitochondrial     (14) Flatworm Mitochondrial
        (15) Blepharisma Macronuclear

  -------------------------------------------------------------
  Available values for organism (Blast2):

      (None)   (DEFAULT; note that the parentheses are required.)
      Arabidopsis thaliana 
      Bacillus subtilis 
      Bos taurus 
      Caenorhabditis elegans 
      Danio rerio 
      Dictyostelium discoideum 
      Drosophila melanogaster 
      Escherichia coli 
      Gallus gallus 
      Homo sapiens 
      Human immunodeficiency virus type 1 
      Mus musculus 
      Oryctolagus cuniculus 
      Oryza sativa 
      Ovis aries 
      Plasmodium falciparum 
      Rattus norvegicus 
      Saccharomyces cerevisiae 
      Schizosaccharomyces pombe 
      Simian immunodeficiency virus 
      Xenopus laevis 
      Zea mays 

  -------------------------------------------------------------
  Available values for align_view (Blast2):

       0             Pairwise  (DEFAULT)
       1             master-slave with identities
       2             master-slave without identities
       3             flat master-slave with identities
       4             flat master-slave without identities
 -------------------------------------------------------------
  Available substitution scoring matrices (NCBI):

      BLAST2 matrices: BLOSUM80, BLOSUM62, BLOSUM45, PAM30, PAM70

      BLAST1 matrices: BLOSUM62, PAM40, PAM120, PAM250, IDENTITY.

      Others members of the BLOSUM and PAM family of matrices
      may be available as well.
      These are exported by this module in the @Blast_matrix_remote array.

      Note that certain combinations of matrices and gap creation/extension
      penalties are disallowed (E.g., PAM250 will work with 12/2 but not 11/1).
 --------------------------------------------------------------
   Limited values for gap creation and extension are supported for 
   blastp, blastx, tblastn.  Some supported and suggested values are:

  Creation     Extension

     10             1
     10             2
     11             1
      8             2
      9             2
  -------------------------------------------------------------
  Available sequence complexity filters:
       SEG, SEG+XNU, XNU, dust, none.


See Also : _set_options(), _adjust_options(), _validate_options(), _blast(), B<Bio::Tools::Blast.pm>

=cut

#-----------------
sub blast_remote {
#-----------------
    my ($bobj, %param) = @_;
    
    unless ($Loaded_LWP) {
      croak ("THE LWP MODULE IS NOT INSTALLED.\n\n".
	     "   LWP (libwww-perl) is now required to run remote Blasts.\n".
	     "   Download it from CPAN: http://www.perl.com/CPAN/.\n");
    }

    unless ($Loaded_IOScalar) {
      croak ("THE IO::Scalar MODULE IS NOT INSTALLED.\n\n".
             "   IO::Scalar is now required to run remote Blasts.\n".
             "   This module is included in the IO-stringy collection\n".
             "   from CPAN: http://www.perl.com/CPAN/.");
  }

    my ($seq_a);
    eval { 
	# _rearrange() is an instance method of Bio::Root::Object.pm and is
	# inherited by the Blast object.
	($seq_a) = $bobj->_rearrange([qw(SEQS)], %param);
    };
    if($@) {
	# Ideally we should use can(), requires Perl 5.004.
	croak "Can't run remote BLAST: failed to obtain parameters using _rearrange().\n".
	    "Make sure that the first argument is a Bio::Tools::Blast.pm object.\n";
    }

    ref($seq_a) eq 'ARRAY' or
	$bobj->throw("Can't get sequences to be Blasted: No sequence data.");

    if( not @queue = @$seq_a) {
	$bobj->throw("Can't Blast: No sequences supplied.");
    }

    @_outFiles = ();

    $_blastObj = $bobj;
    &_set_options(%param);
    &_adjust_options();
    &_validate_options();
    &_blast();

    return @_outFiles;
}

#----------------
sub _set_options {
#----------------
# SAC: Processing function args instead of command-line args
#      Most have the same names as Alex's.
    my (@param) = @_;
    my ($prog, $vers, $db, $exp, $exp_psi, $filt, $mat, $gi, $descr, $aln, 
	$graph, $gap, $gcode, $email, $html, $server, $maxl, $minl, 
	$gap_c, $gap_e, $mmatch, $match, $word, $dir,
	$org, $org_custom, $aln_view, $cut, $in_type, $strnd, $hist, $proxy
       )
	= $_blastObj->_rearrange(
			 [qw(PROG VERSION DATABASE EXPECT EXPECT_PSI FILTER MATRIX 
			  NCBI_GI DESCR ALIGN GRAPH GAP GEN_CODE EMAIL HTML SERVER 
			  MAX_LEN MIN_LEN GAP_CREATE GAP_EXTEND MISMATCH MATCH 
			  WORD OUT_DIR
			  ORGANISM ORGANISM_CUSTOM ALIGN_VIEW CUTOFF 
		          INPUT_TYPE STRAND HISTOGRAM  PROXY_SERVER
			    )], @param);
		      
    if($prog) {
	$program = $prog;
	# SAC: commented out advanced option setting. Current strategy
	#      is to set them as needed by the user. Otherwise, let server
	#      pick the defaults for us.
	if ($prog =~ m/^blastx$/i) {
#	    $advancedOptionG = $_default{'advancedOptionG'} = 10;
#	    $advancedOptionE = $_default{'advancedOptionE'} = 1;
	}
	# reset default if submit protein seq.
	if ($prog =~ m/^blastp$/i || $prog=~m/^tblastn$/i) {
#	    $advancedOptionG = $_default{'advancedOptionG'} = 10;
#	    $advancedOptionE = $_default{'advancedOptionE'} = 1;
	    $minimumLength   = $_default{'minimumLengthProtein'};
	    $maximumLength   = $_default{'maximumLengthProtein'};   
	}
    }
    if($vers)    { $version = $vers;  }
    if($gap)     { $gappedAlignmentFlag = $gap;  }
    if($db)      { $database = lc($db);  }
    if($exp)     { $expect = $exp;  }
    if($exp_psi) { $expectPSI = $exp_psi;  }
    if($filt)    { $filter = $filt;   }
    if($mat)     { $matrix = $mat;   }
    if($gi)      { $NCBI_giFlag = $gi;   }
    if($descr)   { $descriptions = $descr;   }
    if($aln)     { $alignments = $aln;   }
    if($graph)   { $graphicalOverviewFlag = $graph;   }
    if($gcode)   { $queryGeneticCode = $gcode; }
    if($gap_c)   { $advancedOptionG = $gap_c; }
    if($gap_e)   { $advancedOptionE = $gap_e; }
    if($mmatch)  { $advancedOptionQ = $mmatch; }
    if($match)   { $advancedOptionR = $match; }
    if($word)    { $advancedOptionW = $word; }
    if($email)   { $resultReturnProtocol="email";
		   $emailAddress = $email; }
    if($server)  { $blastServerURL = $server;}
           else  { $blastServerURL = '';}
    if($minl)    { $minimumLength = $minl; } # if match minlength
    if($maxl)    { $maximumLength = $maxl; } # if match maxlength
    if($dir)     { $_out_dir = $dir; }

    # sac: New options to support rnc's additions.
    if($org)         { $gi_list = $org; }
    if($org_custom)  { $list_org = $org_custom; }
    if($aln_view)    { $alignment_view = $aln_view; }
    if($cut)         { $cutoff = $cut; }
    if($in_type) { 
  	if($in_type =~ /fasta/i) {
  	  $input_type = 'Sequence in FASTA format'; 
  	} elsif($in_type =~ /gi|acc/i) {
  	  # This is not yet supported.
  	  # $input_type = 'Accession or GI'; 
  	  croak "\nUnsupported sequence input type for Blast: $in_type\n".
  	    "Input type must be 'Fasta'\n";
  	} else {
  	  croak "\nUnsupported sequence input type for Blast: $in_type\n".
  	    "Input type must be 'Fasta' or 'GI' or 'accession'\n";
  	}
    }
    if($strnd)       { $strand = $strnd; }
    if($hist)        { $histogram = $hist; }
    if($proxy)       { $proxyServerURL = $proxy;}

    if($html) {
	# HTML issue is a bit tricky. Blast.pm requires a non-HTML version
	# for parsing. So if the report is to be parsed, we need to make sure
	# to give Blast.pm the correct version of the results.
	# If HTML is 'on', then only the HTML version is saved and we
	# can't parse. (The Blast object may still try to parse it if its mode is
	# set for parsing, but it will fail). 
	if($html !~ /on|both/i) {
	    $html = 'on';
	}
	$resultTypeInHTMLflag = $html; 
    } else { 
	$resultTypeInHTMLflag = 'off'; 
    }
}


#--------------------
sub _adjust_options {
#--------------------
# reset some varibles depending on the flags
# SAC: fixed bug with gappedAlignmentFlag
#      default is gapped. Need to specify ungapped.
    if ($gappedAlignmentFlag=~m/^on$/i) {
	$gappedAlignmentOption="";
    } else {
	# rnc: 	no passing args in $options to postclient.pl any more 
	# 	so only need 'is_set' #$gappedAlignmentOption="UNGAPPED_ALIGNMENT=is_set";
	$gappedAlignmentOption="is_set";
    }
    
    if ($NCBI_giFlag=~m/^on$/i) {
	# rnc: 	no passing args in $options to postclient.pl any more 
	# 	so only need 'is_set'	#$NCBI_giOption="NCBI_GI=is_set";
	$NCBI_giOption="is_set";
    } else {
	$NCBI_giOption="";
    }
    
    if ($graphicalOverviewFlag=~m/^on$/i) {
	# rnc: 	no passing args in $options to postclient.pl any more 
	# 	so only need 'is_set'	#$graphicalOverviewOption="OVERVIEW=is_set";
	$graphicalOverviewOption="is_set";
    } else {
	$graphicalOverviewOption="";
    }
    
    # off means plain text
    if ($resultTypeInHTMLflag=~m/^off$/i) {
	$resultTypeOption="";         # or html
# with html tags
    } else {
	# rnc: 	no passing args in $options to postclient.pl any more 
	# 	so only need 'is_set'	#$resultTypeOption="HTML=is_set";         # or html
	$resultTypeOption="is_set";         # or html
    }
    
    
    if ($resultReturnProtocol ne "www") {
	#check if email is in right format here
	if ($emailAddress !~m/\w+@\w+\.\w+/) {
	    $_blastObj->throw("Email address may be wrong. Please check it and try again.");
	}
    } else {
	$currentDate = $Util->date_format('full');
    }
    
    # SAC: select which Blast server to use depending on $version
    # (if a server has not already been selected).
    if(not $blastServerURL) {
	local($_) = $version;
	SWITCH: {
	    /1/i   && do{ $blastServerURL = $_default{'blastServerURL1'}; last SWITCH;};
	    /2/i   && do{ $blastServerURL = $_default{'blastServerURL2'}; last SWITCH;};
	    /psi/i && do{ $blastServerURL = $_default{'blastServerURLpsi'}; last SWITCH;};
	    /wu|wash/i  && do{ $blastServerURL = $_default{'blastServerURLWashU'}; last SWITCH;}; 
	}
    } else {
	$version = 'custom_server';
    }

    $blastServerURL = 'http://' . $blastServerURL unless $blastServerURL =~ /^http:/i;
}


#----------------------
sub _validate_options {
#----------------------
# some error checks.
# SAC: Throwing exceptions if
#       - can't access postclient.pl
#       - invalid program type
#       - invalid database name

# rnc:  commented out the code for postclient.pl - it is not used anymore
#
# test to see if the program is executable
#    if (!-x "$pathOfPostclient") {
#	## SAC: Added extra info to the die call.
#	$_blastObj->throw("Can't access postclient.pl: $pathOfPostclient.");
#    }
    
# if graphical is on, I need getclient to retrieve the .gif file!
#if ( $graphicalOverviewFlag eq "off" && !-x "/home/ali/script/getclient.pl") {
#    die("\n\nI cannot access getclient.pl, which is needed for graphs. I quit!\a\n\n");
#}

# if program is not of one available program, quit
    # SAC: simplified this conditional
    if (!($program=~m/^t?blast[npx]$/i or $program=~m/^dbest$/i)) {
	$_blastObj->throw("Invalid program name: $program",
			  "It must be one of blastn, blastx, blastp, tblastn, tblastx, or dbest!");
    }
    
# if database is not of one available databases, quit
    # SAC: simplified this conditional and error message
    # SAC: Only checking against nucleotide db's since they are a superset.
    if ( not(grep $database eq $_, @Blast_dbn_remote)) {
	my $valid = join(', ', @Blast_dbn_remote);
	$_blastObj->throw("Invalid Blast database name: $database",
			  "It must be one of $valid");
    }

# if $program is not dbest, use program as output file extension
    if ($program !~m/^dbest$/i) {
	$extension  = $program;  # SAC: note that program case is not changed.
# otherwise use blastn as program, dbest as database name and output file extension
    } else {
	$program="blastn";
	$database="dbest";
	$extension="dbest";
    }
}

#------------
sub _blast {
#------------
# work on all sequences one by one
# SAC: Some alterations in the basic plan of this function.
	my ($seq, $typeDo, $advancedOptions);
	my ($id, $outputFileNameTemp, $outputFileName, $baseFileName, $options);
	my $countSuccessfulBlasting = 0;
	# rnc:	added following
	my $outFile = '';
	my($ua, $sequenceInFastaFormat, $req, $response); 
	my $count = 0;

    	foreach $seq (@queue)  {
	      
	    $count++;
	    # Not catching fatal error if the sequence is the wrong type
	    # (serious error).
	    next unless &_validate_seq($seq);
	    
	    # add >seq_name to make seq in fasta format
            my $sh = new_tie IO::Scalar \$sequenceInFastaFormat;
            my $out = Bio::SeqIO->new ('-fh'   => $sh,   
                                       '-format' =>'Fasta');

            $out->write_seq($seq);
	    $sequenceInFastaFormat =~ s/\n$//;
	    
	    # concatenate advanced options into a string, depending on differnt blast type
	    # options -q and -r must be lower case!
	    # SAC: setting options only if they are defined.
	    $advancedOptions  = $advancedOptionG ?  "-G $advancedOptionG" : '';
	    $advancedOptions .= $advancedOptionE ? " -E $advancedOptionE" : '';
	    $advancedOptions .= $advancedOptionW ? " -W $advancedOptionW" : '';
	    if ($program=~m/^blastn$/i) 
	      {
		  $advancedOptions .= $advancedOptionQ ? " -q $advancedOptionQ" : '';
		  $advancedOptions .= $advancedOptionR ? " -r $advancedOptionR" : '';
	      }
	    
	    # get temp filename and final filename
	    $id = $seq->id;
	    $id =~ s/[|:\/\\<>]/_/g;  # remove special chars that could confuse shell
	    $id =~ s/\.\w+?$//;  # trim off the extension if id is a filename.
	    $baseFileName = $_out_dir. $id . "." . $extension . "$version.$database"; # SAC: new var.
	    $outputFileNameTemp = $baseFileName. ".temp.html";
	    $outputFileName     = $baseFileName. ".html";
	    
	    # test to see if I can write temp file. If not, do nothing:
	    if (open(OUTPUT_FILE_TEMP, ">$outputFileNameTemp")==0 ) 
	      {
		  $MONITOR && do 
		    { 
			print STDERR ("\nI cannot write/overwrite $outputFileNameTemp!\n");
			print STDERR ("I will stop processing this one!\a\n"); 
		    };
		  close(OUTPUT_FILE_TEMP);
	      } elsif (open(OUTPUT_FILE, ">$outputFileName")==0) {
		  $MONITOR && do
		    {
			print STDERR ("\nI cannot write/overwrite $outputFileName!\n");
			print STDERR ("I will stop processing this one!\a\n"); 
		    };
		  close(OUTPUT_FILE);
	      } else {
		  # Okay to POST
		  # close first
		  close(OUTPUT_FILE_TEMP);
		  close(OUTPUT_FILE);
		  
		  $typeDo=$extension;
		  $typeDo=~tr/a-z/A-Z/;
		  # rnc:  $options is no longer needed when using LWP request POST
		  #
		  #	    $options = "-u \"http:\/\/$blastServerURL\" -f PROGRAM=$program DATALIB=$database 
		  #               INPUT_TYPE=\"sequence in fasta format\" SEQUENCE=\"$sequenceInFastaFormat\" 
		  #               EXPECT=$expect FILTER=$filter $NCBI_giOption DESCRIPTIONS=$descriptions 
		  #               ALIGNMENTS=$alignments $graphicalOverviewOption GENETICS_CODE=$queryGeneticCode 
		  #               OTHER_ADVANCED=\"$advancedOptions\" MATRIX=$matrix";
		  
		  
		  # rnc:  Following options added to POST blocks so are no longer required with LWP
		  #
		  #	    $version eq '1'   and $options .= " CUTOFF=default";
		  #	    $version eq '2'   and $options .= " $gappedAlignmentOption";
		  #	    $version eq 'psi' and $options .= " E_THRESH=$expectPSI";
		  
		  #	    print STDERR "\nversion = $version\noptions: \n$options";<STDIN>;
		  
		  $MONITOR && print STDERR "Sending sequence $id to blast server.\n(Prog=$typeDo, Version=$version, DB=$database, Return=$resultReturnProtocol)\n";
		  
		  # rnc:  set up request using LWP
		  $ua = LWP::UserAgent->new();
		  $req = ""; # rnc:    empty to start
		  
		  # rnc:  assuming the user has proxyServerURL entry correct right above
		  #       so only doing a cursory check here that the address starts with
		  #       http: and ends with :somenumbers, maybe should ping the name
		  #       or something as the request will just hang there a long time if
		  #       proxy server 
		  if ($proxyServerURL)
		    {
			# check proxyServerURL starts with http: - this is a very minimal verification
			if ($proxyServerURL =~ m/^(\s+)?http:/i)
			  {
			      ($ua->proxy(['http', 'ftp'] => $proxyServerURL)) 
				&& ($MONITOR && print STDERR "Using Proxy Server: $proxyServerURL\n");
			  } else {
			      $_blastObj->throw("Invalid Proxy Server Name: $proxyServerURL",
						"It must be in form \"http://servername:port\"");
			  }
		    }
		  
		  # rnc:	changing www vs email, post protocol same for both but handle response differently	
		  # if return via www
		  #	    if ($resultReturnProtocol eq "www") {
		  
		  
		  if ($resultReturnProtocol eq "email") 
		    {
			$graphicalOverviewOption = ''; 		# rnc:	no graphic for email
			$emailOption = 'is_set';	
			$MONITOR && print STDERR "Results will be sent via e-mail: $emailAddress.\n";
		    } else {
			$emailAddress = '';	# rnc:	make sure address is empty for www response
			$resultTypeOption = '';	# rnc:	make sure HTML setting for email reply is empty
			$MONITOR && print STDERR "Be patient, especially during peak hours :-)...\n";
		    }
		  
		  
		  # Run the Blast:
		  
		  # rnc:  Commented out postclient related lines, using LWP now
		  #
		  #		system("$pathOfPostclient $options > $outputFileNameTemp") == 0
		  #		  or $_blastObj->throw("postclient.pl terminated abnormally.",
		  #				       "Check the executable $pathOfPostclient\n".
		  #				       "See $outputFileNameTemp for ".
		  #				       "possible error message.");
		  
		  # rnc:  have changed variables $NCBI_giOption etc above so will comment out here now
		  #        - delete this stuff later		
		  #       rnc:  if this works then should really change these further up and remove here
		  #       cut variable  $NCBI_giOption down 
		  #       ($NCBI_giOption eq "NCBI_GI=is_set") && ($NCBI_giOption = "is_set");
		  #       cut variable  $graphicalOverviewOption down 
		  #       ($graphicalOverviewOption eq "OVERVIEW=is_set") && ($graphicalOverviewOption = "is_set");
		  #       cut variable  $graphicalOverviewOption down 
		  #       ($gappedAlignmentOption eq "UNGAPPED_ALIGNMENT=is_set") && ($gappedAlignmentOption = "is_set");
		  
		  # rnc:  added info message to user about Posting data as they do have a wait occasionally
		  $MONITOR && print STDERR "\nPosting data to BLAST server...\n";
		  
		  # rnc:  formulate the request block - using 5 different blocks for blast versions
		  #       -version 1, -version 2, and -version 3 recognised by bioPerl as well as 
		  #       -version washu for default remote WashU Blast2, 
		  #       -version psi for PSI-Blast2 (PHI-Blast2 not yet supported).
		  #	No doubt there is a better way to do this but this way is what my present
		  #	experience with perl dictates I must use. It works this way if a bit ugly. 
		  #	WIll work on improvement later
		  #       MORE IMPORTANT - I am not sure I have ALL the options that the respective
		  #       BLAST servers accept - MUST find out and put in with default values for 
		  #       complete functionality if required
		  
		  # rnc:	-version 1: ungapped blast at NCBI using old blast - runs at NCBI but no longer supported
		  if ($version eq '1') {

		      $req = &_get_request(
					   PROGRAM         =>      $program,
					   DATALIB         =>      $database,
					   INPUT_TYPE      =>      $input_type,
					   SEQUENCE        =>      $sequenceInFastaFormat,
					   EXPECT          =>      $expect,
					   CUTOFF          =>      $cutoff,
					   MATRIX          =>      $matrix,             
					   STRAND          =>      $strand,
					   FILTER          =>      $filter,
					   HISTOGRAM       =>      $histogram,
					   NCBI_GI         =>      $NCBI_giOption,
					   DESCRIPTIONS    =>      $descriptions,
					   ALIGNMENTS      =>      $alignments,
					   ADVANCED  =>      $advancedOptions,				
					  );
		  }
		  
		  # rnc:	NCBI Blast2 server - does both basis blast2 and advanced blast2 
		  if ($version eq '2') {
			
		      $req = &_get_request(
					   PROGRAM                 =>      $program,
					   DATALIB                 =>      $database,
					   UNGAPPED_ALIGNMENT      => $gappedAlignmentOption,
					   INPUT_TYPE              =>      $input_type,
					   # rnc: FSET is in blast2 basic but not blast2 advanced
					   # have not added this anywhere yet, is it needed or will FILTER do it alone
					   #FSET		      =>	$fset, 
					   SEQUENCE                =>      $sequenceInFastaFormat,
					   GI_LIST                 =>      $gi_list,
					   LIST_ORG                =>      $list_org,
					   EXPECT                  =>      $expect,
					   FILTER                  =>      $filter,
					   NCBI_GI                 =>      $NCBI_giOption,
					   OVERVIEW                =>      $graphicalOverviewOption,
					   DESCRIPTIONS            =>      $descriptions,
					   ALIGNMENTS              =>      $alignments,
					   ALIGNMENT_VIEW          =>      $alignment_view,
					   GENETICS_CODE           =>      $queryGeneticCode,
			# rnc:	Blast2 advanced used MAT_PARAM which has value such as 'BLOSUM62	 11	 1'
			# 	This version of Webblast.pm is not dealing this with - will sought out later
			# 	any comments on MAT_PARAM vs MATRIX here ?
			# 	is MATRIX a recognised Blast2 parameter or should it be MAT_PARAM 
			# 	print "MAT_PARAM		=>	$mat_param,
					   MATRIX                  =>      $matrix,             
					   OTHER_ADVANCED          =>      $advancedOptions,
					  );
		  }
		  
		  # NCBI PSI-BLAST2
		  if ($version eq 'psi') {
			# rnc:	make sure that $program = blastp as is it only accepted option
			# 	probably done above, I need to check that
			$program = 'blastp';
		      $req = &_get_request(
					   PROGRAM                 =>      $program,
					   DATALIB                 =>      $database,
					   GAPPED_ALIGNMENT	=>	$gappedAlignmentOption,
					   INPUT_TYPE              =>      $input_type,
					   SEQUENCE                =>      $sequenceInFastaFormat,
					   EXPECT                  =>      $expect,
					   FILTER                  =>      $filter,
					   NCBI_GI                 =>      $NCBI_giOption,
					   GRAPHIC_OVERVIEW	=>	$graphicalOverviewOption,
					   DESCRIPTIONS            =>      $descriptions,
					   ALIGNMENTS              =>      $alignments,
					   E_THRESH		=>	$expectPSI,
			# rnc:	PSI-Blast2 advanced used MAT_PARAM which has value such as 'BLOSUM62	 11	 1'
			# 	This version of Webblast.pm is not dealing this with - will sought out later
			# 	any comments on MAT_PARAM vs MATRIX here ?
			# 	is MATRIX a recognised PSI-Blast2 parameter or should it be MAT_PARAM 
			# 	print FH "MAT_PARAM		=>	$mat_param,
					   MATRIX                  =>      $matrix,             
					   OTHER_ADVANCED          =>      $advancedOptions,
					  );
		    }
		  
		  #rnc:	Default POST block for when user has supplied a www server name
		  if ($version eq 'custom_server') {
		      $req = &_get_request(
					   PROGRAM         =>      $program,
					   DATALIB         =>      $database,
					   INPUT_TYPE      =>      $input_type,
					   SEQUENCE        =>      $sequenceInFastaFormat,
					   EXPECT          =>      $expect,
					   FILTER          =>      $filter,
					   NCBI_GI         =>      $NCBI_giOption,
					   DESCRIPTIONS    =>      $descriptions,
					   ALIGNMENTS      =>      $alignments,
					   OVERVIEW        =>      $graphicalOverviewOption,
					   GENETICS_CODE   =>      $queryGeneticCode,
					   OTHER_ADVANCED  =>      $advancedOptions,
					   MATRIX          =>      $matrix,             
					  );
		  }
		  
		  # rnc: POST block for WashU Blast to default WashU Blast2 remote server
		  #	enables WashU Blast2 specific commands
		  #	NOTE: the params here are not even remotely correct at present!
		  if ($version =~ /wash|wu/i) {
		      $req = &_get_request (
					    PROGRAM         =>      $program,
					    DATALIB         =>      $database,
					    INPUT_TYPE      =>      $input_type,
					    SEQUENCE        =>      $sequenceInFastaFormat,
					    EXPECT          =>      $expect,
					    FILTER          =>      $filter,
					    NCBI_GI         =>      $NCBI_giOption,
					    OVERVIEW        =>      $graphicalOverviewOption,
					    DESCRIPTIONS    =>      $descriptions,
					    ALIGNMENTS      =>      $alignments,
					    GENETICS_CODE   =>      $queryGeneticCode,
					    OTHER_ADVANCED  =>      $advancedOptions,
					    MATRIX          =>      $matrix,             
					   );
		  }
		  
		  # rnc:  deal with response from server
		  $response = $ua->request($req);
		  
		  if ($response->is_success)
		    {
			# rnc:  this prints report to screen but that is a diagnostic for me only at
			#       this point - must diable later and let the parser or calling script
			#       deal with how user wants data to be displayed, ie parsed through 
			#       bioperl or send to netscape for example
			#$MONITOR && print STDERR $response->content;
			
			# rnc:  write out the report to $outputFileNameTemp
			#       The ability to write to this file has already been checked above
			#       so not repeating here, no error capture here - should be!!
			open(FH, ">$outputFileNameTemp");
			print FH $response->content;
			close(FH);
		    } else {
			$MONITOR && do 
			  {
			      # rnc:  something is wrong!
			      #       write error to file
			      open(FH, ">$outputFileNameTemp");
			      print FH $response->error_as_HTML;
			      close(FH);
			      $MONITOR && print STDERR $response->error_as_HTML; 
			  }
		    }
		  
		  #		$MONITOR && print STDERR "Cleaning the file a bit...\n";
		  
		  
		  if ($resultReturnProtocol eq 'www') 
		    {
			
			if (&_removeJunkTagAndText($outputFileNameTemp, $outputFileName)==1) 
			  {
			      unlink $outputFileNameTemp;
			      $outFile = $outputFileName;
			      $MONITOR && print STDERR "$typeDo result saved to file $outputFileName\n";
			      
			      # If html flag is off, then call _removeHTMLtags
			      # If the report is to be parsed, the non-HTML version must be saved
			      # in @_outFiles.
			      if ($resultTypeInHTMLflag=~m/^off$/i || $resultTypeInHTMLflag =~ m/^both$/i ) 
				{
				    $outFile = &_removeHTMLtags($outputFileName, $baseFileName);
				}	
			      $countSuccessfulBlasting++;
			  } else {
			      # SAC: (Note) _removeJunkTagAndText() failed. 
			      #      Going with the temp file. Repercussions?
			      unlink $outputFileName;
			      if(-s $outputFileNameTemp) 
				{
				    $outFile = $outputFileNameTemp;
# Updated for the operation of the NCBI queueing system
# which now always results in the temporary file being saved.
				    $MONITOR && print STDERR "Blast submission result saved to file $outputFileNameTemp\n"; 
# Formerly, when the actual report was returned:
#				    $MONITOR && print STDERR "Possibly bad Blast result saved to file $outputFileNameTemp\n"; 
				} else {
				    unlink $outputFileNameTemp;
				    $MONITOR && print STDERR "No Blast results saved.\n";
				}
			  }
			
			# SAC: Save the name of the Blast output file.
			push @_outFiles, $outFile;
			
		    } else {
			&_removeJunkTagAndTextForEmailResponse($outputFileNameTemp);
			$countSuccessfulBlasting++;
			push @_outFiles, 'email';
		    }
	      } # end of okay to post
	    
	    if ($_errmsg and $count < scalar(@queue)) 
	      {
		  print STDERR "\n$_errmsg\n";
		  $_errmsg = '';
	      }	
	} # end foreach $seq()
	
	if(not $countSuccessfulBlasting) 
	  {
	      $_errmsg ||= "Blast failed: No file created or a error occured upon submission.";
	      $_blastObj->throw($_errmsg);
	  }
	
	if ($resultReturnProtocol eq "www") 
	  {
	      $MONITOR && print STDERR "\n$countSuccessfulBlasting file(s) have been ${typeDo}ed successfully!\n\n";
	  } else {
	      $MONITOR && print STDERR "$countSuccessfulBlasting file(s) have been sent to blast server successfully\n\n";
	  }
	
	#    $MONITOR && print STDERR "All done!\a\n\n";
    }


#------------------
sub _get_request {
#------------------
# SAC: Factoring out common methodology implemented by RNC.

    my (%param) = @_;
    my %extra_param = ('EMAIL' => $emailOption,
		       'PATH'  => $emailAddress,
		       'HTML'  => $resultTypeOption
		      );

    # rnc: write out values to $outputFileNameTemp as diagnostic tool
    # sac: warn if can't write these to file.
    if( open(FH, ">$outputFileNamePost")) {
        local $^W = 0;
	print FH "Posting address = $blastServerURL\n\n";
	foreach (sort keys %param) {
	    printf FH "%-20s %s\n", $_, $param{$_};
	}
	foreach (sort keys %extra_param) {
	    printf FH "%-20s %s\n", $_, $param{$_};
	}
	close FH;
    } else {
	$_blastObj->warn("Can't write Blast posting parameters to $outputFileNamePost: $!");
    }

    %param = (%param, %extra_param) if $resultReturnProtocol eq 'email';

    return POST $blastServerURL, [ %param ];
}


#------------------
sub _validate_seq {
#------------------
    # SAC: sequence error checking using methods on the sequence object.
    # No longer need Alex's various file checking functions since we are
    # obtaining the data directly from sequence objects.
    # Throws  : Exception if the sequence type is incompatible with the 
    #           version of Blast to be run or the database selected.
# Returns : zero if sequence length is too short or too long
#           otherwise returns 1.
    my $seq = shift;
    my $type = $seq->alphabet;
    my $len = $seq->length;

    # Verify that sequence type is correct type for selected program.
    if( (defined $type and $type !~ /unknown/i) 
	and (
	     ($type =~ /[dr]na/i and ($program =~ m/^(blastp|tblastn)$/i))
	     or
	     ($type =~ /amino/i and ($program =~ m/^(blast[nx]|tblastx|dbest)$/i))
	    ))  {
	$_blastObj->throw("The sequence of ${\$seq->id} is wrong type for program $program");
    }

    # Verify that sequence type is correct type for selected database.
    if(($type =~ /[dr]na/i and not (grep $database eq $_, @Blast_dbn_remote))
       or
       ($type =~ /amino/i and not (grep $database eq $_, @Blast_dbp_remote))
       ) {
	$_blastObj->throw("The sequence of ${\$seq->id} is wrong type for database $database");
    }

    # SAC: added more detailed error message.
    if ($len < $minimumLength) {
	$MONITOR && print STDERR ("\nThe sequence length of ${\$seq->id} ($len) too short for blasting: min length = $minimumLength\a\n\n");
	return 0;
    } elsif( $len > $maximumLength) {
	$MONITOR && print STDERR ("\nThe sequence length of ${\$seq->id} ($len) is too long for blasting: max length = $maximumLength\a\n\n");
	return 0;
    }
    1;
}



#--------------------------
sub _removeJunkTagAndText {
#--------------------------
# Purpose : removes junk tags+text from the file, and also retrieve .gif files if graphical is on
# also change relative URL of entries to absolute path.
# return 0 and if any error! or 1 if success
#
# SAC: Added some blocks to detect server errors.
#      Cleaned up some regexps.
#      Reduced the dependency on NCBI output format.
#      Using class-scope private variable to pass error msg back to caller.

   my ($inputFileName, $outputFileName)=@_;
   my ($line, $count);

   # test if files are 0 size!
   if (!-s $inputFileName) {
       $_errmsg = "File $inputFileName is empty! Possible reasons: \nblast server is too busy".
	   " or postclient.pl does not work.";
      return 0;
   }

   # test if files are readable and writable
   if (open(INPUT_FILE_FOR_REMOVE, $inputFileName)==0) {
      $_errmsg = "File $inputFileName is not readable! Can't process it.";
      close(INPUT_FILE_FOR_REMOVE);
      return 0;
   }
   # the following file is writable since I have tested already
   open(OUTPUT_FILE_FOR_REMOVE, ">$outputFileName");

   print OUTPUT_FILE_FOR_REMOVE "<HTML><PRE>\n";

   ## BEGIN EXTRACTING REPORT DATA.

   # ignore header until the name/version of program BLASTN 2.0.3 [Nov-14-1997]
   while (1) {
      $line=<INPUT_FILE_FOR_REMOVE>;

      # SAC: New server error detection block.
      if($line =~ /Error \d+ .+?/ or
	 $line =~ /Internal Server Error/i)  {
	  my $msg = $&;
	  my @remainder = <INPUT_FILE_FOR_REMOVE>;
	  $_errmsg = $1."\n".join(' ',@remainder);
	  return 0;
      }

      # SAC: more general purpose regexp for program line.
      last if $line =~ m/T?BLAST[NPX]\s[\w.-]/;

      if (eof(INPUT_FILE_FOR_REMOVE)) {
          close(INPUT_FILE_FOR_REMOVE);
          close(OUTPUT_FILE_FOR_REMOVE);
          $_errmsg = "Blast file $inputFileName is incomplete: can't find program name.".
	             "\nPossibly an unrecognized report format." ;
          return 0;
      }
   }

   # copy BLASTN 2.0.3 [Nov-14-1997]
   print OUTPUT_FILE_FOR_REMOVE "$line\n";

   # add date/time stamp:
   print OUTPUT_FILE_FOR_REMOVE "\n<b>Done at $currentDate (local time)</b>\n";

   # SAC: The following commented out section is from the original webblast which
   # skips the "Reference:" text. It should be kept for documentation purposes.
   # Also, this code is sensitive to changes in the output format.

   # ignore until Query=   
#   while (1) {
#      $line=<INPUT_FILE_FOR_REMOVE>;
#      if ($line=~m/Query=/) {
#	   last;
#      }
#      if (eof(INPUT_FILE_FOR_REMOVE)) {
#	   close(INPUT_FILE_FOR_REMOVE);
#	   close(OUTPUT_FILE_FOR_REMOVE);
#	   $_errmsg = "Blast file $inputFileName was not completed!".
#		      "\nPossibly an unrecognized report format (2)." ;
#	   return 0;
#      }
#   }
#   # copy Query= line
#   print OUTPUT_FILE_FOR_REMOVE $line;
#
#   # copy next 6 lines, adding time stamp
#   for ($count=0; $count<6; $count++) {
#      $line=<INPUT_FILE_FOR_REMOVE>;
#      if (eof(INPUT_FILE_FOR_REMOVE)) {
#	    close(INPUT_FILE_FOR_REMOVE);
#	    close(OUTPUT_FILE_FOR_REMOVE);
#	    $_errmsg = "Blast file $inputFileName was not completed!".
#		       "\nPossibly an unrecognized report format (3)." ;
#	    return 0;
#      }
#      # add time to this line:
#      $line =~ s/^(Searching.+)$/$1\ndone at $currentDate (local time)/;
#
#      print OUTPUT_FILE_FOR_REMOVE $line;
#   }

   # copy all lines, adding absolute URL stem to all links:
   while (!eof(INPUT_FILE_FOR_REMOVE)) {
      $line=<INPUT_FILE_FOR_REMOVE>;
      # SAC: substituting in one step.
      # rnc: not possible with PSI as it uses ACTION=.. for continuation of posting
      #
      # ACTION="/cgi-bin/BLAST" line added to get PSI to continue
      $line =~ s|ACTION="/cgi-bin/BLAST|ACTION="http://$entryServerURLroot/cgi-bin/BLAST|;
      $line =~ s|<a href="/|<a href="http://$entryServerURLroot/|;
      print OUTPUT_FILE_FOR_REMOVE $line;
   }
   close(INPUT_FILE_FOR_REMOVE);
   close(OUTPUT_FILE_FOR_REMOVE);

   return 1;
}



#--------------------------------------------
sub _removeJunkTagAndTextForEmailResponse {
#--------------------------------------------
# Purpose :b removes junk tags+text from the file, and also retrieve .gif files if graphical is on
# return 0 if any error! or 1 if success

   my ($inputFileName)=@_;
   my ($line);

   # test if files are readable and writable
   if (open(INPUT_FILE_FOR_REMOVE, $inputFileName)==0) {
      close(INPUT_FILE_FOR_REMOVE);
      $_errmsg = "Email response file $inputFileName is not readable! Can't process it.";
      return 0;
   }

   print "\n";
   # copy lines until cgi-bin/imagemap/BLAST/blast_form.map
   while (!eof(INPUT_FILE_FOR_REMOVE)) {
      $line=<INPUT_FILE_FOR_REMOVE>;
      $line=~s/<[^>]+>//g;
      print "$line";
   }

   close(INPUT_FILE_FOR_REMOVE);
   unlink $inputFileName;

   return 1;
}


#----------------------
sub _removeHTMLtags {
#----------------------
# Purpose : removes all html tags!
# return 0 if any error or 1 if success
# SAC:
#   Modifications:
#    - Uses Bio::Tools::Blast::HTML::strip_html() to remove HTML.
#      Stripping now occurs on the whole report in one chunk instead
#      of line-by-line (line-by-line is more error-prone).
#    - Now returns the name of the new file (or zero on failure).
#    - adding time stamp, as was done in _removeJunkTagAndText().
#    - cleaned up the HTML removal code a bit
#    - fixed bug due to failure to close the filehandle
#      (final close calls were missing '_TAG'.
#    - Moved the actual HTML-removing code into Bio::Tools::HTML
#
#    Note that $inputFileName is the ".program.html" file
#    which is deleted by this function if not saving the HTML version. 

   my ($inputFileName, $baseFileName )=@_;
   my ($line, $lineWithTag, $fileNameTagRemoved);
#   $fileNameTagRemoved  = $baseFileName . ".$extension" . ".v2";
   # SAC: don't need '.v2' since it was saved with a .html extension.
   $fileNameTagRemoved   = $baseFileName;

   # test if files are 0 size!
   if (!-s $inputFileName) {
      $_errmsg = "File $inputFileName is empty! Can't remove HTML tags.";
      return 0;
   }

   # test if temp file exist already or not writable
   if (-e $fileNameTagRemoved && !-w $fileNameTagRemoved) {
      $_errmsg = "Can't create temp file $fileNameTagRemoved! Can't process $inputFileName.";
      return 0;
   }

#   $MONITOR && print STDERR ("removing all html tags...\n");
   open(INPUT_FILE_FOR_REMOVE_TAG, "$inputFileName");
   open(OUTPUT_FILE_FOR_REMOVE_TAG, ">$fileNameTagRemoved");

   while (1) {
         $line=<INPUT_FILE_FOR_REMOVE_TAG>;
	 # SAC: more general purpose regexp for program line.
	 if ($line =~ m/T?BLAST[NPX]\s[\w.-]/) {
#         if ($line=~m/[a-z]{5,7} \d+\.\d+\.\d+ \[[a-z]{3}-\d{2}-\d{4}\]/i) {
             # remove tags, print to output, then exit loop
	     $line =~ s/<[^>]+>//g;
             print OUTPUT_FILE_FOR_REMOVE_TAG $line;
             last;
         }

         if (eof(INPUT_FILE_FOR_REMOVE_TAG)) {
             $_errmsg = "Blast file $inputFileName is incomplete: Can't find program name.";
             close(INPUT_FILE_FOR_REMOVE_TAG);
             close(OUTPUT_FILE_FOR_REMOVE_TAG);
             return 0;
         }
    }

   # removing html tags
   # SAC: Slurping the remainder of the file and stripping the whole thing at once
   #      using Bio::Tools::Blast::HTML::strip_html().
   my @lines = <INPUT_FILE_FOR_REMOVE_TAG>;
   $line = join('', @lines);
   &strip_html(\$line); 
   print OUTPUT_FILE_FOR_REMOVE_TAG $line;

   close(INPUT_FILE_FOR_REMOVE_TAG);
   close(OUTPUT_FILE_FOR_REMOVE_TAG);

   $MONITOR && print STDERR ("Output without HTML tags saved to $fileNameTagRemoved\n");

   # now if not ON or BOTH!
   # remove the .html file.
   if ($resultTypeInHTMLflag=~m/^off$/i) {
      unlink $inputFileName; 
      $MONITOR && print STDERR ("Deleted HTML file.\n");
  } 

   
   return $fileNameTagRemoved;
}

######################  Blast Variables  ########################

=head1 APPENDIX 2: Parameter listings

Parameters for Blast (NCBI ungapped, no longer supported by NCBI so
should dicontinue use of ungapped blast), Blast2 (NCBI), PSI-Blast2
(NCBI). WashU-Blast2 has yet to be added as does PHI-Blast2 (NCBI).

These lists of parameters for posting to blast servers were
obtained directly from the respective WWW forms for each server.


=head2 Basic ungapped BLAST Search Server Parameters


PROGRAM 
[default value]:blastn
	blastp tblastn tblastx blastx

DATALIB
[default value]:nr
	month swissprot dbest dbsts pdb vector kabat mito alu epd yeast gss htgs ecoli

INPUT_TYPE
[default value]:Sequence in FASTA format
	Accession or GI

SEQUENCE

EXPECT
[default value]:default
	0.0001 0.01 1 10 100 1000

CUTOFF
[default value]:default
	60 70 80 90 100 110

MATRIX
[default value]:default
	BLOSUM62 PAM40 PAM120 PAM250 IDENTITY

STRAND
[default value]:both
	top bottom

FILTER
[default value]:default
none dust SEG SEG+XNU XNU

HISTOGRAM
[default value]:''
HISTOGRAM

NCBI_GI
[default value]:""
NCBI_GI

DESCRIPTIONS
[default value]:default
0 10 50 100 250 500

ALIGNMENTS
[default value]:default
0 10 50 100 250 500

ADVANCED
[default value]:""

EMAIL
[default value]:''
IS_SET

PATH
[default value]:""

HTML
[default value]:''
HTML


=head2 Basic Blast 2

PROGRAM
[default value]:blastn
blastp blastx tblastn tblastx

DATALIB
[default value]:nr
month swissprot dbest  dbsts est_mouse est_human est_others pdb pat vector kabat mito alu epd yeast ecoli gss htgs

UNGAPPED_ALIGNMENT
[default value]:''
is_set

FSET
[default value]:is_set
''

OVERVIEW
[default value]:is_set
''

INPUT_TYPE
[default value]:Sequence in FASTA format
Accession or GI

SEQUENCE

EMAIL
[default value]:''
IS_SET

PATH
[default value]:""

HTML
[default value]:''
IS_SET


=head2 BLAST2 ADVANCED

PROGRAM
[default value]:blastn
blastp blastx tblastn tblastx 

DATALIB
[default value]:nr
 month swissprot dbest dbsts est_mouse est_human est_others pdb pat vector kabat mito alu epd yeast ecoli gss htgs

UNGAPPED_ALIGNMENT
[default value]:""
is_set

INPUT_TYPE
[default value]:Sequence in FASTA format
Accession or GI

SEQUENCE

GI_LIST
[default value]:(None)
Arabidopsis thaliana   Bacillus subtilis   Bos taurus   Caenorhabditis elegans   Danio rerio   Dictyostelium discoideum   Drosophila melanogaster   Escherichia coli   Gallus gallus   Homo sapiens   Human immunodeficiency virus type 1   Mus musculus   Oryctolagus cuniculus   Oryza sativa   Ovis aries   Plasmodium falciparum   Rattus norvegicus   Saccharomyces cerevisiae   Schizosaccharomyces pombe   
Simian immunodeficiency virus   Xenopus laevis   Zea mays

LIST_ORG

EXPECT
[default value]:10
0.0001 0.01 1 10 100 1000

FILTER
[default value]:default
none


NCBI_GI
[default value]:''
is_set


OVERVIEW
[default value]:is_set
''


DESCRIPTIONS
[default value]:500
0 10 50 100 250 500

ALIGNMENTS
[default value]:500
0 10 50 100 250 500

ALIGNMENT_VIEW
[default value]:0	#Pairwise
1	#master-slave with identities
2	#master-slave without identities
3	#flat master-slave with identities
4	#flat master-slave without identities

GENETIC_CODE
[default value]:Standard (1)
Vertebrate Mitochondrial (2) Yeast Mitochondrial (3) Mold Mitochondrial; ... (4) 
Invertebrate Mitochondrial (5) Ciliate Nuclear; ... (6) Echinoderm Mitochondrial (9) 
Euplotid Nuclear (10) Bacterial (11) Alternative Yeast Nuclear (12) 
Ascidian Mitochondrial (13) Flatworm Mitochondrial (14) Blepharisma Macronuclear (15) 

MAT_PARAM
[default value]:BLOSUM62	 11	 1
PAM30	 9	 1  
PAM70	 10	 1  
BLOSUM80	 10	 1 
BLOSUM62	 11	 1
BLOSUM45	 14	 2
PAM30	 7	 2  
PAM30	 6	 2
PAM30	 5	 2
PAM30	 10	 1
PAM30	 9	 1  #recommended
PAM30	 8	 1  
PAM70	 8	 2  
PAM70	 7	 2  
PAM70	 6	 2  
PAM70	 11	 1  
PAM70	 10	 1  #recommended
PAM70	 9	 1  
BLOSUM80	 8	 2  
BLOSUM80	 7	 2  
BLOSUM80	 6	 2   
BLOSUM80	 11	 1  
BLOSUM80	 10	 1  #recommended
BLOSUM80	 9	 1  
BLOSUM62	 9	 2  
BLOSUM62	 8	 2  
BLOSUM62	 7	 2  
BLOSUM62	 12	 1  
BLOSUM62	 11	 1  #recommended
BLOSUM62	 10	 1  
BLOSUM45	 13	 3  
BLOSUM45	 12	 3  
BLOSUM45	 11	 3  
BLOSUM45	 10	 3  
BLOSUM45	 15	 2  
BLOSUM45	 14	 2  #recommended
BLOSUM45	 13	 2  
BLOSUM45	 12	 2  
BLOSUM45	 19	 1  
BLOSUM45	 18	 1  
BLOSUM45	 17	 1  
BLOSUM45	 16	 1  

OTHER_ADVANCED
[default value]:""

EMAIL
[default value]:''
IS_SET

PATH
[default value]:""

HTML
[default value]:''
IS_SET



=head2 PSI BLAST2


PROGRAM
[default value]:blastp

DATALIB
[default value]:nr
month swissprot pdb kabat alu yeast ecoli

GAPPED_ALIGNMENT
[default value]:is_set
''

INPUT_TYPE
[default value]:Sequence in FASTA format
Accession or GI

SEQUENCE

EXPECT
[default value]:10
0.0001 0.01 1 10  100  1000

FILTER
[default value]:default
none

NCBI_GI
[default value]:''
is_set

GRAPHIC_OVERVIEW
[default value]:is_set
''

DESCRIPTIONS
[default value]:500
0  10 50 100 250 500

ALIGNMENTS
[default value]:500
0  10 50 100 250 500

E_THRESH
[default value]:0.001
#max value is 10

MAT_PARAM
[default value]:BLOSUM62	 11	 1
PAM30	 9	 1
 PAM70	 10	 1
 BLOSUM80	 10	 1
 BLOSUM62	 11	 1
 BLOSUM45	 14	 2
 PAM30	 7	 2
 PAM30	 6	 2
 PAM30	 5	 2
 PAM30	 10	 1
 PAM30	 9	 1
 PAM30	 8	 1
 PAM70	 8	 2
 PAM70	 7	 2
 PAM70	 6	 2
 PAM70	 11	 1
 PAM70	 10	 1
 PAM70	 9	 1
 BLOSUM80	 8	 2
 BLOSUM80	 7	 2
 BLOSUM80	 6	 2
 BLOSUM80	 11	 1
 BLOSUM80	 10	 1
 BLOSUM80	 9	 1
 BLOSUM62	 9	 2
 BLOSUM62	 8	 2
 BLOSUM62	 7	 2
 BLOSUM62	 12	 1
 BLOSUM62	 11	 1
 BLOSUM62	 10	 1
 BLOSUM45	 13	 3
 BLOSUM45	 12	 3
 BLOSUM45	 11	 3
 BLOSUM45	 10	 3
 BLOSUM45	 15	 2
 BLOSUM45	 14	 2
 BLOSUM45	 13	 2 
BLOSUM45	 12	 2
BLOSUM45	 19	 1
BLOSUM45	 18	 1
BLOSUM45	 17	 1
BLOSUM45	 16	 1

OTHER_ADVANCED
[default value]:""

=head2 WashU BLAST2

WU-Blast2 Database Searches
http://www2.ebi.ac.uk/blast2/

email
""

title
Sequence

srchtype
interactive
email

database
swall
swissprot
swnew
trembl 
tremblnew
pdb 
gpcrdb
prints
HLAprot
embl
emnew
est
igvec
emvec
imgt
HLAnuc

program
WU-blastp
WU-blastx
WU-blastn

matrix
blosum62 
blosum30 
blosum35
blosum40 
blosum45 
blosum50 
blosum65 
blosum70
blosum75 
blosum80 
blosum85 
blosum90
blosum100 
GONNET 
pam10 
pam20 
pam30 
pam40
pam50 
pam60 
pam70 
pam80 
pam90 
pam100
pam110 
pam120 
pam130 
pam140 
pam150
pam160 
pam170 
pam180 
pam190 
pam200 
pam210
pam220 
pam230 
pam240 
pam250 
pam260 
pam270
pam280 
pam290 
pam300 
pam310 
pam320 
pam330
pam340 
pam350 
pam360 
pam370 
pam380 
pam390
pam400 
pam410 
pam420 
pam430 
pam440 
pam450
pam460 
pam470 
pam480 
pam490 
pam500 

strand
default
top
bottom

exp
default 
1.0 
10 
100 
1000 

filter
none
seg
xnu
seg+xnu
dust

echofilter
no
yes

histogram
no
yes

stats
sump
poisson

sort
pvalue
count
highscore
totalscore 

scores
default
5
10
20
50 
100 
150 
200 
250 

numal
default
5
10
20
50 
100 
150 
200 
250 

sequence


=cut


1;
__END__

############################################################################
#                             END OF CLASS                              
############################################################################

# MODIFICATION HISTORY :
#  1.24, 15 May 2000 sac:
#     -- Uses SeqIO and IO::Scalar instead of Bio::Seq->layout()
#        to get a Fasta-formatted sequence string.
#     -- Removed message about possibly bad blast report
#        which now always happens given the NCBI queueing system
#
#  1.23, 15 Oct 1999, sac:
#     -- Updated the URLs for the NCBI servers.
#
#  1.22, 26 Jun 1999, sac:
#      -- Bug fix in _blast(): Removing characters from temp file name 
#         that could confuse shell (reported by Bradford Powell).
#      -- Bug fix in _set_options() to reset the blastServerURL between
#         invocations (reported by James Diggans).
#      -- Added general POD for the module itself.
#
#  1.2, 20 Apr 1999, sac:
#      -- Added support for new options introduced by RNC's modifications
#         [changes in blast_remote(), _set_options()].
#      -- Condensed RNC's code [_blast(), created _get_request()].
#      -- Cleaned up comments and variable names.
#
#  1.1, 25 Feb 1999, rnc:
#     Search for comments marked "rnc:"
#     My modifications were as follows 
#       -- commented out references to postclient.pl
#       -- added LWP and HTTP::Request, supporting code for POST
#       -- user needs to enter proxy server address and port or set to ""
#		eg 
#		my $proxyServerURL = "http://proxy.marc.cri.nz:8080";
#		(maybe this is better placed in an "ini" like config file ?)
#		SEARCH FOR "PROXY SERVER NAME MUST BE ADDED HERE"
#		Could have left postclient.pl related code in and used
#		a conditional branch based on value of $proxyServerURL but
#		that would have meant two systems and I considered
#		that could lead to confusion eventually and LWP based 
#		may be better in long term ??
#       -- changed the values of some variables as they are no longer passed to
#		postclient.pl as an options string but are used directly
#		in the UserAgent eg $NCBI_giOption no longer eq "NCBI_GI=is_set" 
#		but $NCBI_giOption = "is_set"
#		and used in the request as: NCBI_GI         =>      $NCBI_giOption,
#       -- added several options to POST for BLAST2 which I am unsure are supported by
#		bioperl (passed by it actually) as yet - help ? (see NCBI blast2 POST block)
#       -- added a default WashU Blast2 server variable and set to 
#		'blastServerURLWashA'	=> 'www2.ebi.ac.uk/cgi-bin/newblast2.pl',
#       -- www2.ebi.ac.uk uses different parameter names in cgi script so probably 
#		should look for a alternate one (an alternative could be 
#		http://dove.embl-heidelberg.de/Blast2 
#		and accessed by -ver dove but I have not set this up yet), 
#		www2.ebi.ac.uk server is selected by passing -ver washa
#		dove.embl-heidelberg.de COULD BE selected by passing -ver 
#		dove (if and when code added)
#	-- if user passes a server name which is not internal to this script they
#		get a default post block at present (yet to be written)
#       PROBLEMS
#       -- 1. ver 3 does not seem to work due to a "400 Bad Request", large because there is no
#		generalised POST request included for a server passed to this script now - need to add this
#       -- 2. the POST blocks are individual in this version, they should be standardised and maybe
#		formualated as strings and eval block ?? - not sure on best way to do this, they are individual
#		blocks at this point just to get them working and tested independently
#       -- 3. code is not really modularised, dont have much experience at this - maybe someone else
#		can assist if this version of Webblast.pm is accepted over the previous one, code size can be 
#		cut down heaps with some more work - as stated above - its linear programming presently
#		so I can check everything is working with out complication of different blast servers
#		having dif reqs not to mention some users wanting email replies - will sought out later
#       -- 4. 	If using -ver psi, the POST works, a return is received but I have only tested with a seq
#		that returned "No hit found". BioPerl handled this by throwing an exception - Not very
#		useful 
#	-- 5.	Blast 2 is sent and returned okay but during parsing is throwing up following errors:
#		Use of uninitialized value at /usr/lib/perl5/site_perl/5.005/Bio/Tools/Blast/Sbjct.pm line 784.
#		Use of uninitialized value at /usr/lib/perl5/site_perl/5.005/Bio/Tools/Blast/Sbjct.pm line 1380.
#		Use of uninitialized value at /usr/lib/perl5/site_perl/5.005/Bio/Tools/Blast/Sbjct.pm line 1870.
#	-- 6.	If -email is used with ./run.pl to test then blast object throws except - this is not nice
#		in a demo, also -noparse does not work as excepted in run.pl demo 
#
#  1.0, 5 Jun 1998, sac: (initial release)
#    A decent amount of Alex's original code (version 0.92) remains, 
#    but has been extensively modified. Search for comments 
#    marked with "SAC:" for details of my modifications. 
#    Generally, my modifications were as follows:
#   
#  	 -- moved code from a stand-alone script into a module.
#  	 -- segmented script code into logical procedures.
#        -- automatically configure the path for postclient.pl.
#  	 -- works with Bio::Seq.pm objects instead of files.
#  	 -- added support for Blast1 and the beginnings of support 
#           for PSI Blast.
#  	 -- changed print() calls to print STDERR 
#           (only if $MONITOR is true).
#  	 -- changed system("rm") calls to unlink().
#  	 -- using 'my' instead of 'local' for subroutine-scope vars.
#  	 -- general cleaning up (got it to run with use strict).
#  	 -- changed `date` calls to $Util->date_format().
#        -- added leading underscores to private methods.
#        -- added advanced option W (word size).
#        -- uses Bio::Tools::Blast::HTML.pm for stripping HTML.
#        -- assorted changes in _removeHTMLtags().
#        -- added miscellaneous comments and POD.
#
# For additional modification notes and the latest version, 
# visit the distribution site:
#    http://bio.perl.org/Projects/Blast
#
#
#   Alex's NOTES FROM webblast 0.92: 
#     changes made to blast files:
#  	 1) delete reference lines and a few junk lines at the
#           beginning of the files
#  	 2) add time/date to search..... ...done line
#  	 3) change the relative URLs of entries to abosulute URLs
#     NEED TO DO:
#     retrieve gif files by getclient--for graphical overview
#     03/05/98: remove .seq and .v2 from output filenames
#     02/10/98: -html option controls more now: -html=ON|off|both
#     02/07/98: add a counter for successful blasting
#     01/24/98: now program can take files in fasta or plain text format
#     01/21/98: add new varible $pathOfPostclient
#  	 so that outside users can easily modify it
#  
