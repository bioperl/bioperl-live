#----------------------------------------------------------------------------
# PACKAGE : Bio::Tools::Blast.pm
# PURPOSE : To encapsulate code for running, parsing, and analyzing 
#           BLAST reports.
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED : March 1996
# REVISION: $Id$
# STATUS  : Alpha
# 
# For the latest version and documentation, visit:
#    http://bio.perl.org/Projects/Blast
#
# To generate documentation, run this module through pod2html
# (preferably from Perl v5.004 or better).
#
# Copyright (c) 1996-2000 Steve A. Chervitz. All Rights Reserved.
#           This module is free software; you can redistribute it and/or 
#           modify it under the same terms as Perl itself.
#----------------------------------------------------------------------------

package Bio::Tools::Blast;
use strict;
use Exporter;

use Bio::Tools::SeqAnal;
use Bio::Root::Global     qw(:std);
use Bio::Root::Utilities  qw(:obj); 

require 5.002;
use Carp;

use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS
            $ID $VERSION $Blast @Blast_programs $Revision $Newline);

@ISA        = qw( Bio::Tools::SeqAnal Exporter);
@EXPORT     = qw();
@EXPORT_OK  = qw($VERSION $Blast);
%EXPORT_TAGS = ( obj => [qw($Blast)],
		 std => [qw($Blast)]);


$ID = 'Bio::Tools::Blast';
$VERSION  = 0.09; 
$Revision = '$Id$';  #'

## Static Blast object. 
$Blast = {};
bless $Blast, $ID;
$Blast->{'_name'} = "Static Blast object";

@Blast_programs  = qw(blastp blastn blastx tblastn tblastx);

use vars qw($DEFAULT_MATRIX $DEFAULT_SIGNIF);
my $DEFAULT_MATRIX   = 'BLOSUM62';
my $DEFAULT_SIGNIF   = 999;# Value used as significance cutoff if none supplied.
my $MAX_HSP_OVERLAP  = 2;  # Used when tiling multiple HSPs.

## POD Documentation:

=head1 NAME

Bio::Tools::Blast.pm - Bioperl BLAST sequence analysis object

=head1 SYNOPSIS

=head2 Parsing Blast reports

Parse an existing Blast report from file:

    use Bio::Tools::Blast;

    $blastObj = Bio::Tools::Blast->new( -file   => '/tmp/blast.out',
					-parse  => 1,  
					-signif => '1e-10',
					);

Parse an existing Blast report from STDIN:

    $blastObj = Bio::Tools::Blast->new( -parse  => 1,  
					-signif => '1e-10',
					);

Then send a Blast report to your script via STDIN.

Full parameters for parsing Blast reports.

 %blastParam = ( 
		-run             => \%runParam,
		-file            => '',
		-parse           => 1,
		-signif          => 1e-5, 
		-filt_func       => \&my_filter,
		-min_len         => 15, 
		-check_all_hits  => 0,
		-strict          => 0,
		-stats           => 1,
		-best            => 0,
		-share           => 0,
		-exec_func       => \&process_blast,
		-save_array      => \@blast_objs,  # not used if -exce_func defined.
	       );

See L<parse>() for a description of parameters and see L<USAGE> for
more examples including how to parse streams containing multiple Blast
reports L<Using the Static $Blast Object>.

See L<Memory Usage Issues> for information about how to make Blast
parsing be more memory efficient.


=head2 Running Blast reports

Run a new Blast2 at NCBI and then parse it:

    %runParam = ( 
		  -method   => 'remote',
		  -prog     => 'blastp',
		  -database => 'swissprot',
		  -seqs     => [ $seq ],  # Bio::Seq.pm objects.
		  );     
 
    $blastObj = Bio::Tools::Blast->new( -run     => \%runParam,
					-parse   => 1,  
					-signif  => '1e-10',
					-strict  => 1,
					);

Full parameters for running Blasts at NCBI using Webblast.pm:

 %runParam = ( 
	      -method   => 'remote',
	      -prog     => 'blastp',
	      -version  => 2,      # BLAST2
	      -database =>'swissprot',
	      -html     => 0,
	      -seqs     => [ $seqObject ],  # Bio::Seq.pm object(s)
	      -descr    => 250,
	      -align    => 250,
	      -expect   => 10,
	      -gap      => 'on',
	      -matrix   => 'PAM250',
	      -email    => undef,  # don't send report via e-mail if parsing.
	      -filter   => undef,  # use default
	      -gap_c    => undef,  # use default
	      -gap_e    => undef,  # use default
	      -word     => undef,  # use default
	      -min_len  => undef,  # use default
	      );     

See L<run>() and L<USAGE> for more information about running Blasts.


=head2 HTML-formatting Blast reports

Print an HTML-formatted version of a Blast report:

    use Bio::Tools::Blast qw(:obj);

    $Blast->to_html($filename);
    $Blast->to_html(-file   => $filename,
		    -header => "<H1>Blast Results</H1>");
    $Blast->to_html(-file   => $filename,
		    -out    => \@array);  # store output
    $Blast->to_html();  # use STDIN

Results are sent directly to STDOUT unless an C<-out =E<gt> array_ref>
parameter is supplied. See L<to_html>() for details.


=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

=head1 DESCRIPTION

The Bio::Tools::Blast.pm module encapsulates data and methods for
running, parsing, and analyzing pre-existing BLAST reports. This
module defines an application programming interface (API) for working
with Blast reports. A Blast object is constructed from raw Blast
output and encapsulates the Blast results which can then be accessed
via the interface defined by the Blast object.

The ways in which researchers use Blast data are many and varied. This
module attempts to be general and flexible enough to accommodate
different uses. The Blast module API is still at an early stage of
evolution and I expect it to continue to evolve as new uses for Blast
data are developed. Your L<FEEDBACK> is welcome.

B<FEATURES:>

=over 2

=item * Supports NCBI Blast1.x, Blast2.x, and WashU-Blast2.x, gapped
and ungapped.

Can parse HTML-formatted as well as non-HTML-formatted reports.

=item * Launch new Blast analyses remotely or locally.

Blast objects can be constructed directly from the results of the
run. See L<run>().

=item * Construct Blast objects from pre-existing files or from a new run.

Build a Blast object from a single file or build multiple Blast
objects from an input stream containing multiple reports. See
L<parse>().

=item * Add hypertext links from a BLAST report.

See L<to_html>().

=item * Generate sequence and sequence alignment objects from HSP
sequences.

If you have Bio::Seq.pm and Bio::UnivAln.pm installed on your system,
they can be used for working with high-scoring segment pair (HSP)
sequences in the Blast alignment.  (A new version of Bio::Seq.pm is
included in the distribution, see L<INSTALLATION>).  For more
information about them, see:

    http://bio.perl.org/Projects/Sequence/
    http://bio.perl.org/Projects/SeqAlign/

=back

A variety of different data can be extracted from the Blast report by
querying the Blast.pm object. Some basic examples are given in the
L<USAGE> section. For some working scripts, see the links provided in
the L<DEMO SCRIPTS> section.

As a part of the incipient Bioperl framework, the Bio::Tools::Blast.pm
module inherits from B<Bio::Tools::SeqAnal.pm>, which provides some
generic functionality for biological sequence analysis. See the
documentation for that module for details (L<Links to related
modules>).


=head2 The BLAST Program

BLAST (Basic Local Alignment Search Tool) is a widely used algorithm
for performing rapid sequence similarity searches between a single DNA
or protein sequence and a large dataset of sequences.  BLAST analyses
are typically performed by dedicated remote servers, such as the ones
at the NCBI. Individual groups may also run the program on local
machines.

The Blast family includes 5 different programs:

              Query Seq        Database
             ------------     ----------
 blastp  --  protein          protein
 blastn  --  nucleotide       nucleotide 
 blastx  --  nucleotide*      protein
 tblastn --  protein          nucleotide*
 tblastx --  nucleotide*      nucleotide*
 
            * = dynamically translated in all reading frames, both strands

See L<References & Information about the BLAST program>.


=head2 Versions Supported

BLAST reports generated by different application front ends are similar
but not exactly the same. Blast reports are not intended to be exchange formats, 
making parsing software susceptible to obsolescence. This module aims to 
support BLAST reports generated by different implementations:

  Implementation    Latest version tested
  --------------    --------------------
  NCBI Blast1       1.4.11   [24-Nov-97] 
  NCBI Blast2       2.0.8    [Jan-5-1999]
  WashU-BLAST2      2.0a19MP [05-Feb-1998]
  GCG               1.4.8    [1-Feb-95]

Support for both gapped and ungapped versions is included. Currently, there
is only rudimentary support for PSI-BLAST in that these reports can be parsed but 
there is no special treatment of separate iteration rounds (they are all
merged together).


=head2 References & Information about the BLAST program

B<WEBSITES:>

   http://www.ncbi.nlm.nih.gov/BLAST/                 - Homepage at NCBI
   http://www.ncbi.nlm.nih.gov/BLAST/blast_help.html  - Help manual
   http://blast.wustl.edu/                            - WashU-Blast2


B<PUBLICATIONS:> (with PubMed links)


     Altschul S.F., Gish W., Miller W., Myers E.W., Lipman D.J. (1990).
     "Basic local alignment search tool", J Mol Biol 215: 403-410.

http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?uid=2231712&form=6&db=m&Dopt=r

     Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer,
     Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997).
     "Gapped BLAST and PSI-BLAST: a new generation of protein database 
     search programs", Nucleic Acids Res. 25:3389-3402.

http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?uid=9254694&form=6&db=m&Dopt=r

     Karlin, Samuel and Stephen F. Altschul (1990).  Methods  for
     assessing the statistical significance of molecular sequence
     features by using general scoring schemes. Proc. Natl. Acad.
     Sci. USA 87:2264-68.

http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?uid=2315319&form=6&db=m&Dopt=b

     Karlin, Samuel and Stephen F. Altschul (1993).  Applications
     and statistics for multiple high-scoring segments in molecu-
     lar sequences. Proc. Natl. Acad. Sci. USA 90:5873-7.

http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?uid=8390686&form=6&db=m&Dopt=b




=head1 USAGE

=head2 Creating Blast objects

A Blast object can be constructed from the contents of a Blast report
using a set of named parameters that specify significance criteria for
parsing.  The report data can be read in from an existing file
specified with the C<-file =E<gt> 'filename'> parameter or from a
STDIN stream containing potentially multiple Blast reports. If the
C<-file> parameter does not contain a valid filename, STDIN will be
used. Separate Blast objects will be created for each report in the
stream.

To parse the report, you must include a C<-parse =E<gt> 1> parameter
in addition to any other parsing parameters
See L<parse>() for a full description of parsing parameters.
To run a new report and then parse it, include a C<-run =E<gt> \%runParams> 
parameter containing a reference to a hash
that hold the parameters required by the L<run>() method.

The constructor for Blast objects is inherited from Bio::Tools::SeqAnal.pm.
See the B<_initialize>() method of that package for general information 
relevant to creating Blast objects. (The B<new>() method, inherited from
B<Bio::Root::Object.pm>, calls B<_initialize>(). See L<Links to related modules>).

The Blast object can read compressed (gzipped) Blast report
files. Compression/decompression uses the gzip or compress programs
that are standard on Unix systems and should not require special
configuration. If you can't or don't want to use gzip as the file 
compression tool, either pre-uncompress your files before parsing with
this module or modify B<Bio::Root::Utilities.pm> to your liking.

Blast objects can be generated either by direct instantiation as in:

 use Bio::Tools::Blast;		 
 $blast = new Bio::Tools::Blast (%parameters);

=head2 Using the Static $Blast Object

 use Bio::Tools::Blast qw(:obj);		 

This exports the static $Blast object into your namespace. "Static"
refers to the fact that it has class scope and there is one of these
created when you use this module. The static $Blast object is
basically an empty object that is provided for convenience and is also
used for various internal chores.

It is exported by this module and can be used for
parsing and running reports as well as HTML-formatting without having
to first create an empty Blast object.

Using the static $Blast object for parsing a STDIN stream of Blast reports:

    use Bio::Tools::Blast qw(:obj);

    sub process_blast {
	my $blastObj = shift;
	print $blastObj->table();
	$blastObj->destroy;
    }

    $Blast->parse( -parse     => 1,
		   -signif    => '1e-10',
		   -exec_func => \&process_blast,
		   );

Then pipe a stream of Blast reports into your script via STDIN.  For
each Blast report extracted from the input stream, the parser will
generate a new Blast object and pass it to the function specified by
C<-exec_func>.  The L<destroy>() call tells Perl to free the memory
associated with the object, important if you are crunching through
many reports. This method is inherited from B<Bio::Root::Object.pm>
(see L<Links to related modules>). See L<parse>() for a full
description of parameters and L<DEMO SCRIPTS> for additional examples.


=head2 Running Blasts

To run a Blast, create a new Blast object with a C<-run =E<gt>
\%runParams> parameter.  Remote Blasts are performed by including a
C<-method =E<gt> 'remote'> parameter; local Blasts are performed by
including a C<-method =E<gt> 'local'> parameter.  See L<Running Blast
reports> as well as the L<DEMO SCRIPTS> for examples.  Note that
running local Blasts is not yet supported, see below.

Note that the C<-seqs =E<gt> [ $seqs ]> run parameter must contain a
reference to an array of B<Bio::Seq.pm> objects (L<Links to related
modules>). Encapsulating the sequence in an object makes sequence
information much easier to handle as it can be supplied in a variety
of formats. Bio::Seq.pm is included with this distribution
(L<INSTALLATION>).

Remote Blasts are implemented using the
B<Bio::Tools::Blast::Run::Webblast.pm> module.  Local Blasts require
that you customize the B<Bio::Tools::Blast::Run::LocalBlast.pm>
module.  The version of LocalBlast.pm included with this distribution
provides the basic framework for running local Blasts. 
See L<Links to related modules>.

=head2 Significance screening

A C<-signif> parameter can be used to screen out all hits with
P-values (or Expect values) above a certain cutoff. For example, to
exclude all hits with Expect values above 1.0e-10: C<-signif =E<gt>
1e-10>. Providing a C<-signif> cutoff can speed up processing
tremendously, since only a small fraction of the report need be
parsed. This is because the C<-signif> value is used to screen hits
based on the data in the "Description" section of the Blast report:

For NCBI BLAST2 reports:

                                                                     Score     E
  Sequences producing significant alignments:                        (bits)  Value
 
  sp|P31376|YAB1_YEAST  HYPOTHETICAL 74.1 KD PROTEIN IN CYS3-MDM10...   957  0.0


For BLAST1 or WashU-BLAST2 reports:

                                                                       Smallest
                                                                         Sum
                                                                High  Probability
  Sequences producing High-scoring Segment Pairs:              Score  P(N)      N
 
  PDB:3PRK_E Proteinase K complexed with inhibitor ...........   504  1.8e-50   1


Thus, the C<-signif> parameter will screen based on Expect values for
BLAST2 reports and based on P-values for BLAST1/WashU-BLAST2 reports.

To screen based on other criteria, you can supply a C<-filt_func>
parameter containing a function reference that takes a
B<Bio::Tools::Sbjct.pm> object as an argument and returns a boolean,
true if the hit is to be screened out. See example below for
L<Screening hits using arbitrary criteria>.


=head2 Get the best hit.

     $hit = $blastObj->hit;  

A "hit" is contained by a B<Bio::Tools::Blast::Sbjct.pm> object. 
 
 
=head2 Get the P-value or Expect value of the most significant hit.

     $p = $blastObj->lowest_p;      
     $e = $blastObj->lowest_expect; 

Alternatively:

     $p = $blastObj->hit->p;      
     $e = $blastObj->hit->expect; 

Note that P-values are not reported in NCBI Blast2 reports.
 
 
=head2 Iterate through all the hits

     foreach $hit ($blastObj->hits) {
	 printf "%s\t %.1e\t %d\t %.2f\t %d\n", 
	                  $hit->name, $hit->expect, $hit->num_hsps,
                          $hit->frac_identical, $hit->gaps;
     }

Refer to the documentation for B<Bio::Tools::Blast::Sbjct.pm> 
for other ways to work with hit objects (L<Links to related modules>).

=head2 Screening hits using arbitrary criteria

    sub filter { $hit=shift;
		 return ($hit->gaps == 0 and 
			 $hit->frac_conserved > 0.5); }
 
     $blastObj = Bio::Tools::Blast->new( -file      => '/tmp/blast.out',
					 -parse     => 1,  
					 -filt_func => \&filter );

While the Blast object is parsing the report, each hit checked by calling
&filter($hit). All hits that generate false return values from &filter
are screened out and will not be added to the Blast object.
Note that the Blast object will normally stop parsing the report after
the first non-significant hit or the first hit that does not pass the
filter function. To force the Blast object to check all hits,
include a C<-check_all_hits =E<gt> 1>  parameter.
Refer to the documentation for B<Bio::Tools::Blast::Sbjct.pm> 
for other ways to work with hit objects.

=over 4

=item Hit start, end coordinates.

      print $sbjct->start('query');
      print $sbjct->end('sbjct');

In array context, you can get information for both query and sbjct with one call:

      ($qstart, $sstart) = $sbjct->start();
      ($qend, $send)     = $sbjct->end();

For important information regarding coordinate information, see 
the L<HSP start, end, and strand> section below.
Also check out documentation for the start and end methods in B<Bio::Tools::Blast::Sbjct.pm>,
which explains what happens if there is more than one HSP.

=back

=head2 Working with HSPs

=over 4 

=item Iterate through all the HSPs of every hit

     foreach $hit ($blastObj->hits) {
	 foreach $hsp ($hit->hsps) {
	 printf "%.1e\t %d\t %.1f\t %.2f\t %.2f\t %d\t %d\n", 
	                  $hsp->expect, $hsp->score, $hsp->bits,
                          $hsp->frac_identical, $hsp->frac_conserved, 
	                  $hsp->gaps('query'), $hsp->gaps('sbjct');
     }

Refer to the documentation for B<Bio::Tools::Blast::HSP.pm> 
for other ways to work with hit objects (L<Links to related modules>).

=back

=over 4

=item Extract HSP sequence data as strings or sequence objects

Get the first HSP of the first hit and the sequences
of the query and sbjct as strings.

      $hsp = $blast_obj->hit->hsp;  
      $query_seq = $hsp->seq_str('query');
      $hsp_seq = $hsp->seq_str('sbjct');

Get the indices of identical and conserved positions in the HSP query seq.

      @query_iden_indices = $hsp->seq_inds('query', 'identical');
      @query_cons_indices = $hsp->seq_inds('query', 'conserved');

Similarly for the sbjct sequence.

      @sbjct_iden_indices = $hsp->seq_inds('sbjct', 'identical');
      @sbjct_cons_indices = $hsp->seq_inds('sbjct', 'conserved');
      	    
      print "Query in Fasta format:\n", $hsp->seq('query')->layout('fasta');
      print "Sbjct in Fasta format:\n", $hsp->seq('sbjct')->layout('fasta');

See the B<Bio::Seq.pm> package for more information about using these sequence objects
(L<Links to related modules>).

=back

=over 4 

=item Create sequence alignment objects using HSP sequences

      $aln = $hsp->get_aln;
      print " consensus:\n", $aln->consensus();
      print $hsp->get_aln->layout('fasta');
 
      $ENV{READSEQ_DIR} = '/home/users/sac/bin/solaris';
      $ENV{READSEQ} = 'readseq';
      print $hsp->get_aln->layout('msf');

MSF formated layout requires Don Gilbert's ReadSeq program (not included). 
See the B<Bio::UnivAln.pm> for more information about using these alignment objects 
(L<Links to related modules>)'.

=back

=over 4

=item HSP start, end, and strand.

To facilitate HSP processing, endpoint data for each HSP sequence are 
normalized so that B<start is always less than end>. This affects TBLASTN 
and TBLASTX HSPs on the reverse complement or "Minus" strand.

Some examples of obtaining start, end coordinates for HSP objects:

      print $hsp->start('query');
      print $hsp->end('sbjct');
      ($qstart, $sstart) = $hsp->start();
      ($qend, $send) = $hsp->end();

Strandedness of the HSP can be assessed using the strand() method 
on the HSP object:

      print $hsp->strand('query');
      print $hsp->strand('sbjct');

These will return 'Minus' or 'Plus'.
Or, to get strand information for both query and sbjct with a single call:

      ($qstrand, $sstrand) = $hsp->strand();

=back

=head2 Report Generation

=over 4

=item Generate a tab-delimited table of all results.

     print $blastObj->table;       
     print $blastObj->table(0);   # don't include hit descriptions.
     print $blastObj->table_tiled; 

The L<table>() method returns data for each B<HSP> of each hit listed one per
line. The L<table_tiled>() method returns data for each B<hit, i.e., Sbjct>
listed one per line; data from multiple HSPs are combined after tiling to
reduce overlaps. See B<Bio::Tools::Blast::Sbjct.pm> for more information about
HSP tiling.  These methods generate stereotypical, tab-delimited data for each
hit of the Blast report. The output is suitable for importation into
spreadsheets or database tables. Feel free to roll your own table function if
you need a custom table.

For either table method, descriptions of each hit can be included if a 
single, true argument is supplied (e.g., $blastObj->table(1)). The description 
will be added as the last field. This will significantly increase the size of
the table. Labels for the table columns can be obtained with L<table_labels>()
and L<table_labels_tiled>().

=back

=over 4

=item Print a summary of the Blast report

     $blastObj->display();     
     $blastObj->display(-show=>'hits');

L<display>() prints various statistics extracted from the Blast report
such as database name, database size, matrix used, etc. The
C<display(-show=E<gt>'hits')> call prints a non-tab-delimited table
attempting to line the data up into more readable columns. The output
generated is similar to L<table_tiled>().

=back

=over 4

=item HTML-format an existing report

     use Bio::Tools::Blast qw(:obj);
 
     # Going straight from a non HTML report file to HTML output using 
     # the static $Blast object exported by Bio::Tools::Blast.pm
     $Blast->to_html(-file   => '/usr/people/me/blast.output.txt',
		     -header => qq|<H1>BLASTP Results</H1><A HREF="home.html">Home</A>|
		     );
 
     # You can also use a specific Blast object created previously.
     $blastObj->to_html;

L<to_html>() will send HTML output, line-by-line, directly to STDOUT
unless an C<-out =E<gt> array_ref> parameter is supplied (e.g., C<-out
=E<gt> \@array>), in which case the HTML will be stored in @array, one
line per array element.  The direct outputting permits faster response
time since Blast reports can be huge. The -header tag can contain a
string containing any HTML that you want to appear at the top of the
Blast report.

=back

=head1 DEMO SCRIPTS

Sample Scripts are included in the central bioperl distribution in the
'examples/blast/' directory (see L<INSTALLATION>). These are also
available at the following URLs (but it would be safer to use the
scripts included with the distribution).

=head2 Handy library for working with Bio::Tools::Blast.pm

   http://bio.perl.org/Core/Examples/blast/blast_config.pl

=head2 Parsing Blast reports one at a time.

   http://bio.perl.org/Core/Examples/blast/parse_blast.pl
   http://bio.perl.org/Core/Examples/blast/parse_blast2.pl
   http://bio.perl.org/Core/Examples/blast/parse_positions.pl

=head2 Parsing sets of Blast reports.

   http://bio.perl.org/Core/Examples/blast/parse_blast.pl
   http://bio.perl.org/Core/Examples/blast/parse_multi.pl

   B<Warning:> See note about L<Memory Usage Issues>.

=head2 Running Blast analyses one at a time.

   http://bio.perl.org/Core/Examples/blast/run_blast_remote.pl

=head2 Running Blast analyses given a set of sequences.

   http://bio.perl.org/Core/Examples/blast/blast_seq.pl

=head2 HTML-formatting Blast reports.

   http://bio.perl.org/Core/Examples/blast/html.pl



=head1 TECHNICAL DETAILS

=head2 Blast Modes

A BLAST object may be created using one of three different modes as
defined by the B<Bio::Tools::SeqAnal.pm> package (See L<Links to
related modules>):

 -- parse - Load a BLAST report and parse it, storing parsed data in
    Blast.pm object.
 -- run      - Run the BLAST program to generate a new report. 
 -- read     - Load a BLAST report into the Blast object without parsing.


B<Run mode support has recently been added>.  The module
B<Bio::Tools::Blast::Run::Webblast.pm> is an modularized adaptation of
the webblast script by Alex Dong Li:

   http://www.genet.sickkids.on.ca/bioinfo_resources/software.html#webblast

for running remote Blast analyses and saving the results locally.  Run
mode can be combined with a parse mode to generate a Blast report and
then build the Blast object from the parsed results of this report
(see L<run>() and L<SYNOPSIS>).

In read mode, the BLAST report is read in by the Blast object but is
not parsed.  This could be used to internalize a Blast report but not
parse it for results (e.g., generating HTML formatted output).



=head2 Significant Hits

This module permits the screening of hits on the basis of
user-specified criteria for significance. Currently, Blast reports can
be screened based on:

   CRITERIA                            PARAMETER       VALUE
   ----------------------------------  ---------      ----------------
  1) the best Expect (or P) value      -signif        float or sci-notation
  2) the length of the query sequence  -min_length    integer
  3) arbitrary criteria                -filt_func     function reference

The parameters are used for construction of the BLAST object or when
running the L<parse>() method on the static $Blast object.  The
-SIGNIF value represents the number listed in the description section
at the top of the Blast report. For Blast2, this is an Expect value,
for Blast1 and WashU-Blast2, this is a P-value.  The idea behind the
C<-filt_func> parameter is that the hit has to pass through a filter
to be considered significant. Refer to the documentation for
B<Bio::Tools::Blast::Sbjct.pm> for ways to work with hit objects.

Using a C<-signif> parameter allows for the following:

=over 2

=item Faster parsing.

Each hit can be screened by examination of the description line alone
without fully parsing the HSP alignment section.

=item Flexibility.

The C<-signif> tag provides a more semantic-free way to specify the
value to be used as a basis for screening hits. Thus, C<-signif> can
be used for screening Blast1 or Blast2 reports. It is up to the user
to understand whether C<-signif> represents a P-value or an Expect
value.

=back

Any hit not meeting the significance criteria will not be added to the
"hit list" of the BLAST object. Also, a BLAST object without any hits
meeting the significance criteria will throw an exception during
object construction (a fatal event).


=head2 Statistical Parameters

There are numerous parameters which define the behavior of the BLAST
program and which are useful for interpreting the search
results. These parameters are extracted from the Blast report:

  filter  --  for masking out low-complexity sequences or short repeats
  matrix  --  name of the substitution scoring matrix (e.g., BLOSUM62)
  E       --  Expect filter (screens out frequent scores)
  S       --  Cutoff score for segment pairs
  W       --  Word length
  T       --  Threshold score for word pairs
  Lambda, --  Karlin-Altschul "sum" statistical parameters dependent on  
   K, H        sequence composition.
  G       --  Gap creation penalty.
  E       --  Gap extension penalty. 

These parameters are not always needed. Extraction may be turned off
explicitly by including a C<-stats =E<gt> 0> parameter during object
construction.  Support for all statistical parameters is not complete.

For more about the meaning of parameters, check out the NCBI URLs given above.


=head2 Module Organization

The modules that comprise this Bioperl Blast distribution are location in the
Bio:: hierarchy as shown in the diagram below.

                            Bio/
                             |
               +--------------------------+
               |                          |
          Bio::Tools                  Bio::Root
               |                          |
    +----------------------+           Object.pm
    |          |           |
 SeqAnal.pm  Blast.pm    Blast/
                           |
            +---------+---------+------------+
            |         |         |            |
          Sbjct.pm   HSP.pm   HTML.pm       Run/
                                             |
                                       +------------+
                                       |            |
                                  Webblast.pm   LocalBlast.pm 


Bio::Tools::Blast.pm is a concrete class that inherits from
B<Bio::Tools::SeqAnal.pm> and relies on other modules for parsing and
managing BLAST data.  Worth mentioning about this hierarchy is the
lack of a "Parse.pm" module.  Since parsing is considered central to
the purpose of the Bioperl Blast module (and Bioperl in general), it
seems somewhat unnatural to segregate out all parsing code. This
segregation could also lead to inefficiencies and harder to maintain
code. I consider this issue still open for debate.

Bio::Tools::Blast.pm, B<Bio::Tools::Blast::Sbjct.pm>, and
B<Bio::Tools::Blast::HSP.pm> are mostly dedicated to parsing and all
can be used to instantiate objects.  Blast.pm is the main "command and
control" module, inheriting some basic behaviors from SeqAnal.pm
(things that are not specific to Blast I<per se>).

B<Bio::Tools::Blast::HTML.pm> contains functions dedicated to
generating HTML-formatted Blast reports and does not generate objects.

=head2 Running Blasts: Details

B<Bio::Tools::Blast::Run::Webblast.pm> contains a set of functions for
running Blast analyses at a remote server and also does not
instantiate objects.  It uses a helper script called postclient.pl,
located in the Run directory.  The proposed LocalBlast.pm module would
be used for running Blast reports on local machines and thus would be
customizable for different sites. It would operate in a parallel
fashion to Webblast.pm (i.e., being a collection of functions, taking
in sequence objects or files, returning result files).

The Run modules are considered experimental. In particular,
Webblast.pm catures an HTML-formatted version of the Blast report from
the NCBI server and strips out the HTML in preparation for parsing. A
more direct approach would be to capture the Blast results directly
from the server using an interface to the NCBI toolkit.  This approach
was recently proposed on the Bioperl mailing list:
http://www.uni-bielefeld.de/mailinglists/BCD/vsns-bcd-perl/9805/0000.html


=head2 Memory Usage Issues

Parsing large numbers of Blast reports (a few thousand or so) with
Bio::Tools::Blast.pm may lead to unacceptable memory usage situations.
This is somewhat dependent of the size and complexity of the reports.

While this problem is under investigation, here are some workarounds
that fix the memory usage problem:

=over 4

=item 1 Don't specify a -signif criterion when calling L<parse>().

The C<-signif> value is used for imposing a upper limit to the expect- or 
P-value for Blast hits to be parsed. For reasons that are still under 
investigation, specifying a value for C<-signif> in the L<parse>() 
method prevents Blast objects from being fully
garbage collected. When using the B<parse_blast.pl> or B<parse_multi.pl> 
scripts in C<examples/blast/> of the bioperl distribution), don't supply 
a C<-signif> command-line parameter.
  

=item 2 If you want to impose a -signif criterion, put it inside a
-filt_func.

For the L<parse>() method, a -signif => 1e-5 parameter is equivalent
to using a filter function parameter of

 -filt_func => sub { my $hit = shift; return $hit->signif <= 1e-5; }

Using the B<examples/blast/parse_multi.pl> script, you can supply a
command-line argument of

 -filt_func '$hit->signif <= 1e-5'

For more information, see L<parse>() and the section L<Screening hits
using arbitrary criteria>.

=back

=head1 TODO

=over 4

=item * Develop a functional, prototype Bio::Tools::Blast::Run::LocalBlast.pm module.

=item * Add support for PSI-BLAST and PHI-BLAST

=item * Parse histogram of expectations and retrieve gif image in
Blast report (if present).

=item * Further investigate memory leak that occurs when parsing Blast
streams whe supplying a -signif parameter to L<parse>().

=item * Access Blast results directly from the NCBI server using a
Perl interface to the NCBI toolkit or XML formated Blast reports (when
available).

=item * Further exploit Bio::UnivAln.pm and multiple-sequence
alignment programs using HSP sequence data. Some of this may best go
into a separate, dedicated module or script as opposed to burdening
Blast.pm, Sbjct.pm, and HSP.pm with additional functionality that is
not always required.

=item * Add an example script for parsing Blast reports containing
HTML formatting.


=back


=head1 VERSION

Bio::Tools::Blast.pm, 0.09


=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

    bioperl-l@bioperl.org             - General discussion
    http://bio.perl.org/MailList.html - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via email
or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           


=head1 AUTHOR

Steve A. Chervitz, sac@genome.stanford.edu

See the L<FEEDBACK> section for where to send bug reports and comments.

=head1 ACKNOWLEDGEMENTS

This module was developed under the auspices of the Saccharomyces Genome
Database:
    http://genome-www.stanford.edu/Saccharomyces

Other contributors include: Alex Dong Li (webblast), Chris Dagdigian
(Seq.pm), Steve Brenner (Seq.pm), Georg Fuellen (Seq.pm, UnivAln.pm),
and untold others who have offered comments (noted in the
Bio/Tools/Blast/CHANGES file of the distribution).

=head1 COPYRIGHT

Copyright (c) 1996-98 Steve A. Chervitz. All Rights Reserved.  This
module is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=head1 SEE ALSO

 Bio::Tools::SeqAnal.pm                  - Sequence analysis object base class.
 Bio::Tools::Blast::Sbjct.pm             - Blast hit object.
 Bio::Tools::Blast::HSP.pm               - Blast HSP object.
 Bio::Tools::Blast::HTML.pm              - Blast HTML-formating utility class.
 Bio::Tools::Blast::Run::Webblast.pm     - Utility module for running Blasts remotely.
 Bio::Tools::Blast::Run::LocalBlast.pm   - Utility module for running Blasts locally.
 Bio::Seq.pm                             - Biosequence object  
 Bio::UnivAln.pm                         - Biosequence alignment object.
 Bio::Root::Object.pm                    - Proposed base class for all Bioperl objects.

=head2 Links to related modules

 Bio::Tools::SeqAnal.pm      
      http://bio.perl.org/Core/POD/Bio/Tools/SeqAnal.html

 Bio::Tools::Blast::Sbjct.pm 
      http://bio.perl.org/Core/POD/Bio/Tools/Blast/Sbjct.html

 Bio::Tools::Blast::HSP.pm   
      http://bio.perl.org/Core/POD/Bio/Tools/Blast/HSP.html

 Bio::Tools::Blast::HTML.pm       
      http://bio.perl.org/Core/POD/Bio/Tools/Blast/HTML.html

 Bio::Tools::Blast::Run::Webblast.pm 
      http://bio.perl.org/Core/POD/Bio/Tools/Blast/Run/Webblast.html

 Bio::Tools::Blast::Run::LocalBlast.pm 
      http://bio.perl.org/Core/POD/Bio/Tools/Blast/Run/LocalBlast.html

 Bio::Seq.pm              
      http://bio.perl.org/Core/POD/Seq.html

 Bio::UnivAln.pm             
      http://bio.perl.org/Projects/SeqAlign/
      Europe:  http://www.techfak.uni-bielefeld.de/bcd/Perl/Bio/#univaln

 Bio::Root::Object.pm        
      http://bio.perl.org/Core/POD/Root/Object.html


 http://bio.perl.org/Projects/modules.html  - Online module documentation
 http://bio.perl.org/Projects/Blast/        - Bioperl Blast Project     
 http://bio.perl.org/                       - Bioperl Project Homepage


L<References & Information about the BLAST program>.

=head1 KNOWN BUGS

There is a memory leak that occurs when parsing parsing streams
containing large numbers of Blast reports (a few thousand or so) and
specifying a -signif parameter to the L<parse>() method. For a
workaround, see L<Memory Usage Issues>.

Not sharing statistical parameters between different Blast objects
when parsing a multi-report stream has not been completely tested and
may be a little buggy.

Documentation inconsistencies or inaccuracies may exist since this
module underwend a fair bit of re-working going from 0.75 to 0.80
(corresponds to versions 0.04.4 to 0.05 of the bioperl distribution).

=cut



#
##
###
#### END of main POD documentation.
###
##
#


=head1 APPENDIX

Methods beginning with a leading underscore are considered private and
are intended for internal use by this module. They are B<not>
considered part of the public interface and are described here for
documentation purposes only.

=cut


##############################################################################
##                          CONSTRUCTOR                                     ##
##############################################################################

## The Blast.pm object relies on the the superclass constructor:
## Bio::Tools::SeqAnal::_initialize(). See that module for details.


#-------------
sub destroy {
#-------------
    my $self=shift; 
    $DEBUG==2 && print STDERR "DESTROYING $self ${\$self->name}";
    if($self->{'_hits'}) {
	foreach($self->hits) { 
	    $_->destroy; 
	    undef $_;
	}
	undef $self->{'_hits'};
        #$self->{'_hits'}->remove_all;  ## When and if this member becomes a vector.
    }

    $self->SUPER::destroy;
}

#####################################################################################
##                                  ACCESSORS                                      ##
#####################################################################################

=head2 run

 Usage     : $object->run( %named_parameters )
 Purpose   : Run a local or remote Blast analysis on one or more sequences.
 Returns   : String containing name of Blast output file if a single Blast 
           : is run.
           :  -- OR --
           : List of Blast objects if multiple Blasts are being run as a group.
 Argument  : Named parameters:  (PARAMETER TAGS CAN BE UPPER OR LOWER CASE).
           :    -METHOD  => 'local' or 'remote' (default = remote),
           :    -PARSE   => boolean, (true if the results are to be parsed after the run)
           :    -STRICT  => boolean, the strict mode to use for the resulting Blast objects.
           :  ADDITIONAL PARAMETERS:
           :      See methods _run_remote() and _run_local() for required
           :      parameters necessary for running the blast report.
 Throws    : Exception if no Blast output file was obtained.
 Comments  : This method is called automatically during construction of a 
           : Blast.pm object when run parameters are sent to the constructor:
           :  $blastObj = new Bio::Tools::Blast (-RUN =>\%runParam,
	   :					 %parseParam );
           :
           : The specific run methods (local or remote) called by run() 
           : must return a list containing  the file name(s) with the Blast output. 
           :
           : The run() method can perform single or multiple Blast runs
           : (analogous to the way parse() works) depending on how many 
           : sequences are submitted. However, the running of multiple
           : Blasts is probably better handled at the script level. See notes in
           : the "TODO" section below.
           :
           : As for what to do with the Blast result file, that decision is 
           : left for the user who can direct the Blast object to delete, compress,
           : or leave it alone.
           :
           : This method does not worry about load balancing, which
           : is probably best handled at the server level.
           :
 TODO:     : Support for running+parsing multiple Blast analyses with a 
           : single run() call is incomplete. One can generate multiple 
           : reports by placing more than one sequence object in the -seqs  
           : reference parameter. This saves some overhead in the code 
           : that executes the Blasts since all options are configured once.
           : (This is analogous to parsing using the static $Blast object 
           : see parse() and _parse_stream()).
           :
           : The trouble is that Blast objects for all runs are constructed,
           : parsed (if necessary), and then returned as a group
           : This can require lots of memory when run+parsing many Blasts
           : but should be fine if you just want to run a bunch Blasts.
           :
           : For now, when running+parsing Blasts, stick to running one 
           : Blast at a time, building the Blast object with the results 
           : of that report, and processing as necessary.
           : 
           : Support for running PSI-Blast is not complete.

See Also:  L<_run_remote>(), L<_run_local>(), L<parse>()

=cut

#---------
sub run {
#---------
    my ($self, %param) = @_;
    my($method, $parse, $strict) = 
	$self->_rearrange([qw(METHOD PARSE STRICT)], %param);

    $strict = $self->strict($strict) if $strict;

    my (@files);
    if($method =~ /loc/i) {
	@files = $self->_run_local(%param);
	
    } else {
	@files = $self->_run_remote(%param);
    }
    
    $self->throw("Run Blast failed: no Blast output created.") if !@files;

    if(scalar(@files) == 1) {
	# If there was just one Blast output file, prepare to incorporate it
	# into the current Blast object. run() is called before parse() in the
	# SeqAnal.pm constructor.
	if($files[0] ne 'email') {
	    $self->file($files[0]);
	} else { 
	    # Can't do anything with the report.
	    $self->throw("Blast report to be sent via e-mail.");
	} 

    } else {
	# If there are multiple report files, build individual Blast objects foreach.
	# In this situation, the static $Blast object is being used to run
	# a set of related Blasts, similar to the way parse() can be used.
	# This strategy is not optimal since all reports are generated first
	# before any are parsed. 
	# Untested.

	my(@objs);
	foreach(@files) {
	    push @objs, new Bio::Tools::Blast(-FILE   => $_,
					      -PARSE  => $parse || 0,
					      -STRICT => $strict,
					      );
	}
	return @objs;
    }
}


=head2 _run_remote

 Usage     : n/a; internal method called by run()
           : $object->_run_remote( %named_parameters )
 Purpose   : Run Blast on a remote server.
 Argument  : Named parameters:
           :   See documentation for function &blast_remote in
           :   Bio::Tools::Blast::Run::Webblast.pm for description 
           :   of parameters.
 Comments  : This method requires the Bio::Tools::Blast::Run::Webblast.pm
           : which conforms to this minimal API:
           :    * export a method called &blast_remote that accepts a 
           :      Bio::Tools::Blast.pm object + named parameters
           :      (specified in the Webblast.pm module).
           :    * return a list of names of files containing the raw Blast reports.
           :      (When building a Blast object, this list would contain a 
           :       single file from which the Blast object is to be constructed).

See Also   : L<run>(), L<_run_local>(), B<Bio::Tools::Blast::Run::Webblast.pm::blast_remote()>, L<Links to related modules>

=cut

#----------------
sub _run_remote {
#----------------
    my ($self, %param) = @_;

    require Bio::Tools::Blast::Run::Webblast; 
    Bio::Tools::Blast::Run::Webblast->import(qw(&blast_remote));

    &blast_remote($self, %param);
}



=head2 _run_local

 Usage     : n/a; internal method called by run()
           : $object->_run_local(%named_parameters)
 Purpose   : Run Blast on a local machine.
 Argument  : Named parameters:
           :   See documentation for function &blast_local in
           :   Bio::Tools::Blast::Run::LocalBlast.pm for description 
           :   of parameters.
 Comments  : This method requires the Bio::Tools::Blast::Run::LocalBlast.pm
           : module which should be customized for your site. This module would 
           : contain all the commands, paths, environment variables, and other 
           : data necessary to run Blast commands on a local machine, but should
           : not contain any semantics for specific query sequences.
           :
           : LocalBlast.pm should also conform to this minimal API:
           :    * export a method called &blast_local that accepts a 
           :       Bio::Tools::Blast.pm object + named parameters
           :      (specified in the LocalBlast.pm module).
           :    * return a list of names of files containing the raw Blast reports.
           :      (When building a Blast object, this list would contain a 
           :       single file from which the Blast object is to be constructed).

See Also   : L<run>(), L<_run_remote>(), B<Bio::Tools::Blast::Run::LocalBlast::blast_local()>, L<Links to related modules>

=cut

#--------------
sub _run_local {
#--------------
    my ($self, %param) = @_;
    
    require Bio::Tools::Blast::Run::Webblast; 
    Bio::Tools::Blast::Run::Webblast->import(qw(&blast_local));

    &blast_local($self, %param);
}


=head2 db_remote

 Usage     : @dbs = $Blast->db_remote( [seq_type] );
 Purpose   : Get a list of available sequence databases for remote Blast analysis.
 Returns   : Array of strings 
 Argument  : seq_type = 'p' or 'n' 
           :  'p' = Gets databases for peptide searches  (default)
           :  'n' = Gets databases for nucleotide searches 
 Throws    : n/a
 Comments  : Peptide databases are a subset of the nucleotide databases.
           : It is convenient to call this method on the static $Blast object
           : as shown in Usage.

See Also   : L<db_local>()

=cut

#----------------
sub db_remote {
#----------------
    my ($self, $type) = @_;
    $type ||= 'p';

    require Bio::Tools::Blast::Run::Webblast; 
    Bio::Tools::Blast::Run::Webblast->import(qw(@Blast_dbp_remote
						 @Blast_dbn_remote));

    # We shouldn't have to fully qualify the Blast_dbX_remote arrays. Hm.

    my(@dbs);
    if( $type =~ /^p|amino/i) {
	@dbs = @Bio::Tools::Blast::Run::Webblast::Blast_dbp_remote;
    } else {
	@dbs = @Bio::Tools::Blast::Run::Webblast::Blast_dbn_remote;
    }
    @dbs;
}



=head2 db_local

 Usage     : @dbs = $Blast->db_local( [seq_type] );
 Purpose   : Get a list of available sequence databases for local Blast analysis.
 Returns   : Array of strings 
 Argument  : seq_type = 'p' or 'n'
           :  'p' = Gets databases for peptide searches  (default)
           :  'n' = Gets databases for nucleotide searches 
 Throws    : n/a
 Comments  : Peptide databases are a subset of the nucleotide databases.
           : It is convenient to call this method on the static $Blast object.
             as shown in Usage.

See Also   : L<db_remote>()

=cut

#----------------
sub db_local {
#----------------
    my ($self, $type) = @_;
    $type ||= 'p';

    require Bio::Tools::Blast::Run::LocalBlast; 
    Bio::Tools::Blast::Run::LocalBlast->import(qw(@Blast_dbp_local
						  @Blast_dbn_local));

    # We shouldn't have to fully qualify the Blast_dbX_local arrays. Hm.

    my(@dbs);
    if( $type =~ /^p|amino/i) {
	@dbs = @Bio::Tools::Blast::Run::LocalBlast::Blast_dbp_local;
    } else {
	@dbs = @Bio::Tools::Blast::Run::LocalBlast::Blast_dbn_local;
    }
    @dbs;
}


=head2 parse

 Usage     : $blast_object->parse( %named_parameters )
 Purpose   : Parse a Blast report from a file or STDIN.
           :   * Parses a raw BLAST data, populating Blast object with report data.
           :   * Sets the significance cutoff.
           :   * Extracts statistical parameters about the BLAST run.
           :   * Handles both single files and streams containing multiple reports.
 Returns   : integer (number of Blast reports parsed)
 Argument  : <named parameters>:  (PARAMETER TAGS CAN BE UPPER OR LOWER CASE).
	   : -FILE       => string (name of file containing raw Blast output. 
           :                         Optional. If a valid file is not supplied, 
	   :		             STDIN will be used).
	   : -SIGNIF     => number (float or scientific notation number to be used
	   :                         as a P- or Expect value cutoff; 
	   :			     default =  $DEFAULT_SIGNIF (999)).
	   : -FILT_FUNC  => func_ref (reference to a function to be used for 
           :                          filtering out hits based on arbitrary criteria. 
           :                          This function should take a
           :                          Bio::Tools::Blast::Sbjct.pm object as its first
           :                          argument and return a boolean value,
	   :                          true if the hit should be filtered out).
           :                          Sample filter function:
           :                          -FILT_FUNC => sub { $hit = shift;
	   :				                  $hit->gaps == 0; },
           : -CHECK_ALL_HITS => boolean (check all hits for significance against
           :                             significance criteria.  Default = false.
	   :			         If false, stops processing hits after the first
           :                             non-significant hit or the first hit that fails
           :                             the filt_func call. This speeds parsing, 
           :                             taking advantage of the fact that the hits 
           :                             are processed in the order they are ranked.)
           : -MIN_LEN     => integer (to be used as a minimum query sequence length
           :                          sequences below this length will not be processed).
           :                          default = no minimum length).
	   : -STATS       => boolean (collect stats for report: matrix, filters, etc.
           :                          default = false).
	   : -BEST        => boolean (only process the best hit of each report; 
           :                          default = false).
           : -OVERLAP     => integer (the amount of overlap to permit between 
           :                          adjacent HSPs when tiling HSPs, 
           :                          Default = $MAX_HSP_OVERLAP (2))
           :
           : PARAMETERS USED WHEN PARSING MULTI-REPORT STREAMS:
           : --------------------------------------------------
	   : -SHARE       => boolean (set this to true if all reports in stream
	   :			      share the same stats. Default = true)
           :                          Must be set to false when parsing both Blast1 and
           :                          Blast2 reports in the same run or if you need
           :                          statistical params for each report, Lambda, K, H).
	   : -STRICT      => boolean (use strict mode for all Blast objects created.
           :                          Increases sensitivity to errors. For single
           :                          Blasts, this is parameter is sent to new().)
           : -EXEC_FUNC   => func_ref (reference to a function for processing each
           :                           Blast object after it is parsed. Should accept a
           :                           Blast object as its sole argument. Return value
           :                           is ignored. If an -EXEC_FUNC parameter is supplied, 
           :                           the -SAVE_ARRAY parameter will be ignored.)
           : -SAVE_ARRAY  =>array_ref, (reference to an array for storing all
           :                            Blast objects as they are created. 
           :                            Experimental. Not recommended.)
           : -SIGNIF_FMT  => boolean   String of 'exp' or 'parts'. Sets the format 
           :                           for reporting P/Expect values. 'exp' reports 
           :                           only the exponent portion. 'parts' reports 
           :                           them as a 2 element list. See signif_fmt()..
           :
 Throws    : Exception if BLAST report contains a FATAL: error.
           : Propagates any exception thrown by read().
           : Propagates any exception thrown by called parsing methods.
 Comments  : This method can be called either directly using the static $Blast object
           : or indirectly (by Bio::Tools::SeqAnal.pm) during constuction of an 
           : individual Blast object.
           :
           : HTML-formatted reports can be parsed as well. No special flag is required
           : since it is detected automatically. The presence of HTML-formatting 
           : will result in slower performace, however, since it must be removed
           : prior to parsing. Parsing HTML-formatted reports is highly
           : error prone and is generally not recommended.
           :               
           : If one has an HTML report, do NOT remove the HTML from it by using the
           : "Save As" option of a web browser to save it as text. This renders the
           : report unparsable.
           : HTML-formatted reports can be parsed after running through the strip_html
           : function of Blast::HTML.pm as in:
           :    require Bio::Tools::Blast::HTML; 
           :    Bio::Tools::Blast::HTML->import(&strip_html); 
           :    &strip_html(\$data);  
           :    # where data contains full contents of an HTML-formatted report.
           : TODO: write a demo script that does this.

See Also   : L<_init_parse_params>(), L<_parse_blast_stream>(), L<overlap>(), L<signif_fmt>(), B<Bio::Root::Object::read()>, B<Bio::Tools::Blast::HTML.pm::strip_html()>, L<Links to related modules>

=cut

#---------
sub parse {
#---------
# $self might be the static $Blast object.
    my ($self, @param) = @_;

    my($signif, $filt_func, $min_len, $check_all, $overlap, $stats, 
       $share, $strict, $best, $signif_fmt, $no_aligns) = 
	$self->_rearrange([qw(SIGNIF FILT_FUNC MIN_LEN CHECK_ALL_HITS 
			      OVERLAP STATS SHARE STRICT 
			      BEST EXPONENT NO_ALIGNS )], @param);

    ## Initialize the static Blast object with parameters that 
    ## apply to all Blast objects within a parsing session.

    &_init_parse_params($share, $filt_func, $check_all,
			$signif, $min_len, $strict,
			$best, $signif_fmt, $stats, $no_aligns
		       );

    my $count = $self->_parse_blast_stream(@param);
	
#    print STDERR "\nDONE PARSING STREAM.\n";

    if($Blast->{'_blast_errs'}) {
      my @errs = @{$Blast->{'_blast_errs'}};
      printf STDERR "\n*** %d BLAST REPORTS HAD FATAL ERRORS:\n", scalar(@errs);
      foreach(@errs) { print STDERR "$_\n"; }
      @{$Blast->{'_blast_errs'}} = ();
    }

    return $count;
}

=head2 _init_parse_params

 Title   : _init_parse_params
 Usage   : n/a; called automatically by parse()
 Purpose : Initializes parameters used during parsing of Blast reports.
         : This is a static method used by the $Blast object.
         : Calls _set_signif().
 Example :
 Returns : n/a
 Args    : Args extracted by parse().

See Also: L<parse>(), L<_set_signif>()

=cut

#----------------------
sub _init_parse_params {
#----------------------
    my ($share, $filt_func, $check_all, 
	$signif, $min_len, $strict,
	$best, $signif_fmt, $stats, $no_aligns) = @_;

    ## Default is to share stats.
    $Blast->{'_share'}  = defined($share) ? $share : 1;
    $Blast->{'_filt_func'} = $filt_func || 0;  
    $Blast->{'_check_all'} = $check_all || 0; 
    $Blast->{'_signif_fmt'} ||= $signif_fmt || ''; 
    $Blast->{'_no_aligns'} = $no_aligns || 0; 

    &_set_signif($signif, $min_len, $filt_func);
    $Blast->strict($strict) if defined $strict;  
    $Blast->best($best) if $best;
    $Blast->{'_blast_count'} = 0;

    ## If $stats is false, miscellaneous statistical and other parameters
    ## are NOT extracted from the Blast report (e.g., matrix name, filter used, etc.).
    ## This can speed processing when crunching tons of Blast reports.
    ## Default is to NOT get stats.
    $Blast->{'_get_stats'} = defined($stats) ? $stats : 0;  

    # Clear any errors from previous parse.
    undef $Blast->{'_blast_errs'};
}



=head2 _set_signif

 Usage     : n/a; called automatically by _init_parse_params()
           : This is now a "static" method used only by $Blast.
           : _set_signif($signif, $min_len, $filt_func); 
 Purpose   : Sets significance criteria for the BLAST object.
 Argument  : Obligatory three arguments:
           :   $signif = float or sci-notation number or undef
           :   $min_len = integer or undef
           :   $filt_func = function reference or undef
           :  
           :   If $signif is undefined, a default value is set 
           :   (see $DEFAULT_SIGNIF; min_length = not set).
 Throws    : Exception if significance value is defined but appears
           :   out of range or invalid.
           : Exception if $filt_func if defined and is not a func ref.
 Comments  : The significance of a BLAST report can be based on
           : the P (or Expect) value and/or the length of the query sequence.
           : P (or Expect) values GREATER than '_significance' are not significant.
           : Query sequence lengths LESS than '_min_length' are not significant.
           :
           : Hits can also be screened using arbitrary significance criteria 
           : as discussed in the parse() method.
           :
           : If no $signif is defined, the '_significance' level is set to 
           : $Bio::Tools::Blast::DEFAULT_SIGNIF (999).

See Also   : L<signif>(), L<min_length>(), L<_init_parse_params>(), parse>()

=cut

#-----------------
sub _set_signif {
#-----------------
    my( $sig, $len, $func ) = @_;
    
    if(defined $sig) {
	$Blast->{'_confirm_significance'} = 1;
	if( $sig =~ /[^\d.e-]/ or $sig <= 0) { 
	    $Blast->throw("Invalid significance value: $sig", 
			 "Must be greater than zero.");
	} 
	$Blast->{'_significance'} = $sig;
    } else {
	$Blast->{'_significance'}   = $DEFAULT_SIGNIF;
	$Blast->{'_check_all'}      = 1 if not $Blast->{'_filt_func'}; 
    }

    if(defined $len) {
	if($len =~ /\D/ or $len <= 0) {
	    $Blast->warn("Invalid minimum length value: $len", 
			"Value must be an integer > 0. Value not set.");
	} else {
	    $Blast->{'_min_length'} = $len;
	} 
    }

    if(defined $func) {
	$Blast->{'_confirm_significance'} = 1;
	if($func and not ref $func eq 'CODE') {
	    $Blast->throw("Not a function reference: $func",
			  "The -filt_func parameter must be function reference.");
	  }
      }
  }


=head2 _parse_blast_stream

 Usage     : n/a. Internal method called by parse()
 Purpose   : Obtains the function to be used during parsing and calls read().
 Returns   : Integer (the number of blast reports read)
 Argument  : Named parameters  (forwarded from parse())
 Throws    : Propagates any exception thrown by _get_parse_blast_func() and read().

See Also   : L<_get_parse_blast_func>(), B<Bio::Root::Object::read()>

=cut

#----------------------
sub _parse_blast_stream {
#----------------------
    my ($self, %param) = @_;

    my $func = $self->_get_parse_blast_func(%param);
#    my $func = sub { my $data = shift; 
#		      printf STDERR "Chunk length = %d\n", length($data);
#		      sleep(3);
#		    };

    # Only setting the newline character once per session.
    $Newline ||= $Util->get_newline(-client => $self, %param);

    $self->read(-REC_SEP  =>"$Newline>", 
		-FUNC     => $func,
		%param);

    return $Blast->{'_blast_count'};
}



=head2 _get_parse_blast_func

 Usage     : n/a; internal method used by _parse_blast_stream()
           : $func_ref = $blast_object->_get_parse_blast_func()
 Purpose   : Generates a function ref to be used as a closure for parsing 
           : raw data as it is being loaded by Bio::Root::IOManager::read().
 Returns   : Function reference (closure).
 Comments  : The the function reference contains a fair bit of logic
           : at present. It could perhaps be split up into separate
           : functions to make it more 'digestible'.

See Also   : L<_parse_blast_stream>()

=cut

#--------------------------
sub _get_parse_blast_func {
#--------------------------
    my ($self, @param) = @_;

    my ($save_a, $exec_func) = 
	$self->_rearrange([qw(SAVE_ARRAY EXEC_FUNC)], @param); 

#    $MONITOR && print STDERR "\nParsing Blast stream (5/dot, 250/line)\n";
    my $count = 0;
    my $strict = $self->strict();

    # Some parameter validation.
    # Remember, all Blast parsing will use this function now.
    # You won't need a exec-func or save_array when just creating a Blast object
    #  as in: $blast = new Bio::Tools::Blast();
    if($exec_func and not ref($exec_func) eq 'CODE') {
	$self->throw("The -EXEC_FUNC parameter must be function reference.",
		    "exec_func = $exec_func");

    } elsif($save_a and not ref($save_a) eq 'ARRAY') {
	$self->throw("The -SAVE_ARRAY parameter must supply an array reference".
		     "when not using an -EXEC_FUNC parameter.");
    }

    ## Might consider breaking this closure up if possible.

     return sub {
	my ($data) = @_;
	## $data should contain one of three possible fragment types
        ## from a Blast report:
        ##   1. Header with description section,
        ##   2. An alignment section for a single hit, or
        ##   3. The final alignment section plus the footer section.
	## (record separator = "Newline>").

#	print STDERR "\n(BLAST) DATA CHUNK: $data\n";

	my ($current_blast, $current_prog, $current_vers, $current_db);
	my $prev_blast;
        my $contains_translation = 0;


### steve --- Wed Mar 15 02:48:07 2000
### In the process of addressing bug PR#95. Tricky.
### Using the $contains_translation to do so. Not complete
### and possibly won't fix. We'll see.

	# Check for header section. Start a new Blast object and 
	# parse the description section.
#        if ($data =~ /\sQuery\s?=/s || ($contains_translation && $data =~ /Database:/s)) {
        if ($data =~ /\sQuery\s?=/s) {
	    $Blast->{'_blast_count'}++;
	    print STDERR ".", $Blast->{'_blast_count'} % 50 ? '' : "\n" if $MONITOR;

            if($data =~ /$Newline\s+Translating/so) {
                print STDERR "\nCONTAINS TRANSLATION\n";
                $contains_translation = 1;
            }

	    # If we're parsing a stream containing multiple reports,
	    # all subsequent header sections will contain the last hit of
	    # the previous report which needs to be parsed and added to that
	    # report if signifcant. It also contains the run parameters
	    # at the bottom of the Blast report.
#	    if($Blast->{'_blast_count'} > 1 || $contains_translation) {
            if($Blast->{'_blast_count'} > 1) {
#	      print STDERR "\nMULTI-BLAST STREAM.\n";
	      $Blast->{'_multi_stream'} = 1;

	      if($data =~ /(.+?)$Newline(<\w+>)?(T?BLAST[NPX])\s+(.+?)$Newline(.+)/so) {
		($current_prog, $current_vers, $data) = ($3, $4, $5);
		# Final chunk containing last hit and last footer.
		$Blast->{'_current_blast'}->_parse_alignment($1);
		$prev_blast = $Blast->{'_current_blast'}; # finalized.
#              }	elsif($contains_translation) {
#                  $data =~ /(T?BLAST[NPX])\s+(.+?)$Newline(.+)/so;
#                  ($current_prog, $current_vers, $data) = ($1, $2, $3);
	      }	else {
		$Blast->throw("Can't determine program type from BLAST report.",
			      "Checked for: @Blast_programs.");
		# This has important implications for how to handle interval
		# information for HSPs. TBLASTN uses nucleotides in query HSP
		# but amino acids in the sbjct HSP sequence.
	      }

	      if($data =~ m/Database:\s+(.+?)$Newline/so ) {
		$current_db = $1;
	      } else {
		# In some reports, the Database is only listed at end.
		#$Blast->warn("Can't determine database name from BLAST report.");
	      }

	      # Incyte_Fix:   Nasty Invisible Bug.
	      # Records in blast report are delimited by '>', but... when
	      #  there are no hits for a query, there won't be a '>'.  That
	      #  causes several blast reports to run together in the data
	      #  passed to this routine.  Need to get rid of non-hits in data
	      if ($data =~ /.+(No hits? found.+Sequences.+)/so) {
		  $data = $1;
	      }
	      # End Incyte_Fix

	    }

	    # Determine if we need to create a new Blast object 
	    # or use the $self object for this method.

	    if($Blast->{'_multi_stream'} or $self->name eq 'Static Blast object') {
	      # Strict mode is not object-specific but may be someday.
#	      print STDERR "\nCreating new Blast object.\n";
	      $current_blast = new Bio::Tools::Blast(-STRICT => $strict);
	    } else {
	      $current_blast = $self;
	    }
	    $Blast->{'_current_blast'} = $current_blast;

	    # If we're not sharing stats, set data on current blast  object.
	    if(defined $current_prog and not $Blast->{'_share'}) {
	      $current_blast->program($current_prog);
	      $current_blast->program_version($current_vers);
	      $current_blast->database($current_db);
	    }

#	    print STDERR "CURRENT BLAST = ", $current_blast->name, "\n";
	    $current_blast->_parse_header($data);
	    
	    # If there were any descriptions in the header,
	    # we know if there are any significant hits.
	    # No longer throwing exception if there were no significant hits
	    # and a -signif parameter was specified. Doing so prevents the
	    # construction of a Blast object, which could still be useful. 
#	    if($current_blast->{'_has_descriptions'} and $Blast->{'_confirm_significance'} and not $current_blast->is_signif) {
#	      $current_blast->throw("No significant BLAST hits for ${\$current_blast->name}");

#	    }

	} # Done parsing header/description section

### For use with $contains_translation - not right - breaks regular report parsing.
# 	elsif(ref $Blast->{'_current_blast'} && $data !~ /\s*\w*\s*/s) {
 	elsif(ref $Blast->{'_current_blast'} ) {
	    # Process an alignment section. 
	    $current_blast = $Blast->{'_current_blast'};
#	    print STDERR "\nCONTINUING PROCESSING ALN WITH ", $current_blast->name, "\n";
#	    print STDERR "DATA: $data\n";
            eval {
                $current_blast->_parse_alignment($data);
            };
            if($@) {
     #           push @{$self->{'_blast_errs'}}, $@;
            }
	}
	
	# If the current Blast object has been completely parsed
	# (occurs with a single Blast stream), or if there is a previous 
	# Blast object (occurs with a multi Blast stream), 
	# execute a supplied function on it or store it in a supplied array.

	if( defined $prev_blast or $current_blast->{'_found_params'}) {
	  my $finished_blast = defined($prev_blast) ? $prev_blast : $current_blast;
	  
	  $finished_blast->_report_errors();
#	  print STDERR "\nNEW BLAST OBJECT: ${\$finished_blast->name}\n";

	  if($exec_func) {
#	    print STDERR "  RUNNING EXEC_FUNC...\n";
	    &$exec_func($finished_blast);  # ignoring any return value.
	    # Report processed, no longer need object.
	    $finished_blast->destroy;
	    undef $finished_blast;
	  } elsif($save_a) {
#	    print STDERR "  SAVING IN ARRAY...\n";
	    # We've already verified that if there is no exec_func
	    # then there must be a $save_array
	    push @$save_a, $finished_blast;
	  }
	}
	1;
      }
  }

=head2 _report_errors

 Title   : _report_errors
 Usage   : n/a; Internal method called by _get_parse_blast_func().
 Purpose : Throw or warn about any errors encountered. 
 Returns : n/a
 Args    : n/a
 Throws  : If all hits generated exceptions, raise exception
         :   (a fatal event for the Blast object.)
         : If some hits were okay but some were bad, generate a warning
         :   (a few bad applies should not spoil the bunch).
         :   This usually indicates a limiting B-value.
         : When the parsing code fails, it is either all or nothing.

=cut

#-------------------
sub _report_errors {
#-------------------
  my $self = shift;

  return unless ref($self->{'_blast_errs'});
#  ref($self->{'_blast_errs'}) || (print STDERR "\nNO ERRORS\n", return );

  my @errs = @{$self->{'_blast_errs'}};

  if(scalar @errs) {
    my ($str);
    @{$self->{'_blast_errs'}} = (); # clear the errs on the object.
    # When there are many errors, in most of the cases, they are
    # caused by the same problem. Only need to see full data for
    # the first one.
    if(scalar @errs > 2) {
      $str = "SHOWING FIRST EXCEPTION ONLY:\n$errs[0]";
      $self->clear_err();  # clearing the existing set of errors (conserve memory).
      # Not necessary, unless the -RECORD_ERR =>1
      # constructor option was used for Blast object.
    } else {
      $str = join("\n",@errs);
    }
    
    if(not $self->{'_num_hits_significant'}) {
      $self->throw(sprintf("Failed to parse any hit data (n=%d).", scalar(@errs)),
		   "\n\nTRAPPED EXCEPTION(S):\n$str\nEND TRAPPED EXCEPTION(S)\n"
		  );
    } else {
      $self->warn(sprintf("Some potential hits were not parsed (n=%d).", scalar(@errs)), 
		  @errs > 2 ? "This may be due to a limiting B value (max alignment listings)." : "",
		  "\n\nTRAPPED EXCEPTION(S):\n$str\nEND TRAPPED EXCEPTION(S)\n"
		 );
    }
  }
}


=head2 _parse_header

 Usage     : n/a; called automatically by the _get_parse_blast_func().
 Purpose   : Parses the header section of a BLAST report.
 Argument  : String containing the header+description section of a BLAST report.
 Throws    : Exception if description data cannot be parsed properly.
           : Exception if there is a 'FATAL' error in the Blast report.
           : Warning if there is a 'WARNING' in the Blast report.
           : Warning if there are no significant hits.
 Comments  : Description section contains a single line for each hit listing
           : the seq id, description, score, Expect or P-value, etc.

See Also   : L<_get_parse_blast_func>()

=cut

#----------------------
sub _parse_header {  
#----------------------
    my( $self, $data ) = @_;
    
#    print STDERR "\n$ID: PARSING HEADER\n"; #$data\n";
    
    $data =~ s/^\s+|\s+>?$//sg;

    if($data =~ /<HTML/i) {
      $self->throw("Can't parse HTML-formatted BLAST reports.",
#		    "Such reports can be parsed with a special parsing \n".
#		    "script included in the examples/blast directory \n".
#		    "of the Bioperl distribution. (TODO)"
		  );
      # This was the old strategy, can't do it with new strategy
      # since we don't have the whole report in one chunk.
      # This could be the basis for the "special parsing script".
#	 require Bio::Tools::Blast::HTML; 
#	 Bio::Tools::Blast::HTML->import(&strip_html); 
#	 &strip_html(\$data);
    }
	
    $data =~ /WARNING: (.+?)$Newline$Newline/so and $self->warn("$1") if $self->strict;
    $data =~ /FATAL: (.+?)$Newline$Newline/so and $self->throw("FATAL BLAST ERROR = $1"); 
    # No longer throwing exception when no hits were found. Still reporting it.
    $data =~ /No hits? found/i and $self->warn("No hits were found.") if $self->strict; 

    # If this is the first Blast, the program, version, and database info
    # pertain to it. Otherwise, they are for the previous report and have
    # already been parsed out.
    # Data is stored in the static Blast object. Data for subsequent reports
    # will be stored in separate objects if the -share parameter is not set.
    # See _get_parse_blast_func().

    if($Blast->{'_blast_count'} == 1) {
      if($data =~ /(<\w+>)?(T?BLAST[NPX])\s+(.+?)$Newline/so) {
	$Blast->program($2);
	$Blast->program_version($3);
      }	else {
	$self->throw("Can't determine program type from BLAST report.",
		     "Checked for: @Blast_programs.");
	# This has important implications for how to handle interval
	# information for HSPs. TBLASTN uses nucleotides in query HSP
	# but amino acids in the sbjct HSP sequence.
      }

      if($data =~ m/Database:\s+(.+?)$Newline/so ) {
	$Blast->database($1);
      } else {
	# In some reports, the Database is only listed at end.
	#$self->warn("Can't determine database name from BLAST report (_parse_header)\n$data\n.");
      }
    }

    my ($header, $descriptions);

    ## For efficiency reasons, we want to to avoid using $' and $`.
    ## Therefore using single-line mode pattern matching.

    if($data =~ /(.+?)\nSequences producing.+?\n(.+)/s ) {
        ($header, $descriptions) = ($1, $2);
	$self->{'_has_descriptions'} = 1;
    } else {
        $header = $data;
	$self->{'_has_descriptions'} = 0;
	# Blast reports can legally lack description section. No need to warn.
	#push @{$self->{'_blast_errs'}}, "Can't parse description data.";
    }

    $self->_set_query($header);  # The name of the sequence will appear in error report.
#    print STDERR "\nQUERY = ", $Blast->{'_current_blast'}->query, "\n";

    $self->_set_date($header) if $Blast->{'_get_stats'};
    $self->_set_length($header);

#    not $Blast->{'_confirm_significance'} and print STDERR "\nNOT PARSING DESCRIPTIONS.\n";

    # Setting the absolute max and min significance levels.
    $self->{'_highestSignif'} = 0;
    $self->{'_lowestSignif'} = $DEFAULT_SIGNIF;

    if ($Blast->{'_confirm_significance'} || $Blast->{'_no_aligns'}) {
      $self->_parse_descriptions($descriptions) if $descriptions;
    } else {
      $self->{'_is_significant'} = 1;
    }
  }


#-----------------------
sub _parse_descriptions {
#-----------------------
  my ($self, $desc) = @_;

    # NOTE: This method will not be called if the report lacks 
    #       a description section.

#    print STDERR "\nPARSING DESCRIPTION DATA\n";

    my @descriptions = split( $Newline, $desc);
    my($line);

    # NOW step through each line parsing out the P/Expect value
    # All we really need to do is check the first one, if it doesn't
    # meet the significance requirement, we can skip the report.
    # BUT: we want to collect data for all hits anyway to get min/max signif.

    my $my_signif = $self->signif;
    my $layout_set = $Blast->{'_layout'} || 0;
    my $layout;
    my $count = 0;
    my $sig;

    desc_loop:
  foreach $line (@descriptions) {
      $count++;
      last desc_loop if $line =~ / NONE |End of List/;
      next desc_loop if $line =~ /^\s*$/;
      next desc_loop if $line =~ /^\.\./;

	## Checking the significance value (P- or Expect value) of the hit
	## in the description line. 

      # These regexps need testing on a variety of reports.
      if ( $line =~ /\d+\s{1,5}[\de.-]+\s*$/) {
	$layout = 2;
      } elsif( $line =~ /\d+\s{1,5}[\de.-]+\s{1,}\d+\s*$/) {
	$layout = 1;
      } else {
	$self->warn("Can't parse significance data in description line $line");
	next desc_loop;
      }
      not $layout_set and ($self->_layout($layout), $layout_set = 1);
      
      $sig = &_parse_signif( $line, $layout );
      
#      print STDERR "  Parsed signif ($layout) = $sig\n"; 

      last desc_loop if ($sig > $my_signif and not $Blast->{'_check_all'});
      $self->_process_significance($sig, $my_signif);
    }

#  printf "\n%d SIGNIFICANT HITS.\nDONE PARSING DESCRIPTIONS.\n", $self->{'_num_hits_significant'};
}


sub _process_significance {
    my($self, $sig, $my_signif) = @_;

    $self->{'_highestSignif'} = ($sig > $self->{'_highestSignif'}) 
   	                        ? $sig : $self->{'_highestSignif'};

    $self->{'_lowestSignif'} = ($sig < $self->{'_lowestSignif'}) 
                                 ? $sig : $self->{'_lowestSignif'};

    # Significance value assessment.
    $sig <= $my_signif and $self->{'_num_hits_significant'}++;
    $self->{'_num_hits'}++;

    $self->{'_is_significant'} = 1 if $self->{'_num_hits_significant'};
}

=head2 _parse_alignment

 Usage     : n/a; called automatically by the _get_parse_blast_func().
 Purpose   : Parses a single alignment section of a BLAST report.
 Argument  : String containing the alignment section.
 Throws    : n/a; All errors are trapped while parsing the hit data
           : and are processed as a group when the report is 
           : completely processed (See _report_errors()).
           :
 Comments  : Alignment section contains all HSPs for a hit.
           : Requires Bio::Tools::Blast::Sbjct.pm.
           : Optionally calls a filter function to screen the hit on arbitrary
           : criteria. If the filter function returns true for a given hit,
           : that hit will be skipped.
           :
           : If the Blast object was created with -check_all_hits set to true,
           : all hits will be checked for significance and processed if necessary.
           : If this field is false, the parsing will stop after the first
           : non-significant hit. 
           : See parse() for description of parsing parameters.

See Also   : L<parse>(), L<_get_parse_blast_func>(), L<_report_errors>(), B<Bio::Tools::Blast::Sbjct()>, L<Links to related modules>

=cut

#----------------------
sub _parse_alignment {  
#----------------------
# This method always needs to check detect if the $data argument
# contains the footer of a Blast report, indicating the last chunk 
# of a single Blast stream.

    my( $self, $data ) = @_;

#    printf STDERR "\nPARSING ALIGNMENT DATA for %s $self.\n", $self->name;

    # NOTE: $self->{'_current_hit'} is an instance variable
    #       The $Blast object will not have this member.

    # If all of the significant hits have been parsed,
    # return if we're not checking all or if we don't need to get 
    # the Blast stats (parameters at footer of report).
    if(defined $self->{'_current_hit'} and 
      defined $self->{'_num_hits_significant'}) {
      return if $self->{'_current_hit'} >= $self->{'_num_hits_significant'} and
	not ($Blast->{'_check_all'} or $Blast->{'_get_stats'});
    }

    # Check for the presence of the Blast footer section.
    # _parse_footer returns the alignment section.
    $data = $self->_parse_footer($data);

    # Return if we're only interested in the best hit.
    # This has to occur after checking for the parameters section 
    # in the footer (since we may still be interested in them).
    return if $Blast->best and ( defined $self->{'_current_hit'} and $self->{'_current_hit'} >=1);


#    print "RETURNED FROM _parse_footer (", $self->to_string, ")";
#    print "\n  --> FOUND PARAMS.\n" if $self->{'_found_params'};
#    print "\n  --> DID NOT FIND PARAMS.\n" unless $self->{'_found_params'};

    require Bio::Tools::Blast::Sbjct;

    $data =~ s/^\s+|\s+>?$//sg;
    $data =~ s/$Newline$Newline/$Newline/sog;  # remove blank lines.
    my @data = split($Newline, $data);
    push @data, 'end';

#    print STDERR "\nALIGNMENT DATA:\n$data\n"; 

    my $prog       = $self->program;
    my $check_all  = $Blast->{'_check_all'};
    my $filt_func  = $Blast->{'_filt_func'} || 0;
    my $signif_fmt = $Blast->{'_signif_fmt'};
    my $my_signif  = $self->signif;
    my $err;

    # Now construct the Sbjct objects from the alignment section

#	debug(1);
    
    $self->{'_current_hit'}++;
    
    # If not confirming significance, _parse_descriptions will not have been run,
    # so we need to count the total number of hits here.
    if( not $Blast->{'_confirm_significance'}) {
      $self->{'_num_hits'}++;
    }

    if($Blast->{'_no_aligns'}) {
#        printf STDERR "\nNOT PARSING ALIGNMENT DATA\n";
        return;
    }

    my $hit;  # Must be my'ed within hit_loop.
    eval {
      $hit = new Bio::Tools::Blast::Sbjct (-DATA      =>\@data, 
					   -PARENT    =>$self, 
					   -NAME      =>$self->{'_current_hit'}, 
					   -RANK      =>$self->{'_current_hit'}, 
					   -RANK_BY   =>'order',
					   -PROGRAM   =>$prog, 
					   -SIGNIF_FMT=>$signif_fmt,
					   -OVERLAP   =>$Blast->{'_overlap'} || $MAX_HSP_OVERLAP,
					  );
#      printf STDERR "NEW HIT: %s, SIGNIFICANCE = %g\n", $hit->name, $hit->expect;  <STDIN>;
      # The BLAST report may have not had a description section.
      if(not $self->{'_has_descriptions'}) {
	  $self->_process_significance($hit->signif, $my_signif);
      }
    };

    if($@) {
      # Throwing lots of errors can slow down the code substantially.
      # Error handling code is not that efficient.
      #print STDERR "\nERROR _parse_alignment: $@\n";
      push @{$self->{'_blast_errs'}}, $@;
      $hit->destroy if ref $hit;
      undef $hit;
    } else {
      # Collect overall signif data if we don't already have it,
      # (as occurs if no -signif parameter is supplied).
      my $hit_signif = $hit->signif;

      if (not $Blast->{'_confirm_significance'} ) {
	$self->{'_highestSignif'} = ($hit_signif > $self->{'_highestSignif'}) 
                                    ? $hit_signif : $self->{'_highestSignif'};

	$self->{'_lowestSignif'} = ($hit_signif < $self->{'_lowestSignif'}) 
                                    ? $hit_signif : $self->{'_lowestSignif'};
      }

      # Test significance using custom function (if supplied)
      if($filt_func) {
	if(&$filt_func($hit)) {
	  push @{$self->{'_hits'}}, $hit;
	} else {
	  $hit->destroy; undef $hit; 
	}
      } elsif($hit_signif <= $my_signif) {
	push @{$self->{'_hits'}}, $hit;
      }
    }
    
  }


=head2 _parse_footer

 Usage     : n/a; internal function. called by _parse_alignment()
 Purpose   : Extracts statistical and other parameters from the BLAST report.
           : Sets various key elements such as the program and version,
           : gapping, and the layout for the report (blast1 or blast2).
 Argument  : Data to be parsed.
 Returns   : String containing an alignment section for processing by
           : _parse_alignment().
 Throws    : Exception if cannot find the parameters section of report.
           : Warning if cannot determine if gapping was used.
           : Warning if cannot determine the scoring matrix used.
 Comments  : This method must always get called, even if the -STATS
           : parse() parameter is false. The reason is that the layout
           : of the report  and the presence of gapping must always be set.
           : The determination whether to set additional stats is made 
           : by methods called by _parse_footer().

See Also   : L<parse>(), L<_parse_alignment>(), L<_set_database>()

=cut

#---------------------
sub _parse_footer {
#---------------------
# Basic strategy:
# 1. figure out if we're supposed to get the stats,
# 2. figure out if the stats are to be shared. some, not all can be shared 
#    (eg., db info and matrix can be shared, karlin altschul params cannot.
#    However, this method assumes they are all sharable.)
# 3. Parse the stats.
# 4. return the block before the parameters section if the supplied data
#    contains a footer parameters section.

    my ($self, $data) = @_;
    my ($client, $last_align, $params);

#    printf STDERR "\nPARSING PARAMETERS for %s $self.\n", $self->name;

    # Should the parameters be shared?
    # If so, set $self to be the static $Blast object and return if 
    # the parameters were already set.
    # Before returning, we need to extract the last alignment section
    # from the parameter section, if any.

    if ($Blast->{'_share'}) {
      $client = $self;
      $self = $Blast if $Blast->{'_share'};
    }

    my $get_stats = $Blast->{'_get_stats'};

    if( $data =~ /(.+?)${Newline}CPU time: (.*)/so) {
	# NCBI-Blast2 format (v2.04).
	($last_align, $params) = ($1, $2);
	return $last_align if $client->{'_found_params'};
	$self->_set_blast2_stats($params);

    } elsif( $data =~ /(.+?)${Newline}Parameters:(.*)/so) {
    # NCBI-Blast1 or WashU-Blast2 format.
    ($last_align, $params) = ($1, $2);
    return $last_align if $client->{'_found_params'};
    $self->_set_blast1_stats($params);

    } elsif( $data =~ /(.+?)$Newline\s+Database:(.*)/so) {
        # Gotta watch out for confusion with the Database: line in the header
        # which will be present in the last hit of an internal Blast report 
        # in a multi-report stream.

	# NCBI-Blast2 format (v2.05).
	($last_align, $params) = ($1, $2);
	return $last_align if $client->{'_found_params'};
        $self->_set_blast2_stats($params);
    
    } elsif( $data =~ /(.+?)$Newline\s*Searching/so) {
	# trying to detect a Searching at the end of a PSI-blast round.
        # Gotta watch out for confusion with the Searching line in the header
        # which will be present in the last hit of an internal Blast report 
        # in a multi-report, non-PSI-blast stream.

	# PSI-Blast format (v2.08).
	($last_align) = ($1);
	return $last_align; # if $client->{'_found_params'};
    }
    
    # If parameter section was found, set a boolean, 
    # otherwise return original data.

    if (defined($params)) {
      $client->{'_found_params'} = 1;
    } else {
      return $data;
    }

    $self->_set_database($params) if $get_stats;

    # The {'_gapped'} member should be set in the _set_blast?_stats() call.
    # This is a last minute attempt to deduce it.
    
    if(!defined($self->{'_gapped'})) {
	if($self->program_version() =~ /^1/) {
	    $self->{'_gapped'} = 0; 
	} else {
	    if($self->strict > 0) {
		$self->warn("Can't determine if gapping was used. Assuming gapped.");
	    }
	    $self->{'_gapped'} = 1; 
	}
    }

    return $last_align;
}


=head2 _set_blast2_stats

 Usage     : n/a; internal function called by _parse_footer()
 Purpose   : Extracts statistical and other parameters from BLAST2 report footer.
           : Stats collected: database release, gapping,
           : posted date, matrix used, filter used, Karlin-Altschul parameters, 
           : E, S, T, X, W.
 Throws    : Exception if cannot get "Parameters" section of Blast report.

See Also   : L<parse>(), L<_parse_footer>(), L<_set_database>(), B<Bio::Tools::SeqAnal::set_date()>,L<Links to related modules>

=cut

#---------------------'
sub _set_blast2_stats {
#---------------------
    my ($self, $data) = (@_);
    
    if($data =~ /$Newline\s*Gapped/so) {
	$self->{'_gapped'} = 1;
    } else {
 	$self->{'_gapped'} = 0;
    }

    # Other stats are not always essential.
    return unless $Blast->{'_get_stats'};

    # Blast2 Doesn't report what filter was used in the parameters section.
    # It just gives a warning that *some* filter was used in the header. 
    # You just have to know the defaults (currently: protein = SEG, nucl = DUST).
    if($data =~ /\bfiltered\b/si) {
	$self->{'_filter'} = 'DEFAULT FILTER';
    } else {
	$self->{'_filter'} = 'NONE';
    }

    if($data =~ /Gapped$Newline\s*Lambda +K +H$Newline +(.+?)$Newline/so) {
	my ($l, $k, $h) = split(/\s+/, $1);
	$self->{'_lambda'} = $l || 'UNKNOWN';
	$self->{'_k'} = $k || 'UNKNOWN';
	$self->{'_h'} = $h || 'UNKNOWN';
    } elsif($data =~ /Lambda +K +H$Newline +(.+?)$Newline/so) {
	my ($l, $k, $h) = split(/\s+/, $1);
	$self->{'_lambda'} = $l || 'UNKNOWN';
	$self->{'_k'} = $k || 'UNKNOWN';
	$self->{'_h'} = $h || 'UNKNOWN';
    }
    
    if($data =~ /$Newline\s*Matrix: (.+?)$Newline/so) {
	$self->{'_matrix'} = $1;
    } else {
	$self->{'_matrix'} = $DEFAULT_MATRIX.'?'; 
	if($self->strict > 0) {
	    $self->warn("Can't determine scoring matrix. Assuming $DEFAULT_MATRIX.");
	}
    }

    if($data =~ /$Newline\s*Gap Penalties: Existence: +(\d+), +Extension: (\d+)$Newline/so) {
	$self->{'_gapCreation'} = $1;
	$self->{'_gapExtension'} = $2;
    }
    if($data =~ /sequences better than (\d+):/s) {
	$self->{'_expect'} = $1;
    }

    if($data =~ /$Newline\s*T: (\d+)/o) { $self->{'_word_size'} = $1; }
    if($data =~ /$Newline\s*A: (\d+)/o) { $self->{'_a'} = $1; }
    if($data =~ /$Newline\s*S1: (\d+)/o) { $self->{'_s'} = $1; }
    if($data =~ /$Newline\s*S2: (\d+)/o) { $self->{'_s'} .= ", $1"; }
    if($data =~ /$Newline\s*X1: (\d+)/o) { $self->{'_x1'} = $1; }
    if($data =~ /$Newline\s*X2: (\d+)/o) { $self->{'_x2'} = $1; }
}


=head2 _set_blast1_stats

 Usage     : n/a; internal function called by _parse_footer()
 Purpose   : Extracts statistical and other parameters from BLAST 1.x style eports.
           : Handles NCBI Blast1 and WashU-Blast2 formats.
           : Stats collected: database release, gapping, 
           : posted date, matrix used, filter used, Karlin-Altschul parameters, 
           : E, S, T, X, W.

See Also   : L<parse>(), L<_parse_footer>(), L<_set_database>(), B<Bio::Tools::SeqAnal::set_date()>,L<Links to related modules>

=cut

#----------------------
sub _set_blast1_stats {
#----------------------
    my ($self, $data) = (@_);
    
    if(!$self->{'_gapped'} and $self->program_version() =~ /^2[\w\-\.]+WashU/) {
	$self->_set_gapping_wu($data);
    } else {
	$self->{'_gapped'} = 0;
    }

    # Other stats are not always essential.
    return unless $Blast->{'_get_stats'};

    if($data =~ /filter=(.+?)$Newline/so) {
	$self->{'_filter'} = $1;
    } elsif($data =~ /filter$Newline +(.+?)$Newline/so) {
	$self->{'_filter'} = $1;
    } else {
	$self->{'_filter'} = 'NONE';
    }
    
    if($data =~ /$Newline\s*E=(\d+)$Newline/so) {  $self->{'_expect'} = $1; }

    if($data =~ /$Newline\s*M=(\w+)$Newline/so) {  $self->{'_matrix'} = $1; }

    if($data =~ /\s*Frame  MatID Matrix name .+?$Newline +(.+?)$Newline/so) {
	## WU-Blast2.
	my ($fr, $mid, $mat, $lu, $ku, $hu, $lc, $kc, $hc) = split(/\s+/,$1);
	$self->{'_matrix'} = $mat || 'UNKNOWN';
	$self->{'_lambda'} = $lu || 'UNKNOWN';
	$self->{'_k'} = $ku || 'UNKNOWN';
	$self->{'_h'} = $hu || 'UNKNOWN';
	
    } elsif($data =~ /Lambda +K +H$Newline +(.+?)$Newline/so) {
	## NCBI-Blast1.
	my ($l, $k, $h) = split(/\s+/, $1);
	$self->{'_lambda'} = $l || 'UNKNOWN';
	$self->{'_k'} = $k || 'UNKNOWN';
	$self->{'_h'} = $h || 'UNKNOWN';
    }

    if($data =~ /E +S +W +T +X.+?$Newline +(.+?)$Newline/so) {
	# WashU-Blast2
	my ($fr, $mid, $len, $elen, $e, $s, $w, $t, $x, $e2, $s2) = split(/\s+/,$1);
	$self->{'_expect'} ||= $e || 'UNKNOWN';
	$self->{'_s'} = $s || 'UNKNOWN';
	$self->{'_word_size'} = $w || 'UNKNOWN';
	$self->{'_t'} = $t || 'UNKNOWN';
	$self->{'_x'} = $x || 'UNKNOWN';
    
    } elsif($data =~ /E +S +T1 +T2 +X1 +X2 +W +Gap$Newline +(.+?)$Newline/so) {
	## NCBI-Blast1.
	my ($e, $s, $t1, $t2, $x1, $x2, $w, $gap) = split(/\s+/,$1);
	$self->{'_expect'} ||= $e || 'UNKNOWN';
	$self->{'_s'} = $s || 'UNKNOWN';
	$self->{'_word_size'} = $w || 'UNKNOWN';
	$self->{'_t1'} = $t1 || 'UNKNOWN';
	$self->{'_t2'} = $t2 || 'UNKNOWN';
	$self->{'_x1'} = $x1 || 'UNKNOWN';
	$self->{'_x2'} = $x2 || 'UNKNOWN';
	$self->{'_gap'} = $gap || 'UNKNOWN';
    }

    if(!$self->{'_matrix'}) {
	$self->{'_matrix'} = $DEFAULT_MATRIX.'?';
	if($self->strict > 0) {
	    $self->warn("Can't determine scoring matrix. Assuming $DEFAULT_MATRIX.");
	}
    }
}


=head2 _set_gapping_wu

 Usage     : n/a; internal function called by _set_blast1_stats()
 Purpose   : Determine if gapping_wu was on for WashU Blast reports.
 Comments  : In earlier versions, gapping was always specified
           : but in the current version (2.0a19MP), gapping is on by default
           : and there is no positive "gapping" indicator in the Parameters
           : section.

See Also   : L<_set_blast1_stats>()

=cut

#--------------------
sub _set_gapping_wu {
#--------------------
    my ($self, $data) = @_;
    
    if($data =~ /gaps?$Newline/so) {
	$self->{'_gapped'} = ($data =~ /nogaps?$Newline/so) ? 0 : 1;
    } else {
	$self->{'_gapped'} = 1;
    }
}


=head2 _set_date

 Usage     : n/a; internal function called by _parse_footer()
 Purpose   : Determine the date on which the Blast analysis was performed.
 Comments  : Date information is not consistently added to Blast output.
           : Uses superclass method set_date() to set date from the file,
           : (if any).

See Also   : L<_parse_footer>(), B<Bio::Tools::SeqAnal::set_date()>,L<Links to related modules>

=cut

#--------------
sub _set_date {
#--------------
    my $self = shift;
    my $data = shift;

    ### Network BLAST reports from NCBI are time stamped as follows:
    #Fri Apr 18 15:55:41 EDT 1997, Up 1 day, 19 mins, 1 user, load: 19.54, 19.13, 17.77    
    if($data =~ /Start:\s+(.+?)\s+End:/s) {
	## Calling superclass method to set the date.
	## If we can't get date from the report, file date is obtained.
	$self->set_date($1);
    } elsif($data =~ /Date:\s+(.*?)$Newline/so) {
	## E-mailed reports have a Date: field
	$self->set_date($1);
    } elsif( $data =~ /done\s+at (.+?)$Newline/so ) {
	$self->set_date($1);  
    } elsif( $data =~ /$Newline([\w:, ]+), Up \d+/so ) {
	$self->set_date($1); 
    } else {
	## Otherwise, let superclass attempt to get the file creation date.
	$self->set_date() if $self->file;
    }
}




=head2 _set_length

 Usage     : n/a; called automatically during Blast report parsing.
 Purpose   : Sets the length of the query sequence (extracted from report).
 Returns   : integer (length of the query sequence)
 Throws    : Exception if cannot determine the query sequence length from
           :           the BLAST report.
           : Exception if the length is below the min_length cutoff (if any).
 Comments  : The logic here is a bit different from the other _set_XXX()
           : methods since the significance of the BLAST report is assessed 
           : if MIN_LENGTH is set.

See Also   : B<Bio::Tools::SeqAnal::length()>, L<Links to related modules>

=cut

#---------------
sub _set_length {
#---------------
    my ($self, $data) = @_;

    my ($length);
    if( $data =~ m/$Newline\s+\(([\d|,]+) letters\)/so ) {
	$length = $1;
	$length =~ s/,//g;
#	printf "Length = $length in BLAST for %s$Newline",$self->name; <STDIN>;
    } else {
	$self->throw("Can't determine sequence length from BLAST report.");
    }

    my($sig_len);
    if(defined($Blast->{'_min_length'})) {
      local $^W = 0;
      if($length < $Blast->{'_min_len'}) {
	$self->throw("Query sequence too short for ${\$self->name} ($length)",
		     "Minimum  length is $Blast->{'_min_len'}");
      }
    }   

    $self->length($length);  # defined in superclass.
}



=head2 _set_database

 Usage     : n/a; called automatically during Blast report parsing.
 Purpose   : Sets the name of the database used by the BLAST analysis.
           : Extracted from raw BLAST report.
 Throws    : Exception if the name of the database cannot be determined.
 Comments  : The database name is used by methods or related objects
           : for database-specific parsing.

See Also   : L<parse>(), B<Bio::Tools::SeqAnal::database()>,B<Bio::Tools::SeqAnal::_set_db_stats()>,L<Links to related modules>

=cut

#------------------
sub _set_database {
#------------------
# This now only sets data base information extracted from the report footer.

    my ($self, $data) = @_;

    my ($name, $date, $lets, $seqs);

    my $strict = $self->strict > 0;

    # This is fail-safe since DB name usually gets set in _parse_header()
    # In some reports, the database is only listed at bottom (NCBI 2.0.8).
    if($data =~ m/Database: +(.+?)$Newline/so ) {
      $name = $1;
    } elsif(not $self->database) {
      $self->warn("Can't determine database name from BLAST report.");
    }

    if($data =~ m/Posted date: +(.+?)$Newline/so ) {
	$date = $1;
    } elsif($data =~ m/Release date: +(.+?)$Newline/so ) {
	$date = $1;
    } elsif($strict) {
	$self->warn("Can't determine database release date.");
    }

    if($data =~ m/letters in database: +([\d,]+)/si ||
       $data =~ m/length of database: +([\d,]+)/si ) {
	$lets = $1;
    } elsif($strict) {
	$self->warn("Can't determine number of letters in database.\n$data\n");
    }

    if($data =~ m/sequences in database: +([\d,]+)/si ||
      $data =~ m/number of sequences: +([\d,]+)/si ) {
	$seqs = $1;
    } elsif($strict) {
	$self->warn("Can't determine number of sequences in database.\n$data\n");
    }

    $self->_set_db_stats( -NAME    => $name,
			  -RELEASE => $date || '',
			  -LETTERS => $lets || '',
			  -SEQS    => $seqs || ''
			  );
}



=head2 _set_query

 Usage     : n/a; called automatically during Blast report parsing.
 Purpose   : Set the name of the query and the query description.
           : Extracted from the raw BLAST report.
 Returns   : String containing name of query extracted from report.
 Throws    : Warning if the name of the query cannont be obtained.

See Also   : B<Bio::Tools::SeqAnal::query_desc()>,L<Links to related modules>

=cut

#---------------
sub _set_query {
#---------------
    my $self = shift;
    my $data = shift;
    
    if($data =~ m/${Newline}Query= *(.+?)$Newline/so ) {
	my $info = $1;
	$info =~ s/TITLE //;
	# Split the query line into two parts.
	# Using \s instead of ' '
	$info =~ /(\S+?)\s(.*)/;
	$self->query_desc($2 || '');
	# set name of Blast object and return.
	$self->name($1 || 'UNKNOWN');
    } else {
	$self->warn("Can't determine query sequence name from BLAST report.");
    }
#    print STDERR "$Newline  NAME = ${\$self->name}$Newline";
}



=head2 _parse_signif

 Usage     : &_parse_signif(string, layout, gapped);
           : This is a class function.
 Purpose   : Extracts the P- or Expect value from a single line of a BLAST description section.
 Example   : &_parse_signif("PDB_UNIQUEP:3HSC_  heat-shock cognate ...   799  4.0e-206  2", 1);
           : &_parse_signif("gi|758803  (U23828) peritrophin-95 precurs   38  0.19", 2);
 Argument  : string = line from BLAST description section
           : layout = integer (1 or 2)
           : gapped = boolean (true if gapped Blast).
 Returns   : Float (0.001 or 1e-03)
 Status    : Static

=cut

#------------------
sub _parse_signif {
#------------------
    my ($line, $layout, $gapped) = @_;

    local $_ = $line;
    my @linedat = split();

    # When processing both Blast1 and Blast2 reports 
    # in the same run, offset needs to be configured each time. 
    
    my $offset  = 0; 
    $offset  = 1 if $layout == 1 or not $gapped;

    my $signif = $linedat[ $#linedat - $offset ];

    # fail-safe check
    if(not $signif =~ /[.-]/) {
	$offset = ($offset == 0 ? 1 : 0);
	$signif = $linedat[ $#linedat - $offset ];
    }

    $signif = "1$signif" if $signif =~ /^e/i;
    return $signif;
}


## 
## BEGIN ACCESSOR METHODS THAT INCORPORATE THE STATIC $Blast OBJECT.
##

sub program { 
## Overridden method to incorporate the BLAST object.
    my $self = shift;  
    return $self->SUPER::program(@_) if @_;          # set
    $self->SUPER::program || $Blast->SUPER::program; # get
}


sub program_version { 
## Overridden method to incorporate the BLAST object.
    my $self = shift;  
    return $self->SUPER::program_version(@_) if @_;                  # set
    $self->SUPER::program_version || $Blast->SUPER::program_version; # get
}

sub database { 
## Overridden method to incorporate the BLAST object.
    my $self = shift;  
    return $self->SUPER::database(@_) if @_;           # set
    $self->SUPER::database || $Blast->SUPER::database; # get
}

sub database_letters { 
## Overridden method to incorporate the BLAST object.
    my $self = shift;  
    return $self->SUPER::database_letters(@_) if @_;           # set
    $self->SUPER::database_letters || $Blast->SUPER::database_letters; # get
}

sub database_release { 
## Overridden method to incorporate the BLAST object.
    my $self = shift;  
    return $self->SUPER::database_release(@_) if @_;           # set
    $self->SUPER::database_release || $Blast->SUPER::database_release; # get
}

sub database_seqs { 
## Overridden method to incorporate the BLAST object.
    my $self = shift;  
    return $self->SUPER::database_seqs(@_) if @_;           # set
    $self->SUPER::database_seqs || $Blast->SUPER::database_seqs; # get
}


sub date { 
## Overridden method to incorporate the BLAST object.
    my $self = shift;  
    return $self->SUPER::date(@_) if @_;       # set
    $self->SUPER::date || $Blast->SUPER::date; # get
}


sub best { 
## Overridden method to incorporate the BLAST object.
    my $self = shift;  
    return $Blast->SUPER::best(@_) if @_;       # set
    $Blast->SUPER::best; # get
}


=head2 signif

 Usage     : $blast->signif();
 Purpose   : Gets the P or Expect value used as significance screening cutoff.
 Returns   : Scientific notation number with this format: 1.0e-05.
 Argument  : n/a
 Comments  : Screening of significant hits uses the data provided on the
           : description line. For Blast1 and WU-Blast2, this data is P-value.
           : for Blast2 it is an Expect value. 
           :
           : Obtains info from the static $Blast object if it has not been set
           : for the current object.

See Also   : L<_set_signif>()

=cut

#-----------
sub signif { 
#-----------
    my $self = shift;  
    my $sig = $self->{'_significance'} || $Blast->{'_significance'};
    sprintf "%.1e", $sig;
}



=head2 is_signif

 Usage     : $blast->is_signif();
 Purpose   : Determine if the BLAST report contains significant hits.
 Returns   : Boolean
 Argument  : n/a
 Comments  : BLAST reports without significant hits but with defined
           : significance criteria will throw exceptions during construction.
           : This obviates the need to check significant() for
           : such objects.

See Also   : L<_set_signif>()

=cut

#------------
sub is_signif { my $self = shift; return $self->{'_is_significant'}; }
#------------

# is_signif() doesn't incorporate the static $Blast object but is included
# here to be with the other 'signif' methods.



=head2 signif_fmt

 Usage     : $blast->signif_fmt( [FMT] );
 Purpose   : Allows retrieval of the P/Expect exponent values only
           : or as a two-element list (mantissa, exponent).
 Usage     : $blast_obj->signif_fmt('exp'); 
           : $blast_obj->signif_fmt('parts');
 Returns   : String or '' if not set.
 Argument  : String, FMT = 'exp' (return the exponent only)
           :             = 'parts'(return exponent + mantissa in 2-elem list)
           :              = undefined (return the raw value)
 Comments  : P/Expect values are still stored internally as the full,
           : scientific notation value. 
           : This method uses the static $Blast object since this issue
           : will pertain to all Blast reports within a given set.
           : This setting is propagated to Bio::Tools::Blast::Sbjct.pm.

=cut

#-------------
sub signif_fmt { 
#-------------
    my $self = shift; 
    if(@_) { $Blast->{'_signif_fmt'} = shift; }
    $Blast->{'_signif_fmt'} || '';
}


=head2 min_length

 Usage     : $blast->min_length();
 Purpose   : Gets the query sequence length used as significance screening criteria.
 Returns   : Integer
 Argument  : n/a
 Comments  : Obtains info from the static $Blast object if it has not been set
           : for the current object.

See Also   : L<_set_signif>(), L<signif>()

=cut

#--------------
sub min_length { 
#--------------
    my $self = shift;  
    $self->{'_min_length'} || $Blast->{'_min_length'};
}




=head2 gapped

 Usage     : $blast->gapped();
 Purpose   : Set/Get boolean indicator for gapped BLAST.
 Returns   : Boolean
 Argument  : n/a
 Comments  : Obtains info from the static $Blast object if it has not been set
           : for the current object.

=cut

#-----------
sub gapped { 
#-----------
    my $self = shift; 
    if(@_) { $self->{'_gapped'} = shift; }
    $self->{'_gapped'} || $Blast->{'_gapped'}; 
}


=head2 _get_stats

 Usage     : n/a; internal method.
 Purpose   : Set/Get indicator for collecting full statistics from report.
 Returns   : Boolean (0 | 1)
 Comments  : Obtains info from the static $Blast object which gets set
           : by _init_parse_params().

=cut

#---------------
sub _get_stats { 
#---------------
    my $self = shift; 
    $Blast->{'_get_stats'};
}


=head2 _layout

 Usage     : n/a; internal method.
 Purpose   : Set/Get indicator for the layout of the report.
 Returns   : Integer (1 | 2)
           : Defaults to 2 if not set.
 Comments  : Blast1 and WashU-Blast2 have a layout = 1.
           : This is intended for internal use by this and closely
           : allied modules like Sbjct.pm and HSP.pm.
           :
           : Obtains info from the static $Blast object if it has not been set
           : for the current object.

=cut

#------------
sub _layout { 
#------------
    my $self = shift; 
    if(@_) { 
      # Optimization if we know all reports share the same stats.
      if($Blast->{'_share'}) {
	$Blast->{'_layout'} = shift;
      } else {
	$self->{'_layout'} = shift; 
      }
    }
    $self->{'_layout'} || $Blast->{'_layout'} || 2; 
}



## 
## END ACCESSOR METHODS THAT INCORPORATE THE STATIC $Blast OBJECT.
##



=head2 hits

 Usage     : $blast->hits();
 Purpose   : Get a list containing all BLAST hit (Sbjct) objects.
           : Get the numbers of significant hits.
 Examples  : @hits       = $blast->hits();
           : $num_signif = $blast->hits();
 Returns   : List context : list of Bio::Tools::Blast::Sbjct.pm objects
           :                or an empty list if there are no hits.
           : Scalar context: integer (number of significant hits)
           :                 or zero if there are no hits.
           :                 (Equivalent to num_hits()).
 Argument  : n/a. Relies on wantarray.
 Throws    : n/a.
           : Not throwing exception because the absence of hits may have
           : resulted from stringent significance criteria, not a failure
           : set the hits.

See Also   : L<hit>(), L<num_hits>(), L<is_signif>(), L<_set_signif>()

=cut

#----------
sub hits {
#----------
    my $self = shift;

    if(wantarray) {
        my @ary = ref($self->{'_hits'}) ? @{$self->{'_hits'}} : ();
        return @ary;
    } else {
        return $self->num_hits();
    }
        
#    my $num = ref($self->{'_hits'}) ? scalar(@{$self->{'_hits'}}) : 0;
#    my @ary = ref($self->{'_hits'}) ? @{$self->{'_hits'}} : ();
#
#    return wantarray 
#        #  returning list containing all hits or empty list.
#        ?  $self->{'_is_significant'} ? @ary : ()
#        #  returning number of hits or 0.
#        :  $self->{'_is_significant'} ? $num : 0;
}


=head2 hit

 Example   : $blast_obj->hit( [class] )
 Purpose   : Get a specific hit object.
           : Provides some syntactic sugar for the hits() method.
 Usage     : $hitObj = $blast->hit();
           : $hitObj = $blast->hit('best');
           : $hitObj = $blast->hit('worst');
           : $hitObj = $blast->hit( $name );
 Returns   : Object reference for a Bio::Tools::Blast::Sbjct.pm object.
           : undef if there are no hit (Sbjct) objects defined.
 Argument  : Class (or no argument).
           :   No argument (default) = highest scoring hit (same as 'best').
           :   'best' or 'first' = highest scoring hit.
           :   'worst' or 'last' = lowest scoring hit.
           :   $name = retrieve a hit by seq id (case-insensitive).
 Throws    : Exception if the Blast object has no significant hits.
           : Exception if a hit cannot be found when supplying a specific
           : hit sequence identifier as an argument.
 Comments  : 'best'  = lowest significance value (P or Expect) among significant hits.
           : 'worst' = highest sigificance value (P or Expect) among significant hits.

See Also   : L<hits>(), L<num_hits>(), L<is_signif>()

=cut

#---------
sub hit {
#---------
    my( $self, $option) = @_;
    $option ||= 'best';
    
    if($Blast->{'_no_aligns'} || ! ref($self->{'_hits'})) {
        return undef;
    }

    $self->{'_is_significant'} or 
	$self->throw("There were no significant hits.",
		     "Use num_hits(), hits(), is_signif() to check.");

    my @hits = @{$self->{'_hits'}};
    
    return $hits[0]      if $option =~ /^(best|first|1)$/i;
    return $hits[$#hits] if $option =~ /^(worst|last)$/i;

    # Get hit by name.	    
    foreach ( @hits ) {
	return $_ if $_->name() =~ /$option/i;
    }

    $self->throw("Can't get hit for: $option");
}



=head2 num_hits

 Usage     : $blast->num_hits( ['total'] );
 Purpose   : Get number of significant hits or number of total hits.
 Examples  : $num_signif = $blast-num_hits;
           : $num_total  = $blast->num_hits('total');
 Returns   : Integer
 Argument  : String = 'total' (or no argument).
           :   No argument (Default) = return number of significant hits.
           :   'total' = number of total hits.
 Throws    : n/a.
           : Not throwing exception because the absence of hits may have
           : resulted from stringent significance criteria, not a failure
           : set the hits.
 Comments  : A significant hit is defined as a hit with an expect value
           : (or P value for WU-Blast) at or below the -signif parameter
           : used when parsing the report. Additionally, if a filter function
           : was supplied, the significant hit must also pass that
           : criteria.
See Also   : L<hits>(), L<hit>(), L<is_signif>(), L<_set_signif>(), L<parse>()

=cut

#-------------
sub num_hits {
#-------------
    my( $self, $option) = @_;
    $option ||= '';

   $option =~ /total/i and return $self->{'_num_hits'} || 0;

    # Default: returning number of significant hits.
#    return $self->{'_num_hits_significant'} || 0;
#    return 0 if not ref $self->{'_hits'};

    if(ref $self->{'_hits'}) {
        return scalar(@{$self->{'_hits'}});
    } else {
        return $self->{'_num_hits_significant'} || 0;
    }
}



=head2 lowest_p

 Usage     : $blast->lowest_p()
 Purpose   : Get the lowest P-value among all hits in a BLAST report.
           : Syntactic sugar for $blast->hit('best')->p().
 Returns   : Float or scientific notation number.
           : Returns -1.0 if lowest_p has not been set.
 Argument  : n/a.
 Throws    : Exception if the Blast report does not report P-values
           : (as is the case for NCBI Blast2).
 Comments  : A value is returned regardless of whether or not there were
           : significant hits ($DEFAULT_SIGNIF, currently  999).

See Also   : L<lowest_expect>(), L<lowest_signif>(), L<highest_p>(), L<signif_fmt>()

=cut

#------------
sub lowest_p {
#------------
    my $self = shift;

    # Layout 2 = NCBI Blast 2.x does not report P-values.
    $self->_layout == 2 and
	$self->throw("Can't get P-value with BLAST2.", 
		     "Use lowest_signif() or lowest_expect()");

    return $self->{'_lowestSignif'} || -1.0;
}



=head2 lowest_expect

 Usage     : $blast->lowest_expect()
 Purpose   : Get the lowest Expect value among all hits in a BLAST report.
           : Syntactic sugar for $blast->hit('best')->expect()
 Returns   : Float or scientific notation number.
           : Returns -1.0 if lowest_expect has not been set.
 Argument  : n/a.
 Throws    : Exception if there were no significant hits and the report
           : does not have Expect values on the description lines
           : (i.e., Blast1, WashU-Blast2).

See Also   : L<lowest_p>(), L<lowest_signif>(), L<highest_expect>(), L<signif_fmt>()

=cut

#------------------
sub lowest_expect {
#------------------
    my $self = shift;
    
    if ($self->_layout == 2) {
      return $self->{'_lowestSignif'} || -1.0;
    }

    if($self->{'_is_significant'}) {
	my $bestHit = $self->{'_hits'}->[0];
	return $bestHit->expect();
    } else {
	$self->throw("Can't get lowest expect value: no significant hits ",
		     "The format of this report requires expect values to be extracted$Newline".
		     "from the hits themselves.");
    }
}


=head2 highest_p

 Example   : $blast->highest_p( ['overall'])
 Purpose   : Get the highest P-value among all hits in a BLAST report.
           : Syntactic sugar for $blast->hit('worst')->p()
           : Can also get the highest P-value overall (not just among signif hits).
 Usage     : $p_signif = $blast->highest_p();
           : $p_all    = $blast->highest_p('overall');
 Returns   : Float or scientific notation number.
           : Returns -1.0 if highest_p has not been set.
 Argument  : String 'overall' or no argument.
           : No argument = get highest P-value among significant hits.
 Throws    : Exception if object is created from a Blast2 report
           : (which does not report P-values).

See Also   : L<highest_signif>(), L<lowest_p>(), L<_set_signif>(), L<signif_fmt>()

=cut

#---------------
sub highest_p {
#---------------
    my ($self, $overall) = @_;
    
    # Layout 2 = NCBI Blast 2.x does not report P-values.
    $self->_layout == 2 and
	$self->throw("Can't get P-value with BLAST2.", 
		     "Use highest_signif() or highest_expect()");

    $overall and  return $self->{'_highestSignif'} || -1.0;
    $self->hit('worst')->p();
}



=head2 highest_expect

 Usage     : $blast_object->highest_expect( ['overall'])
 Purpose   : Get the highest Expect value among all significant hits in a BLAST report.
           : Syntactic sugar for $blast->hit('worst')->expect()
 Examples  : $e_sig = $blast->highest_expect();
           : $e_all = $blast->highest_expect('overall');
 Returns   : Float or scientific notation number.
           : Returns -1.0 if highest_exoect has not been set.
 Argument  : String 'overall' or no argument.
           : No argument = get highest Expect-value among significant hits.
 Throws    : Exception if there were no significant hits and the report
           : does not have Expect values on the description lines
           : (i.e., Blast1, WashU-Blast2).

See Also   : L<lowest_expect>(), L<highest_signif>(), L<signif_fmt>()

=cut

#-------------------
sub highest_expect {
#-------------------
    my ($self, $overall) = @_;
    
    if ( $overall and $self->_layout == 2) {
      return $self->{'_highestSignif'} || -1.0;
    }

    if($self->{'_is_significant'}) {
	return $self->hit('worst')->expect;
    } else {
	$self->throw("Can't get highest expect value: no significant hits ",
		     "The format of this report requires expect values to be extracted$Newline".
		     "from the hits themselves.");
    }
}




=head2 lowest_signif

 Usage     : $blast_obj->lowest_signif();
           : Syntactic sugar for $blast->hit('best')->signif()
 Purpose   : Get the lowest P or Expect value among all hits
           : in a BLAST report.
           : This method is syntactic sugar for $blast->hit('best')->signif()
           : The value returned is the one which is reported in the decription
           : section of the Blast report.
           : For Blast1 and WU-Blast2, this is a P-value, 
           : for NCBI Blast2, it is an Expect value.
 Example   : $blast->lowest_signif();
 Returns   : Float or scientific notation number.
           : Returns -1.0 if lowest_signif has not been set.
 Argument  : n/a.
 Throws    : n/a.
 Status    : Deprecated. Use lowest_expect() or lowest_p().
 Comments  : The signif() method provides a way to deal with the fact that
           : Blast1 and Blast2 formats differ in what is reported in the
           : description lines of each hit in the Blast report. The signif()
           : method frees any client code from having to know if this is a P-value
           : or an Expect value, making it easier to write code that can process 
           : both Blast1 and Blast2 reports. This is not necessarily a good thing, since
           : one should always know when one is working with P-values or
           : Expect values (hence the deprecated status).
           : Use of lowest_expect() is recommended since all hits will have an Expect value.

See Also   : L<lowest_p>(), L<lowest_expect>(), L<signif>(), L<signif_fmt>(), L<_set_signif>()

=cut

#------------------
sub lowest_signif {
#------------------
    my ($self) = @_;

    return $self->{'_lowestSignif'} || -1.0;
}



=head2 highest_signif

 Usage     : $blast_obj->highest_signif('overall');
           : Syntactic sugar for $blast->hit('worst')->signif()
 Purpose   : Get the highest P or Expect value among all hits
           : in a BLAST report.
           : The value returned is the one which is reported in the decription
           : section of the Blast report.
           : For Blast1 and WU-Blast2, this is a P-value, 
           : for NCBI Blast2, it is an Expect value.
 Example   : $blast->highest_signif();
 Returns   : Float or scientific notation number.
           : Returns -1.0 if highest_signif has not been set.
 Argument  : Optional  string 'overall' to get the highest overall significance value.
 Throws    : n/a.
 Status    : Deprecated. Use highest_expect() or highest_p().
 Comments  : Analogous to lowest_signif(), q.v.

See Also   : L<lowest_signif>(), L<lowest_p>(), L<lowest_expect>(), L<signif>(), L<signif_fmt>(), L<_set_signif>()

=cut

#---------------------
sub highest_signif {
#---------------------
    my ($self, $overall) = @_;
    
    $overall and return $self->{'_highestSignif'} || -1.0;

    if($self->{'_is_significant'}) {
        my $worst_hit = $self->hit('worst');
        if(defined $worst_hit) {
            return $worst_hit->signif;
        } else {
            return $self->{'_highestSignif'};
        }
    }
}




=head2 matrix

 Usage     : $blast_object->matrix();
 Purpose   : Get the name of the scoring matrix used.
           : This is extracted from the report.
 Argument  : n/a
 Returns   : string or undef if not defined

=cut

#------------
sub matrix { my $self = shift; $self->{'_matrix'} || $Blast->{'_matrix'}; }
#------------


=head2 filter

 Usage     : $blast_object->filter();
 Purpose   : Get the name of the low-complexity sequence filter used.
           : (SEG, SEG+XNU, DUST, NONE).
           : This is extracted from the report.
 Argument  : n/a
 Returns   : string or undef if not defined

=cut

#----------
sub filter { my $self = shift; $self->{'_filter'}  || $Blast->{'_filter'}; }
#----------



=head2 expect

 Usage     : $blast_object->expect();
 Purpose   : Get the expect parameter (E) used for the Blast analysis.
           : This is extracted from the report.
 Argument  : n/a
 Returns   : string or undef if not defined.

=cut

#-----------
sub expect { my $self = shift; $self->{'_expect'} || $Blast->{'_expect'}; }
#-----------



=head2 karlin_altschul

 Usage     : $blast_object->karlin_altschul();
 Purpose   : Get the Karlin_Altschul sum statistics (Lambda, K, H)
           : These are extracted from the report.
 Argument  : n/a
 Returns   : list of three floats (Lambda, K, H)
           : If not defined, returns list of three zeros)

=cut

#---------------------
sub karlin_altschul { 
#---------------------
    my $self = shift; 
    if(defined($self->{'_lambda'})) {
	($self->{'_lambda'}, $self->{'_k'}, $self->{'_h'});
    } elsif(defined($Blast->{'_lambda'})) {
	($Blast->{'_lambda'}, $Blast->{'_k'}, $Blast->{'_h'});
    } else {
	(0, 0, 0);
    }
}



=head2 word_size

 Usage     : $blast_object->word_size();
 Purpose   : Get the word_size used during the Blast analysis.
           : This is extracted from the report.
 Argument  : n/a
 Returns   : integer or undef if not defined.

=cut

#--------------
sub word_size { 
#--------------
    my $self = shift; 
    $self->{'_word_size'} || $Blast->{'_word_size'}; 
}



=head2 s

 Usage     : $blast_object->s();
 Purpose   : Get the s statistic for the Blast analysis.
           : This is extracted from the report.
 Argument  : n/a
 Returns   : integer or undef if not defined.

=cut

#------
sub s { my $self = shift; $self->{'_s'} || $Blast->{'_s'}; }
#------



=head2 gap_creation

 Usage     : $blast_object->gap_creation();
 Purpose   : Get the gap creation penalty used for a gapped Blast analysis.
           : This is extracted from the report.
 Argument  : n/a
 Returns   : integer or undef if not defined.

See Also   : L<gap_extension>()

=cut

#-----------------
sub gap_creation { 
#-----------------
    my $self = shift; 
    $self->{'_gapCreation'} || $Blast->{'_gapCreation'};
}



=head2 gap_extension

 Usage     : $blast_object->gap_extension();
 Purpose   : Get the gap extension penalty used for a gapped Blast analysis.
           : This is extracted from the report.
 Argument  : n/a
 Returns   : integer or undef if not defined.

See Also   : L<gap_extension>()

=cut

#-------------------
sub gap_extension { 
#-------------------
    my $self = shift; 
    $self->{'_gapExtension'} || $Blast->{'_gapExtension'};
}



=head2 ambiguous_aln

 Usage     : $blast_object->ambiguous_aln();
 Purpose   : Test all hits and determine if any have an ambiguous alignment.
 Example   : print "ambiguous" if $blast->ambiguous_aln();
 Returns   : Boolean (true if ANY significant hit has an ambiguous alignment)
 Argument  : n/a
 Throws    : n/a
 Status    : Experimental
 Comments  : An ambiguous BLAST alignment is defined as one where two or more
           : different HSPs have significantly overlapping sequences such
           : that it is not possible to create a unique alignment
           : by simply concatenating HSPs. This may indicate the presence
           : of multiple domains in one sequence relative to another.
           : This method only indicates the presence of ambiguity in at 
           : least one significant hit. To determine the nature of the
           : ambiguity, each hit must be examined.

See Also   : B<Bio::Tools::Blast::Sbjct::ambiguous_aln()>,L<Links to related modules>

=cut

#----------------
sub ambiguous_aln { 
#----------------
    my $self = shift;
    foreach($self->hits()) {
	return 1 if ($_->ambiguous_aln() ne '-');
    }
    0;
}


=head2 overlap

 Usage     : $blast_object->overlap([integer]);
 Purpose   : Set/Get the number of overlapping residues allowed when tiling multiple HSPs.
           : Delegates to Bio::Tools::Blast::Sbjct::overlap().
 Throws    : Exception if there are no significant hits.
 Status    : Experimental

See Also   : B<Bio::Tools::Blast::Sbjct::overlap()>,L<Links to related modules>

=cut

#------------
sub overlap { 
#------------
    my $self = shift; 
    if(not $self->hits) {
	$self->throw("Can't get overlap data without significant hits.");
    }
    $self->hit->overlap();
}


=head2 homol_data

 Usage     : @data = $blast_object->homo_data( %named_params );
 Purpose   : Gets specific similarity data about each significant hit. 
 Returns   : Array of strings:
           : "Homology data" for each HSP is in the format:
           :  "<integer> <start> <stop>"
           : Data for different HSPs are tab-delimited.
 Argument  : named parameters passed along to the hit objects.
 Throws    : n/a
 Status    : Experimental
 Comments  : This is a very experimental method used for obtaining an 
           : indication of:
           :   1) how many HSPs are in a Blast alignment
           :   2) how strong the similarity is between sequences in the HSP
           :   3) the endpoints of the alignment (sequence monomer numbers)

See Also   : B<Bio::Tools::Blast::Sbjct::homol_data()>,L<Links to related modules>

=cut

#----------------
sub homol_data {
#----------------
    
    my ($self, %param) = @_;
    my @hits = $self->hits();
    my @data = ();
    
    ## Note: Homology data can be either for the query sequence or the hit
    ##       (Sbjct) sequence. Default is for sbjct. This is specifyable via
    ##       $param{-SEQ}='sbjct' || 'query'.

    foreach ( @hits ) {
	push @data, $_->homol_data(%param);
    }
    @data;
}




=head1 REPORT GENERATING METHODS


=head2 table

 Usage     : $blast_obj->table( [get_desc]);
 Purpose   : Output data for each HSP of each hit in tab-delimited format.
 Example   : print $blast->table;
           : print $blast->table(0);  
           : # Call table_labels() to print labels.
 Argument  : get_desc = boolean, if false the description of each hit is not included.
           :            Default: true (if not defined, include description column).
 Returns   : String containing tab-delimited set of data for each HSP
           : of each significant hit. Different HSPs are separated by newlines.
           : Left-to-Right order of fields:
           : 1 QUERY_NAME             # Sequence identifier of the query.
           : 2 QUERY_LENGTH           # Full length of the query sequence.
           : 3 SBJCT_NAME             # Sequence identifier of the sbjct ("hit".
           : 4 SBJCT_LENGTH           # Full length of the sbjct sequence.
           : 5 EXPECT                 # Expect value for the alignment.
           : 6 SCORE                  # Blast score for the alignment.
           : 7 BITS                   # Bit score for the alignment.
           : 8 NUM_HSPS               # Number of HSPs (not the "N" value).
           : 9 HSP_FRAC_IDENTICAL     # fraction of identical substitutions.
           : 10 HSP_FRAC_CONSERVED    # fraction of conserved ("positive") substitutions.
           : 11 HSP_QUERY_ALN_LENGTH  # Length of the aligned portion of the query sequence.
           : 12 HSP_SBJCT_ALN_LENGTH  # Length of the aligned portion of the sbjct sequence.
           : 13 HSP_QUERY_GAPS        # Number of gaps in the aligned query sequence.
           : 14 HSP_SBJCT_GAPS        # Number of gaps in the aligned sbjct sequence.
           : 15 HSP_QUERY_START       # Starting coordinate of the query sequence.
           : 16 HSP_QUERY_END         # Ending coordinate of the query sequence.
           : 17 HSP_SBJCT_START       # Starting coordinate of the sbjct sequence.
           : 18 HSP_SBJCT_END         # Ending coordinate of the sbjct sequence.
           : 19 HSP_QUERY_STRAND      # Strand of the query sequence (TBLASTN/X only)
           : 20 HSP_SBJCT_STRAND      # Strand of the sbjct sequence (TBLASTN/X only)
           : 21 HSP_FRAME             # Frame for the sbjct translation (TBLASTN/X only)
           : 22 SBJCT_DESCRIPTION  (optional)  # Full description of the sbjct sequence from 
           :                                  # the alignment section.
 Throws    : n/a
 Comments  : This method does not collect data based on tiling of the HSPs.
           : The table will contains redundant information since the hit name,
           : id, and other info for the hit are listed for each HSP.
           : If you need more flexibility in the output format than this
           : method provides, design a custom function.

See Also  : L<table_tiled>(), L<table_labels>(), L<_display_hits>()

=cut

#-----------
sub table {
#-----------
    my ($self, $get_desc) = @_;
    my $str = '';

    $get_desc = defined($get_desc) ? $get_desc : 1;
#    $str .= $self->_table_labels($get_desc) unless $self->{'_labels'};

    my $sigfmt = $self->signif_fmt();
    $sigfmt eq 'parts' and $sigfmt = 'exp';  # disallow 'parts' format for this table.
    my $sigprint = $sigfmt eq 'exp' ? 'd' : '.1e';

    my ($hit, $hsp);
    foreach $hit($self->hits) {
	foreach $hsp($hit->hsps) {
	    # Note: range() returns a 2-element list.
	    $str .= sprintf "%s\t%d\t%s\t%d\t%$sigprint\t%d\t%d\t%d\t%.2f\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s$Newline", 
                   $self->name, $self->length, $hit->name, $hit->length, 
	           $hit->expect($sigfmt), $hit->score, $hit->bits,
 	           $hit->num_hsps, $hsp->frac_identical, $hsp->frac_conserved, 
	           $hsp->length('query'), $hsp->length('sbjct'), 
	           $hsp->gaps('list'),
	           $hsp->range('query'), $hsp->range('sbjct'), 
	           $hsp->strand('query'), $hsp->strand('sbjct'), $hsp->frame,
	           ($get_desc ? $hit->desc  : '');
	}
    }
    $str =~ s/\t$Newline/$Newline/gs;
    $str;
}



=head2 table_labels

 Usage     : print $blast_obj->table_labels( [get_desc] );
 Purpose   : Get column labels for table().
 Returns   : String containing column labels. Tab-delimited.
 Argument  : get_desc = boolean, if false the description column is not included.
           : Default: true (if not defined, include description column).
 Throws    : n/a

See Also   : L<table>()

=cut

#----------------
sub table_labels {
#----------------
    my ($self, $get_desc) = @_;
    $get_desc = defined($get_desc) ? $get_desc : 1;
    my $descstr = $get_desc ? 'DESC' : '';
    my $descln = $get_desc ? '-----' : '';

    my $str = sprintf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s$Newline", 
                       'QUERY', 'Q_LEN', 'SBJCT', 'S_LEN', 'EXPCT', 'SCORE', 'BITS', 'HSPS', 
                       'IDEN', 'CONSV', 'Q_ALN', 'S_ALN', 'Q_GAP', 'S_GAP',
                       'Q_BEG', 'Q_END', 'S_BEG', 'S_END', 'Q_STR', 'S_STR', 'FRAM', $descstr;
    $str .= sprintf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s$Newline", 
                       '-----', '-----', '-----', '-----', '-----', '-----', '-----', '-----', 
                       '-----', '-----', '-----', '-----', '-----', '-----', 
                       '-----', '-----', '-----','-----', '-----', '-----','-----', $descln;

    $self->{'_labels'} = 1;    
    $str =~ s/\t$Newline/$Newline/gs;
    $str;
}



=head2 table_tiled

 Purpose   : Get data from tiled HSPs in tab-delimited format.
           : Allows only minimal flexibility in the output format.
           : If you need more flexibility, design a custom function.
 Usage     : $blast_obj->table_tiled( [get_desc]);
 Example   : print $blast->table_tiled;
           : print $blast->table_tiled(0);  
           : # Call table_labels_tiled() if you want labels.
 Argument  : get_desc = boolean, if false the description of each hit is not included.
           :            Default: true (include description).
 Returns   : String containing tab-delimited set of data for each HSP
           : of each significant hit. Multiple hits are separated by newlines.
           : Left-to-Right order of fields:
           : 1 QUERY_NAME           # Sequence identifier of the query.
           : 2 QUERY_LENGTH         # Full length of the query sequence.
           : 3 SBJCT_NAME           # Sequence identifier of the sbjct ("hit".
           : 4 SBJCT_LENGTH         # Full length of the sbjct sequence.
           : 5 EXPECT               # Expect value for the alignment.
           : 6 SCORE                # Blast score for the alignment.
           : 7 BITS                 # Bit score for the alignment.
           : 8 NUM_HSPS             # Number of HSPs (not the "N" value).
           : 9 FRAC_IDENTICAL*      # fraction of identical substitutions.
           : 10 FRAC_CONSERVED*     # fraction of conserved ("positive") substitutions .
           : 11 FRAC_ALN_QUERY*     # fraction of the query sequence that is aligned.
           : 12 FRAC_ALN_SBJCT*     # fraction of the sbjct sequence that is aligned.
           : 13 QUERY_ALN_LENGTH*   # Length of the aligned portion of the query sequence.
           : 14 SBJCT_ALN_LENGTH*   # Length of the aligned portion of the sbjct sequence.
           : 15 QUERY_GAPS*         # Number of gaps in the aligned query sequence.
           : 16 SBJCT_GAPS*         # Number of gaps in the aligned sbjct sequence.
           : 17 QUERY_START*        # Starting coordinate of the query sequence.
           : 18 QUERY_END*          # Ending coordinate of the query sequence.
           : 19 SBJCT_START*        # Starting coordinate of the sbjct sequence.
           : 20 SBJCT_END*          # Ending coordinate of the sbjct sequence.
           : 21 AMBIGUOUS_ALN       # Ambiguous alignment indicator ('qs', 'q', 's').
           : 22 SBJCT_DESCRIPTION  (optional)  # Full description of the sbjct sequence from 
           :                                  # the alignment section.
           :
           : * Items marked with a "*" report data summed across all HSPs
           :   after tiling them to avoid counting data from overlapping regions 
           :   multiple times.
 Throws    : n/a
 Comments  : This function relies on tiling of the HSPs since it calls 
           : frac_identical() etc. on the hit as opposed to each HSP individually.

See Also   : L<table>(), L<table_labels_tiled>(), B<Bio::Tools::Blast::Sbjct::"HSP Tiling and Ambiguous Alignments">, L<Links to related modules>

=cut

#----------------
sub table_tiled {
#----------------
    my ($self, $get_desc) = @_;
    my $str = '';

    $get_desc = defined($get_desc) ? $get_desc : 1;

    my ($hit);
    my $sigfmt = $self->signif_fmt();
    $sigfmt eq 'parts' and $sigfmt = 'exp';  # disallow 'parts' format for this table.
    my $sigprint = $sigfmt eq 'exp' ? 'd' : '.1e';

    foreach $hit($self->hits) {
	$str .= sprintf "%s\t%d\t%s\t%d\t%$sigprint\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s$Newline", 
	               $self->name, $self->length, $hit->name, $hit->length, 
	               $hit->expect($sigfmt), $hit->score, $hit->bits, 
	               $hit->num_hsps, $hit->frac_identical, $hit->frac_conserved, 
	               $hit->frac_aligned_query, $hit->frac_aligned_hit,
	               $hit->length_aln('query'), $hit->length_aln('sbjct'), 
                       $hit->gaps('list'), $hit->range('query'), $hit->range('sbjct'),
	               $hit->ambiguous_aln, ($get_desc ? $hit->desc : '');
    }
    $str =~ s/\t$Newline/$Newline/gs;
    $str;
}


=head2 table_labels_tiled

 Usage     : print $blast_obj->table_labels_tiled( [get_desc] );
 Purpose   : Get column labels for table_tiled().
 Returns   : String containing column labels. Tab-delimited.
 Argument  : get_desc = boolean, if false the description column is not included.
           : Default: true (include description column).
 Throws    : n/a

See Also   : L<table_tiled>()

=cut

#---------------------
sub table_labels_tiled {
#---------------------
    my ($self, $get_desc) = @_;
    my $descstr = $get_desc ? 'DESC' : '';
    my $descln = $get_desc ? '-----' : '';

    my $str = sprintf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s$Newline", 
                       'QUERY', 'Q_LEN', 'SBJCT', 'S_LEN', 'EXPCT', 'SCORE', 'BITS',
                       'HSPS', 'FR_ID', 'FR_CN', 'FR_ALQ', 'FR_ALS', 'Q_ALN', 
                       'S_ALN', 'Q_GAP', 'S_GAP', 'Q_BEG', 'Q_END', 'S_BEG', 'S_END', 
                       'AMBIG', $descstr;
    $str =~ s/\t$Newline/$Newline/;
    $str .= sprintf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s$Newline", 
                       '-----', '-----', '------', '-----', '-----','-----', '-----',
                       '-----', '-----', '-----', '-----', '-----', '-----',
                       '-----', '-----', '-----','-----','-----', '-----', 
                       '-----','-----', $descln;

    $self->{'_labels_tiled'} = 1;    
    $str =~ s/\t$Newline/$Newline/gs;
    $str;
}




=head2 display

 Usage     : $blast_object->display( %named_parameters );
 Purpose   : Display information about Bio::Tools::Blast.pm data members,
           : E.g., parameters of the report, data for each hit., etc.
           : Overrides Bio::Root::Object::display().
 Example   : $object->display(-SHOW=>'stats');
           : $object->display(-SHOW=>'stats,hits');
 Argument  : Named parameters: (TAGS CAN BE UPPER OR LOWER CASE)
           :     -SHOW  => 'file' | 'hits' | 'homol'
           :     -WHERE => filehandle (default = STDOUT)
 Returns   : n/a (print/printf is called)
 Status    : Experimental
 Comments  : For tab-delimited output, see table().

See Also   : L<_display_homol>(), L<_display_hits>(), L<_display_stats>(), L<table>(), B<Bio::Root::Tools::SeqAnal::display()>,L<Links to related modules>, 

=cut

#--------------
sub display {
#--------------
    my( $self, %param ) = @_;
    
    $self->SUPER::display(%param);
    my $OUT = $self->fh();

    $self->show =~ /homol/i and $self->_display_homol($OUT);
    $self->show =~ /hits/i and $self->_display_hits( %param );
    1;
}



=head2 _display_homol

 Usage     : n/a; called automatically by display()
 Purpose   : Print homology data for hits in the BLAST report.
 Example   : n/a
 Argument  : one argument = filehandle object.
 Returns   : printf call.
 Status    : Experimental

See Also   : L<homol_data>(), L<display>()

=cut

#-------------------
sub _display_homol {
#-------------------
    my( $self, $OUT ) = @_;
    
    print $OUT "${Newline}BLAST HOMOLOGY DATA FOR: ${\$self->name()}$Newline";
    print $OUT '-'x40,"$Newline";
    
    foreach ( $self->homol_data()) {
	print $OUT "$_$Newline";
    }
}


=head2 _display_stats

 Usage     : n/a; called automatically by display()
 Purpose   : Display information about the Blast report "meta" data.
           : Overrides Bio::Tools::SeqAnal::_display_stats() calling it first.
 Example   : n/a
 Argument  : one argument = filehandle object.
 Returns   : printf call.
 Status    : Experimental

See Also   : L<display>(), B<Bio::Tools::SeqAnal::_display_stats()>,L<Links to related modules>

=cut

#--------------------
sub _display_stats {
#--------------------
    my( $self, $OUT ) = @_;
    
    $self->SUPER::_display_stats($OUT);
    printf( $OUT "%-15s: %s$Newline", "GAPPED", $self->gapped ? 'YES' : 'NO');
    printf( $OUT "%-15s: %d$Newline", "TOTAL HITS", $self->num_hits('total'));
    printf( $OUT "%-15s: %s$Newline", "CHECKED ALL", $Blast->{'_check_all'} ? 'YES' : 'NO');
    printf( $OUT "%-15s: %s$Newline", "FILT FUNC", $Blast->{'_filt_func'} ? 'YES' : 'NO');
    if($self->min_length) {
	printf( $OUT "%-15s: Length >= %s$Newline", "MIN_LENGTH", $self->min_length);
    }

    my $num_hits =  $self->num_hits;
    my $signif_str = ($self->_layout == 1) ? 'P' : 'EXPECT';

    printf( $OUT "%-15s: %d$Newline", "SIGNIF HITS", $num_hits);
    # Blast1: signif = P-value, Blast2: signif = Expect value.
	
    printf( $OUT "%-15s: %s ($signif_str-VALUE)$Newline", "SIGNIF CUTOFF", $self->signif);
    printf( $OUT "%-15s: %s$Newline", "LOWEST $signif_str", $self->lowest_signif());
    printf( $OUT "%-15s: %s$Newline", "HIGHEST $signif_str", $self->highest_signif());

    printf( $OUT "%-15s: %s (OVERALL)$Newline", "HIGHEST $signif_str", $self->highest_signif('overall'));
    

    if($Blast->_get_stats) {
	my $warn = ($Blast->{'_share'}) ? '(SHARED STATS)' : '';
	printf( $OUT "%-15s: %s$Newline", "MATRIX", $self->matrix() || 'UNKNOWN');
	printf( $OUT "%-15s: %s$Newline", "FILTER", $self->filter() || 'UNKNOWN');
	printf( $OUT "%-15s: %s$Newline", "EXPECT", $self->expect() || 'UNKNOWN');
	printf( $OUT "%-15s: %s, %s, %s %s$Newline", "LAMBDA, K, H", $self->karlin_altschul(), $warn);
	printf( $OUT "%-15s: %s$Newline", "WORD SIZE", $self->word_size() || 'UNKNOWN');
	printf( $OUT "%-15s: %s %s$Newline", "S", $self->s() || 'UNKNOWN', $warn);
	if($self->gapped) {
	    printf( $OUT "%-15s: %s$Newline", "GAP CREATION", $self->gap_creation() || 'UNKNOWN');
	    printf( $OUT "%-15s: %s$Newline", "GAP EXTENSION", $self->gap_extension() || 'UNKNOWN');
	}
    }
    print $OUT "$Newline";
}


=head2 _display_hits

 Usage     : n/a; called automatically by display()
 Purpose   : Display data for each hit. Not tab-delimited.
 Example   : n/a
 Argument  : one argument = filehandle object.
 Returns   : printf call.
 Status    : Experimental
 Comments  : For tab-delimited output, see table().

See Also   : L<display>(), B<Bio::Tools::Blast::Sbjct::display()>, L<table>(),  L<Links to related modules>

=cut

sub _display_hits {
    
    my( $self, %param ) = @_;
    my $OUT = $self->fh();
    my @hits  = $self->hits();
    
    ## You need a wide screen to see this properly.
    # Header.
    print $OUT "${Newline}BLAST HITS FOR: ${\$self->name()} length = ${\$self->length}$Newline";
    print "(This table requires a wide display.)$Newline";
    print $OUT '-'x80,"$Newline";

    print $self->table_labels_tiled(0);
    print $self->table_tiled(0);

    ## Doing this interactively since there is potentially a lot of data here.
    ## Not quite satisfied with this approach.

    if (not $param{-INTERACTIVE}) {
	return 1;
    } else {
	my ($reply);
	print "${Newline}DISPLAY FULL HSP DATA? (y/n): [n] ";
	chomp( $reply = <STDIN> );
	$reply =~ /^y.*/i;
	
	my $count = 0;
	foreach ( @hits ) {
	    $count++;
	    print $OUT "$Newline$Newline",'-'x80,"$Newline";
	    print $OUT "HSP DATA FOR HIT #$count  (hit <RETURN>)";
		print $OUT "$Newline",'-'x80;<STDIN>;
	    $param{-SHOW} = 'hsp';
	    $_->display( %param );
	}
    }
    1;
}


=head2 to_html

 Usage     : $blast_object->to_html( [%named_parameters] )
 Purpose   : To produce an HTML-formatted version of a BLAST report
           : for efficient navigation of the report using a web browser.
 Example   : # Using the static Blast object:
           : # Can read from STDIN or from a designated file:
           :   $Blast->to_html($file); 
           :   $Blast->to_html(-FILE=>$file, -HEADER=>$header);
           :   (if no file is supplied, STDIN will be used).
           : # saving HTML to an array:
           :   $Blast->to_html(-FILE=>$file, -OUT =>\@out);
           : # Using a pre-existing blast object (must have been built from
           : # a file, not STDIN:
           :   $blastObj->to_html();  
 Returns   : n/a, either prints report to STDOUT or saves to a supplied array
           : if an '-OUT' parameter is defined (see below).
 Argument  : %named_parameters: (TAGS ARE AND CASE INSENSITIVE).
           :    -FILE   => string containing name of a file to be processed.
           :               If not a valid file or undefined, STDIN will be used.
           :               Can skip the -FILE tag if supplying a filename 
           :               as a single argument.
           :    -HEADER => string 
           :               This should be an HTML-formatted string to be used 
           :               as a header for the page, typically describing query sequence,
           :               database searched, the date of the analysis, and any
           :               additional links. 
           :               If not supplied, no special header is used.
           :               Regardless of whether a header is supplied, the
           :               standard info at the top of the report is highlighted.
           :               This should include the <HEADER></HEADER> section 
           :               of the page as well.
           :
           :    -IN    => array reference containing a raw Blast report.
           :              each line in a separate element in the array. 
           :              If -IN is not supplied, read() is called
           :              and data is then read either from STDIN or a file.
           :
           :    -OUT   => array reference to hold the HTML output.
           :              If not supplied, output is sent to STDOUT.
 Throws    : Exception is propagated from $HTML::get_html_func()
           : and Bio::Root::Object::read().
 Comments  : The code that does the actual work is located in
           :  Bio::Tools::Blast::HTML::get_html_func().
 Bugs      : Some hypertext links to external databases may not be
           : correct. This due in part to the dynamic nature of
           : the web.
           : Hypertext links are not added to hits without database ids.
 TODO      : Possibly create a function to produce fancy default header
           : using data extracted from the report (requires some parsing).
           : For example, it would be nice to always include a date

See Also   : B<Bio::Tools::Blast::HTML::get_html_func()>, B<Bio::Root::Object::read()>, L<Links to related modules>

=cut

#------------
sub to_html {
#------------
    my ($self, @param) = @_;
    
    # Permits syntax such as: $blast->to_html($filename);
    my ($file, $header_html, $in_aref, $out_aref) = 
	$self->_rearrange([qw(FILE HEADER IN OUT)], @param);

    $self->file($file) if $file;

    # Only setting the newline character once for efficiency.
    $Newline ||= $Util->get_newline(-client => $self, @param);

    $header_html ||= '';  
    (ref($out_aref) eq 'ARRAY') ? push(@$out_aref, $header_html) : print "$header_html$Newline";

    require Bio::Tools::Blast::HTML; 
    Bio::Tools::Blast::HTML->import(qw(&get_html_func)); 

    my ($func);
    eval{ $func = &get_html_func($out_aref);  };
    if($@) {
	my $err = $@; 
	$self->throw($err);
    }

    eval {
	if(!$header_html) {
	    $out_aref ? push(@$out_aref, "<html><body>$Newline") : print "<html><body>$Newline";
	}

	if (ref ($in_aref) =~ /ARRAY/) {
	    # If data is being supplied, process it.
	    foreach(@$in_aref) {
		&$func($_);
	    }
	} else {
	    # Otherwise, read it, processing as we go.
	    
  	    $self->read(-FUNC => $func, @param);
	}    
	$out_aref ? push(@$out_aref, "$Newline</pre></body></html>") : print "$Newline</pre></body></html>";
    };

    if($@) {
	# Check for trivial error (report already HTML formatted).
	if($@ =~ /HTML formatted/) {
	    print STDERR "\a${Newline}Blast report appears to be HTML formatted already.$Newline$Newline";
	} else {
	    my $err = $@; 
	    $self->throw($err);
	}
    }
}



1;
__END__

#####################################################################################
#                                END OF CLASS                                       #
#####################################################################################

=head1 FOR DEVELOPERS ONLY

=head2 Data Members

Information about the various data members of this module is provided for those 
wishing to modify or understand the code. Two things to bear in mind: 

=over 4

=item 1 Do NOT rely on these in any code outside of this module. 

All data members are prefixed with an underscore to signify that they are private.
Always use accessor methods. If the accessor doesn't exist or is inadequate, 
create or modify an accessor (and let me know, too!). (An exception to this might
be for Sbjct.pm or HSP.pm which are more tightly coupled to Blast.pm and
may access Blast data members directly for efficiency purposes, but probably 
should not).

=item 2 This documentation may be incomplete and out of date.

It is easy for these data member descriptions to become obsolete as 
this module is still evolving. Always double check this info and search 
for members not described here.

=back

An instance of Bio::Tools::Blast.pm is a blessed reference to a hash containing
all or some of the following fields:

 FIELD           VALUE
 --------------------------------------------------------------
 _significance    P-value or Expect value cutoff (depends on Blast version:
	          Blast1/WU-Blast2 = P-value; Blast2 = Expect value).
   	          Values GREATER than this are deemed not significant.

 _significant     Boolean. True if the query has one or more significant hit.

 _min_length      Integer. Query sequences less than this will be skipped.

 _confirm_significance  Boolean. True if client has supplied significance criteria.

 _gapped          Boolean. True if BLAST analysis has gapping turned on.

 _hits            List of Sbjct.pm objects. 

 _num_hits        Number of hits obtained from the BLAST report.

 _num_hits_significant Number of significant based on Significant data members.

 _highestSignif   Highest P or Expect value overall (not just what is stored in _hits).

 _lowestSignif    Lowest P or Expect value overall (not just what is stored in _hits).


The static $Blast object has a special set of members:

  _errs
  _share
  _stream
  _get_stats
  _gapped
  _filt_func

 Miscellaneous statistical parameters:
 -------------------------------------
  _filter, _matrix, _word_size, _expect, _gapCreation, _gapExtension, _s,
  _lambda, _k, _h
  

 INHERITED DATA MEMBERS 
 -----------------------
 (See Bio::Tools::SeqAnal.pm for inherited data members.)

=cut

1;
