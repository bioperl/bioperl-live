#-------------------------------------------------------------------------------
# PACKAGE : Bio::Tools::Blast::HTML.pm
# PURPOSE : To encapsulate code for HTML formatting BLAST reports.
# AUTHOR  : Steve Chervitz (sac@bioperl.org)
# CREATED : 28 Apr 1998
# STATUS  : Alpha
# REVISION: $Id$
# 
# For the latest version and documentation, visit the distribution site:
#    http://bio.perl.org/Projects/Blast/
#
# To generate documentation, run this module through pod2html
# (preferably from Perl v5.004 or better).
#
# CUSTOMIZATION NOTE:
#
#   If your Blast reports are not getting marked up correctly, add or
#   modify the regexps in _markup_report() to accomodate the format of
#   your reports.
#
# Copyright (c) 1996-98 Steve Chervitz. All Rights Reserved.
#           This module is free software; you can redistribute it and/or 
#           modify it under the same terms as Perl itself.
#-------------------------------------------------------------------------------

package Bio::Tools::Blast::HTML;
use strict;
use Exporter;

use Bio::Tools::WWW  qw(:obj); 
use Carp;

use vars qw( @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS
             $ID %DbUrl %SGDUrl $Revision
	     $Acc $Pir_acc $Word $Signif $Int $Descrip);

@ISA        = qw(Exporter);
@EXPORT     = qw();
@EXPORT_OK  = qw(&get_html_func &strip_html);
%EXPORT_TAGS = ( std => [qw(&get_html_func  &strip_html)] );

$ID = 'Bio::Tools::Blast::HTML';
$Revision = '$Id$';  #'

my $_set_markup = 0;
my $_gi_link = '';


## POD Documentation:

=head1 NAME

Bio::Tools::Blast::HTML.pm - Bioperl Utility module for HTML
formatting Blast reports

=head1 SYNOPSIS

=head2 Adding HTML-formatting

    use Bio::Tools::Blast::HTML qw(&get_html_func);

    $func = &get_html_func();

    # Now as each line of the report is read, pass it to &$func($line).

See L<get_html_func>() for details.
Also see B<Bio::Tools::Blast::to_html> for an example of usage.


=head2 Removing HTML-formatting

    use Bio::Tools::Blast::HTML qw(&strip_html);

    &strip_html(\$blast_report_string)

See L<strip_html>() for details.


=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

=head1 DESCRIPTION

This module can be used to add HTML formatting to or remove HTML
formatting from a raw Blast sequence analysis report. Hypertext links
to the appropriate database are added for each hit sequence (GenBank,
Swiss-Prot, PIR, PDB, SGD).

This module is intended for use by Bio::Tools::Blast.pm and related modules, 
which provides a front-end to the methods in Bio::Tools::Blast::HTML.pm.

=head1 DEPENDENCIES

Bio::Tools::Blast::HTML.pm does not inherit from any other class
besides Exporter.  It is used by B<Bio::Tools::Blast.pm> only.  This
class relies on B<Bio::Tools::WWW.pm> to provide key URLS for adding
links in the Blast report to specific databases.

The greatest dependency comes from the dynamic state of the web. URLs
are are likely to change in the future, so all links cannot be
guaranteed to work indefinitely.  Feel free to report broken or
incorrect database links (L<FEEDBACK>). Thanks!

=head1 SEE ALSO

 Bio::Tools::Blast.pm    - Blast object.
 Bio::Tools::WWW.pm      - URL repository.

 http://bio.perl.org/Projects/modules.html  - Online module documentation
 http://bio.perl.org/Projects/Blast/        - Bioperl Blast Project     
 http://bio.perl.org/                       - Bioperl Project Homepage

=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

    bioperl-l@bioperl.org          - General discussion
    http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep 
track the bugs and  their resolution. Bug reports can be submitted 
via email or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR

Steve Chervitz, sac@bioperl.org

=head1 COPYRIGHT

Copyright (c) 1998-2000 Steve Chervitz. All Rights Reserved.
This module is free software; you can redistribute it and/or 
modify it under the same terms as Perl itself.


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



=head2 get_html_func

 Usage     : $func_ref = &get_html_func( [array_ref] );
           : This method is exported.
 Purpose   : Provides a function that adds HTML formatting to a
           : raw Blast report line-by-line.
           : Utility method used by to_html() in Bio::Tools::Blast.pm.
 Returns   : Reference to an anonymous function to be used while reading in  
           : the raw report. 
           : The function itself operates on the Blast report line-by-line
           : HTML-ifying it and printing it to STDOUT (or saving in the supplied
           : array ref) as it goes:
           :     foreach( @raw_report ) { &$func_ref($_); }
 Argument  : array ref (optional) for storing the HTML-formatted report.
           : If no argument is supplied, HTML output is sent to STDOUT.
 Throws    : Croaks if an argument is supplied and is not an array ref.
           : The anonymous function returned by this method croaks if 
           : the Blast output appears to be HTML-formatted already.
 Comments  : Adapted from a script by Keith Robison  November 1993 
           : krobison@nucleus.harvard.edu; http://golgi.harvard.edu/gilbert.html
           : Modified extensively by Steve Chervitz and Mike Cherry.
           : Some modifications are customizations for BLAST reports served up
           : by the Saccharomyces Genome Database.
           : Feel free to modify or replace portions of this code as necessary
           : to accomodate new BLAST datasets or changes to the Blast format.

See Also   : B<Bio::Tools::Blast::to_html()>

=cut

#--------------------
sub get_html_func {
#--------------------
    my ($out_aref) = @_;

    ## Key booleans used in parsing.
    my $found_table = 0;  # Located the table at top of report (a.k.a. 'descriptions').
    my $found_data  = 0;  # Nothing is done until this is true
    my $skip        = 0;  # Skipping various items in the report header
    my $ref_skip    = 0;  # so we can include nice HTML versions 
                          # (e.g., references for the BLAST program). 
    my $getNote     = 0; 
    my $getGenBankAlert = 0; 
    my $str = '';
    my $gi_link = \$_gi_link;
    my $prog = '';

    if( defined($out_aref) and not ref($out_aref) eq 'ARRAY') {
	croak("Argument must be an ARRAY ref not a ${\ref $out_aref}.");
    }

    my $refs = &_prog_ref_html;

    &_set_markup_data() if not $_set_markup;

    return sub {
	# $_ contains a single line from a Blast report.
	local $_ = shift;

	croak("Report appears to be HTML formatted already.") if m/<HTML>|<TITLE>|<PRE>/i;

	if(not $found_table) {
	    if($ref_skip) {
		# Replacing an reference data with special HTML.
		$ref_skip = 0 if /^\s+$/;
	    }
	    if($getNote) { 
		## SAC: created this test since we are no longer reading from STDIN.
		$out_aref ? push(@$out_aref, $_) : print $_;
		$getNote = 0 if m/^\s+$/;
	    } elsif( m/(.*), Up \d.*/ or /Date: +(.+)/ or /Start: +(.+?) +End:/ ) {
		### Network BLAST reports from NCBI are time stamped as follows:
		#Fri Apr 18 15:55:41 EDT 1997, Up 1 day, 19 mins, 1 user, load: 19.54, 19.13, 17.77
		my $date = "<b>BLASTed on:</b> $1<p>\n";
		$out_aref ? push(@$out_aref, $date) : print $date;
	    } elsif ( /^(<\w+>)?(T?BLAST[NPX])\s+(.*?)/ ) {
		$found_data = 1;
		local($^W) = 0;
		s#(\S+)\s+(.*)#<P><B>Program:</B> $1 $2 $3<br>#o;
		$out_aref ? push(@$out_aref, $_) : print $_;
		$skip = 1;
		$prog = $2;
		if($prog =~ /BLASTN/) {
		    ## Prevent the error at Entrez when you ask for a nucl 
		    ## entry with a protein GI number.
		    $$gi_link = $DbUrl{'gb_n'};  # nucleotide
		}  else {
		    $$gi_link = $DbUrl{'gb_p'};  # protein
		}
	    } elsif ( m/^Query=/ ) {
		# Keeping the "Query=" format to keep it parsable by Blast.pm
		# (after stripping HTML).
		s#Query= *(.*)#<title>$1</title>\n<p><b>Query=</b> $1#o;
		$out_aref ? push(@$out_aref, $_) : print $_;
		$skip = 1;
	    } elsif ( /Reference:/) {
		$ref_skip = 1;
	    } elsif ( /^Database:/ ) {
		&_markup_database(\$_);
		$out_aref ? push(@$out_aref, $_) : print $_;
		if ( /non-redundant genbank/i and $prog =~ /TBLAST[NX]/i) { 
		    $getGenBankAlert = 1; 
		}
		$skip = 1;
	    } elsif ( /sequences;/ ) {
		$str = "$_<p>";
		$out_aref ? push(@$out_aref, $str) : print $str;
	    } elsif ( /^\s+\(\d+ letters\)\s+/ ) {
		$str = "<br>&nbsp&nbsp&nbsp&nbsp$_";
		$out_aref ? push(@$out_aref, $str) : print $str;
	    } elsif ( /^(WARNING|NOTICE):/i ) {
		s#WARNING: *(.*)#<p><b><font color="red">$1:</font></b> $1#o;
		$out_aref ? push(@$out_aref, $_) : print $_;
		$getNote = 1;
	    } elsif ( /Score +E\s*$/ or /Probability\s*$/ ) {
		# Put the last HTML-formatted lines before the main body of report.
		$found_table = 1;
		$skip = 0;
		$out_aref ? push(@$out_aref, $refs) : print $refs; 
		if($getGenBankAlert) { 
		    $str = &_genbank_alert;
		    $out_aref ? push(@$out_aref, $str) : print $str; 
		}
		$str = "\n<p><pre>";
		$out_aref ? push(@$out_aref, $str) : print $str;
	    }

	} else {
	    &_markup_report(\$_);
	}

	if ($found_data and not($skip or $ref_skip)) {
	    $out_aref ? push(@$out_aref, $_) : print $_;
	}
	1;
    } # end sub {}
}




=head2 _set_markup_data

 Usage     : n/a; utility method used by get_html_func()
 Purpose   : Sets various hashes and regexps used for adding HTML
           : to raw Blast output.
 Returns   : n/a
 Comments  : These items need be set only once. 

See Also   : L<get_html_func>()

=cut

#-------------------
sub _set_markup_data {
#-------------------
    %DbUrl      = $BioWWW->search_url('all');
    %SGDUrl     = $BioWWW->sgd_url('all');

    $Signif  = '[\de.-]{3,}';        # Regexp for a P-value or Expect value. 
    $Int     = ' *\d\d*';            # Regexp for an integer.
    $Descrip = ' +.* {2,}?';         # Regexp for a description line.
    $Acc     = '[A-Z][\d.]+';        # Regexp for GB/EMBL/DDJB/SP accession number
    $Pir_acc = '[A-Z][A-Z0-9]{5,}';  # Regexp for PIR accession number
    $Word    = '[\w_.]+';            # Regexp for a word. Include dot for version.
    
    $_set_markup = 1;
}


=head2 _markup_database

 Usage     : n/a; utility method used by get_html_func()
 Purpose   : Converts a cryptic database ID into a readable name.
 Returns   : n/a
 Comments  : This is used for converting local database IDs into
           : understandable terms. At present, it only recognizes
           : databases used locally at SGD. 

See Also   : L<get_html_func>()

=cut

#---------------------
sub _markup_database {
#---------------------
    my $line_ref = shift;
    local $_ = $$line_ref;

    $_ =~ s#YeastN#<i>S. cerevisiae</i> GenBank Data Set; #;
    $_ =~ s#YeastP#Non-Redundant <i>S. cerevisiae</i> Protein Data Set; #;
    $_ =~ s#genoSC#Complete DNA Sequence for the S. cerevisiae Genome; #;
    $_ =~ s#YeastORF-P#Translation of all Standard S.c. ORFs; #;
    $_ =~ s#YeastORF-N#Coding Sequence of all Standard S.c. ORFs; #;
    s#Database: *(.*)#<p><b>Database:</b> $1#o;

    $$line_ref = $_;
}


=head2 _markup_report

 Usage     : n/a; utility function used by get_html_func()
 Purpose   : Adds HTML links to aid navigation of raw Blast output.
 Returns   : n/a
 Comments  : HTML-formatting is dependent on the Blast server that
           : provided the Blast report. Currently, this function can handle reports
           : produced by NCBI and SGD. Feel free to modify this function
           : to accomodate reports produced by other servers/sites.
           :
           : This function is simply a collection of substitution regexps 
           : that recognize and modify the relevant lines of the Blast report. 
           : All non-header lines of the report are passed through this function,
           : only the ones that match will get modified.
           :
           : The general scheme for adding links is as follows:
           : (Some of the SGD markups do not follow this scheme precisely
           :  but this is the general trend.)
           :
           : For description lines in the summary table at the top of report:
           :
           : DB:SEQUENCE_ID  DESCRIPTION   SIGNIF_VAL
           :        DB          = links to the indicated database (if not Gen/Embl/Ddbj).
           :        SEQUENCE_ID = links to GenBank entry for the sequence.
           :        SIGNIF_VAL  = internal link to relevant alignment section.
           :
           : For the alignment sections in the body of the report:
           :
           : DB:SEQUENCE_ID  (Back | Top) DESCRIPTION 
           :        DB          = links to the indicated database (if not Gen/Embl/Ddbj).
           :        SEQUENCE_ID = links to GenBank entry for the sequence.
           :        SIGNIF_VAL  = internal link to alignment section.
           :        Back        = internal link to description line in summary section.
           :        Top         = internal link to top of page.
           :
           : 'DB' links are created for PDB, PIR, and SwissProt sequences.
           :
           : RE_PARSING HTML-FOMRATTED REPORTS:
           : ----------------------------------
           : HTML-formatted reports generated by this module, as well as reports
           : obtained from the NCBI servers, should be parsable
           : by Bio::Tools::Blast.pm. Parsing HTML-formatted reports is
           : slow, however, since the HTML must be removed prior to parsing.
           : Parsing HTML-formatted reports is dependent on the specific structure
           : of the HTML and is generally not recommended.
           : 
           : Note that since URLs can change without notice, links will need updating.
           : The links are obtained from Bio::Tools::WWW.pm updating that module
           : will update this as well.
           :
 Bugs      : Some links to external databases are incorrect
           : (in particular, for 'bbs' and 'prf' databases on NCBI Blast reports.
           : Some links may fail as a result of the dynamic nature of the web.
           : Hypertext links are not added to hits without database ids.

See Also   : L<get_html_func>(), B<Bio::Tools::WWW.pm>, L<strip_html>()

=cut

#--------------------
sub _markup_report {
#--------------------
    my $line_ref = shift;
    local $_ = $$line_ref;
##
## REGEXPS FOR ALIGNMENT SECTIONS (within the body of the report, 
##                                 the text above the list of HSPs).
##
## If the HSP alignment sections don't start with a '>' we have no way
## of finding them. This occurs with reports saved from HTML-formatted 
## web pages, which we shouldn't be processing here anyway.

## To facilitate parsing of HTML-formatted reports by Bio::Tools::Blast.pm,
## the <a name=...> anchors should be added at the BEGINNING of the HSP 
## alignment section lines and at the END of the description section lines.

    # Removing " ! " addded by GCG. 
    s/ ! / /;

    ### NCBI-specific markups for HSP alignment section lines:

    local($^W) = 0;

  # GenBank/EMBL, DDBJ hits (GenBank Format):
  s@^>(gb|emb|dbj|ref)\|($Word)(\|$Word)?(.*)$@<a name=$2_A></a><b>$1:<a href="$_gi_link$2">$2$3</a></b>$4<br>(<a href="\#$2_H">Back|<a href="\#top">Top</a>)@o;

  s@^>(gb|emb|dbj|ref)\|($Word)(\| \(?$Word\)?)(.*)$@<a name=$2_A></a><b>$1:<a href="$_gi_link$2">$2</a></b>$3$4<br>(<a href="\#$2_H">Back|<a href="\#top">Top</a>)@o;

  # PIR hits
  s@^>pir\|\|($Word)( .*)$@<a name=$1_A></a><b><a href=\"$DbUrl{'pir_acc'}$1\">pir</a>:<a href="$DbUrl{'gb_p'}$1">$1</a></b> $2 <br>(<a href="\#$1_H">Back|<a href="\#top">Top</a>)@o;

  # GI hits (GenBank Format):  using a nested (())
  s@^>(gi)\|($Word)( +\(($Word)\))( .*)$@<a name=$4_A></a><b>$1:<a href="$_gi_link$4">$2</a></b>$3$5<br>(<a href="\#$4_H">Back|<a href="\#top">Top</a>)@o;

  # GNL PID hits (GenBank Format):
  s@^>(gnl)\|($Word)?(\|$Word) +\(($Word)\)( .*)$@<a name=$4_A></a><b>$1:<a href="$_gi_link$4">$2$3</a></b>($4)$5<br>(<a href="\#$4_H">Back|<a href="\#top">Top</a>)@o;

  # BBS and PRF hits (what db?) (GenBank Format):
  s@^>(bbs|prf)\|\|?($Word)( .*)$@<a name=$2_A></a><b>$1:<a href="$_gi_link$2">$2</a></b>$3<br>(<a href="\#$2_H">Back|<a href="\#top">Top</a>)@o;

    # SwissProt hits:
  s@^>sp\|($Word)\|($Word)?( .*)$@<a name=$1_A></a><b><a href="$DbUrl{'swpr'}$1">sp</a>:<a href="$DbUrl{'gb_p'}$1">$1|$2</a></b>$3<br>(<a href="\#$1_H">Back|<a href="\#top">Top</a>)@o;


  ## PDB ids with or without a chain identifier (GenBank format)
  s@^>pdb\|(\d\w{3})\|[\w ] (.*)$@<a name=$1_A></A><b><a href=\"$DbUrl{'3db'}$1\">pdb</A>:<a href="$DbUrl{'gb_struct'}$1">$1</a></b> (<a href="\#$1_H">Back</a>|<a href="\#top">Top</a>)  $2@o;


    ### SGD-specific markups for HSP alignment section lines:

  ## PDB ids without chain identifier
  s@^>PDB_UNIQUEP:(\d\w{3})_ (.*)$@<a name=$1_A></A><b><A HREF="$DbUrl{'3db'}$1">PDB</a>:<A HREF="$DbUrl{'gb_struct'}$1">$1</A></b> (<a href="\#$1_H">Back</a>|<a href="\#top">Top</a>)  $2@o;

  ## PDB ids with chain identifier
  s@^>PDB_UNIQUEP:(\d\w{3})_([\w ]{1})(.*)$@<a name=$1_A></A><b><A HREF="$DbUrl{'3db'}$1">PDB</a>:<A HREF="$DbUrl{'gb_struct'}$1">$1</A></b> Chain:$2, (<a href="\#$1_H">Back</a>|<a href="\#top">Top</a>)  $3@o;

  s@^>($Word)PEPT:GI_(\d+)(.*)$@<a name=$2_A></a><b>$1:<a href="$DbUrl{'gb_p'}$2">GI_$2</a></b> $3 <br>(<a href="\#$2_H">Back|<a href="\#top">Top</a>)@o;

  # The gcg blast dataset generating tools up-case all sbjct sequence IDs.
  # This is fine for yeast but not worm. This is considered a hack here.
  s@WORMPEPT:(\w+\.)(\S+)@WORMPEPT:$1\L$2\E@;

  s@^>WORMPEPT:(\S+)(.*)$@<a name=$1_A></a><b>WORMPEP:<A HREF="$DbUrl{'wormace'}$1">$1</a></b> $2 <br>(<a href="\#$1_H">Back|<a href="\#top">Top</a>)@o;

  s#^>(GB_$Word):($Word) ($Acc) (.*$)#<a name=$2_$3_A></A><a href=\#$2_$3_H>$2|$3</A>$4\t<b>[<A HREF=$_gi_link$3>GenBank</A> / <A HREF=$DbUrl{'embl'}$3>EMBL</A> / <A HREF=\"$SGDUrl{'seq_an'}$2\*\">SGD</A>]</b> #o;

   # Sac's version: ORF name is an external link into SGD:
  s@^>ORFP:(\S*) +([\w-]+)(.*$)@<a name=$1_A></A>ORFP:<a href=\"$SGDUrl{'locus'}$2\">$1 $2</A>$3<br>&nbsp&nbsp&nbsp&nbsp&nbsp<b>[<A HREF=\"$SGDUrl{'seq_an'}$2\">Gene/Sequence Resources</a> / <a href=\"$SGDUrl{'map_orf'}$2\">ORF Map</a></b>] <a href="\#$1_H">Back</a>|<a href="\#top">Top</a>@o;

# Mike's version:
#  s#^>ORFP:(\S*) (.*$)#<a name=$1_A></A><a href=\#$1_H>ORFP:$1</A> $2\t<b>[<A HREF=\"$SGDUrl{'seq_an'}$1\">Gene/Sequence Resources</a> / <a href=\"$SGDUrl{'map_orf'}$1\">ORF Map</a>]</b> #o;

  s#^>ORFN:(\S*) (.*$)#<a name=$1_A></A><a href=\#$1_H>ORFN:$1</A> $2\t<b>[<A HREF=\"$SGDUrl{'seq_an'}$1\">Gene/Sequence Resources</a>] / <a href=\"$SGDUrl{'map_orf'}$1\">ORF Map</a></b> #o;

  s#^>NR_SC:GP-\S* gi\|(\w+)([\w\|]*) (.*$)#<a name=$1_A></A><a href=\#$1_H>GenPept|$1</A> gp|$2 $3\t<b>[<A HREF=$DbUrl{'gb_p'}$1>GenPept</A> / <A HREF=\"$SGDUrl{'gi'}$1\*\">SGD</A>]</b> #o;

  s#^>NR_SC:SW-$Word SW:($Word) ($Acc) (.*$)#<a name=$1_A></A><a href=\#$1_H>SWISS|$1 $2</A> $3\t<b>[<a href=$DbUrl{'swpr'}$2>SwissProt</a> / <A HREF=$DbUrl{'gb_p'}$2>Entrez</A>]</b>#o;

  s#^>NR_SC:PIR-$Word PIR:($Word) (.*$)#<a name=$1_A> </A><a href=\#$1_H>PIR|$1</A> $2\t<b>[<a href=$DbUrl{'pir_uid'}$1>PIR</a> / <A HREF=$DbUrl{'gb_p'}$1>Entrez</A>]</b>#o;

  s#^>CHRS:([A-Z][0-9]*) (.*)$#<a name=$1_A></a><a href=\#$1_H>$1</A> $2:  [<b><a href=$SGDUrl{'seq_an'}$1>Gene/Sequence Resources</A> / <a href=\"$SGDUrl{'map_chr'}$1\">ORF Map</a></b>]#o;

  s#^>NOT:([A-Z]_[0-9]*-[0-9]*)(  *)Chromosome ([0-9]*) from ([0-9]*) to ([0-9]*)$#<a name=$1_A></a><a href=\#$1_H>$1</A> $2Chromosome $3 from $4 to $5 [<b><a href=$SGDUrl{'chr'}$3\&beg=$4\&end=$5>Gene/Sequence Resources</a> / <a href=\"$SGDUrl{'map_chr'}$3\&beg=$4\&end=$5\">ORF Map</a> / <a href=\"$SGDUrl{'chr_old'}$3\&beg=$4\&end=$5\">Retrieve DNA</a></b>]#o;

  s#^>UTR5_SC_[0-9]*:(\S*) 5' untranslated region, chr(\S*) ([0-9]*) - ([0-9]*)(.*$)#<a name=$1_A></A><a href=\#$1_H>UTR5:$1</A> $1 5' untranslated region, chr$2 $3 - $4, $5\t<b>[<A HREF=\"$SGDUrl{'chr'}$2&beg=$3&end=$4\">Gene/Sequence Resources</A> / <a href=\"$SGDUrl{'map_chr'}$2\&beg=$3\&end=$4\">ORF Map</a>]</b>#o;

  # Hits without a db identifier.
  # If any of the previous regexps succeed, the leading '>' will be removed.
  # Otherwise, this regexp could cause trouble.
  s@^>($Word)(.*)$@<a name=$1_A></a>$1 $2<br>(<a href="\#$1_H">Back|<a href="\#top">Top</a>)@o;

##
## REGEXPS FOR SUMMARY TABLE LINES AT TOP OF REPORT (a.k.a. 'descriptions')
## (table of sequence id, description, score, P/Expect value, n)
##
## Not using bold face to highlight the sequence id's since this can throw off
## off formatting of the line when the IDs are different lengths. This lead to 
## the scores and P/Expect values not lining up properly.

    ### NCBI-specific markups for description lines:

  # GenBank/EMBL, DDBJ hits (GenBank Format):
  s@^ ?(gb|emb|dbj|ref)\|($Word)(\|$Word)?($Descrip)($Int +)($Signif)(.*)$@$1:<a href="$_gi_link$2">$2$3</a>$4$5<A href="\#$2_A">$6</a>$7<a name="$2_H"></a>@o;

  s@^ ?(gb|emb|dbj|ref)\|($Word)(\| \(?$Word\)?)($Descrip)($Int +)($Signif)(.*)$@$1:<a href="$_gi_link$2">$2</a>$3$4$5<A href="\#$2_A">$6</a>$7<a name="$2_H"></a>@o;

    # Missing inner ID
  s@^ ?pir\|\|($Word)?($Descrip)($Int)  ($Signif)(.*)$@<a href="$DbUrl{'pir_acc'}$1">pir</a>:<a href="$DbUrl{'gb_p'}$1">$1</a> $2$3  <A href="\#$1_A">$4</a>$5<a name="$1_H"></a>@o;

  # GI hits (GenBank Format):  using a nested (())
  s@^ ?gi\|($Word)( +\(($Word)\))($Descrip)($Int)  ($Signif)(.*)$@gi:<a href="$_gi_link$3">$1</a>$2$4$5  <A href="\#$3_A">$6</a>$7<a name="$3_H"></a>@o;

  s@^ ?(gnl)\|($Word)?(\|$Word +)\(($Word)\)($Descrip)($Int)  ($Signif)(.*)$@$1:<a href="$_gi_link$4">$2$3</a>($4)$5$6  <A href="\#$4_A">$7</a>$8<a name="$4_H"></a>@o;


  s@^ ?(bbs|prf)\|\|?($Word)($Descrip)($Int)  ($Signif)(.*)$@$1:<a href="$_gi_link$2">$2</a> $3$4  <A href="\#$2_A">$5</a>$6<a name="$2_H"></a>@o;


  ## SwissProt accessions (GenBank format)
  s@^ ?sp\|($Word)(\|$Word)?($Descrip)($Int)  ($Signif)(.*)$@<a href="$DbUrl{'swpr'}$1">sp</a>:<a href="$DbUrl{'gb_p'}$1">$1$2</a>$3$4  <a href="\#$1_A">$5</a>$6<a name="$1_H"></a>@o;

  ## PDB ids with or without a chain ID (GenBank format)
  s@^ ?pdb\|($Word)\|($Word)?($Descrip)($Int)  ($Signif)(.*)$@<a href="$DbUrl{'3db'}$1">pdb</a>:<a href="$DbUrl{'gb_struct'}$1">$1_$2</a>$3$4  <a href="\#$1_A">$5</a>$6<a name="$1_H"></a>@o;


 ### SGD-specific markups for description lines:

  ## PDB ids without chain identifier
  s@^ ?PDB_UNIQUEP:(\d\w{3})_($Descrip)($Int)  ($Signif)(.*)$@<a href="$DbUrl{'3db'}$1">PDB</a>:<A HREF="$DbUrl{'gb_struct'}$1">$1</A>       $2$3  <a href="\#$1_A">$4</a>$5<a name="$1_H"></a>@o;


  ## PDB ids with chain identifier
  s@^ ?PDB_UNIQUEP:(\d\w{3})_(\w)($Descrip)($Int)  ($Signif)(.*)$@<a href="$DbUrl{'3db'}$1">PDB</a>:<A HREF="$DbUrl{'gb_struct'}$1">$1</A> Chain:$2$3$4  <a href="\#$1_A">$5</a>$6<a name="$1_H"></a>@o;


  s@^ ?($Word)PEPT:GI_(\d+)($Descrip)($Int)  ($Signif)(.*)$@$1:<A HREF="$DbUrl{'gb_p'}$2">GI_$2</A> $3 $4  <a href="\#$2_A">$5</a> $6<a name="$2_H"></a>@o;

  s@^ *WORMPEPT:(\S+)($Descrip)($Int)  ($Signif)(.*)$@WORMPEP:<A HREF="$DbUrl{'wormace'}$1">$1</a> $2 $3 <a href="\#$1_A">$4</a>$5<a name="$1_H"></a>@o;

    ## Mike Cherry's markups. SAC note: added back database name to allow
    ## the HTML-formatted version to be parsable by Blast.pm.
    
  s#^ ?(GB_$Word:)($Word)( *)($Acc)($Descrip)($Int)  ( *$Signif) ( *\d*)$#GenBank\|<a href="$_gi_link$4">$2</A>\|$4 $3$5$6 <a href="\#$2_$4_A">$7</A> $8<a name="$2_$4_H"></A>#o;

# Mike's version:
#  s#^ ?(ORFP:)(\S*)($Descrip)($Int)  ($Signif) ($Int)$#$1<b>$2</b> $3 $4 <a href="\#$2_A">$5</a> $6<a name="$2_H"></a>#o;

# My modification:
  s@^ ?ORFP:(\S*) +([\w-]+)(.*[ ]{2,3})($Int)  ($Signif) ($Int)$@ORFP:<A HREF=\"$SGDUrl{'locus'}$2\">$1 $2</A>$3$4 <a href="\#$1_A">$5</a> $6<a name="$1_H"></a>@o;

  s#^ ?(ORFN:)(\S*)($Descrip)($Int)  ($Signif) ($Int)$#$1$2 $3 $4 <a href="\#$2_A">$5</a> $6<a name="$2_H"></a>#o;

  s#^ ?(NR_SC:GP-)(\S*) ( *)gi\|(\w+)([\w\|]*)($Descrip)($Int)  ($Signif) ($Int)$#GenPept\|<a href="$DbUrl{'gb_p'}$4">$4</A>$3 gp|$2 $5$6$7  <a href="\#$4_A">$8</A> $9<a name="$4_H"></A>#o;

  s#^ ?(NR_SC:SW-)$Word ( *)SW:($Word) ($Acc)($Descrip)($Int)  ($Signif) ($Int)$#SWISS\|<a href="$DbUrl{'swpr'}$4">$3</A>   SW:$3 $4 $5$6  <a href="\#$3_A">$7</A> $8<a name="$3_H"></A>#o;

  s#^ ?(NR_SC:PIR-)$Word ( *)PIR:($Word)($Descrip)($Int)  ($Signif) ($Int)$#PIR\|<a href="$DbUrl{'pir_uid'}$3">$3</A>   $2   PIR:$3 $4$5  <a href="\#$3_A">$6</A> $7<a name="$3_H"></A>#o;

  s#^ ?(CHRS:)([A-Z][0-9]*)($Descrip)($Int)  ($Signif) ($Int)$#$1Segment:$2 $3  $4  <a href="\#$2_A">$5</a> $6<a name="$2_H"></a>#o;

  s#^ ?(CHR[0-9]*)($Descrip)($Int)  ($Signif) ($Int)$#$1 $2  $3  <a href="\#$1_A">$4</a> $5<a name="$1_H"></a>#o;

  s#^ ?(NOT:)([A-Z]_[0-9]*-[0-9]*)($Descrip)($Int)  ($Signif) ($Int)$#$1$2 $3  $4  <a href="\#$2_A">$5</a> $6<a name="$2_H"></a>#o;

  s#^ ?(UTR5_SC_[0-9]*:)(\S*)($Descrip)($Int)  ($Signif) ($Int)$#UTR5:$2 $3  $4 <a href="\#$2_A">$5</a> $6<a name="$2_H"></a>#o;

  # Hits without a db identifier. 
  s@^ ?($Word)($Descrip)($Int)  ($Signif)(.*)$@$1$2$3  <A href="\#$1_A">$4</a>$5<a name="$1_H"></a>@o;

    $$line_ref = $_;
}




=head2 _prog_ref_html

 Usage     : n/a; utility method used by get_html_func().
 Purpose   : Get a special alert for BLAST reports against all of GenBank/EMBL.
 Returns   : string with HTML

See Also   : L<get_html_func>()

=cut

#------------------
sub _prog_ref_html {
#------------------
    return <<"QQ_REF_QQ";
<p>
<small>
<b>References:</b> 
<ol>
<li>Altschul, Stephen F., Warren Gish, Webb Miller, Eugene W. Myers, and David J. Lipman (1990). 
Basic local alignment search tool.
<a href="http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?uid=2231712&form=6&db=m&Dopt=r">J. Mol. Biol. 215: 403-10</a>.
<li>Altschul et al. (1997), Gapped BLAST and PSI-BLAST: 
a new generation of protein database search programs. 
<a href="http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?uid=9254694&form=6&db=m&Dopt=r">Nucl. Acids Res. 25: 3389-3402</a>.
<li><b>Program Descriptions</b>: 
<a href="http://www.ncbi.nlm.nih.gov/BLAST/newblast.html">BLAST2</a> |
<a href="http://blast.wustl.edu/">WU-BLAST2</a> |
<a href="http://www.ncbi.nlm.nih.gov/BLAST/blast_help.html">Help Manual</a>
</ol>
<small>
HTML formatting provided by the <a href="${\$BioWWW->home_url('bioperl')}Projects/Blast/">Bioperl Blast module</a>.
</small>
</small>
<p>

QQ_REF_QQ

# Not really a reference for the Blast algorithm itself but an interesting usage.
#<li>Gish, Warren, and David J. States (1993). Identification of protein coding regions by database similarity search. 
#<a href="http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?uid=8485583&form=6&db=m&Dopt=r">Nature Genetics 3:266-72</a>.

}


=head2 _genbank_alert

 Usage     : n/a; utility method used by get_html_func().
 Purpose   : Get a special alert for BLAST reports against all of GenBank/EMBL.
 Returns   : string with HTML

See Also   : L<get_html_func>()

=cut

#------------------
sub _genbank_alert {
#------------------
    return << "QQ_GENBANK_QQ";
<p><b><font color="red">CAUTION: Hits reported on this page may be derived from DNA sequences 
         that contain more than one gene. 
         </font>To avoid mis-interpretation, always check database entries
         for any sequence of interest to verify that the similarity 
         occurs within the described sequence. (E.g., A DNA sequence
         for gene X as reported in GenBank may contain a 5' or 3' 
         fragment of coding sequence for a neighboring gene Y, yet will
         be listed as gene X, since gene Y had not yet been identified). </b>
QQ_GENBANK_QQ
}



=head2 strip_html

 Usage     : $boolean = &strip_html( string_ref );
           : This method is exported.
 Purpose   : Removes HTML formatting from a supplied string.
           : Attempts to restore the Blast report to enable
           : parsing by Bio::Tools::Blast.pm.
 Returns   : Boolean: true if string was stripped, false if not.
 Argument  : string_ref = reference to a string containing the whole Blast
           :              report.
 Throws    : Croaks if the argument is not a scalar reference.
 Comments  : Based on code originally written by Alex Dong Li
           : (ali@genet.sickkids.on.ca).
           : This method does some Blast-specific stripping 
           : (adds back a '>' character in front of each HSP 
           : alignment listing).
           :   
           : THIS METHOD IS HIGHLY ERROR-PRONE!
           :
           : Removal of the HTML tags and accurate reconstitution of the
           : non-HTML-formatted report is highly dependent on structure of
           : the HTML-formatted version. For example, it assumes that first 
           : line of each alignment section (HSP listing) starts with a
           : <a name=..> anchor tag. This permits the reconstruction of the 
           : original report in which these lines begin with a ">".
           : This is required for parsing.
           :
           : If the structure of the Blast report itself is not intended to
           : be a standard, the structure of the HTML-formatted version
           : is even less so. Therefore, the use of this method to
           : reconstitute parsable Blast reports from HTML-format versions
           : should be considered a temorary solution.

See Also   : B<Bio::Tools::Blast::parse()>

=cut

#---------------
sub strip_html {
#---------------
      # This may not best way to remove html tags. However, it is simple.
      # it won't work under following conditions:
      #    1) if quoted > appears in a tag  (does this ever happen?)
      #    2) if a tag is split over multiple lines and this method is
      #       used to process one line at a time.
      
    my $string_ref = shift;

    ref $string_ref eq 'SCALAR' or 
	croak ("Can't strip HTML: ".
	       "Argument is should be a SCALAR reference not a ${\ref $string_ref}");

    my $str = $$string_ref;
    my $stripped = 0;

    # Removing "<a name =...>" and adding the '>' character for 
    # HSP alignment listings.
    $str =~ s/(\A|\n)<a name ?=[^>]+> ?/>/sgi and $stripped = 1;

    # Removing all "<>" tags. 
    $str =~ s/<[^>]+>|&nbsp//sgi and $stripped = 1;

    # Re-uniting any lone '>' characters.
    $str =~ s/(\A|\n)>\s+/\n\n>/sgi and $stripped = 1;

    $$string_ref = $str;
    $stripped;
}

1;
__END__

#####################################################################################
#                                END OF CLASS                                       #
#####################################################################################



