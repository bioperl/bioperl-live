#-------------------------------------------------------------------------------
# PACKAGE : Bio::Tools::Blast::Run::LocalBlast.pm
# PURPOSE : Protoype module for running a Blast analysis on a local machine.
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu).
# CREATED : 22 May 1998
# STATUS  : STUB MODULE, MUST BE CUSTOMIZED FOR A LOCAL SITE.
# REVISION: $Id$
#
# For the latest version and documentation, visit the distribution site:
#    http://genome-www.stanford.edu/perlOOP/bioperl/blast/
#
# To generate documentation, run this module through pod2html
# (preferably from Perl v5.004 or better).
#
# This module as-is provides only framework for running Blasts locally.
# Customize it for your site. See comments below for tips.
#
#-------------------------------------------------------------------------------

package Bio::Tools::Blast::Run::LocalBlast;
use strict;

use Bio::Root::Object    qw(&_rearrange);
use Bio::Root::Global    qw(:devel);
use Bio::Root::Utilities qw(:obj);
use Carp;

use Exporter;
use vars qw( @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS );
@ISA        = qw(Exporter);
@EXPORT     = qw();
@EXPORT_OK  = qw(&blast_local @Blast_dbp_local @Blast_dbn_local 
		 @Blast_matrix_local);
%EXPORT_TAGS = ( std => [qw(&blast_local @Blast_dbp_local @Blast_dbn_local 
		 @Blast_matrix_local)] );

use strict;
use vars qw( $ID $VERSION);

$ID      = 'Bio::Tools::Blast::Run::LocalBlast';
$VERSION = '0.01';


## POD Documentation

=head1 NAME

Bio::Tools::Blast::Run::LocalBlast.pm - Bioperl module for running
Blast analyses locally.

=head1 SYNOPSIS

    use Bio::Tools::Blast::Run::LocalBlast qw(&blast_local);

    &blast_local( %named_parameters);

See L<blast_local>() for a description of available parameters.

=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

=head1 DESCRIPTION

Bio::Tools::Blast::Run::LocalBlast.pm contains methods and data
necessary for running Blast sequence analyses on a local machine. This
module must be customized for a specific site.

The basic requirements are that it conform to this minimal API:

=over 4

=item 1 Export a method called L<blast_local>()

that accepts a Bio::Tools::Blast.pm object + named parameters as
specified by L<blast_local>().

=item 2 The L<blast_local>() method should return 

a list of names of files containing the raw Blast reports.

=item 3 Export arrays containing a list of available databases 

in the arrays C<@Blast_dbn_local> and C<@Blast_dbp_local>.

=back

The generic version of this module provides some rudimentary logic,
but feel free to customize as necessary.


=head2 Script Files

Sometimes it is convenient to write an executable shell script for
running a set of Blasts on a local machine. This script can be saved
and re-executed as necessary or saved for documentation purposes. This
module could provide a convenient way to consolidate the logic
necessary for producing such script files or perhaps stubs of script
file that could be further modified for Blast-ing specific datasets.


=head1 DEPENDENCIES

Bio::Tools::Blast::Run::LocalBast.pm is used by B<Bio::Tools::Blast.pm>
The development of this is linked with the Blast.pm module and should
be updated along with that module.

=head1 SEE ALSO

 Bio::Tools::Blast.pm                    - Blast object.
 Bio::Tools::Blast::Run::postclient.pl   - Script for accessing remote server.
 Bio::Tools::Blast::Run::Webblast.pm     - Utility module for running Blasts remotely.
 Bio::Tools::Blast::HTML.pm              - Blast HTML-formating utility class.
 Bio::Seq.pm                             - Biosequence object  


 http://bio.perl.org/Projects/modules.html  - Online module documentation
 http://bio.perl.org/Projects/Blast/        - Bioperl Blast Project     
 http://bio.perl.org/                       - Bioperl Project Homepage

=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

    bioperl-l@bioperl.org          - General discussion
    bioperl-guts-l@bioperl.org     - Technically-oriented discussion
    http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via email
or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR

Steve A. Chervitz, sac@genome.stanford.edu

=head1 VERSION

Bio::Tools::Blast::Run::LocalBlast.pm, 0.01

=head1 COPYRIGHT

Copyright (c) 1998 Steve A. Chervitz. All Rights Reserved.
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


use vars qw(@Blast_dbp_local @Blast_dbn_local @Blast_matrix_local);

# Peptide datasets
@Blast_dbp_local = qw( YOUR LOCAL PEPTIDE BLAST DATASETS HERE );

# Nucleotide datasets
@Blast_dbn_local = qw( YOUR LOCAL NUCLEOTIDE BLAST DATASETS HERE );

# See the blast_local() method below for a list of common matrices.
@Blast_matrix_local  = qw( YOUR LOCAL SUBSTITUTION SCORING MATRICES HERE );

# PRIVATE VARIABLES 
#-------------------
my @seq_queue = ();   # list of sequence objects.
my @file_queue = ();  # list of sequence file names.

my $program                 ='';  # your default blast program 
my $database                ='';  # your default database name
my $gappedAlignmentFlag     ='';  # flag to indicate if gapping is on or off
my $extension               = $program;  # output file name extension base.
my $expect                  = 1;  # default expect value     
my $filter                  ='';   # default filter (e.g., 'seq+xnu').
my $matrix                  ='';   # default scoring matrix
my $descriptions            = 100;  # default V
my $alignments              = 100;  # default B
my $maximumLength           = 100000; # maximum length of query sequence.
my $minimumLength           = 10;     # minimum length of query sequence.

my @_outFiles               = ();
my $_out_dir                = '';    # where to save the Blast output files.
my $_blastObj               = undef;
my $_errmsg                 = '';


######################  BEGIN FUNCTIONS  ########################

=head1 APPENDIX

Methods beginning with a leading underscore are considered private
and are intended for internal use by this module. They are
B<not> considered part of the public interface and are described here
for documentation purposes only.


=head2 blast_local

 Usage     : @files = blast_local($blast_object,  %namedParameters);
           : This method is exported.
 Purpose   : Run a local Blast analysis on one or more sequences.
           : This method defines the API for your LocalBlast.pm module.
 Returns   : Array containing a list of filenames of the Blast reports.
 Argument  : $blast_object = object ref for a Bio::Tools::Blast.pm object.
           : %named parameters: (PARAMETER TAGS CAN BU UPPER OR LOWER CASE)
           : These are some basic parameters. Supply more as desired.
           :
           :   -SEQS       => ref to an array of Bio::Seq.pm objects. 
           :   -SEQ_FILES  => ref to an array of strings containing full-path file names.
           :   -PROG       => name of blast program (blastp, blastx, etc.)
           :   -DATABASE   => name of database (see below.)
           :   -EXPECT     => expect value cutoff
           :   -FILTER     => sequence complexity filter ('default' or 'none')
           :   -MATRIX     => substitution scoring matrix (blast1 only for NCBI server)
           :   -DESCR      => integer, number of on-line descriptions (V, 100)
           :   -ALIGN      => integer, number of alignments (B, 100)
           :   -GAP        => 'on' or 'off'
           :   -OUT_DIR    => output directory to store blast result files
           :   
 Throws    : Exception if:
           :   - Cannot obtain parameters by calling _rearrange() on the
           :     first argument, which should be a Bio::Tools::Blast.pm object ref.
           :   - No sequences are provided (objects or files).
           :   - Sequence type is incompatible with Blast program type.
           :   - Database name is not one of the valid names.
 Comments  :
  -------------------------------------------------------------
  Available programs: blastn, blastx, blastp, tblastn, tblastx

  -------------------------------------------------------------
  Available local databases are: 

   LIST YOUR LOCAL DATABASES HERE. 
    These are exported by this module in the @Blast_dbp_local
    and @Blast_dbn_local arrays.
  -------------------------------------------------------------
  Available substitution scoring matrices: 
    (Here are the standard ones)
    BLOSUM: 100,90,85,80,75,70,65,62,60,55,50,45,40,35,30
    PAM:    500,490,480,470,460,450,440,430,420,410,400,390,380,370,360,350
            340,330,320,310,300,290,280,270,260,250,240,230,220,210,200,190,
            180,170,160,150,140,130,120,110,100,90,80,70,60,50,40,30,20,10
    OTHER: DAYHOFF, GONNET, IDENTITY, MATCH

    These are exported by this module in the @Blast_matrix_local
  -------------------------------------------------------------
  Available sequence complexity filters:
       SEG, SEG+XNU, XNU, dust, none.

See Also : _set_options(), _validate_options(), _blast_seqs(), _blast_files(), B<Bio::Tools::Blast.pm>

=cut

#----------------
sub blast_local {
#----------------
    my ($bobj, %param) = @_;

    my ($seq_a, $file_a);
    eval { 
	# _rearrange() is an instance method of Bio::Root::Object.pm and is
	# inherited by the Blast object.
	($seq_a, $file_a) = $bobj->_rearrange([qw(SEQS SEQ_FILES)], %param);
    };
    if($@) {
	# Ideally we should use can(), requires Perl 5.004.
	croak "Can't run local BLAST: failed to obtain parameters using _rearrange().\n".
	    "Make sure that the first argument is a Bio::Tools::Blast.pm object.\n";
    }

    @seq_queue  = @$seq_a if ref $seq_a eq 'ARRAY';
    @file_queue = @$file_a if ref $file_a eq 'ARRAY';

    if( not (scalar @seq_queue or scalar @file_queue)) {
	$bobj->throw("Can't Blast: No sequences supplied.");
    }

    @_outFiles = ();

    $_blastObj = $bobj;
    &_set_options(%param);
    &_validate_options();
    @seq_queue ? &_blast_seqs() : &_blast_files();

    return @_outFiles;
}

#----------------
sub _set_options {
#----------------
# Not yet supporting all possible parameters, just a few key ones.
    my (%param) = @_;
    my ($prog, $db, $exp, $filt, $mat, $descr, $aln, $gap, 
	$maxl, $minl, $dir)
	= &Bio::Root::Object::_rearrange($_blastObj, 
	  	      [qw(PROG DATABASE EXPECT FILTER MATRIX 
			  DESCR ALIGN GAP MAX_LEN MIN_LEN OUT_DIR)], %param);
		      
    if($prog)    { $program = $prog;  }
    if(!$gap)    { $gappedAlignmentFlag = '-nogap';  }  #default is on
    if($db)      { $database = $db;  }
    if($exp)     { $expect = $exp;  }
    if($filt)    { $filter = $filt;   }
    if($mat)     { $matrix = $mat;   }
    if($aln)     { $alignments = $aln;   }
    if($descr)   { $descriptions = $descr;   }
    if($minl)    { $minimumLength = $minl; }
    if($maxl)    { $maximumLength = $maxl; } 
    if($dir)     { $_out_dir = $dir; }

}

#----------------------
sub _validate_options {
#----------------------

    # if program is not of one available program, quit
    if (not $program =~ m/^t?blast[npx]$/i) {
	$_blastObj->throw("Invalid program name: $program",
			  "It must be one of blastn, blastx, blastp, tblastn, tblastx!");
    }

    my (@db);

    # if database is not of one available databases, quit
    # only need to check nucleotide databases since they are a superset.
    if ( not(grep $database eq $_, @Blast_dbn_local)) {
	my $valid = join(', ', @Blast_dbn_local);
	$_blastObj->throw("Invalid Blast database name: $database",
			  "It must be one of $valid");
    } else {
	$database = $db[0];  # permits case-insensitive db specification
    }
}


#---------------
sub _blast_seqs {
#---------------
# POSSIBLE STRATEGY:
# foreach seq object in @seq_queue:
    # convert sequence to required format (fasta, gcg, etc).
    # create a temp file for sequence.
    # create a path for saving the blast output
    # create argument string ($prog $seqfile $args > $outfile)
    # get environment variables for executing blast (see _get_environment()).
    # write temp script file to execute the blast.
    # execute the temp script file.
    # check for output file to confirm success.

    $_blastObj->throw("Can't Blast sequence objects: Stub method not implemented.");

    my ($seq, $typeDo, $advancedOptions);
    my ($id, $outputFileNameTemp, $outputFileName, $baseFileName, $options);
    my $countSuccessfulBlasting = 0;

    my $count = 0;
    foreach $seq (@seq_queue) {

	$count++;
	# Not catching fatal error if the sequence is the wrong type
	# (serious error).
	next unless &_validate_seq($seq);
	
	# Get sequence in desired format.
#	my $sequenceInFastaFormat= $seq->layout('fasta');
#	my $sequenceInGCGFormat= $seq->layout('gcg');

    }
}


#---------------
sub _blast_files {
#---------------
# Same general strategy as for _blast_seqs() but starting with the existing files.

    my ($file, $typeDo, $advancedOptions);
    my ($id, $outputFileNameTemp, $outputFileName, $baseFileName, $options);
    my $countSuccessfulBlasting = 0;

    $_blastObj->throw("Can't Blast sequence files: Stub method not implemented.");

    my $count = 0;
    foreach $file (@file_queue) {
	$count++;
    }
}


#------------------
sub _validate_seq {
#------------------
    my $seq = shift;
    my $type = $seq->type;
    my $len = $seq->seq_len;

    # Verify that sequence type is correct type for selected program.
    if(($type =~ /[dr]na/i and ($program =~ m/^(blastp|tblastn)$/i))
	  or
       ($type =~ /amino/i and ($program =~ m/^(blast[nx]|tblastx|dbest)$/i))
       ) {
	$_blastObj->throw("The sequence of ${\$seq->id} is in wrong type for $program search");
    }

    # Verify that sequence type is correct type for selected database.
    if(($type =~ /[dr]na/i and not (grep $database eq $_, @Blast_dbn_local))
       or
       ($type =~ /amino/i and not (grep $database eq $_, @Blast_dbp_local))
       ) {
	$_blastObj->throw("The sequence of ${\$seq->id} is wrong type for database $database");
    }

    if ($len < $minimumLength) {
	$MONITOR && print STDERR ("\nThe sequence length of ${\$seq->id} ($len) too short for blasting: min length = $minimumLength\a\n\n");
	return 0;
    } elsif( $len > $maximumLength) {
	$MONITOR && print STDERR ("\nThe sequence length of ${\$seq->id} ($len) is too long for blasting: max length = $maximumLength\a\n\n");
	return 0;
    }
    1;
}

#-------------------
sub _get_environment {
#-------------------
# Returns a string for setting environment variables for local Blast analysis.
<<"QQ_ENV_QQ";

QQ_ENV_QQ
}
    

1;
__END__

#####################################################################################
#                                END OF CLASS                                       #
#####################################################################################


