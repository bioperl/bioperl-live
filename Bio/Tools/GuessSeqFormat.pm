#------------------------------------------------------------------
#
# BioPerl module Bio::Tools::GuessSeqFormat
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Andreas Kähäri, andreas.kahari@ebi.ac.uk
#
# You may distribute this module under the same terms as perl itself
#------------------------------------------------------------------

=encoding utf-8

=head1 NAME

Bio::Tools::GuessSeqFormat - Module for determining the sequence
format of the contents of a file, a string, or through a
filehandle.

=head1 SYNOPSIS

    # To guess the format of a flat file, given a filename:
    my $guesser = Bio::Tools::GuessSeqFormat->new( -file => $filename );
    my $format  = $guesser->guess;

    # To guess the format from an already open filehandle:
    my $guesser = Bio::Tools::GuessSeqFormat->new( -fh => $filehandle );
    my $format  = $guesser->guess;
    # The filehandle will be returned to its original position. Note that this
    # filehandle can be STDIN.

    # To guess the format of one or several lines of text (with
    # embedded newlines):
    my $guesser = Bio::Tools::GuessSeqFormat->new( -text => $linesoftext );
    my $format = $guesser->guess;

    # To create a Bio::Tools::GuessSeqFormat object and set the
    # filename, filehandle, or line to parse afterwards:
    my $guesser = Bio::Tools::GuessSeqFormat->new();
    $guesser->file($filename);
    $guesser->fh($filehandle);
    $guesser->text($linesoftext);

    # To guess in one go, given e.g. a filename:
    my $format = Bio::Tools::GuessSeqFormat->new( -file => $filename )->guess;

=head1 DESCRIPTION

Bio::Tools::GuessSeqFormat tries to guess the format ("swiss",
"pir", "fasta" etc.) of the sequence or MSA in a file, in a
scalar, or through a filehandle.

The guess() method of a Bio::Tools::GuessSeqFormat object will
examine the data, line by line, until it finds a line to which
only one format can be assigned.  If no conclusive guess can be
made, undef is returned.

If the Bio::Tools::GuessSeqFormat object is given a filehandle,
e.g. STDIN, it will be restored to its original position on
return from the guess() method.

=head2 Formats

Tests are currently implemented for the following formats:

=over

=item *

ACeDB ("ace")

=item *

Blast ("blast")

=item *

ClustalW ("clustalw")

=item *

Codata ("codata")

=item *

EMBL ("embl")

=item *

FastA sequence ("fasta")

=item *

FastQ sequence ("fastq")

=item *

FastXY/FastA alignment ("fastxy")

=item *

Game XML ("game")

=item *

GCG ("gcg")

=item *

GCG Blast ("gcgblast")

=item *

GCG FastA ("gcgfasta")

=item *

GDE ("gde")

=item *

Genbank ("genbank")

=item *

Genscan ("genscan")

=item *

GFF ("gff")

=item *

HMMER ("hmmer")

=item *

PAUP/NEXUS ("nexus")

=item *

Phrap assembly file ("phrap")

=item *

NBRF/PIR ("pir")

=item *

Mase ("mase")

=item *

Mega ("mega")

=item *

GCG/MSF ("msf")

=item *

Pfam ("pfam")

=item *

Phylip ("phylip")

=item *

Prodom ("prodom")

=item *

Raw ("raw")

=item *

RSF ("rsf")

=item *

Selex ("selex")

=item *

Stockholm ("stockholm")

=item *

Swissprot ("swiss")

=item *

Tab ("tab")

=item *

Variant Call Format ("vcf")

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and
other Bioperl modules.  Send your comments and suggestions
preferably to one of the Bioperl mailing lists.  Your
participation is much appreciated.

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

Report bugs to the Bioperl bug tracking system to help us
keep track the bugs and their resolution.  Bug reports can be
submitted via the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Andreas KE<228>hE<228>ri, andreas.kahari@ebi.ac.uk

=head1 CONTRIBUTORS

Heikki LehvE<228>slaiho, heikki-at-bioperl-dot-org
Mark A. Jensen, maj-at-fortinbras-dot-us

=cut


package Bio::Tools::GuessSeqFormat;

use strict;
use warnings;


use base qw(Bio::Root::Root);

=head1 METHODS

Methods available to Bio::Tools::GuessSeqFormat objects
are described below.  Methods with names beginning with an
underscore are considered to be internal.

=cut

=head2 new

 Title      : new
 Usage      : $guesser = Bio::Tools::GuessSeqFormat->new( ... );
 Function   : Creates a new object.
 Example    : See SYNOPSIS.
 Returns    : A new object.
 Arguments  : -file The filename of the file whose format is to
                    be guessed, e.g. STDIN, or
              -fh   An already opened filehandle from which a text
                    stream may be read, or
              -text A scalar containing one or several lines of
                    text with embedded newlines.

    If more than one of the above arguments are given, they
    are tested in the order -text, -file, -fh, and the first
    available argument will be used.

=cut

sub new
{
    my $class = shift;
    my @args  = @_;

    my $self = $class->SUPER::new(@args);

    my $attr;
    my $value;

    while (@args) {
        $attr = shift @args;
        $attr = lc $attr;
        $value = shift @args;
        $self->{$attr} = $value;
    }

    return $self;
}

=head2 file

 Title      : file
 Usage      : $guesser->file($filename);
              $filename = $guesser->file;
 Function   : Gets or sets the current filename associated with
              an object.
 Returns    : The new filename.
 Arguments  : The filename of the file whose format is to be
              guessed.

    A call to this method will clear the current filehandle and
    the current lines of text associated with the object.

=cut

sub file
{
    # Sets and/or returns the filename to use.
    my $self = shift;
    my $file = shift;

    if (defined $file) {
        # Set the active filename, and clear the filehandle and
        # text line, if present.
        $self->{-file} = $file;
        $self->{-fh} = $self->{-text} = undef;
    }

    return $self->{-file};
}

=head2 fh

 Title      : fh
 Usage      : $guesser->fh($filehandle);
              $filehandle = $guesser->fh;
 Function   : Gets or sets the current filehandle associated with
              an object.
 Returns    : The new filehandle.
 Arguments  : An already opened filehandle from which a text
              stream may be read.

    A call to this method will clear the current filename and
    the current lines of text associated with the object.

=cut

sub fh
{
    # Sets and/or returns the filehandle to use.
    my $self = shift;
    my $fh = shift;

    if (defined $fh) {
        # Set the active filehandle, and clear the filename and
        # text line, if present.
        $self->{-fh} = $fh;
        $self->{-file} = $self->{-text} = undef;
    }

    return $self->{-fh};
}


=head2 text

 Title      : text
 Usage      : $guesser->text($linesoftext);
              $linesofext = $guesser->text;
 Function   : Gets or sets the current text associated with an
              object.
 Returns    : The new lines of texts.
 Arguments  : A scalar containing one or several lines of text,
              including embedded newlines.

    A call to this method will clear the current filename and
    the current filehandle associated with the object.

=cut

sub text
{
    # Sets and/or returns the text lines to use.
    my $self = shift;
    my $text = shift;

    if (defined $text) {
        # Set the active text lines, and clear the filehandle
        # and filename, if present.
        $self->{-text} = $text;
        $self->{-fh} = $self->{-file} = undef;
    }

    return $self->{-text};
}

=head2 guess

 Title      : guess
 Usage      : $format = $guesser->guess;
              @format = $guesser->guess; # if given a line of text
 Function   : Guesses the format of the data accociated with the
              object.
 Returns    : A format string such as "swiss" or "pir".  If a
              format can not be found, undef is returned.
 Arguments  : None.

    If the object is associated with a filehandle, the position
    of the filehandle will be returned to its original position
    before the method returns.

=cut

our %formats = (
    ace         => { test => \&_possibly_ace        },
    blast       => { test => \&_possibly_blast      },
    bowtie      => { test => \&_possibly_bowtie     },
    clustalw    => { test => \&_possibly_clustalw   },
    codata      => { test => \&_possibly_codata     },
    embl        => { test => \&_possibly_embl       },
    fasta       => { test => \&_possibly_fasta      },
    fastq       => { test => \&_possibly_fastq      },
    fastxy      => { test => \&_possibly_fastxy     },
    game        => { test => \&_possibly_game       },
    gcg         => { test => \&_possibly_gcg        },
    gcgblast    => { test => \&_possibly_gcgblast   },
    gcgfasta    => { test => \&_possibly_gcgfasta   },
    gde         => { test => \&_possibly_gde        },
    genbank     => { test => \&_possibly_genbank    },
    genscan     => { test => \&_possibly_genscan    },
    gff         => { test => \&_possibly_gff        },
    hmmer       => { test => \&_possibly_hmmer      },
    nexus       => { test => \&_possibly_nexus      },
    mase        => { test => \&_possibly_mase       },
    mega        => { test => \&_possibly_mega       },
    msf         => { test => \&_possibly_msf        },
    pfam        => { test => \&_possibly_pfam       },
    phrap       => { test => \&_possibly_phrap      },
    phylip      => { test => \&_possibly_phylip     },
    pir         => { test => \&_possibly_pir        },
    prodom      => { test => \&_possibly_prodom     },
    raw         => { test => \&_possibly_raw        },
    rsf         => { test => \&_possibly_rsf        },
    selex       => { test => \&_possibly_selex      },
    stockholm   => { test => \&_possibly_stockholm  },
    swiss       => { test => \&_possibly_swiss      },
    tab         => { test => \&_possibly_tab        },
    vcf         => { test => \&_possibly_vcf        },
);

sub guess
{
    my $self = shift;

    while (my ($fmt_key) = each (%formats)) {
        $formats{$fmt_key}{fmt_string} = $fmt_key;
    }

    my $fh;
    my $start_pos;
    if (defined $self->{-text}) {
        # Break the text into separate lines.
        my $text = $self->{-text};
        open $fh, '<', \$text or $self->throw("Could not read from string: $!");

    } elsif (defined $self->{-file}) {
        # If given a filename, open the file.
        my $file = $self->{-file};
        open $fh, '<', $file or $self->throw("Could not read file '$file': $!");

    } elsif (defined $self->{-fh}) {
        # If given a filehandle, get the current position in the stream.
        $fh = $self->{-fh};
        if (not seek $fh, 0, 1) { # seek to current position to determine seekability
            # Work around non-seekable filehandles if IO::Scalar is available
            # (adapted from http://www.perlmonks.org/?node_id=33587)
            # IO::Mark may be an option for very large streams?
            $self->throw("Need IO::Scalar to guess from unseekable filehandles")
                if not eval { require IO::Scalar };
            my $data;
            { local $/; $data = <$fh>; $.-- };  # copy raw data from fh
            tie *$fh, 'IO::Scalar', my $s;      # replace fh by scalar-tied fh
            print $fh $data;                    # write raw data to tied fh
            seek $fh, 0, 0;                     # return to start of tied fh
        }
        $start_pos = tell $fh;
    }

    my $done  = 0;
    my $lineno = 0;
    my $guess;
    while (!$done) {
        my $line;       # The next line of the file.
        my $match = 0;  # Number of possible formats of this line.

        last if (!defined($line = <$fh>));
        next if ($line =~ /^\s*$/); # Skip white and empty lines.
        chomp $line;
        $line =~ s/\r$//;   # Fix for DOS files on Unix.
        ++$lineno;

        while (my ($fmt_key, $fmt) = each (%formats)) {
            if ($fmt->{test}($line, $lineno)) {
                ++$match;
                $guess = $fmt->{fmt_string};
            }
        }

        # We're done if there was only one match.
        $done = ($match == 1);
    }

    if (defined $self->{-fh}) {
        # Go back to original position in filehandle
        seek $fh, $start_pos, 0 or $self->throw("Could not reset filehandle $fh: $!");
    } else {
        # Close the filehandle we opened
        close $fh;
    }
    return ($done ? $guess : undef);
}

=head1 HELPER SUBROUTINES

All helper subroutines will, given a line of text and the line
number of the same line, return 1 if the line possibly is from a
file of the type that they perform a test of.

A zero return value does not mean that the line is not part
of a certain type of file, just that the test did not find any
characteristics of that type of file in the line.

=head2 _possibly_ace

From bioperl test data, and from
"http://www.isrec.isb-sib.ch/DEA/module8/B_Stevenson/Practicals/transcriptome_recon/transcriptome_recon.html".

=cut

sub _possibly_ace
{
    my ($line, $lineno) = (shift, shift);
    return ($line =~ /^(?:Sequence|Peptide|DNA|Protein) [":]/);
}

=head2 _possibly_blast

 From various blast results.

=cut

sub _possibly_blast
{
    my ($line, $lineno) = (shift, shift);
    return ($lineno == 1 &&
        $line =~ /^[[:upper:]]*BLAST[[:upper:]]*.*\[.*\]$/);
}

=head2 _possibly_bowtie

Contributed by kortsch.

=cut

sub _possibly_bowtie
{
    my ($line, $lineno) = (shift, shift);
    return ($line =~ /^[[:graph:]]+\t[-+]\t[[:graph:]]+\t\d+\t([[:alpha:]]+)\t([[:graph:]]+)\t\d+\t[[:graph:]]?/)
            && length($1)==length($2);
}

=head2 _possibly_clustalw

From "http://www.ebi.ac.uk/help/formats.html".

=cut

sub _possibly_clustalw
{
    my ($line, $lineno) = (shift, shift);
    return ($lineno == 1 && $line =~ /CLUSTAL/);
}

=head2 _possibly_codata

From "http://www.ebi.ac.uk/help/formats.html".

=cut

sub _possibly_codata
{
    my ($line, $lineno) = (shift, shift);
    return (($lineno == 1 && $line =~ /^ENTRY/) ||
            ($lineno == 2 && $line =~ /^SEQUENCE/) ||
            $line =~ m{^(?:ENTRY|SEQUENCE|///)});
}

=head2 _possibly_embl

From
"http://www.ebi.ac.uk/embl/Documentation/User_manual/usrman.html#3.3".

=cut

sub _possibly_embl
{
    my ($line, $lineno) = (shift, shift);
    return ($lineno == 1 && $line =~ /^ID   / && $line =~ /BP\.$/);
}

=head2 _possibly_fasta

From "http://www.ebi.ac.uk/help/formats.html".

=cut

sub _possibly_fasta
{
    my ($line, $lineno) = (shift, shift);
    return (($lineno != 1 && $line =~ /^[A-IK-NP-Z]+$/i) ||
            $line =~ /^>\s*\w/);
}

=head2 _possibly_fastq

From bioperl test data.

=cut

sub _possibly_fastq
{
    my ($line, $lineno) = (shift, shift);
    return ( ($lineno == 1 && $line =~ /^@/) ||
             ($lineno == 3 && $line =~ /^\+/) );
}

=head2 _possibly_fastxy

From bioperl test data.

=cut

sub _possibly_fastxy
{
    my ($line, $lineno) = (shift, shift);
    return (($lineno == 1 && $line =~ /^ FAST(?:XY|A)/) ||
            ($lineno == 2 && $line =~ /^ version \d/));
}

=head2 _possibly_game

From bioperl testdata.

=cut

sub _possibly_game
{
    my ($line, $lineno) = (shift, shift);
    return ($line =~ /^<!DOCTYPE game/);
}

=head2 _possibly_gcg

From bioperl, Bio::SeqIO::gcg.

=cut

sub _possibly_gcg
{
    my ($line, $lineno) = (shift, shift);
    return ($line =~ /Length: .*Type: .*Check: .*\.\.$/);
}

=head2 _possibly_gcgblast

From bioperl testdata.

=cut

sub _possibly_gcgblast
{
    my ($line, $lineno) = (shift, shift);
    return (($lineno == 1 && $line =~ /^!!SEQUENCE_LIST/) ||
            ($lineno == 2 &&
             $line =~ /^[[:upper:]]*BLAST[[:upper:]]*.*\[.*\]$/));
}

=head2 _possibly_gcgfasta

From bioperl testdata.

=cut

sub _possibly_gcgfasta
{
    my ($line, $lineno) = (shift, shift);
    return (($lineno == 1 && $line =~ /^!!SEQUENCE_LIST/) ||
            ($lineno == 2 && $line =~ /FASTA/));
}

=head2 _possibly_gde

From "http://www.ebi.ac.uk/help/formats.html".

=cut

sub _possibly_gde
{
    my ($line, $lineno) = (shift, shift);
    return ($line =~ /^[{}]$/ ||
            $line =~ /^(?:name|longname|sequence-ID|
                          creation-date|direction|strandedness|
                          type|offset|group-ID|creator|descrip|
                          comment|sequence)/x);
}

=head2 _possibly_genbank

From "http://www.ebi.ac.uk/help/formats.html".
Format of [apparantly optional] file header from
"http://www.umdnj.edu/rcompweb/PA/Notes/GenbankFF.htm". (TODO: dead link)

=cut

sub _possibly_genbank
{
    my ($line, $lineno) = (shift, shift);
    return (($lineno == 1 && $line =~ /GENETIC SEQUENCE DATA BANK/) ||
            ($lineno == 1 && $line =~ /^LOCUS /) ||
            ($lineno == 2 && $line =~ /^DEFINITION /) ||
            ($lineno == 3 && $line =~ /^ACCESSION /));
}

=head2 _possibly_genscan

From bioperl test data.

=cut

sub _possibly_genscan
{
    my ($line, $lineno) = (shift, shift);
    return (($lineno == 1 && $line =~ /^GENSCAN.*Date.*Time/) ||
            ($line =~ /^(?:Sequence\s+\w+|Parameter matrix|Predicted genes)/));
}

=head2 _possibly_gff

From bioperl test data.

=cut

sub _possibly_gff
{
    my ($line, $lineno) = (shift, shift);
    return (($lineno == 1 && $line =~ /^##gff-version/) ||
            ($lineno == 2 && $line =~ /^##date/));
}

=head2 _possibly_hmmer

From bioperl test data.

=cut

sub _possibly_hmmer
{
    my ($line, $lineno) = (shift, shift);
    return (($lineno == 2 && $line =~ /^HMMER/) ||
            ($lineno == 3 &&
             $line =~ /Washington University School of Medicine/));
}

=head2 _possibly_nexus

From "http://paup.csit.fsu.edu/nfiles.html".

=cut

sub _possibly_nexus
{
    my ($line, $lineno) = (shift, shift);
    return ($lineno == 1 && $line =~ /^#NEXUS/);
}

=head2 _possibly_mase

From bioperl test data.
More detail from "http://www.umdnj.edu/rcompweb/PA/Notes/GenbankFF.htm" (TODO: dead link)

=cut

sub _possibly_mase
{
    my ($line, $lineno) = (shift, shift);
    return (($lineno == 1 && $line =~ /^;;/) ||
            ($lineno > 1 && $line =~ /^;[^;]?/));
}

=head2 _possibly_mega

From the ensembl broswer (AlignView data export).

=cut

sub _possibly_mega
{
    my ($line, $lineno) = (shift, shift);
    return ($lineno == 1 && $line =~ /^#mega$/);
}


=head2 _possibly_msf

From "http://www.ebi.ac.uk/help/formats.html".

=cut

sub _possibly_msf
{
    my ($line, $lineno) = (shift, shift);
    return ($line =~ m{^//} ||
            $line =~ /MSF:.*Type:.*Check:|Name:.*Len:/);
}

=head2 _possibly_phrap

From "http://biodata.ccgb.umn.edu/docs/contigimage.html". (TODO: dead link)
From "http://genetics.gene.cwru.edu/gene508/Lec6.htm".    (TODO: dead link)
From bioperl test data ("*.ace.1" files).

=cut

sub _possibly_phrap
{
    my ($line, $lineno) = (shift, shift);
    return ($line =~ /^(?:AS\ |CO\ Contig|BQ|AF\ |BS\ |RD\ |
                          QA\ |DS\ |RT\{)/x);
}

=head2 _possibly_pir

From "http://www.ebi.ac.uk/help/formats.html".
The ".,()" spotted in bioperl test data.

=cut

sub _possibly_pir # "NBRF/PIR" (?)
{
    my ($line, $lineno) = (shift, shift);
    return (($lineno != 1 && $line =~ /^[\sA-IK-NP-Z.,()]+\*?$/i) ||
            $line =~ /^>(?:P1|F1|DL|DC|RL|RC|N3|N1);/);
}

=head2 _possibly_pfam

From bioperl test data.

=cut

sub _possibly_pfam
{
    my ($line, $lineno) = (shift, shift);
    return ($line =~ m{^\w+/\d+-\d+\s+[A-IK-NP-Z.]+}i);
}

=head2 _possibly_phylip

From "http://www.ebi.ac.uk/help/formats.html".  Initial space
allowed on first line (spotted in ensembl AlignView exported
data).

=cut

sub _possibly_phylip
{
    my ($line, $lineno) = (shift, shift);
    return (($lineno == 1 && $line =~ /^\s*\d+\s\d+/) ||
            ($lineno == 2 && $line =~ /^\w\s+[A-IK-NP-Z\s]+/) ||
            ($lineno == 3 && $line =~ /(?:^\w\s+[A-IK-NP-Z\s]+|\s+[A-IK-NP-Z\s]+)/)
           );
}

=head2 _possibly_prodom

From "http://prodom.prabi.fr/prodom/current/documentation/data.php".

=cut

sub _possibly_prodom
{
    my ($line, $lineno) = (shift, shift);
    return ($lineno == 1 && $line =~ /^ID   / && $line =~ /\d+ seq\.$/);
}

=head2 _possibly_raw

From "http://www.ebi.ac.uk/help/formats.html".

=cut

sub _possibly_raw
{
    my ($line, $lineno) = (shift, shift);
    return ($line =~ /^[A-Za-z\s]+$/);
}

=head2 _possibly_rsf

From "http://www.ebi.ac.uk/help/formats.html".

=cut

sub _possibly_rsf
{
    my ($line, $lineno) = (shift, shift);
    return (($lineno == 1 && $line =~ /^!!RICH_SEQUENCE/) ||
            $line =~ /^[{}]$/ ||
            $line =~ /^(?:name|type|longname|
                          checksum|creation-date|strand|sequence)/x);
}

=head2 _possibly_selex

From "http://www.ebc.ee/WWW/hmmer2-html/node27.html".

Assuming presence of Selex file header.  Data exported by
Bioperl on Pfam and Selex formats are identical, but Pfam file
only holds one alignment.

=cut

sub _possibly_selex
{
    my ($line, $lineno) = (shift, shift);
    return (($lineno == 1 && $line =~ /^#=ID /) ||
            ($lineno == 2 && $line =~ /^#=AC /) ||
            ($line =~ /^#=SQ /));
}

=head2 _possibly_stockholm

From bioperl test data.

=cut

sub _possibly_stockholm
{
    my ($line, $lineno) = (shift, shift);
    return (($lineno == 1 && $line =~ /^# STOCKHOLM/) ||
            $line =~ /^#=(?:GF|GS) /);
}



=head2 _possibly_swiss

From "http://ca.expasy.org/sprot/userman.html#entrystruc".

=cut

sub _possibly_swiss
{
    my ($line, $lineno) = (shift, shift);
    return ($lineno == 1 && $line =~ /^ID   / && $line =~ /AA\.$/);
}

=head2 _possibly_tab

Contributed by Heikki.

=cut

sub _possibly_tab
{
    my ($line, $lineno) = (shift, shift);
    return ($lineno == 1 && $line =~ /^[^\t]+\t[^\t]+/) ;
}

=head2 _possibly_vcf

From "http://www.1000genomes.org/wiki/analysis/vcf4.0".

Assumptions made about sanity - format and date lines are line 1 and 2
respectively. This is not specified in the format document.

=cut

sub _possibly_vcf
{
    my ($line, $lineno) = (shift, shift);
    return (($lineno == 1 && $line =~ /##fileformat=VCFv/) ||
            ($lineno == 2 && $line =~ /##fileDate=/));
}



1;
