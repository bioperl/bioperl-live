# Seq.pm 
#
#  $Id$
#
# MODIFICATION NOTES: See bottom of file.

# Copyright (c) 1996 Georg Fuellen, Richard Resnick, Steven E. Brenner,
# Chris Dagdigian, Steve Chervitz, Ewan Birney and others. All Rights Reserved.
# This module is free software; you can redistribute it and/or modify 
# it under the same terms as Perl itself.

package Bio::Seq;

require 5.003;
use Carp;
require Exporter;
use Bio::Root::Object ();  


# Uncomment the next lines if you want to use Autoloading.
# (See Autoloading notes below for more info).
#use AutoLoader;      ## SAC: Added autoloading
#*AUTOLOAD = \&AutoLoader::AUTOLOAD;

@ISA         = qw(Bio::Root::Object Exporter); ## SAC: added Bio::Root::Object to @ISA 
@EXPORT      = qw();
@EXPORT_OK   = qw($VERSION %SeqAlph @SeqAlph %SeqForm @SeqForm);
#%EXPORT_TAGS = (std => [%SeqAlph @SeqAlph %SeqForm @SeqForm]);
$VERSION     = 0.051;
#$NOFILEKEYWORD= "_nofile";


use vars qw(%SeqAlph @SeqAlph %SeqForm @SeqForm $READSEQ_EXISTS %TypeSeq);
# SeqAlph and SeqForm may be used externally

## Try to cleanly determine if Parse.pm is installed
## and configured correctly...

$READSEQ_EXISTS = undef;  # shut strict up

eval { 
    use Bio::Parse;
    my $ok = $Bio::Parse::OK;  # prevents warnings
    if(defined($Bio::Parse::OK)) { $READSEQ_EXISTS = "Y";  }
}; 


## POD-formatted documentation

=head1 NAME

Bio::Seq - bioperl sequence object

=head1 SYNOPSIS

=head2 Object Creation

 $seq = Bio::Seq->new;
 
 $seq = Bio::Seq->new($filename);
 
 $seq = Bio::Seq->new(-seq=>'ACTGTGGCGTCAACTG');
 
 $seq = Bio::Seq->new(-seq=>$sequence_string);
 
 $seq = Bio::Seq->new(-seq=>@character_list);
 
 $seq = Bio::Seq->new(-file=>'seqfile.aa',
		      -desc=>'Sample Bio::Seq sequence',
		      -start=>'1',
		      -type=>'Amino',
		      -ffmt=>'Fasta');
 
 $seq = Bio::Seq->new($file,$seq,$id,$desc,$names,
                     $numbering,$type,$ffmt,$descffmt);


=head2 Object Manipulation

 $seq->[METHOD];

 $result = $seq->[METHOD];
 
 
 
 Accessors
 --------------------------------------------------------
 There are a wide variety of methods designed to give easy
 and flexible access to the contents of sequence objects
 
 The following accessors can be invoked upon a sequence object

 ary()        - access sequence (or slice of sequence) as an array
 str()        - access sequence (or slice of sequence) as a string
 getseq()     - access sequence (or slice) as string or array
 seq_len()    - access sequence length
 id()         - access/change object id 
 desc()       - access/change object description
 names()      - access/change object names
 start()      - access/change start point of the sequence (see note below) 
 end()        - access/change end point of the sequence (see note below)
 numbering()  - access/change sequence numbering offset (deprecated)
 origin()     - access/change sequence origin
 type()       - access/change sequence type
 ffmt()       - access/change default output format
 descffmt()   - access/change description format
 setseq()     - change sequence
 

 Methods
 --------------------------------------------------------
 The following methods can be invoked upon a sequence object

 copy()        - returns an exact copy of an object
 alphabet_ok() - check sequence against genetic alphabet  
 alphabet()    - returns the genetic alphabet currently in use
 layout()      - sequence formatter for output
 revcom()      - reverse complement of sequence
 complement()  - complement of sequence  
 reverse()     - reverse of sequence
 Dna_to_Rna()  - translate Dna seq to Rna
 Rna_to_Dna()  - translate Rna seq to Dna
 translate()   - protein translation of Dna/Rna sequence


=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

=head1 DESCRIPTION

CVS version: $Id$

This module is the generic sequence object which lies at the core of
the bioperl project. It stores Dna, Rna, or Protein sequence
information and annotation. It has associated methods to perform
various manipulations of sequences and support for a reading and
writing sequence data in a variety of file formats.

Bio::Seq has completly superceeded Bio::PreSeq.pm.

The older PreSeq.pm code can be found at Chris Dagdigian's site:
http://www.sonsorol.org/dag/bioperl/top.html


=over 2

=item * BASED ON PreSeq.pm, THIS VERSION OF Seq.pm HAS BEEN INTEGRATED INTO THE BIOPERL FRAMEWORK.

For a complete description of these changes, see the comments
at the top of the source.

=back

=head2 Sequence Types

Currently the following sequence types are recognized:

 Dna
 Rna
 Amino

=head2 Alphabets

This module uses the standard extended single-letter genetic
alphabets to represent nucleotide and amino acid sequences.

In addition to the standard alphabet, the following symbols
are also acceptable in a biosequence:

 ?  (a missing nucleotide or amino acid)
 -  (gap in sequence)

=head2 Extended Dna / Rna alphabet

 (includes symbols for nucleotide ambiguity)
 ------------------------------------------
 Symbol       Meaning      Nucleic Acid
 ------------------------------------------
  A            A           Adenine
  C            C           Cytosine
  G            G           Guanine
  T            T           Thymine
  U            U           Uracil
  M          A or C  
  R          A or G   
  W          A or T    
  S          C or G     
  Y          C or T     
  K          G or T     
  V        A or C or G  
  H        A or C or T  
  D        A or G or T  
  B        C or G or T   
  X      G or A or T or C 
  N      G or A or T or C 

 
 IUPAC-IUB SYMBOLS FOR NUCLEOTIDE NOMENCLATURE:
   Cornish-Bowden (1985) Nucl. Acids Res. 13: 3021-3030.


=head2  Amino Acid alphabet

 ------------------------------------------
 Symbol           Meaning   
 ------------------------------------------
 A        Alanine
 B        Aspartic Acid, Asparagine
 C        Cystine
 D        Aspartic Acid
 E        Glutamic Acid
 F        Phenylalanine
 G        Glycine
 H        Histidine
 I        Isoleucine
 K        Lysine
 L        Leucine
 M        Methionine
 N        Asparagine
 P        Proline
 Q        Glutamine
 R        Arginine
 S        Serine
 T        Threonine
 V        Valine
 W        Tryptophan
 X        Unknown
 Y        Tyrosine
 Z        Glutamic Acid, Glutamine
 *        Terminator

 
 IUPAC-IUP AMINO ACID SYMBOLS:
   Biochem J. 1984 Apr 15; 219(2): 345-373
   Eur J Biochem. 1993 Apr 1; 213(1): 2

=head2 Output Formats

The following output formats are currently supported:
Raw, Fasta, GCG, GenBank, PIR

=head2 Input Formats

In addition to "raw" sequence files, Seq.pm is currently 
only able to read in Fasta and GCG formatted single sequence
files. Support for additional formats is forthcoming.

Seq.pm has the ability to make use of D.G. Gilbert's ReadSeq
program when reading in sequence files. ReadSeq has the ability
to read and interconvert between many different biological
sequence formats.

When readseq is present and Seq.pm has been properly
configured to use it, ReadSeq will be invoked when internal
parsing code fails to recognize the sequence.

Formats which readseq currently understands:

  - IG/Stanford
  - GenBank/GB
  - NBRF
  - EMBL
  - GCG
  - DnaStrider
  - Fitch format
  - Pearson/Fasta
  - Zuker format
  - Olsen format
  - Phylip3.2
  - Phylip
  - Plain/Raw
  * MSF
  * PAUP's multiple sequence (NEXUS) format
  * PIR/CODATA format used by PIR
  * ASN.1 format used by NCBI

  Note: Formats indicated with a '*' allow for multiple
        sequences to be contained within one file. At this
        time, the behaviour of Seq.pm with regard to these
        multiple-sequence files has not been specified.

Readseq is freely distributed and is available in
shell archive (.shar) form via FTP from
ftp.bio.indiana.edu (129.79.224.25) in the
molbio/readseq directory.
(URL) ftp://ftp.bio.indiana.edu/molbio/readseq/

If ReadSeq is not available or Seq.pm is not configured
to use it, internal parsing mechanisms will be used.

Currently supported filetypes for input:
Raw, Fasta, GCG

=head1 USAGE

=head2 Installation 

Seq.pm requires the use of other bioperl modules, particularly
the Bio::Root framework. This module should be installed along
with the rest of the bioperl code.


=head2 Why modules and object-oriented code?

Perl5 is nice in that it allows users to use OO-style programming only in
the situations where they feel like doing so.

=over 4

=item * Simple interfaces to complex tasks.

From the perspective of novice or occasional perl users, objects are useful
because they can offer direct and simple ways to do things that in reality
may be somewhat complex or arcane. Users interact with and manipulate
objects via specific, documented methods and never have to worry about what
is going on "behind the scenes." Many  perl programmers have devoted
significant amounts of time and effort creating easy-to-use "wrappers"
around complex or abstract tasks. Visit the CPAN Module list at 
(URL) http://www.perl.com/perl/CPAN/CPAN.html to see the fruits of their labor.
 
=item * Reusability.

From the prospective of a perl power-user, object-oriented programming
allows programmers to write code that is easily scalable and reusable. This
allows powerful applications to be built rapidly with and with a minimum of
waste or repeated effort.

=back
 
=head2 Using Bio::Seq in your perl programs

Seq.pm is invoked via the perl 'use' command

   use Seq;

=head2 Creating a biosequence object

The "constructor" method in Seq.pm is the L<new>() function.

The proper syntax for accessing the L<new>() function in Seq.pm is as follows:

   $myseq = Bio::Seq->new;

Of course, objects are only useful if they have something in them so you
would probably want to pass along some additional information or arguments
to the constructor. The foundation of any biosequence object is course the
sequence itself.

You can address L<new>() with a sequence directly:

   $myseq = Bio::Seq->new(-seq=>'AACTGGCGTTCGTG');

Or you can pass in a string or a list:

   $myseq = Bio::Seq->new(-seq=>$sequence_string);
   $myseq = Bio::Seq->new(-seq=>@sequence_list);

It is also possible to create a new sequence object based on a sequence
contained in a file. You can tell constructor where to find the sequence
file by passing in the 'file' parameter:

   $myseq  = Bio::Seq->new(-file=>'seqfile.gcg');

Because there are so many different conventions or formats for storing
sequence information in files, it would be polite (although not absolutely
necessary) to tell the constructor what format the sequence file is in. We
can provide that information via the file-format or 'ffmt' field. To create
a sequence object based upon a GCG-formatted sequence file:

   $myseq  = Bio::Seq->new(-file=>'seqfile.gcg',-ffmt=>'GCG');

We've already introduced three different object attributes or arguments
that can be passed to the L<new>() object constructor ('seq','file' and
'ffmt') so now would be a good time to introduce them all:

B<BioSeq Constructor Arguments>

B<file:>
 The "file" argument should be a string value containing path and filename information for a sequence file that is to be read into an object.

B<seq:>
The "seq" argument is for passing in sequence directly instead of reading
in a sequence file. The sequence should consist of RAW info (no whitespace,
newlines or formatting) and can be passed in as either an array/list or
string.

B<id:>
The "id" argument should be a ONE-WORD string value giving a short name for
the sequence.

B<desc:>
The "desc" argument should be a string containing a description of the
sequence. This field is not limited to one word.

B<names:>
The "names" argument should be a hash or reference to a hash that contains
any number of user generated key-value pairs. Various bits of identifying
information can be stored here including name(s), database locations,
accession numbers, URL's, etc.

B<type:>
The "type" argument should be a string value describing the sequence type
eg; "Dna", "Rna" or "Amino".

B<origin:>
The "origin" argument should be a string value describing sequence origin info

B<start:>
The start point, in biological coordinates of the sequence

B<end:>
The end point, in biological coordinates of the last residue in 
the sequence

start/end attributes are not strongly tied to what is actually in 
the sequence (ie, $seq->start()+length($seq->getseq()) doesn't 
necessarily equal $seq->end()-1 - most of the time it should).

This is to allow some oddities to be stored in the Seq object sensibly.
 
The numbering convention is 'biological' coordinates. ie the sequence
ATG would start at 1 (A) and finish at 3 (G). (NB - this is different
from how perl represents ranges in sequences).

numbering() is equivalent to start() (old version). Eventually it
will be removed. numbering() accesses the same attribute as start()

B<numbering:>
(Deprecated) The "numbering" argument should be an integer value containing the sequence
numbering offset value. By default all sequence are numbered starting with
1. 


B<ffmt:>
The "ffmt" argument should be a string describing sequence file-format. If
a sequence is being read from a file via the "file" argument, "ffmt" is
used to invoke the proper parsing code. "ffmt" is also the default format
for sequence output when the layout method is called. See elsewhere in this
documentation for info regarding recognized sequence file-formats.

If most of these arguments were used at once to create a sequence object,
it would look something like this:

   #Set up the name hash
   %names = (
   'CloneID','DB1',
   'Isolate','5',
   'Tissue','Xenopus',
   'Location','/usr2/users/dag/bioperl/sample.tfa'
   );

   $name_ref = \%names;

   #Create the object
   $myseq = new Bio::Seq(-file=>'sample.tfa',
                         -names=>$name_ref,
                         -type=>'Dna',
                         -origin=>'Xenopus mesoderm',
                         -start=>'1',
                         -desc=>'Sample Bio::Seq sequence',
                         -ffmt=>'Fasta');

=head2 Methods

Once an object has been created, there are defined ways to go about
accessing the information -- users are encouraged to poke around "under the
hood" of Seq.pm to see what is going on but it is considered bad form to
bypass the defined accession methods and mess around with the internal
code. Bypassing the defined methods "voids the warrantee" of the module and
can lead to problems down the road. The implied agreement between module
creators and users is that the creators will strive to keep the interface
standard and backwards-compatible while the users will avoid becoming
dependent on bits of internal code that may change or disappear in future
revisions.
 
Detailed information about each method described here can be found in the
Appendix.

=head2 Accessing information 

For each defined way to access information from a biosequence object, there
is a corresponding "method" that is invoked. What follows is a brief
description of each accessor method. For more detailed information see the
individual annotations for each method near the end of this document. 

=over 4

=item * Sequence

The sequence can be accessed in several ways via the L<getseq>() method.
Depending on how it is invoked, it can return either a string or a list
value.

Both examples are appropriate:

   @sequence_list   = $myseq->getseq;
   $sequence_string = $myseq->getseq;

Sequence "slices" can be accessed by passing start and stop integer
position arguments to L<getseq>():

   @slice = $myseq->getseq($start,$stop);
   @slice = $myseq->getseq(1,50);
   @slice = $myseq->getseq(100);

If no stop value is passed in, L<getseq>() will return a slice from the start
position to the end of the sequence. Slices are returned in the context of
the object "start" attribute, not absolute position so be aware of the
objects numbering scheme.

Sequences can also be accessed in with the L<ary>() and L<str>() methods. The
L<ary>() method will always return a list value and L<str>() will always return a
string. Otherwise they are functionally identical to the L<getseq>() method.

   $sequence = $myseq->str;
   @sequence = $myseq->ary;
 
   @slice = $myseq->ary($start,$stop);
   $slice = $myseq->str($start,$stop);


=item * Sequence length

The sequence length can be accessed using the L<seq_len>() method

   $len = $myseq->seq_len;

=item * Sequence ID

The ID field can be accessed using the L<id>() method

   $ID = $myseq->id;

=item * Description

The object description field can be accessed using the L<desc>() method

   $description = $myseq->desc;

=item * Names

The associative array (hash) that contains flexible information regarding
alternative sequence names, database locations, accession numbers, etc. can
be accessed by

   %name_hash = $myseq->names;

=item * Sequence start

The biological position of the first residue in the sequence sequence can be accessed via L<start>()

   $start = $myseq->start;

=item * Sequence end

The biological position of the last residue in the sequence sequence can be accessed via L<end>()

   $end = $myseq->end;

=item * Sequence Origin

The object origin (source organism) field can be accessed via L<origin>()

  $seq_origin = $myseq->origin;

=item * File input format / default output format

The object format field can be accessed using the L<ffmt>() method

   $format = $myseq->ffmt;

=back

=head2 Changing Information in Sequence Objects

In the previous section it was shown how object attributes and values could
be retrieved from a sequence object by calling upon various methods. Many
of the above methods will also allow the user to CHANGE object attributes
by passing in additional arguments. Detailed information on each method can
be found in the L<Appendix>.

=over 4

=item * Changing the sequence

The sequence information for an object can be changed by passing a string
or list value to the L<setseq>() method. Here are some ways that sequence
information can be changed

   $myseq->seqseq($new_sequence_string);
   $myseq->setseq(@new_sequence_list);
   $myseq->setseq("aaccttgcctgc");

The L<setseq>() method checks sequence elements and warns if it finds
non-standard characters. Because of this, arbitrary sequence compositions
are not supported at this time. This method is considered slightly
'insecure' because the 'id','desc' and 'type' fields are not updated
along with the sequence. If necessary, the user must make the appropriate
changes to these fields whenever sequence information is updated or changed.

=item * Changing the sequence ID

The ID field can be changed by passing in a new ID argument to L<id>()

   $myseq->id($new_id);

=item * Changing the object description

The object description field can be changed by passing in a new argument to L<desc>()

   $myseq->desc($new_desc);

=item * Changing the object names hash

The associative array (hash) that contains flexible information regarding
alternative sequence names, database locations, accession numbers, etc. can
be changed by passing in a reference to a new hash to L<names>()

   $hash_ref = \%name_hash;
   $myseq->names($hash_ref);

=item * Changing the sequence start or end

The default numbering offset for the sequence can be changed by passing in
a new value to L<start>() or L<end>()

   $myseq->start(1);
   $myseq->start($new_value);

=item * Sequence Origin

The object origin field can be changed by passing in a new string value to L<origin>()

  $myseq->origin("mitochondrial");
  $myseq->origin($origin_string);

=item * File input format / default output format

The object format field can be accessed by passing in a new value to L<ffmt>()

   $myseq->ffmt("GCG"); 

=back

=head2 Manipulating sequences

Creating, accessing and changing biosequence objects and fields is all well
and good, but eventually you are going to want to actually do some work.

Included with Seq.pm are some commonly used utility methods for
manipulating sequence data. So far Seq.pm contains methods for:

=over 4

=item * Copying a biosequence object

using L<copy>()

    $new_obj = $myseq->copy;

=item * Reversing a sequence 

using L<reverse>()

    $reversed_seq = $myseq->reverse;

=item * Complementing a sequence

The 2nd strand, or "complement" of a biosequence can be obtained by calling
upon the L<complement>() method.

    $comp_seq = $myseq->complement;

=item * Reverse complementing a sequence

using L<revcom>()

    $rev_comp = $myseq->revcom;

=item * Translating Dna to Rna

using L<Dna_to_Rna>()

    $rna_seq = $myseq->Dna_to_Rna;

=item * Translating Rna to Dna

using L<Rna_to_Dna>()

    $dna_seq = $myseq->Rna_to_Dna;

=item * Translating Dna or Rna to protein

using L<translate>()

    $peptide_seq = $myseq->translate;

=item * Checking the sequence alphabet

To check if any nonstandard characters are present in a biosequence, an
L<alphabet_ok>() method is provided. The method returns "1" if everything is
OK, otherwise it returns a "0".

   if($myseq->alphabet_ok) { print "OK!!\n"; }
    else { print "Not OK! \n"; }


To get alphabet itself, use the L<alphabet>() method, which will return
a string containing all characters in the current alphabet.

    $alph = $myseq->alphabet;

To use restrictive alphabets that do not permit ambiguity codes,
include '-strict => 1' in the parameters sent to L<new>().
Or, for any existing sequence object, try:

    $myseq->strict(1); 
    $myseq->alphabet_ok() or die "alphabet not okay.\n";

=back

=head2 Sequence Output

There are several methods for outputting formatted sequences. For your
convenience, a "meta-output" method called L<layout>() also exists.

If L<layout>() is called without any arguments, it calls upon the output
methods as defined by the "ffmt" field.

   print $myseq->layout;

The "ffmt" field is mainly used to describe the format of a sequence
being read in from a file. It is also used as the default format for
all sequence output. If these differ (ie; the format that the 
sequence was read in is not desired as a default output style) then
"ffmt" should be set manually via the L<ffmt>() accessor method. Of course,
after reading the sequence in you are free to change "ffmt" at will.

L<layout>() can also be called with specific formats:

   $gcg_formatted_seq = $myseq->layout("GCG"):
   $fasta_seq = $myseq->layout("Fasta"):

B<Calling output methods directly>

Many output methods accept unique named parameters/arguments that allow a
greater degree of control over output format and style, to take advantage
of these abilities, the formatting methods must be called directly. See the
appendix notes describing each output format for detailed information.

  print $myseq->out_GCG(-date->"10 May 1996",
                        -caps-"up");

Most output methods will return either a string or list value depending
on how they are invoked, check the detailed method  documentation in 
the Appendix to be sure. 

   @formatted_seqlist = $myseq->out_genbank(-id=>'New ID',
                                            -def=>'User defined definition',
                                            -acc=>'User defined accession');
 
   $formatted_seqstring = $myseq->out_genbank(-id=>'New ID',
                                              -def=>'User defined definition',
                                              -acc=>'User defined accession');

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

=head1 ACKNOWLEDGEMENTS

Some pieces of the code were contributed by Steven E. Brenner, 
Steve Chervitz, Ewan Birney, Tim Dudgeon, David Curiel, and other Bioperlers.
Thanks !!!!

=head1 SEE ALSO

  UnivAln.pm - The biosequence alignment object
  Parse.pm   - The perl interface to ReadSeq

=head1 REFERENCES

BioPerl Project Page
http://bio.perl.org/

=head1 VERSION

Bio::Seq.pm, beta 0.051

=head1 COPYRIGHT

 Copyright (c) 1996-1998 Chris Dagdigian, Georg Fuellen, Richard Resnick.
 All Rights Reserved. This module is free software; you can redistribute 
 it and/or modify it under the same terms as Perl itself.

=cut

=head1 Appendix

The following documentation describes the various functions
contained in this module. Some functions are for internal 
use and are not meant to be called by the user; they are 
preceded by an underscore ("_").


=cut

#
##
###
#### END of main POD documentation. Let the code begin..
###
##
#


# List of recognized sequence types
$SeqAlph{'unknown'}  = 0;
$SeqAlph{'dna'}      = 1;
$SeqAlph{'rna'}      = 2;
$SeqAlph{'amino'}    = 3;
$SeqAlph{'otherseq'} = 4; 
$SeqAlph{'aligned'}  = 5; 

# Invert the SeqAlph hash into a %TypeSeq hash
grep {$TypeSeq{$SeqAlph{$_}} = $_} keys %SeqAlph;

# List of recognized file formats 
$SeqForm{'unknown'}  = 0;
$SeqForm{'ig'}       = 1;
$SeqForm{'genbank'}  = 2;
$SeqForm{'nbrf'}     = 3;
$SeqForm{'embl'}     = 4;
$SeqForm{'gcg'}      = 5;
$SeqForm{'strider'}  = 6;
$SeqForm{'fasta'}    = 7;
$SeqForm{'zuker'}    = 8;
$SeqForm{'msf'}      = 9;
$SeqForm{'pir'}      = 10;
$SeqForm{'raw'}      = 11;
$SeqForm{'gcg_seq'}  = 12;   # SAC: added format.
$SeqForm{'gcg_ref'}  = 13;   # SAC: added format.

use vars qw(%FuncParse %FuncOut %Alphabets);
# %FuncParse and %FuncOut are for internal use ONLY
# %FuncParse is an array of (<ffmt>,<parse_meth>), where <ffmt>
# is the file format code (e.g., 0 for "unknown", etc.), and
# <parse_meth> is the method which parses strings in that file
# format. %FuncParse{$i} is mostly set to \&parse_bad right now,
# since we don't have many file formats supported yet.
%FuncParse =
  ($SeqForm{'unknown'} => \&parse_unknown,
   $SeqForm{'fasta'}   => \&parse_fasta,
   $SeqForm{'raw'}     => \&parse_raw,
   $SeqForm{'gcg'}     => \&parse_gcg,
   $SeqForm{'pir'}     => \&parse_bad,
   );
grep {$FuncParse{$_} ||= \&parse_bad} values %SeqForm;

# array of implemented outputting routines, built in analogy to "%FuncParse"
%FuncOut =
  ($SeqForm{'unknown'} => \&out_fasta,   #Fasta is the default format !
   $SeqForm{'fasta'}   => \&out_fasta,
   $SeqForm{'raw'}     => \&out_raw,
   $SeqForm{'gcg'}     => \&out_GCG,
   $SeqForm{'genbank'} => \&out_genbank,
   $SeqForm{'pir'}     => \&out_pir,
   $SeqForm{'nbrf'}    => \&out_nbrf,
   $SeqForm{'gcg_seq'} => \&out_gcgseq,   # SAC: added format.
   $SeqForm{'gcg_ref'} => \&out_gcgref,   # SAC: added format.
   $SeqForm{'ig'}      => \&out_ig, 
   $SeqForm{'strider'} => \&out_strider,
   $SeqForm{'msf'}     => \&out_msf,
   $SeqForm{'zuker'}   => \&out_zuker,
   );

grep {$FuncOut{$_} ||= \&out_bad} values %SeqForm;

my %Alphabets =
    ($SeqAlph{'unknown'}  => [ "A","C","G","T","U","R","Y","M","K","S","W","H","B",
			     "V","D","N","A","R","N","D","C","Q","E","G","H","I",
			     "L","K","M","F","P","S","T","W","X","Y","V","Z","*"], #sac: added Z
     $SeqAlph{'dna'}      => [ "A","C","G","T","R","Y","M","K","S","W","H","B","V","D",
			     "N" ],
     $SeqAlph{'rna'}      => [ "A","C","G","U","R","Y","M","K","S","W","H","B","V","D",
			     "N" ],
     $SeqAlph{'amino'}    => [ "A","R","N","D","C","Q","E","G","H","I","L","K","M","F",
			     "P","S","T","W","X","Y","V","B","Z","*" ], # sac: added B, Z
     $SeqAlph{'aligned'}    => [ "A","R","N","D","C","Q","E","G","H","I","L","K","M","F",
			     "P","S","T","W","X","Y","V","B","Z","*","-","." ], # eb - added for alignments
     $SeqAlph{'otherseq'} => [ ],
     );

# SAC: new strict alphabet: doesn't allow any ambiguity characters.
my %Alphabets_strict =
    ($SeqAlph{'unknown'}  => $Alphabets{'unknown'},
     $SeqAlph{'dna'}      => [ "A","C","G","T" ],
     $SeqAlph{'rna'}      => [ "A","C","G","U" ],
     $SeqAlph{'amino'}    => [ "A","R","N","D","C","Q","E","G","H","I","L","K","M","F",
			     "P","S","T","W","Y","V"], 
     $SeqAlph{'otherseq'} => $Alphabets{'otherseq'},
     );
#create alphabets like ``1Mg''=$SeqAlph{dna}Mg = [ "A","C","G","T","?" ],
#and ``1Gp''=$SeqAlph{dna}Gp = [ "A","C","G","T","-" ],
#where "?" denotes the character for missing nucleotide data, "-" denotes gaps.
#Note that $Alphabets{"1GpMg"} = "A C G T - ?" is defined by this procedure !
grep {$Alphabets{$_."Gp"} ||= [ @{ $Alphabets{$_} },"-" ] } keys %Alphabets;
grep {$Alphabets{$_."Mg"} ||= [ @{ $Alphabets{$_} },"?" ] } keys %Alphabets;
{
local($^W) = 0;
grep {$Alphabets_strict{$_."Gp"} ||= [ @{ $Alphabets_strict{$_} },"-" ] } keys %Alphabets_strict;
grep {$Alphabets_strict{$_."Mg"} ||= [ @{ $Alphabets_strict{$_} },"?" ] } keys %Alphabets_strict;
}

#-----------------------------------------------------------------------

=head2 new

 Title     : new
 Usage     : $mySeq = Bio::Seq->new($file,$seq,$id,$desc,$names,
                         $start,$end,$type,$ffmt,$descffmt);
           :                - or -
           : $mySeq = Bio::Seq->new(-file=>$file,
                                   -seq=>$seq,
                                   -id=>$id,
                                   -desc=>$desc,
                                   -names=>$names,
                                   -start=>$start,
                                   -end=>$end,
                                   -type=>$type,
                                   -origin=>$origin,
                                   -ffmt=>$ffmt,
                                   -descffmt=>$descffmt);
 Function  : The constructor for this class, returns a new object.
 Example   : See usage
 Returns   : Bio::Seq object
 Argument  : $file: file from which the sequence data can be read; all
               the other arguments will overwrite the data read in.
               "_nofile" is recommanded if no file is given.
             $seq: String or array of characters
             $id: String describing the ID the user wishes to assign.
             $desc: String giving a description of the sequence
             $names: A reference to a hash which stores {loc,name}
                     pairs of other database locations and corresponding names
                     where the sequence is located.
             $start: The offset of the sequence, as an integer
             $end: The end point of the sequence, as an integer
             $type: The type of the sequence, see type()
             $origin: The sequence origin
             $ffmt: Sequence format, see ffmt()
             $descffmt: format of $desc, see descffmt()
    

=cut

#-----------------------------------------------------------------------

### SAC: new() is inherited from Bio::Root::Object.

# sub new {
#   my($this) = shift;
#   my($class,$self);
#   # See the ``Perl Module List''
#   $class = ref($this) || $this;
#   $self = {};
#   bless $self, $class;
#   $self->_initialize(@_);
#   return $self;
# }

#############################

=head2 ## Internal methods ##

=cut

#-----------------------------------------------------------------------

=head2 _initialize

 Title     : _initialize
 Usage     : n/a (internal function)
 Function  : Assigns initial parameters to a blessed object.
 Example   : 
 Returns   : 
 Argument  : As Bio::Seq->new, allows for named or listed parameters.
             See ->new for the legal types of these values.

=cut

sub _initialize {
  my($self,@p) = @_;

  # Retaining 'numbering' for backward compatibility (should switch to 'start').
  my($file,$seq,$id,$desc,$names,$numbering,$start,$end,$type,$origin,$ffmt,$descffmt) =
      $self->_rearrange([qw(FILE
			    SEQ
			    ID
			    DESC
			    NAMES
			    NUMBERING
			    START
			    END
			    TYPE
			    ORIGIN
			    FFMT
			    DESCFFMT)],
			@p);

#  printf "Seq: seq length=%d\n%s\nfile %s",length $seq, $seq, $file;<STDIN>; # SAC: tester line.

  my $make = $self->SUPER::_initialize(@p);  ## SAC: Added line.
  
  # Set default values

  $self->{"seq"} = { };
  $self->{"id"} = "No_Id_Given";
  $self->{"desc"} = "No Description Given";
  $self->{"names"} = {"none","none"};
  $self->{"start"} = 1;
  $self->{"type"} = [$SeqAlph{"unknown"},"UnknownOrigin"];
  $self->{"ffmt"} = $SeqForm{"unknown"};
  $self->{"descffmt"} = $SeqForm{"unknown"};

  # Overwrite with values from passed in filepath
  if (defined($file)) {
    $self->_file_read($file,$ffmt);
  }

## Test
#print "IN sub initialize, @type = $type[0]\n";

  # Overwrite with values from @_

  $self->id($id);  # setting id first since error reports use it.

  # SAC: A common error is to forget the type. This can be bad.
  # This helps justify separate subclasses such as AASeq.pm and NASeq.pm.
  # Doing so would also allow the partitioning of NA-specific
  # methods such as complement, DNA_to_RNA, etc. 
  # However, there are some advantages to having one class for both types
  #   -- convenient to use just one class
  #   -- type can be specified at run-time; this makes for simpler code
  #   -- flatter inheritance hierarchy
  # Yet there are also cases when one does not care about the type of sequence
  # (e.g., format conversions).

  if($self->strict) {
      $type or $self->throw("Undefined sequence type for seq $self->{'id'}.", 
			    "Valid types: ".join(',',sort(keys %SeqAlph)) 
			    );
  }
  
  # The intent here is to *set* data. However, if the arguments are undefined,
  # these calls become *gets*. This is known only at runtime and so is a bit
  # dangerous. This is the reason for the separate check for $type above.
  $self->type($type);  
  if(defined($seq)) { $self->_seq($seq); } #in case seq came from a file
  $self->desc($desc);
  $self->names($names);
  $self->start($start || $numbering);  # backwards compatibility (numbering).
  $self->end($end);   
  $self->origin($origin);
  $self->ffmt($ffmt);
  $self->descffmt($descffmt);

  return $make;  # SAC: returning make parameter instead of 1
}

## SAC: _rearrange() is now inherited from Bio::Root::Object.

=head2 _seq

 Title     : _seq()
 Usage     : n/a, internal function
 Function  : called by new() to set sequence field. Checks
           : alphabet before setting.
           :
 Returns   : n/a
 Argument  : sequence string

=cut

sub _seq {
  my($self) = shift;
  my($nseq) = @_;

  my $oseq = $self->{"seq"};
  if (defined($nseq)) {
    if (ref($nseq) eq "ARRAY") {
        $self->{"seq"} = join('', @{$nseq});
    }
    elsif (ref($nseq) eq '') { #???***
        $self->{"seq"} = $nseq;
    }
    else {
        $self->throw("No sequence was assigned since input sequence is not a string, nor an array, but a ${\ref($nseq)} or somesuch.");
    }
     ##Warn if non-standard sequence characters present 
     unless($self->alphabet_ok) {
         $self->warn("Sequence $self->{'id'} contains non-standard alphabet character(s)",
		     "Seq = $self->{'seq'}",
		     "alphabet (${\$self->type}, ${\$self->strict}): ".$self->alphabet,
		     );  # SAC: displaying the problem sequence and alphabet.
     }
  }

}

#_______________________________________________________________________

=head2 _monomer

 Title     : _monomer()
 Usage     : n/a, internal function
 Function  : Returns the internal monomer that represents
           : sequence type.
           :
           : Sequence type is treated internally as a monomer
           : defined by the %SeqAlph hash. The type field
           : is a list of format [monomer,origin]. For any
           : output outside the module, the monomer is resolved
           : back into string form via the %TypeSeq hash.
           :
 Returns   : original type setting [as monomer]
 Argument  : none

=cut

sub _monomer {
  my($self) = shift;
  my $otype = $self->{type}[0];
  return $otype;
}


#_file_read is PRELIMINARY, SeqFileHandle shall be used later
#??? If anyone could write a PRELIMINARY _fileWrite
#along these lines, would be _wonderful_ :-)

=head2 _file_read

 Title     : _file_read()
 Usage     : n/a (Internal Function)
 Function  : _file_read is called whenever the constructor is called 
           : with the name of a sequence to be read from disk.
           :
 Example   : n/a, only called upon by _initialize()
 Returns   : 
 Argument  : 

=cut

sub _file_read {
  my($self, $filename, $ffmt) = @_;
  my($ent);

        ##Read in file and invoke parsing code
        open(Seq::INPUT, $filename);
        $ent = join("\n",<Seq::INPUT>);
        close(Seq::INPUT);

        parse($self, $ent, $ffmt, $filename);

  return 1;
}

######################

=head2 ## ACCESSORS ## 

=cut

#_______________________________________________________________________

=head2 seq_len

 Title       : seq_len()
 Usage       : $len = $myseq->seq_len;
 Function    : Returns a value representing the sequence
             : length
             :
 Example     : see above
 Arguments   : none
 Returns     : integer

=cut
 
sub seq_len {
    my($self) = shift;

    my($seq) = $self->str;
    return length($seq);
}

#_______________________________________________________________________

=head2 ary

 Title     : ary
 Usage     : ary([$start,[$end]])
 Function  : Returns the sequence of the object as an array, or a substring
             of the sequence if $start/$end are defined. If $start is
             defined and $end isn't, the substring is from $start to the
             end of the sequence.
 Example   : @slice = $myObject->ary(3,9);
 Returns   : array of characters
 Argument  : $start,$end (both integers). They are interpreted w.r.t. the
             specific numeration of the sequence!! ($self->{start})

=cut

#------------------------------------'

sub ary {
  my($self,$start,$end) = @_;
  my($firstIndex,$startIndex);
  my($string);

  $firstIndex = $self->start;
  $startIndex = $firstIndex;

  $string = $self->str($start,$end);
  return split('',$string);
}

#_______________________________________________________________________

=head2 str

 Title     : str
 Usage     : str([$start,[$end]])
 Function  : Returns the sequence of the object as a string, or a slice
             of the sequence if $start/$end are defined. If $start is
             defined and $end isn't, the slice is from $start to the
             end of the sequence.
 Example   : $slice = $myObject->str(3,9);
 Returns   : string scalar
 Argument  : $start,$end (both integers). They are interpreted w.r.t. the
             specific numeration of the sequence!! ($self->{start})

=cut

#----------------------------------------------------------------------------'

sub str {
  my($self,$start,$end) = @_;
  my($firstIndex,$startindx);

  $firstIndex = $self->start;
  $startindx = $firstIndex; 
  my $true_end = length($self->{"seq"}) - 1;

  # Make sure $start,$end are in range, and set them to default if they're
  # not supplied
  if (defined($start)) {
      if ($start < $startindx) {
          $self->warn("Requested start location $start of a string starting at $startindx,");
	  $start = $startindx;  # SAC: default to the starindx
# Bug fix suggested by Tim Dudgeon
#      } elsif($start > $true_end) {
      } elsif($start > ($startindx + $true_end)) {
	  ## SAC: a more serious condition.
	  $self->throw("Requested start location $start is beyond end of string ($true_end).");
      }
  }
  else {
      $start = $startindx;
  }
  if (defined($end)) {
      if (($end - $startindx) > $true_end ) {
         $self->warn("Requested end location $end is beyond the end of the string ($true_end).");
	 $end = $true_end;  # SAC: default to the max length of seq.
#	 printf "getting substr() %d - %d\n", $start-$startindx, $end-$start+1;
#	 printf "substr()\n   ", substr($self->{"seq"}, $start-$startindx, $end-$start+1);
#	 <STDIN>;
     } elsif( $end < $startindx) {
	  ## SAC: a more serious condition.
          $self->throw("Requested end location $end is beyond start of string ($startindx).");
      }	 
  }
  else {
      $end = length($self->{"seq"}) + $startindx - 1;
  }

  return substr($self->{"seq"}, $start-$startindx, $end-$start+1);
}

#_______________________________________________________________________

=head2 seq

 Title     : seq
 Usage     : seq([$start,[$end]])
 Function  : Returns the sequence of the object as an array or a char
             string, depending on the value of wantarray. Will rtn a slice
             of the sequence if $start/$end are defined. If $start is
             defined and $end isn't, the slice is from $start to the
             end of the sequence.
 Example   : @slice = $myObject->seq(3,9);
 Returns   : regular array of characters, or a scalar string
 Argument  : $start,$end (both integers). They are interpreted w.r.t. the
             specific numeration of the sequence!! ($self->{start})
 Comments  : 

=cut

sub seq {
  my($self,$start,$end) = @_;
  return wantarray ? $self->ary($start,$end) : $self->str($start,$end);
}


#_______________________________________________________________________'

=head2 getseq

 Title     : getseq
 Usage     : getseq([$start,[$end]])
 Function  : Returns the sequence of the object as an array or a char
             string, depending on the value of wantarray. Will rtn a slice
             of the sequence if $start/$end are defined. If $start is
             defined and $end isn't, the slice is from $start to the
             end of the sequence.
 Example   : @slice = $myObject->seq(3,9);
 Returns   : regular array of characters, or a scalar string
 Throws    : Warning about deprecated method.
 Argument  : $start,$end (both integers). They are interpreted w.r.t. the
             specific numeration of the sequence!! ($self->{start})

=cut

sub getseq {
  my($self,$start,$end) = @_;
  $self->warn("Deprecated method getseq() called.", "Use seq() instead.");
  return wantarray ? $self->ary($start,$end) : $self->str($start,$end);
}



#_______________________________________________________________________'

=head2 id

 Title     : id()
 Usage     : $seq_id = $myseq->id; 
           : $myseq->id($id_string);
           :
 Function  : Sets field if an ID argument string is
           : passed in. If no arguments, returns ID value for
           : object.
           :
 Returns   : original ID value
 Argument  : sequence string

=cut

sub id {
  my($self) = shift;
  my($nid) = @_;

  my $oid = $self->{"id"};
  if (defined($nid)) {
        $self->warn("identifier $nid has illegal whitespace") if $nid =~ /\s/;
        $self->{"id"} = $nid;
        $self->{"id"} = undef if $nid eq "-undef";
    }

  return $oid;
}

#_______________________________________________________________________

=head2 desc

 Title     : desc()
 Usage     : $description = $myseq->desc; 
           : $myseq->desc($desc_string);
           :
 Function  : Sets field if an argument string is
           : passed in. If no arguments, returns original value for
           : object description field.
           :
 Returns   : original value for description
 Argument  : sequence string

=cut

sub desc {
  my($self) = shift;
  my($ndesc) = @_;

  my $odesc = $self->{"desc"};
  $self->{"desc"} = $ndesc if defined($ndesc);

  return $odesc;
}

#_______________________________________________________________________

=head2 names

 Title     : names()
 Usage     : %names = $myseq->names; 
           : $myseq->names($hash_ref);
           :
 Function  : Sets field if a name hash refrence is
           : passed in. If no arguments, returns original 
           : names hash.
           :
 Returns   : hash refrence (associative array)
 Argument  : refrence to a hash (associative array)

=cut

# EB comment - do we really want this sort of thing in this object?
# I think it is better in the heavy-weight object
 
sub names {
  my($self) = shift;
  my($nnames) = @_;

  my $onames = $self->{"names"};
  # Store the hash of $names, which is a
  # reference to key/value pairs. This is 'human-readable'
  # data; each key is a location (whether it be URL, database,
  # database query, etc.) and each value is the id at that location.
  if (defined $nnames) {

      ##Since this method SETS the values, we need to
      ## delete the old names hash contents...
      my($name_ref) = $self->names;
      foreach(keys %$name_ref) { delete $name_ref->{$_}; } 

      ## Set the new values
      foreach (keys %$nnames) {
            $self->{"names"}->{$_} = $nnames->{$_};
      }
  }

  return $onames;
}

#_______________________________________________________________________

=head2 numbering

 Title     : numbering()
 Usage     : $num_start = $myseq->start; 
           : $myseq->start($value);
           :
 Function  : Sets field if an argument is
           : passed in. If no arguments, returns original value.
           :
           : (Deprecated - should switch to start())
 Returns   : original value 
 Argument  : new value

=cut

sub numbering {
    my($self) = shift;
    my($nnums) = @_;

    #chain up to start
    return $self->start(@_);
}

=head2 start

 Title     : start
 Usage     : $start = $myseq->start(); #get
           : $myseq->start($value); #set
 Function  : the set/get for the start position
 Example   :
 Returns   : start value 
 Arguments : new value

=cut

sub start{
  my ($self,$val) = @_;

  if( defined $val ) {
      $self->{'start'} = $val;
  }
  return $self->{'start'};
}

=head2 end

 Title     : end
 Usage     : $end = $myseq->end(); #get
           : $myseq->end($value); #set
 Function  : The set/get for the end position
 Example   :
 Returns   : end value 
 Arguments : new value

=cut

sub end{
  my ($self,$val) = @_;

  if( defined $val ) {
      $self->{'end'} = $val;
  }
  return $self->{'end'};
}


=head2 get_nse

 Title    : get_nse
 Usage    : $tag = $myseq->get_nse() #
 Function : gets a string like "name/start-end". This is likely
          : to be unique in an alignment/database
          : Used alot by SimpleAlign
 Example  :
 Returns  : A string
 Arguments: Two optional arguments - first being the name/ separator, second the
            start-end separator

=cut

sub get_nse{
  my ($self,@args) = @_;
  my $sep1 = shift @args;
  my $sep2 = shift @args;

  if( ! (defined $sep1 )) {
      $sep1 = "/";
  }
  if( ! (defined $sep2)) {
      $sep2 = "-";
  }
  return join('',$self->id(),$sep1,$self->start,$sep2,$self->end);

}


#_______________________________________________________________________

=head2 origin

 Title     : origin()
 Usage     : myseq->origin($value) 
 Function  : Sets the origin field which is actually the second
           : field of the Type list. The {type} field is a 2 value list
           : with a format of ["Monomer","Origin"]
           :
 Returns   : Original value
 Argument  : string
 Comments  : SAC: Consider renaming this method to "organism()" or "species()". 
           : "origin" is ambiguous and can be easily confused with 
           : a coordinate data (0,0).

=cut

sub origin {
  my($self) = shift;
  my($ntype) = @_;

  my $otype = $self->{"type"}[1];
  $self->{"type"}[1] = $ntype if $ntype;

  return $otype;
}


#_______________________________________________________________________

=head2 type

 Title     : type()
 Usage     : myseq->type($value) 
 Function  : Sets the type field which is the first
           : field of the Type list. The {type} field is a 2 value list
           : with a format of ["Monomer","Origin"]
           :
 Returns   : String containing one of the recognized sequence types:
           : 'unknown', 'dna', 'rna', 'amino', 'otherseq', 'aligned'
           : See the %Seq::SeqAlph hash for the current types.
 Argument  : string containing a valid sequence type
           : SAC: case of user-supplied argument does not matter

=cut

sub type {
  my($self) = shift;
  my($ntype) = @_;

  my $otype = $TypeSeq{"$self->{type}[0]"};

  if(defined($ntype)) {  
      $ntype = lc($ntype); # SAC: lower-casing user-supplied string
      if($SeqAlph{$ntype}) {
	  $self->{"type"}[0] = $SeqAlph{$ntype}; 
      } else {
	  # SAC: added this extra test for null types.
	  if(scalar($ntype)) {
	      $self->throw("$ntype is not a supported sequence type for seq $self->{'id'}.", 
			   "Valid types: ".join(',',sort(keys %SeqAlph)));
# SAC: Commented this out. Why store non-information? 
# 	  } else {
#	      $self->{"type"}[0] = $SeqAlph{"unknown"};
	  }
      } 
    }

  return $otype || $TypeSeq{$SeqAlph{"unknown"}};
}


#_______________________________________________________________________

=head2 ffmt

 Title     : ffmt()
 Usage     : $format = $myseq->ffmt;
           : $myseq->ffmt("Fasta");
           : 
 Function  : The file format field is used by the internal
           : sequence parsing code when trying to read 
           : in a sequence file. It is also what is used
           : as a default output format if the layout
           : method is called without an argument.
           :
           : If a sequence object is created without
           : reading in a file, or if the file is read
           : in with the use of the ReadSeq package then
           : the ffmt field can be set to indicate any default
           : output-format preference.
           :
           : If a sequence is read from a file and parsed
           : by internal code (ReadSeq not used) then the ffmt
           : field should describe the format of the sequence
           : file. The ffmt field is used to send the sequence
           : to the correct internal parsing code.
           :
 Returns   : original ffmt value
 Argument  : recognized ffmt string value (see list of recognized 
           : formats) # SAC: What are they?! This list should be obvious.
           : Valid strings: 
           :    RAW, FASTA, GCG, IG, GENBANK, NBRF, EMBL, 
           :    MSF, PIR, GCG_SEQ, GCG_REF, STRIDER, ZUKER,
           : SAC: case of user-supplied argument does not matter

=cut

sub ffmt {
  my($self) = shift;
  my($nffmt) = @_;

  my $offmt = $self->{"ffmt"};
  if (defined($nffmt)) {
      $nffmt = lc($nffmt);  # SAC: lower-casing user-supplied string
      if (defined($SeqForm{$nffmt})) {
	  $self->{"ffmt"} = $SeqForm{$nffmt};
      }
  }
  
  return $offmt;
}


#_______________________________________________________________________

=head2 descffmt

 Title     : descffmt()
 Usage     : $desc = $myseq->descffmt;
           : $myseq->descffmt($new_value); 
 Function  : 
           :
 Returns   : original value
 Argument  : $new_value (one of the formats as defined in $SeqForm).
           : SAC: case of $new_value argument does not matter.

=cut

sub descffmt {
  my($self) = shift;
  my($ndescffmt) = @_;

  my $odescffmt = $self->{"descffmt"};
  if (defined($ndescffmt)) {
      $ndescffmt = lc($ndescffmt);  # SAC: lower-casing user-supplied string
      if (defined($SeqForm{$ndescffmt})) {
	  $self->{"descffmt"} = $SeqForm{$ndescffmt};
      }
  }
  
  return $odescffmt;
}

#_______________________________________________________________________

=head2 setseq

 Title     : setseq()
 Usage     : $self->setseq($new_sequence);
 Function  : Changes the sequence inside a bioseq object
           :
 Returns   : sequence string 
 Argument  : sequence string

=cut

sub setseq {
  my($self) = shift;
  my($nseq) = @_;

  my $oseq = $self->{"seq"};
  if (defined($nseq)) {
    if (ref($nseq) eq "ARRAY") {
        $self->{"seq"} = join('', @{$nseq});
    }
    elsif (ref($nseq) eq '') { #???***
        $self->{"seq"} = $nseq;
    }
    else {
        $self->throw("No sequence was assigned since input sequence is not a string, nor an array, but a
",ref($nseq), "or somesuch.");
    }
     ## Warn if non-standard sequence characters present 
     unless($self->alphabet_ok) {
         $self->warn("Sequence $self->{'id'} contains non-standard alphabet character(s)",
		     "Seq = $self->{'seq'}",
		     "alphabet (${\$self->type}, ${\$self->strict}): ".$self->alphabet,
		     );  # SAC: reporting the problem sequence & alphabet
     }
  }
  $self->{"seq"};
}


##############################################################
# Functions having to do with file formats, parsing,
# formatting, and so on.
##############################################################

#_______________________________________________________________________

=head2 parse

 Title     : parse
 Usage     : parse($ent,[$ffmt]);
 Function  : Invokes the proper parsing code depending on
           : the value of the object 'ffmt' field.
 Example   : $self->parse;
 Returns   : n/a
 Argument  : the prospective sequence to be parsed, 
           : and optionally its format so that it doesn't need to
           : be estimated
           : SAC: case of $ffmt argument does not matter.

=cut

#------------------------------------------------------------------------'

sub parse {
  my($self, $ent, $ffmt, $filename) = @_;

  $ffmt = lc($ffmt) if $ffmt;  # SAC: lower-casing user-supplied string

  if (defined($ffmt) && defined($SeqForm{$ffmt})) { 
    $ffmt=$SeqForm{$ffmt};
    } else {
    $ffmt=$self->{"ffmt"};         
  }

  # We simply need to call the appropriate parsing function, based
  # on the value of $ffmt. Since we've set up the %FuncParse
  # associative array to contain all of the necessary function calls
  # based on the "ffmt" value, this is rather straightforward.

  if (defined($ffmt) && (exists $FuncParse{$ffmt})) {
       return &{$FuncParse{$ffmt}}($self,$ent,$filename);
  }
  else {
        &{$FuncParse{$SeqForm{"unknown"}}}($self,$ent,$filename);
  }

  return 1;
}


#_______________________________________________________________________

=head2 parse_raw

 Title     : parse_raw
 Usage     : parse_raw;
 Function  : parses $ent into the $self->{"seq"} field, using Raw
           : file format.
 Example   : $self->parse_raw;
 Returns   : n/a
 Argument  : n/a

=cut

sub parse_raw {
  my($self,$ent) = @_;

  $self->{"seq"} = join("",split("\n",$ent));
  $self->{"ffmt"} = $SeqForm{"raw"};

  return 1;
}


#_______________________________________________________________________

=head2 parse_fasta

 Title     : parse_fasta
 Usage     : parse_fasta;
 Function  : parses $ent into the "seq" field, using Fasta
           : file format.
           :
 To-do     : use benchmark module to find best/fastest parse
           : method
           :
 Example   : $self->parse_fasta;
 Returns   : n/a
 Argument  : n/a

=cut

sub parse_fasta {
  my ($self) = shift;
  my ($ent) = @_;
  my (@lines, $head);

  @lines = split("\n", $ent);
  $head = shift @lines;
  ($self->{"id"}, $self->{"desc"}) = $head =~ /^>[ \t]*(\S*)[ \t]*(.*)$/;
  $self->{"seq"} = join("",@lines);
  $self->{"descffmt"} = $self->{"ffmt"} = $SeqForm{"fasta"};

  return 1;
}


#_______________________________________________________________________

=head2 parse_gcg

 Title    : parse_gcg
 Usage    : used by internal code
 Function : Parses the sequence out of a gcg-format string and
          : sets the object sequence field accordingly. This is
          : a simple, ineffecient method for grabbing JUST the
          : sequence.
          :
 To-do    : - parse out more info than just sequence 
          : - implement alphabet checking
          : - better regular expressions/efficiency
          : - carp on unexpected / wrong-format situations
          :
 Version  : .01 / 16 Jan 1997 
 Returns  : 1
 Argument : gcg-formatted sequence string

=cut

sub parse_gcg {
  my ($self) = shift;
  my($ent) = @_;

  ## Delete newlines
  $ent =~ s/\n+//g;

  ## Grab everything after ".."
  my($seq) = $ent =~ /\.\.(.*)$/;

  ## Delete numbers and whitespace
  $seq =~ s/\d+//g;
  $seq =~ s/\s+//g;

  $self->{"seq"} = $seq;
  return 1;
}



=head2 ## METHODS FOR FILE FORMAT AND OUTPUT  ##

#_______________________________________________________________________

=head2 layout

  Title    : layout()
 Usage     : layout([$format]);
 Function  : Returns the sequence in whichever format the user specifies,
             or in the "ffmt" field if the user does not specify a format.
 Example   : $fastaFormattedSeq = $myObj->layout("Fasta");
 Returns   : varies
 Argument  : $format (one of the formats as defined in $SeqForm).
           : SAC: case of $ffmt argument does not matter.

=cut

sub layout {
  my($self,$ffmt,$start,$end) = @_;   # SAC: added slicing ability
  
  # We figure out the user-requested format for output, or assign
  # $self->{ffmt} to it.

  $ffmt = lc($ffmt);  # SAC: lower-casing user-supplied string

  if (defined($ffmt) && defined($SeqForm{$ffmt})) { 
      $ffmt=$SeqForm{$ffmt};
  } else {
      $ffmt=$self->{"ffmt"};         
  }
  
  return &{$FuncOut{$ffmt}}($self,$start,$end);
}


#_______________________________________________________________________

=head2 out_raw

 Title     : out_raw
 Usage     : out_raw;
 Function  : Returns the sequence in Raw format.
 Example   : $self->out_raw;
 Returns   : string sequence, in raw format
 Argument  : n/a

=cut

sub out_raw {
  my($self) = @_;

  # The raw format is just the string without any whitespace
  return $self->{"seq"};
}


#_______________________________________________________________________

=head2 out_fasta

 Title     : out_fasta
 Usage     : out_fasta;
 Function  : Returns the sequence as a string in FASTA format.
 Example   : $self->out_fasta;
           :
 To-do     : benchmark code / find fastest method
           :
 Returns   : string sequence in Fasta format
 Argument  : n/a

=cut

sub out_fasta {
  my($self) = @_;
  my($str,$i);

  # First, we have to split up our sequence into lines. We'll do this
  # by sticking a "\n" character every 60 characters in our sequence.
  # ?? Note: this is a really, particualrly costly way of doing this.
  # it would be faster to split into chunks of 60 and then return
  # the string which joins these.
  $str = $self->{"seq"};
  for ($i = 60; $i < length($str); $i += 60+1) {
        substr($str,$i,0) = "\n";
  }

  # Now, we return the result. We can't forget, of course, to put our
  # id on the top.
  return (">$self->{\"id\"} $self->{\"desc\"}\n$str\n");
}


#_______________________________________________________________________

=head2 alphabet_ok

 Title     : alphabet_ok
 Usage     : $myseq->alphabet_ok;
 Function  : Checks the sequence for presence of any characters
           : that are not considered valid members of the genetic
           : alphabet. In addition to the standard genetic alphabet
           : (see documentation), "?" and "-" characters are
           :  considered valid.
           :
 Example   : if($myseq->alphabet_ok) { print "OK!!\n"; }
           :     else { print "Not OK! \n"; }
           :
 Note      : Does not handle '\' characters in sequence robustly
           :
 Returns   : 1 if OK / 0 if not OK
 Argument  : none

=cut

#_______________________________________________________________________'
 
sub alphabet_ok {
    my($self) = @_;

    my($seq) = $self->{"seq"};
    $seq =~ tr/a-z/A-Z/;
    
    ## Make string containing largest possible appropriate alphabet
    ## SAC: Added ability to use a 'strict' alphabet that does not
    ##      allow ambiguity codes. Relying on the strict() method
    ##      inherited from Bio::Root::Object.pm.
    my ($al);
    if($self->strict) {
	$al = join("",@{$Alphabets_strict{$self->_monomer . "GpMg"}});
    } else {
	$al = join("",@{$Alphabets{$self->_monomer . "GpMg"}});
    }

    #CD: This code works ok for me but it fails to deal with the
    #CD: possibility that the sequence may contain the backslash
    #CD: "\" character. That char will then escape the next char
    #CD: in the internal regex that gets built.
    #CD: Example: sequence "AC\GT" will pass the alphabet test

    ##Add backslash escape to the ? and - alphabet characers
    ##(this is needed inside the regular expression)
    $al =~ s/\?/\\?/;
    $al =~ s/\-/\\-/;

    #print "Sequence being checked =$seq=\n";
    #print "Alphabet being checked =$al=\n";

    ## SAC: Corrected error. Was returning 0 if the seq was okay.
    ##Look for non-alphabet characters via regex
    if($seq =~ /[^$al]/i) { return 0; } # not OK  
    else { return 1 ; }                 # OK
 
}


#_______________________________________________________________________

=head2 alphabet

 Title     : alphabet
 Usage     : $myseq->alphabet;
 Function  : Returns the characters in the alphabet in use for the sequence.
 Example   : print "Alphabet: ".$myseq->alphabet;
 Returns   : string containing alphabet characters
 Argument  : none

=cut

sub alphabet {
# SAC: new method.
    my($self) = @_;

    my ($al);
    if($self->strict) {
	$al = join("",@{$Alphabets_strict{$self->_monomer . "GpMg"}});
    } else {
	$al = join("",@{$Alphabets{$self->_monomer . "GpMg"}});
    }

    $al;
}

=head2 GCG_checksum

 Title     : GCG_checksum
 Usage     : $myseq->GCG_checksum;
 Function  : returns a gcg checksum for the sequence
 Example   : 
 Returns   : 
 Argument  : none

=cut
 
sub GCG_checksum {
    my $self = shift;
    my $seq;
    my $index = 0;
    my $checksum = 0;
    my $char;


    $seq = $self->seq();
    $seq =~ tr/a-z/A-Z/;

    foreach $char ( split(/[\.\-]*/, $seq)) {
	$index++;
	$checksum += ($index * (unpack("c",$char)));
	if( $index ==  57 ) {
	    $index = 0;
	}
    }

    return ($checksum % 10000);
}


# Stubs for AUTOLOADED methods

sub copy;
sub revcom;
sub complement;
sub reverse;
sub Dna_to_Rna;
sub Rna_to_Dna;
sub translate;
sub dump ;
sub out_bad;
sub out_GCG;
sub out_zuker;
sub out_msf;
sub out_primer;
sub out_pir;
sub out_genbank;
sub out_nbrf;
sub out_gcgseq;
sub out_gcgref;
sub out_ig;
sub out_strider;
sub parse_unknown;
sub parse_bad;
sub version;

#
# These two lines were used for autoloading.
# (See Autoloading notes above and below for more info).
#1;
#__END__



#######################################################################
#                        AUTOLOADED METHODS  
#######################################################################
## SAC: Autoloading methods that are not always needed.

#### Manipulation Methods
####

#=head2 ## METHODS FOR SEQUENCE MANIPULATION ##

#_______________________________________________________________________

=head2 copy

 Title     : copy
 Usage     : $copyOfObj = $mySeq->copy;
 Function  : Returns an identical copy of the object.
 Example   :
 Returns   : Bio::Seq object ref.
 Argument  : n/a

=cut

sub copy {
  my($self) = @_;
  my(%dup);

  # changes suggested by David Curiel. Done by EB.

  $dup{"seq"}     = $self->{"seq"};
  $dup{"id"}      = $self->{"id"};
  $dup{"desc"}    = $self->{"desc"};
  $dup{"names"}   = {%{$self->{"names"}}}; # copied the name hash EB.
  $dup{'start'}   = $self->{'start'};
  $dup{'end'}     = $self->{'end'};
  $dup{"type"}    = [@{$self->{"type"}}]; # copied the type hash (EB)
  $dup{"ffmt"}    = $self->{"ffmt"};
  $dup{"descffmt"}= $self->{"descffmt"};

  return bless \%dup, ref($self);  #perl magic to support inheritance
}

#_______________________________________________________________________

=head2 revcom

 Title       : revcom
 Usage       : $reverse_complemented_seq = $mySeq->revcom;
 Function    : Returns a char string containing the reverse
             : complement of a nucleotide object sequence
 Example     : $reverse_complemented_seq = $mySeq->revcom;
 Source      : Guts from Jong's <jong@mrc-lmb.cam.ac.uk>
             : library of molbio perl routines
 Note        :
             : The letter codes and compliment translations
             : are those proposed by IUB (Nomenclature Committee,
             : 1985, Eur. J. Biochem. 150; 1-5) and are also
             : used by the GCG package. The IUB/GCG letter codes
             : for nucleotide ambiguity are compatible with
             : EMBL, GenBank and PIR database formats but are
             : *NOT* compatible with Stadem/Sanger ambiguity
             : symbols. Staden/Sanger use different symbols to
             : represent uncertainty and frame abiguity.
             :
             : Currently Staden/Sanger are not recognized
             : sequence types.
             :
             : GCG Documentation on sequence symbols:
 URL         : http://www.neb.com/gcgdoc/GCGdoc/Appendices/appendix_iii.html
             :
 Translation :
             : GCG/IUB    Meaning        Complement
             : ------------------------------------
             :  A            A                T
             :  C            C                G
             :  G            G                C
             :  T            T                A
             :  U            U                A
             :  M          A or C             K
             :  R          A or G             Y
             :  W          A or T             W
             :  S          C or G             S
             :  Y          C or T             R
             :  K          G or T             M
             :  V        A or C or G          B
             :  H        A or C or T          D
             :  D        A or G or T          H
             :  B        C or G or T          V
             :  X      G or A or T or C       X
             :  N      G or A or T or C       N
             :--------------------------------------
 Revision    : 0.01 / 3 Jun 1997
 Returns     : A new sequence object (fixed by eb)
               to get the actual sequence go
               $actual_reversed_sequence = $seq->revcom()->str()
 Argument    : n/a

=cut

#_______________________________________________________________________'

sub revcom {
  my($self,$start,$end)= @_;   # SAC: Added slicing ability 
#  my($self)= @_;  
  my($seq,$revseq);

#  print "revcom: requested range: $start, $end\n";

  # CD: Some type of check be made here to make
  # CD: sure the sequence is nucleotide.

  # SAC: copied code from str() to permit slicing.
  my($firstIndex,$startindx);

  $firstIndex = $self->start;
  $startindx = $firstIndex; 

  # Make sure $start,$end are in range, and set them to default if they're
  # not supplied
  if (defined($start)) {
      if ($start < $startindx) {
          $self->warn("Requested location $start of a string starting at $startindx,");
      }
  }
  else {
      $start = $startindx;
  }
  if (defined($end)) {
      if (($end - $startindx) > length($self->{"seq"}) - 1 ) {
          $self->warn("Requested location $end is beyond the end of the string,");
      }
  }
  else {
      $end = length($self->{"seq"}) + $startindx - 1;
  }

#  printf "revcom slice: %d, len = %d", $start-$startindx, $end-$start+1; <STDIN>;

  $seq = substr($self->{"seq"}, $start-$startindx, $end-$start+1);

#   $seq = $self->{"seq"};   ## SAC: Previous code (no slicing).

  # CD: If a sequence format uses different 
  # CD: symbols for nucleotide ambiguity than 
  # CD: GCG/IUB than this will certainly break.

  $seq =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
  $revseq = CORE::reverse $seq;

  my ($out,$id);
  $id = $self->id();

  ## CD: Added "" chars around id,seq,type params so that perl -w does
  ##     not complain about "ambigiguous use of id resolved to "id" etc. etc.

  $out = Bio::Seq->new('-id'=>"$id.revcom", '-seq'=>$revseq, '-type'=>'Dna' ); 
  return $out;
}


#_______________________________________________________________________

=head2 complement

 Title       : complement
 Usage       : $complemented_seq = $mySeq->compliment;
 Function    : Returns a char string containing 
             : the complementary sequence (eg; other strand)
             : of the original sequence. The translation method
             : is identical to revcom() but the nucleotide order
             : is not reversed. 
             :
             : To be honest *most* of the time you will want
             : to use revcom not this. Be careful!
             :
 Example     :  $complemented_seq = $mySeq->complement;
             :
 Source      : Guts from Jong's <jong@mrc-lmb.cam.ac.uk>
             : library of molbio perl routines
 Note        :
             : The letter codes and complement translations
             : are those proposed by IUB (Nomenclature Committee,
             : 1985, Eur. J. Biochem. 150; 1-5) and are also
             : used by the GCG package. The IUB/GCG letter codes
             : for nucleotide ambiguity are compatible with
             : EMBL, GenBank and PIR database formats but are
             : *NOT* compatible with Stadem/Sanger ambiguity
             : symbols. Staden/Sanger use different symbols to
             : represent uncertainty and frame abiguity.
             :
             : Currently Staden/Sanger are not recognized
             : sequence types.
             :
             : GCG Documentation on sequence symbols:
 URL         : http://www.neb.com/gcgdoc/GCGdoc/Appendices
             : /appendix_iii.html
             :
 Translation :
             : GCG/IUB    Meaning        Complement
             : ------------------------------------
             :  A            A                T
             :  C            C                G
             :  G            G                C
             :  T            T                A
             :  U            U                A
             :  M          A or C             K
             :  R          A or G             Y
             :  W          A or T             W
             :  S          C or G             S
             :  Y          C or T             R
             :  K          G or T             M
             :  V        A or C or G          B
             :  H        A or C or T          D
             :  D        A or G or T          H
             :  B        C or G or T          V
             :  X      G or A or T or C       X
             :  N      G or A or T or C       N
             :--------------------------------------
             :
 Revision    : 0.01 / 6 Dec 1996
 Returns     : char string
 Argument    : n/a

#_______________________________________________________________________'

=cut

sub complement {
  my($self)= @_;
  my($seq);
  
  # CD: See notes in revcom() about checking for nucleotide sequence
  # CD: and dealing with non-GCG/IUB ambiguity symbols.

  $seq = $self->{"seq"};
  $seq =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;

 return $seq;
}



=head2 reverse

 Title     : reverse
 Usage     : $reversed_seq = $mySeq->reverse;
 Function  : Returns a char string containing the
           : reverse of the object sequence
           :
           : Does *NOT* complement it. If you want
           : the other strand, use $mySeq->revcom()
           : 
 Example   :  $reversed_seq = $mySeq->reverse;
           :
 Revision  : 0.01 / 6 Dec 1996
 Returns   : char string
 Argument  : n/a

=cut

sub reverse {
  my($self)= @_;
  my($seq);

  $seq = $self->{"seq"};
  scalar CORE::reverse $seq;
}

#_______________________________________________________________________

=head2 Dna_to_Rna

 Title     : Dna_to_Rna
 Usage     : $translated_seq = $mySeq->Dna_to_Rna;
 Function  : Returns a char string containing the
           : Rna translation of the Dna nucleotide sequence
           : (Replaces T with U)
           : 
 Example   : $translated_seq = $mySeq->Dna_to_Rna;
           :
 Source    : modified from Jong's <jong@mrc-lmb.cam.ac.uk>
           : library of molbio perl routines
           :
 Revision  : 0.01 / 6 Dec 1996
 Returns   : char string
 Argument  : n/a

=cut

#_______________________________________________________________________'

sub Dna_to_Rna {
  my($self)=@_;
  
  ## CD: This is a simple substitution from T -> U
  ## CD: so we shouldnt have to do intricate error checking
  ## CD: here. Right now we only carp if the sequence is
  ## CD: explicitly an amino acid seq.

  $self->throw("Can't translate an amino acid sequence to Rna.") if($self->_monomer eq "3");

  my($seq)  = $self->{"seq"};  # get sequence

  ##Quietly deal with capitalization
  $seq =~ s/T/U/g;             # change (T to U)
  $seq =~ s/t/u/g;             # change (t to u)

  return $seq;
}

#_______________________________________________________________________'

=head2 Rna_to_Dna

 Title     : Rna_to_Dna
 Usage     : $translated_seq = $mySeq->Rna_to_Dna;
 Function  : Returns a char string containing the
           : Dna translation of the Rna nucleotide sequence
           : (Replaces U with T)
           : 
 Example   : $translated_seq = $mySeq->Rna_to_Dna;
           :
 Revision  : 0.01 / 16 MAR 1997
 Returns   : char string
 Argument  : n/a

=cut

sub Rna_to_Dna {
  my($self)=@_;

  ## CD: This is a simple substitution from U -> T
  ## CD: so we shouldnt have to do intricate error checking
  ## CD: here. Right now we only carp if the sequence is
  ## CD: explicitly an amino acid seq.

  $self->throw("Can't translate an amino acid sequence to Dna.") if($self->_monomer eq "3");

  my($seq)  = $self->{"seq"};  # get sequence

  ##Quietly deal with capitalization
  $seq =~ s/U/T/g;             # change (U to T)
  $seq =~ s/U/t/g;             # change (u to t)

  $seq;                        # return value
}

#_______________________________________________________________________

=head2 translate

 Title     : translate
 Usage     : 
 Function  : Returns a new Bio::Seq object with the protein
           : translation from this sequence
           :
           : "*" is the default symbol for a stop codon
           : "X" is the default symbol for an unknown codon
           :
 Example   : $translation = $mySeq->translate;
           :   -or- with user defined stop/unknown codon symbols:
           : $translation = $mySeq->translate($stop_symbol,$unknown_symbol);
           : 
 Source    : modified from Jong's <jong@mrc-lmb.cam.ac.uk>
           : library of molbio perl routines
           :
 To-do     : - allow named parameters (just like new and out_GCG )
           : - allow "frame" parameter to pick translation frame
           :
 Revision  : 0.01 / 6 Dec 1996
 Returns   : new Sequence object. Its id is the original id.trans
 Argument  : n/a

=cut

#_______________________________________________________________________'

sub translate {
  my($self) = shift;
  my($stop,$unknown) = @_;
  my($i,$len,$output) = (0,0,'');
  my($codon)   = "";

  my($seq) = $self->{"seq"};
  
  ## User can pass in symbol for stop and unknown codons
  unless(defined($stop))    { $stop    = "*"; }
  unless(defined($unknown)) { $unknown = "X"; }

  ##Error if monomer is "Amino"
  $self->throw("Can't translate an amino acid sequence.") if($self->_monomer eq "3");


  #If sequence monomer is Dna, we should first pipe it
  #through Dna_to_Rna() to change any T's to U's
  if($self->_monomer eq "1") { $seq = $self->Dna_to_Rna;} 

  for($len=length($seq),$seq =~ tr/a-z/A-Z/,$i=0; $i<($len-2) ; $i+=3) {
    $codon = substr($seq,$i,3);

    # would this be easier with a hash system (?) EB

    if   ($codon =~ /^UC/)     {$output .= 'S'; }       # Serine
    elsif($codon =~ /^UU[UC]/) {$output .= 'F'; }       # Phenylalanine
    elsif($codon =~ /^UU[AG]/) {$output .= 'L'; }       # Leucine
    elsif($codon =~ /^UA[UC]/) {$output .= 'Y'; }       # Tyrosine
    elsif($codon =~ /^UA[AG]/) {$output .= $stop; }     # Stop
    elsif($codon =~ /^UG[UC]/) {$output .= 'C'; }       # Cysteine
    elsif($codon =~ /^UGA/)    {$output .= $stop; }     # Stop
    elsif($codon =~ /^UGG/)    {$output .= 'W'; }       # Tryptophan
    elsif($codon =~ /^CU/)     {$output .= 'L'; }       # Leucine
    elsif($codon =~ /^CC/)     {$output .= 'P'; }       # Proline
    elsif($codon =~ /^CA[UC]/) {$output .= 'H'; }       # Histidine
    elsif($codon =~ /^CA[AG]/) {$output .= 'Q'; }       # Glutamine
    elsif($codon =~ /^CG/)     {$output .= 'R'; }       # Arginine
    elsif($codon =~ /^AU[UCA]/){$output .= 'I'; }       # Isoleucine
    elsif($codon =~ /^AUG/)    {$output .= 'M'; }       # Methionine
    elsif($codon =~ /^AC/)     {$output .= 'T'; }       # Threonine
    elsif($codon =~ /^AA[UC]/) {$output .= 'N'; }       # Asparagine
    elsif($codon =~ /^AA[AG]/) {$output .= 'K'; }       # Lysine
    elsif($codon =~ /^AG[UC]/) {$output .= 'S'; }       # Serine
    elsif($codon =~ /^AG[AG]/) {$output .= 'R'; }       # Arginine
    elsif($codon =~ /^GU/)     {$output .= 'V'; }       # Valine
    elsif($codon =~ /^GC/)     {$output .= 'A'; }       # Alanine
    elsif($codon =~ /^GA[UC]/) {$output .= 'D'; }       # Aspartic Acid
    elsif($codon =~ /^GA[AG]/) {$output .= 'E'; }       # Glutamic Acid
    elsif($codon =~ /^GG/)     {$output .= 'G'; }       # Glycine
    else {$output .= $unknown; }                        # Unknown Codon
  }

  my($out,$id);
  $id = $self->id();

  ## CD: Added "" chars around seq,type,id param so that 'perl -w' does not
  ##     warn about ambiguous usage...

  $out = Bio::Seq->new( '-id' => "$id.trans" , '-seq' => $output, '-type' => 'Amino' );
  return $out;
}


#_______________________________________________________________________

=head2 dump

 Title     : dump
 Usage     : @results = $mySeq->dump; -or- 
           : $results = $mySeq->dump;
           :
 Function  : Returns a formatted array or string (depending on how it
           : is invoked) containing the contents of a 
           : Bio::Seq object. Useful for debugging
           :
           : ***This is used by Chris Dagdigian for debugging ***
           : ***Probably should be removed before distribution***
           :
 Example   :  @results = $mySeq->dump;
           :  foreach(@results){print;}
           :     -or-
           :  print $myseq->dump;
           :
 Returns   : Array or string depending on value of wantarray
 Argument  : n/a

=cut

sub dump {
my($self) = @_;
my(@result,$monomer,$origin,$name_ref);

$name_ref = $self->names;

push(@result,"\nID       : ",$self->id); 
push(@result,"\nType     :",$self->type);
push(@result,"\nMonomer  :",$self->_monomer);
push(@result,"\nOrigin   :",$self->origin);
push(@result,"\nSeq      :",$self->seq);
push(@result,"\nFormat   :",$self->ffmt);
push(@result,"\nStart    :",$self->start);
push(@result,"\nEnd      :",$self->end);
push(@result,"\nDesc     :",$self->desc);
push(@result,"\nDescffmt :",$self->descffmt);
push(@result,"\nNames    -\n");
 foreach(keys %$name_ref) { 
    push(@result," $_  = $name_ref->{$_}\n"); 
  }
push(@result,"\n");

return wantarray ? @result : join("",@result);

}


#_______________________________________________________________________

=head2 out_bad

 Title     : out_bad()
 Usage     : out_bad;
 Function  : Throws a fatal error if we don't know the output format.
 Example   : $self->out_bad;
 Returns   : n/a
 Argument  : n/a

=cut

#----------------------------------------------------------------------'

sub out_bad {
  my($self) = @_;

  $self->throw("Can't write sequence format $self->{\"ffmt\"}"); 

  return 0;
} 


#_______________________________________________________________________

=head2 out_primer

 Title     : out_primer()
 Usage     : $formatted_seq = $myseq->out_primer;
           : @formatted_seq = $myseq->out_primer;
           :
           : print $myseq->out_primer(-id=>'New ID',
           :                          -header=>'This is my header');
           :
 Function  : outputs a sequence in primer format
           :
 Note      : Not a supported output type -  (cant be invoked via layout)
           : Use at your own risk :)
           : 
 Example   : see usage
           :
 Revision  : 0.01 / 20 Dec 1996
 Returns   : string or list, depending on how it is invoked
 Argument  : named list parameters for "id" and "header" are alowed

=cut

sub out_primer {
 my($self,@params) = @_;
 my($i,$j,$len,@out,$in_header);
 my($in_id,$comment,$ID,$seq,$id);

 if(defined(@params)) {
  ($in_header,$in_id) =
      $self->_rearrange([qw(HEADER ID)],@params);
  }

  ## Set default values that may get overwritten
  $seq = $self->{"seq"};
  $id  = $self->{"id"};
  $comment = $self->{"desc"};
  $len = length($seq);

  ## Deal with user arguments (overwrite if they exist)
  if(defined($in_id))     { $id      = $in_id;    }
  if(defined($in_header)) { $comment = $in_header;}

 $out[$i++] = sprintf("*seq: %8s\n",$id);
 $out[$i++] = sprintf("%s", length($comment) > 8 ? "# ".$comment."\n" : "");
 $out[$i++] = sprintf("%s",$id ? "# ". $id . "\n" : "");
 push(@out,"");

  #Format the sequence
  $i = $#out + 1;
  for($j = 0 ; $j < $len ; ) {
    if( $j % 50 == 0) { push(@out,""); } #This shuts strict up ???
    $out[$i] .= sprintf("%s",substr($seq,$j,10));  $j += 10;
     if( $j < $len && $j % 50 != 0 ) {
       $out[$i] .= " ";
     }elsif($j % 50 == 0 ) {
       $out[$i++] .= "\n";
     }                           
  }
  if($j % 50 != 0 ) { $out[$i] .= "\n";  }
  $out[$i] .= "\n";

 return wantarray ? @out : join("",@out);
}

#_______________________________________________________________________

=head2 out_pir

 Title     : out_pir()
 Usage     : $formatted_seq = $myseq->layout("PIR");
           : $formatted_seq = $myseq->out_pir;
           : @formatted_seq = $myseq->out_pir;
           :
           : print $myseq->out_pir(-title=>'New TITLE',
           :                       -entry=>'New ENTRY',
           :                       -acc=>'User defined accession',
           :                       -date=>'User defined date',
           :                       -reference=>'User defined ref info');
           :
 Function  : Returns a string or an array depending on how it
           : is invoked. Can be easily accessed via the layout()
           : method, or if more output control is desired it can
           : be called directly with the folowing named parameters:
           :
           :  -entry      PIR entry
           :  -title      PIR title
           :  -acc        user defined accession number
           :  -reference  user defined reference
           :  -date       user defined date/time info
           :
           : All named parameters will take precedance over any
           : default behavior. When there are no user arguments,
           : the default output is as follows:
           :
           : PIR 'ENTRY'     = sequence object "id" field
           : PIR 'TITLE'     = sequence object "desc" field
           : PIR 'DATE'      = curent date/time
           : PIR 'ACC'       = not used in default output
           : PIR 'REFERENCE' = not used in default output
           :
 Note      : Not tested stringently.
           :
 WARNING   : Does not deal with numbering issue
           :
 To-do     : - Allow user to pass in hash of additional fields/values
           : - Deal with numbering issue
           :
 Example   : see usage
           :
 Revision  : 0.02 / 12 Jan 1997
 Returns   : string or list, depending on how it is invoked
 Argument  : named list parameters are allowed, see above

=cut

sub out_pir {
    my($self,@params) = @_;
    my(@out,$len,$i,$j,$cnt);
    my($in_ent,$in_title,$in_acc,$in_date,$in_ref);
    my($seq,$id,$c_t,$entry,$title,$ref,$acc);

 if($#params >= 0) {
  ($in_ent,$in_title,$in_acc,$in_date,$in_ref) =
      $self->_rearrange([qw(ENTRY
			    TITLE
			    ACC
			    DATE
			    REFERENCE)],
                        @params);
}

  ## Set default values that may get overwritten
  $seq   = $self->{"seq"};
  $entry = $self->{"id"};
  $title = $self->{"desc"};
  $len   = length($seq);

    $c_t = localtime;

  # Deal with user arguments (overwrite default info...)
  if(defined($in_ent))    { $entry   = $in_ent;   }
  if(defined($in_title))  { $title   = $in_title; }
  if(defined($in_acc))    { $acc     = $in_acc;   }
  if(defined($in_ref))    { $ref     = $in_ref;   }  
  if(defined($in_date))   { $c_t     = $in_date;  }

    $out[$i++] = sprintf("ENTRY        %s\n",$entry);
    $out[$i++] = sprintf("TITLE        %s\n",$title);
    if(defined($in_acc))  { $out[$i++] = sprintf("ACCESSION    %s\n",$acc);}
    $out[$i++] = sprintf("DATE         %s\n",$c_t);
    if(defined($in_ref)) { $out[$i++] = sprintf("REFERENCE    %s\n",$ref); }
    $out[$i++] = sprintf("SEQUENCE     %s\n");
    $out[$i++] = sprintf("                5         10        15        20        25        30\n");

    for($j=1; $seq && $j < $len; $j +=30) {
	$out[$i++] = sprintf("%7d ",$j);
	$out[$i++] = sprintf("%s\n", join(" ",split(//,substr($seq, $j-1,length($seq) < 30 ? length($seq) : 30))) );
    }
   $out[$i++] = sprintf("///\n");

   return wantarray ? @out : join("",@out);
}

#_______________________________________________________________________

=head2 out_genbank

 Title     : out_genbank()
 Usage     : $formatted_seq = $myseq->out_genbank;
           : @formatted_seq = $myseq->out_genbank;
           : print $myseq->out_genbank(-id=>'New ID',
           :                           -def=>'User defined definition',
           :                           -acc=>'User defined accession',
           :                           -origin=>'User defined origin info',
           :                           -spacing=>'single',
           :                           -caps=>'up',
           :                           -date=>'DATE GOES HERE',
           :                           -type=>'mRna');
           :   
 Function  : Returns a GenBank formatted sequence array or string
           : depending on the value of wantarray when invoked via layout(). 
           : If more control is desired over output format, out_genbank() 
           : can be addressed directly with the following named parameters:
           :
           : def          - Sequence definition information
           : acc          - Sequence accession number
           : origin       - Sequence origin information
           : id           - short name 
           : date         - new date info
           : type         - sequence type (Dna, mRna, Amino, etc.)
           : spacing      - "single" or "double" sequence line spacing
           : caps         - "up" or "down" sequence capitalization
           :
           : When invoked via layout() or called directly with no 
           : arguments, the following default behaviours apply:
           :  DATE = Current date and time
           :  DEFINITION = object's description field
           :  ID = object's ID field
           :  SPACING = single
           :
           : All named parameters must be strings. Passed in parameters will
           : always take precedence over any fields with default settings.
           :
 Note      : Format not stringently tested for accuracy. Sequence is numbered
           : according to the integer specified in the object 'start' field
           : but the implementation has not been robustly tested.
           :
 To-do     : - allow user hash reference for additional format fields
           :
 Example   : see usage
           :
 Revision  : 0.02 / 12 Jan 1997
 Returns   : string or list, depending on how it is invoked
 Argument  : named list parameters are allowed, see above

=cut

sub out_genbank {
 my($self,@params) = @_;
 my($id,$comment,$len,$type,$seq);
 my($cnt,$sum,$i,$j,$tmp,$offset);
 my($c_t,@out,$origin);
 my($in_def,$in_locus,$in_date,$in_acc,$in_origin,$in_type,$caps,$spacing);
 my($spacer) = "";

 if($#params >= 0) {
  ($in_def,$in_locus,$in_acc,$in_origin,$in_date,$in_type,$caps,$spacing) =
      $self->_rearrange([qw(DEF
			    ID
			    ACC
			    ORIGIN
			    DATE
			    TYPE
			    CAPS
			    SPACING)],
                        @params);
  }

  ## Get current date
  $c_t = localtime; 

  ## Set default values that may get overwritten
  $seq = $self->{"seq"};
  $id  = $self->{"id"};
  $comment = $self->{"desc"};
  $origin = $self->origin;
  $len = length($seq);
  $type = $TypeSeq{$self->_monomer};

  # Deal with user arguments
  if(defined($in_type))   { $type    = $in_type; $type =~ tr/a-z/A-z/;  }
  if(defined($in_locus))  { $id      = $in_locus;    }
  if(defined($in_origin)) { $origin  = $in_origin;   }
  if(defined($in_def))    { $comment = $in_def;}
  if(defined($in_date))   { $c_t     = $in_date;  }

  if(defined($spacing))   { $spacing =~ tr/a-z/A-Z/;
                           if($spacing eq "DOUBLE") { $spacer = "\n"; }
                           else {if($spacing eq "SINGLE") { $spacer = ""; }}
                          }

  if(defined($caps))      { $caps =~ tr/a-z/A-Z/;
                           if($caps eq "UP") { $seq =~ tr/a-z/A-Z/; }
                           else {if($caps eq "DOWN") { $seq =~ tr/A-Z/a-z/; }}
			}

  $offset=0;
  #Set the offset if we have any non-standard numbering going on
  if($self->start < 0)   { $offset = ( 0 + $self->start); }
  if($self->start >= 1)  { $offset = $self->start;}  
  if($self->start == 0)  { $offset = -1;}

  $sum=0; #Need this to shut strict() up

  #Output the sequence header info
  push(@out,"LOCUS\t$id\t$len\t$type\t$c_t\n");                        
  push(@out,"DEFINITION $comment\n");                        
  if(defined($in_acc)) { push(@out,"ACCESSION $in_acc\n"); }                        
  push(@out,"ORIGIN $origin\n");                        

  #Format the sequence
  $len = length($seq);                     
  $i = $#out + 1;
  for($j = 0 ; $j < $len ; ) {
    if( $j % 50 == 0) {
      $out[$i] = sprintf("%8d  ",($j+$offset)); #numbering 
    }
    $out[$i] .= sprintf("%s",substr($seq,$j,10));  $j += 10;
    if( $j < $len && $j % 50 != 0 ) {
      $out[$i] .= " ";
    }elsif($j % 50 == 0 ) {
      $out[$i++] .= "\n";
      if(defined($spacer)) { $out[$i++] = $spacer;} 
    }                           
  }
  local($^W) = 0;
  if($j % 50 != 0 ) { $out[$i] .= "\n"; }
  $out[$i] .= "\n//\n";

return wantarray ? @out : join("",@out);

}


#_______________________________________________________________________

=head2 out_GCG

 Title    : out_GCG
 Usage    : $formatted_seq = $mySeq->layout("GCG"); 
          : @formatted_seq = $mySeq->layout("GCG");
          : 
          : print $myseq->out_GCG(-id=>'New ID',
          :                      -spacing=>'single',
          :                      -caps=>'up',
          :                      -date=>'DATE GOES HERE',
          :                      -header=>'This is a user submitted header',
          :                      -type=>'n');
          :   
 Function : Returns a GCG formatted sequence array or string
          : depending on the value of wantarray when invoked via layout(). 
          : If more control is desired over output format, out_GCG() 
          : can be addressed directly with the following named parameters:
          :
          : header       - first line(s) of formatted sequence
          : id           - short name that appears before 'Length:' field
          : date         - overwrite default date info
          : type         - can be "N" or "P", for nucleotide/protein
          : spacing      - "single" or "double" sequence line spacing
          : caps         - "up" or "down" sequence capitalization
          :
          : When invoked via layout() or called directly with no 
          : arguments, the following default behaviours apply:
          :  DATE = Current date and time
          :  DEFINITION = object's description field
          :  ID = object's ID field
          :  SPACING = single
          :         
          : All named parameters must be strings. Passed in parameters will
          : always take precedence over any fields with default settings.
          :
 Example  :  
 Output   :
          :Sample Bio::Seq sequence
          : sample Length: 240  Wed Nov 27 13:24:28 EST 1996  Type: N Check: 5371  ..
          :
          :       1  aaaacctatg gggtgggctc tcaagctgag accctgtgtg cacagccctc
          :      51  tggctggtgg cagtggagac gggatnnnat gacaagcctg ggggacatga
          :     101  ccccagagaa ggaacgggaa caggatgagt gagaggaggt tctaaattat
          :     151  ccattagcac aggctgccag tggtccttgc ataaatgtat agagcacaca
          :     201  ggtgggggga aagggagaga gagaagaagc cagggtataa
          :
          :
 Note     : GCG formatted sequences contain a "Type:" field.
          : If Type cannot be internally determined and no
          : Type name-parameter is passed in then the Type: 
          : field is not printed.
          :
 Warning  : Unconventional numbering offsets may not
          : be robustly handled
          :
 Revision : 0.06 / 12 Jan 1997
 Source   : Found guts of this code on bionet.gcg, unknown author
 Returns  : Array or String
 Argument : n/a

=cut

sub out_GCG {
 my($self,@params) = @_;
 my($id,$comment,$len,$type,$seq);
 my($cnt,$sum,$i,$j,$tmp,$offset);
 my($c_t,@out);
 my($in_header,$in_id,$in_date,$in_type,$caps,$spacing);
 my($spacer) = "";

 if(defined(@params)) {
  ($in_header,$in_id,$in_date,$in_type,$caps,$spacing) =
      $self->_rearrange([qw(HEADER
			    ID
			    DATE
			    TYPE
			    CAPS
			    SPACING)],
                        @params);
  }


  $c_t = localtime;

  ## Set default values that may get overwritten
  $seq = $self->{"seq"};
  $id  = $self->{"id"};
  $comment = $self->{"desc"};
  $len = length($seq);
  $type = "";


  #Deal with the GCG format Type field 
   if($self->_monomer eq "3") { $type="Type: P";}
   else { if(($self->_monomer eq "1") || ($self->_monomer eq "2")) { $type="Type: N"; }
   } 

  # The default "N" or "P" type has been set if possible. After dealing with any
  # overiding user arguments, we can test it and carp if undefined

  # Deal with user arguments
  if(defined($in_type))   { $in_type =~ tr/a-z/A-Z/; $type    = "Type: $in_type"; }
  if(defined($in_id))     { $id      = $in_id;    }
  if(defined($in_header)) { $comment = $in_header;}
  if(defined($in_date))   { $c_t     = $in_date;  }

  if(defined($spacing))   { $spacing =~ tr/a-z/A-Z/;
                           if($spacing eq "DOUBLE") { $spacer = "\n"; }
                           else {if($spacing eq "SINGLE") { $spacer = ""; }}
                          }

  if(defined($caps))      { $caps =~ tr/a-z/A-Z/;
                           if($caps eq "UP") { $seq =~ tr/a-z/A-Z/; }
                           else {if($caps eq "DOWN") { $seq =~ tr/A-Z/a-z/; }}
                          }

  # Test $type
  unless($type eq "Type: N" || $type eq "Type: P") { 
   $type = "";
   }

  $offset=0;
  #Set the offset if we have any non-standard numbering going on
  if($self->start < 0)   { $offset = ( 0 + $self->start); }
  if($self->start >= 1)  { $offset = $self->start;}  
  if($self->start == 0)  { $offset = -1;}

  $sum=0; #Need this to shut strict() up

  #Generate the GCG Checksum value
  for($i=0; $i<$len ;$i++) {             
    $cnt++;
    $sum += $cnt * ord(substr($seq,$i,1));
    ($cnt == 57) && ($cnt=0);
  }
  $sum %= 10000;

  #Output the sequence header info
  push(@out,"$comment\n");                        
  push(@out," $id Length: $len  $c_t  $type Check: $sum  ..\n\n");

  #Format the sequence
  $len = length($seq);                     
  $i = $#out + 1;
  for($j = 0 ; $j < $len ; ) {
    if( $j % 50 == 0) {
      $out[$i] = sprintf("%8d  ",($j+$offset)); #numbering 
    }
    $out[$i] .= sprintf("%s",substr($seq,$j,10));
    $j += 10;
    if( $j < $len && $j % 50 != 0 ) {
      $out[$i] .= " ";
    }elsif($j % 50 == 0 ) {
      $out[$i++] .= "\n";
      if(defined($spacer)) { $out[$i++] = $spacer;} 
    }                           
  }
  local($^W) = 0;
  if($j % 50 != 0 ) {
    $out[$i] .= "\n";
  }
  $out[$i] .= "\n";

return wantarray ? @out : join("",@out);

} # end of sub 



#_______________________________________________________________________


=head2 out_nbrf

 Title     : out_nbrf()
 Usage     : $self->layout("NBRF") or $self->out_nbrf
           :
 Function  : FORMAT NOT INTERNALLY IMPLEMENTED YET!!!
           :
           : If the ReadSeq wrapper Parse.pm apppears 
           : to be configured properly it is used
           : to generate the output. 
           :
           : If Parse.pm cannot be used then this code
           : carps out with an error message.
           :
 To-do     : write internal output code
           :
 Version   : 1.0 /  16 MAR 1997
 Example   : see Usage
 Returns   : FORMATTED STRING (wantarray is not used here!)
 Argument  : 

=cut

sub out_nbrf {
    my($self) = @_;
    my($formatted_seq) = "";
    my($seq) = $self->{"seq"};

    if(defined($READSEQ_EXISTS)) {
      $formatted_seq = &Bio::Parse::convert_from_raw(-sequence=>$seq,-fmt=>"NBRF");
      return $formatted_seq;
    }
    else {
	$self->throw("NBRF output format is not currently supported.");
    }
}

#_______________________________________________________________________

=head2 out_gcgseq

 Title     : out_gcgseq
 Usage     : out_gcgseq;
 Function  : Returns the sequence as a string in GCG_SEQ format.
 Example   : $self->out_gcgseq;
           :
 Returns   : string sequence in GCG_SEQ format
 Argument  : n/a
 Comments  : SAC: Derived from out_fasta().
           : GCG_SEQ is a format that looks alot like Fasta and is used
           : for building GCG sequence datasets (.seq files).
           : It also has some similarities to NBRF format.

=cut

sub out_gcgseq {
  my $self = shift;
  my($str,$i);

  $str = $self->str(@_);
  # Line length is not limited to 60. Max is probably 500 (as in NBRF).
  # 60 is just more readable.
  for ($i = 60; $i < length($str); $i += 60+1) {
        substr($str,$i,0) = "\n";
  }
  # Terminal '*' is optional.

  return (">>>>$self->{\"id\"}\n$self->{\"desc\"}\n$str\n");
}

#_______________________________________________________________________

=head2 out_gcgref

 Title     : out_gcgref
 Usage     : out_gcgref;
 Function  : Returns the sequence as a string in GCG_REF format.
 Example   : $self->out_gcgref;
           :
 Returns   : string sequence in GCG_REF format
 Argument  : n/a
 Comments  : SAC: Derived from out_gcgseq().
           : GCG_REF is a companion format for GCG_SEQ that is used
           : for building GCG sequence datasets (.ref files).
           : The .ref file is identical to .seq file but without the sequence.

=cut

sub out_gcgref {
  my($self) = @_;
  my($str,$i);

  return (">>>>$self->{\"id\"}\n$self->{\"desc\"}\n");
}

#_______________________________________________________________________


=head2 out_ig

 Title     : out_ig()
 Usage     : $self->layout("IG") or $self->out_ig
           :
 Function  : FORMAT NOT INTERNALLY IMPLEMENTED YET!!!
           :
           : If the ReadSeq wrapper Parse.pm apppears 
           : to be configured properly it is used
           : to generate the output. 
           :
           : If Parse.pm cannot be used then this code
           : carps out with an error message.
           :
 To-do     : write internal output code
           :
 Version   : 1.0 /  16 MAR 1997
 Example   : see Usage
 Returns   : FORMATTED STRING (wantarray is not used here!)
 Argument  : 

=cut

sub out_ig {
    my($self) = @_;
    my($formatted_seq) = "";
    my($seq) = $self->{"seq"};

    if(defined($READSEQ_EXISTS)) {
      $formatted_seq = &Bio::Parse::convert_from_raw(-sequence=>$seq,-fmt=>"IG");
      return $formatted_seq;
    }
    else {
	$self->throw("IG output format is not currently supported.");
    }
}

#_______________________________________________________________________


=head2 out_strider

 Title     : out_strider()
 Usage     : $self->layout("Strider") or $self->out_strider
           :
 Function  : FORMAT NOT INTERNALLY IMPLEMENTED YET!!!
           :
           : If the ReadSeq wrapper Parse.pm apppears 
           : to be configured properly it is used
           : to generate the output. 
           :
           : If Parse.pm cannot be used then this code
           : carps out with an error message.
           :
 To-do     : write internal output code
           :
 Version   : 1.0 /  16 MAR 1997
 Example   : see Usage
 Returns   : FORMATTED STRING (wantarray is not used here!)
 Argument  : 

=cut

sub out_strider {
    my($self) = @_;
    my($formatted_seq) = "";
    my($seq) = $self->{"seq"};

    if(defined($READSEQ_EXISTS)) {
      $formatted_seq = &Bio::Parse::convert_from_raw(-sequence=>$seq,-fmt=>"Strider");
      return $formatted_seq;
    }
    else {
	$self->throw("Strider output format is not currently supported.");
    }
}

#_______________________________________________________________________


=head2 out_zuker

 Title     : out_zuker()
 Usage     : $self->layout("Zuker") or $self->out_zuker
           :
 Function  : FORMAT NOT INTERNALLY IMPLEMENTED YET!!!
           :
           : If the ReadSeq wrapper Parse.pm apppears 
           : to be configured properly it is used
           : to generate the output. 
           :
           : If Parse.pm cannot be used then this code
           : carps out with an error message.
           :
 To-do     : write internal output code
           :
 Version   : 1.0 /  16 MAR 1997
 Example   : see Usage
 Returns   : FORMATTED STRING (wantarray is not used here!)
 Argument  : 

=cut

sub out_zuker {
    my($self) = @_;
    my($formatted_seq) = "";
    my($seq) = $self->{"seq"};

    if(defined($READSEQ_EXISTS)) {
      $formatted_seq = &Bio::Parse::convert_from_raw(-sequence=>$seq,-fmt=>"Zuker");
      return $formatted_seq;
    }
    else {
	$self->throw("Zuker output format is not currently supported.");
    }
}

#_______________________________________________________________________


=head2 out_msf

 Title     : out_msf()
 Usage     : $self->layout("MSF") or $self->out_msf
           :
 Function  : FORMAT NOT INTERNALLY IMPLEMENTED YET!!!
           :
           : If the ReadSeq wrapper Parse.pm apppears 
           : to be configured properly it is used
           : to generate the output. 
           :
           : If Parse.pm cannot be used then this code
           : carps out with an error message.
           :
 To-do     : write internal output code
           :
 Version   : 1.0 /  16 MAR 1997
 Example   : see Usage
 Returns   : FORMATTED STRING (wantarray is not used here!)
 Argument  : 

=cut

sub out_msf {
    my($self) = @_;
    my($formatted_seq) = "";
    my($seq) = $self->{"seq"};

    if(defined($READSEQ_EXISTS)) {
      $formatted_seq = &Bio::Parse::convert_from_raw(-sequence=>$seq,-fmt=>"MSF");
      return $formatted_seq;
    }
    else {
	$self->throw("MSF output format is not currently supported.\n");
    }
}





#_______________________________________________________________________

=head2 parse_unknown

 Title     : parse_unknown
 Usage     : parse_unknown($ent);
 Function  : tries to figure out the format of $ent and then
           : calls the appropriate function to parse it into $self->{"seq"}.
 Example   : $self->parse_unknown;
 Returns   : n/a
 Argument  : $ent : the rough multi-line string to be parsed

=cut

sub parse_unknown {
  my($self, $ent, $filename) = @_;
  if ($ent =~ /^>[ \t]*\S*[ \t]*.*?\n(?:\n|.)*?(?=\n>|\Z)/mg) {
      # >[ \t]*\S*  : The ``>''-sign, space or tab, and the id;
      # [ \t]*.*?\n : space or tab again, the description (arbitrary 
      #               letters until the next newline, so use non-greedy
      #               regexp (``.*?''), ``.*'' would gobble up everthing);
      # (?:\n|.)*?  : now arbitrary text including newlines, non-greedy,
      #               ``(?:...)'' just groups ``...'' together;
      # (?=\n>|\Z)  : up until either a newline followed by the ``>''-sign,
      #               or the end of the multiline string (so-called positive 
      #               lookahead)
    $self->parse_fasta($ent);
  }
  elsif ($ent =~ /.*/mg) {
      # currently in raw format, everything is accepted...
    $self->parse_raw($ent);
  }
  else {
        # some other weird format, so we call parse_bad.
    $self->parse_bad($ent);
  }

  return 1;

}


#_______________________________________________________________________

=head2 parse_bad

 Title     : parse_bad
 Usage     : parse_bad;
 Function  : complains of un-parsable sequence, last-ditch attempt via
           : Parse.pm if sequence is being read from a file.
           :
 Example   : $self->parse_bad;
 Returns   : n/a
 Argument  : n/a

=cut

sub parse_bad
 {
  my($self) = shift;
  my($ent, $filename) = @_;
  my(@lines,$head,$reply);

  if(defined($READSEQ_EXISTS)) {
      ## LAST DITCH ATTEMPT AT PARSING
      ## Use Parse.pm & ReadSeq to parse file

      my($reply) = &Bio::Parse::convert(-sequence=>$ent,
					-location=>$filename,
					-fmt=>'raw');
      $reply=~ s/\n+//g; $reply=~ s/\s+//g;

      unless($reply eq "") { #NO ALPHABET CHECKING HERE! 
                             $self->{"seq"} = "$reply";
                           }
      else { #Give up if no reply
             $self->throw("Can't parse sequence $ent, even tried Parse.pm.\n");
             return 0; 
	   }

   } else {
 
            ## just plain give up.
            $self->throw("Error: Cannot  parse this sequence\n $ent");
            return 0;
          }

  return 1;
}


#=head2 ## Misc. methods  ##

=head2 version

 Title     : version();
 Usage     : $myseq->version;
 Function  : prints Bio::Seq current version number

=cut

sub version {
  my($self) = @_;
  print "Bio::Seq Version is ", $Bio::Seq::VERSION, ".\n";
  return 1;
}


1;
__END__


## =head1 ## END  METHOD DOCS ##

=head1 Bio::Seq Guts

=head2 Sequence Object

 The sequence object is merely a reference to a hash containing
 all or some of the following fields...

 Field         Value
 --------------------------------------------------------------
 seq           the sequence
 
 id            a short identifier for the sequence
 
 desc          a description of the sequence, in descffmt file-format
 
 names         a hash of identifiers that relate to the sequence..
               these could be Database ID's, Accession #'s, URL's,
               pathnames, etc. Currently there is no set format
               for the names hash and no formal definition of databases 
               or names
 
 start         start in bio-coords of the first residue of the sequence

 end           end in bio-coords of the first residue of the sequence
 
 type          the sequence type. Is actually a 2 value list of format
               ["monomer","origin"] where monomer is one of the
               recognized sequence types and origin is a string
               description of the sequences' origin (mitochondrial, etc)
 
 ffmt          file-format for the sequence
 
 descffmt      file-format of the description string

=cut


MODIFICATION NOTES:
--------------------
0.051, 17 Feb 1999, sac:
   * Modified type() so that it always returns a string and does not
     store an 'unknown' string if the type is not set.

### Prior to using a central CVS system on bio.perl.org:
### ++++++++++++++++++++++++++++++++++++++++++++++++++++
###
### VERSION 0.050, 3 Sep 1998:
###
### -- Added start() and end() and deprecated numbering().
###    (Changes made by Ewan Birney).
###    Converted all calls to numbering() to start().
### -- Officially graduated to Seq.pm.
###
### VERSION 0.047, 15 Jul 1998:
###
### -- Bug fixed in str() that caused failure in bounds checking if
###    start does not begin with 1 (suggested by Tim Dudgeon)
###
### VERSION 0.046, 10 June 1998:
###
### -- Added & improved documentation, including internal hyperlinks, for
###    generating docs using pod2html in the Perl 5.004 release.
### -- Not autoloading commonly used/small parsing and outputting
###    methods (parse_raw, parse_fasta, parse_gcg, out_raw, out_fasta).
###
###  VERSION 0.045, 5 June 1998: 
###
###  -- Integrated into the incipient Bioperl framework: 
###       * Inherits from Bio::Root::Object.pm.
###       * Uses _new() and _rearrange() via inheritance.
###       * Changed carps to $self->warn() and croaks to $self->throw().
###  -- revcom() can operate on sequence slices.
###  -- Bug fixed in alphabet_ok().
###  -- Fixed up _rearrange() calls to use non-interpolated lists (qw()).
###  -- Changed &Parse:: calls to &Bio::Parse:: calls.
###  -- Noted a documentation bug referring to a non-existent method seq().
###     Changed docs to refer to getseq(). 
###  -- Added method seq() to the API (see previous note).
###     This is a duplication of getseq().
###     We should settle on one of these and label the other one deprecated.
###     seq() is simpler and easier to remember (eg., was it getseq or get_seq?).
###     However, getseq() is clearer that it is a 'get' not a 'set' call.
###     I favor seq(), since it is simpler and consistent with the other accessors
###     (str(), ffmt(), id()) and the naming of accessors in other modules.
###     Having two methods with different names that do the same thing is
###     a bad idea. For now, getseq() is the official method to keep this
###     version of Bio::Seq in line with previous versions.
###  -- Quoted data members {"seq"} etc. to prevent "Ambiguous use" compiler warnings.
###
###  -- Added AUTOLOADing for seldom-used methods to reduce start-up time.
###     Autoloading the specific 'out_' and 'parse_' methods and the DNA
###     manipulation methods (not needed for protein sequences).
###
###     Autoloading requires that this module be autosplit. If you have installed
###     this module using the standard 'perl Makefile.PL' procedure, autosplitting
###     is done automatically (verify this by checking for an auto/Bio/Seq/ directory 
###     in your perl module library). You can autosplit it manually by first uncommenting
###     the __END__ line and then running the following Perl script: 
###     (substituting the correct path for your perl lib)
###         #!/usr/bin/perl
###         use AutoSplit; 
###         autosplit("/users/me/perl/lib/Bio/Seq.pm", "/users/me/perl/lib/auto", 0, 1, 1);
###
###     To disable autoloading, comment out the __END__ line.
###
###     Autoloading raises the issue that perhaps there should be separate modules
###     such as Bio::Seq::Out.pm and Bio::Seq::DnaTools.pm. This will become
###     more of an issue as more native Perl code for parsing/outputting is
###     added to this module and it becomes large and monolithic.
###
###  -- Added strict version of the basic alphabets that do not allow ambiguity
###     codes. The default alphabets permit ambiguity codes.
###     To use the strict alphabets, include '-strict => 1' in the parameters
###     sent to new(). Or, after constructing the sequence object, try:
###       $seq->strict(1); $seq->alphabet_ok() or die "alphabet not okay.\n";
###     This change affects only alphabet_ok().
###  -- Throws exception during construction in strict mode if type is not defined.
###     This is a common error with occasionally important consequences.
###     (having separate subclasses for amino and nucleic acid seqs would
###     avoid this problem).
###  -- Added methods out_gcgseq() and out_gcgref() for producing GCG-style
###     multiple-sequence files used for building GCG datasets.
###  -- Added method alphabet() to return the alphabet in use for the seq.
###  -- Lower-cased all strings used as hash keys in %SeqForm and %SeqAlph,
###     as discussed in the Bioperl mailing list (see last paragraph):
###     http://www.uni-bielefeld.de/mailinglists/BCD/vsns-bcd-perl/9702/0022.html 
###  -- Enabled case-independence of user-submitted strings
###     for the following methods: 
###     parse(), ffmt(), type(), descffmt(), layout().
###  

