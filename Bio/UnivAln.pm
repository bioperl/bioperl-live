
# WHEN IN DOUBT, ALWAYS CHECK THE FOLLOWING URL FOR THE NEWEST VERSION:
# http://www.techfak.uni-bielefeld.de/bcd/Perl/Bio/welcome.html

package Bio::UnivAln;
use strict;

use vars qw ($VERSION $Revision);

$VERSION    = 1.010; # bio.perl.org Version;
$Revision   = '$Id$';

# Disclaimer from Georg Fuellen:
# UnivAln is now under the CVS system. Georg Fuellen is currently 
# not working on it, nor is he responsible for the version 
# distributed via bio.perl.org. This version is nearly identical 
# to Georg's latest version 1.009.

# Revision History: see the end of the POD, under the header `REVISION HISTORY'
# Last important change in Version 1.008 on 13 May 1998 :
# readseq is now used to parse files that have been processed
# as ``raw'' before; now ``raw'' format is recognized using the
# expression /^[A-Z_0-9$_GAP_SYMBOL$_UNKN_SYMBOL\s]+$/im,
# i.e. the file may only have alphanumerical characters,
# gap and unknown-symbol, and whitespace. If commata, etc, 
# are detected, readseq is used for parsing. Readseq itself 
# seems to be unable to detect ``raw'' format in some cases.

# Copyright (c) 1996, 1997, 1998 Georg Fuellen.
# Some pieces of the code were contributed by Steven E. Brenner, 
# Richard Resnick and Chris Dagdigian. Thanks !!!!
# This module is free software; you can redistribute it and/or modify 
# it under the same terms as Perl itself.

=head1 NAME 

Bio::UnivAln - Bioperl alignment object

=head1 SYNOPSIS

This Perl module is intended to simplify the handling of biosequence alignments.
When in doubt, always check our Homepage for the newest version, contact emails,
help files, etc: http://www.techfak.uni-bielefeld.de/bcd/Perl/Bio/welcome.html

See the L<REVISION HISTORY> for recent bugfixes and enhancements.

=head2 Object Creation in a Nutshell

  use Bio::UnivAln;

  my $aln = Bio::UnivAln->new('t/alnfile.fasta');
  $aln = Bio::UnivAln->new(-file=>'t/alnfile.aa',
                       -desc=>'Sample alignment',
                       -type=>'amino',
                       -ffmt=>'raw'      # 1 line in file -> 1 sequence
                      );
  $aln = Bio::UnivAln->new(-seqs=>"TCCCGCGTCAACTG\nTGGTGCTTCAACCG\nACTTG--TCAACTG");
  $aln = Bio::UnivAln->new(-seqs=>[$sequence_strg,\@character_list,$bioSeqObject]);
  $aln = Bio::UnivAln->new(-seqs=> ['ACCCGCGTCAACTG', 
           ['A','G','G','G','G','C','T','T','C','A','A','C','C','G'], 
           Bio::Seq->new(-seq=>'ACTTG--TCAACTG')
         ]);
  $aln = Bio::UnivAln->new($file,$seqs,$id,$desc,$names,$row_ids,$col_ids,
           $row_descs,$col_descs,$numbering,$type,$ffmt,$descffmt,$inplace);  

=head2 Object Manipulation in a Nutshell

  OUTPUT

  $aln->ffmt('fasta');   # set default output format
  print "\n aln in default format:\n", $aln->layout();
  print "\n aln in raw format:\n", $aln->layout("raw");
  print "\n aln in fasta format:\n", $aln->layout("fasta"), "\n";
  print "\n aln in MSF format (via readseq):\n", $aln->layout("MSF"), "\n";

  SLICING

  my $alnSlice = $aln->seqs(1,3,1,2); # multiline string of rows 1-3, columns 1-2
      print $alnSlice, "\n";
  my @alnSlice = $aln->seqs([1..3], [1,4]); # multidimensional array 
                                            # of rows 1-3, columns 1+4
      for $aref ( @alnSlice ) { print @$aref, "\n"; }
  $alnSlice = $aln->seqs([3,2,3,1], [1,3..5]); # rows 3,2,3,1, cols 1,3..5

  ADVANCED SLICING

      sub has_purine {
        my $str = join "", @{ $_[0] };
        if ($str =~ /[AaGgRr]+/) {return 1;} else {return 0;}
      }
  @alnSlice = $aln->seqs([1,2,3], \&has_purine); # rows 1-3, and from these
      # only the entries of columns for which has_purine returns 1
  $alnSlice = $aln->seqs({ids=>'SeqA SeqB'},{ids=>'ColA ColB ColC'});
  $aln->inplace(1); $aln->seqs({ids=>'A B'},[1..6]); $aln->inplace(0);
      # manipulates the object itself, assigning the slice internally

  MAPPING

  @resSlice = $aln->map_r(\&has_purine, [1..3]); # 1,0,1 if row 1+3 has purine
  @resSlice = $aln->map_c(\&has_purine, [1,4]); # 1,0 if column 1 has purine

  UTILITIES

  $resSlice = $aln->consensus(); # 75% majority needed for consensus letter
  $resSlice = $aln->consensus(0.6, [1,3]); # 60% majority, columns 1+3 only
  $resSlice = $aln->var_sites(); # no columns that are invariable ...
  $indices = $aln->var_inds();   # ... and their indices
  $resSlice = $aln->invar_sites(); # only columns that are invariable
  $indices = $aln->invar_inds();   # (..inds() is available for most utilities)
  $resSlice = $aln->var_sites(0.6); # no columns with >= 60% maj. of 1 letter
  $resSlice = $aln->invar_sites(0.6); # only columns with >= 60% majority
  $resSlice = $aln->gap_free_sites(); 
  $resSlice = $aln->no_allgap_sites(); # exclude sites that only have gaps
  $resSlice = $aln->reverse([1,3]); # reverse of rows 1+3 only
  $resSlice = $aln->complement([1,3]); # dna/rna complement of rows 1+3 only
  $resSlice = $aln->revcom([1,3]); # reverse complement, rows 1+3 only
  $resSlice = $aln->remove_gaps(); # original sequences without gaps

=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

The following are the installation instructions from the original UnivAlign
distribution for installing UnivAlign.pm by itself. These are not necessary if
installing from the central Bioperl distribution. Note that the central
distribution does not currently run the more extensive univaln.t2 
test harness.

The original installation package is available from:

 http://www.techfak.uni-bielefeld.de/bcd/Perl/Bio/#univaln 

To install, untar it and run the following commands in the module directory
created by tar:

        % perl Makefile.PL
        % make
        % make test 
        % make install

For intensive testing, see t/univaln.t2. Run that script via

        % perl t/univaln.t2 > t/my_univaln

and compare the output with the file t/univaln.o .
Expected error messages can be found in t/univaln2_expected_errors .

If you do not have superuser installation privileges, or to install 
in a different directory, you can do one of two things:

1) Either copy the UnivAln.pm module into the ``Bio'' subdirectory of
an accessible Perl module directory (e.g. /my/perl/lib/Bio/). 
(One possible ``accessible Perl module directory'' is your current working 
directory; there you can create a subdirectory named ``Bio'' and place 
UnivAln.pm in that subdirectory. Then, you can start scripts using
Bio::UnivAln from your current working directory.)

2) Or, run the above commands but specify an alternate location for 
the module by supplying a PREFIX argument to the first command:

        % perl Makefile.PL PREFIX=/my/perl/stuff

This will place the module into '/my/perl/stuff/lib/site_lib/'. To
specify a directory for the module file, set the $INSTALLSITELIB 
variable in Makefile.PL (e.g., $INSTALLSITELIB = '/my/perl/lib');

The make install command may report problems with creating documentation
files (pod2man/perllocal.pod); please find the documentation in 
UnivAln.pm.html instead.

If you have Bio::Seq installed and want to test Bio::Seq usage by
Bio::UnivAln, you need enable the line
``use Bio::Seq;'' at the beginning of t/univaln.t2;
Bio::Seq can be found via http://www.techfak.uni-bielefeld.de/bcd/Perl/Bio/ .
Note that the test script will also test error handling; you can expect the
error messages included in the file univaln.t2_expected_errors -- these are OK.

If you wish that the module uses Don Gilbert's readseq package for sequence
format conversion (Version 1 Feb 1993), you can set the environment variable 
`READSEQ_DIR'" appropriately. (Currently, only ``fasta'' and ``raw'' format 
are supported directly by UnivAln.pm.) 
Then, the program detects and uses `readseq' automatically, if it is in the
specified directory (the default directory is ``./''). Modifying the
environment variable `READSEQ' changes the expected name of the executable.
For example, $ENV{READSEQ_DIR} may be ``/vol/biotools/bin/'' and $ENV{READSEQ} 
``readseq2.0''. Readseq will give you support for PIR/CODATA, MSF/GCG and 
PAUP/NEXUS formats; ASN.1 does not seem to work reliably.
(URLs: http://iubio.bio.indiana.edu/IUBio-Software+Data/molbio/readseq/
http://dot.imgen.bcm.tmc.edu:9331/seq-util/Help/readseq.html
http://bimas.dcrt.nih.gov/molbio/readseq/formats.html )

Similar support for conversion from Clustal format, using Clustal as
a converter, is implemented, but not properly tested and documented.
The relevant environment variables are `CLUSTAL_DIR' and `CLUSTAL'.
(URL: http://www-igbmc.u-strasbg.fr/BioInfo/ClustalW/Top.html )

(Thanks to Steve A. Chervitz for his help with bundling the module !)

=head1 DESCRIPTION

This module is the Bio::UnivAln alignment object which is part of 
the Bioperl project. Currently it has some nice methods for accessing
an alignment after reading it in from certain formats, incl. utilities
like consensus and reverse complement. Bio::Seq (single sequences) 
is only needed if you explicitly want to use these.

(Most examples below are taken from the test script(s) that can
be found in directory ``t'' of the Bio::UnivAln distribution.
There you will also find a CGI script producing some graphics, 
which is currently in alpha status: I suspect it needs some
refitting to run on a different server. If you'd like to know more about 
multiple alignments, in theory and practice, check out the tutorial at
http://www.techfak.uni-bielefeld.de/bcd/Curric/MulAli/mulali.html )

=head2 CREATION OF ALIGNMENTS

Alignments can be constructed from files, (multi-line) strings,
arrays and Bio::Seq objects. Files need to be in a standard format,
as described below, under the header L<Alignment Formats>.

  my $aln = Bio::UnivAln->new('t/alnfile.fasta');

The first parameter is regarded as a file name; if you pass
additional parameters, they will overwrite the parameters read in
from the file. You can use named parameters; take a look at
the documentation on the new() method in the appendix for a list of all
parameters, and their names. In the following example, description,
sequence type, and file format are provided. The file format will
relieve Bio::UnivAln from guessing it; however, there are no guarantees if 
you bypass Bio::UnivAln's guessing _and_ provide an incorrect file format.

  $aln = Bio::UnivAln->new(-file=>'t/alnfile.aa',
                       -desc=>'Sample alignment',
                       -type=>'amino',
                       -ffmt=>'raw'      # 1 line in file -> 1 sequence
                      );

If no description (``-desc'') is given, a default one will be based on the 
file name.  The format type is also the default format for output; if both 
differ, you need to specify the input format when you construct the $aln 
object, and then use the accessor ffmt() to set the default output format.
Bio::UnivAln can be passed the aligned sequences directly, using the
named parameter ``-seqs''. It takes a multi-line string, or any mix of 
strings, array references, and Bio::Seq objects:

  $aln = Bio::UnivAln->new(-seqs=>"TCCCGCGTCAACTG\nTGGTGCTTCAACCG\nACTTG--TCAACTG");
  $aln = Bio::UnivAln->new(-seqs=>[$sequence_strg,\@character_list,$bioSeqObject]);
  $aln = Bio::UnivAln->new(-seqs=> ['ACCCGCGTCAACTG',
           ['A','G','G','G','G','C','T','T','C','A','A','C','C','G'],
           Bio::Seq->new(-seq=>'ACTTG--TCAACTG')
         ]);

=head2 ACCESS TO THE DATA, AND MANIPULATION

The layout() method returns the sequence in a specified format;
supported formats are listed under the header L<Alignment Formats>.

  $aln->ffmt('fasta');   # set default output format
  print "\n aln in default format:\n", $aln->layout();
  print "\n aln in raw format:\n", $aln->layout("raw");
  print "\n aln in fasta format:\n", $aln->layout("fasta"), "\n";
  print "\n aln in MSF format (via readseq):\n", $aln->layout("MSF"), "\n";

=head2 Access by Specifying Boundaries

You can calculate slices of alignments in a very flexible way;
interval slices like the intersection of rows 1-3 and columns 1-2 are 
calculated using seqs(). Here, intersection means that those
elements are returned that are both in rows 1-3 and in columns 1-2.

  $alnSlice = $aln->seqs(1,3,1,2);  # rows 1-3, columns 1-2
  $alnSlice = $aln->seqs();  # returns the whole alignment

Here's a diagram illustrating the general case, intersecting rows
$y_lo to $y_hi, and columns $x_lo to $x_hi.

                $x_lo      $x_hi
                 :           :
     .. $y_lo ...:...........:.................
                 :::::::::::::
                 :: SELECTED :
                 ::: PART ::::
                 :::::::::::::
     .. $y_hi ...:::::::::::::.................
                 :           :
                                            Fig.1

Maximal intervals will be assumed if no parameters are provided.
Per default, the first row (sequence) has index 1 (not 0), and the first
column has index 1 (not 0). The latter can be modified using numbering(). 

=head2 Access by Index Lists

If you desire non-consecutive row / column elements,
you can specify the indices as lists, one list of desired row indices 
and one list of desired column indices.

  $alnSlice = $aln->seqs([1..3], [1,4]);  # rows 1-3, columns 1+4

Here, letters in columns 2+3 will not be returned. Another example:

  $alnSlice = $aln->seqs([3,2,3,1], [1,3..5]);  # rows 3,2,3,1, cols 1,3..5

If you specify the empty list (``[]''), all rows/columns will be returned.
The following diagram shows the case where the list of row indices
is [$r1,$r2,$r3,$r4], and the list of column indices is [$c1,$c2,$c3].

                 :    :          :
                 :    :          :
     ............*....*..........*.............  $r1
                 :    :          :
     ............*....*..........*.............  $r2
     ............*....*..........*.............  $r3
                 :    :          :
     ............*....*..........*.............  $r4
                 :    :          :
                 :    :          :
                $c1  $c2        $c3           Fig.2

Again, an element is selected for the slice if and only if it lies
in the intersection of a row and a column which are both desired 
according to the index lists.

=head2 Return Values

In the examples above, a string (scalar) is returned; the standard
sequence accessor seqs() always returns a (multi-line)
string in a scalar context. In a list context, it returns an array;
each element of such an array is a reference to another array holding 
the letters of one sequence, i.e. one single row.

  @alnSlice = $aln->seqs([1..3], [1,4]);  # rows 1-3, columns 1+4
      for $aref ( @alnSlice ) { print @$aref, "\n"; }

If you use the result of an accessor or a utility function in the ``-seqs''
slot of new(), you may need to force that result into a scalar 
context, because the accessor, etc, returns a list in a list context,
and the constructor naturally provides such a list context since it expects a
list of parameters. 

  $aln = new Bio::UnivAln(-seqs=>scalar($aln2->seqs()));

In the example above, the list-context return value of
C<$aln2->seqs()>, i.e. the list of rows of $aln2,
would be fed one by one as additional parameters into the constructor,
if you didn't ``protect'' it by scalar(). You will be warned about the problem
because Bio::UnivAln detects any named parameters that it can't use.

=head2 Access by Id

(The following access method is currently in alpha status, it may need
some revision until the code is fully released.)

Any list of desired row/column indices can be replaced by a hash of desired 
ids, which are recognized if they are in the object's own list of row or column 
ids. You need to pass a reference to a hash that has one key, ``ids'',
and one value, which is a string containing the ids seperated by `` ''(blank) :

  $alnSlice = $aln->seqs({ids=>'SeqA SeqB'},{ids=>'ColA ColB ColC'});

Row (sequence) ids are automatically extracted when reading fasta files and 
Bio::Seq objects. Otherwise, they are set to the default numerical index list,
like (1..20) if the alignment has 20 rows. Since there's currently no way
to extract column (site) ids (none of the supported formats has this feature),
these always hold the default numerical indices. However, both lists
may be set using the accessors row_ids() and col_ids(). Note that arbitrary 
numbering schemes can be supported this way.

=head2 Access by Selector Function

Finally, you can specify a function such that seqs() returns exactly 
those letters which lie in a row (or column) for which your function returns 
true. E.g. has_purine() returns true if the row/column contains A, a, T, t, R, 
or r; the following $alnSlice will contain only columns that have
one of these letters in them.

      sub has_purine {
        my $str = join "", @{ $_[0] };
        if ($str =~ /[AaGgRr]+/) {return 1;} else {return 0;}
      }
  $alnSlice = $aln->seqs([1..3], \&has_purine); # rows 1-3, and from these
      # only the columns for which has_purine returns 1

Similarly, the list of row indices, [1..3], could be replaced by a
function, which is then used to designate the desired rows. (It has the
same role that the expression EXPR has in the Perl code template
B<grep EXPR, LIST>, see the perlfunc manpage. Bio::UnivAln also provides
the equivalent to B<map EXPR, LIST>; this ``Mapping'' will be discussed soon.)

In other words, any list of indices can be replaced by the reference to a 
function r or c. The function then selects those elements (*) that lie in 
a row or column which meets a criterion, returning true if the criterion 
is met, and false otherwise.
If both lists of indices are replaced by functions, the picture is like this:

                 :    :          :
                 :    :          :
     ............*....*..........*.............  r(row) = true
                 :    :          :
     ............*....*..........*.............  r(row) = true
     ............*....*..........*.............  r(row) = true
                 :    :          :
     ............*....*..........*.............  r(row) = true
                 :    :          :
                 :    :          :
           c(column) c(column)  c(column)
            = true    = true     = true       Fig.3

The function out_fasta observes all types of selectors that seqs() does, 
and returns the result in fasta format :

  print "\n aln in fasta format:\n", $aln->out_fasta([1..3],\&has_purine), "\n";

=head2 Mapping functions onto Sequences and Columns

You can map a function onto selected rows or columns, and receive the
results as a list:

  @resSlice = $aln->map_r(\&has_purine, [1..3]); # 1,0,1 if row 1+3 has purine
  @resSlice = $aln->map_c(\&has_purine, [1,4]); # 1,0 if column 1 has purine

As you may expect, maximal index lists (i.e, all rows / all columns) will be 
assumed if no second parameter is provided. In the same way as before,
any one if the index lists may be replaced by a hash of ids, or a selector
function. In the section on L<User-defined Utility Functions>, map_r() and 
map_c() are used to implement user-defined consensus, reverse, complement, etc. 

By the way, the following map is the same as @alnSlice = $aln->seqs([1..3]) 
since mapping the identity function ``sub id { return @_ }''
to desired row/column subsets and collecting the result is just like slicing.

  @alnSlice = $aln->map_r(\&id, [ 1..3 ]);

The number of rows/columns of an alignment can be obtained by using
height() and width() respectively. Ids and descriptions of rows (sequences)
and columns (sites) can be manipulated using row_ids(), col_ids(),
row_descs() and col_descs(). Be B<warned> that these accessors just return
a reference to the array of ids/descriptions; if you'd like to process
the array without changing the object's data, you need to create your
own copy. This decision was taken because especially the result of col_ids()
can be huge, and for a lot of applications a deep copy is unnecessary.

=head2 Inplace Manipulation

Slices can be applied to the object itself, replacing the old alignment
by a new one that is sliced from the old ('inplace' manipulation).  
This is particularly useful is the alignment is large. The row (sequence)
and column (site) ids are taken over from the old alignment. They are
used as the lookup tables for L<Access by Id>, and available via row_ids() 
and col_ids().

The following code sets the ``inplace'' flag, and overwrites the current
alignment with the rows named A and B, and columns 1-6:

  $aln->inplace(1);
  $aln->seqs({ids=>'A B'},[1..6]);
  $aln->inplace(0);

'inplace' manipulation is also available for most utility functions below, 
with the notable exception of consensus(). A full list is given in the 
description of inplace(), see the Appendix.

=head2 UTILITY FUNCTIONS LIKE CONSENSUS, AND REVERSE COMPLEMENT

In the following paragraph, Bio::UnivAln's direct support for utility
functions like consensus, (in)variable sites, gap-free sites, reverse, 
complement, and reverse complement is explained.

  $resSlice = $aln->consensus();

calculates the consensus of the columns, i.e. those columns for which there
exists a letter which has an absolute majority of 75% or more. The letter
'!' designates the case that no consensus is given.
To bypass the default threshold, write

  $resSlice = $aln->consensus(0.6);

Note that for values smaller or equal to 0.5, two letters may have an equal
(relative) majority, and the tie is currently broken arbitrarily.
Thresholds larger than or equal to 1 imply that the site that has a 
consensus letter must be invariant. 

  $resSlice = $aln->consensus(1, [1..10]);

Just like seqs() and map_c()/map_r(), consensus() can be passed a reference 
to a hash of desired column ids, a function that selects the columns,
or the desired column indices themselves. (Here, it's columns 1 to 10.
Internally, map_c() is used to implement consensus().)

The methods var_sites() and invar_sites() are using consensus(), 
just checking whether there is a consensus letter ('invariable'), 
or not ('variable'). Naturally, the default threshold is 1, i.e. 
sites must be truly invariable (100% majority of 1 letter), 
or truly variable (strictly less than 100% majority). Per default,
var_sites() and invar_sites() return a B<multiline string of rows> (sequences),
with elements from the desired columns (sites) only. They do B<not>
return the alignment in a column-by-column fashion ! Also, don't be confused
by the possibilty to specify a list of desired row indices / ids, or
a selector function: It just constrains the output further, without
influencing which columns are selected. (Internally, seqs() is used to 
implement these functions; it is passed a row selector, and a function
that returns true if there is a consensus letter.)

  $resSlice = $aln->var_sites();
  $resSlice = $aln->invar_sites();
  $resSlice = $aln->var_sites(0.6, [1,3]); # no columns with >= 60% majority
                                     # of one letter; also, print rows 1+3 only
  $resSlice = $aln->invar_sites(0.6, [1,3]);  # only columns with >= 60% maj.
  $resSlice = $aln->no_allgap_sites([1,3]); # exclude sites that only have gaps
  $resSlice = $aln->gap_free_sites([1,3]); # exclude sites that have >= 1 gaps

In a similar fashion, the last example selects gap-free columns, and from 
those only prints the elements that happen to be in rows 1+3.
All these utilities support L<Inplace Manipulation>:

  $aln->inplace(1);
  $aln->var_sites(0.6, {ids=>'1 2'});
  $aln->inplace(0);

The utilities reverse(), complement(), and revcom() also allow
for L<Inplace Manipulation>, and again you can
specify the rows which shall form the new alignment by overwriting the 
old, or be returned as reverse'd, complement'ed, or reverse complement'ed :

  $resSlice = $aln->reverse([1,3]);
  $resSlice = $aln->complement([1,3]);
  $resSlice = $aln->revcom([1,3]);

If an inplace manipulation reverses column order, e.g. in the case of 
reverse(), this will be reflected in the column ids available via col_ids().
Note that complement is defined according to the IUPAC code,
using the same substitutions that Bio::Seq uses, such that results
obtained for e.g. Amino Acid sequences are probably nonsensical.

For all functions that have ``sites'' in their name, a corresponding
``inds'' function is available that returns the relevant array of indices 
instead of the sites themselves:

  $indices = $aln->var_inds(); # array of indices of the variable sites
  $indices = $aln->no_allgap_inds([1,3]);

Finally, the original sequences (without gaps) are available via remove_gaps().

  $some_original_sequences = $aln->remove_gaps([1,3]);

In a list context, all these utility functions return an array of array 
references, just like seqs(); the exception is consensus(), which then returns 
a simple array of consensus letters.

=head2 ALIGNMENT FORMATS

The directly supported formats are fasta and raw (1 line in file -> 1 sequence,
where "\n" (newline) is the delimiter.) If available (see L<INSTALLATION>),
readseq is used to parse and write PIR/CODATA, MSF/GCG and PAUP/NEXUS 
formats; ASN.1 does not seem to work reliably. Clustal is used to parse
in clustal format, if available.

=head2 ADVANCED STUFF

An exhaustive list of accessors and methods is given in the appendix.
E.g., Alignments may be copied using copy(), and compared for equality modulo
gaps with equal_nogaps().

Bio::UnivAln doesn't really care whether the sequences passed in are indeed
from an alignment (i.e. have the same length); sequence bags (i.e. multisets
of sequences where the same element may occur more than once) are 
therefore handled, too. However, you need to be careful with some
methods (e.g. the accessor seqs(), and width()) because their default
behavior may depend on inspecting the first sequence only (or the first
requested sequence, if this information is available -- new feature in 1.006) :
If you use seqs() on a sequence bag, and don't provide the number of columns
explicitly, you may be surprised to find out that the length of
the first sequence is used as the default length, and not the length 
of the longest sequence. width() uses the same heuristic.

You can add gap symbols to the elements of a multiset so that
each element takes the length of the longest sequence by calling
equalize_length(); you can even pad more gap characters by specifying
the new width yourself, as an argument to equalize_length().

If you need to have the alignment in an intermediate form, i.e.
neither an array of array references, nor a single multiline string, but
an array of strings, just request the output as a multiline string, and split
it on "\n". For example, to just extract a single column, you can do

  @col = split "\n", $aln->seqs([],[$colindex]);

Checking an alignment for characters that should not be there according 
to the Alignment Type is currently not well-supported; see alphabet_check()
for a _preliminary_ way of doing manual checkup, though.

The hash referenced by $names stores {loc,name} 
pairs of other database locations and corresponding names where 
the alignment is located. 
Currently, loc and name must both be set as text, and must consist
entirely of a string matching the regexp /^[A-Za-z]\w*$/.  That is, they must 
be a single "word" of at least one character, without a leading underscore.
This restriction is not enforced, but code which deviates is subject to
break in future releases.  Note also that the object may place any other
sort of items in the name string, so users of this hash should not rely on
accessing entries conforming the requirements above.  

=head2 Alignment Types

The supported sequence types and corresponding alphabets are the same
as reported in the documentation of the Bio::Seq object: ``dna'', ``rna'',
and ``amino''. They are carried along without much checking, etc.

=head2 User-defined Utility Functions

Here are some user-defined auxiliary functions which can be applied to 
rows and columns of an alignment; this technique can be used if 
no appropriate built-in functions are available, or if you don't want 
to use them.

maj_dominance returns false if all letters in a column/row are the same.

      sub maj_dominance {
        my $str = join "", @{ $_[0] };
        my $first_res = @{ $_[0] }[0];
        $str =~ s/$first_res//g;
        if ($str) {return 0;} else {return 1;}
      }

consensus returns the letter that has the (relative) majority among all
residues, if it exceeds a certain threshold. (For simplicity, the threshold
is hard-wired.)

      sub consensus {
        my @chars = @{ $_[0] };
        my %temp = ();
        my $threshold = 0.34;  # more than 1/3
        my @list = sort { $temp{$a}<=>$temp{$b} } 
                        grep ++$temp{$_} >= $threshold * ($#chars+1), @chars;
          # In case of a tie, it's not specified which residue is in $list[-1]
        return (defined($list[-1]) ? $list[-1] : '!');
      }

reverse_ and complement implement reversing a sequence (note the 
underscore; I'm looking for better ways to implement function passing),
and complementing it according to IUPAC conventions.

      sub reverse_ {
        return [ reverse @{ $_[0] } ];
      }
      sub complement {
        my $str = join "", @{ $_[0] };
        $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
        return [ split "", $str ];
      }

Here, the functions above are applied:

  @resSlice = $aln->map_c(\&maj_dominance);
  print "\nDominated sites:\n", @resSlice, "\n";

  @resSlice = $aln->map_c(\&consensus);
  @resSlice = $aln->map_c(\&consensus, \&has_purine);
  print "\nConsensus of the columns that have purine\n", @resSlice, "\n";

  @resSlice = $aln->map_r(\&reverse_);

  @resSlice = $aln->map_r(\&complement);

You can also parametrize your functions by using ``closures'' (See
``Programming Perl'', 2nd Ed., p.253). Basically, you set up a function
that has the parameter built into it, and pass around the reference to
that function. Here's a function that does the set-up:

      sub _setup_consensus_with_threshold {
        my $threshold = shift;
        return sub {
          my @chars = @{ $_[0] };
          my %temp = ();
          my @list = sort { $temp{$a}<=>$temp{$b} }
                          grep ++$temp{$_} >= $threshold * ($#chars+1), @chars;
            # In case of a tie, it's not specified which residue is in $list[-1]
          return (defined($list[-1]) ? $list[-1] : 'N');
        }
      }

Here, you create an instance of the function, incl. the parameter.
$consensus holds a reference to the new instance.

  my $consensus = _setup_consensus_with_threshold(0.75);

and finally, you pass the reference to the function around, e.g. to map_c().

  @resSlice = $aln->map_c(\&$consensus);

(You can do pretty cute tricks using closures, e.g. you can set a counter
to 1 in the setup function, and increment it in the ``real'' function.)

=head1 Bio::UnivAln Guts

Currently, the object hash has the following keys. This may be subject
to change; in particular the alignment data may at some point be stored
more efficiently in a PDL (Perl Data Language) array.

  $self->{'seqs'}   : An array of array references, each of which holds
                    one sequence of the alignment
  $self->{'id'}     : String specifying the ID; shall be in \w+ (i.e. composed
                    of characters in [a-zA-Z_0-9]; only \S+ is enforced, though)
  $self->{'desc'}   : String giving a description, (later) to be formatted
                    according to $descffmt
  $self->{'names'}  : A reference to a hash which stores {loc,name} pairs of
                    other database locations and corresponding names where
                    the alignment is located.
  $self->{'row_ids'}: A reference to an array which stores row (sequence) ids
  $self->{'col_ids'}: Same as $self->{'row_ids'}, for the columns (sites)
  $self->{'row_descs'}: A reference to an array which stores row (sequence)
                    descriptions
  $self->{'col_descs'}: Same as $self->{'row_descs'}, for the columns (sites)
  $self->{'numbering'}: The offset of the first column, an integer
  $self->{'type'}   : The type of the alignment, concatenated from the molecule 
                    type (see L<Alignment Types>) and a flag that is not 
                    currently used, but intended to flag sequence bags
  $self->{'ffmt'}   : alignment format, see L<Alignment Formats>.
  $self->{'descffmt'}: format of $desc; right now this should be ``raw''
                    or ``fasta'' which just implies that no specific
                    format is being followed, any text is allowed
                    excluding ``\n''(newline). More support is planned.
  $self->{'inplace'}: Flag which is set to true if accessors (and utility
                    functions) should make the modification to the object
                    itself, and just return true on success. (See inplace().)


Some more helpful comments...
In Perl, false is 0 or "", true is everything else.
Not all internal functions have a POD/HTML documentation.
If you read the code of this module, you need to be familiar with
map and grep...

=head1 TO-DO

Soon: 

Better handling of access, copying and slicing of C<$self->{'row_ids'}>, etc.
Add pointer to a UnivAln demo page based on the draft at
http://www.techfak.uni-bielefeld.de/bcd/Perl/Bio/Docs/phylosnapshot.html.
Fix alphabet_check() (current workaround cannot be generalized),
Fix UnivAlnAlph{UnivAlnType{"unknown"}} (current setting is not so nice;
what's the easiest way to set it to all word characters ?),
Support more formats, especially Nexus,
Improved sequence alphabet support; not just providing a function 
alphabet_check() for manual checking.
Then, maybe use strings to represent alphabets.
Test and document aln() for returning a slice as an alignment.

Later: 

Validity marker for correctly initialized/manipulated objects.
Assign _unique_ IDs if none are provided.
Functions like has_seqs(), etc.
Using `undef' during initialization, and functions for regaining such a state.
Use Perl Data Language ?!

=head1 DISCLAIMER

How is this for a maximum of disclaiming warranty ? In short, I'm developing 
this module in my spare time and it's for free, don't sue me :-)

 IN NO EVENT SHALL THE GLOBEWIDE NETWORK ACADEMY, THE VIRTUAL SCHOOL OF NATURAL
 SCIENCES, THE AUTHOR OR THE UNIVERSITY OF BIELEFELD BE LIABLE TO ANY PARTY
 FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING
 OUT OF THE USE OF THIS CODE, EVEN IF THE GLOBEWIDE NETWORK ACADEMY,
 THE VIRTUAL SCHOOL OF NATURAL SCIENCES, THE AUTHOR OR THE UNIVERSITY OF 
 BIELEFELD HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 THE GLOBEWIDE NETWORK ACADEMY, THE VIRTUAL SCHOOL OF NATURAL SCIENCES, THE
 AUTHOR AND THE UNIVERSITY OF BIELEFELD SPECIFICALLY DISCLAIM ANY WARRANTIES,
 INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 AND FITNESS FOR A PARTICULAR PURPOSE.  THE CODE PROVIDED HEREUNDER IS ON AN
 "AS IS" BASIS, AND THERE IS NO OBLIGATION WHATSOEVER TO PROVIDE MAINTENANCEE,
 SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

=head1 COPYRIGHT

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

Copyright (c) 1996, 1997, 1998 Georg Fuellen. All Rights Reserved.
Some pieces of the code were contributed by Steven E. Brenner, 
Richard Resnick and Chris Dagdigian. Thanks !!!!


=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other Bioperl modules.
Send your comments and suggestions preferably to one of the Bioperl mailing lists.
Your participation is much appreciated.

    bioperl-l@bioperl.org          - General discussion
    http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the bugs and 
their resolution. Bug reports can be submitted via email or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR

Georg Fuellen

Technische Fakultaet - AG Praktische Informatik,
Universitaet Bielefeld,
D-33501 Bielefeld,
Germany,
georg.fuellen@uni-bielefeld.de

http://www.techfak.uni-bielefeld.de/~fuellen/

=head1 ACKNOWLEDGEMENTS

Steven E. Brenner, Steve A. Chervitz, Michael Constant, Richard Resnick, 
Chris Dagdigian, Lew Gramer [more to follow]

=head1 SEE ALSO

 Bio::Seq.pm - The biosequence object

 http://bio.perl.org/Projects/modules.html  - Online module documentation
 http://bio.perl.org/Projects/SeqAlign/     - Bioperl sequence alignment project
 http://bio.perl.org/                       - Bioperl Project Homepage

=head1 REFERENCES

If you'd like to acknowledge use of Bio::UnivAln in your work, please cite 

Fuellen, G. (1997). Bio::UnivAln - bioperl alignment object [WWW-Document].
URL http://www.techfak.uni-bielefeld.de/bcd/Perl/Bio/welcome.html

And please drop me a note :-)
An article for the Perl Journal is planned. The page is mirrored at

http://merlin.mbcr.bcm.tmc.edu:8001/bcd/Perl/Bio/welcome.html
http://www.biotech.ist.unige.it/bcd/Perl/Bio/welcome.html

=head1 REVISION HISTORY

Version 1.000 on 12 Feb 1997.  

Version 1.001 on 19 Feb 1997. Fixed a bug that 
triggered _rowbounds() and _colbounds() to use maximal index lists
whenever the first index in an index list was 0. New example in POD:

    $aln = new Bio::UnivAln(-seqs=>scalar($aln2->var_sites()));

Internal: Now avoiding any global parameter passing by using closures,
for the utility functions.

Version 1.002 on 21 Feb 1997. Renamed the module to UnivAln, see the
discussion (Feb 1997) in the vsns-bcd-perl mailing list archive.
Fixed hopefully all problems in out_graph(), all of them triggered by 
bugs/features of PGPLOT.

Version 1.003 on 25 Feb 1997. Added POD on using closures. 
Moved out_graph() into the cgi-script t/univaln.cgi (this is alpha-code !)
so that ``use PGPLOT'' is no longer required.

Version 1.004 on 13 Mar 1997. Fixed bug with reading fasta files.
Changed file format identifiers ``Fasta'' and ``Raw'' to ``fasta'' and ``raw''.
Added support for arbitrary numeration schemes by providing access by name/id 
for rows (sequences) and columns (sites).
Made _arys() and _strs() accept same row/column designations like seqs().
Added 'inplace' manipulation for slicing and most utilities.
Improved speed a little. Deriving default description string from file name, 
if available.

Version 1.005 on 18 Jun 1997. Using simple filename (``basename'') as default 
decription. Changed type identifiers ``Unknown'', ``Dna'', ``Rna'', ``Amino'' 
to ``unknown'', ``dna'', ``rna'', ``amino''. The old acronyms should still 
be supported for the forseeable future.

Version 1.006 on 15 Dec 1997. Added function no_allgap_sites(); out_fasta()
now accepts same row/column designations like seqs(). For advanced users:
the width() of sequence bags can now be controlled by supplying the index of
the row whose width is taken -- by default the width of the first row is used.

Version 1.007 on 15 Mar 1998.                                                 
* Added functions that return indices instead of actual sequence data:
var_inds(), invar_inds(), gap_free_inds(), unknown_free_inds(),
special_free_inds() and no_allgap_inds() return the indices of
(in)variable, gap-free, unknown-free, gap+unknown-free and non-gap-only 
columns; they are the counterparts of var_sites(), etc.
* Added remove_gaps() to return ungapped (original) sequence data
* Added equal_no_gaps() to check for equality ignoring gaps
* Added equalize_lengths() for padding a sequence bag with gaps
such that all rows have equal length (current procedure is slow)
* The default id is now ``_'' (was: ``No_Id_Given'')
* Internal functions _rowbounds() and _colbounds() now return indices 
from the user's perspective, but WITHOUT substracting offset.  Instead, 
all functions subtract offset after applying _rowbounds()/_colbounds()

Version 1.008 on 13 May 1998. 
* Added readseq conversion support, making it possible to read and write
the following formats: MSF, Paup, PIR/Codata, ASN.1. (ASN.1 was not
parsed successfully: readseq seems to be unable to read in its own ASN.1 
output). Technically, readseq is now used to parse files that have been 
processed as ``raw'' before; now ``raw'' format is recognized using the
expression /^[A-Z_0-9$_GAP_SYMBOL$_UNKN_SYMBOL\s]+$/im, i.e. the file 
may only have alphanumerical characters, gap and unknown-symbol, and 
whitespace. If commata, etc, are detected, readseq is used for parsing. 
Readseq itself seems to be unable to detect ``raw'' format in some cases,
causing weird results.
* Added Clustal support for parsing; still to be tested and documented.

Version 1.009 on 25 May 1998.
* The module can now be ``built'' in the standard way (perl Makefile.PL, etc).

=head1 APPENDIX

Please note that functions starting with an underscore (``_'') are
intended for internal use only: Use them only if you know what you're 
doing :-)

=cut

use Exporter;
use vars qw( @ISA  @EXPORT @EXPORT_OK );
@ISA         = qw(Exporter);
@EXPORT      = qw();
@EXPORT_OK   = qw($VERSION %UnivAlnType @UnivAlnType %UnivAlnForm @UnivAlnForm %UnivAlnAlphs @UnivAlnAlphs);
require 5.002;
use Carp;
#if ($] < 5.005) {carp "Not tested for Perl 5.002 or 5.003";}
use POSIX;
use File::Basename;

use vars qw($_NOFILE_FLAG $_GAP_SYMBOL $_UNKN_SYMBOL $_UNKN_SYMBOL2 $_NO_CONSENSUS_SYMBOL $_NaN $_MAX_SIZE_TAXANAME $_READSEQ $_CLUSTAL);
$_NOFILE_FLAG='_undef';
$_GAP_SYMBOL ='-';
$_UNKN_SYMBOL='\?'; # the backslash seems to be needed to be able to write
                    # ($str =~ /$_UNKN_SYMBOL/) ??!!
$_UNKN_SYMBOL2='N'; # ONLY USED IN UNDOCUMENTED FUNCTION special_free_seqs()
$_NO_CONSENSUS_SYMBOL='!';
$_NaN = POSIX::INT_MAX;   #anything better ?
$_MAX_SIZE_TAXANAME=20;
$_READSEQ = "";
$_CLUSTAL = "";
$ENV{READSEQ} = "readseq"; 
 # name of the readseq executable to be found in 
 # $ENV{READSEQ_DIR} or "\.". Set this to "" if readseq is
 # installed, but shall be ignored
$ENV{CLUSTAL} = "clustal";                                         
 # name of the clustal executable to be found in
 # $ENV{CLUSTAL_DIR} or "\.". Set this to "" if clustal is
 # installed, but shall be ignored

use vars qw(%UnivAlnType %TypeUnivAln @TypeUnivAln %UnivAlnForm %FormUnivAln @FormUnivAln);

# List of recognized alignment types
$UnivAlnType{unknown}  = 0;
$UnivAlnType{dna}      = 1;
$UnivAlnType{rna}      = 2;
$UnivAlnType{amino}    = 3;
$UnivAlnType{otherseq} = 4;  # type is explicitly NOT dna/rna/amino/etc, but 
                             # nothing more is known
     # Future support for the following 4 acronyms is doubtful:
$UnivAlnType{Unknown}  = 0;
$UnivAlnType{Dna}      = 1;
$UnivAlnType{Rna}      = 2;
$UnivAlnType{Amino}    = 3;
# %TypeUnivAln = reverse %UnivAlnType; # TypeUnivAln is inverse of %UnivAlnType
                             # i.e. keys and values are exchanged
@TypeUnivAln = ('unknown','dna','rna','amino','otherseq');

# List of recognized file formats 
$UnivAlnForm{unknown}  = 0;
$UnivAlnForm{raw}      = 13;
$UnivAlnForm{raw2}     = 14;
$UnivAlnForm{fasta}    = 7;
     # Future support for the following 3 acronyms is doubtful:
$UnivAlnForm{Unknown}  = 0;
$UnivAlnForm{Raw}      = 13;
$UnivAlnForm{Raw2}     = 14;
$UnivAlnForm{Fasta}    = 7;
# %FormUnivAln = reverse %UnivAlnForm; # FormUnivAln is inverse of %UnivAlnForm
                             # i.e. keys and values are exchanged
@FormUnivAln = ('unknown',undef,undef,undef,undef,undef,undef,
                'raw','raw2',undef,undef,undef,undef,'fasta');

use vars qw(%FuncParse %FuncOut %UnivAlnAlphs);

# %FuncParse and %FuncOut are for internal use ONLY
# %FuncParse is an array of (<ffmt>,<_parse_meth>), where <ffmt>
# is the file format code (e.g., 0 for "unknown", 7 for 'fasta', etc.), and
# <_parse_meth> is the method which parses strings in that file
# format. %FuncParse{$i} is mostly set to \&_parse_bad right now,
# since we don't have many file formats supported yet.

my %FuncParse =
  ($UnivAlnForm{"unknown"} => \&_parse_unknown,
   $UnivAlnForm{'fasta'}   => \&_parse_fasta,
   $UnivAlnForm{'raw'}     => \&_parse_raw,
   );
# Set all other formats to call &_parse_bad, so that currently unsupported 
# formats don't crash the software. 
grep {$FuncParse{$_} ||= \&_parse_bad} values %UnivAlnForm;

# array of implemented outputting routines, built in analogy to "%FuncParse"
my %FuncOut =
  ($UnivAlnForm{"unknown"} => \&out_fasta,  #fasta is the default output !
   $UnivAlnForm{'fasta'}   => \&out_fasta,
   $UnivAlnForm{'raw'}     => \&out_raw,
   $UnivAlnForm{'raw2'}    => \&out_raw2,
   );
grep {$FuncOut{$_} ||= \&out_bad} values %UnivAlnForm;

%UnivAlnAlphs =
 ($UnivAlnType{unknown}  => [ 'A','C','G','T' ], #better alternative ???
  $UnivAlnType{dna}      => [ 'A','C','G','T',$_GAP_SYMBOL ],
  $UnivAlnType{rna}      => [ 'A','C','G','U',$_GAP_SYMBOL ],
  $UnivAlnType{amino}    => [ 'A','R','N','D','C','Q','E','G','H','I','L','K','M','F',
                          'P','S','T','W','Y','V',$_GAP_SYMBOL ],
  $UnivAlnType{otherseq} => [ ], #better alternative ???
   );
#create alphabets like ``1Mg''= $UnivAlnType{dna}Mg 
#                             = [ 'A','C','G','T',$_GAP_SYMBOL,$_UNKN_SYMBOL ],
#where $_GAP_SYMBOL denotes the character for missing nucleotide data.
grep {$UnivAlnAlphs{$_.'Mg'} ||= [ @{ $UnivAlnAlphs{$_} },$_UNKN_SYMBOL ] } keys %UnivAlnAlphs;

# the following PODs should always have a blank at the beginning (i.e. 
# `` Usage'', not ``Usage'' !!) so that pod2html works fine !
=head2 new()

 Usage    :  $myAln = Bio::UnivAln->new($file,$seqs,$id,$desc,$names,
                        $row_ids,$col_ids,$row_descs,$col_descs,$numbering,$type,
                        $ffmt,$descffmt,$inplace);
                           - or -
             $myAln = Bio::UnivAln->new(-file=$file,
                                  -seqs=>$seqs,
                                  -id=>$id,
                                  -desc=>$desc,
                                  -names=>$names,
                                  -row_ids=>$row_ids,
                                  -col_ids=>$col_ids,
                                  -row_descs=>$row_descs,
                                  -col_descs=>$col_descs,
                                  -numbering=>$numbering,
                                  -type=>$type,
                                  -ffmt=>$ffmt,
                                  -descffmt=>$descffmt,
                                  -inplace=>$inplace);
 Function : The constructor for this class, returns a new object.
 Returns  : Bio::UnivAln object
 Argument : $file: file from which the alignment data can be read; all
                   the other arguments will overwrite the data read in. 
                   (see ``Alignment Formats'')
            $seqs: EITHER a reference to a list of either Bio::Seq objects, or
                   arrays of letters, or strings, or any mix of these,
                   OR a single (multi-line) string
            $id: String specifying the ID.
            $desc: String giving a description, (later) to be formatted 
                   according to $descffmt
            $names:A reference to a hash which stores {loc,name} pairs of 
                   other database locations and corresponding names where 
                   the alignment is located. (See L<ADVANCED STUFF>).
            $row_ids  :A reference to an array which stores row (sequence) ids.
            $row_descs:A reference to an array which stores row (sequence)
                    descriptions
            $col_ids:  Same as $self->{'row_ids'}, for the columns (sites)
            $col_descs:Same as $self->{'row_descs'}, for the columns (sites)
            $numbering: The offset of the first column
            $type: The type of the alignment, see ``Alignment Types''.
            $ffmt: alignment format, see ``Alignment Formats''.
            $descffmt: format of $desc; right now this should be ``raw''
                       or ``fasta'' which just implies that no specific
                       format is being followed, any text is allowed
                       excluding ``\n''(newline).
            $inplace: Flag which is set to true if accessors and utility
                      functions should make the modification to the object
                      itself, and just return true on success. See inplace().

=cut

sub new {
  my($this) = shift;
  my($class,$self);
  $class = ref($this) || $this;  # See the ``Perl Module List'' general section
  $self = {};
  bless $self, $class;
  $self->_initialize(@_);
  return $self;
}


=head2 _initialize()

 Usage    : n/a (internal function)
 Function : Assigns initial parameters to a blessed object.
 Returns  : 1 on success
 Argument : As Bio::UnivAln->new, allows for named or listed parameters.
            See ->new for the legal types of these values.

=cut

sub _initialize {
  my($self,@p) = @_;

  my $readseq = $ENV{READSEQ_DIR} || "./";
  $readseq .= ($readseq =~ /\/$/) ? "$ENV{READSEQ}" : "/$ENV{READSEQ}";
  if ((-x $readseq) && !(-d $readseq)) {
    $_READSEQ = $readseq;
    # carp "Using Don Gilbert's readseq program to convert sequence formats";
  } else {
    # carp "readseq not defined-- if installed, you can use its conversion 
    #     abilities by either providing an executable `readseq' in the script's
    #     homedirectory, or by setting the environment variable `READSEQ_DIR'";
  }
    
  my $clustal = $ENV{CLUSTAL_DIR} || "./";
  $clustal .= ($clustal =~ /\/$/) ? "$ENV{CLUSTAL}" : "/$ENV{CLUSTAL}";
  if ((-x $clustal) && !(-d $clustal)) {
    $_CLUSTAL = $clustal;
    # carp "Using clustal to convert clustal sequence formats";
  } else {
    # carp "clustal not defined-- if installed, you can use its conversion 
    #     abilities by either providing an executable `clustal' in the script's
    #     homedirectory, or by setting the environment variable `CLUSTAL_DIR'";
  }
    
    # Set default values that need not be inferred later
  $self->{'seqs'} = ( );
  $self->{'id'} = "_";
  $self->{'desc'} = "No Description Given"; # only preliminary!
  $self->{'names'} = { };
  $self->{'numbering'} = 1;
  $self->{'type'} = [$UnivAlnType{"unknown"},"unknown"];
  $self->{'ffmt'} = $UnivAlnForm{"unknown"};
  $self->{'descffmt'} = $UnivAlnForm{"unknown"};
  $self->{'inplace'} = 0;

  if ( ! @p ) {
    return undef;
  }
    # some operations will not work on ``empty'' placeholder objects

  my($file,$seqs,$id,$desc,$names,$row_ids,$col_ids,$row_descs,$col_descs,$numbering,$type,$ffmt,$descffmt,$inplace) =
    $self->_rearrange(['FILE',
                       'SEQS',
                       'ID',
                       'DESC',
                       'NAMES',
                       'ROW_IDS',
                       'COL_IDS',
                       'ROW_DESCS',
                       'COL_DESCS',
                       'NUMBERING',
                       'TYPE',
                       'FFMT',
                       'DESCFFMT',
                       'INPLACE'],
                       @p);

  $self->{'desc'} = defined($file) ? basename($file) 
                  : "No Description Given"; #cf File::Basename

    # Overwrite with values from file
  if ((defined($file)) && ($file ne $_NOFILE_FLAG) && ($file ne '')) {
    $self->_file_read($file,$ffmt);
  }

    # Overwrite with values from @_
  $self->_seqs($seqs); 
  $self->id($id);
  $self->desc($desc);
  $self->names($names);
  $self->row_ids($row_ids);
  $self->col_ids($col_ids);
  $self->row_descs($row_descs);
  $self->col_descs($col_descs);
  $self->numbering($numbering);
  $self->_type($type);
  $self->ffmt($ffmt);
  $self->descffmt($descffmt);
  $self->inplace($inplace);

    # set default ids for the sequences (rows)
  $self->{'row_ids'} ||= [1..$self->height()];
    # set default ids for the columns
  $self->{'col_ids'} ||= [$self->numbering..$self->width()+$self->numbering-1]
    if (defined($self->width()) && $self->width() > 0);
  $self->{'row_descs'} ||= [];
  $self->{'col_descs'} ||= [];

  return 1;
}



=head2 _rearrange()

 Usage    : n/a (internal function)
 Function : Rearranges named parameters to requested order.
 Returns  : @params - an array of parameters in the requested order.
 Argument : $order : a reference to an array which describes the desired
                     order of the named parameters.
            @param : an array of parameters, either as a list (in
                     which case the function simply returns the list),
                     or as an associative array (in which case the
                     function sorts the values according to @{$order}
                     and returns that new array.

=cut

sub _rearrange {
  # This function was taken from CGI.pm, written by Dr. Lincoln
  # Stein, and adapted for use in Bio::Seq by Richard Resnick.
  my($self,$order,@param) = @_;

  # If there are no parameters, we simply wish to return
  # an empty array which is the size of the @{$order} array.
  return ('') x $#{$order} unless @param;

  # If we've got parameters, we need to check to see whether
  # they are named or simply listed. If they are listed, we
  # can just return them.
  return @param unless (defined($param[0]) && $param[0]=~/^-/);

  # Now we've got to do some work on the named parameters.
  # The next few lines strip out the '-' characters which
  # preceed the keys, and capitalizes them.
  my $i;
  for ($i=0;$i<@param;$i+=2) {
	$param[$i]=~s/^\-//;
	$param[$i]=~tr/a-z/A-Z/;
  }
  
  # Now we'll convert the @params variable into an associative array.
  my(%param) = @param;

  my(@return_array);
  
  # What we intend to do is loop through the @{$order} variable,
  # and for each value, we use that as a key into our associative
  # array, pushing the value at that key onto our return array.
  my($key);

  foreach $key (@{$order}) {
	my($value) = $param{$key};
	delete $param{$key};
	push(@return_array,$value);
  }
  
  # catch user misspellings resulting in unrecognized names
  my(@restkeys) = keys %param;
  if (scalar(@restkeys) > 0) {
       carp("@restkeys not processed in _rearrange(), did you use a
       non-recognized parameter name ? ");
  }

  return (@return_array);
}

# helper functions used by seqs(), etc

=head2 _rowbounds()

 Usage    : $corrected_bounds = $aln->_rowbounds($uncorrected_bounds);
 Function : create default row index list if necessary, 
            create row index list if specified by a hash of ids or 
            by a selector function that acts on rows and returns true/false,
            check row index list for bounds errors, 
            NO LONGER DONE IN VERSION 1.007 AND HIGHER: substract offset of 1.
 Returns  : reference to corrected row index list
 Argument : reference to uncorrected row index list

=cut

sub _rowbounds {
  my($self) = shift;
  my $rowindices = shift;
  my $rowindices2;

  my $firstindx = 1; # currently, it's assumed the first row is row no. 1,
                     # and other conventions aren't supported.

  if ( !defined($rowindices) ) {
      # create default row index list if necessary
    $rowindices2 = [$firstindx .. $self->height() - 1 + $firstindx];
    # cleaner??!! $rowindices2 = [];
  }
  
  if (ref($rowindices) eq 'ARRAY') {
    if (scalar(@$rowindices) == 0) {
      $rowindices2 = [$firstindx .. $self->height() - 1 + $firstindx];
    } else {
      $rowindices2 = [@$rowindices];
    }
  }
 
  if (ref($rowindices) eq 'HASH') {
    my @row_ids = split /\s/, $rowindices->{'ids'};
    my($ii,@rows,$row_id);
    $ii=0;
    my %all_row_ids_hash = map {($_,$ii++)} @{$self->{'row_ids'}};
    foreach $row_id (@row_ids) {
      if (exists $all_row_ids_hash{$row_id}) {
        push @rows, $all_row_ids_hash{$row_id} + $firstindx;
        delete $all_row_ids_hash{$row_id};
      } else {
        carp("Requested row named $row_id not found, or duplicate");
      }
    }
    $rowindices2 = \@rows;
  }

  if (ref($rowindices) eq 'CODE') {
    my $ctr = $firstindx - 1;
    $rowindices2 = [ grep {$_ != $_NaN}
                      # filter out $_NaN-value indexes, see next comment
                    map {($ctr++,$_) ? $ctr : $_NaN} 
                      # convert list of function values into list of indices
                      # in canonical order, except that false values trigger an
                      # index of $_NaN, e.g. (1,$_NaN,3,$_NaN,5) is the result 
                      # of map {($ctr++,$_) ? $ctr : $_NaN} ("A","","B",0,"C")
                      # note that the comma operator evaluates the first
                      # argument and throws the result away, and then 
                      # evaluates the second argument; the side effect
                      # of the first evaluation is the increment of $ctr.
                    $self->map_r($rowindices) ];
  }

  my($ii);
  for $ii (0 .. $#{$rowindices2}) { 
    if ($rowindices2->[$ii] < $firstindx) {
      carp "Requested row $rowindices2->[$ii] of an alignment where $firstindx is the first row";
    }
    if (($rowindices2->[$ii] - $firstindx + 1) > scalar(@{ $self->{'seqs'} })) {
      carp "Requested row $rowindices2->[$ii] of an alignment with fewer rows";
    }
  }
  return $rowindices2;
  #return $self->_rownorm($rowindices2);
}


=head2 _colbounds()

 Usage    : $corrected_bounds = $aln->_colbounds($uncorrected_bounds);
 Function : create default column index list if necessary,
            create column index list if specified by a hash of ids or 
            by a selector function that acts on columns and returns true/false,
            check column index list for bounds errors, 
            NOT IN 1.007 and higher: substract offset (according to numbering scheme).
 Returns  : reference to corrected column index list
 Argument : reference to uncorrected column index list
            reference to row index list: its first row is used to 
            provide the width of the alignment considered, which matters 
            if the alignment is really a sequence bag

=cut

sub _colbounds {
  my($self) = shift;
  my $colindices = shift;
  my $rowindices = shift || [0];
  my $colindices2;
  my $firstindx2 = 1; # currently, it's assumed the first row is row no. 1,
                    # and other conventions aren't supported.
  my $width;
  if (scalar(@$rowindices) != 0) {
    $width = $self->width($rowindices->[0]+$firstindx2);
  } else {
    return [];
  }
  my $firstindx = $self->numbering;

  if ( !defined($colindices) ) {
      # create default col index list if necessary
    $colindices2 = [$firstindx .. $width - 1 + $firstindx]
      if (defined($width) && $width > 0);
    #??!! $colindices2 = [];
  }
  
  if (ref($colindices) eq 'ARRAY') {
    if (scalar(@$colindices) == 0) {
      $colindices2 = [$firstindx .. $width - 1 + $firstindx]
        if (defined($width) && $width > 0);
    } else {
      $colindices2 = [@$colindices];
    }
  }
 
  if (ref($colindices) eq 'HASH') {
    my @col_ids = split /\s/, $colindices->{'ids'};
    my($ii,@cols);
    my %col_ids_hash = map {($_,"")} @col_ids;
    foreach $ii (0..$#{$self->{'col_ids'}}) {
      my $col_id = $self->{'col_ids'}[$ii];
      if (exists $col_ids_hash{$col_id}) {
        push @cols, $ii+$firstindx; 
        delete $col_ids_hash{$col_id};
      }
    }
    my(@restkeys) = keys %col_ids_hash;
    if (scalar(@restkeys) > 0) {
      carp("Requested columns named @restkeys not found");
    }
    $colindices2 = \@cols;
  }
  
  if (ref($colindices) eq 'CODE') {
    my $ctr = $firstindx - 1;
    $colindices2 = [ grep {$_ != $_NaN}                    
                      # see _rowbounds()
                    map {($ctr++,$_) ? $ctr : $_NaN} 
                    $self->map_c($colindices) ];
  }

  my($ii);
  for $ii (0 .. ($#{$colindices2})) {
    if ($colindices2->[$ii] < $firstindx) {
      carp "Requested column $colindices2->[$ii] of an alignment starting at $firstindx";
    }
    if (($colindices2->[$ii] - $firstindx + 1) > $width) {
      carp "Requested column $colindices2->[$ii] is beyond the end of the alignment";
    }
  }
  return $colindices2;
  #return $self->_colnorm($colindices2);
}


    # normalize by substracting the first row index (offset)
sub _rownorm {
  my($self) = shift;
  my $rowindices = [@{$_[0]}];

  my $firstindx = 1; # currently, it's assumed the first row is row no. 1,
                    # and other conventions aren't supported.

  return [ map {$_ -= $firstindx} @$rowindices ];
}

    # normalize by substracting the first column index (offset)
sub _colnorm {
  my($self) = shift;
  my $colindices = [@{$_[0]}];

  my $firstindx = $self->numbering;

  return [ map {$_ -= $firstindx} @$colindices ];
}

=head2 _fixbounds()

 Usage    : ($corrected_rowbounds,$corrected_colbounds) = 
               $aln->_fixbounds($uncorrected_rowbounds,$uncorrected_colbounds);
                    OR
            ($corrected_rowbounds,$corrected_colbounds) = 
               $aln->_fixbounds($firstpos1,$lastpos1,$firstpos2,$lastpos2)
 Function : Convert unfixed index list information into the standard internal
            one, allowing as input either max. 2 references or max. 4 boundary 
            coordinates (2 coord. for row indices and 2 for column indices).
            Call functions to create maximal default index lists if needed, 
            to create index lists if specified by a hash of ids or by a selector
            function that acts on rows/columns and returns true/false,
            to check index lists for bounds errors, and to substract offsets.
 Returns  : 2 references to corrected index lists
 Argument : EITHER max. 2 references to uncorrected index lists, 
            OR max. 4 boundary coordinates (integers)

=cut

sub _fixbounds {
  my($self) = shift;
  my($rrowsel2,$rcolsel2);
  my($ii);

  if (!(defined($_[0])) || (ref($_[0]) =~ /ARRAY|CODE|HASH/)) {
    my($rrowsel,$rcolsel) = @_;
    $rrowsel2 = $self->_rownorm($self->_rowbounds($rrowsel));
    $rcolsel2 = $self->_colnorm($self->_colbounds($rcolsel,$rrowsel2));
  } else {
      # unfortunately, @{[undef..undef]} is ``0'', and not ``undef'', so we 
      # need to create valid bounds right here
    my($firstpos1,$lastpos1,$firstpos2,$lastpos2) = @_;
    my $firstindx1 = 1;
    my $firstindx2 = $self->numbering;
  
    if (!defined($firstpos1)) {
      $firstpos1 = $firstindx1;
    }
    if (!defined($lastpos1)) {
      $lastpos1 = $self->height() - 1 + $firstindx1;
    }
    if (!defined($firstpos2)) {
      $firstpos2 = $firstindx2;
    }
    if (!defined($lastpos2)) {
      $lastpos2 = $self->width() - 1 + $firstindx2;
    }
    $rrowsel2 = $self->_rownorm($self->_rowbounds([$firstpos1..$lastpos1]));
    $rcolsel2 = $self->_colnorm($self->_colbounds([$firstpos2..$lastpos2]));
  }

  return ($rrowsel2,$rcolsel2);
}

=head2 _select()

 Usage    : @alnSlice = $aln->_select($rrowsel2,$rcolsel2);
 Function : select elements from the 2-dimensional array $self->{'seqs'} 
            using index lists
 Returns  : array of references to array of characters
 Argument : 2 index lists, one for the rows and one for the columns,
 Comment  : Here's a diagram of dependencies of some methods relying on _select:

                       seqs()
                      /    |
                    |/_   \|/                                     
                _arys() <-- _strs() 
                  |
                 \|/
               _select()

=cut


sub _select {
  my $self = shift;
  my $array = $self->{'seqs'};
  my $rrowsel = shift; # right row selector; function or index array
  my $rcolsel = shift; # right column selector; function or index array
  my @slice = map [ @$_[ @$rcolsel ] ], map $array->[$_], @$rrowsel ;
  if ($self->{'inplace'}) {
    $self->{'seqs'}    = [@slice]; 
    $self->{'id'}      = $self->{'id'}.'_';
    $self->{'desc'}    = "from ".$self->{'id'}.",".$self->{'desc'};
  
    $self->{'row_ids'} = [ @{ $self->{'row_ids'} }[@$rrowsel] ]
      unless ($#{$self->{'row_ids'}}+1 == 0);
    $self->{'row_descs'} = [ @{ $self->{'row_descs'} }[@$rrowsel] ]
      unless ($#{$self->{'row_descs'}}+1 == 0);
    $self->{'col_ids'} = [ @{ $self->{'col_ids'} }[@$rcolsel] ]
      unless ($#{$self->{'col_ids'}}+1 == 0);
    $self->{'col_descs'} = [ @{ $self->{'col_descs'} }[@$rcolsel] ]
      unless ($#{$self->{'col_descs'}}+1 == 0);
    return (1);
  } else {
    return @slice;
  }
}

=head2 _arys()

 Usage    : @alnSlice = $aln->_arys([$row1,$row2,$rowx],[col1,$col2,$colx])
              (other usages see seqs())
 Function : Same as seqs(), except that an array is returned always
 Returns  : array of references to array of characters
 Argument : Same as seqs()

=cut

# friendly frontend to _select()
sub _arys {
  my $self = shift;
  my($rrowsel,$rcolsel) = $self->_fixbounds(@_);

  return $self->_select($rrowsel,$rcolsel);
}


# convert (array of references to arrays of characters) into (multiline string)
sub _stringify_spaced {
  my ($seqref,$seqs);
  $seqs="";

  foreach $seqref (@_) {
    if (!($seqref && (ref($seqref) eq 'ARRAY'))) {
      carp "Cannot _stringify_spaced $seqref";
      next;
    }
    $seqs.=join(" ",@$seqref);
    $seqs.="\n";
  }
  return $seqs;
}

# convert (array of references to arrays of characters) into (multiline string)
sub _stringify {
  my ($seqref,$seqs);
  $seqs="";

  foreach $seqref (@_) {
    if (!($seqref && (ref($seqref) eq 'ARRAY'))) {     
      carp "Cannot _stringify $seqref";
      next;
    }
    $seqs.=join("",@$seqref);
    $seqs.="\n";
  }
  return $seqs;
}


=head2 _strs()

 Usage    : $alnSlice = $aln->_strs([$row1,$row2,$rowx],[col1,$col2,$colx])
              (other usages see seqs())
 Function : Same as seqs(), except that a multiline string is returned always
 Returns  : multiline string (including newline characters)
 Argument : Same as seqs()

=cut

# friendly frontend to _select(), via _arys()
sub _strs {
  my($self) = shift;
  return _stringify($self->_arys(@_));
}


=head2 seqs()

 Usage    : (1) $alnSlice = $aln->seqs($firstpos1,$lastpos1,$firstpos2,$lastpos2)
            (2) $alnSlice = $aln->seqs([$row1,$row2,$rowx],[col1,$col2,$colx])
            (3) $alnSlice = $aln->seqs(\&function_of_row,\&function_of_column)
            (4) $alnSlice = $aln->seqs({ids=>'r1 r2 r3'},{ids=>'c1 c2 c3'})
               [ (2),(3) and (4) can be intermixed ]
 Function : (1) Returns part of an alignment, from row $firstpos1 to row
            $lastpos1, and from column $firstpos2 to column $lastpos2
            Missing parameters are replaced by default (maximal possible)
            values, where the length of the first row of the _returned_ (*)
            `alignment' determines the default length in case the alignment is
            really a sequence bag. ((*) is a new feature in 1.006.)
            (2) Returns part of an alignment, i.e. those elements that lie
            in a row designated by an index from the first list       
            ([$row1,$row2,$rowx]), and at the same time lie in a column 
            designated in the second list ([col1,$col2,$colx]).
            The empty list (``[]'') is replaced by default (maximal possible)
            values, where the length of the first row of the _returned_
            `alignment' determines the default length in case the alignment is
            really a sequence bag. (Note that the first element
            of rows/columns has index 1 per default.)
            (3) Instead of an index list, a function acting on a row / column
            may be supplied; whenever the function returns true, the
            row / column is designated.
            (4) Instead of an index list, a hash of ids may be supplied;
            the ids are looked up in the alignment's list of row (sequence)
            ids / list of column (site) ids. The former list may be set
            during the construction of the alignment (e.g. it may be read
            from the fasta file, or the Bio::Seq objects), or it may be
            manipulated using row_ids(). col_ids() sets the the latter list.
 Returns  : in a scalar context: multiline string (including newline characters)
            in an array context: array of references to arrays of characters
 Argument : (1) $firstpos1,$lastpos1,$firstpos2,$lastpos2 (all integers; note 
            that the first element of rows/columns has index 1 per default.)
            (2-4) 2 selectors, one for the rows and one for the columns,
            each of which may be a reference to a list of indices, or a
            reference to a hash that has the format ``{ids=>'id1 id2 idx'}''
            where ids is the mandatory key, and the value is a string 
            containing the desired ids, seperated by `` '' (space), or a
            function that acts on a list and returns true/false.

=cut

sub seqs {
  my $self = shift;

  if ($self->{'inplace'}) {
    return $self->_arys(@_);
  } else {
    return wantarray ? $self->_arys(@_) : $self->_strs(@_);   
  }
}

sub aln {
  my($self) = shift;
  my(%copy);

  my($rrowsel,$rcolsel) = $self->_fixbounds(@_);
  $copy{'seqs'}    = [$self->_select($rrowsel,$rcolsel)]; 
  $copy{'id'}      = $self->{'id'}.'_';
  #$copy{'id'}      = $self->{'id'}.'_'.join('_',split(/\s/,_stringify(\@_)));
  $copy{'desc'}    = "sliced copy of ".$self->{'id'}.",".$self->{'desc'};
  foreach (keys %{$self->{'names'}}) {
    $copy{'names'}{$_} = $self->{'names'}{$_};
  }
  if ($#{$self->{'row_ids'}}+1 != 0) {
    $copy{'row_ids'} = [ @{ $self->{'row_ids'} }[@$rrowsel] ]
  } else {
    $copy{'row_ids'} = [];
  }
  if ($#{$self->{'row_descs'}}+1 != 0) {
    $copy{'row_descs'} = [ @{ $self->{'row_descs'} }[@$rrowsel] ]
  } else {
    $copy{'row_descs'} = [];
  }
  if ($#{$self->{'col_ids'}}+1 != 0) {
    $copy{'col_ids'} = [ @{ $self->{'col_ids'} }[@$rcolsel] ]
  } else {
    $copy{'col_ids'} = [];
  }
  if ($#{$self->{'col_descs'}}+1 != 0) {
    $copy{'col_descs'} = [ @{ $self->{'col_descs'} }[@$rcolsel] ]
  } else {
    $copy{'col_descs'} = [];
  }

  $copy{'numbering'} = $self->{'numbering'};
  $copy{'type'}    = [ @{ $self->{'type'} } ];
  $copy{'ffmt'}    = $self->{'ffmt'};
  $copy{'descffmt'}= $self->{'descffmt'};
  $copy{'inplace'}= $self->{'inplace'};

  return bless \%copy, ref($self);
}


=head2 _map_r()

 Usage    : n/a (internal function)
 Function : apply a row function to selected rows, 
            then return the list of all return values
 Returns  : list of all return values
 Argument : $rowf: Reference to the function that is applied to selected rows, 
            its results are put into a list and returned.
            $rrowsel: Reference to the selector designating the rows 
            to which rowf is applied, as in cases (2)-(4) in seqs().

=cut

sub _map_r {
# Developed from an expression from Michael Constant
  my($self) = shift;
  my ($rowf,$rrowsel) = @_;
  my $array = $self->{'seqs'};
  my (@results);
  if (defined($rrowsel)) {
    push @results, map { 
      &$rowf($_)   
    } ref($rrowsel) eq 'CODE'
           ? grep &$rrowsel($_), @{$array}
           : map $array->[$_], @$rrowsel
      ;
  } 
  return @results;
}


=head2 _map_c()

 Usage    : n/a (internal function)
 Function : apply a column function to selected columns,
            then return the list of all return values
 Returns  : list of all return values
 Argument : analogous to _map_r()

=cut

sub _map_c {
# Developed from an expression from Michael Constant
  my($self) = shift;
  my ($colf,$rcolsel) = @_;
  my ($colnr,@colnums, $colnum);
  my $array = $self->{'seqs'};
  my (@results);
  if (defined($rcolsel)) {
    push @results, map { 
      $colnr=$_;
      &$colf([ map { ${$_}[$colnr] } @{ $self->{'seqs'} } ])  
        # Basically, &$colf([ map { ... } ]) is calculated; [ ] creates
        # the reference so that call-by-reference can be done.
        # In the outer map, $colnr=$_ loops thru all indices in @$rcolsel
        # In the inner map, $_ loops thru all rows _and_ from each row
        # the element row[$colnr] is taken. The inner map assembles these
        # elements into a list which is passed to &$colf.
    } ref($rcolsel) eq 'CODE'
        ? @colnums = @colnums
                     ? @colnums
                     : grep {
                             $colnum=$_,
                             &$rcolsel( [ map ${$_}[$colnum], @{$array} ] )
                            } 0..$#{$array->[0]}
        : @$rcolsel;
      # The outmost map processes each desired column, one at time, selected 
      # by either the last grep, or the list of indices, $rcolsel.
      # For each column, if ref($rcolsel) eq 'CODE' DOES hold,
      # desired column numbers are those for which &$rcolsel
      # returns true; @colnums is used as a cache so that they are only
      # calculated once. If the cache is empty, all column indices 
      # (that is, 0..$#{$array[0]}) are processed by a grep which first
      # notes the column index ($colnum memorizes the running variable $_)
      # and then calls $rcolsel with an argument that evaluates to a reference
      # (note the ``[..]'') to the current column. This current column is 
      # obtained by presenting each array row to the innermost map, and 
      # extracting just one element: ${$_}[$colnum], where $_ is the array row,
      # and ${$_}[$colnum] is the row's element belonging to the current column.
      # The single elements are then put into an array, and passed to $rcolsel
  }
  return @results;
}


=head2 map_r()

 Usage    : @resSlice = $aln->map_r($rowf,$rrowsel);
 Function : apply a function to selected rows, 
            then return the list of all return values
 Returns  : list of all return values
 Argument : $rowf: Reference to the function that is applied to selected rows, 
            its results are put into a list and returned.
            $rrowsel: Reference to the selector designating the rows 
            to which rowf is applied, as in cases (2)-(4) in seqs().

=cut

sub map_r {
  my($self,$rowf,$rrowsel) = @_;
  my($rrowsel2);

  if (defined($rowf)) {
    $rrowsel2 = ref($rrowsel) eq 'CODE'
                  ? $rrowsel
                  : $self->_rownorm($self->_rowbounds($rrowsel));
  }

  return $self->_map_r($rowf,$rrowsel2);
}

=head2 map_c()

 Usage    : @resSlice = $aln->map_c($colf,$rcolsel,$rrowsel);
 Function : apply a function to selected columns,
            then return the list of all return values
 Returns  : list of all return values
 Argument : $colf: Reference to the function that is applied to selected 
            columns, its results are put into a list and returned.
            $rcolsel: Reference to the selector designating the columns
            to which colf is applied, as in cases (2)-(4) in seqs().
            $rrowsel: NOT used as a selector, but as a hint for 
            determining the width of a sequence bag if $rcolsel is undef:
            The last index of the first row specified by $rrowsel is 
            taken as the maximum column index.

=cut

sub map_c {
  my($self,$colf,$rcolsel,$rrowsel) = @_;
  my($rcolsel2);

  if (defined($colf)) {
    $rcolsel2 = ref($rcolsel) eq 'CODE'
                  ? $rcolsel
                  : $self->_colnorm($self->_colbounds($rcolsel,$rrowsel));
  }

  return $self->_map_c($colf,$rcolsel2);
}


# Using closures; The following functions are the templates.

sub _c_consensus_of_array {
  my $threshold = shift;
  return sub {
    my %temp = ();
    my @chars = @{ $_[0] };
    my @list = sort { $temp{$a}<=>$temp{$b} }
               #sort { $a cmp $b } #!!!??? Not OK, since next sort scrambles it
               grep ++$temp{$_} >= $threshold * ($#chars+1), 
               @chars;
      # In case of a tie, it's not specified which residue is in $list[-1]
    
    return (defined($list[-1]) ? $list[-1] : $_NO_CONSENSUS_SYMBOL);
  };
}

sub _c_consensus_doesnt_exist {
  my $threshold = shift;
  my $_consensus_of_array = _c_consensus_of_array($threshold);
  return sub {
    return ((&$_consensus_of_array(@_) eq $_NO_CONSENSUS_SYMBOL)) ? 1 : 0;
  }
}

sub _c_consensus_exists {
  my $threshold = shift;
  my $_consensus_of_array = _c_consensus_of_array($threshold);
  return sub {
    return ((&$_consensus_of_array(@_) ne $_NO_CONSENSUS_SYMBOL)) ? 1 : 0;
  }
}

=head2 consensus()

 Usage    : $cons_letters = $aln->consensus($threshold,$rcolsel);
 Function : return the consensus of a (subset of) the columns; the letter '!'
            ($_NO_CONSENSUS_SYMBOL) indicates that no consensus letter exists
 Returns  : in a scalar context: string of consensus letters
            in an array context: array of consensus letters
 Argument : $threshold: A letter is considered consensus of a column
            if the fraction of the letters in the column that form a
            (relative) majority is >= $threshold. Ties between 2 letters with 
            an equal relative majority are broken arbitrarily.
            The default value is 0.75.
            $rcolsel: Reference to the selector designating the columns
            of which the consensus is calculated, as in cases (2)-(4) in seqs().

=cut

sub consensus {
  my $self = shift;
  my($threshold,$rcolsel) = @_;
  $threshold = 0.75 unless defined($threshold);
  my $_consensus_of_array = _c_consensus_of_array($threshold);
  my @consensus = $self->map_c(\&$_consensus_of_array,$rcolsel);
  return wantarray ? @consensus : join("",@consensus);
}

=head2 var_sites()

 Usage    : $resSlice = $aln->var_sites($threshold,$rrowsel);
 Function : return the variable sites of an alignment
 Returns  : in a scalar context: multiline string (including newline characters)
            in an array context: array of references to arrays of characters
 Argument : $threshold: A column is considered variable 
            if the fraction of the letters in the column that form a
            (relative) majority is NOT >= $threshold. 
            The default $threshold is 1, i.e. only constant, INvariable columns 
            are excluded.
            $rrowsel: Reference to the selector designating the rows
            of which in turn the letters in the variable columns are printed,
            as in cases (2)-(4) in seqs(). $rrowsel DOES NOT influence the
            calculation of the variable sites, it just constrains the output
            further !

=cut

sub var_sites {
  my $self = shift;
  my($threshold,$rrowsel) = @_;
  $threshold = 1 unless defined($threshold);
  my $_consensus_doesnt_exist = _c_consensus_doesnt_exist($threshold);
  return $self->seqs($rrowsel,\&$_consensus_doesnt_exist);
}

=head2 var_inds()

 Usage    : $indices = $aln->var_inds($threshold,$rrowsel);
 Function : return the _indices_ of the variable sites of an alignment
 Returns  : reference to array of indices
 Argument : $threshold: A column is considered variable
            if the fraction of the letters in the column that form a
            (relative) majority is NOT >= $threshold.
            The default $threshold is 1, i.e. only constant, INvariable columns
            are excluded.
            $rrowsel: the first row in $rrowsel is used for
            determining the width of a sequence bag. In other words,
            the last index of the first row specified by $rrowsel is
            taken as the maximum column index.

=cut

sub var_inds {    
  my $self = shift;
  my($threshold,$rrowsel) = @_;
  $threshold = 1 unless defined($threshold);
  my $_consensus_doesnt_exist = _c_consensus_doesnt_exist($threshold);
    # here, rrowsel is just used as a hint to compute sequence bag width !!
    #       The last index of the first row specified by $rrowsel is         
    #       taken as the maximum column index.  
  $rrowsel = $self->_rowbounds($rrowsel);
  my $res = $self->_colbounds(\&$_consensus_doesnt_exist,$rrowsel);
  return $res;  #causes weird errors: wantarray ? @$res : join("$,",@$res);
}

=head2 invar_sites()

 Usage    : $resSlice = $aln->invar_sites($threshold,$rrowsel);
 Function : return the INvariable columns of an alignment
 Returns  : in a scalar context: multiline string (including newline characters)
            in an array context: array of references to arrays of characters
 Argument : $threshold: A column is considered INvariable 
            if the fraction of the letters in the column that form a
            (relative) majority is >= $threshold. 
            The default $threshold is 1, i.e. no variability is allowed.
            $rrowsel: see var_sites()

=cut

sub invar_sites {
  my $self = shift;
  my($threshold,$rrowsel) = @_;
  $threshold = 1 unless defined($threshold);
  my $_consensus_exists = _c_consensus_exists($threshold);
  return $self->seqs($rrowsel,\&$_consensus_exists);
}

=head2 invar_inds()

 Comment  : This is the indices-returning version of invar_sites(),
            cf. var_sites() and var_inds().

=cut

sub invar_inds {
  my $self = shift; 
  my($threshold,$rrowsel) = @_;
  $threshold = 1 unless defined($threshold);
  my $_consensus_exists = _c_consensus_exists($threshold);
    # here, rrowsel is just used as a hint to compute sequence bag width !!
    #       The last index of the first row specified by $rrowsel is         
    #       taken as the maximum column index.  
  $rrowsel = $self->_rowbounds($rrowsel);
  my $res = $self->_colbounds(\&$_consensus_exists,$rrowsel);
  return $res;  #causes weird errors: wantarray ? @$res : join("$,",@$res);
} 

sub _gap_free {
  my $str = join "", @{ $_[0] };
  if ($str =~ /$_GAP_SYMBOL/) {return 0;} else {return 1;}
}

=head2 gap_free_sites()

 Usage    : $resSlice = $aln->gap_free_sites($rrowsel);
 Function : return the gap-free columns of an alignment.
 Returns  : in a scalar context: multiline string (including newline characters)
            in an array context: array of references to arrays of characters
 Argument : $rrowsel: see var_sites()

=cut

sub gap_free_sites {
  my $self = shift;
  my($rrowsel) = shift;
  return $self->seqs($rrowsel,\&_gap_free);
}
sub gap_free_cols {
  my $self = shift;
  carp "gap_free_cols deprecated -- use the equivalent gap_free_sites";
  $self->gap_free_sites(@_);
}

=head2 gap_free_inds()

 Comment  : This is the indices-returning version of gap_free_sites(), 
            cf. var_sites() and var_inds().

=cut

sub gap_free_inds {
  my $self = shift;
  my($rrowsel) = shift;
    # here, rrowsel is just used as a hint to compute sequence bag width !!
    #       The last index of the first row specified by $rrowsel is         
    #       taken as the maximum column index.  
  $rrowsel = $self->_rowbounds($rrowsel);
  my $res = $self->_colbounds(\&_gap_free,$rrowsel);
  return $res;  #causes weird errors: wantarray ? @$res : join("$,",@$res);
}

sub _unknown_free {
  my $str = join "", @{ $_[0] };
  if (($str =~ /$_UNKN_SYMBOL/) ||
      ($str =~ /$_UNKN_SYMBOL2/)) {return 0;} else {return 1;}
}

=head2 unknown_free_sites()

 Usage    : $resSlice = $aln->unknown_free_sites($rrowsel);
 Function : return the unknown-free columns of an alignment.
 Returns  : in a scalar context: multiline string (including newline characters)
            in an array context: array of references to arrays of characters
 Argument : $rrowsel: see var_sites()

=cut

sub unknown_free_sites {
  my $self = shift;    
  my($rrowsel) = shift;                                                    
  return $self->seqs($rrowsel,\&_unknown_free);       
}

=head2 unknown_free_inds()

 Comment  : This is the indices-returning version of unknown_free_sites(),
            cf. var_sites() and var_inds().

=cut

sub unknown_free_inds {
  my $self = shift; 
  my($rrowsel) = shift;
    # here, rrowsel is just used as a hint to compute sequence bag width !!
    #       The last index of the first row specified by $rrowsel is         
    #       taken as the maximum column index.  
  $rrowsel = $self->_rowbounds($rrowsel);
  my $res = $self->_colbounds(\&_unknown_free,$rrowsel);
  return $res;  #causes weird errors: wantarray ? @$res : join("$,",@$res);
}

sub _no_allgap {
  my $str = join "", @{ $_[0] };
  if ($str =~ /^$_GAP_SYMBOL*$/) {return 0;} else {return 1;}
}

=head2 no_allgap_sites()

 Usage    : $resSlice = $aln->no_allgap_sites($rrowsel);
 Function : return the columns which do not have gaps only, of an alignment.
 Returns  : in a scalar context: multiline string (including newline characters)
            in an array context: array of references to arrays of characters
 Argument : $rrowsel: see var_sites()

=cut

sub no_allgap_sites {
  my $self = shift;
  my($rrowsel) = shift;
  return $self->seqs($rrowsel,\&_no_allgap);
}

=head2 no_allgap_inds()

 Comment  : This is the indices-returning version of no_allgap_sites(),
            cf. var_sites() and var_inds().

=cut

sub no_allgap_inds {
  my $self = shift;
  my($rrowsel) = shift;
    # here, rrowsel is just used as a hint to compute sequence bag width !!
    #       The last index of the first row specified by $rrowsel is         
    #       taken as the maximum column index.  
  $rrowsel = $self->_rowbounds($rrowsel);
  my $res = $self->_colbounds(\&_no_allgap,$rrowsel);
  return $res;  #causes weird errors: wantarray ? @$res : join("$,",@$res);
}

#sub gap_free_seqs {
  #my $self = shift;
  #my($rcolsel) = shift;
  #return $self->seqs(\&_gap_free,$rcolsel);
#}
#sub gap_free_rows {
  #my $self = shift;
  #$self->gap_free_seqs(@_);
#}

sub _special_free {
  my $str = join "", @{ $_[0] };
  if (($str =~ /$_UNKN_SYMBOL/) ||
      ($str =~ /$_UNKN_SYMBOL2/) ||
      ($str =~ /$_GAP_SYMBOL/)) {return 0;} else {return 1;}
}

=head2 special_free_sites()

 Usage    : $resSlice = $aln->special_free_sites($rrowsel);
 Function : return the special-free (neither gap nor unknown-symbols) columns 
            of an alignment.
 Returns  : in a scalar context: multiline string (including newline characters)
            in an array context: array of references to arrays of characters
 Argument : $rrowsel: see var_sites()

=cut

sub special_free_sites {
  my $self = shift;
  my($rrowsel) = shift;
  return $self->seqs($rrowsel,\&_special_free);
}

=head2 special_free_inds()

 Comment  : This is the indices-returning version of special_free_sites(),
            cf. var_sites() and var_inds().

=cut

sub special_free_inds {
  my $self = shift;
  my($rrowsel) = shift;
    # here, rrowsel is just used as a hint to compute sequence bag width !!
    #       The last index of the first row specified by $rrowsel is         
    #       taken as the maximum column index.  
  $rrowsel = $self->_rowbounds($rrowsel);
  my $res = $self->_colbounds(\&_special_free,$rrowsel);
  return $res;  #causes weird errors: wantarray ? @$res : join("$,",@$res);
}

sub _reverse {
  return [ reverse @{ $_[0] } ];
}

=head2 reverse()

 Usage    : $resSlice = $aln->reverse($rrowsel);
 Function : return the rows of an alignment, in reversed form (right->left) 
 Returns  : in a scalar context: multiline string (including newline characters)
            in an array context: array of references to arrays of characters
 Argument : $rrowsel: Reference to the selector designating the rows
            which are returned in reversed form, as in cases (2)-(4) in seqs().

=cut

sub reverse {
  my $self = shift;
  my($rrowsel) = shift;
  my @rows = $self->map_r(\&_reverse,$rrowsel);
  if ($self->{'inplace'}) {
    my $rrowsel2 = $self->_rownorm($self->_rowbounds($rrowsel));
    $self->{'seqs'}    = [@rows];
    $self->{'row_ids'} = [ @{ $self->{'row_ids'} }[@$rrowsel2] ];
    $self->{'row_descs'} = [ @{ $self->{'row_descs'} }[@$rrowsel2] ];
    $self->{'col_ids'} = [ reverse @{ $self->{'col_ids'} } ];
    $self->{'col_descs'} = [ reverse @{ $self->{'col_descs'} } ];
    return wantarray ? (1) : 1;
  } else {
    return wantarray ? @rows : _stringify(@rows);
  }
}

sub _remove_gaps {
  my $row = join("",@{ $_[0] });
  $row =~ s/$_GAP_SYMBOL//g;
  return [ split("", $row) ];
}

=head2 remove_gaps()

 Usage    : $resSlice = $aln->remove_gaps($rrowsel);
 Function : return the rows of an alignment, gaps removed (i.e. the original 
            sequences)
 Returns  : in a scalar context: multiline string (including newline characters)
            in an array context: array of references to arrays of characters
 Argument : $rrowsel: Reference to the selector designating the rows
            which are returned without gaps, as in cases (2)-(4) in seqs().

=cut

sub remove_gaps {
  my $self = shift;
  my($rrowsel) = shift;
  my @rows = $self->map_r(\&_remove_gaps,$rrowsel);
  if ($self->{'inplace'}) {
    my $rrowsel2 = $self->_rownorm($self->_rowbounds($rrowsel));
    $self->{'seqs'}    = [@rows];
    $self->{'row_ids'} = [ @{ $self->{'row_ids'} }[@$rrowsel2] ];
    $self->{'row_descs'} = [ @{ $self->{'row_descs'} }[@$rrowsel2] ];
    $self->{'col_ids'} = [ ];
    $self->{'col_descs'} = [ ];
    return wantarray ? (1) : 1;
  } else {
    return wantarray ? @rows : _stringify(@rows);
  }
}

#sub _extend_col_ids_if_padded {
#  my $chars = shift;
#  if (join("",@$chars) =~ /^$_GAP_SYMBOL+$/) {
#    return $_GAP_SYMBOL;
#  } else {
#    return shift @{$self->{'col_ids'}};
#  }
#}
#
#sub extend_col_ids_if_padded {
#  my $self = shift;
#  my $col_ids = $self->map_c(\&__extend_col_ids_if_padded);
#  $self->{'col_ids'} = $col_ids;
#}

sub _complement_of_array {
  my $str = join "", @{ $_[0] };
  $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
  return [ split "", $str ];
}

sub _rev_complement_of_array {
  my $str = join "", @{ $_[0] };
  $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
  return [ CORE::reverse split "", $str ];
}

=head2 complement()

 Usage    : $resSlice = $aln->complement($rrowsel);
 Function : return the rows of an alignment, in complemented form 
            In the case of dna/rna, the complement is given according to the 
            IUPAC code; in other cases the result is currently calculated in the
            same way and probably meaningless
 Returns  : in a scalar context: multiline string (including newline characters)
            in an array context: array of references to arrays of characters
 Argument : $rrowsel: see reverse()

=cut

sub complement {
  my $self = shift;
  my($rrowsel) = shift;
  my $type = $self->type();
  if (defined($type) && ($type !~ /rna|dna|unknown|Rna|Dna|Unknown/i)) {
    carp ("Warning: complement'ing an alignment of type $type;
      maybe that's not a good thing to do ?");
  }
  my @rows = $self->map_r(\&_complement_of_array,$rrowsel);
  if ($self->{'inplace'}) {
    my $rrowsel2 = $self->_rownorm($self->_rowbounds($rrowsel));
    $self->{'seqs'}    = [@rows];
    $self->{'row_ids'} = [ @{ $self->{'row_ids'} }[@$rrowsel2] ];
    $self->{'row_descs'} = [ @{ $self->{'row_descs'} }[@$rrowsel2] ];
    return wantarray ? (1) : 1;
  } else {
    return wantarray ? @rows : _stringify(@rows);
  }
}

=head2 revcom()

 Usage    : $aln->revcom($rrowsel);
 Function : return the rows of an alignment, in reversed complement form 
            In the case of dna/rna, the complement is given according to the 
            IUPAC code; otherwise the result is currently calculated in the
            same way and probably nonsensical.
 Returns  : in a scalar context: multiline string (including newline characters)
            in an array context: array of references to arrays of characters
 Argument : $rrowsel: see reverse()

=cut

sub revcom {
  my $self = shift;
  my($rrowsel) = shift;
  my $type = $self->type();
  if (defined($type) && ($type !~ /rna|dna|unknown|Rna|Dna|Unknown/i)) {
    carp ("Warning: complement'ing an alignment of type $type;
      maybe that's not a good thing to do ?");
  }
  my @rows = $self->map_r(\&_rev_complement_of_array,$rrowsel);
  if ($self->{'inplace'}) {
    my $rrowsel2 = $self->_rownorm($self->_rowbounds($rrowsel));
    $self->{'seqs'}    = [@rows];
    $self->{'row_ids'} = [ @{ $self->{'row_ids'} }[@$rrowsel2] ];
    $self->{'row_descs'} = [ @{ $self->{'row_descs'} }[@$rrowsel2] ];
    $self->{'col_ids'} = [ CORE::reverse @{ $self->{'col_ids'} } ];
    $self->{'col_descs'} = [ CORE::reverse @{ $self->{'col_descs'} } ];
    return wantarray ? (1) : 1;
  } else {
    return wantarray ? @rows : _stringify(@rows);
  }
}

=head2 equal_nogaps()

 Usage    : $aln->equal_nogaps($other_aln);
 Function : checks whether two alignments have the same original sequences
            (i.e. gaps are removed and the rows are compared for equality)
 Returns  : 1 if alignments are equal modulo gaps, 0 otherwise
 Argument : $other_aln: the other alignment 

=cut

sub equal_nogaps {
  my $self = shift;
  my($other_aln) = shift;
  if ($self->height() != $other_aln->height()) {
    return 0;
  }
  my $row_index;
  foreach $row_index (1..$self->height()) {
    my $row_gap_free = $self->seqs([$row_index]);
    $row_gap_free =~ s/$_GAP_SYMBOL//g;
    my $other_row_gap_free = $other_aln->seqs([$row_index]);
    $other_row_gap_free =~ s/$_GAP_SYMBOL//g;
    if ($row_gap_free ne $other_row_gap_free) {
      print STDERR "Not equal: in row $row_index, $row_gap_free and $other_row_gap_free";
      return 0;
    }
  }
  return 1;
}

=head2 equalize_lengths()

 Usage    : $aln->equalize_lengths($width);
 Function : modifies the alignment / sequence bag such that all rows
            have length $width; in the case of sequence bags, this
            is a primitive way to obtain a true alignment with
            rows of equal length
            * The current procedure is straightforward but slow *
 Returns  : 1 if alignments are equal modulo gaps, 0 otherwise                  
 Argument : $width: the length until which sequences should be padded;
            per default the length of the longest row is taken

=cut                                                
                        
sub equalize_lengths {
  my $self = shift;
  my $width = shift;
  my $row;
  if (!defined($width)) {
    $width = 0;
    foreach $row (@{$self->{'seqs'}}) {
      $width = ($#$row+1) > $width ? ($#$row+1) : $width;
    }
  }
  foreach $row (@{$self->{'seqs'}}) {
    map {$row->[$_] ||= "$_GAP_SYMBOL"} (0..$width-1); #slow!
  }
}

# Accessor functions

=head2 _seqs()

 Usage    : @oldSeqs = $aln->_seqs(@sequences,$start)
            @oldSeqs = $aln->_seqs($sequences,$start)
 Function : to APPEND/OVERWRITE sequences to an alignment
 Returns  : old list of sequences, i.e. ($self->{'seqs'})
 Argument : 1. EITHER a reference to a list of either Bio::Seq objects, or
               arrays of letters, or strings, or any mix of these,
               OR a single (multi-line) string
            2. $rrowsel: Reference to the selector designating the rows
               which shall be appended/overwritten (feature in alpha status)

=cut

sub _seqs {
  my $self = shift;
  my $nseqs = shift;
  my $rrowsel = shift;
  my $rrowsel2 = [];
  if (defined($rrowsel)) {
    $rrowsel2 = $self->_rownorm($self->_rowbounds($rrowsel));
  } 

  my $oseqs = $self->{'seqs'};
  if (defined $nseqs) {
    my ($seq, @seq, $counter); 
    $counter = 0; my $index;
    if (ref($nseqs) eq 'ARRAY') {
      foreach $seq (@$nseqs) {
        $index = defined($rrowsel2->[$counter]) ? 
                          $rrowsel2->[$counter] : $counter;
        if (ref($seq) eq 'ARRAY') {
                   ##print "\nWriting array[$rrowsel2->[$counter] || $counter] ", @$seq;
           $self->{'seqs'}[$index] = $seq; 
        } elsif (ref($seq) eq "") {  #a string???
                   ##print "\nWriting string[$rrowsel2->[$counter] || $counter] $seq";
           @seq = split("", $seq);
           $self->{'seqs'}[$index] = [ @seq ];
        } elsif (ref($seq) =~ /Seq/) {
               # Currently, there doesn't seem to be another way to find out
               # whether we've been passed a usable object.
               # However, with Perl 5.004, we'll use can() provided 
               # by class ``UNIVERSAL'' to test whether
               # the object passed to us is usable, i.e. has methods ary(), etc
                   ##print "\nWriting sequence[$rrowsel2->[$counter] || $counter] ",$seq->ary;
           $self->{'seqs'}[$index] = [ split(//, $seq->seq()) ];
           $self->{'row_ids'}[$index] = $seq->id();
           $self->{'row_descs'}[$index] = $seq->desc();
        } else {
          carp("sequence nr. $index could not be extracted from input parameter ``-seqs''");
        }
        $counter++;
      }
    } elsif (ref($nseqs) eq "") {  #a string???
      foreach $seq (split("\n", $nseqs)) {
        $index = defined($rrowsel2->[$counter]) ?
                          $rrowsel2->[$counter] : $counter;
        @seq = split("", $seq);
        $self->{'seqs'}[$index] = [ @seq ];
        $counter++;
      } 
    } else {
      $index = defined($rrowsel2->[$counter]) ?
                        $rrowsel2->[$counter] : $counter;
      $self->{'seqs'}[$index] = [ ];
      $counter++;
      carp("alignment sequences could not be extracted from input parameter ``-seqs''");
    }
  }
  return $oseqs;
}




=head2 id()

 Usage    : $aln_id = $aln->id();
            $aln->id($id_string);
 Function : Accessor, also sets field if an ID is passed in. 
 Returns  : (original) ID value
 Argument : sequence string with no whitespace

=cut


sub id {
  my($self) = shift;
  my($nid) = @_;

  my $oid = $self->{'id'};
  if (defined($nid)) {
    $self->{'id'} = $nid;
    $self->{'id'} = undef if $nid eq "-undef";
    if ($nid =~ /\s/) {
      carp("identifier $nid has illegal whitespace");
    }
  }

  return $oid;
}


=head2 desc()

 Usage    : $aln_desc = $aln->desc();
            $aln->desc($desc_string);
 Function : Accessor, also sets field if a description string is passed in. 
 Returns  : (original) description value 
 Argument : sequence string

=cut


sub desc {
  my($self) = shift;
  my($ndesc) = @_;

  my $odesc = $self->{'desc'};
  $self->{'desc'} = $ndesc if defined($ndesc);

  return $odesc;
}



=head2 names()

 Usage    : %names = $aln->names;
            $aln->names($hash_ref)
 Function : Accessor, also sets field if names hash is passed in.
            The names hash is 'human-readable' data; each key is a 
            location (whether it be URL, database, database query, etc.) 
            and each value is the id at that location.
 Returns  : (original) names hash value 
 Argument : reference to a hash 

=cut

sub names {
  my($self) = shift;
  my($nnames) = @_;

  my $onames = $self->{'names'};
  if (defined $nnames) {
    if (ref($nnames) ne 'HASH') {
      carp "Not a reference to a hash";
    } else {
      $self->{'names'} = {};
      foreach (keys %$nnames) {
        $self->{'names'}{$_} = $nnames->{$_};
      }
    }
  }

  return $onames;
}


=head2 row_ids()

 Usage    : $row_ids = $aln->row_ids();
            $aln->row_ids($row_ids);
 Function : Accessor, also sets field.
 Returns  : A reference to the (original) array of row (sequence) ids
 Argument : A reference to an array of row (sequence) ids

=cut

sub row_ids {
  my($self) = shift;
  my($nrow_ids) = @_;

  my $orow_ids = $self->{'row_ids'};
  $self->{'row_ids'} = $nrow_ids if defined($nrow_ids);

  return $orow_ids;
}

=head2 col_ids()

 Usage    : $col_ids = $aln->col_ids();
            $aln->col_ids($col_ids);
 Function : Accessor, also sets field.
 Returns  : A reference to the (original) array of column (site) ids
 Argument : A reference to an array of column (site) ids

=cut

sub col_ids {
  my($self) = shift;
  my($ncol_ids) = @_;

  my $ocol_ids = $self->{'col_ids'};
  $self->{'col_ids'} = $ncol_ids if defined($ncol_ids);

  return $ocol_ids;
}

=head2 row_descs()

 Usage    : $row_descs = $aln->row_descs();
            $aln->row_descs($row_descs);
 Function : Accessor, also sets field.
 Returns  : A reference to the (original) array of row (sequence) descriptions
 Argument : A reference to an array of row (sequence) descriptions

=cut

sub row_descs {
  my($self) = shift;
  my($nrow_descs) = @_;

  my $orow_descs = $self->{'row_descs'};
  $self->{'row_descs'} = $nrow_descs if defined($nrow_descs);

  return $orow_descs;
}

=head2 col_descs()

 Usage    : $col_descs = $aln->col_descs();
            $aln->col_descs($col_descs);
 Function : Accessor, also sets field.
 Returns  : A reference to the (original) array of column (site) descriptions
 Argument : A reference to an array of column (site) descriptions

=cut

sub col_descs {
  my($self) = shift;
  my($ncol_descs) = @_;

  my $ocol_descs = $self->{'col_descs'};
  $self->{'col_descs'} = $ncol_descs if defined($ncol_descs);

  return $ocol_descs;
}

=head2 numbering()

 Usage    : $num_start = $aln->numbering;
            $aln->numbering($value);
 Function : Accessor, also sets field if a new numbering scheme is passed in.
 Returns  : (original) numbering value 
 Argument : number that is used as the offset of the first column

=cut

sub numbering {
  my($self) = shift;
  my($nnumbering) = @_;

  my $onumbering = $self->{'numbering'};
  $self->{'numbering'} = $nnumbering if defined($nnumbering);

  # $self->{'col_ids'} ||= 
  # [$self->numbering()..$self->width()+$self->numbering()-1];

  return $onumbering;
}


sub _type {
  my($self) = shift;
  my($ntype) = @_;

  my $otype = $self->{'type'};
  if (defined($ntype)) {
    my($monomer,$samelength);
    $samelength = "unknown";
    if (defined($UnivAlnType{$ntype})) {
      $monomer = $UnivAlnType{$ntype} ;
    } else {
      $monomer = $UnivAlnType{"unknown"};
      carp("Unrecognized alignment type: $ntype");
    }
    $self->{'type'} = [$monomer,$samelength];
  }

  return @{$otype};
}


sub samelength {
  my($self) = shift;
  my($ntype) = @_;

  my $otype = $self->{'type'}[1];
  $self->{'type'}[1] = $ntype if defined($ntype);

  return $otype;
}


=head2 type()

 Usage    : $aln_type = $aln->type;
            $aln->type($value); # May be dangerous !
 Function : Accessor, also sets field if a new type is passed in.
            The latter is considered dangerous !
 Returns  : (original) type value
 Argument : new type, see the list %UnivAlnTypes in the code
            (currently, ``dna'', ``rna'', ``amino'')

=cut

sub type {
  my($self) = shift;
  my($ntype) = @_;

  my $otype = $self->{'type'}[0];
  $self->{'type'}[0] = $ntype if defined($ntype);

  return (defined($TypeUnivAln[$otype])) ? $TypeUnivAln[$otype] : 'undef';
}

=head2 ffmt()

 Usage    : $ffmt = $aln->ffmt;
            $aln->ffmt($value);
 Function : Accessor, also sets field if a new format acronym is passed in.
            This can be done before reading from a file, so that
            the presumably correct parsing routine is called, or
            the value can be set before writing the object using
            layout(), so that layout uses a specified default format.
 Returns  : (original) format acronym
 Argument : string describing the format, see ``Alignment Formats''

=cut


sub ffmt {
  my($self) = shift;
  my($nffmt) = @_;

  my $offmt = $self->{'ffmt'};
  if (defined($nffmt)) {
    if (defined($UnivAlnForm{$nffmt})) {
      $self->{'ffmt'} = $UnivAlnForm{$nffmt};
    }
    else {
      carp("Can't set to unrecognized sequence format: $nffmt");
    }
  }

  return $FormUnivAln[$offmt];  # return acronym, not number code
}


sub descffmt {
  my($self) = shift;
  my($ndescffmt) = @_;

  my $odescffmt = $self->{'descffmt'};
  if (defined($ndescffmt)) {
    if (defined($UnivAlnForm{$ndescffmt})) {
      $self->{'descffmt'} = $UnivAlnForm{$ndescffmt};
    }
    else {
      carp("Unrecognized description format: $ndescffmt");
    }
  }

  return $FormUnivAln[$odescffmt];  # return acronym, not number code
}


=head2 inplace()

 Usage    : $inplace = $aln->inplace();
            $aln->inplace($value);
 Function : Accessor, also sets field if a value is passed in.
            'inplace' is a flag which is set to true if accessors and
            functions should make the modification to the object itself, 
            and just return true on success. Currently, inplace manipulation
            is supported for seqs(), _arys(), _strs(), remove_gaps(),
            var_sites(), invar_sites(), gap_free_sites(), no_allgap_sites(),
            reverse(), complement(), and revcom() (reverse complement).
 Returns  : (original) inplace value
 Argument : currently 0 or 1

=cut

sub inplace {
  my($self) = shift;
  my($ninplace) = @_;

  my $oinplace = $self->{'inplace'};
  $self->{'inplace'} = $ninplace if defined($ninplace);

  return $oinplace;
}


=head2 width()

 Usage    : $width = $aln->width();
 Function : number of columns (of the first row per default)
 Returns  : number of columns (of the first row per default)
 Argument : row 

=cut

sub width {
  my($self) = shift;
  my($row) = shift || 1;
  my $firstindx = 1; # currently, it's assumed the first row is row no. 1,
                    # and other conventions aren't supported.
  $row -= $firstindx;
  if (defined($self->{'seqs'}[$row])) {
    return scalar( @{ $self->{'seqs'}[$row] });
  }
}


=head2 height()

 Usage    : $height = $aln->height();
 Function : number of rows (sequences)
 Returns  : number of rows (sequences)
 Argument : ./.

=cut

sub height {
  my($self) = shift;
  if (defined($self->{'seqs'}[0])) {  # ``[0]'' on purpose, I hope this is ok.
    return scalar( @{ $self->{'seqs'} });
  }
}


# Functions having to do with file formats, parsing,
# formatting, and so on.

sub _alphabet_check_array {
  my @chars = @{ $_[0] };
  my $chars  = join("", @chars);
  my @alph = @Bio::UnivAln::alph;
  my $letter;

                     #print "Sequence being checked =$chars.\n";
                     #print " Alphabet being checked =";
                     #my $aref;
                     #for $aref (@alph) { print $aref; }; print ".\n";

  for $letter (@alph) { #There must be a better way ???
    if ($letter =~ /\w|-/) {
      $chars =~ s/$letter//gi;
    } else {
      $chars =~ s/\{$letter}//gi;
    }
  }
  return $chars;                
}

=head2 alphabet_check()

 Usage    : @offendig_characters = $aln->alphabet_check($rowsel);
 Function : Check rows of the alignment for ``offending characters'', i.e.
            characters that are not (currently) expected to be found in the
            aligned sequences, because they're not in the default alphabet 
            that belongs to the specified L<Alignment Type>.
            Currently, the Alignment Type can only be set explicitly
            via the constructor, or accessor. Since the default alphabet
            is not the IUPAC code (e.g. just A,C,G,T,-,? in case of DNA),
            this is just a proof of concept.
 Returns  : an array containing the offending characters, one string per row
 Argument : $rrowsel: Reference to the selector designating the rows
            to which the check is applied, as in cases (2)-(4) in seqs().

=cut

sub alphabet_check {
  my $self = shift;
  my $rrowsel = shift;
  
  @Bio::UnivAln::alph = @{$UnivAlnAlphs{$UnivAlnType{$self->type()}."Mg"}};
  return $self->map_r(\&_alphabet_check_array,$rrowsel);  
}


=head2 _file_read()

 Usage    : n/a (internal function)
 Function : Read data from file, and call the parsing system to store data
            into the object fields
 Returns  : 1 on success
 Argument : filename, and (optionally) the format if known.

=cut

sub _file_read {
  my($self, $filename, $ffmt) = @_;
  my($ent);
  if (!(-e $filename)) {
    carp("File $filename not found");
  }
  if (!(-r $filename)) {
    carp("File $filename not readable");
  }
  open(UnivAln::INPUT, $filename);
  $ent = join("",<UnivAln::INPUT>);
  close(UnivAln::INPUT);
  _parse($self, $ent, $ffmt);

  return 1;
}

sub parse {
  my($self) = shift;
  carp "You've called a parsing routine that may be removed in future versions
        of Bio::UnivAln; currently the only long-term supported method for
        reading in files is via the new() constructor; revision of the parsing
        system should be finalized in mid-1997 ";
  $self->_parse(@_);
}
sub parse_unknown {
  my($self) = shift;
  carp "You've called a parsing routine that may be removed in future versions
        of Bio::UnivAln; currently the only long-term supported method for
        reading in files is via the new() constructor; revision of the parsing
        system should be finalized in mid-1997 ";
  $self->_parse_unknown(@_);
}
sub parse_fasta {
  my($self) = shift;
  carp "You've called a parsing routine that may be removed in future versions
        of Bio::UnivAln; currently the only long-term supported method for
        reading in files is via the new() constructor; revision of the parsing
        system should be finalized in mid-1997 ";
  $self->_parse_fasta(@_);
}
sub parse_raw {
  my($self) = shift;
  carp "You've called a parsing routine that may be removed in future versions
        of Bio::UnivAln; currently the only long-term supported method for
        reading in files is via the new() constructor; revision of the parsing
        system should be finalized in mid-1997 ";
  $self->_parse_raw(@_);
}
sub parse_bad {
  my($self) = shift;
  carp "You've called a parsing routine that may be removed in future versions
        of Bio::UnivAln; currently the only long-term supported method for
        reading in files is via the new() constructor; revision of the parsing
        system should be finalized in mid-1997 ";
  $self->_parse_bad(@_);
}


=head2 _parse()

 Usage    : $aln->_parse($ent,[$ffmt]);
 Function : Parses $ent into the object fields, according to
            $ffmt or $self->{'ffmt'}.
 Returns  : 1 on success
 Argument : the prospective alignment to be parsed, 
            and optionally its format so that it doesn't need to be estimated
            Note that ``raw'' is not estimated in a reliable way by readseq
            if readseq is installed; I presume it is then estimated only if
            the sequence(s) are longer than 111 basepairs


=cut

sub _parse {
  my($self, $ent, $ffmt) = @_;

  if (defined($ffmt) && defined($UnivAlnForm{$ffmt})) {
    $ffmt=$UnivAlnForm{$ffmt};
    } else {
      #warn("Unrecognized or undefined alignment format, using built-in "); 
      # no carp->focussed debugging
      # We later assume that self->{'ffmt'} in a key in $UnivAlnForm if it is defined
    $ffmt=$self->{'ffmt'};
  }

  # We simply need to call the appropriate parsing function, based
  # on the value of $ffmt. Since we've set up the %FuncParse
  # associative array to contain all of the necessary function calls
  # based on the 'ffmt' value, this is rather straightforward.

  if (defined($ffmt) && (exists $FuncParse{$ffmt})) {
    &{$FuncParse{$ffmt}}($self,$ent);
  } else {
    &{$FuncParse{$UnivAlnForm{"unknown"}}}($self,$ent);
  }

  return 1;
}


=head2 _parse_unknown()

 Usage    : $aln->_parse_unknown($ent);
 Function : tries to figure out the format of $ent and then
            calls the appropriate function to parse it into $self->{'seqs'}.
 Returns  : 1 on success
 Argument : $ent : the rough multi-line string to be parsed

=cut

sub _parse_unknown {
  my($self, $ent) = @_;

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
    $self->_parse_fasta($ent);
  } elsif ($ent =~ /^[A-Z_0-9$_GAP_SYMBOL$_UNKN_SYMBOL\s]+$/ig) {
      # currently in raw format, all letters+numbers+gap+missing 
      # symbol is accepted...
    $self->_parse_raw($ent);
  } elsif (($ent =~ /\s*^CLUSTAL/g) && ($_CLUSTAL ne "")) {
    $self->_parse_clustal($ent);
  } elsif ($_READSEQ ne "") { 
    $self->_parse_readseq($ent);
  } else {
      # some other weird format, so we call _parse_bad.
    $self->_parse_bad($ent);
  }

  return 1;
}



=head2 _parse_bad()

 Usage    : $aln->_parse_bad;
 Function : Carp on the bad data that the user gave us.
 Returns  : undef
 Argument : (multiline) string that cannot be parsed

=cut

sub _parse_bad {
  my($self,$ent) = @_;

  carp "Can't parse unknown alignment format:\n$ent";
  return undef
}

=head2 _parse_readseq()

                             
 Usage    : $aln->_parse_readseq;                          
 Function : Try readseq to parse data.
 Returns  : 1 on success
 Argument : $ent : the rough multi-line string to be parsed
    

=cut

sub _parse_readseq {      
  my($self,$ent) = @_;
                                                     
  my $tmp1 = POSIX::tmpnam;                                                 
  open(TMP1,">$tmp1") or carp "Unable to sucessfully open file for readseq";
  print TMP1 $ent;                                                  
  close TMP1;                                                       
  open(READSEQ,"$_READSEQ $tmp1 -pipe -format=Fasta |") or
    carp "Unable to sucessfully open pipe to readseq due to $!";
  $ent = join("\n",<READSEQ>);
  close(READSEQ);
    #going thru fasta format to preserve id/desc
  $self->_parse_fasta($ent);

  return 1;
}

sub _parse_clustal {
  my($self,$ent) = @_;

  my $tmp1 = POSIX::tmpnam;
  open(TMP1,">$tmp1") or carp "Unable to sucessfully open file for CLUSTAL";
  print TMP1 $ent;
  close TMP1;
  my $tmp2 = POSIX::tmpnam;
  open(CLUSTAL,"$_CLUSTAL -infile=$tmp1 -outfile=$tmp2 -convert -output=GDE|") or
    carp "Unable to sucessfully open pipe to CLUSTAL";
  my $msg = join("\n",<CLUSTAL>);
  close(CLUSTAL);
  open(TMP2,"$tmp2") or carp "Unable to sucessfully open file for CLUSTAL";
  $ent = join("\n",<TMP2>);
  close TMP2;
  $ent =~ s/^#/>/mg;
  $ent =~ s/^%/>/mg;
  $self->_parse_fasta($ent);

  return $msg;
}


=head2 _parse_raw()

 Usage    : $aln->_parse_raw;
 Function : parses $ent into the $self->{'seqs'} field, using raw
            file format.
 Returns  : 1 on success
 Argument : (multiline) string to be parsed

=cut

sub _parse_raw {
  my ($self) = shift;
  my ($seqs) = @_;

  $self->_seqs($seqs);
  $self->ffmt('raw');
  $self->descffmt('raw');

  return 1;
}



=head2 _parse_fasta()

 Usage    : $aln->_parse_fasta;
 Function : parses $ent into the 'seqs' field, using fasta
            file format.
 Returns  : 1 on success
 Argument : (multiline) string to be parsed

=cut

sub _parse_fasta {
  my ($self) = shift;
  my ($ent) = @_;
  my (@lines,$line,@seqs);

  @lines = split("\n", $ent);

  my $counter = -1;
  foreach $line (@lines) {
    if ($line =~ /^>[ \t]*(\S*)[ \t]*(.*)$/) {
      $counter++;
      if (defined($1)) {
        my $id = $1;
        $id =~ s/,//;
        $self->{'row_ids'}[$counter] = $id;
      }
      if (defined($2)) {
        $self->{'row_descs'}[$counter] = $2;
      }
      $seqs[$counter] = "";
    }
    elsif ($line =~ /^(.+)$/) {
      if (defined($1)) {
        $seqs[$counter] .= $1;
      }
    }
    else{ next; }
  }

  $self->_seqs([ @seqs ]);
  $self->ffmt('fasta');
  $self->descffmt('fasta');
  return 1;
}

# Functions having to do with file formats, and OUTPUT


=head2 copy()

 Usage    : $copyOfObj = $myUnivAln->copy;
 Function : Returns an identical copy of the object.
 Returns  : Bio::UnivAln
 Argument : n/a

=cut

sub copy {
  my($self) = shift;
  my(%copy);
  my($ii);

  $copy{'seqs'}    = [$self->_arys()]; #SLOW!!
  $copy{'id'}      = $self->{'id'};
  $copy{'desc'}    = $self->{'desc'};
  foreach (keys %{$self->{'names'}}) {
    $copy{'names'}{$_} = $self->{'names'}{$_};
  }
  $copy{'row_ids'} = [ @{ $self->{'row_ids'} } ];
  $copy{'row_descs'} = [ @{ $self->{'row_descs'} } ];
  $copy{'col_ids'} = [ @{ $self->{'col_ids'} } ];
  $copy{'col_descs'} = [ @{ $self->{'col_descs'} } ];
  $copy{'numbering'} = $self->{'numbering'};
  $copy{'type'}    = [ @{ $self->{'type'} } ];
  $copy{'ffmt'}    = $self->{'ffmt'};
  $copy{'descffmt'}= $self->{'descffmt'};
  $copy{'inplace'}= $self->{'inplace'};

  return bless \%copy, ref($self);
}

=head2 layout()

 Usage    : $aln->layout($format);
 Function : Returns the alignment in whichever format the user specifies,
            or according to the "ffmt" field if the user does not specify 
            a format.
 Returns  : varies; "" if unsuccessful
 Argument : $format (one of the formats as defined in $UnivAlnForm). 

=cut

sub layout {
  my($self,$ffmt) = @_;

  $ffmt ||= $self->{"ffmt"};

  # replace internally handled formats by number code
  if (defined($ffmt) && defined($UnivAlnForm{$ffmt})) {
    $ffmt = $UnivAlnForm{$ffmt};
  } 

  if (defined($ffmt) && (exists $FuncOut{$ffmt})) {
    return &{$FuncOut{$ffmt}}($self);
  } elsif ($_READSEQ ne "") {
    return $self->out_readseq($ffmt);
  } else {
    carp "Can't find routine to write format $ffmt";
    return "";
  }
}



=head2 out_bad()

 Usage    : $aln->out_bad;
 Function : Carp if we don't know the output format.
 Returns  : undef
 Argument : n/a

=cut

sub out_bad {
  my($self) = shift;

  carp "Can't write alignment format $self->{\"ffmt\"}";
  return undef;
}




=head2 out_raw()

 Usage    : $aln->out_raw;
 Function : Returns the alignment in raw format.
 Returns  : multiline string
 Argument : n/a

=cut

sub out_readseq {
  my $self = shift;
  my $ffmt = shift;
  my $seqs = $self->out_fasta();
    # The following sometimes causes seg-faults if $ffmt is PAUP or PIR,
    # but seem to work nicely for MSF :-(
    # Example:
    # % cat fastafile | /vol/biotools/bin/readseq -pipe -format=MSF > AAA
    # % cat fastafile | /vol/biotools/bin/readseq -pipe -format=PAUP > AAA
    # Segmentation fault
    # % /vol/biotools/bin/readseq fastafile -pipe -format=PAUP > AAA
    # open(READSEQ,"echo \"$seqs\" | $_READSEQ -pipe -format=$ffmt |") ||
    #   carp "Unable to sucessfully open pipe to readseq due to $!"; 

    # Therefore, better use tmp files:
  my $tmp1 = POSIX::tmpnam;
  open(TMP1,">$tmp1") or carp "Unable to sucessfully open file for readseq";
  print TMP1 $seqs;
  close TMP1;
  open(READSEQ,"$_READSEQ $tmp1 -pipe -format=$ffmt |") ||
    carp "Unable to sucessfully open pipe to readseq due to $!";
  $seqs = join("",<READSEQ>);                        
  close(READSEQ);            
  return $seqs;
}

sub out_raw {
  my($self) = shift;
  # The raw format is just a multiline string 
  return $self->_strs();
}

# UNTESTED AND UNDOCUMENTED PRETTY-PRINTING ROUTINE
sub out_raw2 {
  my($self) = shift;
  my($seq,$row,$id,$desc);

  my($rrowsel,$rcolsel) = $self->_fixbounds(@_);
  my $seqs = "Aln ".$self->id()." (#=".scalar(@{$self->no_allgap_inds()})
    .",var=".scalar(@{$self->var_inds()})."), ".$self->desc()."\n";

  foreach $row (@$rrowsel) {
    $id = $self->{'row_ids'}[$row] || $row;
    my $format = "%".$_MAX_SIZE_TAXANAME."s";
    $seqs .= sprintf("$format",$id);
    $seqs.=" ";
    my $tmp = $self->{'inplace'};
    $self->{'inplace'} = 0; # The user does not want/expect inplace to be acting
                           # in the out_fasta function !
    $seq = [$self->_select([$row],$rcolsel)];
    $self->{'inplace'} = $tmp;
    $seqs.=join("",@{$seq->[0]}); #one row only, it's no. 0 !!!
    $seqs.="\n";
  }
  #carp "? matched" if ($seqs =~ /\?/); #??!!
  return $seqs;
}


=head2 out_fasta()

 Usage    : $aln->out_fasta;
 Function : Returns the alignment as a string in fasta format.
 Returns  : multiline string
 Argument : Same as seqs()
 Comment  : The old ``out_fasta()'' function with no arguments
            is slightly faster and simpler; still available as out_fasta2().

=cut

sub out_fasta2 {
  my($self) = shift;
  my $seqs = "";
  my($seq,$id,$desc);

  my $counter = -1;
  foreach $seq (@{ $self->{'seqs'} }) {
    $counter++;
    $seqs.='> ';
    $id = $self->{'row_ids'}[$counter];
    if (defined($id)) {
     $seqs.=$id;
    }
    $seqs.=" ";
    $desc = $self->{'row_descs'}[$counter];
    if (defined($desc)) {
     $seqs.=$desc;
    }
    $seqs.="\n";
    $seqs.=join("",@$seq);
    $seqs.="\n";
  }
  return $seqs;
}

sub out_fasta {
  my($self) = shift;
  my($seq,$row,$id,$desc);

  my($rrowsel,$rcolsel) = $self->_fixbounds(@_);
  my $seqs = "";

  foreach $row (@$rrowsel) {
    $seqs.='> ';
    $id = $self->{'row_ids'}[$row];
    if (defined($id)) {
      $seqs.=$id;
    }
    $seqs.=" ";
    $desc = $self->{'row_descs'}[$row];
    if (defined($desc)) {
      $seqs.=$desc;
    }
    $seqs.="\n";
    my $tmp = $self->{'inplace'};
    $self->{'inplace'} = 0; # The user does not want/expect inplace to be acting
                           # in the out_fasta function !
    $seq = [$self->_select([$row],$rcolsel)];
    $self->{'inplace'} = $tmp;;
    $seqs.=join("",@{$seq->[0]}); #one row only, it's no. 0 !!!
    $seqs.="\n";
  }
  #carp "? matched" if ($seqs =~ /\?/); #??!!
  return $seqs;
}


sub catch {
  print "Caught internal error";
}

1;
__END__

# The following routine is an undocumented shortcut
sub col {
  my $self = shift;
  my $rcolsel =shift;  # should be a valid column index

  my $array = $self->{'seqs'};
  my @col = map { @$_[$rcolsel] } @$array;

  return wantarray ? @col : _stringify(@col);
}




+increment version number + revision history
diff t/univaln.t2 ~/Bcd/Perl/Bio/Bio/t/univaln.t
diff UnivAln.pm ~/Bcd/Perl/Bio/Bio/UnivAln.pm
cleanup directories
Check pod2html UnivAln.pm > UnivAln.pm.html
rm -fr Makefile blib
perl Makefile.PL
make 
make test
Test w/ t/univaln.t2
rm Makefile t/my_univaln blib
(Create error log file + diff)
Test w/ Adm.pm
cd ~/Perl/Bio/
tar cvf UnivAln-1.009.tar UnivAln-1.009
cp UnivAln-1.009.tar ~/Bcd/Perl/Bio
cd ~/Bcd/Perl/Bio
  tar cvf UnivAln1.008.tar UnivAln
  mv UnivAln1.008.tar ~/Perl/Bio/Docus/
tar xvf UnivAln-1.009.tar
cp UnivAln-1.009.tar UnivAln.tarUnivAln-1.009.tar
gzip UnivAln-1.009.tar
mv UnivAln.tarUnivAln-1.009.tar UnivAln-1.009.tar 
Test at luft
update welcome.html
update publications

possibly into sub _c_consensus_of_array {
    return (defined($list[-1]) 
      ? (($list[-1] =~ /$_GAP_SYMBOL/)&&(_no_allgap(\@chars)==0)) 
        ? $_NO_CONSENSUS_SYMBOL 
        : $list[-1])
      : $_NO_CONSENSUS_SYMBOL);
