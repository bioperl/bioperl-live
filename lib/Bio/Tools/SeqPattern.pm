#
# bioperl module for Bio::Tools::SeqPattern
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by  Steve Chervitz  (sac-at-bioperl.org)
#
# Copyright  Steve Chervitz
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::SeqPattern - represent a sequence pattern or motif

=head1 SYNOPSIS

 use Bio::Tools::SeqPattern;

 my $pat1     = 'T[GA]AA...TAAT';
 my $pattern1 = Bio::Tools::SeqPattern->new(-SEQ =>$pat1, -TYPE =>'Dna');

 my $pat2     = '[VILM]R(GXX){3,2}...[^PG]';
 my $pattern2 = Bio::Tools::SeqPattern->new(-SEQ =>$pat2, -TYPE =>'Amino');

=head1 DESCRIPTION

L<Bio::Tools::SeqPattern> module encapsulates generic data and
methods for manipulating regular expressions describing nucleic or
amino acid sequence patterns (a.k.a, "motifs"), such as the ones produced by
L<Bio::Tools::IUPAC>.

L<Bio::Tools::SeqPattern> is a concrete class that inherits from L<Bio::Seq>.

This class grew out of a need to have a standard module for doing routine
tasks with sequence patterns such as:

  -- Forming a reverse-complement version of a nucleotide sequence pattern
  -- Expanding patterns containing ambiguity codes
  -- Checking for invalid regexp characters
  -- Untainting yet preserving special characters in the pattern

Other features to look for in the future:

  -- Full pattern syntax checking
  -- Conversion between expanded and condensed forms of the pattern

=head1 MOTIVATIONS

A key motivation for L<Bio::Tools::SeqPattern> is to have a way to
generate a reverse complement of a nucleotide sequence pattern.
This makes possible simultaneous pattern matching on both sense and
anti-sense strands of a query sequence.

In principle, one could do such a search more inefficiently by testing
against both sense and anti-sense versions of a sequence.
It is entirely equivalent to test a regexp containing both sense and
anti-sense versions of the *pattern* against one copy of the sequence.
The latter approach is much more efficient since:

   1) You need only one copy of the sequence.
   2) Only one regexp is executed.
   3) Regexp patterns are typically much smaller than sequences.

Patterns can be quite complex and it is often difficult to
generate the reverse complement pattern. The Bioperl SeqPattern.pm
addresses this problem, providing a convenient set of tools
for working with biological sequence regular expressions.

Not all patterns have been tested. If you discover a pattern that
is not handled properly by Bio::Tools::SeqPattern.pm, please
send me some email (sac@bioperl.org). Thanks.

=head1 OTHER FEATURES

=head2 Extended Alphabet Support

This module supports the same set of ambiguity codes for nucleotide
sequences as supported by L<Bio::Seq>. These ambiguity codes
define the behavior or the L<expand> method.

 ------------------------------------------
 Symbol       Meaning      Nucleic Acid
 ------------------------------------------
  A            A           (A)denine
  C            C           (C)ytosine
  G            G           (G)uanine
  T            T           (T)hymine
  U            U           (U)racil
  M          A or C        a(M)ino group
  R          A or G        pu(R)ine
  W          A or T        (W)eak bond
  S          C or G        (S)trong bond
  Y          C or T        p(Y)rimidine
  K          G or T        (K)eto group
  V        A or C or G
  H        A or C or T
  D        A or G or T
  B        C or G or T
  X      G or A or T or C
  N      G or A or T or C
  .      G or A or T or C



 ------------------------------------------
 Symbol           Meaning
 ------------------------------------------
 A        Alanine
 C        Cysteine
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
 Y        Tyrosine

 B        Aspartic Acid, Asparagine
 Z        Glutamic Acid, Glutamine
 X        Any amino acid
 .        Any amino acid


=head2 Multiple Format Support

Ultimately, this module should be able to build SeqPattern.pm objects
using a variety of pattern formats such as ProSite, Blocks, Prints, GCG, etc.
Currently, this module only supports patterns using a grep-like syntax.

=head1 USAGE

A simple demo script called seq_pattern.pl is included in the examples/
directory of the central Bioperl distribution.

=head1 SEE ALSO

L<Bio::Seq> - Lightweight sequence object.

L<Bio::Tools::IUPAC> - The IUPAC code for degenerate residues and their
conversion to a regular expression.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules.  Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Steve Chervitz, sac-at-bioperl.org

=head1 COPYRIGHT

Copyright (c) 1997-8 Steve Chervitz. All Rights Reserved.
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
# CREATED : 28 Aug 1997


package Bio::Tools::SeqPattern;

use base qw(Bio::Root::Root);
use strict;
use vars qw ($ID);
$ID  = 'Bio::Tools::SeqPattern';

## These constants may be more appropriate in a Bio::Dictionary.pm
## type of class.
my $PURINES      = 'AG';
my $PYRIMIDINES  = 'CT';
my $BEE      = 'DN';
my $ZED      = 'EQ';
my $Regexp_chars = '\w,.\*()\[\]<>\{\}^\$';  # quoted for use in regexps

## Package variables used in reverse complementing.
my (%Processed_braces, %Processed_asterics);

#####################################################################################
##                                 CONSTRUCTOR                                     ##
#####################################################################################

=head1 new

 Title     : new
 Usage     : my $seqpat = Bio::Tools::SeqPattern->new();
 Purpose   : Verifies that the type is correct for superclass (Bio::Seq.pm)
           : and calls superclass constructor last.
 Returns   : n/a
 Argument  : Parameters passed to new()
 Throws    : Exception if the pattern string (seq) is empty.
 Comments  : The process of creating a new SeqPattern.pm object
           : ensures that the pattern string is untained.

See Also   : L<Bio::Root::Root::new>,
             L<Bio::Seq::_initialize>

=cut

#----------------
sub new {
#----------------
    my($class, %param) = @_;

    my $self = $class->SUPER::new(%param);
    my ($seq,$type) = $self->_rearrange([qw(SEQ TYPE)], %param);

    $seq || $self->throw("Empty pattern.");
    my $t;
    # Get the type ready for Bio::Seq.pm
    if ($type =~ /nuc|[dr]na/i) {
	$t = 'Dna';
    } elsif ($type =~ /amino|pep|prot/i) {
	$t = 'Amino';
    }
    $seq =~ tr/a-z/A-Z/;  #ps 8/8/00 Canonicalize to upper case
    $self->str($seq);
    $self->type($t);

    return $self;
}


=head1 alphabet_ok

 Title     : alphabet_ok
 Usage     : $mypat->alphabet_ok;
 Purpose   : Checks for invalid regexp characters.
           : Overrides Bio::Seq::alphabet_ok() to allow
           : additional regexp characters ,.*()[]<>{}^$
           : in addition to the standard genetic alphabet.
           : Also untaints the pattern and sets the sequence
           : object's sequence to the untained string.
 Returns   : Boolean (1 | 0)
 Argument  : n/a
 Throws    : Exception if the pattern contains invalid characters.
 Comments  : Does not call the superclass method.
           : Actually permits any alphanumeric, not just the
           : standard genetic alphabet.

=cut

#----------------'
sub alphabet_ok {
#----------------
    my( $self) = @_;

    return 1 if $self->{'_alphabet_checked'};

    $self->{'_alphabet_checked'} = 1;

    my $pat = $self->seq();

    if($pat =~ /[^$Regexp_chars]/io) {
	$self->throw("Pattern contains invalid characters: $pat",
		     'Legal characters: a-z,A-Z,0-9,,.*()[]<>{}^$ ');
    }

    # Untaint pattern (makes code taint-safe).
    $pat  =~ /([$Regexp_chars]+)/io;
    $self->setseq(uc($1));
#    print STDERR "\npattern ok: $pat\n";
    1;
}

=head1 expand

 Title     : expand
 Usage     : $seqpat_object->expand();
 Purpose   : Expands the sequence pattern using special ambiguity codes.
 Example   : $pat = $seq_pat->expand();
 Returns   : String containing fully expanded sequence pattern
 Argument  : n/a
 Throws    : Exception if sequence type is not recognized
           : (i.e., is not one of [DR]NA, Amino)

See Also   : L<Extended Alphabet Support>, L<_expand_pep>(), L<_expand_nuc>()

=cut

#----------
sub expand {
#----------
    my $self = shift;

    if($self->type =~ /[DR]na/i) { $self->_expand_nuc(); }
    elsif($self->type =~ /Amino/i) { $self->_expand_pep(); }
    else{
	$self->throw("Don't know how to expand ${\$self->type} patterns.\n");
    }
}


=head1 _expand_pep

 Title     : _expand_pep
 Usage     : n/a; automatically called by expand()
 Purpose   : Expands peptide patterns
 Returns   : String (the expanded pattern)
 Argument  : String (the unexpanded pattern)
 Throws    : n/a

See Also   : L<expand>(), L<_expand_nuc>()

=cut

#----------------
sub _expand_pep {
#----------------
    my ($self,$pat) = @_;
    $pat ||= $self->str;
    $pat =~ s/X/./g;
    $pat =~ s/^</\^/;
    $pat =~ s/>$/\$/;

    ## Avoid nested situations: [bmnq] --/--> [[$ZED]mnq]
    ## Yet correctly deal with: fze[bmnq] ---> f[$BEE]e[$ZEDmnq]
    if($pat =~ /\[\w*[BZ]\w*\]/) {
	$pat =~ s/\[(\w*)B(\w*)\]/\[$1$ZED$2\]/g;
	$pat =~ s/\[(\w*)Z(\w*)\]/\[$1$BEE$2\]/g;
	$pat =~ s/B/\[$ZED\]/g;
	$pat =~ s/Z/\[$BEE\]/g;
    } else {
	$pat =~ s/B/\[$ZED\]/g;
	$pat =~ s/Z/\[$BEE\]/g;
    }
    $pat =~ s/\((.)\)/$1/g;  ## Doing these last since:
    $pat =~ s/\[(.)\]/$1/g;  ## Pattern could contain [B] (for example)

    return $pat;
}



=head1 _expand_nuc

 Title     : _expand_nuc
 Purpose   : Expands nucleotide patterns
 Returns   : String (the expanded pattern)
 Argument  : String (the unexpanded pattern)
 Throws    : n/a

See Also   : L<expand>(), L<_expand_pep>()

=cut

#---------------
sub _expand_nuc {
#---------------
    my ($self,$pat) = @_;

    $pat ||= $self->str;
    $pat =~ s/N|X/./g;
    $pat =~ s/pu/R/ig;
    $pat =~ s/py/Y/ig;
    $pat =~ s/U/T/g;
    $pat =~ s/^</\^/;
    $pat =~ s/>$/\$/;

    ## Avoid nested situations: [ya] --/--> [[ct]a]
    ## Yet correctly deal with: sg[ya] ---> [gc]g[cta]
    if($pat =~ /\[\w*[RYSWMK]\w*\]/) {
	$pat =~ s/\[(\w*)R(\w*)\]/\[$1$PURINES$2\]/g;
	$pat =~ s/\[(\w*)Y(\w*)\]/\[$1$PYRIMIDINES$2\]/g;
	$pat =~ s/\[(\w*)S(\w*)\]/\[$1GC$2\]/g;
	$pat =~ s/\[(\w*)W(\w*)\]/\[$1AT$2\]/g;
	$pat =~ s/\[(\w*)M(\w*)\]/\[$1AC$2\]/g;
	$pat =~ s/\[(\w*)K(\w*)\]/\[$1GT$2\]/g;
	$pat =~ s/\[(\w*)V(\w*)\]/\[$1ACG$2\]/g;
	$pat =~ s/\[(\w*)H(\w*)\]/\[$1ACT$2\]/g;
	$pat =~ s/\[(\w*)D(\w*)\]/\[$1AGT$2\]/g;
	$pat =~ s/\[(\w*)B(\w*)\]/\[$1CGT$2\]/g;
	$pat =~ s/R/\[$PURINES\]/g;
	$pat =~ s/Y/\[$PYRIMIDINES\]/g;
	$pat =~ s/S/\[GC\]/g;
	$pat =~ s/W/\[AT\]/g;
	$pat =~ s/M/\[AC\]/g;
	$pat =~ s/K/\[GT\]/g;
	$pat =~ s/V/\[ACG\]/g;
	$pat =~ s/H/\[ACT\]/g;
	$pat =~ s/D/\[AGT\]/g;
	$pat =~ s/B/\[CGT\]/g;
    } else {
	$pat =~ s/R/\[$PURINES\]/g;
	$pat =~ s/Y/\[$PYRIMIDINES\]/g;
	$pat =~ s/S/\[GC\]/g;
	$pat =~ s/W/\[AT\]/g;
	$pat =~ s/M/\[AC\]/g;
	$pat =~ s/K/\[GT\]/g;
	$pat =~ s/V/\[ACG\]/g;
	$pat =~ s/H/\[ACT\]/g;
	$pat =~ s/D/\[AGT\]/g;
	$pat =~ s/B/\[CGT\]/g;
    }
    $pat =~ s/\((.)\)/$1/g;  ## Doing thses last since:
    $pat =~ s/\[(.)\]/$1/g;  ## Pattern could contain [y] (for example)

    return $pat;
}



=head1 revcom

 Title     : revcom
 Usage     : revcom([1]);
 Purpose   : Forms a pattern capable of recognizing the reverse complement
           : version of a nucleotide sequence pattern.
 Example   : $pattern_object->revcom();
           : $pattern_object->revcom(1); ## returns expanded rev complement pattern.
 Returns   : Object reference for a new Bio::Tools::SeqPattern containing
           : the revcom of the current pattern as its sequence.
 Argument  : (1) boolean (optional) (default= false)
           :     true : expand the pattern before rev-complementing.
           :     false: don't expand pattern before or after rev-complementing.
 Throws    : Exception if called for amino acid sequence pattern.
 Comments  : This method permits the simultaneous searching of both
           : sense and anti-sense versions of a nucleotide pattern
           : by means of a grep-type of functionality in which any
           : number of patterns may be or-ed into the recognition
           : pattern.
           : Overrides Bio::Seq::revcom() and calls it first thing.
           : The order of _fixpat() calls is critical.

See Also   : L<Bio::Seq::revcom>, L</_fixpat_1>, L</_fixpat_2>, L</_fixpat_3>, L</_fixpat_4>, L</_fixpat_5>

=cut

#-----------'
sub revcom {
#-----------
    my($self,$expand) = @_;

    if ($self->type !~ /Dna|Rna/i) {
	$self->throw("Can't get revcom for ${\$self->type} sequence types.\n");
    }
#    return $self->{'_rev'} if defined $self->{'_rev'};

    $expand ||= 0;
    my $str = $self->str;
    $str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    my $rev = CORE::reverse $str;
    $rev    =~ tr/[](){}<>/][)(}{></;

    if($expand) {
	$rev = $self->_expand_nuc($rev);
#	print "\nExpanded: $rev\n";
    }

    %Processed_braces = ();
    %Processed_asterics = ();

    my $fixrev = _fixpat_1($rev);
#   print "FIX 1: $fixrev";<STDIN>;

     $fixrev = _fixpat_2($fixrev);
#   print "FIX 2: $fixrev";<STDIN>;

     $fixrev = _fixpat_3($fixrev);
#    print "FIX 3: $fixrev";<STDIN>;

     $fixrev = _fixpat_4($fixrev);
#    print "FIX 4: $fixrev";<STDIN>;

     $fixrev = _fixpat_5($fixrev);
#    print "FIX 5: $fixrev";<STDIN>;

##### Added by ps 8/7/00 to allow non-greedy matching
     $fixrev = _fixpat_6($fixrev);
#    print "FIX 6: $fixrev";<STDIN>;

#    $self->{'_rev'} = $fixrev;

     return Bio::Tools::SeqPattern->new(-seq =>$fixrev, -type =>$self->type);
}

=head1 backtranslate

 Title     : backtranslate
 Usage     : backtranslate();
 Purpose   : Produce a degenerate oligonucleotide whose translation would produce
           : the original protein motif.
 Example   : $pattern_object->backtranslate();
 Returns   : Object reference for a new Bio::Tools::SeqPattern containing
           : the reverse translation of the current pattern as its sequence.
 Throws    : Exception if called for nucleotide sequence pattern.

=cut

sub backtranslate {
    my $self = shift;
    
    # _load_module loads dynamically, caches call if successful
    $self->_load_module('Bio::Tools::SeqPattern::Backtranslate');
    Bio::Tools::SeqPattern::Backtranslate->import("_reverse_translate_motif");

    if ($self->type ne 'Amino') {
        $self->throw(
            "Can't get backtranslate for ${\$self->type} sequence types.\n"
        );
    }

    return __PACKAGE__->new(
        -SEQ  => _reverse_translate_motif($self->str),
        -TYPE => 'Dna',
    );
}

=head1 _fixpat_1

 Title     : _fixpat_1
 Usage     : n/a; called automatically by revcom()
 Purpose   : Utility method for revcom()
           : Converts all {7,5} --> {5,7}     (Part I)
           :           and [T^] --> [^T]      (Part II)
           :           and *N   --> N*        (Part III)
 Returns   : String (the new, partially reversed pattern)
 Argument  : String (the expanded pattern)
 Throws    : n/a

See Also   : L<revcom>()

=cut

#--------------
sub _fixpat_1 {
#--------------
    my $pat = shift;

    ## Part I:
    my (@done,@parts);
    while(1) {
	$pat =~ /(.*)\{(\S+?)\}(.*)/ or do{ push @done, $pat; last; };
	$pat = $1.'#{'.reverse($2).'}'.$3;
#	print "1: $1\n2: $2\n3: $3\n";
#	print "modified pat: $pat";<STDIN>;
	@parts = split '#', $pat;
	push @done, $parts[1];
	$pat = $parts[0];
#	print "done: $parts[1]<---\nnew pat: $pat<---";<STDIN>;
	last if not $pat;
    }
    $pat = join('', reverse @done);

    ## Part II:
    @done = ();
    while(1) {
	$pat =~ /(.*)\[(\S+?)\](.*)/ or do{ push @done, $pat; last; };
	$pat = $1.'#['.reverse($2).']'.$3;
#	print "1: $1\n2: $2\n3: $3\n";
#	print "modified pat: $pat";<STDIN>;
	@parts = split '#', $pat;
	push @done, $parts[1];
	$pat = $parts[0];
#	print "done: $parts[1]<---\nnew pat: $pat<---";<STDIN>;
	last if not $pat;
    }
    $pat = join('', reverse @done);

    ## Part III:
    @done = ();
    while(1) {
	$pat =~ /(.*)\*([\w.])(.*)/ or do{ push @done, $pat; last; };
	$pat = $1.'#'.$2.'*'.$3;
	$Processed_asterics{$2}++;
#	print "1: $1\n2: $2\n3: $3\n";
#	print "modified pat: $pat";<STDIN>;
	@parts = split '#', $pat;
	push @done, $parts[1];
	$pat = $parts[0];
#	print "done: $parts[1]<---\nnew pat: $pat<---";<STDIN>;
	last if not $pat;
    }
    return join('', reverse @done);

}


=head1 _fixpat_2

 Title     : _fixpat_2
 Usage     : n/a; called automatically by revcom()
 Purpose   : Utility method for revcom()
           : Converts all {5,7}Y ---> Y{5,7}
           :          and {10,}. ---> .{10,}
 Returns   : String (the new, partially reversed pattern)
 Argument  : String (the expanded, partially reversed pattern)
 Throws    : n/a

See Also   : L<revcom>()

=cut

#--------------
sub _fixpat_2 {
#--------------
    my $pat = shift;

    local($^W) = 0;
    my (@done,@parts,$braces);
    while(1) {
#	$pat =~ s/(.*)([^])])(\{\S+?\})([\w.])(.*)/$1$2#$4$3$5/ or do{ push @done, $pat; last; };
	$pat =~ s/(.*)(\{\S+?\})([\w.])(.*)/$1#$3$2$4/ or do{ push @done, $pat; last; };
	$braces = $2;
	$braces =~ s/[{}]//g;
	$Processed_braces{"$3$braces"}++;
#	print "modified pat: $pat";<STDIN>;
	@parts = split '#', $pat;
	push @done, $parts[1];
	$pat = $parts[0];
#	print "done: $parts[1]<---\nnew pat: $pat<---";<STDIN>;
	last if not $pat;
    }
    return join('', reverse @done);
}


=head1 _fixpat_3

 Title     : _fixpat_3
 Usage     : n/a; called automatically by revcom()
 Purpose   : Utility method for revcom()
           : Converts all {5,7}(XXX) ---> (XXX){5,7}
 Returns   : String (the new, partially reversed pattern)
 Argument  : String (the expanded, partially reversed pattern)
 Throws    : n/a

See Also   : L<revcom>()

=cut

#-------------
sub _fixpat_3 {
#-------------
    my $pat = shift;

    my (@done,@parts,$braces,$newpat,$oldpat);
    while(1) {
#	$pat =~ s/(.+)(\{\S+\})(\(\w+\))(.*)/$1#$3$2$4/ or do{ push @done, $pat; last; };
	if( $pat =~ /(.*)(.)(\{\S+\})(\(\w+\))(.*)/) {
	    $newpat = "$1#$2$4$3$5";
##ps	    $oldpat = "$1#$2$3$4$5";
#	    print "1: $1\n2: $2\n3: $3\n4: $4\n5: $5\n";
##ps	    $braces = $3;
##ps	    $braces =~ s/[{}]//g;
##ps	    if( exists $Processed_braces{"$2$braces"} || exists $Processed_asterics{$2}) {
##ps		$pat = $oldpat;  # Don't change it. Already processed.
#		print "saved pat: $pat";<STDIN>;
##ps	    } else {
#		print "new pat: $newpat";<STDIN>;
		$pat = $newpat;  # Change it.
##ps	    }
	} elsif( $pat =~ /^(\{\S+\})(\(\w+\))(.*)/) {
	    $pat = "#$2$1$3";
	} else {
	    push @done, $pat; last;
	}
	@parts = split '#', $pat;
	push @done, $parts[1];
	$pat = $parts[0];
#	print "done: $parts[1]<---\nnew pat: $pat<---";<STDIN>;
	last if not $pat;
    }
    return join('', reverse @done);
}


=head1 _fixpat_4

 Title     : _fixpat_4
 Usage     : n/a; called automatically by revcom()
 Purpose   : Utility method for revcom()
           : Converts all {5,7}[XXX] ---> [XXX]{5,7}
 Returns   : String (the new, partially reversed pattern)
 Argument  : String (the expanded, partially reversed  pattern)
 Throws    : n/a

See Also   : L<revcom>()

=cut

#---------------
sub _fixpat_4 {
#---------------
    my $pat = shift;

    my (@done,@parts,$braces,$newpat,$oldpat);
    while(1) {
#	$pat =~ s/(.*)(\{\S+\})(\[\w+\])(.*)/$1#$3$2$4/ or do{ push @done, $pat; last; };
#	$pat =~ s/(.*)([^\w.])(\{\S+\})(\[\w+\])(.*)/$1$2#$4$3$5/ or do{ push @done, $pat; last; };
	if( $pat =~ /(.*)(.)(\{\S+\})(\[\w+\])(.*)/) {
	    $newpat = "$1#$2$4$3$5";
	    $oldpat = "$1#$2$3$4$5";
#	    print "1: $1\n2: $2\n3: $3\n4: $4\n5: $5\n";
	    $braces = $3;
	    $braces =~ s/[{}]//g;
	    if( (defined $braces and defined $2) and
		exists $Processed_braces{"$2$braces"} || exists $Processed_asterics{$2}) {
		$pat = $oldpat;  # Don't change it. Already processed.
#		print "saved pat: $pat";<STDIN>;
	    } else {
		$pat = $newpat;  # Change it.
#		print "new pat: $pat";<STDIN>;
	    }
	} elsif( $pat =~ /^(\{\S+\})(\[\w+\])(.*)/) {
	    $pat = "#$2$1$3";
	} else {
	    push @done, $pat; last;
	}

	@parts = split '#', $pat;
	push @done, $parts[1];
	$pat = $parts[0];
#	print "done: $parts[1]<---\nnew pat: $pat<---";<STDIN>;
	last if not $pat;
    }
    return join('', reverse @done);
}


=head1 _fixpat_5

 Title     : _fixpat_5
 Usage     : n/a; called automatically by revcom()
 Purpose   : Utility method for revcom()
           : Converts all *[XXX]  ---> [XXX]*
           :          and *(XXX)  ---> (XXX)*
 Returns   : String (the new, partially reversed pattern)
 Argument  : String (the expanded, partially reversed pattern)
 Throws    : n/a

See Also   : L<revcom>()

=cut

#--------------
sub _fixpat_5 {
#--------------
    my $pat = shift;

    my (@done,@parts,$newpat,$oldpat);
    while(1) {
#	$pat =~ s/(.*)(\{\S+\})(\[\w+\])(.*)/$1#$3$2$4/ or do{ push @done, $pat; last; };
#	$pat =~ s/(.*)([^\w.])(\{\S+\})(\[\w+\])(.*)/$1$2#$4$3$5/ or do{ push @done, $pat; last; };
	if( $pat =~ /(.*)(.)\*(\[\w+\]|\(\w+\))(.*)/) {
	    $newpat = "$1#$2$3*$4";
	    $oldpat = "$1#$2*$3$4";
#	    print "1: $1\n2: $2\n3: $3\n4: $4\n";
	    if( exists $Processed_asterics{$2}) {
		$pat = $oldpat;  # Don't change it. Already processed.
#		print "saved pat: $pat";<STDIN>;
	    } else {
		$pat = $newpat;  # Change it.
#		print "new pat: $pat";<STDIN>;
	    }
	} elsif( $pat =~ /^\*(\[\w+\]|\(\w+\))(.*)/) {
	    $pat = "#$1*$3";
	} else {
	    push @done, $pat; last;
	}

	@parts = split '#', $pat;
	push @done, $parts[1];
	$pat = $parts[0];
#	print "done: $parts[1]<---\nnew pat: $pat<---";<STDIN>;
	last if not $pat;
    }
    return join('', reverse @done);
}





############################
#
#  PS: Added 8/7/00 to allow non-greedy matching patterns
#
######################################

=head1 _fixpat_6

 Title     : _fixpat_6
 Usage     : n/a; called automatically by revcom()
 Purpose   : Utility method for revcom()
           : Converts all ?Y{5,7}  ---> Y{5,7}?
           :          and ?(XXX){5,7}  ---> (XXX){5,7}?
           :          and ?[XYZ]{5,7}  ---> [XYZ]{5,7}?
 Returns   : String (the new, partially reversed pattern)
 Argument  : String (the expanded, partially reversed pattern)
 Throws    : n/a

See Also   : L<revcom>()

=cut

#--------------
sub _fixpat_6 {
#--------------
    my $pat = shift;
    my (@done,@parts);

   @done = ();
    while(1) {
	$pat =~   /(.*)\?(\[\w+\]|\(\w+\)|\w)(\{\S+?\})?(.*)/ or do{ push @done, $pat; last; };
     my $quantifier = $3 ? $3 : ""; # Shut up warning if no explicit quantifier
 	$pat = $1.'#'.$2.$quantifier.'?'.$4;
#	$pat = $1.'#'.$2.$3.'?'.$4;

#	print "1: $1\n2: $2\n3: $3\n";
#	print "modified pat: $pat";<STDIN>;
	@parts = split '#', $pat;
	push @done, $parts[1];
	$pat = $parts[0];
#	print "done: $parts[1]<---\nnew pat: $pat<---";<STDIN>;
	last if not $pat;
    }
    return join('', reverse @done);

 }

=head2 str

 Title   : str
 Usage   : $obj->str($newval)
 Function:
 Returns : value of str
 Args    : newvalue (optional)


=cut

sub str{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'str'} = $value;
    }
    return $obj->{'str'};

}

=head2 type

 Title   : type
 Usage   : $obj->type($newval)
 Function:
 Returns : value of type
 Args    : newvalue (optional)


=cut

sub type{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'type'} = $value;
    }
    return $obj->{'type'};

}

1;

__END__

#########################################################################
#  End of class
#########################################################################

=head1 FOR DEVELOPERS ONLY

=head2 Data Members

Information about the various data members of this module is provided
for those wishing to modify or understand the code. Two things to bear
in mind:

=over 2

=item 1 Do NOT rely on these in any code outside of this module.

All data members are prefixed with an underscore to signify that they
are private.  Always use accessor methods. If the accessor doesn't
exist or is inadequate, create or modify an accessor (and let me know,
too!).

=item 2 This documentation may be incomplete and out of date.

It is easy for this documentation to become obsolete as this module is
still evolving.  Always double check this info and search for members
not described here.

=back

An instance of Bio::Tools::RestrictionEnzyme.pm is a blessed reference
to a hash containing all or some of the following fields:

 FIELD          VALUE
 ------------------------------------------------------------------------
 _rev     : The corrected reverse complement of the fully expanded pattern.

 INHERITED DATA MEMBERS:

 _seq     : (From Bio::Seq.pm) The original, unexpanded input sequence after untainting.
 _type    : (From Bio::Seq.pm) 'Dna' or 'Amino'


=cut
