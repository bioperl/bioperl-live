#-----------------------------------------------------------------------------
# PACKAGE : Bio::Tools::SeqPattern.pm
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED : 28 Aug 1997
# REVISION: $Id$
#            
# Copyright (c) 1997-8 Steve A. Chervitz. All Rights Reserved.
#           This module is free software; you can redistribute it and/or 
#           modify it under the same terms as Perl itself.
#-----------------------------------------------------------------------------

package Bio::Tools::SeqPattern;

use Bio::Root::Global qw(:devel);
use Bio::Seq ();

@ISA = qw(Bio::Seq);
use strict;
use vars qw ($ID $VERSION);
$ID  = 'Bio::Tools::SeqPattern';
$VERSION = 0.011;

## These constants may be more appropriate in a Bio::Dictionary.pm 
## type of class.
my $PURINES      = 'AG';
my $PYRIMIDINES  = 'CT';
my $PHILICS      = 'TSHEDQNKR';
my $PHOBICS      = 'IFVLWMAGCY';
my $Regexp_chars = '\w,.\*()\[\]<>\{\}^\$';  # quoted for use in regexps

## Package variables used in reverse complementing.
my (%Processed_braces, %Processed_asterics);

## POD Documentation:

=head1 NAME

Bio::Tools::SeqPattern.pm - Bioperl object for a sequence pattern or motif

=head1 SYNOPSIS

=head2 Object Creation

    use Bio::Tools::SeqPattern ();

    $pat1     = 'T[GA]AA...TAAT';
    $pattern1 = new Bio::Tools::SeqPattern(-SEQ =>$pattern, -TYPE =>'Dna'); 

    $pat2     = '[VILM]R(GXX){3,2}...[^PG]';
    $pattern2 = new Bio::Tools::SeqPattern(-SEQ =>$pattern, -TYPE =>'Amino'); 

=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

=head1 DESCRIPTION

The Bio::Tools::SeqPattern.pm module encapsulates generic data and methods for manipulating 
regular expressions describing nucleic or amino acid sequence patterns (a.k.a, "motifs").

Bio::Tools::SeqPattern.pm is a concrete class that inherits from B<Bio::Seq.pm>.

This class grew out of a need to have a standard module for doing routine
tasks with sequence patterns such as:

  -- Forming a reverse-complement version of a nucleotide sequence pattern
  -- Expanding patterns containing ambiguity codes
  -- Checking for invalid regexp characters
  -- Untainting yet preserving special characters in the pattern

Other features to look for in the future:
 
  -- Full pattern syntax checking
  -- Conversion between expanded and ondensed forms of the pattern

=head1 MOTIVATIONS

A key motivation for Bio::Tools::SeqPattern.pm is to have a way to
generate a reverse complement of a nucleotide sequence pattern.
This makes possible simultaneous pattern matching on both sense and 
anti-sense strands of a query sequence. 
 
In principle, one could do such a search more inefficiently by testing 
ainst both sense and anti-sense versions of a sequence. 
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
send me some email (sac@genome.stanford.edu). Thanks.

=head1 OTHER FEATURES

=head2 Extended Alphabet Support

This module supports the same set of ambiguity codes for nucleotide 
sequences as supported by B<Bio::Seq.pm>. These ambiguity codes
define the behavior or the expand() method.
Amino acid alphabet support is different from that of Seq.pm (see below).

 ------------------------------------------
 Symbol       Meaning      Nucleic Acid
 ------------------------------------------
  A            A           Adenine
  C            C           Cytosine
  G            G           Guanine
  T            T           Thymine
  U            U           Uracil
  M          A or C  
  R          A or G        Any purine
  W          A or T    
  S          C or G     
  Y          C or T        Any pyrimidine
  K          G or T     
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

 B        Any hydrophobic: IFVLWMAGCY
 Z        Any hydrophilic: TSHEDQNKR
 X        Any amino acid
 .        Any amino acid


=head2   Multiple Format Support

Ultimately, this module should be able to build SeqPattern.pm objects objects 
using a variety of pattern formats such as ProSite, Blocks, Prints, GCG, etc.
Currently, this module only supports patterns using a grep-like syntax. 

=head1 USAGE

A simple demo script is included with the central Bioperl distribution (L<INSTALLATION>)
and is also available from:

    http://bio.perl.org/Core/Examples/seq_pattern.pl

=head1 SEE ALSO

  Bio::Root::Object.pm    - Base class.
  Bio::Seq.pm             - Lightweight sequence object.

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

See the L<FEEDBACK> section for where to send bug reports and comments.

=head1 VERSION

Bio::Tools::SeqPattern.pm, 0.011

=head1 COPYRIGHT

Copyright (c) 1997-8 Steve A. Chervitz. All Rights Reserved.
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



#####################################################################################
##                                 CONSTRUCTOR                                     ##
#####################################################################################


=head1 _initialize

 Title     : _initialize
 Usage     : n/a; automatically called by Bio::Root::Object::new()
 Purpose   : Verifies that the type is correct for superclass (Bio::Seq.pm)
           : and calls superclass constructor last.
 Returns   : n/a
 Argument  : Parameters passed to new()
 Throws    : Exception if the pattern string (seq) is empty.
 Comments  : The process of creating a new SeqPattern.pm object
           : ensures that the pattern string is untained.

See Also   : L<_untaint_pat>(), B<Bio::Root::Object::new()>, B<Bio::Seq::_initialize()>

=cut

#----------------
sub _initialize {
#----------------
    my($self, %param) = @_;
    
    my($seq, $type) = $self->_rearrange([qw(SEQ TYPE)], %param);

    $seq || $self->throw("Empty pattern.");

    # Get the type ready for Bio::Seq.pm
    if ($type =~ /nuc|[dr]na/i) {
	$param{-TYPE} = 'Dna';
    } elsif ($type =~ /amino|pep|prot/i) {
	$param{-TYPE} = 'Amino';
    }

    my $make = $self->SUPER::_initialize(%param);
#	 if($self->make() =~ /^nuc/i) {
#	     $param{-SEQ} = $self->_expand_nuc($param{-SEQ}); 
#	 } else {
##	    $param{-SEQ} = $self->_expand_pep($param{-SEQ});  #unimplemented here.
#	 }

    $make;
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
           : Actaully permits any alphanumeric, not just the
           : standard genetic alphabet.

See Also   : B<Bio::Seq::alphabet_ok()>, L<_initialize>()

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
    $pat  =~ /[$Regexp_chars]+/io; 
    $self->setseq(uc($&));
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

See Also   : B<Extended Alphabet Support>, L<_expand_pep>(), L<_exapand_nuc>()

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

    ## Avoid nested situations: [bmnq] --/--> [[$PHOBICS]mnq]
    ## Yet correctly deal with: fze[bmnq] ---> f[$PHILICS]e[$PHOBICSmnq]
    if($pat =~ /\[\w*[BZ]\w*\]/) {
	$pat =~ s/\[(\w*)B(\w*)\]/\[$1$PHOBICS$2\]/g;
	$pat =~ s/\[(\w*)Z(\w*)\]/\[$1$PHILICS$2\]/g;
	$pat =~ s/B/\[$PHOBICS\]/g;
	$pat =~ s/Z/\[$PHILICS\]/g;
    } else {
	$pat =~ s/B/\[$PHOBICS\]/g;
	$pat =~ s/Z/\[$PHILICS\]/g;
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

See Also   : B<Bio::Seq::revcom()>, L<_fixpat_1>(), L<_fixpat_2>(), L<_fixpat_3>(), L<_fixpat_4>(), L<_fixpat_5>()

=cut
#'

#-----------
sub revcom {
#-----------
    my($self,$expand) = @_;
    
    if ($self->type !~ /Dna|Rna/i) {
	$self->throw("Can't get revcom for ${\$self->type} sequence types.\n");
    }
#    return $self->{'_rev'} if defined $self->{'_rev'};

    $expand ||= 0;
    my $rev = $self->SUPER::revcom->str;  ## Invoke Bio::Seq::revcom()
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
    
#    $self->{'_rev'} = $fixrev;

     return new Bio::Tools::SeqPattern(-seq =>$fixrev, -type =>$self->type);
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
	    $oldpat = "$1#$2$3$4$5";
#	    print "1: $1\n2: $2\n3: $3\n4: $4\n5: $5\n";
	    $braces = $3;
	    $braces =~ s/[{}]//g;
	    if( exists $Processed_braces{"$2$braces"} || exists $Processed_asterics{$2}) {
		$pat = $oldpat;  # Don't change it. Already processed.
#		print "saved pat: $pat";<STDIN>;
	    } else {
#		print "new pat: $newpat";<STDIN>;
		$pat = $newpat;  # Change it.
	    }
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

1;
__END__

#########################################################################
#  End of class 
#########################################################################

=head1 FOR DEVELOPERS ONLY

=head2 Data Members

Information about the various data members of this module is provided for those 
wishing to modify or understand the code. Two things to bear in mind: 

=over 4

=item 1 Do NOT rely on these in any code outside of this module. 

All data members are prefixed with an underscore to signify that they are private.
Always use accessor methods. If the accessor doesn't exist or is inadequate, 
create or modify an accessor (and let me know, too!). 

=item 2 This documentation may be incomplete and out of date.

It is easy for this documentation to become obsolete as this module is still evolving. 
Always double check this info and search for members not described here.

=back

An instance of Bio::Tools::RestrictionEnzyme.pm is a blessed reference to a hash
containing all or some of the following fields:

 FIELD          VALUE
 ------------------------------------------------------------------------
 _rev     : The corrected reverse complement of the fully expanded pattern.

 INHERITED DATA MEMBERS:

 _seq     : (From Bio::Seq.pm) The original, unexpanded input sequence after untainting.
 _type    : (From Bio::Seq.pm) 'Dna' or 'Amino' 


=cut



