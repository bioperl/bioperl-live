#-----------------------------------------------------------------------------
# PACKAGE : Bio::Tools::Sigcleave.pm
# AUTHOR  : Chris Dagdigian, dag@sonsorol.org
# CREATED : Jan 28 1999
# REVISION: $Id$
#            
# Copyright (c) 1997-9 Chris Dagdigian and others. All Rights Reserved.
#           This module is free software; you can redistribute it and/or 
#           modify it under the same terms as Perl itself.
#
# _History_
#
# Object framework ripped from Steve Chervits's SeqPattern.pm
# 
# Core EGCG Sigcleave emulation from perl code developed by
# Danh Nguyen & Kamalakar Gulukota which itself was based 
# loosely on Colgrove's signal.c program.
#
# The overall idea is to replicate the output of the sigcleave
# program which was distributed with the EGCG extension to the GCG sequence
# analysis package. There is also an accessor method for just getting at
# the raw results.
#
#-----------------------------------------------------------------------------

package Bio::Tools::Sigcleave;

use Bio::Root::Global qw(:devel);
use Bio::Seq ();

@ISA = qw(Bio::Seq);
use strict;
use vars qw ($ID $VERSION %WeightTable);
$ID  = 'Bio::Tools::Sigcleave';
$VERSION = 0.01;

=head1 NAME

Bio::Tools::Sigcleave.pm - Bioperl object for sigcleave analysis

=head1 SYNOPSIS

=head2 Object Creation

    use Bio::Tools::Sigcleave ();

    $sigcleave_object = new Bio::Tools::Sigcleave(-file=>'sigtest.aa',
                                                  -desc=>'test sigcleave protein seq',
                                                  -type=>'AMINO',
                                                  -threshold=>'3.5',
                                                 );

Sigcleave objects can be created via the same methods as Bio::Seq objects. The
one additional parameter is "-threshold" which sets the score reporting limit
for the algorithim. The above exmple shows a sigcleave object being created
from a protein sequence file. See the Bio::Seq documention to see the other ways
that objects can be created.

=head2 Object Methods & Accessors

     %raw_results      = $sigcleave_object->signals;

     $formatted_output = $sigcleave_object->pretty_print;

=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

=head1 DESCRIPTION

"Sigcleave" was a program distributed as part of the free EGCG add-on to
earlier versions of the GCG Sequence Analysis package. 

From the EGCG documentation:
  SigCleave uses the von Heijne method to locate signal sequences, and to identify 
  the cleavage site. The method is 95% accurate in resolving signal sequences from 
  non-signal sequences with a cutoff score of 3.5, and 75-80% accurate in identifying 
  the cleavage site. The program reports all hits above a minimum value. 

The EGCG Sigcleave program was written by Peter Rice 
(E-mail: pmr@sanger.ac.uk Post: Informatics Division, The Sanger Centre,
Wellcome Trust Genome Campus, Hinxton, Cambs, CB10 1SA, UK). 

Since EGCG is no longer distributed for the latest versions of GCG, this code
was developed to emulate the output of the original program as much as possible for
those who lost access to sigcleave when upgrading to newer versions of GCG.

The EGCG website can be found at: http://www.sanger.ac.uk/Software/EGCG/

There are 2 accessor methods for this object. "signals" will return a perl
associative array containing the sigcleave scores keyed by amino acid position.
"pretty_print" returns a formatted string similar to the output of the original
sigcleave utility.

In both cases, the "threshold" setting controls the score reporting level. If no
value for threshold is passed in by the user, the code defaults to a reporting value
of 3.5. 

Example pretty_print output:

	SIGCLEAVE of sigtest from: 1 to 146

	Report scores over 3.5
	Maximum score 4.9 at residue 131

	 Sequence:  FVILAAMSIQGSA-NLQTQWKSTASLALET
        	    | (signal)    | (mature peptide)
          	118            131

	 Other entries above 3.5

	Maximum score 3.7 at residue 112

	 Sequence:  CSRQLFGWLFCKV-HPGAIVFVILAAMSIQGSANLQTQWKSTASLALET
         	   | (signal)    | (mature peptide)
           	99            112


=head1 USAGE

No warranty implied or expressed. Use at your own risk :) Users unfamiliar
with the original Sigcleave application should read the von Heijne papers. 

The weight matrix in thos code is for eukaryote signal sequences.

Please see the example script located in the bioperl distribution
to see how this code can be used.

=head1 FEEDBACK

When updating and maintaining a module, it helps to know that people
are actually using it. Let us know if you find a bug, think this code
is useful or have any improvements/features to suggest.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the bugs and 
their resolution. Bug reports can be submitted via email or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR

Chris Dagdigian, dag@sonsorol.org  & others

=head1 VERSION

Bio::Tools::Sigcleave.pm, $Id$

=head1 COPYRIGHT

Copyright (c) 1999 Chris Dagdigian & others. All Rights Reserved.
This module is free software; you can redistribute it and/or 
modify it under the same terms as Perl itself.

=head1 REFERENCES / SEE ALSO

von Heijne G. (1986) "A new method for predicting signal sequences cleavage sites." 
Nucleic Acids Res. 14, 4683-4690. 

von Heijne G. (1987) in "Sequence Analysis in Molecular Biology: Treasure Trove or Trivial Pursuit" 
(Acad. Press, (1987), 113-117). 

EGCG website: http://www.sanger.ac.uk/Software/EGCG/

=head1 APPENDIX

The following documentation describes the various functions
contained in this module. Some functions are for internal 
use and are not meant to be called by the user; they are 
preceded by an underscore ("_").


=cut

#
##
###
#### END of main POD documentation.
###
##
#


%Bio::Tools::Sigcleave::WeightTable = (
A => [ -0.109199, -0.035091, 0.033902, 0.321584, 0.216223, 0.216223, 0.159065, 0.544727, 0.033902, 1.175999, -0.882389, 1.707878, 0.216223, -0.882389, 0.000000, 0.000000], 
C => [ 0.287682, 0.693147, 0.441833, 0.693147, 1.134980, 0.287682, 0.575364, 0.105361, 0.287682, 1.440362, -0.405465, 0.693147, 0.575364, -0.405465, 0.000000, 0.000000], 
D => [ -2.186051, -2.186051, -2.186051, -2.186051, -2.186051, -2.186051, -2.186051, -0.576613, -1.087439, -25.211902, -0.576613, -25.211902, 0.116534, 0.211844, 0.000000, 0.000000], 
E => [ -2.302585, -2.302585, -2.302585, -2.302585, -2.302585, -2.302585, -2.302585, -1.203973, -0.356675, -25.328436, -0.356675, -25.328436, 0.262364, 0.336472, 0.000000, 0.000000], 
F => [ 0.474458, 0.675129, 0.675129, 0.068993, 0.223144, 1.167605, 0.842183, -0.336472, -0.113329, -24.748618, 0.842183, -24.748618, 0.068993, -0.336472, 0.000000, 0.000000], 
G => [ -1.106911, -1.394593, -0.701446, -1.394593, 0.071744, -1.394593, -1.800058, 0.451234, 1.033155, -0.883768, -0.547295, 1.170356, -0.190620, -0.547295, 0.000000, 0.000000], 
H => [ -1.223775, -1.223775, -1.223775, -1.223775, -1.223775, -1.223775, -1.223775, 0.385662, -1.223775, -24.249626, 0.567984, -24.249626, 0.162519, -0.530628, 0.000000, 0.000000], 
I => [ 0.706570, 0.077962, -0.209721, 0.396415, -0.392042, -0.615186, 0.077962, -0.392042, -2.001480, 0.301105, -0.392042, -25.027331, 0.077962, -0.055570, 0.000000, 0.000000], 
K => [ -2.424803, -2.424803, -2.424803, -2.424803, -2.424803, -2.424803, -2.424803, -2.424803, -1.038508, -25.450654, -1.731656, -25.450654, -0.026907, -0.227578, 0.000000, 0.000000], 
L => [ 1.726302, 1.783461, 1.876242, 1.863503, 1.313457, 1.665678, 1.398615, -0.190620, 0.642289, -0.413764, 0.502527, -2.493205, -0.413764, -1.106911, 0.000000, 0.000000], 
M => [ 0.105361, 0.952658, 0.393043, -0.993252, 0.798508, -0.300105, -0.300105, -0.993252, -0.993252, -24.019103, -0.993252, -24.019103, -0.993252, -0.300105, 0.000000, 0.000000], 
N => [ -1.960095, -1.960095, -1.960095, -1.960095, -1.960095, -1.960095, -1.960095, -0.861482, -0.861482, -24.985946, 0.342490, -24.985946, -0.573800, -0.014185, 0.000000, 0.000000], 
P => [ -2.001480, -1.308333, -2.001480, -2.001480, -0.615186, -2.001480, 0.077962, 0.994252, 0.637577, -25.027331, -2.001480, -0.902868, -2.001480, 1.089562, 0.000000, 0.000000], 
Q => [ -1.840550, -1.840550, -1.840550, -1.840550, -0.048790, -1.840550, -1.840550, 0.462035, 0.238892, -24.866401, 1.049822, -0.741937, 1.103889, 0.462035, 0.000000, 0.000000], 
R => [ -2.028148, -2.028148, -2.028148, -2.028148, -2.028148, -2.028148, -2.028148, -0.082238, -0.641854, -25.053999, 0.679902, -25.053999, 0.456758, 0.169076, 0.000000, 0.000000], 
S => [ -1.335001, -0.354172, -0.641854, 0.131336, -0.131028, 0.274437, 0.338975, 0.824483, -0.035718, 0.701881, 0.399600, 0.562119, 0.274437, -0.131028, 0.000000, 0.000000], 
T => [ 0.030459, -0.662688, -0.885832, -0.662688, 0.292823, -0.326216, -0.326216, 0.212781, -0.480366, 0.561087, -0.192684, -0.480366, -1.173514, 0.030459, 0.000000, 0.000000], 
V => [ 0.811931, 0.301105, 0.483427, 0.158004, 0.301105, -0.009050, 0.888892, -2.406945, 0.077962, 1.058791, -1.308333, -25.432796, -0.327504, 0.426268, 0.000000, 0.000000], 
W => [ 0.510826, 0.510826, -0.587787, -0.587787, 0.105361, 1.203973, 0.510826, -0.587787, 0.510826, -23.613638, 1.609438, -23.613638, 0.105361, -0.587787, 0.000000, 0.000000], 
Y => [ -1.722767, -0.336472, -1.722767, -1.722767, -1.722767, -0.624154, -1.722767, -1.722767, -1.029619, -24.748618, -0.113329, -24.748618, -1.722767, 0.223144, 0.000000, 0.000000] 
);


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
 Throws    : Exception if given type is not protein.
 Comments  : 
See Also   : B<Bio::Root::Object::new()>, B<Bio::Seq::_initialize()>

=cut

#----------------
sub _initialize {
#----------------
    my($self, %param) = @_;
    
    my($threshold,$type) = $self->_rearrange([qw(THRESHOLD
                                                 TYPE)], %param);

    my $make = $self->SUPER::_initialize(%param);

    # complain if not protein
    if ($type =~ /nuc|[dr]na/i) {
	$self->throw("Sigcleave.pm only supports protein sequences.");
    } elsif ($type =~ /amino|pep|prot/i) {
	$self->{type} = 'amino';
    }

    # set threshold if supplied, otherwise default to 3.5
    if (defined $threshold) {
	$self->{threshold} = $threshold;
    } else {
      $self->{threshold} = 3.5;
    }

    $self->_Analyze;
    $make;
}


=head1 _Analyze

 Title     : _Analyze
 Usage     : N/A This is an internal method. Not meant to be called from outside
           : the package
           :
 Purpose   : calculates sigcleave score and amino acid position for the
           : given protein sequence. The score reporting threshold can
           : be adjusted by passing in the "threshold" parameter during
           : object construction. If no threshold is passed in, the code
           : defaults to reporting any scores equal to or above 3.5
           :
 Returns   : nothing. results are added to the object
 Argument  : none.
 Throws    : nothing.
 Comments  : nothing.
See Also   : n/a

=cut

#----------------
sub _Analyze {
#----------------
my($self) = @_;

## need to shut strict() up

my %signals;
my @hitWeight = ();
my @hitsort   = ();
my @hitpos    = ();
my $maxSite   = "";
my $seqPos    = "";
my $istart    = "";
my $iend      = "";
my $icol      = "";
my $i         = "";
my $weight    = "";
my $k         = 0;
my $c         = 0;
my $seqBegin  = 0; 
my $pVal      = -13; 
my $nVal      = 2;
my $nHits     = 0;
my $seqEnd    = $self->seq_len;
my $pep       = $self->str;
my $minWeight = $self->threshold;

 ## The weight table is keyed by UPPERCASE letters so we uppercase
 ## the pep string because we don't want to alter the actual object
 ## sequence.

 $pep =~ tr/a-z/A-Z/;

for($seqPos = $seqBegin; $seqPos < $seqEnd; $seqPos++)
  {
    $istart = (0 > $seqPos + $pVal)? 0 : $seqPos + $pVal;
    $iend = ($seqPos + $nVal - 1 < $seqEnd)? $seqPos + $nVal - 1 : $seqEnd;
    $icol= $iend - $istart + 1; 
    $weight = 0.00;
    for ($k=0; $k<$icol; $k++)
     {
	$c = substr($pep, $istart + $k, 1);

        ## CD: The if(defined) stuff was put in here because Sigcleave.pm
        ## CD: kept getting warnings about undefined vals during 'make test' ... 

	if(defined $WeightTable{$c}[$k]) { $weight += $WeightTable{$c}[$k]; }

     }

    if ($weight >= $minWeight)
        {  $signals{$seqPos+1} = sprintf ("%.1f", $weight);   }
  }

  $self->{"signal_scores"} = { %signals };
}



=head1 threshold

 Title     : threshold
 Usage     : $value = $self->threshold 
           :
 Purpose   : Accessor method sigcleave score reporting threshold.
           : 
 Returns   : float.
           : 
 Argument  : none. 
 Throws    : none.
 Comments  : none.
See Also   : n/a

=cut

#----------------
sub threshold {
#----------------
my $self = shift;
return $self->{threshold};
}



=head1 signals

 Title     : signals
 Usage     : %sigcleave_results = $sigcleave_object->signals;
           :
 Purpose   : Accessor method for sigcleave results
           : 
 Returns   : Associative array. The key value represents the amino acid position
           : and the value represents the score. Only scores that
           : are greater than or equal to the THRESHOLD value are reported.
           : 
 Argument  : none. 
 Throws    : none.
 Comments  : none.
See Also   : THRESHOLD

=cut

#----------------
sub signals {
#----------------
my $self = shift;
my %results;
my $position;

 foreach $position ( sort keys %{ $self->{signal_scores} } ) {
     $results{$position} = $self->{signal_scores}{$position};    
 }

return %results;
}



=head1 pretty_print

 Title     : pretty_print
 Usage     : $output = $sigcleave_object->pretty_print;
           : print $sigcleave_object->pretty_print;
           :
 Purpose   : Emulates the output of the EGCG Sigcleave
           : utility.  
           : 
 Returns   : A formatted string.
 Argument  : none.
 Throws    : none.
 Comments  : none.
See Also   : n/a

=cut

#----------------
sub pretty_print {
#----------------
my $self = shift;
my $pos;
my $output;
my $cnt = 1;
my %results  = $self->signals;
my @hits     = keys %results;
my $hitcount = $#hits; $hitcount++;
my $thresh   = $self->threshold;
my $seqlen   = $self->seq_len;
my $name     = $self->id;
my $pep      = $self->str;
   $pep      =~ tr/a-z/A-Z/;

$output = "SIGCLEAVE of $name from: 1 to $seqlen\n\n";

if($hitcount > 0) {
   $output .= "Report scores over $thresh\n";
	foreach $pos ((sort { $results{$b} cmp $results{$a} } keys %results)) {
 	  my $start = $pos - 13;
          $output .= sprintf ("Maximum score %1.1f at residue %3d\n",$results{$pos},$pos + 1);
    	  $output .= "\n";
    	  $output .= " Sequence:  ";
    	  $output .= substr($pep,$start,13);
    	  $output .= "-";
    	  $output .= substr($pep,$pos,50);
    	  $output .= "\n";
    	  $output .= "            | \(signal\)    | \(mature peptide\)\n";
    	  $output .= sprintf("          %3d            %3d\n\n",$start+1,$pos+1);

   		if(($hitcount > 1) && ($cnt == 1)) { 
    		$output .= " Other entries above $thresh\n\n";
   		} 
    	$cnt++;
 	}
}
$output;
}


1;
__END__


#########################################################################
#  End of class 
#########################################################################








