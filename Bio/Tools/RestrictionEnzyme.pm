#-----------------------------------------------------------------------------
# PACKAGE : Bio::Tools::RestrictionEnzyme.pm
# AUTHOR  : Steve A. Chervitz (sac@genome.stanford.edu)
# CREATED : 3 June 1997
# REVISION: $Id$
# STATUS  : Alpha
#            
# MODIFIED: 
#  sac --- Tue Nov 24 04:39:59 1998
#    * Removed -terse setting in constructor (deprecated).
#  sac --- Mon Jul 13 15:35:08 1998
#    * Fixed parameter lists (added hyphens to tags).
#  sac --- Mon Nov 17 13:57:24 1997
#    * Added a few new restriction enzymes to the %RE hash.
#    * Fixed bug in _make_custom() (now strips white space from RE site).
#
# Copyright (c) 1997-8 Steve A. Chervitz. All Rights Reserved.
#           This module is free software; you can redistribute it and/or 
#           modify it under the same terms as Perl itself.
#-----------------------------------------------------------------------------

package Bio::Tools::RestrictionEnzyme;
use strict;

use Bio::Root::RootI;
use Exporter;

use vars qw (@ISA @EXPORT_OK %EXPORT_TAGS $ID $version @RE_available $Revision);

@ISA         = qw(Bio::Root::RootI Exporter);
@EXPORT_OK   = qw(@RE_available);
%EXPORT_TAGS = ( std => [qw(@RE_available)] );

$ID = 'Bio::Tools::RestrictionEnzyme';
$version = 0.04;
$Revision = '$Id$';  #'

# Generated from REBASE version 802 (strider format), dated Jan 29 98
# by rebase2perl.pl (JA Feb 98). Merged with previous list by Ewan, Nov 1998
# Syntax: RE-name => 'SITE CUTS-AT' where SITE and CUTS-AT are separated by a space.

my %RE = (
          'AatII'   => 'GACGTC 5',
          'AccI'    => 'GTMKAC 2',
          'AclI'    => 'AACGTT 2',
          'AcyI'    => 'GRCGYC 2',
          'AflII'   => 'CTTAAG 1',
          'AflIII'  => 'ACRYGT 1',
          'AgeI'    => 'ACCGGT 1',
          'AhaIII'  => 'TTTAAA 3',
	  'AhdI'    => 'GACNNNNNGTC 6',
          'AluI'    => 'AGCT 2',
          'AlwNI'   => 'CAGNNNCTG 6',
          'ApaBI'   => 'GCANNNNNTGC 8',
          'ApaI'    => 'GGGCCC 5',
          'ApaLI'   => 'GTGCAC 1',
          'ApoI'    => 'RAATTY 1',
          'AscI'    => 'GGCGCGCC 2',
          'AsuI'    => 'GGNCC 1',
          'AsuII'   => 'TTCGAA 2',
          'AvaI'    => 'CYCGRG 1',
          'AvaII'   => 'GGWCC 1',
          'AvrII'   => 'CCTAGG 1',
          'BalI'    => 'TGGCCA 3',
          'BamHI'   => 'GGATCC 1',
          'BclI'    => 'TGATCA 1',
          'BetI'    => 'WCCGGW 1',
          'BglI'    => 'GCCNNNNNGGC 7',
          'BglII'   => 'AGATCT 1',
          'BsaAI'   => 'YACGTR 3',
          'BsaBI'   => 'GATNNNNATC 5',
          'BsePI'   => 'GCGCGC 1',
          'BsiYI'   => 'CCNNNNNNNGG 7',
          'Bsp1407I'=> 'TGTACA 1',
          'BspHI'   => 'TCATGA 1',
          'BspLU11I'=> 'ACATGT 1',
          'BspMII'  => 'TCCGGA 1',
          'BstEII'  => 'GGTNACC 1',
          'BstXI'   => 'CCANNNNNNTGG 8',
          'Cac8I'   => 'GCNNGC 3',
          'CauII'   => 'CCSGG 2',
          'Cfr10I'  => 'RCCGGY 1',
          'CfrI'    => 'YGGCCR 1',
          'ClaI'    => 'ATCGAT 2',
          'CviJI'   => 'RGCY 2',
          'CviRI'   => 'TGCA 2',
          'DdeI'    => 'CTNAG 1',
          'DpnI'    => 'GATC 2',
	  'DraI'    => 'TTTAAA 3',
          'DraII'   => 'RGGNCCY 2',
          'DraIII'  => 'CACNNNGTG 6',
          'DrdI'    => 'GACNNNNNNGTC 7',
          'DsaI'    => 'CCRYGG 1',
          'Eam1105I'=> 'GACNNNNNGTC 6',
          'Eco47III'=> 'AGCGCT 3',
          'EcoNI'   => 'CCTNNNNNAGG 5',
          'EcoRI'   => 'GAATTC 1',
          'EcoRII'  => 'CCWGG 0',
          'EcoRV'   => 'GATATC 3',
          'EspI'    => 'GCTNAGC 2',
          'Fnu4HI'  => 'GCNGC 2',
          'FnuDII'  => 'CGCG 2',
          'FseI'    => 'GGCCGGCC 6',
          'HaeI'    => 'WGGCCW 3',
          'HaeII'   => 'RGCGCY 5',
          'HaeIII'  => 'GGCC 2',
          'HgiAI'   => 'GWGCWC 5',
          'HgiCI'   => 'GGYRCC 1',
          'HgiJII'  => 'GRGCYC 5',
          'HhaI'    => 'GCGC 3',
	  'HincII'  => 'GTYRAC 3',
          'HindII'  => 'GTYRAC 3',
          'HindIII' => 'AAGCTT 1',
          'HinfI'   => 'GANTC 1',
          'HpaI'    => 'GTTAAC 3',
          'HpaII'   => 'CCGG 1',
          'KpnI'    => 'GGTACC 5',
          'MaeI'    => 'CTAG 1',
          'MaeII'   => 'ACGT 1',
          'MaeIII'  => 'GTNAC 0',
          'MboI'    => 'GATC 0',
          'McrI'    => 'CGRYCG 4',
          'MfeI'    => 'CAATTG 1',
          'MluI'    => 'ACGCGT 1',
          'MseI'    => 'TTAA 1',
          'MslI'    => 'CAYNNNNRTG 5',
          'MstI'    => 'TGCGCA 3',
          'MwoI'    => 'GCNNNNNNNGC 7',
          'NaeI'    => 'GCCGGC 3',
          'NarI'    => 'GGCGCC 2',
          'NcoI'    => 'CCATGG 1',
          'NdeI'    => 'CATATG 2',
          'NheI'    => 'GCTAGC 1',
          'NlaIII'  => 'CATG 4',
          'NlaIV'   => 'GGNNCC 3',
          'NotI'    => 'GCGGCCGC 2',
          'NruI'    => 'TCGCGA 3',
          'NspBII'  => 'CMGCKG 3',
          'NspI'    => 'RCATGY 5',
          'PacI'    => 'TTAATTAA 5',
          'PflMI'   => 'CCANNNNNTGG 7',
          'PmaCI'   => 'CACGTG 3',
          'PmeI'    => 'GTTTAAAC 4',
          'PpuMI'   => 'RGGWCCY 2',
          'PshAI'   => 'GACNNNNGTC 5',
          'PstI'    => 'CTGCAG 5',
          'PvuI'    => 'CGATCG 4',
          'PvuII'   => 'CAGCTG 3',
          'RsaI'    => 'GTAC 2',
          'RsrII'   => 'CGGWCCG 2',
          'SacI'    => 'GAGCTC 5',
          'SacII'   => 'CCGCGG 4',
	  'Sau96I'  => 'GGNCC 1',
          'SalI'    => 'GTCGAC 1',
          'SanDI'   => 'GGGWCCC 2',
          'SauI'    => 'CCTNAGG 2',
	  'SbfI'    => 'CCTGCAGG 6',
          'ScaI'    => 'AGTACT 3',
          'ScrFI'   => 'CCNGG 2',
          'SduI'    => 'GDGCHC 5',
          'SecI'    => 'CCNNGG 1',
          'SexAI'   => 'ACCWGGT 1',
          'SfeI'    => 'CTRYAG 1',
          'SfiI'    => 'GGCCNNNNNGGCC 8',
          'SgfI'    => 'GCGATCGC 5',
          'SgrAI'   => 'CRCCGGYG 2',
          'SmaI'    => 'CCCGGG 3',
          'SmlI'    => 'CTYRAG 1',
          'SnaBI'   => 'TACGTA 3',
          'SpeI'    => 'ACTAGT 1',
          'SphI'    => 'GCATGC 5',
          'SplI'    => 'CGTACG 1',
          'SrfI'    => 'GCCCGGGC 4',
          'Sse8387I'=> 'CCTGCAGG 6',
          'Sse8647I'=> 'AGGWCCT 2',
          'SspI'    => 'AATATT 3',
          'StuI'    => 'AGGCCT 3',
          'StyI'    => 'CCWWGG 1',
          'SwaI'    => 'ATTTAAAT 4',
          'TaqI'    => 'TCGA 1',
          'TatI'    => 'WGTACW 1',
          'TfiI'    => 'GAWTC 1',
          'TseI'    => 'GCWGC 1',
          'Tsp45I'  => 'GTSAC 0',
          'Tsp4CI'  => 'ACNGT 3',
          'TspEI'   => 'AATT 0',
          'TspRI'   => 'CASTGNN 7',
          'Tth111I' => 'GACNNNGTC 4',
          'VspI'    => 'ATTAAT 2',
          'XbaI'    => 'TCTAGA 1',
          'XcmI'    => 'CCANNNNNNNNNTGG 8',
          'XhoI'    => 'CTCGAG 1',
          'XhoII'   => 'RGATCY 1',
          'XmaIII'  => 'CGGCCG 1',
          'XmnI'    => 'GAANNNNTTC 5',
         );

@RE_available = sort keys %RE;

## POD Documentation:

=head1 NAME

Bio::Tools::RestrictionEnzyme.pm - Bioperl object for a restriction endonuclease object.

=head1 SYNOPSIS

=head2 Object Creation

    require Bio::Tools::RestrictionEnzyme;

    ## Create a new object by name.

    $re1 = new Bio::Tools::RestrictionEnzyme(-NAME =>'EcoRI');

    ## Create a new object using special syntax
    ## which specifies the enzyme name, recognition site, and cut position.
    ## Used for enzymes not known to this module.

    $re2 = new Bio::Tools::RestrictionEnzyme(-NAME =>'EcoRV--GAT^ATC', 
				  	     -MAKE =>'custom');


=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

=head1 DESCRIPTION

The Bio::Tools::RestrictionEnzyme.pm module encapsulates generic data and 
methods for using restriction endonucleases for in silico restriction
analysis of DNA sequences.

=head2 Considerations

This module is a precursor for a more full featured version that may do such things as
download data from online databases such as REBase http://www.neb.com/rebase/.
Thus, there is currently no functionality for obtaining data about commercial
availability for a restriction enzyme.

At some point in the future, it may be best to derive RestrictionEnzymes from
a class such as Bio::Enzyme.pm or Bio::Prot::Protein.pm so that more data about 
the enzyme and related information can be easily obtained.

This module is currently in use at 

 http://genome-www.stanford.edu/Sacch3D/analysis/

B<This module is at an early stage of development and is not yet ready for general use. API documentation is presently incomplete.>


=head1 DEPENDENCIES 

Bio::Tools::RestrictionEnzyme.pm is a concrete class that inherits from B<Bio::Root::Object.pm>
and uses by delegation B<Bio::Seq.pm>.

=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other Bioperl modules.
Send your comments and suggestions preferably to one of the Bioperl mailing lists.
Your participation is much appreciated.

   bioperl-l@bioperl.org             - General discussion
   http://bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the bugs and 
their resolution. Bug reports can be submitted via email or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR

Steve A. Chervitz, sac@genome.stanford.edu

=head1 VERSION

Bio::Tools::RestrictionEnzyme.pm, 0.04

=head1 COPYRIGHT

Copyright (c) 1997-8 Steve A. Chervitz. All Rights Reserved.
This module is free software; you can redistribute it and/or 
modify it under the same terms as Perl itself.

=head1 SEE ALSO

  Bio::Root::Object.pm    - Base class.
  Bio::Seq.pm             - Lightweight sequence object.

  http://bio.perl.org/Projects/modules.html  - Online module documentation
  http://bio.perl.org/Projects/Blast/        - Bioperl Blast Project     
  http://bio.perl.org/                       - Bioperl Project Homepage

=cut

#
##
###
#### END of main POD documentation.
###
##
#'


=head1 APPENDIX

Methods beginning with a leading underscore are considered private
and are intended for internal use by this module. They are
B<not> considered part of the public interface and are described here
for documentation purposes only.

=cut


#######################################################
#               CONSTRUCTOR/DESTRUCTOR                #
#######################################################


=head1 new

 Title     : new
 Purpose   : Initializes the RestrictionEnzyme object and calls
           : superclass constructor last (Bio:Seq.pm).
 Returns   : n/a
 Argument  : Parameters passed to new()
 Comments  : The process of creating a new SeqPattern.pm object
           : ensures that the pattern string is untained.

See Also   : L<_make_custom>(), L<_make_standard>(), B<Bio::Seq.pm::_initialize()>

=cut

#---------------
sub new {
#---------------
    my($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    my ($name,$make) = $self->_rearrange([qw(NAME MAKE)],@args);
    $name && $self->name($name);
    my %data;
    if(defined $make && $make eq 'custom') {
	%data = $self->_make_custom($name); 
    } else {
	%data = $self->_make_standard($name);
    }
    $self->{'_seq'} = new Bio::Seq(%data, 
				   -VERBOSE =>$self->verbose,
				   -moltype => 'dna',
				   );
    return $self;
}


=head1 _make_standard

 Title     : _make_standard
 Usage     : n/a; automatically called by _initialize()
 Purpose   : Permits custom RE object construction from name.
           : 'EcoRI'.
 Returns   : Hash containing named parameters for Bio::Seq.pm constructor.
 Argument  : String containing string with special syntax.
 Throws    : Exception if the requested enzyme name is unavailable.
           : NOTE: Case sensitive.

See Also   : L<_initialize>(), L<_make_custom>()

=cut

#------------------
sub _make_standard {
#------------------
    my($self, $name) = @_;

    $name =~ s/^\s+|\s+$//g;
 
    $self->is_available($name) || 
	$self->throw("Unavailable or undefined enzyme: $name (Note: CASE SENSITIVE)",
		     "Currently available enzymes: \n@RE_available\n");

    my @data = split( ' ', $RE{$name});
    my (%dat);
    $dat{-SEQ} = $data[0];
    $dat{-NAME} = $dat{-ID}= $name;
    $self->{'_cuts_after'} = $data[1];

    return %dat;
}


=head1 _make_custom

 Title     : _make_custom
 Usage     : n/a; automatically called by _initialize()
 Purpose   : Permits custom RE object construction from strings 
           : such as 'EcoRI--G^AATTC' as the name of the enzyme.
 Returns   : Hash containing named parameters for Bio::Seq.pm constructor.
 Argument  : String containing string with special syntax.
 Throws    : Exception if the string has bad syntax.
           : Warning if the string did not specify cut position.
           :         Places cut site after 5'-most position.

See Also   : L<_initialize>()

=cut

#'
#-----------------
sub _make_custom {
#-----------------
    my($self, $name) = @_;

    $name =~ s/\s+//g;
    my @parts  = split '--', $name;
    my (%dat);
    $dat{-NAME} = $dat{-ID} = $parts[0];
    $self->name($parts[0]);  ## Reset name

    $parts[1] || return $self->throw("Undefined recognition site for $parts[0].",
				      "Use this syntax: EcoRV--GAT^ATC");
    ## Determine the cuts_after point.
    my $cut_index = index $parts[1], '^';
    if( $cut_index <0) { $cut_index = 0;
			 $self->warn("Unknown cut position for $parts[0]. Assuming position 0",
				     "Use carat to specify cut position (e.g., G^AATTC)"); }
    $self->{'_cuts_after'} =  $cut_index;

    ## Save the recognition sequence after removing the '^'
    $parts[1] =~ s/\^//g;
    $dat{-SEQ} = $parts[1];
    return %dat;
}
    

=head1 cuts_after

 Title     : cuts_after
 Usage     : $re->cuts_after();
 Purpose   : Sets/Gets the position of cleavage relative to the 5' end.
 Example   : $num = $re->cuts_after() 
 Returns   : Integer
 Argument  : Integer (optional)
 Throws    : Exception if argument is non-numeric.
 Access    : Public
 Comments  : This method is only needed to change the cuts at
           : position. This data is automatically set during
           : construction.

See Also   : L<_make_standard>(), L<_make_custom>()

=cut

#'
#---------------
sub cuts_after { 
#---------------
    my $self = shift; 
    if(@_) { my $num = shift;
	     if($num == 0 and $num ne '0') {
		 $self->throw("Bad number: $num", 
			      "The cuts_after position be an integer.");
	     }
	     $self->{'_cuts_after'} = $num;
	 }
    $self->{'_cuts_after'}; 
}



=head1 site

 Title     : site
 Usage     : $re->site();
 Purpose   : Gets the recognition sequence for the enzyme. 
 Example   : $seq_string = $re->site();
 Returns   : String containing recognition sequence indicating 
           : cleavage site as in  'G^AATTC'.
 Argument  : n/a
 Throws    : n/a

=cut

#---------
sub site {
#---------
    my $self = shift;
    my $seq = $self->seq;
    my $cuts_after = $self->cuts_after;
    if    (!$cuts_after)                { return '^' . $seq->subseq(1, $seq->length); }
    elsif ($cuts_after == $seq->length) { return $seq->subseq(1, $seq->length) . '^'; }
    else                                { return $seq->subseq(1, $cuts_after) . '^' . $seq->subseq($cuts_after + 1, $seq->length); }
}
    



=head1 seq

 Title     : seq
 Usage     : $re->seq();
 Purpose   : Get the Bio::Seq.pm-derived object representing 
           : the recognition sequence
 Returns   : String
 Argument  : n/a
 Throws    : n/a

See Also   : L<string>(), L<revcom>()

=cut

#---------
sub seq    {  my $self = shift; $self->{'_seq'}; }
#---------



=head1 string

 Title     : string
 Usage     : $re->string();
 Purpose   : Get a string representing the recognition sequence.
 Returns   : String
 Argument  : n/a
 Throws    : n/a
 Comments  : Delegates to the Bio::Seq.pm-derived object.

See Also   : L<seq>(), L<revcom>()

=cut

#-----------
sub string {  my $self = shift; $self->{'_seq'}->seq; }
#-----------



=head1 revcom

 Title     : revcom
 Usage     : $re->revcom();
 Purpose   : Get a string representing the reverse complement of
           : the recognition sequence.
 Returns   : String
 Argument  : n/a
 Throws    : n/a
 Comments  : Delegates to the Bio::Seq.pm-derived object, but needs to get
             out the string from it, as now Bio::Seq->revcom makes a Bio::Seq
             object

See Also   : L<seq>(), L<string>()

=cut

#-----------
sub revcom {  my $self = shift; $self->{'_seq'}->revcom->seq(); }
#-----------



=head1 cut_seq

 Title     : cut_seq
 Usage     : $re->cut_seq(<sequence object>);
 Purpose   : Conceptually cut or "digest" a DNA sequence with the given enzyme.
 Example   : $string = $re->cut_seq(<sequence object>); 
 Returns   : List of strings containing the resulting fragments.
 Argument  : Reference to a Bio::Seq.pm-derived object.
 Throws    : Exception if argument is not an object.
           : (Does not yet verify that it is derived from Bio::Seq.pm.)
 Comments  : Strategy relies on Perl's built-in split() function.
           : Since split removes the recognition pattern, the resulting
           : fragments must be repaired after split()-ing.
           : There is currently no support for partial digestions.

=cut

#'
#-------------
sub cut_seq {
#-------------
    my( $self, $seqObj) = @_;

    # Could check that $seqObj is derived from Seq (Perl 5.004).
    ref $seqObj || $self->throw( "Can't cut sequence. Missing or invalid object",
				 "seqObj: $seqObj");
    
#    print "$ID: cutting sequence.\n";

    my $cuts_after = $self->{'_cuts_after'};
    my ($site_3prime_seq, $site_5prime_seq);
    my $reSeq = $self->seq;
    if($cuts_after == 0) {
	$site_3prime_seq = '';
	$site_5prime_seq = $reSeq->seq();
    } elsif($cuts_after == $reSeq->length) {
	$site_3prime_seq = $reSeq->seq();
	$site_5prime_seq = '';
    } else {
	$site_3prime_seq = $reSeq->subseq(1, $self->{'_cuts_after'});
	$site_5prime_seq = $reSeq->subseq($self->{'_cuts_after'}+1, $reSeq->length);
    }

#    print "3' site: $site_3prime_seq\n5' site: $site_5prime_seq";<STDIN>;

    my(@re_frags);
    my $seq = $reSeq->seq;
    $seq =~ s/N/\./g;
    $seq =~ s/R/\[AG\]/g;
    $seq =~ s/Y/\[CT\]/g;
    $seq =~ s/S/\[GC\]/g;
    $seq =~ s/W/\[AT\]/g;
    if(!$self->palindromic) {
	my $revseq = $reSeq->revcom;
	$revseq =~ s/N/\./g;
	$revseq =~ s/R/\[AG\]/g;
	$revseq =~ s/Y/\[CT\]/g;
	$revseq =~ s/S/\[GC\]/g;
	$revseq =~ s/W/\[AT\]/g;
	$seq .= '|'.$revseq;
    }
#    printf "$ID: site seq: %s\n\n", $seq;
#    printf "$ID: splitting %s\n\n",$reSeq->str;
    @re_frags = split(/$seq/, $seqObj->seq);

#    print "$ID: cut_seq, ",scalar @re_frags, " fragments.\n";

    ## Re-attach the split recognition site back to the frags
    ## since perl zapped them in the split() call.
    my($i);
    my $numFrags = scalar @re_frags;
    for($i=0; $i<$numFrags; $i++) {
	$i < $#re_frags  and $re_frags[$i] = $re_frags[$i].$site_3prime_seq;
	$i > 0           and $re_frags[$i] = $site_5prime_seq.$re_frags[$i];
    }

    @re_frags;
}

=head1 cut_locations

 Title     : cut_locations
 Usage     : my $locations = $re->cut_locations(<sequence_object>);
 Purpose   : Report the location of the recognition site(s) within
           : an input sequence. 
 Example   : my $locations = $re->annotate_seq($seqObj);
 Returns   : Arrayref of starting locations where enzyme would cut 
 Argument  : Reference to a Bio::SeqI-derived sequence object.
 Throws    : n/a
 Comments  : 

=cut

#-----------------
sub cut_locations {
#-----------------
    my($self, $seqobj) = @_;

    my $site = $self->string;
    my $seq = $seqobj->seq;
    study($seq);
    $site =~ s/N|X/\./g;
    $site =~ s/R/\[AG\]/g;
    $site =~ s/Y/\[CT\]/g;
    $site =~ s/S/\[GC\]/g;
    $site =~ s/W/\[AT\]/g;
    my @locations;
    while( $seq =~ /($site)/g ) {
        # $` is preceding string before pattern so length returns position
	push @locations, length($`); 	
    }
    return \@locations;
}    


=head1 annotate_seq

 Title     : annotate_seq
 Usage     : $re->annotate_seq(<sequence_object>);
 Purpose   : Identify the location of the recognition site(s) within
           : an input sequence. Uses HTML.
 Example   : $annot_seq = $re->annotate_seq($seqObj);
 Returns   : String containing the annotated sequence.
 Argument  : Reference to a Bio::Seq.pm-derived sequence object.
 Throws    : n/a
 Comments  : The annotated sequence must be viewed with a web
           : browser to see the location(s) of the recognition site(s).

=cut

#-----------------
sub annotate_seq {
#-----------------
    my($self, $seqObj) = @_;

    my $site = $self->string;
    my $seq = $seqObj->seq;

    $site =~ s/N|X/\./g;
    $site =~ s/R/\[AG\]/g;
    $site =~ s/Y/\[CT\]/g;
    $site =~ s/S/\[GC\]/g;
    $site =~ s/W/\[AT\]/g;

    $seq =~ s|$site|<b>$site</b>|g;
    return $seq;
}    


=head1 palindromic

 Title     : palindromic
 Usage     : $re->palindromic();
 Purpose   : Determines if the recognition sequence is palindromic
           : for the current restriction enzyme.
 Returns   : Boolean
 Argument  : n/a
 Throws    : n/a
 Access    : Public 
 Comments  : A palindromic site (EcoRI): 5-GAATTC-3
           :                             3-CTTAAG-5

=cut

#----------------
sub palindromic {
#----------------
    my $self = shift;
    $self->string eq $self->revcom;
}



=head1 is_available

 Title     : is_available
 Usage     : $re->is_available(<string containing name of enzyme>);
 Purpose   : Determine if an enzyme is available (to this module).
           : (see the package lexical %RE).
 Example   : $re->is_available('EcoRI');
           : &Bio::Tools::RestrictionEnzyme::is_available($object,'EcoRI');
 Returns   : Boolean
 Argument  : String
 Throws    : n/a
 Comments  : This method does NOT give information about
           : commercial availability (yet). 
           : Enzyme names are CASE SENSITIVE.

See Also   : L<available_list>()

=cut

#----------------
sub is_available {
#----------------
    my($self,$name) = @_;
    exists $RE{$name};
}

#--------------
sub available {
#--------------
    my($self,$name) = @_;
    print STDERR "\nDeprecated method: $ID:: available(); ".
	"use is_available() instead.\n";
    $self->is_available($name);
}


=head2 name

 Title   : name
 Usage   : $obj->name($newval)
 Function: 
 Example : 
 Returns : value of name
 Args    : newvalue (optional)


=cut

sub name{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'name'} = $value;
    }
    return $obj->{'name'};

}

=head1 available_list

 Title     : available_list
 Usage     : $re->available_list([<integer>]);
 Purpose   : Retrieve a list of currently available enzymes.
 Example   : @all = $re->available_list();  ## All enzymes
           : @six_cutters = $re->available_list(6);  ## All 6-cutters
 Returns   : List of strings
 Argument  : Integer (optional)
 Throws    : n/a
 Comments  : This method may be more appropriate for a REData.pm class.

See Also   : L<is_available>()

=cut

#-------------------
sub available_list {
#-------------------
    my($self,$size) = @_;
    $size ||= 'all';

    $size eq 'all' and return @RE_available;

    my(@data, @names);
    foreach (@RE_available) {
	@data = split /\s/, $RE{$_};
	if(length $data[0] == $size) {
	    push @names, $_;
	}
    }
    @names;
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
 _seq         : A Bio::Seq.pm-derived object.
              :
 _site        : String containing the recognition sequence.
              :
 _cuts_after  : Integer indicating the cleavage position relative to the 
              : 5' end of the recognition sequence.

 INHERITED DATA MEMBERS:

 _name      : (From Bio::Bio::Root::Object.pm) String containing name of the enzyme.

=cut
