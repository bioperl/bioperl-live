#------------------------------------------------------------------
# $Id$
#
# BioPerl module Bio::Tools::RestrictionEnzyme
#
# Cared for by Steve Chervitz <sac@bioperl.org>
#
# You may distribute this module under the same terms as perl itself
#------------------------------------------------------------------

## POD Documentation:

=head1 NAME

Bio::Tools::RestrictionEnzyme - Bioperl object for a restriction endonuclease
(cuts DNA at specific locations)

DEPRECATED. Please use the Bio::Restriction modules instead.

=head1 SYNOPSIS

  use Bio::Tools::RestrictionEnzyme;

  ## Create a new object by name.

  $re1 = Bio::Tools::RestrictionEnzyme->new(-NAME =>'EcoRI');

  ## Create a new object using special syntax
  ## which specifies the enzyme name, recognition site, and cut position.
  ## Used for enzymes not known to this module.

  $re2 = Bio::Tools::RestrictionEnzyme->new(-NAME =>'EcoRV--GAT^ATC', 
                                            -MAKE =>'custom');

  ## Get a list of the resulting fragments when a sequence is cut with
  ## the given enzyme. The method expects a Bio::Seq object.

  @fragments = $re2->cut_seq($seqobj);

  ## Get a list of names of all available restriction enzymes 
  ## known to this module.

  @all = $re->available_list();

  ## Get the names of restriction enzymes that have 6 bp 
  ## recognition sequences.

  @sixcutters = $re->available_list(6);

=head1 DESCRIPTION

The Bio::Tools::RestrictionEnzyme.pm module encapsulates generic data and 
methods for using restriction endonucleases for in silico restriction
analysis of DNA sequences.

=head2 Considerations

This module is a precursor for a more full featured version that may do such
things as download data from online databases such as REBase L<http://rebase.neb.com/>.
Thus, there is currently no functionality for obtaining data regarding commercial
availability of a restriction enzyme.

At some point in the future, it may make sense to derive RestrictionEnzymes
from a class such as Bio::Enzyme or Bio::Prot::Protein (neither of which now
exist) so that more data about the enzyme and related information can be 
easily obtained.

This module is currently in use at L<http://genome-www.stanford.edu/Sacch3D/analysis/>.

=head2 Digesting on Runs of N

To digest a sequence on runs of N's in the sequence. Here's what you can do:

    $re_n  = Bio::Tools::RestrictionEnzyme->new(-name=>'N--NNNNN',
                                                -make=>'custom');

Specify the number of N's you want to match in the -name parameter. 
So the above example will recognize and cut at runs of 5 Ns.
If you wanted to cut at runs of 10 N's, you would use 

     -name => 'N--NNNNNNNNNN'

Note that you must use a specific number of N's, you cannot use a regexp to
digest at N+ for example, because the actual number of N's at each site are
not recorded when the sequence is analyzed. So cut_locations( ) wouldn't be 
correct. 

=head1 EXAMPLES

See the script C<examples/restriction.pl> in the Bioperl distribution.

=head1 DEPENDENCIES 

Bio::Tools::RestrictionEnzyme is a concrete class that inherits from 
L<Bio::Root::Root> and uses by delegation L<Bio::PrimarySeq>.

=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other Bioperl
modules. Send your comments and suggestions preferably to one of the Bioperl
mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR

Steve Chervitz, E<lt>sac-at-bioperl.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 1997-2002 Steve A. Chervitz. All Rights Reserved.
This module is free software; you can redistribute it and/or 
modify it under the same terms as Perl itself.

=head1 SEE ALSO

L<Bio::Root::Root>, 
L<Bio::PrimarySeq>.

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


package Bio::Tools::RestrictionEnzyme;
use strict;

use vars qw (@EXPORT_OK %EXPORT_TAGS @RE_available);

use base qw(Bio::Root::Root Exporter);
@EXPORT_OK   = qw(@RE_available);
%EXPORT_TAGS = ( std => [qw(@RE_available)] );

# Generated from REBASE version 208 (strider format), dated Aug 1 2002
# using scripts/contributed/rebase2list.pl
# Syntax: RE-name => 'SITE CUTS-AT' where SITE and CUTS-AT are separated 
# by a space.

my %RE = (
 'AasI'	=> 'GACNNNNNNGTC 7',
 'AatI'	=> 'AGGCCT 3',
 'AatII'	=> 'GACGTC 5',
 'AauI'	=> 'TGTACA 1',
 'AccI'	=> 'GTMKAC 2',
 'AccII'	=> 'CGCG 2',
 'AccIII'	=> 'TCCGGA 1',
 'Acc16I'	=> 'TGCGCA 3',
 'Acc65I'	=> 'GGTACC 1',
 'Acc113I'	=> 'AGTACT 3',
 'AccB1I'	=> 'GGYRCC 1',
 'AccB7I'	=> 'CCANNNNNTGG 7',
 'AclI'	=> 'AACGTT 2',
 'AcsI'	=> 'RAATTY 1',
 'AcvI'	=> 'CACGTG 3',
 'AcyI'	=> 'GRCGYC 2',
 'AdeI'	=> 'CACNNNGTG 6',
 'AfaI'	=> 'GTAC 2',
 'AfeI'	=> 'AGCGCT 3',
 'AflI'	=> 'GGWCC 1',
 'AflII'	=> 'CTTAAG 1',
 'AflIII'	=> 'ACRYGT 1',
 'AgeI'	=> 'ACCGGT 1',
 'AhaIII'	=> 'TTTAAA 3',
 'AhdI'	=> 'GACNNNNNGTC 6',
 'AhlI'	=> 'ACTAGT 1',
 'AleI'	=> 'CACNNNNGTG 5',
 'AluI'	=> 'AGCT 2',
 'Alw21I'	=> 'GWGCWC 5',
 'Alw44I'	=> 'GTGCAC 1',
 'AlwNI'	=> 'CAGNNNCTG 6',
 'Ama87I'	=> 'CYCGRG 1',
 'AocI'	=> 'CCTNAGG 2',
 'Aor51HI'	=> 'AGCGCT 3',
 'ApaI'	=> 'GGGCCC 5',
 'ApaBI'	=> 'GCANNNNNTGC 8',
 'ApaLI'	=> 'GTGCAC 1',
 'ApoI'	=> 'RAATTY 1',
 'AscI'	=> 'GGCGCGCC 2',
 'AseI'	=> 'ATTAAT 2',
 'AsiAI'	=> 'ACCGGT 1',
 'AsiSI'	=> 'GCGATCGC 5',
 'AsnI'	=> 'ATTAAT 2',
 'AspI'	=> 'GACNNNGTC 4',
 'Asp700I'	=> 'GAANNNNTTC 5',
 'Asp718I'	=> 'GGTACC 1',
 'AspEI'	=> 'GACNNNNNGTC 6',
 'AspHI'	=> 'GWGCWC 5',
 'AspLEI'	=> 'GCGC 3',
 'AspS9I'	=> 'GGNCC 1',
 'AsuI'	=> 'GGNCC 1',
 'AsuII'	=> 'TTCGAA 2',
 'AsuC2I'	=> 'CCSGG 2',
 'AsuNHI'	=> 'GCTAGC 1',
 'AvaI'	=> 'CYCGRG 1',
 'AvaII'	=> 'GGWCC 1',
 'AviII'	=> 'TGCGCA 3',
 'AvrII'	=> 'CCTAGG 1',
 'AxyI'	=> 'CCTNAGG 2',
 'BalI'	=> 'TGGCCA 3',
 'BamHI'	=> 'GGATCC 1',
 'BanI'	=> 'GGYRCC 1',
 'BanII'	=> 'GRGCYC 5',
 'BanIII'	=> 'ATCGAT 2',
 'BbeI'	=> 'GGCGCC 5',
 'BbrPI'	=> 'CACGTG 3',
 'BbuI'	=> 'GCATGC 5',
 'Bbv12I'	=> 'GWGCWC 5',
 'BclI'	=> 'TGATCA 1',
 'BcnI'	=> 'CCSGG 2',
 'BcoI'	=> 'CYCGRG 1',
 'BcuI'	=> 'ACTAGT 1',
 'BetI'	=> 'WCCGGW 1',
 'BfaI'	=> 'CTAG 1',
 'BfmI'	=> 'CTRYAG 1',
 'BfrI'	=> 'CTTAAG 1',
 'BfrBI'	=> 'ATGCAT 3',
 'BfuCI'	=> 'GATC 0',
 'BglI'	=> 'GCCNNNNNGGC 7',
 'BglII'	=> 'AGATCT 1',
 'BlnI'	=> 'CCTAGG 1',
 'BloHII'	=> 'CTGCAG 5',
 'BlpI'	=> 'GCTNAGC 2',
 'Bme18I'	=> 'GGWCC 1',
 'Bme1390I'	=> 'CCNGG 2',
 'Bme1580I'	=> 'GKGCMC 5',
 'BmtI'	=> 'GCTAGC 5',
 'BmyI'	=> 'GDGCHC 5',
 'BoxI'	=> 'GACNNNNGTC 5',
 'Bpu14I'	=> 'TTCGAA 2',
 'Bpu1102I'	=> 'GCTNAGC 2',
 'Bsa29I'	=> 'ATCGAT 2',
 'BsaAI'	=> 'YACGTR 3',
 'BsaBI'	=> 'GATNNNNATC 5',
 'BsaHI'	=> 'GRCGYC 2',
 'BsaJI'	=> 'CCNNGG 1',
 'BsaOI'	=> 'CGRYCG 4',
 'BsaWI'	=> 'WCCGGW 1',
 'BscI'	=> 'ATCGAT 2',
 'Bsc4I'	=> 'CCNNNNNNNGG 7',
 'BscBI'	=> 'GGNNCC 3',
 'BscFI'	=> 'GATC 0',
 'Bse8I'	=> 'GATNNNNATC 5',
 'Bse21I'	=> 'CCTNAGG 2',
 'Bse118I'	=> 'RCCGGY 1',
 'BseAI'	=> 'TCCGGA 1',
 'BseBI'	=> 'CCWGG 2',
 'BseCI'	=> 'ATCGAT 2',
 'BseDI'	=> 'CCNNGG 1',
 'BseJI'	=> 'GATNNNNATC 5',
 'BseLI'	=> 'CCNNNNNNNGG 7',
 'BsePI'	=> 'GCGCGC 1',
 'BseSI'	=> 'GKGCMC 5',
 'BseX3I'	=> 'CGGCCG 1',
 'BshI'	=> 'GGCC 2',
 'Bsh1236I'	=> 'CGCG 2',
 'Bsh1285I'	=> 'CGRYCG 4',
 'BshFI'	=> 'GGCC 2',
 'BshNI'	=> 'GGYRCC 1',
 'BshTI'	=> 'ACCGGT 1',
 'BsiBI'	=> 'GATNNNNATC 5',
 'BsiCI'	=> 'TTCGAA 2',
 'BsiEI'	=> 'CGRYCG 4',
 'BsiHKAI'	=> 'GWGCWC 5',
 'BsiHKCI'	=> 'CYCGRG 1',
 'BsiLI'	=> 'CCWGG 2',
 'BsiMI'	=> 'TCCGGA 1',
 'BsiQI'	=> 'TGATCA 1',
 'BsiSI'	=> 'CCGG 1',
 'BsiWI'	=> 'CGTACG 1',
 'BsiXI'	=> 'ATCGAT 2',
 'BsiYI'	=> 'CCNNNNNNNGG 7',
 'BsiZI'	=> 'GGNCC 1',
 'BslI'	=> 'CCNNNNNNNGG 7',
 'BsoBI'	=> 'CYCGRG 1',
 'Bsp13I'	=> 'TCCGGA 1',
 'Bsp19I'	=> 'CCATGG 1',
 'Bsp68I'	=> 'TCGCGA 3',
 'Bsp106I'	=> 'ATCGAT 2',
 'Bsp119I'	=> 'TTCGAA 2',
 'Bsp120I'	=> 'GGGCCC 1',
 'Bsp143I'	=> 'GATC 0',
 'Bsp143II'	=> 'RGCGCY 5',
 'Bsp1286I'	=> 'GDGCHC 5',
 'Bsp1407I'	=> 'TGTACA 1',
 'Bsp1720I'	=> 'GCTNAGC 2',
 'BspA2I'	=> 'CCTAGG 1',
 'BspCI'	=> 'CGATCG 4',
 'BspDI'	=> 'ATCGAT 2',
 'BspEI'	=> 'TCCGGA 1',
 'BspHI'	=> 'TCATGA 1',
 'BspLI'	=> 'GGNNCC 3',
 'BspLU11I'	=> 'ACATGT 1',
 'BspMII'	=> 'TCCGGA 1',
 'BspTI'	=> 'CTTAAG 1',
 'BspT104I'	=> 'TTCGAA 2',
 'BspT107I'	=> 'GGYRCC 1',
 'BspXI'	=> 'ATCGAT 2',
 'BsrBRI'	=> 'GATNNNNATC 5',
 'BsrFI'	=> 'RCCGGY 1',
 'BsrGI'	=> 'TGTACA 1',
 'BssAI'	=> 'RCCGGY 1',
 'BssECI'	=> 'CCNNGG 1',
 'BssHI'	=> 'CTCGAG 1',
 'BssHII'	=> 'GCGCGC 1',
 'BssKI'	=> 'CCNGG 0',
 'BssNAI'	=> 'GTATAC 3',
 'BssT1I'	=> 'CCWWGG 1',
 'Bst98I'	=> 'CTTAAG 1',
 'Bst1107I'	=> 'GTATAC 3',
 'BstACI'	=> 'GRCGYC 2',
 'BstAPI'	=> 'GCANNNNNTGC 7',
 'BstBI'	=> 'TTCGAA 2',
 'BstBAI'	=> 'YACGTR 3',
 'Bst4CI'	=> 'ACNGT 3',
 'BstC8I'	=> 'GCNNGC 3',
 'BstDEI'	=> 'CTNAG 1',
 'BstDSI'	=> 'CCRYGG 1',
 'BstEII'	=> 'GGTNACC 1',
 'BstENI'	=> 'CCTNNNNNAGG 5',
 'BstENII'	=> 'GATC 0',
 'BstFNI'	=> 'CGCG 2',
 'BstH2I'	=> 'RGCGCY 5',
 'BstHHI'	=> 'GCGC 3',
 'BstHPI'	=> 'GTTAAC 3',
 'BstKTI'	=> 'GATC 3',
 'BstMAI'	=> 'CTGCAG 5',
 'BstMCI'	=> 'CGRYCG 4',
 'BstMWI'	=> 'GCNNNNNNNGC 7',
 'BstNI'	=> 'CCWGG 2',
 'BstNSI'	=> 'RCATGY 5',
 'BstOI'	=> 'CCWGG 2',
 'BstPI'	=> 'GGTNACC 1',
 'BstPAI'	=> 'GACNNNNGTC 5',
 'BstSCI'	=> 'CCNGG 0',
 'BstSFI'	=> 'CTRYAG 1',
 'BstSNI'	=> 'TACGTA 3',
 'BstUI'	=> 'CGCG 2',
 'Bst2UI'	=> 'CCWGG 2',
 'BstXI'	=> 'CCANNNNNNTGG 8',
 'BstX2I'	=> 'RGATCY 1',
 'BstYI'	=> 'RGATCY 1',
 'BstZI'	=> 'CGGCCG 1',
 'BstZ17I'	=> 'GTATAC 3',
 'Bsu15I'	=> 'ATCGAT 2',
 'Bsu36I'	=> 'CCTNAGG 2',
 'BsuRI'	=> 'GGCC 2',
 'BsuTUI'	=> 'ATCGAT 2',
 'BtgI'	=> 'CCRYGG 1',
 'BthCI'	=> 'GCNGC 4',
 'Cac8I'	=> 'GCNNGC 3',
 'CaiI'	=> 'CAGNNNCTG 6',
 'CauII'	=> 'CCSGG 2',
 'CciNI'	=> 'GCGGCCGC 2',
 'CelII'	=> 'GCTNAGC 2',
 'CfoI'	=> 'GCGC 3',
 'CfrI'	=> 'YGGCCR 1',
 'Cfr9I'	=> 'CCCGGG 1',
 'Cfr10I'	=> 'RCCGGY 1',
 'Cfr13I'	=> 'GGNCC 1',
 'Cfr42I'	=> 'CCGCGG 4',
 'ChaI'	=> 'GATC 4',
 'ClaI'	=> 'ATCGAT 2',
 'CpoI'	=> 'CGGWCCG 2',
 'CspI'	=> 'CGGWCCG 2',
 'Csp6I'	=> 'GTAC 1',
 'Csp45I'	=> 'TTCGAA 2',
 'CspAI'	=> 'ACCGGT 1',
 'CviAII'	=> 'CATG 1',
 'CviJI'	=> 'RGCY 2',
 'CviRI'	=> 'TGCA 2',
 'CviTI'	=> 'RGCY 2',
 'CvnI'	=> 'CCTNAGG 2',
 'DdeI'	=> 'CTNAG 1',
 'DpnI'	=> 'GATC 2',
 'DpnII'	=> 'GATC 0',
 'DraI'	=> 'TTTAAA 3',
 'DraII'	=> 'RGGNCCY 2',
 'DraIII'	=> 'CACNNNGTG 6',
 'DrdI'	=> 'GACNNNNNNGTC 7',
 'DsaI'	=> 'CCRYGG 1',
 'DseDI'	=> 'GACNNNNNNGTC 7',
 'EaeI'	=> 'YGGCCR 1',
 'EagI'	=> 'CGGCCG 1',
 'Eam1105I'	=> 'GACNNNNNGTC 6',
 'Ecl136II'	=> 'GAGCTC 3',
 'EclHKI'	=> 'GACNNNNNGTC 6',
 'EclXI'	=> 'CGGCCG 1',
 'Eco24I'	=> 'GRGCYC 5',
 'Eco32I'	=> 'GATATC 3',
 'Eco47I'	=> 'GGWCC 1',
 'Eco47III'	=> 'AGCGCT 3',
 'Eco52I'	=> 'CGGCCG 1',
 'Eco72I'	=> 'CACGTG 3',
 'Eco81I'	=> 'CCTNAGG 2',
 'Eco88I'	=> 'CYCGRG 1',
 'Eco91I'	=> 'GGTNACC 1',
 'Eco105I'	=> 'TACGTA 3',
 'Eco130I'	=> 'CCWWGG 1',
 'Eco147I'	=> 'AGGCCT 3',
 'EcoHI'	=> 'CCSGG 0',
 'EcoICRI'	=> 'GAGCTC 3',
 'EcoNI'	=> 'CCTNNNNNAGG 5',
 'EcoO65I'	=> 'GGTNACC 1',
 'EcoO109I'	=> 'RGGNCCY 2',
 'EcoRI'	=> 'GAATTC 1',
 'EcoRII'	=> 'CCWGG 0',
 'EcoRV'	=> 'GATATC 3',
 'EcoT14I'	=> 'CCWWGG 1',
 'EcoT22I'	=> 'ATGCAT 5',
 'EcoT38I'	=> 'GRGCYC 5',
 'EgeI'	=> 'GGCGCC 3',
 'EheI'	=> 'GGCGCC 3',
 'ErhI'	=> 'CCWWGG 1',
 'EsaBC3I'	=> 'TCGA 2',
 'EspI'	=> 'GCTNAGC 2',
 'FatI'	=> 'CATG 0',
 'FauNDI'	=> 'CATATG 2',
 'FbaI'	=> 'TGATCA 1',
 'FblI'	=> 'GTMKAC 2',
 'FmuI'	=> 'GGNCC 4',
 'FnuDII'	=> 'CGCG 2',
 'Fnu4HI'	=> 'GCNGC 2',
 'FriOI'	=> 'GRGCYC 5',
 'FseI'	=> 'GGCCGGCC 6',
 'FspI'	=> 'TGCGCA 3',
 'FspAI'	=> 'RTGCGCAY 4',
 'Fsp4HI'	=> 'GCNGC 2',
 'FunI'	=> 'AGCGCT 3',
 'FunII'	=> 'GAATTC 1',
 'HaeI'	=> 'WGGCCW 3',
 'HaeII'	=> 'RGCGCY 5',
 'HaeIII'	=> 'GGCC 2',
 'HapII'	=> 'CCGG 1',
 'HgiAI'	=> 'GWGCWC 5',
 'HgiCI'	=> 'GGYRCC 1',
 'HgiJII'	=> 'GRGCYC 5',
 'HhaI'	=> 'GCGC 3',
 'Hin1I'	=> 'GRCGYC 2',
 'Hin6I'	=> 'GCGC 1',
 'HinP1I'	=> 'GCGC 1',
 'HincII'	=> 'GTYRAC 3',
 'HindII'	=> 'GTYRAC 3',
 'HindIII'	=> 'AAGCTT 1',
 'HinfI'	=> 'GANTC 1',
 'HpaI'	=> 'GTTAAC 3',
 'HpaII'	=> 'CCGG 1',
 'Hpy8I'	=> 'GTNNAC 3',
 'Hpy99I'	=> 'CGWCG 5',
 'Hpy178III'	=> 'TCNNGA 2',
 'Hpy188I'	=> 'TCNGA 3',
 'Hpy188III'	=> 'TCNNGA 2',
 'HpyCH4I'	=> 'CATG 3',
 'HpyCH4III'	=> 'ACNGT 3',
 'HpyCH4IV'	=> 'ACGT 1',
 'HpyCH4V'	=> 'TGCA 2',
 'HpyF10VI'	=> 'GCNNNNNNNGC 8',
 'Hsp92I'	=> 'GRCGYC 2',
 'Hsp92II'	=> 'CATG 4',
 'HspAI'	=> 'GCGC 1',
 'ItaI'	=> 'GCNGC 2',
 'KasI'	=> 'GGCGCC 1',
 'KpnI'	=> 'GGTACC 5',
 'Kpn2I'	=> 'TCCGGA 1',
 'KspI'	=> 'CCGCGG 4',
 'Ksp22I'	=> 'TGATCA 1',
 'KspAI'	=> 'GTTAAC 3',
 'Kzo9I'	=> 'GATC 0',
 'LpnI'	=> 'RGCGCY 3',
 'LspI'	=> 'TTCGAA 2',
 'MabI'	=> 'ACCWGGT 1',
 'MaeI'	=> 'CTAG 1',
 'MaeII'	=> 'ACGT 1',
 'MaeIII'	=> 'GTNAC 0',
 'MamI'	=> 'GATNNNNATC 5',
 'MboI'	=> 'GATC 0',
 'McrI'	=> 'CGRYCG 4',
 'MfeI'	=> 'CAATTG 1',
 'MflI'	=> 'RGATCY 1',
 'MhlI'	=> 'GDGCHC 5',
 'MlsI'	=> 'TGGCCA 3',
 'MluI'	=> 'ACGCGT 1',
 'MluNI'	=> 'TGGCCA 3',
 'Mly113I'	=> 'GGCGCC 2',
 'Mph1103I'	=> 'ATGCAT 5',
 'MroI'	=> 'TCCGGA 1',
 'MroNI'	=> 'GCCGGC 1',
 'MroXI'	=> 'GAANNNNTTC 5',
 'MscI'	=> 'TGGCCA 3',
 'MseI'	=> 'TTAA 1',
 'MslI'	=> 'CAYNNNNRTG 5',
 'MspI'	=> 'CCGG 1',
 'Msp20I'	=> 'TGGCCA 3',
 'MspA1I'	=> 'CMGCKG 3',
 'MspCI'	=> 'CTTAAG 1',
 'MspR9I'	=> 'CCNGG 2',
 'MssI'	=> 'GTTTAAAC 4',
 'MstI'	=> 'TGCGCA 3',
 'MunI'	=> 'CAATTG 1',
 'MvaI'	=> 'CCWGG 2',
 'MvnI'	=> 'CGCG 2',
 'MwoI'	=> 'GCNNNNNNNGC 7',
 'NaeI'	=> 'GCCGGC 3',
 'NarI'	=> 'GGCGCC 2',
 'NciI'	=> 'CCSGG 2',
 'NcoI'	=> 'CCATGG 1',
 'NdeI'	=> 'CATATG 2',
 'NdeII'	=> 'GATC 0',
 'NgoAIV'	=> 'GCCGGC 1',
 'NgoMIV'	=> 'GCCGGC 1',
 'NheI'	=> 'GCTAGC 1',
 'NlaIII'	=> 'CATG 4',
 'NlaIV'	=> 'GGNNCC 3',
 'Nli3877I'	=> 'CYCGRG 5',
 'NmuCI'	=> 'GTSAC 0',
 'NotI'	=> 'GCGGCCGC 2',
 'NruI'	=> 'TCGCGA 3',
 'NruGI'	=> 'GACNNNNNGTC 6',
 'NsbI'	=> 'TGCGCA 3',
 'NsiI'	=> 'ATGCAT 5',
 'NspI'	=> 'RCATGY 5',
 'NspIII'	=> 'CYCGRG 1',
 'NspV'	=> 'TTCGAA 2',
 'NspBII'	=> 'CMGCKG 3',
 'OliI'	=> 'CACNNNNGTG 5',
 'PacI'	=> 'TTAATTAA 5',
 'PaeI'	=> 'GCATGC 5',
 'PaeR7I'	=> 'CTCGAG 1',
 'PagI'	=> 'TCATGA 1',
 'PalI'	=> 'GGCC 2',
 'PauI'	=> 'GCGCGC 1',
 'PceI'	=> 'AGGCCT 3',
 'PciI'	=> 'ACATGT 1',
 'PdiI'	=> 'GCCGGC 3',
 'PdmI'	=> 'GAANNNNTTC 5',
 'Pfl23II'	=> 'CGTACG 1',
 'PflBI'	=> 'CCANNNNNTGG 7',
 'PflFI'	=> 'GACNNNGTC 4',
 'PflMI'	=> 'CCANNNNNTGG 7',
 'PfoI'	=> 'TCCNGGA 1',
 'PinAI'	=> 'ACCGGT 1',
 'Ple19I'	=> 'CGATCG 4',
 'PmaCI'	=> 'CACGTG 3',
 'PmeI'	=> 'GTTTAAAC 4',
 'PmlI'	=> 'CACGTG 3',
 'Ppu10I'	=> 'ATGCAT 1',
 'PpuMI'	=> 'RGGWCCY 2',
 'PpuXI'	=> 'RGGWCCY 2',
 'PshAI'	=> 'GACNNNNGTC 5',
 'PshBI'	=> 'ATTAAT 2',
 'PsiI'	=> 'TTATAA 3',
 'Psp03I'	=> 'GGWCC 4',
 'Psp5II'	=> 'RGGWCCY 2',
 'Psp6I'	=> 'CCWGG 0',
 'Psp1406I'	=> 'AACGTT 2',
 'PspAI'	=> 'CCCGGG 1',
 'Psp124BI'	=> 'GAGCTC 5',
 'PspEI'	=> 'GGTNACC 1',
 'PspGI'	=> 'CCWGG 0',
 'PspLI'	=> 'CGTACG 1',
 'PspN4I'	=> 'GGNNCC 3',
 'PspOMI'	=> 'GGGCCC 1',
 'PspPI'	=> 'GGNCC 1',
 'PspPPI'	=> 'RGGWCCY 2',
 'PssI'	=> 'RGGNCCY 5',
 'PstI'	=> 'CTGCAG 5',
 'PsuI'	=> 'RGATCY 1',
 'PsyI'	=> 'GACNNNGTC 4',
 'PvuI'	=> 'CGATCG 4',
 'PvuII'	=> 'CAGCTG 3',
 'RcaI'	=> 'TCATGA 1',
 'RsaI'	=> 'GTAC 2',
 'RsrII'	=> 'CGGWCCG 2',
 'Rsr2I'	=> 'CGGWCCG 2',
 'SacI'	=> 'GAGCTC 5',
 'SacII'	=> 'CCGCGG 4',
 'SalI'	=> 'GTCGAC 1',
 'SanDI'	=> 'GGGWCCC 2',
 'SatI'	=> 'GCNGC 2',
 'SauI'	=> 'CCTNAGG 2',
 'Sau96I'	=> 'GGNCC 1',
 'Sau3AI'	=> 'GATC 0',
 'SbfI'	=> 'CCTGCAGG 6',
 'ScaI'	=> 'AGTACT 3',
 'SciI'	=> 'CTCGAG 3',
 'ScrFI'	=> 'CCNGG 2',
 'SdaI'	=> 'CCTGCAGG 6',
 'SduI'	=> 'GDGCHC 5',
 'SecI'	=> 'CCNNGG 1',
 'SelI'	=> 'CGCG 0',
 'SexAI'	=> 'ACCWGGT 1',
 'SfcI'	=> 'CTRYAG 1',
 'SfeI'	=> 'CTRYAG 1',
 'SfiI'	=> 'GGCCNNNNNGGCC 8',
 'SfoI'	=> 'GGCGCC 3',
 'Sfr274I'	=> 'CTCGAG 1',
 'Sfr303I'	=> 'CCGCGG 4',
 'SfuI'	=> 'TTCGAA 2',
 'SgfI'	=> 'GCGATCGC 5',
 'SgrAI'	=> 'CRCCGGYG 2',
 'SgrBI'	=> 'CCGCGG 4',
 'SinI'	=> 'GGWCC 1',
 'SlaI'	=> 'CTCGAG 1',
 'SmaI'	=> 'CCCGGG 3',
 'SmiI'	=> 'ATTTAAAT 4',
 'SmiMI'	=> 'CAYNNNNRTG 5',
 'SmlI'	=> 'CTYRAG 1',
 'SnaBI'	=> 'TACGTA 3',
 'SpaHI'	=> 'GCATGC 5',
 'SpeI'	=> 'ACTAGT 1',
 'SphI'	=> 'GCATGC 5',
 'SplI'	=> 'CGTACG 1',
 'SrfI'	=> 'GCCCGGGC 4',
 'Sse9I'	=> 'AATT 0',
 'Sse232I'	=> 'CGCCGGCG 2',
 'Sse8387I'	=> 'CCTGCAGG 6',
 'Sse8647I'	=> 'AGGWCCT 2',
 'SseBI'	=> 'AGGCCT 3',
 'SspI'	=> 'AATATT 3',
 'SspBI'	=> 'TGTACA 1',
 'SstI'	=> 'GAGCTC 5',
 'SstII'	=> 'CCGCGG 4',
 'StuI'	=> 'AGGCCT 3',
 'StyI'	=> 'CCWWGG 1',
 'SunI'	=> 'CGTACG 1',
 'SwaI'	=> 'ATTTAAAT 4',
 'TaaI'	=> 'ACNGT 3',
 'TaiI'	=> 'ACGT 4',
 'TaqI'	=> 'TCGA 1',
 'TasI'	=> 'AATT 0',
 'TatI'	=> 'WGTACW 1',
 'TauI'	=> 'GCSGC 4',
 'TelI'	=> 'GACNNNGTC 4',
 'TfiI'	=> 'GAWTC 1',
 'ThaI'	=> 'CGCG 2',
 'TliI'	=> 'CTCGAG 1',
 'Tru1I'	=> 'TTAA 1',
 'Tru9I'	=> 'TTAA 1',
 'TscI'	=> 'ACGT 4',
 'TseI'	=> 'GCWGC 1',
 'Tsp45I'	=> 'GTSAC 0',
 'Tsp509I'	=> 'AATT 0',
 'Tsp4CI'	=> 'ACNGT 3',
 'TspEI'	=> 'AATT 0',
 'Tth111I'	=> 'GACNNNGTC 4',
 'TthHB8I'	=> 'TCGA 1',
 'UnbI'	=> 'GGNCC 0',
 'Van91I'	=> 'CCANNNNNTGG 7',
 'Vha464I'	=> 'CTTAAG 1',
 'VneI'	=> 'GTGCAC 1',
 'VpaK11AI'	=> 'GGWCC 0',
 'VpaK11BI'	=> 'GGWCC 1',
 'VspI'	=> 'ATTAAT 2',
 'XagI'	=> 'CCTNNNNNAGG 5',
 'XapI'	=> 'RAATTY 1',
 'XbaI'	=> 'TCTAGA 1',
 'XceI'	=> 'RCATGY 5',
 'XcmI'	=> 'CCANNNNNNNNNTGG 8',
 'XhoI'	=> 'CTCGAG 1',
 'XhoII'	=> 'RGATCY 1',
 'XmaI'	=> 'CCCGGG 1',
 'XmaIII'	=> 'CGGCCG 1',
 'XmaCI'	=> 'CCCGGG 1',
 'XmaJI'	=> 'CCTAGG 1',
 'XmiI'	=> 'GTMKAC 2',
 'XmnI'	=> 'GAANNNNTTC 5',
 'XspI'	=> 'CTAG 1',
 'ZhoI'	=> 'ATCGAT 2',
 'ZraI'	=> 'GACGTC 3',
 'Zsp2I'	=> 'ATGCAT 5',
);

@RE_available = sort keys %RE;


=head1 new

 Title     : new
 Purpose   : Initializes the RestrictionEnzyme object and calls
           : superclass constructor last (Bio:Seq.pm).
 Returns   : n/a
 Argument  : Parameters passed to new()
 Comments  : A RestrictionEnzyme object manages its recognition sequence
           : as a Bio::PrimarySeq object.

See Also   : L<_make_custom>(), L<_make_standard>(), L<Bio::PrimarySeq.pm::_initialize()>

=cut

#---------------
sub new {
#---------------
    my($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    $self->warning("Use of Bio::Tools::RestrictionEnzyme is deprecated".
                   "Use Bio::Restriction classes instead");
    my ($name,$make) = $self->_rearrange([qw(NAME MAKE)],@args);

    $name && $self->name($name);
    my %data;
    if(defined $make && $make eq 'custom') {
	%data = $self->_make_custom($name); 
    } else {
	%data = $self->_make_standard($name);
    }
    $self->{'_seq'} = new Bio::PrimarySeq(%data, 
				   -VERBOSE =>$self->verbose,
 				   -alphabet => 'dna',
				   );
    return $self;
}


#=head1 _make_standard
#
# Title     : _make_standard
# Usage     : n/a; automatically called by _initialize()
# Purpose   : Permits custom RE object construction from name.
#	    : 'EcoRI'.
# Returns   : Hash containing named parameters for Bio::PrimarySeq.pm constructor.
# Argument  : String containing string with special syntax.
# Throws    : Exception if the requested enzyme name is unavailable.
#	    : NOTE: Case sensitive.
#
#See Also   : L<Bio::PrimarySeq::_initialize()|Bio::PrimarySeq>, L<_make_custom()|_make_custom>
#
#=cut

#------------------
sub _make_standard {
#------------------
    my($self, $name) = @_;

    $name =~ s/^\s+|\s+$//g;
 
    $self->is_available($name) || 
	$self->throw("Unavailable or undefined enzyme: $name (Note: CASE SENSITIVE)\n" .
		     "Currently available enzymes: \n@RE_available\n");

    my @data = split( ' ', $RE{$name});
    my (%dat);
    $dat{-SEQ} = $data[0];
    $dat{-NAME} = $dat{-ID}= $name;    
    $self->{'_cuts_after'} = $data[1];

    return %dat;
}


#=head1 _make_custom
#
# Title     : _make_custom
# Usage     : n/a; automatically called by _initialize()
# Purpose   : Permits custom RE object construction from strings 
#	    : such as 'EcoRI--G^AATTC' as the name of the enzyme.
# Returns   : Hash containing named parameters for Bio::PrimarySeq.pm constructor.
# Argument  : String containing string with special syntax.
# Throws    : Exception if the string has bad syntax.
#	    : Warning if the string did not specify cut position.
#	    :         Places cut site after 5'-most position.
#
#See Also   : L<Bio::PrimarySeq::_initialize()|Bio::PrimarySeq>
#
#=cut

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
			 $self->warn("Unknown cut position for $parts[0]. Assuming position 0\n" . 
				     "Use carat to specify cut position (e.g., G^AATTC)"); }
    $self->{'_cuts_after'} =  $cut_index;

    ## Save the recognition sequence after removing the '^'
    $parts[1] =~ s/\^//g;
    $dat{-SEQ} = $parts[1];
    return %dat;
}
    

=head1 cuts_after

 Title     : cuts_after
 Usage     : $num = $re->cuts_after();
 Purpose   : Sets/Gets an integer indicating the position of cleavage 
           : relative to the 5' end of the recognition sequence.
 Returns   : Integer
 Argument  : Integer (optional)
 Throws    : Exception if argument is non-numeric.
 Access    : Public
 Comments  : This method is only needed to change the cuts at
           : position. This data is automatically set during
           : construction.

See Also   : L<_make_standard()|_make_standard>, L<_make_custom()|_make_custom>

=cut

#'
#---------------
sub cuts_after { 
#---------------
    my $self = shift; 
    if(@_) { my $num = shift;
	     if($num == 0 and $num ne '0') {
		 $self->throw("The cuts_after position be an integer ($num)");
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
 Comments  : If you want a simple string representing the site without 
             any '^', use the string() method.

See Also   : L<string()|string>

=cut

#---------
sub site {
#---------
    my $self = shift;
    my $seq = $self->seq;
    my $cuts_after = $self->cuts_after;
    if($cuts_after > 0) {
	if( $cuts_after >= $seq->length) {
	    return $seq->seq.'^';
	} else { 
	    return $seq->subseq(1, $self->cuts_after).'^'.$seq->subseq($self->cuts_after+1, $seq->length); 
	}
    } else {
        return $seq->seq;
    }
}


=head1 seq

 Title     : seq
 Usage     : $re->seq();
 Purpose   : Get the Bio::PrimarySeq.pm-derived object representing 
           : the recognition sequence
 Returns   : String
 Argument  : n/a
 Throws    : n/a

See Also   : L<string()|string>, L<revcom()|revcom>

=cut

#---------
sub seq    {  my $self = shift; $self->{'_seq'}; }
#---------



=head1 string

 Title     : string
 Usage     : $re->string();
 Purpose   : Get a string representing the recognition sequence.
 Returns   : String. Does NOT contain a  '^' representing the cut location
             as returned by the site() method
 Argument  : n/a
 Throws    : n/a
 Comments  : Delegates to the Bio::PrimarySeq-derived object.

See Also   : L<seq()|seq>, L<site()|site>, L<revcom()|revcom>

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
 Comments  : Delegates to the Bio::PrimarySeq.pm-derived object, but needs to
             get out the string from it, as now Bio::PrimarySeq->revcom makes a
             Bio::PrimarySeq object

See Also   : L<seq()|seq>, L<string()|string>

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
 Argument  : Reference to a Bio::PrimarySeq.pm-derived object.
 Throws    : Exception if argument is not an object.
           : (Does not yet verify that it is derived from Bio::PrimarySeq.pm.)
 Comments  : Strategy relies on Perl's built-in split() function.
           : Since split removes the recognition pattern, the resulting
           : fragments are repaired after split()-ing.
           : A side-effect of this is that for sites with ambiguous
           : recognition sequence (i.e., containing N), the fragments
           : will contain ambiguity characters instead of AGCT.
           :
           : There is currently no support for partial digestions.
           : There is currently no support for circular sequences.
           : (This should just involve merging the first and last frag
           : if $seqObj->is_circular returns true).

=cut

#'
#-------------
sub cut_seq {
#-------------
    my( $self, $seqObj) = @_;
    if( !ref($seqObj) || 
	! $seqObj->isa('Bio::PrimarySeqI') ) {
	$self->throw( "Can't cut sequence. Missing or invalid object".
		      "seqObj: $seqObj");
    }

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

    $self->debug("3' site: $site_3prime_seq\n5' site: $site_5prime_seq\n");

    my(@re_frags);
    my $seq = uc $self->_expanded_string;

    if(!$self->palindromic and $self->name ne 'N') {
	my $revseq = $self->_expanded_string( $reSeq->revcom->seq() );
	$seq .= '|'.uc($revseq);
    }
    $self->debug(__PACKAGE__, ": site seq: $seq\n");
    $self->debug(__PACKAGE__, ": splitting ", $reSeq->seq, "\n");
    @re_frags = split(/$seq/i, $seqObj->seq);

    $self->debug(__PACKAGE__, ": cut_seq, ", scalar(@re_frags), " fragments.\n");

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
 Argument  : Reference to a Bio::PrimarySeqI-derived sequence object.
 Throws    : n/a
 Comments  : 

=cut

#-----------------
sub cut_locations {
#-----------------
    my($self, $seqobj) = @_;

    my $site = $self->_expanded_string;
    my $seq = $seqobj->seq;
    study($seq);
    my @locations;
    while( $seq =~ /($site)/ig ) {
        # $` is preceding string before pattern so length returns position
	push @locations, length($`); 	
    }
    return \@locations;
}    

# Purpose : Expand nucleotide ambiguity codes to their representative letters
# Argument: (optional) the string to be expanded. If not supplied, used
#           the string returned by $self->string().
# Returns : String
sub _expanded_string {
    my ($self, $str) = @_;
    
    $str ||= $self->string;

    if( $self->name ne 'N' ) {
        $str =~ s/N|X/\./g;
        $str =~ s/R/\[AG\]/g;
        $str =~ s/Y/\[CT\]/g;
        $str =~ s/S/\[GC\]/g;
        $str =~ s/W/\[AT\]/g;
        $str =~ s/M/\[AC\]/g;
        $str =~ s/K/\[TG\]/g;
        $str =~ s/B/\[CGT\]/g;
        $str =~ s/D/\[AGT\]/g;
        $str =~ s/H/\[ACT\]/g;
        $str =~ s/V/\[ACG\]/g;
    }
    return $str;
}


=head1 annotate_seq

 Title     : annotate_seq
 Usage     : $re->annotate_seq(<sequence_object>);
 Purpose   : Identify the location of the recognition site(s) within
           : an input sequence. Uses HTML.
 Example   : $annot_seq = $re->annotate_seq($seqObj);
 Returns   : String containing the annotated sequence.
 Argument  : Reference to a Bio::PrimarySeq.pm-derived sequence object.
 Throws    : n/a
 Comments  : The annotated sequence must be viewed with a web
           : browser to see the location(s) of the recognition site(s).

=cut

#-----------------
sub annotate_seq {
#-----------------
    my($self, $seqObj) = @_;

    my $site = $self->_expanded_string;
    my $seq = $seqObj->seq;

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

See Also   : L<available_list()|available_list>

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
    $self->warn('available() is deprecated; use is_available() instead');
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

See Also   : L<is_available()|is_available>

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
