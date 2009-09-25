package Bio::Tools::SeqPattern::Backtranslate;
use strict;
use warnings;

use base qw(Bio::Root::Root);
use base qw(Exporter);

=head1 NAME

Bio::Tools::SeqPattern::Backtranslate

=head1 DESCRIPTION

This module should not be used directly. It provides helper methods to
Bio::Tools::SeqPattern to reverse translate protein patterns.

=cut

use Bio::Seq;
use Bio::Tools::CodonTable;

use List::MoreUtils qw(uniq);
use Carp qw(croak);

our @EXPORT_OK = qw(_reverse_translate_motif);

our @EXPORT = @EXPORT_OK;

sub _reverse_translate_motif {
   # Main subroutine. It takes a Profam-like motif and returns its
   # reverse translation using degenerate codons.

   # Steps:
   # 1. Tokenize, then parse tokens.
   # 2. Reverse translate each token type.
   # 3. Join tokens in original order. Return the resulting string.

   my $motif = shift;

   $motif =~ s/\./X/g;
   $motif = uc $motif;

   ### 1. Tokenize, parse the motif.
   my ( $ordered, $classified ) = _parse_motif($motif);

   ### 2. Reverse translate each token type.
   # Reverse translate the plain (unambiguous) tokens.
   my $ct = Bio::Tools::CodonTable->new;
   foreach my $seq ( @{ $classified->{plain} } ) {
      my $seqO
          = Bio::Seq->new( -seq => $$seq, -alphabet => 'protein' );
      $$seq = $ct->reverse_translate_all($seqO);
   }

   # Reverse translate the ambiguous tokens.
   foreach my $token ( @{ $classified->{ambiguous} } ) {
      my ($aas) = $$token =~ m(([A-Za-z\.]+));
      my @codons_to_contract;

      foreach my $residue ( split '', $aas ) {
         push @codons_to_contract, $ct->revtranslate($residue);
      }

      my $ambiguous_codon = _contract_codons(@codons_to_contract);
      $$token             = $ambiguous_codon;
   }

   # Reverse translate the negated residues.
   foreach my $token ( @{ $classified->{negated} } ) {
      my ($aas) = $$token =~ m(([A-Za-z\.]+));
      my $ambiguous_codon = _negated_aas_to_codon($aas);
      $$token = $ambiguous_codon;
   }

   ### 3. Join the profile back from its tokens.
   return join '', map {$$_} @{$ordered};

}

sub _parse_motif {
   # Profam-like motif parser. It takes the pattern as a string, and
   # returns two data structures that contain the tokens, organized
   # by order of appearance in the pattern (first return value) and by
   # category (second return value).

   my $motif = shift;
   my $parser  = _tokenize_motif($motif);
   my ( %tokens, @tokens );

   while ( my $token = $parser->() ) {
      croak ("Unknown syntax token: <", $token->[1], ">")
          if ( $token->[0] eq 'UNKNOWN' );
      push @{ $tokens{ $token->[0] } }, \$token->[1];
      push @tokens, \$token->[1];
   }
   return ( \@tokens, \%tokens );
}

sub _tokenize_motif {

   # Return a tokenizer iterator that sequentially recognizes and
   # returns each token in the input pattern.
   # Examples of each token type:

   # ambiguous: a position with more than one possible residue.
   #     eg. [ALEP]
   # negated: a position in which some residues are excluded.
   #     eg. [^WY]
   # plain: a common sequence of residues. One position, one residue.
   #     eg. MAAEIK
   # open_par, close_par: tags surrounding a motif that is repeated
   #     a certain number of times.
   #     eg. (...){3}

   my $target = shift;
   return sub {
      return [ 'ambiguous', $1 ]
          if $target =~ /\G (\[[A-Za-z\.]+\])     /gcx;
      return [ 'negated',   $1 ]
          if $target =~ /\G (\[\^[A-Za-z\.]+\])   /gcx;
      return [ 'plain',     $1 ]
          if $target =~ /\G ([A-Za-z\.]+)         /gcx;
      return [ 'open_par',  $1 ]
          if $target =~ /\G (\()                  /gcx;
      return [ 'close_par', $1 ]
          if $target =~ /\G (\)[\{\d+[,\d+]*\}]*) /gcx;
      return [ 'UNKNOWN',   $1 ]
          if $target =~ /\G (.)                   /gcx;
      return;
   };
}

sub _contract_codons {

   # Take a list of codons, return an ambiguous codon.
   my @codons = map { uc $_ } @_;

   my @by_letter = ( [], [], [], );
   my $ambiguous_codon;
   foreach my $codon (@codons) {
      my @letters = split '', $codon;
      for my $i ( 0 .. 2 ) {
         push @{ $by_letter[$i] }, $letters[$i];
      }
   }
   for my $i ( 0 .. 2 ) {
      $ambiguous_codon
          .= _convert( 'dna', _uniq_string( @{ $by_letter[$i] } ) );
   }
   return $ambiguous_codon;
}

sub _expand_codon {

   # Given a degenerate codon, return a list with all its
   # constituents. Takes a three-letter string (codon) as
   # input, returns a list with three-letter scalars.

   my $codon = shift;
   die "Wrong codon length!\n" if length $codon != 3;


   my ( @codons, @return_bases );
   my @orig_bases = split '', $codon;

   for my $i ( 0 .. 2 ) {

      # from each redundant base, create a list with all their
      # components (e.g., N -> (A, C, G, T) );
      my @components = split '', _convert('dna', $orig_bases[$i] );
      $orig_bases[$i] = [@components];
   }

   # Combine all the bases of each of the three positions of the
   # codons, and build the return list.
   for my $i ( @{ $orig_bases[0] } ) {
      for my $j ( @{ $orig_bases[1] } ) {
         for my $k ( @{ $orig_bases[2] } ) {
            push @return_bases, $i . $j . $k;
         }
      }
   }
   return @return_bases;
}

{
   my %convert;

   sub _convert {
      # Interconvert between redundant and non-redundant protein and
      # dna alphabets. Takes an alphabet (protein or dna) and a string
      # with the letter, and returns its equivalent in
      # redundant/non-redundant alphabet. Example ACTG -> N.

      my ($alphabet, $letter) = @_;
      unless (
         $alphabet and $alphabet =~ /^dna$|^protein$/i
         and $letter and length $letter <= 4
      ) { croak "Wrong arguments!\n"; }

      unless (%convert) {
         %convert = (
            'dna' => {
               qw(N ACGT B CGT D AGT H ACT V ACG K GT
                   M AC R AG S CG W AT Y CT A A C C T T G G)
            },
            'protein' => {
               '.' => 'ACDEFGHIJKLMNOPQRSTUVWY',
               X   => 'ACDEFGHIJKLMNOPQRSTUVWY',
               Z   => 'QE',
               B   => 'ND',
            },
         );

         # Make %convert hash key/value agnostic.
         foreach my $alphabet ( keys %convert ) {
            map { $convert{$alphabet}->{ $convert{$alphabet}{$_} } = $_ }
                keys %{ $convert{$alphabet} };
         }
      }

      return $convert{$alphabet}{$letter};
   }

}

sub _uniq_string {
   # Takes a list of letters and returns an alphabetically sorted
   # list with unique elements.

   my @letters = @_;
   return join '', sort { $a cmp $b } uniq @letters;
}

{
   my ( @codon_library, $ct );

   sub _negated_aas_to_codon {

      # Given a string of residues, returns a degenerate codon that will
      # not be translated into any of them, while maximizing degeneracy
      # (ie, it tries to also translate into as many residues as possible).

      # This functionality is required for reverse translating profiles
      # that contain negative patterns: [^X]. This means that the current
      # position should not contain aminoacid X, but can have any of the
      # others. The reverse translated nucleotide sequence should
      # reflect this.

      # Approach: construct a list of all possible codons, incluiding all
      # degenerate bases. This is an array of 15x15x15 = 3375 elements.
      # Order them by descendent "degeneracy".
      # Return the first one whose expansion in 4-lettered codons
      # doesn't contain a codon that translates into any of the
      # non-wanted residues.
      
      # * Since this takes some time, I presorted them and saved them.
      # Reading them from a file takes a fraction of the time that it taes
      # to re-sort them every time the application is launched.

      my $aas_to_avoid = shift;

      # Initialize reusable variables if it's the first time the sub
      # is called.
      unless (@codon_library) {
         while (<DATA>) { chomp; push @codon_library, split ' ', $_ }
      }
      unless ($ct) { $ct = Bio::Tools::CodonTable->new; }

      # Reverse translate the unwanted aminoacids to unwanted codons.
      my @unwanted_codons;
      foreach my $aa ( split '', $aas_to_avoid ) {
         push @unwanted_codons, $ct->revtranslate($aa);
      }

      foreach my $degenerate_codon (@codon_library) {
         my @codons  = _expand_codon($degenerate_codon);
         my $success = 1;

         foreach my $unwanted (@unwanted_codons) {
            if ( grep { uc $unwanted eq $_ } @codons ) {
               $success = 0;
            }
         }

         if ($success) { return $degenerate_codon }
      }
   }

}

1;

=head1 COPYRIGHT & LICENSE

Copyright 2009 Bruno Vecchi, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

__DATA__
NNN NNB NND NNH NNV NBN NDN NHN NVN BNN DNN HNN VNN NBB NBD NBH NBV NDB NDD
NDH NDV NHB NHD NHH NHV NVB NVD NVH NVV BNB BND BNH BNV BBN BDN BHN BVN DNB
DND DNH DNV DBN DDN DHN DVN HNB HND HNH HNV HBN HDN HHN HVN VNB VND VNH VNV
VBN VDN VHN VVN NNK NNM NNR NNS NNW NNY NKN NMN NRN NSN NWN NYN KNN MNN RNN
SNN WNN YNN BBB BBD BBH BBV BDB BDD BDH BDV BHB BHD BHH BHV BVB BVD BVH BVV
DBB DBD DBH DBV DDB DDD DDH DDV DHB DHD DHH DHV DVB DVD DVH DVV HBB HBD HBH
HBV HDB HDD HDH HDV HHB HHD HHH HHV HVB HVD HVH HVV VBB VBD VBH VBV VDB VDD
VDH VDV VHB VHD VHH VHV VVB VVD VVH VVV NBK NBM NBR NBS NBW NBY NDK NDM NDR
NDS NDW NDY NHK NHM NHR NHS NHW NHY NVK NVM NVR NVS NVW NVY NKB NKD NKH NKV
NMB NMD NMH NMV NRB NRD NRH NRV NSB NSD NSH NSV NWB NWD NWH NWV NYB NYD NYH
NYV BNK BNM BNR BNS BNW BNY BKN BMN BRN BSN BWN BYN DNK DNM DNR DNS DNW DNY
DKN DMN DRN DSN DWN DYN HNK HNM HNR HNS HNW HNY HKN HMN HRN HSN HWN HYN VNK
VNM VNR VNS VNW VNY VKN VMN VRN VSN VWN VYN KNB KND KNH KNV KBN KDN KHN KVN
MNB MND MNH MNV MBN MDN MHN MVN RNB RND RNH RNV RBN RDN RHN RVN SNB SND SNH
SNV SBN SDN SHN SVN WNB WND WNH WNV WBN WDN WHN WVN YNB YND YNH YNV YBN YDN
YHN YVN BBK BBM BBR BBS BBW BBY BDK BDM BDR BDS BDW BDY BHK BHM BHR BHS BHW
BHY BVK BVM BVR BVS BVW BVY BKB BKD BKH BKV BMB BMD BMH BMV BRB BRD BRH BRV
BSB BSD BSH BSV BWB BWD BWH BWV BYB BYD BYH BYV DBK DBM DBR DBS DBW DBY DDK
DDM DDR DDS DDW DDY DHK DHM DHR DHS DHW DHY DVK DVM DVR DVS DVW DVY DKB DKD
DKH DKV DMB DMD DMH DMV DRB DRD DRH DRV DSB DSD DSH DSV DWB DWD DWH DWV DYB
DYD DYH DYV HBK HBM HBR HBS HBW HBY HDK HDM HDR HDS HDW HDY HHK HHM HHR HHS
HHW HHY HVK HVM HVR HVS HVW HVY HKB HKD HKH HKV HMB HMD HMH HMV HRB HRD HRH
HRV HSB HSD HSH HSV HWB HWD HWH HWV HYB HYD HYH HYV VBK VBM VBR VBS VBW VBY
VDK VDM VDR VDS VDW VDY VHK VHM VHR VHS VHW VHY VVK VVM VVR VVS VVW VVY VKB
VKD VKH VKV VMB VMD VMH VMV VRB VRD VRH VRV VSB VSD VSH VSV VWB VWD VWH VWV
VYB VYD VYH VYV KBB KBD KBH KBV KDB KDD KDH KDV KHB KHD KHH KHV KVB KVD KVH
KVV MBB MBD MBH MBV MDB MDD MDH MDV MHB MHD MHH MHV MVB MVD MVH MVV RBB RBD
RBH RBV RDB RDD RDH RDV RHB RHD RHH RHV RVB RVD RVH RVV SBB SBD SBH SBV SDB
SDD SDH SDV SHB SHD SHH SHV SVB SVD SVH SVV WBB WBD WBH WBV WDB WDD WDH WDV
WHB WHD WHH WHV WVB WVD WVH WVV YBB YBD YBH YBV YDB YDD YDH YDV YHB YHD YHH
YHV YVB YVD YVH YVV NNA NNC NNT NNG NKK NKM NKR NKS NKW NKY NMK NMM NMR NMS
NMW NMY NRK NRM NRR NRS NRW NRY NSK NSM NSR NSS NSW NSY NWK NWM NWR NWS NWW
NWY NYK NYM NYR NYS NYW NYY NAN NCN NTN NGN KNK KNM KNR KNS KNW KNY KKN KMN
KRN KSN KWN KYN MNK MNM MNR MNS MNW MNY MKN MMN MRN MSN MWN MYN RNK RNM RNR
RNS RNW RNY RKN RMN RRN RSN RWN RYN SNK SNM SNR SNS SNW SNY SKN SMN SRN SSN
SWN SYN WNK WNM WNR WNS WNW WNY WKN WMN WRN WSN WWN WYN YNK YNM YNR YNS YNW
YNY YKN YMN YRN YSN YWN YYN ANN CNN TNN GNN NBA NBC NBT NBG NDA NDC NDT NDG
NHA NHC NHT NHG NVA NVC NVT NVG NAB NAD NAH NAV NCB NCD NCH NCV NTB NTD NTH
NTV NGB NGD NGH NGV BNA BNC BNT BNG BKK BKM BKR BKS BKW BKY BMK BMM BMR BMS
BMW BMY BRK BRM BRR BRS BRW BRY BSK BSM BSR BSS BSW BSY BWK BWM BWR BWS BWW
BWY BYK BYM BYR BYS BYW BYY BAN BCN BTN BGN DNA DNC DNT DNG DKK DKM DKR DKS
DKW DKY DMK DMM DMR DMS DMW DMY DRK DRM DRR DRS DRW DRY DSK DSM DSR DSS DSW
DSY DWK DWM DWR DWS DWW DWY DYK DYM DYR DYS DYW DYY DAN DCN DTN DGN HNA HNC
HNT HNG HKK HKM HKR HKS HKW HKY HMK HMM HMR HMS HMW HMY HRK HRM HRR HRS HRW
HRY HSK HSM HSR HSS HSW HSY HWK HWM HWR HWS HWW HWY HYK HYM HYR HYS HYW HYY
HAN HCN HTN HGN VNA VNC VNT VNG VKK VKM VKR VKS VKW VKY VMK VMM VMR VMS VMW
VMY VRK VRM VRR VRS VRW VRY VSK VSM VSR VSS VSW VSY VWK VWM VWR VWS VWW VWY
VYK VYM VYR VYS VYW VYY VAN VCN VTN VGN KBK KBM KBR KBS KBW KBY KDK KDM KDR
KDS KDW KDY KHK KHM KHR KHS KHW KHY KVK KVM KVR KVS KVW KVY KKB KKD KKH KKV
KMB KMD KMH KMV KRB KRD KRH KRV KSB KSD KSH KSV KWB KWD KWH KWV KYB KYD KYH
KYV MBK MBM MBR MBS MBW MBY MDK MDM MDR MDS MDW MDY MHK MHM MHR MHS MHW MHY
MVK MVM MVR MVS MVW MVY MKB MKD MKH MKV MMB MMD MMH MMV MRB MRD MRH MRV MSB
MSD MSH MSV MWB MWD MWH MWV MYB MYD MYH MYV RBK RBM RBR RBS RBW RBY RDK RDM
RDR RDS RDW RDY RHK RHM RHR RHS RHW RHY RVK RVM RVR RVS RVW RVY RKB RKD RKH
RKV RMB RMD RMH RMV RRB RRD RRH RRV RSB RSD RSH RSV RWB RWD RWH RWV RYB RYD
RYH RYV SBK SBM SBR SBS SBW SBY SDK SDM SDR SDS SDW SDY SHK SHM SHR SHS SHW
SHY SVK SVM SVR SVS SVW SVY SKB SKD SKH SKV SMB SMD SMH SMV SRB SRD SRH SRV
SSB SSD SSH SSV SWB SWD SWH SWV SYB SYD SYH SYV WBK WBM WBR WBS WBW WBY WDK
WDM WDR WDS WDW WDY WHK WHM WHR WHS WHW WHY WVK WVM WVR WVS WVW WVY WKB WKD
WKH WKV WMB WMD WMH WMV WRB WRD WRH WRV WSB WSD WSH WSV WWB WWD WWH WWV WYB
WYD WYH WYV YBK YBM YBR YBS YBW YBY YDK YDM YDR YDS YDW YDY YHK YHM YHR YHS
YHW YHY YVK YVM YVR YVS YVW YVY YKB YKD YKH YKV YMB YMD YMH YMV YRB YRD YRH
YRV YSB YSD YSH YSV YWB YWD YWH YWV YYB YYD YYH YYV ANB AND ANH ANV ABN ADN
AHN AVN CNB CND CNH CNV CBN CDN CHN CVN TNB TND TNH TNV TBN TDN THN TVN GNB
GND GNH GNV GBN GDN GHN GVN BBA BBC BBT BBG BDA BDC BDT BDG BHA BHC BHT BHG
BVA BVC BVT BVG BAB BAD BAH BAV BCB BCD BCH BCV BTB BTD BTH BTV BGB BGD BGH
BGV DBA DBC DBT DBG DDA DDC DDT DDG DHA DHC DHT DHG DVA DVC DVT DVG DAB DAD
DAH DAV DCB DCD DCH DCV DTB DTD DTH DTV DGB DGD DGH DGV HBA HBC HBT HBG HDA
HDC HDT HDG HHA HHC HHT HHG HVA HVC HVT HVG HAB HAD HAH HAV HCB HCD HCH HCV
HTB HTD HTH HTV HGB HGD HGH HGV VBA VBC VBT VBG VDA VDC VDT VDG VHA VHC VHT
VHG VVA VVC VVT VVG VAB VAD VAH VAV VCB VCD VCH VCV VTB VTD VTH VTV VGB VGD
VGH VGV ABB ABD ABH ABV ADB ADD ADH ADV AHB AHD AHH AHV AVB AVD AVH AVV CBB
CBD CBH CBV CDB CDD CDH CDV CHB CHD CHH CHV CVB CVD CVH CVV TBB TBD TBH TBV
TDB TDD TDH TDV THB THD THH THV TVB TVD TVH TVV GBB GBD GBH GBV GDB GDD GDH
GDV GHB GHD GHH GHV GVB GVD GVH GVV NKA NKC NKT NKG NMA NMC NMT NMG NRA NRC
NRT NRG NSA NSC NST NSG NWA NWC NWT NWG NYA NYC NYT NYG NAK NAM NAR NAS NAW
NAY NCK NCM NCR NCS NCW NCY NTK NTM NTR NTS NTW NTY NGK NGM NGR NGS NGW NGY
KNA KNC KNT KNG KKK KKM KKR KKS KKW KKY KMK KMM KMR KMS KMW KMY KRK KRM KRR
KRS KRW KRY KSK KSM KSR KSS KSW KSY KWK KWM KWR KWS KWW KWY KYK KYM KYR KYS
KYW KYY KAN KCN KTN KGN MNA MNC MNT MNG MKK MKM MKR MKS MKW MKY MMK MMM MMR
MMS MMW MMY MRK MRM MRR MRS MRW MRY MSK MSM MSR MSS MSW MSY MWK MWM MWR MWS
MWW MWY MYK MYM MYR MYS MYW MYY MAN MCN MTN MGN RNA RNC RNT RNG RKK RKM RKR
RKS RKW RKY RMK RMM RMR RMS RMW RMY RRK RRM RRR RRS RRW RRY RSK RSM RSR RSS
RSW RSY RWK RWM RWR RWS RWW RWY RYK RYM RYR RYS RYW RYY RAN RCN RTN RGN SNA
SNC SNT SNG SKK SKM SKR SKS SKW SKY SMK SMM SMR SMS SMW SMY SRK SRM SRR SRS
SRW SRY SSK SSM SSR SSS SSW SSY SWK SWM SWR SWS SWW SWY SYK SYM SYR SYS SYW
SYY SAN SCN STN SGN WNA WNC WNT WNG WKK WKM WKR WKS WKW WKY WMK WMM WMR WMS
WMW WMY WRK WRM WRR WRS WRW WRY WSK WSM WSR WSS WSW WSY WWK WWM WWR WWS WWW
WWY WYK WYM WYR WYS WYW WYY WAN WCN WTN WGN YNA YNC YNT YNG YKK YKM YKR YKS
YKW YKY YMK YMM YMR YMS YMW YMY YRK YRM YRR YRS YRW YRY YSK YSM YSR YSS YSW
YSY YWK YWM YWR YWS YWW YWY YYK YYM YYR YYS YYW YYY YAN YCN YTN YGN ANK ANM
ANR ANS ANW ANY AKN AMN ARN ASN AWN AYN CNK CNM CNR CNS CNW CNY CKN CMN CRN
CSN CWN CYN TNK TNM TNR TNS TNW TNY TKN TMN TRN TSN TWN TYN GNK GNM GNR GNS
GNW GNY GKN GMN GRN GSN GWN GYN BKA BKC BKT BKG BMA BMC BMT BMG BRA BRC BRT
BRG BSA BSC BST BSG BWA BWC BWT BWG BYA BYC BYT BYG BAK BAM BAR BAS BAW BAY
BCK BCM BCR BCS BCW BCY BTK BTM BTR BTS BTW BTY BGK BGM BGR BGS BGW BGY DKA
DKC DKT DKG DMA DMC DMT DMG DRA DRC DRT DRG DSA DSC DST DSG DWA DWC DWT DWG
DYA DYC DYT DYG DAK DAM DAR DAS DAW DAY DCK DCM DCR DCS DCW DCY DTK DTM DTR
DTS DTW DTY DGK DGM DGR DGS DGW DGY HKA HKC HKT HKG HMA HMC HMT HMG HRA HRC
HRT HRG HSA HSC HST HSG HWA HWC HWT HWG HYA HYC HYT HYG HAK HAM HAR HAS HAW
HAY HCK HCM HCR HCS HCW HCY HTK HTM HTR HTS HTW HTY HGK HGM HGR HGS HGW HGY
VKA VKC VKT VKG VMA VMC VMT VMG VRA VRC VRT VRG VSA VSC VST VSG VWA VWC VWT
VWG VYA VYC VYT VYG VAK VAM VAR VAS VAW VAY VCK VCM VCR VCS VCW VCY VTK VTM
VTR VTS VTW VTY VGK VGM VGR VGS VGW VGY KBA KBC KBT KBG KDA KDC KDT KDG KHA
KHC KHT KHG KVA KVC KVT KVG KAB KAD KAH KAV KCB KCD KCH KCV KTB KTD KTH KTV
KGB KGD KGH KGV MBA MBC MBT MBG MDA MDC MDT MDG MHA MHC MHT MHG MVA MVC MVT
MVG MAB MAD MAH MAV MCB MCD MCH MCV MTB MTD MTH MTV MGB MGD MGH MGV RBA RBC
RBT RBG RDA RDC RDT RDG RHA RHC RHT RHG RVA RVC RVT RVG RAB RAD RAH RAV RCB
RCD RCH RCV RTB RTD RTH RTV RGB RGD RGH RGV SBA SBC SBT SBG SDA SDC SDT SDG
SHA SHC SHT SHG SVA SVC SVT SVG SAB SAD SAH SAV SCB SCD SCH SCV STB STD STH
STV SGB SGD SGH SGV WBA WBC WBT WBG WDA WDC WDT WDG WHA WHC WHT WHG WVA WVC
WVT WVG WAB WAD WAH WAV WCB WCD WCH WCV WTB WTD WTH WTV WGB WGD WGH WGV YBA
YBC YBT YBG YDA YDC YDT YDG YHA YHC YHT YHG YVA YVC YVT YVG YAB YAD YAH YAV
YCB YCD YCH YCV YTB YTD YTH YTV YGB YGD YGH YGV ABK ABM ABR ABS ABW ABY ADK
ADM ADR ADS ADW ADY AHK AHM AHR AHS AHW AHY AVK AVM AVR AVS AVW AVY AKB AKD
AKH AKV AMB AMD AMH AMV ARB ARD ARH ARV ASB ASD ASH ASV AWB AWD AWH AWV AYB
AYD AYH AYV CBK CBM CBR CBS CBW CBY CDK CDM CDR CDS CDW CDY CHK CHM CHR CHS
CHW CHY CVK CVM CVR CVS CVW CVY CKB CKD CKH CKV CMB CMD CMH CMV CRB CRD CRH
CRV CSB CSD CSH CSV CWB CWD CWH CWV CYB CYD CYH CYV TBK TBM TBR TBS TBW TBY
TDK TDM TDR TDS TDW TDY THK THM THR THS THW THY TVK TVM TVR TVS TVW TVY TKB
TKD TKH TKV TMB TMD TMH TMV TRB TRD TRH TRV TSB TSD TSH TSV TWB TWD TWH TWV
TYB TYD TYH TYV GBK GBM GBR GBS GBW GBY GDK GDM GDR GDS GDW GDY GHK GHM GHR
GHS GHW GHY GVK GVM GVR GVS GVW GVY GKB GKD GKH GKV GMB GMD GMH GMV GRB GRD
GRH GRV GSB GSD GSH GSV GWB GWD GWH GWV GYB GYD GYH GYV NAA NAC NAT NAG NCA
NCC NCT NCG NTA NTC NTT NTG NGA NGC NGT NGG KKA KKC KKT KKG KMA KMC KMT KMG
KRA KRC KRT KRG KSA KSC KST KSG KWA KWC KWT KWG KYA KYC KYT KYG KAK KAM KAR
KAS KAW KAY KCK KCM KCR KCS KCW KCY KTK KTM KTR KTS KTW KTY KGK KGM KGR KGS
KGW KGY MKA MKC MKT MKG MMA MMC MMT MMG MRA MRC MRT MRG MSA MSC MST MSG MWA
MWC MWT MWG MYA MYC MYT MYG MAK MAM MAR MAS MAW MAY MCK MCM MCR MCS MCW MCY
MTK MTM MTR MTS MTW MTY MGK MGM MGR MGS MGW MGY RKA RKC RKT RKG RMA RMC RMT
RMG RRA RRC RRT RRG RSA RSC RST RSG RWA RWC RWT RWG RYA RYC RYT RYG RAK RAM
RAR RAS RAW RAY RCK RCM RCR RCS RCW RCY RTK RTM RTR RTS RTW RTY RGK RGM RGR
RGS RGW RGY SKA SKC SKT SKG SMA SMC SMT SMG SRA SRC SRT SRG SSA SSC SST SSG
SWA SWC SWT SWG SYA SYC SYT SYG SAK SAM SAR SAS SAW SAY SCK SCM SCR SCS SCW
SCY STK STM STR STS STW STY SGK SGM SGR SGS SGW SGY WKA WKC WKT WKG WMA WMC
WMT WMG WRA WRC WRT WRG WSA WSC WST WSG WWA WWC WWT WWG WYA WYC WYT WYG WAK
WAM WAR WAS WAW WAY WCK WCM WCR WCS WCW WCY WTK WTM WTR WTS WTW WTY WGK WGM
WGR WGS WGW WGY YKA YKC YKT YKG YMA YMC YMT YMG YRA YRC YRT YRG YSA YSC YST
YSG YWA YWC YWT YWG YYA YYC YYT YYG YAK YAM YAR YAS YAW YAY YCK YCM YCR YCS
YCW YCY YTK YTM YTR YTS YTW YTY YGK YGM YGR YGS YGW YGY ANA ANC ANT ANG AKK
AKM AKR AKS AKW AKY AMK AMM AMR AMS AMW AMY ARK ARM ARR ARS ARW ARY ASK ASM
ASR ASS ASW ASY AWK AWM AWR AWS AWW AWY AYK AYM AYR AYS AYW AYY AAN ACN ATN
AGN CNA CNC CNT CNG CKK CKM CKR CKS CKW CKY CMK CMM CMR CMS CMW CMY CRK CRM
CRR CRS CRW CRY CSK CSM CSR CSS CSW CSY CWK CWM CWR CWS CWW CWY CYK CYM CYR
CYS CYW CYY CAN CCN CTN CGN TNA TNC TNT TNG TKK TKM TKR TKS TKW TKY TMK TMM
TMR TMS TMW TMY TRK TRM TRR TRS TRW TRY TSK TSM TSR TSS TSW TSY TWK TWM TWR
TWS TWW TWY TYK TYM TYR TYS TYW TYY TAN TCN TTN TGN GNA GNC GNT GNG GKK GKM
GKR GKS GKW GKY GMK GMM GMR GMS GMW GMY GRK GRM GRR GRS GRW GRY GSK GSM GSR
GSS GSW GSY GWK GWM GWR GWS GWW GWY GYK GYM GYR GYS GYW GYY GAN GCN GTN GGN
BAA BAC BAT BAG BCA BCC BCT BCG BTA BTC BTT BTG BGA BGC BGT BGG DAA DAC DAT
DAG DCA DCC DCT DCG DTA DTC DTT DTG DGA DGC DGT DGG HAA HAC HAT HAG HCA HCC
HCT HCG HTA HTC HTT HTG HGA HGC HGT HGG VAA VAC VAT VAG VCA VCC VCT VCG VTA
VTC VTT VTG VGA VGC VGT VGG ABA ABC ABT ABG ADA ADC ADT ADG AHA AHC AHT AHG
AVA AVC AVT AVG AAB AAD AAH AAV ACB ACD ACH ACV ATB ATD ATH ATV AGB AGD AGH
AGV CBA CBC CBT CBG CDA CDC CDT CDG CHA CHC CHT CHG CVA CVC CVT CVG CAB CAD
CAH CAV CCB CCD CCH CCV CTB CTD CTH CTV CGB CGD CGH CGV TBA TBC TBT TBG TDA
TDC TDT TDG THA THC THT THG TVA TVC TVT TVG TAB TAD TAH TAV TCB TCD TCH TCV
TTB TTD TTH TTV TGB TGD TGH TGV GBA GBC GBT GBG GDA GDC GDT GDG GHA GHC GHT
GHG GVA GVC GVT GVG GAB GAD GAH GAV GCB GCD GCH GCV GTB GTD GTH GTV GGB GGD
GGH GGV KAA KAC KAT KAG KCA KCC KCT KCG KTA KTC KTT KTG KGA KGC KGT KGG MAA
MAC MAT MAG MCA MCC MCT MCG MTA MTC MTT MTG MGA MGC MGT MGG RAA RAC RAT RAG
RCA RCC RCT RCG RTA RTC RTT RTG RGA RGC RGT RGG SAA SAC SAT SAG SCA SCC SCT
SCG STA STC STT STG SGA SGC SGT SGG WAA WAC WAT WAG WCA WCC WCT WCG WTA WTC
WTT WTG WGA WGC WGT WGG YAA YAC YAT YAG YCA YCC YCT YCG YTA YTC YTT YTG YGA
YGC YGT YGG AKA AKC AKT AKG AMA AMC AMT AMG ARA ARC ART ARG ASA ASC AST ASG
AWA AWC AWT AWG AYA AYC AYT AYG AAK AAM AAR AAS AAW AAY ACK ACM ACR ACS ACW
ACY ATK ATM ATR ATS ATW ATY AGK AGM AGR AGS AGW AGY CKA CKC CKT CKG CMA CMC
CMT CMG CRA CRC CRT CRG CSA CSC CST CSG CWA CWC CWT CWG CYA CYC CYT CYG CAK
CAM CAR CAS CAW CAY CCK CCM CCR CCS CCW CCY CTK CTM CTR CTS CTW CTY CGK CGM
CGR CGS CGW CGY TKA TKC TKT TKG TMA TMC TMT TMG TRA TRC TRT TRG TSA TSC TST
TSG TWA TWC TWT TWG TYA TYC TYT TYG TAK TAM TAR TAS TAW TAY TCK TCM TCR TCS
TCW TCY TTK TTM TTR TTS TTW TTY TGK TGM TGR TGS TGW TGY GKA GKC GKT GKG GMA
GMC GMT GMG GRA GRC GRT GRG GSA GSC GST GSG GWA GWC GWT GWG GYA GYC GYT GYG
GAK GAM GAR GAS GAW GAY GCK GCM GCR GCS GCW GCY GTK GTM GTR GTS GTW GTY GGK
GGM GGR GGS GGW GGY AAA AAC AAT AAG ACA ACC ACT ACG ATA ATC ATT ATG AGA AGC
AGT AGG CAA CAC CAT CAG CCA CCC CCT CCG CTA CTC CTT CTG CGA CGC CGT CGG TAA
TAC TAT TAG TCA TCC TCT TCG TTA TTC TTT TTG TGA TGC TGT TGG GAA GAC GAT GAG
GCA GCC GCT GCG GTA GTC GTT GTG GGA GGC GGT GGG
