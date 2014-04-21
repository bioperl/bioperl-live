#
# BioPerl module for Bio::Restriction::IO::base
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Rob Edwards <redwards@utmem.edu>
#
# Copyright Rob Edwards
#
# You may distribute this module under the same terms as perl itself
#

# POD documentation - main docs before the code

=head1 NAME

Bio::Restriction::IO::base - base enzyme set

=head1 SYNOPSIS

Do not use this module directly.  Use it via the Bio::Restriction::IO class.

=head1 DESCRIPTION


This class defines some base methods for restriction enzyme input and
at the same time gives a base list of common enzymes.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Rob Edwards, redwards@utmem.edu

=head1 CONTRIBUTORS

Heikki Lehvaslaiho, heikki-at-bioperl-dot-org
Mark A. Jensen, maj-at-fortinbras-dot-us

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Restriction::IO::base;

use strict;

use Bio::Restriction::Enzyme;
use Bio::Restriction::EnzymeCollection;
use Bio::Restriction::Enzyme::MultiCut;
use Bio::Restriction::Enzyme::MultiSite;

use base qw(Bio::Restriction::IO);

my $offset; # class variable

sub new {
    my($class, @args) = @_;
    $class = ref $class ? ref $class : $class;
    my $self = bless {}, $class;
    $self->_initialize(@args);
    return $self;
}


{
    
my %FILE_FORMAT = (
            #'itype2'    => 'itype2', # itype2 format doesn't work with 'current'
            #'8'         => 'itype2',
            'withrefm'  => 'withrefm',
            '31'        => 'withrefm',
            #'bairoch'   => 'bairoch', # bairoch format doesn't work with 'current'
            #'19'        => 'bairoch',
            #'macvector' => 'bairoch',
            #'vectorNTI' => 'bairoch',
            'neo'       => 'neos',
            'prototype' => 'proto'
);

sub _initialize {
    my($self,@args) = @_;
    my ($current, $url, $file, $fh, $format, $verbose) =
        $self->_rearrange([qw(CURRENT URL FILE FH FORMAT VERBOSE)],@args);
    $verbose || 0;
    $self->verbose($verbose);
    if ($current && $format) {
        $self->throw("Can't use -current with file, fh, or url set")  if ($url || $file || $fh);
        $self->throw("Format $format not retrievable using 'current'") if (!exists $FILE_FORMAT{$format});
        my $io = $self->new(-url => 'ftp://ftp.neb.com/pub/rebase/VERSION');
        chomp (my $version = $io->_readline);
        push @args, (-url => "ftp://ftp.neb.com/pub/rebase/$FILE_FORMAT{$format}.$version", -retries => 1);
    }

    $self->_companies;
    return unless $self->SUPER::_initialize(@args);
}

}

=head2 read

 Title   : read
 Usage   : $renzs = $stream->read
 Function: reads all the restrction enzymes from the stream
 Returns : a Bio::Restriction::Restriction object
 Args    : none

=cut

sub read {
    my $self = shift;

    my $renzs = Bio::Restriction::EnzymeCollection->new(-empty => 1);
    seek DATA,($offset||=tell DATA), 0;
    while (<DATA>) {
        chomp;
        next if /^\s*$/;
        my ($name, $site, $cut) = split /\s+/;
        my $re = Bio::Restriction::Enzyme->new(-name => $name,
                                              -site => $site,
                                              -cut => $cut);
        $renzs->enzymes($re);
    }
    return $renzs;
}

=head2 _xln_sub

 Title   : _xln_sub
 Function: Translates withrefm coords to Bio::Restriction coords
 Args    : Bio::Restriction::Enzyme object, scalar integer (cut posn)
 Note    : Used internally; pass as a coderef to the B:R::Enzyme 
           constructor
 Note    : It is convenient for each format module to have its own 
           version of this; not currently demanded by the interface.
=cut

sub _xln_sub { # for base.pm, a no-op
    my ($z,$c) = @_; 
    return $c;
}


=head2 write

 Title   : write
 Usage   : $stream->write($renzs)
 Function: writes restriction enzymes into the stream
 Returns : 1 for success and 0 for error
 Args    : a Bio::Restriction::Enzyme
           or a Bio::Restriction::EnzymeCollection object

=cut

sub write {
    my $self = shift;
    foreach (@_) {
        map { printf "%s\t%s\t%s\n", $_->name, $_->string, $_->cut
          } sort {$a->name cmp $b->name} $_->each_enzyme
            if $_->isa('Bio::Restriction::EnzymeCollection');
        printf "%s\t%s\t%s\n", $_->name, $_->string, $_->cut 
            if $_->isa('Bio::Restriction::Enzyme');
    }
}

=head2 verify_prototype

 Title     : verify_prototype
 Purpose   : checks enzyme against current prototype list (retrieved remotely)
 Returns   : returns TRUE if enzyme is prototype
 Argument  : Bio::Restriction::EnzymeI
 Comments  : This is an auxiliary method to retrieve and check an enzyme
             as a prototype.  It retrieves the current list, stores it
             as a singleton instance, then uses it to check the prototype
             and modify is_prototype() to true or false.  Use as follows:

             my $col = $io->read;
             for my $enz ($col->each_enzyme) {
                 print $enz->name.":".$enz->site."\n";
                 print "\t".$io->verify_prototype($enz)."\n";
             }

=cut

my $protodb;

sub verify_prototype {
    my ($self, $enz) = @_;
    $self->throw("Must pass a Bio::Restriction::EnzymeI") unless
        $enz && ref $enz && $enz->isa("Bio::Restriction::EnzymeI");
    if (!defined $protodb) {
        my $io = Bio::Restriction::IO->new(-format => 'prototype',
                      -current => 1);
        $protodb = $io->read;
    }
    if ($protodb->get_enzyme($enz->name)) {
        $enz->is_prototype(1);
    } else {
        $enz->is_prototype(0);
    }
    $enz->is_prototype;
}

=head2 Common REBASE parsing methods

The rest of the methods in this file are to be used by other REBASE
parsers. They are not to be used outside subclasses of this base
class. (They are 'protected' in the sense the word is used in Java.)

=cut

=head2 _cuts_from_site

 Title   : _cuts_from_site
 Usage   : ($site, $cut, $comp_cut) = _cuts_from_site('ACGCGT(4/5)');
 Function: Separates cut positions from a single site string.
           Does nothing to site if it does not have the cut string
 Returns : array of site_string, forward_cut_position, reverse_cut_position
 Args    : recognition site string
 Note    : Not used in withrefm refactor/maj

=cut

sub _cuts_from_site {
    my ($self, $site) = @_;
    my ($cut, $comp_cut) = $site =~ /\((-?\d+)\/(-?\d+)\)/;
    $site =~ s/\(.*\)$//;
    return ($site, $cut, $comp_cut);
}


=head2 _meth

 Title   : _meth
 Usage   : ($pos, $meth) = $self->_meth('2(5)');
 Function: Separates methylation postion and coce from a string.
           Adjusts the postion depending on enzyme site length
           and symmetry 
 Returns : array of position and methylation code
 Args    : 1. reference to Enzyme object
           2. methylation description string

=cut

sub _meth {
    my ($self, $re, $meth) = @_;

    $meth =~ /(\S+)\((\d+)\)/;
    my ($pos, $m) = ($1, $2);
    $pos = 0 if $pos eq '?';
    $pos = $re->seq->length + $pos if $pos and ! $re->palindromic;
    return ($pos, $m);

    $self->warn("Unknown methylation format [$meth]") if $self->verbose >0;
}


=head2 _coordinate_shift_to_cut

 Title   : _coordinate_shift_to_cut
 Usage   : $cut = $self->_coordinate_shift_to_cut($oricut, offset);
 Function: Adjust cut position coordinates to start from the 
           first nucleotides of site
 Returns : Cut position in correct coordinates
 Args    : 1. Original cut position
           2. Length of the recognition site
 Note    : Not used in withrefm.pm refactor/maj

=cut

sub _coordinate_shift_to_cut {
    my ($self, $cut, $site_length) = @_;
    return $cut + $site_length;
}


=head2 _make_multisites

 Title   : _make_multisites
 Usage   : $self->_make_multisites($first_enzyme, \@sites, \@mets)
 Function: Bless a Bio::Restriction::Enzyme  into
           Bio::Restriction::Enzyme::MultiSite and clone it as many
           times as there are alternative sites.
 Returns : nothing, does in place editing
 Args    : 1. a Bio::Restriction::Enzyme
           2. reference to an array of recognition site strings
           3. reference to an array of methylation code strings, optional

=cut

# removed the enzyme collection from arg list /maj

sub _make_multisites {
    my ($self, $re, $sites, $meths, $xln_sub) = @_;

    bless $re, 'Bio::Restriction::Enzyme::MultiSite';

    my $count = 0;
    while ($count < scalar @{$sites}) {
	# this should probably be refactored to use the constructor
	# too, rather than the clone/accessor method /maj
#        my $re2 = $re->clone;
#	my $re2;

        my $site = @{$sites}[$count];
	my ($precut, $recog, $postcut) = ( $site =~ m/^(?:\((\w+\/\w+)\))?([\w^]+)(?:\((\w+\/\w+)\))?/ );
	
	# set the site attribute
#	$re2->site($recog);

	# set the recog attribute (which will make the regexp transformation
	# if necessary:
#	$re2->recog($recog);
#	$recog = $re2->string;
	
# 	no warnings; # avoid 'uninitialized value' warning against $postcut
#         my ($cut, $comp_cut) = ( $postcut =~  /(-?\d+)\/(-?\d+)/ );
# 	use warnings;
	
	# note the following hard codes the coordinate transformation
	# used for rebase/itype2 : this method will break on the 
	# base.pm format. 
#         if ($cut) {
#             $re2->cut($cut + length $recog);
#             $re2->complementary_cut($comp_cut + length $recog);
# 	}
	
	my $re2 = Bio::Restriction::Enzyme::MultiSite->new(
	    -name     => $re->name,
	    -site     => $recog,
	    -recog    => $recog,
	    -precut   => $precut,
	    -postcut  => $postcut,
	    -xln_sub  => $xln_sub
	    );

        if ($meths and @$meths) {
            $re2->purge_methylation_sites;
            $re2->methylation_sites($self->_meth($re2, @{$meths}[$count]));
        }

        $re->others($re2);
        $count++;
    }

    foreach my $enz ($re->others) {
        $enz->others($re, grep {$_ ne $enz} $re->others);
    }

    1;
}

=head2 _make_multicuts

 Title   : _make_multicuts
 Usage   : $self->_make_multicuts($first_enzyme, $precuts)
 Function: 

           Bless a Bio::Restriction::Enzyme into
           Bio::Restriction::Enzyme::MultiCut and clone it. The precut
           string is processed to replase the cut sites in the cloned
           object. Both objects refer to each other through others() method.

 Returns : nothing, does in place editing
 Args    : 1. a Bio::Restriction::Enzyme
           2. precut string, e.g. '12/7'


The examples we have of multiply cutting enzymes cut only four
times. This protected method deals only with a string of two
integers separated with a slash, e.g. '12/7'. The numbers represent the postions
BEFORE the start of the recognition site, i.e. negative positions.

=cut

# removed the enzyme collection from arg list /maj

sub _make_multicuts {
    my ($self, $re, $precut) = @_;

    bless $re, 'Bio::Restriction::Enzyme::MultiCut';
 
    my ($cut, $comp_cut) = $precut =~ /(-?\d+)\/(-?\d+)/;
    
    my $re2 = $re->clone;

    $re2->cut("-$cut");
    $re2->complementary_cut("-$comp_cut");

    $re->others($re2);

    1;
}

=head2 _companies

 Title     : _companies
 Purpose   : Defines the companies that we know about
 Returns   : A hash
 Argument  : Nothing
 Comments  : An internal method to define the companies that we know about
             REBASE uses a code, and this converts the code to the real name
	     (e.g. A = Amersham Pharmacia Biotech)

=cut

sub _companies {
    # this is just so it is easy to set up the codes that REBASE uses
    my $self=shift;
    my %companies=(
                   'A'=>'Amersham Pharmacia Biotech (1/03)',
                   'C'=>'Minotech Biotechnology (6/01)',
                   'E'=>'Stratagene (1/03)',
                   'F'=>'Fermentas AB (1/03)',
                   'G'=>'Qbiogene (1/03)',
                   'H'=>'American Allied Biochemical, Inc. (10/98)',
                   'I'=>'SibEnzyme Ltd. (1/03)',
                   'J'=>'Nippon Gene Co., Ltd. (6/00)',
                   'K'=>'Takara Shuzo Co. Ltd. (1/03)',
                   'M'=>'Roche Applied Science (1/03)',
                   'N'=>'New England Biolabs (1/03)',
                   'O'=>'Toyobo Biochemicals (11/98)',
                   'P'=>'Megabase Research Products (5/99)',
                   'Q'=>'CHIMERx (1/03)',
                   'R'=>'Promega Corporation (1/03)',
                   'S'=>'Sigma Chemical Corporation (1/03)',
                   'U'=>'Bangalore Genei (1/03)',
                   'V'=>'MRC-Holland (1/03)',
                   'X'=>'EURx Ltd. (1/03)');
    $self->{company}=\%companies;
}

1;

__DATA__
AasI	GACNNNNNNGTC	7
AatI	AGGCCT	3
AccII	CGCG	2
AatII	GACGTC	5
AauI	TGTACA	1
Acc113I	AGTACT	3
Acc16I	TGCGCA	3
Acc65I	GGTACC	1
AccB1I	GGYRCC	1
AccB7I	CCANNNNNTGG	7
AccI	GTMKAC	2
AccIII	TCCGGA	1
AciI	CCGC	1
AclI	AACGTT	2
AcsI	RAATTY	1
AcvI	CACGTG	3
AcyI	GRCGYC	2
AdeI	CACNNNGTG	6
AfaI	GTAC	2
AfeI	AGCGCT	3
AflI	GGWCC	1
AflII	CTTAAG	1
AflIII	ACRYGT	1
AgeI	ACCGGT	1
AhaIII	TTTAAA	3
AhdI	GACNNNNNGTC	6
AhlI	ACTAGT	1
AleI	CACNNNNGTG	5
AluI	AGCT	2
Alw21I	GWGCWC	5
Alw44I	GTGCAC	1
AlwNI	CAGNNNCTG	6
Ama87I	CYCGRG	1
AocI	CCTNAGG	2
Aor51HI	AGCGCT	3
ApaBI	GCANNNNNTGC	8
ApaI	GGGCCC	5
ApaLI	GTGCAC	1
ApoI	RAATTY	1
AscI	GGCGCGCC	2
AseI	ATTAAT	2
AsiAI	ACCGGT	1
AsiSI	GCGATCGC	5
AsnI	ATTAAT	2
Asp700I	GAANNNNTTC	5
Asp718I	GGTACC	1
AspEI	GACNNNNNGTC	6
AspHI	GWGCWC	5
AspI	GACNNNGTC	4
AspLEI	GCGC	3
AspS9I	GGNCC	1
AsuC2I	CCSGG	2
AsuI	GGNCC	1
AsuII	TTCGAA	2
AsuNHI	GCTAGC	1
AvaI	CYCGRG	1
AvaII	GGWCC	1
AviII	TGCGCA	3
AvrII	CCTAGG	1
AxyI	CCTNAGG	2
BalI	TGGCCA	3
BamHI	GGATCC	1
BanI	GGYRCC	1
BanII	GRGCYC	5
BanIII	ATCGAT	2
BbeI	GGCGCC	5
BbrPI	CACGTG	3
BbuI	GCATGC	5
Bbv12I	GWGCWC	5
BclI	TGATCA	1
BcnI	CCSGG	2
BcoI	CYCGRG	1
BcuI	ACTAGT	1
BetI	WCCGGW	1
BfaI	CTAG	1
BfmI	CTRYAG	1
BfrBI	ATGCAT	3
BfrI	CTTAAG	1
BfuCI	GATC	0
BglI	GCCNNNNNGGC	7
BglII	AGATCT	1
BlnI	CCTAGG	1
BloHII	CTGCAG	5
BlpI	GCTNAGC	2
Bme1390I	CCNGG	2
Bme1580I	GKGCMC	5
Bme18I	GGWCC	1
BmtI	GCTAGC	5
BmyI	GDGCHC	5
BoxI	GACNNNNGTC	5
Bpu1102I	GCTNAGC	2
Bpu14I	TTCGAA	2
Bsa29I	ATCGAT	2
BsaAI	YACGTR	3
BsaBI	GATNNNNATC	5
BsaHI	GRCGYC	2
BsaJI	CCNNGG	1
BsaOI	CGRYCG	4
BsaWI	WCCGGW	1
Bsc4I	CCNNNNNNNGG	7
BscBI	GGNNCC	3
BscFI	GATC	0
BscI	ATCGAT	2
Bse118I	RCCGGY	1
Bse21I	CCTNAGG	2
Bse8I	GATNNNNATC	5
BseAI	TCCGGA	1
BseBI	CCWGG	2
BseCI	ATCGAT	2
BseDI	CCNNGG	1
BseJI	GATNNNNATC	5
BseLI	CCNNNNNNNGG	7
BsePI	GCGCGC	1
BseSI	GKGCMC	5
BseX3I	CGGCCG	1
Bsh1236I	CGCG	2
Bsh1285I	CGRYCG	4
BshFI	GGCC	2
BshI	GGCC	2
BshNI	GGYRCC	1
BshTI	ACCGGT	1
BsiBI	GATNNNNATC	5
BsiCI	TTCGAA	2
BsiEI	CGRYCG	4
BsiHKAI	GWGCWC	5
BsiHKCI	CYCGRG	1
BsiLI	CCWGG	2
BsiMI	TCCGGA	1
BsiQI	TGATCA	1
BsiSI	CCGG	1
BsiWI	CGTACG	1
BsiXI	ATCGAT	2
BsiYI	CCNNNNNNNGG	7
BsiZI	GGNCC	1
BslI	CCNNNNNNNGG	7
BsoBI	CYCGRG	1
Bsp106I	ATCGAT	2
Bsp119I	TTCGAA	2
Bsp120I	GGGCCC	1
Bsp1286I	GDGCHC	5
Bsp13I	TCCGGA	1
Bsp1407I	TGTACA	1
Bsp143I	GATC	0
Bsp143II	RGCGCY	5
Bsp1720I	GCTNAGC	2
Bsp19I	CCATGG	1
Bsp68I	TCGCGA	3
BspA2I	CCTAGG	1
BspCI	CGATCG	4
BspDI	ATCGAT	2
BspEI	TCCGGA	1
BspHI	TCATGA	1
BspLI	GGNNCC	3
BspLU11I	ACATGT	1
BspMII	TCCGGA	1
BspT104I	TTCGAA	2
BspT107I	GGYRCC	1
BspTI	CTTAAG	1
BspXI	ATCGAT	2
BsrBRI	GATNNNNATC	5
BsrFI	RCCGGY	1
BsrGI	TGTACA	1
BssAI	RCCGGY	1
BssECI	CCNNGG	1
BssHI	CTCGAG	1
BssHII	GCGCGC	1
BssKI	CCNGG	0
BssNAI	GTATAC	3
BssT1I	CCWWGG	1
Bst1107I	GTATAC	3
Bst2UI	CCWGG	2
Bst4CI	ACNGT	3
Bst98I	CTTAAG	1
BstACI	GRCGYC	2
BstAPI	GCANNNNNTGC	7
BstBAI	YACGTR	3
BstBI	TTCGAA	2
BstC8I	GCNNGC	3
BstDEI	CTNAG	1
BstDSI	CCRYGG	1
BstEII	GGTNACC	1
BstENI	CCTNNNNNAGG	5
BstENII	GATC	0
BstFNI	CGCG	2
BstH2I	RGCGCY	5
BstHHI	GCGC	3
BstHPI	GTTAAC	3
BstKTI	GATC	3
BstMAI	CTGCAG	5
BstMCI	CGRYCG	4
BstMWI	GCNNNNNNNGC	7
BstNI	CCWGG	2
BstNSI	RCATGY	5
BstOI	CCWGG	2
BstPAI	GACNNNNGTC	5
BstPI	GGTNACC	1
BstSCI	CCNGG	0
BstSFI	CTRYAG	1
BstSNI	TACGTA	3
BstUI	CGCG	2
BstX2I	RGATCY	1
BstXI	CCANNNNNNTGG	8
BstYI	RGATCY	1
BstZ17I	GTATAC	3
BstZI	CGGCCG	1
Bsu15I	ATCGAT	2
Bsu36I	CCTNAGG	2
BsuRI	GGCC	2
BsuTUI	ATCGAT	2
BtgI	CCRYGG	1
BthCI	GCNGC	4
Cac8I	GCNNGC	3
CaiI	CAGNNNCTG	6
CauII	CCSGG	2
CciNI	GCGGCCGC	2
CelII	GCTNAGC	2
CfoI	GCGC	3
Cfr10I	RCCGGY	1
Cfr13I	GGNCC	1
Cfr42I	CCGCGG	4
Cfr9I	CCCGGG	1
CfrI	YGGCCR	1
ChaI	GATC	4
ClaI	ATCGAT	2
CpoI	CGGWCCG	2
Csp45I	TTCGAA	2
Csp6I	GTAC	1
CspAI	ACCGGT	1
CspI	CGGWCCG	2
CviAII	CATG	1
CviJI	RGCY	2
CviRI	TGCA	2
CviTI	RGCY	2
CvnI	CCTNAGG	2
DdeI	CTNAG	1
DpnI	GATC	2
DpnII	GATC	0
DraI	TTTAAA	3
DraII	RGGNCCY	2
DraIII	CACNNNGTG	6
DrdI	GACNNNNNNGTC	7
DsaI	CCRYGG	1
DseDI	GACNNNNNNGTC	7
EaeI	YGGCCR	1
EagI	CGGCCG	1
Eam1105I	GACNNNNNGTC	6
Ecl136II	GAGCTC	3
EclHKI	GACNNNNNGTC	6
EclXI	CGGCCG	1
Eco105I	TACGTA	3
Eco130I	CCWWGG	1
Eco147I	AGGCCT	3
Eco24I	GRGCYC	5
Eco32I	GATATC	3
Eco47I	GGWCC	1
Eco47III	AGCGCT	3
Eco52I	CGGCCG	1
Eco72I	CACGTG	3
Eco81I	CCTNAGG	2
Eco88I	CYCGRG	1
Eco91I	GGTNACC	1
EcoHI	CCSGG	0
EcoICRI	GAGCTC	3
EcoNI	CCTNNNNNAGG	5
EcoO109I	RGGNCCY	2
EcoO65I	GGTNACC	1
EcoRI	GAATTC	1
EcoRII	CCWGG	0
EcoRV	GATATC	3
EcoT14I	CCWWGG	1
EcoT22I	ATGCAT	5
EcoT38I	GRGCYC	5
EgeI	GGCGCC	3
EheI	GGCGCC	3
ErhI	CCWWGG	1
EsaBC3I	TCGA	2
EspI	GCTNAGC	2
FatI	CATG	0
FauNDI	CATATG	2
FbaI	TGATCA	1
FblI	GTMKAC	2
FmuI	GGNCC	4
Fnu4HI	GCNGC	2
FnuDII	CGCG	2
FriOI	GRGCYC	5
FseI	GGCCGGCC	6
Fsp4HI	GCNGC	2
FspAI	RTGCGCAY	4
FspI	TGCGCA	3
FunI	AGCGCT	3
FunII	GAATTC	1
HaeI	WGGCCW	3
HaeII	RGCGCY	5
HaeIII	GGCC	2
HapII	CCGG	1
HgiAI	GWGCWC	5
HgiCI	GGYRCC	1
HgiJII	GRGCYC	5
HhaI	GCGC	3
Hin1I	GRCGYC	2
Hin6I	GCGC	1
HinP1I	GCGC	1
HincII	GTYRAC	3
HindI	CAC	2
HindII	GTYRAC	3
HindIII	AAGCTT	1
HinfI	GANTC	1
HpaI	GTTAAC	3
HpaII	CCGG	1
Hpy178III	TCNNGA	2
Hpy188I	TCNGA	3
Hpy188III	TCNNGA	2
Hpy8I	GTNNAC	3
Hpy99I	CGWCG	5
HpyCH4I	CATG	3
HpyCH4III	ACNGT	3
HpyCH4IV	ACGT	1
HpyCH4V	TGCA	2
HpyF10VI	GCNNNNNNNGC	8
Hsp92I	GRCGYC	2
Hsp92II	CATG	4
HspAI	GCGC	1
ItaI	GCNGC	2
KasI	GGCGCC	1
Kpn2I	TCCGGA	1
KpnI	GGTACC	5
Ksp22I	TGATCA	1
KspAI	GTTAAC	3
KspI	CCGCGG	4
Kzo9I	GATC	0
LpnI	RGCGCY	3
LspI	TTCGAA	2
MabI	ACCWGGT	1
MaeI	CTAG	1
MaeII	ACGT	1
MaeIII	GTNAC	0
MamI	GATNNNNATC	5
MboI	GATC	0
McrI	CGRYCG	4
MfeI	CAATTG	1
MflI	RGATCY	1
MhlI	GDGCHC	5
MlsI	TGGCCA	3
MluI	ACGCGT	1
MluNI	TGGCCA	3
Mly113I	GGCGCC	2
Mph1103I	ATGCAT	5
MroI	TCCGGA	1
MroNI	GCCGGC	1
MroXI	GAANNNNTTC	5
MscI	TGGCCA	3
MseI	TTAA	1
MslI	CAYNNNNRTG	5
Msp20I	TGGCCA	3
MspA1I	CMGCKG	3
MspCI	CTTAAG	1
MspI	CCGG	1
MspR9I	CCNGG	2
MssI	GTTTAAAC	4
MstI	TGCGCA	3
MunI	CAATTG	1
MvaI	CCWGG	2
MvnI	CGCG	2
MwoI	GCNNNNNNNGC	7
NaeI	GCCGGC	3
NarI	GGCGCC	2
NciI	CCSGG	2
NcoI	CCATGG	1
NdeI	CATATG	2
NdeII	GATC	0
NgoAIV	GCCGGC	1
NgoMIV	GCCGGC	1
NheI	GCTAGC	1
NlaIII	CATG	4
NlaIV	GGNNCC	3
Nli3877I	CYCGRG	5
NmuCI	GTSAC	0
NotI	GCGGCCGC	2
NruGI	GACNNNNNGTC	6
NruI	TCGCGA	3
NsbI	TGCGCA	3
NsiI	ATGCAT	5
NspBII	CMGCKG	3
NspI	RCATGY	5
NspIII	CYCGRG	1
NspV	TTCGAA	2
OliI	CACNNNNGTG	5
PacI	TTAATTAA	5
PaeI	GCATGC	5
PaeR7I	CTCGAG	1
PagI	TCATGA	1
PalI	GGCC	2
PauI	GCGCGC	1
PceI	AGGCCT	3
PciI	ACATGT	1
PdiI	GCCGGC	3
PdmI	GAANNNNTTC	5
Pfl23II	CGTACG	1
PflBI	CCANNNNNTGG	7
PflFI	GACNNNGTC	4
PflMI	CCANNNNNTGG	7
PfoI	TCCNGGA	1
PinAI	ACCGGT	1
Ple19I	CGATCG	4
PmaCI	CACGTG	3
PmeI	GTTTAAAC	4
PmlI	CACGTG	3
Ppu10I	ATGCAT	1
PpuMI	RGGWCCY	2
PpuXI	RGGWCCY	2
PshAI	GACNNNNGTC	5
PshBI	ATTAAT	2
PsiI	TTATAA	3
Psp03I	GGWCC	4
Psp124BI	GAGCTC	5
Psp1406I	AACGTT	2
Psp5II	RGGWCCY	2
Psp6I	CCWGG	0
PspAI	CCCGGG	1
PspEI	GGTNACC	1
PspGI	CCWGG	0
PspLI	CGTACG	1
PspN4I	GGNNCC	3
PspOMI	GGGCCC	1
PspPI	GGNCC	1
PspPPI	RGGWCCY	2
PssI	RGGNCCY	5
PstI	CTGCAG	5
PsuI	RGATCY	1
PsyI	GACNNNGTC	4
PvuI	CGATCG	4
PvuII	CAGCTG	3
RcaI	TCATGA	1
RsaI	GTAC	2
Rsr2I	CGGWCCG	2
RsrII	CGGWCCG	2
SacI	GAGCTC	5
SacII	CCGCGG	4
SalI	GTCGAC	1
SanDI	GGGWCCC	2
SatI	GCNGC	2
Sau3AI	GATC	0
Sau96I	GGNCC	1
SauI	CCTNAGG	2
SbfI	CCTGCAGG	6
ScaI	AGTACT	3
SciI	CTCGAG	3
ScrFI	CCNGG	2
SdaI	CCTGCAGG	6
SduI	GDGCHC	5
SecI	CCNNGG	1
SelI	CGCG	0
SexAI	ACCWGGT	1
SfcI	CTRYAG	1
SfeI	CTRYAG	1
SfiI	GGCCNNNNNGGCC	8
SfoI	GGCGCC	3
Sfr274I	CTCGAG	1
Sfr303I	CCGCGG	4
SfuI	TTCGAA	2
SgfI	GCGATCGC	5
SgrAI	CRCCGGYG	2
SgrBI	CCGCGG	4
SinI	GGWCC	1
SlaI	CTCGAG	1
SmaI	CCCGGG	3
SmiI	ATTTAAAT	4
SmiMI	CAYNNNNRTG	5
SmlI	CTYRAG	1
SnaBI	TACGTA	3
SpaHI	GCATGC	5
SpeI	ACTAGT	1
SphI	GCATGC	5
SplI	CGTACG	1
SrfI	GCCCGGGC	4
Sse232I	CGCCGGCG	2
Sse8387I	CCTGCAGG	6
Sse8647I	AGGWCCT	2
Sse9I	AATT	0
SseBI	AGGCCT	3
SspBI	TGTACA	1
SspI	AATATT	3
SstI	GAGCTC	5
SstII	CCGCGG	4
StuI	AGGCCT	3
StyI	CCWWGG	1
SunI	CGTACG	1
SwaI	ATTTAAAT	4
TaaI	ACNGT	3
TaiI	ACGT	4
TaqI	TCGA	1
TasI	AATT	0
TatI	WGTACW	1
TauI	GCSGC	4
TelI	GACNNNGTC	4
TfiI	GAWTC	1
ThaI	CGCG	2
TliI	CTCGAG	1
Tru1I	TTAA	1
Tru9I	TTAA	1
TscI	ACGT	4
TseI	GCWGC	1
Tsp45I	GTSAC	0
Tsp4CI	ACNGT	3
Tsp509I	AATT	0
TspEI	AATT	0
Tth111I	GACNNNGTC	4
TthHB8I	TCGA	1
UnbI	GGNCC	0
Van91I	CCANNNNNTGG	7
Vha464I	CTTAAG	1
VneI	GTGCAC	1
VpaK11AI	GGWCC	0
VpaK11BI	GGWCC	1
VspI	ATTAAT	2
XagI	CCTNNNNNAGG	5
XapI	RAATTY	1
XbaI	TCTAGA	1
XceI	RCATGY	5
XcmI	CCANNNNNNNNNTGG	8
XhoI	CTCGAG	1
XhoII	RGATCY	1
XmaCI	CCCGGG	1
XmaI	CCCGGG	1
XmaIII	CGGCCG	1
XmaJI	CCTAGG	1
XmiI	GTMKAC	2
XmnI	GAANNNNTTC	5
XspI	CTAG	1
ZhoI	ATCGAT	2
ZraI	GACGTC	3
Zsp2I	ATGCAT	5
