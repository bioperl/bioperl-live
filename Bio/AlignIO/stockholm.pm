#
# BioPerl module for Bio::AlignIO::stockholm
#
#   Based on the Bio::SeqIO::stockholm module
#       by Ewan Birney <birney@ebi.ac.uk>
#       and Lincoln Stein  <lstein@cshl.org>
#
#       and the SimpleAlign.pm module of Ewan Birney
#
# Copyright Peter Schattner, Chris Fields
#
# You may distribute this module under the same terms as perl itself
# _history
# September 5, 2000
# November 6, 2006 - completely refactor read_aln(), add write_aln()
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::stockholm - stockholm sequence input/output stream

=head1 SYNOPSIS

  # Do not use this module directly.  Use it via the L<Bio::AlignIO> class.

  use Bio::AlignIO;
  use strict;

  my $in = Bio::AlignIO->new(-format => 'stockholm',
                             -file   => 't/data/testaln.stockholm');
  while( my $aln = $in->next_aln ) {

  }

=head1 DESCRIPTION

This object can transform L<Bio::Align::AlignI> objects to and from
stockholm flat file databases.  This has been completely refactored
from the original stockholm parser to handle annotation data and now
includes a write_aln() method for (almost) complete stockholm
format output.

Stockholm alignment records normally contain additional sequence-based
and alignment-based annotation

  GF Lines (alignment feature/annotation):
  #=GF <featurename> <Generic per-file annotation, free text>
  Placed above the alignment

  GC Lines (Alignment consensus)
  #=GC <featurename> <Generic per-column annotation, exactly 1
       character per column>
  Placed below the alignment

  GS Lines (Sequence annotations)
  #=GS <seqname> <featurename> <Generic per-sequence annotation, free
       text>

  GR Lines (Sequence meta data)
  #=GR <seqname> <featurename> <Generic per-sequence AND per-column
       mark up, exactly 1 character per column>

Currently, sequence annotations (those designated with GS tags) are
parsed only for accession numbers and descriptions.  It is intended that
full parsing will be added at some point in the near future along with
a builder option for optionally parsing alignment annotation and meta data.

The following methods/tags are currently used for storing and writing
the alignment annotation data.

    Tag        SimpleAlign
                 Method
    ----------------------------------------------------------------------
     AC        accession
     ID        id
     DE        description
    ----------------------------------------------------------------------

    Tag        Bio::Annotation   TagName                    Parameters
               Class
    ----------------------------------------------------------------------
     AU        SimpleValue       record_authors             value
     SE        SimpleValue       seed_source                value
     GA        SimpleValue       gathering_threshold        value
     NC        SimpleValue       noise_cutoff               value
     TC        SimpleValue       trusted_cutoff             value
     TP        SimpleValue       entry_type                 value
     SQ        SimpleValue       num_sequences              value
     PI        SimpleValue       previous_ids               value
     DC        Comment           database_comment           comment
     CC        Comment           alignment_comment          comment
     DR        Target            dblink                     database
                                                            primary_id
                                                            comment
     AM        SimpleValue       build_method               value
     NE        SimpleValue       pfam_family_accession      value
     NL        SimpleValue       sequence_start_stop        value
     SS        SimpleValue       sec_structure_source       value
     BM        SimpleValue       build_model                value
     RN        Reference         reference                  *
     RC        Reference         reference                  comment
     RM        Reference         reference                  pubmed
     RT        Reference         reference                  title
     RA        Reference         reference                  authors
     RL        Reference         reference                  location
    ----------------------------------------------------------------------
  * RN is generated based on the number of Bio::Annotation::Reference objects

=head2 Custom annotation

Some users may want to add custom annotation beyond those mapped above.
Currently there are two methods to do so; however, the methods used for adding
such annotation may change in the future, particularly if alignment Writer
classes are introduced. In particular, do not rely on changing the global
variables @WRITEORDER or %WRITEMAP as these may be made private at some point.

1) Use (and abuse) the 'custom' tag.  The tagname for the object can differ
from the tagname used to store the object in the AnnotationCollection.

    # AnnotationCollection from the SimpleAlign object
    my $coll = $aln->annotation;
    my $factory = Bio::Annotation::AnnotationFactory->new(-type =>
        Bio::Annotation::SimpleValue');
    my $rfann = $factory->create_object(-value => $str,
                                        -tagname => 'mytag');
    $coll->add_Annotation('custom', $rfann);
    $rfann = $factory->create_object(-value => 'foo',
                                    -tagname => 'bar');
    $coll->add_Annotation('custom', $rfann);

OUTPUT:

    # STOCKHOLM 1.0

    #=GF ID myID12345
    #=GF mytag katnayygqelggvnhdyddlakfyfgaglealdffnnkeaaakiinwvaEDTTRGKIQDLV??
    #=GF mytag TPtd~????LDPETQALLV???????????????????????NAIYFKGRWE?????????~??
    #=GF mytag ??HEF?A?EMDTKPY??DFQH?TNen?????GRI??????V???KVAM??MF?????????N??
    #=GF mytag ???DD?VFGYAEL????DE???????L??D??????A??TALELAY??????????????????
    #=GF mytag ?????????????KG??????Sa???TSMLILLP???????????????D??????????????
    #=GF mytag ???????????EGTr?????AGLGKLLQ??QL????????SREef??DLNK??L???AH????R
    #=GF mytag ????????????L????????????????????????????????????????R?????????R
    #=GF mytag ??QQ???????V???????AVRLPKFSFefefdlkeplknlgmhqafdpnsdvfklmdqavlvi
    #=GF mytag gdlqhayafkvd????????????????????????????????????????????????????
    #=GF mytag ????????????????????????????????????????????????????????????????
    #=GF mytag ????????????????????????????????????????????????????????????????
    #=GF mytag ????????????????????????????????????????????????????????????????
    #=GF mytag ?????????????INVDEAG?TEAAAATAAKFVPLSLppkt??????????????????PIEFV
    #=GF mytag ADRPFAFAIR??????E?PAT?G????SILFIGHVEDPTP?msv?
    #=GF bar foo
    ...

2) Modify the global @WRITEORDER and %WRITEMAP.

    # AnnotationCollection from the SimpleAlign object
    my $coll = $aln->annotation;

    # add to WRITEORDER
    my @order = @Bio::AlignIO::stockholm::WRITEORDER;
    push @order, 'my_stuff';
    @Bio::AlignIO::stockholm::WRITEORDER = @order;

    # make sure new tag maps to something
    $Bio::AlignIO::stockholm::WRITEMAP{my_stuff} = 'Hobbit/SimpleValue';

    my $rfann = $factory->create_object(-value => 'Frodo',
                                        -tagname => 'Hobbit');
    $coll->add_Annotation('my_stuff', $rfann);
    $rfann = $factory->create_object(-value => 'Bilbo',
                                     -tagname => 'Hobbit');
    $coll->add_Annotation('my_stuff', $rfann);

OUTPUT:

    # STOCKHOLM 1.0

    #=GF ID myID12345
    #=GF Hobbit Frodo
    #=GF Hobbit Bilbo
    ....

=head1 FEEDBACK

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

=head1 AUTHORS - Chris Fields, Peter Schattner

Email: cjfields-at-uiuc-dot-edu, schattner@alum.mit.edu

=head1 CONTRIBUTORS

Andreas Kahari, ak-at-ebi.ac.uk
Jason Stajich, jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::stockholm;
use strict;

use Bio::Seq::Meta;
use Bio::AlignIO::Handler::GenericAlignHandler;
use Text::Wrap qw(wrap);

use base qw(Bio::AlignIO);

my $STKVERSION = 'STOCKHOLM 1.0';

# This maps the two-letter annotation key to a Annotation/parameter/tagname
# combination.  Some data is stored using get/set methods ('Methods')  The rest
# is mapped to Annotation objects using the parameter for the parsed data
# and the tagname for, well, the Annotation tagname.  A few are treated differently
# based on the type of data stored (Reference data in particular).

my %MAPPING = (
    'AC'    =>  'ACCESSION',
    'ID'    =>  'ID',
    'DE'    =>  ['DESCRIPTION' => 'DESCRIPTION'],
    'AU'    =>  ['RECORD_AUTHORS' => 'RECORD_AUTHORS'],
    'SE'    =>  'SEED_SOURCE',
    'BM'    =>  'BUILD_COMMAND',
    'GA'    =>  'GATHERING_THRESHOLD',
    'NC'    =>  'NOISE_CUTOFF',
    'TC'    =>  'TRUSTED_CUTOFF',
    'TP'    =>  'ENTRY_TYPE',
    'SQ'    =>  'NUM_SEQUENCES',
    'PI'    =>  'PREVIOUS_IDS',
    'DC'    =>  ['DATABASE_COMMENT' => 'DATABASE_COMMENT'],
    'DR'    =>  'DBLINK',
    'RN'    =>  ['REFERENCE' => 'REFERENCE'],
    'RC'    =>  ['REFERENCE' => 'COMMENT'],
    'RM'    =>  ['REFERENCE' => 'PUBMED'],
    'RT'    =>  ['REFERENCE' => 'TITLE'],
    'RA'    =>  ['REFERENCE' => 'AUTHORS'],
    'RL'    =>  ['REFERENCE' => 'JOURNAL'],
    'CC'    =>  ['ALIGNMENT_COMMENT' => 'ALIGNMENT_COMMENT'],
    #Pfam-specific
    'AM'    =>  'BUILD_METHOD',
    'NE'    =>  'PFAM_FAMILY_ACCESSION',
    'NL'    =>  'SEQ_START_STOP',
    # Rfam-specific GF lines
    #'SS'    =>  'SEC_STRUCTURE_SOURCE',
    'SEQUENCE' => 'SEQUENCE'
);

# this is the order that annotations are written
our @WRITEORDER = qw(accession
  id
  description
  previous_ids
  record_authors
  seed_source
  sec_structure_source
  gathering_threshold
  trusted_cutoff
  noise_cutoff
  entry_type
  build_command
  build_method
  pfam_family_accession
  seq_start_stop
  reference
  database_comment
  custom
  dblink
  alignment_comment
  num_sequences
  seq_annotation
  );

# This maps the tagname back to a tagname-annotation value combination.
# Some data is stored using get/set methods ('Methods'), others
# are mapped b/c of more complex annotation types.

our %WRITEMAP = (
    'accession'             =>  'AC/Method',
    'id'                    =>  'ID/Method',
    'description'           =>  'DE/Method',
    'record_authors'        =>  'AU/SimpleValue',
    'seed_source'           =>  'SE/SimpleValue',
    'build_command'         =>  'BM/SimpleValue',
    'gathering_threshold'   =>  'GA/SimpleValue',
    'noise_cutoff'          =>  'NC/SimpleValue',
    'trusted_cutoff'        =>  'TC/SimpleValue',
    'entry_type'            =>  'TP/SimpleValue',
    'num_sequences'         =>  'SQ/SimpleValue',
    'previous_ids'          =>  'PI/SimpleValue',
    'database_comment'      =>  'DC/SimpleValue',
    'dblink'                =>  'DR/DBLink',
    'reference'             =>  'RX/Reference',
    'ref_number'            =>  'RN/number',
    'ref_comment'           =>  'RC/comment',
    'ref_pubmed'            =>  'RM/pubmed',
    'ref_title'             =>  'RT/title',
    'ref_authors'           =>  'RA/authors',
    'ref_location'          =>  'RL/location',
    'alignment_comment'     =>  'CC/Comment',
    'seq_annotation'        =>  'DR/Collection',
    #Pfam-specific
    'build_method'          =>  'AM/SimpleValue',
    'pfam_family_accession' =>  'NE/SimpleValue',
    'seq_start_stop'        =>  'NL/SimpleValue',
    # Rfam-specific GF lines
    'sec_structure_source'  =>  'SS/SimpleValue',
    # custom; this is used to carry over anything from the input alignment
    # not mapped to LocatableSeqs or SimpleAlign in a meaningful way
    'custom'                =>  'XX/SimpleValue'
);

# This maps the tagname back to a tagname-annotation value combination.
# Some data is stored using get/set methods ('Methods'), others
# are mapped b/c of more complex annotation types.

=head2 new

 Title   : new
 Usage   : my $alignio = Bio::AlignIO->new(-format => 'stockholm'
					  -file   => '>file');
 Function: Initialize a new L<Bio::AlignIO::stockholm> reader or writer
 Returns : L<Bio::AlignIO> object
 Args    : -line_length :  length of the line for the alignment block
           -alphabet    :  symbol alphabet to set the sequences to.  If not set,
                           the parser will try to guess based on the alignment
                           accession (if present), defaulting to 'dna'.
           -spaces      :  (optional, def = 1) boolean to add a space in between
                           the "# STOCKHOLM 1.0" header and the annotation and
                           the annotation and the alignment.

=cut

sub _initialize {
    my ( $self, @args ) = @_;
    $self->SUPER::_initialize(@args);
    my ($handler, $linelength, $spaces) = $self->_rearrange([qw(HANDLER LINE_LENGTH SPACES)],@args);
    $spaces = defined $spaces ? $spaces : 1;
    $self->spaces($spaces);
    # hash for functions for decoding keys.
    $handler ? $self->alignhandler($handler) :
    $self->alignhandler(Bio::AlignIO::Handler::GenericAlignHandler->new(
                    -format => 'stockholm',
                    -verbose => $self->verbose,
                    ));
    $linelength && $self->line_length($linelength);
}

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
 Returns : L<Bio::Align::AlignI> object
 Args    : NONE

=cut

sub next_aln {
    my $self = shift;

    my $handler = $self->alignhandler;
    # advance to alignment header
    while( defined(my $line = $self->_readline) ) {
        if ($line =~ m{^\#\s*STOCKHOLM\s+}xmso) {
            last;
        }
    }

    $self->{block_line} = 0;
    # go into main body of alignment
    my ($data_chunk, $isa_primary, $name, $alphabet);
    my $last_feat = '';
    while( defined(my $line = $self->_readline) ) {
        # only blank lines are in between blocks, so reset block line
        my ($primary_tag, $secondary_tag, $data, $nse, $feat, $align, $concat);
        if ($line =~ m{^\s*$}xmso) {
            $self->{block_line} &&= 0;
            next;
        }

        # End of Record
        if (index($line, '//') == 0) {
            # fencepost
            $handler->data_handler($data_chunk);
            undef $data_chunk;
            $handler->data_handler({ALIGNMENT => 1,
                                    NAME => 'ALPHABET',
                                    DATA => $self->alphabet})
                if $self->alphabet;
            last;
        }
        elsif ($line =~ m{^\#=([A-Z]{2})\s+([^\n]+?)\s*$}xmso) {
            ($primary_tag, $data) = ($1, $2);
            if ($primary_tag eq 'GS' || $primary_tag eq 'GR') {
                ($nse, $feat, $data) = split(/\s+/, $data, 3);
            } else {
                ($feat, $data) = split(/\s+/, $data, 2);
            }
            $align = ($primary_tag eq 'GF' || $primary_tag eq 'GR') ? 1 : 0;
        }
        elsif ($line =~ m{^(\S+)\s+([^\s]+)\s*}) {
            $self->{block_line}++;
            ($feat, $nse, $data) = ('SEQUENCE', $1, $2);
        }
        else {
            $self->debug("Missed line : $line\n");
        }
        $primary_tag ||= ''; # when no #= line is present
        $align ||= 0;

        # array refs where the two values are equal indicate the start of a
        # primary chunk of data, otherwise it is to be folded into the last
        # data chunk under a secondary tag.  These are also concatenated
        # to previous values if the

        if (exists($MAPPING{$feat}) && ref $MAPPING{$feat} eq 'ARRAY') {
            ($name, $secondary_tag, $isa_primary) = ( $MAPPING{$feat}->[0] eq $MAPPING{$feat}->[1] ) ?
                ($MAPPING{$feat}->[0], 'DATA', 1) :
                (@{ $MAPPING{$feat} }, 0) ;
            $concat = $last_feat eq $feat ? 1 : 0;
        } elsif (exists($MAPPING{$feat})) {
            ($name, $secondary_tag, $isa_primary) = ($MAPPING{$feat}, 'DATA', 1);
            # catch alphabet here if possible
            if ($align && $name eq 'ACCESSION' && !$self->alphabet) {
                if ($data =~ m{^(P|R)F}) {
                    $self->alphabet($1 eq 'R' ? 'rna' : $1 eq 'P' ? 'protein' : undef );
                }
            }
        } else {
            $name = ($primary_tag eq 'GR') ? 'NAMED_META' :
                    ($primary_tag eq 'GC') ? 'CONSENSUS_META' :
                    'CUSTOM';
            ($secondary_tag, $isa_primary) = ('DATA', 1);
        }

        # Since we can't determine whether data should be passed into the
        # Handler until the next round (due to concatenation and combining
        # data), we always check for the presence of the last chunk when the
        # occasion calls for it (i.e. when the current data string needs to go
        # into a new data chunk). If the data needs to be concatenated it is
        # flagged above and checked below (and passed by if the conditions
        # warrant it).

        # We run into a bit of a fencepost problem, (one chunk left over at
        # the end); that is taken care of above when the end of the record is
        # found.

        if ($isa_primary && defined $data_chunk && !$concat) {
            $handler->data_handler($data_chunk);
            undef $data_chunk;
        }
        $data_chunk->{NAME} = $name;       # used for the handler
        $data_chunk->{ALIGNMENT} = $align; # flag that determines chunk destination
        $data_chunk->{$secondary_tag} .= (defined($data_chunk->{$secondary_tag})) ?
            ' '.$data : $data;
        $data_chunk->{NSE} = $nse if $nse;
        if ($name eq 'SEQUENCE' || $name eq 'NAMED_META' || $name eq 'CONSENSUS_META') {
            $data_chunk->{BLOCK_LINE} = $self->{block_line};
            $data_chunk->{META_TAG} = $feat if ($name ne 'SEQUENCE');
        }
        $last_feat = $feat;
    }

    my $aln = $handler->build_alignment;
    $handler->reset_parameters;
    return $aln;
}

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in stockholm format
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object

=cut

{
    my %LINK_CB = (
        'PDB' => sub {join('; ',($_[0]->database,
                                 $_[0]->primary_id.' '.
                                 ($_[0]->optional_id || ''),
                                 $_[0]->start,
                                 $_[0]->end)).';'},
        'SCOP' => sub {join('; ',($_[0]->database,
                                 $_[0]->primary_id || '',
                                 $_[0]->optional_id)).';'},
        '_DEFAULT_' => sub {join('; ',($_[0]->database,
                                 $_[0]->primary_id)).';'},
    );

sub write_aln {
    # enable array of SimpleAlign objects as well (see clustalw write_aln())
    my ($self, @aln) = @_;
    for my $aln (@aln) {
    $self->throw('Need Bio::Align::AlignI object')
          if (!$aln || !($aln->isa('Bio::Align::AlignI')));

    my $coll = $aln->annotation;
    my ($aln_ann, $seq_ann) =
       ('#=GF ', '#=GS ');
    $self->_print("# $STKVERSION\n") || return 0;
    $self->spaces && $self->_print("\n");
    # annotations first

    #=GF XX ....
    for my $param (@WRITEORDER) {
        my @anns;
        # no point in going through this if there is no annotation!
        last if !$coll;
        # alignment annotations
        my $ct = 1;
        $self->throw("Bad parameter: $param") if !exists $WRITEMAP{$param};
        # get the data, act on it based on the tag
        my ($tag, $key) = split q(/), $WRITEMAP{$param};
        if ($key eq 'Method') {
            push @anns, $aln->$param;
        } else {
            @anns = $coll->get_Annotations($param);
        }
        my $rn = 1;
        ANNOTATIONS:
        for my $ann (@anns) {
            # using Text::Wrap::wrap() for word wrap
            my ($text, $alntag, $data);
            if ($tag eq 'RX') {
                REFS:
                for my $rkey (qw(ref_comment ref_number ref_pubmed
                              ref_title ref_authors ref_location)) {
                    my ($newtag, $method) = split q(/), $WRITEMAP{$rkey};
                    $alntag = sprintf('%-10s',$aln_ann.$newtag);
                    if ($rkey eq 'ref_number') {
                        $data = "[$rn]";
                    } else {
                        $data = $ann->$method;
                    }
                    next REFS unless $data;
                    $text = wrap($alntag, $alntag, $data);
                    $self->_print("$text\n") or return 0;
                }
                $rn++;
                next ANNOTATIONS;
            }
            elsif ($tag eq 'XX') { # custom
                my $newtag = $ann->tagname;
                my $tmp = $aln_ann.$newtag;
                $alntag = sprintf('%-*s',length($tmp) + 1, $tmp);
                $data = $ann->display_text;
            }
            elsif ($tag eq 'SQ') {
                # use the actual number, not the stored Annotation data
                my $tmp = $aln_ann.$tag;
                $alntag = sprintf('%-*s',length($tmp) + 1, $tmp);
                $data = $aln->num_sequences;
            }
            elsif ($tag eq 'DR') {
                my $tmp = $aln_ann.$tag;
                $alntag = sprintf('%-*s',length($tmp) + 1, $tmp);
                my $db = uc $ann->database;
                my $cb = exists $LINK_CB{$db} ? $LINK_CB{$db} : $LINK_CB{_DEFAULT_};
                $data = $ann->display_text($cb);
            }
            else {
                my $tmp = $aln_ann.$tag;
                $alntag = sprintf('%-*s',length($tmp) + 1, $tmp);
                $data = ref $ann ? $ann->display_text : $ann;
            }
            next unless $data;
            $text = wrap($alntag, $alntag, $data);
            $self->_print("$text\n") || return 0;
        }
    }

    #=GS <seq-id> AC xxxxxx
    my $tag = 'AC';
    for my $seq ($aln->each_seq) {
        if (my $acc = $seq->accession_number) {
	    my $text = sprintf("%-4s%-22s%-3s%s\n",$seq_ann,
			       $aln->displayname($seq->get_nse), $tag, $acc);
	    $self->_print($text) || return 0;
        }
    }

    #=GS <seq-id> DR xxxxxx
    $tag = 'DR';
    for my $sf ($aln->get_SeqFeatures) {
        if (my @links = $sf->annotation->get_Annotations('dblink')) {
            for my $link (@links) {
                my $db = uc $link->database;
                my $cb = exists $LINK_CB{$db} ? $LINK_CB{$db} : $LINK_CB{_DEFAULT_};
                my $text = sprintf("%-4s%-22s%-3s%s\n",$seq_ann,
                                   $aln->displayname($sf->entire_seq->get_nse),
                                   $tag,
                                   $link->display_text($cb));
                $self->_print($text) || return 0;
            }
        }
    }

    $self->spaces && $self->_print("\n");
    # now the sequences...

    my $blocklen = $self->line_length;
    my $maxlen = $aln->maxdisplayname_length() + 3;
    my $metalen = $aln->max_metaname_length() || 0;
    if ($blocklen) {
        my $blockstart = 1;
        my $alnlen = $aln->length;
        while ($blockstart < $alnlen) {
            my $subaln = $aln->slice($blockstart, $blockstart+$blocklen-1 ,1);
            $self->_print_seqs($subaln,$maxlen,$metalen);
            $blockstart += $blocklen;
            $self->_print("\n") unless $blockstart >= $alnlen;
        }
    } else {
        $self->_print_seqs($aln,$maxlen,$metalen);
    }

    $self->_print("//\n") || return 0;
    }
    $self->flush() if $self->_flush_on_write && defined $self->_fh;

    return 1;
}

}

=head2 line_length

 Title   : line_length
 Usage   : $obj->line_length($newval)
 Function: Set the alignment output line length
 Returns : value of line_length
 Args    : newvalue (optional)

=cut

sub line_length {
    my ( $self, $value ) = @_;
    if ( defined $value ) {
        $self->{'_line_length'} = $value;
    }
    return $self->{'_line_length'};
}

=head2 spaces

 Title   : spaces
 Usage   : $obj->spaces(1)
 Function: Set the 'spaces' flag, which prints extra newlines between the
           header and the annotation and the annotation and the alignment
 Returns : sequence data type
 Args    : newvalue (optional)

=cut

sub spaces {
    my $self = shift;
    return $self->{'_spaces'} = shift if @_;
    return $self->{'_spaces'};
};

=head2 alignhandler

 Title   : alignhandler
 Usage   : $stream->alignhandler($handler)
 Function: Get/Set the Bio::HandlerBaseI object
 Returns : Bio::HandlerBaseI
 Args    : Bio::HandlerBaseI

=cut

sub alignhandler {
    my ($self, $handler) = @_;
    if ($handler) {
        $self->throw("Not a Bio::HandlerBaseI") unless
        ref($handler) && $handler->isa("Bio::HandlerBaseI");
        $self->{'_alignhandler'} = $handler;
    }
    return $self->{'_alignhandler'};
}

############# PRIVATE INIT/HANDLER METHODS #############

sub _print_seqs {
    my ($self, $aln, $maxlen, $metalen) = @_;

    my ($seq_meta, $aln_meta) = ('#=GR','#=GC');
    # modified (significantly) from AlignIO::pfam

    my ($namestr,$seq,$add);

    # pad extra for meta lines

    for $seq ( $aln->each_seq() ) {
        my ($s, $e, $str) = ($seq->start, $seq->end, $seq->strand);
        $namestr = $aln->displayname($seq->get_nse());
        $self->_print(sprintf("%-*s%s\n",$maxlen+$metalen,
                              $namestr,
                              $seq->seq())) || return 0;
        if ($seq->isa('Bio::Seq::MetaI')) {
            for my $mname ($seq->meta_names) {
                 $self->_print(sprintf("%-*s%s\n",$maxlen+$metalen,
                                       $seq_meta.' '.$namestr.' '.$mname,
                                       $seq->named_meta($mname))) || return 0;
            }
        }
    }
    # alignment consensus
    my $ameta = $aln->consensus_meta;
    if ($ameta) {
        for my $mname ($ameta->meta_names) {
            $self->_print(sprintf("%-*s%s\n",$maxlen+$metalen,
                                  $aln_meta.' '.$mname,
                                  $ameta->named_meta($mname))) || return 0;
        }
    }
}

1;
