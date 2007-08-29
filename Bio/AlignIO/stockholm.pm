# $Id$
#
# BioPerl module for Bio::AlignIO::stockholm

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
     DR        DBLink            aln_dblink                 database
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

=head1 FEEDBACK

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

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
use Bio::Annotation::AnnotationFactory;
use Data::Dumper;
use Text::Wrap qw(wrap);

use base qw(Bio::AlignIO);

our $STKVERSION = 'STOCKHOLM 1.0';

# This maps the two-letter annotation key to a Annotation/parameter/tagname
# combination.  Some data is stored using get/set methods ('Methods')  The rest 
# is mapped to Annotation objects using the parameter for the parsed data
# and the tagname for, well, the Annotation tagname.  A few are treated differently
# based on the type of data stored (Reference data in particular).

our %READMAP = (
            'AC'   => 'Method/accession', 
            'ID'   => 'Method/id', 
            'DE'   => 'Method/description',
            'AU'   => 'SimpleValue/-value/record_authors',
            'SE'   => 'SimpleValue/-value/seed_source', 
            'GA'   => 'SimpleValue/-value/gathering_threshold',
            'NC'   => 'SimpleValue/-value/noise_cutoff', 
            'TC'   => 'SimpleValue/-value/trusted_cutoff', 
            'TP'   => 'SimpleValue/-value/entry_type', 
            'SQ'   => 'SimpleValue/-value/num_sequences', 
            'PI'   => 'SimpleValue/-value/previous_ids', 
            'DC'   => 'Comment/-text/database_comment',
            'CC'   => 'Comment/-text/alignment_comment',
            # DBLink, treated differently
            'DR'   => 'DBLink/-value/aln_dblink',
            # Pfam-specific
            'AM'   => 'SimpleValue/-value/build_method', 
            'NE'   => 'SimpleValue/-value/pfam_family_accession',
            'NL'   => 'SimpleValue/-value/sequence_start_stop',
            # Rfam-specific GF lines
            'SS'   => 'SimpleValue/-value/sec_structure_source',
            # Reference objects mapped differently
            'RN'   => '-number',  # reference number is dumped
            'RC'   => '-comment',
            'RM'   => '-pubmed', 
            'RT'   => '-title', 
            'RA'   => '-authors',
            'RL'   => '-location',
            # Build model mapped differently
            'BM'   => '-value',            
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
  aln_dblink
  alignment_comment
  num_sequences
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
            'aln_dblink'            =>  'DR/DBLink',
            'reference'             =>  'RX/Reference',
            'ref_number'            =>  'RN/number',
            'ref_comment'           =>  'RC/comment',
            'ref_pubmed'            =>  'RM/pubmed',
            'ref_title'             =>  'RT/title',
            'ref_authors'           =>  'RA/authors',
            'ref_location'          =>  'RL/location',
            'alignment_comment'     =>  'CC/Comment',
            #Pfam-specific 
            'build_method'          =>  'AM/SimpleValue',
            'pfam_family_accession' =>  'NE/SimpleValue',
            'seq_start_stop'        =>  'NL/SimpleValue',
            # Rfam-specific GF lines
            'sec_structure_source'  =>  'SS/SimpleValue',
            # custom
            'custom'                =>  'XX/SimpleValue'
            );

sub _initialize {
    my ( $self, @args ) = @_;
    $self->SUPER::_initialize(@args);
    # add arguments to handle build object, interleaved format
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
    my $line;

    my ($id, $name, $seqname, $seq, $count, $tag, $data);
    my $start = my $end = 0;
    my $seen_rc;
    my ($refct, $bct, $lnkct) = (0,0,0);
    my @c2name;
    my (%align, %accession, %desc, %seq_meta, %aln_meta, %annotation);

    # in stockholm format, every non-blank line that does not start
    # with '#=' is an alignment segment; the '#=' lines are mark up lines.
    # Of particular interest are the '#=GF <name/st-ed> AC <accession>'
    # lines, which give accession numbers for each segment

    my $aln =  Bio::SimpleAlign->new(-source => 'stockholm');
    while( defined($line = $self->_readline) ) {
        next if $line =~ /^\s*$/;
        if ($line =~ /^#\s*STOCKHOLM\s+/) {
            last;
        } else {
            $line = $self->_readline;
            # commented this out since some programs have headers with extra output
            # $self->throw("Not Stockholm format: Expecting \"# STOCKHOLM\"; Found \"$line\"");
        }
    }
    
    READLINE:
    while( defined($line = $self->_readline) ) {
        #skip empty lines
        next if $line =~ /^\s+$/;
        
        # Double slash (//) signals end of file.
        last if $line =~ m{^//};
        
        # GF/GS lines, by convention, should be at the top of the alignment
        if ($line =~ m{^\#=GF\s+(\S+?)\s+([^\n]*)$}xms) {

            # alignment annotation
            ($tag, $data) = ($1, $2);
            if (exists $READMAP{$tag}) {

                # reference data (multi line)
                if (index($tag, 'R') == 0) {
                    # comments come before numbering, tricky
                    $refct++ if ( ($tag eq 'RN' && !$seen_rc) || $tag eq 'RC');
                    $seen_rc = 1 if $tag eq 'RC';
                    # Don't need
                    next READLINE if $tag eq 'RN';
                    #                           # of ref       parameter     
                    $annotation{ 'reference' }->[$refct-1]->{ $READMAP{$tag} } .= $data.' ';

                # Build commands (single line)
                } elsif ($tag eq 'BM') {
                    #                            # build cmd    parameter     
                    $annotation{ 'build_command' }->[$bct]->{ $READMAP{$tag} } = $data;
                    $bct++;
                    
                # DBLinks (single line)
                } elsif ($tag eq 'DR') {
                    my ($dbase, $uid, $extra) = split /\s*;\s*/ , $data, 3;
                    my $ref;
                    $ref->{'-database'} = $dbase;
                    $ref->{'-primary_id'} = ($dbase eq 'URL') ? $uid : uc $uid;
                    $ref->{'-comment'} = $extra if $extra;
                    #                       # dblink       parameter list    
                    $annotation{ 'aln_dblink' }->[$lnkct] = $ref;
                    $lnkct++;
                    
                # Everything else (single and multi line)
                } else {
                    #       # param/-value/tagname 
                    $annotation{ $READMAP{$tag} } .= $data.' ';
                }

            } else {
                # unknown or custom data treated with simplevalue objects
                #$self->debug("Unknown tag: $tag:\t$data\n");
                $annotation{ 'custom' }->{ $tag } .= $data.' ';
            }

        } elsif( $line =~ m{^\#=GS\s+(\S+)\s+(\w{2})\s+(\S+)}xms ) {
            # sequence annotation and data
            ($id, $tag, $data) = ($1, $2, $3);
            if ($tag eq 'AC') {
                $accession{$id} .= $data;
            } elsif ($tag eq 'DE') {
                $desc{$id} .= $data;
            }
            # Bio::Seq::Meta is not AnnotationI, so can't add seq-based
            # Annotations yet; uncomment to see what is passed by
            #else {
            #    $self->debug("Missed data: $entry");
            #}
        } elsif( $line =~ m{^\#=GR\s+(\S+)\s+(\S+)\s+([^\n]+)} ) {
            # meta strings per sequence
            ($name, $tag, $data) = ($1, $2, $3);
            $seq_meta{$name}->{$tag} .= $data;
        } elsif( $line =~ m{^\#=GC\s+(\S+)\s+([^\n]+)}xms ) {
            # meta strings per alignment
            ($tag, $data) = ($1, $2);
            $aln_meta{$tag} .= $data;
        } elsif( $line =~ m{^([^\#]\S+)\s+([A-Za-z.\-\*]+)\s*}xms ) {
            ($name,$seq) = ($1,$2);
            if( ! exists $align{$name}  ) {
                push @c2name, $name;
            }
            $align{$name} .= $seq;
        } else {
            # debugging to catch missed data; uncomment to turn on
            #$self->debug("Missed Data: $line");
        }
    }
    
    # ok... now we can make the sequences
    
    for my $name ( @c2name ) {
        if( $name =~ m{(\S+)/(\d+)-(\d+)}xms ) {
            ($seqname, $start, $end) = ($1, $2, $3);
        } else {
            $seqname=$name;
            $start = 1;
            $end = length($align{$name});
        }
        $seq = Bio::Seq::Meta->new
            ('-seq'              => $align{$name},
             '-display_id'       => $seqname,
             '-start'            => $start,
             '-end'              => $end,
             '-description'      => $desc{$name},
             '-accession_number' => $accession{$name}
             );
        if (exists $seq_meta{$name}) {
            for my $tag (sort keys %{ $seq_meta{$name} }) {
                $seq->named_meta($tag, $seq_meta{$name}->{$tag});
            }
        }
        $aln->add_seq($seq);
    }
    
    # add meta strings w/o sequence for consensus meta data
    my $ameta = Bio::Seq::Meta->new();
    for my $tag (sort keys %aln_meta) {
        $ameta->named_meta($tag, $aln_meta{$tag});
    }
    
    $aln->consensus_meta($ameta);
    
    # Make the annotation collection...
    
    my $coll = Bio::Annotation::Collection->new();
    my $factory = Bio::Annotation::AnnotationFactory->new();
    for my $tag (sort keys %annotation) {
        
        # most annotations
        if (!ref($annotation{$tag})) {
            my ($atype, $aparam, $tagname) = split q(/), $tag;
            # remove trailing newline, convert internal newlines to spaces
            $annotation{$tag} =~ s{\s+$}{}g;
            # split the READTYPE map to determine Annotation type, parameters, etc.
            if ($atype eq 'Method') {
                $aln->$aparam($annotation{$tag});
            } else {
                $factory->type("Bio::Annotation::$atype");
                $coll->add_Annotation
                ($tagname, $factory->create_object($aparam  => $annotation{$tag}));
            }
            
        } elsif ($tag eq 'custom') {
            
            for my $key (sort keys %{ $annotation{$tag} }) {
                $factory->type("Bio::Annotation::SimpleValue");
                $coll->add_Annotation(
                    $tag, $factory->create_object(-tagname => $key,
                                                  -value => $annotation{$tag}->{$key}));
            }
        
        # more complex annotations
        
        } else {
            my $atype = #($tag eq 'custom')          ? 'SimpleValue'   :
                        ($tag eq 'reference')       ? 'Reference'   :
                        ($tag eq 'aln_dblink')      ? 'DBLink'   :
                        ($tag eq 'build_command')   ? 'SimpleValue' :
                        'BadValue'; # this will cause the factory to choke
            $self->throw("Bad tag value : $tag.") if $atype eq 'BadValue';
            $factory->type("Bio::Annotation::$atype");                
            while (my $data = shift @{ $annotation{$tag} }) {
                next unless $data;
                # remove trailing spaces for concatenated data
                my %clean_data = map {
                    $data->{$_} =~ s{\s+$}{}g;
                    $_ => $data->{$_};
                    } keys %{ $data };
                my $ann = $factory->create_object(%clean_data);
                $coll->add_Annotation($tag, $ann);
                $refct++;
            }
        }
    }

    # add annotations
    $aln->annotation($coll); 
    
    #  If $end <= 0, we have either reached the end of
    #  file in <fh> or we have encountered some other error
    return if ($end <= 0);
    return $aln;
}


=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in stockholm format
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object

=cut

sub write_aln {
    # enable array of SimpleAlign objects as well (see clustalw write_aln())
    my ($self, @aln) = @_;
    for my $aln (@aln) {
    $self->throw('Need Bio::Align::AlignI object')
          if (!$aln || !($aln->isa('Bio::Align::AlignI')));

    my @anns;
    my $coll = $aln->annotation;
    my ($aln_ann, $seq_ann, $aln_meta, $seq_meta) =
       ('#=GF ', '#=GS ', '#=GC ', '#=GR' );
    $self->_print("# $STKVERSION\n\n") or return 0;
    
    # annotations first
    
    for my $param (@WRITEORDER) {
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
            } elsif ($tag eq 'XX') { # custom
                my $newtag = $ann->tagname;
                $alntag = sprintf('%-10s',$aln_ann.$newtag);
                $data = $ann->display_text;
            } elsif ($tag eq 'SQ') {
                # use the actual number, not the stored Annotation data
                $alntag = sprintf('%-10s',$aln_ann.$tag);
                $data = $aln->no_sequences;
            } else {
                $alntag = sprintf('%-10s',$aln_ann.$tag);
                $data = ref $ann ? $ann->display_text : $ann;
            }
            $text = wrap($alntag, $alntag, $data);
            $self->_print("$text\n") or return 0;
        }
    }
    
    $self->_print("\n");
    
    # now the sequences...
    
    # modified (significantly) from AlignIO::pfam
    
    my ($namestr,$seq,$add);
    
    # pad extra for meta lines
    my $maxlen = $aln->maxdisplayname_length() + 5;
    my $metalen = $aln->max_metaname_length() || 0;
    
    for $seq ( $aln->each_seq() ) {
        $namestr = $aln->displayname($seq->get_nse());
        $self->_print(sprintf("%-*s  %s\n",$maxlen+$metalen, $namestr, $seq->seq())) or return 0;
        if ($seq->isa('Bio::Seq::MetaI')) {
            for my $mname ($seq->meta_names) {
                 $self->_print(sprintf("%-*s%*s  %s\n",$maxlen, $seq_meta.' '.$namestr, $metalen,
                                       $mname, $seq->named_meta($mname))) or return 0;
            }
        }
    }
    # alignment consensus
    my $ameta = $aln->consensus_meta;
    if ($ameta) {
        for my $mname ($ameta->meta_names) {
            $self->_print(sprintf("%-*s%*s  %s\n",$maxlen, $aln_meta, $metalen,
                                  $mname, $ameta->named_meta($mname))) or return 0; 
        }
    }
    $self->_print("//\n") or return 0;
    }
    $self->flush() if $self->_flush_on_write && defined $self->_fh;
    
    return 1;
}

1;