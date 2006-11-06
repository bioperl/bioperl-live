# $Id$
#
# BioPerl module for Bio::AlignIO::stockholm

#   Based on the Bio::SeqIO::stockholm module
#       by Ewan Birney <birney@ebi.ac.uk>
#       and Lincoln Stein  <lstein@cshl.org>
#
#       and the SimpleAlign.pm module of Ewan Birney
#
# Copyright Peter Schattner
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
stockholm flat file databases.

Note that many Stockholm sequence files contain additional sequence-based
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

This module is currently being refactored to incorporate Meta data for
sequences and alignments.  Annotations are also now added for alignments.

Note that sequence names in the alignment are cut off at 

=head1 FEEDBACK

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHORS - Peter Schattner, Chris Fields

Email: schattner@alum.mit.edu, cjfields-at-uiuc-dot-edu

=head1 CONTRIBUTORS

Chris Fields, cjfields-at-uiuc-dot-edu
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
our $SEQNAMELENGTH = 30;
our $INTERVEAVED = 80;

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
            'DR'   => 'SimpleValue/-value/database_reference',
            'CC'   => 'Comment/-text/alignment_comment',
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
            'database_reference'    =>  'DR/SimpleValue',
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
            'sec_structure_source'  =>  'SS/SimpleValue' 
            );

# Add mapping hash do deal with alignment annotations
# GS lines => sequence data
# GC lines are consensus-like, need new method for meta string to be stored in SimpleAlign
# GF lines => Bio::Annotation objects in a Bio::Annotation::Collection
# GR lines => Meta data stored in Bio::Seq::Meta objects 

sub _initialize {
    my ( $self, @args ) = @_;
    $self->SUPER::_initialize(@args);
    my ( $build_anno ) = $self->_rearrange( [qw(BUILD_ANNO)], @args );
    $build_anno     && $self->build_annotation($build_anno);
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

    my ($start, $end, $id, $name, $seqname, $seq, $count, $tag, $data);
    my $seen_rc;
    my $refct = 0;
    my $bct = 0;
    my @c2name;
    my (%align, %accession, %desc, %seq_meta, %aln_meta, %annotation);

    # in stockholm format, every non-blank line that does not start
    # with '#=' is an alignment segment; the '#=' lines are mark up lines.
    # Of particular interest are the '#=GF <name/st-ed> AC <accession>'
    # lines, which give accession numbers for each segment

    my $aln =  Bio::SimpleAlign->new(-source => 'stockholm');
    while( defined($line = $self->_readline) ) {
        next unless $line =~ /\w+/;
        if ($line =~ /^#\s*STOCKHOLM\s+/) {
            last;
        } else {
            $self->throw("Not Stockholm format: Expecting \"# STOCKHOLM 1.0\"; Found \"$_\"");
        }
    }
    
    READLINE:
    while( defined($line = $self->_readline) ) {
        #skip empty lines
        next if $line =~ /^\s+$/;
        
        # Double slash (//) signals end of file.
        last if $line =~ m{^//};
        
        # GF/GS lines, by convention, should be at the top of the alignment
        if ($line =~ m{^\#=GF\s+(\w{2})\s+([^\n]*)$}xms) {

            # alignment annotation
            ($tag, $data) = ($1, $2);
            if (exists $READMAP{$tag}) {

                # reference data
                if (index($tag, 'R') == 0) {
                    # comments come before numbering, tricky
                    $refct++ if ( ($tag eq 'RN' && !$seen_rc) || $tag eq 'RC');
                    $seen_rc = 1 if $tag eq 'RC';
                    # Don't need
                    next READLINE if $tag eq 'RN';
                    #                           # of ref       parameter     
                    $annotation{ 'reference' }->[$refct]->{ $READMAP{$tag} } .= $data.' ';

                # Build commands 
                } elsif ($tag eq 'BM') {
                    $bct++;
                    #                            # build cmd    parameter     
                    $annotation{ 'build_command' }->[$bct]->{ $READMAP{$tag} } .= $data.' ';
                    
                # Everything else (for now)    
                } else {
                    #       # param/-value/tagname 
                    $annotation{ $READMAP{$tag} } .= $data.' ';
                }

            } else {
                $self->debug("Unknown tag: $tag:\t$data");
            }

        } elsif( $line =~ m{^\#=GS\s+(\S+)\s+(\w{2})\s+(\S+)}xms ) {
            # sequence annotation and data
            ($id, $tag, $data) = ($1, $2, $3);
            if ($tag eq 'AC') {
                $accession{$id} .= $data;
            } elsif ($tag eq 'DE') {
                $desc{$id} .= $data;
            }
            # Does not catch other GS-based tags (DR, OS, OC, LO) yet
            #else {
            #    $self->debug("Missed data: $entry");
            #}
        } elsif( $line =~ m{^\#=GR\s+(\S+)\s+(\w+)\s+([^\n]+)} ) {
            # meta strings per sequence
            ($name, $tag, $data) = ($1, $2, $3);
            $seq_meta{$name}->{$tag} .= $data;
        } elsif( $line =~ m{^\#=GC\s+(\S+)\s+(\S+)}xms ) {
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
            $self->debug("Missed Data: $line");
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
    
    my $aln_meta;
    
    # Make the annotation collection...
    
    my $coll = Bio::Annotation::Collection->new();

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
                my $factory = Bio::Annotation::AnnotationFactory->new(
                    -type => "Bio::Annotation::$atype");
                my $ann = $factory->create_object(-tagname    => $tagname,
                                                  $aparam  => $annotation{$tag});
                $coll->add_Annotation($ann);
            }
        }
        else {
            my $data;
            my $atype = ($tag eq 'reference')       ? 'Reference'   :
                        ($tag eq 'build_command')   ? 'SimpleValue' :
                        'BadValue'; # this will cause the factory to choke
            $self->throw("Bad tag value : $tag.") if $atype eq 'BadValue';
            my $factory = Bio::Annotation::AnnotationFactory->new(
                -type => "Bio::Annotation::$atype");                
            while(@{ $annotation{$tag} }) {
                $data = shift @{ $annotation{$tag} };
                next unless $data;
                # remove trailing spaces
                my %clean_data = map {
                    $data->{$_} =~ s{\s+$}{}g;
                    $_ => $data->{$_};
                    } keys %{ $data };
                my $ann = $factory->create_object(%clean_data);
                $coll->add_Annotation($tag,$ann);
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
 Function: writes the $aln object into the stream in stockholm format  ###Not yet implemented!###
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object


=cut

sub write_aln {
    # enable array of SimpleAlign objects as well (see clustalw write_aln())
    my ($self, $aln) = @_;
    $self->throw('Need Bio::Align::AlignI object')
          if (!$aln || !($aln->isa('Bio::Align::AlignI')));

    my (@anns, @ann_objs);
    my $coll = $aln->annotation;
    my ($aln_ann, $seq_ann, $aln_meta, $seq_meta) =
       ('#=GF ', '#=GS ', '#=GC ', '#=GR' );
    $self->_print("# $STKVERSION\n\n") or return;
    
    # annotations first
    
    for my $param (@WRITEORDER) {
        # no point in going through this if there is no annotation!
        last if !$coll;
        # alignment annotations
        my $ct = 1;
        $self->throw("Bad parameter: $param") if !exists $WRITEMAP{$param};
        my ($tag, $key) = split q(/), $WRITEMAP{$param};
        if ($key eq 'Method') {
            push @anns, $aln->$param;
        } else {
            @anns = $coll->get_Annotations($param);
        }
        my $rn = 1;
        while (my $ann = shift @anns) {
            # using Text::Wrap::wrap() for word wrap
            # references
            if ($tag eq 'RX') {
                for my $rkey (qw(ref_comment ref_number ref_pubmed
                              ref_title ref_authors ref_location)) {
                    my $text;
                    my ($newtag, $method) = split q(/), $WRITEMAP{$rkey};
                    if ($rkey eq 'ref_number') {
                        $text = wrap($aln_ann.$newtag.' ', $aln_ann.$newtag.' ', "[$rn]");
                    } else {
                        $text = wrap($aln_ann.$newtag.' ', $aln_ann.$newtag.' ', $ann->$method)
                            if $ann->$method;
                    }
                    #$self->_print("REFERENCE\n");
                    $self->_print("$text\n") or return if $text;
                }
                $rn++;
            } else {
                my $text = wrap($aln_ann.$tag.' ', $aln_ann.$tag.' ', $ann);
                $self->_print("$text\n") or return;
            }
        }
    }
    
    # now the sequences...
    
    # copied from pfam.pm for now; will modify to interleave sequences, print
    # meta data, etc.
    
    my ($namestr,$seq,$add);
    
    # 
    my $maxlen = $aln->maxdisplayname_length() + 10;
    
    foreach $seq ( $aln->each_seq() ) {
        $namestr = $aln->displayname($seq->get_nse());
        $self->_print(sprintf("%-*s  %s\n",$maxlen, $namestr, $seq->seq())) or return;
        if ($seq->isa('Bio::Seq::MetaI')) {
            for my $mname ($seq->meta_names) {
                 $self->_print(sprintf("%-*s%5s  %s\n",$maxlen-5, $seq_meta.' '.$namestr, $mname, $seq->named_meta($mname))) or return;
            }
        }
    }
    $self->flush() if $self->_flush_on_write && defined $self->_fh;
    
    $self->_print("//\n") or return;
    return;
}

1;