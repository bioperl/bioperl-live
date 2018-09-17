# Let the code begin...

package Bio::AlignIO::Handler::GenericAlignHandler;

use strict;
use warnings;

use Bio::Annotation::Collection;
use Bio::Annotation::Comment;
use Bio::Annotation::SimpleValue;
use Bio::Annotation::Target;
use Bio::Annotation::DBLink;
use Bio::Annotation::Reference;
use Bio::SimpleAlign;
use Data::Dumper;

use base qw(Bio::Root::Root Bio::HandlerBaseI);

# only stockholm is defined for now...
my %HANDLERS = (
    # stockholm has sequence and alignment specific annotation; this
    'stockholm'   => {
        'CONSENSUS_META'    => \&_generic_consensus_meta,
        'SEQUENCE'          => \&_generic_metaseq,
        'NAMED_META'        => \&_generic_metaseq,
        'ACCESSION'         => \&_generic_store,
        'ALPHABET'          => \&_generic_store,
        'ID'                => \&_generic_store,
        'DESCRIPTION'       => \&_generic_store,
        'REFERENCE'         => \&_generic_reference,
        'DBLINK'            => \&_stockholm_target,
        'DATABASE_COMMENT'  => \&_generic_comment,
        'ALIGNMENT_COMMENT' => \&_generic_comment,
        '_DEFAULT_'         => \&_generic_simplevalue
        },
    );

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($format, $verbose) = $self->_rearrange([qw(FORMAT VERBOSE)], @args);
    $self->throw("Must define alignment record format") if !$format;
    $verbose   && $self->verbose($verbose);
    $self->format($format);
    $self->handler_methods();
    # if we intend at a later point we can add a Builder 
    #$builder  &&  $self->alignbuilder($builder);
    return $self;
}

sub handler_methods {
    my $self = shift;
    if (!($self->{'handlers'})) {
        $self->throw("No handlers defined for alignment format ",$self->format)
            unless exists $HANDLERS{$self->format};
        $self->{'handlers'} = $HANDLERS{$self->format};
    }
    return ($self->{'handlers'});
}

sub data_handler {
    my ($self, $data) = @_;
    my $nm = $data->{NAME} || $self->throw("No name tag defined!");
    # this should handle data on the fly w/o caching; any caching should be 
    # done in the driver!
    my $method = (exists $self->{'handlers'}->{$nm}) ? ($self->{'handlers'}->{$nm}) :
                (exists $self->{'handlers'}->{'_DEFAULT_'}) ? ($self->{'handlers'}->{'_DEFAULT_'}) :
                undef;
    if (!$method) {
        $self->debug("No handler defined for $nm\n");
        return;
    };
    $self->$method($data);
}

sub reset_parameters {
    my $self = shift;
    $self->{'_params'} = undef;
    $self->{'_nse_cache'} = undef;
    $self->{'_features'} = undef;
}

sub format {
    my $self = shift;
    if (@_) {
        my $format = lc shift;
        $self->throw("Format $format not supported") unless exists $HANDLERS{$format};
        $self->{'_alignformat'} = $format;
    };
    return $self->{'_alignformat'};
}

sub get_params {
    my ($self, @ids) = @_;
    my $data;
    if (scalar(@ids)) {
        for my $id (@ids) {
            if (!index($id, '-')==0) {
                $id = '-'.$id ;
            }
            $data->{$id} = $self->{'_params'}->{$id} if (exists $self->{'_params'}->{$id});
        }
        $data ||= {};
    } else {
        $data = $self->{'_params'};
    }
    return $data;
}

sub set_params {
    shift->throw('Not implemented yet!');
}

sub build_alignment {
    my $self = shift;
    my %init;
    $self->process_seqs;
    my $param = $self->get_params;
    if (defined $param->{-seqs}) {
        return Bio::SimpleAlign->new(%$param, -source => $self->format);
    }
    return;
}

sub annotation_collection {
    my ($self, $coll) = @_;
    if ($coll) {
        $self->throw("Must have Bio::AnnotationCollectionI ".
                     "when explicitly setting annotation_collection()")
            unless (ref($coll) && $coll->isa('Bio::AnnotationCollectionI'));
        $self->{'_params'}->{'-annotation'} = $coll;
    } elsif (!exists($self->{'_params'}->{'-annotation'})) {
        $self->{'_params'}->{'-annotation'} = Bio::Annotation::Collection->new()
    }
    return $self->{'_params'}->{'-annotation'};
}

sub seq_annotation_collection {
    my ($self, $coll) = @_;
    if ($coll) {
        $self->throw("Must have Bio::AnnotationCollectionI ".
                     "when explicitly setting seq_annotation_collection()")
            unless (ref($coll) && $coll->isa('Bio::AnnotationCollectionI'));
        $self->{'_params'}->{'-seq_annotation'} = $coll;
    } elsif (!exists($self->{'_params'}->{'-seq_annotation'})) {
        $self->{'_params'}->{'-seq_annotation'} = Bio::Annotation::Collection->new()
    }
    return $self->{'_params'}->{'-seq_annotation'};
}

sub process_seqs {
    my $self = shift;

    my $data = $self->get_params(qw(-seqs -seq_class -consensus_meta));
    my $class = $data->{-seq_class} || 'Bio::LocatableSeq';
    # cache classes loaded already
    if (!exists($self->{'_loaded_modules'}->{$class})) {
        $self->_load_module($class);
        $self->{'_loaded_modules'}->{$class}++;
    }
    # process any meta sequence data
    if ( $data->{-consensus_meta} && !UNIVERSAL::isa($data->{-consensus_meta},'Bio::Seq::Meta')) {
        my $ref = $data->{-consensus_meta};
        if (!exists($self->{'_loaded_modules'}->{'Bio::Seq::Meta'})) {
            $self->_load_module('Bio::Seq::Meta');
            $self->{'_loaded_modules'}->{'Bio::Seq::Meta'}++;
        }
        my $ms = Bio::Seq::Meta->new();
        for my $tag (sort keys %{$ref}) {
            $ms->named_meta($tag, $ref->{$tag});
        }
        $self->{'_params'}->{'-consensus_meta'} = $ms;
    }
    # this should always be an array ref!
    for my $seq (@{$data->{-seqs}}) {
        next if (UNIVERSAL::isa($seq,'Bio::LocatableI'));
        # process anything else
        $self->_from_nse($seq) if $seq->{NSE};
        if (UNIVERSAL::isa($seq,'HASH')) {
            my %param;
            for my $p (keys %$seq) {
                $param{'-'.lc $p} = $seq->{$p} if exists $seq->{$p};
            }
            my $ls = $class->new(%param);
            # a little switcheroo to attach the sequence
            # (though using it to get seq() doesn't work correctly yet!)
            if (defined $seq->{NSE} &&
                exists $self->{'_features'} &&
                exists $self->{'_features'}->{ $seq->{NSE} }) {
                for my $feat (@{ $self->{'_features'}->{ $seq->{NSE} } }) {
                    push @{ $self->{'_params'}->{'-features'} }, $feat;
                    $feat->attach_seq($ls);
                }
            }
            $seq = $ls;
        }
    }
}

####################### SEQUENCE HANDLERS #######################

# any sequence data for a Bio::Seq::Meta
sub _generic_metaseq {
    my ($self, $data) = @_;
    return unless $data;
    $self->throw("No alignment position passed") if !exists($data->{BLOCK_LINE});
    $self->throw("Alignment position must be an index greater than 0") if $data->{BLOCK_LINE} < 1;
    $self->{'_params'}->{'-seq_class'} = 'Bio::Seq::Meta';
    my $index = $data->{BLOCK_LINE} - 1;
    if (my $nse = $self->{'_params'}->{'-seqs'}->[$index]->{NSE}) {
        $self->throw("NSE in passed data doesn't match stored data in same position: $nse") unless $nse eq $data->{NSE};
    } else {
        $self->{'_params'}->{'-seqs'}->[$index]->{NSE} = $data->{NSE};
    }
    if ($data->{NAME} eq 'SEQUENCE') {
        $self->{'_params'}->{'-seqs'}->[$index]->{SEQ} .= $data->{DATA};
    } elsif ($data->{NAME} eq 'NAMED_META') {
        $self->{'_params'}->{'-seqs'}->[$index]->{NAMED_META}->{$data->{META_TAG}} .= $data->{DATA};
    }
}

sub _generic_consensus_meta {
    my ($self, $data) = @_;
    return unless $data;
    if ($data->{NAME} eq 'CONSENSUS_META') {
        $self->{'_params'}->{'-consensus_meta'}->{$data->{META_TAG}} .= $data->{DATA};
    }
}

# any sequence data for a Bio::LocatableSeq
sub _generic_locatableseq {
    my ($self, $data) = @_;
    return unless $data;
    $self->throw("No alignment position passed") if !exists($data->{BLOCK_LINE});
    $self->throw("Alignment position must be an index greater than 0") if $data->{BLOCK_LINE} < 1;
    my $index = $data->{BLOCK_LINE} - 1;
    if (my $nse = $self->{'_params'}->{'-seqs'}->[$index]->{NSE}) {
        $self->throw("NSE in passed data doesn't match stored data in same position: $nse") if $nse ne $data->{NSE};
    } else {
        $self->{'_params'}->{'-seqs'}->[$index]->{NSE} = $data->{NSE};
    }
    if ($data->{NAME} eq 'SEQUENCE') {
        $self->{'_params'}->{'-seqs'}->[$index]->{SEQ} .= $data->{DATA};
    }
}

####################### RAW DATA HANDLERS #######################

# store by data name (ACCESSION, ID, etc), which can be mapped to the
# appropriate alignment or sequence parameter
sub _generic_store {
    my ($self, $data) = @_;
    return unless $data;
    if ($data->{ALIGNMENT}) {
        $self->{'_params'}->{'-'.lc $data->{NAME}} = $data->{DATA};
    } else {
        $self->{'_params'}->{'-seq_'.lc $data->{NAME}}->{$data->{NSE}} = $data->{DATA}
    }
}

sub _generic_reference {
    my ($self, $data) = @_;
    my $ref = Bio::Annotation::Reference->new(-title => $data->{TITLE},
                                              -authors => $data->{AUTHORS},
                                              -pubmed => $data->{PUBMED},
                                              -location => $data->{JOURNAL},
                                              -tagname  => lc $data->{NAME});
    $self->annotation_collection->add_Annotation($ref);
}

sub _generic_simplevalue {
    my ($self, $data) = @_;
    my $sv = Bio::Annotation::SimpleValue->new(-value => $data->{DATA},
                                            -tagname  => lc $data->{NAME});
    $self->annotation_collection->add_Annotation($sv);
}

sub _generic_comment {
    my ($self, $data) = @_;
    my $comment = Bio::Annotation::Comment->new(-type => lc $data->{NAME},
                                                -text => $data->{DATA},
                                                -tagname  => lc $data->{NAME});
    $self->annotation_collection->add_Annotation($comment);
}

# Some DBLinks in Stockholm format are unique, so a unique handler for them
sub _stockholm_target {
    my ($self, $data) = @_;
    # process database info
    $self->_from_stk_dblink($data);
    my $comment;
    # Bio::Annotation::Target is now a DBLink, but has additional (RangeI) 
    # capabilities (for PDB data)
    my $dblink = Bio::Annotation::Target->new(
        -database => $data->{DBLINK_DB},
        -primary_id => $data->{DBLINK_ACC},
        -optional_id => $data->{DBLINK_OPT},
        -start => $data->{DBLINK_START},
        -end => $data->{DBLINK_END},
        -strand => $data->{DBLINK_STRAND},
        -comment => $comment,
        -tagname => 'dblink',
    );
    if ($data->{ALIGNMENT}) {
        # Alignment-specific DBLinks
        $self->annotation_collection->add_Annotation($dblink);
    } else {
        # Sequence-specific DBLinks
        # These should come with identifying information of some sort
        # (ID/START/END/STRAND).  Make into a SeqFeature (SimpleAlign is
        # FeatureHolderI) spanning the length acc. to the NSE. Add the DBLink as
        # Annotation specific to that SeqFeature, store in an internal hash by
        # NSE so we can tie the LocatableSeq to the proper Features
        $self->_from_nse($data) if $data->{NSE};
        $self->throw("Must supply an sequence DISPLAY_ID or NSE for sequence-related
            DBLinks") unless $data->{ACCESSION_NUMBER} || $data->{DISPLAY_ID};
        my $sf = Bio::SeqFeature::Generic->new(-seq_id => $data->{DISPLAY_ID},
                                               -accession_number => $data->{ACCESSION_NUMBER},
                                               -start => $data->{START},
                                               -end => $data->{END},
                                               -strand => $data->{STRAND}
                                               );
        $sf->annotation->add_Annotation($dblink);
        # index by NSE
        push @{ $self->{'_features'}->{ $data->{NSE} } }, $sf;
        #$self->seq_annotation_collection->add_Annotation($dblink);
    }
}

####################### HELPER METHODS #######################

# returns ACCESSION VERSION START END STRAND ALPHABET
# cached for multiple lookups, should reset in between uses
sub _from_nse {
    my ($self, $data) = @_;
    return unless my $nse = $data->{NSE};
    $data->{ALPHABET} =  $self->get_params('-alphabet')->{'-alphabet'} || 'protein';
    # grab any accessions if present, switch out with ACCESSION from NSE
    # (move that to primary_id)
    my $new_acc;
    if (exists $self->{'_params'}->{'-seq_accession'}) {
        $new_acc = $self->{'_params'}->{'-seq_accession'}->{$data->{NSE}};
    }
    if ($nse =~ m{(\S+?)(?:\.(\d+))?/(\d+)-(\d+)}xmso) {
        my $strand = $data->{ALPHABET} eq 'dna' || $data->{ALPHABET} eq 'rna' ? 1 : undef;
        my ($start, $end) = ($3, $4);
        if ($start > $end) {
            ($start, $end, $strand) = ($end, $start, -1);
        }
        $data->{ACCESSION_NUMBER} = $new_acc || $1;
        $data->{DISPLAY_ID} = $1;
        $data->{VERSION} = $2;
        $data->{START} = $start;
        $data->{END} = $end;
        $data->{STRAND} = $strand;
    } else {
        # we can parse for version here if needed
        $data->{DISPLAY_ID} = $data->{NSE};
    }
}

# this will probably be split up into subhandlers based on Record/DB 
sub _from_stk_dblink {
    my ($self, $data) = @_;
    return unless my $raw = $data->{DATA};
    my @rawdata = split(m{\s*;\s*}, $raw);
    my %dblink_data;
    if ($rawdata[0] eq 'PDB') {
        # fix for older Stockholm PDB range format
        if (scalar(@rawdata) == 3 && $rawdata[2] =~ m{-}) {
            @rawdata[2,3] = split('-',$rawdata[2],2);
        }
        $self->throw("Not standard PDB form: ".$data->{DATA}) if scalar(@rawdata) != 4;
        my ($main, $chain) = split(m{\s+}, $rawdata[1]);
        %dblink_data = (
            DBLINK_DB => $rawdata[0],
            DBLINK_ACC => $main,
            DBLINK_OPT => $chain || '',
            DBLINK_START => $rawdata[2],
            DBLINK_END => $rawdata[3]
        );
    } elsif ($rawdata[0] eq 'SCOP') {
        $self->throw("Not standard SCOP form: ".$data->{DATA}) if scalar(@rawdata) != 3;
        %dblink_data = (
            DBLINK_DB => $rawdata[0],
            DBLINK_ACC => $rawdata[1],
            DBLINK_OPT => $rawdata[2],
        );        
    } else {
        $self->warn("Some data missed: ".$data->{DATA}) if scalar(@rawdata) > 2;
        %dblink_data = (
            DBLINK_DB => $rawdata[0],
            DBLINK_ACC => $rawdata[1],
        );        
    }
    while (my ($k, $v) = each %dblink_data) {
        $data->{$k} = $v if $v;
    }    
}

1;

__END__

# $Id: GenericAlignHandler.pm 14816 2008-08-21 16:00:12Z cjfields $
#
# BioPerl module for Bio::AlignIO::Handler::GenericAlignHandler
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chris Fields
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
#
# Documentation after the __END__ marker

=head1 NAME

Bio::AlignIO::Handler::GenericAlignHandler - Bio::HandlerI-based
generic data handler class for alignment-based data

=head1 SYNOPSIS

  # MyHandler is a GenericAlignHandler object.
  # inside a parser (driver) constructor....

  $self->alignhandler($handler || MyHandler->new(-format => 'stockholm'));

  # in next_aln() in driver...

  $hobj = $self->alignhandler();

  # roll data up into hashref chunks, pass off into Handler for processing...

  $hobj->data_handler($data);

  # or retrieve Handler methods and pass data directly to Handler methods...

  my $hmeth = $hobj->handler_methods;

  if ($hmeth->{ $data->{NAME} }) {
      my $mth = $hmeth->{ $data->{NAME} };
      $hobj->$mth($data);
  }

=head1 DESCRIPTION

This is an experimental implementation of a alignment-based HandlerBaseI parser
and may change over time. It is possible that the way handler methods are set up
will change over development to allow more flexibility. 

Standard Developer caveats:

Here thar be dragoons...

Consider yourself warned!

=head2 NOTES

As in the SeqIO Handler object (still in development), data is passed in as
chunks. The Annotation and SeqFeatures are essentially the same as the SeqIO
parser; the significant difference is that data hash being passed could pertain
to either the alignment or to a specific sequence, so an extra tag may be needed
to disambiguate between the two in some cases. Here I use the ALIGNMENT tag as a
boolean flag: it must be present and set to 0 for the data to be tagged for
Bio::LocatableSeq or similar (in all other cases it is assumed to be for the
alignment). In some cases this will not matter (the actual sequence data, for
instance) but it is highly recommmended adding this tag in to prevent possible
ambiguities.

This is the current Annotation data chunk (via Data::Dumper):

  $VAR1 = {
            'NAME' => 'REFERENCE',
            'DATA' => '1  (bases 1 to 10001)'
            'AUTHORS' => 'International Human Genome Sequencing Consortium.'
            'TITLE' => 'The DNA sequence of Homo sapiens'
            'JOURNAL' => 'Unpublished (2003)'
            'ALIGNMENT' => 1,
          };

In the case of LocatableSeqs, one can pass them in as follows for simplicity
(note the block line):

  $VAR1 = {
            'NAME' => 'SEQUENCE', 
            'BLOCK_LINE' => 0,
            'NSE' => 'Q7WNI7_BORBR/113-292',
            'ALPHABET' => 'protein',
            'DATA' => 'VALILGVYRRL...CYVNREM..RAG....QW',
            'ALIGNMENT' => 0            
          };

This can be done as the parser parses each block instead of parsing all the
blocks and then passing them in one at a time; the handler will store the
sequence data by the block line in an internal hash, concatenating them along
the way.  This behaviour is b/c the alignment building step requires that
the sequence be checked for start/end/strand, possible meta sequence, optional
accession, etc.

Similarly, a Meta sequence line can be passed in as follows:

  $VAR1 = {
            'NAME' => 'NAMED_META',
            'BLOCK_LINE' => 0,
            'NSE' => 'Q7WNI7_BORBR/113-292',
            'META_KEY' => 'pAS',
            'DATA' => '................................',
            'ALIGNMENT' => 0
          };

The meta sequence will be checked against the NSE for the block position and
stored based on the meta tag. A meta sequence does not have to correspond to a
real sequence. At this time, unique meta sequence tags must be used for each
sequence or they will be overwritten (this may change).

An alignment consensus string: 

  $VAR1 = {
            'NAME' => 'CONSENSUS',
            'DATA' => 'VALILGVYRRL...CYVNREM..RAG....QW',
            'ALIGNMENT' => 1
          };

A consensus meta sequence:

  $VAR1 = {
            'NAME' => 'CONSENSUS_META',
            'META_KEY' => 'pAS',
            'DATA' => '................................',
            'ALIGNMENT' => 1
          };

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Chris Fields

Email cjfields at bioperl dot org

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut

=head2 new

 Title   :  new
 Usage   :  
 Function:  
 Returns :  
 Args    :  -format    Sequence format to be mapped for handler methods
            -builder   Bio::Seq::SeqBuilder object (normally defined in
                       SequenceStreamI object implementation constructor)
 Throws  :  On undefined '-format' sequence format parameter
 Note    :  Still under heavy development

=cut

=head1 L<Bio::HandlerBaseI> implementing methods

=head2 handler_methods

 Title   :  handler_methods
 Usage   :  $handler->handler_methods('GenBank')
            %handlers = $handler->handler_methods();
 Function:  Retrieve the handler methods used for the current format() in
            the handler.  This assumes the handler methods are already
            described in the HandlerI-implementing class.
 Returns :  a hash reference with the data type handled and the code ref
            associated with it.
 Args    :  [optional] String representing the sequence format.  If set here
            this will also set sequence_format()
 Throws  :  On unimplemented sequence format in %HANDLERS

=cut

=head2 data_handler

 Title   :  data_handler
 Usage   :  $handler->data_handler($data)
 Function:  Centralized method which accepts all data chunks, then distributes
            to the appropriate methods for processing based on the chunk name
            from within the HandlerBaseI object.

            One can also use 
 Returns :  None
 Args    :  an hash ref containing a data chunk.  

=cut

=head2 reset_parameters

 Title   :  reset_parameters
 Usage   :  $handler->reset_parameters()
 Function:  Resets the internal cache of data (normally object parameters for
            a builder or factory)
 Returns :  None
 Args    :  None

=cut

=head2 format

 Title   :  format
 Usage   :  $handler->format('GenBank')
 Function:  Get/Set the format for the report/record being parsed. This can be
            used to set handlers in classes which are capable of processing
            similar data chunks from multiple driver modules.
 Returns :  String with the sequence format
 Args    :  [optional] String with the sequence format
 Note    :  The format may be used to set the handlers (as in the
            current GenericRichSeqHandler implementation)

=cut

=head2 get_params

 Title   :  get_params
 Usage   :  $handler->get_params('-species')
 Function:  Convenience method used to retrieve the specified
            parameters from the internal parameter cache
 Returns :  Hash ref containing parameters requested and data as
            key-value pairs.  Note that some parameter values may be 
            objects, arrays, etc.
 Args    :  List (array) representing the parameters requested

=cut

=head2 set_params

 Title   :  set_params
 Usage   :  $handler->set_param({'-seqs' => $seqs})
 Function:  Convenience method used to set specific parameters
 Returns :  None
 Args    :  Hash ref containing the data to be passed as key-value pairs

=cut

=head1 Methods unique to this implementation

=head2 build_alignment

 Title   :  build_alignment
 Usage   :  
 Function:  
 Returns :  a Bio::SimpleAlign
 Args    :
 Throws  :
 Note    :  This may be replaced by a Builder object at some point 

=cut

=head2 annotation_collection

 Title   :  annotation_collection
 Usage   :  
 Function:  
 Returns :  
 Args    :
 Throws  :
 Note    :  

=cut

=head2 seq_annotation_collection

 Title   :  seq_annotation_collection
 Usage   :  
 Function:  
 Returns :  
 Args    :
 Throws  :
 Note    :  

=cut

=head2 process_seqs

 Title   :  process_seqs
 Usage   :  $handler->process_seqs;
 Function:  checks internal sequences to ensure they are converted over
            to the proper Bio::AlignI-compatible sequence class
 Returns :  1 if successful
 Args    :  none

=cut
