#
# BioPerl module for Bio::Tools::ERPIN
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chris Fields <cjfields-at-uiuc-dot-edu>
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::ERPIN -  a parser for ERPIN output

=head1 SYNOPSIS

  use Bio::Tools::ERPIN;
  my $parser = Bio::Tools::ERPIN->new( -file => $rna_output,
                                      -motiftag => 'protein_bind'
                                      -desctag => 'TRAP_binding');
  #parse the results
  while( my $motif = $parser->next_prediction) {
    # do something here
  }

=head1 DESCRIPTION

Parses raw ERPIN output.

This module is not currently complete.  As is, it will parse raw
ERPIN long format output and pack information into
Bio::SeqFeature::Generic objects.  

Several values have also been added in the 'tag' hash.  These can be
accessed using the following syntax:

  my ($entry) = $feature->get_Annotations('SecStructure');

Added tags are : 
   tset         - training set used for the sequence
   tsetdesc     - training set description line
   cutoff       - cutoff value used
   database     - name of database
   dbdesc       - description of database
   dbratios     - nucleotide ratios of database (used to calculate evalue)
   descline     - entire description line (in case the regex used for
                  sequence ID doesn't adequately catch the name
   accession    - accession number of sequence (if present)
   logodds      - logodds score value
   sequence     - sequence from hit, separated based on training set

See t/ERPIN.t for example usage.

At some point a more complicated feature object may be used to support
this data rather than forcing most of the information into tag/value
pairs in a SeqFeature::Generic.  This will hopefully allow for more
flexible analysis of data (specifically RNA secondary structural
data).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Chris Fields

Email cjfields-at-uiuc-dot-edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::ERPIN;
use strict;

use Bio::SeqFeature::Generic;

use base qw(Bio::Tools::AnalysisResult);

our($MotifTag,$SrcTag,$DescTag) = qw(misc_binding ERPIN erpin);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::ERPIN->new();
 Function: Builds a new Bio::Tools::ERPIN object 
 Returns : an instance of Bio::Tools::ERPIN
 Args    : -fh/-file for input filename
           -motiftag => primary tag used in gene features (default 'misc_binding')
           -desctag => tag used for display_name name (default 'erpin')
           -srctag  => source tag used in all features (default 'ERPIN')

=cut

sub _initialize {
    my($self,@args) = @_;
    $self->warn('Use of this module is deprecated.  Use Bio::SearchIO::erpin instead');    
    $self->SUPER::_initialize(@args);
    my ($motiftag,$desctag,$srctag) =  $self->SUPER::_rearrange([qw(MOTIFTAG
                                                                    DESCTAG
                                                                    SRCTAG
                                   )],
                                    @args);
    $self->motif_tag(defined $motiftag ? $motiftag : $MotifTag);
    $self->source_tag(defined $srctag ? $srctag : $SrcTag);
    $self->desc_tag(defined $desctag ? $desctag : $DescTag);
    foreach (qw(_tset _tset_desc _cutoff _db _db_desc
               _db_ratios _eval_cutoff _seqid _secacc _seqdesc )) {
        $self->{$_}='';
    }
}

=head2 motif_tag

 Title   : motiftag
 Usage   : $obj->motiftag($newval)
 Function: Get/Set the value used for 'motif_tag', which is used for setting the
           primary_tag.
           Default is 'misc_binding' as set by the global $MotifTag.
           'misc_binding' is used here because a conserved RNA motif is capable
           of binding proteins (regulatory proteins), antisense RNA (siRNA),
           small molecules (riboswitches), or nothing at all (tRNA,
           terminators, etc.).  It is recommended that this be changed to other
           tags ('misc_RNA', 'protein_binding', 'tRNA', etc.) where appropriate.
           For more information, see:
           http://www.ncbi.nlm.nih.gov/collab/FT/index.html
 Returns : value of motif_tag (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub motif_tag{
    my $self = shift;

    return $self->{'motif_tag'} = shift if @_;
    return $self->{'motif_tag'};
}

=head2 source_tag

 Title   : source_tag
 Usage   : $obj->source_tag($newval)
 Function: Get/Set the value used for the 'source_tag'.
           Default is 'ERPIN' as set by the global $SrcTag
 Returns : value of source_tag (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub source_tag{
    my $self = shift;

    return $self->{'source_tag'} = shift if @_;
    return $self->{'source_tag'};
}

=head2 desc_tag

 Title   : desc_tag
 Usage   : $obj->desc_tag($newval)
 Function: Get/Set the value used for the query motif.  This will be placed in
           the tag '-display_name'.  Default is 'erpin' as set by the global
           $DescTag.  Use this to manually set the descriptor (motif searched for).
           Since there is no way for this module to tell what the motif is from the
           name of the descriptor file or the ERPIN output, this should
           be set every time an ERPIN object is instantiated for clarity
 Returns : value of exon_tag (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub desc_tag{
    my $self = shift;

    return $self->{'desc_tag'} = shift if @_;
    return $self->{'desc_tag'};
}

=head2 analysis_method

 Usage     : $obj->analysis_method();
 Purpose   : Inherited method. Overridden to ensure that the name matches
             /ERPIN/i.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method { 
#-------------
    my ($self, $method) = @_;  
    if($method && ($method !~ /ERPIN/i)) {
    $self->throw("method $method not supported in " . ref($self));
    }
    return $self->SUPER::analysis_method($method);
}

=head2 next_feature

 Title   : next_feature
 Usage   : while($gene = $obj->next_feature()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the ERPIN result
           file. Call this method repeatedly until FALSE is returned.
           The returned object is actually a SeqFeatureI implementing object.
           This method is required for classes implementing the
           SeqAnalysisParserI interface, and is merely an alias for 
           next_prediction() at present.
 Returns : A Bio::Tools::Prediction::Gene object.
 Args    : None (at present)

=cut

sub next_feature {
    my ($self,@args) = @_;
    # even though next_prediction doesn't expect any args (and this method
    # does neither), we pass on args in order to be prepared if this changes
    # ever
    return $self->next_prediction(@args);
}

=head2 next_prediction

 Title   : next_prediction
 Usage   : while($gene = $obj->next_prediction()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the ERPIN result
           file. Call this method repeatedly until FALSE is returned.
 Returns : A Bio::Tools::Prediction::Gene object.
 Args    : None (at present)

=cut

sub next_prediction {
    my ($self) = @_;
    my ($motiftag,$srctag,$desctag) = ( $self->motif_tag,
                       $self->source_tag,
                       $self->desc_tag);
    # hit vars
    my ($strand, $start, $end, $sequence, $logodds, $score)=0;
    while($_ = $self->_readline) {
        #skip blank lines
        next if /^\s+$/;
        # parse header; there's probably a better way to do this, perhaps by
        # mapping, but this works for now...
        if(/^Training set:\s+\"(.*)\":$/) {
            $self->{'_tset'}=$1;
        }
        elsif(/\s+(\d+ sequences of length \d+)/){
            $self->{'_tset_descr'}=$1;
        }
        elsif(/^Cutoff:\s+(\S+)\s+$/) {
            $self->{'_cutoff'}=$1;
        }
        elsif(/^Database:\s+\"(.*)\"$/) {
            $self->{'_db'}=$1;
        }
        elsif(/^\s+(\d+ nucleotides to be processed in \d+ sequence)$/) {
            $self->{'_db_desc'}=$1;
        }
        elsif(/^\s+ATGC ratios:\s(\d.\d+)\s+(\d.\d+)\s+(\d.\d+)\s+(\d.\d+)$/) {
            my $atgc=sprintf("A=%0.3f T=%0.3f G=%0.3f C=%0.3f", $1, $2, $3, $4);
            $self->{'_db_ratios'}=$atgc;
        }
        elsif(/^E-value at cutoff \S+ for \S+(?:G|M|k)?b double strand data: (\S+)/) {
            $self->{'_eval_cutoff'}=$1;
        }
        # catch hit, store in private hash keys
        elsif (/^>(.*)/) {
            $self->{_seq_desc} = $1;
            if($self->{_seq_desc} =~
               /(?:P<db>gb|gi|emb|dbj|sp|pdb|bbs|ref|lcl)\|(\d+)((?:\:|\|)\w+\|(\S*.\d+)\|)?/) { 
                $self->{_seqid} = $1; # pulls out gid
                $self->{_seq_acc} = $3;
            } else {
                $self->{_seqid} = $self->{_seq_desc};
                $self->{_seq_acc} = '';
            }
        }
        # parse next hit
        elsif (/^(FW|RC)\s+\d+\s+(\d+)..(\d+)\s+(\d+.\d+)\s+(.*)/) {
            ($strand, $start, $end, $logodds, $score)=($1, $2, $3, $4, $5);
            $score =~ s/^e/1e/i;
            chomp ($sequence = $self->_readline); # grab next line, which is the sequence hit
            my $gene = Bio::SeqFeature::Generic->new(-seq_id => $self->{_seqid},
                                                      -start  => $start,
                                                      -end    => $end,
                                                      -strand => $strand eq 'FW' ? 1 : -1,
                                                      -score  => $score,
                                                      -primary_tag => $motiftag,
                                                      -source_tag  => $srctag,
                                                      -display_name => $desctag,
                                                      -tag     => {
                                                        'tset'          => $self->{_tset},
                                                        'tsetdesc'      => $self->{_tset_descr},
                                                        'cutoff'        => $self->{_cutoff},
                                                        'database'      => $self->{_db},
                                                        'dbdesc'        => $self->{_db_desc},
                                                        'dbratios'      => $self->{_db_ratios},
                                                        'descline'      => $self->{_seq_desc},
                                                        'accession'     => $self->{_seq_acc},
                                                        'logodds'       => $logodds,
                                                        'sequence'      => $sequence}
                                                    );
            return $gene;
        }
        #else {
        #    $self->debug("unrecognized line: $_");
        #}
    }
}

1;