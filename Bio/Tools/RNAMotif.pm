#
# BioPerl module for Bio::Tools::RNAMotif
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

Bio::Tools::RNAMotif - A parser for RNAMotif output

=head1 SYNOPSIS

  use Bio::Tools::RNAMotif;
  my $parser = Bio::Tools::RNAMotif->new(-file => $rna_output,
                                        -motiftag => 'protein_bind'
                                        -desctag => 'TRAP_binding');
  #parse the results
  while( my $motif = $parser->next_prediction) {
    # do something here
  }

=head1 DESCRIPTION

Parses raw RNAMotif output.  RNAMotif uses a RNA profile, consisting
of sequence and structural elements stored in a descriptor file, to
search for potential motifs in a DNA sequence file.  For more
information, see:

Macke TJ, Ecker DJ, Gutell RR, Gautheret D, Case DA, Sampath R. 
RNAMotif, an RNA secondary structure definition and search algorithm.
Nucleic Acids Res. 2001 Nov 15;29(22):4724-35. 
http://www.scripps.edu/mb/case/casegr-sh-3.5.html.

This module is not currently complete.  As is, it will parse raw
RNAMotif output (i.e. information not passed through the secondary
programs rmfmt or rm2ct) and pack information into
Bio::SeqFeature::Generic objects.  Currently, parsing extra output
utilized by the sprintf() function in an RNAMotif descriptor is not
implemented; this information is instead packed into the score tag,
which can be accessed by using the following:

  my ($score) = $feature->score; 

If the score contains anything besides a digit, it will throw a
warning that sprintf() may have been used.
Several values have also been added in the 'tag' hash.  These can be
accessed using the following syntax:

  my ($entry) = $feature->get_Annotations('secstructure');

Added tags are : 

   descline     - entire description line (in case the regex used for
                  sequence ID doesn't adequately catch the name
   descfile     - name of the descriptor file (may include path to file)
   secstrucure  - contains structural information from the descriptor
                  used as a query
   sequence     - sequence of motif, separated by spaces according to
                  matches to the structure in the descriptor (in
                  SecStructure).

See t/RNAMotif.t for example usage.

The clean_features method can also be used to return a list of seqfeatures (in a
Bio::SeqFeature::Collection object) that are within a particular region.   RNAMotif
is prone with some descriptors to returning redundant hits; an attempt to rectify
this problem is attempted with RNAMotif's companion program rmprune, which returns
the structure with the longest helices (and theoretically the best scoring structure).
However, this doesn't take into account alternative foldings which may score better.
This method adds a bit more flexibility, giving the user the ability to screen folds
based on where the feature is found and the score.  Passing a positive integer x
screens SeqFeatures based on the highest score within x bp, while a negative integer
screens based on the lowest score. So, to return the highest scoring values within
20 bp (likely using an arbitrary scroing system in the SCORE section of a descriptor
file), one could use:

  $list = $obj->clean_features(20); 

... and returning the lowest scoring structures within the same region (when the
score is based on calculated free energies from efn2) can be accomplished
by the following:

  $list = $obj->clean_features(-20);

If you wanted the best feature in a sequence, you could set this to a large number,
preferrably on that exceeds the bases in a sequence

  $list = $obj->clean_features(10000000);

Each seqfeature in the collection can then be acted upon:

  @sf = $list->get_all_features;
  for my $f (@sf) {
    # do crazy things here
  }

At some point a more complicated feature object may be used to support
this data rather than forcing most of the information into tag/value
pairs in a SeqFeature::Generic.  This will hopefully allow for more
flexible analysis of data (specifically RNA secondary structural
data).  It works for now...

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

package Bio::Tools::RNAMotif;
use strict;

use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Collection;

use base qw(Bio::Tools::AnalysisResult);

our($MotifTag,$SrcTag,$DescTag) = qw(misc_binding RNAMotif rnamotif);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::RNAMotif->new();
 Function: Builds a new Bio::Tools::RNAMotif object 
 Returns : an instance of Bio::Tools::RNAMotif
 Args    : -fh/-file for input filename
           -motiftag => primary tag used in gene features (default 'misc_binding')
           -desctag => tag used for display_name name (default 'rnamotif')
           -srctag  => source tag used in all features (default 'RNAMotif')

=cut

sub _initialize {
  my($self,@args) = @_;
  $self->warn('Use of this module is deprecated.  Use Bio::SearchIO::rnamotif instead');
  $self->SUPER::_initialize(@args);
  my ($motiftag,$desctag,$srctag) =  $self->SUPER::_rearrange([qw(MOTIFTAG
                                                                  DESCTAG
                                                                  SRCTAG
                                 )],
                                  @args);
  $self->motif_tag(defined $motiftag ? $motiftag : $MotifTag);
  $self->source_tag(defined $srctag ? $srctag : $SrcTag);
  $self->desc_tag(defined $desctag ? $desctag : $DescTag);
  $self->{'_sec_structure' => '',
          '_dfile' => ''};
}

=head2 motif_tag

 Title   : motif_tag
 Usage   : $obj->motif_tag($newval)
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
           Default is 'RNAMotif' as set by the global $SrcTag
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
           the tag '-display_name'.  Default is 'rnamotif' as set by the global
           $DescTag.  Use this to manually set the descriptor (motif searched for).
           Since there is no way for this module to tell what the motif is from the
           name of the descriptor file or the RNAMotif output, this should
           be set every time an RNAMotif object is instantiated for clarity
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
             /RNAMotif/i.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method { 
#-------------
    my ($self, $method) = @_;  
    if($method && ($method !~ /RNAMotif/i)) {
    $self->throw("method $method not supported in " . ref($self));
    }
    return $self->SUPER::analysis_method($method);
}

=head2 next_feature

 Title   : next_feature
 Usage   : while($gene = $obj->next_feature()) {
                  # do something
           }
 Function: Returns the next gene structure prediction of the RNAMotif result
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
 Function: Returns the next gene structure prediction of the RNAMotif result
           file. Call this method repeatedly until FALSE is returned.
 Returns : A Bio::SeqFeature::Generic object
 Args    : None (at present)

=cut

sub next_prediction {
    my ($self) = @_;
    my ($motiftag,$srctag,$desctag) = ( $self->motif_tag,
                       $self->source_tag,
                       $self->desc_tag);
    my ($score, $strand, $start, $length, $sequence, $end, $seqid, $description)=0;
    while($_ = $self->_readline) {
        while($_ =~ /^#RM/) { # header line
            if(/^#RM\sdescr\s(.*)$/) { # contains sec structure
                $self->{'_sec_structure'}=$1;
            }
            if(/^#RM\sdfile\s(.*)$/) { # contains dfile
                $self->{'_dfile'}=$1;
            }
            $_ = $self->_readline;
        }
        if(m/^>((\S*)\s.*)$/) {
            $seqid = $2;
            $description = $1; # contains entire description line if needed
            if($seqid =~  /(gb|emb|dbj|sp|pdb|bbs|ref|lcl)\|(.*)\|/) {
                $seqid = $2; # pulls out gid
            }
        }
        # start pulling out hit information...
        # regex is m/^\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s(.*)$/
        # m/^\S+\s+     # seqID, not needed
        # (.+)\s+       # score, or extra info using sprintf()
        # (\d+)\s+      # strand
        # (\d+)\s+      # start
        # (\d+)\s       # length
        # (.*)$/        # sequence, divided according to descriptor
        if(m/^\S+\s+(.+)\s+(\d+)\s+(\d+)\s+(\d+)\s(.*)$/) {
            ($score, $strand, $start, $length, $sequence, $end)=
                ($1, $2, $3, $4, $5, 0);
            if( $strand==0 ) {
                $end = $start + $length -1;
                $strand = 1;
            } else {
                $end = $start - $length + 1;
                ($start, $end, $strand) = ($end, $start, -1);
            }
            my $gene = Bio::SeqFeature::Generic->new(-seq_id => $seqid,
                                                      -start  => $start,
                                                      -end    => $end,
                                                      -strand => $strand,
                                                      -score  => $score,
                                                      -primary_tag => $motiftag,
                                                      -source_tag  => $srctag,
                                                      -display_name => $desctag,
                                                      -tag     => {
                                                        'descline'       => $description,
                                                        'descfile'      => $self->{'_dfile'},
                                                        'secstructure'  => $self->{'_sec_structure'},
                                                        'sequence'       => $sequence});
            return $gene;
        }
    }
}

=head2 clean_features

 Title   : next_prediction
 Usage   : @list = $obj->clean_features(-10);
 Function: Cleans (reduces redundant hits) based on score, position
 Returns : a Bio::SeqFeature::Collection object
 Args    : Pos./Neg. integer (for highest/lowest scoring seqfeature within x bp).
 Throws  : Error unless digit is entered.  

=cut

sub clean_features {
    my $self = shift;
    my $bp = shift;
    $self->throw("No arg, need pos. or neg. integer") if !$bp;
    $self->throw("Need pos. or neg. integer") if ($bp !~ /\-?\d/ || $bp =~ /\./);
    my ($b, $sf2);
    my @list = ();
    my @features = ();
    while (my $pred = $self->next_prediction) {
        push @features, $pred;
    }
    while (@features) {
        $b = shift @features if !defined($b);
        $sf2 = shift @features;
        # from same sequence?
        if ($sf2) { # if there is no feature, then...
            if ($b->seq_id == $sf2->seq_id) {
                # close starts (probable redundant hit)?
                if(abs(($b->start)-($sf2->start)) <= abs($bp)) {
                    # which score is better?
                    if( (($bp < 0) && ($b->score > $sf2->score))  ||  # lowest score
                        (($bp > 0) && ($b->score < $sf2->score)) ){   # highest score
                        $b = $sf2;
                        next;
                    } else {
                        next;
                    }
                }
                push @list, $b;
                $b = $sf2;
            }
        }
    }
    push @list, $b if $b;
    my $col = Bio::SeqFeature::Collection->new;
    $col->add_features(\@list);
    return $col;
}

1;
