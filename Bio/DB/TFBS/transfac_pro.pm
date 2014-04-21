# $Id: transfac_pro.pm,v 1.15 2006/08/12 11:00:03 sendu Exp $
#
# BioPerl module for Bio::DB::TFBS::transfac_pro
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright Sendu Bala
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::DB::TFBS::transfac_pro - An implementation of Bio::DB::TFBS
which uses local flat files for transfac pro

=head1 SYNOPSIS

  use Bio::DB::Taxonomy;

  my $db = new Bio::DB::Taxonomy(-source => 'transfac_pro'
                                 -dat_dir => $directory);

  # we're interested in the gene P5
  my ($gene_id) = $db->get_gene_ids(-name => 'P5'); # G000001

  # we want all the transcription factors that bind to our gene
  my @factor_ids = $db->get_factor_ids(-gene => $gene_id);

  # get info about those TFs
  foreach my $factor_id (@factor_ids) {
    my $factor = $db->get_factor($factor_id);
    my $name = $factor->universal_name;
    # etc. - see Bio::Map::TranscriptionFactor, eg. find out where it binds
  }

  # get a matrix
  my $matrix = $db->get_matrix('M00001');

  # get a binding site sequence
  my $seq = $db->get_site('R00001');

=head1 DESCRIPTION

This is an implementation which uses local flat files and the DB_File
module RECNO data structures to manage a local copy of the Transfac Pro TFBS
database.

Required database files require a license which can be obtained via
http://www.biobase-international.com/pages/index.php?id=170

Within the linux installation tarball you will find a cgibin tar ball, and
inside that is a data directory containing the .dat files needed by this
module. Point to that data directory with -dat_dir

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 CONTRIBUTORS

Based on Bio::DB::Taxonomy::flatfile by Jason Stajich

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::DB::TFBS::transfac_pro;
use strict;
use Bio::Annotation::Reference;
use Bio::Annotation::SimpleValue;
use Bio::LocatableSeq;
use Bio::SimpleAlign;
use Bio::Matrix::PSM::SiteMatrix;
use Bio::AlignIO;
use Bio::Map::GeneMap;
use Bio::Map::TranscriptionFactor;
use Bio::Map::Position;
use Bio::Map::Relative;
use DB_File;

use constant SEPARATOR => ':!:';
use constant INTERNAL_SEPARATOR => '!:!';

$DB_BTREE->{'flags'} = R_DUP; # allow duplicate values in DB_File BTREEs

use base qw(Bio::DB::TFBS);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::DB::TFBS::transfac_pro();
 Function: Builds a new Bio::DB::TFBS::transfac_pro object 
 Returns : an instance of Bio::DB::TTFBS::transfac_pro
 Args    : -dat_dir   => name of directory where Transfac Pro .dat files
                         (required to initially build indexes)
           -tax_db    => Bio::DB::Taxonomy object, used when initially building
                         indexes, gives better results for species information
                         but not required.
           -index_dir => name of directory where index files should be created
                         or already exist. (defaults to -dat_dir, required if
                         -dat_dir not supplied)
           -force     => 1 replace current indexes even if they exist

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    my ($dat_dir, $index_dir, $tax_db, $force) = $self->_rearrange([qw(DAT_DIR INDEX_DIR TAX_DB FORCE)], @args);
    $self->throw("At least one of -dat_dir and -index_dir must be supplied") unless ($dat_dir || $index_dir);
    
    $self->index_directory($index_dir || $dat_dir);
    $self->{_tax_db} = $tax_db if $tax_db;
    
    if ($dat_dir) {
        $self->_build_index($dat_dir, $force);
    }
    
    $self->_db_connect;
    return $self;
}

=head2 Bio::DB::TFBS Interface implementation

=cut

sub _get_ids {
    my ($self, $dat, @args) = @_;
    @args % 2 == 0 || $self->throw("Must provide key => value pairs");
    my $hash = $self->{$dat} || $self->throw("Unknown .dat type '$dat'");
    
    if (@args) {
        # get a subset corresponding to args
        my @final;
        my %args = @args;
        my $multiple = 0;
        while (my ($type, $value) = each %args) {
            unless ($value) {
                $self->warn("Arguement '$type' has no value, ignored");
                next;
            }
            $type =~ s/-//;
            $type = lc($type);
            my $converter = $hash->{$type};
            unless ($converter) {
                $self->warn("Unknown search type '$type' for .dat type '$dat'");
                next;
            }
            
            my @ids = $converter->get_dup($value);
            unless (@ids) {
                @ids = $converter->get_dup(lc($value));
            }
            
            if ($multiple) {
                # we can have multiple types given at once, find the ids that
                # satisfy all criteria
                @final || return;
                my %final = map { $_ => 1 } @final;
                @final = grep { $final{$_} } @ids;
            }
            else {
                @final = @ids;
                $multiple++;
            }
        }
        
        return @final;
    }
    else {
        # get them all
        my $db_file_hash = $self->{$dat}->{id};
        
        my ($key, $prev_key, $value) = ('_!_', '!_!');
        my @ids;
        while (1) {
            $db_file_hash->seq($key, $value, R_NEXT);
            last if $prev_key eq $key;
            push(@ids, $value); # confusing? when creating objects we store
                                # $value as accession and $key as id, but from
                                # this method we return $value as id given $id!
            $prev_key = $key;
        }
        
        return @ids;
    }
}

=head2 get_reference

 Title   : get_reference
 Usage   : my $ref = $obj->get_reference($id);
 Function: Get a literature reference.
 Returns : Bio::Annotation::Reference
 Args    : string - a reference id ('RE...')

=cut

sub get_reference {
    my ($self, $id) = @_;
    $id || return;
    my $data = $self->{reference}->{data}->{$id} || return;
    my @data = split(SEPARATOR, $data);
    
    return Bio::Annotation::Reference->new(-pubmed   => $data[0],
                                           -authors  => $data[1],
                                           -title    => $data[2],
                                           -location => $data[3] );
}

=head2 get_genemap

 Title   : get_genemap
 Usage   : my $map = $obj->get_genemap($id);
 Function: Get a GeneMap for a gene.
 Returns : Bio::Map::GeneMap
 Args    : string - a gene id ('G...'), and optionally int (number of bp
           upstream)

=cut

sub get_genemap {
    my ($self, $id, $upstream) = @_;
    $id || return;
    return $self->{got_map}->{$id} if defined $self->{got_map}->{$id};
    $upstream ||= 1000;
    my $data = $self->{gene}->{data}->{$id} || return;
    my @data = split(SEPARATOR, $data);
    
    # accession = id name description species_tax_id_or_raw_string
    my $taxon = $self->{_tax_db} ? $self->{_tax_db}->get_taxon($data[3]) || $data[3] : $data[3];
    my $map = Bio::Map::GeneMap->get(-uid => $id,
                                     -gene => $data[1],
                                     -species => $taxon,
                                     -description => $data[2],
                                     -upstream => $upstream);
    $self->{got_map}->{$id} = $map; # prevents infinite recurse when we call get_factor below
    
    # spawn all the factors that belong on this gene map
    # get_factor_ids(-gene => ...) only works for genes that encode factors;
    # have to go via sites
    foreach my $sid ($self->get_site_ids(-gene => $id)) {
        foreach my $fid ($self->get_factor_ids(-site => $sid)) {
            # it is quite deliberate that we deeply recurse to arrive at the
            # correct answer, which involves pulling in most of the database
            no warnings "recursion";
            $self->get_factor($fid);
        }
    }
    
    return $map;
}

=head2 get_seq

 Title   : get_seq
 Usage   : my $seq = $obj->get_seq($id);
 Function: Get the sequence of a site. The sequence will be annotated with the
           the tags 'relative_start', 'relative_end', 'relative_type' and
           'relative_to'.
 Returns : Bio::Seq
 Args    : string - a site id ('R...')

=cut

sub get_seq {
    my ($self, $id) = @_;
    $id || return;
    my $data = $self->{site}->{data}->{$id} || return;
    my @data = split(SEPARATOR, $data);
    
    my $seq = Bio::Seq->new(-seq              => $data[2],
                            -accession_number => $id,
                            -description      => $data[6] ? 'Genomic sequence' : 'Consensus or artificial sequence',
                            -id               => $data[0],
                            -strand           => 1,
                            -alphabet         => $data[7] || 'dna',
                            -species          => $data[6]);
    
    my $annot = $seq->annotation;
    my $sv = Bio::Annotation::SimpleValue->new(-tagname => 'relative_start', -value => $data[4] || 1);
    $annot->add_Annotation($sv);
    $sv = Bio::Annotation::SimpleValue->new(-tagname => 'relative_end', -value => $data[5] || ($data[4] || 1 + length($data[2]) - 1));
    $annot->add_Annotation($sv);
    $sv = Bio::Annotation::SimpleValue->new(-tagname => 'relative_type', -value => $data[3] || 'artificial');
    $annot->add_Annotation($sv);
    $sv = Bio::Annotation::SimpleValue->new(-tagname => 'relative_to', -value => $data[1]);
    $annot->add_Annotation($sv);
    
    return $seq;
}

=head2 get_fragment

 Title   : get_fragment
 Usage   : my $seq = $obj->get_fragment($id);
 Function: Get the sequence of a fragment.
 Returns : Bio::Seq
 Args    : string - a site id ('FR...')

=cut

sub get_fragment {
    my ($self, $id) = @_;
    $id || return;
    my $data = $self->{fragment}->{data}->{$id} || return;
    my @data = split(SEPARATOR, $data);
    
    # accession = id gene_id1 gene_id2 species_tax_id_or_raw_string sequence source
    return new Bio::Seq( -seq              => $data[4],
                         -accession_number => $id,
                         -description      => 'Between genes '.$data[1].' and '.$data[2],
                         -species          => $data[3],
                         -id               => $data[0],
                         -alphabet         => 'dna' );
}

=head2 get_matrix

 Title   : get_matrix
 Usage   : my $matrix = $obj->get_matrix($id);
 Function: Get a matrix that describes a binding site.
 Returns : Bio::Matrix::PSM::SiteMatrix
 Args    : string - a matrix id ('M...'), optionally a sequence string from
           which base frequencies will be calcualted for the matrix model
           (default 0.25 each)

=cut

sub get_matrix {
    my ($self, $id, $seq) = @_;
    $id || return;
    $seq ||= 'atgc';
    $seq = lc($seq);
    my $data = $self->{matrix}->{data}->{$id} || return;
    my @data = split(SEPARATOR, $data);
    $data[4] || $self->throw("Matrix data missing for $id");
    
    my ($a, $c, $g, $t);
    foreach my $position (split(INTERNAL_SEPARATOR, $data[4])) {
        my ($a_count, $c_count, $g_count, $t_count) = split("\t", $position);
        push(@{$a}, $a_count);
        push(@{$c}, $c_count);
        push(@{$g}, $g_count);
        push(@{$t}, $t_count);
    }
    
    # our psms include a simple background model so we can use
    # sequence_match_weight() if desired
    my $a_freq = ($seq =~ tr/a//) / length($seq);
    my $c_freq = ($seq =~ tr/c//) / length($seq);
    my $g_freq = ($seq =~ tr/g//) / length($seq);
    my $t_freq = ($seq =~ tr/t//) / length($seq);
    
    my $psm = Bio::Matrix::PSM::SiteMatrix->new(-pA => $a,
                                                -pC => $c,
                                                -pG => $g,
                                                -pT => $t,
                                                -id => $data[0],
                                                -accession_number => $id,
                                                -sites => $data[3],
                                                -width => scalar(@{$a}),
                                                -correction => 1,
                                                -model => { A => $a_freq, C => $c_freq, G => $g_freq, T => $t_freq } );
    
    #*** used to make a Bio::Matrix::PSM::Psm and add references, but it
    #    didn't seem worth it. You can get references from the database by:
    #foreach my $ref_id ($db->get_reference_ids(-matrix => $id)) {
    #    my $ref = $db->get_reference($ref_id);
    #}
    
    return $psm;
}

=head2 get_aln

 Title   : get_aln
 Usage   : my $aln = $obj->get_aln($id);
 Function: Get the alignment that was used to generate a matrix. Each sequence
           in the alignment will have an accession_number corresponding to the
           Transfac site id, and id() based on that but unique within the
           alignment.
 Returns : Bio::SimpleAlign
 Args    : string - a matrix id ('M...'), optionally true to, when a matrix
           lists no sequences, search for sequences via the matrix's factors,
           picking the sites that best match the matrix

=cut

my %VALID_STRAND = map {$_ => 1} qw(-1 0 1);

sub get_aln {
    my ($self, $id, $via_factors) = @_;
    $id || return;
    my $data = $self->{matrix}->{data}->{$id} || $self->throw("matrix '$id' had no data in DB_File");
    my @data = split(SEPARATOR, $data);

    if (! $data[5] && $via_factors) {
        # This is a matrix with no site sequences given in matrix.dat.
        # Find some matching site sequences via factors.
        
        # First, check its factors for sites
        my %site_seqs;
        my %factor_ids;
        foreach my $factor_id ($self->get_factor_ids(-matrix => $id)) {
            $factor_ids{$factor_id} = 1;
            foreach my $site_id ($self->get_site_ids(-factor => $factor_id)) {
                next if defined $site_seqs{$site_id};
                my $seq = $self->get_seq($site_id);
                
                # skip sites that have no sequence, or have IUPAC symbols in
                # their sequence (most probably the 'consensus' sequence itself
                # that was used to make and exactly corresponds to the matrix)
                my $seq_str = $seq->seq || next;
                $seq_str =~ /[MRWSYKVHDB]/ and next;
                
                $site_seqs{$site_id} = $seq;
            }
        }
        my @seqs = values %site_seqs;
        
        if (@seqs > 1) {
            # pick the sub-seqs that match to the matrix with the best scores
            my $matrix = $self->get_matrix($id);
            my $desired_sequences = $matrix->sites;
            return if @seqs < $desired_sequences;
            
            my $desired_length = $matrix->width;
            my %best_seqs;
            foreach my $seq (@seqs) {
                my $for_str = $seq->seq;
                next if length($for_str) < $desired_length;
                my $rev_str = $seq->revcom->seq;
                
                my $best_score = 0;
                my $best_subseq = '';
                my $best_i = 0;
                my $best_subseq_caps = 0;
                my $best_revcom;
                my $revcom = 0;
                foreach my $seq_str ($for_str, $rev_str) {
                    for my $i (0..(length($seq_str) - $desired_length)) {
                        my $subseq = substr($seq_str, $i, $desired_length);
                        $subseq =~ s/[^ACGTacgt]//g; # can only score atcg
                        next unless length($subseq) == $desired_length; # short or 0-length seqs could get the highest scores!
                        my $score = $matrix->sequence_match_weight($subseq);
                        
                        # caps represent the author-chosen bit of a site
                        # sequence so we would prefer to choose a subseq that
                        # contains it
                        my $caps = $subseq =~ tr/ACGT//;
                        
                        #*** (don't know why numeric == fails for comparing
                        #     scores, when the string eq works)
                        if ($score > $best_score || ("$score" eq "$best_score" && $caps > $best_subseq_caps)) {
                            $best_score = $score;
                            $best_subseq_caps = $caps;
                            $best_subseq = $subseq;
                            $best_i = $i;
                            $best_revcom = $revcom;
                        }
                    }
                    $revcom++;
                }
                
                if ($best_score) {
                    $best_seqs{$seq->accession_number} = [$best_subseq, $seq->accession_number, ($best_i + 1), $revcom ? -1 : 1, $best_score];
                }
            }
            my @sorted = sort { $best_seqs{$b}->[-1] <=> $best_seqs{$a}->[-1] } keys %best_seqs;
            return if @sorted < $desired_sequences;
            splice(@sorted, $desired_sequences);
            my %wanted = map { $_ => 1 } @sorted;
            
            my @site_data;
            foreach my $seq (@seqs) {
                next unless exists $wanted{$seq->accession_number};
                my @data = @{$best_seqs{$seq->accession_number}};
                pop(@data);
                push(@site_data, join('_', @data));
            }
            
            $data[5] = join(INTERNAL_SEPARATOR, @site_data);
            $self->{matrix}->{data}->{$id} = join(SEPARATOR, @data);
        }
    }
    $data[5] || return;
    
    my @blocks = split(INTERNAL_SEPARATOR, $data[5]);
    
    # append gap chars to all sequences to make them the same length
    # (applies to sequences found via factors, presumably, since we already
    # do this for matrix alignments in transfac_pro.pm)
    my $longest = 0;
    foreach (@blocks) {
        my ($seq) = split('_', $_);
        my $length = length($seq);
        if ($length > $longest) {
            $longest = $length;
        }
    }
    foreach my $i (0..$#blocks) {
        my $block = $blocks[$i];
        my ($seq, $seq_id) = split('_', $block);
        my $length = length($seq);
        if ($length < $longest) {
            my $orig_seq = $seq;
            $seq .= '-'x($longest - $length);
            $block =~ s/^${orig_seq}_/${seq}_/;
            $blocks[$i] = $block;
        }
    }
    
    # build the alignment
    my $aln = Bio::SimpleAlign->new(-source => 'transfac_pro');
    my %done_ids;
    foreach (@blocks) {
        my ($seq, $seq_acc, $start, $strand) = split('_', $_);
        
        $self->throw("Invalid strand $strand found in block $_")
            unless exists $VALID_STRAND{$strand};
        # we can get back multiple different subparts of the same site (sequence),
        # so $seq_acc isn't unique across this loop. Can't use it as the seq id
        # of the alignment (ids must be unique in SimpleAlign), so we
        # uniquify the id and store the original id as the accession_number
        my $seq_id;
        $done_ids{$seq_acc}++;
        if ($done_ids{$seq_acc} > 1) {
            $seq_id = $seq_acc.'_'.$done_ids{$seq_acc};
        }
        else {
            $seq_id = $seq_acc;
        }
        
        my $gaps = $seq =~ tr/-//;
        my $length = length($seq) - $gaps;
        $self->throw("seq '$seq_id' for matrix '$id' had seq '$seq'") unless $length;
        $aln->add_seq(Bio::LocatableSeq->new(-seq    => $seq,
                                             -id     => $seq_id,
                                             -accession_number => $seq_acc,
                                             -start  => $start,
                                             -end    => $start + $length - 1,
                                             -strand => $strand));
    }
    $aln->id($id);
    # could also store score? of?
    
    return $aln;
}

=head2 get_factor

 Title   : get_factor
 Usage   : my $factor = $obj->get_factor($id);
 Function: Get the details of a transcription factor.
 Returns : Bio::Map::TranscriptionFactor
 Args    : string - a factor id ('T...')

=cut

sub get_factor {
    my ($self, $id) = @_;
    $id || return;
    return $self->{got_factor}->{$id} if defined $self->{got_factor}->{$id};
    my $data = $self->{factor}->{data}->{$id} || return;
    my @data = split(SEPARATOR, $data);
    
    # accession = id name species sequence
    my $tf = Bio::Map::TranscriptionFactor->get(-id => $id,
                                                -universal_name => $data[1]);
    #*** not sure what to do with species and sequence, since we don't want to
    # confuse the idea that a TF is a general thing that could bind to any
    # species... then again, you might want to model species-specific variants
    # of a TF with different binding abilities...
    #*** idea of having inclusion and exclusion species so you can prevent/
    # ignore a tf that binds to the wrong species (a species that doesn't even
    # have the tf), and associating sequence with each species/tf combo so you
    # can see how diverged the tf is and make assumptions about site difference
    # allowance
    
    # place it on all its genemaps
    foreach my $sid ($self->get_site_ids(-factor => $id)) {
        my $s_data = $self->{site}->{data}->{$sid} || next;
        my @s_data = split(SEPARATOR, $s_data);
        
        # accession = id gene_id sequence relative_to first_position last_position species_tax_id_or_raw_string
        $s_data[1] || next; # site isn't relative to a gene, meaningless
        $s_data[4] || next; # don't know where its supposed to be, can't model it
        $s_data[5] ||= $s_data[4] + ($s_data[2] ? length($s_data[2]) - 1 : 0);
        
        # it is quite deliberate that we deeply recurse to arrive at the
        # correct answer, which involves pulling in most of the database
        no warnings "recursion";
        my $gene_map = $self->get_genemap($s_data[1]) || next;
        return $self->{got_factor}->{$id} if defined $self->{got_factor}->{$id};
        
        #*** not always relative to gene start...
        #    we need Bio::Map::Gene s to have some default tss and atg positions
        #    that we can be relative to
        my $rel = Bio::Map::Relative->new(-element => $gene_map->gene, -description => $s_data[3]);
        Bio::Map::Position->new(-map => $gene_map, -element => $tf, -start => $s_data[4], -end => $s_data[5], -relative => $rel);
    }
    
    $self->{got_factor}->{$id} = $tf;
    return $tf;
}

# since get_factor() is uncertain, just have direct access methods to factor
# information
sub get_factor_name {
    my ($self, $id) = @_;
    my $details = $self->_get_factor_details($id) || return;
    return $details->{name};
}
sub get_factor_species {
    my ($self, $id) = @_;
    my $details = $self->_get_factor_details($id) || return;
    return $details->{species};
}
sub get_factor_sequence {
    my ($self, $id) = @_;
    my $details = $self->_get_factor_details($id) || return;
    return $details->{sequence};
}
sub _get_factor_details {
    my ($self, $id) = @_;
    $id || return;
    
    return $self->{factor_details}->{$id} if defined $self->{factor_details}->{$id};
    
    my $data = $self->{factor}->{data}->{$id} || return;
    my @data = split(SEPARATOR, $data);
    
    # accession = id name species sequence
    
    my %details = (name => $data[1], species => $data[2], sequence => $data[3]);
    $self->{factor_details}->{$id} = \%details;
    
    return \%details;
}

=head2 get_reference_ids

 Title   : get_reference_ids
 Usage   : my @ids = $obj->get_reference_ids(-key => $value);
 Function: Get all the reference ids that are associated with the supplied
           args.
 Returns : list of strings (ids)
 Args    : -key => value, where value is a string id, and key is one of:
           -pubmed -site -gene -matrix -factor

=cut

sub get_reference_ids {
    my $self = shift;
    return $self->_get_ids('reference', @_);
}

# -id -name -species -site -factor -reference
sub get_gene_ids {
    my $self = shift;
    return $self->_get_ids('gene', @_);
}

=head2 get_site_ids

 Title   : get_site_ids
 Usage   : my @ids = $obj->get_site_ids(-key => $value);
 Function: Get all the site ids that are associated with the supplied
           args.
 Returns : list of strings (ids)
 Args    : -key => value, where value is a string id, and key is one of:
           -id -species -gene -matrix -factor -reference

=cut

sub get_site_ids {
    my $self = shift;
    return $self->_get_ids('site', @_);
}

=head2 get_matrix_ids

 Title   : get_matrix_ids
 Usage   : my @ids = $obj->get_matrix_ids(-key => $value);
 Function: Get all the matrix ids that are associated with the supplied
           args.
 Returns : list of strings (ids)
 Args    : -key => value, where value is a string id, and key is one of:
           -id -name -site -factor -reference

=cut

sub get_matrix_ids {
    my $self = shift;
    return $self->_get_ids('matrix', @_);
}

=head2 get_factor_ids

 Title   : get_factor_ids
 Usage   : my @ids = $obj->get_factor_ids(-key => $value);
 Function: Get all the factor ids that are associated with the supplied
           args.
 Returns : list of strings (ids)
 Args    : -key => value, where value is a string id, and key is one of:
           -id -name -species -interactors -gene -matrix -site -reference
           NB: -gene only gets factor ids for genes that encode factors

=cut

sub get_factor_ids {
    my $self = shift;
    return $self->_get_ids('factor', @_);
}

=head2 get_fragment_ids

 Title   : get_fragment_ids
 Usage   : my @ids = $obj->get_fragment_ids(-key => $value);
 Function: Get all the fragment ids that are associated with the supplied
           args.
 Returns : list of strings (ids)
 Args    : -key => value, where value is a string id, and key is one of:
           -id -species -gene -factor -reference

=cut

sub get_fragment_ids {
    my $self = shift;
    return $self->_get_ids('fragment', @_);
}

=head2 Helper methods 

=cut

# internal method which does the indexing
sub _build_index {
    my ($self, $dat_dir, $force) = @_;
    
    # MLDBM would give us transparent complex data structures with DB_File,
    # allowing just one index file, but its yet another requirment and we
    # don't strictly need it
    
    my $index_dir = $self->index_directory;
    my $gene_index      = "$index_dir/gene.dat.index";
    my $reference_index = "$index_dir/reference.dat.index";
    my $matrix_index    = "$index_dir/matrix.dat.index";
    my $factor_index    = "$index_dir/factor.dat.index";
    my $fragment_index  = "$index_dir/fragment.dat.index";
    my $site_index      = "$index_dir/site.dat.index";
    
    my $reference_dat = "$dat_dir/reference.dat";
    if (! -e $reference_index || $force) {
        open my $REF, '<', $reference_dat or $self->throw("Could not read reference file '$reference_dat': $!");
        
        my %references;
        unlink $reference_index;
        my $ref = tie(%references, 'DB_File', $reference_index, O_RDWR|O_CREAT, 0644, $DB_HASH)
            or $self->throw("CCould not open file '$reference_index': $!");
        
        my %pubmed;
        my $reference_pubmed = $reference_index.'.pubmed';
        unlink $reference_pubmed;
        my $pub = tie(%pubmed, 'DB_File', $reference_pubmed, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$reference_pubmed': $!");
        
        my %gene;
        my $reference_gene = $gene_index.'.reference';
        unlink $reference_gene;
        my $gene = tie(%gene, 'DB_File', $reference_gene, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$reference_gene': $!");
        
        my %site;
        my $reference_site = $site_index.'.reference';
        unlink $reference_site;
        my $site = tie(%site, 'DB_File', $reference_site, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$reference_site': $!");
        
        my %fragment;
        my $reference_fragment = $fragment_index.'.reference';
        unlink $reference_fragment;
        my $fragment = tie(%fragment, 'DB_File', $reference_fragment, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$reference_fragment': $!");
        
        my %factor;
        my $reference_factor = $factor_index.'.reference';
        unlink $reference_factor;
        my $factor = tie(%factor, 'DB_File', $reference_factor, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$reference_factor': $!");
        
        my %matrix;
        my $reference_matrix = $matrix_index.'.reference';
        unlink $reference_matrix;
        my $matrix = tie(%matrix, 'DB_File', $reference_matrix, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$reference_matrix': $!");
        
        # skip the first three header lines
        <$REF>; <$REF>; <$REF>;
        
        my @data;
        while (<$REF>) {
            if (/^AC  (\S+)/) {
                $data[0] = $1;
            }
            elsif (/^RX  PUBMED: (\d+)/) {
                $data[1] = $1;
                $pub->put("$1", $data[0]);
            }
            elsif (/^RA  (.+)\n$/) {
                $data[2] = $1;
            }
            elsif (/^RT  (.+?)\.?\n$/) {
                $data[3] = $1;
            }
            elsif (/^RL  (.+?)\.?\n$/) {
                $data[4] = $1;
            }
            elsif (/^GE  TRANSFAC: (\w\d+)/) {
                $gene->put($data[0], "$1");
            }
            elsif (/^BS  TRANSFAC: (\w\d+)/) {
                $site->put($data[0], "$1");
            }
            elsif (/^FA  TRANSFAC: (\w\d+)/) {
                $factor->put($data[0], "$1");
            }
            elsif (/^FR  TRANSFAC: (FR\d+)/) {
                $fragment->put($data[0], "$1");
            }
            elsif (/^MX  TRANSFAC: (\w\d+)/) {
                $matrix->put($data[0], "$1");
            }
            elsif (/^\/\//) {
                # end of a record, store previous data and reset
                
                # accession = pubmed authors title location
                $references{$data[0]} = join(SEPARATOR, ($data[1] || '',
                                                         $data[2] || '',
                                                         $data[3] || '',
                                                         $data[4] || ''));
                
                @data = ();
            }
        }
        close $REF;
        
        $ref = $pub = $gene = $site = $fragment = $factor = $matrix = undef;
        untie %references;
        untie %pubmed;
        untie %gene;
        untie %site;
        untie %fragment;
        untie %factor;
        untie %matrix;
    }
    
    my $gene_dat = "$dat_dir/gene.dat";
    if (! -e $gene_index || $force) {
        open my $GEN, '<', $gene_dat or $self->throw("Could not read gene file '$gene_dat': $!");
        
        my %genes;
        unlink $gene_index;
        my $gene = tie(%genes, 'DB_File', $gene_index, O_RDWR|O_CREAT, 0644, $DB_HASH)
            or $self->throw("Could not open file '$gene_index': $!");
        
        my %id;
        my $gene_id = $gene_index.'.id';
        unlink $gene_id;
        my $id = tie(%id, 'DB_File', $gene_id, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$gene_id': $!");
        
        my %name;
        my $gene_name = $gene_index.'.name';
        unlink $gene_name;
        my $name = tie(%name, 'DB_File', $gene_name, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$gene_name': $!");
        
        my %species;
        my $gene_species = $gene_index.'.species';
        unlink $gene_species;
        my $species = tie(%species, 'DB_File', $gene_species, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$gene_species': $!");
        
        my %site;
        my $gene_site = $site_index.'.gene';
        unlink $gene_site;
        my $site = tie(%site, 'DB_File', $gene_site, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$gene_site': $!");
        
        my %factor;
        my $gene_factor = $factor_index.'.gene';
        unlink $gene_factor;
        my $factor = tie(%factor, 'DB_File', $gene_factor, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$gene_factor': $!");
        
        my %fragment;
        my $gene_fragment = $fragment_index.'.gene';
        unlink $gene_fragment;
        my $fragment = tie(%fragment, 'DB_File', $gene_fragment, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$gene_fragment': $!");
        
        my %reference;
        my $gene_reference = $reference_index.'.gene';
        unlink $gene_reference;
        my $reference = tie(%reference, 'DB_File', $gene_reference, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$gene_reference': $!");
        
        # skip the first three header lines
        <$GEN>; <$GEN>; <$GEN>;
        
        my @data;
        while (<$GEN>) {
            if (/^AC  (\S+)/) {
                $data[0] = $1;
            }
            elsif (/^ID  (\S+)/) {
                $data[1] = $1;
                $id->put("$1", $data[0]);
            }
            elsif (/^SD  (.+)$/) {
                $data[2] = lc("$1");
                $name->put(lc("$1"), $data[0]);
            }
            elsif (/^SY  (.+)\.$/) {
                foreach (split('; ', lc("$1"))) {
                    $name->put($_, $data[0]);
                }
            }
            elsif (/^DE  (.+)$/) {
                $data[3] = $1;
            }
            elsif (/^OS  (.+)$/) {
                my $raw_species = $1;
                my $taxid = $self->_species_to_taxid($raw_species);
                $data[4] = $taxid || $raw_species;
                $species->put($data[4], $data[0]);
            }
            elsif (/^RN  .+?(RE\d+)/) {
                $reference->put($data[0], "$1");
            }
            elsif (/^BS  .+?(R\d+)/) {
                $site->put($data[0], "$1");
            }
            elsif (/^FA  (T\d+)/) {
                $factor->put($data[0], "$1");
            }
            elsif (/^BR  (FR\d+)/) {
                $fragment->put($data[0], "$1");
            }
            elsif (/^\/\//) {
                # end of a record, store previous data and reset
                
                # accession = id name description species_tax_id_or_raw_string
                $genes{$data[0]} = join(SEPARATOR, ($data[1] || '',
                                                    $data[2] || '',
                                                    $data[3] || '',
                                                    $data[4] || ''));
                
                @data = ();
            }
        }
        close $GEN;
        
        $gene = $id = $name = $species = $site = $factor = $reference = undef;
        untie %genes;
        untie %id;
        untie %name;
        untie %species;
        untie %site;
        untie %factor;
        untie %reference;
    }
    
    my $site_dat = "$dat_dir/site.dat";
    if (! -e $site_index || $force) {
        open my $SIT, '<', $site_dat or $self->throw("Could not read site file '$site_dat': $!");
        
        my %sites;
        unlink $site_index;
        my $site = tie(%sites, 'DB_File', $site_index, O_RDWR|O_CREAT, 0644, $DB_HASH)
            or $self->throw("Could not open file '$site_index': $!");
        
        my %id;
        my $site_id = $site_index.'.id';
        unlink $site_id;
        my $id = tie(%id, 'DB_File', $site_id, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$site_id': $!");
        
        my %species;
        my $site_species = $site_index.'.species';
        unlink $site_species;
        my $species = tie(%species, 'DB_File', $site_species, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$site_species': $!");
        
        my %qualities;
        my $site_qualities = $site_index.'.qual';
        unlink $site_qualities;
        my $quality = tie(%qualities, 'DB_File', $site_qualities, O_RDWR|O_CREAT, 0644, $DB_HASH)
            or $self->throw("Could not open file '$site_qualities': $!");
        
        my %gene;
        my $site_gene = $gene_index.'.site';
        unlink $site_gene;
        my $gene = tie(%gene, 'DB_File', $site_gene, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$site_gene': $!");
        
        my %matrix;
        my $site_matrix = $matrix_index.'.site';
        unlink $site_matrix;
        my $matrix = tie(%matrix, 'DB_File', $site_matrix, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$site_matrix': $!");
        
        my %factor;
        my $site_factor = $factor_index.'.site';
        unlink $site_factor;
        my $factor = tie(%factor, 'DB_File', $site_factor, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$site_factor': $!");
        
        my %reference;
        my $site_reference = $reference_index.'.site';
        unlink $site_reference;
        my $reference = tie(%reference, 'DB_File', $site_reference, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$site_reference': $!");
        
        # skip the first three header lines
        <$SIT>; <$SIT>; <$SIT>;
        
        my @data;
        while (<$SIT>) {
            if (/^AC  (\S+)/) {
                $data[0] = $1;
            }
            elsif (/^ID  (\S+)/) {
                $data[1] = $1;
                $id->put("$1", $data[0]);
            }
            elsif (/^TY  (.+)$/) {
                $data[8] = $1;
            }
            elsif (/^DE  .*Gene: (G\d+)/) {
                $data[2] = $1;
                $gene->put($data[0], "$1");
                
                # if it has no gene it is an artificial sequence, unless it
                # has a species (OS line), in which case it is unassigned
                # genomic; either way we won't be able to make a
                # Bio::Map::PositionI later on, so such sites won't be
                # on any MapI.
            }
            elsif (/^OS  (.+)$/) {
                # Since not all sites in site.dat with a species have a gene,
                # (small handful are unassigned 'genomic') can't delegate to
                # gene.dat and must parse species here (effectively again)
                my $raw_species = $1;
                my $taxid = $self->_species_to_taxid($raw_species);
                $data[7] = $taxid || $raw_species;
                $species->put($data[7], $data[0]);
            }
            elsif (/^SQ  (.+)\.$/) {
                $data[3] = $1;
                # there can actually be more than one SQ line, seemingly with
                # variations of the sequence (not a long sequence split over
                # two lines); not sure what to do with data; currently we end
                # up storing only the last variant.
            }
            elsif (/^S1  (.+)$/) {
                $data[4] = $1;
                # if S1 not present, means transcriptional start
            }
            elsif (/^SF  (.+)$/) {
                $data[5] = $1;
            }
            elsif (/^ST  (.+)$/) {
                $data[6] = $1;
            }
            elsif (/^RN  .+?(RE\d+)/) {
                $reference->put($data[0], "$1");
            }
            elsif (/^MX  (M\d+)/) {
                $matrix->put($data[0], "$1");
            }
            elsif (/^BF  (T\d+); .+?; Quality: (\d)/) {
                $factor->put($data[0], "$1");
                $qualities{$data[0].SEPARATOR.$1} = $2;
            }
            elsif (/^\/\//) {
                # end of a record, store previous data and reset
                
                # accession = id gene_id sequence relative_to first_position last_position species_tax_id_or_raw_string type
                $sites{$data[0]} = join(SEPARATOR, ($data[1] || '',
                                                    $data[2] || '',
                                                    $data[3] || '',
                                                    $data[4] || 'TSS',
                                                    $data[5] || '',
                                                    $data[6] || '',
                                                    $data[7] || '',
                                                    $data[8] || ''));
                
                @data = ();
            }
        }
        close $SIT;
        
        $site = $id = $species = $quality = $gene = $matrix = $factor = $reference = undef;
        untie %sites;
        untie %id;
        untie %species;
        untie %qualities;
        untie %gene;
        untie %matrix;
        untie %factor;
        untie %reference;
    }
    
    my $matrix_dat = "$dat_dir/matrix.dat";
    if (! -e $matrix_index || $force) {
        open my $MAT, '<', $matrix_dat or $self->throw("Could not read matrix file '$matrix_dat': $!");
        
        my %matrices;
        unlink $matrix_index;
        my $matrix = tie(%matrices, 'DB_File', $matrix_index, O_RDWR|O_CREAT, 0644, $DB_HASH)
            or $self->throw("Could not open file '$matrix_index': $!");
        
        my %id;
        my $matrix_id = $matrix_index.'.id';
        unlink $matrix_id;
        my $id = tie(%id, 'DB_File', $matrix_id, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$matrix_id': $!");
        
        my %name;
        my $matrix_name = $matrix_index.'.name';
        unlink $matrix_name;
        my $name = tie(%name, 'DB_File', $matrix_name, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$matrix_name': $!");
        
        my %site;
        my $matrix_site = $site_index.'.matrix';
        unlink $matrix_site;
        my $site = tie(%site, 'DB_File', $matrix_site, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$matrix_site': $!");
        
        my %factor;
        my $matrix_factor = $factor_index.'.matrix';
        unlink $matrix_factor;
        my $factor = tie(%factor, 'DB_File', $matrix_factor, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$matrix_factor': $!");
        
        my %reference;
        my $matrix_reference = $reference_index.'.matrix';
        unlink $matrix_reference;
        my $reference = tie(%reference, 'DB_File', $matrix_reference, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$matrix_reference': $!");
        
        # skip the first three header lines
        <$MAT>; <$MAT>; <$MAT>;
        
        my @data;
        my @matrix_data;
        my @site_data;
        while (<$MAT>) {
            if (/^AC  (\S+)/) {
                $data[0] = $1;
            }
            elsif (/^ID  (\S+)/) {
                $data[1] = $1;
                $id->put("$1", $data[0]);
            }
            elsif (/^NA  (.+)$/) {
                $data[2] = $1;
                $name->put("$1", $data[0]);
            }
            elsif (/^DE  (.+)$/) {
                $data[3] = $1;
            }
            elsif (/^\d\d  \s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
                # a, c, g, t counts/weights
                push(@matrix_data, join("\t", ($1, $2, $3, $4)));
                
                # Work out the number of sites as the largest number of
                # sites amongst all positions in the sequences. (The BA
                # line isn't reliable for telling us the correct number of
                # sites all the time)
                my $num = $1 + $2 + $3 + $4;
                $data[4] ||= 0;
                if ($num > $data[4]) {
                    $data[4] = $num;
                }
            }
            elsif (/^BS  ([\sa-zA-Z]+); (.+?); (-?\d+); \d+;.*; ([np])/) {
                # sequence id start strand
                push(@site_data, join('_', ($1, $2, $3, $4 eq 'p' ? 1 : -1)));
                $site->put($data[0], $2);
            }
            elsif (/^BF  (T\d+)/) {
                $factor->put($data[0], "$1");
            }
            elsif (/^RN  .+?(RE\d+)/) {
                $reference->put($data[0], "$1");
            }
            elsif (/^\/\//) {
                # end of a record, store previous data and reset
                my $matrix_data = join(INTERNAL_SEPARATOR, @matrix_data) || '';
                
                # sites of a matrix are pre-aligned but padded with spaces on
                # the left and no padding on the right; pad with -s both sides
                my $longest_seq = 0;
                
                # For all the work, does anything meaningful actually get passed
                # on here? Commenting out fixes the latest crashes on trunk.
                # 5-10-10 cjfields
                
                #foreach my $site_seq (map {my ($seq) = split("_", $_ ,2); $seq;} @site_data) {
                #    $site_seq =~ s/ /-/g;
                #    my $length = length($site_seq);
                #    if ($length > $longest_seq) {
                #        $longest_seq = $length;
                #    }
                #}
                #foreach my $site (@site_data) {
                #    my ($site_seq) = split("_", $site ,2);
                #    my $length = length($site_seq);
                #    if ($length < $longest_seq) {
                #        $site_seq .= '-' x ($longest_seq - $length);
                #    }
                #}

                my $site_data = join(INTERNAL_SEPARATOR, @site_data) || '';
                
                # accession = id name description num_of_sites matrix_data site_data
                $matrices{$data[0]} = join(SEPARATOR, ($data[1] || '',
                                                       $data[2] || '',
                                                       $data[3] || '',
                                                       $data[4],
                                                       $matrix_data,
                                                       $site_data));
                
                @data = @matrix_data = @site_data = ();
            }
        }
        close $MAT;
        
        $matrix = $id = $name = $site = $factor = $reference = undef;
        untie %matrices;
        untie %id;
        untie %name;
        untie %site;
        untie %factor;
        untie %reference;
    }
    
    my $factor_dat = "$dat_dir/factor.dat";
    if (! -e $factor_index || $force) {
        open my $FAC, '<', $factor_dat or $self->throw("Could not read factor file '$factor_dat': $!");
        
        my %factors;
        unlink $factor_index;
        my $factor = tie(%factors, 'DB_File', $factor_index, O_RDWR|O_CREAT, 0644, $DB_HASH)
            or $self->throw("Could not open file '$factor_index': $!");
        
        my %id;
        my $factor_id = $factor_index.'.id';
        unlink $factor_id;
        my $id = tie(%id, 'DB_File', $factor_id, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$factor_id': $!");
        
        my %name;
        my $factor_name = $factor_index.'.name';
        unlink $factor_name;
        my $name = tie(%name, 'DB_File', $factor_name, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$factor_name': $!");
        
        my %species;
        my $factor_species = $factor_index.'.species';
        unlink $factor_species;
        my $species = tie(%species, 'DB_File', $factor_species, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$factor_species': $!");
        
        my %interactors;
        my $factor_interactors = $factor_index.'.interactors';
        unlink $factor_interactors;
        my $interact = tie(%interactors, 'DB_File', $factor_interactors, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$factor_interactors': $!");
        
        my %gene;
        my $factor_gene = $gene_index.'.factor';
        unlink $factor_gene;
        my $gene = tie(%gene, 'DB_File', $factor_gene, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$factor_gene': $!");
        
        my %matrix;
        my $factor_matrix = $matrix_index.'.factor';
        unlink $factor_matrix;
        my $matrix = tie(%matrix, 'DB_File', $factor_matrix, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$factor_matrix': $!");
        
        my %site;
        my $factor_site = $site_index.'.factor';
        unlink $factor_site;
        my $site = tie(%site, 'DB_File', $factor_site, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$factor_site': $!");
        
        my %fragment;
        my $factor_fragment = $fragment_index.'.factor';
        unlink $factor_fragment;
        my $fragment = tie(%fragment, 'DB_File', $factor_fragment, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$factor_fragment': $!");
        
        my %reference;
        my $factor_reference = $reference_index.'.factor';
        unlink $factor_reference;
        my $reference = tie(%reference, 'DB_File', $factor_reference, O_RDWR|O_CREAT, 0644, $DB_BTREE)
            or $self->throw("Could not open file '$factor_reference': $!");
        
        # skip the first three header lines
        <$FAC>; <$FAC>; <$FAC>;
        
        my @data;
        my $sequence = '';
        while (<$FAC>) {
            if (/^AC  (\S+)/) {
                $data[0] = $1;
            }
            elsif (/^ID  (\S+)/) {
                # IDs are always the same as AC? Is this needed?
                $data[1] = $1;
                $id->put("$1", $data[0]);
            }
            elsif (/^FA  (.+)$/) {
                $data[2] = $1;
                $name->put("$1", $data[0]);
            }
            elsif (/^OS  (.+)$/) {
                # This is the species the actual factor came from, which may
                # differ from the species of any sequences it is described as
                # binding to. Not all factors that have a species have a gene,
                # so can't delegate species to a gene lookup.
                my $raw_species = $1;
                my $taxid = $self->_species_to_taxid($raw_species);
                $data[3] = $taxid || $raw_species;
                $species->put($data[3], $data[0]);
            }
            elsif (/^GE  (G\d+)/) {
                $gene->put($data[0], "$1");
            }
            elsif (/^SQ  (.+)$/) {
                $sequence .= $1;
            }
            elsif (/^IN  (T\d+)/) {
                $interact->put($data[0], "$1");
            }
            elsif (/^MX  (M\d+)/) {
                $matrix->put($data[0], "$1");
            }
            elsif (/^BS  (R\d+)/) {
                $site->put($data[0], "$1");
            }
            elsif (/^BR  (FR\d+)/) {
                $fragment->put($data[0], "$1");
            }
            elsif (/^RN  .+?(RE\d+)/) {
                $reference->put($data[0], "$1");
            }
            elsif (/^\/\//) {
                # end of a record, store previous data and reset
                
                # accession = id name species sequence
                $factors{$data[0]} = join(SEPARATOR, ($data[1] || '',
                                                      $data[2] || '',
                                                      $data[3] || '',
                                                      $sequence));
                
                @data = ();
                $sequence = '';
            }
        }
        close $FAC;
        
        $factor = $id = $name = $species = $interact = $gene = $matrix = $site = $fragment = $reference = undef;
        untie %factors;
        untie %id;
        untie %name;
        untie %species;
        untie %interactors;
        untie %gene;
        untie %matrix;
        untie %site;
        untie %fragment;
        untie %reference;
    }
    
    my $fragment_dat = "$dat_dir/fragment.dat";
    if (! -e $fragment_index || $force) {
        if (open my $FRA, '<', $fragment_dat) {
            my %fragments;
            unlink $fragment_index;
            my $fragment = tie(%fragments, 'DB_File', $fragment_index, O_RDWR|O_CREAT, 0644, $DB_HASH)
                or $self->throw("Could not open file '$fragment_index': $!");
            
            my %id;
            my $fragment_id = $fragment_index.'.id';
            unlink $fragment_id;
            my $id = tie(%id, 'DB_File', $fragment_id, O_RDWR|O_CREAT, 0644, $DB_BTREE)
                or $self->throw("Could not open file '$fragment_id': $!");
            
            my %qualities;
            my $fragment_qualities = $fragment_index.'.qual';
            unlink $fragment_qualities;
            my $quality = tie(%qualities, 'DB_File', $fragment_qualities, O_RDWR|O_CREAT, 0644, $DB_HASH)
                or $self->throw("Could not open file '$fragment_qualities': $!");
            
            my %species;
            my $fragment_species = $fragment_index.'.species';
            unlink $fragment_species;
            my $species = tie(%species, 'DB_File', $fragment_species, O_RDWR|O_CREAT, 0644, $DB_BTREE)
                or $self->throw("Could not open file '$fragment_species': $!");
            
            my %gene;
            my $fragment_gene = $gene_index.'.fragment';
            unlink $fragment_gene;
            my $gene = tie(%gene, 'DB_File', $fragment_gene, O_RDWR|O_CREAT, 0644, $DB_BTREE)
                or $self->throw("Could not open file '$fragment_gene': $!");
            
            my %factor;
            my $fragment_factor = $factor_index.'.fragment';
            unlink $fragment_factor;
            my $factor = tie(%factor, 'DB_File', $fragment_factor, O_RDWR|O_CREAT, 0644, $DB_BTREE)
                or $self->throw("Could not open file '$fragment_factor': $!");
            
            my %reference;
            my $fragment_reference = $reference_index.'.fragment';
            unlink $fragment_reference;
            my $reference = tie(%reference, 'DB_File', $fragment_reference, O_RDWR|O_CREAT, 0644, $DB_BTREE)
                or $self->throw("Could not open file '$fragment_reference': $!");
            
            # skip the first three header lines
            <$FRA>; <$FRA>; <$FRA>;
            
            my @data;
            while (<$FRA>) {
                if (/^AC  (\S+)/) {
                    $data[0] = $1;
                }
                elsif (/^ID  (\S+)/) {
                    # IDs are always the same as AC? Is this needed?
                    $data[1] = $1;
                    $id->put("$1", $data[0]);
                }
                elsif (/^DE  Gene: (G\d+)(?:.+Gene: (G\d+))?/) {
                    my ($gene1, $gene2) = ($1, $2);
                    $data[2] = $gene1;
                    $data[3] = $gene2; # could be undef
                    $gene->put($data[0], $gene1);
                    $gene->put($data[0], $gene2) if $gene2;
                }
                elsif (/^OS  (.+)$/) {
                    # As per the site.dat parsing
                    my $raw_species = $1;
                    my $taxid = $self->_species_to_taxid($raw_species);
                    $data[4] = $taxid || $raw_species;
                    $species->put($data[4], $data[0]);
                }
                elsif (/^SQ  [atcgn]*([ATCGN]+)[atcgn]*/) {
                    $data[5] .= $1;
                    # there can be (usually are) multiple SQ lines with a single
                    # long seq split over them. The 'real' sequence is in caps
                }
                elsif (/^SC  Build (\S+):$/) {
                    $data[6] = $1;
                    # maybe parse it out a little more? We have build,
                    # chromosomal coords and strand, eg.
                    # SC  Build HSA_May2004: Chr.2 43976692..43978487 (FORWARD).
                }
                elsif (/^RN  .+?(RE\d+)/) {
                    $reference->put($data[0], "$1");
                }
                elsif (/^BF  (T\d+); .+?; Quality: (\d)/) {
                    $factor->put($data[0], "$1");
                    $qualities{$data[0].SEPARATOR.$1} = $2;
                }
                elsif (/^\/\//) {
                    # end of a record, store previous data and reset
                    
                    # accession = id gene_id1 gene_id2 species_tax_id_or_raw_string sequence source
                    $fragments{$data[0]} = join(SEPARATOR, ($data[1] || '',
                                                            $data[2] || '',
                                                            $data[3] || '',
                                                            $data[4] || '',
                                                            $data[5] || '',
                                                            $data[6] || ''));
                    
                    @data = ();
                }
            }
            close $FRA;
            
            $fragment = $id = $species = $quality = $gene = $factor = $reference = undef;
            untie %fragments;
            untie %id;
            untie %species;
            untie %qualities;
            untie %gene;
            untie %factor;
            untie %reference;
        }
        else {
            $self->warn("Could not read fragment file '$fragment_dat', assuming you have an old version of Transfac Pro with no fragment.dat file");
        }
    }
}

# connect the internal db handle
sub _db_connect {
    my $self = shift;
    return if $self->{'_initialized'};
    
    my $index_dir = $self->index_directory;
    my $gene_index = "$index_dir/gene.dat.index";
    my $reference_index = "$index_dir/reference.dat.index";
    my $matrix_index = "$index_dir/matrix.dat.index";
    my $factor_index = "$index_dir/factor.dat.index";
    my $site_index = "$index_dir/site.dat.index";
    my $fragment_index = "$index_dir/fragment.dat.index";
    
    foreach ($gene_index, $reference_index, $matrix_index, $factor_index, $site_index, $fragment_index) {
        if (! -e $_) {
            #$self->warn("Index files have not been created");
            #return 0;
        }
    }
    
    # reference
    {
        $self->{reference}->{data} = {};
        tie (%{$self->{reference}->{data}}, 'DB_File', $reference_index, O_RDWR, undef, $DB_HASH) || $self->throw("Cannot open file '$reference_index': $!");
        
        my $reference_pubmed = $reference_index.'.pubmed';
        $self->{reference}->{pubmed} = tie (%{$self->{reference}->{pubmed}}, 'DB_File', $reference_pubmed, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$reference_pubmed': $!");
        
        my $reference_gene = $gene_index.'.reference';
        $self->{gene}->{reference} = tie (%{$self->{gene}->{reference}}, 'DB_File', $reference_gene, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$reference_gene': $!");
        
        my $reference_site = $site_index.'.reference';
        $self->{site}->{reference} = tie (%{$self->{site}->{reference}}, 'DB_File', $reference_site, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$reference_site': $!");
        
        my $reference_fragment = $fragment_index.'.reference';
        $self->{fragment}->{reference} = tie (%{$self->{fragment}->{reference}}, 'DB_File', $reference_fragment, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$reference_fragment': $!");
        
        my $reference_factor = $factor_index.'.reference';
        $self->{factor}->{reference} = tie (%{$self->{factor}->{reference}}, 'DB_File', $reference_factor, undef, 0644, $DB_BTREE) || $self->throw("Cannot open file '$reference_factor': $!");
        
        my $reference_matrix = $matrix_index.'.reference';
        $self->{matrix}->{reference} = tie (%{$self->{matrix}->{reference}}, 'DB_File', $reference_matrix, undef, 0644, $DB_BTREE) || $self->throw("Cannot open file '$reference_matrix': $!");
    }
    
    # gene
    {
        $self->{gene}->{data} = {};
        tie (%{$self->{gene}->{data}}, 'DB_File', $gene_index, O_RDWR, undef, $DB_HASH) || $self->throw("Cannot open file '$gene_index': $!");
        
        my $gene_id = $gene_index.'.id';
        $self->{gene}->{id} = tie(%{$self->{gene}->{id}}, 'DB_File', $gene_id, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$gene_id': $!");
        
        my $gene_name = $gene_index.'.name';
        $self->{gene}->{name} = tie(%{$self->{gene}->{name}}, 'DB_File', $gene_name, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$gene_name': $!");
        
        my $gene_species = $gene_index.'.species';
        $self->{gene}->{species} = tie(%{$self->{gene}->{species}}, 'DB_File', $gene_species, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$gene_species': $!");
        
        my $gene_site = $site_index.'.gene';
        $self->{site}->{gene} = tie(%{$self->{site}->{gene}}, 'DB_File', $gene_site, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$gene_site': $!");
        
        my $gene_fragment = $fragment_index.'.gene';
        $self->{fragment}->{gene} = tie(%{$self->{fragment}->{gene}}, 'DB_File', $gene_fragment, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$gene_fragment': $!");
        
        my $gene_factor = $factor_index.'.gene';
        $self->{factor}->{gene} = tie(%{$self->{factor}->{gene}}, 'DB_File', $gene_factor, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$gene_factor': $!");
        
        my $gene_reference = $reference_index.'.gene';
        $self->{reference}->{gene} = tie(%{$self->{reference}->{gene}}, 'DB_File', $gene_reference, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$gene_reference': $!");
    }
    
    # site
    {
        $self->{site}->{data} = {};
        tie (%{$self->{site}->{data}}, 'DB_File', $site_index, O_RDWR, undef, $DB_HASH) || $self->throw("Cannot open file '$site_index': $!");
        
        my $site_id = $site_index.'.id';
        $self->{site}->{id} = tie(%{$self->{site}->{id}}, 'DB_File', $site_id, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$site_id': $!");
        
        my $site_species = $site_index.'.species';
        $self->{site}->{species} = tie(%{$self->{site}->{species}}, 'DB_File', $site_species, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file $site_species': $!");
        
        #*** quality not actually used by anything (yet)
        my $site_qualities = $site_index.'.qual';
        $self->{quality} = {};
        tie(%{$self->{quality}}, 'DB_File', $site_qualities, O_RDWR, undef, $DB_HASH) || $self->throw("Cannot open file '$site_qualities': $!");
        
        my $site_gene = $gene_index.'.site';
        $self->{gene}->{site} = tie(%{$self->{gene}->{site}}, 'DB_File', $site_gene, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$site_gene': $!");
        
        my $site_matrix = $matrix_index.'.site';
        $self->{matrix}->{site} = tie(%{$self->{matrix}->{site}}, 'DB_File', $site_matrix, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$site_matrix': $!");
        
        my $site_factor = $factor_index.'.site';
        $self->{factor}->{site} = tie(%{$self->{factor}->{site}}, 'DB_File', $site_factor, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$site_factor': $!");
        
        my $site_reference = $reference_index.'.site';
        $self->{reference}->{site} = tie(%{$self->{reference}->{site}}, 'DB_File', $site_reference, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$site_reference': $!");
    }
    
    # fragment (may not be in older databases)
    if (-e $fragment_index) {
        $self->{fragment}->{data} = {};
        tie (%{$self->{fragment}->{data}}, 'DB_File', $fragment_index, O_RDWR, undef, $DB_HASH) || $self->throw("Cannot open file '$fragment_index': $!");
        
        my $fragment_id = $fragment_index.'.id';
        $self->{fragment}->{id} = tie(%{$self->{fragment}->{id}}, 'DB_File', $fragment_id, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$fragment_id': $!");
        
        my $fragment_species = $fragment_index.'.species';
        $self->{fragment}->{species} = tie(%{$self->{fragment}->{species}}, 'DB_File', $fragment_species, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file $fragment_species': $!");
        
        #*** quality not actually used by anything (yet)
        my $fragment_qualities = $fragment_index.'.qual';
        $self->{fragment_quality} = {};
        tie(%{$self->{fragment_quality}}, 'DB_File', $fragment_qualities, O_RDWR, undef, $DB_HASH) || $self->throw("Cannot open file '$fragment_qualities': $!");
        
        my $fragment_gene = $gene_index.'.fragment';
        $self->{gene}->{fragment} = tie(%{$self->{gene}->{fragment}}, 'DB_File', $fragment_gene, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$fragment_gene': $!");
        
        my $fragment_factor = $factor_index.'.fragment';
        $self->{factor}->{fragment} = tie(%{$self->{factor}->{fragment}}, 'DB_File', $fragment_factor, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$fragment_factor': $!");
        
        my $fragment_reference = $reference_index.'.fragment';
        $self->{reference}->{fragment} = tie(%{$self->{reference}->{fragment}}, 'DB_File', $fragment_reference, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$fragment_reference': $!");
    }
    else {
        die "no fragment_index at '$fragment_index'\n";
    }
    
    # matrix
    {
        $self->{matrix}->{data} = {};
        tie (%{$self->{matrix}->{data}}, 'DB_File', $matrix_index, O_RDWR, undef, $DB_HASH) || $self->throw("Cannot open file '$matrix_index': $!");
        
        my $matrix_id = $matrix_index.'.id';
        $self->{matrix}->{id} = tie(%{$self->{matrix}->{id}}, 'DB_File', $matrix_id, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$matrix_id': $!");
        
        my $matrix_name = $matrix_index.'.name';
        $self->{matrix}->{name} = tie(%{$self->{matrix}->{name}}, 'DB_File', $matrix_name, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$matrix_name': $!");
        
        my $matrix_site = $site_index.'.matrix';
        $self->{site}->{matrix} = tie(%{$self->{site}->{matrix}}, 'DB_File', $matrix_site, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$matrix_site': $!");
        
        my $matrix_factor = $factor_index.'.matrix';
        $self->{factor}->{matrix} = tie(%{$self->{factor}->{matrix}}, 'DB_File', $matrix_factor, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$matrix_factor': $!");
        
        my $matrix_reference = $reference_index.'.matrix';
        $self->{reference}->{matrix} = tie(%{$self->{reference}->{matrix}}, 'DB_File', $matrix_reference, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$matrix_reference': $!");
    }
    
    # factor
    {
        $self->{factor}->{data} = {};
        tie (%{$self->{factor}->{data}}, 'DB_File', $factor_index, O_RDWR, undef, $DB_HASH) || $self->throw("Cannot open file '$factor_index': $!");
        
        my $factor_id = $factor_index.'.id';
        $self->{factor}->{id} = tie(%{$self->{factor}->{id}}, 'DB_File', $factor_id, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file 'factor_id': $!");
        
        my $factor_name = $factor_index.'.name';
        $self->{factor}->{name} = tie(%{$self->{factor}->{name}}, 'DB_File', $factor_name, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$factor_name': $!");
        
        my $factor_species = $factor_index.'.species';
        $self->{factor}->{species} = tie(%{$self->{factor}->{species}}, 'DB_File', $factor_species, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$factor_species': $!");
        
        my $factor_interactors = $factor_index.'.interactors';
        $self->{factor}->{interactors} = tie(%{$self->{factor}->{interactors}}, 'DB_File', $factor_interactors, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$factor_interactors': $!");
        
        my $factor_gene = $gene_index.'.factor';
        $self->{gene}->{factor} = tie(%{$self->{gene}->{factor}}, 'DB_File', $factor_gene, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$factor_gene': $!");
        
        my $factor_matrix = $matrix_index.'.factor';
        $self->{matrix}->{factor} = tie(%{$self->{matrix}->{factor}}, 'DB_File', $factor_matrix, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$factor_matrix': $!");
        
        my $factor_site = $site_index.'.factor';
        $self->{site}->{factor} = tie(%{$self->{site}->{factor}}, 'DB_File', $factor_site, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$factor_site': $!");
        
        my $factor_fragment = $fragment_index.'.factor';
        $self->{fragment}->{factor} = tie(%{$self->{fragment}->{factor}}, 'DB_File', $factor_fragment, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$factor_fragment': $!");
        
        my $factor_reference = $reference_index.'.factor';
        $self->{reference}->{factor} = tie(%{$self->{reference}->{factor}}, 'DB_File', $factor_reference, O_RDWR, undef, $DB_BTREE) || $self->throw("Cannot open file '$factor_reference': $!");
    }
    
    $self->{'_initialized'}  = 1;
}

=head2 index_directory

 Title   : index_directory
 Funtion : Get/set the location that index files are stored. (this module
           will index the supplied database)
 Usage   : $obj->index_directory($newval)
 Returns : value of index_directory (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub index_directory {
    my $self = shift;
    return $self->{'index_directory'} = shift if @_;
    return $self->{'index_directory'};
}

# resolve a transfac species string into an ncbi taxid
sub _species_to_taxid {
    my ($self, $raw_species) = @_;
    $raw_species or return;
    
    my $species_string;
    my @split = split(', ', $raw_species);
    (@split > 1) ? ($species_string = $split[1]) : ($species_string = $split[0]);
    
    my $ncbi_taxid;
    if ($species_string =~ /^[A-Z]\S+ \S+$/) {
        SWITCH: for ($species_string) {
            # some species don't classify so custom handling
            /^Darnel ryegrass/ && do { $ncbi_taxid = 34176; last; };
            /^Coix lacryma/ && do { $ncbi_taxid = 4505; last; };
            /^Rattus spec/ && do { $ncbi_taxid = 10116; last; };
            /^Mus spec/ && do { $ncbi_taxid = 10090; last; };
            /^Equus spec/ && do { $ncbi_taxid = 9796; last; };
            /^Cavia sp/ && do { $ncbi_taxid = 10141; last; };
            /^Marsh marigold/ && do { $ncbi_taxid = 3449; last; };
            /^Phalaenopsis sp/ && do { $ncbi_taxid = 36900; last; };
            /^Anthirrhinum majus/ && do { $ncbi_taxid = 4151; last; };
            /^Equus spec/ && do { $ncbi_taxid = 9796; last; };
            /^Lycopodium spec/ && do { $ncbi_taxid = 13840; last; };
            /^Autographa californica/ && do { $ncbi_taxid = 307456; last; };
            /^E26 AEV/ && do { $ncbi_taxid = 31920; last; };
            /^Pseudocentrotus miliaris/ && do { $ncbi_taxid = 7677; last; }; # the genus is 7677 but this species isn't there
            /^SL3-3 (?:retro)?virus/ && do { $ncbi_taxid = 53454; last; }; # 53454 is unclassified MLV-related, SL3-3 a variant of that?
            /^Petunia sp/ && do { $ncbi_taxid = 4104; last; };
        }
        if (! $ncbi_taxid && defined $self->{_tax_db}) {
            ($ncbi_taxid) = $self->{_tax_db}->get_taxonids($species_string);
        }
    }
    else {
        # some species lines are poorly formated so custom handling
        SWITCH: for ($raw_species) {
            # for speed, go by common first letters
            my $first_letter = substr($raw_species, 0, 1);
            
            $first_letter eq 'A' && do {
                /^Adiantum raddianum/ && do { $ncbi_taxid = 32168; last; };
                /^Avian sarcoma virus \(strain 17\)/ && do { $ncbi_taxid = 11877; last; };
                /^AMV/ && do { $ncbi_taxid = 11866; last; };
                /^AEV/ && do { $ncbi_taxid = 11861; last; };
                /^AS42|^Avian musculoaponeurotic/ && do { $ncbi_taxid = 11873; last; };
                /^Avian myelocytomatosis/ && do { $ncbi_taxid = 11869; last; };
                /^ASV 31/ && do { $ncbi_taxid = 35270; last; };
                /^A-MuLV/ && do { $ncbi_taxid = 188539; last; };
                /^Asparagus officinalis/ && do { $ncbi_taxid = 4686; last; };
                /^Agrobacterium tumefaciens/ && do { $ncbi_taxid = 358; last; };
                /^ALV/ && do { $ncbi_taxid = 11864; last; };
                /^AAV/ && do { $ncbi_taxid = 272636; last; };
                /^AKV MLV/ && do { $ncbi_taxid = 11791; last; };
                last;
            };
            
            $first_letter eq 'B' && do {
                /^BPV-1/ && do { $ncbi_taxid = 10559; last; };
                /^BKV/ && do { $ncbi_taxid = 10629; last; };
                /^Bolivian squirrel monkey/ && do { $ncbi_taxid = 39432; last; };
                last;
            };
            
            $first_letter eq 'C' && do {
                /^Cauliflower/ && do { $ncbi_taxid = 3715; last; };
                /^Chamek/ && do { $ncbi_taxid = 118643; last; };
                /^Candida albicans/ && do { $ncbi_taxid = 5476; last; };
                /^CaMV/ && do { $ncbi_taxid = 10641; last; };
                last;
            };
            
            $first_letter eq 'E' && do {
                /^Eucalyptus gunnii/ && do { $ncbi_taxid = 3933; last; };
                /^EBV, Epstein-Barr virus/ && do { $ncbi_taxid = 10376; last; };
                /^Eucalyptus globulus subsp. bicostata/ && do { $ncbi_taxid = 71272; last; };
                /^Eucalyptus globulus subsp. globulus/ && do { $ncbi_taxid = 71271; last; };
                last;
            };
            
            $first_letter eq 'F' && do {
                /^FBR MuLV/ && do { $ncbi_taxid = 11806; last; };
                /^FBJ MuLV/ && do { $ncbi_taxid = 11805; last; };
                /^FeLV|Feline leukemia/ && do { $ncbi_taxid = 11923; last; };
                /^Flaveria trinervia/ && do { $ncbi_taxid = 4227; last; };
                /^FSV/ && do { $ncbi_taxid = 11885; last; };
                /^F-MuLV/ && do { $ncbi_taxid = 11795; last; };
                last;
            };
            
            $first_letter eq 'H' && do {
                /^HSV-1/ && do { $ncbi_taxid = 10298; last; };
                /^HTLV-I/ && do { $ncbi_taxid = 11908; last; };
                /^HIV-1/ && do { $ncbi_taxid = 11676; last; };
                /^HPV-16/ && do { $ncbi_taxid = 333760; last; };
                /^HBV/ && do { $ncbi_taxid = 10407; last; };
                /^HBI/ && do { $ncbi_taxid = 11867; last; };
                /^HPV-8/ && do { $ncbi_taxid = 10579; last; };
                /^HPV-11/ && do { $ncbi_taxid = 10580; last; };
                /^HPV-18/ && do { $ncbi_taxid = 333761; last; };
                /^HCMV/ && do { $ncbi_taxid = 10359; last; };
                /^HSV/ && do { $ncbi_taxid = 126283; last; };
                /^HSV-2/ && do { $ncbi_taxid = 10310; last; };
                /^HCV/ && do { $ncbi_taxid = 11108; last; };
                /^HIV-2/ && do { $ncbi_taxid = 11709; last; };
                last;
            };
            
            $first_letter eq 'M' && do {
                /^MMTV/ && do { $ncbi_taxid = 11757; last; };
                /^Mo-MuLV/ && do { $ncbi_taxid = 11801; last; };
                /^MuLV/ && do { $ncbi_taxid = 11786; last; };
                /^MSV/ && do { $ncbi_taxid = 11802; last; };
                /^MC29/ && do { $ncbi_taxid = 11868; last; };
                /^MVM/ && do { $ncbi_taxid = 10794; last; };
                /^MH2E21/ && do { $ncbi_taxid = 11955; last; }; # 11955 is a species, presumably MH2E21 is the strain
                last;
            };
            
            $first_letter eq 'R' && do {
                /^Raphanus sativus/ && do { $ncbi_taxid = 3726; last; };
                /^REV-T/ && do { $ncbi_taxid = 11636; last; };
                /^RAV-0/ && do { $ncbi_taxid = 11867; last; }; # should be rous-associated virus 0 variant
                /^RSV/ && do { $ncbi_taxid = 11886; last; };
                /^RadLV/ && do { $ncbi_taxid = 31689; last; };
                /^RTBV/ && do { $ncbi_taxid = 10654; last; };
                last;
            };
            
            $first_letter eq 'S' && do {
                /^SV40/ && do { $ncbi_taxid = 10633; last; };
                /^Sesbania rostrata/ && do { $ncbi_taxid = 3895; last; };
                /^SIV/ && do { $ncbi_taxid = 11723; last; };
                /^Spinacia oleracea/ && do { $ncbi_taxid = 3562; last; };
                /^SCMV/ && do { $ncbi_taxid = 10364; last; }; # supposed to be AGM isolate
                last;
            };
            
            # and lower case
            $first_letter eq 'a' && do {
                /^adenovirus type 5/ && do { $ncbi_taxid = 28285; last; };
                /^adenovirus type 2/ && do { $ncbi_taxid = 10515; last; };
                /^adenovirus/ && do { $ncbi_taxid = 189831; last; }; # 189831 ('unclassified Adenoviridae') is the closest I can get, but this has no genus and is not a species
                last;
            };
            
            $first_letter eq 'b' && do {
                /^bell pepper/ && do { $ncbi_taxid = 4072; last; };
                /^baculovirus, Autographa californica/ && do { $ncbi_taxid = 46015; last; };
                /^broccoli/ && do { $ncbi_taxid = 36774; last; };
                /^barley/ && do { $ncbi_taxid = 112509; last; };
                last;
            };
            
            $first_letter eq 'c' && do {
                /^clawed frog/ && do { $ncbi_taxid = 8355; last; };
                /^chipmunk/ && do { $ncbi_taxid = 64680; last; };
                /^common tree shrew/ && do { $ncbi_taxid = 37347; last; };
                /^cat/ && do { $ncbi_taxid = 9685; last; };
                last;
            };
            
            # and misc
            /^NK24/ && do { $ncbi_taxid = 11955; last; };
            /^OK10/ && do { $ncbi_taxid = 11871; last; };
            /^Dendrobium grex/ && do { $ncbi_taxid = 84618; last; };
            /^KSHV/ && do { $ncbi_taxid = 37296; last; };
            /^Oncidium/ && do { $ncbi_taxid = 96474; last; };
            /^Japanese quail/ && do { $ncbi_taxid = 93934; last; };
            /^Nile tilapia/ && do { $ncbi_taxid = 8128; last; };
            /^GALV/ && do { $ncbi_taxid = 11840; last; };
            /^JCV/ && do { $ncbi_taxid = 10632; last; };
            /^LPV/ && do { $ncbi_taxid = 10574; last; };
            /^Py,/ && do { $ncbi_taxid = 36362; last; };
            /^DHBV/ && do { $ncbi_taxid = 12639; last; };
            /^VZV/ && do { $ncbi_taxid = 10335; last; };
            /^Vicia faba/ && do { $ncbi_taxid = 3906; last; };
            
            /^hamster/ && do { $ncbi_taxid = 10029; last; };
            /^sea urchin/ && do { $ncbi_taxid = 7668; last; };
            /^fruit fly/ && do { $ncbi_taxid = 7227; last; };
            /^halibut/ && do { $ncbi_taxid = 8267; last; };
            /^vaccinia virus/ && do { $ncbi_taxid = 10245; last; };
            /^taxonomic class Mammalia/ && do { $ncbi_taxid = 40674; last; }; # not a species
            /^taxonomic class Vertebrata/ && do { $ncbi_taxid = 7742; last; }; # not a species
            /^dog/ && do { $ncbi_taxid = 9615; last; };
            /^parsley/ && do { $ncbi_taxid = 4043; last; };
            /^mouse, Mus domesticus Torino/ && do { $ncbi_taxid = 10092; last; }; # 10092 is domesticus subspecies, but not the Torino strain
            /^lemur, Eulemur fulvus collaris/ && do { $ncbi_taxid = 47178; last; };
            /^red sea bream/ && do { $ncbi_taxid = 143350; last; };
            /^zebra finch/ && do { $ncbi_taxid = 59729; last; };
            /^mung bean/ && do { $ncbi_taxid = 3916; last; };
            /^soybean/ && do { $ncbi_taxid = 3847; last; };
            /^oat/ && do { $ncbi_taxid = 4498; last; };
            /^pseudorabies virus/ && do { $ncbi_taxid = 10345; last; };
        }
    }
    
    $self->warn("Didn't know what species '$raw_species' was, unable to classify") unless $ncbi_taxid;
    return $ncbi_taxid;
}

sub DESTROY {
    my $self = shift;
    # Destroy tied references to close filehandles
    # and allow proper temporary files deletion
    undef $self->{_tax_db}->{'_nodes'};
    undef $self->{_tax_db}->{'_id2name'};
    undef $self->{_tax_db}->{'_name2id'};
    undef $self->{_tax_db}->{'_parent2children'};
    undef $self->{_tax_db}->{'_parentbtree'};
}

1;
