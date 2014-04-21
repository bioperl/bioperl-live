#
# BioPerl module for Bio::Assembly::IO::sam
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Mark A. Jensen <maj -at- fortinbras -dot- us>
#
# Copyright Mark A. Jensen
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Assembly::IO::sam - An IO module for assemblies in Sam format *BETA*

=head1 SYNOPSIS

 $aio = Bio::Assembly::IO( -file => "mysam.bam",
                           -refdb => "myrefseqs.fas");
 $assy = $aio->next_assembly;

=head1 DESCRIPTION

This is a (currently) read-only IO module designed to convert
Sequence/Alignment Map (SAM; L<http://samtools.sourceforge.net/>)
formatted alignments to L<Bio::Assembly::Scaffold> representations,
containing .L<Bio::Assembly::Contig> and L<Bio::Assembly::Singlet>
objects. It uses lstein's L<Bio::DB::Sam> to parse binary formatted SAM
(.bam) files guided by a reference sequence fasta database.

B<NB>: C<Bio::DB::Sam> is not a BioPerl module; it can be obtained via
CPAN. It in turn requires the C<libbam> library; source can be
downloaded at L<http://samtools.sourceforge.net/>. 

=head1 DETAILS

=over

=item * Required files

A binary SAM (C<.bam>) alignment and a reference sequence database in
FASTA format are required. Various required indexes (C<.fai>, C<.bai>)
will be created as necessary (via L<Bio::DB::Sam>).

=item * Compressed files

...can be specified directly , if L<IO::Uncompress::Gunzip> is
installed. Get it from your local CPAN mirror.

=item * BAM vs. SAM

The input alignment should be in (possibly gzipped) binary SAM
(C<.bam>) format. If it isn't, you will get a message explaining how
to convert it, viz.:

 $ samtools view -Sb mysam.sam > mysam.bam

The bam file must also be sorted on coordinates: do

 $ samtools sort mysam.unsorted.bam > mysam.bam

=item * Contigs

Contigs are calculated by this module, using the 'coverage' feature of
the L<Bio::DB::Sam> object. A contig represents a contiguous portion
of a reference sequence having non-zero coverage at each base.

The bwa assembler (L<http://bio-bwa.sourceforge.net/>) can assign read
sequences to multiple reference sequence locations. The present
implementation currently assigns such reads only to the first contig
in which they appear.

=item * Consensus sequences

Consensus sequence and quality objects are calculated by this module,
using the C<pileup> callback feature of C<Bio::DB::Sam>. The consensus
is (currently) simply the residue at a position that has the maximum
sum of quality values. The consensus quality is the integer portion of
the simple average of quality values for the consensus residue.

=item * SeqFeatures

Read sequences stored in contigs are accompanied by the following
features:
 
 contig : name of associated contig
 cigar  : CIGAR string for this read

If the read is paired with a successfully mapped mate, these features
will also be available:

 mate_start  : coordinate of to which the mate was aligned
 mate_len    : length of mate read
 mate_strand : strand of mate (-1 or 1)
 insert_size : size of insert spanned by the mate pair

These features are obtained as follows:

 @ids = $contig->get_seq_ids;
 $an_id = $id[0]; # or whatever
 $seq = $contig->get_seq_by_name($an_id);
 # Bio::LocatableSeq's aren't SeqFeature containers, so...
 $feat = $contig->get_seq_feat_by_tag( 
            $seq, "_aligned_coord:".$s->id
         );
 ($cigar) = $feat->get_tag_values('cigar');
 # etc.

=back

=head1 TODO

=over

=item * Supporting both text SAM (TAM) and binary SAM (BAM)

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support

Please direct usage questions or support issues to the mailing list:

L<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and
reponsive experts will be able look at the problem and quickly
address it. Please include a thorough description of the problem
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Mark A. Jensen

Email maj -at- fortinbras -dot- us

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::Assembly::IO::sam;
use strict;
use warnings;

# Object preamble - inherits from Bio::Root::Root

use Bio::Seq::Quality;
use Bio::Seq::PrimaryQual;
use Bio::LocatableSeq;
use Bio::Assembly::IO;
use Bio::Assembly::Scaffold;
use Bio::Assembly::Contig;
use Bio::Assembly::Singlet;
use Bio::SeqIO;
use File::Spec;
use File::Basename;
use File::Temp qw(tempfile);
use Carp;
use Bio::Root::Root;
use base qw(Bio::Root::Root Bio::Assembly::IO Bio::Root::IO);

our $HAVE_IO_UNCOMPRESS;
BEGIN {
# check requirements
    unless ( eval "require Bio::DB::Sam; 1" ) {
        Bio::Root::Root->throw(__PACKAGE__.' requires installation of samtools (libbam) and Bio::DB::Sam (available on CPAN; not part of BioPerl)');
    }
    unless ( eval "require IO::Uncompress::Gunzip; \$HAVE_IO_UNCOMPRESS = 1") {
        Bio::Root::Root->warn("IO::Uncompress::Gunzip is not available; you'll have to do your decompression by hand.");
    }
}

my $progname = 'sam';

sub new {
    my $class = shift;
    my @args = @_;
    my $self = $class->SUPER::new(@args);
    my ($file, $refdb, $format) = $self->_rearrange([qw(FILE REFDB FORMAT)], @args);
    $self->file($file);
    $refdb && $self->refdb($refdb);
    $self->_init_sam() or croak( "Sam initialization failed" );
    $self->{_assigned} = {};
    return $self;
}


=head1 Bio::Assembly::IO compliance

=head2 next_assembly()

    Title   : next_assembly
    Usage   : my $scaffold = $asmio->next_assembly();
    Function: return the next assembly in the sam-formatted stream
    Returns : Bio::Assembly::Scaffold object
    Args    : none

=cut

sub next_assembly {
    my $self = shift;
    my @contig_set;
    # get Bio::DB::Sam object 
    # could add a refdb or fasfile attribute to contain the reference db name
    # iterate through seq_ids...
    my @refseq_ids = $self->sam->seq_ids;
    my $assembly = Bio::Assembly::Scaffold->new( -progname => $progname );

    foreach my $id (@refseq_ids) {
        #### this is choice 1: all refseqs into one assy...###
        $self->_current_refseq_id( $id );
        # Load contigs and singlets in the scaffold
        while ( my $obj = $self->next_contig()) {
            # Add contig /singlet to assembly
            if ($obj->isa('Bio::Assembly::Singlet')) { # a singlet
                $assembly->add_singlet($obj);
            } else { # a contig
                $assembly->add_contig($obj);
            }
        }
    }
    return $assembly;
}

=head2 next_contig()

    Title   : next_contig
    Usage   : my $contig = $asmio->next_contig();
    Function: return the next contig or singlet from the sam stream
    Returns : Bio::Assembly::Contig or Bio::Assembly::Singlet
    Args    : none

=cut

sub next_contig {
    my $self = shift;
    if (!defined $self->_current_contig_seg_idx) {
        $self->_current_contig_seg_idx(0);
    }
    else {
        $self->_current_contig_seg_idx( 1+$self->_current_contig_seg_idx );
    }
    unless ($self->_current_refseq_id) {
        croak("No current refseq id set");
    }
    my $contig_segs = $self->_segset($self->_current_refseq_id);
    unless ($contig_segs && @$contig_segs) {
        croak("No contig segset for current id '".$self->_current_refseq_id."'")
    }
    # each segment in @$contig_segs represents a contig within the
    # current refseq
    my $contig_seg = $$contig_segs[$self->_current_contig_seg_idx];
    return if (!defined $contig_seg); # iterator finished
    # each 'feature' in $contig_seg represents a read;
    # $seqio lets us iterate efficiently over the reads:
    my $seqio = $contig_seg->features(-iterator => 1);


    # Contig and read related
    my $contigobj = $self->_store_contig($contig_seg);
    my $numseq = 0;

    while ( my $read = $seqio->next_seq ) {
        if ($self->{_assigned}->{$read->name}) {
            next;
        }
        $self->{_assigned}->{$read->name}=1;
        $self->_store_read($read, $contigobj);
        $numseq++;
    }
    if ($numseq == 1) { # ooh! a singlet!
        $contigobj = $self->_store_singlet($contigobj);
    }
    return $contigobj;
}

=head2 write_assembly()

 Title   : write_assembly
 Usage   : 
 Function: not implemented (module currrently read-only)
 Returns : 
 Args    : 

=cut

sub write_assembly {
    my $self = shift;
    $self->throw( "This module is currently read-only" );
}

=head1 Internal

=head2  _store_contig()

    Title   : _store_contig
    Usage   : my $contigobj = $self->_store_contig(\%contiginfo);
    Function: create and load a contig object
    Returns : Bio::Assembly::Contig object
    Args    : Bio::DB::Sam::Segment object

=cut

sub _store_contig {
    my ($self, $contig_seg) = @_;

    # Create a contig 
    my $contigobj = Bio::Assembly::Contig->new(
        -id     => 'sam_assy['.$self->_basename.':'.$self->_current_refseq_id.']/'.$contig_seg->start.'-'.$contig_seg->end,
        -source => $progname,
        -strand => 1
    );
    my $consobj = $self->_calc_consensus($contig_seg);
    my $consensus_seq = Bio::LocatableSeq->new(
        -id    => $contigobj->id,
        -seq   => $consobj->seq,
        -start => 1,
    );
    $contigobj->set_consensus_sequence($consensus_seq);
    my $consensus_qual = Bio::Seq::PrimaryQual->new(
        -id    => $contigobj->id,
        -qual  => $consobj->qual,
        -start => 1,
    );
    $contigobj->set_consensus_quality($consensus_qual);

    # Add other misc contig information as subsequence feature
    #my @other = grep !/asmbl_id|end|qualobj|start/, keys %$contiginfo;
    #my %other;
    #@other{@other} = @$contiginfo{@other};
    #my $contigtags = Bio::SeqFeature::Generic->new(
    #    -primary     => '_main_contig_feature',
    #    -source      => $$contiginfo{'asmbl_id'},
    #    -start       => 1,
    #    -end         => $contig_seg->length,
    #    -strand      => 1,
    #    # dumping ground:
    #    -tag         => \%other
    #);
    #$contigobj->add_features([ $contigtags ], 1);

    return $contigobj;
}

=head2 _store_read()

    Title   : _store_read
    Usage   : my $readobj = $self->_store_read($readobj, $contigobj);
    Function: store information of a read belonging to a contig in a contig object
    Returns : Bio::LocatableSeq
    Args    : Bio::DB::Bam::AlignWrapper, Bio::Assembly::Contig

=cut

sub _store_read {
    my $self = shift;
    my ($read, $contigobj) = @_;
    my $readseq = Bio::LocatableSeq->new(
        -display_id => $read->name,
        -primary_id => $read->name,
        -seq        => $read->dna,
        -start      => 1,
        -strand     => $read->strand,
        -alphabet   => 'dna'
    );
    my $qual = Bio::Seq::PrimaryQual->new(
        -id         => $read->name."_qual",
        -qual       => [$read->qscore]
    );

    # add pair information
    my @pair_info;
    if ($read->proper_pair) { # mate also aligned
        @pair_info = (
            mate_start => $read->mate_start,
            mate_len   => $read->mate_len,
            mate_strand => $read->mstrand,
            insert_size => $read->isize
        );
    }

    my $alncoord = Bio::SeqFeature::Generic->new(
        -primary    => $read->name,
        -start      => $read->start,
        -end        => $read->end,
        -strand     => $read->strand,
        -qual       => join(' ',$read->qscore),
        -tag        => { 'contig' => $contigobj->id,
                         'cigar'  => $read->cigar_str,
                         @pair_info }
    );

    $contigobj->set_seq_coord($alncoord, $readseq);
    $contigobj->set_seq_qual( $readseq, $qual );

    #add other misc read info as subsequence feature:
    #my @other = grep !/aln_(?:end|start)|seq(?:str)?|strand/, keys %$readinfo;
    #my %other;
    #@other{@other} = @$readinfo{@other};
    #my $readtags = Bio::SeqFeature::Generic->new(
    #    -primary     => '_main_read_feature',
    #    -source      => $readobj->id,
    #    -start       => $$readinfo{'aln_start'},
    #    -end         => $$readinfo{'aln_end'},
    #    -strand      => $$readinfo{'strand'},
    #    # dumping ground:
    #    -tag         => \%other
    #);
    #$contigobj->get_features_collection->add_features([$readtags]);
    #$contigobj->get_features_collection->add_SeqFeature($alncoord, $readtags);

    return $readseq;
}

=head2  _store_singlet()

    Title   : _store_singlet
    Usage   : my $singletobj = $self->_store_singlet($contigobj);
    Function: convert a contig object containing a single read into
              a singlet object
    Returns : Bio::Assembly::Singlet
    Args    : Bio::Assembly::Contig (previously loaded with only one seq)

=cut

sub _store_singlet {
    my $self = shift;
    my ($contigobj) = @_;

    my ($readseq) = $contigobj->each_seq;

    my $singletobj = Bio::Assembly::Singlet->new( -id => $contigobj->id,
                                                  -seqref => $readseq );

# may want to attach this someday
#    my $qual = $contigobj->get_qual_by_name($readseq->id);    

    return $singletobj;
}

=head1 REALLY Internal

=head2 _init_sam()

 Title   : _init_sam
 Usage   : $self->_init_sam($fasfile)
 Function: obtain a Bio::DB::Sam parsing of the associated sam file
 Returns : true on success
 Args    : [optional] name of the fasta reference db (scalar string)
 Note    : The associated file can be plain text (.sam) or binary (.bam);
           If the fasta file is not specified, and no file is contained in 
           the refdb() attribute, a .fas file with the same
           basename as the sam file will be searched for.
           
=cut

sub _init_sam {
    my $self = shift;
    my $fasfile = shift;
    my $file = $self->file;
    my $sam;
    $fasfile ||= $self->refdb;
    $file =~ s/^[<>+]*//; # byebye parasitic mode chars
    my ($fname, $dir, $suf) = fileparse($file, ".sam", ".bam");
    $self->_set_from_args({ '_basename' => $fname }, 
                          -methods => [qw( _basename)],
                          -create => 1);
    if (!defined $fasfile) {
        for (qw( fas fa fasta fas.gz fa.gz fasta.gz )) {
            $fasfile = File::Spec->catdir($dir, $self->_basename.$_);
            last if -e $fasfile;
            undef $fasfile;
        }
    }
    unless (-e $fasfile) {
        croak( "Can't find associated reference fasta db" );
    }
    !$self->refdb && $self->refdb($fasfile);
    # compression
    if ($fasfile =~ /\.gz[^.]*$/) {
        unless ($HAVE_IO_UNCOMPRESS) {
            croak( "IO::Uncompress::Gunzip not available; can't decompress on the fly");
        }
        my ($tfh, $tf) = tempfile( UNLINK => 1);
        my $z = IO::Uncompress::Gunzip->new($fasfile) or croak("Can't expand: $@");
        while (<$z>) { print $tfh $_ }
        close $tfh;
        $fasfile = $tf;
    }
    if ($file =~ /\.gz[^.]*$/) {
        unless ($HAVE_IO_UNCOMPRESS) {
            croak( "IO::Uncompress::Gunzip not available; can't decompress on the fly");
        }
        my ($tfh, $tf) = tempfile( UNLINK => 1);
        my $z = IO::Uncompress::Gunzip->new($file) or croak("Can't expand: $@");
        while (<$z>) { 
            print $tfh $_;
        }
        close $tfh;
        $file = $tf;
    }
    # sam conversion : just barf for now
    if (-T $file) {
        my $bam = $file;
        $bam =~ s/\.sam/\.bam/;
        croak( "'$file' looks like a text file.\n\tTo convert to the required .bam (binary SAM) format, run\n\t\$ samtools view -Sb $file > $bam\n");
    }

    $sam = Bio::DB::Sam->new( -bam => $file, 
                              -fasta => $fasfile,
                              -autoindex  => 1,
                              -expand_flags => 1);
    unless (defined $sam) {
        croak( "Couldn't create the Bio::DB::Sam object" );
    }
    $self->{sam} = $sam;
    # now produce the contig segments for each seq_id...
    for ($sam->seq_ids) {
        my $seg = $sam->segment(-seq_id=>$_, -start=>1, -end=>$sam->length($_));
        ${$self->{_segset}}{$_} = [$self->_get_contig_segs_from_coverage($seg)];
    }
 
    return 1;

}

=head2 _get_contig_segs_from_coverage()

 Title   : _get_contig_segs_from_coverage
 Usage   : 
 Function: calculates separate contigs using coverage info 
           in the segment
 Returns : array of Bio::DB::Sam::Segment objects, representing
           each contig
 Args    : Bio::DB::Sam::Segment object

=cut

sub _get_contig_segs_from_coverage {
    my $self = shift;
    my $segment = shift;
    unless ($self->sam) {
        croak("Sam object not yet initialized (call _init_sam)");
    }
    unless ( ref($segment) =~ /Bio::DB::Sam::Segment/ ) {
        croak("Requires Bio::DB::Sam::Segment object at arg 1");
    }
    my ($cov, @covdata, @rngs, @segs);
    ($cov) = $segment->features('coverage');
    unless ($cov) {
        croak("No coverage data!");
    }
    @covdata = $cov->coverage;
    
    # calculate contigs: 
    my $in_contig;
    my ($lf_end,$rt_end);
    for (0..$#covdata) {
        if ($covdata[$_]) {
            if ($in_contig) {
                $rt_end = $_+1;
                next;
            }
            else {
                $in_contig = 1;
                # push previous range
                if (defined $lf_end && defined $rt_end) {
                    push @rngs, [$lf_end, $rt_end];
                }
                $lf_end = $_+1;
            } 
        }
        else {
            $in_contig = 0;
        }
    }
    # last one
    push @rngs, [$lf_end, $rt_end] if (defined $lf_end and 
                                       defined $rt_end and
                                       $lf_end <= $rt_end);
    unless (@rngs) {
        carp ("No coverage for this segment!");
        return;
    }
    for (@rngs) {
        push @segs, $self->sam->segment(-seq_id=>$segment->seq_id,
                                        -start=>$$_[0], 
                                        -end=>$$_[1]);
    }
    return @segs;
}

=head2 _calc_consensus_quality()

 Title   : _calc_consensus_quality
 Usage   : @qual = $aio->_calc_consensus_quality( $contig_seg );
 Function: calculate an average or other data-reduced quality
           over all sites represented by the features contained
           in a Bio::DB::Sam::Segment
 Returns : 
 Args    : a Bio::DB::Sam::Segment object

=cut

sub _calc_consensus_quality {
    # just an average over sites for now...
    my $self = shift;
    my $seg = shift;

    my @quals;
    my $region = $seg->seq_id.':'.$seg->start.'..'.$seg->end;
    my $qual_averager = sub {
        my ($seqid, $pos, $p) = @_;
        return unless ($seg->start <= $pos and $pos <= $seg->end);
        my $acc = 0;
        my $n = 0;
        for my $pileup (@$p) {
            my $b = $pileup->alignment;
            $acc += $b->qscore->[$pileup->qpos];
            $n++;
        }
        push @quals, int($acc/$n);
    };
    $self->sam->pileup($region, $qual_averager);
    return @quals;
}

=head2 _calc_consensus()

 Title   : _calc_consensus
 Usage   : @qual = $aio->_calc_consensus( $contig_seg );
 Function: calculate a simple quality-weighted consensus sequence
           for the segment
 Returns : a SeqWithQuality object
 Args    : a Bio::DB::Sam::Segment object

=cut

sub _calc_consensus {
    # just an average over sites for now...
    my $self = shift;
    my $seg = shift;

    my @quals;
    my $conseq ='';
    my $region = $seg->seq_id.':'.$seg->start.'-'.$seg->end;

    my $weighted_consensus = sub {
        my ($seqid, $pos, $p) = @_;
        return unless ($seg->start <= $pos and $pos <= $seg->end);
        my %wt_tbl;
        my %n;
        for my $pileup (@$p) {
            my $b = $pileup->alignment;
            my $res = substr($b->qseq,$pileup->qpos,1);
            $wt_tbl{$res} += $b->qscore->[$pileup->qpos] || 0;
            $n{$res} ||= 0;
            $n{$res}++;
        }
        # really simple
        my $c = (sort { $wt_tbl{$b}<=>$wt_tbl{$a} } keys %wt_tbl)[0];
        $conseq .= $c;
        push @quals, int($wt_tbl{$c}/$n{$c});
    };

    $self->sam->pileup($region, $weighted_consensus);
    return Bio::Seq::Quality->new(
        -qual => join(' ', @quals ),
        -seq => $conseq,
        -id => $region
    );
}

=head2 refdb()

 Title   : refdb
 Usage   : $obj->refdb($newval)
 Function: the name of the reference db fasta file
 Example : 
 Returns : value of refdb (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub refdb {
    my $self = shift;
    
    return $self->{'refdb'} = shift if @_;
    return $self->{'refdb'};
}

=head2 _segset()

 Title   : _segset
 Usage   : $segset_hashref = $self->_segset()
 Function: hash container for the Bio::DB::Sam::Segment objects that
           represent each set of contigs for each seq_id
           { $seq_id => [@contig_segments], ... }
 Example : 
 Returns : value of _segset (a hashref) if no arg, 
           or the arrayref of contig segments, if arg == a seq id
 Args    : [optional] seq id (scalar string)
 Note    : readonly; hash elt set in _init_sam()

=cut

sub _segset {
    my $self = shift;
    return $self->{'_segset'} unless @_;
    return ${$self->{'_segset'}}{shift()};
}

=head2 _current_refseq_id()

 Title   : _current_refseq_id
 Usage   : $obj->_current_refseq_id($newval)
 Function: the "current" reference sequence id
 Example : 
 Returns : value of _current_refseq (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub _current_refseq_id {
    my $self = shift;
    
    return $self->{'_current_refseq_id'} = shift if @_;
    return $self->{'_current_refseq_id'};
}

=head2 _current_contig_seg_idx()

 Title   : current_contig_seg_idx
 Usage   : $obj->current_contig_seg_idx($newval)
 Function: the "current" segment index in the "current" refseq
 Example : 
 Returns : value of current_contig_seg_idx (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub _current_contig_seg_idx {
    my $self = shift;
    
    return $self->{'_current_contig_seg_idx'} = shift if @_;
    return $self->{'_current_contig_seg_idx'};
}


=head2 sam()

 Title   : sam
 Usage   : $obj->sam($newval)
 Function: store the associated Bio::DB::Sam object
 Example : 
 Returns : value of sam (a Bio::DB::Sam object)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub sam {
    my $self = shift;
    return $self->{'sam'};
}

sub DESTROY {
    my $self = shift;
    undef $self->{'sam'};
    delete $self->{_segset}->{$_} foreach (keys %{$self->{_segset}});
}

1;
