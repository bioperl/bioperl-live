#
# BioPerl module for Bio::Assembly::IO::sam
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Dan Kortschak <dan.kortschak adelaide.edu.au>
#
# Copyright Dan Kortschak
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Assembly::IO::bowtie - An IO module for assemblies in Bowtie format *BETA*

=head1 SYNOPSIS

 $aio = Bio::Assembly::IO( -file   => "mybowtie.bowtie",
                           -index  => "myindex",
                           -format => "bowtie");
 $assy = $aio->next_assembly;

=head1 DESCRIPTION

This is a read-only IO module designed to convert Bowtie
(L<http://bowtie-bio.sourceforge.net/>) formatted alignments to
L<Bio::Assembly::Scaffold> representations, containing 
L<Bio::Assembly::Contig> and L<Bio::Assembly::Singlet> objects.
It is a wrapper that converts the Bowtie format to BAM format taken
by the Bio::Assembly::IO::sam module which in turn uses lstein's
L<Bio::DB::Sam> to parse binary formatted SAM (.bam) files guided by a
reference sequence fasta database.

Some information is lost in conversion from bowtie format to SAM/BAM format
that is provided by Bowtie using the SAM output option and the conversion
to SAM format from bowtie format is slower than using bowtie's SAM option.
If you plan to use SAM/BAM format it is preferable to use this Bowtie
option rather than convert the format after the fact.

See the Bio::Assembly::IO::sam documentation for relevant details.

=head1 DETAILS

=over

=item * Required files

A bowtie (C<.bowtie>) alignment and the bowtie index or fasta
file used to generate the alignment are required.

=item * Compressed files

...can be specified directly , if L<IO::Uncompress::Gunzip> is
installed. Get it from your local CPAN mirror.

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

=head1 AUTHOR - Dan Kortschak

Email dan.kortschak adelaide.edu.au

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::Assembly::IO::bowtie;
use strict;
use warnings;

# Object preamble - inherits from Bio::Root::Root

use Bio::SeqIO;
use Bio::Tools::Run::Samtools;
use Bio::Assembly::IO;

use base qw( Bio::Assembly::IO );

our $HD = "\@HD\tVN:1.0\tSO:unsorted\n";
our $PG = "\@PG\tID=Bowtie\n";

our $HAVE_IO_UNCOMPRESS;
our $HAVE_BOWTIE;

BEGIN {
# check requirements
    eval "require Bio::Tools::Run::Bowtie; \$HAVE_BOWTIE = 1";
    unless ( eval "require IO::Uncompress::Gunzip; \$HAVE_IO_UNCOMPRESS = 1") {
        Bio::Root::Root->warn("IO::Uncompress::Gunzip is not available; you'll have to do your decompression by hand.");
    }
}

=head2 new()

 Title   : new
 Usage   : my $obj = new Bio::Assembly::IO::bowtie();
 Function: Builds a new Bio::Assembly::IO object
 Returns : an instance of Bio::Assembly::IO
 Args    : hash of options:

            -file    => bowtie_output_file
            -index   => bowtie_index or fasta_file used to create index
            -no_head => boolean skip SAM header
            -no_sq   => boolean skip SQ lines of SAM header

 Note    : bowtie_output and fasta files may be gzipped

=cut

sub new {
    my $class = shift;
    my @args = @_;
    my $self = $class->SUPER::new(@args);
    $self->_initialize(@args);
    $self->{'_tempdir'} = $self->tempdir(CLEANUP=>1);
    my ($file, $index, $no_head, $no_sq) = $self->_rearrange([qw(FILE INDEX NO_HEAD NO_SQ)], @args);
    $file =~ s/^<//;
    $self->{'_no_head'} = $no_head;
    $self->{'_no_sq'} = $no_sq;

    # get the sequence so Bio::DB::Sam can work with it
    my $refdb;
    my $inspector;
    if (-e $index && -r _ ) {
        $refdb = ($index =~ m/\.gz[^.]*$/) ? $self->_uncompress($index) : $index;
        my $guesser = Bio::Tools::GuessSeqFormat->new(-file=>$refdb);
        $self->throw("'$index' is not a fasta file.")
            unless $guesser->guess =~ m/^fasta$/;
    } elsif ($HAVE_BOWTIE) {
        $inspector = Bio::Tools::Run::Bowtie->new( -command => 'inspect' );
        $refdb = $inspector->run($index);
    } else {
        $self->throw("Bio::Tools::Run::Bowtie is not available - cannot extract refdb from index.");
    }

    my $bam_file = $self->_make_bam($self->_bowtie_to_sam($file, $refdb));
    my $sam = Bio::Assembly::IO->new( -file => "<$bam_file", -refdb => $refdb , -format => 'sam' );

    return $sam;
}

sub _bowtie_to_sam {
    my ($self, $file, $refdb) = @_;

    $self->throw("'$file' does not exist or is not readable.")
        unless ( -e $file && -r _ );

    if ($file =~ m/\.gz[^.]*$/) {
        $file = $self->_uncompress($file);
        $self->close;
        open my $fh, '<', $file or $self->throw("Could not read file '$file': $!");
        $self->file($file);
        $self->_fh($fh);
    }

    my $guesser = Bio::Tools::GuessSeqFormat->new(-file=>$file);
    $self->throw("'$file' is not a bowtie formatted file.") unless $guesser->guess =~ m/^bowtie$/;

    my %SQ;
    my $mapq = 255;
    my $in_pair;
    my @mate_line;
    my $mlen;

    # create temp file for working
    my ($sam_tmp_h, $sam_tmp_f) = $self->tempfile( -dir => $self->{'_tempdir'}, -suffix => '.sam' );
    
    while (my $line = $self->_readline) {
        chomp($line);
        my ($qname,$strand,$rname,$pos,$seq,$qual,$m,$details)=split("\cI",$line);
        $SQ{$rname} = 1;
        
        my $paired_f =  ($qname =~ m#/[12]#) ? 0x03 : 0;
        my $strand_f = ($strand eq '-') ? 0x10 : 0;
        my $op_strand_f = ($strand eq '+' && $paired_f) ? 0x20 : 0;
        my $first_f =  ($qname =~ m#/1#) ? 0x40 : 0;
        my $second_f = ($qname =~ m#/2#) ? 0x80 : 0;
        my $flag = $paired_f | $strand_f | $op_strand_f | $first_f | $second_f;

        $pos++;
        my $len = length $seq;
        die unless $len == length $qual;
        my $cigar = $len.'M';
        my @detail = split(',',$details);
        my $dist = 'NM:i:'.scalar @detail;

        my @mismatch;
        my $last_pos = 0;
        for (@detail) {
            m/(\d+):(\w)>\w/;
            my $err = ($1-$last_pos);
            $last_pos = $1+1;
            push @mismatch,($err,$2);
        }
        push @mismatch, $len-$last_pos;
        @mismatch = reverse @mismatch if $strand eq '-';
        my $mismatch = join('',('MD:Z:',@mismatch));

        if ($paired_f) {
            my $mrnm = '=';
            if ($in_pair) {
                my $mpos = $mate_line[3];
                $mate_line[7] = $pos;
                my $isize = $mpos-$pos-$len;
                $mate_line[8] = -$isize;
                print $sam_tmp_h join("\t",@mate_line),"\n";
                print $sam_tmp_h join("\t",$qname, $flag, $rname, $pos, $mapq, $cigar, $mrnm, $mpos, $isize, $seq, $qual, $mismatch, $dist),"\n";
                $in_pair = 0;
            } else {
                $mlen = $len;
                @mate_line = ($qname, $flag, $rname, $pos, $mapq, $cigar, $mrnm, undef, undef, $seq, $qual, $mismatch, $dist);
                $in_pair = 1;
            }
        } else {
            my $mrnm = '*';
            my $mpos = 0;
            my $isize = 0;
            print $sam_tmp_h join("\t",$qname, $flag, $rname, $pos, $mapq, $cigar, $mrnm, $mpos, $isize, $seq, $qual, $mismatch, $dist),"\n";
        }
    }

    $sam_tmp_h->close;

    return $sam_tmp_f if $self->{'_no_head'};

    my ($samh, $samf) = $self->tempfile( -dir => $self->{'_tempdir'}, -suffix => '.sam' );

    # print header
    print $samh $HD;

    # print sequence dictionary
    unless ($self->{'_no_sq'}) {
        my $db  = Bio::SeqIO->new( -file => $refdb, -format => 'fasta' );
        while ( my $seq = $db->next_seq() ) {
            $SQ{$seq->id} = $seq->length if $SQ{$seq->id};
        }
    
        map { print $samh join("\t", ('@SQ', "SN:$_", "LN:$SQ{$_}")), "\n" } keys %SQ;
    }
    
    # print program
    print $samh $PG;
    
    # print alignments
    open $sam_tmp_h, '<', $sam_tmp_f or
        $self->throw("Could not read file '$sam_tmp_f': $!");

    print $samh $_ while (<$sam_tmp_h>);
    
    close($sam_tmp_h);
    $samh->close;
    
    return $samf;
}

sub _uncompress {
    my ($self, $file) = @_;

    unless ($HAVE_IO_UNCOMPRESS) {
        croak( "IO::Uncompress::Gunzip not available, can't expand '$file'" );
    }
    my ($tfh, $tf) = $self->tempfile( -dir => $self->{'_tempdir'} );
    my $z = IO::Uncompress::Gunzip->new($file);
    while (my $block = $z->getline) { print $tfh $block }
    close $tfh;
    $file = $tf;

    return $file
}

sub _make_bam {
    my ($self, $file) = @_;
    
    $self->throw("'$file' does not exist or is not readable")
        unless ( -e $file && -r _ );

    # make a sorted bam file from a sam file input
    my ($bamh, $bamf) = $self->tempfile( -dir => $self->{'_tempdir'}, -suffix => '.bam' );
    my ($srth, $srtf) = $self->tempfile( -dir => $self->{'_tempdir'}, -suffix => '.srt' );
    $_->close for ($bamh, $srth);
    
    my $samt = Bio::Tools::Run::Samtools->new( -command => 'view',
                                               -sam_input => 1,
                                               -bam_output => 1 );

    $samt->run( -bam => $file, -out => $bamf );

    $samt = Bio::Tools::Run::Samtools->new( -command => 'sort' );

    $samt->run( -bam => $bamf, -pfx => $srtf);

    return $srtf.'.bam'
}

1;
