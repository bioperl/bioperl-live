#
# BioPerl module for Bio::Assembly::IO::maq
#
# Copyright by Mark A. Jensen
#
# You may distribute this module under the same terms as Perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Assembly::IO::maq - Driver to read assembly files in maq format *BETA*

=head1 SYNOPSIS

    # convert the native maq map format to plain text
    $ maq mapview all.map > all.maq

    # Building an input stream
    use Bio::Assembly::IO;

    # Assembly loading methods
    my $asmio = Bio::Assembly::IO->new( -file   => 'all.maq',
                                        -format => 'maq' );
    my $scaffold = $asmio->next_assembly;


=head1 DESCRIPTION

This package loads and writes map information in/from C<maq> map files converted by the C<maq mapview> utility. This module is a driver module for
Bio::Assembly::IO input/output.

Parsing is based on Heng Li's description of C<maq mapview> output, found 
at the C<maq> manpage: L<http://maq.sourceforge.net/maq-manpage.shtml>. 

The basic C<maq> workflow is: map reads to a reference sequence (with
C<maq map>), then create a consensus from the map (with C<maq
assemble>). To read a complete assembly with this module, the
following files need to be available:

 [basename].maq : created by maq mapview [basename].map > [basename].maq
 [basename].cns.fastq : created as follows
    $ maq assemble [basename].cns [refseq].bfa [basename].map
    $ maq cns2fq [basename].cns > [basename].cns.fastq

C<maq> produces only one "contig"; all reads map to the reference
sequence, which covers everything. This module breaks the reads into
contigs by dividing the C<maq> consensus into pieces for which there
are contiguous non-zero quality values.

The module C<Bio::Tools::Run::Maq> will help in this process (eventually).

This module has no write capability.

=head2 Implementation

Assemblies are loaded into Bio::Assembly::Scaffold objects composed of
Bio::Assembly::Contig and Bio::Assembly::Singlet objects. Contigs are 
not explicitly specified in C<map> files; the division of the map into 
contigs is calculated in this module.

Additional assembly information is stored as features. Contig objects have
SeqFeature information associated with the primary_tag:

    _main_contig_feature:$contig_id -> misc contig information

Read objects have sub_seqFeature information associated with the
primary_tag:

    _main_read_feature:$read_id     -> misc read information

Singlets are contigs of a single sequence, as calculated within this module. 
They are cataloged separately, as specified in L<Bio::Assembly::Scaffold>.

=head1 TODO

=over

=item *

Add pod descriptions of maq descriptive data (currently SeqFeatures
added to each contig component)

=item *

Add features describing the aggregate status of reads and contigs based
on the maq "paired flag"

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

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

Report bugs to the BioPerl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Mark A. Jensen

Email maj -at- fortinbras -dot- us


=head1 CONTRIBUTORS

Further improvements by Florent Angly 
(florent dot angly at gmail dot com) 

=head1 ACKNOWLEDGEMENT

Code and some POD text ripped liberally from Florent Angly's
L<Bio::Assembly::IO::tigr>.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a "_".

=cut

package Bio::Assembly::IO::maq;

use strict;
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

use base qw(Bio::Assembly::IO);

# paired flag constants
use constant { 
    FF => 1, FR => 2, RF => 4, RR => 8, PE => 16, 
    XC => 32, UN => 64, CP => 18
};

my $progname = 'maq';

=head2 next_assembly

 Title   : next_assembly
 Usage   : $scaffold = $stream->next_assembly()
 Function: return the assembly defined by the map and cns files
 Returns : Bio::Assembly::Scaffold object
 Args    : none

=cut

sub next_assembly {
    my $self = shift;

    my $assembly = Bio::Assembly::Scaffold->new( -progname => $progname );

    # Load contigs and singlets in the scaffold
    while ( my $obj = $self->next_contig()) {
        # Add contig /singlet to assembly
        if ($obj->isa('Bio::Assembly::Singlet')) { # a singlet
            $assembly->add_singlet($obj);
        } else { # a contig
            $assembly->add_contig($obj);
        }
    }

    return $assembly;
}


=head2 next_contig

 Title   : next_contig
 Usage   : $scaffold = $stream->next_contig()
 Function: Returns the next contig or singlet in the ACE stream.
 Returns : a Bio::Assembly::Contig or Bio::Assembly::Single object
 Args    : none

=cut

sub next_contig {
    my $self = shift; # object reference

    # Read the file of consensus sequences if it has not already been done for
    # this Bio:::Assembly::IO stream already
    if (not defined $self->_cons) {
        $self->_parse_cns_file or
            $self->throw("Associated maq consensus file is not available");
    }

    # Contig and read related
    my $contigobj;
    my %contiginfo;

    # Loop over all assembly file lines
    while ($_ = $self->_readline) {
        chomp;
        next if /^$/;

        # mapview format parsing ; every line is a read...
        my %readinfo;
        @readinfo{ qw(read_name chr posn strand insert_size
            paired_flag map_qual se_map_qual alt_map_qual
            num_mm_best_hit sum_qual_mm_best_hit zero_mm_hits
            one_mm_hits read_len seqstr qualstr) } = split(/\s+/);

        # sanger conversion
        my @qual = map { ord($_)-33 } split('', $readinfo{qualstr});
        $readinfo{seq} = Bio::Seq::Quality->new(
            -id   => $readinfo{read_name},
            -seq  => $readinfo{seqstr},
            -qual => \@qual
            );

        if ( not defined $contiginfo{start} ) {
            # First read of new contig or singlet
            $contiginfo{'seqnum'}   = 1;
            $contiginfo{'qualobj'}  = $self->_next_cons;
            $contiginfo{'start'}    = $contiginfo{'qualobj'}->start;
            $contiginfo{'end'}      = $contiginfo{'qualobj'}->end;
            $contiginfo{'asmbl_id'} = 'maq_assy['.$self->_basename.']/'.$contiginfo{start}.'-'.$contiginfo{end};
            # It may be a singlet, but assume it's a contig for now
            $contigobj = $self->_init_contig(\%contiginfo);
            $self->_store_read(\%readinfo, $contigobj);
        } else {
            if ( $readinfo{'posn'} <= $contiginfo{end} ) {
                # Add read to existing contig
                $contiginfo{'seqnum'}++;
                $self->_store_read(\%readinfo, $contigobj);
            } else {
                # Read belongs in a new contig
                if ($contiginfo{'seqnum'} > 1) {
                    $self->_store_contig(\%contiginfo, $contigobj);
                }
                else { # singlet
                    # Create a new singlet object from the read info
                    $contigobj = $self->_store_singlet(\%contiginfo, $contigobj);
                }
                # do a pushback
                $self->_pushback($_);
                last;
            }
        }
    }

    return $contigobj;
}

=head2 _init_contig()

    Title   : _init_contig
    Usage   : my $contigobj; $contigobj = $self->_init_contig(
              \%contiginfo, $scaffoldobj);
    Function: store information of a contig belonging to a scaffold in the
              appropriate object
    Returns : Bio::Assembly::Contig object
    Args    : hash, Bio::Assembly::Scaffold

=cut

sub _init_contig {
    my ($self, $contiginfo) = @_;
    # Create a contig and attach it to scaffold
    my $contigobj = Bio::Assembly::Contig->new(
        -id     => $$contiginfo{'asmbl_id'},
        -source => $progname,
        -strand => 1
    );
    return $contigobj;
}

=head2 _store_contig()

    Title   : _store_contig
    Usage   : my $contigobj; $contigobj = $self->_store_contig(
              \%contiginfo, $contigobj);
    Function: store information of a contig belonging to a scaffold
              in the appropriate object
    Returns : Bio::Assembly::Contig object
    Args    : hash, Bio::Assembly::Contig

=cut

sub _store_contig {
    my ($self, $contiginfo, $contigobj) = @_;

    $self->throw("Contig object must be defined") unless $contigobj;

    my $consensus_seq = Bio::LocatableSeq->new(
        -id    => $$contiginfo{'asmbl_id'},
        -seq   => $$contiginfo{'qualobj'}->seq,
        -start => 1,
    );
    $contigobj->set_consensus_sequence($consensus_seq);
    my $consensus_qual = Bio::Seq::PrimaryQual->new(
        -id    => $$contiginfo{'asmbl_id'},
        -qual  => $$contiginfo{'qualobj'}->qual,
        -start => 1,
    );
    $contigobj->set_consensus_quality($consensus_qual);

    # Add other misc contig information as features of the contig
    # Add other misc read information as subsequence feature
    my @other = grep !/asmbl_id|end|qualobj|start/, keys %$contiginfo;
    my %other;
    @other{@other} = @$contiginfo{@other};
    my $contigtags = Bio::SeqFeature::Generic->new(
        -primary     => '_main_contig_feature',
        -source      => $$contiginfo{'asmbl_id'},
        -start       => 1,
        -end         => $contigobj->get_consensus_length(),
        -strand      => 1,
        # dumping ground:
        -tag         => \%other
    );
    $contigobj->add_features([ $contigtags ], 1);

    return $contigobj;
}

=head2 _parse_cns_file()

 Title   : _parse_cns_file
 Usage   : $self->_parse_cns_file
 Function: parse the .cns.fastq (consensus) file
           associated with the present map;
           set the objects cns attribute
 Returns : true on success; nil if file dne
 Args    : none

=cut

sub _parse_cns_file {
    my ($self) = @_;
    my @cons;

    $self->{'_cns_parsed'} = 1;
    my $file = $self->file;
    $file =~ s/^[<>+]*//; # byebye parasitic mode chars
    my ($fname, $dir, $suf) = fileparse($file, ".maq");
    my $cnsf = File::Spec->catdir($dir, "$fname.cns.fastq");
    return unless (-e $cnsf );
    my $fqio = Bio::SeqIO->new( -file => $cnsf );
    my $cns = $fqio->next_seq;
    # now, infer the contigs on the basis of quality values
    # - assuming quality of zero => no coverage
    my $qual = $cns->qual;
    # covered sites
    my @sites = grep { $$qual[$_] > 0 } (0..$#$qual);
    my @ranges = ($sites[0]+1);
    for my $i (1..$#sites) {
        if ($sites[$i]-$sites[$i-1]>1) {
            push @ranges, $sites[$i-1]+1, $sites[$i]+1;
        }
    }
    push @ranges, $sites[-1];
    for (my $i = 0; $i<$#ranges; $i+=2) {
        push @cons, Bio::Seq::Quality->new(
            -display_id => "${fname}/".$ranges[$i]."-".$ranges[$i+1],
            -start => $ranges[$i],
            -end => $ranges[$i+1],
            -seq => $cns->subseq($ranges[$i], $ranges[$i+1]),
            -qual => [@{$cns->qual}[$ranges[$i]-1..$ranges[$i+1]-1]]
            );
    }

    $self->{'_cons'} = \@cons;
    return 1;
}


=head2 _cons()

 Title   : _cons
 Usage   : @cons = $self->_cons
 Function: get the array of consensus fastq Bio::Seq::Quality objects
 Returns : array of Bio::Seq::Quality objects
 Args    : none

=cut

sub _cons {
    my $self = shift;
    my $cons = undef;
    if (defined $self->{'_cons'}) {
      $cons = $self->{'_cons'};
    }
    return $cons;
}


=head2 _next_cons()

=cut

sub _next_cons { shift(@{shift->{'_cons'}}) }


=head2 _store_read()

    Title   : _store_read
    Usage   : my $readobj = $self->_store_read(\%readinfo, $contigobj);
    Function: store information of a read belonging to a contig 
              in the appropriate object
    Returns : a Bio::LocatableSeq object
    Args    : hash, Bio::Assembly::Contig

=cut

#       @readinfo{ qw(read_name chr posn strand insert_size,
#           paired_flag map_qual se_map_qual alt_map_qual,
#           num_mm_best_hit sum_qual_mm_best_hit zero_mm_hits,
#           one_mm_hits seqstr qualstr) } = split(/\s+/);

sub _store_read {
   my ($self, $readinfo, $contigobj) = @_;

   # Create an aligned read object

   $$readinfo{'strand'} = ($$readinfo{strand} eq '+' ? 1 : -1);

   my $readobj = Bio::LocatableSeq->new(
       -display_id => $$readinfo{'read_name'},
       -primary_id => $$readinfo{'read_name'},
       -seq        => $$readinfo{'seqstr'},
       -start      => 1,
       -strand     => $$readinfo{'strand'},
       -alphabet   => 'dna'
   );

   # Add read location and sequence to contig (in 'gapped consensus' coordinates)
   $$readinfo{'aln_start'} = $$readinfo{'posn'};
   $$readinfo{'aln_end'}   = $$readinfo{'posn'} + length($$readinfo{'seqstr'})-1;

   my $alncoord = Bio::SeqFeature::Generic->new(
       -primary     => $readobj->id,
       -start       => $$readinfo{'aln_start'},
       -end         => $$readinfo{'aln_end'},
       -strand      => $$readinfo{'strand'},
       -qual        => join(' ', $$readinfo{seq}->qual),
       # check here, need to create contigs "by hand"...
       -tag         => { 'contig' => $contigobj->id() }
   );

   $contigobj->set_seq_coord($alncoord, $readobj);

   # Add other misc read information as subsequence feature
   my @other = grep !/aln_(?:end|start)|seq(?:str)?|strand/, keys %$readinfo;
   my %other;
   @other{@other} = @$readinfo{@other};
   my $readtags = Bio::SeqFeature::Generic->new(
       -primary     => '_main_read_feature',
       -source      => $readobj->id,
       -start       => $$readinfo{'aln_start'},
       -end         => $$readinfo{'aln_end'},
       -strand      => $$readinfo{'strand'},
       # dumping ground:
       -tag         => \%other
   );
   $contigobj->get_features_collection->add_features([$readtags]);
   $contigobj->get_features_collection->add_SeqFeature($alncoord, $readtags);

   return $readobj;
}

#### revamp for maq

=head2 _store_singlet()

    Title   : _store_singlet
    Usage   : my $singletobj = $self->_store_read(\%readinfo, \%contiginfo);
    Function: store information of a singlet belonging to a scaffold in a singlet object
    Returns : Bio::Assembly::Singlet
    Args    : hash, hash

=cut

sub _store_singlet {
    my ($self, $contiginfo, $contigobj) = @_;

    my $contigid = $$contiginfo{'asmbl_id'};
    my $seqref = ($contigobj->each_seq())[0];
    my $singletobj = Bio::Assembly::Singlet->new( -id     => $contigid,
                                                  -seqref => $seqref   );

    # Add other misc contig information as features of the contig
    # Add other misc read information as subsequence feature
    #my @other = grep !/_sfc|_assembly|_elem/, keys %$contiginfo; # remove the objects; _elem contains a code ref and can't be frozen. Just shooting blind here.
    #my %other;
    #@other{@other} = @$contiginfo{@other};
    #my $contigtags = Bio::SeqFeature::Generic->new(
    #    -primary     => '_main_contig_feature',
    #    -source      => $$contiginfo{asmbl_id},
    #    -start       => 1,
    #    -end         => $singletobj->get_consensus_length(),
    #    -strand      => 1,
    #    # dumping ground:
    #    -tag         => \%other
    #);
    #$singletobj->add_features([ $contigtags ], 1);

    #$$readinfo{'aln_start'} = $$readinfo{'start'};
    #$$readinfo{'aln_end'} = $$readinfo{'end'};
    #$$readinfo{'strand'} = ($$readinfo{strand} eq '+' ? 1 : -1);
    #my $alncoord = Bio::SeqFeature::Generic->new(
    #    -primary     => '_aligned_coord',
    #    -source      => $$readinfo{read_name},
    #    -start       => $$readinfo{'start'},
    #    -end         => $$readinfo{'end'},
    #    -strand      => $$readinfo{'strand'},
    #    -tag         => { 'contig' => $$contiginfo{asmbl_id} }
    #    );
    #$alncoord->attach_seq($singletobj->seqref);
    #$singletobj->add_features([ $alncoord ], 0);

    # Add other misc read information as subsequence feature
    #my @other = grep !/seqstr|strand/, keys %$readinfo;
    #my %other;
    #@other{@other} = @$readinfo{@other};
    #my $readtags = Bio::SeqFeature::Generic->new(
    #    -primary     => '_main_read_feature',
    #    -source      => $$readinfo{read_name},
    #    -start       => $$readinfo{'aln_start'},
    #    -end         => $$readinfo{'aln_end'},
    #    -strand      => $$readinfo{'strand'},
    #    # dumping ground:
    #    -tag         => \%other
    #    );
    #$singletobj->get_features_collection->add_features([$readtags]);
    #$singletobj->get_features_collection->add_SeqFeature($alncoord, $readtags);

    return $singletobj;
}

###### writes -- need them??

=head2 write_assembly()

    Title   : write_assembly
    Usage   : 
    Function: not currently available for maq assemblies
    Returns : throw
    Args    : 

=cut

sub write_assembly {
    my ($self,@args) = @_;
    $self->throw("Writes not currently available for maq assemblies. Complain to author.")
}



=head2 _basename()

 Title   : _basename
 Usage   : $self->_basename
 Function: return the basename of the associate IO file
 Returns : scalar string
 Args    : none

=cut

sub _basename {
    my $self = shift;
    return (fileparse($self->file, ".maq"))[0];
}

1;
__END__
