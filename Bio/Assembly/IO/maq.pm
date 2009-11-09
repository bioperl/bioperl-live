# $Id$
#
# BioPerl module for Bio::Assembly::IO::maq
#
# Copyright by Mark A. Jensen
#
# You may distribute this module under the same terms as Perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::Assembly::IO::maq - Driver to read assembly files in maq format *ALPHA!*

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
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Mark A. Jensen

Email maj -at- fortinbras -dot- us

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
 Usage   : my $scaffold = $asmio->next_assembly()
 Function: return the assembly defined by the map and cns files
 Returns : Bio::Assembly::Scaffold object
 Args    : none

=cut

sub next_assembly {
    my $self = shift; # object reference
    
    # Create a new scaffold to hold the contigs
    my $scaffoldobj = Bio::Assembly::Scaffold->new(-source => $progname);
    
    # Contig and read related
    my $contigobj;
    my $iscontig = 1;
    my %contiginfo;
    my $isread = 0;
    my %readinfo;
    
    # Loop over all assembly file lines
    if ( $self->_parse_cns_file ) {
	$contiginfo{qualobj} = $self->_next_cons;
	$contiginfo{start} = $contiginfo{qualobj}->start;
	$contiginfo{end} = $contiginfo{qualobj}->end;
    }
    else { # just choke
	$self->throw("Associated maq consensus file is not available");
    }

    while ($_ = $self->_readline) {
        chomp;
	# mapview format parsing ; every line is a read...
	undef %readinfo;
	@readinfo{ qw(read_name chr posn strand insert_size,
	    paired_flag map_qual se_map_qual alt_map_qual,
	    num_mm_best_hit sum_qual_mm_best_hit zero_mm_hits,
	    one_mm_hits read_len seqstr qualstr) } = split(/\s+/);
	# sanger conversion
	my @qual = map { ord($_)-33 } split('', $readinfo{qualstr});
	$readinfo{seq} = Bio::Seq::Quality->new(
	    -id => $readinfo{read_name},
	    -seq => $readinfo{seqstr},
	    -qual => join(' ', @qual)
	    );

	if (!defined $contiginfo{start} || 
	    $readinfo{'posn'} > $contiginfo{end} ) {
	    # new contig
	    # close old contigobj, if nec.
	    if (defined $contiginfo{start}) {
		# store old contig object
		if ($contiginfo{'seqnum'} > 1) {
		    $self->_store_contig(\%contiginfo, $contigobj, $scaffoldobj);
		}
		else { #singlet
		    $self->_store_singlet(\%contiginfo, $contigobj, $scaffoldobj);
		    1;
		}
	    }
	    # create new contigobj
	    undef %contiginfo;
	    $contiginfo{'seqnum'} = 1;
	    $contiginfo{'asmbl_id'} = 'maq_assy';
	    $contiginfo{'qualobj'} = $self->_next_cons;
	    $contigobj = $self->_init_contig(\%contiginfo, $scaffoldobj);
	    # reset the contig trackers
	    $contiginfo{start} = $contiginfo{'qualobj'}->start;
	    $contiginfo{end} = $contiginfo{'qualobj'}->end;
	}
	else {
	    # update this contig's info...
	    $contiginfo{'seqnum'}++;
	}

    }
    # Store read info for last read
    if (defined $contiginfo{'seqnum'}) {
        if ($contiginfo{'seqnum'} > 1) {
            # This is a read in a contig
            my $readobj = $self->_store_read(\%readinfo, $contigobj);
        } elsif ($contiginfo{'seqnum'} == 1) {
            # This is a singlet
            my $singletobj = $self->_store_singlet(\%readinfo, \%contiginfo,
                $scaffoldobj);
        } else {
            # this shouldn't happen
            $self->throw("Unhandled exception");
        }
    }
    
    $scaffoldobj->update_seq_list();
    
    return $scaffoldobj;
}

=head2 _init_contig

    Title   : _init_contig
    Usage   : my $contigobj; $contigobj = $self->_init_contig(
              \%contiginfo, $scaffoldobj);
    Function: store information of a contig belonging to a scaffold in the
              appropriate object
    Returns : Bio::Assembly::Contig object
    Args    : hash, Bio::Assembly::Scaffold

=cut

sub _init_contig {
    my ($contiginfo, $scaffoldobj) = @_;
    # Create a contig and attach it to scaffold
    my $contigobj = Bio::Assembly::Contig->new(
        -id     => $$contiginfo{'asmbl_id'},
        -source => $progname,
        -strand => 1
    );
    $scaffoldobj->add_contig($contigobj);
    return $contigobj;
}

=head2 _store_contig

    Title   : _store_contig
    Usage   : my $contigobj; $contigobj = $self->_store_contig(
              \%contiginfo, $contigobj, $scaffoldobj);
    Function: store information of a contig belonging to a scaffold in the
              appropriate object
    Returns : Bio::Assembly::Contig object
    Args    : hash, Bio::Assembly::Contig, Bio::Assembly::Scaffold

=cut

sub _store_contig {
    my ($self, $contiginfo, $contigobj, $scaffoldobj) = @_;

    $self->throw("Contig object must be defined") unless $contigobj;

    $contigobj->set_consensus_quality($$contiginfo{qualobj});
    my $consensus = Bio::LocatableSeq->new(
        -id    => $$contiginfo{'asmbl_id'},
        -seq   => $$contiginfo{'qualobject'}->seq,
        -start => 1,
    );
    $contigobj->set_consensus_sequence($consensus);


    # Add other misc contig information as features of the contig
   # Add other misc read information as subsequence feature
   my @other = grep !/asmbl_id|qualobj/, keys %$contiginfo;
   my %other;
   @other{@other} = @$contiginfo{@other};
    my $contigtags = Bio::SeqFeature::Generic->new(
        -primary_tag => "_main_contig_feature:$$contiginfo{'asmbl_id'}",
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
    my ($fname, $dir, $suf) = fileparse($self->file);
    my $cnsf = File::Spec->catdir($dir, $fname, ".cns.fastq");
    return unless (-e $cnsf );
    my $fqio = Bio::SeqIO->new( -file => $cnsf );
    my $cns = $fqio->next_seq;
    # now, infer the contigs on the basis of quality values
    # - assuming quality of zero => no coverage
    my @qual = $cns->qual;
    # covered sites
    my @sites = grep { $qual[$_] > 0 } (0..$#qual);
    my @ranges = ($sites[0]+1);
    for my $i (1..$#sites) {
	if ($sites[$i]-$sites[$i-1]>1) {
	    push @ranges, $sites[$i-1]+1, $sites[$i]+1;
	}
    }
    push @ranges, $sites[-1];
    for (my $i = 0; $i<$#ranges; $i+=2) {
	push @cons, Bio::Seq::Quality->new(
	    -id => "${fname}/".$ranges[$i]."-".$ranges[$i+1],
	    -start => $ranges[$i],
	    -end => $ranges[$i+1],
	    -seq => $cns->subseq($ranges[$i], $ranges[$i+1]),
	    -qual => @{$cns->qual}[$ranges[$i]-1..$ranges[$i+1]-1]
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

sub _cons { @{shift->{'_cons'}} };


=head2 _next_cons()

=cut

sub _next_cons() { shift(@{shift->{'_cons'}}) }

=head2 _store_read

    Title   : _store_read
    Usage   : my $readobj = $self->_store_read(\%readinfo, $contigobj);
    Function: store information of a read belonging to a contig 
              in the appropriate object
    Returns : a Bio::LocatableSeq object
    Args    : hash, Bio::Assembly::Contig

=cut

# 	@readinfo{ qw(read_name chr posn strand insert_size,
# 	    paired_flag map_qual se_map_qual alt_map_qual,
# 	    num_mm_best_hit sum_qual_mm_best_hit zero_mm_hits,
# 	    one_mm_hits seqstr qualstr) } = split(/\s+/);

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
   $$readinfo{'aln_end'} = $$readinfo{'posn'} + length($$readinfo{'seqstr'})-1;

   my $alncoord = Bio::SeqFeature::Generic->new(
       -primary_tag => $readobj->id,
       -start       => $$readinfo{'aln_start'},
       -end         => $$readinfo{'aln_end'},
       -strand      => $$readinfo{'strand'},
       -qual        => join(' ', $$readinfo{seq}->qual),
       # check here, need to create contigs "by hand"...
       -tag         => { 'contig' => $contigobj->id() }
   );

   $contigobj->set_seq_coord($alncoord, $readobj);

   # Add other misc read information as subsequence feature
   my @other = grep !/seqstr|strand/, keys %$readinfo;
   my %other;
   @other{@other} = @$readinfo{@other};
   my $readtags = Bio::SeqFeature::Generic->new(
       -primary_tag => '_main_read_feature:'.$readobj->id,
       -start       => $$readinfo{'aln_start'},
       -end         => $$readinfo{'aln_end'},
       -strand      => $$readinfo{'strand'},
       # dumping ground:
       -tag         => \%other
   );
   $alncoord->add_sub_SeqFeature($readtags);

   return $readobj;
}

#### revamp for maq

=head2 _store_singlet

    Title   : _store_singlet
    Usage   : my $singletobj = $self->_store_read(\%readinfo, \%contiginfo,
                  $scaffoldobj);
    Function: store information of a singlet belonging to a scaffold in the appropriate object
    Returns : Bio::Assembly::Singlet
    Args    : hash, hash, Bio::Assembly::Scaffold

=cut

sub _store_singlet {
    my ($self, $readinfo, $contiginfo, $scaffoldobj) = @_;

    my $singletobj = Bio::Assembly::Singlet->new( -seqref => $$readinfo{seq} );
    $scaffoldobj->add_singlet($singletobj);

    # Add other misc contig information as features of the contig
   # Add other misc read information as subsequence feature
    my @other = grep !/asmbl_id|qualobj/, keys %$contiginfo;
    my %other;
    @other{@other} = @$contiginfo{@other};
    my $contigtags = Bio::SeqFeature::Generic->new(
        -primary_tag => "_main_contig_feature:$$contiginfo{asmbl_id}",
        -start       => 1,
         -end         => $singletobj->get_consensus_length(),
        -strand      => 1,
	# dumping ground:
	-tag         => \%other
    );
    $singletobj->add_features([ $contigtags ], 1);

    $$readinfo{'aln_start'} = $$readinfo{'posn'};
    $$readinfo{'aln_end'} = $$readinfo{'posn'} + length($$readinfo{'seqstr'})-1;
    $$readinfo{'strand'} = ($$readinfo{strand} eq '+' ? 1 : -1);
    my $alncoord = Bio::SeqFeature::Generic->new(
	-primary_tag => "_aligned_coord:$$readinfo{read_name}",
	-start       => $$readinfo{'aln_start'},
	-end         => $$readinfo{'aln_end'},
	-strand      => $$readinfo{'strand'},
	-tag         => { 'contig' => $$contiginfo{asmbl_id} }
	);
    $alncoord->attach_seq($singletobj->seqref);
    $singletobj->add_features([ $alncoord ], 0);

    # Add other misc read information as subsequence feature
    my @other = grep !/seqstr|strand/, keys %$readinfo;
    my %other;
    @other{@other} = @$readinfo{@other};
    my $readtags = Bio::SeqFeature::Generic->new(
	-primary_tag => "_main_read_feature:$$readinfo{read_name}",
	-start       => $$readinfo{'aln_start'},
	-end         => $$readinfo{'aln_end'},
	-strand      => $$readinfo{'strand'},
	# dumping ground:
	-tag         => \%other
	);
    $alncoord->add_sub_SeqFeature($readtags);

    return $singletobj;
}

###### writes -- need them??

=head2 write_assembly

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

1;
__END__
