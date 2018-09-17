#
# BioPerl module for Bio::AlignIO::arp
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::AlignIO::arp - ARP MSA Sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::AlignIO> 
class.

=head1 DESCRIPTION

This object can create L<Bio::SimpleAlign> objects from
ARP flat files.  These are typically configuration-like data files
for the program Arlequin.  For more information, see:

  http://lgb.unige.ch/arlequin/

For the moment, this retains the allele sequence data in the DATA section and
inserts them into SimpleAlign objects. ARP files that contain other data (RFLP,
etc.) are not expected to parse properly.  Also, if the DNA data is actually SNP
data, then the LocatableSeq object instantiation will throw an error.

This is now set up as a generic parser (i.e. it parses everything) and
collects as much data as possible into the SimpleAlign object.  The following
in a general mapping of where data can be found:

    Tag        SimpleAlign
               Method  
    ----------------------------------------------------------------------
    Title      description
    SampleName id  
    ----------------------------------------------------------------------

    Tag        Bio::Annotation   TagName                    Bio::Annotation
               Class                                        Parameters
    ----------------------------------------------------------------------
     NE        SimpleValue       pfam_family_accession      value
     NL        SimpleValue       sequence_start_stop        value
     SS        SimpleValue       sec_structure_source       value
     BM        SimpleValue       build_model                value
     RN        Reference         reference                  *
    ----------------------------------------------------------------------
  * RN is generated based on the number of Bio::Annotation::Reference objects

In addition, the number of samples found in the alignment is retained in a
Bio::Annotation::TagTree object in the annotation collection and is accessible
via:

  ($samples) = $aln->annotation->get_Annotations('Samples');
  say $samples->display_text;
  # or use other relevant TagTree methods to retrieve data

=head1 FEEDBACK

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

=head1 AUTHORS

Chris Fields (cjfields)

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::arp;
use strict;
use base qw(Bio::AlignIO);

use Data::Dumper;
use Bio::Annotation::AnnotationFactory;

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln
 Function: returns the next alignment in the stream.
 Returns : Bio::Align::AlignI object - returns 0 on end of file
            or on error
 Args    : -width => optional argument to specify the width sequence
           will be written (60 chars by default)

See L<Bio::Align::AlignI>

=cut

sub next_aln {
    my $self = shift;
    my $aln = Bio::SimpleAlign->new(-source => 'arp');
    my ($data, $cur_block, $cur_type, $cur_data);
    SCAN:
    while (defined ($data = $self->_readline) ) {
        next if $data =~ m{^\s*$}xms;
        if ($data =~ m{\[{1,2}(\w+)\]{1,2}}xms) {
            $self->{state}->{current_block} = $1;
            next SCAN;
        }
        elsif ($data =~ m{^\s*(\w+)=\s?(\S[^\n]*$)}xms) {
            ($cur_type, $cur_data) = ($1, $2);
            if ($cur_data =~ m{^\s*\{\s*$}) {
                $self->throw("Curly block must be embedded in a named Block")
                    if !exists($self->{state}->{current_block});
                $self->{state}->{in_curly_block} = 1;
                next SCAN;
            }
            $cur_data =~ s{[\"\']}{}g;
            $cur_data =~ s{\s*$}{};
            # per alignment annotation data (i.e. Sample Blocks) or
            # annotation data retained for each alignment?
            $self->{state}->{current_block} eq 'Samples' ?
                push @{$self->{state}->{SampleAnnotation}->{$cur_type}}, $cur_data :
                push @{$self->{state}->{Annotation}->{$cur_type}}, $cur_data;
        }
        elsif ($data =~ m{^\s*\}\s*$}xms) {
            $self->throw("Unmatched bracket in ARP file:\n$data") if
                !exists($self->{state}->{in_curly_block});
            if ($self->{state}->{current_block} eq 'Samples') {;
                my $ac = $self->_process_annotation($aln);
                delete $self->{state}->{SampleAnnotation};
            } else {
                # process other data at a later point
            }
            delete $self->{state}->{blockdata};
            $self->{state}->{in_curly_block} = 0;
            last SCAN;
        }
        else {
            # all other data should be in a curly block and have a block title
            $self->throw("Data found outside of proper block:\n$data") if
                !exists($self->{state}->{current_block}) && !$self->{state}->{in_curly_block};
            # bypass commented stuff (but we may want to process it at a later
            # point, so turn back here)
            next if $data =~ m{^\s*\#}xms;
            if ($self->{state}->{current_block} eq 'Samples') {
                chomp $data;
                # we have two possible ways to deal with sample number, either
                # clone the LocatableSeq (in which case we need to deal with ID
                # duplication), or store as annotation data. I chose the latter
                # route using a Bio::Annotation::TagTree. YMMV - cjfields 10-15-08
                my ($ls, $samples) = $self->_process_sequence($data);
                my $id = $ls->id;
                push @{ $self->{state}->{SampleAnnotation}->{Samples} }, [$id => $samples];
                $aln->add_seq($ls);
            } else {
                # add elsif's for further processing
                #$self->debug('Unmatched data in block '.
                #             $self->{state}->{current_block}.
                #             ":\n$data\n");
                $self->{state}->{blockdata} .= $data;
            }
        }
    }
    # alignments only returned if they contain sequences
    return $aln if $aln->num_sequences;
    return;
}

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in xmfa format
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object

See L<Bio::Align::AlignI>

=cut

sub write_aln {
    my ($self,@aln) = @_;
    $self->throw_not_implemented;
}

################ PRIVATE SUBS ################ 

sub _process_sequence {
    my ($self, $raw) = @_;
    return unless defined $raw;
    $raw =~ s{(?:^\s+|\s+$)}{}g;
    my ($id, $samples, $seq) = split(' ', $raw);
    my $ls = Bio::LocatableSeq->new('-seq'        => $seq,
                                    '-start'      => 1,
                                    '-display_id' => $id,
				    '-alphabet'   => $self->alphabet);
    return($ls, $samples);
}

sub _process_annotation {
    my ($self, $aln) = @_;
    my $coll = Bio::Annotation::Collection->new();
    my $factory = Bio::Annotation::AnnotationFactory->new(-type => 'Bio::Annotation::SimpleValue');
    for my $anntype (qw(SampleAnnotation Annotation)) {
        for my $key (keys %{ $self->{state}->{$anntype} }) {
            if ($key eq 'Title') {
                $aln->description($self->{state}->{$anntype}->{$key}[0]);
            } elsif ($key eq 'Samples') {
                $factory->type('Bio::Annotation::TagTree');
                $coll->add_Annotation($key, $factory->create_object(
                    -value => [$key => $self->{state}->{$anntype}->{$key}]));
                $factory->type('Bio::Annotation::SimpleValue');
            } elsif ($key eq 'SampleName') {
                $aln->id($self->{state}->{$anntype}->{$key}[0]);
            } else {
                $self->throw('Expecting an array reference') unless
                    ref $self->{state}->{$anntype}->{$key} eq 'ARRAY';
                for my $a (@{ $self->{state}->{$anntype}->{$key} }) {
                    $coll->add_Annotation($key, $factory->create_object(
                        -value => $a) );
                }
            }
        }
    }
    #$self->debug("Collection:".Dumper($coll)."\n");
    $aln->annotation($coll);
}

1;
