#
# BioPerl module for Bio::Tools::Phylo::Gumby
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

Bio::Tools::Phylo::Gumby - Parses output from gumby

=head1 SYNOPSIS

  #!/usr/bin/perl -Tw
  use strict;

  use Bio::Tools::Phylo::Gumby;

  my $parser = Bio::Tools::Phylo::Gumby->new(-file => "out.align");
  my @features = $parser->next_result();

=head1 DESCRIPTION

This module is used to parse the output from 'gumby' by Shyam Prabhakar. You
can get details here: http://pga.lbl.gov/gumby/

It works on the .align files produced.

The result is a list of feature objects.

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

=head1 AUTHOR - Sendu Bala

Email bix@sendu.me.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::Phylo::Gumby;
use strict;

use Bio::SeqFeature::Annotated;
use Bio::Annotation::SimpleValue;

use base qw(Bio::Root::Root Bio::Root::IO);


=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Phylo::Gumby->new();
 Function: Builds a new Bio::Tools::Phylo::Gumby object
 Returns : Bio::Tools::Phylo::Gumby
 Args    : -file (or -fh) should contain the contents of a gumby .align output
           file

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    
    $self->_initialize_io(@args);
    
    return $self;
}

=head2 next_result

 Title   : next_result
 Usage   : $result = $obj->next_result();
 Function: Returns the next set of results available from the input, or undef if
           there are no more results.
 Returns : list of Bio::SeqFeature::Annotated (one per sequence). Features are
           annotated with tags for pvalue and 'kind' (holding 'all', 'exon', or
           'nonexon').

           NB: Gumby ignores sequence coordinates in input alignments, treating
           each sequence as if it started at position 1. If you're running this
           manually (ie. not via the Bio::Tools::Run::Phylo::Gumby) you will
           have to adjust the coordinates to match up with your input alignment
           and sequences.
 Args    : none

=cut

sub next_result {
    my ($self) = @_;
    
    my $line = $self->_readline || return;
    while ($line !~ /^start/) {
        $line = $self->_readline || return;
        
        if ($line =~ /^(all|exon|nonexon):/) {
            $self->{_kind} = $1;
        }
    }
    
    my ($score, $pvalue) = $line =~ /score (\d+), pvalue (\S+)/;
    
    my @feats;
    while ($line = $self->_readline) {
        $line =~ /^$/ && last;
        $line || last;
        my ($seq_id, $start, $end) = split(/\s+/, $line);
        
        my $feature = Bio::SeqFeature::Annotated->new(-seq_id => $seq_id,
                                                      -start  => $start,
                                                      -end    => $end,
                                                      -score  => $score,
                                                      -strand => 1,
                                                      -source => 'gumby');
        my $sv = Bio::Annotation::SimpleValue->new(-tagname => 'pvalue', -value => $pvalue);
        $feature->annotation->add_Annotation($sv);
        $sv = Bio::Annotation::SimpleValue->new(-tagname => 'kind', -value => $self->{_kind});
        $feature->annotation->add_Annotation($sv);
        $sv = Bio::Annotation::SimpleValue->new(-tagname => 'predicted', -value => 1);
        $feature->annotation->add_Annotation($sv);
        
        push(@feats, $feature);
    }
    
    return @feats;
}

1;
