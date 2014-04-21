# $Id: Gumby.pm,v 1.2 2007/06/14 18:01:52 nathan Exp $
#
# BioPerl module for Bio::Tools::Phylo::Gerp
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

Bio::Tools::Phylo::Gerp - Parses output from GERP

=head1 SYNOPSIS

  use strict;

  use Bio::Tools::Phylo::Gerp;

  my $parser = Bio::Tools::Phylo::Gerp->new(-file => "alignment.rates.elems");

  while (my $feat = $parser->next_result) {
    my $start = $feat->start;
    my $end = $feat->end;
    my $rs_score = $feat->score;
    my $p_value = ($feat->annotation->get_Annotations('p-value'))[0]->value;
  }

=head1 DESCRIPTION

This module is used to parse the output from 'GERP' (v2) by Eugene Davydov
(originally Gregory M. Cooper et al.). You can get details here:
http://mendel.stanford.edu/sidowlab/

It works on the .elems files produced by gerpelem.

Each result is a Bio::SeqFeature::Annotated representing a single constrained
element.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list. Your participation is much appreciated.

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

package Bio::Tools::Phylo::Gerp;
use strict;

use Bio::SeqFeature::Generic;
use Bio::Annotation::SimpleValue;

use base qw(Bio::Root::Root Bio::Root::IO);


=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Phylo::Gerp->new();
 Function: Builds a new Bio::Tools::Phylo::Gerp object
 Returns : Bio::Tools::Phylo::Gerp
 Args    : -file (or -fh) should contain the contents of a gerpelem .elems file

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
 Function: Returns the next result available from the input, or undef if there
           are no more results.
 Returns : Bio::SeqFeature::Annotated object. Features are annotated with a tag
           for 'pvalue', and a 'predicted' tag. They have no sequence id unless
           the input GERP file is non-standard, with the seq id as the 6th
           column.

           NB: feature coordinates are alignment columns of the alignment
           used to create the result file.
 Args    : none

=cut

sub next_result {
    my ($self) = @_;
    
    my $line = $self->_readline || return;
    
    while ($line !~ /^\d+\s+\d+\s+\d+\s+\S+\s+\S+\s*(?:\S+\s*)?$/) {
        $line = $self->_readline || return;
    }
    
    #start   end     length  RS-score        p-value
    # code elsewhere adds seq_id on the end (not valid GERP), so we capture that
    # if present
    my ($start, $end, undef, $rs_score, $p_value, $seq_id) = split(/\s+/, $line);
    my $feat = Bio::SeqFeature::Generic->new(
        $seq_id ? (-seq_id => $seq_id) : (),
        -start        => $start, 
        -end          => $end,
        -strand       => 1,
        -score        => $rs_score,
        #-type         => 'conserved_region', ***causes 740x increase in SeqFeatureDB storage requirments!
        -source       => 'GERP');
    
    my $sv = Bio::Annotation::SimpleValue->new(-tagname => 'predicted', -value => 1);
    $feat->annotation->add_Annotation($sv);
    $sv = Bio::Annotation::SimpleValue->new(-tagname => 'pvalue', -value => $p_value);
    $feat->annotation->add_Annotation($sv);
    
    return $feat;
}

1;
