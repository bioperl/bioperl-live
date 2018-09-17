#
# BioPerl module for Bio::Tools::isPcr
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sheldon McKay <mckays@cshl.edu>
#
# Copyright Sheldon McKay
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::isPcr - Parse isPcr output and make features

=head1 SYNOPSIS

    # A simple annotation pipeline wrapper for isPcr data
    # assuming isPcr data is already generated in file seq1.isPcr
    # and sequence data is in fasta format in file called seq1.fa

    # Note: this parser is meant for the default fasta output from
    # isPcr.  bed and psl output formats are not supported.

    use Bio::Tools::IsPcr;
    use Bio::SeqIO;
    my $parser = Bio::Tools::isPcr->new(-file => 'seq1.isPcr');
    my $seqio = Bio::SeqIO->new(-format => 'fasta', -file => 'seq1.fa');
    my $seq = $seqio->next_seq || die("cannot get a seq object from SeqIO");

    while( my $feat = $parser->next_feature ) {
	# add isPcr annotation to a sequence
	$seq->add_SeqFeature($feat);
    }
    my $seqout = Bio::SeqIO->new(-format => 'embl');
    $seqout->write_seq($seq);


=head1 DESCRIPTION

This object serves as a parser for isPcr data (in the default fasta
format), creating a Bio::SeqFeatureI for each isPcr hit.  
These can be processed or added as annotation to an existing
Bio::SeqI object for the purposes of automated annotation.

This module is adapted from the Bio::Tools::EPCR module
written by Jason Stajich (jason-at-bioperl.org).

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

=head1 AUTHOR - Sheldon McKay

Email mckays@cshl.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::isPcr;
use strict;

use Bio::SeqIO;
use Bio::SeqFeature::Generic;


use base qw(Bio::Root::Root);


=head2 new

 Title   : new
 Usage   : my $ispcr = Bio::Tools::isPcr->new( -file => $file,
					      -primary => $fprimary, 
					      -source => $fsource,
					      -groupclass => $fgroupclass);

 Function: Initializes a new isPcr parser
 Returns : Bio::Tools::isPcr
 Args    : -fh   => filehandle
           OR
           -file => filename

           -primary => a string to be used as the common value for
                       each features '-primary' tag.  Defaults to
                       the sequence ontology term 'PCR_product'.  
                       (This in turn maps to the GFF 'type'
                       tag (aka 'method')).

            -source => a string to be used as the common value for
                       each features '-source' tag.  Defaults to
                       'isPcr'. (This in turn maps to the GFF 'source'
                       tag)

            -groupclass => a string to be used as the name of the tag
                           which will hold the sts marker namefirst
                           attribute.  Defaults to 'name'.

=cut


sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($primary,$source,$groupclass) = 
      $self->_rearrange([qw/PRIMARY SOURCE GROUPCLASS/],@args);
  $self->primary(defined $primary ? $primary : 'PCR_product');
  $self->source(defined $source ? $source : 'isPcr');
  $self->groupclass(defined $groupclass ? $groupclass : 'name');

  # default output for isPcr is fasta format
  $self->{io} = Bio::SeqIO->new(-format => 'fasta', @args);

  return $self;
}

=head2 next_feature

 Title   : next_feature
 Usage   : $seqfeature = $obj->next_feature();
 Function: Returns the next feature available in the analysis result, or
           undef if there are no more features.
 Example :
 Returns : A Bio::SeqFeatureI implementing object, or undef if there are no
           more features.
 Args    : none    

=cut

sub next_feature {
    my ($self) = @_;
    my $result = $self->{io}->next_seq;
    return unless defined $result;
    
    my ($seqname,$location)    = split ':', $result->primary_id;
    my ($pcrname,$left,$right) = split /\s+/, $result->desc;
    my ($start,$strand,$end)   = $location =~ /^(\d+)([-+])(\d+)$/;
    my $amplicon = $result->seq;

    # if there are multiple hits, increment the name for
    # the groupclass
    if (++$self->{seen}->{$pcrname} > 1) {
	$pcrname .= "\.$self->{seen}->{$pcrname}";
    }

    my $tags = {
	$self->groupclass => $pcrname,
	amplicon          => $amplicon,
	left_primer       => $left,
	right_primer       => $right
	};

    my $markerfeature = Bio::SeqFeature::Generic->new( 
						       '-start'   => $start,
						       '-end'     => $end,
						       '-strand'  => $strand,
						       '-source'  => $self->source,
						       '-primary' => $self->primary,
						       '-seq_id'  => $seqname,
						       '-tag'     => $tags
						       );

    return $markerfeature;
}

=head2 source

 Title   : source
 Usage   : $obj->source($newval)
 Function: 
 Example : 
 Returns : value of source (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub source{
    my $self = shift;
    return $self->{'_source'} = shift if @_;
    return $self->{'_source'};
}

=head2 primary

 Title   : primary
 Usage   : $obj->primary($newval)
 Function: 
 Example : 
 Returns : value of primary (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub primary{
    my $self = shift;
    return $self->{'_primary'} = shift if @_;
    return $self->{'_primary'};
}

=head2 groupclass

 Title   : groupclass
 Usage   : $obj->groupclass($newval)
 Function: 
 Example : 
 Returns : value of groupclass (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub groupclass{
    my $self = shift;

    return $self->{'_groupclass'} = shift if @_;
    return $self->{'_groupclass'};
}

1;
