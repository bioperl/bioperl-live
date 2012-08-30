#
# BioPerl module for Bio::Tools::HMMER::Domain
#
# Please direct questions and support issues to <bioperl-l@bioperl.org>
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::HMMER::Domain - One particular domain hit from HMMER

=head1 SYNOPSIS

Read the Bio::Tools::HMMER::Results docs

=head1 DESCRIPTION

A particular domain score. We reuse the Homol SeqFeature system
here, so this inherits off Homol SeqFeature. As this code
originally came from a separate project, there are some backward
compatibility stuff provided to keep this working with old code.

Don't forget this inherits off Bio::SeqFeature, so all your usual
nice start/end/score stuff is ready for use.

=head1 CONTACT

Ewan Birney, birney@ebi.ac.uk

=head1 CONTRIBUTORS

Jason Stajich, jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#'
package Bio::Tools::HMMER::Domain;

use Bio::SeqFeature::Generic;
use strict;


use base qw(Bio::SeqFeature::FeaturePair);

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->{'alignlines'} = [];

  my $hmmf1 = Bio::SeqFeature::Generic->new(@args);
  my $hmmf2 = Bio::SeqFeature::Generic->new(@args);

  $self->feature1($hmmf1);
  $self->feature2($hmmf2);

  return $self;
}

=head2 add_alignment_line

 Title   : add_alignment_line
 Usage   : $domain->add_alignment_line($line_from_hmmer_output);
 Function: add an alignment line to this Domain object
 Returns : Nothing
 Args    : scalar

 Adds an alignment line, mainly for storing the HMMER alignments
as flat text which can be reguritated. You're right. This is *not
nice* and not the right way to do it.  C'est la vie.

=cut

sub add_alignment_line {
    my $self = shift;
    my $line = shift;
    push(@{$self->{'alignlines'}},$line);
}

=head2 each_alignment_line

 Title   : each_alignment_line
 Usage   : foreach $line ( $domain->each_alignment_line )
 Function: reguritates the alignment lines as they were fed in.
           only useful realistically for printing.
 Example :
 Returns :
 Args    : None


=cut

sub each_alignment_line {
    my $self = shift;
    return @{$self->{'alignlines'}};
}

=head2 get_nse

 Title   : get_nse
 Usage   : $domain->get_nse()
 Function: Provides a seqname/start-end format, useful
           for unique keys. nse stands for name-start-end
           It is used a lot in Pfam
 Example :
 Returns : A string
 Args    : Optional separator 1 and separator 2 (default / and -)


=cut



sub get_nse {
    my $self = shift;
    my $sep1 = shift;
    my $sep2 = shift;

    if( !defined $sep2 ) {
	$sep2 = "-";
    }
    if( !defined $sep1 ) {
	$sep1 = "/";
    }

    return sprintf("%s%s%d%s%d",$self->seq_id,$sep1,$self->start,$sep2,$self->end);
}


#  =head2 start_seq

#   Title   : start_seq
#   Usage   : Backward compatibility with old HMMER modules.
#             should use $domain->start
#   Function:
#   Example :
#   Returns :
#   Args    :

#  =cut

sub start_seq {
    my $self = shift;
    my $start = shift;

    $self->warn("Using old domain->start_seq. Should use domain->start");
    return $self->start($start);
}

#  =head2 end_seq

#   Title   : end_seq
#   Usage   : Backward compatibility with old HMMER modules.
#             should use $domain->end
#   Function:
#   Example :
#   Returns :
#   Args    :

#  =cut

sub end_seq {
    my $self = shift;
    my $end = shift;

    $self->warn("Using old domain->end_seq. Should use domain->end");
    return $self->end($end);
}

#  =head2 start_hmm

#   Title   : start_hmm
#   Usage   : Backward compatibility with old HMMER modules, and
#             for convience. Equivalent to $self->homol_SeqFeature->start
#   Function:
#   Example :
#   Returns :
#   Args    :

#  =cut

sub start_hmm {
    my $self = shift;
    my $start = shift;
    $self->warn("Using old domain->start_hmm. Should use domain->hstart");
    return $self->hstart($start);
}

#  =head2 end_hmm

#   Title   : end_hmm
#   Usage   : Backward compatibility with old HMMER modules, and
#             for convience. Equivalent to $self->homol_SeqFeature->start
#   Function:
#   Example :
#   Returns :
#   Args    :

#  =cut

sub end_hmm {
    my $self = shift;
    my $end = shift;

    $self->warn("Using old domain->end_hmm. Should use domain->hend");
    return $self->hend($end);
}

=head2 hmmacc

 Title   : hmmacc
 Usage   : $domain->hmmacc($newacc)
 Function: set get for HMM accession number. This is placed in the homol
           feature of the HMM
 Example :
 Returns :
 Args    :


=cut

sub hmmacc{
   my ($self,$acc) = @_;
   if( defined $acc ) {
       $self->feature2->add_tag_value('accession',$acc);
   }
   my @vals = $self->feature2->each_tag_value('accession');
   return shift @vals;
}

=head2 hmmname

 Title   : hmmname
 Usage   : $domain->hmmname($newname)
 Function: set get for HMM accession number. This is placed in the homol
           feature of the HMM
 Example :
 Returns :
 Args    :

=cut

sub hmmname {
    return shift->hseq_id(@_);
}

=head2 bits

 Title   : bits
 Usage   :
 Function: backward compatibility. Same as score
 Example :
 Returns :
 Args    :

=cut

sub bits{
   return shift->score(@_);
}

=head2 evalue

 Title   : evalue
 Usage   :
 Function: $domain->evalue($value);
 Example :
 Returns :
 Args    :

=cut

sub evalue{
    return shift->_tag_value('evalue',@_);
}

=head2 seqbits

 Title   : seqbits
 Usage   :
 Function: $domain->seqbits($value);
 Example :
 Returns :
 Args    :

=cut

sub seqbits {
    return shift->_tag_value('seqbits',@_);
}

=head2 seq_range

 Title   : seq_range
 Usage   :
 Function: Throws an exception to catch scripts which need to upgrade
 Example :
 Returns :
 Args    :

=cut

sub seq_range{
   my ($self,@args) = @_;

   $self->throw("You have accessed an old method. Please recode your script to the new bioperl HMMER module");
}

=head2 hmm_range

 Title   : hmm_range
 Usage   :
 Function: Throws an exception to catch scripts which need to upgrade
 Example :
 Returns :
 Args    :


=cut

sub hmm_range{
   my ($self,@args) = @_;

   $self->throw("You have accessed an old method. Please recode your script to the new bioperl HMMER module");
}

1;  # says use was ok
__END__
