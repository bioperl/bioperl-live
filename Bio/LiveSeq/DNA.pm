#
# bioperl module for Bio::LiveSeq::DNA
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Joseph Insana <insana@ebi.ac.uk> <jinsana@gmx.net>
#
# Copyright Joseph Insana
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::LiveSeq::DNA - DNA object for LiveSeq

=head1 SYNOPSIS

  # documentation needed

=head1 DESCRIPTION

This holds the DNA sequence (or the RNA in the case of cDNA entries)
and is accessed by exons, genes, transcripts... objects

=head1 AUTHOR - Joseph A.L. Insana

Email:  Insana@ebi.ac.uk, jinsana@gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::LiveSeq::DNA;

use strict;
use base qw(Bio::LiveSeq::SeqI);

=head2 new

  Title   : new
  Usage   : $dna = Bio::LiveSeq::DNA->new(-seq => "atcgaccaatggacctca",
					  -offset => 3 );

  Function: generates a new Bio::LiveSeq::DNA
  Returns : reference to a new object of class DNA
  Errorcode -1
  Args    : a string
        AND an optional offset to create nucleotide labels (default is 1, i.e.
            starting the count of labels from "1") -> do not bother using it ->
            it could be used by alternative loaders !EMBL format
  NOTE    : strand of DNA is set to 1 by default

=cut

sub new {
  my ($thing, %args) = @_;
  my $class = ref($thing) || $thing;
  my (%empty,$obj);

  if ($args{-seq}) {
    $obj = $thing->string2chain($args{-seq},$args{-offset}); # inherited from ChainI
    $obj = bless $obj, $class;
  } else {
    $obj=\%empty;
    $obj = bless $obj, $class;
    $obj->throw("$class not initialized properly");
  }

  $obj->{'alphabet'}='dna'; # set alphabet default
  $obj->{'strand'}=1; # set strand default = 1
  $obj->{'seq'}=$obj; # set seq field to itself

  return $obj;
}

# START method
# it has to be redefined here because default from SeqI accesses field "start"
sub start {
  my $self = shift;
  return $self->{'begin'}; # the chain's start is called begin
}

# it is overridden to provide faster output
sub length {
  my $self=shift;
  return $self->chain_length();
}

# it is overridden to provide MUCH faster output
sub valid {
  my $self=shift(@_);
  return $self->label_exists(@_);
}

1;
