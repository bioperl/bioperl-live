# $Id$
#
# bioperl module for Bio::LiveSeq::DNA
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


=head1 DESCRIPTION

This holds the DNA sequence (or the RNA in the case of cDNA entries)
and is accessed by exons, genes, transcripts... objects

=head1 AUTHOR - Joseph A.L. Insana

Email:  Insana@ebi.ac.uk, jinsana@gmx.net

Address: 

     EMBL Outstation, European Bioinformatics Institute
     Wellcome Trust Genome Campus, Hinxton
     Cambs. CB10 1SD, United Kingdom 

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::LiveSeq::DNA;
$VERSION=1.4;

# Version history:
# Mon Mar 20 19:21:22 GMT 2000 v.1.0 begun
# Tue Mar 21 14:20:30 GMT 2000 v.1.1 new() is now here, not inherited
# Wed Mar 22 19:43:20 GMT 2000 v.1.2 length override
# Thu Jun 22 20:02:39 BST 2000 v 1.3 valid() from SeqI now moved here, as override
# Wed Mar 28 17:01:59 BST 2001 v 1.4 changed croaks into throw

use strict;
use vars qw($VERSION @ISA);
use Bio::LiveSeq::SeqI 3.2; # uses SeqI, inherits from it
@ISA=qw(Bio::LiveSeq::SeqI);

=head1 new

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

  $obj->{'moltype'}='dna'; # set moltype default
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
