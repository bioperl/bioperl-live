#
# Module for Bio::PhyloNetwork::muVector
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Gabriel Cardona <gabriel(dot)cardona(at)uib(dot)es>
#
# Copyright Gabriel Cardona
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PhyloNetwork::muVector - Module to compute with vectors of arbitrary
dimension

=head1 SYNOPSIS

 use strict;
 use warnings;

 use Bio::PhyloNetwork::muVector;

 my $vec1=Bio::PhyloNetwork::muVector->new(4);
 my $vec2=Bio::PhyloNetwork::muVector->new([1,2,3,4]);
 my $vec3=Bio::PhyloNetwork::muVector->new([10,20,30,40]);

 my $vec4=$vec3-10*$vec2;
 if (($vec4 cmp $vec1) == 0) {
   print "$vec4 is zero\n";
 }

 my $vec5=Bio::PhyloNetwork::muVector->new([8,2,2,4]);
 my $vec6=Bio::PhyloNetwork::muVector->new([1,2,3,4]);

 print "Test poset $vec5 > $vec6: ".$vec5->geq_poset($vec6)."\n";
 print "Test lex $vec5 > $vec6: ".($vec5 cmp $vec6)."\n";

=head1 DESCRIPTION

This is a module to work with vectors. It creates
vectors of arbitrary length, defines its basic arithmetic operations,
its lexicographic ordering and the natural structure of poset.

=head1 AUTHOR

Gabriel Cardona, gabriel(dot)cardona(at)uib(dot)es

=head1 APPENDIX

The rest of the documentation details each of the object methods.

=cut

package Bio::PhyloNetwork::muVector;

use strict;
use warnings;

use base qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : my $mu = new Bio::PhyloNetwork::muVector();
 Function: Creates a new Bio::PhyloNetwork::muVector object
 Returns : Bio::PhyloNetwork::muVector
 Args    : integer or (reference to) an array

If given an integer as argument, returns a Bio::PhyloNetwork::muVector
object with dimension the integer given and initialized to zero.
If it is an anonimous array, then the vector is initialized with the values
in the array and with the corresponding dimension.

=cut

sub new {
  my ($pkg,$cont)=@_;
  my $self=$pkg->SUPER::new();
  my @arr=();
  if (!ref($cont)) {
    #$cont is a number; initialize to a zero-vector
    for (my $i=0; $i < $cont; $i++) {
      $arr[$i]=0;
    }
    $self->{arr}=\@arr;
  } else {
    #$cont points to an array
    @arr=@{$cont};
  }
  $self->{dim}=scalar @arr;
  $self->{arr}=\@arr;
  bless($self,$pkg);
  return $self;
}

sub dim {
  return shift->{dim}
}

use overload
  "+" => \&add,
  "-" => \&substract,
  "*" => \&scalarproduct,
  "<=>" => \&comparelex,
  "cmp" => \&comparelex,
  '""' => \&display,
  '@{}' => \&as_array;

sub as_array {
  return shift->{arr};
}

=head2 display

 Title   : display
 Usage   : my $str=$mu->display()
 Function: returns an string displaying its contents
 Returns : string
 Args    : none

This function is also overloaded to the "" operator.

=cut

sub display {
  my ($self)=@_;
  my @arr=@{$self->{arr}};
  return "(@arr)";
}

=head2 add

 Title   : add
 Usage   : $mu->add($mu2)
 Function: returns the sum of $mu and $mu2
 Returns : Bio::PhyloNetwork::muVector
 Args    : Bio::PhyloNetwork::muVector

This function is also overloaded to the + operator.

=cut

sub add {
  my ($v1,$v2)=@_;

  $v1->throw("Vectors not the same size") unless ($v1->{dim} == $v2->{dim});
  my $dim=$v1->{dim};
  my @sum=();
  for (my $i=0; $i<$dim; $i++) {
    $sum[$i]=$v1->[$i]+$v2->[$i];
  }
  my $result=Bio::PhyloNetwork::muVector->new(\@sum);
  return $result;
}

=head2 substract

 Title   : substract
 Usage   : $mu->substract($mu2)
 Function: returns the difference of $mu and $mu2
 Returns : Bio::PhyloNetwork::muVector
 Args    : Bio::PhyloNetwork::muVector

This function is also overloaded to the - operator.

=cut

sub substract {
  my ($v1,$v2)=@_;

  $v1->throw("Vectors not the same size") unless ($v1->{dim} == $v2->{dim});
  my $dim=$v1->{dim};
  my @sum=();
  for (my $i=0; $i<$dim; $i++) {
    $sum[$i]=$v1->{arr}->[$i]-$v2->{arr}->[$i];
  }
  my $result=Bio::PhyloNetwork::muVector->new(\@sum);
  return $result;
}

=head2 scalarproduct

 Title   : scalarproduct
 Usage   : $mu->scalarproduct($ct)
 Function: returns the scalar product of $ct and $mu
 Returns : Bio::PhyloNetwork::muVector
 Args    : scalar

This function is also overloaded to the * operator.

=cut

sub scalarproduct {
  my ($v1,$num,$swapped)=@_;

  my $dim=$v1->{dim};
  my @sum=();
  for (my $i=0; $i<$dim; $i++) {
    $sum[$i]=$num*$v1->{arr}->[$i];
  }
  my $result=Bio::PhyloNetwork::muVector->new(\@sum);
  return $result;
  return $result;
}

=head2 comparelex

 Title   : comparelex
 Usage   : $mu1->comparelex($mu2)
 Function: compares $mu and $mu2 w.r.t. the lexicographic ordering
 Returns : scalar (-1 if $mu1<$mu2, 0 if $mu1=$mu2, 1 if $mu1>$mu2)
 Args    : Bio::PhyloNetwork::muVector

This function is also overloaded to the E<lt>=E<gt> and cmp operator.

=cut

sub comparelex {
  my ($v1,$v2)=@_;

  $v1->throw("Vectors not the same size") unless ($v1->{dim} == $v2->{dim});
  my $dim=$v1->{dim};
  for (my $i=0; $i<$dim; $i++) {
    return -1 if $v1->{arr}->[$i] < $v2->{arr}->[$i];
    return 1 if $v1->{arr}->[$i] > $v2->{arr}->[$i];
  }
  return 0;
}

=head2 geq_poset

 Title   : geq_poset
 Usage   : $mu1->geq_poset($mu2)
 Function: compares $mu and $mu2 w.r.t. the natural partial ordering
 Returns : boolean (1 if $mu >= $mu2, 0 otherwise)
 Args    : Bio::PhyloNetwork::muVector

=cut

sub geq_poset {
  my ($v1,$v2)=@_;

  $v1->throw("Vectors not the same size") unless ($v1->{dim} == $v2->{dim});
  my $dim=$v1->{dim};
  for (my $i=0; $i<$dim; $i++) {
    return 0 unless $v1->[$i] >= $v2->[$i];
  }
  return 1;
}

=head2 is_positive

 Title   : is_positive
 Usage   : $mu->is_positive()
 Function: tests if all components of $mu are positive (or zero)
 Returns : boolean
 Args    : none

=cut

sub is_positive {
  my ($v1)=@_;

  my $dim=$v1->{dim};
  for (my $i=0; $i<$dim; $i++) {
    return 0 unless $v1->[$i] >= 0;
  }
  return 1;
}

=head2 hamming

 Title   : hamming
 Usage   : $mu1->hamming($mu2)
 Function: returns the Hamming distance between $mu1 and $mu2
 Returns : scalar
 Args    : Bio::PhyloNetwork::muVector

=cut

sub hamming {
  my ($v1,$v2)=@_;

  $v1->throw("Vectors not the same size") unless ($v1->{dim} == $v2->{dim});
  my $dim=$v1->{dim};
  my $w=0;
  for (my $i=0; $i<$dim; $i++) {
    $w++ unless $v1->[$i] == $v2->[$i];
  }
  return $w;
}

=head2 manhattan

 Title   : manhattan
 Usage   : $mu1->manhattan($mu2)
 Function: returns the Manhattan distance between $mu1 and $mu2
 Returns : scalar
 Args    : Bio::PhyloNetwork::muVector

=cut

sub manhattan {
  my ($v1,$v2)=@_;

  $v1->throw("Vectors not the same size") unless ($v1->{dim} == $v2->{dim});
  my $dim=$v1->{dim};
  my $w=0;
  for (my $i=0; $i<$dim; $i++) {
    $w+= abs($v1->[$i] - $v2->[$i]);
  }
  return $w;
}

1;
