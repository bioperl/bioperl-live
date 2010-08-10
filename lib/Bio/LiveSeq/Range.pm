#
# bioperl module for Bio::LiveSeq::Range
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

Bio::LiveSeq::Range - Range abstract class for LiveSeq

=head1 SYNOPSIS

  # documentation needed

=head1 DESCRIPTION

This is used as parent for exon and intron classes.

=head1 AUTHOR - Joseph A.L. Insana

Email:  Insana@ebi.ac.uk, jinsana@gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::LiveSeq::Range;
use strict;
use base qw(Bio::LiveSeq::SeqI);

=head2 new

  Title   : new
  Usage   : $range1 = Bio::LiveSeq::Range->new(-seq => $obj_ref,
                                               -start => $beginlabel,
                                               -end => $endlabel, -strand => 1);

  Function: generates a new Bio::LiveSeq::Range
  Returns : reference to a new object of class Range
  Errorcode -1
  Args    : two labels, an obj_ref and an integer
            strand 1=forward strand, strand -1=reverse strand
            if strand not specified, it defaults to 1
            the -seq argument must point to the underlying DNA LiveSeq object

=cut

sub new {
  my ($thing, %args) = @_;
  my $class = ref($thing) || $thing;
  my ($obj,%range);

  my ($seq,$start,$end,$strand)=($args{-seq},$args{-start},$args{-end},$args{-strand});

  $obj = \%range;
  $obj = bless $obj, $class;

  unless ($seq->valid($start)) {
    $obj->warn("$class not initialised because start label not valid");
    return (-1);
  }
  unless ($seq->valid($end)) {
    $obj->warn("$class not initialised because end label not valid");
    return (-1);
  }
  unless (defined $strand) {
    $strand = 1;
  }
  if (($strand != 1)&&($strand != -1)) {
    $obj->warn("$class not initialised because strand identifier not valid. Use 1 (forward strand) or -1 (reverse strand).");
    return (-1);
  }
  if ($start eq $end) {
    $obj->warn("$class reports: start and end label are the same....");
  } else {
    unless ($seq->follows($start,$end,$strand)==1) {
      $obj->warn("Fatal: end label $end doesn't follow start label $start for strand $strand!");
      return (-1);
    }
  }
  #if ($strand == 1) {
  #  unless ($seq->is_downstream($start,$end)==1) {
  #    croak "Fatal: end label not downstream of start label for forward strand!";
  #  }
  #} else {
  #  unless ($seq->is_upstream($start,$end)==1) {
  #    croak "Fatal: end label not upstream of start label for reverse strand!";
  #  }
  #}
  $obj->{'seq'}=$seq;
  $obj->{'start'}=$start;
  $obj->{'end'}=$end;
  $obj->{'strand'}=$strand;
  return $obj;
}

=head2 valid

  Title   : valid
  Usage   : $boolean = $obj->valid($label)
  Function: tests if a label exists AND is part of the object
  Returns : boolean
  Args    : label

=cut

1;
