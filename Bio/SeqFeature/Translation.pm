
#
# BioPerl module for Bio::SeqFeature::Translation
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Translation - Represents one coding path on genomic DNA

=head1 SYNOPSIS

   # Really a coordinating object with Gene. Read Bio::SeqFeature::Gene
   # documentation

    foreach $tls ( $gene->each_Translation ) {
	print "Translation start at ",$tls->start," end  ",$tls->end,"\n";

	# assumming this gene has been attached to an annseq
	# $pep is a Bio::Seq object
	$pep = $tls->seq();
    }
   

=head1 DESCRIPTION

Describe the object here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Translation;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::SeqFeature::Generic;


@ISA = qw(Bio::SeqFeature::Generic);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;
  $self->primary_tag("Translation");
  $self->source_tag("Bioperl");

  $self->{'_internal_exons'} = [];
# set stuff in self from @args
  return $make; # success - we hope!
}

=head2 seq

 Title   : seq
 Usage   : $seq = $tls->seq()
 Function: Gives translated sequence of gene, as a Bio::Seq object
 Example :
 Returns : 
 Args    :


=cut

sub seq{
   my ($self,@args) = @_;
   my ($str,$seq);

   $seq = $self->first_exon->seq();
   $str = $seq->str();
   foreach my $ex ( $self->_each_internal_exon() ) {
       $seq = $ex->seq();
       $str .= $seq->str();
   }
   $seq = $self->last_exon->seq();
   $str .= $seq->str();

   $seq = Bio::Seq->new(-seq => $str,-type => 'Dna');

   my $out = $seq->translate();
   
   return $out;
}

=head2 _add_internal_exon

 Title   : _add_internal_exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _add_internal_exon{
   my ($self,@args) = @_;
   push(@{$self->{'_internal_exons'}},@args);
}

=head2 _each_internal_exon

 Title   : _each_internal_exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _each_internal_exon{
   my ($self,@args) = @_;
   return @{$self->{'_internal_exons'}};
}

=head2 first_exon

 Title   : first_exon
 Usage   : $obj->first_exon($newval)
 Function: 
 Example : 
 Returns : value of first_exon
 Args    : newvalue (optional)


=cut

sub first_exon{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'first_exon'} = $value;
    }
    return $obj->{'first_exon'};

}

=head2 last_exon

 Title   : last_exon
 Usage   : $obj->last_exon($newval)
 Function: 
 Example : 
 Returns : value of last_exon
 Args    : newvalue (optional)


=cut

sub last_exon{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'last_exon'} = $value;
    }
    return $obj->{'last_exon'};

}

=head2 _attach_seq

 Title   : _attach_seq
 Usage   : used internally to attach sequences to seqfeatures
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _attach_seq{
   my ($self,$seq) = @_;

   $self->first_exon->attach_seq($seq);
   $self->last_exon->attach_seq($seq);
}

1;
