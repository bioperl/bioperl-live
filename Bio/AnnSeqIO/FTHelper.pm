
#
# BioPerl module for Bio::AnnSeqIO::FTHelper
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::AnnSeqIO::FTHelper - Helper class for Embl/Genbank feature tables

=head1 SYNOPSIS

Used by Bio::AnnSeqIO::EMBL to help process the Feature Table

=head1 DESCRIPTION

Represents one particular Feature with the following fields

      key - the key of the feature
      loc - the location string of the feature
      <other fields> - other fields

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


package Bio::AnnSeqIO::FTHelper;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object Exporter);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;
  $self->{'_field'} = {};
# set stuff in self from @args
  return $make; # success - we hope!
}

=head2 from_SeqFeature

 Title   : from_SeqFeature
 Usage   : @fthelperlist = Bio::AnnSeqIO::FTHelper::from_SeqFeature($sf);
 Function: constructor of fthelpers from SeqFeatures
 Returns : an array of FThelpers
 Args    : seq features


=cut

sub from_SeqFeature {
    my ($sf) = shift;
    my @ret;
    my $key;

    # if the parent homogenous flag is set, build things from the
    # sub level
    my $loc;
    if( $sf->_parse->{'parent_homogenous'} == 1 ) {
	$key = $sf->primary_tag();
	$key =~ s/_span//g;
	$loc = "join("; 
	my $touch = 0;
	foreach my $sub ( $sf->sub_SeqFeature() ) {
	    if( $touch == 1 ) {
		$loc .= ",";
	    } else {
		$touch = 1;
	    }

	    $loc .= $sub->start() . ".." . $sub->end();
	}
	$loc .= ")";
	if( $sf->strand == -1 ) {
	    $loc = "complement($loc)";
	} 
    } else {
	$loc = $sf->start() . ".." . $sf->end();
	$key = $sf->primary_tag();
	if( $sf->strand == -1 ) {
	    $loc = "complement($loc)";
	}
	# going into sub features
	foreach my $sub ( $sf->sub_SeqFeature() ) {
	    my $subfth = &Bio::AnnSeqIO::FTHelper::from_SeqFeature($sub);
	    push(@ret,$subfth);
	}
    }

    my $fth = new Bio::AnnSeqIO::FTHelper();

    $fth->loc($loc);
    $fth->key($key);
    foreach my $tag ( $sf->all_tags ) {
	if( !defined $fth->field->{$tag} ) {
	    $fth->field->{$tag} = [];
	}
	foreach my $val ( $sf->each_tag_value($tag) ) {
	    push(@{$fth->field->{$tag}},$val);
	}
    }

    push(@ret,$fth);
    return @ret;

}

=head2 key

 Title   : key
 Usage   : $obj->key($newval)
 Function: 
 Example : 
 Returns : value of key
 Args    : newvalue (optional)


=cut

sub key{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'key'} = $value;
    }
    return $obj->{'key'};

}

=head2 loc

 Title   : loc
 Usage   : $obj->loc($newval)
 Function: 
 Example : 
 Returns : value of loc
 Args    : newvalue (optional)


=cut

sub loc{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'loc'} = $value;
    }
    return $obj->{'loc'};

}


=head2 field

 Title   : field
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub field{
   my ($self) = @_;

   return $self->{'_field'};
}


