
#
# BioPerl module for Bio::Tools::Sim4::Exon
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Sim4::Exon - A single 

=head1 SYNOPSIS

Give standard usage here

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


package Bio::Tools::Sim4::Exon;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::SeqFeature::FeaturePair;


@ISA = qw(Bio::SeqFeature::FeaturePair);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  $self->primary_tag('exon'); # set 
  $self->source_tag('Sim4');
  $self->strand(0);
# set stuff in self from @args
  return $make; # success - we hope!
}

#
# Everything else is just inherited from SeqFeature::Generic. Cool.
#

=head2 percentage_id

 Title   : percentage_id
 Usage   : $obj->percentage_id($newval)
 Function: 
 Returns : value of percentage_id
 Args    : newvalue (optional)


=cut

sub percentage_id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'percentage_id'} = $value;
    }
    return $obj->{'percentage_id'};

}

1;


