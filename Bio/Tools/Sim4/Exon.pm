
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

See Bio::Tools::Sim4::Results for a description of the context.     

=head1 DESCRIPTION

This class inherits from Bio::SeqFeature::FeaturePair. Feature1 refers
to the 'exon' on the genomic sequence, so that $exon->start(), $exon->end()
etc will always return what you expect. To get the match to this exon (as
aligned by sim4) as a Generic feature, use
     $esthit = $exon->est_hit()
which is the same as referring to $exon->feature2().

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
Bug fixes and enhancements (doc & code) by Hilmar Lapp <hlapp@gmx.net> or
<hilmar.lapp@pharma.novartis.com>.

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
use Bio::SeqFeature::Generic;

@ISA = qw(Bio::SeqFeature::FeaturePair);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
    my ($self,@args) = @_;
    my $make = $self->SUPER::_initialize(@args);

    # If the caller didn't specify a -feature1 object (and he/she isn't
    # supposed to know about that), we create it here, passing all arguments
    # to it for proper initialization (i.e., we pretend to be a Generic
    # feature to allow for arguments like -start, -end, etc). Apart from that,
    # because calling things like primary_tag() actually refer to the feature1
    # object (as implemented by FeaturePair), we need to have it defined
    # anyway.
    if(! defined($self->feature1())) {
	my $fea1 = Bio::SeqFeature::Generic->new(@args);
	$self->feature1($fea1);
    }
    $self->primary_tag('exon'); # set 
    $self->source_tag('Sim4');
    $self->strand(0) unless defined($self->strand());
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

sub percentage_id {
    my $obj = shift;
    if( @_ ) {
	my $value = shift;
	$obj->{'percentage_id'} = $value;
    }
    return $obj->{'percentage_id'};

}

=head2 est_hit

 Title   : est_hit
 Usage   : $est_feature = $obj->est_hit();
 Function: Returns the EST hit pointing to (i.e., aligned to by Sim4) this
           exon (i.e., genomic region). At present, merely a synonym for
           $obj->feature2().
 Returns : An Bio::SeqFeatureI implementing object.
 Args    : 


=cut

sub est_hit {
    my $self = shift;
    return $self->feature2();
}

1;


