
#
# BioPerl module for Bio::Tools::Prediction::Exon
#
# Cared for by Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Prediction::Exon

=head1 SYNOPSIS

See documentation of methods.

=head1 DESCRIPTION

A feature representing a predicted exon. This class actually
inherits off Bio::SeqFeature::Generic and therefore has all that
functionality, plus a few methods supporting predicted features, like various
scores and a significance. Even though these were inspired by GenScan results,
at least a subset should be generally useable for exon prediction results.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Hilmar Lapp

Email hlapp@gmx.net or hilmar.lapp@pharma.novartis.com

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Prediction::Exon;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::SeqFeature::Generic;


@ISA = qw(Bio::SeqFeature::Generic);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
    my($self,@args) = @_;
    
    my $make = $self->SUPER::_initialize(@args);

    my ($primary) =
	$self->_rearrange([qw(PRIMARY
			      )],@args);

    $primary = 'predicted_exon' unless $primary;
    $self->primary_tag($primary);
    # set stuff in self from @args
    return $make; # success - we hope!
}

#
# Everything else is just inherited from SeqFeature::Generic.
#

=head2 _tag_value

 Title   : _tag_value
 Usage   : 
 Function: 
 Returns : 
 Args    : 


=cut

sub _tag_value {
    my ($self, $tag, $value) = @_;

    if(defined($value) || (! $self->has_tag($tag))) {
	$self->remove_tag($tag) if($self->has_tag($tag));
	$self->add_tag_value($tag, $value);
    }
    return ($self->each_tag_value($tag))[0];
}

=head2 significance

 Title   : significance
 Usage   : $evalue = $obj->significance();
           $obj->significance($evalue);
 Function: 
 Returns : 
 Args    : 


=cut

sub significance {
    my ($self, $value) = @_;

    return $self->_tag_value('signif', $value);
}

=head2 start_signal_score

 Title   : start_signal_score
 Usage   : $sc = $obj->start_signal_score();
           $obj->start_signal_score($evalue);
 Function: Get/Set a score for the exon start signal (acceptor splice site
           or initiation signal).
 Returns : 
 Args    : 


=cut

sub start_signal_score {
    my ($self, $value) = @_;

    return $self->_tag_value('AccScore', $value);
}

=head2 end_signal_score

 Title   : end_signal_score
 Usage   : $sc = $obj->end_signal_score();
           $obj->end_signal_score($evalue);
 Function: Get/Set a score for the exon end signal (donor splice site
           or termination signal).
 Returns : 
 Args    : 


=cut

sub end_signal_score {
    my ($self, $value) = @_;

    return $self->_tag_value('DonScore', $value);
}

=head2 coding_signal_score

 Title   : coding_signal_score
 Usage   : $sc = $obj->coding_signal_score();
           $obj->coding_signal_score($evalue);
 Function: Get/Set a score for the exon coding signal (e.g., coding potential).
 Returns : 
 Args    : 


=cut

sub coding_signal_score {
    my ($self, $value) = @_;

    return $self->_tag_value('CodScore', $value);
}

1;

