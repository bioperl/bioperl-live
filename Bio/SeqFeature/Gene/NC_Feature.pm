#
# BioPerl module for Bio::SeqFeature::Gene::NC_Feature.pm
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by David Block <dblock@gene.pbi.nrc.ca>
#
# Copyright David Block
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Gene::NC_Feature.pm - superclass for non-coding features

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - David Block

Email dblock@gnf.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::SeqFeature::Gene::NC_Feature;
use strict;

# Object preamble - inherits from Bio::Root::Root


use base qw(Bio::SeqFeature::Generic);

sub new {
    my($class,@args) = @_;
    
    my $self = $class->SUPER::new(@args);

    my ($is_coding) =
	$self->_rearrange([qw(IS_CODING)],@args);
    # default is non-coding
    $self->is_coding(defined($is_coding) ? $is_coding : 0);

    return $self;
}



=head2 is_coding

 Title   : is_coding
 Usage   : if ($feature->is_coding()) {
                     #do something
            }
 Function: Whether or not the feature codes for amino acid.
 Returns : FALSE
 Args    : none

=cut

sub is_coding{
    my $self = shift;

    return $self->{'is_coding'} = shift if @_;
    return $self->{'is_coding'};
}

=head2 cds

 Title   : cds
 Usage   : $cds=$feature->cds();
 Function: get the coding sequence of this feature
 Returns : undef
 Args    : none


=cut

sub cds {
   my ($self,@args) = @_;
   return;

}


1;
