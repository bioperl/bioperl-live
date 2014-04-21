#
# BioPerl module for Bio::SeqFeature::Gene::UTR
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

Bio::SeqFeature::Gene::UTR - A feature representing an untranslated region
          that is part of a transcriptional unit

=head1 SYNOPSIS

See documentation of methods

=head1 DESCRIPTION

A UTR is a Bio::SeqFeature::Gene::ExonI compliant object that is
non-coding, and can be either 5' or 3' in a transcript.

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

Email dblock@gene.pbi.nrc.ca

=head1 CONTRIBUTORS

This is based on the Gene Structure scaffolding erected by Hilmar Lapp
(hlapp@gmx.net).

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Gene::UTR;
use strict;

# Object preamble - inherits from Bio::Root::Root


use base qw(Bio::SeqFeature::Gene::Exon);

=head2 new

 Title   : new
 Usage   :
 Function: We override the constructor here to set is_coding to false
           unless explicitly overridden.

 Example :
 Returns : 
 Args    :


=cut

sub new{
    my ($caller, @args) = @_;

    if(! grep { lc($_) eq '-is_coding'; } @args) {
	push(@args, '-is_coding', 0);
    }
    my $self = $caller->SUPER::new(@args);

    my ($primary, $prim) = 
	$self->_rearrange([qw(PRIMARY PRIMARY_TAG)],@args);

    $self->primary_tag('utr') unless $primary || $prim;

    return $self;
}

=head2 primary_tag

 Title   : primary_tag
 Usage   : $tag = $feat->primary_tag()
 Function: Returns the primary tag for a feature,
           eg 'utr5prime'.  This method insures that 5prime/3prime information
           is uniformly stored
 Returns : a string 
 Args    : none

=cut

sub primary_tag{
    my $self = shift;
    if(@_ && defined($_[0])) {
	my $val = shift;
	if ($val =~ /(3|5)/ ) { 
	    $val = "utr$1prime";
	} else { 
	    $self->warn("Primary tag should indicate if this is 3 or 5'. ".
			"Preferred text is 'utr3prime' or 'utr5prime'.");
	}
	unshift(@_,$val);
    }
    return $self->SUPER::primary_tag(@_);
}

1;
