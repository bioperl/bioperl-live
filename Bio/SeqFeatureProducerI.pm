#
# BioPerl module for Bio::SeqFeatureProducerI
#
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>,
# and Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Jason Stajich, Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeatureProducerI - interface describing objects capable of annotating
              sequences with features

=head1 SYNOPSIS

    # get a SeqFeatureProducerI somehow
    $seqfeatureproducer->add_features($seq,$seqanalparser);

=head1 DESCRIPTION

SeqFeatureProducerI is an interface for describing objects that will
add Sequence Features to sequences.  

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                - General discussion
  http://bio.perl.org/MailList.html    - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Hilmar Lapp, Jason Stajich

Email Hilmar Lapp <hlapp@gmx.net>, Jason Stajich <jason@chg.mc.duke.edu>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::SeqFeatureProducerI;
use strict;

use Bio::Root::RootI;
use Bio::SeqFeatureI;
use Carp;

use vars qw(@ISA );
@ISA = qw(Bio::Root::RootI); 

=head2 add_features

 Title   : add_features
 Usage   : $obj->add_features($seq,$analfeatParser)
 Function: adds features to a sequence based on the analysis parser
 Example :
 Returns : void
 Args    : B<seq>    - a sequence object to receive new features
	   B<analFeatureParser> - a SeqAnalysisParserI compliant object    

=cut

sub add_features {
    my ($self, $seq, $analfeatparser) = @_;
    $self->_abstractDeath();
}

sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  
  confess "Abstract method '$caller' defined in interface Bio::SeqFeatureProducerI not implemented by pacakge $package. Not your fault - author of $package should be blamed!";
}

1;
