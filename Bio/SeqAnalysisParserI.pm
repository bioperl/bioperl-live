#
# BioPerl module for Bio::SeqAnalysisParserI
#
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>,
# and Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Jason Stajich, Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqAnalysisParserI - Sequence analysis output parser interface  

=head1 SYNOPSIS

    # get a SeqAnalysisParserI somehow
    $seqanalparser->parse('inputfile');
    while( my $feature = $seqanalparser->next_feature() ) {
	print "Feature from ", $feature->start, " to ", $feature->end, "\n";
    }

=head1 DESCRIPTION

SeqAnalysisParserI is an interface for describing AnalysisResult parsers.  These parsers will produce features that can be used to annotate a sequence object.

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

package Bio::SeqAnalysisParserI;
use strict;

use Bio::Root::RootI;
use Bio::SeqFeatureI;
use Carp;

use vars qw(@ISA );
@ISA = qw(Bio::Root::RootI); 

=head2 parse

 Title   : parse
 Usage   : $obj->parse(-input=>$inputobj, 
           $obj->parse(-input=>$inputobj, -params=>[@params],
		       [ -method => $method ] )
 Function: sets up parsing for feature retrieval from an analysis file, 
           or object
 Example :
 Returns : void
 Args    : B<input>  - object/file where analysis are coming from,
	   B<params> - parameter to use when parsing/running analysis
	   B<method> - method of analysis (optional)
    
=cut

sub parse {
    my ( $self, @args) = @_;

    $self->_abstractDeath();
}

=head2 next_feature

 Title   : next_feature
 Usage   : $seqfeature = $obj->next_feature();
 Function: returns the next SeqFeatureI object available
 Example :
 Returns : SeqFeatureI
 Args    : none    
=cut

sub next_feature {
    my ($self);
    $self->_abstractDeath();
}

sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  
  confess "Abstract method '$caller' defined in interface Bio::SeqAnalsysiParserI not implemented by pacakge $package. Not your fault - author of $package should be blamed!";
}

1;
