# $Id$
#
# BioPerl module for Bio::Factory::SeqAnalysisParserFactoryI
#
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>,
# and Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Jason Stajich, Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Factory::SeqAnalysisParserFactoryI - interface describing objects capable
              of creating SeqAnalysisParserI compliant parsers

=head1 SYNOPSIS

    # initialize an object implementing this interface, e.g.
    $factory = Bio::Factory::SeqAnalysisParserFactory->new();
    # obtain a parser object
    $parser = $factory->get_parser(-input=>$inputobj,
                                   -params=>[@params],
		                   -method => $method);
    # $parser is an object implementing Bio::SeqAnalysisParserI
    # annotate sequence with features produced by parser
    while(my $feat = $parser->next_feature()) {
	$seq->add_SeqFeature($feat);
    }

=head1 DESCRIPTION

This is an interface for factory classes capable of instantiating
SeqAnalysisParserI implementing parsers.

The concept behind the interface is a generic analysis result parsing
in high-throughput automated sequence annotation pipelines. See
Bio::SeqAnalysisParserI for more documentation of this concept.

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

Email Hilmar Lapp E<lt>hlapp@gmx.netE<gt>, Jason Stajich E<lt>jason@chg.mc.duke.eduE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Factory::SeqAnalysisParserFactoryI;
use strict;

#use Bio::Root::Root;
use Carp;

use vars qw(@ISA );
#@ISA = qw(Bio::Root::Root); 

=head2 get_parser

 Title   : get_parser
 Usage   : $factory->get_parser(-input=>$inputobj, 
                                [ -params=>[@params] ],
		                -method => $method)
 Function: Creates and returns a parser object for the given input and method.
           The type of input which is suitable depends on the implementation,
           but a good-style implementation should allow both file names and
           streams (filehandles).

           A particular implementation may not be able to create a parser for
           the requested method. In this case it shall return undef.

           Parameters (-params argument) are passed on to the parser object
           and therefore are specific to the parser to be created. An
           implementation of this interface should make this argument optional.
 Example :
 Returns : A Bio::SeqAnalysisParserI implementing object.
 Args    : B<input>  - object/file where analysis results are coming from,
	   B<params> - parameter to use when parsing/running analysis
	   B<method> - method of analysis

=cut

sub get_parser {
    my ( $self, @args) = @_;
    $self->_abstractDeath();
}

sub _abstractDeath {
  my $self = shift;
  my $package = ref $self;
  my $caller = (caller)[1];
  
  confess "Abstract method '$caller' defined in interface Bio::Factory::SeqAnalysisParserFactoryI not implemented by pacakge $package. Not your fault - author of $package should be blamed!";
}

1;
