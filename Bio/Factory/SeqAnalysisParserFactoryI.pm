#
# BioPerl module for Bio::Factory::SeqAnalysisParserFactoryI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>,
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
L<Bio::SeqAnalysisParserI> for more documentation of this concept.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Hilmar Lapp, Jason Stajich

Email Hilmar Lapp E<lt>hlapp@gmx.netE<gt>, Jason Stajich E<lt>jason@bioperl.orgE<gt>

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Factory::SeqAnalysisParserFactoryI;
use strict;

use Carp;

use base qw(Bio::Root::RootI);

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
    $self->throw_not_implemented();
}

1;
