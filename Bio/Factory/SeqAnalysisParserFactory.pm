# $Id$
#
# BioPerl module for Bio::Factory::SeqAnalysisParserFactory
#
# Cared for by Jason Stajich <jason@chg.mc.duke.edu>,
# and Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Jason Stajich, Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Factory::SeqAnalysisParserFactoryI - class capable of
              creating SeqAnalysisParserI compliant parsers

=head1 SYNOPSIS

    # initialize an object implementing this interface, e.g.
    $factory = Bio::Factory::SeqAnalysisParserFactory->new();
    # find out the methods it knows about
    print "registered methods: ",
          join(', ', keys($factory->driver_table())), "\n";
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

This is a factory class capable of instantiating SeqAnalysisParserI 
implementing parsers.

The concept behind this class and the interface it implements
(Bio::Factory::SeqAnalysisParserFactoryI) is a generic analysis result parsing
in high-throughput automated sequence annotation pipelines. See
Bio::SeqAnalysisParserI for more documentation of this concept.

You can always find out the methods an instance of this class knows
about by the way given in the SYNOPSIS section. By default, and
assuming that the documentation is up-to-date, this will comprise of
genscan, mzef, estscan, blast, hmmer, gff, and sim4 (all case-insensitive).

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

package Bio::Factory::SeqAnalysisParserFactory;
use strict;

use Bio::Factory::SeqAnalysisParserFactoryI;
use Bio::Factory::DriverFactory;

use vars qw(@ISA);
@ISA = qw(Bio::Factory::DriverFactory Bio::Factory::SeqAnalysisParserFactoryI);

BEGIN {
    Bio::Factory::DriverFactory->register_driver(
			 "genscan" => "Bio::Tools::Genscan",
			 "mzef" => "Bio::Tools::MZEF",
			 "estscan" => "Bio::Tools::ESTScan",
			 "bplite" => "Bio::Tools::BPlite",
			 "blast" => "Bio::Tools::BPlite",
                         "hmmer" => "Bio::Tools::HMMER::Result",
			 "gff" => "Bio::Tools::GFF",
			 "sim4" => "Bio::Tools::Sim4::Results",
						 );
}

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    # no per-object initialization right now - registration of default drivers
    # is only done once when the module is loaded
    return $self;
}

=head2 get_parser

 Title   : get_parser
 Usage   : $factory->get_parser(-input=>$inputobj, 
                                [ -params=>[@params] ],
		                -method => $method)
 Function: Creates and returns a parser object for the given input and method.
           Both file names and streams (filehandles) are allowed.

           Parameters (-params argument) are passed on to the parser object
           and therefore are specific to the parser to be created.
 Example :
 Returns : A Bio::SeqAnalysisParserI implementing object. Exception if
           creation of the parser object fails.
 Args    : B<input>  - object/file where analysis results are coming from,
	   B<params> - parameter to use when parsing/running analysis
	   B<method> - method of analysis
    
=cut

sub get_parser {
    my ($self, @args) = @_;
    my $parser;
    my $module;
    
    my ($input, $params, $method) =
	$self->_rearrange([qw(INPUT PARAMS METHOD)], @args);
    
    # retrieve module name for requested method
    $method = lc $method; # method is case-insensitive
    $module = $self->get_driver($method);
    if(! defined($module)) {
	$self->throw("Analysis parser driver for method $method not registered.");
    }
    # load module
    $self->_load_module($module); # throws an exception on failure to load
    # make sure parameters is not undef
    $params = [] if( !defined $params );
    # figure out input method (file or stream)
    my $inputmethod = '-file';
    if( ref($input) =~ /GLOB/i ) {
	$inputmethod = '-fh';
    }
    # instantiate parser and return the result
    $parser = $module->new($inputmethod => $input, @$params);
    if(! $parser->isa('Bio::SeqAnalysisParserI')) {
	$self->throw("Driver $module registered for method $method does not implement Bio::SeqAnalyisParserI. How come?");
    }
    return $parser;
}

1;
