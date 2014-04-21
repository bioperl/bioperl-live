#
# BioPerl module for Bio::Factory::SeqAnalysisParserFactory
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

Bio::Factory::SeqAnalysisParserFactory - class capable of creating
SeqAnalysisParserI compliant parsers

=head1 SYNOPSIS

    # initialize an object implementing this interface, e.g.
    $factory = Bio::Factory::SeqAnalysisParserFactory->new();
    # find out the methods it knows about
    print "registered methods: ",
          join(', ', keys %{$factory->driver_table}), "\n";
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
(Bio::Factory::SeqAnalysisParserFactoryI) is a generic analysis result
parsing in high-throughput automated sequence annotation
pipelines. See Bio::SeqAnalysisParserI for more documentation of this
concept.

You can always find out the methods an instance of this class knows
about by the way given in the SYNOPSIS section. By default, and
assuming that the documentation is up-to-date, this will comprise of
genscan, mzef, estscan, blast, hmmer, gff, and sim4 (all
case-insensitive).

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

package Bio::Factory::SeqAnalysisParserFactory;
use strict;


use base qw(Bio::Factory::DriverFactory Bio::Factory::SeqAnalysisParserFactoryI);

BEGIN {
    Bio::Factory::DriverFactory->register_driver
	(
	 "genscan"   => "Bio::Tools::Genscan",
	 "mzef"      => "Bio::Tools::MZEF",
	 "estscan"   => "Bio::Tools::ESTScan",
	 "hmmer"     => "Bio::Tools::HMMER::Result",
	 "gff"       => "Bio::Tools::GFF",
	 "sim4"      => "Bio::Tools::Sim4::Results",
	 "epcr"      => "Bio::Tools::EPCR", 
	 "exonerate" => "Bio::Tools::Exonerate",
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
	$self->throw("Driver $module registered for method $method does not ".
                     "implement Bio::SeqAnalyisParserI. How come?");
    }
    return $parser;
}


=head2 register_driver

 Title   : register_driver
 Usage   : $factory->register_driver("genscan", "Bio::Tools::Genscan");
 Function: Registers a driver a factory class should be able to instantiate.

           This method can be called both as an instance and as a
           class method.

 Returns : 
 Args    : Key of the driver (string) and the module implementing the driver
           (string).

=cut

=head2 driver_table

 Title   : driver_table
 Usage   : $table = $factory->driver_table();
 Function: Returns a reference to the hash table storing associations of
           methods with drivers.

           You use this table to look up registered methods (keys) and
           drivers (values).

           In this implementation the table is class-specific and
           therefore shared by all instances. You can override this in
           a derived class, but note that this method can be called
           both as an instance and a class method.

           This will be the table used by the object internally. You
           should definitely know what you're doing if you modify the
           table's contents.  Modifications are shared by _all_
           instances, those present and those yet to be created.

 Returns : A reference to a hash table.
 Args    : 


=cut


1;
