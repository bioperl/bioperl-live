#
# BioPerl module for Bio::HandlerI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chris Fields
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::HandlerBaseI - Interface class for handler methods which interact with any
event-driven parsers (drivers).

=head1 SYNOPSIS

  # MyHandler is a Bio::HandlerBaseI-derived class for dealing with GenBank
  # sequence data, derived from a GenBank event-driven parser

  # inside a parser (driver) constructor

  $self->seqhandler($handler || MyHandler->new(-format => 'genbank'));

  # in the driver parsing method ( such as next_seq() ) ...

  $handler = $self->seqhandler();

  # roll data up into hashref chunks, pass off into Handler for processing...

  $hobj->data_handler($data);

  # or retrieve Handler methods and pass data directly to Handler methods

  my $hmeth = $hobj->handler_methods;

  if ($hmeth->{ $data->{NAME} }) {
      my $mth = $hmeth->{ $data->{NAME} }; # code ref
      $hobj->$mth($data);
  }

=head1 DESCRIPTION

This interface describes simple class methods used for processing data from an
event-based parser (a driver). This is similar in theme to an XML SAX-based
driver but differs in that one can optionally pass related data
semi-intelligently as chunks (defined in a hash reference) vs. passing as single
data elements in a stream. For instance, any reference-related and
species-related data as well as individual sequence features could be passed as
chunks of data to be processed in part or as a whole (from Data::Dumper output):

Annotation Data (References):

  $VAR1 = {
          'NAME' => 'REFERENCE',
          'DATA' => '1  (bases 1 to 10001)'
          'AUTHORS' => 'International Human Genome Sequencing Consortium.'
          'TITLE' => 'The DNA sequence of Homo sapiens'
          'JOURNAL' => 'Unpublished (2003)'
          };

Sequence features (source seqfeature):

  $VAR1 = {
          'mol_type' => 'genomic DNA',
          'LOCATION' => '<1..>10001',
          'NAME' => 'FEATURES',
          'FEATURE_KEY' => 'source',
          'note' => 'Accession AL451081 sequenced by The Sanger Centre',
          'db_xref' => 'taxon:9606',
          'clone' => 'RP11-302I18',
          'organism' => 'Homo sapiens'
          };

These would be 'handled' accordingly by methods specified in a
HandlerI-based class. The data in a chunk is intentionally left vague
here since this may vary from implementation to implementation and can
be somewhat open to interpretation. A data chunk in a sequence record,
for instance, will be different than a data chunk in a BLAST
report. This also allows one the flexibility to pass data as more
XML-like small bits, as huge chunks, or even as indexed locations in a
file (such as when using a "pull" parser, like a Bio::PullParserI).

For an sequence-based implementation see
Bio::SeqIO::RichSeq::GenericRichSeqHandler, which handles any GenBank,
UniProt, and EMBL data from their respective driver modules
(Bio::SeqIO::gbdriver, Bio::SeqIO::swissdriver, and
Bio::SeqIO::embldriver).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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

=head1 AUTHOR - Chris Fields

Email cjfields at bioperl dot org

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::HandlerBaseI;
use strict;
use warnings;

use base qw(Bio::Root::RootI);

my %HANDLERS = ('foo' => \&noop);

=head2 data_handler

 Title   :  data_handler
 Usage   :  $handler->data_handler($data)
 Function:  Centralized method which accepts all data chunks, then distributes
            to the appropriate methods for processing based on the chunk name
            from within the HandlerBaseI object.

            One can also use
 Returns :  None
 Args    :  an hash ref containing a data chunk.

=cut

sub data_handler {
    shift->throw_not_implemented
}

=head2 handler_methods

 Title   :  handler_methods
 Usage   :  $handler->handler_methods('GenBank')
            %handlers = $handler->handler_methods();
 Function:  Retrieve the handler methods used for the current format() in
            the handler.  This assumes the handler methods are already
            described in the HandlerI-implementing class.
 Returns :  a hash reference with the data type handled and the code ref
            associated with it.
 Args    :  [optional] String representing the sequence format.  If set here
            this will also set sequence_format()
 Throws  :  On unimplemented sequence format in %HANDLERS

=cut

sub handler_methods {
    shift->throw_not_implemented
}

=head2 format

 Title   :  format
 Usage   :  $handler->format('GenBank')
            $handler->format('BLAST')
 Function:  Get/Set the format for the report/record being parsed. This can be
            used to set handlers in classes which are capable of processing
            similar data chunks from multiple driver modules.
 Returns :  String with the sequence format
 Args    :  [optional] String with the sequence format
 Note    :  The format may be used to set the handlers (as in the
            current GenericRichSeqHandler implementation)

=cut

sub format {
    shift->throw_not_implemented
}

=head2 get_params

 Title   :  get_params
 Usage   :  $handler->get_params('-species')
 Function:  Convenience method used to retrieve the specified
            parameters from the internal parameter cache
 Returns :  Hash ref containing parameters requested and data as
            key-value pairs.  Note that some parameter values may be
            objects, arrays, etc.
 Args    :  List (array) representing the parameters requested

=cut

sub get_params {
    shift->throw_not_implemented
}

=head2 set_params

 Title   :  set_params
 Usage   :  $handler->set_params({
                                '-species' => $species,
                                '-accession_number' => $acc
                                });
 Function:  Convenience method used to set specific parameters
 Returns :  None
 Args    :  Hash ref containing the data to be passed as key-value pairs

=cut

sub set_params {
    shift->throw_not_implemented
}

=head2 reset_parameters

 Title   :  reset_parameters
 Usage   :  $handler->reset_parameters()
 Function:  Resets the internal cache of data (normally object parameters for
            a builder or factory)
 Returns :  None
 Args    :  None

=cut

sub reset_parameters {
    shift->throw_not_implemented
}

1;

