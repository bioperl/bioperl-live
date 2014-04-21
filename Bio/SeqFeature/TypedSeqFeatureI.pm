#
# BioPerl module for Bio::SeqFeature::TypedSeqFeatureI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::TypedSeqFeatureI - a strongly typed SeqFeature

=head1 SYNOPSIS


   # get Sequence Features in some manner, eg
   # from a Sequence object

    foreach $sf ( $seq->get_SeqFeatures() ) {
        # all sequence features must have primary_tag() return a string
        $type_as_string = $sf->primary_tag();
	if( $sf->isa("Bio::SeqFeature::TypedSeqFeatureI") ) {
            $ot = $sf->ontology_term();
            print "Ontology identifier:",$ot->identifier(),
                  " name:",$ot->name(),
                  " Description:",$ot->description(),"\n";

        } else {
            print "Sequence Feature does not have an ontology type\n";
	}

    }

=head1 DESCRIPTION

This interface describes the extension of SeqFeatureI 
to being a strongly typed SeqFeature.

Bio::SeqFeature::TypedSeqFeatureI extends the Bio::SeqFeatureI
interface (ie, a TypedSeqFeatureI feature must also implement
all the Bio::SeqFeatureI interface as well). 

It is suggested that the primary_tag() method of SeqFeatureI
return the same as the ontology_term()-E<gt>name() of the OntologyTypedI
(ie, the "string" name of the ontology type is used as the primary
tag), but this should not be assummed by client code as they
are scenarios where one would like to maintain the difference.


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

=head1 AUTHOR - Ewan Birney

Email - please email the BioPerl mailing list above.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::TypedSeqFeatureI;

use strict;
use Bio::Root::RootI;

use base qw(Bio::SeqFeatureI);


=head2 ontology_term

  Title   : ontology_term
  Usage   : my $ot = $seqfeature->ontology_term()
  Returns : a Bio::Ontology::TermI compliant object
  Args    : none
  Status  : public

This method returns the ontology term for a 
strongly typed sequence feature. 

=cut

sub ontology_term {
    shift->throw_not_implemented();
}

1;
