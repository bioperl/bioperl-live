# $Id$
#
# bioperl module for Bio::SeqFeature::Tools::TypeMapper
#
# Cared for by Chris Mungall <cjm@fruitfly.org>
#
# Copyright Chris Mungall
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Tools::TypeMapper - maps $seq_feature->primary_tag

=head1 SYNOPSIS

  use Bio::SeqIO;
  use Bio::SeqFeature::Tools::TypeMapper;

  # first fetch a genbank SeqI object
  $seqio =
    Bio::SeqIO->new(-file=>'AE003644.gbk',
                    -format=>'GenBank');
  $seq = $seqio->next_seq();
  
  $tm = Bio::SeqFeature::Tools::TypeMapper->new;

  # map all the types in the sequence
  $tm->map_types(-seq=>$seq,
		 {CDS=>'ORF',
		  variation=>sub {
		      my $f = shift;
		      $f->length > 1 ?
			'variation' : 'SNP'
		  },
		 });

   # alternatively, use the hardcoded SO mapping
   $tm->map_types_to_SO(-seq=>$seq);

=head1 DESCRIPTION

This class implements an object for mapping between types; for
example, the types in a genbank feature table, and the types specified
in the Sequence Ontology.

You can specify your own mapping, either as a simple hash index, or by
providing your own subroutines.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                         - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Chris Mungall

Email:  cjm@fruitfly.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::SeqFeature::Tools::TypeMapper;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root
use Bio::Root::Root;

@ISA = qw(Bio::Root::Root);

=head2 new

 Title   : new
 Usage   : $unflattener = Bio::SeqFeature::Tools::TypeMapper->new();
 Function: constructor
 Example : 
 Returns : a new Bio::SeqFeature::Tools::TypeMapper
 Args    : see below


=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my($typemap) =
	$self->_rearrange([qw(TYPEMAP
			     )],
                          @args);

    $typemap  && $self->typemap($typemap);
    return $self; # success - we hope!
}

=head2 typemap

 Title   : typemap
 Usage   : $obj->typemap($newval)
 Function: 
 Example : 
 Returns : value of typemap (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub typemap{
    my $self = shift;

    return $self->{'typemap'} = shift if @_;
    return $self->{'typemap'};
}

=head2 map_types

 Title   : map_types
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub map_types{
   my ($self,@args) = @_;

   my($sf, $seq, $type_map) =
     $self->_rearrange([qw(FEATURE
                           SEQ
			   TYPE_MAP
                          )],
                          @args);
   if (!$sf && !$seq) {
       $self->throw("you need to pass in either -feature or -seq");
   }

   my @sfs = ($sf);
   if ($seq) {
       $seq->isa("Bio::SeqI") || $self->throw("$seq NOT A SeqI");
       @sfs = $seq->get_all_SeqFeatures;
   }
   $type_map = $type_map || $self->type_map;
   foreach my $sf (@sfs) {

       $sf->isa("Bio::SeqFeatureI") || $self->throw("$sf NOT A SeqFeatureI");
       $sf->isa("Bio::FeatureHolderI") || $self->throw("$sf NOT A FeatureHolderI");

       my $type = $sf->primary_tag;
       my $mtype = $type_map->{$type};
       if ($mtype) {
	   if (ref($mtype)) {
	       if (ref($mtype) eq 'CODE') {
		   $mtype = $mtype->($sf);
	       }
	       else {
		   $self->throw('must be scalar or CODE ref');
	       }
	   }
	   $sf->primary_tag($mtype);
       }
   }
   return;
}

=head2 map_types_to_SO

 Title   : map_types_to_SO
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

hardcodes the genbank to SO mapping

Based on revision 1.22 of SO

Please see the actual code for the mappings

!!!NOT COMPLETE!!!

=cut

sub map_types_to_SO{
   my ($self,@args) = @_;

   push(@args,
	(-type_map=>{

		     # this is the most generic form for RNAs;
		     # we always represent the processed form of
		     # the transcript
		     misc_RNA=>'processed_transcript',

		     misc_feature=>'located_sequence_feature',

		     # not sure about this one...
		     source=>'database_entry',

		     LTR=>'LTR_retrotransposon',

		     rep_origin=>'origin_of_replication',
		     
		    }));
   return $self->map_types(@args);

}


1;
