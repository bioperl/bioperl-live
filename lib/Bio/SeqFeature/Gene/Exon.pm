#
# BioPerl module for Bio::SeqFeature::Gene::Exon
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Hilmar Lapp <hlapp@gmx.net>
#
# Copyright Hilmar Lapp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Gene::Exon - a feature representing an exon

=head1 SYNOPSIS

    # obtain an exon instance $exon somehow
    print "exon from ", $exon->start(), " to ", $exon->end(),
          " on seq ", $exon->seq_id(), ", strand ", $exon->strand(),
          ", encodes the peptide sequence ", 
          $exon->cds()->translate()->seq(), "\n";

=head1 DESCRIPTION

This module implements a feature representing an exon by implementing
the Bio::SeqFeature::Gene::ExonI interface. By default an Exon is
coding. Supply -is_coding =E<gt> 0 to the constructor or call
$exon-E<gt>is_coding(0) otherwise.

Apart from that, this class also implements Bio::SeqFeatureI by
inheriting off Bio::SeqFeature::Generic.

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

=head1 AUTHOR - Hilmar Lapp

Email hlapp@gmx.net

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Gene::Exon;
use strict;


use base qw(Bio::SeqFeature::Generic Bio::SeqFeature::Gene::ExonI);

#
# A list of allowed exon types. See primary_tag().
#
my @valid_exon_types = ('initial', 'internal', 'terminal');

sub new {
    my ($caller, @args) = @_;
    my $self = $caller->SUPER::new(@args);

    my ($is_coding) =
	$self->_rearrange([qw(IS_CODING)],@args);
    $self->primary_tag('exon') unless $self->primary_tag();
    $self->is_coding(defined($is_coding) ? $is_coding : 1);
    $self->strand(0) if(! defined($self->strand()));
    return $self;
}


=head2 is_coding

 Title   : is_coding
 Usage   : if($exon->is_coding()) {
                   # do something
           }
           if($is_utr) {
               $exon->is_coding(0);
           }
 Function: Get/set whether or not the exon codes for amino acid.
 Returns : TRUE if the object represents a feature translated into protein,
           and FALSE otherwise.
 Args    : A boolean value on set.


=cut

sub is_coding {
    my ($self,$val) = @_;

    if(defined($val)) {
	$self->{'_iscoding'} = $val;
    }
    return $self->{'_iscoding'};
}

=head2 primary_tag

 Title   : primary_tag
 Usage   : $tag = $feat->primary_tag()
           $feat->primary_tag('exon')
 Function: Get/set the primary tag for the exon feature.

           This method is overridden here in order to allow only for
           tag values following a certain convention. For consistency reasons,
           the tag value must either contain the string 'exon' or the string
           'utr' (both case-insensitive). In the case of 'exon', a string
           describing the type of exon may be appended or prefixed. Presently,
           the following types are allowed: initial, internal, and terminal
           (all case-insensitive). 

           If the supplied tag value matches 'utr' (case-insensitive),
           is_coding() will automatically be set to FALSE, and to TRUE
           otherwise.

 Returns : A string.
 Args    : A string on set.


=cut

# sub primary_tag {
#    my ($self,$value) = @_;

#    if(defined($value)) {
#        if((lc($value) =~ /utr/i) || (lc($value) eq "exon") ||
# 	  ((lc($value) =~ /exon/i) &&
# 	   (grep { $value =~ /$_/i; } @valid_exon_types))) {
# 	   $self->is_coding($value =~ /utr/i ? 0 : 1);
#        } else {
# 	   $self->throw("primary tag $value is invalid for object of class ".
# 			ref($self));
#        }
#    }
#    return $self->SUPER::primary_tag($value);
# }

=head2 location

 Title   : location
 Usage   : my $location = $exon->location()
 Function: Returns a location object suitable for identifying the location 
	   of the exon on the sequence or parent feature.

           This method is overridden here to restrict allowed location types
           to non-compound locations.

 Returns : Bio::LocationI object
 Args    : none


=cut

sub location {
   my ($self,$value) = @_;  

   if(defined($value) && $value->isa('Bio::Location::SplitLocationI')) {
       $self->throw("split or compound location is not allowed ".
		    "for an object of type " . ref($self));
   }
   return $self->SUPER::location($value);
}

=head2 cds

 Title   : cds()
 Usage   : $cds = $exon->cds();
 Function: Get the coding sequence of the exon as a sequence object.

           The sequence of the returned object is prefixed by Ns (lower case)
           if the frame of the exon is defined and different from zero. The
           result is that the first base starts a codon (frame 0).

           This implementation returns undef if the particular exon is
           not translated to protein, i.e., is_coding() returns FALSE. Undef
           will also be returned if no sequence is attached to this exon
           feature.

 Returns : A Bio::PrimarySeqI implementing object.
 Args    : 


=cut

sub cds {
    my ($self) = @_;

    # UTR is not translated
    return if(! $self->is_coding());

    my $seq = $self->seq();
    if(defined($seq) && defined($self->frame()) && ($self->frame() != 0)) {
	my $prefix = "n" x $self->frame();
	$seq->seq($prefix . $seq->seq());
    }
    return $seq;
}

1;
