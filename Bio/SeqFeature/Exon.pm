#
# $Id$
#

=head1 NAME

Bio::SeqFeature::Exon

=head1 SYNOPSIS

  use Bio::SeqFeature::Exon;
  my $Exon = Bio::SeqFeature::Exon->new( # takes same args as SeqFeature::Generic
      -start => 10,
      -end => 100,
      -frame => '.',
      -primary => 'high_GC_exon',
      -strand => '1',
      -source => 'cDNA_alignment',
      -score => '100',
      type => 'internal');   # Descr. of the Exon type; *not* Cont.Vocab.


=head1 DESCRIPTION

Creates Intron type sequence features.  These are essentially SeaFeature::Generic
features, but report themselves as "Bio::SeqFeature::Intron" when you query them
with a $Feature-E<gt>isa.

=head1 AUTHORS - Hilmar Lapp, Mark Wilkinson

Based on original code and concept from Hilmar Lapp (hlapp@gmx.net)

Mark Wilkinson (mwilkinson@gene.pbi.nrc.ca)
Plant Biotechnology Institute, National Research Council of Canada.
Copyright (c) National Research Council of Canada, April, 2001.

=head1 DISCLAIMER

Anyone who intends to use and uses this software and code acknowledges
and agrees to the following: The National Research Council of Canada
(herein "NRC") disclaims any warranties, expressed, implied, or
statutory, of any kind or nature with respect to the software,
including without limitation any warranty or merchantability or
fitness for a particular purpose.  NRC shall not be liable in any
event for any damages, whether direct or indirect, consequential or
incidental, arising from the use of the software.

=head1 CONTACT

Mark Wilkinson (mwilkinson@gene.pbi.nrc.ca)

=head1 METHODS

identical to SeqFeature::Generic except for:

=head2 type

  Usage:  $Exon->type($type);
  Args:   optional string indicating new type
  Returns: current or newly set type


=cut


package Bio::SeqFeature::Exon;

use strict;
use Carp;
use vars qw(@ISA $AUTOLOAD);
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Gene::ExonI;

@ISA = qw(Bio::SeqFeature::Generic Bio::SeqFeature::Gene::ExonI);

{
	#class data	
	#___________________________________________________________
	#ATTRIBUTES (default)
    my %_attr_data = #     				DEFAULT    	ACCESSIBILITY
                  (
					valid_types  		=> 	[['initial', 'internal', 'terminal', 'single'],    'read/write'],  # should be something like 'init', 'intr', 'term', 'sngl' (GenScan exon types)					                  	
                    is_coding			=>  [1, 		'read/write'],   # set to "true" as a default
                    type				=>  ['internal', 		'read/write'],  # set default type as internal (should think about this.. perhaps undef is better?)

                    );
   #_____________________________________________________________
   #Class attribs and methods

    # Is a specified object attribute accessible in a given mode
    sub _accessible  {
    my ($self, $attr, $mode) = @_;
    $_attr_data{$attr}[1] =~ /$mode/
    }

    # Classwide default value for a specified object attribute
    sub _default_for {
    my ($self, $attr) = @_;
    $_attr_data{$attr}[0];
    }

    # List of names of all specified object attributes
    sub _standard_keys {
    keys %_attr_data;
	}
	
	sub _valid_types {
		@{$_attr_data{valid_types}[0]};  # return as a list
	}
			

}

sub new {
	my ($caller, %args) = @_;
	
	#create a generic feature based on %args
	# is this a call to duplicate an object?
	my $caller_is_obj = ref($caller);
    my $class = $caller_is_obj || $caller;

    my $self = $caller->SUPER::new(%args);
    my $primary = $args{primary}?$args{primary}:'exon'; # if primary tag is not set, then set to exon
    $self->primary_tag($primary);
    $self->strand(0) unless (defined ($self->strand));

    foreach my $attrname ( $self->_standard_keys ) {   # SET OBJECT PROPERTIES
    	if (exists $args{$attrname}) {                 # based on initiating args
			$self->{$attrname} = $args{$attrname} }
		elsif ($caller_is_obj) {                       # OR based on existing object
			$self->{$attrname} = $caller->{$attrname} }
    	else {
			$self->{$attrname} = $self->_default_for($attrname) } # OR based on the encapsulated default values
    }   	


    # we need to ensure that 'type' set above is valid based on the current valid_types list
    if ($self->type($self->type)){  # simply re-send it it's current value and have a look at what is returned
    	return $self;
    } else {warn "type ".$args{type}." not allowed.  Valid exon type is required - Exon not created";
    	return 0;
    }

}
=head2 valid_types

 Title   : valid_types
 Usage   : $type_list_ref = $exon->valid_types()
           $exon->valid_types(\@types)
 Function: Get/set the valid types for the exon feature.

           Creates a controlled vocab. for the type. Defaults to
           ["initial", "terminal" "internal", "single"]

 Returns : A list ref.
 Args    : A list ref on set.

=cut

=head2 type

 Title   : type
 Usage   : $tag = $exon->type()
           $exon->type('internal')
 Function: Get/set the type for the exon feature.

           Conrolled Vocab. based on $exon->valid_types

 Returns : A string.
 Args    : A string on set.

=cut

sub type {
	my ($self, $type) = @_;
	unless ($type){return $self->{type}}  # if this is just a call to get the value, then return the value
	
	if (grep {$type =~ /$_/i} $self->_valid_types){   # grep it against the valid types.  If it matches
		$self->{type} = $type;                        # then set
		return 1;  # valid type                       # and return true
	} else {return 0} # invalid type                  # otherwise return false
}

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

sub location {  # copied directly from Hilmars Exon object
   my ($self,$value) = @_;

   if(defined($value) && $value->isa('Bio::Location::SplitLocationI')) {
       $self->throw("split or compound location is not allowed ".
		    "for an object of type " . ref($self));
   }
   return $self->SUPER::location($value);
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

sub is_coding { # this can be set independantly of exon type, since not all exons are coding
	my ($self, $coding) = @_;
	unless (defined $coding) {return $self->{is_coding}}
	if ($coding == 1){$self->{is_coding} = 1; return 1}
	elsif ($coding == 0){$self->{is_coding} = 0; return 1}
	else {warn "$coding is an invalid value for is_coding.  Must be '1' or '0'"; return 0}
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
    return undef unless ($self->is_coding());

    my $seq = $self->seq();
    if(defined($seq) && defined($self->frame()) && ($self->frame() != 0)) {
		my $prefix = "n" x $self->frame();
		$seq->seq($prefix . $seq->seq());
    }
    return $seq;
}


sub AUTOLOAD {
    no strict "refs";
    my ($self, $newval) = @_;

    $AUTOLOAD =~ /.*::(\w+)/;

    my $attr=$1;
    if ($self->_accessible($attr,'write')) {

	*{$AUTOLOAD} = sub {
	    if (defined $_[1]) { $_[0]->{$attr} = $_[1] }
	    return $_[0]->{$attr};
	};    ### end of created subroutine

###  this is called first time only
	if (defined $newval) {
	    $self->{$attr} = $newval
	}
	return $self->{$attr};

    } elsif ($self->_accessible($attr,'read')) {

	*{$AUTOLOAD} = sub {
	    return $_[0]->{$attr} }; ### end of created subroutine
	return $self->{$attr}  }


    # Must have been a mistake then...
    croak "No such method: $AUTOLOAD";
}

1;


