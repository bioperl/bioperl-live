#
# BioPerl module for Bio::Seq::RichSeq
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

Bio::Seq::RichSeq - Module implementing a sequence created from a rich
sequence database entry

=head1 SYNOPSIS

See L<Bio::Seq::RichSeqI> and documentation of methods.

=head1 DESCRIPTION

This module implements Bio::Seq::RichSeqI, an interface for sequences
created from or created for entries from/of rich sequence databanks,
like EMBL, GenBank, and SwissProt. Methods added to the Bio::SeqI
interface therefore focus on databank-specific information. Note that
not every rich databank format may use all of the properties provided.

For more information, please see the relevant 

=head1 Implemented Interfaces

This class implementes the following interfaces.

=over 4

=item L<Bio::Seq::RichSeqI>

Note that this includes implementing L<Bio::PrimarySeqI> and L<Bio::SeqI>,
specifically via L<Bio::Seq> and L<Bio::PrimarySeq>. Please review the
documentation for those modules on implementation details relevant to those
interfaces, as well as the ones below.

=item L<Bio::IdentifiableI>

=item L<Bio::DescribableI>

=item L<Bio::AnnotatableI>

=back

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

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::RichSeq;
use vars qw($AUTOLOAD);
use strict;



use base qw(Bio::Seq Bio::Seq::RichSeqI);


=head2 new

 Title   : new
 Usage   : $seq    = Bio::Seq::RichSeq->new( -seq => 'ATGGGGGTGGTGGTACCCT',
                                             -id  => 'human_id',
				             -accession_number => 'AL000012',
				            );

 Function: Returns a new seq object from
           basic constructors, being a string for the sequence
           and strings for id and accession_number
 Returns : a new Bio::Seq::RichSeq object

=cut

sub new {
    # standard new call..
    my($caller,@args) = @_;
    my $self = $caller->SUPER::new(@args);
    
    $self->{'_dates'} = [];
    $self->{'_secondary_accession'} = [];

    my ($dates, $xtra, $sv,
	$keywords, $pid, $mol, 
	$division ) = $self->_rearrange([qw(DATES 
					    SECONDARY_ACCESSIONS
					    SEQ_VERSION 
					    KEYWORDS
					    PID
					    MOLECULE
					    DIVISION
					    )],
					@args);
    defined $division && $self->division($division);
    defined $mol && $self->molecule($mol);
    if(defined($keywords)) {
	if(ref($keywords) && (ref($keywords) eq "ARRAY")) {
	    $self->add_keyword(@$keywords);
	} else {
	    # got a string - use the old API
	    $self->keywords($keywords);
	}
    }
    defined $sv && $self->seq_version($sv);
    defined $pid && $self->pid($pid);

    if( defined $dates ) {
	if( ref($dates) eq "ARRAY" ) {
	    foreach ( @$dates) {
		$self->add_date($_);
	    } 
	} else { 
	    $self->add_date($dates);
	}
    }

    if( defined $xtra ) {
	if( ref($xtra) eq "ARRAY" ) {
	    foreach ( @$xtra) {
		$self->add_secondary_accession($_);
	    } 
	} else { 
	    $self->add_secondary_accession($xtra);
	}
    }
    
    return $self;
}


=head2 division

 Title   : division
 Usage   : $obj->division($newval)
 Function: 
 Returns : value of division
 Args    : newvalue (optional)


=cut

sub division {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_division'} = $value;
    }
    return $obj->{'_division'};

}

=head2 molecule

 Title   : molecule
 Usage   : $obj->molecule($newval)
 Function: 
 Returns : type of molecule (DNA, mRNA)
 Args    : newvalue (optional)


=cut

sub molecule {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_molecule'} = $value;
    }
    return $obj->{'_molecule'};

}

=head2 add_date

 Title   : add_date
 Usage   : $self->add_date($datestr)
 Function: adds one or more dates

           This implementation stores dates as keyed annotation, the
           key being 'date_changed'. You can take advantage of this
           fact when accessing the annotation collection directly.

 Example :
 Returns : 
 Args    : a date string or an array of such strings


=cut

sub add_date {
    return shift->_add_annotation_value('date_changed',@_);
}

=head2 get_dates

 Title   : get_dates
 Usage   : my @dates = $seq->get_dates;
 Function: Get the dates of the sequence (usually, when it was created and
           changed.
 Returns : an array of date strings
 Args    :


=cut

sub get_dates{
    return shift->_get_annotation_values('date_changed');
}


=head2 pid

 Title   : pid
 Usage   : my $pid = $seq->pid();
 Function: Get (and set, depending on the implementation) the PID property
           for the sequence.
 Returns : a string
 Args    :


=cut

sub pid{
    my $self = shift;

    return $self->{'_pid'} = shift if @_;
    return $self->{'_pid'};
}


=head2 accession

 Title   : accession
 Usage   : $obj->accession($newval)
 Function: Whilst the underlying sequence object does not 
           have an accession, so we need one here.

           In this implementation this is merely a synonym for
           accession_number().
 Example : 
 Returns : value of accession
 Args    : newvalue (optional)


=cut

sub accession {
   my ($obj,@args) = @_;
   return $obj->accession_number(@args);
}

=head2 add_secondary_accession

 Title   : add_secondary_accession
 Usage   : $self->add_domment($ref)
 Function: adds a secondary_accession

           This implementation stores secondary accession numbers as
           keyed annotation, the key being 'secondary_accession'. You
           can take advantage of this fact when accessing the
           annotation collection directly.

 Example :
 Returns : 
 Args    : a string or an array of strings


=cut

sub add_secondary_accession {
    return shift->_add_annotation_value('secondary_accession',@_);
}

=head2 get_secondary_accessions

 Title   : get_secondary_accessions
 Usage   : my @acc = $seq->get_secondary_accessions();
 Function: Get the secondary accession numbers as strings.
 Returns : An array of strings
 Args    : none


=cut

sub get_secondary_accessions{
    return shift->_get_annotation_values('secondary_accession');
}

=head2 seq_version

 Title   : seq_version
 Usage   : $obj->seq_version($newval)
 Function: Get/set the sequence version
 Returns : value of seq_version (a scalar)
 Args    : on set, new value (a scalar or undef, optional)
 Note    : this differs from Bio::PrimarySeq version() in that this explicitly
           refers to the sequence record version one would find in a typical
           sequence file.  

=cut

sub seq_version{
    my $self = shift;

    return $self->{'_seq_version'} = shift if @_;
    return $self->{'_seq_version'};
}


=head2 add_keyword

 Title   : add_keyword
 Usage   : $obj->add_keyword($newval)
 Function: Add a new keyword to the annotation of the sequence.

           This implementation stores keywords as keyed annotation,
           the key being 'keyword'. You can take advantage of this
           fact when accessing the annotation collection directly.

 Returns : 
 Args    : value to be added (optional) (a string)


=cut

sub add_keyword {
    return shift->_add_annotation_value('keyword',@_);
}

=head2 get_keywords

 Title   : get_keywords
 Usage   : $obj->get_keywords($newval)
 Function: Get the keywords for this sequence as an array of strings.
 Returns : an array of strings
 Args    : 


=cut

sub get_keywords {
    return shift->_get_annotation_values('keyword');
}

=head1 Private methods and synonyms for backward compatibility

=cut

=head2 _add_annotation_value

 Title   : _add_annotation_value
 Usage   :
 Function: Adds a value to the annotation collection under the specified
           key. Note that this is not a public method.
 Returns : 
 Args    : key (a string), value(s) (one or more scalars)


=cut

sub _add_annotation_value{
    my $self = shift;
    my $key  = shift;

    foreach my $val (@_) {
	$self->annotation->add_Annotation(
			Bio::Annotation::SimpleValue->new(-tagname => $key,
							  -value => $val)
					  );
    }
}

=head2 _get_annotation_values

 Title   : _get_annotation_values
 Usage   :
 Function: Gets the values of a specific annotation as identified by the
           key from the annotation collection. Note that this is not a
           public method.
 Example :
 Returns : an array of strings
 Args    : the key (a string)


=cut

sub _get_annotation_values{
    my $self = shift;

    return map { $_->value(); } $self->annotation->get_Annotations(shift);
}

#
##
### Deprecated methods kept for ease of transition
##
#

sub keywords {
    my $self = shift;

    # have we been called in set mode?
    if(@_) {
	# yes; translate to the new API
	foreach my $kwd (@_) {
	    $self->add_keyword(split(/\s*;\s*/,$kwd));
	}
    } else {
	# no; translate read-only to the new API
	return join("; ",$self->get_keywords());
    }
}

sub each_date {
   my ($self) = @_;
   $self->warn("Deprecated method... please use get_dates");
   return $self->get_dates;
}


sub each_secondary_accession {
   my ($self) = @_;
   $self->warn("each_secondary_accession - deprecated method. use get_secondary_accessions");
   return $self->get_secondary_accessions;

}

sub sv {
   my ($obj,$value) = @_;
   $obj->warn("sv - deprecated method. use seq_version");
   $obj->seq_version($value);
}


1;
