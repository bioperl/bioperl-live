
#
# BioPerl module for Bio::Seq::RichSeq
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Seq::RichSeq - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Seq::RichSeq;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Seq;
use Bio::Seq::RichSeqI

@ISA = qw(Bio::Seq Bio::Seq::RichSeqI);


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

    $self->{'date'} = [];
    $self->{'secondary_accession'} = [];

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
      $obj->{'division'} = $value;
    }
    return $obj->{'division'};

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
      $obj->{'molecule'} = $value;
    }
    return $obj->{'molecule'};

}

=head2 add_date

 Title   : add_date
 Usage   : $self->add_date($ref)
 Function: adds a date
 Example :
 Returns : 
 Args    :


=cut

sub add_date {
   my ($self) = shift;
   foreach my $dt ( @_ ) {
       push(@{$self->{'date'}},$dt);
   }
}

=head2 get_dates

 Title   : get_dates
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_dates{
   my ($self,@args) = @_;
   return @{$self->{'date'}}; 

}


=head2 accession

 Title   : accession
 Usage   : $obj->accession($newval)
 Function: Whilst the underlying sequence object does not 
           have an accession, so we need one here. Won't stay
           when we do the reimplementation.
 Example : 
 Returns : value of accession
 Args    : newvalue (optional)


=cut
#'

sub accession {
   my ($obj,$value) = @_;
   return $obj->accession_number($value);
}

=head2 add_secondary_accession

 Title   : add_secondary_accession
 Usage   : $self->add_domment($ref)
 Function: adds a secondary_accession
 Example :
 Returns : 
 Args    :


=cut

sub add_secondary_accession {
   my ($self) = shift;
   foreach my $dt ( @_ ) {
       push(@{$self->{'secondary_accession'}},$dt);
   }
}

=head2 get_secondary_accessions

 Title   : get_secondary_accessions
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_secondary_accessions{
   my ($self,@args) = @_;
   return @{$self->{'secondary_accession'}}; 
}

=head2 seq_version

 Title   : seq_version
 Usage   : $obj->seq_version($newval)
 Function: 
 Example : 
 Returns : value of seq_version
 Args    : newvalue (optional)


=cut

sub seq_version{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'seq_version'} = $value;
    }
    return $obj->{'seq_version'};

}


=head2 keywords

 Title   : keywords
 Usage   : $obj->keywords($newval)
 Function: 
 Returns : value of keywords
 Args    : newvalue (optional)


=cut

sub keywords {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'keywords'} = $value;
    }
    return $obj->{'keywords'};

}

=head2 Deprecated methods kept for ease of transtion

=head2 each_date

 Title   : each_date
 Usage   : foreach $dt ( $self->each_date() )
 Function: gets an array of dates
 Example :
 Returns : 
 Args    :


=cut

sub each_date {
   my ($self) = @_;
   $self->warn("Deprecated method... please use get_dates");
   return $self->get_dates;
}


=head2 each_secondary_accession

 Title   : each_secondary_accession
 Usage   : foreach $dt ( $self->each_secondary_accession() )
 Function: gets an array of secondary_accessions
 Example :
 Returns : 
 Args    :


=cut

sub each_secondary_accession {
   my ($self) = @_;
   $self->warn("each_secondary_accession - deprecated method. use get_secondary_accessions");
   return $self->get_secondary_accessions;

}


=head2 sv

 Title   : sv
 Usage   : $obj->sv($newval)
 Function: 
 Returns : value of sv
 Args    : newvalue (optional)


=cut

sub sv {
   my ($obj,$value) = @_;
   $obj->warn("sv - deprecated method. use seq_version");
   $obj->seq_version($value);
}


1;






