#
# BioPerl module for Bio::Tools::HMMER::Set
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::HMMER::Set - Set of identical domains from HMMER matches

=head1 SYNOPSIS

    # get a Set object probably from the results object
    print "Bits score over set ",$set->bits," evalue ",$set->evalue,"\n";

    foreach $domain ( $set->each_Domain ) {
	print "Domain start ",$domain->start," end ",$domain->end,"\n";
    }

=head1 DESCRIPTION

Represents a set of HMMER domains hitting one sequence. HMMER reports two
different scores, a per sequence total score (and evalue) and a per
domain score and evalue. This object represents a collection of the same
domain with the sequence bits score and evalue. (these attributes are also
on the per domain scores, which you can get there).

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
the bugs and their resolution.Bug reports can be submitted via the
web: 

 https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Ewan Birney

Email birney-at-ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::HMMER::Set;
use strict;

use Bio::Tools::HMMER::Domain;

use base qw(Bio::Root::Root);

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
    my ($name,$acc,$desc) = $self->_rearrange([qw(NAME ACCESSION DESC)],
					      @args);
    $name && $self->name($name);
    $acc  && $self->accession($acc);
    $desc && $self->desc($desc);


    $self->{'domains'} = [];
    $self->{'domainnames'} = {};
    return $self;
}

=head2 add_Domain

 Title   : add_Domain
 Usage   : $set->add_Domain($domain)
 Function: adds the domain to the list
 Returns : nothing
 Args    : A Bio::Tools::HMMER::Domain object

=cut

sub add_Domain{
   my ($self,$domain) = @_;


   if( ! defined $domain || ! $domain->isa("Bio::Tools::HMMER::Domain") ) {
       $self->throw("[$domain] is not a Bio::Tools::HMMER::Domain. aborting");
   }
   return if $self->{'domainnames'}->{$domain->get_nse}++;
   push(@{$self->{'domains'}},$domain);

}

=head2 each_Domain

 Title   : each_Domain
 Usage   : foreach $domain ( $set->each_Domain() )
 Function: returns an array of domain objects in this set
 Returns : array
 Args    : none


=cut

sub each_Domain{
   my ($self,@args) = @_;

   return @{$self->{'domains'}};
}

=head2 name

 Title   : name
 Usage   : $obj->name($newval)
 Function:
 Example :
 Returns : value of name
 Args    : newvalue (optional)


=cut

sub name{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'name'} = $value;
    }
    return $obj->{'name'};

}

=head2 desc

 Title   : desc
 Usage   : $obj->desc($newval)
 Function:
 Example :
 Returns : value of desc
 Args    : newvalue (optional)

=cut

sub desc{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'desc'} = $value;
    }
    return $self->{'desc'};

}

=head2 accession

 Title   : accession
 Usage   : $obj->accession($newval)
 Function:
 Example :
 Returns : value of accession
 Args    : newvalue (optional)


=cut

sub accession{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'accession'} = $value;
    }
    return $self->{'accession'};
}


=head2 bits

 Title   : bits
 Usage   : $obj->bits($newval)
 Function:
 Example :
 Returns : value of bits
 Args    : newvalue (optional)


=cut

sub bits{
   my ($obj,$value) = @_;

   if( defined $value) {
      $obj->{'bits'} = $value;
    }
    return $obj->{'bits'};

}

=head2 evalue

 Title   : evalue
 Usage   : $obj->evalue($newval)
 Function:
 Example :
 Returns : value of evalue
 Args    : newvalue (optional)


=cut

sub evalue{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'evalue'} = $value;
    }
    return $obj->{'evalue'};

}


sub addHMMUnit {
    my $self = shift;
    my $unit = shift;

    $self->warn("Using old addHMMUnit call on Bio::Tools::HMMER::Set. Should replace with add_Domain");
    return $self->add_Domain($unit);
}

sub eachHMMUnit {
    my $self = shift;
    $self->warn("Using old eachHMMUnit call on Bio::Tools::HMMER::Set. Should replace with each_Domain");
    return $self->each_Domain();
}

1;  # says use was ok
__END__

