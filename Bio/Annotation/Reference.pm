#
# BioPerl module for Bio::Annotation::Reference
#
# Cared for by Ewan Birney <pfam@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Annotation::Reference - Literature Reference Block

=head1 SYNOPSIS

    $annotation = $seq->annotation();
    foreach $ref ( $annotation->each_Reference() ) {
	print "Title: ",$ref->title,"\n";
        print "Authors:",$ref->authors,"\n";
        print "Location:",$ref->location,"\n";
    }

    

=head1 DESCRIPTION

Object which presents a literature reference. This is considered to be
a specialised form of database link. The literature reference contains
the title, authors and location lines as text: no attempt is made to
parse these into more computationally accessibly aspects, as it
becomes a nightmare. A slot for the medline accession number is also
provided.

As most people will realise this object is really designed to
supported literature references made in flat files. For more complex
representations, including abstract and possibly full text a richer
object that can somehow be easily connected (? inheritence or
composition) to this object would be ideal. Please air your views on
the guts list and if possible, nominate yourself as someone who can do
this.

=head1 CONTACT

Ewan Birney <birney@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Annotation::Reference;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Annotation::DBLink;

@ISA = qw(Bio::Annotation::DBLink);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);

# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 authors

 Title   : authors
 Usage   : $self->authors($newval)
 Function: Gives the author line. No attempt is made to parse the author line
 Example : 
 Returns : value of authors
 Args    : newvalue (optional)

=cut

sub authors{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'authors'} = $value;
    }
    return $self->{'authors'};

}

=head2 location

 Title   : location
 Usage   : $self->location($newval)
 Function: Gives the location line. No attempt is made to parse the location line
 Example : 
 Returns : value of location
 Args    : newvalue (optional)

=cut

sub location{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'location'} = $value;
    }
    return $self->{'location'};

}

=head2 title

 Title   : title
 Usage   : $self->title($newval)
 Function: Gives the title line (if exists)
 Example : 
 Returns : value of title
 Args    : newvalue (optional)

=cut

sub title{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'title'} = $value;
    }
    return $self->{'title'};

}

=head2 medline

 Title   : medline
 Usage   : $self->medline($newval)
 Function: Gives the medline number
 Example : 
 Returns : value of medline
 Args    : newvalue (optional)

=cut

sub medline{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'medline'} = $value;
    }
    return $self->{'medline'};

}

=head2 database

 Title   : database
 Usage   :
 Function: Overrides DBLink database to be hard coded to 'MEDLINE'
 Example :
 Returns : 
 Args    :

=cut

sub database{
   my ($self) = @_;

   return 'MEDLINE';
}

=head2 primary_id

 Title   : primary_id
 Usage   :
 Function: Overrides DBLink primary_id to provide medline number
 Example :
 Returns : 
 Args    :

=cut

sub primary_id{
   my ($self, @args) = @_;

   return $self->medline(@args);
}

=head2 start

 Title   : start
 Usage   : $self->start($newval)
 Function: Gives the reference start base
 Example : 
 Returns : value of start
 Args    : newvalue (optional)

=cut

sub start{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'start'} = $value;
    }
    return $self->{'start'};

}
=head2 end

 Title   : end
 Usage   : $self->end($newval)
 Function: Gives the reference end base 
 Example : 
 Returns : value of end
 Args    : newvalue (optional)

=cut

sub end{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'end'} = $value;
    }
    return $self->{'end'};

}

=head2 rp

 Title   : rp
 Usage   : $self->rp($newval)
 Function: Gives the RP line. No attempt is made to parse this line.
 Example : 
 Returns : value of rp
 Args    : newvalue (optional)

=cut

sub rp{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'rp'} = $value;
    }
    return $self->{'rp'};

}

1;

