#
# BioPerl module for Bio::SeqFeature::Primer
#
# Cared for by Chad Matsalla
#chad
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature::Primer - Primer Generic SeqFeature

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Chad Matsalla

Chad Matsalla <lt>bioinformatics1@dieselwurks.comE<gt>

=head1 DEVELOPERS

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature::Primer;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::SeqFeatureI;
use Bio::SeqFeature::Generic;
     # I am really not sure which of these I need so I included them at this time
use Bio::Annotation;
use Bio::Location::Simple;
use Bio::Seq;
use Dumpvalue qw(dumpValue);

my $dumper = new Dumpvalue();


@ISA = qw(Bio::Root::Root Bio::SeqFeature::Generic);



=head2 new()

 Title   : new()
 Usage   :
 Function:
 Example :
 Returns :
 Args    : 
Devel notes: I think that I want to accept a hash

=cut


sub new {
     my ($class, @args) = @_;  
     my %arguments = @args;
     my $self = $class->SUPER::new(@args);
          # these are from generic.pm, with which i started
     $self->{'_parse_h'}       = {};
     $self->{'_gsf_tag_hash'}  = {}; 
          # things that belong with the primer
     my ($sequence, $primer_sequence_id, $id) = 
     $self->{'_gsf_seqname'} = $self->{primer_sequence_id};
          # i am going to keep an array of the things that have been passed
          # into the object on construction. this will aid retrieval of these
          # things later
     foreach my $argument (sort keys %arguments) {
          if ($argument eq "-SEQUENCE" || $argument eq "-sequence") {
               if (ref(%arguments->{$argument}) eq "Bio::Seq") {
                    $self->{seq} = %arguments->{$argument};
               }
               else {
                    $self->{seq} = new Bio::Seq( -seq => %arguments->{$argument},
                                                  -id => %arguments->{-id});
                    # print("Created new bioseq with seq ".%arguments->{$argument}." and name ".%arguments->{-id}."\n");
               }
               push @{$self->{arguments}},$argument;
          }
          else {
               (my $fixed = $argument) =~ s/-//;
               $self->{$fixed} = %arguments->{$argument};
               push @{$self->{arguments}},$fixed;
          }
     }
     if (!$self->{seq}) {
          $self->throw("You must pass in a sequence to construct this object.");
     }
               # a bunch of things now need to be set for this SeqFeature
               # things like:
                    # TARGET=513,26
                    # PRIMER_FIRST_BASE_INDEX=1
                    # PRIMER_LEFT=484,20
          return $self;
}


=head2 seq()

 Title   : seq()
 Usage   : $seq = $primer->seq();
 Function: Return the _entire_ sequence associated with this Primer. 
 Returns : A Bio::Seq object
 Args    : None.
Develper Note: Do you want to be able to set the sequence associated with this
     SeqFeature?

=cut

sub seq {
     my $self = shift;
     return $self->{seq};
}



=head2 all_tags()

 Title   : all_tags()
 Usage   : @tags = $primer->all_tags();
 Function: Return a list of tag names for this Primer.
 Returns : An array of strings representing the names of tags in this Primer
 Args    : None.
 Notes   : When the Bio::SeqFeature::Primer object is created, the user can
     pass in an arbitrary hash containing key->value pairs. This is allowed
     because I didn't want to assume that the user was trying to model a
     primer3 construct. 

=cut

sub all_tags {
     my $self = shift;
     my @tags = sort keys %{$self->{tags}};
     return @tags;
}


=head2 primary_tag()

 Title   : primary_tag()
 Usage   : $tag = $feature->primary_tag();
 Function: Returns the string "Primer"
 Returns : A string.
 Args    : None.

=cut

sub primary_tag {
     return "Primer";
}

=head2 source_tag()

 Title   : source_tag()
 Usage   : $tag = $feature->source_tag();
 Function: Returns the source of this tag.
 Returns : A string.
 Args    : If an argument is provided, the source of this SeqFeature
     is set to that argument.

=cut

sub source_tag {
     my ($self,$insource) = @_;
     if ($insource) { $self->{source} = $insource; }
     return $self->{source};
}

=head2 has_tag()

 Title   : has_tag()
 Usage   : $true_or_false = $feature->has_tag('MELTING_TEMPERATURE');
 Function: Does this SeqFeature have this tag?
 Returns : TRUE or FALSE
 Args    : A string.

=cut

sub has_tag {
     my ($self,$tagname) = @_;
     if ($self->{tags}->{$tagname}) { return "TRUE"; }
     return { "FALSE" };
}

=head2 each_tag_value()

 Title   : each_tag_value()
 Usage   : $tag = $feature->each_tag_value('MELTING_TEMPERATURE');
 Function: Returns the value of this tag.
 Returns : Unknown. Whatever the value of the given tag was.
 Args    : None.

=cut

sub each_tag_value {
     my ($self,$tagname) = @_;
     return $self->{tags}->{$tagname};
}

=head2 location()

 Title   : location()
 Usage   : $tag = $feature->location();
 Function: Returns a location object suitable for identifying location of
     feature on sequence or parent feature  
 Returns : A Bio::LocationI object.
 Args    : None.

=cut

sub location {
     my $self = shift;
     $self->warn("Chad has not written the code for this yet.");
}



sub start {
}

sub end {
}

sub strand {
}

sub display_id {

}



sub overlaps {

}

sub contains {
}

sub equals {

}

sub intersection {

}

sub union {

}




1;
