# BioPerl module for Bio::Tools::Phylo::Phylip::ProtDist
#
# Cared for by Shawn Hoon <shawnh@fugu-sg.org>
#
# Copyright Shawn Hoon
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Phylo::Phylip::ProtDist - DESCRIPTION of Object

=head1 SYNOPSIS

    use Bio::Tools::Phylo::Phylip::ProtDist;
    my $parser = new Bio::Tools::Phylo::Phylip::ProtDist(-file => 'outfile');
    while( my $result = $parser->next_matrix) {

    }

=head1 DESCRIPTION

A parser for ProtDist  output into a Bio::Matrix::PhylipDist object

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Shawn Hoon 

Email shawnh@fugu-sg.org 

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Phylo::Phylip::ProtDist;
use vars qw(@ISA);
use strict;

use Bio::Root::Root;
use Bio::Matrix::PhylipDist;

use Bio::Root::IO;

@ISA = qw(Bio::Root::Root Bio::Root::IO );

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::Phylo::Phylip::ProtDist();
 Function: Builds a new Bio::Tools::Phylo::Phylip::ProtDist object 
 Returns : Bio::Tools::ProtDist
 Args    : -fh/-file => $val, # for initing input, see Bio::Root::IO


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->_initialize_io(@args);

  return $self;
}

=head2 next_result

 Title   : next_result
 Usage   : my $matrix = $parser->next_result
 Function: Get the next result set from parser data
 Returns : L<Bio::Matrix::PhylipDist>
 Args    : none


=cut

sub next_matrix{
   my ($self) = @_;
   my @names;
   my @values;
   my $entry;
   while ($entry=$self->_readline) {
       if($#names >=0 && $entry =~/^\s+\d+$/){
           last;
       }
       elsif($entry=~/^\s+\d+\n$/){
           next;
       }
       my ($n,@line) = split( /\s+/,$entry);
       push @names, $n;
       push @values, [@line];
    }
    $#names>=0 || return;
    my %dist;
    my $i=0;
    foreach my $name(@names){
      my $j=0;
      foreach my $n(@names) {
        $dist{$name}{$n} = [$i,$j];
        $j++;
      }
     $i++;
    }
    my $matrix = Bio::Matrix::PhylipDist->new(-matrix=>\%dist,
                                              -names =>\@names,
                                              -values=>\@values);
    return $matrix;
}

1;
