# BioPerl module for Bio::Tools::Phylo::Phylip::ProtDist
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Shawn Hoon <shawnh@fugu-sg.org>
#
# Copyright Shawn Hoon
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Phylo::Phylip::ProtDist - parser for ProtDist output

=head1 SYNOPSIS

    use Bio::Tools::Phylo::Phylip::ProtDist;
    my $parser = Bio::Tools::Phylo::Phylip::ProtDist->new(-file => 'outfile');
    while( my $result = $parser->next_matrix) {
      # do something with it
    }

=head1 DESCRIPTION

A parser for ProtDist output into a L<Bio::Matrix::PhylipDist> object.
See also L<Bio::Matrix::IO::phylip> this module may go away.

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

=head1 AUTHOR - Shawn Hoon 

Email shawnh@fugu-sg.org 

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::Phylo::Phylip::ProtDist;
use strict;

use Bio::Matrix::PhylipDist;


use base qw(Bio::Root::Root Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Phylo::Phylip::ProtDist->new();
 Function: Builds a new Bio::Tools::Phylo::Phylip::ProtDist object 
 Returns : Bio::Tools::ProtDist
 Args    : -fh/-file => $val, # for initing input, see Bio::Root::IO
           -program  => 'programname' # name of the program

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->_initialize_io(@args);
  my ($prog) = $self->_rearrange([qw(PROGRAM)], @args);
  $self->{'_program'} = $prog;
  return $self;
}

=head2 next_matrix

 Title   : next_matrix
 Usage   : my $matrix = $parser->next_matrix
 Function: Get the next result set from parser data
 Returns : L<Bio::Matrix::PhylipDist>
 Args    : none


=cut

sub next_matrix{
   my ($self) = @_;
   my @names;
   my @values;
   my $entry;
   my $size = 0;
   while ($entry=$self->_readline) {
       if($#names >=0 && $entry =~/^\s+\d+\n$/){
	   $self->_pushback($entry);
	   last;
       } elsif($entry=~/^\s+(\d+)\n$/){	   
	   $size = $1;
	   next;
       } elsif( $entry =~ s/^\s+(\-?\d+\.\d+)/$1/ ) {
	   my (@line) = split( /\s+/,$entry);
	   push @{$values[-1]}, @line;
	   next;
       }
       my ($n,@line) = split( /\s+/,$entry);
       
       push @names, $n;
       push @values, [@line];
   }
   if( scalar @names != $size ) {
       $self->warn("The number of entries ".(scalar @names).
		   " is not the same $size");
   }
   $#names>=0 || return;
   my %dist;
   my $i=0;
   for my $name (@names){
       my $j=0;
       for my $n (@names) {
	   $dist{$name}{$n} = [$i,$j];
	   $j++;
       }
       $i++;
   }
   return Bio::Matrix::PhylipDist->new(-program => $self->{'_program'},
					     -matrix  => \%dist,
					     -names   => \@names,
					     -values  => \@values);
}

1;
