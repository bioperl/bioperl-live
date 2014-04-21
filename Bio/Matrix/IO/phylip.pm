#
# BioPerl module for Bio::Matrix::IO::phylip
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Matrix::IO::phylip - A parser for PHYLIP distance matricies

=head1 SYNOPSIS

  use Bio::Matrix::IO;
  my $parser = Bio::Matrix::IO->new(-format   => 'phylip',
                                   -file     => 't/data/phylipdist.out');
  my $matrix = $parser->next_matrix;

=head1 DESCRIPTION

This is a parser for PHYLIP distance matrix output.

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
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Matrix::IO::phylip;
use vars qw($DEFAULTPROGRAM);
use strict;

$DEFAULTPROGRAM = 'phylipdist';

use Bio::Matrix::PhylipDist;

use base qw(Bio::Matrix::IO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Matrix::IO::phylip->new();
 Function: Builds a new Bio::Matrix::IO::phylip object 
 Returns : an instance of Bio::Matrix::IO::phylip
 Args    :


=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  my ($prog) = $self->_rearrange([qw(PROGRAM)], @args);
  $self->{'_program'} = $prog || $DEFAULTPROGRAM;
  return $self;
}


=head2 next_matrix

 Title   : next_matrix
 Usage   : my $matrix = $parser->next_matrix
 Function: Get the next result set from parser data
 Returns : L<Bio::Matrix::PhylipDist>
 Args    : none


=cut

sub next_matrix {
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
   foreach my $name(@names){
       my $j=0;
       foreach my $n(@names) {
	   $dist{$name}{$n} = [$i,$j];
	   $j++;
       }
       $i++;
   }
   my $matrix = Bio::Matrix::PhylipDist->new
       (-matrix_name => $self->{'_program'},
	-matrix      => \%dist,
	-names       => \@names,
	-values      => \@values);
    return $matrix;
}

=head2 write_matrix

 Title   : write_matrix
 Usage   : $matio->write_matrix($matrix)
 Function: Write out a matrix in the phylip distance format
 Returns : none
 Args    : L<Bio::Matrix::PhylipDist>


=cut

sub write_matrix {
    my ($self,@matricies) = @_;
    foreach my $matrix ( @matricies ) {
	my @names = @{$matrix->names};
	my @values = @{$matrix->_values};
	my %matrix = %{$matrix->_matrix};
	my $str;
	$str.= (" "x 4). scalar(@names)."\n";
	foreach my $name (@names){
	    my $newname = $name. (" " x (15-length($name)));
	    if( length($name) >= 15 ) { $newname .= " " }
	    $str.=$newname;
	    my $count = 0;
	    foreach my $n (@names){
		my ($i,$j) = @{$matrix{$name}{$n}};
		if($count < $#names){
		    $str.= $values[$i][$j]. "  ";
		}
		else {
		    $str.= $values[$i][$j];
		}
		$count++;
	    }
	    $str.="\n";
	}
	$self->_print($str);
    }
}


1;
