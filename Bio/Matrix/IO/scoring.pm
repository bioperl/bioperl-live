#
# BioPerl module for Bio::Matrix::IO::scoring
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Matrix::IO::scoring - A parser for PAM/BLOSUM matricies

=head1 SYNOPSIS

  use Bio::Matrix::IO;
  my $parser = Bio::Matrix::IO->new(-format => 'scoring',
                                   -file   => 'BLOSUM50');
  my $matrix = $parser->next_matrix;

=head1 DESCRIPTION

Describe the object here

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

Email jason-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Matrix::IO::scoring;
use strict;

use Bio::Matrix::Scoring;
use base qw(Bio::Matrix::IO);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Matrix::IO::scoring->new();
 Function: Builds a new Bio::Matrix::IO::scoring object 
 Returns : an instance of Bio::Matrix::IO::scoring
 Args    :


=cut

=head2 next_matrix

 Title   : next_matrix
 Usage   : my $matrux = $parser->next_matrix
 Function: parses a scoring matrix (BLOSUM,PAM styles) 
 Returns : L<Bio::Matrix::Scoring>
 Args    : none


=cut

sub next_matrix{
   my ($self) = @_;
   local ($_);
   my (@matrix,@cols,@rows,%extras,$inmatrix);
   while( defined ( $_ = $self->_readline ) ) {
       next if ( /^\s*$/);
       if( /^\#/ ) {
	   if( $inmatrix ) { 
	       $self->_pushback($_);
	       last;
	   }
	   if( m/Entropy\s+\=\s+(\S+)\,\s+
	       Expected\s+\=\s+(\S+)/ox ) {
	       $extras{'-entropy'} = $1;
	       $extras{'-expected'} = $2;
	   } elsif ( m/Expected\s+score\s+\=\s+(\S+)\,
		     \s+Entropy\s+\=\s+(\S+)/xo ){
	       $extras{'-entropy'} = $2;
	       $extras{'-expected'} = $1;
	   } elsif( m/(PAM\s+\d+)\s+substitution.+
		    scale\s+\=\s+(\S+)\s+\=\s+(\S+)/ox ) {
	       $extras{'-matrix_name'} = $1;
	       $extras{'-scale'}       = $2;	       
	       $extras{'-scale_value'} = $3;
	   } elsif( /Blocks Database\s+\=\s+(\S+)/o ) {
	       $extras{'-database'} = $1;
	   } elsif( m/(\S+)\s+Bit\s+Units/ox ) {
	       $extras{'-scale'} = $1;
	   } elsif( m/Lowest score\s+\=\s+(\S+)\,\s+
		    Highest score\s+\=\s+(\S+)/ox ) {
	       $extras{'-lowest_score'} = $1;
	       $extras{'-highest_score'} = $2;
	   } elsif( m/(Lambda)\s+\=\s+(\S+)\s+bits\,
		    \s+(H)\s+\=\s+(\S+)/ox ) {
	       # This is a DNA matrix
	       $extras{$1} = $2;
	       $extras{$3} = $4;
	   }	       
       } elsif( s/^\s+(\S+)/$1/ ) {
	   @cols = split;
	   if( $cols[0] ne 'A' ) {
	       $self->warn("Unrecognized first line of matrix, we might not have parsed it correctly");
	   }
	   $inmatrix = 1;
       } elsif( $inmatrix ) {
	   if( ! /^(\S+)/ ) { $inmatrix = 0; next }
	   my ($rowname,@row) = split;
	   push @rows, $rowname;
	   push @matrix, [@row];
       } else { 
	   print;
       }
   }
   my $matrix = Bio::Matrix::Scoring->new(-values     => \@matrix,
					 -rownames   => \@rows,
					 -colnames   => \@cols,
					 %extras);
}

=head2 write_matrix

 Title   : write_matrix
 Usage   : $matio->write_matrix($matrix)
 Function: Write out a matrix in the BLOSUM/PAM format
 Returns : none
 Args    : L<Bio::Matrix::Scoring>


=cut

sub write_matrix{
   my ($self,@args) = @_;
   $self->warn("cannot actually use this function yet - it isn't finished");
   return;
}


1;
