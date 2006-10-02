# $Id$
#
# BioPerl module for Bio::Symbol::DNAAlphabet
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Symbol::DNAAlphabet - A ready made DNA alphabet

=head1 SYNOPSIS

    use Bio::Symbol::DNAAlphabet;
    my $alpha = new Bio::Symbol::DNAAlphabet();
    foreach my $symbol ( $alpha->symbols ) {
	print "symbol is $symbol\n";
    }

=head1 DESCRIPTION

This object builds an Alphabet with DNA symbols.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Symbol::DNAAlphabet;
use strict;

use Bio::Symbol::Symbol;
use Bio::Tools::IUPAC;

use base qw(Bio::Symbol::Alphabet);

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Symbol::DNAAlphabet();
 Function: Builds a new Bio::Symbol::DNAAlphabet object 
 Returns : Bio::Symbol::DNAAlphabet
 Args    :


=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);  
  my %alphabet = Bio::Tools::IUPAC::iupac_iub();
  my %symbols;
  foreach my $let ( keys %alphabet ) {
      next unless @{$alphabet{$let}} == 1 || $let eq 'U';
      $symbols{$let} = new Bio::Symbol::Symbol(-name => $let,
					       -token => $let);      
  }
  
  foreach my $let ( keys %alphabet ) {
      next if( $symbols{$let} || $let eq 'U');
      my @subsymbols;
      
      foreach my $sublet ( @{$alphabet{$let}} ) {
	  push @subsymbols, $symbols{$sublet};
      }
      my $alpha = new Bio::Symbol::Alphabet(-symbols => \@subsymbols);
      $symbols{$let} = new Bio::Symbol::Symbol(-name    => $let,
					       -token   => $let,
					       -matches => $alpha,  
					       -symbols => \@subsymbols); 
  }
  
  $self->symbols(values %symbols); 
  return $self;
}


1;
