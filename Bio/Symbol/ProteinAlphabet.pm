#
# BioPerl module for Bio::Symbol::ProteinAlphabet
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason@bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Symbol::ProteinAlphabet - A ready made Protein alphabet

=head1 SYNOPSIS

    use Bio::Symbol::ProteinAlphabet;
    my $alpha = Bio::Symbol::ProteinAlphabet->new();
    foreach my $symbol ( $alpha->symbols ) {
	print "symbol is $symbol\n";
    }

=head1 DESCRIPTION

This object builds an Alphabet with Protein symbols.

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

=head1 AUTHOR - Jason Stajich

Email jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Symbol::ProteinAlphabet;
use strict;

use Bio::Symbol::Symbol;
use Bio::Tools::IUPAC;
use Bio::SeqUtils;

use base qw(Bio::Symbol::Alphabet);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Symbol::ProteinAlphabet->new();
 Function: Builds a new Bio::Symbol::ProteinAlphabet object 
 Returns : Bio::Symbol::ProteinAlphabet
 Args    :


=cut

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);  
  my %aa = Bio::SeqUtils->valid_aa(2);
  my %codes = Bio::Tools::IUPAC->iupac_iup();
  my %symbols;
  my @left;
  
  foreach my $let ( keys %codes  ) {  
      if( scalar @{$codes{$let}} != 1) { push @left, $let; next; }
      $symbols{$let} = Bio::Symbol::Symbol->new(-name => $aa{$let},
					       -token => $let);      
  }
  foreach my $l ( @left ) {
      my @subsym;
      foreach my $sym ( @{$codes{$l}} ) {
	  push @subsym, $symbols{$sym};
      }
      my $alpha = Bio::Symbol::Alphabet->new(-symbols => \@subsym);
      $symbols{$l} = Bio::Symbol::Symbol->new(-name => $aa{$l},
					       -token => $l,
					       -matches => $alpha,
					       -symbols => \@subsym);
  }
  
  $self->symbols(values %symbols); 
  return $self;
}


1;
