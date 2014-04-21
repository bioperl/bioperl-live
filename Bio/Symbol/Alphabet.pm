#
# BioPerl module for Bio::Symbol::Alphabet
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

Bio::Symbol::Alphabet - BSANE/BioCORBA compliant symbol list alphabet

=head1 SYNOPSIS

  {
      my $alphabet = Bio::Symbols::Alphabet->new(-symbols => [ @s ],
  						-subalphabets => [ @alphas ] );

      my @symbols = $alphabet->symbols;
      my @subalphas = $alphabet->alphabets;
      if( $alphabet->contains($symbol) ) {
  	  # do something
      }
  }

=head1 DESCRIPTION

Alphabet contains set of symbols, which can be concatenated to
form symbol lists. Sequence string, for example, is stringified
representation of the symbol list (tokens of symbols).

This module was implemented for the purposes of meeting the
BSANE/BioCORBA spec 0.3 only.

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


package Bio::Symbol::Alphabet;
use strict;

# Object preamble - inherits from Bio::Root::Root


use base qw(Bio::Root::Root Bio::Symbol::AlphabetI);

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Symbol::Alphabet->new();
 Function: Builds a new Bio::Symbol::Alphabet object 
 Returns : Bio::Symbol::Alphabet
 Args    : -symbols  => Array ref of Bio::Symbol::SymbolI objects
           -subalphas=> Array ref of Bio::Symbol::AlphabetI objects 
                        representing sub alphabets

=cut

sub new {
    my($class,@args) = @_;

    my $self = $class->SUPER::new(@args);
    $self->{'_symbols'} = [];
    $self->{'_alphabets'} = [];
    my ($symbols, $subalphas) = $self->_rearrange([qw(SYMBOLS SUBALPHAS)],
						  @args);

    defined $symbols && ref($symbols) =~ /array/i && $self->symbols(@$symbols);
    defined $subalphas && ref($subalphas) =~ /array/i && $self->alphabets(@$subalphas);
    return $self;
}

=head2 AlphabetI Interface methods

=cut

=head2 symbols

 Title   : symbols
 Usage   : my @symbols = $alphabet->symbols();
 Function: Get/Set Symbol list for an alphabet
           List of symbols, which make up this alphabet.
 Returns : Array of Bio::Symbol::SymbolI objects
 Args    : (optionalalphabets) Array of Bio::Symbol::SymbolI objects

=cut

sub symbols {
    my ($self,@args) = @_;
    if( @args ) { 
	$self->{'_symbols'} = [];
	foreach my $symbol ( @args ) {
	    if( ! defined $symbol || ! ref($symbol) || 
		! $symbol->isa('Bio::Symbol::SymbolI') ) {
		$self->warn("Did not provide a proper Bio::Symbol::SymbolI to method 'symbols' (got $symbol)");
	    } else { 
		push @{$self->{'_symbols'}}, $symbol;
	    }
	}
    }
    return @{$self->{'_symbols'}};
}

=head2 alphabets

 Title   : alphabets
 Usage   : my @alphabets = $alphabet->alphabets();
 Function: Get/Set Sub Alphabet list for an alphabet 
           Sub-alphabets. E.g. codons made from DNAxDNAxDNA alphabets
 Returns : Array of Bio::Symbol::AlphabetI objects
 Args    : (optional) Array of Bio::Symbol::AlphabetI objects

=cut

sub alphabets {
    my ($self,@args) = @_;
   if( @args ) { 
       $self->{'_alphabets'} = [];
       foreach my $alpha ( @args ) {
	   if( ! $alpha->isa('Bio::Symbol::AlphabetI') ) {
	       $self->warn("Did not provide a proper Bio::Symbol::AlphabetI to method 'alphabets' (got $alpha)");
	   } else { 
	       push @{$self->{'_alphabets'}}, $alpha;
	   }
       }
   }
    return @{$self->{'_alphabets'}};
}

=head2 contains

 Title   : contains
 Usage   : if($alphabet->contains($symbol)) { }
 Function: Tests of Symbol is contained in this alphabet
 Returns : Boolean
 Args    : Bio::Symbol::SymbolI

=cut

sub contains{
   my ($self,$testsymbol) = @_;
   foreach my $symbol ( $self->symbols ) {
       return 1 if( $symbol->equals($testsymbol) );
   }
   return 0;
}

1;
