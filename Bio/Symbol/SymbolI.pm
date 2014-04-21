#
# BioPerl module for Bio::Symbol::SymbolI
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

Bio::Symbol::SymbolI - Interface for a Symbol

=head1 SYNOPSIS

    # get a Bio::Symbol::SymbolI object somehow

    my ($name,$token) = ($symbol->name, $symbol->token);
    my @symbols       = $symbol->symbols;
    my $matches       = $symbol->matches;

=head1 DESCRIPTION

Symbol represents a single token in the sequence. Symbol can have
multiple synonyms or matches within the same Alphabet, which
makes possible to represent ambiguity codes and gaps.

Symbols can be also composed from ordered list other symbols. For
example, codons can be represented by single Symbol using a
compound Alphabet made from three DNA Alphabets.

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


package Bio::Symbol::SymbolI;
use strict;
use base qw(Bio::Root::RootI);

=head2 Bio::Symbol::SymbolI interface methods

=cut

=head2 name

 Title   : name
 Usage   : my $name = $symbol->name();
 Function: Get/Set Descriptive name for Symbol
 Returns : string
 Args    : (optional) string

=cut

sub name{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 token

 Title   : token
 Usage   : my $token = $self->token();
 Function: Get/Set token for this symbol
 Example : Letter A,C,G,or T for a DNA alphabet Symbol
 Returns : string
 Args    : (optional) string

=cut

sub token{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 symbols

 Title   : symbols
 Usage   : my @symbols = $self->symbols();
 Function: Get/Set Symbols this Symbol is composed from
 Example : A codon is composed of 3 DNA symbols
 Returns : Array of Bio::Symbol::SymbolI objects
 Args    : (optional) Array of Bio::Symbol::SymbolI objects


=cut

sub symbols{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 matches

 Title   : matches
 Usage   : my $matchalphabet = $symbol->matches();
 Function: Get/Set (Sub) alphabet of symbols matched by this symbol
           including the symbol itself (i.e. if symbol is DNA
           ambiguity code W then the matches contains symbols for W
           and T)
 Returns : Bio::Symbol::AlphabetI
 Args    : (optional) Bio::Symbol::AlphabetI

=cut

sub matches{
   my ($self,@args) = @_;
   $self->throw_not_implemented();   
}

=head2 equals

 Title   : equals
 Usage   : if( $symbol->equals($symbol2) ) { }
 Function: Tests if a symbol is equal to another 
 Returns : Boolean
 Args    : Bio::Symbol::SymbolI

=cut

sub equals{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

1;
