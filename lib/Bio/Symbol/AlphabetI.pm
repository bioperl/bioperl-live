#
# BioPerl module for Bio::Symbol::AlphabetI
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

Bio::Symbol::AlphabetI - A Symbol Alphabet

=head1 SYNOPSIS

    # get a Bio::Symbol::AlphabetI object somehow
    my @symbols = $alphabet->symbols;
    my @subalphas = $alphabet->alphabets;
    if( $alphabet->contains($symbol) ) {
	# do something
    }

=head1 DESCRIPTION

Alphabet contains set of symbols, which can be concatenated to form
symbol lists. Sequence string, for example, is stringified
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

package Bio::Symbol::AlphabetI;
use strict;
use Bio::Root::RootI;

=head2 AlphabetI Interface methods

=cut

=head2 symbols

 Title   : symbols
 Usage   : my @symbols = $alphabet->symbols();
 Function: Get/Set Symbol list for an alphabet
           List of symbols, which make up this alphabet.
 Returns : Array of L<Bio::Symbol::SymbolI> objects
 Args    : (optional) Array of L<Bio::Symbol::SymbolI> objects

=cut

sub symbols{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

=head2 alphabets

 Title   : alphabets
 Usage   : my @alphabets = $alphabet->alphabets();
 Function: Get/Set Sub Alphabet list for an alphabet 
           Sub-alphabets. E.g. codons made from DNAxDNAxDNA alphabets
 Returns : Array of L<Bio::Symbol::AlphabetI> objects
 Args    : (optional) Array of L<Bio::Symbol::AlphabetI> objects

=cut

sub alphabets{
    my ($self,@args) = @_;
    $self->throw_not_implemented();
}

=head2 contains

 Title   : contains
 Usage   : if($alphabet->contains($symbol)) { }
 Function: Tests of Symbol is contained in this alphabet
 Returns : Boolean
 Args    : L<Bio::Symbol::SymbolI>

=cut

sub contains{
   my ($self,@args) = @_;
   $self->throw_not_implemented();
}

# Other methods from BSANE - not sure if we will implement here or only in
# BioCORBA implementation

# Resolve symbols from the token string.
#    SymbolList to_symbol(in string tokens) raises ( IllegalSymbolException) ;

# Convinience method, which returns gap symbol that do not
# match with any other symbols in the alphabet.
#   Symbol get_gap_symbol() raises ( DoesNotExist) ;


# Returns a ambiguity symbol, which represent list of
# symbols. All symbols in a list must be members of
# this alphabet otherwise IllegalSymbolException is
# thrown.
# Symbol get_ambiguity( in SymbolList symbols) raises( IllegalSymbolException);


#  Returns a Symbol, which represents ordered list of symbols
#  given as a parameter. Each symbol in the list must be member of
#  different sub-alphabet in the order defined by the alphabets
#  attribute. For example, codons can be represented by a compound
#  Alphabet of three DNA Alphabets, in which case the get_symbol(
#  SymbolList[ a,g,t]) method of the Alphabet returns Symbol for
#  the codon agt.<p>

#  IllegalSymbolException is raised if members of symbols
#  are not Symbols over the alphabet defined by
#  get_alphabets()-method
#  Symbol get_symbol(in SymbolList symbols) raises(IllegalSymbolException) ;

1;
