
#
# BioPerl module for Bio::SeqUtils
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqUtils - Additional methods for PrimarySeq objects

=head1 SYNOPSIS

    # get a Bio::PrimarySeqI compliant object, $seq, somehow
    $util = new Bio::SeqUtils;
    $poplypeptide_3char = $util->seq3($seq);
    #or
    $poplypeptide_3char = Bio::SeqUtils->seq3($seq);

    #set the sequence string (stored in one char code in the object)
    Bio::SeqUtils->seq3($seq, $poplypeptide_3char);

=head1 DESCRIPTION

This class is a holder of methods that work on L<Bio::PrimarySeqI>
compliant sequence objects, e.g. L<Bio::PrimarySeq> and
L<Bio::Seq>. These methods are not part of the L<Bio::PrimarySeqI>
interface and should in general not essential to the primary function
of sequence objects. If you are thinking of adding essential
functions, it might be better to create your own sequence class.

The methods take as their first argument a sequence object. It is
possible to use methods without first creating a SeqUtils object,
i.e. use it as an anonymous hash.

The first two methodsgive out or read in protein sequences coded in
three letter IUPAC amino acid codes.

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

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki@ebi.ac.uk
Address:

     EMBL Outstation, European Bioinformatics Institute
     Wellcome Trust Genome Campus, Hinxton
     Cambs. CB10 1SD, United Kingdom

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqUtils;
use vars qw(@ISA);
use strict;
use Carp;
use Bio::Tools::CodonTable;

@ISA = qw(Bio::Root::RootI);
# new inherited from RootI

{

my  %onecode =
    ('Ala' => 'A', 'Asx' => 'B', 'Cys' => 'C', 'Asp' => 'D',
     'Glu' => 'E', 'Phe' => 'F', 'Gly' => 'G', 'His' => 'H',
     'Ile' => 'I', 'Lys' => 'K', 'Leu' => 'L', 'Met' => 'M',
     'Asn' => 'N', 'Pro' => 'P', 'Gln' => 'Q', 'Arg' => 'R',
     'Ser' => 'S', 'Thr' => 'T', 'Val' => 'V', 'Trp' => 'W',
     'Xaa' => 'X', 'Tyr' => 'Y', 'Glx' => 'Z', 'Ter' => '*',
     'Sel' => 'U'
     );

my  %threecode =
    ('A' => 'Ala', 'B' => 'Asx', 'C' => 'Cys', 'D' => 'Asp',
     'E' => 'Glu', 'F' => 'Phe', 'G' => 'Gly', 'H' => 'His',
     'I' => 'Ile', 'K' => 'Lys', 'L' => 'Leu', 'M' => 'Met',
     'N' => 'Asn', 'P' => 'Pro', 'Q' => 'Gln', 'R' => 'Arg',
     'S' => 'Ser', 'T' => 'Thr', 'V' => 'Val', 'W' => 'Trp',
     'Y' => 'Tyr', 'Z' => 'Glx', 'X' => 'Xaa', '*' => 'Ter',
     'U' => 'Sel'
     );

=head2 seq3

 Title   : seq3
 Usage   : $string = Bio::SeqUtils->seq3($seq)
 Function:

           Read only method that returns the amino acid sequence as a
           string of three letter codes. alphabet has to be
           'protein'. Output follows the IUPAC standard plus 'Ter' for
           terminator. Any unknown character, including the default
           unknown character 'X', is changed into 'Xaa'. A noncoded
           aminoacid selenocystein is recognized (Sel, U).

 Returns : A scalar
 Args    : character used for stop in the protein seqence optional,
           defaults to '*' string used to separate the output amino
           acid codes, optional, defaults to ''

=cut

sub seq3 {
   my ($self, $seq, $stop, $sep ) = @_;

   $seq->isa('Bio::PrimarySeqI') ||
       $self->throw('Not a Bio::PrimarySeqI object but [$self]');
   $seq->alphabet eq 'protein' ||
       $self->throw('Not a protein sequence');

   if (defined $stop) {
       length $stop != 1 and $self->throw('One character stop needed, not [$stop]');
       $threecode{$stop} = "Ter";
   }
   $sep ||= '';

   my $aa3s;
   foreach my $aa  (split //, uc $seq->seq) {
       $threecode{$aa} and $aa3s .= $threecode{$aa}. $sep, next;
       $aa3s .= 'Xaa'. $sep;
   }
   $sep and substr($aa3s, -(length $sep), length $sep) = '' ;
   return $aa3s;
}

=head2 seq3in

 Title   : seq3in
 Usage   : $string = Bio::SeqUtils->seq3in($seq, 'MetGlyTer')
 Function:

           Read only method that returns the amino acid sequence as a
           string of three letter codes. alphabet has to be
           'protein'. Output follows the IUPAC standard plus 'Ter' for
           terminator. Any unknown character, including the default
           unknown character 'X', is changed into 'Xaa'

 Returns : Bio::PrimarySeq object;
 Args    : character to be used for stop in the protein seqence,
              optional, defaults to '*'
           character to be used for unknown in the protein seqence,
              optional, defaults to 'X'
           string used to separate the output amino acid codes,
              optional, defaults to ''

=cut

sub seq3in {
   my ($self, $seq, $string, $stop, $unknown) = @_;

   $seq->isa('Bio::PrimarySeqI') ||
       $self->throw('Not a Bio::PrimarySeqI object but [$self]');
   $seq->alphabet eq 'protein' ||
       $self->throw('Not a protein sequence');

   if (defined $stop) {
       length $stop != 1 and $self->throw('One character stop needed, not [$stop]');
       $onecode{'Ter'} = $stop;
   }
   if (defined $unknown) {
       length $unknown != 1 and $self->throw('One character stop needed, not [$unknown]');
       $onecode{'Xaa'} = $unknown;
   }

   my ($aas, $aa3);
   my $length = (length $string) - 2;
   for (my $i = 0 ; $i < $length ; $i += 3)  {
       $aa3 = substr($string, $i, 3);
       $onecode{$aa3} and $aas .= $onecode{$aa3}, next;
       warn("Unknown three letter amino acid code [$aa3] ignored");
   }
   $seq->seq($aas);
   return $seq;
}

}

1;




