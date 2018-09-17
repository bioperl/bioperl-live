#---------------------------------------------------------

=head1 NAME

Bio::Matrix::PSM::IO::psiblast - PSM psiblast parser

=head1 SYNOPSIS

See Bio::Matrix::PSM::IO for documentation

=head1 DESCRIPTION

Parser for ASCII matrices from PSI-BLAST (blastpgp program in
BLAST distribution).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - James Thompson

Email tex@biosysadmin.com

=head1 APPENDIX

=cut


# Let the code begin...
package Bio::Matrix::PSM::IO::psiblast;
use Bio::Matrix::PSM::Psm;
use Bio::Matrix::PSM::ProtMatrix;
use strict;

use base qw(Bio::Matrix::PSM::PsmHeader Bio::Matrix::PSM::IO);

# define the order in which amino acids are listed in the psiblast matrix file
our @ordered_alphabet = qw/A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V/;

=head2 new

 Title   : new
 Usage   : my $psmIO  =  Bio::Matrix::PSM::IO->new(-format=>'psiblast', 
                                                  -file=>$file);
 Function: Associates a file with the appropriate parser
 Throws  :
 Example :
 Args    :
 Returns : Bio::Matrix::PSM::ProtMatrix->new(@args);

=cut

sub new {
   my ($class,@args)=@_;
   my $line;

   my $self   = $class->SUPER::new(@args);
   my ($file) = $self->_rearrange(['FILE'], @args);
   $self->_initialize_io(@args) || warn "Did you intend to use STDIN?"; # Read only for now
   $self->_initialize;

   $self->{_ordered_alphabet} = \@ordered_alphabet;
   return $self;
}

=head2 next_psm

 Title   : next_psm
 Usage   : my $psm = $psmIO->next_psm();
 Function: Reads the next PSM from the input file, associated with this object
 Throws  : None
 Returns : Bio::Matrix::PSM::ProtPsm object
 Args    : none

=cut

sub next_psm {
   my $self = shift;
   my $line;

   return if ($self->{_end});

   my %args;
   my @ordered_alphabet = @{$self->{_ordered_alphabet}};

   while ( defined( $line = $self->_readline) ) {
      # remove leading and trailing whitespace
      chomp $line;
      $line =~ s/^\s+//g;
      $line =~ s/\s+$//g;
   
      if ( $line =~ /^(\d+)\s+(\w{1})/ ) {  # match reference aa and position number
         my @elements = split /\s+/, $line;
         
         my $position = shift @elements;
         my $letter   = shift @elements;
   
         my $ratio = pop @elements;
         my $ic    = pop @elements;
   
         # put the next 20 elements into the correct array in %args
         for ( 0 .. 19 ) { push @{$args{'l'.$ordered_alphabet[$_]}}, shift @elements; }
         for ( 0 .. 19 ) { push @{$args{'p'.$ordered_alphabet[$_]}}, shift @elements; }
         
         push @{$args{'ic'}}, $ic;
      }
   }

   $self->{_end} = 1; # psiblast matrix files currently only hold one PSM per file

   my $psm = Bio::Matrix::PSM::ProtMatrix->new( %args );
   return $psm;
}

sub DESTROY {
   my $self=shift;
   $self->close;
}

1;
