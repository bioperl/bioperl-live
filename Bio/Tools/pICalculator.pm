#
# BioPerl module for Bio::Tools::pICalculator
#
# Copyright (c) 2002, Merck & Co. Inc. All Rights Reserved.
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::pICalculator - calculate the isoelectric point of a protein

=head1 DESCRIPTION

Calculates the isoelectric point of a protein, the pH at which there
is no overall charge on the protein. Calculates the charge on a protein
at a given pH. Can use built-in sets of pK values or custom pK sets.

=head1 SYNOPSIS

  use Bio::Tools::pICalculator;
  use Bio::SeqIO;

  my $in = Bio::SeqIO->new( -fh => \*STDIN ,
                            -format => 'Fasta' );

  my $calc = Bio::Tools::pICalculator->new(-places => 2,
                                           -pKset => 'EMBOSS');

  while ( my $seq = $in->next_seq ) {
     $calc->seq($seq);
     my $iep = $calc->iep;
     print sprintf( "%s\t%s\t%.2f\n",
                    $seq->id,
                    $iep,
                    $calc->charge_at_pH($iep) );

     for( my $i = 0; $i <= 14; $i += 0.5 ){
        print sprintf( "pH = %.2f\tCharge = %.2f\n",
                       $i,
                       $calc->charge_at_pH($i) );
     }
  }

=head1 SEE ALSO

http://fields.scripps.edu/DTASelect/20010710-pI-Algorithm.pdf
http://emboss.sourceforge.net/apps/cvs/emboss/apps/iep.html
http://us.expasy.org/tools/pi_tool.html

=head1 LIMITATIONS

There are various sources for the pK values of the amino acids. 
The set of pK values chosen will affect the pI reported.

The charge state of each residue is assumed to be independent of 
the others. Protein modifications (such as a phosphate group) that 
have a charge are ignored.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions 
preferably to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via the 
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

Mark Southern (mark_southern@merck.com). From an algorithm by David 
Tabb found at http://fields.scripps.edu/DTASelect/20010710-pI-Algorithm.pdf.
Modification for Bioperl, additional documentation by Brian Osborne.

=head1 COPYRIGHT

Copyright (c) 2002, Merck & Co. Inc. All Rights Reserved. This module is
free software. It may be used, redistributed and/or modified under the terms
of the Perl Artistic License (see http://www.perl.com/perl/misc/Artistic.html)

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Private methods are usually preceded by a _.

=cut

# Let the code begin...

package Bio::Tools::pICalculator;
use strict;


use base qw(Bio::Root::Root);

# pK values from the DTASelect program from Scripps
# http://fields.scripps.edu/DTASelect
my $DTASelect_pK = {  N_term   =>  8.0,
                      K        => 10.0, # Lys
                      R        => 12.0, # Arg
                      H        =>  6.5, # His
                      D        =>  4.4, # Asp
                      E        =>  4.4, # Glu
                      C        =>  8.5, # Cys
                      Y        => 10.0, # Tyr
                      C_term   =>  3.1
                    };

# pK values from the iep program from EMBOSS
# http://emboss.sourceforge.net/apps/cvs/emboss/apps/iep.html
my $Emboss_pK  = { N_term   =>  8.6,
                   K        => 10.8, # Lys
                   R        => 12.5, # Arg
                   H        =>  6.5, # His
                   D        =>  3.9, # Asp
                   E        =>  4.1, # Glu
                   C        =>  8.5, # Cys
                   Y        => 10.1, # Tyr
                   C_term   =>  3.6
                 };

=head2 desc

 Title   : new
 Usage   : Bio::Tools::pICalculator->new
 Function: Instantiates the Bio::Tools::pICalculator object
 Example : $calc = Bio::Tools::pICalculator->new( -pKset => \%pKvalues,
                                                  # a Bio::Seq object
                                                  -seq => $seq,
                                                  -places => 2 );
           or:

           $calc = Bio::Tools::pICalculator->new( -pKset => 'string',
                                                  # a Bio::Seq object
                                                  -seq => $seq,
                                                  -places => 1 );

           Constructs a new pICalculator. Arguments are a flattened hash.
           Valid, optional keys are:

           pKset - A reference to a hash with key value pairs for the 
                   pK values of the charged amino acids. Required keys
                   are:

                   N_term   C_term   K   R   H   D   E   C   Y

           pKset - A string ( 'DTASelect' or 'EMBOSS' ) that will 
                   specify an internal set of pK values to be used. The 
                   default is 'EMBOSS'

           seq - A Bio::Seq sequence object to analyze

           places - The number of decimal places to use in the
                    isoelectric point calculation. The default is 2.

 Returns : The description
 Args    : The description or none

=cut

sub new {
   my( $class, %opts ) = @_;
   my $self = $class->SUPER::new(%opts);
   $self = bless {}, ref $self || $self;
   $self->seq( $opts{-seq} ) if exists $opts{-seq};
   $self->pKset( $opts{-pKset} || 'EMBOSS' );
   exists $opts{-places} ? $self->places( $opts{-places} ) :
     $self->places(2);
   return $self;
}

=head2 seq

 Title   : seq
 Usage   : $calc->seq($seqobj)
 Function: Sets or returns the Bio::Seq used in the calculation
 Example : $seqobj = Bio::Seq->new(-seq=>"gghhhmmm",-id=>"GHM");
           $calc = Bio::Tools::pICalculator->new;
           $calc->seq($seqobj);
 Returns : Bio::Seq object
 Args    : Bio::Seq object or none

=cut

sub seq {
   my( $this, $seq ) = @_;
   unless( defined $seq && UNIVERSAL::isa($seq,'Bio::Seq') ){
      $this->throw("$seq is not a valid Bio::Seq object");
   }
   $this->{-seq} = $seq;
   $this->{-count} = count_charged_residues( $seq );
   return $this->{-seq};
}

=head2 pKset

 Title   : pKset
 Usage   : $pkSet = $calc->pKSet(\%pKSet)
 Function: Sets or returns the hash of pK values used in the calculation
 Example : $calc->pKset('emboss')
 Returns : reference to pKset hash
 Args    : The reference to a pKset hash, a string, or none. Examples:

           pKset - A reference to a hash with key value pairs for the
                   pK values of the charged amino acids. Required keys
                   are:

                   N_term   C_term   K   R   H   D   E   C   Y

           pKset - A valid string ( 'DTASelect' or 'EMBOSS' ) that will 
                   specify an internal set of pK values to be used. The 
                   default is 'EMBOSS'

=cut

sub pKset {
   my ( $this, $pKset ) = @_;
   if( ref $pKset eq 'HASH' ){         # user defined pK values
      $this->{-pKset} = $pKset;
   }elsif( $pKset =~ /^emboss$/i ){    # from EMBOSS's iep program
      $this->{-pKset} = $Emboss_pK;
   }elsif( $pKset =~ /^dtaselect$/i ){ # from DTASelect (scripps)
      $this->{-pKset} = $DTASelect_pK;
   }else{                              # default to EMBOSS
      $this->{-pKset} = $Emboss_pK;
   }
   return $this->{-pKset};
}

sub places {
   my $this = shift;
   $this->{-places} = shift if @_;
   return $this->{-places};
}

=head2 iep

 Title   : iep
 Usage   : $calc->iep
 Function: Returns the isoelectric point
 Example : $calc = Bio::Tools::pICalculator->new(-places => 2);
           $calc->seq($seqobj);
           $iep = $calc->iep;
 Returns : The isoelectric point of the sequence in the Bio::Seq object
 Args    : None

=cut

sub iep {
   my $this = shift;
   return _calculate_iep($this->{-pKset},
	    		 $this->{-places},
			 $this->{-seq},
			 $this->{-count}
		        );
}

=head2 charge_at_pH

 Title   : charge_at_pH
 Usage   : $charge = $calc->charge_at_pH($pH)
 Function: Sets or gets the description of the sequence
 Example : $calc = Bio::Tools::pICalculator->new(-places => 2);
           $calc->seq($seqobj);
           $charge = $calc->charge_at_ph("7");
 Returns : The predicted charge at the given pH
 Args    : pH

=cut

sub charge_at_pH {
   my $this = shift;
   return _calculate_charge_at_pH( shift, $this->{-pKset},
				  $this->{-count} );
}

sub count_charged_residues {
   my $seq = shift;
   my $sequence = $seq->seq;
   my $count;
   for ( qw( K R H D E C Y ) ){ # charged AA's
      $count->{$_}++ while $sequence =~ /$_/ig;
   }
   return $count;
}

sub _calculate_iep {
    my( $pK, $places, $seq, $count ) = @_;
    my $pH = 7.0;
    my $step = 3.5;
    my $last_charge = 0.0;
    my $format = "%.${places}f";

    unless( defined $count ){
       $count = count_charged_residues($seq);
    }
    while(1){
       my $charge = _calculate_charge_at_pH( $pH, $pK, $count );
       last if sprintf($format,$charge) == 
	 sprintf($format,$last_charge);
       $charge > 0 ? ( $pH += $step ) : ( $pH -= $step );
       $step /= 2.0;
       $last_charge = $charge;
    }
    return sprintf( $format, $pH );
}

# it's the sum of all the partial charges for the
# termini and all of the charged aa's!
sub _calculate_charge_at_pH {
   no warnings; # don't complain if a given key doesn't exist
   my( $pH, $pK, $count ) = @_;
   my $charge =    _partial_charge( $pK->{N_term}, $pH )
   + $count->{K} * _partial_charge( $pK->{K},      $pH )
   + $count->{R} * _partial_charge( $pK->{R},      $pH )
   + $count->{H} * _partial_charge( $pK->{H},      $pH )
   - $count->{D} * _partial_charge( $pH,      $pK->{D} )
   - $count->{E} * _partial_charge( $pH,      $pK->{E} )
   - $count->{C} * _partial_charge( $pH,      $pK->{C} )
   - $count->{Y} * _partial_charge( $pH,      $pK->{Y} )
   -               _partial_charge( $pH, $pK->{C_term} );
   return $charge;
}

# Concentration Ratio is 10**(pK - pH) for positive groups
# and 10**(pH - pK) for negative groups
sub _partial_charge {
   my $cr = 10 ** ( $_[0] - $_[1] );
   return $cr / ( $cr + 1 );
}

1;

__END__
