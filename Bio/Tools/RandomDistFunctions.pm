#
# BioPerl module for Bio::Tools::RandomDistFunctions
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Jason Stajich <jason-at-bioperl.org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::RandomDistFunctions - A set of routines useful for
generating random data in different distributions

=head1 SYNOPSIS

  use Bio::Tools::RandomDistFunctions;
  my $dist = Bio::Tools::RandomDistFunctions->new();
  for my $v ( 1..1000 ) { 
    my $birth_dist = $dist->rand_birth_distribution($lambda);
    # ... do something with the variable
  }

=head1 DESCRIPTION

Most of the code is based on the C implementation of these routines in
Mike Sanderson's r8s's package.  See http://loco.biosci.arizona.edu/r8s/ for
information on his software.

=for comment
This code tries to be fast and use available faster BigInt and GMP
library methods when those modules are available.

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

Email jason-at-bioperl.org

=head1 CONTRIBUTORS

Thanks to Mike Sanderson for assistance in the getting this
implementation together.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Tools::RandomDistFunctions;
require Exporter;
use vars qw(%LOADED @EXPORT_OK); use strict;

#use Math::BigFloat lib => 'GMP,Bit::Vector';
#use Math::BigInt   lib => 'GMP,Bit::Vector';
use POSIX;

use base qw(Bio::Root::Root);

=head2 birth_distribution

 Title   : rand_birth_distribution
 Usage   : my $randvar = $dist->
            rand_birth_distribution($lambda);
 Function: Returns a random number from a birth process waiting 
           time with a fixed interval
           1.0.  Times are measured from 0=present,1=root;
 Returns : floating point number
 Args    : $lambda ( > 0 )
 References : This is based on code by Mike Sanders in r8s.
              Ross, Stochastic Processes, p. 145 for the density

=cut

sub rand_birth_distribution{
   my ($self,$lambda) = @_;
   if( ! ref($self) && 
       $self !~ /RandomDistFunctions/ ) { 
       $lambda = $self;
   }
   unless( $lambda ) { 
       $self->throw("Cannot call birth_distribution without a valid lambda value (>0)");
   }
   return 1 - (log(rand(1) * (exp($lambda) - 1)+1)/ $lambda); 
}


=head2 rand_geometric_distribution

 Title   : rand_geometric_distribution
 Usage   : my $randvar = $dist->rand_geometric_distribution($param);
 Function: Returns a random geometric variate distributed with 
           parameter $param, according to
           c.d.f. 1 - ( 1- param) ^ n 
 Returns : integer
 Args    : $param ( 0 > $param < 1 )


=cut

sub rand_geometric_distribution{
   my ($self,$param) = @_;
   if( ! ref($self) && 
       $self !~ /RandomDistFunctions/ ) { 
       $param = $self;
   }
   unless( $param ) { 
       $self->throw("Cannot call rand_geometric_distribution without a valid param value (>0)");
   }

   my $den;
   if( $param < 1e-8) { 
       $den = (-1 * $param) - ( $param * $param ) / 2;
   } else { 
       $den = log(1 - $param);
   }
   my $z = log(1 - rand(1)) / $den;
   return POSIX::floor($z) + 1;
   # MSanderson comments from r8s code
   # Is this the right truncation of the real-valued expression above?
   # YES
   # Checked by reference to the expected mean of the distribution in
   # 100,000 replicates
   # EX = 1/param Var = (1-param)/param^2 See Olkin, Gleser, and
   # Derman, p. 193ff. Probability Models and Applications, 1980.
}

=head2 rand_exponentional_distribution

 Title   : rand_exponentional_distribution
 Usage   : my $var = $dist->rand_exponentional_distribution($param);
 Function: Returns a random exponential variate distributed with parameter
           $param, according to c.d.f 1 - e^(-param * x)
 Returns : floating point number 
 Args    : $param ( > 0 )


=cut

sub rand_exponentional_distribution {
   my ($self,$param) = @_;
      if( ! ref($self) && 
       $self !~ /RandomDistFunctions/ ) { 
       $param = $self;
   }
   unless( $param ) { 
       $self->throw("Cannot call rand_exponentional_distribution without a valid param value (>0)");
   }
   return log( 1- rand(1)) / $param;
}

=head2 rand_normal_distribution

 Title   : rand_normal_distribution
 Usage   : my $var = $dist->rand_normal_distribution()
 Function: Returns a random normal (gaussian) variate distributed 
 Returns : floating point number
 Args    : none


=cut

sub rand_normal_distribution{
    my $gset;
    my ($rsq,$v1,$v2) = ( 0,0,0);
    do { 
	$v1 = 2 * rand(1) - 1;
	$v2 = 2 * rand(1) - 1;
	$rsq= $v1**2 + $v2 ** 2;
    } while( $rsq >= 1 || $rsq == 0);
    my $fac = sqrt(-2.0 * log($rsq) / $rsq );
    $gset = $v1 * $fac;
    return $v2 * $fac;
}

1;
