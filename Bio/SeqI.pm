
#
# BioPerl module for Bio::SeqI
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqI - Interface definition for a Bio::SeqI

=head1 SYNOPSIS

    # get a Bio::SeqI compliant object somehow

    # to test this is a seq object

    $obj->isa("Bio::SeqI") || $obj->throw("$obj does not implement the Bio::SeqI interface");

    # accessors

    $string = $obj->seq();
    $substring = $obj->seq(12,50);
    $name   = $obj->name();
    $id     = $obj->identifier(); # unique idenitifer assign to the sequence by the system
    
    # object manipulation
   
    eval {
	$rev    = $obj->revcom();
    };
    if( $@ ) {
	$obj->throw("Could not reverse complement. Probably not DNA. Actual exception\n$@\n");
    }

    $trunc = $obj->trunc(12,50);
    
    # $rev and $trunc are Bio::SeqI compliant objects



=head1 DESCRIPTION

 This object defines an abstract interface to sequences. There is a
pure perl implementation of this in Bio::Seq. If you just want to use
Bio::Seq objects, then please read that module first.

This interface defines what bioperl consideres necessary to "be" a sequence,
without 

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@sanger.ac.ul

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqI;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


use AutoLoader;
@ISA = qw(Bio::Root::Object Exporter);
@EXPORT_OK = qw();
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

# set stuff in self from @args
 return $make; # success - we hope!
}
