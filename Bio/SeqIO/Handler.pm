
#
# BioPerl module for Bio::SeqIO::Handler
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::Handler - Handler of the SeqIO classes for tieing to filehandles

=head1 SYNOPSIS

    $stream = Bio::SeqIO->new(-file => "filename", -format => 'Fasta');

    tie *SEQ, 'Bio::SeqIO::Handler', $stream;

    while $seq ( <SEQ> ) {
	# $seq is Bio::Seq object
    }

=head1 DESCRIPTION

This object wraps over the Bio::SeqIO:: format classes, providing
the correct glue for the tie mechanism in to filehandles. This 
allows normal looking perl filehandles can be used to provide a 
streams for converting sequences to and from ascii formats

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

Email birney@sanger.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqIO::Handler;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object Exporter);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,$obj) = @_;

  my $make = $self->SUPER::_initialize;
  $self->_streamobj($obj);

# set stuff in self from @args
 return $make; # success - we hope!
}

sub TIEHANDLE {
    my $class = shift;
    my $obj  = shift;
    my $self;

    $self = $class->new($obj);
    
    return $self;
}

sub READLINE {
    my $self = shift;
    my $obj = $self->_streamobj;

    return $obj->next_seq();
}

sub PRINT {
    my ($self,@args) = @_;
    my $obj = $self->_streamobj;
    my ($seq,$count);
    
    foreach $seq ( @args ) {
	if( ! $seq->isa("Bio::Seq") ) {
	    $self->warn("$seq does not seem to be a sequence object. No passing it onto the Bio::SeqIO stream");
	} else {
	    $obj->write_seq($seq);
	    $count++;
	}
    }

    return $count;
}

sub WRITE {
    my ($self,@args) = @_;
    
    $self->warn("Using a write mechanism on Bio::SeqIO::Handler. Should use print HANDLE \$seq, not write. Forwarding onto print now");
    $self->PRINT(@args);
}

sub PRINTF {
    my $self = shift;

    $self->throw("Attempting to use a Bio::SeqIO::Handler in an inappropiate context !");
}

sub READ {
    my $self = shift;

    $self->throw("Attempting to use a Bio::SeqIO::Handler in an inappropiate context !");
}

sub GETC {
    my $self = shift;

    $self->throw("Attempting to use a Bio::SeqIO::Handler in an inappropiate context !");
}

sub CLOSE {
    my $self = shift;

    $self->warn("Attempting to close a Bio::SeqIO::Handler. No need - just reassign it to a different file or use untie");
}


=head2 _streamobj

 Title   : _streamobj
 Usage   : $obj->_streamobj($newval)
 Function: 
 Example : 
 Returns : value of _streamobj
 Args    : newvalue (optional)


=cut

sub _streamobj{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_streamobj'} = $value;
    }
    return $obj->{'_streamobj'};

}

sub DESTROY {
    my $self=shift;

    $self->{'_streamobj'} = '';
}





