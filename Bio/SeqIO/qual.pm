# $Id$
#
# Copyright (c) 1997-9 bioperl, Chad Matsalla. All Rights Reserved.
#           This module is free software; you can redistribute it and/or
#           modify it under the same terms as Perl itself.
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::Qual - .qual file input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::SeqIO> class.

=head1 DESCRIPTION

This object can transform .qual (similar to fasta) objects to and from
?objects? <fix this>

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

=head1 AUTHOR Chad Matsalla

Chad Matsalla
bioinformatics@dieselwurks.com

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::qual;
use vars qw(@ISA);
use strict;
use Bio::SeqIO;
use Bio::Seq::PrimaryQual;
use Bio::Seq::SeqWithQuality;
require 'dumpvar.pl';

@ISA = qw(Bio::SeqIO);

=head2 next_qual

 Title   : next_qual
 Usage   : $seq = $stream->next_qual()
 Function: returns the next sequence in the stream
 Returns : Bio::Seq::PrimaryQual object
 Args    : NONE

=cut

sub next_qual {
		# print("CSM next_qual\n");
		# print("CSM Invoking next_primary_qual with $_[0] and 1\n");
		# ::dumpValue($_[0]);
		# print("\n");
	my $something = next_primary_qual( $_[0], 1 );
		# print("CSM The result is:\n");
		# ::dumpValue($something);
		# print("\n");
	return $something;    
		# return next_primary_qual( $_[0], 1 );
}

=head2 next_primary_qual()

 Title   : next_primary_qual()
 Usage   : $seq = $stream->next_primary_qual()
 Function: returns the next sequence in the stream
 Returns : Bio::PrimaryQual object
 Args    : NONE

=cut

sub next_primary_qual {
	# print("CSM next_primary_qual!\n");
  my( $self, $as_next_qual ) = @_;
  my ($qual,$seq);
  my $moltype;
  local $/ = "\n>";

  return unless my $entry = $self->_readline;

  if ($entry eq '>')  {  # very first one
    return unless $entry = $self->_readline;
  }

  	# original: my ($top,$sequence) = $entry =~ /^(.+?)\n([^>]*)/s
  my ($top,$sequence) = $entry =~ /^(.+?)\n([^>]*)/s
    or $self->throw("Can't parse entry [$entry]");
  my ($id,$fulldesc) = $top =~ /^\s*(\S+)\s*(.*)/
    or $self->throw("Can't parse fasta header");
  $id =~ s/^>//;
  	# $sequence =~ s/\s//g; # Remove whitespace
  # for empty sequences we need to know the mol.type
	# no we don't, not for PrimaryQuals because... well just because.
	# $moltype = $self->moltype();
	# print("CSM \$moltype is $moltype\n");
  if(length($sequence) == 0) {
      if(! defined($moltype)) {
          # let's default to dna
		# lets not.
		# $moltype = "dna";
      }
  } else {
      # we don't need it really, so disable
	# you bet we don't need it because PrimaryQual doesn't pay it any mind anyway
	# $moltype = undef;
  }

  # create the seq object
	$sequence =~ s/\n//g;
  if ($as_next_qual) {
	# print("CSM qual.pm: creating a primaryqual object with $sequence\n");
    $qual = Bio::Seq::PrimaryQual->new(-qual        => $sequence,
                                -id         => $id,
                                -primary_id => $id,
				-display_id => $id,
                                -desc       => $fulldesc
                                );
  }
  # if there wasn't one before, set the guessed type
	# no, don't.
  	# $self->moltype($qual->moltype());
  	# print("CSM next_primary_qual: returning $qual.\n");
  return $qual;
}

=head2 _initialize()

 Title   : _initialize()
 Usage   :
 Function: Bioperl initialize.
 Returns :
 Args    :
 Notes   :

=cut

sub _initialize {
  my($self,@args) = @_;
  return unless my $make = $self->SUPER::_initialize(@args);
}

=head2 write_qual(-source => $source, -header => "some information")

 Title   : write_qual(-source => $source, -header => "some information")
 Usage   : $obj->write_qual(	-source => $source,
				-header => "some information");
 Function: Write out an list of quality values to a fasta-style file.
 Returns : Nothing.
 Args    : Requires: a reference to a SeqWithQuality object or a
	PrimaryQual object as the -source. Optional: information
	for the header.
 Notes   : If no -header is provided, $obj->id() will be used where
	$obj is a reference to either a SeqWithQuality object or a
	PrimaryQual object. If $source->id() fails, ">unknown" will be
	the header. If the SeqWithQuality object has $source->length() of
	"DIFFERENT" (read the pod, luke), write_qual will use the length
	of the PrimaryQual object within the SeqWithQuality object.

=cut

sub write_qual {
	my ($self,%args) = @_;
	my $source = $args{-source};
	if (!$source || ( ref($source) ne "Bio::Seq::SeqWithQuality" && ref($source) ne "Bio::Seq::PrimaryQual")) {
		$self->throw("You must pass a Bio::Seq::SeqWithQuality or a Bio::Seq::PrimaryQual object to write_qual as a parameter named \"source\"");
	}
	my $header = $source->id();
	if (!$header) { $header = "unknown"; }
	my @quals = $source->qual();
	# ::dumpValue(\@quals);
	my $fh = $self->_filehandle();
	print $fh (">$header \n");
	my (@slice,$max,$length);
	$length = $source->length();
	if ($length eq "DIFFERENT") {
		$self->warn("You passed a SeqWithQuality object that contains a sequence and quality of differing lengths. Using the length of the PrimaryQual component of the SeqWithQuality object.");
		$length = $source->qual_obj()->length();
	}
		# print("Printing $header to a file.\n");
	for (my $count = 1; $count<$length; $count+= 50) {
		if ($count+50 > $length) { $max = $length; }
		else { $max = $count+49; }
		my @slice = @{$source->subqual($count,$max)};
		print $fh join(' ',@slice);
		print $fh " \n";
	}
}


1;
__END__
