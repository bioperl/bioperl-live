# $Id$
#
# Copyright (c) 1997-2001 bioperl, Chad Matsalla. All Rights Reserved.
#           This module is free software; you can redistribute it and/or
#           modify it under the same terms as Perl itself.
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::phd - .phd file input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::SeqIO> class.

=head1 DESCRIPTION

This object can transform .phd files (from Phil Green\'s phred basecaller)
to and from Bio::Seq::SeqWithQuality objects

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://www.bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR Chad Matsalla

Chad Matsalla
bioinformatics@dieselwurks.com

=head1 CONTRIBUTORS

Jason Stajich, jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::phd;
use vars qw(@ISA);
use strict;
use Bio::SeqIO;
use Bio::Seq::PrimaryQual;
use Bio::PrimarySeq;
use Bio::Seq::SeqWithQuality;

@ISA = qw(Bio::SeqIO);

=head2 next_phd()

 Title   : next_phd()
 Usage   : $seq = $stream->next_phd()
 Function: returns the next phred sequence in the stream
 Returns : Bio::Seq::SeqWithQuality object
 Args    : NONE
 Notes   : This is really redundant because AFAIK there is no such thing as
	a .phd file that contains more then one sequence. It is included as
	an interface thing and for consistency.

=cut

sub next_phd {
    # print("CSM next_phd\n");
    # print("CSM Invoking next_primary_phd with $_[0] and 1\n");
    # ::dumpValue($_[0]);
    # print("\n");
    my $something = next_primary_phd( $_[0], 1 );
    # print("CSM The result is:\n");
    # ::dumpValue($something);
    # print("\n");
    return $something;    
    # return next_primary_phd( $_[0], 1 );
}

=head2 next_primary_phd()

 Title   : next_primary_phd()
 Usage   : $seq = $stream->next_primary_phd()
 Function: returns the next phred sequence in the stream
 Returns : Bio::Seq::SeqWithQuality object
 Args    : NONE (huh?)

=cut

sub next_primary_phd {
		# print("CSM next_primary_phd!\n");
	my( $self, $as_next_qual ) = @_;
	my $entry;
	# if (!($entry = $self->_readline)) { return; }
	my ($qual,$seq);
        my $in_dna = 0;
        my $base_number = 0;
	my $done;
	my (@lines, @bases, @qualities,$id);
	while ($entry = $self->_readline) {
		return if (!$entry);
		chomp($entry);
		if ($entry =~ /^BEGIN_SEQUENCE\s+(\S+)/) {
			$self->warn("Setting id to $1\n") if($self->verbose > 0);
			$id = $1;
			$entry = $self->_readline();
		}
		if ($entry =~ /^BEGIN_DNA/) {
			$entry =~ /^BEGIN_DNA/;
			$in_dna = 1;
			$entry = $self->_readline();
		}
		if ($entry =~ /^END_DNA/) {
			$in_dna = 0;
		}
		if ($entry =~ /^END_SEQUENCE/) {
		}
		if (!$in_dna) { next;  }
		$entry =~ /(\S+)\s+(\S+)/;
		push @bases,$1;
		push @qualities,$2;
		push(@lines,$entry);
	}
	$self->warn("Creating objects with id = $id\n") if( $self->verbose > 0);
	my $swq = Bio::Seq::SeqWithQuality->new(
					-seq	=>	join('',@bases),
					-qual        => \@qualities,
					-id         => $id,
					-primary_id => $id,
					-display_id => $id,
                              );
	return $swq;
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


=head2 write_phd

 Title   : write_phd(-SeqWithQuality => $swq, <comments>)
 Usage   : $obj->write_swq(     -SeqWithQuality => $swq,);
 Function: Write out an scf.
 Returns : Nothing.
 Args    : Requires: a reference to a SeqWithQuality object to form the
        basis for the scf. Any other arguments are assumed to be comments
        and are put into the comments section of the scf. Read the
        specifications for scf to decide what might be good to put in here.
 Notes   : These are the comments that reside in the header of a phd file
	at the present time. If not provided in the parameter list for
	write_phd(), the following default values will be used:
	CHROMAT_FILE: $swq->id()
	ABI_THUMBPRINT: 0
	PHRED_VERSION: 0.980904.e
	CALL_METHOD: phred
	QUALITY_LEVELS: 99
	TIME: <current time>
	TRACE_ARRAY_MIN_INDEX: 0
	TRACE_ARRAY_MAX_INDEX: unknown
	CHEM: unknown
	DYE: unknown


=cut

sub write_phd {
	my ($self,%args) = @_;
        my %comments;
	my @phredstack;
        my ($label,$arg);
        my $swq = $args{-SeqWithQuality};
	unless (ref($swq) eq "Bio::Seq::SeqWithQuality") {
		$self->throw("You must pass a Bio::Seq::SeqWithQuality object to write_scf as a parameter named \"SeqWithQuality\"");
	}
	my $id = $swq->id();
	if (!$id) { $id = "UNDEFINED in SeqWithQuality Object"; }
	push @phredstack,("BEGIN_SEQUENCE $id","","BEGIN_COMMENT","");
	if ($args{-CHROMAT_FILE}) { push @phredstack,"CHROMAT_FILE: ".$args{-CHROMAT_FILE}; }
	else { push @phredstack,"CHROMAT_FILE: undefined in write_phd"; }
	if ($args{-ABI_THUMBPRINT}) { push @phredstack,"ABI_THUMBPRINT: ".$args{-ABI_THUMBPRINT}; }
	else { push @phredstack,"ABI_THUMBPRINT: 0"; }
	if ($args{-PHRED_VERSION}) { push @phredstack,"PHRED_VERSION: ".$args{-PHRED_VERSION}; }
	else { push @phredstack,"PHRED_VERSION: 0.980904.e"; }
	if ($args{-CALL_METHOD}) { push @phredstack,"CALL_METHOD: ".$args{-CALL_METHOD}; }
	else { push @phredstack,"CALL_METHOD: phred"; }
	if ($args{-QUALITY_LEVELS}) { push @phredstack,"".$args{-QUALITY_LEVELS}; }
	else { push @phredstack,"QUALITY_LEVELS: 99"; }
	my $localtime = localtime();
	if ($args{-TIME}) { push @phredstack,"TIME: ".$args{-TIME}; }
	else { push @phredstack,"TIME: $localtime"; }
	if ($args{-TRACE_ARRAY_MIN_INDEX}) { push @phredstack,"TRACE_ARRAY_MIN_INDEX: ".$args{-TRACE_ARRAY_MIN_INDEX}; }
	else { push @phredstack, "TRACE_ARRAY_MIN_INDEX: 0"; }
	if ($args{-TRACE_ARRAY_MAX_INDEX}) { push @phredstack,"TRACE_ARRAY_MAX_INDEX: ".$args{-TRACE_ARRAY_MAX_INDEX}; }
	else { push @phredstack, "TRACE_ARRAY_MAX_INDEX: unknown";}
	if ($args{-CHEM}) { push @phredstack,"CHEM: ".$args{-CHEM}; }
	else { push @phredstack, "CHEM: unknown"; }
	if ($args{-DYE} ) { push @phredstack,"DYE: ".$args{-DYE}; }
	else { push @phredstack, "DYE: unknown"; }
	push @phredstack,("END_COMMENT","","BEGIN_DNA");
	my $fh = $self->_filehandle();
	foreach (@phredstack) { print $fh ($_."\n"); }
	my $length = $swq->length();
	if ($length eq "DIFFERENT") {
		$self->throw("Can't create the phd because the sequence and the quality in the SeqWithQuality object are of different lengths.");
	}
	for (my $curr = 1; $curr<=$length; $curr++) {
		print $fh (uc($swq->baseat($curr))." ".$swq->qualat($curr)."\n");
	}
	print $fh ("END_DNA\n\nEND_SEQUENCE\n");
}

1;
__END__
