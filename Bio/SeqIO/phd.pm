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

This object can transform .phd files (from Phil Green's phred basecaller)
to and from Bio::Seq::Quality objects

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.
Bug reports can be submitted via the web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR Chad Matsalla

Chad Matsalla
bioinformatics@dieselwurks.com

=head1 CONTRIBUTORS

Jason Stajich, jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# 'Let the code begin...

package Bio::SeqIO::phd;
use strict;
use Bio::Seq::SeqFactory;
use Dumpvalue();
my $dumper = new Dumpvalue();

use base qw(Bio::SeqIO);

sub _initialize {
  my($self,@args) = @_;
  $self->SUPER::_initialize(@args);
  if( ! defined $self->sequence_factory ) {
      $self->sequence_factory(Bio::Seq::SeqFactory->new
			      (-verbose => $self->verbose(),
			       -type => 'Bio::Seq::Quality'));
  }
}

=head2 next_seq()

 Title   : next_seq()
 Usage   : $swq = $stream->next_seq()
 Function: returns the next phred sequence in the stream
 Returns : Bio::Seq::Quality object
 Args    : NONE
 Notes   : This is really redundant because AFAIK there is no such thing as
  	   a .phd file that contains more then one sequence. It is included as
	   an interface thing and for consistency.

=cut

sub next_seq {
    my ($self,@args) = @_;
    my ($entry,$done,$qual,$seq);
    my ($id,@lines, @bases, @qualities, @trace_indices) = ('');
    if (!($entry = $self->_readline)) { return; }
	if ($entry =~ /^BEGIN_SEQUENCE\s+(\S+)/) {
          $id = $1;
     }
     my $in_comments = 0;
    my $in_dna = 0;
    my $base_number = 0;
     my $comments = {};
    while ($entry = $self->_readline) {
	return if (!$entry);
	chomp($entry);
     if ($entry =~ /^BEGIN_COMMENT/) {
          $in_comments = 1;
          while ($in_comments == 1) {
              $entry = $self->_readline();
               chomp($entry);
              if ($entry) {
                    if ($entry =~ /^END_COMMENT/) {
                         $in_comments = 0;
                    }
                    else {
                         my ($name,$content) = split(/:/,$entry);
                         if ($content) { $content =~ s/^\s//g; }
                         $comments->{$name} = $content;
                    }
               }
          }
     }
	if ($entry =~ /^BEGIN_CHROMAT:\s+(\S+)/) {
	     # this is where I used to grab the ID
          if (!$id) {
               $id = $1;
          }
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
	$entry =~ /(\S+)\s+(\S+)(?:\s+(\S+))?/;
	push @bases,$1;
	push @qualities,$2;
	# Trace index values are required for phd file
	push(@trace_indices,$3) if defined $3;
	push(@lines,$entry);
    }
     # $self->debug("csmCreating objects with id = $id\n");
    my $swq = $self->sequence_factory->create
	(-seq        => join('',@bases),
	 -qual       => \@qualities,
	 -trace      => \@trace_indices,
	 -id         => $id,
	 -primary_id => $id,
	 -display_id => $id,
	 );
     # this should be an actual object to assist in serialization
     # but I don't have time for this now.
     if ($comments) { $swq->{comments} = $comments; }
    return $swq;
}

=head2 write_seq

 Title   : write_seq(-Quality => $swq, <comments>)
 Usage   : $obj->write_seq(     -Quality => $swq,);
 Function: Write out an scf.
 Returns : Nothing.
 Args    : Requires: a reference to a Quality object to form the
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
     IMPORTANT: This method does not write the trace index where this
          call was made. All base calls are placed at index 1.


=cut

sub write_seq {
    my ($self,@args) = @_;
    my @phredstack;
    my ($label,$arg);

    my ($swq, $swq2, $chromatfile, $abithumb,
	$phredversion, $callmethod,
	$qualitylevels,$time,
	$trace_min_index,
	$trace_max_index,
	$chem, $dye
	) = $self->_rearrange([qw(QUALITY
                                  SEQWITHQUALITY
				  CHROMAT_FILE
				  ABI_THUMBPRINT
				  PHRED_VERSION
				  CALL_METHOD
				  QUALITY_LEVELS
				  TIME
				  TRACE_ARRAY_MIN_INDEX
				  TRACE_ARRAY_MAX_INDEX
				  CHEM
				  DYE
				  )], @args);

    $swq = $swq2 if not $swq and $swq2;
    unless (ref($swq) eq "Bio::Seq::Quality") {
	$self->throw("You must pass a Bio::Seq::Quality object to write_scf as a parameter named \"Quality\"");
    }
    my $id = $swq->id();
    if (!$id) { $id = "UNDEFINED in Quality Object"; }
    push @phredstack,("BEGIN_SEQUENCE $id","","BEGIN_COMMENT","");

    $chromatfile = 'undefined in write_phd' unless defined $chromatfile;
    push @phredstack,"CHROMAT_FILE: $chromatfile";

    $abithumb = 0 unless defined $abithumb;
    push @phredstack,"ABI_THUMBPRINT: $abithumb";

    $phredversion = "0.980904.e" unless defined $phredversion;
    push @phredstack,"PHRED_VERSION: $phredversion";

    $callmethod = 'phred' unless defined $callmethod;
    push @phredstack,"CALL_METHOD: $callmethod";

    $qualitylevels = 99 unless defined $qualitylevels;
    push @phredstack,"QUALITY_LEVELS: $qualitylevels";

    $time = localtime() unless defined $time;
    push @phredstack,"TIME: $time";

    $trace_min_index = 0 unless defined $trace_min_index;
    push @phredstack,"TRACE_ARRAY_MIN_INDEX: $trace_min_index";

    $trace_max_index = '10000'  unless defined $trace_max_index;
    push @phredstack,"TRACE_ARRAY_MAX_INDEX: $trace_max_index";

    $chem = 'unknown' unless defined $chem;
    push @phredstack,"CHEM: $chem";

    $dye = 'unknown' unless defined $dye;
    push @phredstack, "DYE: $dye";

    push @phredstack,("END_COMMENT","","BEGIN_DNA");

    foreach (@phredstack) {  $self->_print($_."\n"); }

    my $length = $swq->length();
    if ($length eq "DIFFERENT") {
	$self->throw("Can't create the phd because the sequence and the quality in the Quality object are of different lengths.");
    }
    for (my $curr = 1; $curr<=$length; $curr++) {
	$self->_print (sprintf("%s %s %s\n",
			       uc($swq->baseat($curr)),
			       $swq->qualat($curr),
			       $swq->trace_index_at($curr)));
    }
    $self->_print ("END_DNA\n\nEND_SEQUENCE\n");

    $self->flush if $self->_flush_on_write && defined $self->_fh;
    return 1;
}

1;
__END__
