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

Bio::SeqIO::scf - .scf file input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::SeqIO> class.

=head1 DESCRIPTION

This object can transform .scf files to and from Bio::Seq::SeqWithQuality objects.
Mechanisms are present to retrieve trace data from scf files.

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

package Bio::SeqIO::csmscf;
use vars qw(@ISA);
use strict;
use Bio::SeqIO;
use Bio::Seq::PrimaryQual;
use Bio::PrimarySeq;
use Bio::Seq::SeqWithQuality;

require 'dumpvar.pl';

@ISA = qw(Bio::SeqIO);

=head2 next_scf()
                                                                                                                                  
Title   : next_scf()                                                            
Usage   : $scf = $stream->next_scf()                                            
Function: returns the next scf sequence in the stream                           
Returns : Bio::Seq::SeqWithQuality object                                       
Args    : NONE (huh?)
Notes   : The SCF specification does not provide for having more then
        one sequence in a given scf. So once the filehandle has been open
        and passed to SeqIO don't expect to run this function more then
        once on a given scf unless you embraced and extended the SCF
	standard. (But that's just C R A Z Y talk, isn't it.)

=cut

sub next_scf {
	my( $self, @args ) = @_;
	my ($seq, $seqc, $fh, $buffer, $offset, $length, $read_bytes, @read, %names);
	$fh = $self->_filehandle();
	unless ($fh) {  # simulate the <> function
		if ( !fileno(ARGV) or eof(ARGV) ) {
			return unless my $ARGV = shift;
			open(ARGV,$ARGV) or
				$self->throw("Could not open $ARGV for SCF stream reading $!");
		}
		$fh = \*ARGV;
	}
	binmode $fh; # for the Win32/Mac crowds
	return unless read $fh, $buffer, 128;  # no exception; probably end of file
	$self->_parse_header($buffer);
	print("After getting the header information, the read position is ".tell($fh)."\n");
		# gather the trace information
	$length = $self->{samples}*$self->{sample_size}*4;
	print("Going to read trace sample information. Current file position is: ".tell($fh).". Read length is: $length\n");
	read $fh,$buffer,$length;
	unless (length($buffer) == $length) {
			$self->throw('Unexpected end of file while reading from SCF file');
	}
	@read = unpack "n$length",$buffer;
	$self->_split_traces(\@read);
		# do more sample stuff here

		# now go and get the base information
	$offset = $self->{bases_offset};
	$length = ($self->{bases} * 12);
	print("Reading the bases: Going to $offset to read $length bytes which is ".$self->{bases}." times 12.\n");
	seek $fh,$offset,0;
	print("Just before the read, the curent filehandle position is: ".tell($fh)."\n");
	read $fh, $buffer, $length;
	unless (length($buffer) == $length) {
			$self->throw('Unexpected end of file while reading from SCF file');
	}
	$self->_parse_bases($buffer);

		# now go and get the base information
	$offset = $self->{comments_offset};
	$length = $self->{comment_size};
	print("Reading the comments: Going to $offset to read $length bytes.\n");
	seek $fh,$offset,0;
	print("Just before the read, the curent filehandle position is: ".tell($fh)."\n");
	read $fh, $buffer, $length;
	unless (length($buffer) == $length) {
			$self->throw('Unexpected end of file while reading from SCF file');
	}
	$self->_parse_comments($buffer);
		# my @qualities = @{$self->{qualities}};
		# print("The ref for selfqualities is ".ref($self->{qualities})." while the ref for \@qualities is ".ref(@qualities)."\n");
	my $swq = Bio::Seq::SeqWithQuality->new(	-seq	=>	$self->{parsed}->{sequence},
							-qual	=>	$self->{parsed}->{qualities},
							-id	=>	$self->{comments}->{NAME});
	return $swq;
}

=head2 _parse_comments($buffer)                                                                                                            
                      
Title   : _parse_comments($buffer)
Usage   : $self->_parse_comments($buffer);
Function: Gather the comments section from the scf and parse it into its
	components.
Returns : Nothing. Modifies $self.
Args    : The buffer. It is expected that the buffer contains a binary
	string for the comments section of an scf file according to the
	scf file specifications.
Notes   : None. Works like Jello.

=cut

sub _parse_comments {
	my ($self,$buffer) = @_;
	my $size = length($buffer);
	my $comments_retrieved = unpack "a$size",$buffer;
	$comments_retrieved =~ s/\0//;
	my @comments_split = split/\n/,$comments_retrieved;
	foreach (@comments_split) {
		/(\w+)=(.*)/;
		$self->{comments}->{$1} = $2;
	}

	return;
}

=head2 _parse_header()

 Title   : _parse_header($buffer)
 Usage   : $self->_parse_header($buffer);
 Function: Gather the header section from the scf and parse it into its
        components.
 Returns : Nothing. Modifies $self.
 Args    : The buffer. It is expected that the buffer contains a binary
        string for the header section of an scf file according to the
        scf file specifications.
 Notes   : None.

=cut

sub _parse_header {
	my ($self,$buffer) = @_;
	($self->{scf},
		$self->{samples},
		$self->{sample_offset},
		$self->{bases},
		$self->{bases_left_clip},
		$self->{bases_right_clip},
		$self->{bases_offset},
		$self->{comment_size},
		$self->{comments_offset},
		$self->{version},
		$self->{sample_size},
		$self->{code_set},
		@{$self->{header_spare}} ) = unpack "a4 NNNNNNNN a4 NN N20", $buffer;
	return;
}

=head2 _parse_bases($buffer)

 Title   : _parse_bases($buffer)
 Usage   : $self->_parse_bases($buffer);
 Function: Gather the bases section from the scf and parse it into its
        components.
 Returns : Nothing. Modifies $self.
 Args    : The buffer. It is expected that the buffer contains a binary
        string for the bases section of an scf file according to the
        scf file specifications.
 Notes   : None.

=cut

sub _parse_bases {
	my ($self,$buffer) = @_;
	my $length = length($buffer);
	my ($offset2,$currbuff,$currbase,$currqual,$sequence,@qualities);
	my @read;
	for ($offset2=0;$offset2<$length;$offset2+=12) {
		@read = unpack "N C C C C a C3", substr($buffer,$offset2,$length);
		# print("\@read is @read\n");
		$currbase = uc($read[5]);
		if ($currbase eq "A") { $currqual = $read[1]; }
		if ($currbase eq "C") { $currqual = $read[2]; }
		if ($currbase eq "G") { $currqual = $read[3]; }
		if ($currbase eq "T") { $currqual = $read[4]; }
		$sequence .= $currbase;
		push @qualities,$currqual;
	}

	print("\$sequence is $sequence\n\@qualities are @qualities\n");
	$self->{parsed}->{sequence} = $sequence;
		# $self->{parsed}->{qualities} = \@qualities;
		# for debugging
	$self->{parsed}->{qualities} = join(' ',@qualities);
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

=head2 _split_traces(\@traces_array)

 Title   : _split_traces(\@traces_array)
 Usage   : $self->_split_traces(\@traces_array);
 Function: Parses an scf Version2 trace array into its base components.
 Returns : Nothing. Modifies $self.
 Args    : A reference to an array of the unpacked traces section of an
        scf version2 file.
 Notes   : There really should be a _split_traces_3 (for scf Version3) but
        I haven't run into any here so I didn't bother to write it.

=cut

sub _split_traces {
	my ($self,$rread) = @_;
	my @read = @$rread;
	my $array = 0;
	for (my $offset2 = 0; $offset2< scalar(@read); $offset2+=4) {
		# print("\$offset2 is $offset2\n");
		# print(int($offset2/4)."\t$read[$offset2] $read[$offset2+1] $read[$offset2+2] $read[$offset2+3]\n");
		if ($array) {
			push @{$self->{traces}->{A}},$read[$offset2];
			push @{$self->{traces}->{C}},$read[$offset2+1];
			push @{$self->{traces}->{G}},$read[$offset2+3];
			push @{$self->{traces}->{T}},$read[$offset2+2];
		}
		else {
			$self->{traces}->{A} .= " ".$read[$offset2];
			$self->{traces}->{C} .= " ".$read[$offset2+1];
			$self->{traces}->{G} .= " ".$read[$offset2+2];
			$self->{traces}->{T} .= " ".$read[$offset2+3];
		}
	}
	# my @ta = split(' ',$self->{traces}->{A});
	# my @tc = split(' ',$self->{traces}->{C});
	# my @tg = split(' ',$self->{traces}->{G});
	# my @tt = split(' ',$self->{traces}->{T});
	# print("Total number of samples is ".(scalar(@ta)+scalar(@tc)+scalar(@tg)+scalar(@tt))."\n");
	# &_dump_traces(\@ta,\@tc,\@tg,\@tt);
	return;
}

=head2 get_trace($base_channel)

 Title   : get_trace($base_channel)
 Usage   : @a_trace = @{$obj->get_trace("A")};
 Function: Return the trace data for the given base.
 Returns : A reference to an array containing the trace data for the given
        base.
 Args    : A,C,G, or T. Any other input throws.
 Notes   :

=cut

sub get_trace {
	my ($self,$base_channel) = @_;
	$base_channel =~ tr/a-z/A-Z/;
	print("get_trace: you asked for the colour channel for $base_channel\n");
	if ($base_channel !~ /A|T|G|C/) {
		$self->throw("You tried to ask for a base channel that wasn't A,T,G, or C. Ask for one of those next time.");
	}
	elsif ($base_channel) {
		return $self->{traces}->{$base_channel};
	}
}

=head2 _dump_traces($ra,$rc,$rg,$rt)

 Title   : _dump_traces($ra,$rc,$rg,$rt)
 Usage   : &_dump_traces($ra,$rc,$rg,$rt);
 Function: Used in debugging. Prints all traces one beside each other.
 Returns : Nothing.
 Args    : References to the arrays containing the traces for A,C,G, and
        T.
 Notes   : Beats using dumpValue, I'll tell ya. Much better then using
        join' ' too.

=cut

sub _dump_traces {
	my ($ra,$rc,$rg,$rt) = @_;
	my @as = @$ra; my @cs = @$rc; my @gs = @$rg; my @ts = @$rt;
	print("Count\ta\tc\tg\tt\n");
	for (my $curr=0; $curr < scalar(@as); $curr++) {
		print("$curr\t$as[$curr]\t$cs[$curr]\t$gs[$curr]\t$ts[$curr]\n");
	}
	return;
}

=head2 write_scf(-SeqWithQuality => $swq,<comments>)

 Title   : write_scf(-SeqWithQuality => $swq, <comments>)
 Usage   : $obj->write_scf(	-SeqWithQuality => $swq,
				-CONV => "Bioperl-Chads Mighty SCF writer.");
 Function: Write out an scf.
 Returns : Nothing.
 Args    : Requires: a reference to a SeqWithQuality object to form the
        basis for the scf. Any other arguments are assumed to be comments
        and are put into the comments section of the scf. Read the
	specifications for scf to decide what might be good to put in here.
 Notes   :
        Someday: (All of this stuff is easy easy easy I just don't have
                the requirement or the time.)
                - Change the peak scaling factor?
                - Change the width of the peak?
                - Change the overlap between peaks?

=cut

sub write_scf {
	# modeled a bit after SeqIO::fasta.pm, but...
	my ($self,%args) = @_;
	my %comments;
	my ($label,$arg);
	my $swq = $args{-SeqWithQuality};
	unless (ref($swq) eq "Bio::Seq::SeqWithQuality") {
		$self->throw("You must pass a Bio::Seq::SeqWithQuality object to write_scf as a parameter named \"SeqWithQuality\"");
	}
		# print("write_scf!!! Woowoo. Received a swq object. It is $swq\n");
		# print("Here are the args:\n");
		# all of the rest of the arguments are comments, at the moment
	foreach $arg (sort keys %args) {
		next if ($arg =~ /SeqWithQuality/);
		($label = $arg) =~ s/^\-//;
		$comments{$label} = %args->{$arg};
	}
	if (!$comments{NAME}) { $comments{NAME} = $swq->id(); }
		# HA! Bwahahahaha.
	$comments{CONV} = "Bioperl-Chads Mighty SCF writer.";

		# set a few things in the header
	$self->{header}->{magic} = ".scf";
	$self->{header}->{sample_size} = "2";
	$self->{header}->{bases} = length($swq->seq());
	$self->{header}->{bases_left_clip} = "0";
	$self->{header}->{bases_right_clip} = "0";
	$self->{header}->{version} = "2.00";
	$self->{header}->{sample_size} = "2";
	$self->{header}->{code_set} = "9";
	@{$self->{header}->{spare}} = qw(0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0);
 
	my $b_comments = $self->_give_comments_binary(\%comments);
	my ($b_traces,$b_bases);
	($b_traces,$b_bases,$self->{header}->{samples}) = $self->_give_tracesbases_binary($swq->seq(),$swq->qual());
	$self->{header}->{samples_offset} = "128";
	my $samples_size = $self->{header}->{samples} * 4 * $self->{header}->{sample_size};
	my $bases_size = length($swq->seq()) * 12;
	$self->{header}->{bases_offset} = 128 + length($b_traces);
	$self->{header}->{comments_offset} = 128 + length($b_traces) + length($b_bases); 
	$self->{header}->{comments_size} = length($b_comments);
	$self->{header}->{private_size} = "0";
	$self->{header}->{private_offset} = 128 + $samples_size + $bases_size + $self->{header}->{comments_size};

	my $b_header = $self->_give_header_binary();

	print("Lengths:\n");
	print("Header  : ".length($b_header)."\n");
	print("Traces  : ".length($b_traces)."\n");
	print("Bases   : ".length($b_bases)."\n");
	print("Comments: ".length($b_comments)."\n");
	my $fh = $self->_filehandle();
	print $fh ($b_header) or return;
	print $fh ($b_traces) or return;
	print $fh ($b_bases) or return;
	print $fh ($b_comments) or return;

}

=head2 _give_header_binary();

 Title   : _give_header_binary();
 Usage   : $self->_give_header_binary();
 Function: Provide the binary string that will be used as the header for a
	scfv2 document.
 Returns : A binary string.
 Args    : None. Uses the entries in the $self->{header} hash. These are set
	on construction of the object (hopefully correctly!).
 Notes   : 

=cut


sub _give_header_binary {
	my ($self,$r_header) = @_;
	my $binary = pack "a4 NNNNNNNN a4 NN N20", (
	        $self->{header}->{magic},
	        $self->{header}->{samples},
	        $self->{header}->{samples_offset},
	        $self->{header}->{bases},
	        $self->{header}->{bases_left_clip},
		$self->{header}->{bases_right_clip},
       		$self->{header}->{bases_offset},
        	$self->{header}->{comments_size},
        	$self->{header}->{comments_offset},
        	$self->{header}->{version},
        	$self->{header}->{sample_size},
        	$self->{header}->{code_set},
        	@{$self->{header}->{spare}});
	return $binary;
}

=head2 _give_tracesbases_binary($sequence,\@wqualities);

 Title   : _give_tracesbases_binary($sequence,\@wqualities)
 Usage   : $self->_give_tracesbases_binary($sequence,\@wqualities);
 Function: Return the binary strings for the samples and the bases sections
	of an scfv2 file.
 Returns : 	1. A binary string for the traces section of the scfv2 file
		2. A binary string for the bases section of the scfv2 file
		3. The length of the trace string.
 Args    : The sequence and a reference to an array of quality values.
 Notes   : Get the args for this routine from the SeqWithQuality object to
	make sure that they are valid. Don't call this outside of new().
	Really.

=cut

sub _give_tracesbases_binary {
	my ($self,$sequence,$rqual) = @_;
	my $samples = {};
	$sequence =~ tr/a-z/A-Z/;
        $samples->{sequence} = $sequence;
        $samples->{sequence_length} = length($sequence);	
	my @quals = @$rqual;
		# build the ramp for the first base.
		# a ramp looks like this "1 4 13 29 51 71 80 71 51 29 13 4 1" times the quality score.
		# REMEMBER: A C G T
		# note to self-> smooth this thing out a bit later
        @{$samples->{ramp}} = qw( 1 4 13 29 51 75 80 75 51 29 13 4 1 );
		# the width of the ramp
        $samples->{ramp_width} = scalar(@{$samples->{ramp}});
		# how far should the peaks overlap?
	$samples->{overlap} = 1;
		# where should the peaks be located?
	$samples->{peak_at} = 7;
	$samples->{total_length} = $samples->{sequence_length} * $samples->{ramp_width} - $samples->{sequence_length} * $samples->{overlap};
		# create some empty arrays
		# my (@sam_a,@sam_c,@sam_g,@sam_t,$pos);
	my $pos;
	my $total_length = $samples->{total_length};
	for ($pos=0;$pos<$total_length;$pos++) {
		$samples->{arrays}->{sam_a}[$pos] = $samples->{arrays}->{sam_c}[$pos] = $samples->{arrays}->{sam_g}[$pos] = $samples->{arrays}->{sam_t}[$pos] = "0";
	}
		# now populate them
	my ($current_base,$place_base_at,$peak_quality,$ramp_counter,$current_ramp,$ramp_position);
		# print("sequence is $samples->{sequence}\n");
	my $sequence_length = $samples->{sequence_length};
	my $half_ramp = int($samples->{ramp_width}/2);
	for ($pos = 0; $pos<$sequence_length;$pos++) {
		$current_base = substr($samples->{sequence},$pos,1);
			# print("Current base is $current_base\n");
			# where should the peak for this base be placed? Modeled after a mktrace scf
		$place_base_at = ($pos * $samples->{ramp_width}) - ($pos * $samples->{overlap}) - $half_ramp + $samples->{ramp_width} - 1;
			# print("Placing base $current_base at $place_base_at.\n");
			# next;
		$peak_quality = $quals[$pos];
			# print("$current_base has a quality $peak_quality\n");
		if ($current_base eq "A") {
			$ramp_position = $place_base_at - $half_ramp;
			for ($current_ramp = 0; $current_ramp < $samples->{ramp_width}; $current_ramp++) {
				$samples->{arrays}->{sam_a}[$ramp_position+$current_ramp] = $peak_quality * $samples->{ramp}[$current_ramp];
			}
				# print(($pos+12)."\t".$peak_quality."\t0\t0\t0\t".$current_base."\t0\t0\t0\n");
			push @{$samples->{arrays}->{all_bases}},($place_base_at+1,$peak_quality,0,0,0,$current_base,0,0,0);
		}
		elsif ($current_base eq "C") {
			$ramp_position = $place_base_at - $half_ramp;
			for ($current_ramp = 0; $current_ramp < $samples->{ramp_width}; $current_ramp++) {
				$samples->{arrays}->{sam_c}[$ramp_position+$current_ramp] = $peak_quality * $samples->{ramp}[$current_ramp];
			}
				# print(($pos+12)."\t0\t".$peak_quality."\t0\t0\t".$current_base."\t0\t0\t0\n");
			push @{$samples->{arrays}->{all_bases}},($place_base_at+1,0,$peak_quality,0,0,$current_base,0,0,0);
		}
		elsif ($current_base eq "G") {
			$ramp_position = $place_base_at - $half_ramp;
			for ($current_ramp = 0; $current_ramp < $samples->{ramp_width}; $current_ramp++) {
				$samples->{arrays}->{sam_g}[$ramp_position+$current_ramp] = $peak_quality * $samples->{ramp}[$current_ramp];
			}
				# print(($pos+12)."\t0\t0\t".$peak_quality."\t0\t".$current_base."\t0\t0\t0\n");
			push @{$samples->{arrays}->{all_bases}},($place_base_at+1,0,0,$peak_quality,0,$current_base,0,0,0);
		}
		elsif ($current_base eq "T") {
			$ramp_position = $place_base_at - $half_ramp;
			for ($current_ramp = 0; $current_ramp < $samples->{ramp_width}; $current_ramp++) {
				$samples->{arrays}->{sam_t}[$ramp_position+$current_ramp] = $peak_quality * $samples->{ramp}[$current_ramp];
			}
				# print(($pos+12)."\t0\t0\t0\t".$peak_quality."\t".$current_base."\t0\t0\t0\n");
			push @{$samples->{arrays}->{all_bases}},($place_base_at+1,0,0,0,$peak_quality,$current_base,0,0,0);
		}
		else {
			print("The current base is not a base. Hmmm.\n");
		}
	}
		# dumpValue($samples);
		# &dump_traces(\@sam_a,\@sam_c,\@sam_g,\@sam_t);
		# return $trace_string,\@traces_view,$length,\@traces;
	($samples->{strings}->{trace_string}->{string},$samples->{arrays}->{samples_view},$samples->{strings}->{trace_string}->{length},$samples->{arrays}->{traces}) = &_make_trace_string($samples->{arrays}->{sam_a},$samples->{arrays}->{sam_c},$samples->{arrays}->{sam_g},$samples->{arrays}->{sam_t});
	$samples->{strings}->{bases} = join(' ',@{$samples->{arrays}->{all_bases}});
		# $self->{strings}->{bases} = $samples->{strings}->{bases};
		# my @traces_to_bin = split(' ',$trace_string);
	@{$samples->{arrays}->{all_traces}} = split(' ',$samples->{strings}->{trace_string}->{string});
		# print("There are ".scalar(@{$samples->{bases}})." bases in the base string.\n")
	my ($packstring,@pack_array,$pos2,$tester,@unpacked);
	for ($pos = 0; $pos<$sequence_length;$pos++) {
		my @pack_array = @{$samples->{arrays}->{all_bases}}[$pos*9..$pos*9+8];
			# print("Trying to pack @pack_array\n");
			# $tester = pack "N C C C C a C3",@pack_array;
			# print("Now trying to _un_pack the same string.\n");
			# @unpacked = unpack "N C C C C a C3",$tester;
			# print("The unpacked string is @unpacked\n");

		$samples->{binaries}->{bases} .= pack "N C C C C a C3",@pack_array;
	}
		print("binary bases string has length ".length($samples->{binaries}->{bases})."\n");
		# $samples->{bases} = undef;
	my $trace_pack_length = $samples->{strings}->{trace_string}->{length} * 4 * $self->{header}->{sample_size};
	print("\$length of trace pack is ".$trace_pack_length."\n");
		# $samples->{binaries}->{traces} = pack "n
	$samples->{binaries}->{traces} .= pack "n$trace_pack_length",@{$samples->{arrays}->{all_traces}};
		# print("\$trace_string is $trace_string.\n");	
	print("Length of binary bases is ".length($samples->{binaries}->{bases})."\n");
	return ($samples->{binaries}->{traces},$samples->{binaries}->{bases},$samples->{strings}->{trace_string}->{length});

}

=head2 _make_trace_string(\@as,\@cs,\@gs,\@ts)

 Title   : _make_trace_string(\@as,\@cs,\@gs,\@ts)
 Usage   : $self->_make_trace_string(\@as,\@cs,\@gs,\@ts)
 Function: Merges trace data for the four bases to produce an scf version2
	trace string.
 Returns :
	1. A string containing a join' ' of the array of trace data
	2. A reference to an array containing trace data that can be
		easily viewed
	3. The length of the trace string.
	4. A reference to the array containing the scfv2 trace data.
 Args    : References to four arrays containing trace data.
 Notes   : This was a nightmare to create. The explains the gratuitous debug
	information. Note to self-> remove it.

=cut

sub _make_trace_string {
	my ($ra,$rc,$rg,$rt) = @_;
	my @traces;
	my @traces_view;
	my @as = @$ra; my @cs = @$rc; my @gs = @$rg; my @ts = @$rt;
        for (my $curr=0; $curr < scalar(@as); $curr++) {
		push @traces,($as[$curr],$cs[$curr],$gs[$curr],$ts[$curr]);
		push @traces_view,($as[$curr]." ".$cs[$curr]." ".$gs[$curr]." ".$ts[$curr]);
	}
	my $length = scalar(@traces)/4;
	my $trace_string = join' ',@traces;
	return $trace_string,\@traces_view,$length,\@traces;
}	

=head2 _give_comments_binary(\@comments);

 Title   : _give_comments_binary(\@comments)
 Usage   : $self->_give_comments_binary(\@comments);
 Function: Provide a binary string that will be the comments section of the
	scf file. See the scf specifications for detailed specifications for
	the comments section of an scf file. Hint:
		CODE=something\nBODE=something\n\0
 Returns : The binary string for the coments section.
 Args    : A reference to an array containing comments.
 Notes   : None.

=cut

sub _give_comments_binary {
	my ($self,$rcomments) = @_;
	my %comments = %$rcomments;
	my $comments_string;
	foreach (sort keys %comments) {
		$comments_string .= "$_=".$comments{$_}."\n";
	}
	$comments_string .= "\n\0";
		# print("Comments string is $comments_string\n");
	my $length = length($comments_string);
	my $binary .= pack "A$length",$comments_string;
	return $binary;
}

1;
__END__

