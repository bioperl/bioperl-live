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

Do not use this module directly. Use it via the Bio::SeqIO class, see
L<Bio::SeqIO> for more information.

=head1 DESCRIPTION

This object can transform .scf files to and from
Bio::Seq::SeqWithQuality objects.  Mechanisms are present to retrieve
trace data from scf files.

=head1 FEEDBACK

=head2 Mailing Lists


User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://www.bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

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

package Bio::SeqIO::scf;
use vars qw(@ISA $DEFAULT_QUALITY);
use strict;
use Bio::SeqIO;
use Bio::Seq::SeqFactory;
require 'dumpvar.pl';

BEGIN { 
    $DEFAULT_QUALITY= 10;
}

@ISA = qw(Bio::SeqIO);

sub _initialize {
  my($self,@args) = @_;
  $self->SUPER::_initialize(@args);    
  if( ! defined $self->sequence_factory ) {
      $self->sequence_factory(new Bio::Seq::SeqFactory
			      (-verbose => $self->verbose(), 
			       -type => 'Bio::Seq::SeqWithQuality'));      
  }
}

=head2 next_seq()

 Title   : next_seq()
 Usage   : $scf = $stream->next_seq()
 Function: returns the next scf sequence in the stream
 Returns : Bio::Seq::SeqWithQuality object
 Args    : NONE
 Notes   : Fills the interface specification for SeqIO.
	   The SCF specification does not provide for having more then
           one sequence in a given scf. So once the filehandle has been open
           and passed to SeqIO don't expect to run this function more then
           once on a given scf unless you embraced and extended the SCF
  	   standard. (But that's just C R A Z Y talk, isn't it.)

=cut

#'
sub next_seq {
    my ($self) = @_;
    my ($seq, $seqc, $fh, $buffer, $offset, $length, $read_bytes, @read,
	%names);
    # set up a filehandle to read in the scf    
    $fh = $self->_filehandle();
    unless ($fh) {		# simulate the <> function
	if ( !fileno(ARGV) or eof(ARGV) ) {
	    return unless my $ARGV = shift;
	    open(ARGV,$ARGV) or
		$self->throw("Could not open $ARGV for SCF stream reading $!");
	}
	$fh = \*ARGV;
    }
    binmode $fh;		# for the Win32/Mac crowds
    return unless read $fh, $buffer, 128; # no exception; probably end of file
    # the first thing to do is parse the header. This is common
    # among all versions of scf.
    $self->_set_header($buffer);
    # the rest of the the information is different between the
    # the different versions of scf.
    if ($self->{'version'} lt "3.00") {
	# first gather the trace information
	$length = $self->{'samples'}*$self->{sample_size}*4;
	$buffer = $self->read_from_buffer($fh,$buffer,$length);
	@read = unpack "n$length",$buffer;
	# these traces need to be split
	$self->_set_v2_traces(\@read);
	# now go and get the base information
	$offset = $self->{bases_offset};
	$length = ($self->{bases} * 12);
	seek $fh,$offset,0;
	$buffer = $self->read_from_buffer($fh,$buffer,$length);
	# now distill the information into its fractions.
	$self->_set_v2_bases($buffer);
    }
    else {
	my $transformed_read;
	foreach (qw(A C G T)) {
	    $length = $self->{'samples'}*$self->{sample_size};
	    $buffer = $self->read_from_buffer($fh,$buffer,$length);
	    @read = unpack "n$length",$buffer;
				# this little spurt of nonsense is because
				# the trace values are given in the binary
				# file as unsigned shorts but they really
				# are signed. 30000 is an arbitrary number
				# (will there be any traces with a given
				# point greater then 30000? I hope not.
				# once the read is read, it must be changed
				# from relative 
	    for (my $element=0; $element < scalar(@read); $element++) {
		if ($read[$element] > 30000) {
		    $read[$element] = $read[$element] - 65536;
		}
	    }
	    $transformed_read = $self->_delta(\@read,"backward");
	    $self->{'traces'}->{$_} = join(' ',@{$transformed_read});
	}
	# now go and get the peak index information
	$offset = $self->{bases_offset};
	$length = ($self->{bases} * 4);
	seek $fh,$offset,0;
	$buffer = $self->read_from_buffer($fh,$buffer,$length);
	$self->_set_v3_peak_indices($buffer);
	# now go and get the accuracy information
	$buffer = $self->read_from_buffer($fh,$buffer,$length);
	$self->_set_v3_base_accuracies($buffer);
	# OK, now go and get the base information.
	$length = $self->{bases};
	$buffer = $self->read_from_buffer($fh,$buffer,$length);
	$self->{'parsed'}->{'sequence'} = unpack("a$length",$buffer);
	# now, finally, extract the calls from the accuracy information.
	$self->_set_v3_quality($self);
    }
    # now go and get the comment information
    $length = $self->{comment_size};
    $buffer = $self->read_from_buffer($fh,$buffer,$length);
    $self->_set_comments($buffer);
    return $self->sequence_factory->create
	(-seq  =>	$self->{'parsed'}->{'sequence'},
	 -qual =>	$self->{'parsed'}->{'qualities'},
	 -id   =>	$self->{'comments'}->{'NAME'}
	 );
}

=head2 _set_v3_quality()

 Title   : _set_v3_quality()
 Usage   : $self->_set_v3_quality()
 Function: Set the base qualities from version3 scf's
 Returns : Nothing. Alters $self.
 Args    : None.
 Notes   : 

=cut

#'
sub _set_v3_quality {
    my $self = shift;
    my @bases = split//,$self->{'parsed'}->{'sequence'};
    my (@qualities,$currbase,$currqual,$counter);
    for ($counter=0; $counter <= $#bases ; $counter++) {
	$currbase = uc($bases[$counter]);
	if ($currbase eq "A") { $currqual = $self->{'parsed'}->{'base_accuracies'}->{'A'}->[$counter]; }
	elsif ($currbase eq "C") { $currqual = $self->{'parsed'}->{'base_accuracies'}->{'C'}->[$counter]; }
	elsif ($currbase eq "G") { $currqual = $self->{'parsed'}->{'base_accuracies'}->{'G'}->[$counter]; }
	elsif ($currbase eq "T") { $currqual = $self->{'parsed'}->{'base_accuracies'}->{'T'}->[$counter]; }
	else { $currqual = "unknown"; }
	push @qualities,$currqual;
    }
    $self->{'parsed'}->{'qualities'} = \@qualities;
}

=head2 _set_v3_peak_indices($buffer)

 Title   : _set_v3_peak_indices($buffer)
 Usage   : $self->_set_v3_peak_indices($buffer);
 Function: Unpacks the base accuracies for version3 scf
 Returns : Nothing. Alters $self
 Args    : A scalar containing binary data.
 Notes   : 

=cut

sub _set_v3_peak_indices {
    my ($self,$buffer) = @_;
    my $length = length($buffer);
    my ($offset,@read,@positions);
    @read = unpack "N$length",$buffer;
    $self->{'parsed'}->{'peak_indices'} = join(' ',@read);
}

=head2 _set_v3_base_accuracies($buffer)

 Title   : _set_v3_base_accuracies($buffer)
 Usage   : $self->_set_v3_base_accuracies($buffer)
 Function: Set the base accuracies for version 3 scf's
 Returns : Nothing. Alters $self.
 Args    : A scalar containing binary data.
 Notes   : 

=cut

#'
sub _set_v3_base_accuracies {
    my ($self,$buffer) = @_;
    my $length = length($buffer);
    my $qlength = $length/4;
    my $offset = 0;
    my (@qualities,@sorter,$counter,$round,$last_base);
    foreach (qw(A C G T)) {
	my @read;
	$last_base = $offset + $qlength;
	for (;$offset < $last_base; $offset += $qlength) {
	    @read = unpack "c$qlength", substr($buffer,$offset,$qlength);
	    $self->{'parsed'}->{'base_accuracies'}->{"$_"} = \@read;
	}
    }
}


=head2 _set_comments($buffer)

 Title   : _set_comments($buffer)
 Usage   : $self->_set_comments($buffer);
 Function: Gather the comments section from the scf and parse it into its
	   components.
 Returns : Nothing. Modifies $self.
 Args    : The buffer. It is expected that the buffer contains a binary
	   string for the comments section of an scf file according to
	   the scf file specifications.
 Notes   : None. Works like Jello.

=cut

sub _set_comments {
    my ($self,$buffer) = @_;
    my $size = length($buffer);
    my $comments_retrieved = unpack "a$size",$buffer;
    $comments_retrieved =~ s/\0//;
    my @comments_split = split/\n/,$comments_retrieved;
    if (@comments_split) {
	foreach (@comments_split) {
	    /(\w+)=(.*)/;
	    if ($1 && $2) {
		$self->{'comments'}->{$1} = $2;
	    }
	}
    }
    return;
}

=head2 _set_header()

 Title   : _set_header($buffer)
 Usage   : $self->_set_header($buffer);
 Function: Gather the header section from the scf and parse it into its
           components.
 Returns : Nothing. Modifies $self.
 Args    : The buffer. It is expected that the buffer contains a binary
           string for the header section of an scf file according to the
           scf file specifications.
 Notes   : None.

=cut

sub _set_header {
    my ($self,$buffer) = @_;
    ($self->{'scf'},
     $self->{'samples'},
     $self->{'sample_offset'},
     $self->{'bases'},
     $self->{'bases_left_clip'},
     $self->{'bases_right_clip'},
     $self->{'bases_offset'},
     $self->{'comment_size'},
     $self->{'comments_offset'},
     $self->{'version'},
     $self->{'sample_size'},
     $self->{'code_set'},
     @{$self->{'header_spare'}} ) = unpack "a4 NNNNNNNN a4 NN N20", $buffer;
    return;

}

=head2 _set_v2_bases($buffer)

 Title   : _set_v2_bases($buffer)
 Usage   : $self->_set_v2_bases($buffer);
 Function: Gather the bases section from the scf and parse it into its
           components.
 Returns : Nothing. Modifies $self.
 Args    : The buffer. It is expected that the buffer contains a binary
           string for the bases section of an scf file according to the
           scf file specifications.
 Notes   : None.

=cut

sub _set_v2_bases {
    my ($self,$buffer) = @_;
    my $length = length($buffer);
    my ($offset2,$currbuff,$currbase,$currqual,$sequence,@qualities,@indices);
    my @read;
    for ($offset2=0;$offset2<$length;$offset2+=12) {
	@read = unpack "N C C C C a C3", substr($buffer,$offset2,$length);
	push @indices,$read[0];
	$currbase = uc($read[5]);
	if ($currbase eq "A") { $currqual = $read[1]; }
	elsif ($currbase eq "C") { $currqual = $read[2]; }
	elsif ($currbase eq "G") { $currqual = $read[3]; }
	elsif ($currbase eq "T") { $currqual = $read[4]; }
	else { $currqual = "UNKNOWN"; }
	$sequence .= $currbase;
	push @qualities,$currqual;
    }
    unless (!@indices) {
	$self->{'parsed'}->{'peak_indices'} = join(' ',@indices);
    }
    $self->{'parsed'}->{'sequence'} = $sequence;
    unless (!@qualities) {
	$self->{'parsed'}->{'qualities'} = join(' ',@qualities);
    }
}

=head2 _set_v2_traces(\@traces_array)

 Title   : _set_v2_traces(\@traces_array)
 Usage   : $self->_set_v2_traces(\@traces_array);
 Function: Parses an scf Version2 trace array into its base components.
 Returns : Nothing. Modifies $self.
 Args    : A reference to an array of the unpacked traces section of an
           scf version2 file.

=cut

sub _set_v2_traces {
    my ($self,$rread) = @_;
    my @read = @$rread;
    my $array = 0;
    for (my $offset2 = 0; $offset2< scalar(@read); $offset2+=4) {
	if ($array) {
	    push @{$self->{'traces'}->{'A'}},$read[$offset2];
	    push @{$self->{'traces'}->{'C'}},$read[$offset2+1];
	    push @{$self->{'traces'}->{'G'}},$read[$offset2+3];
	    push @{$self->{'traces'}->{'T'}},$read[$offset2+2];
	} else {
	    $self->{'traces'}->{'A'} .= " ".$read[$offset2];
	    $self->{'traces'}->{'C'} .= " ".$read[$offset2+1];
	    $self->{'traces'}->{'G'} .= " ".$read[$offset2+2];
	    $self->{'traces'}->{'T'} .= " ".$read[$offset2+3];
	}
    }
    return;
}

=head2 get_trace($base_channel)

 Title   : get_trace($base_channel)
 Usage   : @a_trace = @{$obj->get_trace("A")};
 Function: Return the trace data for the given base.
 Returns : A reference to an array containing the trace data for the
	   given base.
 Args    : A,C,G, or T. Any other input throws.
 Notes   :

=cut

sub get_trace {
    my ($self,$base_channel) = @_;
    $base_channel =~ tr/a-z/A-Z/;
    if ($base_channel !~ /A|T|G|C/) {
	$self->throw("You tried to ask for a base channel that wasn't A,T,G, or C. Ask for one of those next time.");
    } elsif ($base_channel) {
	my @temp = split(' ',$self->{'traces'}->{$base_channel});
	return \@temp;
    }
}

=head2 get_peak_indices()

 Title   : get_peak_indices()
 Usage   : @a_trace = @{$obj->get_peak_indices()};
 Function: Return the peak indices for this scf.
 Returns : A reference to an array containing the peak indices for this scf. 
 Args    : None.
 Notes   :

=cut

sub get_peak_indices {
    my ($self) = shift;
    my @temp = split(' ',$self->{'parsed'}->{'peak_indices'});
    return \@temp;
}


=head2 get_header()

 Title   : get_header()
 Usage   : %header = %{$obj->get_header()};
 Function: Return the header for this scf.
 Returns : A reference to a hash containing the header for this scf. 
 Args    : None.
 Notes   :

=cut

sub get_header {
    my ($self) = shift;
    my %header;
    foreach (qw(scf samples sample_offset bases bases_left_clip 
		bases_right_clip bases_offset comment_size comments_offset 
		version sample_size code_set peak_indices)) {
	%header->{"$_"} = $self->{"$_"};
    }
    return \%header;
}

=head2 _dump_traces_incoming($transformed)

 Title   : _dump_traces_incoming("transformed")
 Usage   : &_dump_traces($ra,$rc,$rg,$rt);
 Function: Used in debugging. Prints all traces one beside each other.
 Returns : Nothing.
 Args    : References to the arrays containing the traces for A,C,G,T.
 Notes   : Beats using dumpValue, I'll tell ya. Much better then using
           join' ' too.
	- if a scalar is included as an argument (any scalar), this
	procedure will dump the _delta'd trace. If you don't know what
	that means you should not be using this.

=cut

#'
sub _dump_traces_incoming {
    my ($self) = @_;
    my (@sA,@sT,@sG,@sC);
    # @sA = @{$self->{'traces'}->{'A'}};
    # @sC = @{$self->{'traces'}->{'C'}};
    # @sG = @{$self->{'traces'}->{'G'}};
    # @sT = @{$self->{'traces'}->{'T'}};
    @sA = @{$self->get_trace('A')};
    @sC = @{$self->get_trace('C')};
    @sG = @{$self->get_trace('G')};
    @sT = @{$self->get_trace('t')};
    print ("Count\ta\tc\tg\tt\n");
    for (my $curr=0; $curr < scalar(@sG); $curr++) {
	print("$curr\t$sA[$curr]\t$sC[$curr]\t$sG[$curr]\t$sT[$curr]\n");
    }
    return;
}

=head2 _dump_traces_outgoing($transformed)

 Title   : _dump_traces_outgoing("transformed")
 Usage   : &_dump_traces_outgoing(($ra,$rc,$rg,$rt);
 Function: Used in debugging. Prints all traces one beside each other.
 Returns : Nothing.
 Args    : References to the arrays containing the traces for A,C,G,T.
 Notes   : Beats using dumpValue, I\'ll tell ya. Much better then using
           join' ' too.
	- if a scalar is included as an argument (any scalar), this
	procedur will dump the _delta'd trace. If you don't know what
	that means you should not be using this.

=cut

sub _dump_traces_outgoing {
    my ($self,$transformed) = @_;
    my (@sA,@sT,@sG,@sC);
    if ($transformed) {
	@sA = @{$self->{'text'}->{'t_samples_a'}};
	@sC = @{$self->{'text'}->{'t_samples_c'}};
	@sG = @{$self->{'text'}->{'t_samples_g'}};
	@sT = @{$self->{'text'}->{'t_samples_t'}};
    }
    else {
	@sA = @{$self->{'text'}->{'samples_a'}};
	@sC = @{$self->{'text'}->{'samples_c'}};
	@sG = @{$self->{'text'}->{'samples_g'}};
	@sT = @{$self->{'text'}->{'samples_t'}};
    }
    print ("Count\ta\tc\tg\tt\n");
    for (my $curr=0; $curr < scalar(@sG); $curr++) {
	print("$curr\t$sA[$curr]\t$sC[$curr]\t$sG[$curr]\t$sT[$curr]\n");
    }
    return;
}

=head2 write_seq

 Title   : write_seq(-SeqWithQuality => $swq, <comments>)
 Usage   : $obj->write_seq(	-SeqWithQuality => $swq,
			-version => 2,
			-CONV => "Bioperl-Chads Mighty SCF writer.");
 Function: Write out an scf.
 Returns : Nothing.
 Args    : Requires: a reference to a SeqWithQuality object to form the
           basis for the scf. 
	   if -version is provided, it should be "2" or "3". A SCF of that
	   version will be written.
	   Any other arguments are assumed to be comments and are put into
	   the comments section of the scf. Read the specifications for scf
	   to decide what might be good to put in here.

 Notes   :
        Someday: (All of this stuff is easy easy easy I just don't have
                the requirement or the time.)
                - Change the peak scaling factor?
                - Change the width of the peak?
                - Change the overlap between peaks?

=cut

#'
sub write_seq {
    my ($self,%args) = @_;
    my %comments;
    my ($label,$arg);

    my ($swq) = $self->_rearrange([qw(SEQWITHQUALITY)], %args);
    unless (ref($swq) eq "Bio::Seq::SeqWithQuality") {
	$self->throw("You must pass a Bio::Seq::SeqWithQuality object to write_seq as a parameter named \"SeqWithQuality\"");
    }
    # verify that there is some sequence or some qualities
    # If the $swq with quality has no qualities, set all qualities to 0.
    # If the $swq has no sequence, set the sequence to N\'s.
    $self->_fill_missing_data($swq);

    # all of the rest of the arguments are comments for the scf
    foreach $arg (sort keys %args) {
	next if ($arg =~ /SeqWithQuality/i);
	($label = $arg) =~ s/^\-//;
	$comments{$label} = %args->{$arg};
    }
    if (!$comments{'NAME'}) { $comments{'NAME'} = $swq->id(); }
    # HA! Bwahahahaha.
    $comments{'CONV'} = "Bioperl-Chads Mighty SCF writer." unless defined $comments{'CONV'};
    # now deal with the version of scf they want to write
    if (%comments->{version}) {
	if (%comments->{version} != 2 && %comments->{version} != 3) {
	    $self->warn("This module can only write version 2.0 or 3.0 scf's. Writing a version 2.0 scf by default.");
	    %comments->{version} = "2.00";
	}
	if (%comments->{'version'} > 2) {
	    %comments->{'version'} = "3.00";
	}
    }
    else {
	%comments->{'version'} = "2.00";
    }



    # set a few things in the header
    $self->{'header'}->{'magic'} = ".scf";
    $self->{'header'}->{'sample_size'} = "2";
    $self->{'header'}->{'bases'} = length($swq->seq());
    $self->{'header'}->{'bases_left_clip'} = "0";
    $self->{'header'}->{'bases_right_clip'} = "0";
    $self->{'header'}->{'version'} = %comments->{'version'};
    $self->{'header'}->{'sample_size'} = "2";
    $self->{'header'}->{'code_set'} = "9";
    @{$self->{'header'}->{'spare'}} = qw(0 0 0 0 0 0 0 0 0 0 
					 0 0 0 0 0 0 0 0 0 0);

    # create the binary for the comments and file it in $self->{'binaries'}->{'comments'}
    $self->_set_binary_comments(\%comments);
    # create the binary and the strings for the traces, bases, offsets (if necessary), and accuracies (if necessary)
    $self->_set_binary_tracesbases(%comments->{'version'},$swq->seq(),$swq->qual());

    # now set more things in the header
    $self->{'header'}->{'samples_offset'} = "128";

    my ($b_base_offsets,$b_base_accuracies,$samples_size,$bases_size);
    #
    # version 2
    #
    if ($self->{'header'}->{'version'} == 2) {
	$samples_size = $self->{'header'}->{'samples'} * 4 * 
	    $self->{'header'}->{'sample_size'};
	$bases_size = length($swq->seq()) * 12;
	$self->{'header'}->{'bases_offset'} = 128 + length($self->{'binaries'}->{'samples_all'});
	$self->{'header'}->{'comments_offset'} = 128 + length($self->{'binaries'}->{'samples_all'}) + length($self->{'binaries'}->{'v2_bases'}); 
	$self->{'header'}->{'comments_size'} = length($self->{'binaries'}->{'comments'});
	$self->{'header'}->{'private_size'} = "0";
	$self->{'header'}->{'private_offset'} = 128 + $samples_size + 
	    $bases_size + $self->{'header'}->{'comments_size'};
    }
    else {
	$self->{'header'}->{'bases_offset'} = 128 + length($self->{'binaries'}->{'samples_all'});
	$self->{'header'}->{'comments_size'} = length($self->{'binaries'}->{'comments'});
	# this is:
	# bases_offset + base_offsets + accuracies + called_bases + reserved
	$self->{'header'}->{'comments_offset'} = $self->{'header'}->{'bases_offset'} + 4*$self->{header}->{'bases'} + 4*$self->{header}->{'bases'} + $self->{header}->{'bases'} + 3*$self->{header}->{'bases'}; 
	$self->{'header'}->{'private_size'} = "0";
	$self->{'header'}->{'private_offset'} = $self->{'header'}->{'comments_offset'} + $self->{'header'}->{'comments_size'};
    }

    $self->_set_binary_header();

    # should something better be done rather then returning after
    # writing? I don't do any exception trapping here
    if (%comments->{'version'} == 2) {
	# print ("Lengths:\n");
	# print("Header  : ".length($self->{'binaries'}->{'header'})."\n");
	# print("Traces  : ".length($self->{'binaries'}->{'samples_all'})."\n");
	# print("Bases   : ".length($self->{'binaries'}->{'v2_bases'})."\n");
	# print("Comments: ".length($self->{'binaries'}->{'comments'})."\n");
	$self->_print ($self->{'binaries'}->{'header'}) or return;
	$self->_print ($self->{'binaries'}->{'samples_all'}) or return;
	$self->_print ($self->{'binaries'}->{'v2_bases'}) or return;
	$self->_print ($self->{'binaries'}->{'comments'}) or return;
    }
    elsif (%comments->{'version'} ==3) {
	# print ("Lengths:\n");
	# print("Header  : ".length($self->{'binaries'}->{'header'})."\n");
	# print("Traces  : ".length($self->{'binaries'}->{'samples_all'})."\n");
	# print("Offsets : ".length($self->{'binaries'}->{'v3_peak_offsets'})."\n");
	# print("Accuracy: ".length($self->{'binaries'}->{'v3_accuracies_all'})."\n");
	# print("Bases   : ".length($self->{'binaries'}->{'v3_called_bases'})."\n");
	# print("Reserved: ".length($self->{'binaries'}->{'v3_reserved'})."\n");
	# print("Comments: ".length($self->{'binaries'}->{'comments'})."\n");
	$self->{'header'}->{'comments_offset'} = 
	    128+length($self->{'binaries'}->{'samples_all'})+
		length($self->{'binaries'}->{'v3_peak_offsets'})+
		    length($self->{'binaries'}->{'v3_accuracies_all'})+
			length($self->{'binaries'}->{'v3_called_bases'})+
			    length($self->{'binaries'}->{'v3_reserved'});
	$self->{'header'}->{'spare'}->[1] = 
	    $self->{'header'}->{'comments_offset'} +
		length($self->{'binaries'}->{'comments'});
	$self->_set_binary_header();
	$self->_print ($self->{'binaries'}->{'header'}) or print("Couldn't write header\n");
	$self->_print ($self->{'binaries'}->{'samples_all'}) or print("Couldn't write samples\n");
	$self->_print ($self->{'binaries'}->{'v3_peak_offsets'}) or print("Couldn't write peak offsets\n");
	$self->_print ($self->{'binaries'}->{'v3_accuracies_all'}) or print("Couldn't write accuracies\n");
	$self->_print ($self->{'binaries'}->{'v3_called_bases'}) or print("Couldn't write called_bases\n");
	$self->_print ($self->{'binaries'}->{'v3_reserved'}) or print("Couldn't write reserved\n");
	$self->_print ($self->{'binaries'}->{'comments'}) or print ("Couldn't write comments\n");
    }

    # kinda unnecessary, given the close() below, but maybe that'll go
    # away someday.
    $self->_fh->flush if $self->_flush_on_write && defined $self->_fh;

    $self->close();
}

=head2 _set_binary_header()

 Title   : _set_binary_header();
 Usage   : $self->_set_binary_header();
 Function: Provide the binary string that will be used as the header for
	   a scfv2 document.
 Returns : A binary string.
 Args    : None. Uses the entries in the $self->{'header'} hash. These
	   are set on construction of the object (hopefully correctly!).
 Notes   : 

=cut

sub _set_binary_header {
    my ($self) = shift;
    my $binary = pack "a4 NNNNNNNN a4 NN N20", 
    (
     $self->{'header'}->{'magic'},
     $self->{'header'}->{'samples'},
     $self->{'header'}->{'samples_offset'},
     $self->{'header'}->{'bases'},
     $self->{'header'}->{'bases_left_clip'},
     $self->{'header'}->{'bases_right_clip'},
     $self->{'header'}->{'bases_offset'},
     $self->{'header'}->{'comments_size'},
     $self->{'header'}->{'comments_offset'},
     $self->{'header'}->{'version'},
     $self->{'header'}->{'sample_size'},
     $self->{'header'}->{'code_set'},
     @{$self->{'header'}->{'spare'}});
	$self->{'binaries'}->{'header'} = $binary;
}

=head2 _set_binary_tracesbases($version,$sequence,$ref_quality)

 Title   : _set_binary_tracesbases($version,$sequence,$ref_quality)
 Usage   : $self->_set_binary_tracesbases($version,$sequence,
	   $ref_quality);
 Function: Constructs the trace and base strings for all scfs
 Returns : Nothing. Alters self.
 Args    : $version - "2" or "3"
	   $sequence - a scalar containing arbitrary sequence data
	   $ref_quality - a reference to an array containing quality
	   values
 Notes   : This is a really complicated thing.

=cut

sub _set_binary_tracesbases {
    my ($self,$version,$sequence,$rqual) = @_;
    $sequence =~ tr/a-z/A-Z/;
    $self->{'info'}->{'sequence'} = $sequence;
    $self->{'info'}->{'sequence_length'} = length($sequence);	
    my @quals = @$rqual;
	    # build the ramp for the first base.
	    # a ramp looks like this "1 4 13 29 51 71 80 71 51 29 13 4 1" times the quality score.
	    # REMEMBER: A C G T
	    # note to self-> smooth this thing out a bit later
    @{$self->{'text'}->{'ramp'}} = qw( 1 4 13 29 51 75 80 75 51 29 13 4 1 );
	    # the width of the ramp
    $self->{'text'}->{'ramp_width'} = scalar(@{$self->{'text'}->{'ramp'}});
	    # how far should the peaks overlap?
    $self->{'text'}->{'ramp_overlap'} = 1;
    # where should the peaks be located?
    $self->{'text'}->{'peak_at'} = 7;
    $self->{'text'}->{'ramp_total_length'} =
		$self->{'info'}->{'sequence_length'} * $self->{'text'}->{'ramp_width'}
		- $self->{'info'}->{'sequence_length'} * $self->{'text'}->{'ramp_overlap'};
    # create some empty arrays
    # my (@sam_a,@sam_c,@sam_g,@sam_t,$pos);
    my $pos;
    my $total_length = $self->{'text'}->{ramp_total_length};
    for ($pos=0;$pos<=$total_length;$pos++) {
	$self->{'text'}->{'samples_a'}[$pos] = $self->{'text'}->{'samples_c'}[$pos] 
		= $self->{'text'}->{'samples_g'}[$pos] = $self->{'text'}->{'samples_t'}[$pos] = "0";
    }
	# $self->_dump_traces();
	    # now populate them
    my ($current_base,$place_base_at,$peak_quality,$ramp_counter,$current_ramp,$ramp_position);
    my $sequence_length = $self->{'info'}->{'sequence_length'};
    my $half_ramp = int($self->{'text'}->{'ramp_width'}/2);
    for ($pos = 0; $pos<$sequence_length;$pos++) {
	$current_base = substr($self->{'info'}->{'sequence'},$pos,1);
		# where should the peak for this base be placed? Modeled after a mktrace scf
	$place_base_at = ($pos * $self->{'text'}->{'ramp_width'}) - 
	                 ($pos * $self->{'text'}->{'ramp_overlap'}) - 
		         $half_ramp + $self->{'text'}->{'ramp_width'} - 1;
	push @{$self->{'text'}->{'v3_peak_offsets'}},$place_base_at;
	$peak_quality = $quals[$pos];
	if ($current_base eq "A") {
		$ramp_position = $place_base_at - $half_ramp;
		for ($current_ramp = 0; $current_ramp < $self->{'text'}->{'ramp_width'};  $current_ramp++) {
			$self->{'text'}->{'samples_a'}[$ramp_position+$current_ramp] = $peak_quality * $self->{'text'}->{'ramp'}[$current_ramp];
		}
		push @{$self->{'text'}->{'v2_bases'}},($place_base_at+1,$peak_quality,0,0,0,$current_base,0,0,0);
		push @{$self->{'text'}->{'v3_base_accuracy_a'}},$peak_quality;
		foreach (qw(g c t)) {
			push @{$self->{'text'}->{"v3_base_accuracy_$_"}},0;
		}
	}
	elsif ($current_base eq "C") {
		$ramp_position = $place_base_at - $half_ramp;
		for ($current_ramp = 0; $current_ramp < $self->{'text'}->{'ramp_width'}; $current_ramp++) {
			$self->{'text'}->{'samples_c'}[$ramp_position+$current_ramp] = $peak_quality * $self->{'text'}->{'ramp'}[$current_ramp];
		}
		push @{$self->{'text'}->{'v2_bases'}},($place_base_at+1,0,$peak_quality,0,0,$current_base,0,0,0);
		push @{$self->{'text'}->{'v3_base_accuracy_c'}},$peak_quality;
		foreach (qw(g a t)) {
			push @{$self->{'text'}->{"v3_base_accuracy_$_"}},0;
		}
	} elsif ($current_base eq "G") {
		$ramp_position = $place_base_at - $half_ramp;
		for ($current_ramp = 0; 
			$current_ramp < $self->{'text'}->{'ramp_width'}; 
			$current_ramp++) {
			$self->{'text'}->{'samples_g'}[$ramp_position+$current_ramp] = $peak_quality * $self->{'text'}->{'ramp'}[$current_ramp];
		}
		push @{$self->{'text'}->{'v2_bases'}},($place_base_at+1,0,0,$peak_quality,0,$current_base,0,0,0);
		push @{$self->{'text'}->{"v3_base_accuracy_g"}},$peak_quality;
		foreach (qw(a c t)) {
			push @{$self->{'text'}->{"v3_base_accuracy_$_"}},0;
		}
	}
	elsif( $current_base eq "T" ) { 
		$ramp_position = $place_base_at - $half_ramp;
		for ($current_ramp = 0; $current_ramp < $self->{'text'}->{'ramp_width'}; $current_ramp++) {
			$self->{'text'}->{'samples_t'}[$ramp_position+$current_ramp] = $peak_quality * $self->{'text'}->{'ramp'}[$current_ramp];
		}
		push @{$self->{'text'}->{'v2_bases'}},($place_base_at+1,0,0,0,$peak_quality,$current_base,0,0,0);
		push @{$self->{'text'}->{'v3_base_accuracy_t'}},$peak_quality;
		foreach (qw(g c a)) {
			push @{$self->{'text'}->{"v3_base_accuracy_$_"}},0;
		}
	} elsif ($current_base eq "N") {
	    $ramp_position = $place_base_at - $half_ramp;
	    for ($current_ramp = 0; 
		 $current_ramp < $self->{'text'}->{'ramp_width'}; 
		 $current_ramp++) {
		$self->{'text'}->{'samples_a'}[$ramp_position+$current_ramp] = $peak_quality * $self->{'text'}->{'ramp'}[$current_ramp];
	    }
	    push @{$self->{'text'}->{'v2_bases'}},($place_base_at+1,$peak_quality,
							 $peak_quality,$peak_quality,$peak_quality,
							 $current_base,0,0,0);
		foreach (qw(a c g t)) {
			push @{$self->{'text'}->{"v3_base_accuracy_$_"}},0;
		}
	}
	else {
		# don't print this.
		# print ("The current base ($current_base) is not a base. Hmmm.\n");
	}
    }
	foreach (qw(a c g t)) {
		pop @{$self->{'text'}->{"samples_$_"}};
	}

			# set the samples in the header
		$self->{'header'}->{'samples'} = scalar(@{$self->{'text'}->{'samples_a'}});

			# create the final trace string (this is version dependent)
		$self->_make_trace_string($version);
			# create the binary for v2 bases
		if ($self->{'header'}->{'version'} == 2) {
			my ($packstring,@pack_array,$pos2,$tester,@unpacked);
			for ($pos = 0; $pos<$sequence_length;$pos++) {
				my @pack_array = @{$self->{'text'}->{'v2_bases'}}[$pos*9..$pos*9+8];
				$self->{'binaries'}->{'v2_bases'} .= pack "N C C C C a C3",@pack_array;
			}
			# now create the binary for the traces
			my $trace_pack_length = scalar(@{$self->{'text'}->{'samples_all'}});
	    		$self->{'binaries'}->{'samples_all'} .= pack "n$trace_pack_length",@{$self->{'text'}->{'samples_all'}};
		}
		else {
				# now for the version 3 stuff!
				# delta the trace data
			my @temp;
			foreach (qw(a c g t)) {
				$self->{'text'}->{"t_samples_$_"} = $self->_delta($self->{'text'}->{"samples_$_"},"forward");
					if ($_ eq 'a') {
						@temp = @{$self->{'text'}->{"t_samples_a"}};
						@{$self->{'text'}->{'samples_all'}} = @{$self->{'text'}->{"t_samples_a"}};
					}
					else {
						push @{$self->{'text'}->{'samples_all'}},@{$self->{'text'}->{"t_samples_$_"}};
					}
			}
				# now create the binary for the traces
			my $trace_pack_length = scalar(@{$self->{'text'}->{'samples_all'}});

			$self->{'binaries'}->{'samples_all'} .= pack "n$trace_pack_length",@{$self->{'text'}->{'samples_all'}};

				# peak offsets
			my $length = scalar(@{$self->{'text'}->{'v3_peak_offsets'}});
			$self->{'binaries'}->{'v3_peak_offsets'} = pack "N$length",@{$self->{'text'}->{'v3_peak_offsets'}};
				# base accuracies
			@{$self->{'text'}->{'v3_accuracies_all'}} = @{$self->{'text'}->{"v3_base_accuracy_a"}};
			foreach (qw(c g t)) {
				@{$self->{'text'}->{'v3_accuracies_all'}} = (@{$self->{'text'}->{'v3_accuracies_all'}},@{$self->{'text'}->{"v3_base_accuracy_$_"}});
			}
			$length = scalar(@{$self->{'text'}->{'v3_accuracies_all'}});
			
			$self->{'binaries'}->{'v3_accuracies_all'} = pack "c$length",@{$self->{'text'}->{'v3_accuracies_all'}};
				# called bases
			$length = length($self->{'info'}->{'sequence'});
			my @seq = split(//,$self->{'info'}->{'sequence'});
				# pack the string
			$self->{'binaries'}->{'v3_called_bases'} = $self->{'info'}->{'sequence'};
				# finally, reserved for future use
			$length = $self->{'info'}->{'sequence_length'};
			for (my $counter=0; $counter < $length; $counter++) {
				push @temp,0;
			}
			$self->{'binaries'}->{'v3_reserved'} = pack "N$length",@temp;
		}
}

=head2 _make_trace_string($version)

 Title   : _make_trace_string($version)
 Usage   : $self->_make_trace_string($version)
 Function: Merges trace data for the four bases to produce an scf
	   trace string. _requires_ $version
 Returns : Nothing. Alters $self.
 Args    : $version - a version number. "2" or "3"
 Notes   : 

=cut

sub _make_trace_string {
	my ($self,$version) = @_;
	my @traces;
	my @traces_view;
	my @as = @{$self->{'text'}->{'samples_a'}};
	my @cs = @{$self->{'text'}->{'samples_c'}};
	my @gs = @{$self->{'text'}->{'samples_g'}};
	my @ts = @{$self->{'text'}->{'samples_t'}};
	if ($version == 2) {
	    for (my $curr=0; $curr < scalar(@as); $curr++) {
		$as[$curr] = $DEFAULT_QUALITY unless defined $as[$curr];
		$cs[$curr] = $DEFAULT_QUALITY unless defined $cs[$curr];
		$gs[$curr] = $DEFAULT_QUALITY unless defined $gs[$curr];
		$ts[$curr] = $DEFAULT_QUALITY unless defined $ts[$curr];
		push @traces,($as[$curr],$cs[$curr],$gs[$curr],$ts[$curr]);
	    }
	}
	elsif ($version == 3) {
		@traces = (@as,@cs,@gs,@ts);
	}
	else {
		$self->throw("No idea what version required to make traces here. You gave #$version#  Bailing.");
	}
	my $length = scalar(@traces);
	$self->{'text'}->{'samples_all'} = \@traces;

}	

=head2 _set_binary_comments(\@comments)

 Title   : _set_binary_comments(\@comments)
 Usage   : $self->_set_binary_comments(\@comments);
 Function: Provide a binary string that will be the comments section of
	   the scf file. See the scf specifications for detailed
	   specifications for the comments section of an scf file. Hint:
	   CODE=something\nBODE=something\n\0
 Returns : Nothing. Alters self.
 Args    : A reference to an array containing comments.
 Notes   : None.

=cut

sub _set_binary_comments {
    my ($self,$rcomments) = @_;
    my $comments_string = '';
    my %comments = %$rcomments;
    foreach my $key (sort keys %comments) {
	$comments{$key} ||= '';
	$comments_string .= "$key=$comments{$key}\n";
    }
    $comments_string .= "\n\0";
	$self->{'header'}->{'comments'} = $comments_string;
    my $length = length($comments_string);
    $self->{'binaries'}->{'comments'} = pack "A$length",$comments_string;
	$self->{'header'}->{'comments'} = $comments_string;
}

=head2 _fill_missing_data($swq)

 Title   : _fill_missing_data($swq)
 Usage   : $self->_fill_missing_data($swq);
 Function: If the $swq with quality has no qualities, set all qualities
	   to 0.
	   If the $swq has no sequence, set the sequence to N's.
 Returns : Nothing. Modifies the SeqWithQuality that was passed as an
	   argument.
 Args    : A reference to a Bio::Seq::SeqWithQuality
 Notes   : None.

=cut

#'
sub _fill_missing_data {
    my ($self,$swq) = @_;
    my $qual_obj = $swq->qual_obj();
    my $seq_obj = $swq->seq_obj();
    if ($qual_obj->length() == 0 && $seq_obj->length() != 0) {
	my $fake_qualities = ("$DEFAULT_QUALITY ")x$seq_obj->length();
	$swq->qual($fake_qualities);
    }
    if ($seq_obj->length() == 0 && $qual_obj->length != 0) {
	my $sequence = ("N")x$qual_obj->length();
	$swq->seq($sequence);
    }
}

=head2 _delta(\@trace_data,$direction)

 Title   : _delta(\@trace_data,$direction)
 Usage   : $self->_delta(\@trace_data,$direction);
 Function: 
 Returns : A reference to an array containing modified trace values.
 Args    : A reference to an array containing trace data and a string
	   indicating the direction of conversion. ("forward" or
	   "backward").
 Notes   : This code is taken from the specification for SCF3.2.
	   http://www.mrc-lmb.cam.ac.uk/pubseq/manual/formats_unix_4.html

=cut


sub _delta {
	my ($self,$rsamples,$direction) = @_;
	my @samples = @$rsamples;
		# /* If job == DELTA_IT:
		# *  change a series of sample points to a series of delta delta values:
		# *  ie change them in two steps:
		# *  first: delta = current_value - previous_value
		# *  then: delta_delta = delta - previous_delta
		# * else
		# *  do the reverse
		# */
		# int i;
		# uint_2 p_delta, p_sample;

	my ($i,$num_samples,$p_delta,$p_sample,@samples_converted);

		# c-programmers are funny people with their single-letter variables

	if ( $direction eq "forward" ) {
		$p_delta  = 0;
		for ($i=0; $i < scalar(@samples); $i++) {
			$p_sample = $samples[$i];
			$samples[$i] = $samples[$i] - $p_delta;
			$p_delta  = $p_sample;
		}
		$p_delta  = 0;
		for ($i=0; $i < scalar(@samples); $i++) {
			$p_sample = $samples[$i];
			$samples[$i] = $samples[$i] - $p_delta;
			$p_delta  = $p_sample;
		}
	}
	elsif ($direction eq "backward") {
		$p_sample = 0;
		for ($i=0; $i < scalar(@samples); $i++) {
			$samples[$i] = $samples[$i] + $p_sample;
			$p_sample = $samples[$i];
		}
		$p_sample = 0;
		for ($i=0; $i < scalar(@samples); $i++) {
			$samples[$i] = $samples[$i] + $p_sample;
			$p_sample = $samples[$i];
		}
	}
	else {
		$self->warn("Bad direction. Use \"forward\" or \"backward\".");
	}
	return \@samples;
}

=head2 _unpack_magik($buffer)

 Title   : _unpack_magik($buffer)
 Usage   : $self->_unpack_magik($buffer)
 Function: What unpack specification should be used? Try them all.
 Returns : Nothing.
 Args    : A buffer containing arbitrary binary data.
 Notes   : Eliminate the ambiguity and the guesswork. Used in the
	   adaptation of _delta(), mostly.

=cut

sub _unpack_magik {
	my ($self,$buffer) = @_;
	my $length = length($buffer);
	my (@read,$counter);
	foreach (qw(c C s S i I l L n N v V)) {
		@read = unpack "$_$length", $buffer;
		print ("----- Unpacked with $_\n");
		for ($counter=0; $counter < 20; $counter++) {
			print("$read[$counter]\n");
		}
	}
}

=head2 read_from_buffer($filehandle,$buffer,$length)

 Title   : read_from_buffer($filehandle,$buffer,$length)
 Usage   : $self->read_from_buffer($filehandle,$buffer,$length);
 Function: Read from the buffer.
 Returns : $buffer, containing a read of $length
 Args    : a filehandle, a buffer, and a read length
 Notes   : I just got tired of typing
	   "unless (length($buffer) == $length)" so I put it here.

=cut

sub read_from_buffer {
	my ($self,$fh,$buffer,$length) = @_;
	read $fh, $buffer, $length;
	unless (length($buffer) == $length) {
		$self->warn("The read was incomplete! Trying harder.");
		my $missing_length = $length - length($buffer);
		my $buffer2;
		read $fh,$buffer2,$missing_length;
		$buffer .= $buffer2;
		if (length($buffer) != $length) {	
			$self->throw("Unexpected end of file while reading from SCF file. I should have read $length but instead got ".length($buffer)."! Current file position is ".tell($fh).".");
		}
	}
	
	return $buffer;
}

=head2 _dump_keys()

 Title   : _dump_keys()
 Usage   : &_dump_keys($a_reference_to_some_hash)
 Function: Dump out the keys in a hash.
 Returns : Nothing.
 Args    : A reference to a hash.
 Notes   : A debugging method.

=cut

sub _dump_keys {
	my $rhash = shift;
	if ($rhash !~ /HASH/) {
		print("_dump_keys: that was not a hash.\nIt was #$rhash# which was this reference:".ref($rhash)."\n");
		return;
	}
	print("_dump_keys: The keys for $rhash are:\n");
	foreach (sort keys %$rhash) {
		print("$_\n");
	}
}

=head2 _dump_base_accuracies()

 Title   : _dump_base_accuracies()
 Usage   : $self->_dump_base_accuracies();
 Function: Dump out the v3 base accuracies in an easy to read format.
 Returns : Nothing.
 Args    : None.
 Notes   : A debugging method.

=cut

sub _dump_base_accuracies {
	my $self = shift;
	print("Dumping base accuracies! for v3\n");
	print("There are this many elements in a,c,g,t:\n");
	print(scalar(@{$self->{'text'}->{'v3_base_accuracy_a'}}).",".scalar(@{$self->{'text'}->{'v3_base_accuracy_c'}}).",".scalar(@{$self->{'text'}->{'v3_base_accuracy_g'}}).",".scalar(@{$self->{'text'}->{'v3_base_accuracy_t'}})."\n");
	my $number_traces = scalar(@{$self->{'text'}->{'v3_base_accuracy_a'}});
	for (my $counter=0; $counter < $number_traces; $counter++ ) {
		print("$counter\t");
		print $self->{'text'}->{'v3_base_accuracy_a'}->[$counter]."\t";
		print $self->{'text'}->{'v3_base_accuracy_c'}->[$counter]."\t";
		print $self->{'text'}->{'v3_base_accuracy_g'}->[$counter]."\t";
		print $self->{'text'}->{'v3_base_accuracy_t'}->[$counter]."\t";
		print("\n");
	}
}

=head2 _dump_peak_indices_incoming()

 Title   : _dump_peak_indices_incoming()
 Usage   : $self->_dump_peak_indices_incoming();
 Function: Dump out the v3 peak indices in an easy to read format.
 Returns : Nothing.
 Args    : None.
 Notes   : A debugging method.

=cut

sub _dump_peak_indices_incoming {
	my $self = shift;
	print("Dump peak indices incoming!\n");
	my $length = $self->{'bases'};
	print("The length is $length\n");
	for (my $count=0; $count < $length; $count++) {
		print("$count\t$self->{parsed}->{peak_indices}->[$count]\n");
	}
}

=head2 _dump_base_accuracies_incoming()

 Title   : _dump_base_accuracies_incoming()
 Usage   : $self->_dump_base_accuracies_incoming();
 Function: Dump out the v3 base accuracies in an easy to read format.
 Returns : Nothing.
 Args    : None.
 Notes   : A debugging method.

=cut

sub _dump_base_accuracies_incoming {
	my $self = shift;
	print("Dumping base accuracies! for v3\n");
		# print("There are this many elements in a,c,g,t:\n");
		# print(scalar(@{$self->{'parsed'}->{'v3_base_accuracy_a'}}).",".scalar(@{$self->{'text'}->{'v3_base_accuracy_c'}}).",".scalar(@{$self->{'text'}->{'v3_base_accuracy_g'}}).",".scalar(@{$self->{'text'}->{'v3_base_accuracy_t'}})."\n");
	my $number_traces = $self->{'bases'};
	for (my $counter=0; $counter < $number_traces; $counter++ ) {
		print("$counter\t");
		foreach (qw(A T G C)) {
			print $self->{'parsed'}->{'base_accuracies'}->{$_}->[$counter]."\t";
		}
		print("\n");
	}
}



1;
__END__


