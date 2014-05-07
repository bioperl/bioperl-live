#
# BioPerl module for Bio::AlignIO::nexus
#
# Copyright Heikki Lehvaslaiho
#

=head1 NAME

Bio::AlignIO::nexus - NEXUS format sequence input/output stream

=head1 SYNOPSIS

Do not use this module directly.  Use it via the L<Bio::AlignIO> class.

    use Bio::AlignIO;

    my $in = Bio::AlignIO->new(-format => 'nexus',
                              -file   => 'aln.nexus');
    while( my $aln = $in->next_aln ) {
        # do something with the alignment
    }

=head1 DESCRIPTION

This object can transform L<Bio::Align::AlignI> objects to and from NEXUS
data blocks. See method documentation for supported NEXUS features.

=head1 ACKNOWLEDGEMENTS

Will Fisher has written an excellent standalone NEXUS format parser in
Perl, readnexus. A number of tricks were adapted from it.

=head1 FEEDBACK

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

=head1 AUTHORS - Heikki Lehvaslaiho

Email: heikki-at-bioperl-dot-org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::nexus;
use vars qw(%valid_type);
use strict;
no strict "refs";


use base qw(Bio::AlignIO);

BEGIN {
    %valid_type = map {$_, 1} qw( dna rna protein standard );
    # standard throws error: inherited from Bio::PrimarySeq
}

=head2 new

 Title   : new
 Usage   : $alignio = Bio::AlignIO->new(-format => 'nexus', -file => 'filename');
 Function: returns a new Bio::AlignIO object to handle clustalw files
 Returns : Bio::AlignIO::clustalw object
 Args    : -verbose => verbosity setting (-1,0,1,2)
           -file    => name of file to read in or with ">" - writeout
           -fh      => alternative to -file param - provide a filehandle 
                       to read from/write to 
           -format  => type of Alignment Format to process or produce

           Customization of nexus flavor output

           -show_symbols => print the symbols="ATGC" in the data definition
                            (MrBayes does not like this)
                            boolean [default is 1] 
           -show_endblock => print an 'endblock;' at the end of the data
                            (MyBayes does not like this)
                            boolean [default is 1] 

=cut

sub _initialize {
    my ($self, @args) = @_;
    $self->SUPER::_initialize(@args);
    my ($show_symbols, $endblock) = 
	$self->_rearrange([qw(SHOW_SYMBOLS SHOW_ENDBLOCK)], @args);
    my @names = qw(symbols endblock);
    for my $v ( $show_symbols, $endblock ) {
	$v = 1 unless defined $v; # default value is 1
	my $n = shift @names;
	$self->flag($n, $v);
    }
}


=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: Returns the next alignment in the stream.

           Supports the following NEXUS format features:
           - The file has to start with '#NEXUS'
           - Reads in the name of the alignment from a comment
             (anything after 'TITLE: ') .
           - Sequence names can be given in a taxa block, too.
           - If matchchar notation is used, converts
             them back to sequence characters.
           - Does character conversions specified in the
             NEXUS equate command.
           - Sequence names of type 'Homo sapiens' and
             Homo_sapiens are treated identically.

 Returns : L<Bio::Align::AlignI> object
 Args    :


=cut


sub next_aln {
    my $self = shift;
    my $entry;
    my ($aln_name, $seqcount, $residuecount, %hash, $alphabet,
	$match, $gap, $missing, $equate, $interleave,
	$name,$str,@names,$seqname,$start,$end,$count,$seq);
    local $Bio::LocatableSeq::OTHER_SYMBOLS = '\*\?\.';
    local $Bio::LocatableSeq::GAP_SYMBOLS = '\-';
    my $aln =  Bio::SimpleAlign->new(-source => 'nexus');

    # file starts with '#NEXUS' but we allow white space only lines before it
    $entry = $self->_readline;
    $entry = $self->_readline while defined $entry && $entry =~ /^\s+$/;

    return unless $entry;
    $self->throw("Not a valid interleaved NEXUS file! [#NEXUS] not starting the file\n$entry")
	unless ($entry && $entry =~ /^#NEXUS/i);

    # skip anything before either the taxa or data block
    # but read in the optional title in a comment
    while (defined($entry = $self->_readline)) {
	local ($_) = $entry;
	/\[TITLE. *([^\]]+)]\s+/i and $aln_name = $1;
	last if /^begin +data/i || /^begin +taxa/i;
    }
    $aln_name =~ s/\s/_/g and $aln->id($aln_name) if $aln_name;

    # data and taxa blocks
    my $incomment;
    while (defined ($entry = $self->_readline)) {
	local ($_) =  $entry;
	next if s/\[[^\]]+\]//g; # remove comments
	if( s/\[[^\]]+$// ) {
	    $incomment = 1;
		 # skip line if it is now empty or contains only whitespace
	    next if /^\s*$/;
	} elsif($incomment) {
	    if( s/^[^\]]*\]// ) {
			 $incomment = 0;
	    } else {
			 next;
	    }
	} elsif( /taxlabels/i ) {
	    # doesn't deal with taxlabels adequately and can mess things up!
	    # @names = $self->_read_taxlabels;
	} else {

	    /ntax\s*=\s*(\d+)/i        and $seqcount = $1;
	    /nchar\s*=\s*(\d+)/i       and $residuecount = $1;
	    /matchchar\s*=\s*(.)/i     and $match = $1;
	    /gap\s*=\s*(.)/i           and $gap = $1;
	    /missing\s*=\s*(.)/i       and $missing = $1;
	    /equate\s*=\s*\"([^\"]+)/i and $equate = $1;  # "e.g. equate="T=C G=A";
	    /datatype\s*=\s*(\w+)/i    and $alphabet = lc $1;
	    /interleave/i              and $interleave = 1 ;
	    last if /matrix/io;
	}
    }
    $self->throw("Not a valid NEXUS sequence file. Datatype not specified.")
	unless $alphabet;
    $self->throw("Not a valid NEXUS sequence file. Datatype should not be [$alphabet]")
	unless $valid_type{$alphabet};
    $self->throw("\"$gap\" is not a valid gap character. For compatability, gap char can not be one of: ()[]{}/\,;:=*'`\"<>^")
    	if $gap && $gap =~ /[\(\)\[\]\{\}\/\\\,\;\:\=\*\'\`\<\>\^]/;
    $self->throw("\"$missing\" is not a valid missing character. For compatability, missing char can not be one of: ()[]{}/\,;:=*'`\"<>^")
    	if $missing && $missing =~ /[\(\)\[\]\{\}\/\\\,\;\:\=\*\'\`\<\>\^]/;

    $aln->gap_char($gap);
    $aln->missing_char($missing);

    #
    # if data is not right after the matrix line
    #  read the empty lines out
    #
    while ($entry = $self->_readline) {
	unless ($entry =~ /^\s+$/) {
	    $self->_pushback($entry);
	    last;
	}
    }

    #
    # matrix command
    #
    # first alignment section
    if (@names == 0) {		# taxa block did not exist
	while ($entry = $self->_readline) {
		local ($_) =  $entry;
		if( s/\[[^\]]+\]//g ) { #]  remove comments
			next if /^\s*$/; 
			# skip line if it is now empty or contains only whitespace
		}
		if ($interleave && defined$count && ($count <= $seqcount)) {
			/^\s+$/ and last;
		} else {
			/^\s+$/ and next;
		}
		/^\s*;/ and last;	# stop if colon at end of matrix is on it's own line
		#/^\s*;\s*$/ and last;
		if ( /^\s*([\"\'](.+?)[\"\']|(\S+))\s+(.*)\s*$/ ) {	
			# get single and double quoted names, or all the first 
         # nonwhite word as the name, and remained is seq
			#if (/^\s*('([^']*?)'|([^']\S*))\s+(.*)$/) { #'
			$name = ($2 || $3);
			if  ($4) {
				# seq is on same line as name
				# this is the usual NEXUS format
				$str = $4;
			} else {
				# otherwise get seq from following lines. No comments allowed
				# a less common matrix format, usually used for very long seqs
				$str='';
				while (local ($_) = $self->_readline) {
					my $str_tmp = $_;
					$str_tmp =~ s/[\s;]//g;
					$str .= $str_tmp;
					last if length$str == $residuecount;
				}
			}
			$name =~ s/ /_/g;
			push @names, $name;

			$str =~ s/[\s;]//g;
			$count =  @names;
			$hash{$count} = $str;
		}
		$self->throw("Not a valid interleaved NEXUS file! seqcount [$count] > predeclared [$seqcount] in the first section") if $count > $seqcount;
		/;/ and last;	# stop if colon at end of matrix is on the same line as the last seq
	}
}

    # interleaved sections
    $count = 0;
    if ( $interleave ) {	# only read next section if file is interleaved
	while( $entry = $self->_readline) {
	    local ($_) =  $entry;
	    if( s/\[[^\]]+\]//g ) { #]  remove comments
		next if /^\s*$/; # skip line if it is now empty or contains only whitespace
	    }
	    /^\s*;/ and last;		# stop if colon at end of matrix is on it's own line
	    $count = 0, next if $entry =~ /^\s*$/;
	    if (/^\s*(\'([^\']*?)\'|([^\']\S*))\s+(.*)$/) { 
		$str = $4;
		$str =~ s/[\s;]//g;
		$count++;
		$hash{$count} .= $str;
	    };
	    $self->throw("Not a valid interleaved NEXUS file!
    		seqcount [$count] > predeclared [$seqcount] ") if $count > $seqcount;
	    /;/ and last;	# stop if colon at end of matrix is on the same line as the last seq
	}
    }
    
    return if @names < 1;
    
    # sequence creation
    $count = 0;
    foreach $name ( @names ) {
	$count++;
	if( $name =~ /(\S+)\/(\d+)-(\d+)/ ) {
	    ($seqname,$start,$end) = ($1,$2,$3);
	} else {
	    ($seqname,$start,$str) = ($name,1,$hash{$count});
	    $str =~ s/[$Bio::LocatableSeq::GAP_SYMBOLS]//g;
	    $end = length($str);
	}
	
	# consistency test
	$self->throw("Length of sequence [$seqname] is not [$residuecount]; got".CORE::length($hash{$count}))
	    unless CORE::length($hash{$count}) == $residuecount;
	
	$seq = Bio::LocatableSeq->new('-seq'        => $hash{$count},
				      '-display_id' => $seqname,
				      '-start'      => $start,
				      '-end'        => $end,
				     '-alphabet'    => $alphabet
				      );
	$aln->add_seq($seq);
    }

    # if matchchar is used
    $aln->unmatch($match) if $match;

    # if equate ( e.g. equate="T=C G=A") is used
    if ($equate) {
	$aln->map_chars($1, $2) while $equate =~ /(\S)=(\S)/g;
    }

    while  (defined $entry &&
	    $entry !~ /endblock/i) {
        $entry = $self->_readline;
    }

    return $aln if $aln->num_sequences;
	return;
}

sub _read_taxlabels {
    my ($self) = @_;
    my ($name, @names);
    while (my $entry = $self->_readline) {
	last if $entry =~ m/^\s*(END)?;/i;
	if( $entry =~ m/\s*(\S+)\s+/ ) {
	    ($name) = ($1);
	    $name =~ s/\[[^\[]+\]//g;
	    $name =~ s/\W/_/g;
	    push @names, $name;
	}
    }
    return @names;
}

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: Writes the $aln object into the stream in interleaved NEXUS
           format. Everything is written into a data block.
           SimpleAlign methods match_char, missing_char and gap_char must be set
           if you want to see them in the output.
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object

=cut

sub write_aln {
    my ($self,@aln) = @_;
    my $count = 0;
    my $wrapped = 0;
    my $maxname;
    my ($length,$date,$name,$seq,$miss,$pad,%hash,@arr,$tempcount,$index );
    my ($match, $missing, $gap,$symbols) = ('', '', '','');

    foreach my $aln (@aln) {
	if( ! $aln || ! $aln->isa('Bio::Align::AlignI')  ) {
	    $self->warn("Must provide a Bio::Align::AlignI object when calling write_aln");
	    next;
	}
	$self->throw("All sequences in the alignment must be the same length")
	    unless $aln->is_flush($self->verbose);

	$length  = $aln->length();

	$self->_print (sprintf("#NEXUS\n[TITLE: %s]\n\nbegin data;\ndimensions ntax=%s nchar=%s;\n",
			       $aln->id, $aln->num_sequences, $length));
	$match = "match=". $aln->match_char if $aln->match_char;
	$missing = "missing=". $aln->missing_char if $aln->missing_char;
	$gap = "gap=". $aln->gap_char if $aln->gap_char;

	$symbols = 'symbols="'.join('',$aln->symbol_chars). '"' 
	    if( $self->flag('symbols') && $aln->symbol_chars);
	$self->_print 
	    (sprintf("format interleave datatype=%s %s %s %s %s;\n\nmatrix\n",
		     $aln->get_seq_by_pos(1)->alphabet, $match, 
		     $missing, $gap, $symbols));

                     # account for single quotes round names
	my $indent = $aln->maxdisplayname_length+2;

	$aln->set_displayname_flat();
	foreach $seq ( $aln->each_seq() ) {
	    my $nmid = $aln->displayname($seq->get_nse());
	    if( $nmid =~ /[^\w\d\.]/ ) {
              # put name in single quotes incase it contains any of
              # the following chars: ()[]{}/\,;:=*'"`+-<> that are not
              # allowed in PAUP* and possible other software

		$name = sprintf("%-${indent}s", "\'" . $nmid . "\'");
	    } else { 
		$name = sprintf("%-${indent}s", $nmid);
	    }
	    $hash{$name} = $seq->seq;
	    push(@arr,$name);
	}

	while( $count < $length ) {
	    # there is another block to go!
	    foreach $name ( @arr ) {
		my $dispname = $name;
#		 $dispname = '' if $wrapped;
		$self->_print (sprintf("%${indent}s  ",$dispname));
		$tempcount = $count;
		$index = 0;
		while( ($tempcount + 10 < $length) && ($index < 5)  ) {
		    $self->_print (sprintf("%s ",substr($hash{$name},$tempcount,10)));
		    $tempcount += 10;
		    $index++;
		}
		# last
		if( $index < 5) {
		    # space to print!
		    $self->_print (sprintf("%s ",substr($hash{$name},$tempcount)));
		    $tempcount += 10;
		}
		$self->_print ("\n");
	    }
	    $self->_print ("\n\n");
	    $count = $tempcount;
	    $wrapped = 1;
	}
	if( $self->flag('endblock') ) {
	    $self->_print (";\n\nendblock;\n");
	} else { 
	    $self->_print (";\n\nend;\n");
	}
    }
    $self->flush if $self->_flush_on_write && defined $self->_fh;
    return 1;
}

=head2 flag

 Title   : flag
 Usage   : $obj->flag($name,$value)
 Function: Get/Set a flag value
 Returns : value of flag (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub flag{
    my ($self,$name,$val) = @_;
    return $self->{'flag'}->{$name} = $val if defined $val;
    return $self->{'flag'}->{$name};
}

1;
