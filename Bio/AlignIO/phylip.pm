#
# BioPerl module for Bio::AlignIO::phylip
#
# Copyright Heikki Lehvaslaiho
#

=head1 NAME

Bio::AlignIO::phylip - PHYLIP format sequence input/output stream

=head1 SYNOPSIS

# Do not use this module directly.  Use it via the Bio::AlignIO class.

    use Bio::AlignIO;
    use Bio::SimpleAlign;
    #you can set the name length to something other than the default 10
    #if you use a version of phylip (hacked) that accepts ids > 10
    my $phylipstream = Bio::AlignIO->new(-format  => 'phylip',
                                        -fh      => \*STDOUT,
                                        -idlength=>30);
    # convert data from one format to another
    my $gcgstream     =  Bio::AlignIO->new(-format => 'msf',
                                          -file   => 't/data/cysprot1a.msf');

    while( my $aln = $gcgstream->next_aln ) {
        $phylipstream->write_aln($aln);
    }

    # do it again with phylip sequential format format
    $phylipstream->interleaved(0);
    # can also initialize the object like this
    $phylipstream = Bio::AlignIO->new(-interleaved => 0,
                                     -format => 'phylip',
                                     -fh   => \*STDOUT,
                                     -idlength=>10);
    $gcgstream     =  Bio::AlignIO->new(-format => 'msf',
                                       -file   => 't/data/cysprot1a.msf');

    while( my $aln = $gcgstream->next_aln ) {
        $phylipstream->write_aln($aln);
    }

=head1 DESCRIPTION

This object can transform Bio::SimpleAlign objects to and from PHYLIP
format. By default it works with the interleaved format. By specifying
the flag -interleaved =E<gt> 0 in the initialization the module can
read or write data in sequential format.

Long IDs up to 50 characters are supported by flag -longid =E<gt>
1. ID strings can be surrounded by single quoted. They are mandatory
only if the IDs contain spaces.

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
the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHORS - Heikki Lehvaslaiho and Jason Stajich

Email: heikki at ebi.ac.uk
Email: jason at bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::phylip;
use vars qw($DEFAULTIDLENGTH $DEFAULTLINELEN $DEFAULTTAGLEN);
use strict;

use Bio::SimpleAlign;
use POSIX; # for the rounding call

use base qw(Bio::AlignIO);

BEGIN {
    $DEFAULTIDLENGTH = 10;
    $DEFAULTLINELEN = 60;
    $DEFAULTTAGLEN = 10;
}

=head2 new

 Title   : new
 Usage   : my $alignio = Bio::AlignIO->new(-format => 'phylip'
					  -file   => '>file',
					  -idlength => 10,
					  -idlinebreak => 1);
 Function: Initialize a new L<Bio::AlignIO::phylip> reader or writer
 Returns : L<Bio::AlignIO> object
 Args    : [specific for writing of phylip format files]
           -idlength => integer - length of the id (will pad w/
						    spaces if needed)
           -interleaved => boolean - whether interleaved
                                     or sequential format required
           -line_length  => integer of how long a sequence lines should be
           -idlinebreak => insert a line break after the sequence id
                           so that sequence starts on the next line
           -flag_SI => whether or not write a "S" or "I" just after
                       the num.seq. and line len., in the first line
           -tag_length => integer of how long the tags have to be in
                         each line between the space separator. set it
                         to 0 to have 1 tag only.
           -wrap_sequential => boolean for whether or not sequential
                                   format should be broken up or a single line
                                   default is false (single line)
           -longid => boolean for allowing arbitrary long IDs (default is false)

=cut

sub _initialize {
  my($self,@args) = @_;
  $self->SUPER::_initialize(@args);

  my ($interleave,$linelen,$idlinebreak,
      $idlength, $flag_SI, $tag_length,$ws, $longid) =
          $self->_rearrange([qw(INTERLEAVED
                                LINE_LENGTH
                                IDLINEBREAK
                                IDLENGTH
                                FLAG_SI
                                TAG_LENGTH
				WRAP_SEQUENTIAL
                                LONGID)],@args);
  $self->interleaved($interleave ? 1 : 0) if defined $interleave;
  $self->idlength($idlength || $DEFAULTIDLENGTH);
  $self->id_linebreak(1) if( $idlinebreak );
  $self->line_length($linelen) if defined $linelen && $linelen > 0;
  $self->flag_SI(1) if ( $flag_SI );
  $self->tag_length($tag_length) if ( $tag_length || $DEFAULTTAGLEN );
  $self->wrap_sequential($ws ? 1 : 0);
  $self->longid($longid ? 1 : 0);
  1;
}

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
           Throws an exception if trying to read in PHYLIP
           sequential format.
 Returns : L<Bio::SimpleAlign> object
 Args    :

=cut

sub next_aln {
    my $self = shift;
    my $entry;
    my ($seqcount, $residuecount, %hash, $name,$str,
	@names,$seqname,$start,$end,$count,$seq);

    my $aln =  Bio::SimpleAlign->new(-source => 'phylip');

    # First, parse up through the header.
    # If we see a non-blank line that isn't the seqcount and residuecount line
    # then bail out of next_aln (return)
    while ($entry = $self->_readline) {
        if ($entry =~ /^\s?$/) {
            next;
        } elsif ($entry =~ /\s*(\d+)\s+(\d+)/) {
            ($seqcount, $residuecount) = ($1, $2);
            last;
        } else {
        	$self->warn ("Failed to parse PHYLIP: Did not see a sequence count and residue count.");
            return;
        }
    }
    return unless $seqcount and $residuecount;
    
    # First alignment section. We expect to see a name and (part of) a sequence.
    my $idlen = $self->idlength;
    $count = 0;

    while ($entry = $self->_readline) {
        if ($entry =~ /^\s?$/) { # eat the newlines
            next;
        }

        # Names can be in a few different formats:
        # 1. they can be traditional phylip: 10 chars long, period. If this is the case, that name can have spaces.
        # 2. they can be hacked with a long ID, as passed in with the flag -longid.
        # 3. if there is a long ID, the name can have spaces as long as it is wrapped in single quotes.
        if ($self->longid()) { # 2 or 3
            if ($entry =~ /^'(.+)'\s+(.+)$/) { # 3. name has single quotes.
                $name = $1;
                $str = $2;
            } else {    # 2. name does not have single quotes, so should not have spaces.
                # therefore, the first part of the line is the name and the rest is the seq.
                # make sure that the line does not lead with extra spaces.
                $entry =~ s/^\s+//;
                ($name, $str) = split (/\s+/,$entry, 2);
            }
        } else { # 1. traditional phylip.
            $entry =~ /^(.{10})\s+(.+)$/;
            $name = $1;
            $str = $2;
            $name =~ s/\s+$//; # eat any trailing spaces
            $name =~ s/\s+/_/g;
        }
        push @names, $name;
        #clean sequence of spaces:
        $str =~ s/\s+//g;

        # are we sequential? If so, we should keep adding to the sequence until we've got all the residues.
        if (($self->interleaved) == 0) {
            while (length($str) < $residuecount) {
                $entry = $self->_readline;
                $str .= $entry;
                $str =~ s/\s+//g;
                if ($entry =~ /^\s*$/) { # we ran into a newline before we got a complete sequence: bail!
                    $self->warn("Failed to parse PHYLIP: Sequence $name was shorter than expected: " . length($str) . " instead of $residuecount.");
                    last;
                }
            }
        }
        $hash{$count} = $str;

        $count++;
        # if we've read as many seqs as we're supposed to, move on.
        if ($count == $seqcount) {
            last;
        }
    }

    # if we are interleaved, we're going to keep seeing chunks of sequence until we get all of it.
    if ($self->interleaved) {
        while (length($hash{$seqcount-1}) < $residuecount) {
            $count = 0;
            while ($entry = $self->_readline) {
                if ($entry =~ /^\s*$/) { # eat newlines
                    if ($count != 0) { # there was a newline at an unexpected place!
                        $self->warn("Failed to parse PHYLIP: Interleaved file is missing a segment: saw $count, expected $seqcount.");
                        return;
                    }
                    next;
                } else { # start taking in chunks
                    $entry =~ s/\s//g;
                    $hash{$count} .= $entry;
                    $count++;
                }
                if ($count >= $seqcount) { # we've read all of the sequences for this chunk, so move on.
                    last;
                }
            }
        }
    }
    if ((scalar @names) != $seqcount) {
        $self->warn("Failed to parse PHYLIP: Did not see the correct number of seqs: saw " . scalar(@names) . ", expected $seqcount.");
        return;
    }
    for ($count=0; $count<$seqcount; $count++) {
        $str = $hash{$count};
        my $seqname = $names[$count];
        if (length($str) != $residuecount) {
            $self->warn("Failed to parse PHYLIP: Sequence $seqname was the wrong length: " . length($str) . " instead of $residuecount.");
        }
        $seq = Bio::LocatableSeq->new('-seq'  => $hash{$count},
                          '-display_id'    => $seqname);
        $aln->add_seq($seq);
    }
    return $aln;
}

=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in phylip format
 Returns : 1 for success and 0 for error
 Args    : L<Bio::Align::AlignI> object

=cut

sub write_aln {
    my ($self,@aln) = @_;
    my $count = 0;
    my $wrapped = 0;
    my $maxname;
    my $width = $self->line_length();
    my ($length,$date,$name,$seq,$miss,$pad,
	%hash,@arr,$tempcount,$index,$idlength,$flag_SI,$line_length, $tag_length);

    foreach my $aln (@aln) {
	if( ! $aln || ! $aln->isa('Bio::Align::AlignI')  ) {
	    $self->warn("Must provide a Bio::Align::AlignI object when calling write_aln");
	    next;
	}
	$self->throw("All sequences in the alignment must be the same length")
	    unless $aln->is_flush(1) ;

        $flag_SI = $self->flag_SI();
	$aln->set_displayname_flat(); # plain
	$length  = $aln->length();
        if ($flag_SI) {
            if ($self->interleaved() ) {
                $self->_print (sprintf(" %s %s I\n", $aln->num_sequences, $aln->length));
            } else {
                $self->_print (sprintf(" %s %s S\n", $aln->num_sequences, $aln->length));
            }
        } else {
            $self->_print (sprintf(" %s %s\n", $aln->num_sequences, $aln->length));
        }

	$idlength = $self->idlength();
	$line_length = $self->line_length();
	$tag_length = $self->tag_length();
	foreach $seq ( $aln->each_seq() ) {
	    $name = $aln->displayname($seq->get_nse);
	    if ($self->longid) {
		$self->warn("The length of the name is over 50 chars long [$name]")
		    if length($name) > 50;
		$name = "'$name'  "
	    } else {
		$name = substr($name, 0, $idlength) if length($name) > $idlength;
		$name = sprintf("%-".$idlength."s",$name);
		if( $self->interleaved() ) {
		    $name .= '   ' ;
		} elsif( $self->id_linebreak) {
		    $name .= "\n";
		}
	    }
	    #phylip needs dashes not dots
	    my $seq = $seq->seq();
	    $seq =~ s/\./-/g;
	    $hash{$name} = $seq;
	    push(@arr,$name);
	}

	if( $self->interleaved() ) {
            my $numtags;
            if ($tag_length <= $line_length) {
                $numtags = floor($line_length/$tag_length);
                $line_length = $tag_length*$numtags;
            } else {
                $numtags = 1;
            }
	    while( $count < $length ) {

		# there is another block to go!
		foreach $name ( @arr ) {
		    my $dispname = $name;
		    $dispname = '' if $wrapped;
		    $self->_print (sprintf("%".($idlength+3)."s",$dispname));
		    $tempcount = $count;
                    $index = 0;
                    $self->debug("residue count: $count\n") if ($count%100000 == 0);
		    while( ($tempcount + $tag_length < $length) &&
			   ($index < $numtags)  ) {
			$self->_print (sprintf("%s ",substr($hash{$name},
							    $tempcount,
							    $tag_length)));
			$tempcount += $tag_length;
			$index++;
		    }
		    # last
		    if( $index < $numtags) {
			# space to print!
			$self->_print (sprintf("%s",substr($hash{$name},
							    $tempcount)));
			$tempcount += $tag_length;
		    }
		    $self->_print ("\n");
		}
		$self->_print ("\n");
		$count = $tempcount;
		$wrapped = 1;
	    }
	} else {
	    foreach $name ( @arr ) {
		my $dispname = $name;
		my $line = sprintf("%s%s\n",$dispname,$hash{$name});
		if( $self->wrap_sequential ) {
		    $line =~ s/(.{1,$width})/$1\n/g;
		}
		$self->_print ($line);
	    }
	}
    }
    $self->flush if $self->_flush_on_write && defined $self->_fh;
    return 1;
}

=head2 interleaved

 Title   : interleaved
 Usage   : my $interleaved = $obj->interleaved
 Function: Get/Set Interleaved status
 Returns : boolean
 Args    : boolean


=cut

sub interleaved {
   my ($self,$value) = @_;
   if( defined $value ) {
       if ($value) {$self->{'_interleaved'} = 1 }
       else {$self->{'_interleaved'} = 0 }
   }
   return 1 unless defined $self->{'_interleaved'};
   return $self->{'_interleaved'};
}

=head2 flag_SI

 Title   : flag_SI
 Usage   : my $flag = $obj->flag_SI
 Function: Get/Set if the Sequential/Interleaved flag has to be shown
           after the number of sequences and sequence length
 Example :
 Returns : boolean
 Args    : boolean


=cut

sub flag_SI{
   my ($self,$value) = @_;
   my $previous = $self->{'_flag_SI'};
   if( defined $value ) {
       $self->{'_flag_SI'} = $value;
   }
   return $previous;
}

=head2 idlength

 Title   : idlength
 Usage   : my $idlength = $obj->idlength
 Function: Get/Set value of id length
 Returns : string
 Args    : string


=cut

sub idlength {
	my($self,$value) = @_;
	if (defined $value){
	   $self->{'_idlength'} = $value;
	}
	return $self->{'_idlength'};
}

=head2 line_length

 Title   : line_length
 Usage   : $obj->line_length($newval)
 Function:
 Returns : value of line_length
 Args    : newvalue (optional)


=cut

sub line_length{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_line_length'} = $value;
    }
    return $self->{'_line_length'} || $DEFAULTLINELEN;

}

=head2 tag_length

 Title   : tag_length
 Usage   : $obj->tag_length($newval)
 Function:
 Example : my $tag_length = $obj->tag_length
 Returns : value of the length for each space-separated tag in a line
 Args    : newvalue (optional) - set to zero to have one tag per line


=cut

sub tag_length{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_tag_length'} = $value;
    }
    return $self->{'_tag_length'} || $DEFAULTTAGLEN;
}


=head2 id_linebreak

 Title   : id_linebreak
 Usage   : $obj->id_linebreak($newval)
 Function:
 Returns : value of id_linebreak
 Args    : newvalue (optional)


=cut

sub id_linebreak{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_id_linebreak'} = $value;
    }
    return $self->{'_id_linebreak'} || 0;
}


=head2 wrap_sequential

 Title   : wrap_sequential
 Usage   : $obj->wrap_sequential($newval)
 Function:
 Returns : value of wrap_sequential
 Args    : newvalue (optional)


=cut

sub wrap_sequential{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_wrap_sequential'} = $value;
    }
    return $self->{'_wrap_sequential'} || 0;
}

=head2 longid

 Title   : longid
 Usage   : $obj->longid($newval)
 Function:
 Returns : value of longid
 Args    : newvalue (optional)


=cut

sub longid{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_longid'} = $value;
    }
    return $self->{'_longid'} || 0;
}

1;
