# $Id$
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
    my $phylipstream = new Bio::AlignIO(-format  => 'phylip',
                                        -fh      => \*STDOUT,
                                        -idlength=>30);
    # convert data from one format to another
    my $gcgstream     =  new Bio::AlignIO(-format => 'msf',
                                          -file   => 't/data/cysprot1a.msf');

    while( my $aln = $gcgstream->next_aln ) {
        $phylipstream->write_aln($aln);
    }

    # do it again with phylip sequential format format
    $phylipstream->interleaved(0);
    # can also initialize the object like this
    $phylipstream = new Bio::AlignIO(-interleaved => 0,
                                     -format => 'phylip',
                                     -fh   => \*STDOUT,
                                     -idlength=>10);
    $gcgstream     =  new Bio::AlignIO(-format => 'msf',
                                       -file   => 't/data/cysprot1a.msf');

    while( my $aln = $gcgstream->next_aln ) {
        $phylipstream->write_aln($aln);
    }

=head1 DESCRIPTION

This object can transform Bio::SimpleAlign objects to and from PHYLIP
interleaved format. It will not work with PHYLIP sequencial format.

This module will output PHYLIP sequential format.  By specifying the
flag -interleaved =E<gt> 0 in the initialization the module can output
data in interleaved format.

=head1 FEEDBACK

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

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
 Usage   : my $alignio = new Bio::AlignIO(-format => 'phylip'
					  -file   => '>file',
					  -idlength => 10,
					  -idlinebreak => 1);
 Function: Initialize a new L<Bio::AlignIO::phylip> reader or writer
 Returns : L<Bio::AlignIO> object
 Args    : [specific for writing of phylip format files]
           -idlength => integer - length of the id (will pad w/
						    spaces if needed)
           -interleaved => boolean - whether or not write as interleaved
                                     or sequential format
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

=cut

sub _initialize {
  my($self,@args) = @_;
  $self->SUPER::_initialize(@args);

  my ($interleave,$linelen,$idlinebreak,
      $idlength, $flag_SI, $tag_length,$ws) =
          $self->_rearrange([qw(INTERLEAVED
                                LINE_LENGTH
                                IDLINEBREAK
                                IDLENGTH
                                FLAG_SI
                                TAG_LENGTH
				WRAP_SEQUENTIAL)],@args);
  $self->interleaved(1) if( $interleave || ! defined $interleave);
  $self->idlength($idlength || $DEFAULTIDLENGTH);
  $self->id_linebreak(1) if( $idlinebreak );
  $self->line_length($linelen) if defined $linelen && $linelen > 0;
  $self->flag_SI(1) if ( $flag_SI );
  $self->tag_length($tag_length) if ( $tag_length || $DEFAULTTAGLEN );
  $self->wrap_sequential($ws ? 1 : 0);
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
    $entry = $self->_readline and
        ($seqcount, $residuecount) = $entry =~ /\s*(\d+)\s+(\d+)/;
    return 0 unless $seqcount and $residuecount;

    # first alignment section
    my $idlen = $self->idlength;
    $count = 0;
    my $iter = 1;
    my $interleaved = $self->interleaved;
    while( $entry = $self->_readline) {
	last if( $entry =~ /^\s?$/ && $interleaved );

	if( $entry =~ /^\s+(\d+)\s+(\d+)\s*$/) {
	    $self->_pushback($entry);
	    last;
	}
	if( $entry =~ /^\s+(.+)$/ ) {
	    $interleaved = 0;
	    $str = $1;
	    $str =~ s/\s//g;
	    $count = scalar @names;
	    $hash{$count} .= $str;

       	} elsif( $entry =~ /^(.{$idlen})\s+(.*)\s$/ ||
		 $entry =~ /^(.{$idlen})(\S{$idlen}\s+.+)\s$/ # Handle weirdnes s when id is too long
		 ) {
	    $name = $1;
	    $str = $2;
	    $name =~ s/[\s\/]/_/g;
	    $name =~ s/_+$//; # remove any trailing _'s

	    push @names, $name;
	    $str =~ s/\s//g;
	    $count = scalar @names;
	    $hash{$count} = $str;
	} elsif( $interleaved ) {
	    if( $entry =~ /^(\S+)\s+(.+)/ ||
		$entry =~ /^(.{$idlen})(.*)\s$/ ) {
		$name = $1;
		$str = $2;
		$name =~ s/[\s\/]/_/g;
		$name =~ s/_+$//; # remove any trailing _'s
		push @names, $name;
		$str =~ s/\s//g;
		$count = scalar @names;
		$hash{$count} = $str;
	    } else {
		$self->debug("unmatched line: $entry");
	    }
	}
	$self->throw("Not a valid interleaved PHYLIP file!") if $count > $seqcount;
    }

    if( $interleaved ) {
	# interleaved sections
	$count = 0;
	while( $entry = $self->_readline) {
            # finish current entry
	    if($entry =~/\s*\d+\s+\d+/){
		$self->_pushback($entry);
		last;
	    }
	    $count = 0, next if $entry =~ /^\s$/;
	    $entry =~ /\s*(.*)$/ && do {
		$str = $1;
		$str =~ s/\s//g;
		$count++;
		$hash{$count} .= $str;
	    };
	    $self->throw("Not a valid interleaved PHYLIP file! [$count,$seqcount] ($entry)") if $count > $seqcount;
	}
    }
    return 0 if scalar @names < 1;

    # sequence creation
    $count = 0;
    foreach $name ( @names ) {
	$count++;
	if( $name =~ /(\S+)\/(\d+)-(\d+)/ ) {
	    $seqname = $1;
	    $start = $2;
	    $end = $3;
	} else {
	    $seqname=$name;
	    $start = 1;
	    $str = $hash{$count};
	    $str =~ s/[^A-Za-z]//g;
	    $end = length($str);
	}
	# consistency test
	$self->throw("Length of sequence [$seqname] is not [$residuecount] it is ".CORE::length($hash{$count})."! ")
	    unless CORE::length($hash{$count}) == $residuecount;

       $seq = new Bio::LocatableSeq('-seq'=>$hash{$count},
				    '-id'=>$seqname,
				    '-start'=>$start,
				    '-end'=>$end,
				   );
	$aln->add_seq($seq);

   }
   return $aln;
}


=head2 write_aln

 Title   : write_aln
 Usage   : $stream->write_aln(@aln)
 Function: writes the $aln object into the stream in MSF format
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
                $self->_print (sprintf(" %s %s I\n", $aln->no_sequences, $aln->length));
            } else {
                $self->_print (sprintf(" %s %s S\n", $aln->no_sequences, $aln->length));
            }
        } else {
            $self->_print (sprintf(" %s %s\n", $aln->no_sequences, $aln->length));
        }

	$idlength = $self->idlength();
	$line_length = $self->line_length();
	$tag_length = $self->tag_length();
	foreach $seq ( $aln->each_seq() ) {
	    $name = $aln->displayname($seq->get_nse);
	    $name = substr($name, 0, $idlength) if length($name) > $idlength;
	    $name = sprintf("%-".$idlength."s",$name);
	    if( $self->interleaved() ) {
		$name .= '   ' ;
	    } elsif( $self->id_linebreak) {
		$name .= "\n";
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
			$self->_print (sprintf("%s ",substr($hash{$name},
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

sub interleaved{
   my ($self,$value) = @_;
   my $previous = $self->{'_interleaved'};
   if( defined $value ) {
       $self->{'_interleaved'} = $value;
   }
   return $previous;
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

1;
