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
    my $phylipstream = new Bio::AlignIO(-format => 'phylip',
					-fh   => \*STDOUT);
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
				     -fh   => \*STDOUT);
    $gcgstream     =  new Bio::AlignIO(-format => 'msf',
				       -file   => 't/data/cysprot1a.msf');    

    while( my $aln = $gcgstream->next_aln ) {
	$phylipstream->write_aln($aln);
    }

=head1 DESCRIPTION

This object can transform Bio::SimpleAlign objects to and from PHYLIP
interleaved format. It will not work with PHYLIP sequencial format.

This module will output PHYLIP sequential format.  By specifying the flag
-interleaved => 0 in the initialization the module can output data in interleaved format.

=head1 FEEDBACK

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHORS - Heikki Lehvaslaiho and Jason Stajich

Email: heikki@ebi.ac.uk
Email: jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::AlignIO::phylip;
use vars qw(@ISA);
use strict;

use Bio::AlignIO;

@ISA = qw(Bio::AlignIO);


sub _initialize {
  my($self,@args) = @_;
  $self->SUPER::_initialize(@args);

  my ($interleave) = $self->_rearrange([qw(INTERLEAVED)],@args);
  if( ! defined $interleave ) { $interleave = 1 }  # this is the default
  $self->interleaved(1) if( $interleave);

  1;
}

=head2 next_aln

 Title   : next_aln
 Usage   : $aln = $stream->next_aln()
 Function: returns the next alignment in the stream.
           Throws an exception if trying to read in PHYLIP
           sequential format.
 Returns : SimpleAlign object
 Args    : 

=cut

sub next_aln {
    my $self = shift;
    my $entry;
    my ($seqcount, $residuecount, %hash, $name,$str,
	@names,$seqname,$start,$end,$count,$seq);
    
    my $aln =  Bio::SimpleAlign->new();
    $entry = $self->_readline and 
        ($seqcount, $residuecount) = $entry =~ /\s*(\d+)\s+(\d+)/;
    return 0 unless $seqcount and $residuecount;
    
    # first alignment section
    while( $entry = $self->_readline) {
	$entry =~ /^\s$/ and last;
	$entry =~ /^(.{10})\s+(.*)\s$/ && do {
	    $name = $1;
	    $str = $2;
	    $name =~ s/[\s\/]/_/g;
	    push @names, $name;
	    
	    $str =~ s/\s//g;
	    $count = scalar @names;
	    $hash{$count} = $str;
	};
	$self->throw("Not a valid interleaved PHYLIP file!") if $count > $seqcount; 
    }
    
    # interleaved sections
    $count = 0;
    while( $entry = $self->_readline) {
	$count = 0, next if $entry =~ /^\s$/;
	$entry =~ /\s*(.*)$/ && do {
	    $str = $1;
	    $str =~ s/\s//g;
	    $count++;
	    $hash{$count} .= $str;
	};
	$self->throw("Not a valid interleaved PHYLIP file!") if $count > $seqcount; 
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
	$self->throw("Length of sequence [$seqname] is not [$residuecount]! ") 
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
 Args    : Bio::SimpleAlign object

=cut

sub write_aln {
    my ($self,@aln) = @_;
    my $count = 0;
    my $wrapped = 0;
    my $maxname;
    my ($length,$date,$name,$seq,$miss,$pad,%hash,@arr,$tempcount,$index);
    
    foreach my $aln (@aln) {
	$self->throw("All sequences in the alignment must be the same length") 
	    unless $aln->is_flush ;

	$length  = $aln->length();
	$self->_print (sprintf(" %s %s\n", $aln->no_sequences, $aln->length));
	
	$aln->set_displayname_flat();
	foreach $seq ( $aln->each_seq() ) {
	    $name = $aln->displayname($seq->get_nse());	     
	    ($name) = substr($name,0,10);
	    $name = sprintf("%-10s",$name);
	    $name .= '   ' if( $self->interleaved());
	    $hash{$name} = $seq->seq();
	    push(@arr,$name);
	}

	if( $self->interleaved() ) {
	    while( $count < $length ) {	
		
		# there is another block to go!
		foreach $name ( @arr ) {
		    my $dispname = $name;
		    $dispname = '' if $wrapped;
		    $self->_print (sprintf("%13s  ",$dispname));
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
	} else {
	    foreach $name ( @arr ) {
		my $dispname = $name;
		$dispname = '' if $wrapped;
		$self->_print (sprintf("%s%s\n",$dispname,$hash{$name}));
	    }	
	}
    }
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

1;
