
#
# BioPerl module for SimpleAlign
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

SimpleAlign - Multiple alignments held as a set of sequences

=head1 SYNOPSIS

    $aln = new Bio::SimpleAlign;
   
    $aln->read_MSF(\*STDIN);
 
    $aln->write_fasta(\*STDOUT);

=head1 INSTALLATION

This module is included with the central Bioperl distribution:

   http://bio.perl.org/Core/Latest
   ftp://bio.perl.org/pub/DIST

Follow the installation instructions included in the README file.

=head1 DESCRIPTION

SimpleAlign handles multiple alignments of sequences. It is very permissive
of types (it wont insist on things being all same length etc): really
it is a SequenceSet explicitly held in memory with a whole series of
built in manipulations and especially file format systems for
read/writing alignments.

SimpleAlign basically views an alignment as an immutable block of text.
SimpleAlign *is not* the object to be using if you want to manipulate an
alignment (eg, truncate an alignment or remove columns that are all gaps).
These functions are much better done by UnivAln by Georg Fuellen. 

However for lightweight display/formatting - this is the one to use.

Tricky concepts. SimpleAlign expects name,start,end to be 'unique' in
the alignment, and this is the key for the internal hashes.
(name,start,end is abreviated nse in the code). However, in many cases
people don't want the name/start-end to be displayed: either multiple
names in an alignment or names specific to the alignment
(ROA1_HUMAN_1, ROA1_HUMAN_2 etc). These names are called
'displayname', and generally is what is used to print out the
alignment. They default to name/start-end

The SimpleAlign Module came from Ewan Birney's Align module

=head1 PROGRESS 

SimpleAlign is being slowly converted to bioperl coding standards,
mainly by Ewan.

=over

=item Use Bio::Root::Object - done

=item Use proper exceptions - done

=item Use hashed constructor - not done!

=back

=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other Bioperl modules.
Send your comments and suggestions preferably to one of the Bioperl mailing lists.
Your participation is much appreciated.

    vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
    vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
    http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track the bugs and 
their resolution. Bug reports can be submitted via email or the web:

    bioperl-bugs@bio.perl.org                   
    http://bio.perl.org/bioperl-bugs/           

=head1 AUTHOR

Ewan Birney, birney@sanger.ac.uk

=head1 SEE ALSO

 Bio::Seq.pm - The biosequence object

 http://bio.perl.org/Projects/modules.html  - Online module documentation
 http://bio.perl.org/Projects/SeqAlign/     - Bioperl sequence alignment project
 http://bio.perl.org/                       - Bioperl Project Homepage

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SimpleAlign;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::Seq;         # uses Seq's as list 


@ISA = qw(Bio::Root::Object Exporter);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  # we need to set up internal hashs first!

  $self->{'seq'} = {};
  $self->{'order'} = {};
  $self->{'dis_name'} = {};
  $self->{'id'} = 'NoName';

  # maybe we should automatically read in from args. Hmmm...

# set stuff in self from @args
  return $make; # success - we hope!
}




=head2 id

 Title     : id
 Usage     : $myalign->id("Ig")
 Function  : Gets/sets the id field of the alignment
           :
 Returns   : An id string
 Argument  : An id string (optional)

=cut

sub id {
    my ($self, $name) = @_;

    if (defined( $name )) {
	$self->{'id'} = $name;
    }
    
    return $self->{'id'};
}




=head2 addSeq

 Title     : addSeq
 Usage     : $myalign->addSeq($newseq);
           : 
           :
 Function  : Adds another sequence to the alignment
           : *doesn't* align it - just adds it to the
           : hashes
           :
 Returns   : nothing
 Argument  : 

=cut


sub addSeq {
    my $self = shift;
    my $seq  = shift;
    my $order = shift;
    my ($name,$id,$start,$end);

    $id = $seq->id();
    $start = $seq->start();
    $end  = $seq->end(); 


    if( !defined $order ) {
	$order = keys %{$self->{'seq'}};
    }

    $name = sprintf("%s-%d-%d",$id,$start,$end);

    if( $self->{'seq'}->{$name} ) {
	$self->warn("Replacing one sequence [$name]\n");
    }

    $self->{'seq'}->{$name} = $seq;
    $self->{'order'}->{$order} = $name;

}


=head2 removeSeq

 Title     : removeSeq
 Usage     : $aln->removeSeq($seq);
 Function  : removes a single sequence from an alignment

=cut



sub removeSeq {
    my $self = shift;
    my $seq = shift;
    my ($name,$id,$start,$end);

    $seq->out_fasta(\*STDOUT);
    $id = $seq->id();
    $start = $seq->start();
    $end  = $seq->end(); 
    $name = sprintf("%s-%d-%d",$id,$start,$end);

    if( !exists $self->{'seq'}->{$name} ) {
	$self->throw("Sequence $name does not exist in the alignment to remove!");
    }

    delete $self->{'seq'}->{$name};

    # we can't do anything about the order hash but that is ok
    # because eachSeq will handle it
}

    
=head2 eachSeq

 Title     : eachSeq
 Usage     : foreach $seq ( $align->eachSeq() ) 
           : 
           :
 Function  : gets an array of Seq objects from the
           : alignment
           : 
           :
 Returns   : an array
 Argument  : nothing

=cut

sub eachSeq {
    my $self = shift;
    my (@arr,$order);

    foreach $order ( sort { $a <=> $b } keys %{$self->{'order'}} ) {
	if( exists $self->{'seq'}->{$self->{'order'}->{$order}} ) {
	    #print STDOUT sprintf("Putting %s\n", $self->{'order'}->{$order});
	    push(@arr,$self->{'seq'}->{$self->{'order'}->{$order}});
	}
    }

    return @arr;
}



=head2 consensus_string

 Title     : consensus_string
 Usage     : $str = $ali->consensus_string()
           : 
           :
 Function  : Makes a consensus
           : 
           : 
           :
 Returns   : 
 Argument  : 

=cut

sub consensus_string {
    my $self = shift;
    my $len;
    my ($out,$count);

    $out = "";

    $len = $self->length_aln();

    foreach $count ( 0 .. $len ) {
	$out .= $self->consensus_aa($count);
    }

    return $out;
    
}

sub consensus_aa {
    my $self = shift;
    my $point = shift;
    my ($seq,%hash,$count,$letter,$key);


    foreach $seq ( $self->eachSeq() ) {
	$letter = substr($seq->{'seq'},$point,1);
	($letter =~ /\./) && next;
	# print "Looking at $letter\n";
	$hash{$letter}++;
    }

    $count = -1;
    $letter = '?';

    foreach $key ( keys %hash ) {
	# print "Now at $key $hash{$key}\n";
	if( $hash{$key} > $count ) {
	    $letter = $key;
	    $count = $hash{$key};
	}
    }
    return $letter;
    
}



=head2 read_MSF

 Title   : read_MSF
 Usage   : $al->read_MSF(\*STDIN);
 Function: reads MSF formatted files. Tries to read *all* MSF
          It reads all non whitespace characters in the alignment
          area. For MSFs with weird gaps (eg ~~~) map them by using
          $al->map_chars('~','-');
 Example :
 Returns : 
 Args    : filehandle


=cut

sub read_MSF{
   my ($self,$fh) = @_;
   my (%hash,$name,$str,@names,$seqname,$start,$end,$count,$seq);

   # read in the name section

   while( <$fh> ) {
       /\/\// && last; # move to alignment section
       /Name:\s+(\S+)/ && do { $name = $1;
			       $hash{$name} = ""; # blank line
			       push(@names,$name); # we need it ordered!
			   };
       # otherwise - skip
   }

   # alignment section

   while( <$fh> ) {
       /^\s+(\S+)\s+(.*)$/ && do {
	   $name = $1;
	   $str = $2;
	   if( ! exists $hash{$name} ) {
	       $self->throw("$name exists as an alignment line but not in the header. Not confident of what is going on!");
	   }
	   $str =~ s/\s//g;
	   $hash{$name} .= $str;
       };
   }

   # now got this as a name - sequence hash. Lets make some sequences!

   $count = 0;

   foreach $name ( @names ) {
       if( $name =~ /(\S+)\/(\d+)-(\d+)/ ) {
	   $seqname = $1;
	   $start = $2;
	   $end = $3;
       } else {
	   $seqname=$name;
	   $start = 1;
	   $str = $hash{$name};
	   $str =~ s/[^A-Za-z]//g;
	   $end = length($str);
       }
    
       $seq = new Bio::Seq('-seq'=>$hash{$name},
			   '-id'=>$seqname,
			   '-start'=>$start,
			   '-end'=>$end, 
			   '-type'=>'aligned');
    
       $self->addSeq($seq);
    
       $count++;
   }

   return $count;
}

=head2 write_MSF

 Title     : write_MSF
 Usage     : $ali->write_MSF(\*FH)
           : 
           :
 Function  : writes MSF format output
           : 
           : 
           :
 Returns   : 
 Argument  : 

=cut

sub write_MSF {
    my $self = shift;
    my $file = shift;
    my $msftag;
    my $type;
    my $count = 0;
    my $maxname;
    my ($length,$date,$name,$seq,$miss,$pad,%hash,@arr,$tempcount,$index);
    
    $date = localtime(time);
    $msftag = "MSF";
    $type = "P";
    $maxname = $self->maxnse_length();
    $length  = $self->length_aln();
    $name = $self->id();
    if( !defined $name ) {
	$name = "Align";
    }
    
    print $file sprintf("\n%s   MSF: %d  Type: P  %s  Check: 00 ..\n\n",$name,$self->no_sequences,$date);

    foreach $seq ( $self->eachSeq() ) {
	$name = $self->get_displayname($seq->get_nse());
	$miss = $maxname - length ($name);
	$miss += 2;
	$pad  = " " x $miss;

	print $file sprintf(" Name: %s%sLen:    %d  Check:  %d  Weight:  1.00\n",$name,$pad,length $seq->str(),$seq->GCG_checksum());
	
	$hash{$name} = $seq->str();
	push(@arr,$name);
    }

    #
    # ok - heavy handed, but there you go.
    #
    print "\n//\n";

    while( $count < $length ) {
	
	# there is another block to go!
    
	foreach $name  ( @arr ) {
	    print $file sprintf("%20s  ",$name);
	    
	    $tempcount = $count;
	    $index = 0;
	    while( ($tempcount + 10 < $length) && ($index < 5)  ) {
		
		print $file sprintf("%s ",substr($hash{$name},$tempcount,10));
				    
		$tempcount += 10;
		$index++;
	    }

	    #
	    # ok, could be the very last guy ;)
	    #

	    if( $index < 5) {
		
		#
		# space to print! 
		#

		print $file sprintf("%s ",substr($hash{$name},$tempcount));
		$tempcount += 10;
	    }

	    print $file "\n";
	} # end of each sequence

	print "\n\n";

	$count = $tempcount;
    }				
}


sub maxname_length {
    my $self = shift;
    my $maxname = (-1);
    my ($seq,$len);

    foreach $seq ( $self->eachSeq() ) {
	$len = length $seq->id();

	if( $len > $maxname ) {
	    $maxname = $len;
	}
    }

    return $maxname;
}

sub maxnse_length {
    my $self = shift;
    my $maxname = (-1);
    my ($seq,$len);

    foreach $seq ( $self->eachSeq() ) {
	$len = length $seq->get_nse();

	if( $len > $maxname ) {
	    $maxname = $len;
	}
    }

    return $maxname;
}

sub maxdisplayname_length {
    my $self = shift;
    my $maxname = (-1);
    my ($seq,$len);

    foreach $seq ( $self->eachSeq() ) {
	$len = length $self->get_displayname($seq->get_nse());

	if( $len > $maxname ) {
	    $maxname = $len;
	}
    }

    return $maxname;
}


=head2 length_aln

 Title     : length_aln()
 Usage     : $len = $ali->length_aln() 
           : 
           :
 Function  : returns the maximum length of the alignment.
           : To be sure the alignment is a block, use is_flush
           : 
           :
 Returns   : 
 Argument  : 

=cut

sub length_aln {
    my $self = shift;
    my $seq;
    my $length = (-1);
    my ($temp,$len);

    foreach $seq ( $self->eachSeq() ) {
	$temp = length($seq->str());
	if( $temp > $length ) {
	    $length = $temp;
	}
    }

    return $length;
}


=head2 is_flush

 Title     : is_flush
 Usage     : if( $ali->is_flush() )  
           : 
           :
 Function  : Tells you whether the alignment 
           : is flush, ie all of the same length
           : 
           :
 Returns   : 1 or 0
 Argument  : 

=cut

sub is_flush {
    my $self = shift;
    my $seq;
    my $length = (-1);
    my $temp;
    
    foreach $seq ( $self->eachSeq() ) {
	if( $length == (-1) ) {
	    $length = length($seq->str());
	    next;
	}

	$temp = length($seq->str());

	if( $temp != $length ) {
	    return 0;
	}
    }

    return 1;
}


=head2 read_fasta

 Title     : read_fasta
 Usage     : $ali->read_fasta(\*INPUT)
           : 
           :
 Function  : reads in a fasta formatted
           : file for an alignment
           : 
           :
 Returns   : 
 Argument  : 

=cut

sub read_fasta {
    my $self = shift;
    my $in = shift;
    my $count = 0;
    my ($start,$end,$name,$seqname,$seq,$seqchar,$tempname,%align);

    while( <$in> ) {
	if( /^>(\S+)/ ) {
	    $tempname = $1;
	    if( defined $name ) {
		# put away last name and sequence 

		if( $name =~ /(\S+)\/(\d+)-(\d+)/ ) {
		    $seqname = $1;
		    $start = $2;
		    $end = $3;
		} else {
		    $seqname=$name;
		    $start = 1;
		    $end = length($align{$name});
		}
		
		$seq = new Bio::Seq('-seq'=>$seqchar,
				    '-id'=>$seqname,
				    '-start'=>$start,
				    '-end'=>$end, 
				    '-type'=>'aligned');

		$self->addSeq($seq);

		$count++;
	    }
	    $name = $tempname;
	    $seqchar  = "";
	    next;
	}
	s/[^A-Za-z\.\-]//g;
	$seqchar .= $_;

    }
    # put away last name and sequence 
    
    if( $name =~ /(\S+)\/(\d+)-(\d+)/ ) {
	$seqname = $1;
	$start = $2;
	$end = $3;
    } else {
	$seqname=$name;
	$start = 1;
	$end = length($align{$name});
    }
    
    $seq = new Bio::Seq('-seq'=>$seqchar,
			'-id'=>$seqname,
			'-start'=>$start,
			'-end'=>$end, 
			'-type'=>'aligned');
    
    $self->addSeq($seq);
    
    $count++;

    return $count;
}
	    

=head2 read_selex

 Title     : read_selex
 Usage     : $ali->read_selex(\*INPUT) 
           : 
           :
 Function  : reads selex (hmmer) format
           : alignments
           : 
           :
 Returns   : 
 Argument  : 

=cut

sub read_selex {
    my $self = shift;
    my $in = shift;
    my ($start,$end,%align,$name,$seqname,$seq,$count,%hash,%c2name,$no);
    
    # not quite sure 'what' selex format is, but here we go!

    while( <$in> ) {
	#/^#=SQ/ && last;
	/^#=RF/ && last;
    }

    while( <$in> ) {
	/^#/ && next;
	/^\s/ && next;
	/\w/ && last;
    }

    if( eof($in) ) {
	warn("Exhausted file without reading selex format");
	return undef;
    }
    # first line starting with word
    chop;

    $count =0;

    if( !/^(\S+)\s+([A-Za-z\.\-]+)\s*/ ) {
	warn("Ok. I think [$_] does not look like selex format. Not sure...");
	return (-1);
    }

    $name = $1;
    $seq = $2;

    $align{$name} .= $seq;
    $c2name{$count} = $name;
    $count++;

    while( <$in> ) {
	!/^([^\#]\S+)\s+([A-Za-z\.\-]+)\s*/ && next;
	
	$name = $1;
	$seq = $2;

	if( ! defined $align{$name}  ) {
	    $c2name{$count} = $name;
	    $count++;
	}

	$align{$name} .= $seq;
    }

    # ok... now we can make the sequences

    $count = 0;
    foreach $no ( sort { $a <=> $b } keys %c2name ) {
	$name = $c2name{$no};

	if( $name =~ /(\S+)\/(\d+)-(\d+)/ ) {
	    $seqname = $1;
	    $start = $2;
	    $end = $3;
	} else {
	    $seqname=$name;
	    $start = 1;
	    $end = length($align{$name});
	}

	$seq = new Bio::Seq('-seq'=>$align{$name},
			    '-id'=>$seqname,
			    '-start'=>$start,
			    '-end'=>$end, 
			    '-type'=>'aligned'
			    );

	$self->addSeq($seq);

	$count++;
    }

    return $count;
}


=head2 read_mase

 Title     : read_mase
 Usage     : $ali->read_mase(\*INPUT)
           : 
           :
 Function  : reads mase (seaview) 
           : formatted alignments
           : 
           :
 Returns   : 
 Argument  : 

=cut
	
sub read_mase {
    my $self = shift;
    my $in = shift;
    my $name;
    my $start;
    my $end;
    my $seq;
    my $add;
    my $count = 0;
	
	
    while( <$in> ) {
	/^;/ && next;
	if( /^(\S+)\/(\d+)-(\d+)/ ) {
	    $name = $1;
	    $start = $2;
	    $end = $3;
	} else {
	    s/\s//g;
	    $name = $_;
	    $end = -1;
	}

	$seq = "";

	while( <$in> ) {
	    /^;/ && last;
	    s/[^A-Za-z\.\-]//g;
	    $seq .= $_;
	}
	if( $end == -1) {
	    $start = 1;
	    $end = length($seq);
	}

	$add = new Bio::Seq('-seq'=>$seq,
			    '-id'=>$name,
			    '-start'=>$start,
			    '-end'=>$end, 
			    '-type'=>'aligned');


	
	$self->addSeq($add);

	$count++;
    }

    return $count;
}

	

=head2 read_Pfam_file

 Title     : read_Pfam_file
 Usage     : $ali->read_Pfam_file("thisfile");
           : 
 Function  : opens a filename, reads
           : a Pfam (mul) formatted alignment
           :
           : 
           :
 Returns   : 
 Argument  : 

=cut

sub read_Pfam_file {
    my $self = shift;
    my $file = shift;
    my $out;

    if( open(_TEMPFILE,$file) == 0 ) {
	$self->throw("Could not open $file for reading Pfam style alignment");
	return 0;
    }

    $out = read_Pfam($self,\*_TEMPFILE);
    
    close(_TEMPFILE);

    return $out;
}


=head2 read_Pfam

 Title     : read_Pfam
 Usage     : $ali->read_Pfam(\*INPUT)
           : 
           :
 Function  : reads a Pfam formatted
           : Alignment (Mul format).
           : - this is the format used by Belvu
           :
 Returns   : 
 Argument  : 

=cut

sub read_Pfam {
    my $self = shift;
    my $in = shift;
    my $name;
    my $start;
    my $end;
    my $seq;
    my $add;
    my $acc;
    my %names;
    my $count = 0;

    while( <$in> ) {
	chop;
	/^\/\// && last;
	if( !/^(\S+)\/(\d+)-(\d+)\s+(\S+)\s*/ ) { 
	    $self->throw("Found a bad line [$_] in the pfam format alignment");
	    next;
	}

	$name = $1;
	$start = $2;
	$end = $3;
	$seq = $4;


	# there may be an accession number at the end of the line

	($acc) = ($' =~ /\s*(\S+)\s*/);
	$names{'acc'} = $acc;

	$add = new Bio::Seq('-seq'=>$seq,
			    '-id'=>$name,
			    '-names'=>\%names,
			    '-start'=>$start,
			    '-end'=>$end, 
			    '-type'=>'aligned');
     
	$self->addSeq($add);
	
	$count++;
    }

    return $count;
}

sub write_Pfam_link {
    my $self = shift;
    my $link = shift;
    my $out  = shift;
    my $len  = shift;
    my $acc   = shift;
    my ($namestr,$linkstr,$seq,$name,$start,$end,$disname,$add,$place);

    if( !defined $len ) {
	$len = 22;
    }

    foreach $seq ( $self->eachSeq() ) {
	$name = $seq->id();
	$disname = $self->get_displayname($seq->get_nse());

	$add = $len - length($disname);
	$place = " " x $add;
	
	
#	print $out "Going to eval \$linkstr = $link\n";
	eval ("\$linkstr = \"$link\"");
	if( defined $acc ) { 
	    print $out sprintf("%s%s%s %s\n",$linkstr,$place,$seq->str(),$seq->names()->{'acc'});
	} else {
	    print $out sprintf("%s%s%s\n",$linkstr,$place,$seq->str());
	}
    }
}


=head2 write_Pfam

 Title     : write_Pfam
 Usage     : $ali->write_Pfam(\*OUTPUT) 
           : 
           :
 Function  : writes a Pfam/Mul formatted
           : file
           : 
           :
 Returns   : 
 Argument  : 

=cut

sub write_Pfam {
    my $self = shift;
    my $out  = shift;
    my $acc  = shift;
    my ($namestr,$seq,$add);
    my ($maxn);

    $maxn = $self->maxdisplayname_length();
    
    foreach $seq ( $self->eachSeq() ) {
	$namestr = $self->get_displayname($seq->get_nse());
	$add = $maxn - length($namestr) + 2;
	$namestr .= " " x $add;
	if( defined $acc && $acc == 1) {
	    print $out sprintf("%s  %s %s\n",$namestr,$seq->str(),$seq->names()->{'acc'});
	} else {
	    print $out sprintf("%s  %s\n",$namestr,$seq->str());
	}
    }
}


=head2 write_clustalw

 Title     : write_clustalw
 Usage     : $ali->write_clustalw 
           : 
           :
 Function  : writes a clustalw formatted
           : (.aln) file
           : 
           :
 Returns   : 
 Argument  : 

=cut

sub write_clustalw {
    my $self = shift;
    my $file = shift;
    my ($count,$length,$seq,@seq,$tempcount);

    print $file "CLUSTAL W(1.4) multiple sequence alignment\n\n\n";
    
    $length = $self->length_aln();
    $count = $tempcount = 0;
    @seq = $self->eachSeq();

    while( $count < $length ) {
	foreach $seq ( @seq ) {
	    print $file sprintf("%-22s %s\n",$self->get_displayname($seq->get_nse()),substr($seq->str(),$count,50));
	}
	print $file "\n\n";
	$count += 50;
    }

    
}


=head2 write_fasta

 Title     : write_fasta
 Usage     : $ali->write_fasta(\*OUTPUT) 
           : 
           :
 Function  : writes a fasta formatted alignment
           : 
 Returns   : 
 Argument  : reference-to-glob to file or filehandle object 

=cut

sub write_fasta {
    my $self = shift;
    my $file  = shift;
    my ($seq,$rseq,$name,$count,$length,$seqsub);

    foreach $rseq ( $self->eachSeq() ) {
	$name = $self->get_displayname($rseq->get_nse());
	$seq  = $rseq->str();
	
	print $file ">$name\n";
	
	$count =0;
	$length = length($seq);
	while( ($count * 60 ) < $length ) {
	    $seqsub = substr($seq,$count*60,60);
	    print $file "$seqsub\n";
	    $count++;
	}
    }
}



sub set_displayname {
    my $self = shift;
    my $name = shift;
    my $disname = shift;

    # print "Setting $name to $disname\n";
    $self->{dis_name}->{$name} = $disname;
}

sub get_displayname {
    my $self = shift;
    my $name = shift;

    if( defined $self->{dis_name}->{$name} ) {
	return  $self->{dis_name}->{$name};
    } else {
	return $name;
    }
}

=head2 set_displayname_flat

 Title     : set_displayname_flat
 Usage     : $ali->set_displayname_flat() 
           : 
           :
 Function  : Makes all the sequences be displayed
           : as just their name, not name/start-end
           : 
           :
 Returns   : 
 Argument  : 

=cut

sub set_displayname_flat {
    my $self = shift;
    my ($nse,$seq);

    foreach $seq ( $self->eachSeq() ) {
	$nse = $seq->get_nse();
	$self->set_displayname($nse,$seq->id());
    }
}

=head2 set_displayname_normal

 Title     : set_displayname_normal
 Usage     : $ali->set_displayname_normal() 
           : 
           :
 Function  : Makes all the sequences be displayed
           : as name/start-end
           : 
           :
 Returns   : 
 Argument  : 

=cut

sub set_displayname_normal {
    my $self = shift;
    my ($nse,$seq);

    foreach $seq ( $self->eachSeq() ) {
	$nse = $seq->get_nse();
	$self->set_displayname($nse,sprintf("%s/%d-%d",$seq->id(),$seq->start(),$seq->end()));
    }
}

=head2 set_displayname_count

 Title     : set_displayname_count
 Usage     : $ali->set_displayname_count 
           : 
           :
 Function  : sets the names to be name_#
           : where # is the number of times this
           : name has been used. 
           :
 Returns   : 
 Argument  : 

=cut


sub set_displayname_count {
    my $self= shift;
    my (@arr,$name,$seq,$count,$temp,$nse);

    foreach $seq ( $self->each_alphabetically() ) {
	$nse = $seq->get_nse();

	#name will be set when this is the second
	#time (or greater) is has been seen

	if( $name eq ($seq->id()) ) {
	    $temp = sprintf("%s_%s",$name,$count);
	    $self->set_displayname($nse,$temp);
	    $count++;
	} else {
	    $count = 1;
	    $name = $seq->id();
	    $temp = sprintf("%s_%s",$name,$count);
	    $self->set_displayname($nse,$temp);
	    $count++;
	}
    }

}

=head2 each_alphabetically

 Title     : each_alphabetically
 Usage     : foreach $seq ( $ali->each_alphabetically() )
           : 
           :
 Function  : returns an array of sequence object sorted
           : alphabetically by name and then by start point
           : 
           : Does not change the order of the alignment
 Returns   : 
 Argument  : 

=cut

	
sub each_alphabetically {
    my $self = shift;
    my ($seq,$nse,@arr,%hash,$count);

    foreach $seq ( $self->eachSeq() ) {
	$nse = $seq->get_nse("-","-");
	$hash{$nse} = $seq;
    }


    foreach $nse ( sort alpha_startend keys %hash) {
	push(@arr,$hash{$nse});
    }

    return @arr;

}

   
=head2 sort_alphabetically

 Title     : sort_alphabetically
 Usage     : $ali->sort_alphabetically
           : 
           :
 Function  : changes the order of the alignemnt
           : to alphabetical on name followed by
           : numerical by number
           :
 Returns   : 
 Argument  : 

=cut
 
sub sort_alphabetically {
    my $self = shift;
    my ($seq,$nse,@arr,%hash,$count);


    foreach $seq ( $self->eachSeq() ) {
	$nse = $seq->get_nse("-","-");
	$hash{$nse} = $seq;
    }

    $count = 0;

    %{$self->{'order'}} = (); # reset the hash;

    foreach $nse ( sort alpha_startend keys %hash) {
	$self->{'order'}->{$count} = $nse;

	$count++;
    }

}

=head2 map_chars

 Title     : map_chars
 Usage     : $ali->map_chars('\.','-')
           : 
           :
 Function  : does a s/$arg1/$arg2/ on 
           : the sequences. Useful for
           : gap characters
           :
           : Notice that the from (arg1) is interpretted 
           : as a regex, so be careful about quoting meta
           : characters (eg $ali->map_chars('.','-') wont
           : do what you want)
 Returns   : 
 Argument  : 

=cut


sub map_chars {
    my $self = shift;
    my $from = shift;
    my $to   = shift;
    my ($seq,$temp);

    foreach $seq ( $self->eachSeq() ) {
	$temp = $seq->str();
	$temp =~ s/$from/$to/g;
	$seq->setseq($temp);
    }
}

=head2 uppercase

 Title     : uppercase()
 Usage     : $ali->uppercase()
           : 
           :
 Function  : Sets all the sequences
           : to uppercase
           : 
           :
 Returns   : 
 Argument  : 

=cut

sub uppercase {
    my $self = shift;
    my $seq;
    my $temp;

    foreach $seq ( $self->eachSeq() ) {
      $temp = $seq->str();
      $temp =~ tr/[a-z]/[A-Z]/;

      $seq->setseq($temp);
    }
}


sub alpha_startend {
    my ($aname,$astart,$bname,$bstart);
    ($aname,$astart) = split (/-/,$a);
    ($bname,$bstart) = split (/-/,$b);

    if( $aname eq $bname ) {
	return $astart <=> $bstart;
    }
    else { 
	return $aname cmp $bname;
    }

}

=head2 no_sequences

 Title     : no_sequences
 Usage     : $depth = $ali->no_sequences
           : 
           :
 Function  : number of sequence in the
           : sequence alignment
           : 
           :
 Returns   : 
 Argument  : 

=cut
    
sub no_sequences {
    my $self = shift;

    return keys %{$self->{'seq'}};
}

=head2 no_residues

 Title     : no_residues
 Usage     : $no = $ali->no_residues
           : 
           :
 Function  : number of residues in total
           : in the alignment
           : 
           :
 Returns   : 
 Argument  : 

=cut


sub no_residues {
    my $self = shift;
    my $count = 0;
    my ($key,$temp);

    foreach $key ( keys %{$self->{'seq'}} ) {

	$temp = $self->{'seq'}->{$key}->str();


	$temp =~ s/[^A-Za-z]//g;
	
	$count += length($temp);
    }

    return $count;
}



=head2 purge

 Title   : purge
 Usage   : $aln->purge(0.7);
 Function: removes sequences above whatever %id
 Example :
 Returns : An array of the removed sequences
 Arguments

 This function will grind on large alignments. Beware!

 (perhaps not ideally implemented)

=cut

sub purge{
  my ($self,$perc) = @_;
  my (@seqs,$seq,%removed,$i,$j,$count,@one,@two,$seq2,$k,$res,$ratio,@ret);
  
 # $self->write_Pfam(\*STDOUT);

  @seqs = $self->eachSeq();
  
  #$self->write_Pfam(\*STDOUT);

#  foreach $seq ( @seqs ) {
#      printf("$seq %s %s\n",$seq->get_nse(),join(' ',$seq->dump()));
#  }

  for($i=0;$i< @seqs;$i++ ) {
      $seq = $seqs[$i];
      #printf "%s\n", $seq->out_fasta();
      #print "\n\nDone\n\n";

      # if it has already been removed, skip

      if( $removed{$seq->get_nse()} == 1 ) {
	  next;
      }

      # if not ... look at the other sequences
      
      # make the first array once

      @one = $seq->seq();
      for($j=$i+1;$j < @seqs;$j++) {
	  $seq2 = $seqs[$j];
	  if ( $removed{$seq2->get_nse()} == 1 ) {
	      next;
	  }
	  @two = $seq2->seq();
	  $count = 0;
	  $res = 0;
	  for($k=0;$k<@one;$k++) {
	      if( $one[$k] ne '.' && $one[$k] ne '-' && $one[$k] eq $two[$k]) {
		  $count++;
	      }
	      if( $one[$k] ne '.' && $one[$k] ne '-' && $two[$k] ne '.' && $two[$k] ne '-' ) {
		  $res++;
	      }
	  }
	  if( $res == 0 ) {
	      $ratio = 0;
	  } else {
	      $ratio = $count/$res;
	  }

	  if( $ratio > $perc ) {
	      $removed{$seq2->get_nse()} = 1;
	      $self->removeSeq($seq2);
	      push(@ret,$seq2);
	  } else {
	      # could put a comment here!
	  }
      }
  }
  
  return @ret;
}



=head2 percentage_identity

 Title   : percentage_identity
 Usage   : $id = $align->percentage_identity
 Function:
    The function uses a fast method to calculate the average percentage identity of the alignment
 Returns : The average percentage identity of the alignment
 Args    : None

=cut

sub percentage_identity{
   my ($self,@args) = @_;

   my @alphabet = ('A','B','C','D','E','F','G','H','I','J','K','L','M',
                   'N','O','P','Q','R','S','T','U','V','W','X','Y','Z');

   my ($len, $total, $subtotal, $divisor, $subdivisor, @seqs, @countHashes);

   if (! $self->is_flush()) {
       $self->throw("All sequences in the alignment must be the same length");
   }

   @seqs = $self->eachSeq();
   $len = $self->length_aln();

   # load the each hash with correct keys for existence checks
   for( my $index=0; $index < $len; $index++) { 
       foreach my $letter (@alphabet) {
	   $countHashes[$index]->{$letter} = 0;
       }
   }


   foreach my $seq (@seqs)  {
       my @seqChars = $seq->ary();
       for( my $column=0; $column < @seqChars; $column++ ) {
	   my $char = uc($seqChars[$column]);
	   if (exists $countHashes[$column]->{$char}) {
	       $countHashes[$column]->{$char}++;
	   }
       }
   }

   $total = 0;
   $divisor = 0;
   for(my $column =0; $column < $len; $column++) {
       my %hash = %{$countHashes[$column]};
       $subdivisor = 0;
       foreach my $res (keys %hash) {
	   $total += $hash{$res}*($hash{$res} - 1);
	   $subdivisor += $hash{$res};
       }
       $divisor += $subdivisor * ($subdivisor - 1);
   }
   return ($total / $divisor )*100.0;
}


=head2 read_Prodom

 Title   : read_Prodom
 Usage   : $ali->read_Prodom( $file )
 Function: Reads in a Prodom format alignment
 Returns : 
    Args    : A filehandle glob or ref. to a filehandle object

=cut

sub read_Prodom{
   my $self = shift;
   my $file = shift;

   my ($acc, $fake_id, $start, $end, $seq, $add, %names);

   while (<$file>){
       if (/^AL\s+(\S+)\|(\S+)\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)$/){
	   $acc=$1;
	   $fake_id=$2;  # Accessions have _species appended
	   $start=$3;
	   $end=$4;
	   $seq=$5;
	   
	   $names{'fake_id'} = $fake_id;

	   $add = new Bio::Seq('-seq'=>$seq,
			       '-id'=>$acc,
			       '-names'=>\%names,
			       '-start'=>$start,
			       '-end'=>$end, 
			       '-type'=>'aligned');
	   
	   $self->addSeq($add);
        }
    }
}


1;

