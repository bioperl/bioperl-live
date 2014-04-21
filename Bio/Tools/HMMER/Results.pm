#
# Perl Module for HMMResults
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
#Copyright Genome Research Limited (1997).

=head1 NAME

Bio::Tools::HMMER::Results - Object representing HMMER output results

=head1 SYNOPSIS

   # parse a hmmsearch file (can also parse a hmmpfam file)
   $res = Bio::Tools::HMMER::Results->new( -file => 'output.hmm' ,
                                          -type => 'hmmsearch');

   # print out the results for each sequence
   foreach $seq ( $res->each_Set ) {
       print "Sequence bit score is",$seq->bits,"\n";
       foreach $domain ( $seq->each_Domain ) {
           print " Domain start ",$domain->start," end ",$domain->end,
	   " score ",$domain->bits,"\n";
       }
   }

   # new result object on a sequence/domain cutoff of
   # 25 bits sequence, 15 bits domain
   $newresult = $res->filter_on_cutoff(25,15);

   # alternative way of getting out all domains directly
   foreach $domain ( $res->each_Domain ) {
       print "Domain on ",$domain->seq_id," with score ",
       $domain->bits," evalue ",$domain->evalue,"\n";
   }

=head1 DESCRIPTION

This object represents HMMER output, either from hmmsearch or
hmmpfam. For hmmsearch, a series of HMMER::Set objects are made, one
for each sequence, which have the the bits score for the object. For
hmmpfam searches, only one Set object is made.


These objects come from the original HMMResults modules used
internally in Pfam, written by Ewan Birney. Ewan then converted them to
BioPerl objects in 1999. That conversion is meant to be backwardly
compatible, but may not be (caveat emptor).

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

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

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

=head1 CONTRIBUTORS

Jason Stajich, jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::Tools::HMMER::Results;

use strict;

use Bio::Tools::HMMER::Domain;
use Bio::Tools::HMMER::Set;
use Symbol;

use base qw(Bio::Root::Root Bio::Root::IO Bio::SeqAnalysisParserI);

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);

  $self->{'domain'} = []; # array of HMMUnits
  $self->{'seq'} = {};

  my ($parsetype) = $self->_rearrange([qw(TYPE)],@args);
  $self->_initialize_io(@args);
  if( !defined $parsetype ) {
      $self->throw("No parse type provided. should be hmmsearch or hmmpfam");
  }
  $self->parsetype($parsetype);
  if( defined $self->_fh() ) {
      if( $parsetype eq 'hmmsearch' ) {
	  $self->_parse_hmmsearch($self->_fh());
      } elsif ( $parsetype eq 'hmmpfam' ) {
	  $self->_parse_hmmpfam($self->_fh());
      } else {
	  $self->throw("Did not recoginise type $parsetype");
      }
  }

  return $self; # success - we hope!
}


=head2 next_feature

 Title   : next_feature
 Usage   : while( my $feat = $res->next_feature ) { # do something }
 Function: SeqAnalysisParserI implementing function
 Example :
 Returns : A Bio::SeqFeatureI compliant object, in this case,
           each DomainUnit object, ie, flattening the Sequence
           aspect of this.
 Args    : None


=cut

sub next_feature{
   my ($self) = @_;

   if( $self->{'_started_next_feature'} == 1 ) {
       return shift @{$self->{'_next_feature_array'}};
   } else {
       $self->{'_started_next_feature'} = 1;
       my @array;
       foreach my $seq ( $self->each_Set() ) {
	   foreach my $unit ( $seq->each_Domain() ) {
	       push(@array,$unit);
	   }
       }
       my $res = shift @array;
       $self->{'_next_feature_array'} = \@array;
       return $res;
   }

   $self->throw("Should not reach here! Error!");
}


=head2 number

 Title   : number
 Usage   : print "There are ",$res->number," domains hit\n";
 Function: provides the number of domains in the HMMER report

=cut

sub number {
    my $self = shift;
    my @val;
    my $ref;
    $ref = $self->{'domain'};


    @val = @{$self->{'domain'}};
    return scalar @val;
}

=head2 seqfile

 Title   : seqfile
 Usage   : $obj->seqfile($newval)
 Function:
 Example :
 Returns : value of seqfile
 Args    : newvalue (optional)


=cut

sub seqfile{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'seqfile'} = $value;
    }
    return $self->{'seqfile'};

}

=head2 hmmfile

 Title   : hmmfile
 Usage   : $obj->hmmfile($newval)
 Function:
 Example :
 Returns : value of hmmfile
 Args    : newvalue (optional)


=cut

sub hmmfile{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'hmmfile'} = $value;
    }
    return $self->{'hmmfile'};

}

=head2 add_Domain

 Title   : add_Domain
 Usage   : $res->add_Domain($unit)
 Function: adds a domain to the results array. Mainly used internally.
 Args    : A Bio::Tools::HMMER::Domain


=cut

sub add_Domain {
    my $self = shift;
    my $unit = shift;
    my $name;

    $name = $unit->seq_id();

    if( ! exists $self->{'seq'}->{$name} ) {
	$self->warn("Adding a domain of $name but with no HMMSequence. Will be kept in domain array but not added to a HMMSequence");
    } else {
	$self->{'seq'}->{$name}->add_Domain($unit);
    }
    push(@{$self->{'domain'}},$unit);
}


=head2 each_Domain

 Title   : each_Domain
 Usage   : foreach $domain ( $res->each_Domain() )
 Function: array of Domain units which are held in this report
 Returns : array
 Args    : none


=cut

sub each_Domain {
    my $self = shift;
    my (@arr,$u);

    foreach $u ( @{$self->{'domain'}} ) {
	push(@arr,$u);
    }

    return @arr;
}


=head2 domain_bits_cutoff_from_evalue

 Title   : domain_bits_cutoff_from_evalue
 Usage   : $cutoff = domain_bits_cutoff_from_evalue(0.01);
 Function: return a bits cutoff from an evalue using the
           scores here. Somewhat interesting logic:
            Find the two bit score which straddle the evalue
            if( 25 is between these two points) return 25
            else return the midpoint.

           This logic tries to ensure that with large signal to
           noise separation one still has sensible 25 bit cutoff
 Returns :
 Args    :

=cut

sub domain_bits_cutoff_from_evalue {
    my $self = shift;
    my $eval = shift;
    my ($dom,$prev,@doms,$cutoff,$sep,$seen);

    @doms = $self->each_Domain;


    @doms = map { $_->[0] }
            sort { $b->[1] <=> $a->[1] }
            map { [ $_, $_->bits] } @doms;
    $seen = 0;
    foreach $_ ( @doms ) {
	if( $_->evalue > $eval ) {
	    $seen = 1;
	    $dom = $_;
	    last;
	}
	$prev = $_;
    }

    if( ! defined $prev || $seen == 0) {
	$self->throw("Evalue is either above or below the list...");
	return;
    }

    $sep = $prev->bits - $dom->bits ;

    if( $sep < 1 ) {
	return $prev->bits();
    }
    if( $dom->bits < 25 && $prev->bits > 25 ) {
	return 25;
    }

    return int( $dom->bits + $sep/2 ) ;

}


sub dictate_hmm_acc {
    my $self = shift;
    my $acc = shift;
    my ($unit);


    foreach $unit ( $self->eachHMMUnit() ) {
	$unit->hmmacc($acc);
    }
}

=head2 write_FT_output

 Title   : write_FT_output
 Usage   : $res->write_FT_output(\*STDOUT,'DOMAIN')
 Function: writes feature table output ala swissprot
 Returns :
 Args    :


=cut

sub write_FT_output {
    my $self = shift;
    my $file = shift;
    my $idt  = shift;
    my ($seq,$unit);

    if( !defined $idt ) {
	$idt = "DOMAIN";
    }

    foreach $seq ( $self->each_Set() ) {
	print $file sprintf("ID   %s\n",$seq->name());
	foreach $unit ( $seq->each_Domain() ) {
	    print $file sprintf("FT   %s   %d %d %s\n",$idt,
				$unit->start,$unit->end,$unit->hmmname);
	}
	print $file "//\n";
    }
}

=head2 filter_on_cutoff

 Title   : filter_on_cutoff
 Usage   : $newresults = $results->filter_on_cutoff(25,15);
 Function: Produces a new HMMER::Results module which has
           been trimmed at the cutoff.
 Returns : a Bio::Tools::HMMER::Results module
 Args    : sequence cutoff and domain cutoff. in bits score
           if you want one cutoff, simply use same number both places

=cut

sub filter_on_cutoff {
    my $self = shift;
    my $seqthr = shift;
    my $domthr = shift;
    my ($new,$seq,$unit,@array,@narray);

    if( !defined $domthr ) {
       $self->throw("hmmresults filter on cutoff needs two arguments");
    }

    $new = Bio::Tools::HMMER::Results->new(-type => $self->parsetype);

    foreach $seq ( $self->each_Set()) {
	next if( $seq->bits() < $seqthr );
	$new->add_Set($seq);
	foreach $unit ( $seq->each_Domain() ) {
	    next if( $unit->bits() < $domthr );
	    $new->add_Domain($unit);
	}
    }
    $new;
}

=head2 write_ascii_out

 Title   : write_ascii_out
 Usage   : $res->write_ascii_out(\*STDOUT)
 Function: writes as
           seq seq_start seq_end model-acc model_start model_end model_name
 Returns :
 Args    :

  FIXME: Now that we have no modelacc, this is probably a bad thing.

=cut

# writes as seq sstart send modelacc hstart hend modelname

sub write_ascii_out {
    my $self = shift;
    my $fh = shift;
    my ($unit,$seq);

    if( !defined $fh) {
	$fh = \*STDOUT;
    }


    foreach $seq ( $self->each_Set()) {
	foreach $unit ( $seq->each_Domain()) {
	    print $fh sprintf("%s %4d %4d %s %4d %4d %4.2f %4.2g %s\n",
			      $unit->seq_id(),$unit->start(),$unit->end(),
			      $unit->hmmacc,$unit->hstart,$unit->hend,
			      $unit->bits,$unit->evalue,$unit->hmmname);
	}
    }

}

=head2 write_GDF_bits

 Title   : write_GDF_bits
 Usage   : $res->write_GDF_bits(25,15,\*STDOUT)
 Function: writes GDF format with a sequence,domain threshold
 Returns :
 Args    :

=cut

sub write_GDF_bits {
    my $self = shift;
    my $seqt = shift;
    my $domt = shift;
    my $file = shift;
    my $seq;
    my $unit;
    my (@array,@narray);

    if( !defined $file ) {
	$self->throw("Attempting to use write_GDF_bits without passing in correct arguments!");
	return;
    }

    foreach $seq ( $self->each_Set()) {

	if( $seq->bits() < $seqt ) {
	    next;
	}

	foreach $unit ( $seq->each_Domain() ) {
	    if( $unit->bits() < $domt ) {
		next;
	    }
	    push(@array,$unit);
	}

    }

    @narray = sort { my ($aa,$bb,$st_a,$st_b);
		     $aa = $a->seq_id();
		     $bb = $b->seq_id();
		     if ( $aa eq $bb) {
			 $st_a = $a->start();
			 $st_b = $b->start();
			 return $st_a <=> $st_b;
			 }
		     else {
			 return $aa cmp $bb;
		     } } @array;

    foreach $unit ( @narray ) {
	print $file sprintf("%-24s\t%6d\t%6d\t%15s\t%.1f\t%g\n",$unit->get_nse(),$unit->start(),$unit->end(),$unit->seq_id(),$unit->bits(),$unit->evalue);
    }

}

sub write_scores_bits {
    my $self = shift;
    my $seqt = shift;
    my $domt = shift;
    my $file = shift;
    my $seq;
    my $unit;
    my (@array,@narray);

    if( !defined $file ) {
	$self->warn("Attempting to use write_scores_bits without passing in correct arguments!");
	return;
    }

    foreach $seq ( $self->eachHMMSequence()) {

	if( $seq->bits() < $seqt ) {
	    next;
	}

	foreach $unit ( $seq->eachHMMUnit() ) {
	    if( $unit->bits() < $domt ) {
		next;
	    }
	    push(@array,$unit);
	}

    }

    @narray = sort { my ($aa,$bb,$st_a,$st_b);
		     $aa = $a->bits();
		     $bb = $b->bits();
		     return $aa <=> $bb;
		     } @array;

    foreach $unit ( @narray ) {
	print $file sprintf("%4.2f %s\n",$unit->bits(),$unit->get_nse());
    }

}

sub write_GDF {
    my $self = shift;
    my $file = shift;
    my $unit;

    if( !defined $file ) {
	$file = \*STDOUT;
    }


    foreach $unit ( $self->eachHMMUnit() ) {
	print $file sprintf("%-24s\t%6d\t%6d\t%15s\t%.1f\t%g\n",$unit->get_nse(),$unit->start(),$unit->end(),$unit->seq_id(),$unit->bits(),$unit->evalue);
    }

}

sub highest_noise {
    my $self = shift;
    my $seqt = shift;
    my $domt = shift;
    my ($seq,$unit,$hseq,$hdom,$noiseseq,$noisedom);

    $hseq = $hdom = -100000;

    foreach $seq ( $self->eachHMMSequence()) {
	if( $seq->bits() < $seqt && $seq->bits() > $hseq  ) {
	    $hseq = $seq->bits();
	    $noiseseq = $seq;
	}
	foreach $unit ( $seq->eachHMMUnit() ) {
	    if( (($seq->bits() < $seqt) || ($seq->bits() > $seqt && $unit->bits < $domt)) && $unit->bits() > $hdom ) {
		$hdom  = $unit->bits();
		$noisedom = $unit;
	    }
	}
    }


    return ($noiseseq,$noisedom);

}


sub lowest_true {
    my $self = shift;
    my $seqt = shift;
    my $domt = shift;
    my ($seq,$unit,$lowseq,$lowdom,$trueseq,$truedom);

    if( ! defined $domt ) {
	$self->warn("lowest true needs at least a domain threshold cut-off");
	return (0,0);
    }

    $lowseq = $lowdom = 100000;

    foreach $seq ( $self->eachHMMSequence()) {

	if( $seq->bits() >= $seqt && $seq->bits() < $lowseq  ) {
	    $lowseq = $seq->bits();
	    $trueseq = $seq;
	}
	if( $seq->bits() < $seqt ) {
	    next;
	}

	foreach $unit ( $seq->eachHMMUnit() ) {
	    if( $unit->bits() >= $domt && $unit->bits() < $lowdom ) {
		$lowdom  = $unit->bits();
		$truedom = $unit;
	    }
	}
    }


    return ($trueseq,$truedom);

}



=head2 add_Set

 Title   : add_Set
 Usage   : Mainly internal function
 Function:
 Returns :
 Args    :


=cut

sub add_Set {
    my $self = shift;
    my $seq  = shift;
    my $name;

    $name = $seq->name();

    if( exists $self->{'seq'}->{$name} ) {
	$self->throw("You alredy have $name in HMMResults!");
    }
    $self->{'seq'}->{$name} = $seq;
}


=head2 each_Set

 Title   : each_Set
 Usage   :
 Function:
 Returns :
 Args    :


=cut

sub each_Set {
    my $self = shift;
    my (@array,$name);


    foreach $name ( keys %{$self->{'seq'}} ) {
	push(@array,$self->{'seq'}->{$name});
    }
    return @array;
}


=head2 get_Set

 Title   : get_Set
 Usage   : $set = $res->get_Set('sequence-name');
 Function: returns the Set for a particular sequence
 Returns : a HMMER::Set object
 Args    : name of the sequence


=cut

sub get_Set {
    my $self = shift;
    my $name = shift;

    return $self->{'seq'}->{$name};
}


=head2 _parse_hmmpfam

 Title   : _parse_hmmpfam
 Usage   : $res->_parse_hmmpfam($filehandle)
 Function:
 Returns :
 Args    :


=cut

sub _parse_hmmpfam {
    my $self = shift;
    my $file = shift;

    my ($id,$sqfrom,$sqto,$hmmf,$hmmt,$sc,$ev,
	$unit,$nd,$seq,$name,$seqname,$from,
	$to,%hash,%acc,$acc);
    my $count = 0;

    while(<$file>) {
        if( /^HMM file:\s+(\S+)/ ) { $self->hmmfile($1); next; }
	elsif( /^Sequence file:\s+(\S+)/ ) { $self->seqfile($1); next }
	elsif( /^Query(\s+sequence)?:\s+(\S+)/ ) {

	    $seqname = $2;

	    $seq     = Bio::Tools::HMMER::Set->new();

	    $seq ->name($seqname);
	    $self->add_Set($seq);
	    %hash = ();

	    while(<$file>){

		if( /Accession:\s+(\S+)/ ) { $seq->accession($1); next }
		elsif( s/^Description:\s+// ) { chomp; $seq->desc($_); next }
		/^Parsed for domains/ && last;

		# This is to parse out the accession numbers in old Pfam format.
		# now not support due to changes in HMMER.

		if( (($id,$acc, $sc, $ev, $nd) = /^\s*(\S+)\s+(\S+).+?\s(\S+)\s+(\S+)\s+(\d+)\s*$/)) {
		    $hash{$id} = $sc; # we need this for the sequence
		                      # core of the domains below!
		    $acc {$id} = $acc;

		    # this is the more common parsing routine
		} elsif ( (($id,$sc, $ev, $nd) =
			   /^\s*(\S+).+?\s(\S+)\s+(\S+)\s+(\d+)\s*$/) ) {

		    $hash{$id} = $sc; # we need this for the
		                      # sequence score of hte domains below!

		}
	    }

	    while(<$file>) {
		/^Align/ && last;
		m{^//} && last;
		# this is meant to match

		#Sequence Domain  seq-f seq-t    hmm-f hmm-t      score  E-value
		#-------- ------- ----- -----    ----- -----      -----  -------
		#PF00621    1/1     198   372 ..     1   207 []   281.6    1e-80

		if( (($id, $sqfrom, $sqto, $hmmf,$hmmt,$sc, $ev) =
		     /(\S+)\s+\S+\s+(\d+)\s+(\d+).+?(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)\s*$/)) {
		    $unit = Bio::Tools::HMMER::Domain->new();
		    $unit->seq_id  ($seqname);
		    $unit->hmmname  ($id);
		    $unit->start    ($sqfrom);
		    $unit->end      ($sqto);
		    $unit->hstart($hmmf);
		    $unit->hend  ($hmmt);
		    $unit->bits     ($sc);
		    $unit->evalue   ($ev);

		    if( !exists($hash{$id}) ) {
			$self->throw("HMMResults parsing error in hmmpfam for $id - can't find sequecne score");
		    }

		    $unit->seqbits($hash{$id});

		    if( defined $acc{$id} ) {
			$unit->hmmacc($acc{$id});
		    }

		    # this should find it's own sequence!
		    $self->add_Domain($unit);
		}
	    }
	    if( m{^//} ) { next; }

	    $_ = <$file>;
	    # parses alignment lines. Icky as we have to break on the same line
	    # that we need to read to place the alignment lines with the unit.

	    while(1) {
		(!defined $_ || m{^//}) && last;

		# matches:
		# PF00621: domain 1 of 1, from 198 to 372
		if( /^\s*(\S+):.*from\s+(\d+)\s+to\s+(\d+)/ ) {

		    $name = $1;
		    $from = $2;
		    $to   = $3;

		    # find the HMMUnit which this alignment is from

		    $unit = $self->get_unit_nse($seqname,$name,$from,$to);
		    if( !defined $unit ) {
			$self->warn("Could not find $name $from $to unit even though I am reading it in. ugh!");
			$_ = <$file>;
			next;
		    }
		    while(<$file>) {
			m{^//} && last;
			/^\s*\S+:.*from\s+\d+\s+to\s+\d+/ && last;
			$unit->add_alignment_line($_);
		    }
		} else {
		    $_ = <$file>;
		}
	    }

	    # back to main 'Query:' loop
	}
    }
}

# mainly internal function

sub get_unit_nse {
    my $self    = shift;
    my $seqname = shift;
    my $domname = shift;
    my $start   = shift;
    my $end     = shift;

    my($seq,$unit);

    $seq = $self->get_Set($seqname);

    if( !defined $seq ) {
	$self->throw("Could not get sequence name $seqname - so can't get its unit");
    }

    foreach $unit ( $seq->each_Domain() ) {
	if( $unit->hmmname() eq $domname && $unit->start() == $start &&  $unit->end() == $end ) {
	    return $unit;
	}
    }

    return;
}


=head2 _parse_hmmsearch

 Title   : _parse_hmmsearch
 Usage   : $res->_parse_hmmsearch($filehandle)
 Function:
 Returns :
 Args    :


=cut

sub _parse_hmmsearch {
    my $self = shift;
    my $file = shift;
    my ($id,$sqfrom,$sqto,$sc,$ev,$unit,$nd,$seq,$hmmf,$hmmt,
	$hmmfname,$hmmacc, $hmmid, %seqh);
    my $count = 0;

    while(<$file>) {
        /^HMM file:\s+(\S+)/ and do { $self->hmmfile($1); $hmmfname = $1 };
	/^Accession:\s+(\S+)/ and do { $hmmacc = $1 };
	/^Query HMM:\s+(\S+)/ and do { $hmmid = $1 };
	/^Sequence database:\s+(\S+)/ and do { $self->seqfile($1) };
        /^Scores for complete sequences/ && last;
    }

    $hmmfname = "given" if not $hmmfname;

    while(<$file>) {
	/^Parsed for domains/ && last;
	if( (($id, $sc, $ev, $nd) = /(\S+).+?\s(\S+)\s+(\S+)\s+(\d+)\s*$/)) {
	    $seq = Bio::Tools::HMMER::Set->new();
	    $seq->name($id);
	    $seq->bits($sc);
	    $seqh{$id} = $sc;
	    $seq->evalue($ev);
	    $self->add_Set($seq);
	    $seq->accession($hmmacc);
	}
    }

    while(<$file>) {
	/^Alignments of top-scoring domains/ && last;
	if( (($id, $sqfrom, $sqto, $hmmf, $hmmt, $sc, $ev) = /(\S+)\s+\S+\s+(\d+)\s+(\d+).+?(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)\s*$/)) {
	    $unit = Bio::Tools::HMMER::Domain->new();

	    $unit->seq_id($id);
	    $unit->hmmname($hmmfname);
	    $unit->start($sqfrom);
	    $unit->end($sqto);
	    $unit->bits($sc);
	    $unit->hstart($hmmf);
	    $unit->hend($hmmt);
	    $unit->evalue($ev);
	    $unit->seqbits($seqh{$id});
	    $self->add_Domain($unit);
	    $count++;
	}
    }

    $_ = <$file>;

    ## Recognize and store domain alignments

    while(1) {
	if( !defined $_ ) {
	    last;
	}
        /^Histogram of all scores/ && last;

        # matches:
        # PF00621: domain 1 of 1, from 198 to 372
        if( /^\s*(\S+):.*from\s+(\d+)\s+to\s+(\d+)/ ) {
            my $name = $1;
            my $from = $2;
            my $to   = $3;

            # find the HMMUnit which this alignment is from
            $unit = $self->get_unit_nse($name,$hmmfname,$from,$to);

            if( !defined $unit ) {
                $self->warn("Could not find $name $from $to unit even though I am reading it in. ugh!");
                next;
            }
            while(<$file>) {
                /^Histogram of all scores/ && last;
                /^\s*\S+:.*from\s+\d+\s+to\s+\d+/ && last;
                $unit->add_alignment_line($_);
            }
        }
        else {
            $_ = <$file>;
        }
    }

    return $count;
}

=head2 parsetype

 Title   : parsetype
 Usage   : $obj->parsetype($newval)
 Function:
 Returns : value of parsetype
 Args    : newvalue (optional)


=cut

sub parsetype{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_parsetype'} = $value;
    }
    return $self->{'_parsetype'};
}

1;  # says use was ok
__END__


