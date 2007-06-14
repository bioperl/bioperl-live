# $Id$
# Bioperl module Bio::Tools::BPlite::Iteration
#	based closely on the Bio::Tools::BPlite modules
#	Ian Korf (ikorf@sapiens.wustl.edu, http://sapiens.wustl.edu/~ikorf), 
#	Lorenz Pollak (lorenz@ist.org, bioperl port)
#
# Copyright Peter Schattner
#
# You may distribute this module under the same terms as perl itself
# _history
# October 20, 2000
# POD documentation - main docs before the code
#
# Added to get a simple_align object for a psiblast run with the -m 6 flag /AE
# 

=head1 NAME

Bio::Tools::BPlite::Iteration - object for parsing single iteration
of a PSIBLAST report

=head1 SYNOPSIS

   use Bio::Tools::BPpsilite;

   open my $FH, "t/psiblastreport.out";
   $report = Bio::Tools::BPpsilite->new(-fh=>\*FH);

   # determine number of iterations executed by psiblast
   $total_iterations = $report->number_of_iterations;
   $last_iteration = $report->round($total_iterations);

   # Process only hits found in last iteration ...
   $oldhitarray_ref = $last_iteration->oldhits;
   HIT: while($sbjct = $last_iteration->nextSbjct) {
       $id = $sbjct->name;
       $is_old =  grep  { $_ eq $id } @$oldhitarray_ref;
       if ($is_old ){next HIT;}
   #  do something with new hit...
   }

=head2 ALIGNMENTS

  # This assumed that you have $db pointing to a database, $out to an output file
  # $slxdir to a directory and $psiout    
  # note the alignments can only be obtained if the flag "-m 6" is run.
  # It might also be necessary to use the flag -v to get all alignments
  # 
    my @psiparams = ('database' => $db , 'output' => $out, 'j' => 3, 'm' => 6,
		     'h' => 1.e-3 , 'F' => 'T' , 'Q' => $psiout ); 
    my $factory = Bio::Tools::Run::StandAloneBlast->new(@psiparams);
    my $report = $factory->blastpgp($seq);
    my $total_iterations = $report->number_of_iterations();
    my $last_iteration = $report->round($total_iterations);
    my $align=$last_iteration->Align;
    my $slxfile=$slxdir.$id.".slx";
    my $slx = Bio::AlignIO->new('-format' => 'selex','-file' => ">".$slxfile );
    $slx->write_aln($align);

=head1 DESCRIPTION

See the documentation for BPpsilite.pm for a description of the
Iteration.pm module.

=head1 AUTHORS - Peter Schattner

Email: schattner@alum.mit.edu

=head1 CONTRIBUTORS

Jason Stajich, jason@cgt.mc.duke.edu

=head1 ACKNOWLEDGEMENTS

Based on work of:
Ian Korf (ikorf@sapiens.wustl.edu, http://sapiens.wustl.edu/~ikorf), 
Lorenz Pollak (lorenz@ist.org, bioperl port)

=head1 COPYRIGHT

BPlite.pm is copyright (C) 1999 by Ian Korf. 

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut

package Bio::Tools::BPlite::Iteration;

use strict;
use Bio::Tools::BPlite; #
use Bio::Tools::BPlite::Sbjct;

use base qw(Bio::Root::Root);

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    ($self->{'PARENT'},$self->{'ROUND'}) =
	$self->_rearrange([qw(PARENT
			      ROUND
			      )],@args);
    
    $self->{'QUERY'} = $self->{'PARENT'}->{'QUERY'};
    $self->{'LENGTH'} = $self->{'PARENT'}->{'LENGTH'};

    if($self->_parseHeader) {$self->{'REPORT_DONE'} = 0} # there are alignments
    else                    {$self->{'REPORT_DONE'} = 1} # empty report
  
    return $self; # success - we hope!
}

=head2 query

 Title    : query
 Usage    : $query = $obj->query();
 Function : returns the query object
 Example  :
 Returns  : query object
 Args     :

=cut

sub query    {shift->{'QUERY'}}

=head2 qlength

 Title    : qlength
 Usage    : $len = $obj->qlength();
 Returns  : length of query
 Args     : none

=cut

sub qlength  {shift->{'LENGTH'}}

=head2 newhits

 Title    :  newhits
 Usage    : $newhits = $obj->newhits();
 Returns  : reference to an array listing all the hits 
            from the current iteration which were not identified 
            in the previous iteration
 Args     : none

=cut

sub  newhits  {shift->{'NEWHITS'}}

=head2 oldhits

 Title    :  oldhits
 Usage    : $oldhits = $obj->oldhits();
 Returns  : reference to an array listing all the hits from 
            the current iteration which were identified and 
            above threshold in the previous iteration
 Args     : none

=cut

sub  oldhits  {shift->{'OLDHITS'}}


=head2 nextSbjct

 Title    : nextSbjct
 Usage    : $sbjct = $obj->nextSbjct();
 Function : Method of iterating through all the Sbjct retrieved
            from parsing the report 
#Example  : while ( my $sbjct = $obj->nextSbjct ) {}
 Returns  : next Sbjct object or undef if finished
 Args     :

=cut

sub nextSbjct {
    my ($self) = @_;
    $self->_fastForward or return;

    #######################
    # get all sbjct lines #
    #######################
    my $def = $self->_readline();

    while(defined ($_ = $self->_readline) )  {
	if    ($_ !~ /\w/)            {next}
	elsif ($_ =~ /Strand HSP/)    {next} # WU-BLAST non-data
	elsif ($_ =~ /^\s{0,2}Score/) {$self->_pushback( $_); last}
	elsif ($_ =~ /^(\d+) .* \d+$/) { # This is not correct at all
	    $self->_pushback($_); # 1: HSP does not work for -m 6 flag
	    $def = $1;		  # 2: length/name are incorrect     
	    my $length = undef;	  # 3: Names are repeated many times.
	    my $sbjct = Bio::Tools::BPlite::Sbjct->new('-name'=>$def,
						      '-length'=>$length,
						      '-parent'=>$self);
	    return $sbjct;
	}			# m-6 
	elsif ($_ =~ /^Parameters|^\s+Database:|^\s+Posted date:/) {
	    $self->_pushback( $_); 
	    last;
	} else {$def .= $_}
    }
    $def = '' unless defined $def;
    $def =~ s/\s+/ /g;
    $def =~ s/\s+$//g;
    $def =~ s/Length = ([\d,]+)$//g;
    my $length = $1;
    return 0 unless $def =~ /^>/;
    $def =~ s/^>//;

    ####################
    # the Sbjct object #
    ####################
    my $sbjct = Bio::Tools::BPlite::Sbjct->new('-name'=>$def,
					      '-length'=>$length,
					      '-parent'=>$self);
    return $sbjct;
}


# This is added by /AE

=head2 Align

 Title    : Align
 Usage    : $SimpleAlign = $obj->Align();
 Function : Method to obtain a simpleAlign object from psiblast
 Example  : $SimpleAlign = $obj->Align();
 Returns  : SimpleAlign object or undef if not found.
 BUG      : Only works if psiblast has been run with m 6 flag
 Args     :

=cut

sub Align {
    use Bio::SimpleAlign;
    my ($self) = @_;
    $self->_fastForward or return;
    my $lastline = $self->_readline();
    return unless $lastline =~ /^QUERY/; # If psiblast not run correctly
    my (%sequence,%first,%last,$num);

    if ( $lastline =~ /^QUERY\s+(\d*)\s*([-\w]+)\s*(\d*)\s*$/){
	my $name='QUERY';
	my $start=$1; 
	my $seq=$2; 
	my $stop=$3; 
	$seq =~ s/-/\./g; 
	$start =~ s/ //g; 
	$stop =~ s/ //g; 
	$sequence{$name} .= $seq; 
	if ($first{$name} eq ''){$first{$name}=$start;} 
	if ($stop ne ''){$last{$name}=$stop;} 
#     print "FOUND:\t$seq\t$start\t$stop\n"; 
	$num=0;
    } 
    while(defined($_ = $self->_readline()) ){
	chomp($_);
	if ( $_ =~ /^QUERY\s+(\d+)\s*([\-A-Z]+)\s*(\+)\s*$/){
	    my $name='QUERY';
	    my $start=$1; 
	    my $seq=$2; 
	    my $stop=$3; 
	    $seq   =~ s/-/\./g; 
	    $start =~ s/ //g; 
	    $stop  =~ s/ //g; 
	    $sequence{$name} .= $seq; 
	    if ($first{$name} eq '') { $first{$name} = $start;} 
	    if ($stop ne '') { $last{$name}=$stop;} 
	    $num=0;
	} elsif ( $_ =~ /^(\d+)\s+(\d+)\s*([\-A-Z]+)\s*(\d+)\s*$/ ){
	    my $name=$1.".".$num;
	    my $start=$2;
	    my $seq=$3;
	    my $stop=$4;
	    $seq =~ s/-/\./g;
	    $start =~ s/ //g;
	    $stop =~ s/ //g;
	    $sequence{$name} .= $seq;
	    if ($first{$name} eq ''){$first{$name}=$start;}
	    if ($stop ne ''){$last{$name}=$stop;}
	    $num++;
	} 
    } 
    my $align = Bio::SimpleAlign->new();
    my @keys=sort keys(%sequence);
    foreach my $name (@keys){
	my $nse = $name."/".$first{$name}."-".$last{$name};
	my $seqobj = Bio::LocatableSeq->new( -seq => $sequence{$name},
					     -id  => $name,
					     -name  => $nse,
					     -start  => $first{$name},
					     -end  => $last{$name}
					     );

	$align->add_seq($seqobj);
    }
    return $align;
}

# Start of internal subroutines.

sub _parseHeader {
  my ($self) = @_;
  my (@old_hits, @new_hits);

  my $newhits_true = ($self->{'ROUND'} < 2) ? 1  : 0 ;
  while(defined($_ = $self->_readline()) ) {
    if ($_ =~ /(\w\w|.*|\w+.*)\s\s+(\d+)\s+([-\.e\d]+)$/)    {
	my $id = $1;
	my $score= $2;	#not used currently
	my $evalue= $3; 	#not used currently
    	if ($newhits_true) { push ( @new_hits, $id);}
    	else { push (@old_hits, $id);}
    }
    elsif ($_ =~ /^Sequences not found previously/)  {$newhits_true = 1 ;}
# This is changed for "-m 6" option /AE
    elsif ($_ =~ /^>/ || $_ =~ /^QUERY/)
    {
	$self->_pushback($_);
	$self->{'OLDHITS'} = \@old_hits;
	$self->{'NEWHITS'} = \@new_hits;
	return 1;
    }
    elsif ($_ =~ /^Parameters|^\s+Database:|^\s*Results from round\s+(d+)/) {
      	$self->_pushback($_);
      	return 0; #  no sequences found in this iteration
    }
  }
  return 0; # no sequences found in this iteration
}

sub _fastForward {
  my ($self) = @_;
  return 0 if $self->{'REPORT_DONE'}; # empty report

  while(defined($_ = $self->_readline()) ) {
      if( $_ =~ /^>/ ||
	  $_ =~ /^QUERY|^\d+ .* \d+$/ ) { # Changed to also handle "-m 6" /AE
	  $self->_pushback($_);
	  return 1;
      }
#    print "FASTFORWARD",$_,"\n";
      if ($_ =~ /^>|^Parameters|^\s+Database:/) {
	  $self->_pushback($_);
	  return 1;
      }
  }
  $self->warn("Possible error (2) while parsing BLAST report!");
}


=head2 _readline

 Title   : _readline
 Usage   : $obj->_readline
 Function: Reads a line of input.

           Note that this method implicitely uses the value of $/ that is
           in effect when called.

           Note also that the current implementation does not handle pushed
           back input correctly unless the pushed back input ends with the
           value of $/.
 Example :
 Returns : 

=cut

sub _readline{
   my ($self) = @_;
   return $self->{'PARENT'}->_readline();
}

=head2 _pushback

 Title   : _pushback
 Usage   : $obj->_pushback($newvalue)
 Function: puts a line previously read with _readline back into a buffer
 Example :
 Returns :
 Args    : newvalue

=cut

sub _pushback {
   my ($self, $arg) = @_;   
   return $self->{'PARENT'}->_pushback($arg);    
}
1;
__END__
