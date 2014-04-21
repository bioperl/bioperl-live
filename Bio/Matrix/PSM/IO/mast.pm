
=head1 NAME

Bio::Matrix::PSM::IO::mast - PSM mast parser implementation

=head1 SYNOPSIS

See Bio::Matrix::PSM::IO for detailed documentation on how to 
use PSM parsers

=head1 DESCRIPTION

Parser for mast. This driver unlike meme or transfac for example is
dedicated more to PSM sequence matches, than to PSM themselves.

=head1 TO DO

Section III should be parsed too, otherwise no real sequence is
available, so we supply 'NNNNN....' as a seq which is not right.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
to one of the Bioperl mailing lists. Your participation is much appreciated.

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

=head1 AUTHOR - Stefan Kirov

Email skirov@utk.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::Matrix::PSM::IO::mast;
use Bio::Matrix::PSM::InstanceSite;
use Bio::Matrix::PSM::Psm;
use Bio::Root::Root;
use strict;

use base qw(Bio::Matrix::PSM::PsmHeader Bio::Matrix::PSM::IO);

=head2 new

 Title   : new
 Usage   : my $psmIO =  Bio::Matrix::PSM::IO->new(-format=>'mast', 
						                               -file=>$file);
 Function: Associates a file with the appropriate parser
 Throws  : Throws if the file passed is in HTML format or if 
           some criteria for the file
           format are not met.
 Example :
 Returns : psm object, associated with a file with matrix file
 Args    : hash
 return  : "Bio::Matrix::PSM::$format"->new(@args);

=cut


sub new {
    my($class, @args)=@_;
    my $self = $class->SUPER::new(@args);
    my (%instances,@header,$n);
    my ($file)=$self->_rearrange(['FILE'], @args);
    $self->{file} = $file;
    $self->{_factor}=1;
    $self->_initialize_io(@args) || warn "Did you intend to use STDIN?"; #Read only for now
    $self->{_end}=0;
    undef $self->{hid};
    return $self if ($file=~/^>/);#Just writing
    my $buf=$self->_readline;
	$self->throw('Cannot parse HTML format yet') if ($buf =~/^<HTML>/); 
    # this should probably be moved to its own function
    while ( defined($buf=$self->_readline)) {
	chomp($buf);
	if ($buf=~/DATABASE AND MOTIFS/) {
		while ($buf=$self->_readline) {
			if ($buf=~/DATABASE/) {
					$buf=~s/^[\s\t]+//;
					chomp $buf;
					($n,$self->{_dbname},$self->{_dbtype})=split(/\s/,$buf);
					$self->{_dbtype}=~s/[\(\)]//g;
			}
			if ($buf=~/MOTIFS/) {
					$buf=~s/^[\s\t]+//;
					chomp $buf;
					($n,$self->{_mrsc},$self->{_msrctype})=split(/\s/,$buf);
					$self->{_msrctype}=~s/[\(\)]//g;
					last;
			}
		}
		if ($self->{_msrctype} ne $self->{_dbtype}) {#Assume we have protein motifs, nuc DB (not handling opp.)
			$self->{_factor}=3;
			$self->{_mixquery}=1;
		}
	}
	if ($buf=~m/MOTIF WIDTH BEST POSSIBLE MATCH/) {
	    $self->_readline;
	    while (defined($buf=$self->_readline)) {
		last if ($buf!~/\w/);
		$buf=~s/\t+//g;
		$buf=~s/^\s+//g;
		my ($id,$width,$seq)=split(/\s+/,$buf);
		push @{$self->{hid}},$id;
		$self->{length}->{$id}=$width;
		$self->{seq}->{$id}=$seq;
	    }
	    next;
	}
	if ($buf=~m/section i:/i) {
	    $self->_readline;
	    $self->_readline;
	    $self->_readline;
	    %instances=_get_genes($self);
	    $self->{instances}=\%instances;
      	if (!(%instances)) {
        	$self->warn ("Your MAST analysis did not find any matches satisfying the current threshold.\nSee MAST documentation for more information.\n");
        	return $self; #The header might be useful so we return the object, not undef
      	}
	    next;
	}
	if ($buf=~m/section ii:/i) {
	    $self->_readline;
	    $self->_readline;
	    $self->_readline;
	    last;
	}
	$buf=~s/[\t+\s+]/ /g;
	push @header,$buf unless (($buf=~/\*{10,}/)||($buf!~/\w/));
    }
    $self->throw('Could not read Section I, probably wrong format, make sure it is not HTML, giving up...') if !(%instances);
    $self->warn( "This file might be an unreadable version, proceed with caution!\n") if (!grep(/\s+MAST\s+version\s+3/,@header));

    $self->{unstructured} = \@header;
    $self->_initialize;
    return $self;
}


# Get the file header and put store it as a hash, which later we'll use to create
# the header for each Psm. See Bio::Matrix::PSM::PsmI for header function.
sub _get_genes {
	my $self=shift;
	my %llid;
	my $ok=0;
	my $i=0;
	my %instances;
	while (my $line=$self->_readline) {
		last if ($line=~/^[\s\t*]/); # Well, ids can be nearly anything...???
		chomp($line);
		$i++;
		next if ($line eq '');
		$line=~s/\s+/,/g;
		my ($id,$key,$eval,$len)=split(/,/,$line);
		unless ($len) {
			warn "Malformed data found: $line\n";
			next;
		}
		$instances{$id}=Bio::Matrix::PSM::InstanceSite->new(-id=>$id,
																			  -desc=>$key,
																			  -score=>$eval, 
																			  -width=>$len,
																			  -seq=>'ACGT');
	}
	return %instances;
}


=head2 next_psm

 Title   : next_psm
 Usage   : my $psm=$psmIO->next_psm();
 Function: Reads the next PSM from the input file, associated with this object
 Throws  : Throws if there ara format violations in the input file (checking is not
            very strict with all drivers).
 Example :
 Returns : Bio::Matrix::PSM::Psm object
 Args    : none

=cut


sub next_psm {
    my $self=shift;
    return if ($self->{_end}==1);
    my (@lmotifsm,%index,$eval,$scheme,$sid);
    %index= %{$self->{length}};
    my (@instances,%instances);
    my $line=$self->_readline;
    $line=~s/[\t\n]//;
    if ($line =~ /\*{10,}/) { #Endo of Section II if we do only section II
        $self->{_end}=1;
        return ;
    }
    do {
	if ($line!~/^\s/) {
	    ($sid,$eval,$scheme)=split(/\s+/,$line,3);
	}
	else
	{ $scheme .=$line; }
	$line=$self->_readline;
	$line=~s/[\t\n]//;
    } until ($line!~/^\s/);
    my $pos=1;
    $scheme=~s/\s+//g;
    $scheme=~s/\n//g;
    my @motifs=split(/_/,$scheme);
    while (@motifs) {
	my $next=shift(@motifs);
	if (!($next=~/\D/)) {
	    last if (!@motifs);
	    $pos+=$next;
	    next;
	}
        my $id=$next;
	my $score= $id=~m/\[/ ? 'strong' : 'weak' ;
	my $frame;
	my $strand = $id =~ m/\-\d/ ? -1 : 1 ;
	if ($self->{_mixquery}) {
		$frame = 0 if $id =~ m/\d+a/ ;
		$frame = 1 if $id =~ m/\d+b/ ;
		$frame = 2 if $id =~ m/\d+c/ ;
	}
	$id=~s/\D+//g;

	my @s;
	my $width=$index{$id};
    #We don't know the sequence, but we know the length
	my $seq='N' x ($width*$self->{_factor}); #Future version will have to parse Section tree nad get the real seq
	my $instance=Bio::Matrix::PSM::InstanceSite->new 
	    ( -id=>"$id\@$sid", 
	      -mid=>$id, 
	      -accession_number=>$sid,
	      -desc=>"Motif $id occurrance in $sid",
	      -score=>$score, 
	      -seq=>$seq,
		  -alphabet => 'dna', 
	      -start=>$pos,
	      -strand=>$strand);
	  $instance->frame($frame) if ($self->{_mixquery});
	push @instances,$instance;
	$pos+=$index{$id}*$self->{_factor};
    }
    my $psm= Bio::Matrix::PSM::Psm->new(-instances=> \@instances, 
					-e_val    => $eval, 
					-id       => $sid);
    $self->_pushback($line);
    return $psm;
}


=head2 write_psm

 Title   : write_psm
 Usage   : #Get SiteMatrix object somehow (see Bio::Matrix::PSM::SiteMatrix)
            my $matrix=$psmin->next_matrix;
            #Create the stream
            my $psmio=new(-file=>">psms.mast",-format=>'mast');
            $psmio->write_psm($matrix);
            #Will warn if only PFM data is contained in $matrix, recalculate the PWM
            #based on normal distribution (A=>0.25, C=>0.25, etc)
 Function: writes pwm in mast format
 Throws  :
 Example : 
 Args    : SiteMatrix object
 Returns : 

=cut

sub write_psm {
    my ($self,$matrix)=@_;
#    my $idline=">". $matrix->id . "\n";
    my $w=$matrix->width;
    my $header="ALPHABET= ACGT\nlog-odds matrix: alength= 4 w= $w\n";
    $self->_print($header);
    unless ($matrix->get_logs_array('A')) {
        warn "No log-odds data, available, using normal distribution to recalculate the PWM";
        $matrix->calc_weight({A=>0.25, C=>0.25, G=>0.25,T=>0.25});
    }
    while (my %h=$matrix->next_pos) {
	$self->_print (join("\t",$h{lA},$h{lC},$h{lG},$h{lT},"\n"));
    }
}

1;
