#---------------------------------------------------------

=head1 NAME

Bio::Matrix::PSM::IO::meme - PSM meme parser implementation

=head1 SYNOPSIS

See Bio::Matrix::PSM::IO for detailed documentation on how to use PSM parsers

=head1 DESCRIPTION

Parser for meme.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

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

=cut


# Let the code begin...
package Bio::Matrix::PSM::IO::meme;
use Bio::Matrix::PSM::InstanceSite;
use Bio::Matrix::PSM::SiteMatrix;
use Bio::Matrix::PSM::Psm;
use vars qw(@HEADER);
use strict;

use base qw(Bio::Matrix::PSM::PsmHeader Bio::Matrix::PSM::IO);

@Bio::Matrix::PSM::IO::meme::HEADER = qw(e_val sites IC width);

=head2 new

 Title   : new
 Usage   : my $psmIO =  Bio::Matrix::PSM::IO->new(-format=>'meme', 
						 -file=>$file);
 Function: Associates a file with the appropriate parser
 Throws  : Throws if the file passed is in HTML format or 
           if the MEME header cannot be found.
 Example :
 Args    : hash
 Returns : "Bio::Matrix::PSM::$format"->new(@args);

=cut

sub new {
    my($class, @args)=@_;
    my $self = $class->SUPER::new(@args);
    my ($file)=$self->_rearrange(['FILE'], @args);
    my ($query,$tr1)=split(/\./,$file,2);
    $self->{file} = $file;
    $self->{query}= $query;
    $self->{end}  = 0;
    $self->{_strand}=0; #This we'll need to see if revcom option is used
    $self->_initialize_io(@args) || warn "Did you intend to use STDIN?"; #Read only for now
    #Skip header
    my $line;
    while (my $line=$self->_readline) {
	$self->throw('Cannot parse HTML, please use text output\n') if ($line=~/<HEAD>/); #Should start parsing HTML output, not a bug deal
	chomp($line);
	if ($line=~"^ALPHABET") {
	    $self=_parse_coordinates($self);
	    last;
	}
	push @{$self->{unstructured}},$line unless (($line=~/\*{10,}/) || ($line eq ''));
    }
    $self->_initialize;
    return $self;
}

=head2 _parse_coordinates

 Title   : _parse_coordinates
 Usage   :
 Function:
 Throws  :
 Example : Internal stuff
 Returns :
 Args    :

=cut

sub _parse_coordinates {
    my $self=shift;
    $self->_readline;
    $self->_readline;
    my $line=$self->_readline;
    while ($line !~ /^\*{10,}/ ) {
	chomp $line;
	$line =~ s/\s+/,/g;
	my ($id1,$w1,$l1,$id2,$w2,$l2)=split(/,/,$line);
	push @{$self->{hid}},$id1;
	$self->{weight}->{$id1}=$w1;
	$self->{length}->{$id1}=$l1;
	if ($id2) {
	    push @{$self->{hid}},$id2;
	    $self->{weight}->{$id2}=$w2;
	    $self->{length}->{$id2}=$l2;
	}
	$line=$self->_readline;
    }
    return $self;
}

=head2 header

 Title   : header
 Usage   :  my %header=$psmIO->header;
 Function:  Returns the header for the MEME file
 Throws  :
 Example : Fetching all the sequences included in the MEME analysis, 
           being parsed
           my %header=$psmIO->header;
            foreach my $seqid (@{$header{instances}}) {
               my $seq=$db->get_Seq_by_acc($id);
               #Do something with the sequence
            }
            where $db might be Bio::DB:GenBank object, see
 Returns : Hash with three keys: instances, weights and lengths, which
           should be self-explenatory. Each value is an array
           reference. Each array element corresponds to the same
           element in the other two arrays. So $header{instances}->[$i]
           will refer to the same sequence in the motif file as
           $header{weights}->[$i] and $header{lengths}->[$i]
 Args    :  none
 Notes   :  OBSOLETE!

=cut

sub header {
    my $self=shift;
    my @instances=@{$self->{_inst_name}};
    my @weights=@{$self->{_inst_weight}};
    my @lengths=@{$self->{_inst_coord}};
    return (instances=>\@instances,weights=>\@weights,lengths=>\@lengths);
}

=head2 next_psm

 Title   : next_psm
 Usage   : my $psm=$psmIO->next_psm();
 Function: Reads the next PSM from the input file, associated with this object
 Throws  : Throws if the format is inconsistent with the rules for MEME 3.0.4:
            no SUMMARY Section present or some keywords are missing/altered.
 Example :
 Returns : Bio::Matrix::PSM::Psm object
 Args    : none

=cut

sub next_psm {
    #Parses the next prediction and returns a psm objects
    my $self=shift;
    return if ($self->{end});
    my ($endm,$line,$instances,$tr,$width,$motif_id,$sites,$e_val,$id,$ic,$lA,$lC,$lG,$lT);
    while (defined( $line = $self->_readline) ) {
#Check if revcom is enabled, not very original check....
  $self->{_strand}=1 if (($line=~/^Sequence name/) && ($line=~/Strand/));
	if ($line=~ m/\sSite\s/) {
	    $instances= $self->_parseInstance;
	}
	#Here starts the next motif
	if ( ($line=~/width/) && ($line=~/sites/)) {
	    chomp($line);
	    $line=~s/[\t\s=]+/,/g;
	    $line=~s/\t/,/g;
	    #Parsing the general information for this prediction
	    ($tr,$motif_id,$tr,$width,$tr,$sites,
	     $tr,$tr,$tr,$e_val)=split(/,/,$line);
	    $self->{id}=$self->{query} . $motif_id;
	}
	if ($line =~ /content/i) {
	    $line=$self->_readline;
	    chomp($line);
	    $line=~s/[\)\(]//g;
	    ($ic)=split(/\s/,$line);
	}
        #Last info-prob matrix data
	if ($line=~/position-specific\s+scoring matrix/) {
		($lA,$lC,$lG,$lT)=_parse_logs($self);
	}
	if ($line=~/^letter-probability\smatrix/) {
	    my %matrix_dat=$self->_parseMatrix($motif_id);
	    my $psm= Bio::Matrix::PSM::Psm->new(%matrix_dat, 
					       -instances=>$instances, 
					       -e_val=>$e_val,
					       -IC=>$ic, 
					       -width=>$width, 
					       -sites=>$sites,
						   -lA=>$lA,
						   -lC=>$lC,
						   -lG=>$lG,
						   -lT=>$lT,
						   );
	    return $psm;
	}
	if ($line=~"SUMMARY OF MOTIFS") {
	    $self->{end}=1;
	    return;
	}
	$endm=1 if ($line=~/^Time\s/); 
    }
	if ($endm) { #End of file found, end of current motif too, but not all predictions were made as requested (No summary)
	    $self->{end}=1;
            warn "This MEME analysis was terminated prematurely, you may have less motifs than you requested\n";
	    return;
	}
    $self->throw("Wrong format\n"); # Multiple keywords not found, probably wrong format
}

=head2 _parseMatrix

 Title   : _parseMatrix
 Usage   :
 Function: Parses the next site matrix information in the meme file
 Throws  :
 Example :  Internal stuff
 Returns :  hash as for constructing a SiteMatrix object (see SiteMatrixI)
 Args    :  string

=cut

sub _parseMatrix {
    my ($self,$id)=@_;
    my (@pA,@pC,@pG,@pT);
    my $i=0;
    my $line = $self->_readline;
    #Most important part- the probability matrix
    do {
	chomp $line;
	last if ($line eq '');
  $line=~s/^\s+//;
	$line=~s/\s+/,/g;
	($pA[$i],$pC[$i],$pG[$i],$pT[$i])=split(/,/,$line);
	$i++;
	$line=$self->_readline;
    } until $line =~ /\-{10,}/;
    return (-pA=>\@pA,-pC=>\@pC,-pG=>\@pG,-pT=>\@pT,-id=>$id);
}

=head2 _parse_logs

 Title   : _parse_logs
 Usage   :
 Function: Parses the next site matrix log values in the meme file
 Throws  :
 Example :  Internal stuff
 Returns :  array of array refs
 Args    :  string

=cut

sub _parse_logs {
    my $self=shift;
    my (@lA,@lC,@lG,@lT);
    my $i=0;
    $self->_readline;   $self->_readline;
    my $line = $self->_readline;
    #Most important part- the probability matrix
    do {
	chomp $line;
	last if ($line eq '');
  $line=~s/^\s+//;
	$line=~s/\s+/,/g;
	($lA[$i],$lC[$i],$lG[$i],$lT[$i])=split(/,/,$line);
	$i++;
	$line=$self->_readline;
    } until $line =~ /\-{10,}/;
    
    return (\@lA,\@lC,\@lG,\@lT);
}

=head2 _parseInstance

 Title   : _parseInstance
 Usage   :
 Function:  Parses the next sites instances from the meme file
 Throws  :
 Example :  Internal stuff
 Returns :  Bio::Matrix::PSM::InstanceSite object
 Args    :  none

=cut

sub _parseInstance {
    my $self = shift;
    my $i=0;
    $self->_readline;
    my ($line,@instance);
    while (defined($line=$self->_readline) ) {
	last if ($line =~ /\-{5}/ );
	chomp($line);
	my @comp=split(/\s+/,$line);
	my ($id,$start,$score,$strand,$s1,$s2,$s3);
	if ( $self->{_strand}) {
	    ($id,$strand,$start,$score,$s1,$s2,$s3)=@comp;
	} else {
	    ($id,$start,$score,$s1,$s2,$s3)=@comp;
	    $strand=1;
	}
  	my $seq= $s1.$s2.$s3;
	if ($seq =~ /[^ACGTacgtNnXx\-\.]/) {
            my $col=$#comp;
	    $self->throw("I have not been able to parse the correct instance sequence: $seq, $col columns\n");
	}
	my $sid = $self->{id} . '@' . $id;
	$instance[$i] = Bio::Matrix::PSM::InstanceSite->new
	    (-mid      => $self->{id}, 
	     -start    => $start, 
	     -score    => $score,
	     -seq      => $seq, 
	     -strand   => $strand,
	     -accession_number => $id, 
	     -primary_id => $sid, 
	     -desc => 'Bioperl MEME parser object' );
	$i++;
    }
    $self->{instances} = \@instance;
    return \@instance;
}

				
			

1;
