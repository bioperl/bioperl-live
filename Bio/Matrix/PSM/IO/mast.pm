#---------------------------------------------------------
# $Id$

=head1 NAME

Bio::Matrix::PSM::IO::mast - PSM mast parser implementation

=head1 SYNOPSIS

See Bio::Matrix::PSM::IO for detailed documentation on how to use PSM parsers

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
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://bio.perl.org/MailList.html     - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

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
use Bio::Matrix::PSM::IO;
use Bio::Matrix::PSM::PsmHeader;
use Bio::Root::Root;
use strict;
use vars qw(@ISA);

@ISA=qw(Bio::Matrix::PSM::PsmHeader Bio::Root::Root Bio::Matrix::PSM::IO);

=head2 new

 Title   : new
 Usage   : my $psmIO =  new Bio::Matrix::PSM::IO(-format=>'mast', 
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
    my (%instances,@header);
    my ($file)=$self->_rearrange(['FILE'], @args);
    $self->{file} = $file;
    $self->_initialize_io(@args) || warn "Did you intend to use STDIN?"; #Read only for now
    $self->{_end}=0;
    undef $self->{hid};
    my $buf;
    # this should probably be moved to its own function
    while ( defined($buf=$self->_readline)) {
	chomp($buf);
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
        	$self->warn ("Your MAST analysis did not find any matches satisfying the current thershold.\nSee MAST documentation for more information.\n");
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


#Get the file header and put store it as a hash, which later we'll use to create
#the header for each Psm. See Bio::Matrix::PSM::PsmI for header function.
sub _get_genes {
    my $self=shift;
    my %llid;
    my $ok=0;
    my $i=0;
    my %instances;
    while (my $line=$self->_readline) {
	last if ($line=~/^\D{10,}/);
	chomp($line);
	$i++;
	next if ($line eq '');
	$line=~s/\s+/,/g;
	my ($id,$key,$eval,$len)=split(/,/,$line);
	$instances{$id}=new Bio::Matrix::PSM::InstanceSite ( -id=>$id,
							     -desc=>$key,-score=>$eval, -width=>$len,-seq=>'ACGT');
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
    return undef if ($self->{_end}==1);
    my (@lmotifsm,%index,$eval,$scheme,$sid);
    my $i=0;
    %index= %{$self->{length}};
    my (@instances,%instances);
    my $line=$self->_readline;
    $line=~s/[\t\n]//;
    if ($line =~ /\*{10,}/) { #Endo of Section II if we do only section II
        $self->{_end}=1;
        return undef;
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
    my $pos=0;
    $scheme=~s/\s+//g;
    $scheme=~s/\n//g;
    my @motifs=split(/_/,$scheme);
    $i++;
    while (@motifs) {
	my $next=shift(@motifs);
	if (!($next=~/\D/)) {
	    last if (!@motifs);
	    $pos+=$next;
	    next;
	}
        my $id=$next;
	my $score= $id=~m/\[/ ? 'strong' : 'weak' ;
	$id=~s/\D+//g;
	my @s;
	my $width=$index{$id};
	foreach (1..$width) {push @s,'N';} #We don't know the sequence, but we know the length
	my $seq=join('N',@s); #Future version will have to parse Section tree nad get the real seq
	my $instance=new Bio::Matrix::PSM::InstanceSite 
	    ( -id=>"$id\@$sid", 
	      -mid=>$id, 
	      -accession_number=>$sid,
	      -desc=>"Motif $id occurrance in $sid",
	      -score=>$score, 
	      -seq=>$seq,
		  -alphabet => 'dna', 
	      -start=>$pos);
	push @instances,$instance;
	$pos+=$index{$id};
    }
    my $psm= new Bio::Matrix::PSM::Psm (-instances=> \@instances, 
					-e_val    => $eval, 
					-id       => $sid);
    $self->_pushback($line);
    return $psm;
}


1;
