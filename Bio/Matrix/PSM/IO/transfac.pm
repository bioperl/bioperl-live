#---------------------------------------------------------
# $Id$

=head1 NAME

Bio::Matrix::PSM::transfac - PSM transfac parser

=head1 SYNOPSIS

See Bio::Matrix::PSM::IO for documentation

=head1 DESCRIPTION

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Stefan Kirov

Email skirov@utk.edu

=head1 APPENDIX

=cut


# Let the code begin...
package Bio::Matrix::PSM::IO::transfac;
use Bio::Matrix::PSM::Psm;
use Bio::Matrix::PSM::IO;
use Bio::Matrix::PSM::PsmHeader;
use Bio::Root::Root;
use vars qw(@ISA);
use strict;

@ISA=qw(Bio::Matrix::PSM::PsmHeader Bio::Root::Root Bio::Matrix::PSM::IO);

=head2 new

 Title   : new
 Usage   : my $psmIO =  new Bio::Matrix::PSM::IO(-format=>'transfac', 
						 -file=>$file);
 Function: Associates a file with the appropriate parser
 Throws  :
 Example :
 Args    :
 Returns : "Bio::Matrix::PSM::$format"->new(@args);

=cut

sub new {
    my ($class,@args)=@_;
    my $line;
    my $self = $class->SUPER::new(@args);
    my ($file)=$self->_rearrange(['file'], @args);
    $self->_initialize_io("<$file") || warn "Did you intend to use STDIN?"; #Read only for now
    #Remove header
    do {
	$line=$self->_readline;
	chomp $line;
	push @{$self->{unstructured}},$line if (length($line)>2); } until ($line =~ /^\/\//) || (!defined($line)); #Unstructured header
    $self->_initialize;
    return $self;
}
   

=head2 next_psm

 Title   : next_psm
 Usage   : my $psm=$psmIO->next_psm();
 Function: Reads the next PSM from the input file, associated with this object
 Throws  :
 Returns : Bio::Matrix::PSM::Psm object
 Args    : none

=cut

sub next_psm {
    my $self=shift;
    my $line;
    return undef if ($self->{end});
    my (@a,@c,@g,@t, $id, $tr1, $accn, $bf, $sites);
    my $i=0;
    while (defined( $line=$self->_readline)) {
	chomp($line);
	if ($line=~/^\d{2}/) {	#Begining of the frequency data
	    ($a[$i],$c[$i],$g[$i],$t[$i])=_parse_matrix($line);
	    $i++;
	}
	($tr1,$accn)=split(/\s{2}/,$line) if ($line=~/^AC\s/);
	($tr1,$bf)=split(/\s{2}/,$line) if ($line=~/^BF\s/);
	($tr1,$id)=split(/\s{2}/,$line) if ($line=~/^ID\s/);
	last if (($line=~/^XX/) && ($i>0));
    }
    if (!(defined($id) && defined($accn))) {
	$self->{end}=1;
	return undef;
    }
    while (defined( $line=$self->_readline)) {	#How many sites?
	if ($line=~/^BA\s/) {
	    my ($tr1,$ba)=split(/\s{2}/,$line);
	    ($sites)=split(/\s/,$ba);
	    last;
	}
	last if ($line=~/^\/\//);
    }
    # We have the frequencies, let's create a SiteMatrix object
    my %matrix = &_make_matrix(\@a,\@c,\@g,\@t,$id, $accn);
    $matrix{-sites}=$sites if ($sites);
    $matrix{-width}=@a;
    my $psm=new Bio::Matrix::PSM::Psm(%matrix);
    return $psm;
}

=head2 _parseMatrix

 Title   : _parseMatrix
 Usage   :
 Function: Parses a line
 Throws  :
 Example :  Internal stuff
 Returns :  array (frequencies for A,C,G,T in this order).
 Args    :  string

=cut

sub _parse_matrix {
    my $line=shift;
    $line=~s/\s+/,/g;
    my ($tr,$a,$c,$g,$t)=split(/,/,$line);
    return $a,$c,$g,$t;
}


=head2 _make_matrix

 Title   : _make_matrix
 Usage   :
 Function:
 Throws  :
 Example :  Internal stuff
 Returns :
 Args    :

=cut

sub _make_matrix {
    my ($a, $c, $g, $t, @fa, @fc,@fg, @ft, @a,@c,@g,@t);
    my $ave=0;
    my ($cA,$cC,$cG,$cT, $id, $accn)= @_;
    my $len = @{$cA}+1;

    for (my $i=0; $i <= $len;$i++) {
	# reset to 0 if value is undefined - this might need to be set properly somewhere else
	map { ${$_}[$i] ||= 0 } ($cA, $cG, $cC, $cT);
	
	if ( (${$cA}[$i] + ${$cC}[$i] + 
	      ${$cG}[$i] + ${$cT}[$i] ) ==0 ) {
	    push @a,$ave;
	    push @c,$ave;
	    push @g,$ave;
	    push @t,$ave;
	}
	else {
	    push @a,${$cA}[$i];
	    push @c,${$cC}[$i];
	    push @g,${$cG}[$i];
	    push @t,${$cT}[$i];
	    $ave = ((${$cA}[$i]+${$cC}[$i]+
		     ${$cG}[$i]+${$cT}[$i]) / 4 +$ave)/2;
	}
    }
    warn("a is $#a\n");
    for (my $i=0; $i<$#a;$i++) {
	my $zero=($a[$i]+$c[$i]+$g[$i]+$t[$i]);
	next if ($zero==0);
	push @fa, $a[$i];
	push @fc, $c[$i];
	push @fg, $g[$i];
	push @ft, $t[$i];
    }
    return (-pA=>\@fa,-pC=>\@fc,-pG=>\@fg,-pT=>\@ft, -id=>$id, -accession_number=>$accn)
    }


sub DESTROY {
    my $self=shift;
    $self->close;
}

1;
  
