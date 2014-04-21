#---------------------------------------------------------

=head1 NAME

Bio::Matrix::PSM::IO::masta - motif fasta format parser

=head1 SYNOPSIS 

MASTA is a position frequency matrix format similar to
fasta. It contains one ID row just like fasta and then the actual
data, which is tab delimited:

  0.1	0.62	.017	0.11
  0.22	0.13	0.54	0.11

Or A,C,G and T could be horizontally positioned (positioning is
automatically detected).  Please note masta will parse only DNA at the
moment.

It will also convert a set of aligned sequences:
ACATGCAT
ACAGGGAT
ACAGGCAT
ACCGGCAT

to a PFM (SiteMatrix object). When writing if you supply SEQ it will
write 10 random instances, which represent correctly the frequency and
can be used as an input for weblogo creation purposes.

See Bio::Matrix::PSM::IO for detailed documentation on how to use masta parser

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
package Bio::Matrix::PSM::IO::masta;
use Bio::Matrix::PSM::SiteMatrix;
use vars qw(@HEADER);
use strict;

use base qw(Bio::Matrix::PSM::IO Bio::Root::Root);



=head2 new

 Title   : new
 Usage   : my $psmIO =  Bio::Matrix::PSM::IO->new(-format=> 'masta',
						 -file  => $file, 
                                                 -mtype => 'PWM');
 Function: Associates a file with the appropriate parser
 Throws  :
 Example :
 Args    : hash
 Returns : "Bio::Matrix::PSM::$format"->new(@args);

=cut

sub new {
    my($class, @args)=@_;
    my $self = $class->SUPER::new(@args);
    my ($file)=$self->_rearrange(['FILE'], @args);
    my ($query,$tr1)=split(/\./,$file,2);
    $self->{file}  = $file;
    $self->{_end}  = 0;
    $self->{mtype} = uc($self->_rearrange(['MTYPE'], @args) || "PFM");
    $self->_initialize_io(@args) || $self->warn("Did you intend to use STDIN?"); #Read only for now
    return $self;
}

=head2 write_psm

 Title   : write_psm
 Usage   : 
 Function: writes a pfm/pwm/raw sequence in a simple masta format
 Throws  :
 Example : 
 Args    : SiteMatrix object, type (optional string: PWM, SEQ or PFM)
 Returns : 

=cut

sub write_psm {
    my ($self,$matrix,$type)=@_;
    $self->{mtype} = uc($type) if ($type);
    my $idline=">". $matrix->id . "\n";
    $self->_print($idline);
    unless ($self->{mtype} eq 'SEQ') {
	while (my %h=$matrix->next_pos) {
	    my $row=$self->{mtype} eq 'PWM' ? join("\t",$h{lA},$h{lC},$h{lG},$h{lT},"\n"):join("\t",$h{pA},$h{pC},$h{pG},$h{pT},"\n");
	    $self->_print ($row);
	}
    } else {
	my @seq;
	while (my %h=$matrix->next_pos) {
	    my ($a,$c,$g,$t)=_freq_to_count(\%h);
	    $self->throw("Could not convert from frequency to count\n") if (($a+$c+$g+$t) !=10);
	    for my $i (0..$a-1) {$seq[$i].='A';}
	    my $m=$a+$c;
	    for my $i ($a..$m-1) {$seq[$i].='C';}
	    my $n=$a+$c+$g;
	    for my $i ($m..$n-1) {$seq[$i].='G';}
	    for my $i ($n..9) {$seq[$i].='T';}
	}	
	foreach my $s (@seq) {
	    $s.="\n";
	    $self->_print ($s);
	}
    }
}

=head2 next_matrix

  Title   : next_matrix
  Usage   : my $matrix = $psmio->next_matrix;
  Function: Alias of next_psm function

=cut

sub next_matrix { 
    shift->next_psm(@_);
}

=head2 next_psm

 Title   : next_psm
 Usage   : my $matrix=$psmio->next_psm;
 Function: returns the next matrix in the stream
 Throws  : If there is you mix different types, for example weights and
           frequencies occur in the same entry You can mix weights, but these
           should be designated by different ID lines
 Example :
 Args    :
 Returns : Bio::Matrix::PSM::SiteMatrix

=cut

sub next_psm {
    my $self=shift;
    return if ($self->{_end});
    my $line=$self->_readline;
    $self->throw("No ID line- wrong format\n") unless ($line=~/^>/);
    my ($id,$desc)=split(/[\t\s]+/,$line,2);
    $id=~s/>//;
    my ($mtype,$format,@mdata,$len);
    $self->{_mtype} = 0;
    while ($line=$self->_readline) {
	next if $line =~ /^\s+$/;# There should not be empty lines, but just in case...
	chomp $line;
	if ($line =~ /^>/) {
	    $self->_pushback($line);
	    last;
	}

	if ($line !~ /[^ACGTacgt]/g) {
	    # This is a set of aligned sequences
	    $self->throw("Mixing between types is not allowed or a parsing error occured\n") 
		if (($self->{_mtype} != 3) && ($mtype)) ;
	    $self->throw("Bad sequence- different length: $line\n") 
		if (($len) && ($len!=length($line)));
	    $len=length($line) unless ($len);
	    push @mdata,$line;
	    $self->{_mtype}=3;
	} else {
		# do not strip 'e's since they are part of number notation for small/big numbers
	    $line=~s/[a-df-zA-DF-Z]//g; #Well we may wanna do a hash and auto check for letter order if there is a really boring talk...
	    $line=~s/^[\s\t]+//;
	    $line=~s/[\s\t]+/\t/g;
	    my @data=split(/[\s\t]+/,$line);
	    if ($#data==3) {
		$self->throw("Mixing between types is not allowed or a parsing error occured\n") if (($mtype)&&($self->{_mtype} !=1)) ;
		$self->{_mtype}=1;
		$mtype=1;
	    }
	    else   {
		$self->throw("Mixing between types is not allowedor a parsing error occured\n") if (($mtype)&&($self->{_mtype} !=2)) ;
		$self->{_mtype}=2;
		$mtype=1;
	    }
	    push @mdata,\@data;
	}
    }
    $self->{_end} = 1 if (!defined $line || $line !~ /^>/);
    return _make_matrix(\@mdata,$self->{_mtype},$id,$desc);
}

sub _make_matrix {
    my ($mdata,$type,$id,$desc)=@_;
    if ($type==1) {
	my @rearr=_rearrange_matrix($mdata); 
	$mdata=\@rearr;
    }
#Auto recognition for what type is this entry (PFM, PWM or simple count)
#A bit dangerous, I hate too much auto stuff, but I want to be able to mix different
#types in a single file
    my $mformat='count';
    my ($a,$c,$g,$t);
    if ($type == 3 ) {
	($a,$c,$g,$t)= &_count_positions($mdata);
    } else {
	($a,$c,$g,$t)=@{$mdata};	
	my $k=$a->[0]+$c->[0]+$g->[0]+$t->[0];
	my $l= ($a->[0]+$c->[0]+$g->[0]+$t->[0]) - 
	    (abs($a->[0])+abs($c->[0])+abs($g->[0])+abs($t->[0]));
	$mformat='freq' if (($k==1) && ($l==0));
	$mformat='pwm' if ($l!=0);
    }
    my (@fa,@fc,@fg,@ft,%mparam);

    if ($mformat eq 'pwm') {
	foreach my $i (0..$#{$a}) {
	    my $ca=exp $a->[$i];
	    my $cc=exp $c->[$i];
	    my $cg=exp $g->[$i];
	    my $ct=exp $t->[$i];
	    my $all=$ca+$cc+$cg+$ct;
	    push @fa,($ca/$all)*100;
	    push @fc,($cc/$all)*100;
	    push @fg,($cg/$all)*100;
	    push @ft,($ct/$all)*100;
	}
    }
    $desc.=", source is $mformat";
    if ($mformat eq 'pwm') {
	$desc=~s/^pwm//;
	%mparam=(-pA=>\@fa,-pC=>\@fc,-pG=>\@fg,-pT=>\@ft,-id=>$id,-desc=>$desc,
		 -lA=>$a,-lC=>$c,-lG=>$g,-lT=>$t);
    }
    else {
	%mparam=(-pA=>$a,-pC=>$c,-pG=>$g,-pT=>$t,-id=>$id,-desc=>$desc);
    }
    return new Bio::Matrix::PSM::SiteMatrix(%mparam);
}

sub _rearrange_matrix {
    my $mdata=shift;
    my (@a,@c,@g,@t);
    foreach my $entry (@{$mdata}) {
	my ($a,$c,$g,$t)=@$entry;
	push @a,$a;
	push @c,$c;
	push @g,$g;
	push @t,$t;
    }
    return \@a,\@c,\@g,\@t;
}


sub _count_positions {
    my $seq=shift;
    my %pos;
    my $l=length($seq->[0])-1;
    for( my $i = 0; $i <= $l; $i++ ) {
	for ( qw(A C G T) ) {
	    $pos{$_}->[$i] = 0;
	}
    }
    foreach my $sequence (@{$seq}) {
	my @let= split(//,$sequence);
	for my $i (0..$#let) {
	    $pos{uc($let[$i])}->[$i]++;
	}
    }
    return $pos{A},$pos{C},$pos{G},$pos{T};
}


sub _freq_to_count {
    my $h=shift;
    my $a=int(10*$h->{pA}+0.5);
    my $c=int(10*$h->{pC}+0.5);
    my $g=int(10*$h->{pG}+0.5);
    my $t=int(10*$h->{pT}+0.5);
    return ($a,$c,$g,$t);
}

1;
