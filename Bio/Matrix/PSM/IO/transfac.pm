#---------------------------------------------------------

=head1 NAME

Bio::Matrix::PSM::IO::transfac - PSM transfac parser

=head1 SYNOPSIS

See Bio::Matrix::PSM::IO for documentation

=head1 DESCRIPTION

#

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

=head1 AUTHOR - Stefan Kirov

Email skirov@utk.edu

=head1 APPENDIX

=cut


# Let the code begin...
package Bio::Matrix::PSM::IO::transfac;
use Bio::Matrix::PSM::Psm;
use Bio::Root::Root;
use Bio::Annotation::Reference;
use Bio::Annotation::Comment;
use Bio::Annotation::DBLink;
use strict;

use base qw(Bio::Matrix::PSM::PsmHeader Bio::Matrix::PSM::IO);

=head2 new

 Title   : new
 Usage   : my $psmIO =  Bio::Matrix::PSM::IO->new(-format=>'transfac', 
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
    my ($file)=$self->_rearrange(['FILE'], @args);
    $self->_initialize_io(@args) || warn "Did you intend to use STDIN?"; #Read only for now
    #Remove header
    do {
	$line=$self->_readline;
	chomp $line;
	push @{$self->{unstructured}},$line if (length($line)>2); } until ($line =~ m{^//}) || (!defined($line)); #Unstructured header
    $self->_initialize;
    return $self;
}


=head2 next_psm

 Title   : next_psm
 Usage   : my $psm=$psmIO->next_psm();
 Function: Reads the next PSM from the input file, associated with this object
 Throws  : Upon finding a line, defining the matrix, where one or more positions
            are not defined, see _make_matrix
 Returns : Bio::Matrix::PSM::Psm object
 Args    : none

=cut

sub next_psm {
    my $self=shift;
    my $line;
    return if ($self->{end});
    my (@a,@c,@g,@t, $id, $tr1, @refs,$accn, $bf, $sites);
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
	return;
    }
    while (defined( $line=$self->_readline)) {	#How many sites?
	if ($line=~/^BA\s/) {
	    my ($tr1,$ba)=split(/\s{2}/,$line);
	    ($sites)=split(/\s/,$ba);
	}
   if ($line=~/^RN/) { #Adding a reference as Bio::Annotation object (self)
    # not interested in RN line itself, since has only transfac-specific
    # reference id? - no push back of line
    my $ref=_parse_ref($self);
    push @refs,$ref
  }
	last if ($line=~m{^//});
    }
    # We have the frequencies, let's create a SiteMatrix object
    my %matrix = &_make_matrix($self,\@a,\@c,\@g,\@t,$id, $accn);
    $matrix{-sites}=$sites if ($sites);
    $matrix{-width}=@a;
    my $psm=Bio::Matrix::PSM::Psm->new(%matrix);
    foreach my $ref (@refs) { $psm->add_Annotation('reference',$ref); }
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
 Throws  :  If a position is undefined, for example if you have line like this
            in the file you are parsing: 08  4,7,,9
 Example :  Internal stuff
 Returns :
 Args    :

=cut

sub _make_matrix {
    my ($a, $c, $g, $t, @fa, @fc,@fg, @ft, @a,@c,@g,@t);
    my $ave=0;
    my ($self,$cA,$cC,$cG,$cT, $id, $accn)= @_;

    for (my $i=0; $i < @{$cA};$i++) {
	#No value can be undefined -throw an exception, since setting to 0 probably would be wrong
  #If this happens it would indicate most probably that the file, being parsed is in a different format
	map {  $self->throw('Parsing error, a position is not defined') unless  defined(${$_}[$i]) } ($cA, $cG, $cC, $cT);
	
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

    for (my $i=0; $i<@a;$i++) {
	my $zero=($a[$i]+$c[$i]+$g[$i]+$t[$i]);
	next if ($zero==0);
	push @fa, $a[$i];
	push @fc, $c[$i];
	push @fg, $g[$i];
	push @ft, $t[$i];
    }
    return (-pA=>\@fa,-pC=>\@fc,-pG=>\@fg,-pT=>\@ft, -id=>$id, -accession_number=>$accn)
    }

sub _parse_ref {
my $self=shift;
my ($authors,$title,$loc,@refs,$tr,$db,$dbid);
    while (my $refline=$self->_readline) { #Poorely designed, should go through an array with fields
      chomp $refline;
      my ($field,$arg)=split(/\s+/,$refline,2);
      last if ($field=~/XX/);
      $field.=' ';
      REF: {
          if ($field=~/RX/) {  #DB Reference
              $refline=~s/[;\.]//g;
              ($tr, $db, $dbid)=split(/\s+/,$refline);
              last REF;
          }
         if ($field=~/RT/) {   #Title
            $title .= $arg;
            last REF;
          }
          if ($field=~/RA/) {  #Author
            $authors .= $arg;
            last REF;
          }
          if ($field=~/RL/) {  #Journal
            $loc .= $arg;
            last REF;
          }
        }
     }
     my $reference=Bio::Annotation::Reference->new(-authors=>$authors, -title=>$title,
                                                    -location=>$loc);
     if ($db eq 'MEDLINE') {
        # does it ever equal medline?
        $reference->medline($dbid);
     }
     elsif ($dbid) {
        $reference->pubmed($dbid);
     }
     return $reference;
}

sub DESTROY {
    my $self=shift;
    $self->close;
}

1;
  
