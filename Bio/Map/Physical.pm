#
# BioPerl module for Bio::Map::Physical
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Sendu Bala <bix@sendu.me.uk>
#
# Copyright AGCoL
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Map::Physical - A class for handling a Physical Map (such as FPC)

=head1 SYNOPSIS

    use Bio::MapIO;

    # accquire a Bio::Map::Physical using Bio::MapIO::fpc
    my $mapio = Bio::MapIO->new(-format => "fpc",-file => "rice.fpc",
                               -readcor => 0);

    my $physical = $mapio->next_map();

    # get all the markers ids
    foreach my $marker ( $physical->each_markerid() ) {
      print "Marker $marker\n";

      # acquire the marker object using Bio::Map::FPCMarker
      my $markerobj = $physical->get_markerobj($marker);

      # get all the clones hit by this marker
      foreach my $clone ($markerobj->each_cloneid() ) {
          print " +++$clone\n";
      }
  }

=head1 DESCRIPTION

This class is basically a continer class for a collection of Contig maps and
other physical map information.

Bio::Map::Physical has been tailored to work for FPC physical maps, but
could probably be used for others as well (with the appropriate MapIO
module).

This class also has some methods with specific functionalities:

  print_gffstyle()     : Generates GFF; either Contigwise[Default] or
                         Groupwise

  print_contiglist()   : Prints the list of Contigs, markers that hit the
                         contig, the global position and whether the marker
                         is a placement (<P>) or a Framework (<F>) marker.

  print_markerlist()   : Prints the markers list; contig and corresponding
                         number of clones.

  matching_bands()     : Given two clones [and tolerence], this method
                         calculates how many matching bands do they have.

  coincidence_score()  : Given two clones [,tolerence and gellen], this
                         method calculates the Sulston Coincidence score.

For faster access and better optimization, the data is stored internally in
hashes. The corresponding objects are created on request.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

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
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Gaurav Gupta

Email gaurav@genome.arizona.edu

=head1 CONTRIBUTORS

Sendu Bala  bix@sendu.me.uk

=head1 PROJECT LEADERS

Jamie Hatfield      jamie@genome.arizona.edu
Dr. Cari Soderlund  cari@genome.arizona.edu

=head1 PROJECT DESCRIPTION

The project was done in Arizona Genomics Computational Laboratory (AGCoL)
at University of Arizona.

This work was funded by USDA-IFAFS grant #11180 titled "Web Resources for 
the Computation and Display of Physical Mapping Data".

For more information on this project, please refer: 
  http://www.genome.arizona.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Map::Physical;
use vars qw($MAPCOUNT);
use strict;
use POSIX;

use Bio::Map::Clone;
use Bio::Map::Contig;
use Bio::Map::FPCMarker;

use base qw(Bio::Map::SimpleMap);
BEGIN { $MAPCOUNT = 1; }

=head1 Access Methods

These methods let you get and set the member variables

=head2 version

 Title   : version
 Usage   : my $version = $map->version();
 Function: Get/set the version of the program used to
           generate this map
 Returns : scalar representing the version
 Args    : none to get, OR string to set

=cut

sub version {
    my ($self,$value) = @_;
    if (defined($value)) {
	$self->{'_version'} = $value;
    }
    return $self->{'_version'};
}

=head2 modification_user

 Title   : modification_user
 Usage   : my $modification_user = $map->modification_user();
 Function: Get/set the name of the user who last modified this map
 Returns : scalar representing the username
 Args    : none to get, OR string to set

=cut

sub modification_user {
    my ($self,$value) = @_;
    if (defined($value)) {
	$self->{'_modification_user'} = $value;
    }
    return $self->{'_modification_user'};
}

=head2 group_type

 Title   : group_type
 Usage   : $map->group_type($grptype);
	       my $grptype = $map->group_type();
 Function: Get/set the group type of this map
 Returns : scalar representing the group type
 Args    : none to get, OR string to set

=cut

sub group_type {
    my ($self,$value) = @_;
    if (defined($value)) {
	$self->{'_grouptype'} = $value;
    }
    return $self->{'_grouptype'};
}

=head2 group_abbr

 Title   : group_abbr
 Usage   : $map->group_abbr($grpabbr);
	       my $grpabbr = $map->group_abbr();
 Function: get/set the group abbrev of this map
 Returns : string representing the group abbrev
 Args    : none to get, OR string to set

=cut

sub group_abbr {
    my ($self,$value) = @_;
    if (defined($value)) {
	$self->{'_groupabbr'} = $value;
    }
    return $self->{'_groupabbr'};
}

=head2 core_exists

 Title   : core_exists
 Usage   : my $core_exists = $map->core_exists();
 Function: Get/set if the FPC file is accompanied by COR file
 Returns : boolean
 Args    : none to get, OR 1|0 to set

=cut

sub core_exists {
    my ($self,$value) = @_;
    if (defined($value)) {
	$self->{'_corexists'} = $value ? 1 : 0;
    }
    return $self->{'_corexists'};
}

=head2 each_cloneid

 Title   : each_cloneid
 Usage   : my @clones = $map->each_cloneid();
 Function: returns an array of clone names
 Returns : list of clone names
 Args    : none

=cut

sub each_cloneid {
    my ($self) = @_;
    return keys %{$self->{'_clones'}};
}

=head2 get_cloneobj

 Title   : get_cloneobj
 Usage   : my $cloneobj = $map->get_cloneobj('CLONEA');
 Function: returns an object of the clone given in the argument
 Returns : object of the clone
 Args    : scalar representing the clone name

=cut

sub get_cloneobj {
    my ($self,$clone) = @_;

    return 0     if(!defined($clone));
    return if($clone eq "");
    return if(!exists($self->{'_clones'}{$clone}));

    my ($type,$contig,$bands,$gel,$group,$remark,$fp_number);
    my ($sequence_type,$sequence_status,$fpc_remark,@amatch,@pmatch,@ematch,
        $startrange,$endrange);
    my %clones = %{$self->{'_clones'}{$clone}};
    my @markers;

    if (ref($clones{'clone'}) eq 'Bio::Map::Clone') {
	return $clones{'clone'};
    }

    $type    = $clones{'type'}              if (exists($clones{'type'}));
    @markers = (keys %{$clones{'markers'}}) if (exists($clones{'markers'}));
    $contig  =  $clones{'contig'}           if (exists($clones{'contig'}));
    $bands   =  $clones{'bands'}            if (exists($clones{'bands'}));
    $gel     =  $clones{'gel'}              if (exists($clones{'gel'}));
    $group   =  $clones{'group'}            if (exists($clones{'group'}));
    $remark  =  $clones{'remark'}           if (exists($clones{'remark'}));

    $fp_number  =  $clones{'fp_number'}  if (exists($clones{'fp_number'}));
    $fpc_remark =  $clones{'fpc_remark'} if (exists($clones{'fpc_remark'}));

    $sequence_type   =  $clones{'sequence_type'}
        if (exists($clones{'sequence_type'}));
    $sequence_status =  $clones{'sequence_status'}
        if (exists($clones{'sequence_status'} ));

    @amatch  =  (keys %{$clones{'matcha'}})  if (exists($clones{'matcha'}));
    @ematch  =  (keys %{$clones{'matche'}})  if (exists($clones{'matche'}));
    @pmatch  =  (keys %{$clones{'matchp'}})  if (exists($clones{'matchp'}));

    $startrange =  $clones{'range'}{'start'}
        if (exists($clones{'range'}{'start'}));
    $endrange   =  $clones{'range'}{'end'}
        if (exists($clones{'range'}{'end'}));

    #*** why doesn't it call Bio::Map::Clone->new ? Seems dangerous...
    my $cloneobj = bless( {
	_name       => $clone,
	_markers    => \@markers,
	_contig     => $contig,
	_type       => $type,
	_bands      => $bands,
	_gel        => $gel,
	_group      => $group,
	_remark     => $remark,
	_fpnumber   => $fp_number,
	_sequencetype   => $sequence_type,
	_sequencestatus => $sequence_status,
	_fpcremark      => $fpc_remark,
	_matche     => \@ematch, 		
	_matcha     => \@amatch,
	_matchp     => \@pmatch,
	_range      => Bio::Range->new(-start => $startrange,
				       -end   => $endrange),	
    }, 'Bio::Map::Clone'); 		

    $self->{'_clones'}{$clone}{'clone'} = $cloneobj;
    return $cloneobj;
}

=head2 each_markerid

 Title   : each_markerid
 Usage   : my @markers = $map->each_markerid();
 Function: returns list of marker names
 Returns : list of marker names
 Args    : none

=cut

sub each_markerid {
   my ($self) = @_;
   return keys (%{$self->{'_markers'}});
}

=head2 get_markerobj

 Title   : get_markerobj
 Usage   : my $markerobj = $map->get_markerobj('MARKERA');
 Function: returns an object of the marker given in the argument
 Returns : object of the marker
 Args    : scalar representing the marker name

=cut

sub get_markerobj {
    my ($self,$marker) = @_;

    return 0 if(!defined($marker));
    return if($marker eq "");
    return if(!exists($self->{'_markers'}{$marker}));

    my ($global,$framework,$group,$anchor,$remark,$type,$linkage,$subgroup);
    my %mkr = %{$self->{'_markers'}{$marker}};

    return $mkr{'marker'} if (ref($mkr{'marker'}) eq 'Bio::Map::FPCMarker');

    $type       = $mkr{'type'}       if(exists($mkr{'type'}));
    $global     = $mkr{'global'}     if(exists($mkr{'global'} ));
    $framework  = $mkr{'framework'}  if(exists($mkr{'framework'}));
    $anchor     = $mkr{'anchor'}     if(exists($mkr{'anchor'}));
    $group      = $mkr{'group'}      if(exists($mkr{'group'}));
    $subgroup   =  $mkr{'subgroup'}  if(exists($mkr{'subgroup'}));
    $remark     =  $mkr{'remark'}    if(exists($mkr{'remark'}));

    my %clones  = %{$mkr{'clones'}};
    my %contigs = %{$mkr{'contigs'}};

    my %markerpos = %{$mkr{'posincontig'}} if(exists($mkr{'posincontig'}));

    #*** why doesn't it call Bio::Map::FPCMarker->new ? Seems dangerous...
    my $markerobj = bless( {
	_name    => $marker,
	_type    => $type,
	_global  => $global,
	_frame   => $framework,
    _group   => $group,
	_subgroup   => $subgroup,
	_anchor     => $anchor,
    _remark     => $remark,
	_clones     => \%clones,
	_contigs    => \%contigs,
	_position   => \%markerpos,	
    }, 'Bio::Map::FPCMarker');

    $self->{'_markers'}{$marker}{'marker'} = $markerobj;
    return $markerobj;
}

=head2 each_contigid

 Title   : each_contigid
 Usage   : my @contigs = $map->each_contigid();
 Function: returns a list of contigs (numbers)
 Returns : list of contigs
 Args    : none

=cut

sub each_contigid {
    my ($self) = @_;
    return keys (%{$self->{'_contigs'}});
}

=head2 get_contigobj

 Title   : get_contigobj
 Usage   : my $contigobj = $map->get_contigobj('CONTIG1');
 Function: returns an object of the contig given in the argument
 Returns : object of the contig
 Args    : scalar representing the contig number

=cut

sub get_contigobj {
    my ($self,$contig) = @_;

    return 0     if(!defined($contig));
    return if($contig eq "");
    return if(!exists($self->{'_contigs'}{$contig}));

    my ($group,$anchor,$uremark,$tremark,$cremark,$startrange,$endrange,
	$linkage,$subgroup);
    my %ctg = %{$self->{'_contigs'}{$contig}};
    my (%position, %pos);

    return $ctg{'contig'} if (ref($ctg{'contig'}) eq 'Bio::Map::Contig');

    $group        =  $ctg{'group'}        if (exists($ctg{'group'}));
    $subgroup     =  $ctg{'subgroup'}     if (exists($ctg{'subgroup'}));
    $anchor       =  $ctg{'anchor'}       if (exists($ctg{'anchor'}));
    $cremark      =  $ctg{'chr_remark'}   if (exists($ctg{'chr_remark'}));
    $uremark      =  $ctg{'usr_remark'}   if (exists($ctg{'usr_remark'}));
    $tremark      =  $ctg{'trace_remark'} if (exists($ctg{'trace_remark'}));

    $startrange =  $ctg{'range'}{'start'}
        if (exists($ctg{'range'}{'start'}));
    $endrange   =  $ctg{'range'}{'end'}
        if (exists($ctg{'range'}{'end'}));

    my %clones    =  %{$ctg{'clones'}}     if (exists($ctg{'clones'}));
    my %markers   =  %{$ctg{'markers'}}    if (exists($ctg{'markers'}));

    my $pos       =  $ctg{'position'};

    #*** why doesn't it call Bio::Map::Contig->new ? Seems dangerous...
    my $contigobj = bless( {
	_group      => $group,
	_subgroup   => $subgroup,
	_anchor     => $anchor,
	_markers    => \%markers,
	_clones     => \%clones,
	_name       => $contig,
	_cremark    => $cremark,
	_uremark    => $uremark,
	_tremark    => $tremark,
	_position   => $pos,
	_range      => Bio::Range->new(-start => $startrange,
				       -end => $endrange),	
    }, 'Bio::Map::Contig');

    $self->{'_contigs'}{$contig}{'contig'} = $contigobj;
    return $contigobj;
}

=head2 matching_bands

 Title   : matching_bands
 Usage   : $self->matching_bands('cloneA','cloneB',[$tol]);
 Function: given two clones [and tolerence], this method calculates how many
           matching bands do they have.
           (this method is ported directly from FPC)
 Returns : scalar representing the number of matching bands
 Args    : names of the clones ('cloneA', 'cloneB') [Default tolerence=7]

=cut

sub matching_bands {
    my($self,$cloneA,$cloneB,$tol) = @_;
    my($lstart,$kband,$match,$diff,$i,$j);

    return 0 if(!defined($cloneA) || !defined($cloneB) ||
		!($self->core_exists()));

    $tol = 7 if (!defined($tol));

    my %_clones  = %{$self->{'_clones'}};

    my @bandsA = @{$_clones{$cloneA}{'bands'}};
    my @bandsB = @{$_clones{$cloneB}{'bands'}};

    $match  = 0;
    $lstart = 0;

    for ($i=0; $i<scalar(@bandsA);$i++) {
	$kband = $bandsA[$i];
	for ($j = $lstart; $j<scalar(@bandsB); $j++) {
	    $diff = $kband - $bandsB[$j];
	    if (abs($diff)  <= $tol ) {
		$match++;
		$lstart = $j+1;
		last;
	    }
	    elsif ($diff < 0) {
		$lstart = $j;
		last;
	    }
	}
    }
    return $match;
}

=head2 coincidence_score

 Title   : coincidence_score
 Usage   : $self->coincidence_score('cloneA','cloneB'[,$tol,$gellen]);
 Function: given two clones [,tolerence and gellen], this method calculates
           the Sulston Coincidence score.
           (this method is ported directly from FPC)
 Returns : scalar representing the Sulston coincidence score.
 Args    : names of the clones ('cloneA', 'cloneB')
           [Default tol=7 gellen=3300.0]

=cut

sub coincidence_score {
    my($self,$cloneA,$cloneB,$tol,$gellen) = @_;

    return 0 if(!defined($cloneA) || !defined($cloneB) ||
		!($self->core_exists()));

    my %_clones  = %{$self->{'_clones'}};

    my $numbandsA = scalar(@{$_clones{$cloneA}{'bands'}});
    my $numbandsB = scalar(@{$_clones{$cloneB}{'bands'}});

    my ($nL,$nH,$m,$i,$psmn,$pp,$pa,$pb,$t,$c,$a,$n);
    my @logfact;
    my $score;

    $gellen = 3300.0 if (!defined($gellen));
    $tol    = 7      if (!defined($tol));

    if ($numbandsA > $numbandsB) {
	$nH = $numbandsA;
	$nL = $numbandsB;
    }
    else {
	$nH = $numbandsB;
	$nL = $numbandsA;
    }

    $m = $self->matching_bands($cloneA, $cloneB,$tol);

    $logfact[0] = 0.0;
    $logfact[1] = 0.0;
    for ($i=2; $i<=$nL; $i++) {
	$logfact[$i] = $logfact[$i - 1] + log($i);
    }

    $psmn = 1.0 - ((2*$tol)/$gellen);

    $pp = $psmn ** $nH;
    $pa = log($pp);
    $pb = log(1 - $pp);
    $t  = 1e-37;

    for ($n = $m; $n <= $nL; $n++)  {
	$c = $logfact[$nL] - $logfact[$nL - $n] - $logfact[$n];
	$a = exp($c + ($n * $pb) + (($nL - $n) * $pa));
	$t += $a;
    }

    $score = sprintf("%.e",$t);
    return $score;
}

=head2 print_contiglist

 Title   : print_contiglist
 Usage   : $map->print_contiglist([showall]); #[Default 0]
 Function: prints the list of contigs, markers that hit the contig, the
           global position and whether the marker is a placement (P) or
           a Framework (F) marker.
 Returns : none
 Args    : [showall] [Default 0], 1 includes all the discrepant markers

=cut

sub print_contiglist{
    my ($self,$showall) = @_;
    my $pos;

    $showall = 0 if (!defined($showall));
    my %_contigs = %{$self->{'_contigs'}};
    my %_markers = %{$self->{'_markers'}};
    my %_clones  = %{$self->{'_clones'}};

    my @contigs       = $self->each_contigid();
    my @sortedcontigs = sort {$a <=> $b } @contigs;

    print "\n\nContig List\n\n";
    foreach my $contig (@sortedcontigs) {
        my %list;
	my %alist;
	
	my $ctgAnchor  = $_contigs{$contig}{'anchor'};
	my $ctgGroup   = $_contigs{$contig}{'group'};	
	
	my @mkr = keys ( %{$_contigs{$contig}{'markers'}} );
	
	foreach my $marker (@mkr)  {	
	    my $mrkGroup       = $_markers{$marker}{'group'};
	    my $mrkGlobal      = $_markers{$marker}{'global'};
	    my $mrkFramework   = $_markers{$marker}{'framework'};
	    my $mrkAnchor      = $_markers{$marker}{'anchor'}; 	    	

	    if($ctgGroup =~ /\d+|\w/ && $ctgGroup != 0)  {		
		if ($mrkGroup eq $ctgGroup) {
		    if ($mrkFramework == 0)  {		
			$pos = $mrkGlobal."P";
		    }
		    else {
			$pos = $mrkGlobal."F";
		    }		
		    $list{$marker} = $pos;
		}
		elsif ($showall == 1) {			
		    my $chr = $self->group_abbr().$mrkGroup;
		    $alist{$marker} = $chr;
		} 	
	    }
	    elsif ($showall == 1 &&  $ctgGroup !~ /\d+/) {
		my $chr = $self->group_abbr().$mrkGroup;
		$alist{$marker} = $chr;
	    }
	}
	
	my $chr = $ctgGroup;
	$chr = $self->group_abbr().$ctgGroup if ($ctgGroup =~ /\d+|\w/);
	
	if ($showall == 1 ) {
	   	
	    print "   ctg$contig  ", $chr, "  "
		if ($_contigs{$contig}{'group'} !~ /\d+|\w/);  		
        }
	elsif ($ctgGroup =~ /\d+|\w/ && $ctgGroup ne 0){
	        print "   ctg",$contig, "  ",$chr, "  ";
	}  	
	
	while (my ($k,$v) = each %list) {
            print "$k/$v  ";		
	}
	
	print "\n" if ($showall == 0 && $ctgGroup =~ /\d+|\w/ &&
		       $ctgGroup ne 0 );
	
	if ($showall == 1) {
            while (my ($k,$v) = each %alist) {
		print "$k/$v  ";		
            }  		
	    print "\n";
        }
    }
}

=head2 print_markerlist

 Title    : print_markerlist
 Usage    : $map->print_markerlist();
 Function : prints the marker list; contig and corresponding number of
            clones for each marker.
 Returns  : none
 Args     : none

=cut

sub print_markerlist {
    my ($self) = @_;

    my %_contigs = %{$self->{'_contigs'}};
    my %_markers = %{$self->{'_markers'}};
    my %_clones  = %{$self->{'_clones'}};

    print "Marker List\n\n";

    foreach my $marker ($self->each_markerid()) {
        print "  ",$marker, "  ";
	
	my %list;
	my %mclones = %{$_markers{$marker}{'clones'}};
	
	foreach my $clone (%mclones) {
	    if (exists($_clones{$clone}{'contig'}) ) {
		my $ctg = $_clones{$clone}{'contig'};
		
		if (exists($list{$ctg})) {
		    my $clonehits = $list{$ctg};
		    $clonehits++;
		    $list{$ctg} = $clonehits;
		}
		else {
		    $list{$ctg} = 1;
		}
	    }
	}
	while (my ($k,$v) = each %list) {
	    print "$k/$v  ";
        }
        print "\n";
    }
}

=head2 print_gffstyle

 Title    : print_gffstyle
 Usage    : $map->print_gffstyle([style]);
 Function : prints GFF; either Contigwise (default) or Groupwise
 Returns  : none
 Args     : [style] default = 0 contigwise, else
                              1 groupwise (chromosome-wise).

=cut

sub print_gffstyle {
    my ($self,$style) = @_;

    $style = 0 if(!defined($style));

    my %_contigs = %{$self->{'_contigs'}};
    my %_markers = %{$self->{'_markers'}};
    my %_clones  = %{$self->{'_clones'}};

    my $i;
    my ($depth, $save_depth);
    my ($x, $y);
    my @stack;
    my ($k, $j, $s);
    my $pos;
    my $contig;

    # Calculate the position for the marker in the contig

    my @contigs       = $self->each_contigid();
    my @sortedcontigs = sort {$a <=> $b } @contigs;
    my $offset = 0;
    my %gffclones;
    my %gffcontigs;
    my %gffmarkers;
    my $basepair = 4096;

    foreach my $contig (@sortedcontigs) {
        if($_contigs{$contig}{'range'} ) {	
	    $offset =  $_contigs{$contig}{'range'}{'start'};	
	
	    if ($offset <= 0){
	        $offset = $offset * -1;	
		$gffcontigs{$contig}{'start'} = 1;
		$gffcontigs{$contig}{'end'}   =
		    ($_contigs{$contig}{'range'}{'end'} +
		     $offset ) * $basepair + 1;				
	    }
	    else {
	        $offset = 0;
		$gffcontigs{$contig}{'start'} =
		    $_contigs{$contig}{'range'}{'start'} * $basepair;
		$gffcontigs{$contig}{'end'}   =
		    $_contigs{$contig}{'range'}{'end'} * $basepair;
	    }	    		
	}
	else {
	    $gffcontigs{$contig}{'start'} = 1;
            $gffcontigs{$contig}{'end'}   = 1;		
	} 	
	
	my @clones  =  keys %{$_contigs{$contig}{'clones'}};	
	foreach my $clone (@clones) {
	    if(exists ($_clones{$clone}{'range'}) ) {
	        my $gffclone = $clone;
		
		$gffclone =~ s/sd1$//;
		
		$gffclones{$gffclone}{'start'} =
		    (($_clones{$clone}{'range'}{'start'} + $offset) *
		     $basepair + 1);

		$gffclones{$gffclone}{'end'}   =
		    (($_clones{$clone}{'range'}{'end'}
		      + $offset) * $basepair + 1);
	    }
	
	    if(!$contig) {	
	        my %markers = %{$_clones{$clone}{'markers'}}
		if (exists($_clones{$clone}{'markers'}));

	        while (my ($k,$v) = each %markers) {
		    $gffmarkers{$contig}{$k} =
		    ( ( $_clones{$clone}{'range'}{'start'} +
			$_clones{$clone}{'range'}{'end'} ) / 2 ) *
			$basepair + 1 ;
		}	
	    }
	}	
	
	if($contig) {
	    my %markers = %{$_contigs{$contig}{'markers'}}
	        if (exists($_contigs{$contig}{'markers'}));

	    while (my ($k,$v) = each %markers) {
	        $gffmarkers{$contig}{$k} = ($v + $offset) * $basepair + 1;
	    }
	}
    }

    if (!$style) {
	foreach my $contig (@sortedcontigs) {
	   	
	    if(exists ($_contigs{$contig}{'range'} )  ) {	
		print join("\t","ctg$contig","assembly","contig",
			   $gffcontigs{$contig}{'start'},
			   $gffcontigs{$contig}{'end'},".",".",".",
			   "Sequence \"ctg$contig\"; Name \"ctg$contig\"\n"
                          );
	    }
	
	    my @clones = (keys %{$_contigs{$contig}{'clones'}} );
	
	    foreach my $clone (@clones) {
		if(exists ($_clones{$clone}{'range'}) ) {	
		    print join("\t","ctg$contig","FPC");
		
		    my $type = $_clones{$clone}{'type'};
		
		    if($clone =~ /sd1$/) {
			$clone =~ s/sd1$//;
   		        $type  = "sequenced";
		    }		
		    print join ("\t","\t$type",$gffclones{$clone}{'start'},
				$gffclones{$clone}{'end'},".",".",".",
				"$type \"$clone\"; Name \"$clone\"");

		    my @markers = keys %{$_clones{$clone}{'markers'}};
		    print "; Marker_hit" if (scalar(@markers));
		
		    foreach my $mkr(@markers) {
			if (exists($_markers{$mkr}{'framework'})) {
			    print " \"$mkr ",$_markers{$mkr}{'group'}," ",
				   $_markers{$mkr}{'global'},"\"";
			}
			else {
			    print " \"$mkr 0 0\"";
			}
		    }	
		    print "; Contig_hit \"",$_clones{$clone}{'contig'},"\" "
		        if (defined($_clones{$clone}{'contig'}));
		}
		print "\n";
	    }
	
	    if (exists ($_contigs{$contig}{'markers'}) ) {	
		my %list = %{$_contigs{$contig}{'markers'}};
		
		while (my ($k,$v) = each %list) {
		    print "ctg", $contig, "\tFPC\t";
		    my $position = $gffmarkers{$contig}{$k};
		
		    my $type = "marker";
		
		    $type = "electronicmarker"
		         if ($_markers{$k}{'type'} eq "eMRK");
		
		    if( exists($_markers{$k}{'framework'})) {
			$type = "frameworkmarker"
			    if($_markers{$k}{'framework'} == 1);
			
			$type = "placementmarker"
			    if($_markers{$k}{'framework'} == 0);
		    }	
		
		    print join ("\t","$type",$position,$position,".",".",
                                ".","$type \"$k\"; Name \"$k\"");
		
	            my @clonelist;
		    my @clones  = keys %{$_markers{$k}{'clones'}};
		
		    foreach my $cl (@clones) {
			push (@clonelist, $cl)
			    if($_clones{$cl}{'contig'} == $contig);
		    }
		
		    $" = " ";
		    print("; Contig_hit \"ctg$contig - ",scalar(@clonelist),
			  "\" (@clonelist)\n");
		}
	    }  		   	
	}
    }
    else {
	my %_groups;
	my $margin       = 2 * $basepair;
	my $displacement = 0;
	my @grouplist;
	
	foreach my $contig (@sortedcontigs) {
	    my $recordchr;
            my $chr = $_contigs{$contig}{'group'};		
	    $chr = 0 if ($chr !~ /\d+|\w+/);
	
            $recordchr->{group}      = $chr;
	    $recordchr->{contig}     = $contig;
	    $recordchr->{position}   = $_contigs{$contig}{'position'};

	    push @grouplist, $recordchr;	
	}
	
	my @chr       = keys (%{$_groups{'group'}});
	my @sortedchr;
	
	if ($self->group_type eq 'Chromosome') {
	    @sortedchr = sort { $a->{'group'} <=> $b->{'group'}
				               ||
				$a->{'contig'} <=> $b->{'contig'}
                              } @grouplist;
	}
	else {
	    @sortedchr = sort { $a->{'group'}  cmp $b->{'group'} 	
				                ||
				$a->{'contig'} cmp $b->{'contig'}
                              } @grouplist;
	}
	my $lastchr   = -1;
	my $chrend    = 0;

	foreach my $chr (@sortedchr) {
	    my $chrname = $self->group_abbr().$chr->{'group'};	
	
	    if ($lastchr eq -1 || $chr->{'group'} ne $lastchr ) {
		$lastchr = $chr->{'group'} if ($lastchr eq -1);		
		$displacement = 0;	
		
		# caluclate the end position of the contig		
		my $ctgcount = 0;
		my $prevchr  = 0;		
		$chrend = 0;
		
		if ($chr->{contig} != 0) {		
		    foreach my $ch (@sortedchr) {
			if ($ch->{'group'} eq $chr->{'group'}) {
			    if($ch->{'contig'} != 0) {	
				my $ctg  = $ch->{'contig'}
				    if($ch->{'contig'} != 0);

				$chrend += $gffcontigs{$ctg}->{'end'};
				++$ctgcount;
			    }			    	
			}
		    }	
		    $chrend += ($ctgcount-1) * $margin;
		}
		else {
		    $chrend  = $gffcontigs{'0'}->{'end'};
		}
		
		$chrname    = $self->group_abbr()."ctg0"
		if ($chr->{'contig'} == 0);
		
		print join ("\t", $chrname,"assembly","Chromosome",1,
			    "$chrend",".",".",".",
			    "Sequence \"$chrname\"; Name \"$chrname\"\n");
	    }
	
	    print join ("\t", $chrname,"assembly","Chromosome",1,
			"$chrend",".",".",".",
			"Sequence \"$chrname\"; Name \"$chrname\"\n")
	        if ($chr->{'group'} ne $lastchr && $chr->{'group'} eq 0 );
	
	    $lastchr = $chr->{'group'};
	    $lastchr = -1 if ($chr->{'contig'} == 0);	
	
	    my $contig = $chr->{'contig'};
	    	
	    if(exists ($_contigs{$contig}{'range'} )  ) {
		
		print join ("\t",$chrname, "FPC","contig",
			    $gffcontigs{$contig}{'start'}+$displacement,
		            $gffcontigs{$contig}{'end'}+$displacement,
			    ".",".",".",
			    "contig \"ctg$contig\"; Name \"ctg$contig\"\n");
	    }
	
	    my @clones = (keys %{$_contigs{$contig}{'clones'}} );
	    foreach my $clone (@clones) {
		if(exists ($_clones{$clone}{'range'}) ) {	
		    print join ("\t",$chrname,"FPC");
		    my $type = $_clones{$clone}{'type'};
		
		    if ($clone =~ /sd1$/) {
			$clone =~ s/sd1$//;
			$type  = "sequenced";
		    }
		
		    print join ("\t","\t$type",$gffclones{$clone}{'start'}
				+$displacement,$gffclones{$clone}{'end'}
				+$displacement,".",".",".",
				"$type \"$clone\"; Name \"$clone\"");
		
		    my @markers = keys %{$_clones{$clone}{'markers'}};
		    print "; Marker_hit" if (scalar(@markers));
		    		
		    foreach my $mkr(@markers) {
			if (exists($_markers{$mkr}{'framework'})) {
			    print " \"$mkr ",$_markers{$mkr}{'group'}," ",
				   $_markers{$mkr}{'global'},"\"";
			}
			else {
			    print (" \"$mkr 0 0\"");
			}
		    }	
		    print "; Contig_hit \"",$_clones{$clone}{'contig'},"\" "
		        if (defined($_clones{$clone}{'contig'}));
		}
		print "\n";
	    }
	
	    if (exists ($_contigs{$contig}{'markers'}) ) {	
		my %list = %{$_contigs{$contig}{'markers'}};
		
		while (my ($k,$v) = each %list) {
		    print join ("\t",$chrname,"FPC");
		    my $type = "marker";
		
		    $type = "electronicmarker"
		        if ($_markers{$k}{'type'} eq "eMRK");
		
		    if( exists($_markers{$k}{'framework'})) {
			$type = "frameworkmarker"
			    if($_markers{$k}{'framework'} == 1);
			
			$type = "placementmarker"
			    if($_markers{$k}{'framework'} == 0);	
		    }	
		    		    		    	
		    print join ("\t","\t$type",$gffmarkers{$contig}{$k}
				+ $displacement,$gffmarkers{$contig}{$k}
				+ $displacement,".",".",".",
				"$type \"$k\"; Name \"$k\"");

		    my @clonelist;
		    my @clones  = keys %{$_markers{$k}{'clones'}};
		
		    foreach my $cl (@clones) {
			push (@clonelist, $cl)
			    if($_clones{$cl}{'contig'} == $contig);
		    }
		
		    $" = " ";		
		    print("; Contig_hit \"ctg$contig - ",
			  scalar(@clonelist),"\" (@clonelist)\n");
		}
	    }
	    $displacement += $margin + $gffcontigs{$contig}{'end'};
	}
    }
}

=head2 _calc_markerposition

 Title   : _calc_markerposition
 Usage   : $map->_calc_markerposition();
 Function: Calculates the position of the marker in the contig
 Returns : none
 Args    : none

=cut

sub _calc_markerposition {
    my ($self) = @_;
    my %_contigs = %{$self->{'_contigs'}};
    my %_markers = %{$self->{'_markers'}};
    my %_clones  = %{$self->{'_clones'}};

    my $i;
    my ($depth, $save_depth);
    my ($x, $y);
    my @stack;
    my ($k, $j, $s);
    my $pos;
    my $contig;

    # Calculate the position for the marker in the contig

    my @contigs       = $self->each_contigid();
    my @sortedcontigs = sort {$a <=> $b } @contigs;
    my $offset;
    my %gffclones;
    my %gffcontigs;

    foreach my $marker ($self->each_markerid()) {
        my (@ctgmarker, @sortedctgmarker);
	
	my @clones = (keys %{$_markers{$marker}{'clones'}})
	    if (exists ($_markers{$marker}{'clones'} ));
	
        foreach my $clone (@clones) {
	    my $record;
	    $record->{contig} = $_clones{$clone}{'contig'};		
	    $record->{start}  = $_clones{$clone}{'range'}{'start'};
	    $record->{end}    = $_clones{$clone}{'range'}{'end'};
	    push @ctgmarker,$record;
	}
	
	# sorting by contig and left position
	@sortedctgmarker = sort { $a->{'contig'} <=> $b->{'contig'}
				                  ||
				  $b->{'start'}  <=> $a->{'start'}
				                  ||
				  $b->{'end'}    <=> $a->{'end'}
			        } @ctgmarker;
				
	my $ctg = -1;
	
	for ($i=0; $i < scalar(@sortedctgmarker); $i++) {
	    if ($ctg != $sortedctgmarker[$i]->{'contig'}) {
		if ($ctg == -1) {
		    $ctg = $sortedctgmarker[$i]->{'contig'};
		}
		else  {	
		    if ($depth > $save_depth){
			$pos = ($x + $y) >> 1;
			$_contigs{$ctg}{'markers'}{$marker}      = $pos;
			$_markers{$marker}{'posincontig'}{$ctg}  = $pos;
		    }
		}
		
		$ctg      = $sortedctgmarker[$i]->{'contig'};
		$x        = $sortedctgmarker[$i]->{'start'};
		$y        = $sortedctgmarker[$i]->{'end'};
		$stack[0] = $y;
		
		$pos = ($x + $y) >> 1;
		$_contigs{$ctg}{'markers'}{$marker}     = $pos;
		$_markers{$marker}{'posincontig'}{$ctg} = $pos;
		
		$depth = $save_depth = 1;
	    }
	    elsif ($sortedctgmarker[$i]->{'end'} <= $y) {
		$stack[$depth++] = $sortedctgmarker[$i]->{'end'};
		# MAX
		if ($x < $sortedctgmarker[$i]->{'start'} ) {
		    $x = $sortedctgmarker[$i]->{'start'};
		}
		# MIN
		if ($y > $sortedctgmarker[$i]->{'end'}) {
		    $y = $sortedctgmarker[$i]->{'end'};
		}	
	    }
	    else {
		if ($depth > $save_depth) {
		    $save_depth = $depth;
		    $pos = ($x + $y) >> 1;
		    $_contigs{$ctg}{'markers'}{$marker}     = $pos;
		    $_markers{$marker}{'posincontig'}{$ctg} = $pos;
		}
		
		$x               = $sortedctgmarker[$i]->{'start'};
		$y               = $sortedctgmarker[$i]->{'end'};
		$stack[$depth++] = $y;
		
		for($j=-1, $k=0, $s=0; $s<$depth; $s++) {
		    if ($stack[$s] <$x) {
			$stack[$s] = -1;
			$j = $s if ($j == -1);
		    }
		    else {
			$k++;
			# MIN
			$y = $stack[$s] if ($y > $stack[$s]);
			if ($stack[$j] == -1) {
			    $stack[$j] = $stack[$s];
			    $stack[$s] = -1;
			    while ($stack[$j] != -1) {$j++;}
			}
			else {
			    $j = $s;
			}
		    }
		    $depth = $k;
		}	
	    }
	    if ($depth > $save_depth) {
		$pos = ($x + $y) >> 1;
		$_contigs{$ctg}{'markers'}{$marker}     = $pos;
		$_markers{$marker}{'posincontig'}{$ctg} = $pos;
	    }
	}	
    }
}

=head2 _calc_contigposition

 Title   : _calc_contigposition
 Usage   : $map->_calc_contigposition();
 Function: calculates the position of the contig in the group
 Returns : none
 Args    : none

=cut

sub _calc_contigposition{
    my ($self) = @_;

    my %_contigs = %{$self->{'_contigs'}};
    my %_markers = %{$self->{'_markers'}};
    my %_clones  = %{$self->{'_clones'}};

    my @contigs       = $self->each_contigid();
    my @sortedcontigs = sort {$a <=> $b } @contigs;

    foreach my $contig (@sortedcontigs) {
		my $position = 0;
	my $group;
	
	if (exists($_contigs{$contig}{'group'}) ) {		
	
	    my %weightedmarkers;
	    my @mkrs = keys (%{$_contigs{$contig}{'markers'}})
	        if (exists($_contigs{$contig}{'markers'})) ;

	    my $chr = $_contigs{$contig}{'group'};
	    $chr = 0 if ($_contigs{$contig}{'group'} =~ /\?/);	

	    foreach my $mkr (@mkrs) {
		if (exists($_markers{$mkr}{'group'})) {
		    if ( $_markers{$mkr}{'group'} == $chr ) {
			my @mkrclones = keys( %{$_markers{$mkr}{'clones'}});
			my $clonescount = 0;
			foreach my $clone (@mkrclones) {
			    ++$clonescount
			        if ($_clones{$clone}{'contig'} == $contig);
			}
			$weightedmarkers{$_markers{$mkr}{'global'}} =
			    $clonescount;			
		    }
		}
	    }
	
	    my $weightedctgsum = 0;
	    my $totalhits      = 0;

	    while (my ($mpos,$hits) = each %weightedmarkers) {
		$weightedctgsum += ($mpos * $hits);
		$totalhits      += $hits;
	    }
	
	    $position = sprintf("%.2f",$weightedctgsum / $totalhits)
	        if ($totalhits != 0);	
	
	    $_contigs{$contig}{'position'} = $position;	
	}
    }
}

=head2 _calc_contiggroup

 Title   : _calc_contiggroup
 Usage   : $map->_calc_contiggroup();
 Function: calculates the group of the contig
 Returns : none
 Args    : none

=cut

sub _calc_contiggroup {
    my ($self)  = @_;
    my %_contig = %{$self->{'_contigs'}};
    my @contigs = $self->each_contigid();

    foreach my $ctg (@contigs) {
        my $chr = floor($ctg/1000);
		$_contig{$ctg}{'group'} = $chr;
    }
}

=head2 _setI<E<lt>TypeE<gt>>Ref

 Title   : _set<Type>Ref
 Usage   : These are used for initializing the reference of the hash in
           Bio::MapIO (fpc.pm) to the corresponding hash in Bio::Map
           (physical.pm). Should be used only from Bio::MapIO System.
               $map->setCloneRef(\%_clones);
               $map->setMarkerRef(\%_markers);
               $map->setContigRef(\%_contigs);
 Function: sets the hash references to the corresponding hashes
 Returns : none
 Args    : reference of the hash.

=cut

sub _setCloneRef {
    my ($self, $ref)    = @_;
    %{$self->{'_clones'}} = %{$ref};
}

sub _setMarkerRef {
    my ($self, $ref)     = @_;
    %{$self->{'_markers'}} = %{$ref};
}

sub _setContigRef {
    my ($self, $ref)    = @_;
    %{$self->{'_contigs'}} = %{$ref};
}

1;
