# fpc.pm,v 1.2.2.1 2005/10/09 15:16:27 jason Exp
#
# BioPerl module for Bio::MapIO::fpc
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Gaurav Gupta <gaurav@genome.arizona.edu>
#
# Copyright AGCoL
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::MapIO::fpc - A FPC Map reader

=head1 SYNOPSIS

# do not use this object directly it is accessed through the Bio::MapIO system

    use Bio::MapIO;

     -format  : specifies the format of the file format is "fpc",
     -file    : specifies the name of the .fpc file
     -readcor : boolean argument, indicating if .cor is to be read
                 or not. It looks for the .cor file in the same path
                 as .fpc file.
                 0 : doesn't read .cor file
                 1 : reads the .cor file
                 [default 0]
     -verbose : indicates the process of loading of fpc file
    my $mapio = Bio::MapIO->new(-format  => "fpc",
                               -file    => "rice.fpc",
                               -readcor => 0,
                               -verbose => 0);

    my $map = $mapio->next_map();

    foreach my $marker ( $map->each_markerid() ) {
         # loop through the markers associated with the map
         # likewise for contigs, clones, etc.
    }


=head1 DESCRIPTION

This object contains code for parsing and processing FPC files and creating
L<Bio::Map::Physical> object from it.

For faster access and better optimization, the data is stored internally in
hashes. The corresponding objects are created on request.

We handle reading of the FPC ourselves, since MapIO module of Bioperl adds
too much overhead.

=cut

# Let the code begin...

package Bio::MapIO::fpc;
use strict;
use POSIX;

use Bio::Map::Physical;
use Bio::Map::Clone;
use Bio::Map::Contig;
use Bio::Map::FPCMarker;
use Bio::Range;

use base qw(Bio::MapIO);

my $_readcor;

=head1 Initializer

=head2 _initialize

 Title   : _initialize
 Usage   : called implicitly
 Function: calls the SUPER::_initialize
 Returns : nothing
 Args    : species, readcor

=cut

sub _initialize{
    my ($self,@args) = @_;
    my $species;
    $self->SUPER::_initialize(@args);
    ($species,$_readcor) = $self->_rearrange([qw(SPECIES READCOR)], @args);
    $_readcor = 0 unless (defined($_readcor));
}

=head1 Access Methods

These methods let you get and set the member variables

=head2 next_map

 Title   : next_map
 Usage   : my $fpcmap = $mapio->next_map();
 Function: gets the fpcmap from MapIO
 Returns : object of type L<Bio::Map::MapI>
 Args    : none

=cut

sub next_map{

    my ($self) = @_;

    my $line;
    my ($name,$fpcver,$moddate,$moduser,$contigcnt,$clonecnt,$markerscnt,
        $bandcnt,$marker,$seqclone);
    my ($corfile,$corindex,$BUFFER);
    my @cordata;
    my %fpcmarker;
    my ($contig, $contigNumber);
    my $curClone  = 0;
    my $curMarker = 0;
    my $curContig = 0;
    my %_clones;
    my %_markers;
    my %_contigs;
    my $ctgzeropos = 1;

    my $map = Bio::Map::Physical->new('-units' => 'CB',
                                     '-type'  => 'physical');

    my $filename = $self->file();
    my $fh = $self->{'_filehandle'};

    if (defined($_readcor)) {
        $map->core_exists($_readcor);
    }
    else {
        $map->core_exists(0);
    }

    if ($map->core_exists()) {
        $corfile = substr($filename,0,length($filename)-3)."cor";
        if (open my $CORE, '<', $corfile) {
            while( read($CORE, $BUFFER, 2) ) {
                push @cordata, unpack('n*', $BUFFER);
            }
        }
        else {
            $map->core_exists(0);
        }
    }

    ## Read in the header
    while (defined($line = <$fh>)) {
        chomp($line);

        if ($line =~ m{^//\s+fpc\s+project\s+(.+)}) { $map->name($1); }
        if ($line =~ m{^//\s+([\d.]+)}) {
            my $version = $1;
            $version =~ /((\d+)\.(\d+))(.*)/;
            $map->version($1);
            if ($line =~ /User:\s+(.+)/) { $map->modification_user($1); }
        }

        if ($line =~ m{^//\s+Framework\s+(\w+)\s+(\w+)\s+([-\w]+)\s+(\w+)\s+(\w+)\s+(.+)$})
        {
            $map->group_type($3) if ($2 eq "Label");
            $map->group_abbr($5) if ($4 eq "Abbrev");
        }

        last unless ($line =~ m{^//});
    }

    if (!defined($map->group_type()) || !defined($map->group_abbr()) ) {
        $map->group_type("Chromosome");
        $map->group_abbr("Chr");
    }

    $_contigs{0}{'range'}{'end'}   = 0;
    $_contigs{0}{'range'}{'start'} = 0;

    ## Read in the clone data
    while (defined($line = <$fh>)) {
        $marker = 0;
        $contig = 0;
        $seqclone = 0;
        $contigNumber = 0;

        my ($type,$name);
        my (@amatch,@pmatch,@ematch);

        my $bandsread = 0;

        last if ($line =~ /^Markerdata/);


        $line =~ /^(\w+)\s+:\s+"(.+)"/;

        ## these will be set if we did find the clone line
        ($type, $name) = ($1, $2);

        if ($name =~ /sd1/) {
            $seqclone = 1;
        }

        $_clones{$name}{'type'} = $type;
        $_clones{$name}{'contig'} = 0;
        $_contigs{'0'}{'clones'}{$name}  = 0;

        my $temp;

        ## Loop through the following lines, getting attributes for clone
        while (defined($line = <$fh>) && $line !~ /^\s*\n$/)  {

            if ($line =~ /^Map "ctg(\d+)" Ends (Left|Right) ([-\d]+)/)  {
                $_clones{$name}{'contig'} = $1;
                $_contigs{$1}{'clones'}{$name} = 0;

                delete($_contigs{'0'}{'clones'}{$name});

                $temp = $3;
                $contigNumber = $1;
                $line = <$fh>;
                $line =~ /^Map "ctg(\d+)" Ends (Left|Right) ([\d]+)/;
                $_clones{$name}{'range'}{'start'} = $temp;

                $_contigs{$contigNumber}{'range'}{'start'} = $temp
                    if (!exists($_contigs{$contigNumber}{'range'}{'start'})
                        || $_contigs{$contigNumber}{'range'}{'start'}
                        >  $temp );

                $_clones{$name}{'range'}{'end'} = $3;

                $_contigs{$contigNumber}{'range'}{'end'} = $3
                    if (!exists($_contigs{$contigNumber}{'range'}{'end'})
                        || $_contigs{$contigNumber}{'range'}{'end'} < $3 );

            }
            elsif ($line =~ /^([a-zA-Z]+)_match_to_\w+\s+"(.+)"/) {
                my $matchtype = "match" . lc(substr($1, 0, 1));
                $_clones{$name}{$matchtype}{$2} = 0;
            }
            elsif ($line =~ /^Positive_(\w+)\s+"(.+)"/) {
                $_clones{$name}{'markers'}{$2} = 0;
                $_markers{$2}{'clones'}{$name} = 0;
                $_markers{$2}{'type'} = $1;
                $_markers{$2}{'contigs'}{$contigNumber} = 0;
                $_contigs{$contigNumber}{'markers'}{$2} = 0;
            }
            elsif ($line =~ /^Bands\s+(\d+)\s+(\d+)/ && !$bandsread) {
                my $i = 0;
                my @numbands;
                $bandsread = 1;

                if ($map->core_exists()) {
                    while($i<$2){
                        push(@numbands,$cordata[($1-1)+$i]);
                        $i++;
                    }
                    $_clones{$name}{'bands'} = \@numbands;
                }
                else {
                    push(@numbands,$1,$2);
                    $_clones{$name}{'bands'} = \@numbands;
                }
                if (exists($_contigs{0}{'clones'}{$name})) {
                    $_clones{$name}{'range'}{'start'} = $ctgzeropos;
                    $_clones{$name}{'range'}{'end'} = $ctgzeropos + $2;
                    $_contigs{0}{'range'}{'end'} = $ctgzeropos + $2;
                    $ctgzeropos += $2;
                }
            }
            elsif ($line =~ /^Gel_number\s+(.+)/) {
                $_clones{$name}{'gel'} = $1;
            }
            elsif ($line =~ /^Remark\s+"(.+)"/)  {
                $_clones{$name}{'remark'} .= $1;
                $_clones{$name}{'remark'} .= "\n";
                if($seqclone == 1 ) {
                    if( $1 =~ /\,\s+Chr(\d+)\s+/){
                        $_clones{$name}{'group'} = $1;
                    }
                }
            }
            elsif ($line =~ /^Fp_number\s+"(.+)"/) {
                $_clones{$name}{'fp_number'} = $1;
            }
            elsif ($line =~ /^Shotgun\s+(\w+)\s+(\w+)/) {
                $_clones{$name}{'sequence_type'} = $1;
                $_clones{$name}{'sequence_status'} = $2;
            }
            elsif ($line =~ /^Fpc_remark\s+"(.+)"/) {
                $_clones{$name}{'fpc_remark'} .= $1;
                $_clones{$name}{'fpc_remark'} .= "\n";
            }
        }

        $curClone++;
        print "Adding clone $curClone...\n\r"
            if ($self->verbose()  && $curClone % 1000 == 0);
    }

    $map->_setCloneRef(\%_clones);
    $line = <$fh>;

    while (defined($line = <$fh>) && $line !~ /Contigdata/) {
        my ($type,$name);

        last if ($line !~ /^Marker_(\w+)\s+:\s+"(.+)"/);

        ($type, $name) = ($1, $2);

        $_markers{$name}{'type'}   = $type;
        $_markers{$name}{'group'}  = 0;
        $_markers{$name}{'global'} = 0;
        $_markers{$name}{'anchor'} = 0;

        while (defined($line = <$fh>) && $line !~ /^\s*\n$/)  {
            if ($line =~ /^Global_position\s+([\d.]+)\s*(Frame)?/)  {
                my $position = $1 - floor($1/1000)*1000;
                $position = sprintf("%.2f",$position);

                $_markers{$name}{'global'} = $position;
                $_markers{$name}{'group'}  = floor($1/1000);
                $_markers{$name}{'anchor'} = 1;

                if(defined($2)) {
                    $_markers{$name}{'framework'} = 1;
                }
                else {
                    $_markers{$name}{'framework'} = 0;
                }
            }
            elsif ($line =~ /^Anchor_bin\s+"([\w\d.]+)"/) {
                my $grpmatch = $1;
                my $grptype  = $map->group_type();

                $grpmatch =~ /(\d+|\w)(.*)/;

                my ($group,$subgroup);
                $group    = $1;
                $subgroup = $2;

                $subgroup = substr($subgroup,1) if ($subgroup =~ /^\./);

                $_markers{$name}{'group'}      = $group;
                $_markers{$name}{'subgroup'}   = $subgroup;
            }
            elsif ($line =~ /^Anchor_pos\s+([\d.]+)\s+(F|P)?/){
                $_markers{$name}{'global'}  = $1;
                $_markers{$name}{'anchor'}  = 1;

                if ($2 eq 'F') {
                    $_markers{$name}{'framework'} = 1;
                }
                else {
                    $_markers{$name}{'framework'} = 0;
                }
            }
            elsif ($line =~ /^anchor$/) {
                $_markers{$name}{'anchor'} = 1;
            }
            elsif ($line =~ /^Remark\s+"(.+)"/)  {
                $_markers{$name}{'remark'} .= $1;
                $_markers{$name}{'remark'} .= "\n";
            }
        }
        $curMarker++;
        print "Adding Marker $curMarker...\n"
            if ($self->verbose() && $curMarker % 1000 == 0);
    }

    $map->_setMarkerRef(\%_markers);

    my $ctgname;
    my $grpabbr = $map->group_abbr();
    my $chr_remark;

    $_contigs{0}{'group'} = 0;

    while (defined($line = <$fh>)) {

        if ($line =~ /^Ctg(\d+)/) {
            $ctgname = $1;
            $_contigs{$ctgname}{'group'}      = 0;
            $_contigs{$ctgname}{'anchor'}     = 0;
            $_contigs{$ctgname}{'position'}   = 0;

            if ($line =~ /#\w*(.*)\w*$/) {
                $_contigs{$ctgname}{'remark'} = $1;
                if ($line =~ /#\s+Chr(\d+)\s+/) {
                    $_contigs{$ctgname}{'group'}  = $1;
                    $_contigs{$ctgname}{'anchor'} = 1;
                }
            }
        }
        elsif ($line =~ /^Chr_remark\s+"(-|\+|Chr(\d+))\s+(.+)"$/) {

            $_contigs{$ctgname}{'anchor'}     = 1;
            $_contigs{$ctgname}{'chr_remark'} = $3 if(defined($3));

            if (defined($2)) {
                $_contigs{$ctgname}{'group'}  = $2;
            }
            else {
                $_contigs{$ctgname}{'group'}  = "?";
            }
        }
        elsif ($line =~ /^User_remark\s+"(.+)"/) {
            $_contigs{$ctgname}{'usr_remark'} = $1;
        }
        elsif ($line =~ /^Trace_remark\s+"(.+)"/) {
            $_contigs{$ctgname}{'trace_remark'} = $1;
        }
        elsif ($grpabbr && $line =~ /^Chr_remark\s+"(\W|$grpabbr((\d+)|(\w+)|([.\w\d]+)))\s*(\{(.*)\}|\[(.*)\])?"\s+(Pos\s+((\d.)+|NaN))(NOEDIT)?/)
        {
            my $grpmatch = $2;
            my $pos = $10;
            if ($pos eq "NaN") {
                $pos = 0;
                print "Warning: Nan encountered for Contig position \n";
            }
            $_contigs{$ctgname}{'chr_remark'}   = $6;
            $_contigs{$ctgname}{'position'} = $pos;
            $_contigs{$ctgname}{'subgroup'} = 0;

            if (defined($grpmatch)) {
                $_contigs{$ctgname}{'anchor'} = 1;

                if ($grpmatch =~ /((\d+)((\D\d.\d+)|(.\d+)))|((\w+)(\.\d+))/) {

                    my ($group,$subgroup);
                    $group    = $2 if($grpabbr eq "Chr");
                    $subgroup = $3 if($grpabbr eq "Chr");

                    $group    = $7 if($grpabbr eq "Lg");
                    $subgroup = $8 if($grpabbr eq "Lg");

                    $subgroup = substr($subgroup,1) if ($subgroup =~ /^\./);
                    $_contigs{$ctgname}{'group'}     = $group;
                    $_contigs{$ctgname}{'subgroup'}  = $subgroup;

                }
                else {
                    $_contigs{$ctgname}{'group'} = $grpmatch;
                }
            }
            else {
                $_contigs{$ctgname}{'anchor'} = 1;
                $_contigs{$ctgname}{'group'}  = "?";
            }
        }
        $curContig++;
        print "Adding Contig $curContig...\n"
            if ($self->verbose() && $curContig % 100 == 0);
    }

    $map->_setContigRef(\%_contigs);
    $map->_calc_markerposition();
    $map->_calc_contigposition() if ($map->version() < 7.0);
    $map->_calc_contiggroup() if ($map->version() == 4.6);

    return $map;
}


=head2 write_map

 Title   : write_map
 Usage   : $mapio->write_map($map);
 Function: Write a map out
 Returns : none
 Args    : Bio::Map::MapI

=cut

sub write_map{
    my ($self,@args) = @_;
    $self->throw_not_implemented();
}

1;

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

=head1 PROJECT LEADERS

Jamie Hatfield            jamie@genome.arizona.edu

Dr. Cari Soderlund        cari@genome.arizona.edu

=head1 PROJECT DESCRIPTION

The project was done in Arizona Genomics Computational Laboratory
(AGCoL) at University of Arizona.

This work was funded by USDA-IFAFS grant #11180 titled "Web Resources
for the Computation and Display of Physical Mapping Data".

For more information on this project, please refer:
  http://www.genome.arizona.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut
