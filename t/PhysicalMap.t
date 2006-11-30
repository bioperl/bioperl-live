# -*-Perl-*-
## Bioperl Test Harness Script for Modules
## $Id$
#

use strict;

BEGIN {
    use vars qw($DEBUG);
    $DEBUG = $ENV{'BIOPERLDEBUG'};
    # to handle systems with no installed Test module
    # we include the t dir (where a copy of Test.pm is located)
    # as a fallback
    eval { require Test; };
    if( $@ ) {
        use lib 't';
    }
    use Test;
    plan tests => 38;
}

use Bio::Map::Clone;
use Bio::Map::Contig;
use Bio::Map::FPCMarker;
use Bio::Map::OrderedPositionWithDistance;

use Bio::Map::Physical;
ok 1;

ok my $phm = new Bio::Map::Physical;
ok $phm->version(2), 2;
ok $phm->version(), 2;
ok $phm->modification_user('me'), 'me';
ok $phm->modification_user(), 'me';

ok $phm->group_type('xx'), 'xx';
ok $phm->group_type(), 'xx';

ok $phm->group_abbr('xx'), 'xx';
ok $phm->group_abbr(), 'xx';

ok $phm->core_exists, undef, 'code holds and returns a string, definition requires a boolean';

ok $phm->core_exists(3), 1, 'code holds and returns a string, definition requires a boolean';

ok $phm->core_exists(1), 1;
ok $phm->core_exists(), 1;


use Bio::MapIO::fpc;

my $fpcpath = Bio::Root::IO->catfile('t','data','biofpc.fpc');

my $mapio = new Bio::MapIO(-format => "fpc", -species => 'demo', -readcor => 1, -file => $fpcpath);
my $fobj = $mapio->next_map();

ok $fobj->group_abbr(), "Chr";
ok $fobj->core_exists(), 1;

test_clones($fobj);
test_contigs($fobj);
test_markers($fobj);

#########################################################

sub test_markers
{
    my $nmrk = 0;
    my $nrem = 0;
    my %types;
    my $nanch = 0;
    my $nfrm = 0;
    my %grps;
    my $pos = 0;
    my $ctgpos = 0;

    my $f = shift;
    foreach my $mid ($f->each_markerid())
    {
        $nmrk++;
        my $mobj = $f->get_markerobj($mid);
        if (not defined $mobj)
        {
            ok 1, 0;
            next;
        }
        my @remarks = split /\n/, $mobj->remark();
        $nrem += scalar(@remarks);
        $types{$mobj->type()} = 1;
        if ($mobj->anchor())
        {
            $nanch++;
            $grps{$mobj->group()} = 1;
            $pos += $mobj->global();
        }
        if ($mobj->framework())
        {
            $nfrm++;
        }
        foreach my $ctgid ($f->each_contigid())
        {
            $ctgpos += $mobj->position($ctgid);
        }
    }
    ok $nmrk, 15;
    ok $nrem, 17;
    ok scalar(keys %types), 2;
    ok $nanch, 9;
    ok $nfrm, 7;
    ok scalar (keys %grps), 4;
    ok $pos, 36;
    ok $ctgpos, 1249;
}

#########################################################

sub test_contigs
{
    my $f = shift;
    my $nchr = 0;
    my $nuser = 0;
    my $ntrace = 0;
    my $nctg = 0;
    my $ncb = 0;
    my $psum = 0;
    my %grps;
    
    foreach my $cid ($f->each_contigid())
    {
        $nctg++;
        my $cobj = $f->get_contigobj($cid);
        if (not defined $cobj)
        {
            ok 1, 0;
            next;
        }
        if ($cobj->chr_remark() ne "")
        {
            $nchr++;
        }
        if ($cobj->user_remark() eq "test")
        {
            $nuser++;
        }
        if ($cobj->trace_remark() eq "test")
        {
            $ntrace++;
        }
        if ($cid > 0)
        {
            $ncb += ($cobj->range()->end() - $cobj->range()->start() + 1);
        }
        if ($cobj->anchor())
        {
            $psum += $cobj->position(); 
            $grps{$cobj->group()} = 1;
        }
    }
    ok $nctg, 11;
    ok $nchr, 3;
    ok $nuser, 1;
    ok $ntrace, 1;
    ok $ncb, 880; 
    ok $psum, 15.55;
    ok scalar(keys %grps), 3;
}

#########################################################

sub test_clones
{
    my $f = shift;
    my $nclones = 0;
    my $nbands = 0;
    my $nrem = 0;
    my %ctgs;
    my $nmrkhits = 0;
    my $nfprem = 0;
    my %stati;
    foreach my $cid ($f->each_cloneid())
    {
        $nclones++;
        my $cobj = $f->get_cloneobj($cid);
        if (not defined $cobj)
        {
            ok 1, 0;
            next;
        }
        my $pbands = $cobj->bands();
        $nbands += scalar(@$pbands);
        $ctgs{$cobj->contigid()} = 1;
        if ($cobj->contigid() > 0)
        {
            if (not defined $cobj->range()->start() or 
                not defined $cobj->range()->end() or
                $cobj->range()->end() < $cobj->range()->start())
            {
                ok 1, 0;
            }
        }
        foreach my $mid ($cobj->each_markerid())
        {
            $nmrkhits++;
        }
        my @remarks;
        if ($cobj->remark) {
            @remarks = split /\n/, $cobj->remark();
            $nrem += scalar(@remarks);
        }
        if ($cobj->fpc_remark) {
            @remarks = split /\n/, $cobj->fpc_remark();
            $nfprem += scalar(@remarks);
        }
        $stati{$cobj->sequence_status()} = 1 if $cobj->sequence_status;
    }
    ok $nclones, 355;
    ok $nbands, 9772;
    ok scalar(keys %ctgs), 11;
    ok $nmrkhits, 46;
    ok $nrem, 12;
    ok $nfprem, 162;
    ok scalar(keys %stati), 5;
}