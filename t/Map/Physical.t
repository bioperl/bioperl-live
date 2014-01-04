# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 40);
	
    use_ok('Bio::Map::Physical');
    use_ok('Bio::MapIO');
}

ok my $phm = Bio::Map::Physical->new();
is $phm->version(2), 2;
is $phm->version(), 2;
is $phm->modification_user('me'), 'me';
is $phm->modification_user(), 'me';

is $phm->group_type('xx'), 'xx';
is $phm->group_type(), 'xx';

is $phm->group_abbr('xx'), 'xx';
is $phm->group_abbr(), 'xx';

is $phm->core_exists, undef, 'code holds and returns a string, definition requires a boolean';

is $phm->core_exists(3), 1, 'code holds and returns a string, definition requires a boolean';

is $phm->core_exists(1), 1;
is $phm->core_exists(), 1;

my $fpcpath = test_input_file('biofpc.fpc');

# TODO? get Bio::MapIO::fpc to load from a Bio::MapIO call
my $mapio = Bio::MapIO->new(-format => "fpc", -species => 'demo', -readcor => 1, -file => $fpcpath);
my $fobj = $mapio->next_map();

is $fobj->group_abbr(), "Chr";
is $fobj->core_exists(), 1;

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
    my @ctgpos;

    my $f = shift;
    foreach my $mid ($f->each_markerid())
    {
        $nmrk++;
        my $mobj = $f->get_markerobj($mid);
        if (not defined $mobj)
        {
            is 1, 0;
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
            push @ctgpos, $mobj->position($ctgid);
        }
    }
    is $nmrk, 15;
    is $nrem, 17;
    is scalar(keys %types), 2;
    is $nanch, 9;
    is $nfrm, 7;
    is scalar (keys %grps), 4;
    is $pos, 36;
    is @ctgpos, 165;
    my $sum = 0;
    $sum += $_ for @ctgpos;
    is $sum, 1177;
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
            is 1, 0;
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
    is $nctg, 11;
    is $nchr, 3;
    is $nuser, 1;
    is $ntrace, 1;
    is $ncb, 880; 
    is $psum, 15.55;
    is scalar(keys %grps), 3;
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
            is 1, 0;
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
                is 1, 0;
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
    is $nclones, 355;
    is $nbands, 9772;
    is scalar(keys %ctgs), 11;
    is $nmrkhits, 46;
    is $nrem, 12;
    is $nfprem, 162;
    is scalar(keys %stati), 5;
}
