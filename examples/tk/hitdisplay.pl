#!/usr/local/bin/perl
#
# PROGRAM  : hitdisplay.pl
# PURPOSE  : Demonstrate Bio::Tk::HitDisplay
# AUTHOR   : Keith James kdj@sanger.ac.uk
# CREATED  : Nov 1 2000
#
# Requires Bio::Tk::HitDisplay
#
# To use, just pipe Blast output into this script. Try clicking on
# the blue Subject ids with the left button to activate a callback
# or with the right button to show text describing the hit.
# 

use strict;
use Text::Wrap qw(wrap $columns);
use Bio::Tools::BPlite;
BEGIN { 
	print STDERR "This example uses deprecated BioPerl code; feel free to refactor as needed\n";
	exit;
    eval { 
	require 'Tk.pm';
	require 'Bio/Tk/HitDisplay.pm'; 
    };
    if( $@ ) {
	print STDERR "Must have bioperl-gui and Tk installed to run this test, see bioperl website www.bioperl.org for instructions on how to installed bioperl-gui modules\n";    
	exit;
    }

}
use Tk;
    $columns = 80;

my $report = Bio::Tools::BPlite->new(-fh => \*STDIN);

# Normally the code ref below is in a separate package and I do 
# something like:
#
# my $adapter = Bio::PSU::IO::Blast::HitAdapter->new;
#
# while (my $hit = $result->next_hit)
# {
#     my $text     = " ... ";
#     my $callback = sub { ... };
#     push(@hits, $adapter->($sbjct, $text, $callback));
# }
#
# It's easy to roll your own for Fasta, or whatever.

my $adapter = sub
{
    my ($sbjct, $text, $callback) = @_;

    my (@data, $expect, $percent, $length);
    my ($q_id, $s_id, $q_len, $s_len);

    while (my $hsp = $sbjct->nextHSP)
    {
	$q_id ||= $hsp->query->seqname;
	$s_id ||= $hsp->subject->seqname;

	$q_len ||= $hsp->query->seqlength;
	$s_len ||= $hsp->subject->seqlength;

	my $q_x1 = $hsp->query->start;
	my $q_x2 = $hsp->query->end;

	my $s_x1 = $hsp->subject->start;
	my $s_x2 = $hsp->subject->end;

	push(@data, [$q_x1, $q_x2,
		     $s_x1, $s_x2]);

	if (defined $expect)
	{
	    if ($hsp->P < $expect)
	    {
		$expect  = $hsp->P;
		$percent = $hsp->percent;
		$length  = $hsp->length;
	    }
	}
	else
	{
	    $expect  = $hsp->P;
	    $percent = $hsp->percent;
	    $length  = $hsp->length;
	}
    }

    return { q_id     => $q_id,
	     s_id     => $s_id,
	     expect   => $expect,
	     score    => $percent,
	     overlap  => $length,
	     q_len    => $q_len,
	     s_len    => $s_len,
	     data     => \@data,
	     text     => $text,
	     callback => $callback }

};

my @hits;

while (my $sbjct = $report->nextSbjct)
{
    # Make some text to show when the left button is clicked
    my $text = wrap("", "", "Blast hit to: ", $sbjct->name, "\n");

    # Make a callback to actiavte when the right button is clicked
    my $callback = sub { print "Blast hit to ", $sbjct->name, "\n" };

    # Convert Subjct, text and callback into hash
    push(@hits, $adapter->($sbjct, $text, $callback));
}

# Create the main window and HitDisplay
my $mw = MainWindow->new;
my $hds = $mw->Scrolled('HitDisplay',
			-borderwidth => 5,
			-scrollbars  => 'ose',
			-width       => 600,
			-height      => 300,
			-background  => 'white',
			-hitcolours  => {
					 10 => 'pink',
					 20 => 'purple',
					 40 => 'yellow',
					 60 => 'gold',
					 70 => 'orange',
					 90 => 'red'
					},
			-interval    => 15,
			-hitdata     => \@hits);

$hds->pack(-side => 'top', -fill => 'both', -expand => 1);
$hds->waitVisibility;
$hds->configure(-height => 900);
$hds->configure(-scrollregion => [$hds->bbox("all")]);

MainLoop;
