#!/usr/bin/perl

=head1 NAME

waba2gff3.pl - convert waba output into GFF3 suitable for Gbrowse

=head1 DESCRIPTION

This script turns WABA output into GFF3 output for the query sequence.
If you need to get this where the Hit sequence is the reference
sequence you'll want to flip-flop the code to use hit instead of
query.  I didn't try and make it that general yet.

I don't (yet) know how the 'score' field is calculate by Wormbase
folks for WABA data in their GFF dumps.  I'm checking on that but it
shouldn't make a difference for Gbrowse.

=head1 AUTHOR

Jason Stajich, jason-at-bioperl-dot-org
Duke University, 

=head1 LICENSE

This script is available under the Perl Artistic License meaning you
can do with it what you wish.  

Please do tell me about bugs or improvements so I can roll those back
in for other people to use.


=cut

 

use strict;
use Bio::SearchIO;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
use Getopt::Long;

my %States = ('1' => 'coding',
	      '2' => 'coding',
	      '3' => 'coding',
	      'L'  => 'weak',
	      'H'  => 'strong',
	     );
my ($infile,$outfile,$verbose,$version);
$version = 3;
my $ptag = 'nucleotide_match';
GetOptions( 
    'i|input:s'  => \$infile,
    'o|output:s' => \$outfile,
    'v|verbose'  => \$verbose,
    'version'    => \$version,
    'p|primary|primary_tag:s' => \$ptag,
);
$infile = shift unless $infile;

my $in;

if( $infile ) {
    $in = new Bio::SearchIO(-verbose => $verbose,
			    -format  => 'waba',
			    -file    => $infile);
} else {
    $in = new Bio::SearchIO(-verbose => $verbose,
			    -format  => 'waba',
			    -fh      => \*ARGV);
}

my $out;
if( defined $outfile) {
    $out = new Bio::Tools::GFF(-gff_version => $version,
			       -file => ">$outfile",
			       -verbose => $verbose);
} else {
    $out = new Bio::Tools::GFF(-gff_version => $version,
			       -verbose     => $verbose);
}

while( my $r = $in->next_result ) {
    while( my $hit = $r->next_hit ) {
	while( my $hsp = $hit->next_hsp ) {
	    # now split this HSP up into pieces
	    my ($qs,$qe,$hs,$he)= ($hsp->query->start,
				   $hsp->query->end,
				   $hsp->hit->start,
				   $hsp->hit->end);
	    my $i = 0;
	    # grab the HMM states from Jim's WABA output
	    my $stateseq = $hsp->hmmstate_string;
	    my $state_len = length($stateseq);
	    my ($piece,$gap,@pieces);
	    $piece = {'length'   => 0,
		      'str'      => '',
		      'start'    => $i};
	    $gap = 0;
	    
	    # parse the state string, finding the gaps (Q and T states)
	    # runs of Non Q or T letters indicate a 'piece'
	    while($i <  $state_len ) {
		my $char = substr($stateseq,$i,1);
		if($char =~ /[QT]/ ) {
		    $gap++;
		} elsif( $gap ) {
		    # just finished a gap
		    $piece->{'length'} = length($piece->{'str'});
		    push @pieces, $piece;
		    $piece = {'length' => 0,
			      'str'    => '',
			      'start'  => $i };
		    $gap = 0;
		} else {
		    $piece->{'str'} .= $char;
		}
		$i++;
	    }
	    # for each piece, this could be made up of things either 
	    # as H,L, or 123 state. 
	    # In retrospect this could all probably be contained in a 
	    # single loop, but now I'm feeling lazy. I had just converted this
	    # from using 'split' in the first place if you want to know
	    # why it is structured this way....
	    for my $piece ( @pieces ) {
		
		my $len = $piece->{length};
		my $start = $piece->{start};
		my $end = $start + $len;
		my ($j) = 0;
		my $state = substr($piece->{str},$j++,1);
		warn("start is $start end is $end len is $len\n") if $verbose;
		my ($set,@sets) = ($state);
		while( $j < $len ) {
		    my $char = substr($piece->{str},$j++,1);
		    next unless( $char);
		    if( ($char =~ /\d/ && $state =~ /\d/) ||
			($char =~ /\w/ && $char eq $state) ) {
			$set .= $char;
		    } else {
			push @sets, $set;
			$set = $state = $char;
		    }		    
		}
		push @sets, $set;
		for my $set (@sets ) {
		    my $c = substr($set,0,1);
		    if( ! $c ) {
			warn("no char for '$set'\n") if $verbose;
			next;
		    }
		    my $type ='waba_'.$States{$c};
		    my $f = Bio::SeqFeature::Generic->new(
			-start => $qs + $start,
			-end   => $qs + $start + length($set),
			-strand=> $hsp->query->strand,
			-seq_id=> $hsp->query->seq_id,
			-score => $hsp->query->score,
			-primary_tag => $ptag,
			-source_tag  => $type,
			-tag    => {
			    'ID' => $hsp->hit->seq_id
			    });
		    $f->add_tag_value('ID',$hs+$start,$hs+$start+$f->length);
		    $out->write_feature($f);
		    $start += $f->length+1;
		}
	    }
	}
    }
}
