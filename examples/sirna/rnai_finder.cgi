#!/usr/bin/perl

=pod

=head1 NAME

rnai_finder.cgi

=head1 DESCRIPTION

A CGI script using the Bio::Tools::SiRNA package to design RNAi reagents.  
Retrieves sequences from NCBI and generates output in graphic and tabular form.

=head1 INSTALLATION

To use this script, place it in an appropriate cgi-bin directory on a web server.  
The script needs to write its graphic maps to a temporary directory.  Please update
$TMPDIR and $TMPURL to suit your local configuation.

=head1 AUTHOR

Donald Jackson (donald.jackson@bms.com)

=head1 SEE ALSO

L<Bio::Tools::SIRNA>, L<Bio::Graphics::Panel>, L<Bio::DB::NCBIHelper>, L<CGI>

=cut

use Bio::Tools::SiRNA;

use Bio::Graphics::Panel;
use Bio::DB::NCBIHelper;
use Bio::Seq::RichSeq; # for hand-entry
use Bio::SeqFeature::Generic; 

use GD::Text::Align;
use Clone qw(clone);

use CGI;
use CGI::Carp qw (fatalsToBrowser carpout);

my $q = CGI->new;



# define a bunch of constants
my %COLORRANKS = ( 1 => 'red',
		   2 => 'orchid',
		   3 => 'blue' );
my $TMPDIR = '/var/www/htdocs/tmp/';
my $TMPURL = '/tmp/';
my $ATGPAD = 75; # how far from start do we wait?
my $NOLIGOS = 3;

my $log = $TMPDIR . 'RNAiFinder.log';
open (LOG, ">>$log") or die $!;
carpout(LOG);

print $q->header,
    $q->start_html;
    
print $q->h1('RNAi Finder');


if ($q->param('Design')) {
    if ($q->param('accession') and !$q->param('seq')) {
	$target = get_target();
    }
    else { 
	$target = make_target();
    }
    get_rnai($target);
}
else {
    get_settings();
}


sub get_settings {
    print <<EOM1;
<P>Oligos are designed as described on the <A HREF="http://www.mpibpc.gwdg.de/abteilungen/100/105/sirna.html" TARGET="Tuschl">Tuschl lab web page</A> and are ranked as follows: 
<UL>
<LI><B>New:</B> Selecting 'Pol3-compatible targets' looks for oligos with the 
pattern NAR(N17)YNN which can be synthesized or expressed from a Pol3 promoter.
<br>This selection <b>overrides</b> the 'Cutoff' rank.
<LI>Oligos with Rank = 1 (best) match the AAN(19)TT rule.
<LI>Oligos with Rank = 2 match the AAN(21) rule 
<LI>Oligos with Rank = 3 match the NAN(21) rule.
</UL>
<P>If percent GC and specificity are similar, Rank 1 oligos are better. All 3 prime overhangs are converted to TT; the rest of the sequence is transcribed into RNA</P>

<h3>Modifications to published rules:</h3>
<ul>
<li>
Runs of 3 or more consecutive Gs on either strand are skipped - these can cause problems in synthesis. 
<li>Users may choose to exclude oligos that overlap single nucleotide polymorphisms (ON by default).  SNP data comes from the NCBI dbSNP database.  
<li>'Low-complexity' regions (such as runs of a single nucleotide) are also excluded.
</ul>

EOM1
    
    print $q->start_form;
    print $q->h2('Enter your sequence and other parameters:'), "\n";
    print $q->p('The values already here are DEFAULTS - you should change them to suit YOUR sequence');
    print $q->start_table();
    print $q->TR( $q->td({-align=> 'left'}, 
		       [
			$q->textfield(-name 	=> 'mingc', -default => '0.40'),
			$q->textfield(-name 	=> 'maxgc', -default => '0.60'),
			]
		       ),
		      $q->td({-align=> 'left'},
			     $q->popup_menu(-name	 => 'worstrank', 
				      -values 	=> [1,2,3], 
				      -default 	=> 2,
				      ),
			     $q->b('OR'),
			     $q->checkbox(-name 		=> 'pol3',
				    -label		=> 'Pol3 compatible',
				    -default	=> 0,
				    ),
			     ),
		      );	
    print    $q->TR( $q->th({-align=> 'left'}, 'Exclude oligos with SNPs?'),
		     $q->td($q->radio_group(-name => 'avoid_snps', 
					    -values => [1,0],
					    -default => 1,
					    -labels => {1 => 'Yes', 0 => 'No'}
					    )),
		     );

    print $q->TR( $q->th({-align=> 'left'}, 'Sequence Name:'),
		  $q->td({-align=> 'left'},$q->textfield('accession')),
		  $q->td({-align=> 'left'}, 
			 $q->em( q(Enter an accession and you won't have to enter the <br>sequence or start/stop. Use accessions beginning with NM_ if possible.))),
		  );

    print $q->TR( $q->th({-align=> 'left'}, ['Position of initiator ATG:', 
					     'NT after start to exclude:',
					     'Position of Stop codon:' ]));
    print $q->TR( $q->td({-align=> 'left'}, 
					   [$q->textfield(-name => 'cdstart', -default => 1),
					    $q->textfield(-name => 'atgpad', -default => $ATGPAD),
					    $q->textfield('cdend'), ]));
    print $q->TR( $q->th({-align=> 'left'}, ['Minimum Fraction GC:', 
					     'Maximum Fraction GC:', 
					     'Rank cutoff',
					     ]));
    print $q->TR($q->th({-align=> 'left', -colspan=>2},'cDNA Sequence in plain text or FASTA format'),
		 $q->td( $q->a({-href =>'Fasta_format.html', -target => 'Fasta_desc'}, 'What is FASTA format?')),
		 );
    print $q->TR($q->td({-align => 'left', -colspan=>3},
		  $q->textarea( -name =>'seq',
				-rows => 4,
				-columns => 80, 
				-wrap => 'virtual',
				)));
    print $q->TR( $q->th({-align => 'left', -colspan=>3},
			 'Output options: '));
    print $q->TR( $q->td({-align=> 'left'},
		   [ $q->checkbox(-name => 'Graphic', -checked => 'checked'), 
		     $q->checkbox(-name =>  'Table',  -checked => 'checked'), 
						 ]));			       		   		
    print $q->TR($q->td({-align=> 'left', -colspan=>3}, $q->submit('Design')));
    print $q->end_table();		
    print $q->end_form;

}

sub get_rnai {
    # design and output RNAi reagents
    my ($gene) = @_;

    my $factory = Bio::Tools::SiRNA->new( -target 	=> $gene, 
					  -tmpdir	=> $TMPDIR,
					  -cutoff 	=> $q->param('worstrank') || 2,
					  -avoid_snps	=> $q->param('avoid_snps') || 1, 
					  -min_gc	=> $q->param('min_gc') || 0.40,
					  -max_gc	=> $q->param('max_gc') || 0.60,
					  -pol3		=> $q->param('pol3') || 0,
					  );

    print $q->p('Designing Pol3-compatible oligos') if ($q->param('pol3'));

    my @pairs = $factory->design;

    draw_gene($gene) if ($q->param('Graphic'));
    print_table($gene->accession, \@pairs) if ($q->param('Table'));
    print_text($gene->accession, \@pairs) if ($q->param('Text'));
}

sub get_target {
    my ($acc) = $q->param('accession');
    my $gb = Bio::DB::NCBIHelper->new();
    my $seq = $gb->get_Seq_by_acc($acc);

    if ($seq) { 
	return $seq;
    }
    else {
	print_error("Unable to retrieve sequence from GenBank using accession $acc");
	return;
    }

}

sub make_target {
      # sanity chex - do we have the necessary info?
      $q->param('seq') or print_error("Please supply a sequence", 1);
      my $seq = $q->param('seq');
      my $name;

      # is sequence in fasta format?
      if ($seq =~ /^>/) {
	  my ($head, $realseq) = split (/\n/, $seq, 2);
	  $head =~ /^>(.+?) /;
	  $name = $1;
	  $realseq =~ s/[\n|\r|\s]//g;
	  $seq = $realseq;
      }
      elsif ($q->param('accession')) {
  	$name = $q->param('accession');
	$seq =~ s/[\n|\r|\s]//g;
      }
      else {
  	print_error('Please supply a sequence name!');
  	return;
      }

      $cds_start = $q->param('cds_start') || 1;
      $cds_end = $q->param('cds_end') || length($seq);

      # create a new Bio::Seq::RichSeq object from parameters 
      my $seqobj = Bio::Seq::RichSeq->new( -seq 		=> $seq,
					   -accession_number	=> $name,
					   -molecule		=> 'DNA',
					   
  				     );
      my $cds = Bio::SeqFeature::Generic->new( -start 	=> $cds_start,
					       -end	=> $cds_end,
  					     );
      $cds->primary_tag('CDS');
      $seqobj->add_SeqFeature($cds);
      return $seqobj;
     
}
sub draw_gene {
# now draw a pretty picture
    my ($gene) = @_;

    my $panel = Bio::Graphics::Panel->new( -segment 	=> $gene,
					   -width 	=> 600,
					   -pad_top	=> 100,
					   -pad_bottom  => 20,
					   -pad_left	=> 50,
					   -pad_right	=> 50,
					   -fontcolor	=> 'black',
					   -fontcolor2  => 'black',
					   -key_color	=> 'white',
					   -grid	=> 1,
					   -key_style	=> 'between',
					   #-gridcolor	=> 'lightgray',
					   );

    my $genefeat = Bio::SeqFeature::Generic->new( -start	=> 1,
						  -end 	       	=> $gene->length);

    $panel->add_track( arrow	=> $genefeat,
		       -bump	=> 0,
		       -tick	=> 2,
		       -label 	=> 1,
		       );

    my %feature_classes;

    foreach $feat($gene->top_SeqFeatures) {
	$feature_classes{ $feat->primary_tag } ||= [];

	push(@{ $feature_classes{ $feat->primary_tag } }, $feat);
    }

# for some reason, Bio::Graphics insists on drawing subfeatures for SiRNA::Pair objects...
    $cleanpairs = cleanup_feature($feature_classes{'SiRNA::Pair'});

# draw
    $panel->add_track( transcript	=> $feature_classes{'gene'},
		       -bgcolor	=> 'green',
		       -fgcolor	=> 'black',
		       -fontcolor2  => 'black',
		       -key		=> 'Gene',
		       -bump	=> +1,
		       -height	=> 8,
		       -label	=> \&feature_label,
		       -description	=> 1,
		       );

    $panel->add_track( transcript2	=> $feature_classes{'CDS'},
		       -bgcolor		=> 'blue',
		       -fontcolor2  => 'black',
		       -fgcolor		=> 'black',
		       -key		=> 'CDS',
		       -bump		=> +1,
		       -height		=> 8,
		       -label		=> \&feature_label,
		       -description		=> \&feature_desc,
		       );

    $panel->add_track( $feature_classes{'variation'},
		       -bgcolor	=> 'black',
		       -fgcolor	=> 'black',
		       -fontcolor2  => 'black',
		       -key	=> 'SNPs',
		       -bump	=> +1,
		       -height	=> 8,
		       -label	=> \&snp_label,
		       #-glyph	=> 'triangle',
		       -glyph	=> 'diamond',
		       -description		=> \&feature_desc,
		       );

    $panel->add_track( generic	=> $feature_classes{'Excluded'},
		       -bgcolor	=> 'silver',
		       -fgcolor	=> 'black',
		       -fontcolor  => 'black',
		       -fontcolor2  => 'black',
		       -key	=> 'Excluded Regions',
		       -bump	=> +1,
		       -height	=> 6,
		       -label	=> \&feature_label,
		       -description		=> \&feature_desc,
		       );

    $panel->add_track( 
		       generic => $cleanpairs,
		       -bgcolor	=> \&feature_color,
		       -fgcolor	=> \&feature_color,
		       -fontcolor  => 'black',
		       -fontcolor2  => 'black',
		       -key	=> 'SiRNA Reagents',
		       -bump	=> +1,
		       -height	=> 8,
		       -label	=> \&feature_label,
		       -glyph	=> 'generic',
		       -description		=> \&feature_desc,
		       );

    my $gd = $panel->gd;
    my $black = $gd->colorAllocate(0,0,0);
    my $txt = GD::Text::Align->new($gd);    
    $txt->set( valign => 'center', align => 'center', color => $black);
    #$txt->set_font(['/usr/share/fonts/truetype/VERDANA.TTF',gdGiantFont ], 10);
    $txt->set_font(gdGiantFont);
    $txt->set_text("RNAi Reagents for ".$gene->accession );
    $txt->draw(200, 50, 0);

    my $pngfile = $TMPDIR . $gene->accession . '.png';
    my $pngurl = $TMPURL . $gene->accession . '.png';
    open (IMG, ">$pngfile") or die $!;
    binmode IMG;
    print IMG $gd->png;
    close IMG;

    # also get the imagemap boxes
    my @pairboxes = extract_pairs($panel->boxes);

    print $q->img({-src => $pngurl, -usemap=>"#MAP"});
    print $q->p('Oligos are color coded: rank 1 in ', 
		$q->font({-color => $COLORRANKS{1}}, $COLORRANKS{1}),
		', rank 2 in ',
		$q->font({-color => $COLORRANKS{2}}, $COLORRANKS{2}),
		' and rank 3 in ',
		$q->font({-color => $COLORRANKS{3}}, $COLORRANKS{3}),
		'. Click on an oligo to bring it up in the table below');

    print_imagemap(@pairboxes);

}

sub feature_label {
    my ($feature) = @_;
    my (@notes, @label);
    #$label = ucfirst($feature->primary_tag);
    foreach (qw(note name product gene)) {
	if ($feature->has_tag($_)) {
	    @notes = $feature->each_tag_value($_);
	    #$label .= ': ' . $notes[0];
	    push(@label, $notes[0]);
	    last;
	}
    }
    return join(': ', @label);
    #return $label;
}

sub feature_color {
    my ($feature) = @_;
    my ($rank) = $feature->each_tag_value('rank');
    #print STDERR "Feature rank: $rank COLOR $COLORRANKS{$rank}\n";
    return $COLORRANKS{$rank};
    #return 'red';
}


sub print_table {
    my ($accession, $pairs) = @_;

    print $q->h2("RNAi Reagents for $accession");
    print $q->start_table({-border => 1, -cellpadding => 2});
    print $q->TR( $q->th(['Reagent #', 'Start', 'Stop', 'Rank', 'Fxn GC', 'Sense Oligo', 'Antisense Oligo', 'Target' ]) ), "\n";


    my $i = 1;

    foreach $pair ( sort { $a->start <=> $b->start } @$pairs ) {
	my $sense = $pair->sense;
	my $anti = $pair->antisense;
	my $color = feature_color($pair);

#  	my $blasturl = "http://nunu.hpw.pri.bms.com/biocgi/versablast.pl?p=blastn&sequence=";
#  	$blasturl .= $pair->seq->seq;
#  	$blasturl .= "&action=Nucleotide Databases";

	print 
	    $q->TR( $q->td( [ $q->a({-name => 'RNAi' . $pair->start}) . $i,
			      $pair->start, 
			      $pair->end,
			      $q->font({-color => $color},$pair->rank), 
			      $pair->fxGC,
			      $q->tt($sense->seq), 	
			      $q->tt($anti->seq),
			      $q->tt($pair->seq->seq),
#  			      $q->a({-href=>$blasturl,
#  				     -target=>"blastn"},
#  				    "BLAST this target"),
			      ] ) ),
	"\n";
	$i++;
    }
    print $q->end_table;
}





sub print_text {
    my ($accession,  $pairs ) = @_;
    my ($pair);

    print "RNAi reagents for $accession \n";

    print join("\t", qw(Start Stop Rank Sense Antisense)), "\n";
    foreach $pair (@$pairs ) {
	my $sense = $pair->sense;
	my $anti = $pair->antisense;

	print join("\t", $pair->start, $pair->end, $pair->rank, $sense->seq, $anti->seq), "\n";


    }


}

sub cleanup_feature {
    my ($flist) = @_;

    my ($feat, @clean, $cfeat);

    foreach $feat(@$flist) {
	$cfeat = clone($feat);
#	$cfeat = $feat->clone;
	$cfeat->flush_sub_SeqFeature;
	push (@clean, $cfeat); # will they 
    }
    return \@clean;
}


sub extract_pairs {
    # get SiRNA::Pair features ONLY for imagemap
    return ( grep {ref($_->[0]) eq "Bio::SeqFeature::SiRNA::Pair"} @_ );
}

sub print_imagemap {
    my @items = @_;

    print q(<MAP NAME="MAP">), "\n";

    my $i = 1;
    
    foreach $item (@items) {
	my ($feature, $x1, $y1, $x2, $y2) = @$item;
	my $fstart = $feature->start; # should be unique
	my $text = 'RNAi #' . $i. ' Start=' . $feature->start . ' Rank='.$feature->rank;
	print qq(<AREA SHAPE="RECT" COORDS="$x1,$y1,$x2,$y2" TITLE="$text" HREF="#RNAi$fstart">), "\n";
	warn "Mouseover text: $text\n";

	$i++;
    }
    print "</MAP>\n";
}


sub print_error {
    # print error messages in big red type. Provide more graceful die/warn to end user
    my ($msg, $fatal) = @_;
    print $q->h3($q->font({-color=>'RED'}, $msg));
    
    if ($fatal) {
	print $q->end_html;
	die "$msg \n";
    }
    else {
	warn $msg;
    }
}

sub dump {
    print $q->start_ul;

    foreach ($q->param) {
	print $q->li($_),
	$q->ul($q->li([ $q->param($_) ]));
    }
}
    
sub snp_label {
    # special format for SNPs
    my ($feature) = @_;
    my $label;

    if ( $feature->has_tag('db_xref') ) {
	my @notes = $feature->each_tag_value('db_xref');
	$label .= $notes[0];
	$label .= ' ';
    }
    if ( $feature->has_tag('allele') ) {
	my ($nt1, $nt2) = $feature->each_tag_value('allele');
	$label .=  $nt1 . '->' . $nt2;
    }
    return $label;
}

sub feature_desc {
    my ($feature) = @_;
    my $desc = $feature->start;
    $desc .= '-' . $feature->end unless ($feature->start == $feature->end);
    return $desc;
}
