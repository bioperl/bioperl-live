#!/usr/local/bin/perl


=head1 NAME

gff2ps - you will want to change this script

=head2 SYNOPSIS

   perl gff2ps < file.gff > file.ps

=head2 DESCRIPTION

This script provides GFF to postscript handling. Due to the ... ummm
... potential for flexible reinterpretation that is GFF, this script
will almost certainly need modifying for anyone elses use (basically,
you need to know what you want to get out of the GFF file and how to
draw it). But it does include code to draw the most challenging thing
out there - genes - and should give you a good example of where to
start

=head2 AUTHOR

Ewan Birney

=cut


use Bio::Tools::GFF;

my $font      = 8;
my $scale     = 200;
my $rotate    = 1;
my $feature_off = 0;

use Getopt::Long;

&GetOptions(
	    "scale=i"   => \$scale,
	    "font=i"    => \$font,
	    "rotate=i"  => \$rotate,
	    "start=i"   => \$feature_off
	    );


my $gffio = Bio::Tools::GFF->new(-fh => \*STDIN, -gff_version => 1);
my $feature;

use Data::Dumper;

my %set;

# loop over the input stream
while( my $f = $gffio->next_feature()) {
    $f->start($f->start - $feature_off);
    $f->end  ($f->end   - $feature_off);

    if( $f->start < 0 ) {
	next;
    }

    if( $f->start > $scale*1000 ) {
      next;
    }
    
    
    if( $f->primary_tag =~ /coding_exon/ ) {
	#print STDERR "Seen ",$f->start," ",$f->end,"\n";
	($group) = $f->each_tag_value('group');
	$group =~ s/\s+//g;

	

	if( !defined $set{$group} ) {
	    $set{$group} = Bio::SeqFeature::Generic->new();
	    $set{$group}->seqname($f->seqname);
	    $set{$group}->primary_tag('transcript');
	    $set{$group}->source_tag($f->source_tag);
	    $set{$group}->add_tag_value('id',$group);
	}
	$set{$group}->add_sub_SeqFeature($f,'EXPAND');
	$set{$group}->strand($f->strand);
    }
}
$gffio->close();


#foreach my $set ( values %set ) {
#    print $set->gff_string,"\n";
#    foreach $sub ( $set->sub_SeqFeature ) {
#	print $sub->gff_string,"\n";
#    }
#}


# sort into forward and reverse strands

my @forward;
my @reverse;

$max = 0;

foreach my $set ( values %set ) {
    if( $set->end > $max ) {
	$max = $set->end;
    }

    if( $set->strand == -1 ) {
	push(@reverse,$set);
    } else {
	push(@forward,$set);
    }
}

@forward = sort { $a->start <=> $b->start } @forward;
@reverse = sort { $a->start <=> $b->start } @reverse;

&print_header(\*STDOUT);

if( $rotate ) {
   print "0 700 translate\n";
   print "-90 rotate\n";
}

print "0 200 moveto 900 200 lineto stroke\n";

my $bp_max = $scale*900;

for(my $bp = 0;$bp < $bp_max ;$bp = $bp + 5000) {
    print STDOUT $bp/$scale," 200 moveto ",$bp/$scale," 197 lineto\n";
    $text = int( $feature_off + ($bp/1000));
    print STDOUT $bp/$scale," 195 moveto ($text) show\n";
}


&draw_gene(\@forward,1,$scale,220,\*STDOUT);
&draw_gene(\@reverse,-1,$scale,180,\*STDOUT);

print "showpage\n";




sub draw_gene {
    my ($gene_array,$strand,$scale,$offset,$fh) = @_;


    my @bump_array;
    my $bump_row_max = 1;
    my $bump_end = int $max/$scale;

    $bump_array[0] = '0' x $bump_end;


    foreach my $f ( @$gene_array ) {

	#
	# Bump it baby!
	#

	# We keep an array of strings for currently draw areas. Do this in pixel
	# coordinates to save on space. If the region has all 0's then we know we
	# can draw here. If so, we set it to 1's. If not, we go up a row and see if
	# we can fit it in there. If we exhausted the rows we make a new row. 

	$bump_start = (int $f->start/$scale)-1;
	$bump_len   = int(($f->end - $f->start)/$scale) +1;

	# text might be longer than gene. Mystic number 5 looks good for 8 point helvetica
	# you will have to change this otherwise.

	my ($gene_id) = $f->each_tag_value('id');
	if( (length $gene_id)*5 > $bump_len ) {
	    $bump_len = (length $gene_id)*5; 
	}

	# figure out the first place to fit in this gene;
	for($i=0;$i<$bump_row_max;$i++) {
	    #print STDERR "Seeing $bump_start $bump_len $i ",substr($bump_array[$i],$bump_start,$bump_len),"\n";
	    
	    if( substr($bump_array[$i],$bump_start,$bump_len) !~ /1/ ) {
		#print STDERR "Going to break with $i\n";
		last;
	    }
	}
	#print STDERR "i is $i\n";
	# if $i == bump_row_max then we need a new bump row
	if( $i == $bump_row_max ) {
	    $bump_array[$bump_row_max] = '0' x $bump_end;
	    $bump_row_max++;
	}
	
	# now blank out this bump row to 1's
	
	substr($bump_array[$i],$bump_start,$bump_len) = '1' x $bump_len;
	
	# now print it out ;)

	
	#
	# Need to be portable between strands. Gene hats go the
	# other way up on reverse strand, but not the text. 
	#

	if( $strand == 1 ) {
	    $text   = $offset+($i*20)+1;
	    $bottom = $offset+($i*20)+10;
	    $top    = $offset+($i*20)+20;
	    $mid    = $offset+($i*20)+15;
	} else {
	    $text   = $offset-($i*20)-19;
	    $bottom = $offset-($i*20);
	    $top    = $offset-($i*20)-10;
	    $mid    = $offset-($i*20)-5;
	}

	print $fh $f->start/$scale," ",$text," moveto\n";
	print $fh "($gene_id) show\n";

	my $prev = undef;
	    
	foreach $exon ( $f->sub_SeqFeature ) {
	    print $fh $exon->start/$scale," ",$bottom," moveto\n";
	    print $fh $exon->end/$scale," ",$bottom," lineto\n";
	    print $fh $exon->end/$scale," ",$top, " lineto\n";
	    print $fh $exon->start/$scale," ",$top," lineto\n";
	    print $fh "closepath stroke\n";

	    # draw the intron hat
	    if( defined $prev ) {
		print $prev->end/$scale," ",$mid," moveto\n";
		my $intron_len = $exon->start - $prev->end;
		
		print $fh ($prev->end+($intron_len/2))/$scale," ",$top," lineto\n";
		print $fh $exon->start/$scale," ",$mid," lineto stroke\n";
	    }
		
	    $prev = $exon;
	}
	
    }
}


sub print_header {
    my $fh = shift;

    print $fh <<EOF;
%!PS-Adobe-2.0
% Created by Genome2ps. Ewan Birney <birney\@ebi.ac.uk>
0.5 setlinewidth
/Helvetica findfont $font scalefont setfont
EOF

}
