# -*-Perl-*- Test Harness script for Bioperl
# $Id: Sim4.t 11525 2007-06-27 10:16:38Z sendu $

use strict;
BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 102);
	
	use_ok('Bio::SearchIO');
}

# parse align format 0
my $parser = Bio::SearchIO->new(-format => 'sim4',
			       -file   => test_input_file('crypto.sim4-0')
			       );
my $r = $parser->next_result;
is ($r->query_name, 'cn416');
is ($r->query_length, 630);

my $hit = $r->next_hit;
is ($hit->name, 'Contig147');
is ($hit->description, 'Contig147.fa');
is ($hit->length, 1086);

my $hsp = $hit->next_hsp;
is ($hsp->query->start, 36);
is ($hsp->query->end, 132);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 191);
is ($hsp->hit->end, 286);
is ($hsp->hit->strand, 1);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 133);
is ($hsp->query->end, 191);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 343);
is ($hsp->hit->end, 401);
is ($hsp->hit->strand, 1);

# parse align format 3
$parser = Bio::SearchIO->new(-format => 'sim4',
			    -file   => test_input_file('crypto.sim4-3')
			    );
$r = $parser->next_result;
is ($r->query_name, 'cn416');
is ($r->query_length, 630);
$hit = $r->next_hit;
is ($hit->name, 'Contig147');
is ($hit->description, 'Contig147.fa');
is ($hit->length, 1086);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 36);
is ($hsp->query->end, 132);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 191);
is ($hsp->hit->end, 286);
is ($hsp->hit->strand, 1);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 133);
is ($hsp->query->end, 191);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 343);
is ($hsp->hit->end, 401);
is ($hsp->hit->strand, 1);

# parse align format 4
$parser = Bio::SearchIO->new(-format => 'sim4',
			    -file   => test_input_file('crypto.sim4-4')
			    );
$r = $parser->next_result;
is ($r->query_name, 'cn416');
is ($r->query_length, 630);

$hit = $r->next_hit;
is ($hit->name, 'Contig147');
is ($hit->length, 1086);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 36);
is ($hsp->query->end, 132);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 191);
is ($hsp->hit->end, 286);
is ($hsp->hit->strand, 1);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 133);
is ($hsp->query->end, 191);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 343);
is ($hsp->hit->end, 401);
is ($hsp->hit->strand, 1);


# do the other sim4 files
$parser = Bio::SearchIO->new(-format => 'sim4',
			    -file   => test_input_file('sim4.rev')
			    );
$r = $parser->next_result;
is ($r->query_name, '/nfs/disk21/birney/prog/wise2/example/human.rev');
is ($r->query_length, 5368);
$hit = $r->next_hit;
is ($hit->name, 'HSHNCPA1');
is ($hit->description, 'temp.cdna');
is ($hit->length, 1198);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 486);
is ($hsp->query->end, 503);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 10);
is ($hsp->hit->end, 27);
is ($hsp->hit->strand, -1);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 1048);
is ($hsp->query->end, 1117);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 194);
is ($hsp->hit->end, 265);
is ($hsp->hit->strand, -1);

# do the other sim4 files fwd
$parser = Bio::SearchIO->new(-format => 'sim4',
			    -file   => test_input_file('sim4.for.for')
			    );
$r = $parser->next_result;
is ($r->query_name, 'human.genomic');
is ($r->query_length, 5368);
$hit = $r->next_hit;
is ($hit->name, 'hs_est');
is ($hit->description, 'est.for');
is ($hit->length, 479);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 695);
is ($hsp->query->end, 813);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 1);
is ($hsp->hit->end, 119);
is ($hsp->hit->strand, 1);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 1377);
is ($hsp->query->end, 1500);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 120);
is ($hsp->hit->end, 243);
is ($hsp->hit->strand, 1);

# do the other sim4 files fwd rev
$parser = Bio::SearchIO->new(-format => 'sim4',
			    -file   => test_input_file('sim4.for.rev')
			    );
$r = $parser->next_result;
is ($r->query_name, 'human.genomic');
is ($r->query_length, 5368);
$hit = $r->next_hit;
is ($hit->name, 'REVCOMP');
is ($hit->description, 'hn_est.rev');
is ($hit->length, 479);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 695);
is ($hsp->query->end, 813);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 1);
is ($hsp->hit->end, 119);
is ($hsp->hit->strand, -1);

$hsp = $hit->next_hsp;
is ($hsp->query->start, 1377);
is ($hsp->query->end, 1500);
is ($hsp->query->strand, 1);
is ($hsp->hit->start, 120);
is ($hsp->hit->end, 243);
is ($hsp->hit->strand, -1);
