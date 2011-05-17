#!usr/bin/perl
use strict;

use Bio::Root::Test;
test_begin( -tests => 41,
            -requires_modules => [qw(GD)]); 

#Check if module and all its methods can be loaded
use_ok('Bio::Align::Graphics');
require_ok('Bio::Align::Graphics');
can_ok('Bio::Align::Graphics', qw(new draw height width aln_length aln_format no_sequences));
    
 
#Get an alignment file
my $file = Bio::Root::IO->catfile("t","data","pep-266.aln");
ok($file, 'input is defined');
	
#Create an AlignI object using AlignIO
my $in=new Bio::AlignIO(-file=>$file, -format=>'clustalw');
ok(defined $in, 'AlignIO object is defined');
isa_ok($in, 'Bio::AlignIO');

#Read the alignment
my $aln=$in->next_aln();
ok(defined $aln, 'alignment is there and defined');

#Create some domains for highlighting
my @domain_start = ( 25, 50, 80 );
my @domain_end = ( 40 , 65 , 100 );
my @domain_color = ( 'red' , 'cyan' , 'green' );
ok(exists $domain_start[2], 'all starts are present');
ok(exists $domain_end[2],'all ends are present');
ok(exists $domain_color[2], 'all colors are present');
cmp_ok($domain_start[0], '<=', $domain_end[0],'first end is further than first start');   #Some  
cmp_ok($domain_start[1], '<=', $domain_end[1],'second end is further than second start');	#logical
cmp_ok($domain_start[2], '<=', $domain_end[2],'third end is further than third start');   #tests 

#Create Labels for the domains
my @dml = ("CARD", "Proline Rich", "Transmembrane");
my @dml_start = (25, 50, 80);
my @dml_end = (40, 65, 100);
my @dml_color = ("lightpink", "lightblue", "lightgreen");
ok(exists $dml[2], 'domain labels are present');
ok(exists $dml_start[2], 'domain starts are present');
ok(exists $dml_end[2], 'domain ends are present');
ok(exists $dml_color[2], 'domain colors are present');

#Some logical tests 
cmp_ok($dml_start[0], '<=', $dml_end[0],'label - first end is further than first start');         
cmp_ok($dml_start[1], '<=', $dml_end[1],'label - second end is further than second start');	
cmp_ok($dml_start[2], '<=', $dml_end[2],'label - third end is further than third start');          
cmp_ok($domain_start[0], '>=', $dml_start[0],'first label start is within domain range');        
cmp_ok($domain_start[1], '>=',$dml_start[1],'second label start is within domain range');        
cmp_ok($domain_start[2], '>=',$dml_start[2],'third label start is within domain range');        
cmp_ok($domain_end[0], '>=',$dml_end[0],'first label end is within domain range');        
cmp_ok($domain_end[1], '>=',$dml_end[1],'second label end is within domain range');        
cmp_ok($domain_end[2], '>=',$dml_end[2],'third label end is within domain range');        
 
#Create individual labels
my %labels = ( 145 => "Hep-c target");
ok(exists $labels{145}, 'individual labels work');

#my $output_file = test_output_file();
 
my $print_align = Bio::Align::Graphics->new( align => $aln,
					pad_bottom => 5,
					dm_start => \@domain_start,
					dm_end => \@domain_end,
					dm_color => \@domain_color,
					dm_labels => \@dml,
					dml_start => \@dml_start,
					dml_end => \@dml_end,
					dml_color => \@dml_color,
				 	labels => \%labels,
					out_format => "png",
#					output=>$output_file,
					wrap=>80);

isa_ok($print_align, 'Bio::Align::Graphics');
ok( defined $print_align, 'new object is defined');
is($print_align->{pad_bottom}, 5, '  pad_bottom is right');
is($print_align->{pad_top}, 5, '  default pad_top is right');
is_deeply($print_align->{domain_start}, \@domain_start,'  start point loaded');
is_deeply($print_align->{domain_end}, \@domain_end,'  end point loaded');
is_deeply($print_align->{domain_color}, \@domain_color,'  color of domain loaded');
is_deeply($print_align->{dm_labels}, \@dml, '  domain labels loaded');
is_deeply($print_align->{dm_label_start}, \@dml_start, '  label starts loaded');
is_deeply($print_align->{dm_label_end}, \@dml_end, '  label ends loaded');
is_deeply($print_align->{dm_label_color}, \@dml_color, '  label colors loaded');
is_deeply($print_align->{labels}, \%labels, '  labels loaded');
is($print_align->{out_format}, 'png', '  output file is png');
isnt($print_align->{wrapping}, 0, '  wrapping length is not zero');

exit;
