# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 26,
			   -requires_modules => ['Text::Wrap 98', 'XML::Writer']);
	
	use_ok('Bio::Variation::IO');
}

sub io {
    my ($t_file, $o_file, $out_format) = @_; 
    my $res;
	
    my ($o_ext) = $out_format eq 'flat' ? 'dat' : 'xml';
    my ($o_format) = $out_format;
    my ($t_name) = $t_file =~ /(.*)....$/;
	
    my( $before );
    {
        local $/ = undef;
        open my $BEFORE, '<', "$t_name.$o_ext" or die "Could not read file '$t_name.$o_ext': $!\n";
        $before = <$BEFORE>;
        close $BEFORE;
    }

    ok $before;#"Error in reading input file [$t_name.$o_ext]";

    my $in = Bio::Variation::IO->new( -file => $t_file);
    my  @entries ;
    while (my $e = $in->next) {
        push @entries, $e;
    }
    my $count = scalar @entries;
    cmp_ok @entries, '>', 0;# "No SeqDiff objects [$count]";

    my $out = Bio::Variation::IO->new( -FILE => ">$o_file", 
				       -FORMAT => $o_format);
    my $out_ok = 1;
    foreach my $e (@entries) {
        $out->write($e) or $out_ok = 0;
    }
    undef($out);  # Flush to disk
    ok $out_ok;#  "error writing variants";

    my( $after );
    {
        local $/ = undef;
        open my $AFTER, '<', $o_file or die "Could not read file '$o_file': $!\n";
        $after = <$AFTER>;
        close $AFTER;
    }

    ok $after;# "Error in reading in again the output file [$o_file]";
    is $before, $after, "test output file compared to input";
    print STDERR `diff $t_file $o_file` if $before ne $after;
}

io  (test_input_file('mutations.dat'), 
     test_output_file(), 'flat'); #1..5
io  (test_input_file('polymorphism.dat'), 
     test_output_file(), 'flat'); #6..10

SKIP: {
	test_skip(-tests => 15, -requires_modules => [qw(XML::Twig
												     XML::Writer
												     IO::String)]);

	eval {
		if( $XML::Writer::VERSION >= 0.5 ) { 
		io  (test_input_file('mutations.xml'), 
			 test_output_file(), 'xml'); #10..12
		} else { 
		io  (test_input_file('mutations.old.xml'), 
			 test_output_file(), 'xml'); #10..12
		}
	};
	
	eval {
		if( $XML::Writer::VERSION >= 0.5 ) { 
		io  (test_input_file('polymorphism.xml'), 
			 test_output_file(), 'xml'); #13..14
		} else { 
		io  (test_input_file('polymorphism.old.xml'), 
			 test_output_file(), 'xml'); #13..14
	
		}
	};
	
	eval { 
		if( $XML::Writer::VERSION >= 0.5 ) { 
		io  (test_input_file('mutations.dat'), 
			 test_output_file(), 'xml'); #15..25
		} else { 
		io  (test_input_file('mutations.old.dat'), 
			 test_output_file(), 'xml'); #15..25
		}
	};
}
