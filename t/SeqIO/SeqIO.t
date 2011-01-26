# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN {     
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 5);
	
	use_ok('Bio::SeqIO');
}

# simple tests specific to Bio::SeqIO interface (applicable to all SeqIO
# modules)

############ EXCEPTION HANDLING ############

throws_ok {
    Bio::SeqIO->new();
} qr/No file or fh argument provided/, 'Must pass a file or file handle';

throws_ok {
    Bio::SeqIO->new(-fh => undef);
} qr/fh argument provided, but with an undefined value/,
    'Must pass a file or file handle';

throws_ok {
    Bio::SeqIO->new(-file => undef);
} qr/file argument provided, but with an undefined value/,
    'Must pass a file or file handle';

throws_ok {
    Bio::SeqIO->new(-file => 'foo.bar');
} qr/Can not open 'foo.bar' for reading: No such file or directory/,
    'Must pass a real file';

