#!perl
# $Id: einfo.PLS 15088 2008-12-04 02:49:09Z bosborne $
use strict;
use warnings;

=head1 NAME 

einfo.pl - query einfo to find all available databases or information about a
            specific database (field information or links to other NCBI
            databases)

=head1 SYNOPSIS

 einfo [-d database] [-f Field Code] [-l Link Name] [-o outfile]

=head1 DESCRIPTION

Command line options:

  -e/--email
        Valid email (required by NCBI policy)

  -d/--db/--database
        NCBI database to query
        (default = none, which shows available databases)
        
  -f/--field
        print out information about a specific field code
        (default = none)
        
  -l/--link
        print out information about a specific link name
        (default = none)
        
  -o/--out
        outfile
        (default = STDOUT)
        
  -h/--help
        show this documentation

As per NCBI's policy regarding eutils access, a valid email is required.  This
is not enforced here (if one is provided you will get a standard warning), but
don't be surprised if this doesn't work after June 1, 2010 unless one is
supplied.

If -d is not specified, field and link arguments are ignored and all available
databases are printed instead.

If either link names or field codes (or both) are specified, nothing else is
printed out (only the info requested).  You can specify as many fields and/or
links as you want by using multiple -f/-l E<lt>ARGE<gt> on the command line.

=head1 AUTHOR - Chris Fields

Chris Fields cjfields at bioperl dot org

=cut

use Getopt::Long;
use Bio::DB::EUtilities;

my ($db, @fields, @links, $outfile, $email);
GetOptions(
       'e|email:s'          => \$email,
       'd|db|database:s'    => \$db,
       'f|field:s'          => \@fields,
       'l|link:s'           => \@links,
       'o|out|outfile:s'    => \$outfile,
       'h|help'             => sub { exec('perldoc',$0); exit; }
       );

my $outfh;
if( $outfile ) { 
    open($outfh, ">$outfile") || die("$outfile: $!");
} else { 
    $outfh = \*STDOUT; 
}

if (!defined $db) {
    my $eutil = Bio::DB::EUtilities->new(-eutil => 'einfo',
                                         -email => $email);
    print $outfh join("\n",$eutil->get_available_databases);
    exit;
} else {
    my $eutil = Bio::DB::EUtilities->new(-eutil => 'einfo',
                                         -db => $db,
                                         -email => $email);
    if (@links || @fields) {
        for my $fi ($eutil->get_FieldInfo) {
            my $code = $fi->get_field_code;
            if (grep {$_ eq $code} @fields) {
                print $outfh $fi->to_string."\n";
            }
        }
        for my $li ($eutil->get_LinkInfo) {
            my $nm = $li->get_link_name;
            if (grep {$_ eq $nm} @links) {
                print $outfh $li->to_string."\n";
            }
        }
    } else {
        $eutil->print_FieldInfo;
        $eutil->print_LinkInfo;
    }
}
