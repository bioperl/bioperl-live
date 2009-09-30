#$Id$
package Bio::DB::AnyDBMImporter;
use strict;

=head1 Bio::DB::AnyDBMImporter - A hack to import DBM package symbols when using AnyDBM_File

=head1 SYNOPSIS

 BEGIN {
    @AnyDBM_File::ISA = qw( DB_File SDBM_File );
 }
 use AnyDBM_File;
 use Bio::DB::AnyDBMImporter; 
 use Fcntl; # often required; do if you get an error involving 'O_SVWST'

 my %db;
 $DB_BTREE->{'flags'} = R_DUP;
 tie( %db, 'AnyDBM_File', O_CREAT | O_RDWR, 0644, $DB_BTREE);

=head1 DESCRIPTION

This is a hack that allows symbols (like $DB_HASH, R_DUP, etc.) to be
imported into the caller's namespace when using the L<AnyDBM_File> DBM
auto-selection package. L<AnyDBM_File> includes its auto-selected module
by using C<require>, which unlike C<use> does not export symbols in
the required packages C<@EXPORT> array.

The hack is so-called because it relies on L<AnyDBM_File> internal
behavior. Specifically, at the time of DBM module selection,
C<AnyDBM_File> sets its C<@ISA> to a length 1 array containing the
package name of the selected DBM module.

=head1 AUTHOR - Mark A. Jensen

 Email: maj -at- fortinbras -dot- us

=cut

my ($pkg, @goob) = caller;
if (!@AnyDBM_File::ISA) {
    die "No packages specified for AnyDBM_File (have you forgotten to include AnyDBM_File?)"
}
elsif (@AnyDBM_File::ISA > 1) {
    warn "AnyDBM_File has not selected a single DBM package; returning..."
}
else {
    my @export = eval "\@$AnyDBM_File::ISA[0]::EXPORT";
    for (@export) {
	m/^\$(.*)/ && do {
	    eval "\$${pkg}::$1 = \$$AnyDBM_File::ISA[0]\::$1";
	};
	m/^\@(.*)/ && do {
	    eval "\@${pkg}::$1 = \@$AnyDBM_File::ISA[0]\::$1";
	};
	m/^\%(.*)/ && do {
	    eval "\%${pkg}::$1 = \%$AnyDBM_File::ISA[0]\::$1";
	};
	m/^[^\$@%]/ && do {
	    eval "*{${pkg}::$_\} = \\\&$AnyDBM_File::ISA[0]\::$_";
	};
	die $@ if $@;
    }
}
1;
