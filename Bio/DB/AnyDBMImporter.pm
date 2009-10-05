#$Id$
package Bio::DB::AnyDBMImporter;
use strict;

=head1 Bio::DB::AnyDBMImporter - A hack to import DBM package symbols when using AnyDBM_File

=head1 SYNOPSIS

 BEGIN {
    @AnyDBM_File::ISA = qw( DB_File SDBM_File ) unless @AnyDBM_File::ISA;
 }
 use AnyDBM_File;
 use vars qw( $DB_BTREE &R_DUP); # must declare the globals you expect to use
 use Bio::DB::AnyDBMImporter qw(:bdb); # an import tag is REQUIRED

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

=head1 USAGE NOTES

Use of AnyDBMImporter within module code requires a kludge.
Symbols of imported variables or constants need to
be declared globals, as in the SYNOPSIS above. This is not necessary
when AnyDBMImporter is used in package main.

AnyDBMImporter consists entirely of an import function. To import the
symbols, a tag must be given. More than one tag can be supplied. Symbols
cannot be individually specified at the moment.

 :bdb    DB_File (BDB) symbols ($DB_*, R_*, O_*)
 :db     $DB_* type hashrefs
 :R      R_* constants (R_DUP, R_FIRST, etc)
 :O      O_* constants (O_CREAT, O_RDWR, etc)
 :other  Exportable symbols not in any of the above groups
 :all    All exportable symbols
 
=head1 AUTHOR - Mark A. Jensen

 Email: maj -at- fortinbras -dot- us

=cut

use constant { R_CONST => 1, O_CONST => 2, DB_TYPES => 4, OTHER => 8 };

sub import {
    my ($class, @args) = @_;
    my ($pkg, $fn, $ln) = caller;
    my $flags = 0;
    for (@args) {
	/^:all$/ && do {
	    $flags |= (R_CONST | O_CONST |  DB_TYPES | OTHER );
	    next;
	};
	/^:other$/ && do {
	    $flags |= OTHER;
	    next;
	};
	/^:bdb/ && do {
	    $flags |= (R_CONST | O_CONST |  DB_TYPES );
	    next;
	};
	/^:db$/ && do {
	    $flags |= DB_TYPES;
	    next;
	};
	/^:R$/ && do {
	    $flags |= R_CONST;
	    next;
	};
	/^:O$/ && do {
	    $flags |= O_CONST;
	    next;
	};
	do {
	    die "Tag '$_' not recognized";
	};
    }
    die "No symbols exported" unless $flags;
    
    if (!@AnyDBM_File::ISA) {
	die "No packages specified for AnyDBM_File (have you forgotten to include AnyDBM_File?)"
    }
    elsif (@AnyDBM_File::ISA > 1) {
	warn "AnyDBM_File has not yet selected a single DBM package; returning..."
    }
    else {
	my @export = eval "(\@$AnyDBM_File::ISA[0]::EXPORT, \@$AnyDBM_File::ISA[0]::EXPORT_OK)";
	my $ref;
	for (@export) {
	    # kludge: ignore gnu perl undefined symbols
	    my $qm = quotemeta;
	    next if grep(/^$qm$/, qw(R_NOKEY R_SNAPSHOT O_ALIAS O_ASYNC O_DEFER O_DIRECTORY O_EXLOCK O_LARGEFILE O_RANDOM O_RAW O_RSRC O_SEQUENTIAL O_SHLOCK O_TEMPORARY ));
	    m/^\$(.*)/ && do {
		$_ = substr $_, 1;
		eval "\$ref = *${pkg}::$_\{SCALAR}";
		die $@ if $@;
		if ( ($flags & DB_TYPES and ($1 =~ /^DB_/)) ||
		     ($flags & OTHER and ($1 !~ /^DB_/)) ) {
		    $$ref = eval "\$$AnyDBM_File::ISA[0]\::$_";
		}
		next;
	    };
	    m/^\@(.*)/ && do {
		$_ = substr $_, 1;
		eval "\$ref = *${pkg}::$_\{ARRAY}";
		die $@ if $@;
		if  ($flags & OTHER) {
		    $$ref = eval "\@$AnyDBM_File::ISA[0]\::$1";
		}
		next;
	    };
	    m/^\%(.*)/ && do {
		$_ = substr $_, 1;
		eval "\$ref = *${pkg}::$_\{HASH}";
		die $@ if $@;
		if  ($flags & OTHER) {
		    $$ref = eval "\%$AnyDBM_File::ISA[0]\::$1";
		}
		next;
	    };
	    m/^[^\$@%]/ && do {
		eval "*{${pkg}::$_} = \\\&$AnyDBM_File::ISA[0]\::$_" if 
		   ( ($flags & R_CONST and /^R_/) ||
		    ($flags & O_CONST and /^O_/) ||
		     ($flags & OTHER and /^[RO]_/) );

		next;
	    };
	}
    }
}

1;
