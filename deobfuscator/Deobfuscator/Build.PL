use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Deobfuscator',
    license             => 'perl',
    dist_author         => 'Dave Messina <dave-pause@davemessina.net>',
    dist_version_from   => 'lib/Deobfuscator.pm',
    script_files        => ['bin/deob_index.pl',],
    requires => {
        'Test::More'       => 0,
        'version'          => 0,
		'Class::Inspector' => 0,
		'DB_File'          => 0,
        'CGI'              => 0,
    },
    add_to_cleanup      => [ 'Deobfuscator-*' ],
);

$builder->create_build_script();
