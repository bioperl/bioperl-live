use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME                => 'Deobfuscator',
    AUTHOR              => 'Dave Messina <dave-pause@davemessina.net>',
    VERSION_FROM        => 'lib/Deobfuscator.pm',
    ABSTRACT_FROM       => 'lib/Deobfuscator.pm',
    EXE_FILES           => ['bin/deob_index.pl',],
    PL_FILES            => {},
    PREREQ_PM => {
        'Test::More'       => 0,
        'version'          => 0,
		'Class::Inspector' => 0,
		'DB_File'          => 0,
        'CGI'              => 0,
    },
    dist                => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean               => { FILES => 'Deobfuscator-*' },
);
