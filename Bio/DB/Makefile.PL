use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'		=> 'Bio::DB::GFF',
    'VERSION_FROM'	=> 'GFF.pm', # finds $VERSION
    'PREREQ_PM'		=> { 'DBI'            => 1.0,
			     'Bio::DB::Fasta' => 1.0,
			   }, # e.g., Module::Name => 1.1
    'LIBS'		=> [''], # e.g., '-lm'
    'DEFINE'		=> '', # e.g., '-DHAVE_SOMETHING'
    'INC'		=> '', # e.g., '-I/usr/include/other'
);
