package Dummy::Config;
use strict;
use warnings;
no warnings qw(qw);
use Exporter;
our (@ISA, @EXPORT, @EXPORT_OK);
push @ISA, 'Exporter';
@EXPORT = qw(
             $program_name
             $program_dir
             $use_dash
             $join
             @program_commands
             %command_prefixes
             @program_params
             @program_switches
             %param_translation
             %command_files
             %composite_commands
            );

@EXPORT_OK = qw();

our $program_name = '*flurb';
#our $program_dir = 'C:\cygwin\usr\local\bin';
our $use_dash = 'mixed';
our $join = ' ';


our @program_commands = qw(
 rpsblast
 find
 goob
 blorb
 multiglob
  );

our %command_prefixes = (
    blastp => 'blp',
    tblastn => 'tbn',
    goob => 'g',
    blorb => 'b',
    multiglob => 'm'
    );

our @program_params = qw(
    command
    g|narf
    g|schlurb
    b|scroob
    b|frelb
    m|trud
);

our @program_switches = qw(
    g|freen
    b|klep
);

our %param_translation = (
    'g|narf'     => 'n',
    'g|schlurb'  => 'schlurb',
    'g|freen'    => 'f',
    'b|scroob'   => 's',
    'b|frelb'    => 'frelb'
    );

our %command_files = (
    'goob'       => [qw( fas faq )],
    );

our %composite_commands = (
    'multiglob' => [qw( blorb goob )]
    );
1;
