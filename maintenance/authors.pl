#!/usr/bin/perl
#
=head1 NAME

authors.pl - check modules and scripts for authors not in AUTHORS file

=head1 SYNOPSIS

B<authors.pl> [B<-d|--dir> path ] [B<-v|--verbose>] B<-a|--authorsfile>
    [B<-?|-h|--help>]

=head1 DESCRIPTION

Checks Plain Old Documentation (POD) of all bioperl live modules for
AUTHORS and CONTRIBUTORS tags and prints out any emails missing from
the AUTHORS file

=cut

use Data::Dumper;
use File::Find;
use Getopt::Long;
use strict;

sub findauthors;

#
# command line options
#

my ($verbose, $dir, $authorsfile, $help) = (0, undef, "../AUTHORS", undef);
GetOptions(
           'v|verbose' => \$verbose,
           'dir:s' => \$dir,
           'authorsfile:s' => \$authorsfile,
	   'h|help|?' => sub{ exec('perldoc',$0); exit(0) }
	   );

#
# global variables
#

# known authors from the AUTHORS file are read into
# the hash which is initialized with known synonymes
our %AUTHORS = map {$_=>1} qw{
                              birney@sanger.ac.uk
                              jinsana@gmx.net
                              Insana@ebi.ac.uk
                              fugui@worf.fugu-sg.org
                              cjm@fruitfly.bdgp.berkeley.edu
                              elia@tll.org.sg
                              heikki-at-bioperl-dot-org
                              bioinformatics@dieselwurks.com
                              bioinformatics1@dieselwurks.com
                              bioperl-l@bio.perl.org
                              paul@systemsbiology.org
                              gattiker@isb-sib.ch
                              elia@fugu-sg.org
                              jason@cgt.mc.duke.edu
                              jason@chg.mc.duke.edu
                              jason@open-bio.org
                              hilmar.lapp@pharma.novartis.com
                              richard.adams@ed.ac.uk
                              dblock@gene.pbi.nrc.ca
                              ak@ebi.ac.uk
                              day@cshl.org
                              bala@tll.org.sg
                              mrp@sanger.ac.uk
                              m.w.e.j.fiers@plant.wag-ur.nl
                              cmzmasek@yahoo.com
                              fuguteam@fugu-sg.org
                              shawnh@gmx.net
                          };
our %NEWAUTHORS;     # new authors
our %FIND_OPTIONS = ( wanted => \&findauthors, no_chdir => 1 );


# Directories to check
my @dirs = qw( ../Bio/ ../scripts . );

#print Dumper \%AUTHORS;

#
# Read the AUTHORS file
#


open my $F, '<', $authorsfile or die "Could not read file '$authorsfile': $!\n";


while (my $line = <$F>) {
    my ($email) = ($line =~/([\.\w_-]+ at [\.\w_-]+)/);
    next unless $email;
    #print $email, "\n";
    $email =~ s/ at /@/;
    $AUTHORS{$email} = 1;
}
close $F;


#
# run
#

if ($dir) {
    find \%FIND_OPTIONS, $dir;
} else {
    find \%FIND_OPTIONS, @dirs;
}

#
# results
#
print Dumper \%NEWAUTHORS;


#
##
### end main
##
#

#
# this is where the action is
#
sub findauthors {
    my ($filename) = @_;
    return unless ($filename =~ /\.PLS$/ or $filename =~ /\.p[ml]$/);
    return unless -e $filename;

    print "$filename\n" if $verbose;
    #local $/=undef;
    open my $F, '<', $filename or die "Could not read file '$filename': $!\n";
    while (my $line = <$F>) {
        #print $line;
        last if $line =~ /=head1 +AUTHOR/;
    }
    my $authorblock;
    while (my $line = <$F>) {
        last if $line =~ /=head/ and $line !~ /CONTRIBUTORS/;
        $authorblock .= $line;
    }
    close $F;
    return unless $authorblock;

    while ( $authorblock =~ /([\.\w_-]+@[\.a-z_-]+)/g) {
        #my $email = $1;
        #$email =~ //
        next if $AUTHORS{$1};
        #print "$filename\t$1\n";

        push @{$NEWAUTHORS{$1}}, $filename;
    }
}



=head1 OPTIONS

=over 3

=item B<-d | --dir> path

Overides the default directories to check by one directory 'path' and
all its subdirectories.

=item B<-a | --authorsfile>

path from working directory the AUTHORS file.

Redundant as this information could be had from --dir option butI am
feeling too lazy to change the code.

=cut

sub blankline {
    my ($file) = @_;
    return unless ($file =~ /\.PLS$/ or $file =~ /\.p[ml]$/);
    return unless -e $file;

    open my $F, '<', $file or warn "Could not read file '$file': $!\n" && return;
    local $/ = "";
    while (my $line = <$F>) {
        print "$file: +|$1|\n"  if $line =~ /[ \t]\n(=[a-z][^\n]+$)/m and $verbose;
        print "$file: ++|$1|\n" if $line =~ /\w\n(=[a-z][^\n]+$)/m    and $verbose;
        print "$file:|$1|+\n"   if $line =~ /(^=[a-z][^\n]+)\n[\t ]/m;
        #print "$file:|$1|++\n" if /(^=[^\n]+)\n\w/m;
    }
    close $F;
}

__END__

=item B<-v | --verbose>

Show the progress through files during the POD checking.

=item B<-? | -h  | --help>

This help text.

=back

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Heikki Lehvaslaiho

Email heikki-at-bioperl-dot-org

=cut


# find . -name '*.pm' -print | xargs  perl -e '$/=""; while (<>) {$n = $1 if /^package\s+([\w:]+)/; print "$n:|$1|"  if  /(\s\s^=[^\n]+$)/m ; }'  ;

# find . -name '*.pm' -print | xargs  perl -e '$/=""; while (<>) {$n = $1 if /^package\s+([\w:]+)/; print "$n:|$1|\n"  if /(^=[^\n]+\n[\t ])/m ; }'  ;
