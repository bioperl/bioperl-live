# BioPerl module for Bio::Mapmaker
#
# Cared for by Jason Stajich <bioinformatics1@dieselwurks.com>
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Mapmaker - Parser for Mapmaker Linkage Map output files

=head1 SYNOPSIS

=head1 DESCRIPTION

Parse Mapmaker files.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  bioperl-bugs@bioperl.org
  http://bioperl.org/bioperl-bugs/

=head1 AUTHOR - Chad Matsalla

Email bioinformatics1@dieselwurks.com

=head1 CONTRIBUTORS

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::Map::Mapmaker;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Root

use Bio::Root::Root;
use Bio::Root::IO;
use Bio::Map::Marker;
use Bio::Map::OrderedPositionWithDistance;

@ISA = qw(Bio::Root::Root Bio::Root::IO);

=head2 new(-file => $filename)

 Title   : new(-file => $filename)
 Usage   : my $obj = new Bio::Map::Mapmaker(-file => $filename);
 Function: Parses a Mapmaker object
 Returns : A Bio::Map::Mapmaker object
 Args    : -file => the name of a file

=cut

sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->_initialize_io(@args);
    ($self->{'type'}) = $self->_rearrange([qw(TYPE)], @args);
    return $self;
}


=head2 next_marker()

 Title   : next_marker()
 Usage   : 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub next_marker {
    my $self = shift;
    if ($self->{'done_markers'}) { return; }
    my $line;
    if (!$self->{'processed_top'}) {
	until (($line = $self->_readline()) =~ /^\s+Markers\s+Distance/) { }
	$self->{'processed_top'} = 1;
    }
    $line = $self->_readline();
    chomp $line;
    if ($line =~ /-{5,}/) { # terminator is ------- 
	$self->{'done_markers'} = 1;
	return undef;
    }
    $line =~ s/^\s+//;
    my ($marker_number,$marker_name,$marker_distance) = split(/\s+/,$line);

    my $pos = new Bio::Map::OrderedPositionWithDistance
	(-positions => $marker_number,
	 -distance => $marker_distance
	 );
    my $o_marker = new Bio::Map::Marker('-name'=> $marker_name,
					'-position' => $pos);
    return $o_marker;
}

=head2 summary_info()

 Title   : summary_info()
 Usage   : $o_mapmaker->summary_info();
 Function:
 Example :
 Returns :
 Args    :


=cut

sub summary_info {
    my $self = shift;
    if (!$self->{'done_markers'}) {
	$self->warn("You can't get the summary info until all markers are read from the mapmaker file. Sorry.\n");
	return;
    }
    my $line = $self->_readline();
    chomp $line;
    return $line;
}

sub DESTROY {
    my $self = shift;

    $self->close();
}

1;
