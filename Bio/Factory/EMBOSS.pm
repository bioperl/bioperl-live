# $Id$
#
# BioPerl module for Bio::Factory::EMBOSS
#
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Factory::EMBOSS - EMBOSS appliaction factory class

=head1 SYNOPSIS

  # get an EMBOSS factory
  use Bio::Factory::EMBOSS;
  $f = Bio::Factory::EMBOSS -> new();
  # get an EMBOSS application  object from the factory
  $water = $f->program('water');

  # here is an example of running the application
  # water can compare 1 seq against 1->many sequences
  # in a database using Smith-Waterman
  my $seq_to_test; # this would have a seq here
  my @seqs_to_check; # this would be a list of seqs to compare 
                     # (could be just 1)
  my $wateroutfile = 'out.water';
  $water->run({ '-sequencea' => $seq_to_test,
              '-seqall'    => \@seqs_to_check,
              '-gapopen'   => '10.0',
              '-gapextend' => '0.5',
              '-outfile'   => $wateroutfile});
  # now you might want to get the alignment
  use Bio::AlignIO;
  my $alnin = new Bio::AlignIO(-format => 'emboss',
			       -file   => $wateroutfile);

  while( my $aln = $alnin->next_aln ) {
      # process the alignment -- these will be Bio::SimpleAlign objects
  }

=head1 DESCRIPTION

The EMBOSS factory class encapsulates access to EMBOSS programs.  A
factory object allows creation of only known applications.

If you want to check command line options before sending them to the
program set $factory-E<gt>verbose to positive integer. The value is
passed on to programs objects and the ADC description of the available
command line options is parsed and compared to input.

See also L<Bio::Tools::Run::EMBOSSApplication> and
L<Bio::Tools::Run::EMBOSSacd>.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                         - General discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.  Bug reports can be submitted via
 email or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHOR - Heikki Lehvaslaiho

Email:  heikki@ebi.ac.uk
Address:

     EMBL Outstation, European Bioinformatics Institute
     Wellcome Trust Genome Campus, Hinxton
     Cambs. CB10 1SD, United Kingdom

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Factory::EMBOSS;
use vars qw(@ISA $EMBOSSVERSION);
use strict;

use Bio::Root::Root;
use Bio::Tools::Run::EMBOSSApplication;
use Bio::Factory::ApplicationFactoryI;
@ISA = qw(Bio::Root::Root Bio::Factory::ApplicationFactoryI );

$EMBOSSVERSION = "2.0.0";

sub new {
  my($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  # set up defaults

  my($location) =
      $self->_rearrange([qw(LOCATION )],
			@args);

  $self->{ '_programs' } = {};
  $self->{ '_programgroup' } = {};
  $self->{ '_groups' } = {};

  $self->location($location) if $location;

  $self->_program_list; # retrieve info about available programs

  return $self;

}

=head2 location

 Title   : location
 Usage   : $embossfactory->location
 Function: get/set the location of EMBOSS programs.
           Valid values are 'local' and 'novella'.
 Returns : string, defaults to 'local'
 Args    : string 

=cut

sub location {
    my ($self, $value) = @_;
    my %location = ('local' => '1',
		    'novella' => '1'
		    );
    if (defined $value) {
	$value = lc $value;
	if ($location{$value}) {
	    $self->{'_location'} = $value;
	} else {
	    $self->warn("Value [$value] not a valid value for ".
			"location(). Defaulting to [local]");
	    $self->{'_location'} = 'local';
	}
    }
    $self->{'_location'} ||= 'local';
    return $self->{'_location'};
}


=head2 program

 Title   : program
 Usage   : $embossfactory->program('program_name')
 Function: Creates a representation of a single EMBOSS program
 Returns : Bio::Tools::Run::EMBOSSApplication object
 Args    : string, program name

=cut

sub program {
    my ($self, $value) = @_;

    unless( $self->{'_programs'}->{$value} ) {
	$self->warn("Application [$value] is not available!");
	return undef;
    }
    my $attr = {};
    $attr->{name} = $value;
    $attr->{verbose} = $self->verbose;

    my $appl = Bio::Tools::Run::EMBOSSApplication ->new($attr);
    return $appl;
}

=head2 version

 Title   : $self->version
 Usage   : $embossfactory->version()
 Function: gets the version of EMBOSS programs
 Throws  : if EMBOSS suite is not accessible
 Returns : version value
 Args    : None

=cut

sub version {
    my ($self) = @_;
    my ($version);
    eval {
	$version = `embossversion -auto`;
    };
    $self->throw("EMBOSS suite of programs is not available \n\n$@")
	if $@;
    chop $version;

    # compare versions
    my ($thisv, $embossv);
    $version =~ /(\d+)\.(\d+)\.(\d+)/;
    $thisv = "$1.$2$3";
    $EMBOSSVERSION =~ /(\d+)\.(\d+)\.(\d+)/;
    $embossv = "$1.$2$3";
    $self->throw("EMBOSS has to be at least version $EMBOSSVERSION\n")
	if $thisv < $embossv;

    return $version;
}


=head2 Programs

These methods allow the programmer to query the EMBOSS suite and find
out which program names can be used and what arguments can be used.

=head2 program_info

 Title   : program_info
 Usage   : $embossfactory->program_info('emma')
 Function: Finds out if the program is available.
 Returns : definition string of the program, undef if program name not known
 Args    : string, prgramname

=cut

sub program_info {
    my ($self, $value) = @_;
    return $self->{'_programs'}->{$value};
}


=head2 Internal methods

Do not call these methods directly

=head2 _program_list

 Title   : _program_list
 Usage   : $embossfactory->_program_list()
 Function: Finds out what programs are available.
           Writes the names into an internal hash.
 Returns : true if successful
 Args    : None

=cut

sub _program_list {
    my ($self) = @_;
    if( $^O =~ /MSWIN/i ||
	$^O =~ /Mac/i ) { return; }
    {
	local * SAVERR;
	open SAVERR, ">&STDERR";
	open STDERR, ">/dev/null";
	open(WOSSOUT, "wossname -auto |") || return;
	open STDERR, ">&SAVERR";

    }
    local $/ = "\n\n";
    while(<WOSSOUT> ) {	
	my ($groupname) = (/^([A-Z][A-Z0-9 ]+)$/m);
	#print $groupname, "\n" if $groupname;
	$self->{'_groups'}->{$groupname} = [] if $groupname;
	while ( /^([a-z]\w+) +(.+)$/mg ) {
	    #print "$1\t$2 \n" if $1;
	    $self->{'_programs'}->{$1} = $2 if $1; 
	    $self->{'_programgroup'}->{$1} = $groupname if $1;
	    push @{$self->{'_groups'}->{$groupname}}, $1 if $1;
	}
    }
    close(WOSSOUT);

}

1;
