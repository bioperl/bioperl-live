# $Id$
#
# BioPerl module for Bio::Tools::Run::EMBOSSApplication
#
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::EMBOSSApplication -  class for EMBOSS Applications

=head1 SYNOPSIS

#

=head1 DESCRIPTION

#

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
  http://bio.perl.org/bioperl-bugs/

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

package Bio::Tools::Run::EMBOSSApplication;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;
@ISA = qw(Bio::Root::RootI);


sub new {
  my($class, $args) = @_;
  my $self = $class->SUPER::new();
  #print join( " ", @args), "\n";
  $self->{ '_attributes' } = $args;

  $self->name($self->{ '_attributes' }->{'name'});
  delete $self->{ '_attributes' }->{'name'};

  $self->descr($self->{ '_attributes' }->{'documentation'});
  delete $self->{ '_attributes' }->{'documentation'};
#reorganize into group and subgroup
  $self->groups($self->{ '_attributes' }->{'groups'});
  delete $self->{ '_attributes' }->{'groups'};

#  $self-> test_input;

  return $self;
}

=head2 run

 Title   : run
 Usage   : $embossapplication->run
 Function: Runs the EMBOSS program.
 Returns : string only, will return objects!
 Args    : hash of input to the program

=cut

sub run {
    my ($self, $input) = @_;
    # test input
    use Data::Dumper;
    print Dumper($input);

    # match acd attributes against the input
    foreach my $attr (keys %{$self->{'_attributes'}}) {
	print $attr, "\n";
	 my $attr_name = substr($attr, 1) if substr($attr, 0, 1) =~ /\W/;
	 my $input_value = '';
	 $input_value = %{$input}->{$attr} if defined %{$input}->{$attr};
	 $self->throw('Attribute [$attr] not set') if 
	     defined %{$self}->{'_attributes'}->{$attr_name}->{optional} and $input_value ;

	 $self->throw('Attribute [$attr] not set')
	     if ( %{$self}->{'_attributes'}->{$attr}->{optional} eq 'N' and defined %{$input}->{$attr_name} and %{$input}->{$attr_name} eq '' ) or 

		( %{$self}->{'_attributes'}->{$attr}->{required} eq 'Y' and not exists %{$input}->{$attr_name})  ;
    }

    # collect the options into a string 
    my $option_string = '';
    foreach my $attr (keys %{$input}) {
	my $attr_name = substr($attr, 1) if substr($attr, 0, 1) =~ /\W/;

	# validate the values against acd

	print $attr_name, " ", %{$input}->{$attr}, "\n";
	$option_string .= $attr;
	$option_string .= " ". %{$input}->{$attr} 
	   if %{$input}->{$attr};
    }
#    print $option_string, "\n";
    my $runstring = join (' ', $self->name, $option_string, '-auto');
    print STDERR "Command line: ", $runstring, "\n" if $self->verbose > 0;
    return `$runstring`;
}

=head2 name

 Title   : name
 Usage   : $embossprogram->name
 Function: sets/gets the name of the EMBOSS program
           Setting is done by the EMBOSSFactory object,
           you should only get it.
 Throws  : 
 Returns : name string
 Args    : None

=cut

sub name {
    my ($self,$value) = @_;
    if (defined $value) {
	$self->{'name'} = $value;
    }
    return $self->{'name'};
}


=head2 descr

 Title   : descr
 Usage   : $embossprogram->descr
 Function: sets/gets the descr of the EMBOSS program
           Setting is done by the EMBOSSFactory object,
           you should only get it.
 Throws  : 
 Returns : description string
 Args    : None

=cut

sub descr {
    my ($self,$value) = @_;
    if (defined $value) {
	$self->{'descr'} = $value;
    }
    return $self->{'descr'};
}


=head2 groups

 Title   : groups
 Usage   : $embossprogram->groups
 Function: sets/gets the groups of the EMBOSS program
           Setting is done by the EMBOSSFactory object,
           you should only get it.

           There can be more than one group in which case 
           names are separated by ':'character.
 Throws  : 
 Returns : groups string
 Args    : None

=cut

sub groups {
    my ($self,$value) = @_;
    if (defined $value) {
	$self->{'groups'} = $value;
    }
    return $self->{'groups'};
}


1;
