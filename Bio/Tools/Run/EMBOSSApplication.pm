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

  use Bio::Factory::EMBOSS;
  # get an EMBOSS application object from the EMBOSS factory
  $factory = new Bio::Factory::EMBOSS
  $application = $factory->program('embossversion');
  # run the application with an optional hash containing parameters
  $result = $application->run(); # returns a string or creates a file
  print $result . "\n";

  $water = $factory->program('water');

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

The EMBOSSApplication class can represent EMBOSS any program. It is
created by a Bio::Factory::EMBOSS object which primes it by reading
in the ADC description of the command line options. See also
L<Bio::Factory::EMBOSS>.

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

=head2 CONTRIBUTORS

Email: jason@bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::Run::EMBOSSApplication;
use vars qw(@ISA);
use strict;
use Data::Dumper;
use Bio::Root::Root;
use Bio::Root::IO;

@ISA = qw(Bio::Root::Root);


sub new {
  my($class, $args) = @_;
  my $self = $class->SUPER::new();  
  $self->{ '_attributes' } = $args;

  $self->name($self->{ '_attributes' }->{'name'});
  delete $self->{ '_attributes' }->{'name'};

  $self->descr($self->{ '_attributes' }->{'documentation'});
  delete $self->{ '_attributes' }->{'documentation'};

  $self->group($self->{ '_attributes' }->{'groups'});
  delete $self->{ '_attributes' }->{'groups'};

  $self->{'_io'} = new Bio::Root::IO('-verbose' => $self->verbose);
  return $self;
}

=head2 run

 Title   : run
 Usage   : $embossapplication->run($attribute_hash)
 Function: Runs the EMBOSS program.
 Returns : string or creates files for now; will return objects!
 Args    : hash of input to the program

=cut

sub run {
    my ($self, $input) = @_;
    $self->{'_io'}->_io_cleanup();
    # test input
    print Dumper($input) if $self->verbose > 0;
    $self->_acd2input($input);

    # collect the options into a string 
    my $option_string = '';
    foreach my $attr (keys %{$input}) {
	my $attr_name = substr($attr, 1) if substr($attr, 0, 1) =~ /\W/;
	
	my $array = 0;
	
	if( defined $input->{$attr} && ref($input->{$attr}) ) {

	    my (@pieces);
	    
	    if( $array = (ref($input->{$attr}) =~ /array/i) ) {
		foreach my $s ( @{$input->{$attr}} ) {
		    @pieces = @{$input->{$attr}};
		}
	    } else {
		@pieces = ($input->{$attr});
	    }
	    if( $pieces[0]->isa('Bio::PrimarySeqI') ) {
		require Bio::SeqIO;
		my ($tfh,$tempfile) = $self->{'_io'}->tempfile();
		my $out = new Bio::SeqIO(-format => 'fasta',
					 -fh     => $tfh);
		foreach my $seq ( @pieces ) {
			$out->write_seq($seq);
		    }
		$out->close();
		$input->{$attr} = $tempfile;
	    } elsif( $pieces[0]->isa('Bio::Align::AlignI') ) {
		require Bio::AlignIO;
		my ($tfh,$tempfile) = $self->{'_io'}->tempfile();
		my $out = new Bio::AlignIO(-format => 'msf',
					   -fh     => $tfh);
		foreach my $p ( @pieces ) {
		    $out->write_aln($p);
		}
		$input->{$attr} = $tempfile;
	    }
	}
	# ADD: validate the values against acd
	
	$self->debug("Input attr: ". $attr_name. " => ". 
		     %{$input}->{$attr}. "\n"); 
	$option_string .= " " . $attr;
	$option_string .= " ". %{$input}->{$attr} 
	   if %{$input}->{$attr};
    }

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
	$self->{'_name'} = $value;
    }
    return $self->{'_name'};
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
	$self->{'_descr'} = $value;
    }
    return $self->{'_descr'};
}


=head2 group

 Title   : group
 Usage   : $embossprogram->group
 Function: sets/gets the group of the EMBOSS program
           Setting is done by the EMBOSSFactory object,
           you should only get it.

           If the application is assigned into a subgroup
           use l<subgroup> to get it.

 Throws  : 
 Returns : string, group name
 Args    : group string

=cut

sub group {
    my ($self,$value) = @_;
    if (defined $value) {
	my ($group, $subgroup) = split ':', $value;
	$self->{'_group'} = $group;
	$self->{'_subgroup'} = $subgroup;
    }
    return $self->{'_group'};
}


=head2 subgroup

 Title   : subgroup
 Usage   : $embossprogram->subgroup
 Function: sets/gets the subgroup of the EMBOSS program
           Setting is done by the EMBOSSFactory object,
           you should only get it.
 Throws  : 
 Returns : string, subgroup name; undef if not defined
 Args    : None

=cut

sub subgroup {
    my ($self) = @_;
    return $self->{'_subgroup'};
}




=head2 Internal methods

Do not call these methods directly

=head2 _acd2input

 Title   : _acd2input
 Usage   : $embossfactory->_acd2input()
 Function: compares ACD file requirements to input hash
 Returns : 
 Throws  : if requirements are not met 
 Args    : None

=cut

sub _acd2input($input) {
    my ($self, $input) = @_;
    # match acd attributes against the input
    foreach my $attr (keys %{$self->{'_attributes'}}) {
	$self->debug( "ACD Attr: |". $attr. "|\n");
	my $input_value = '';
	my $input_key = '';
	$input_key = 1, $input_value = %{$input}->{"-$attr"} if defined %{$input}->{"-$attr"};
#	my $input_name = substr($input_value, 1) if substr($input_value, 0, 1) =~ /\W/;
	$self->debug("--------------$input_value, $attr\n") if $input_key;

#	 $self->throw('Attribute [$attr] not set') if 
#	     defined %{$self}->{'_attributes'}->{$attr_name}->{optional} and $input_value ;
#	 
#	 $self->throw('Attribute [$attr] not set')
#	     if ( %{$self}->{'_attributes'}->{$attr}->{optional} eq 'N' and defined %{$input}->{$attr_name} and %{$input}->{$attr_name} eq '' ) or 
#		 
##
#		 ( %{$self}->{'_attributes'}->{$attr}->{required} eq 'Y' and not exists %{$input}->{$attr_name})  ;
    }
}


1;
