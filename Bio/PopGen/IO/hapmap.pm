# $Id$
#
# BioPerl module for Bio::PopGen::IO::hapmap
#
# Cared for by Rich  <dobbo@thevillas.eclipse.co.uk>
#
# Copyright Rich 
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PopGen::IO::hapmap - A parser for HapMap output data

=head1 SYNOPSIS

# Do not use directly, use through the Bio::PopGen::IO driver

  use Bio::PopGen::IO;
  my $io = new Bio::PopGen::IO(-format => 'hapmap',
                               -file   => 'data.hapmap');

  # Some IO might support reading in a population at a time

  my @population;
  while( my $ind = $io->next_individual ) {
      push @population, $ind;
  }

=head1 DESCRIPTION

Describe the object here

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
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Rich 

Email dobbo@thevillas.eclipse.co.uk

=head1 CONTRIBUTORS

Jason Stajich, jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::PopGen::IO::hapmap;
use vars qw(@ISA $FieldDelim $AlleleDelim $NoHeader $InsertMissing);
use strict;

($FieldDelim,$AlleleDelim,$NoHeader,$InsertMissing) =( '\s+','',0,0);

use Bio::PopGen::IO;
use Bio::PopGen::Individual;
use Bio::PopGen::Population;
use Bio::PopGen::Genotype;

@ISA = qw(Bio::PopGen::IO );


=head2 new

 Title   : new
 Usage   : my $obj = new Bio::PopGen::IO::hapmap();
 Function: Builds a new Bio::PopGen::IO::hapmap object 
 Returns : an instance of Bio::PopGen::IO::hapmap
 Args    : [optional, these are the current defaults] 
           -field_delimiter => ','
           -allele_delimiter=> '\s+'
           -no_header       => 0,


=cut


sub _initialize  {

    my($self, @args) = @_;

    $Bio::PopGen::Genotype::BlankAlleles='';

    my ($fieldsep,$all_sep, 
	$noheader,$insertmissing) = $self->_rearrange([qw(FIELD_DELIMITER
							  ALLELE_DELIMITER
							  NO_HEADER)],@args);


    $self->flag('no_header', defined $noheader ? $noheader : $NoHeader);
    $self->flag('field_delimiter',defined $fieldsep ? $fieldsep : $FieldDelim);
    $self->flag('allele_delimiter',defined $all_sep ? $all_sep : $AlleleDelim);


    $self->{'_header'} = undef;
    return 1;

}

=head2 flag

 Title   : flag
 Usage   : $obj->flag($flagname,$newval)
 Function: Get/Set the flag value
 Returns : value of a flag (a boolean)
 Args    : A flag name, currently we expect 
           'no_header', 'field_delimiter', or 'allele_delimiter' 
           on set, new value (a boolean or undef, optional)


=cut

sub flag  {

    my $self = shift;
    my $fieldname = shift;
    return unless defined $fieldname;
    return $self->{'_flag'}->{$fieldname} = shift if @_;
    return $self->{'_flag'}->{$fieldname};

}

sub _pivot {
    my ($self) = @_;

    my (@cols,@rows,@idheader);
    while ($_ = $self->_readline){
	chomp($_);
	next if( /^\s*\#/ || /^\s+$/ || ! length($_) );
	if( /^rs\#\s+alleles\s+chrom\s+pos\s+strand/ ) {
	    @idheader = split $self->flag('field_delimiter');
	} else { 
	    push @cols, [split $self->flag('field_delimiter')];
	}
    }
    my $startingcol = 10;

    $self->{'_header'} = [ map { $_->[0] } @cols];
    for my $n ($startingcol.. $#{ $cols[ 0 ]}) { 
	my $column = [ $idheader[$n],
		       map{ $_->[ $n ] } @cols ];	
	push (@rows, $column); 
    }
    $self->{'_pivot'} = [@rows];
    $self->{'_i'} = 0;
}


=head2 next_individual

 Title   : next_individual
 Usage   : my $ind = $popgenio->next_individual;
 Function: Retrieve the next individual from a dataset
 Returns : L<Bio::PopGen::IndividualI> object
 Args    : none


=cut

sub next_individual  {
    my ($self) = @_;
    unless($self->{'_pivot'}){
	#if it's the first time then pivot the table and store.
	#Lines will now be read from the stored pivot version of the input file
	$self->_pivot;
    }

    $_ = $self->{'_pivot'}->[$self->{'_i'}++];

    return unless defined $_;

    # Store all the marker related info. Now that the pivot has taken
    # place this is in the first few lines of the file Maybe this
    # should be put in a marker object. Doesn't seem to fit too well
    # though

    my ($samp,@marker_results) = @$_;	

    # at some point use all this info
    my $i = 1;
    foreach my $m ( @marker_results ) {
	$m =~ s/^\s+//;
	$m =~ s/\s+$//;
	my $markername;
	if( defined $self->{'_header'} ) {
	    $markername = $self->{'_header'}->[$i-1];
	} else { 
	    $markername = "Marker$i";
	}

	my @alleles = split($self->flag('allele_delimiter'), $m);
	if( @alleles != 2 ) { 
	    warn("$m for $samp\n");
	} else { 
	    $m = new Bio::PopGen::Genotype(-alleles      => \@alleles,
					   -marker_name  => $markername,
					   -individual_id=> $samp);
	}
	$i++; 
    }

    return new Bio::PopGen::Individual(-unique_id => $samp,
				       -genotypes => \@marker_results);

}

=head2 next_population

 Title   : next_population
 Usage   : my $ind = $popgenio->next_population;
 Function: Retrieve the next population from a dataset
 Returns : L<Bio::PopGen::PopulationI> object
 Args    : none
 Note    : Many implementation will not implement this

=cut

sub next_population {
    my ($self) = @_;
    my @inds;
    while( my $ind = $self->next_individual ) {
	push @inds, $ind;
    }
    Bio::PopGen::Population->new(-individuals => \@inds);
}

=head2 write_individual

 Title   : write_individual
 Usage   : $popgenio->write_individual($ind);
 Function: Write an individual out in the file format
           NOT SUPPORTED  BY hapmap format
 Returns : none
 Args    : L<Bio::PopGen::PopulationI> object(s)

=cut

sub write_individual {
    my ($self,@inds) = @_;

    # data from hapmap is output, not input, so 
    # we don't need a method for writing and input file

    $self->throw('writing to hapmap format is not implemented');
}

=head2 write_population

 Title   : write_population
 Usage   : $popgenio->write_population($pop);
 Function: Write a population out in the file format
           NOT SUPPORTED  BY hapmap format
 Returns : none
 Args    : L<Bio::PopGen::PopulationI> object(s)
 Note    : Many implementation will not implement this

=cut

sub write_population {
    my ($self,@inds) = @_;
    $self->throw('writing to hapmap format is not implemented');
}

1;
