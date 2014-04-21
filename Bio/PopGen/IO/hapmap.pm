#
# BioPerl module for Bio::PopGen::IO::hapmap
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Rich Dobson <r.j.dobson-at-qmul.ac.uk>
#
# Copyright Rich Dobson
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::PopGen::IO::hapmap - A parser for HapMap output data

=head1 SYNOPSIS

  # Do not use directly, use through the Bio::PopGen::IO driver

  use Bio::PopGen::IO;
  my $io = Bio::PopGen::IO->new(-format => 'hapmap',
                               -file   => 'data.hapmap');

  # Some IO might support reading in a population at a time

  my @population;
  while( my $ind = $io->next_individual ) {
      push @population, $ind;
  }

=head1 DESCRIPTION

A driver module for Bio::PopGen::IO for parsing hapmap data.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Rich Dobson

Email r.j.dobson-at-qmul.ac.uk

=head1 CONTRIBUTORS

Jason Stajich, jason-at-bioperl.org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::PopGen::IO::hapmap;
use vars qw($FieldDelim $AlleleDelim $NoHeader $StartingCol);
use strict;

($FieldDelim,$AlleleDelim,$NoHeader,$StartingCol) =( '\s+','',0,11);

use Bio::PopGen::Individual;
use Bio::PopGen::Population;
use Bio::PopGen::Genotype;

use base qw(Bio::PopGen::IO);


=head2 new

 Title   : new
 Usage   : my $obj = Bio::PopGen::IO::hapmap->new();
 Function: Builds a new Bio::PopGen::IO::hapmap object 
 Returns : an instance of Bio::PopGen::IO::hapmap
 Args    : [optional, these are the current defaults] 
           -field_delimiter => ','
           -allele_delimiter=> '\s+'
           -no_header       => 0,
           -starting_column => 11

=cut


sub _initialize  {

    my($self, @args) = @_;

    $Bio::PopGen::Genotype::BlankAlleles='';

    my ($fieldsep,$all_sep, 
	$noheader, $start_col) = $self->_rearrange([qw(FIELD_DELIMITER
						       ALLELE_DELIMITER
						       NO_HEADER
						       STARTING_COLUMN)],
						   @args);

    $self->flag('no_header', defined $noheader ? $noheader : $NoHeader);
    $self->flag('field_delimiter',defined $fieldsep ? $fieldsep : $FieldDelim);
    $self->flag('allele_delimiter',defined $all_sep ? $all_sep : $AlleleDelim);
    $self->starting_column(defined $start_col ? $start_col : $StartingCol );

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
    my $startingcol = $self->starting_column;

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
 Returns : A Bio::PopGen::IndividualI object
 Args    : none

See L<Bio::PopGen::IndividualI>

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
	    $self->warn("$m for $samp\n");
	} else { 
	    $m = Bio::PopGen::Genotype->new(-alleles       => \@alleles,
					    -marker_name   => $markername,
					    -marker_type   => 'S',          # Guess hapmap only has SNP data
					    -individual_id => $samp);
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
 Returns : Bio::PopGen::PopulationI object
 Args    : none
 Note    : Many implementation will not implement this

See L<Bio::PopGen::PopulationI>

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
 Args    : Bio::PopGen::PopulationI object(s)

See L<Bio::PopGen::PopulationI>

=cut

sub write_individual {
    my ($self,@inds) = @_;

    # data from hapmap is output, not input, so 
    # we don't need a method for writing and input file

    $self->throw_not_implemented();
}

=head2 write_population

 Title   : write_population
 Usage   : $popgenio->write_population($pop);
 Function: Write a population out in the file format
           NOT SUPPORTED  BY hapmap format
 Returns : none
 Args    : Bio::PopGen::PopulationI object(s)
 Note    : Many implementation will not implement this

See L<Bio::PopGen::PopulationI>

=cut

sub write_population {
    my ($self,@inds) = @_;
    $self->throw_not_implemented();
}


=head2 starting_column

 Title   : starting_column
 Usage   : $obj->starting_column($newval)
 Function: Column where data starts
 Example : 
 Returns : value of starting_column (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub starting_column{
    my $self = shift;

    return $self->{'starting_column'} = shift if @_;
    return $self->{'starting_column'};
}

1;
