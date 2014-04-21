#
# BioPerl module for Bio::PopGen::IO::phase
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

Bio::PopGen::IO::phase - A parser for Phase format data

=head1 SYNOPSIS

# Do not use directly, use through the Bio::PopGen::IO driver

  use Bio::PopGen::IO;
  my $io = Bio::PopGen::IO->new(-format => 'phase',
                               -file   => 'data.phase');

  # Some IO might support reading in a population at a time

  my @population;
  while( my $ind = $io->next_individual ) {
      push @population, $ind;
  }


=head1 DESCRIPTION

A driver module for Bio::PopGen::IO for parsing phase data.

PHASE is defined in http://www.stat.washington.edu/stephens/instruct2.1.pdf

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


package Bio::PopGen::IO::phase;
use vars qw($FieldDelim $AlleleDelim $NoHeader);
use strict;

($FieldDelim, $AlleleDelim, $NoHeader) = (' ', '\s+', 1);




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
           -field_delimiter => ' '
           -allele_delimiter=> '\s+'
           -no_header       => 0,


=cut


sub _initialize  {

    my($self, @args) = @_;

    $Bio::PopGen::Genotype::BlankAlleles='';

    my ($fieldsep,$all_sep, 
	$noheader) = $self->_rearrange([qw(FIELD_DELIMITER
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

=head2 next_individual

 Title   : next_individual
 Usage   : my $ind = $popgenio->next_individual;
 Function: Retrieve the next individual from a dataset
 Returns : L<Bio::PopGen::IndividualI> object
 Args    : none


=cut

sub next_individual  {
    my ($self) = @_;

    my ($sam,@marker_results,$number_of_ids,$number_of_markers,
	$marker_positions,$micro_snp);

    while( defined( $_ = $self->_readline) ) {
	chomp;
	next if( /^\s+$/ || ! length($_) );
	last;
    }
    
    return unless defined $_; 
    if( $self->flag('no_header') || defined $self->{'_header'} ) {

	####### sometimes there is some marker info at the start of a phase input file
	####### we collect it in the next few lines if there is. Should this info be held in a marker object?

	if(!$self->{'_count'} && /^\s*\d+$/){
	    $self->flag('number_of_ids',$_);
	    #print "number_of_ids : $number_of_ids\n";
	    $self->{'_count'}++;
	    return $self->next_individual;
	} elsif($self->{'_count'} == 1 && /^\s*\d+$/){
	    $self->flag('number_of_markers',$_);
	    #print "number_of_markers : $number_of_markers\n";
	    $self->{'_count'}++;
	    return $self->next_individual;
	} elsif($self->{'_count'} == 2 && /^\s*P\s\d/){
	    $self->flag('marker_positions',$_);
	    #print "marker_position : $marker_positions\n";
	    $self->{'_count'}++;
	    return $self->next_individual;
	} elsif($self->{'_count'} == 3 && /^\s*(M|S)+\s*$/i){
	    $self->flag('micro_snp',$_);
	    #print "microsat or snp : $micro_snp\n";
	    $self->{'_count'}++;
	    return $self->next_individual;
	} elsif(/^\s*\#/){
	    ($self->{'_sam'}) = /^\s*\#(.+)/;
	    #print "sample : $self->{'_sam'}\n";
	    $self->{'_count'}++;
	    return $self->next_individual;
	} else {
	    if( $self->{'_row1'} ) {
		# if we are looking at the 2nd row of alleles for this id

		@{$self->{'_second_row'}} = 
		    split($self->flag('field_delimiter'),$_);

		for my $i(0 .. $#{$self->{'_first_row'}}){

		    push(@{$self->{'_marker_results'}},
			 $self->{'_first_row'}->[$i].
			 $self->flag('field_delimiter').
			 $self->{'_second_row'}->[$i]);
		}
		$self->{'_row1'} = 0;
	    } else {
		# if we are looking at the first row of alleles for this id
		@{$self->{'_marker_results'}} = ();
		@{$self->{'_first_row'}} = split($self->flag('field_delimiter'),$_);
		$self->{'_row1'} = 1;
		return $self->next_individual;
	    }
	}

	my $i = 1;
	foreach my $m ( @{$self->{'_marker_results'}} ) {
	    $m =~ s/^\s+//;
	    $m =~ s/\s+$//;
	    my $markername;
	    if( defined($self->flag('marker_positions')) ) {
		$markername = (split($self->flag('field_delimiter'), $self->flag('marker_positions')))[$i];
	    } elsif( defined $self->{'_header'} ) {
		$markername = $self->{'_header'}->[$i] || "$i";
	    } else { 
		$markername = "$i";
	    }

	    my $markertype;
	    if( defined($self->flag('marker_positions')) ) {
		$markertype = (split('', $self->flag('micro_snp')))[$i-1];
	    } else {
		$markertype = "S";
	    }

	    $self->debug( "markername is $markername alleles are $m\n");
	    my @alleles = split($self->flag('allele_delimiter'), $m);	

	    $m = Bio::PopGen::Genotype->new(-alleles       =>\@alleles,
					    -marker_name   => $markername,
					    -marker_type   => $markertype,
					    -individual_id => $self->{'_sam'}); 
	    $i++; 
	}
	return Bio::PopGen::Individual->new(-unique_id => $self->{'_sam'},
					   -genotypes =>\@{$self->{'_marker_results'}},
					   );

    } else {
	$self->{'_header'} = [split($self->flag('field_delimiter'),$_)];
	return $self->next_individual; # rerun loop again
    }
    return;
}

=head2 next_population

 Title   : next_population
 Usage   : my $ind = $popgenio->next_population;
 Function: Retrieve the next population from a dataset
 Returns : L<Bio::PopGen::PopulationI> object
 Args    : none
 Note    : Many implementation will not implement this

=cut

sub next_population{
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
 Returns : none
 Args    : L<Bio::PopGen::PopulationI> object(s)

=cut


sub write_individual {
    my ($self,@inds) = @_;
    my $fielddelim  = $self->flag('field_delimiter');
    my $alleledelim = $self->flag('allele_delimiter');

    # For now capture print_header flag from @inds
    my $header = 1;
    $header = pop(@inds) if($inds[-1] =~ m/^[01]$/);

    foreach my $ind ( @inds ) {
	if (! ref($ind) || ! $ind->isa('Bio::PopGen::IndividualI') ) {
	    $self->warn("Cannot write an object that is not a Bio::PopGen::IndividualI object ($ind)");
	    next;
	}

	# sort lexically until we have a better way to insure a consistent order
	my @marker_names = sort $ind->get_marker_names;

	if ($header) {
	    my $n_markers = scalar(@marker_names);
	    $self->_print( "1\n");
	    $self->_print( $n_markers, "\n");
	    if( $self->flag('no_header') && 
		! $self->flag('header_written') ) {
		$self->_print(join($fielddelim, ('P', @marker_names)), "\n");
		$self->flag('header_written',1);
	    }
	    foreach my $geno ($ind->get_Genotypes()) {
		$self->_print($geno->marker_type);
	    }
	    $self->_print("\n");
	}
	
	my(@row1,@row2);
	for (@marker_names){
	    my $geno = $ind->get_Genotypes(-marker => $_);
	    my @alleles = $geno->get_Alleles(1);
	    push(@row1,$alleles[0]);
	    push(@row2,$alleles[1]);
	}
	$self->_print("#",$ind->unique_id,"\n",
		      join($fielddelim,@row1),"\n",
		      join($fielddelim,@row2),"\n");
    }
}

=head2 write_population

 Title   : write_population
 Usage   : $popgenio->write_population($pop);
 Function: Write a population out in the file format
 Returns : none
 Args    : L<Bio::PopGen::PopulationI> object(s)
 Note    : Many implementation will not implement this

=cut


sub write_population {
    my ($self,@pops) = @_;
    my $fielddelim  = $self->flag('field_delimiter');
    my $alleledelim = $self->flag('allele_delimiter');

    foreach my $pop ( @pops ) {
	if (! ref($pop) || ! $pop->isa('Bio::PopGen::PopulationI') ) {
	    $self->warn("Cannot write an object that is not a Bio::PopGen::PopulationI object");
	    next;
	}
	# sort lexically until we have a better way to insure a consistent order
	my @marker_names = sort $pop->get_marker_names;
	my $n_markers = scalar(@marker_names);
	$self->_print( $pop->get_number_individuals, "\n");
	$self->_print( $n_markers, "\n");
	if( $self->flag('no_header') && 
	    ! $self->flag('header_written') ) {
	    $self->_print( join($fielddelim, ('P', @marker_names)), "\n");
	    $self->flag('header_written',1);
	}

	foreach (@marker_names) {
	    $self->_print(($pop->get_Genotypes($_))[0]->marker_type);
	}
	$self->_print("\n");
	
	$self->write_individual( $pop->get_Individuals, 0 );
    }
}

1;
