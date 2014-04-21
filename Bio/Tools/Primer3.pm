#
# BioPerl module for Bio::Tools::Primer3
#
# Copyright (c) 2003 bioperl, Rob Edwards. All Rights Reserved.
#           This module is free software; you can redistribute it and/or
#           modify it under the same terms as Perl itself.
#
# Copyright Rob Edwards
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Primer3 - Create input for and work with the output from
the program primer3

=head1 SYNOPSIS

 # parse primer3 output to get some data
 # this is also called from Bio::Tools::Run::Primer3
 use Bio::Tools::Primer3;

 # read a primer3 output file
 my $p3 = Bio::Tools::Primer3->new(-file=>"data/primer3_output.txt");

 # how many results were there?
 my $num = $p3->number_of_results;
 print "There were $num results\n";

 # get all the results
 my $all_results = $p3->all_results;
 print "ALL the results\n";
 foreach my $key (keys %{$all_results}) {
    print "$key\t${$all_results}{$key}\n";
 }

 # get specific results
 my $result1 = $p3->primer_results(1);
 print "The first primer is\n";
 foreach my $key (keys %{$result1}) {
    print "$key\t${$result1}{$key}\n";
 }

 # get the results as a Bio::Seq::PrimedSeq stream
 my $primer = $p3->next_primer;
 print "The left primer in the stream is ",
   $primer->get_primer('-left_primer')->seq->seq, "\n";

=head1 DESCRIPTION

Bio::Tools::Primer3 creates the input files needed to design primers using
primer3 and provides mechanisms to access data in the primer3 output files.

This module provides a bioperl interface to the program primer3. See
http://www-genome.wi.mit.edu/genome_software/other/primer3.html
for details and to download the software.

This module is based on one written by Chad Matsalla
(bioinformatics1@dieselwurks.com)

I have ripped some of his code, and added a lot of my own. I hope he
is not mad at me!

This is probably best run in one of the two following ways:

  i. To parse the output from Bio::Tools::Run::Primer3.
     You will most likely just use next_primer to get the results from
     Bio::Tools::Run::Primer3.
  ii. To parse the output of primer3 handed to it as a file name.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

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
the bugs and their resolution.  Bug reports can be submitted via the web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR -

  Rob Edwards

  redwards@utmem.edu

  Based heavily on work of

  Chad Matsalla

  bioinformatics1@dieselwurks.com

=head1 CONTRIBUTORS

  Brian Osborne bosborne at alum.mit.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::Primer3;

use strict;
use Bio::Seq;
use Bio::Seq::PrimedSeq;
use Bio::SeqFeature::Primer;

use vars qw($AUTOLOAD @PRIMER3_PARAMS %OK_FIELD $ID);

BEGIN {
    @PRIMER3_PARAMS = qw(results seqobject);
    foreach my $attr (@PRIMER3_PARAMS) {$OK_FIELD{$attr}++}
}


use base qw(Bio::Root::Root Bio::Root::IO);


sub AUTOLOAD {
    my $self = shift;
    my $attr = $AUTOLOAD;
    $attr =~ s/.*:://;
    $self->throw("Unallowed parameter: $attr !") unless $OK_FIELD{$attr};
    $self->{$attr} = shift if @_;
    return $self->{$attr};
}


=head2 new

  Title   : new()
  Usage   : my $primer3 = Bio::Tools::Primer3->new(-file=>$file);
  Function: Parse primer3 output
  Returns : Does not return anything. If called with a filename will
            allow you to retrieve the results
  Args    : -file (optional) file of primer3 results to parse -verbose
            (optional) set verbose output
  Notes   :

=cut

sub new {
    my($class,%args) = @_;
    my $self = $class->SUPER::new(%args);
    if ($args{'-file'}) {
        $self->_readfile($args{'-file'});
    }
    if ($args{'-verbose'}) {
        $self->{'verbose'} = 1;
    }
    return $self;
}


=head2 number_of_results

  Title   : number_of_results()
  Usage   : my $count = $primer3->number_of_results();
  Function: Retrieve the number of primers returned from Primer3.
  Returns : A scalar
  Args    : None
  Notes   : This returns the count of the primers returned by Primer3
             (aka how many of them there are).
             This is one more than the maximum offset into the zero
             based list of primers that is accessed by primer_results().

=cut

sub number_of_results {
    my $self = shift;
    return $self->{'maximum_primers_returned'} + 1;
}


=head2 all_results

  Title   : all_results()
  Usage   : my $results = $primer3->all_results();
               or
            my $results = $primer3->all_results('primer3 result name', 'other results');
  Function: Retrieve the results returned from Primer3.
  Returns : A reference to a hash
  Args    : Optional array of specific results to retrieve

=cut

sub all_results {
    my ($self, @results) = @_;
    my %hash;
    if (@results) {
        # we only want a few things
        foreach my $result (@results) {
            $hash{$result} = $self->{'results'}->$result;
        }
    } else {
        foreach my $result (keys %{$self->{'results'}}) {
            $hash{$result}=$self->{'results'}->{$result};
        }
    }

    return \%hash;
}


=head2 primer_results

  Title   : primer_results()
  Usage   : my $results = $primer3->primer_results(2); # results for third primer
  Function: Retrieve the results returned from Primer3 for specific primer pairs.
  Returns : A reference to a hash
  Args    : A number between 0 and the maximum number of primers to retrieve

=cut

sub primer_results {
    my ($self, $toget) = @_;
    if ($toget > $self->{'maximum_primers_returned'}) {
        $self->warn("Didn't get any results for $toget");
        return 0;
    } else {
        return \%{$self->{'results_by_number'}->{$toget}};
    }
}


=head2 _readfile

  Title   : _readfile()
  Usage   : $self->_readfile();
  Function: An internal function that reads a file and sets up the results
  Returns : Nothing.
  Args    : None
  Notes   :

=cut

sub _readfile {
    my ($self, $file) = @_;
    $self->_initialize_io(-file=>$file);
    my $line;
    my $id='primer 3 parsed results'; # hopefully we'll get this, but we can set a temp id in case not.
    while (defined($line = $self->_readline()) ) {
        chomp $line;
        next unless ($line);
        my ($return, $value) = split /=/, $line;
        if (uc($return) eq "SEQUENCE") {
            $self->{seqobject} = Bio::Seq->new(-seq=>$value, $id=>$id);
            next;
        }
        if (uc($return) eq "PRIMER_SEQUENCE_ID") {
            if ($self->{seqobject}) {$self->{seqobject}->id($value)} else {$id=$value}
        }

        $self->{'results'}->{$return} = $value;
    }

    # convert the results to individual results
    $self->_separate();
}


=head2 next_primer

  Title   : next_primer()
  Usage   : while (my $primed_seq  = $primer3->next_primer()) {
  Function: Retrieve the primed sequence and a primer pair, one at a time
  Returns : Returns a Bio::Seq::PrimedSeq object, one at a time
  Args    : None
  Notes   : Use $primed_seq->annotated_seq to get an annotated sequence
            object you can write out.

=cut

sub next_primer {
    my $self = shift;
    # here we are going to convert the primers to Bio::SeqFeature::Primer objects
    # and the primer/sequence to Bio::Seq::PrimedSeq objects
    # the problem at the moment is that PrimedSeq can only take one sequence/primer pair, and
    # yet for each sequence we can have lots of primer pairs. We need a way to overcome this.
    # at the moment we can do this as a stream, I guess.

    $self->warn("No primers were found for: ".$self->{'seqobject'}->{'primary_id'})
      if (! $self->number_of_results);

    $self->{'next_to_return'} = 0 unless ($self->{'next_to_return'});
    return if ($self->{'next_to_return'} >= $self->number_of_results);
    my $results = $self->primer_results($self->{'next_to_return'});

    $self->throw("No left primer sequence") unless (${$results}{'PRIMER_LEFT_SEQUENCE'});
    $self->throw("No right primer sequence") unless (${$results}{'PRIMER_RIGHT_SEQUENCE'});
    $self->throw("No target sequence") unless ($self->{'seqobject'});

    my $left_seq = Bio::SeqFeature::Primer->new(
         -id         => 'left_primer',
         -seq        => ${$results}{'PRIMER_LEFT_SEQUENCE'},
         -display_id => ($self->{'next_to_return'} + 1),
    );

    my $right_seq = Bio::SeqFeature::Primer->new(
             -id         => "right_primer",
             -seq        => ${$results}{'PRIMER_RIGHT_SEQUENCE'},
             -display_id => ($self->{'next_to_return'} + 1) );

    # add data to the Primer objects
    for my $key (%$results) {
        # skip the primer sequence data, already added above
        next if ($key =~ /PRIMER_(LEFT|RIGHT)_SEQUENCE/i );
        if ($key =~ /PRIMER_LEFT/i) {
            $left_seq->add_tag_value($key, $$results{$key});
        } elsif ($key =~ /PRIMER_RIGHT/i) {
            $right_seq->add_tag_value($key, $$results{$key});
        }
    }

    my $primed_seq = Bio::Seq::PrimedSeq->new(
         -target_sequence => $self->{'seqobject'}->clone,
         -left_primer     => $left_seq,
         -right_primer    => $right_seq,
    );

    # add data to the the PrimedSeq object that's not specific to the Primers
    for my $key (%$results) {
        next if ($key =~ /PRIMER_(LEFT|RIGHT)/i );
            $primed_seq->add_tag_value($key, $$results{$key});
    }

    $self->{'next_to_return'}++;
    return $primed_seq;
}


=head2 primer_stream

  Title   : primer_stream()
  Usage   : while (my $primed_seq  = $primer3->primer_stream()) {
  Function: Retrieve the primer/sequences one at a time
  Returns : Returns a Bio::Seq::PrimedSeq object, one at a time
  Args    : None
  Notes   : Deprecated, just a link to next_primer

=cut

sub primer_stream {
    my $self = shift;
    my $primedseq = $self->next_primer;
    return $primedseq;
}


=head2 _separate

  Title   : _separate()
  Usage   : $self->_separate();
  Function: An internal function that groups the results by number
            (e.g. primer pair 1, etc)
  Returns : Nothing.
  Args    : None
  Notes   :

=cut

sub _separate {
    my $self = shift;
    my %results; # the results that we find
    my $maxlocation = -1; # the maximum number of primers returned
    foreach my $key (keys %{$self->{'results'}}) {
        next if (${$self->{'input_options'}}{$key}); # don't process it if it is an input key

        my $location; # the number of the primer pair
        # names will have values like
        # PRIMER_RIGHT_SEQUENCE, PRIMER_RIGHT_2_SEQUENCE, PRIMER_PRODUCT_SIZE, and
        # PRIMER_PRODUCT_SIZE_3 hence we need to find and remove the number
        my $tempkey = $key;
        if ($tempkey =~ s/_(\d+)//) {
            $location = $1;
            if ($location > $maxlocation) {$maxlocation = $location}
        } elsif ( $tempkey =~ /PRIMER_(RIGHT|LEFT)_SEQUENCE/ ) {
            # first primers reported without a number, therefore set $location to 0
            $location = 0;
            if ($location > $maxlocation) {$maxlocation = $location}
        } else {
            $location = 0;
        }
        # we will hash the results by number, and then by name
        ${$results{$location}}{$tempkey}=${$self->{'results'}}{$key};
    }
    $self->{'results_by_number'} = \%results;
    $self->{'maximum_primers_returned'} = $maxlocation;
}


=head2 _set_variable

  Title   : _set_variable()
  Usage   : $self->_set_variable('variable name', 'value');
  Function: An internal function that sets a variable
  Returns : Nothing.
  Args    : None
  Notes   : Used to set $self->{results} and $self->seqobject

=cut

sub _set_variable {
    my ($self, $name, $value) = @_;
    next unless ($name);
    $self->{$name} = $value;
}

1;

__END__

