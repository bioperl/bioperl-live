# $Id$
#
# bioperl module for Bio::Coordinate::GeneMapper
#
# Cared for by Heikki Lehvaslaiho <heikki@ebi.ac.uk>
#
# Copyright Heikki Lehvaslaiho
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Coordinate::GeneMapper - transformations between gene related coordinate systems

=head1 SYNOPSIS

  # to use
  use Bio::Coordinate::GeneMapper;

  # defaults to ID 1 "Standard"
  $gmap = Bio::Coordinate::GeneMapper->new();
  $gmap = Bio::Coordinate::GeneMapper -> new (-in => 'chr',
                                              -out=> 'cds',
                                              -strict => 1,
                                              -nozero => 0 );


=head1 DESCRIPTION


Bio::Coordinate::GeneMapper is a module for simple mapping of
coodinate locations between various gene related locations in human
genetics. It uses a relaxed form of Bio::Coordinate::Pair called
Bio::Coordinate::ExtrapolatingPair which in addition to freely
extrapolaiting values beyond boundaries disallows the use of location
zero.

It understands by name the following coordinate systems and mapping
between them:

                         peptide (peptide length)
                            ^
                            | -peptide_offset
                            |
                  (frame) propeptide (propeptide length)
                       ^    ^
                        \   |
            translate    \  |
                          \ |
                           cds  (transcript start and end)
                            ^
        negative_introns    | \
                    ^       |  \
`                    \      |   \
        inex  =  introns  exons  \
                       ^    ^    /
                        \   |   /
               splice    \  |  /
                          \ | /
                          gene (gene length)
                            ^
                            | - gene_offset
                            |
                           chr (or entry)


Of these, two operations are special cases, translate and splice.
Translating and reverse translating are implemented as internal
methods that do the simple 1E<lt>-E<gt>3 conversion. Splicing needs
additional information that in BioPerl is represented by
Bio::SeqFeature::Gene::GeneStructureI modules. Splicing depends on
method exons() which takes in a more general array of Bio::LocationI
objects.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to the
Bioperl mailing lists  Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://bio.perl.org/MailList.html     - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

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

package Bio::Coordinate::GeneMapper;
use vars qw(@ISA %COORDINATE_SYSTEMS  %COORDINATE_INTS $TRANSLATION);
use strict;

# Object preamble - inherits from Bio::Root::Root
use Bio::Root::Root;
use Bio::Coordinate::Result;

@ISA = qw(Bio::Root::Root);

# first set internal values for all translation tables

%COORDINATE_SYSTEMS = (
		       peptide          => 7,
		       propeptide       => 6,
		       cds              => 5,
		       negative_introns => 4, #
		       inex             => 3, #
		       gene             => 2,
		       chr              => 1
		      );

%COORDINATE_INTS = (
		    7 => 'peptide',
		    6 => 'propeptide',
		    5 => 'cds',
		    4 => 'negative_introns',
		    3 => 'inex',
		    2 => 'gene',
		    1 => 'chr'
		   );

$TRANSLATION =  $COORDINATE_SYSTEMS{'cds'}. $COORDINATE_SYSTEMS{'propeptide'};


sub new {
    my($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my($in, $out, $peptide_offset, $transcript, $exons, $cds, $strict) =
	$self->_rearrange([qw(IN
                              OUT
                              PEPTIDE_OFFSET
                              TRANSCRIPT
                              EXONS
                              CDS
                              STRICT
			     )],
			 @args);

    # direction of mapping
    $self->{_direction} = 1;

    $in  && $self->in($in);
    $out  && $self->out($out);
    $cds && $self->cds($cds);
    $exons  && $self->exons($exons);
    $transcript  && $self->transcript($transcript);
    $peptide_offset && $self->peptide_offset($peptide_offset);

    $strict && $self->strict($strict);

    return $self; # success - we hope!
}

=head2 in

 Title   : in
 Usage   : $obj->in('peptide');
 Function: Set and read the input coordinate system.
 Example :
 Returns : value of input system
 Args    : new value (optional)

=cut

sub in {
   my ($self,$value) = @_;
   if( defined $value) {
       $self->throw("Not a valid input coordinate system name [$value]\n".
		    "Valid values are ". join(", ", keys %COORDINATE_SYSTEMS ))
	   unless defined $COORDINATE_SYSTEMS{$value};

       $self->{'_in'} = $COORDINATE_SYSTEMS{$value};
   }
   #return $self->{'_in'};
   return $COORDINATE_INTS{ $self->{'_in'} };
}


=head2 out

 Title   : out
 Usage   : $obj->out('peptide');
 Function: Set and read the output coordinate system.
 Example :
 Returns : value of output system
 Args    : new value (optional)

=cut

sub out {
   my ($self,$value) = @_;
   if( defined $value) {
       $self->throw("Not a valid input coordinate system name [$value]\n".
		    "Valid values are ". join(", ", keys %COORDINATE_SYSTEMS ))
	   unless defined $COORDINATE_SYSTEMS{$value};

       $self->{'_out'} = $COORDINATE_SYSTEMS{$value};
   }
   #return $self->{'_out'};
   return $COORDINATE_INTS{ $self->{'_out'} };
}

=head2 strict

 Title   : strict
 Usage   : $obj->strict('peptide');
 Function: Set and read weather strict boundaried of coordinate
           systems are enforced.
           When strict is on, the end of the coordinate range must be defined.
 Example :
 Returns : boolean
 Args    : boolean (optional)

=cut

sub strict {
   my ($self,$value) = @_;
   if( defined $value) {
       $value ? ( $self->{'_strict'} = 1 ) : ( $self->{'_strict'} = 0 );
       ## update in each mapper !!
   }
   return $self->{'_strict'} || 0 ;
}

=head2 peptide

 Title   : peptide
 Usage   : $obj->peptide_offset($peptide_coord);
 Function: Read and write the offset of peptide from the start of propeptide
           and peptide length
 Returns : a Bio::Location::Simple object
 Args    : a Bio::LocationI object

=cut

sub peptide {
   my ($self, $value) = @_;
   if( defined $value) {
       $self->throw("I need a Bio::LocationI, not  [". $value. "]")
	   unless $value->isa('Bio::LocationI');

       $self->throw("Peptide start not defined")
	   unless defined $value->start;
       $self->{'_peptide_offset'} = $value->start - 1;

       $self->throw("Peptide end not defined")
	   unless defined $value->end;
       $self->{'_peptide_length'} = $value->end - $self->{'_peptide_offset'};


       my $a = $self->_create_pair
	   ('propeptide', 'peptide', $self->strict, 
	    $self->{'_peptide_offset'}, $self->{'_peptide_length'} );
       my $mapper =  $COORDINATE_SYSTEMS{'propeptide'}. + $COORDINATE_SYSTEMS{'peptide'};
       $self->{'_mappers'}->{$mapper} = $a;
   }
   return  Bio::Location::Simple->new
       (-seq_id => 'propeptide', 
	-start => $self->{'_peptide_offset'} + 1 ,
	-end => $self->{'_peptide_length'} + $self->{'_peptide_offset'},
	-strand => 1
       );
}

=head2 peptide_offset

 Title   : peptide_offset
 Usage   : $obj->peptide_offset(20);
 Function: Set and read the offset of peptide from the start of propeptide
 Returns : set value or 0
 Args    : new value (optional)

=cut

sub peptide_offset {
   my ($self,$offset, $len) = @_;
   if( defined $offset) {
       $self->throw("I need an integer, not [$offset]")
	   unless $offset =~ /^[+-]?\d+$/;
       $self->{'_peptide_offset'} = $offset;

       if (defined $len) {
	   $self->throw("I need an integer, not [$len]")
	       unless $len =~ /^[+-]?\d+$/;
	   $self->{'_peptide_length'} = $len;
       }

       my $a = $self->_create_pair
	   ('propeptide', 'peptide', $self->strict, $offset, $self->{'_peptide_length'} );
       my $mapper =  $COORDINATE_SYSTEMS{'propeptide'}. + $COORDINATE_SYSTEMS{'peptide'};
       $self->{'_mappers'}->{$mapper} = $a;
   }
   return $self->{'_peptide_offset'} || 0;
}

=head2 peptide_length

 Title   : peptide_length
 Usage   : $obj->peptide_length(20);
 Function: Set and read the offset of peptide from the start of propeptide
 Returns : set value or 0
 Args    : new value (optional)

=cut


sub peptide_length {
   my ($self, $len) = @_;
   if( defined $len) {
       $self->throw("I need an integer, not [$len]")
	   if defined $len && $len !~ /^[+-]?\d+$/;
       $self->{'_peptide_length'} = $len;
   }
   return $self->{'_peptide_length'};
}


=head2 transcript

 Title   : transcript
 Usage   : $obj->transcript(20);
 Function: Set and read the offset of CDS from the start of transcipt


 Returns : set value or 0
 Args    : new value (optional)

=cut

sub transcript {
   my ($self,$value) = @_;
   if( defined $value) {
       $self->throw("I need an integer, not [$value]")
	   unless $value =~ /^[+-]?\d+$/;
       $self->{'_transcript'} = $value;

       my $a = $self->_create_pair('cds', 'cds', 0, $value );
       my $mapper =  $COORDINATE_SYSTEMS{'cds'}. + $COORDINATE_SYSTEMS{'cds'};
       $self->{'_mappers'}->{$mapper} = $a;
   }
   return $self->{'_transcript'} || 0;
}


=head2 exons

 Title   : exons
 Usage   : $obj->exons(\@exons);
 Function: Set and read the offset of CDS from the start of transcipt
           You do not have to sort the exons before calling this method as
           they will be sorted automatically
 Returns : array of Bio::LocationI exons in genome coordinates or 0
 Args    : array of Bio::LocationI exons in genome (or entry) coordinates

=cut

sub exons {
   my ($self,@value) = @_;
   my $mapper =  $COORDINATE_SYSTEMS{'gene'}. + $COORDINATE_SYSTEMS{'cds'};
   my $intron_exon_mapper =
       $COORDINATE_SYSTEMS{'gene'}. + $COORDINATE_SYSTEMS{'inex'};
   my $negative_intron_mapper =
       $COORDINATE_SYSTEMS{'gene'}. + $COORDINATE_SYSTEMS{'negative_introns'};

   if(@value) {
       $self->throw("I need an array , not [@value]")
	   unless ref \@value eq 'ARRAY';
       $self->throw("I need an  reference to an array of Bio::LocationIs, not to [".
		    $value[0]. "]")
	   unless ref $value[0] and $value[0]->isa('Bio::LocationI');

       # sort the input array
       if (defined $value[0]->strand && $value[0]->strand == - 1) {  #reverse strand
	   @value = map { $_->[0] }
	               sort { $b->[1] <=> $a->[1] }
                       map { [ $_, $_->start] }
                       @value;
       } else {               #undef or forward strand
	   @value = map { $_->[0] }
	               sort { $a->[1] <=> $b->[1] }
                       map { [ $_, $_->start] }
                       @value;
       }

       $self->{'_chr_exons'} = \@value;


       # transform exons from chromosome to gene coordinates
       my $tmp_in = $self->in;
       my $tmp_out = $self->out;

       $self->in('chr');
       $self->out('gene');
       my @exons = map {$self->map($_)} @value;

       $self->in($tmp_in);
       $self->out($tmp_out);
       #

       my $map = Bio::Coordinate::Collection->new;
       my $exon_map = Bio::Coordinate::Collection->new;
       my $negative_intron_map = Bio::Coordinate::Collection->new;

       my $tr_end = 0;
       my $coffset;
       my $exon_counter;
       my $prev_exon_end;

       for my $exon ( @exons ) {

	   $exon_counter++;

	   #
	   # gene -> cds
	   #

	   unless (defined $tr_end) {
	       $tr_end = $exon->start - 1 ;
	   }
	   #print "--------------------------------->$tr_end\n";
	   my $match1 = Bio::Location::Simple->new
	       (-seq_id =>'gene' ,
		-start => $exon->start,
		-end => $exon->end, -strand=>1 );
	   my $match2 = Bio::Location::Simple->new
	       (-seq_id => 'cds',
		-start => $exon->start - $tr_end ,
		-end => $exon->end -  $tr_end,
		-strand=>$exon->strand );

	   my $pair = Bio::Coordinate::Pair->new(-in => $match1,
						 -out => $match2,
						);
	   $map->add_mapper($pair);

	   if ($exon->start <= 1 and $exon->end >= 1) {
	       $coffset = $tr_end ;
	   }

#	   if ($exon->start <= 1) {
#	       $coffset += $exon->start - $prev_exon_end;
#	   }

	   $tr_end = $tr_end  + $exon->end - $exon->start + 1;

	   #
	   # gene -> intron
	   #

	   if (defined $prev_exon_end) {
	       my $match3 = Bio::Location::Simple->new
		   (-seq_id =>'gene',
		    -start => $prev_exon_end + 1,
		    -end => $exon->start -1, -strand=>1 );

	       my $match4 = Bio::Location::Simple->new
		   (-seq_id => 'intron'. ($exon_counter -1),
		    -start => 1,
		    -end => $exon->start - 1 - $prev_exon_end,
		    -strand=>$exon->strand );

	       # negetive intron coordinates
	       my $match6 = Bio::Location::Simple->new
		   (-seq_id => 'intron'. ($exon_counter -1),
		    -start => -1 * ($exon->start - 2 - $prev_exon_end),
		    -end => 0,
		    -strand=>$exon->strand );

	       my $pair = Bio::Coordinate::Pair->new(-in => $match3,
						     -out => $match4
						    );
	       my $neg_pair = Bio::Coordinate::Pair->new(-in => $match3,
							 -out => $match6
							);

	       $exon_map->add_mapper($pair);
	       $negative_intron_map->add_mapper($neg_pair);

	   }

	   # store the value
	   $prev_exon_end = $exon->end;

	   #
	   # gene -> exon
	   #
	   my $match5 = Bio::Location::Simple->new
	       (-seq_id => 'exon'. $exon_counter,
		-start => 1,
		-end => $exon->end - $exon->start +1,
		-strand=> $exon->strand );

	   my $pair2 = Bio::Coordinate::Pair->new(-in => $match1,
						  -out => $match5
						 );
	   $exon_map->add_mapper($pair2);

       }

       #$coffset = $tr_offset unless $coffset;
       # move coordinate start if exons have negative values
       if ($coffset) {
	   foreach my $m ($map->each_mapper) {
	       $m->out->start($m->out->start + $coffset);
	       $m->out->end($m->out->end + $coffset);
	   }

       }

       $self->{'_mappers'}->{$mapper} = $map;
       $self->{'_mappers'}->{$intron_exon_mapper} = $exon_map;
       $self->{'_mappers'}->{$negative_intron_mapper} = $negative_intron_map;
   }
   return  @{$self->{'_chr_exons'}}  || 0;
}


=head2 cds

 Title   : cds
 Usage   : $obj->cds(20);
 Function: Set and read the offset of CDS from the start of transcipt

           Simple input can be an integer which gives the start of the
           coding region in genomic coordinate. If you want to provide
           the end of the coding region or indicate the use of the
           opposite strand, you have to pass a Bio::Location::Simple
           object to this method.

 Returns : set value or 0
 Args    : new value (optional)

=cut

sub cds {
   my ($self,$value) = @_;
   if( defined $value) {
       my ($start, $end);

       if ($value =~ /^[+-]?\d+$/ ) {
	   my $loc = Bio::Location::Simple->new(-start=>$value);
	   $self->{'_cds'} = $loc;
       }
       elsif (ref $value &&  $value->isa('Bio::LocationI') ) {
	   $self->{'_cds'} = $value;
       } else {
	   $self->throw("I need an integer or Bio::LocationI, not [$value]")
       }
       # strand !!
       my $len;

       $len = $self->{'_cds'}->end - $self->{'_cds'}->start +1
	   if defined $self->{'_cds'}->end;

       my $a = $self->_create_pair
	   ('chr', 'gene', 0,
	    $self->{'_cds'}->start-1,
	    $len,
	    $self->{'_cds'}->strand);
       my $mapper =  $COORDINATE_SYSTEMS{'chr'}. + $COORDINATE_SYSTEMS{'gene'};
       $self->{'_mappers'}->{$mapper} = $a;

       # recalculate exon-based mappers
       if ( defined $self->{'_chr_exons'} ) {
	   $self->exons(@{$self->{'_chr_exons'}});
       }

   }
   return $self->{'_cds'} || 0;
}


=head2 map

 Title   : map
 Usage   : $newpos = $obj->map(5);
 Function: Map the location from the input coordinate system
           to a new value in the output coordinate system.
 Example :
 Returns : new value in the output coordiante system
 Args    : a Bio::Location::Simple

=cut

sub map {
   my ($self,$value) = @_;
   my ($res);

   $self->throw("Need to pass me a Bio::Location::Simple")
       unless defined $value;
   $self->throw("Need to pass me a Bio::Location::Simple, not [".
		ref($value). "]")
       unless defined $value->isa('Bio::Location::Simple');
   $self->throw("Input coordinate system not set")
       unless $self->{'_in'};
   $self->throw("Output coordinate system not set")
       unless $self->{'_out'};
   $self->throw("Do not be silly. Input and output coordinate ".
		"systems are the same!")
       unless $self->{'_in'} != $self->{'_out'};

   $self->_check_direction();

   print STDERR "=== Start location: ". $value->start. ",".
       $value->end. " (". $value->strand. ")\n" if $self->verbose > 0;

   my $counter = $self->{'_in'};
   while ($counter != $self->{'_out'}) {

       my $mapper;
       if ($self->direction == 1 ) {
	   if ($COORDINATE_INTS{$counter} eq 'gene') {
	       if ( $self->{'_out'} > ($counter+2) ) {
		   $mapper = "$counter". ($counter+3);
	       }
	       elsif ( $self->{'_out'} == ($counter+2) ) {
		   $mapper = "$counter". ($counter+2);
	       } else {
		   $mapper = "$counter". ($counter+1);
	       }
	   } else { #default
	       $mapper = "$counter". ($counter+1);
	   }
       } else { # direction = -1
	   if ($COORDINATE_INTS{$counter} eq 'cds') {
	       if ( $self->{'_out'} < ($counter-2) ) {
		   $mapper = ($counter-3). "$counter";
	       }
	       elsif ( $self->{'_out'} == ($counter-2) ) {
		   $mapper = ($counter-2). "$counter";
	       } else { # cds -> negetive_intron
		   $self->throw(" 'cds' -> 'negative_intron' is not possible mapping!");
	       }
	   } else {
	       $mapper = ($counter-1). "$counter";
	   }
       }
       print "=================================>$counter:$mapper(", $self->{'_out'},")\n" if $self->verbose > 0;
       # handle exception : translation
       #  cds -> propeptide
       if ($mapper == $TRANSLATION) {
	   if ($self->direction == 1) {
	       $value = $self->_translate($value);
	       print STDERR "+   cds -> propeptide (translate) "
		   if $self->verbose > 0;
	   } else {
	       $value = $self->_reverse_translate($value);
	       print STDERR "+   propeptide -> cds (reverse translate) "
		   if $self->verbose > 0;
	   }
	   next;
       }

       # keep the start and end values, and go on to next iteration
       #  if this mapper is not set
       unless ( defined $self->{'_mappers'}->{$mapper} ) {
	   # update mapper name
	   $value->seq_id($COORDINATE_INTS{$counter});
	   print STDERR "-   ". $COORDINATE_INTS{$counter}. " -> ".
	       $COORDINATE_INTS{$counter+ $self->direction}. " " if $self->verbose > 0;
	   next;
       }

       # generic mapping
       #       print "counter = $counter|$mapper\n";
       my $res = $self->{'_mappers'}->{$mapper}->map($value);

       my @matches = $res->each_match;
       if (@matches > 1) {
	   $self->throw('Multiple matches in different coordinate systems not handled')
	       if $matches[0]->seq_id ne $matches[-1]->seq_id;
	   $value = $res->match;
	   $value->start($matches[0]->start);

       } else {
	   $value = $res->match;
       }

#       use Data::Dumper;
#       print Dumper $res;
       return undef  unless  $value;
       print STDERR "+   ". $COORDINATE_INTS{$counter}. " -> ".
	   $COORDINATE_INTS{$counter+ $self->direction}. " " if $self->verbose > 0;


   } continue {
       # move counter
       print STDERR "    ". $value->start. ",".
	   $value->end. " (". $value->strand. ")\n" if $self->verbose > 0;


       $counter = $counter + 2 
	   if $COORDINATE_INTS{$counter} eq 'gene'
	      and $self->direction == 1
	      and $self->{'_out'} > ($counter+2);

       $self->direction == 1 ? $counter++ : $counter-- ;

   }

   return $value;

}

=head2 direction

 Title   : direction
 Usage   : $obj->direction('peptide');
 Function: Read-only method for the direction of mapping deduced from
           predefined input and output coordinate names.
 Example :
 Returns : 1 or -1, mapping direction
 Args    : new value (optional)

=cut

sub direction {
   my ($self) = @_;
   return $self->{'_direction'};
}


=head2 swap

 Title   : swap
 Usage   : $obj->swap;
 Function: Swap the direction of transformation
           (input <-> output)
 Example :
 Returns : 1
 Args    : 

=cut

sub swap {
   my ($self,$value) = @_;
   my ($tmp);

   $tmp = $self->{'_out'};;
   $self->{'_out'} = $self->{'_in'};
   $self->{'_in'} = $tmp;
   foreach my $map (keys %{$self->{'_mappers'}}){
       $self->{'_mappers'}->{$map}->swap;
   }

   # record the changed direction;
   $self->{_direction} *= -1;

   return 1;
}


=head2 to_string

 Title   : to_string
 Usage   : $newpos = $obj->to_string(5);
 Function: Dump the mapping values into human readable format
 Example :
 Returns : string
 Args    : 

=cut

sub to_string {
   my ($self) = shift;

   print "-" x 40, "\n";
   my $counter = 1;
   while (defined $COORDINATE_INTS{$counter} && $COORDINATE_INTS{$counter+1}) {
       my $in = $COORDINATE_INTS{$counter};
       my $out = $COORDINATE_INTS{$counter+1};
       my $mapper = $counter. ($counter+1);


       printf "\n%12s -> %-12s (%s)\n", $in, $out, $mapper
	   unless $mapper eq '34' or $mapper eq '45';


       if ($mapper eq $TRANSLATION) {
	   printf "%9s%-12s\n", "", '"translate"';
       }
       elsif ($mapper eq '34' or $mapper eq '45') {
	   next;
       }
       elsif ($mapper eq '23') {
		   next unless $self->{'_mappers'}->{$mapper};
	   	   printf "%10s%-12s\n", "", '"splice"';

		   my $i = 1;
		   my $ex = 1;

		   my $inex;
		   foreach my $pair ( $self->{'_mappers'}->{$mapper}->each_mapper ) {
		       $i%2 == 1 ? ($inex = 'exon') : ($inex = 'intron');
		       printf "%2s :%8s -> %-12s%-6s\n", $i, $pair->in->start, $pair->out->start, "$inex$ex" ;
		       printf "%2s :%8s -> %-12s\n", '', $pair->in->end, $pair->out->end ;
		       $i++; 
		       $ex++ if $inex eq 'exon';
		   }

		   printf "\n%12s -> %-12s (%s)\n", $in, 'cds', '25';
	   	   printf "%10s%-12s\n", "", '"splice"';
		   $i = 1;
		   foreach my $pair ( $self->{'_mappers'}->{'25'}->each_mapper ) {
		       printf "%2s :%8s -> %-12s\n", $i, $pair->in->start, $pair->out->start ;
		       printf "%2s :%8s -> %-12s\n", '', $pair->in->end, $pair->out->end ;
		       $i++;
		   }


       }
       elsif (not defined $self->{'_mappers'}->{$mapper}) {
	   printf "%12s%-12s\n", "", 'undef';
       } else {
	   if ($mapper eq '12') {
	       my $end= 'undef';
	       $end = $self->cds->end -1 if defined $self->cds->end;
	       printf "%16s%s: %s (%s)\n", ' ', 'gene offset', $self->cds->start-1 , $end;
	       printf "%16s%s: %s\n", ' ', 'gene strand', $self->cds->strand;

	   }
	   elsif ($mapper eq '34') {
	       printf "%16s%s: %s\n", ' ', "transcript (transcipt offset)", $self->transcript;
	   }
	   elsif ($mapper eq '56') {
	       printf "%16s%s: %s\n", ' ', "peptide offset", $self->peptide_offset;
	   }


	   printf "%12s -> %-12s\n",
	       $self->{'_mappers'}->{$mapper}->{'_in'}->start, 
		   $self->{'_mappers'}->{$mapper}->{'_out'}->start;
       }
   } continue {
       $counter++;
   }

   print "\nin : ", $self->in, "\n";
   print "out: ", $self->out, "\n";
   my $dir;
   $self->direction ? ($dir='forward') : ($dir='reverse');
   printf "direction: %-8s(%s)\n",  $dir, $self->direction;
   print "\n", "-" x 40, "\n";

}


=head2 _create_pair

 Title   : _create_pair
 Usage   : $mapper = $obj->_create_pair('chr', 'gene', 0, 2555, 10000, -1);
 Function: Internal helper method to create a mapper between 
           two coordinate systems
 Returns : a Bio::Coordinate::Pair object
 Args    : string, input coordinate system name,
           string, output coordinate system name,
           boolean, strict mapping
           positive integer, offset
           positive integer, length
           1 || -1 , strand

=cut

sub _create_pair {
   my ($self, $in, $out, $strict, $offset, $length, $strand ) = @_;
   $strict ||=0;
   $strand ||=1;
   $length ||= 20;

   my $match1 = Bio::Location::Simple->new
       (-seq_id => $in,
	-start => $offset+1,
	-end => $offset+$length, -strand=>1 );

   my $match2 = Bio::Location::Simple->new
       (-seq_id => $out,
	-start => 1,
	-end => $length, -strand=>$strand );

   my $pair = Bio::Coordinate::ExtrapolatingPair->
       new(-in => $match1,
	   -out => $match2,
	   -strict => $strict
	  );

   return $pair;

}


=head2 _translate

 Title   : _translate
 Usage   : $newpos = $obj->_translate(5);
 Function: Translate the location from the CDS coordinate system
           to a new value in the propeptide coordinate system.
 Example :
 Returns : new location
 Args    : a Bio::Location::Simple

=cut

sub _translate {
   my ($self,$value) = @_;

   $self->throw("Need to pass me a Bio::Location::Simple, not [".
		ref($value). "]")
       unless defined $value->isa('Bio::Location::Simple');

   my $loc = new Bio::Location::Simple;
   $loc->start(int($value->start / 3 )+1);
   $loc->end(int($value->end / 3 )+1);
   $loc->seq_id('propeptide');
   $loc->strand(1);
   return $loc;
}

sub _frame {
   my ($self,$value) = @_;

   $self->throw("Need to pass me a Bio::Location::Simple, not [".
		ref($value). "]")
       unless defined $value->isa('Bio::Location::Simple');

   my $loc = new Bio::Location::Simple;
   $loc->start( ($value->start-1) % 3 +1);
   $loc->end( ($value->end-1) % 3 +1);
   $loc->seq_id('frame');
   $loc->strand(1);
   return $loc;
}


=head2 _reverse_translate

 Title   : _reverse_translate
 Usage   : $newpos = $obj->_reverse_translate(5);
 Function: Reverse translate the location from the propeptide
           coordinate system to a new value in the CSD.
           Note that a single peptide location expands to cover
           the codon triplet
 Example :
 Returns : new location in the CDS coordinate system
 Args    : a Bio::Location::Simple

=cut

sub _reverse_translate {
   my ($self,$value) = @_;

   $self->throw("Need to pass me a Bio::Location::Simple, not [".
		ref($value). "]")
       unless defined $value->isa('Bio::Location::Simple');

   my $loc = new Bio::Location::Simple;
   $loc->start($value->start * 3 - 2);
   $loc->end($value->end * 3 );
   $loc->seq_id('cds');
   $loc->strand(1);
   return $loc;

}


=head2 _check_direction

 Title   : _check_direction
 Usage   : $obj->_check_direction();
 Function: Check and swap when needed the direction the location
           mapping Pairs based on input and output values
 Example :
 Returns : new location
 Args    : a Bio::Location::Simple

=cut

sub _check_direction {
   my ($self) = @_;

   my $new_direction = 1;
   $new_direction = -1 if $self->{'_in'} > $self->{'_out'};

   unless ($new_direction == $self->{_direction} ) {
       foreach my $map (keys %{$self->{'_mappers'}}){
	   $self->{'_mappers'}->{$map}->swap;
       }
       # record the changed direction;
       $self->{_direction} *= -1;
   }
   1;
}


1;
