#
# Copyright (c) 1997-2001 bioperl, Chad Matsalla. All Rights Reserved.
#           This module is free software; you can redistribute it and/or
#           modify it under the same terms as Perl itself. 
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::SeqIO::Primer3 - Create input for and work with the output from the 
program primer3

=head1 SYNOPSIS

Do not use this module directly. Use it via the Bio::SeqIO class, see
L<Bio::SeqIO> for more information.

=head1 DESCRIPTION

Bio::SeqIO::Primer3 creates the input files needed to design primers using
primer3 and provides mechanisms to access data in the primer3 output files.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org          - General discussion
  http://www.bioperl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/


=head1 AUTHOR - Chad Matsalla

bioinformatics1@dieselwurks.com

=head1 CONTRIBUTORS

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::SeqIO::primer3;

use vars qw(@ISA);
use strict;
use Bio::SeqIO;
# use Bio::SeqFeature::Primer;
use Bio::Seq::SeqFactory;
use Dumpvalue;

@ISA = qw(Bio::SeqIO);

my $dumper = new Dumpvalue;


sub _initialize {
  my($self,@args) = @_;
  $self->SUPER::_initialize(@args);
  if( ! defined $self->sequence_factory ) {
      $self->sequence_factory(new Bio::Seq::SeqFactory
                     (-verbose => $self->verbose(),
                      -type => 'Bio::Seq'));
  }
}





=head2 next_seq()

 Title   : next_seq()
 Usage   : $primer3 = $stream->next_seq()
 Function: returns the next scf sequence in the stream
 Returns : Bio::SeqFeature::Primer object
 Args    : NONE (huh?)
 Notes   : Fills the interface specification for SeqIO.

=cut

sub next_seq {
     my $self = shift;
     my $fh = $self->_filehandle();
     my ($line,%primer);
          # first, read in the next set of primers
     while ($line = $self->_readline()) {
          chomp ($line);
          last if ($line =~ /^=/);
          $line =~ m/(^.*)\=(.*$)/;
          %primer->{$1} = $2; 
     }
          # then, get the primers as SeqFeature::Primer objects
     # chad! fix this!!!
     # my ($left,$right) = &_create_primer_features(\%primer);
          # then, create the sequence to place them on
     my $sequence = Bio::Seq->new(-seq => %primer->{SEQUENCE},
                                         -display_id => %primer->{PRIMER_SEQUENCE_ID});

          # now, place the primers on the sequence
     # $sequence->add_SeqFeature($left,$right);
     # $dumper->dumpValue($sequence);

     return $sequence;
}



sub _create_primer_features {
     my $rdat = shift;
     my (%left,%right,$updir,$downdir,$var,$trunc);
     my @variables = qw(
        PRIMER_DIRECTION
        PRIMER_DIRECTION_END_STABILITY
        PRIMER_DIRECTION_EXPLAIN
        PRIMER_DIRECTION_GC_PERCENT
        PRIMER_DIRECTION_PENALTY
        PRIMER_DIRECTION_SELF_ANY
        PRIMER_DIRECTION_SELF_END
        PRIMER_DIRECTION_SEQUENCE
        PRIMER_DIRECTION_TM
         PRIMER_FIRST_BASE_INDEX
     );
          # create the hash to pass into the creation routine
     foreach $updir (qw(LEFT RIGHT)) {
            my %dat;
             foreach (@variables) {
               ($var = $_) =~ s/DIRECTION/$updir/e;
                 if (/^PRIMER_DIRECTION$/) {
                      $trunc = "colen";
                 }
                 elsif (/^PRIMER_FIRST_BASE_INDEX/) {
                      $trunc = "first_base_index";    
                 }
                 else {
                      ($trunc = $_) =~ s/PRIMER_DIRECTION_//;
                      $trunc =~ tr/[A-Z]/[a-z]/;
                 }
               %dat->{-$trunc} = $rdat->{$var};
             }
               if ($updir eq "LEFT") {
                    %left = %dat;
                    %left->{-primer_sequence_id} = $rdat->{PRIMER_SEQUENCE_ID}."-left";
               }
               else {
                    %right = %dat;
                    %right->{-primer_sequence_id} = $rdat->{PRIMER_SEQUENCE_ID}."-right";
               }
     }
     my $primer_left = new Bio::SeqFeature::Primer(%left);
     my $primer_right = new Bio::SeqFeature::Primer(%right);
     return($primer_left,$primer_right);
}









=head2 get_amplified_region()

 Title   : get_amplified_region()
 Usage   : 
 Function: 
 Returns : 
 Args    : 
 Notes   : 

=cut

sub get_amplified_region {
	my $primer = $_[1];
	return $Primer3::primers{$primer}{amplified};
} # end get_amplified_region

=head2 get_amplification_error()

 Title   : get_amplification_error()
 Usage   : 
 Function: 
 Returns : 
 Args    : 
 Notes   : 

=cut

sub get_amplification_error {
	my $primer = $_[1];
	my $error = $Primer3::primers{$primer}{PRIMER_ERROR};
	if ($error) { return $error; }
	else { return "Some error that primer3 didn't define.\n"; }
}

=head2 _set_target()

 Title   : _set_target()
 Usage   : &_set_target($self);
 Function: 
 Returns : 
 Args    : 
 Notes   : 

=cut

sub _set_target {
	my $self = shift;
	my ($sequence,$primer,$primer_left,$primer_right,$position_left,$position_right,$boggle);
	$boggle = 1;
	foreach $primer (sort keys %{$self->{primers}}) {
		$sequence = $self->{primers}{$primer}{SEQUENCE};
		$primer_left = $self->{primers}{$primer}{PRIMER_LEFT};
		$primer_right = $self->{primers}{$primer}{PRIMER_RIGHT};
		if (!$primer_left) {
			$self->{primers}{$primer}{design_failed} = "1";
		}
		else {
			$primer_left =~ m/(.*)\,(.*)/;
				# print("left: \$1 is $1 and \$2 is $2\n");
			$position_left = $1+$2-1;
			$primer_right =~ m/(.*)\,(.*)/;
			$position_right = $1-$2;
			$self->{primers}{$primer}{left} = $position_left;
			$self->{primers}{$primer}{right} = $position_right;
			$self->{primers}{$primer}{amplified} = substr($sequence,$position_left,$position_right-$position_left);
		}
	}
}

=head2 _read_file($self,$filename)

 Title   : _read_file($self,$filename)
 Usage   : 
 Function: 
 Returns : A scalar containing the contents of $filename
 Args    : $self and the name of a file to parse.
 Notes   : 

=cut

sub _read_file {
	my ($self,$filename) = @_;
		# print("in _parse_report, \$self is $self\n");
		# set this to keep track of things....
	$self->{outfilename} = $filename;
	my $fh = new FileHandle;
	open($fh,$filename) or die "I can't open the primer report ($filename) : $!\n";
		# _parse_report();
		# my %Primer3::primers;
	my ($output,$line);
	while ($line=<$fh>) {
			# print("Adding $line\n");
		$output .= $line;
	} # end while
		# print("\$output is $output\n");
	return $output;
}





=head2 _parse_report()

 Title   : _parse_report()
 Usage   : &_parse_report($self,$filename);
 Function: Parse a primer3 outfile and place everything into an object under
	{primers} with PRIMER_SEQUENCE_ID being the name of the keys for the
	{primers} hash.
 Returns : Nothing.
 Args    : $self and the name of a file to parse.
 Notes   : 

=cut

sub _parse_report {
		# old
		# my ($self,$filename) = @_;
	my ($self,$outputs) = @_;
		# print("\$self is $self, \$outputs are $outputs\n");
		# print("in _parse_report, \$self is $self\n");
		# set this to keep track of things....
	my ($sequence_name,$line,$counter,$variable_name,$variable_value);
	my @output = split/\n/,$outputs;	
	foreach $line (@output) {
			# print("Reading line $line\n");
		next if ($line =~ /^\=/);
		if ($line =~ m/^PRIMER_SEQUENCE_ID/) {
			$line =~ m/(\S+)=(.*$)/;
			$variable_name = $1;
			$sequence_name = $2;
			$variable_value = $2;
		}
		else {
			$line =~ m/(\S+)=(.*$)/;
			$variable_name = $1;
			$variable_value = $2;
		}
			# print("$sequence_name\t$variable_name\t$variable_value\n");
		$self->{primers}{$sequence_name}{$variable_name} = $variable_value;
	} # end while <>
} # end parse_report

=head2 _construct_empty()

 Title   : _construct_empty()
 Usage   : &_construct_empty($self);
 Function: Construct an empty object that will be used to construct a primer3
	input "file" so that it can be run.
 Returns : 
 Args    : 
 Notes   : 

=cut

sub _construct_empty {
	my $self = shift;
	$self->{inputs} = {};
	return;
}

=head2 add_target(%stuff)

 Title   : add_target(%stuff)
 Usage   : $o_primer->add_target(%stuff);
 Function: Add an target to the infile constructor.
 Returns : 
 Args    : A hash. Looks something like this:
	$o_primer2->add_target(
		-PRIMER_SEQUENCE_ID     =>      "sN11902",
		-PRIMER_COMMENT         =>      "3831",
		-SEQUENCE               =>      "some_sequence",
		-TARGET                 =>      "513,26",
		-PRIMER_PRODUCT_SIZE_RANGE      =>      "100-500",
		-PRIMER_FILE_FLAG       =>      "0",
		-PRIMER_LIBERAL_BASE    =>      "1",
		-PRIMER_NUM_RETURN      =>      "1",
		-PRIMER_FIRST_BASE_INDEX        =>      "1",
		-PRIMER_EXPLAIN_FLAG    =>      "1");
	The add_target() method does not validate the things you put into
	this parameter hash. Read the docs for Primer3 to see which fields
	do what and how they should be used.
 Notes   : To design primers, first create a new CSM::Primer3 object with the
	-construct_infile parameter. Then, add targets using this method
	(add_target()) with the target hash as above in the Args: section.
	Be careful. No validation will be done here. All of those parameters
	will be fed straight into primer3.
	Once you are done adding targets, invoke the function run_primer3().
	Then retrieve the results using something like a loop around the array
	from get_primer_sequence_IDs();

=cut


sub add_target {
	my ($self,%args) = @_;
	my ($currkey,$renamed,$sequence_id,$value);
	if (!$args{-PRIMER_SEQUENCE_ID}) {
		print("You cannot add an element to the primer3 infile without specifying the PRIMER_SEQUENCE_ID. Sorry.\n");
	}
	else {
		$sequence_id = $args{-PRIMER_SEQUENCE_ID};
		foreach $currkey (keys %args) {
			print("\$currkey is $currkey\n");
			next if ($currkey eq "-PRIMER_SEQUENCE_ID");
			($renamed = $currkey) =~ s/-//;
				# print("Adding $renamed to the hash under $sequence_id\n");
			$value = $args{$currkey};
				# print("\$value is $value\n");
			if ($renamed eq "SEQUENCE") { $value =~ s/\n//g; }
			$self->{infile}{$sequence_id}{$renamed} = $value;
		}
	}
}

=head2 get_primer_sequence_IDs()

 Title   : get_primer_sequence_IDs()
 Usage   : $o_phred->get_primer_sequence_IDs();
 Function: Return the primer sequence ID's. These normally correspond to
	the name of a sequence in a database but can be whatever was used when
	the primer3 infile was constructed.
 Returns : An array containing the names of the primer sequence ID's
 Args    : None.
 Notes   : This would be used as the basis for an iterator to loop around each
	primer that was designed.

=cut

sub get_primer_sequence_IDs {
	my $self = shift;
	return sort keys %{$self->{primers}};
} # end get keys

=head2 dump_hash()

 Title   : dump_hash()
 Usage   : $o_primer->dump_hash();
 Function: Dump out the CSM::Primer3 object.
 Returns : Nothing.
 Args    : None.
 Notes   : Used extensively in debugging.

=cut

sub dump_hash {
	my $self = shift;
	my $dumper = new Dumpvalue;
	$dumper->dumpValue($self);
} # end dump_hash

=head2 dump_infile_hash()

 Title   : dump_infile_hash()
 Usage   : $o_primer->dump_infile_hash();
 Function: Dump out the contents of the infile hash.
 Returns : Nothing.
 Args    : None.
 Notes   : Used for debugging the construction of the infile.

=cut

sub dump_infile_hash {
	my $self = shift;
	my $dumper = new Dumpvalue;
	$dumper->dumpValue($self->{infile});
}

=head2 run_primer3()

 Title   : run_primer3()
 Usage   : $o_consed->run_primer3();
 Function: Run primer3 on the stuff in {infile}
 Returns : 
 Args    : 
 Notes   : 

=cut

sub run_primer3() {
	my $self = shift;
	print("Running primer3\n");
	my $primer3_cli = "/usr/bin/primer3";
		# deprecated
		# open PRIMER, "| primer3 > $results_file" or die "can't fork: $!";
	open2(\*README, \*WRITEME, "$primer3_cli");
		# print WRITEME $sequence;
	my ($seq_id,$parm,$result);
	foreach $seq_id (sort keys %{$self->{infile}}) {
		print WRITEME ("PRIMER_SEQUENCE_ID=$seq_id\n");
		foreach $parm (sort keys %{$self->{infile}{$seq_id}}) {
			print WRITEME ("$parm=$self->{infile}{$seq_id}{$parm}\n");
		}
	}
	print WRITEME ("=\n");
	close(WRITEME);
	my $outfh = select(README);
	$| = 1;
	select($outfh);
	while(<README>) {
		$result .= $_;
	}
	_parse_report($self,$result);
	_set_target($self);
}

=head2 _new_deprecated()

 Title   : new()
 Usage   : $o_primer = CSM::Primer3->new( -file => $filename);
	_or_
           $o_primer = CSM::Primer3->new( -construct_infile	=>	1);
 Function: If constructed with -file, reads in qualities and sequence from a
	primer3 output file.
	If constructed with -construct_infile, prepares to create a file that
	will be sent into primer3.
 Returns : A new CSM::Primer3 object
 Args    : Hash, with one field: -file _or_ -construct_infile
 Notes   : I know, modality sucks.

=cut

sub _new_deprecated {
	my ($class,$filename,%args);
        ($class,%args) = @_;
	my $self = {};
		# print("In new(), the self is $self\n");
	bless ($self,$class);
	if ($args{-file}) {
			# we are now going to stay in parse mode
		$self->{mode} = "parse";
		print("Preparing to parse file called $args{-file}...\n");
		my $primer3_output = _read_file($self,$args{-file});
			# print("\$primer3_output is $primer3_output\n");
		_parse_report($self,$primer3_output);
		_set_target($self);
	}
	elsif ($args{-construct_infile}) {
			# we are now going to enter into construct_infile mode
		$self->{mode} = "construct";
		_construct_empty($self);
	}
	else {
			# uh, oh.
		die "Use -file or -construct_infile when constructing a CSM::Primer3 object. Read the pod, luke.\n";
	}
		# set_target();
	return $self;
} # end new



1;
__END__

=head2 null

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 
 Notes   : 

=cut

=head1 SEE ALSO

perl(1).

=cut
