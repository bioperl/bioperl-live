# Bio::Tools::Alignment::Consed
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chad Matsalla
#
# Copyright Chad Matsalla
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Alignment::Consed - A module to work with objects from consed .ace files

=head1 SYNOPSIS

  # a report for sequencing stuff
  my $o_consed = Bio::Tools::Alignment::Consed->new( 
      -acefile => "/path/to/an/acefile.ace.1",
      -verbose => 1);
  my $foo = $o_consed->set_reverse_designator("r");
  my $bar = $o_consed->set_forward_designator("f");

  # get the contig numbers
  my @keys = $o_consed->get_contigs();

  # construct the doublets
  my $setter_doublets = $o_consed->choose_doublets();

  # get the doublets
  my @doublets = $o_consed->get_doublets();

=head1 DESCRIPTION

L<Bio::Tools::Alignment::Consed> provides methods and objects to deal
with the output from the Consed software suite. Specifically,
takes an C<.ace> file and provides objects for the results.

A word about doublets: This module was written to accommodate a large
EST sequencing operation. In this case, EST's were sequenced from the
3' and from the 5' end of the EST. The objective was to find a
consensus sequence for these two reads.  Thus, a contig of two is what
we wanted, and this contig should consist of the forward and reverse
reads of a getn clone. For example, for a forward designator of "F"
and a reverse designator of "R", if the two reads chad1F and chad1R
were in a single contig (for example Contig 5) it will be determined
that the consensus sequence for Contig 5 will be the sequence for
clone chad1.

Doublets are good!

This module parses C<.ace> and related files. A detailed list of methods
can be found at the end of this document.

I wrote a detailed rationale for design that may explain the reasons
why some things were done the way they were done. That document is
beyond the scope of this pod and can probably be found in the
directory from which this module was 'made' or at
L<http://www.dieselwurks.com/bioinformatics/consedpm_documentation.pdf>.

Note that the POD in that document might be old but the original
rationale still stands.

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
the bugs and their resolution.  Bug reports can be submitted via the
web:

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR - Chad Matsalla

Email chad-at-dieselwurks.com

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#' 

package Bio::Tools::Alignment::Consed;

use strict;

use FileHandle;
use Dumpvalue qw(dumpValue);
use Bio::Tools::Alignment::Trim;
use File::Spec;

use base qw(Bio::Root::Root Bio::Root::IO);

our %DEFAULTS = ( 'f_designator' => 'f',
		  'r_designator' => 'r');

=head2 new()

 Title   : new(-acefile => $path_to_some_acefile, -verbose => "1")
 Usage   : $o_consed = Bio::Tools::Alignment::Consed->
              new(-acefile => $path_to_some_acefile, -verbose => "1");
 Function: Construct the Bio::Tools::Alignment::Consed object. Sets
	   verbosity for the following procedures, if necessary:
	   1. Construct a new Bio::Tools::Alignment::Trim object, to
	   handle quality trimming 2. Read in the acefile and parse it

 Returns : A reference to a Bio::Tools::Alignment::Consed object.
 Args    : A hash. (-acefile) is the filename of an acefile. If a full path
	   is not specified "./" is prepended to the filename and used from
	   instantiation until destruction. If you want 
           Bio::Tools::Alignment::Consed to be noisy during parsing of
           the acefile, specify some value for (-verbose).

=cut

sub new {
    my ($class,%args) = @_;
    my $self = $class->SUPER::new(%args);

    $self->{'filename'} = $args{'-acefile'};

    # this is special to UNIX and should probably use catfile : DONE!
#    if (!($self->{'filename'} =~ m{/})) { 
#	$self->{'filename'} = "./".$self->{'filename'}; 
#    }
#    $self->{'filename'} =~ s#\\#\/#g if $^O =~ m/mswin/i;
#    $self->{'filename'} =~ m/(.*\/)(.*)ace.*$/;
#    $self->{'path'} = $1;

    # this is more generic and should work on most systems   
    (undef, $self->{'path'}, undef) = File::Spec->splitpath($self->{'filename'});

    $self->_initialize_io('-file'=>$self->{'filename'});
    $self->{'o_trim'} = Bio::Tools::Alignment::Trim->new(-verbose => $self->verbose());
    $self->set_forward_designator($DEFAULTS{'f_designator'});
    $self->set_reverse_designator($DEFAULTS{'r_designator'});

    $self->_read_file();
    return $self;
}

=head2 set_verbose()

 Title   : set_verbose()
 Usage   : $o_consed->set_verbose(1);
 Function: Set the verbosity level for debugging messages. On instantiation
	   of the Bio::Tools::Alignment::Consed object the verbosity level
           is set to 0 (quiet).
 Returns : 1 or 0.
 Args    : The verbosity levels are:
	      0 - quiet
	      1 - noisy
	      2 - noisier
	      3 - annoyingly noisy

This method for setting verbosity has largely been superseded by a
sub-by-sub way, where for every sub you can provide a (-verbose)
switch. I am doing converting this bit-by-bit so do not be surprised
if some subs do not honour this.

=cut

# from RootI

# backwards compat
sub set_verbose { (shift)->verbose(@_) }

=head2 get_filename()

 Title   : get_filename()
 Usage   : $o_consed->get_filename();
 Function: Returns the name of the acefile being used by the
	   Bio::Tools::Alignment::Consed object.
 Returns : A scalar containing the name of a file.
 Args    : None.

=cut


sub get_filename {
    my $self = shift;
    return $self->{'filename'};
}

=head2 count_sequences_with_grep()

 Title   : count_sequences_with_grep()
 Usage   : $o_consed->count_sequences_with_grep();
 Function: Use /bin/grep to scan through the files in the ace project dir
	   and count sequences in those files. I used this method in the
	   development of this module to verify that I was getting all of the
	   sequences. It works, but it is (I think) unix-like platform
	   dependent.
 Returns : A scalar containing the number of sequences in the ace project
	   directory.
 Args    : None.

If you are on a non-UNIX platform, you really do not have to use
this. It is more of a debugging routine designed to address very
specific problems.

This method was reimplemented to be platform independent with a pure
perl implementation.  The above note can be ignored.

=cut

sub count_sequences_with_grep {
    my $self = shift;
    my ($working_dir,$grep_cli,@total_grep_sequences);
    # this should be migrated to a pure perl implementation ala
    # Tom Christiansen's 'tcgrep'
    # http://www.cpan.org/modules/by-authors/id/TOMC/scripts/tcgrep.gz

    open my $FILE, '<', $self->{'filename'} or do {
        $self->warn("Could not read file '$self->{'filename'}' for grepping: $!");
        return
    };
    my $counter =0;
    while(<$FILE>) { $counter++ if(/^AF/); }
    close $FILE;

    opendir my $SINGLETS, $self->{'path'};
    foreach my $f ( readdir($SINGLETS) ) {
        next unless ($f =~ /\.singlets$/);

        my $singlet_file = File::Spec->catfile($self->{'path'}, $f);
        open my $S_FILE, '<', $singlet_file or do {
            $self->warn("Could not read file '$singlet_file': $!");
            next
        };
        while(<$S_FILE>) { $counter++ if(/^>/) }
        close $S_FILE;
    }
    closedir $SINGLETS;
    return $counter;
}

=head2 get_path()

 Title   : get_path()
 Usage   : $o_consed->get_path();
 Function: Returns the path to the acefile this object is working with.
 Returns : Scalar. The path to the working acefile.
 Args    : None.

=cut

sub get_path {
    my $self = shift;
    return $self->{'path'};
}

=head2 get_contigs()

 Title   : get_contigs()
 Usage   : $o_consed->get_contigs();
 Function: Return the keys to the Bio::Tools::Alignment::Consed object.
 Returns : An array containing the keynames in the
           Bio::Tools::Alignment::Consed object.
 Args    : None.

This would normally be used to get the keynames for some sort of
iterator. These keys are worthless in general day-to-day use because
in the Consed acefile they are simply Contig1, Contig2, ...

=cut

sub get_contigs {
    my ($self,$contig) = @_;
    my @contigs = sort keys %{$self->{'contigs'}};
    return @contigs;
}

=head2 get_class($contig_keyname)

 Title   : get_class($contig_keyname)
 Usage   : $o_consed->get_class($contig_keyname);
 Function: Return the class name for this contig
 Returns : A scalar representing the class of this contig.
 Args    : None.
 Notes   : 

=cut

sub get_class {
    my ($self,$contig) = @_;
    return $self->{contigs}->{$contig}->{class};
}

=head2 get_quality_array($contig_keyname)

 Title   : get_quality_array($contig_keyname)
 Usage   : $o_consed->get_quality_array($contig_keyname);
 Function: Returns the quality for the consensus sequence for the given
	   contig as an array. See get_quality_scalar to get this as a scalar.
 Returns : An array containing the quality for the consensus sequence with
	   the given keyname.
 Args    : The keyname of a contig. Note: This is a keyname. The key would
	   normally come from get_contigs.

Returns an array, not a reference. Is this a bug? I<thinking> No.
Well, maybe.  Why was this developed like this? I was using FreezeThaw
for object persistence, and when it froze out these arrays it took a
long time to thaw it. Much better as a scalar.

See L<get_quality_scalar()|get_quality_scalar>

=cut

sub get_quality_array {
    my ($self,$contig) = @_;
    return split ' ', $self->{contigs}->{$contig}->{quality};
}

=head2 get_quality_scalar($contig_keyname)

 Title   : get_quality_scalar($contig_keyname)
 Usage   : $o_consed->get_quality_scalar($contig_keyname);
 Function: Returns the quality for the consensus sequence for the given
	   contig as a scalar. See get_quality_array to get this as an array.
 Returns : An scalar containing the quality for the consensus sequence with
           the given keyname.
 Args    : The keyname of a contig. Note this is a _keyname_. The key would
	   normally come from get_contigs.

Why was this developed like this? I was using FreezeThaw for object
persistence, and when it froze out these arrays it took a coon's age
to thaw it. Much better as a scalar.

See L<get_quality_array()|get_quality_array>

=cut

#'
sub get_quality_scalar {
    my ($self,$contig) = @_;
    return $self->{'contigs'}->{$contig}->{'quality'};
}

=head2 freeze_hash()

 Title   : freeze_hash()
 Usage   : $o_consed->freeze_hash();

 Function: Use Ilya's FreezeThaw module to create a persistent data
	   object for this Bio::Tools::Alignment::Consed data
	   structure. In the case of AAFC, we use
	   Bio::Tools::Alignment::Consed to pre-process bunches of
	   sequences, freeze the structures, and send in a harvesting
	   robot later to do database stuff.
 Returns : 0 or 1;
 Args    : None.

This procedure was removed so Consed.pm won't require FreezeThaw.

=cut

#'
sub freeze_hash {
    my $self = shift;
    $self->warn("This method (freeze_hash) was removed ".
                "from the bioperl consed.pm. Sorry.\n");
    if (1==2) {
        $self->debug("Bio::Tools::Alignment::Consed::freeze_hash:".
                     " \$self->{path} is $self->{path}\n");
        my $filename = $self->{'path'}."frozen";
        my %contigs = %{$self->{'contigs'}};
        my $frozen = freeze(%contigs);
        umask 0001;
        open my $FREEZE, '>', $filename or do {
            $self->warn( "Bio::Tools::Alignment::Consed could not ".
                         "freeze the contig hash because the file ".
                         "($filename) could not be opened: $!");
            return 1;
        };
        print $FREEZE $frozen;
        return 0;
    }
}

=head2 get_members($contig_keyname)

 Title   : get_members($contig_keyname)
 Usage   : $o_consed->get_members($contig_keyname);
 Function: Return the _names_ of the reads in this contig.
 Returns : An array containing the names of the reads in this contig.
 Args    : The keyname of a contig. Note this is a keyname. The keyname
	   would normally come from get_contigs.

See L<get_contigs()|get_contigs>

=cut

sub get_members {
    my ($self,$contig) = @_;
    if (!$contig) {
	$self->warn("You need to provide the name of a contig to ".
                    "use Bio::Tools::Alignment::Consed::get_members!\n");
	return;
    }
    return @{$self->{'contigs'}->{$contig}->{'member_array'}};
}

=head2 get_members_by_name($some_arbitrary_name)

 Title   : get_members_by_name($some_arbitrary_name)
 Usage   : $o_consed->get_members_by_name($some_arbitrary_name);
 Function: Return the names of the reads in a contig. This is the name given
	   to $contig{key} based on what is in the contig. This is different
	   from the keys retrieved through get_contigs().
 Returns : An array containing the names of the reads in the contig with this
	   name.
 Args    : The name of a contig. Not a key, but a name.

Highly inefficient. use some other method if possible.
See L<get_contigs()|get_contigs>

=cut

sub get_members_by_name {
    my ($self,$name) = @_;
    # build a list to try to screen for redundancy
    my @contigs_with_that_name;
    foreach my $currkey ( sort keys %{$self->{'contigs'}} ) {
	next if (!$self->{'contigs'}->{$currkey}->{'name'});
	if ($self->{'contigs'}->{$currkey}->{'name'} eq "$name") {
	    push @contigs_with_that_name,$currkey;
	}
    }
    my $count = @contigs_with_that_name;
    if ($count == 1) {
	my $contig_num = $contigs_with_that_name[0];
	return @{$self->{'contigs'}->{$contig_num}->{'member_array'}};
    }
}

=head2 get_contig_number_by_name($some_arbitrary_name)

 Title   : get_contig_number_by_name($some_arbitrary_name)
 Usage   : $o_consed->get_contig_number_by_name($some_arbitrary_name);
 Function: Return the names of the reads in a contig. This is the name given
	   to $contig{key} based on what is in the contig. This is different
	   from the keys retrieved through get_contigs().
 Returns : An array containing the names of the reads in the contig with this
	   name.
 Args    : The name of a contig. Not a key, but a name.

See L<get_contigs()|get_contigs>

=cut

sub get_contig_number_by_name {
    my ($self,$name) = @_;
    foreach my $currkey (sort keys %{$self->{'contigs'}}) {
	if ($self->{'contigs'}->{$currkey}->{'name'} && 
	    $self->{'contigs'}->{$currkey}->{'name'} eq "$name") {
	    return $currkey;
	}
    }
}	

=head2 get_sequence($contig_keyname)

 Title   : get_sequence($contig_keyname)
 Usage   : $o_consed->get_sequence($contig_keyname); 
 Function: Returns the consensus sequence for a given contig.
 Returns : A scalar containing a sequence.
 Args    : The keyname of a contig. Note this is a key. The key would
	   normally come from get_contigs.

See L<get_contigs()|get_contigs>

=cut

sub get_sequence {
    my ($self,$contig) = @_;
    return $self->{'contigs'}->{$contig}->{'consensus'};
}

=head2 set_final_sequence($some_sequence)

 Title   : set_final_sequence($name,$some_sequence)
 Usage   : $o_consed->set_final_sequence($name,$some_sequence);
 Function: Provides a manual way to set the sequence for a given key in the
	   contig hash. Rarely used.
 Returns : 0 or 1;
 Args    : The name (not the keyname) of a contig and an arbitrary string.

A method with a questionable and somewhat mysterious origin. May raise
the dead or something like that.

=cut

sub set_final_sequence {
    my ($self,$name,$sequence) = @_;
    if (!$self->{'contigs'}->{$name}) {
	$self->warn("You cannot set the final sequence for ".
                    "$name because it doesn't exist!\n");
	return 1;
    }
    else {
	$self->{'contigs'}->{$name}->{'final_sequence'} = $sequence;
    }
    return 0;
}

=head2  _read_file()

 Title   : _read_file();
 Usage   : _read_file();
 Function: An internal subroutine used to read in an acefile and parse it
	   into a Bio::Tools::Alignment::Consed object.
 Returns : 0 or 1.
 Args    : Nothing.

This routine creates and saves the filhandle for reading the files in
{fh}

=cut

sub _read_file {
    my ($self) = @_;
    my ($line,$in_contig,$in_quality,$contig_number,$top);
    # make it easier to type $fhl
    while (defined($line=$self->_readline()) ) {
	chomp $line;
	# check if there is anything on this line
	# if not, you can stop gathering consensus sequence
	if (!$line) {
	    # if the line is blank you are no longer to gather consensus 
	    # sequence or quality values
	    $in_contig = 0;
	    $in_quality = 0;
	}
	# you are currently gathering consensus sequence
	elsif ($in_contig) {
	    if ($in_contig == 1) {
		$self->debug("Adding $line to consensus of contig number $contig_number.\n");
		$self->{'contigs'}->{$contig_number}->{'consensus'} .= $line;
	    }
	}
	elsif ($in_quality) {
	    if (!$line) {
		$in_quality = undef;
	    }
	    else {

		# I wrote this in here because acefiles produced by
		# cap3 do not have a leading space like the acefiles
		# produced by phrap and there is the potential to have
		# concatenated quality values like this: 2020 rather
		# then 20 20 whre lines collide. Thanks Andrew for
		# noticing.

		if ($self->{'contigs'}->{$contig_number}->{'quality'} &&
                    !($self->{'contigs'}->{$contig_number}->{'quality'} =~ m/\ $/)) {
		    $self->{'contigs'}->{$contig_number}->{'quality'} .= " ";
		}
		$self->{'contigs'}->{$contig_number}->{'quality'} .= $line;
	    }
	}
	elsif ($line =~ /^BQ/) {
	    $in_quality = 1;
	}

	# the line /^CO/ like this:
	# CO Contig1 796 1 1 U
	# can be broken down as follows:
	# CO - Contig!
	# Contig1 - the name of this contig
	# 796 - Number of bases in this contig
	# 1 - Number of reads in this contig
	# 1 - number of base segments in this contig
	# U - Uncomplemented

	elsif ($line =~ /^CO/) {
	    $line =~ m/^CO\ Contig(\d+)\ \d+\ \d+\ \d+\ (\w)/;
	    $contig_number = $1;
	    if ($2 eq "C") {
		$self->debug("Contig $contig_number is complemented!\n");
	    }
	    $self->{'contigs'}->{$contig_number}->{'member_array'} = [];
	    $self->{'contigs'}->{$contig_number}->{'contig_direction'} = "$2";
	    $in_contig = 1;
	}

	# 000713
	# this BS is deprecated, I think.
	# haha, I am really witty. <ew>

	elsif ($line =~ /^BSDEPRECATED/) {
	    $line =~ m/^BS\s+\d+\s+\d+\s+(.+)/;
	    my $member = $1;
	    $self->{'contigs'}->{$contig_number}->{$member}++;
	}
	# the members of the contigs are determined by the AF line in the ace file
	elsif ($line =~ /^AF/) {
	    $self->debug("I see an AF line here.\n");
	    $line =~ /^AF\ (\S+)\ (\w)\ (\S+)/;

            # push the name of the current read onto the member array for this contig
	    push @{$self->{'contigs'}->{$contig_number}->{'member_array'}},$1;

            # the first read in the contig will be named the "top" read
	    if (!$top) {
		$self->debug("\$top is not set.\n");
		if ($self->{'contigs'}->{$contig_number}->{'contig_direction'} eq "C") {
		    $self->debug("Reversing the order of the reads. The bottom will be $1\n");

		    # if the contig sequence is marked as the
		    # complement, the top becomes the bottom and$
		    $self->{'contigs'}->{$contig_number}->{'bottom_name'} = $1;
		    $self->{'contigs'}->{$contig_number}->{'bottom_complement'} = $2;
		    $self->{'contigs'}->{$contig_number}->{'bottom_start'} = $3;
		}
		else {
		    $self->debug("NOT reversing the order of the reads. ".
                                 "The top_name will be $1\n");
		    # if the contig sequence is marked as the
		    # complement, the top becomes the bottom and$
		    $self->{'contigs'}->{$contig_number}->{'top_name'} = $1;
		    $self->{'contigs'}->{$contig_number}->{'top_complement'} = $2;
		    $self->{'contigs'}->{$contig_number}->{'top_start'} = $3;
		}
		$top = 1;
	    }
	    else {

		# if the contig sequence is marked as the complement,
		# the top becomes the bottom and the bottom becomes
		# the top
		if ($self->{'contigs'}->{$contig_number}->{'contig_direction'} eq "C") {
		    $self->debug("Reversing the order of the reads. The top will be $1\n");
		    $self->{'contigs'}->{$contig_number}->{'top_name'} = $1;
		    $self->{'contigs'}->{$contig_number}->{'top_complement'} = $2;
		    $self->{'contigs'}->{$contig_number}->{'top_start'} = $3;
		}
		else {
		    $self->debug("NOT reversing the order of the reads. The bottom will be $1\n");
		    $self->{'contigs'}->{$contig_number}->{'bottom_name'} = $1;
		    $self->{'contigs'}->{$contig_number}->{'bottom_complement'} = $2;
		    $self->{'contigs'}->{$contig_number}->{'bottom_start'} = $3;
		}
		$top = undef;
	    }
	}
    }
    return 0;
}

=head2 set_reverse_designator($some_string)

 Title   : set_reverse_designator($some_string)
 Usage   : $o_consed->set_reverse_designator($some_string);
 Function: Set the designator for the reverse read of contigs in this
	   Bio::Tools::Alignment::Consed object. Used to determine if
           contigs containing two reads can be named.
 Returns : The value of $o_consed->{reverse_designator} so you can check
	   to see that it was set properly.
 Args    : An arbitrary string.

May be useful only to me. I<shrug>

=cut

sub set_reverse_designator {
    my ($self,$reverse_designator) = @_;
    $self->{'reverse_designator'} = $reverse_designator;
    $self->{'o_trim'}->set_reverse_designator($reverse_designator);
    return $self->{'reverse_designator'};
}				# end set_reverse_designator

=head2 set_forward_designator($some_string)

 Title   : set_forward_designator($some_string)
 Usage   : $o_consed->set_forward_designator($some_string);
 Function: Set the designator for the forward read of contigs in this
	   Bio::Tools::Alignment::Consed object. Used to determine if
           contigs containing two reads can be named.
 Returns : The value of $o_consed->{forward_designator} so you can check
	   to see that it was set properly.
 Args    : An arbitrary string.

May be useful only to me. I<shrug>

=cut

sub set_forward_designator {
    my ($self,$forward_designator) = @_;
    $self->{'forward_designator'} = $forward_designator;
    $self->{'o_trim'}->set_forward_designator($forward_designator);
    return $self->{'forward_designator'};
}				# end set_forward_designator

=head2 set_designator_ignore_case("yes")

 Title   : set_designator_ignore_case("yes")
 Usage   : $o_consed->set_designator_ignore_case("yes");
 Function: Deprecated.
 Returns : Deprecated.
 Args    : Deprecated.

Deprecated. Really. Trust me.

=cut

sub set_designator_ignore_case {
    my ($self,$ignore_case) = @_;
    if ($ignore_case eq "yes") {
	$self->{'designator_ignore_case'} = 1;
    }
    return $self->{'designator_ignore_case'};
}				# end set_designator_ignore_case

=head2 set_trim_points_singlets_and_singletons()

 Title   : set_trim_points_singlets_and_singletons()
 Usage   : $o_consed->set_trim_points_singlets_and_singletons();
 Function: Set the trim points for singlets and singletons based on
	   quality.  Uses the Bio::Tools::Alignment::Trim object. Use
	   at your own risk because the Bio::Tools::Alignment::Trim
	   object was designed specifically for me and is mysterious
	   in its ways. Every time somebody other then me uses it a
	   swarm of locusts decends on a small Central American
	   village so do not say you weren't warned.
 Returns : Nothing.
 Args    : None.

Working on exceptions and warnings here.

See L<Bio::Tools::Alignment::Trim> for more information

=cut

#' to make my emacs happy

sub set_trim_points_singlets_and_singletons {
    my ($self) = @_;
    $self->debug("Consed.pm : \$self is $self\n");
    my (@points,$trimmed_sequence);
    if (!$self->{'doublets_set'}) {
        $self->debug("You need to set the doublets before you use ".
                     "set_trim_points_singlets_and_doublets. Doing that now.");
	$self->set_doublets();
    }
    foreach (sort keys %{$self->{'contigs'}}) {
	if ($self->{'contigs'}->{$_}->{'class'} eq "singlet") {
	    $self->debug("Singlet $_\n");
	    # this is what Warehouse wants
	    #         my ($self,$sequence,$quality,$name) = @_;
	    # this is what Bio::Tools::Alignment::Trim::trim_singlet wants:
	    # my ($self,$sequence,$quality,$name,$class) = @_;
	    # the following several lines are to make the parameter passing legible.
	    my ($sequence,$quality,$name,$class);
	    $sequence = $self->{'contigs'}->{$_}->{'consensus'};
	    if (!$self->{'contigs'}->{$_}->{'quality'}) { $quality = "unset"; }
	    else { $quality = $self->{'contigs'}->{$_}->{'quality'}; }
	    $name = $self->{'contigs'}->{$_}->{'name'};
	    $class = $self->{'contigs'}->{$_}->{'class'};
	    @points = @{$self->{'o_trim'}->trim_singlet($sequence,$quality,$name,$class)};
	    $self->{'contigs'}->{$_}->{'start_point'} = $points[0];
	    $self->{'contigs'}->{$_}->{'end_point'} = $points[1];
	    $self->{'contigs'}->{$_}->{'sequence_trimmed'} = 
                substr($self->{contigs}->{$_}->{'consensus'},$points[0],$points[1]-$points[0]);
	}
    }
    $self->debug("Bio::Tools::Alignment::Consed::set_trim_points_singlets".
                 "_and_singletons: Done setting the quality trimpoints.\n");
    return;
}  # end set_trim_points_singlet

=head2 set_trim_points_doublets()

 Title   : set_trim_points_doublets()
 Usage   : $o_consed->set_trim_points_doublets();
 Function: Set the trim points for doublets based on quality. Uses the
	   Bio::Tools::Alignment::Trim object. Use at your own risk because
           the Bio::Tools::Alignment::Trim object was designed specifically
           for me and is mysterious in its ways. Every time somebody other
           then me uses it you risk a biblical plague being loosed on your
           city.
 Returns : Nothing.
 Args    : None.
 Notes   : Working on exceptions here.

See L<Bio::Tools::Alignment::Trim> for more information

=cut

sub set_trim_points_doublets {
    my $self = shift;
    my @points;
    $self->debug("Bio::Tools::Alignment::Consed::set_trim_points_doublets:".
                 " Restoring zeros for doublets.\n");
    # &show_missing_sequence($self);
    $self->debug("Bio::Tools::Alignment::Consed::set_trim_points_doublets:".
                 " Setting doublet trim points.\n");
    foreach (sort keys %{$self->{'contigs'}}) {
	if ($self->{'contigs'}->{$_}->{'class'} eq "doublet") {
            # my ($self,$sequence,$quality,$name,$class) = @_;
            my @quals = split(' ',$self->{'contigs'}->{$_}->{'quality'});

	    @points = $self->{o_trim}->trim_doublet
                ($self->{'contigs'}->{$_}->{'consensus'},
                 $self->{'contigs'}->{$_}->{'quality'},
                 $self->{'contigs'}->{$_}->{name},
                 $self->{'contigs'}->{$_}->{'class'});
	    $self->{'contigs'}->{$_}->{'start_point'} = $points[0];
	    $self->{'contigs'}->{$_}->{'end_point'} = $points[1];
            # now set this
	    $self->{'contigs'}->{$_}->{'sequence_trimmed'} =
                substr($self->{contigs}->{$_}->{'consensus'},
                       $points[0],$points[1]-$points[0]);
	    # 010102 the deprecated way to do things:
	}
    }
    $self->debug("Bio::Tools::Alignment::Consed::set_trim_points_doublets:".
                 " Done setting doublet trim points.\n"); 
    return;
} # end set_trim_points_doublets

=head2 get_trimmed_sequence_by_name($name)

 Title   : get_trimmed_sequence_by_name($name)
 Usage   : $o_consed->get_trimmed_sequence_by_name($name);
 Function: Returns the trimmed_sequence of a contig with {name} eq $name.
 Returns : A scalar- the trimmed sequence.
 Args    : The {name} of a contig.
 Notes   : 

=cut

sub get_trimmed_sequence_by_name {
    my ($self,$name) = @_;
    my $trimmed_sequence;
    my $contigname = &get_contig_number_by_name($self,$name);
    my $class = $self->{'contigs'}->{$contigname}->{'class'};
    # what is this business and who was smoking crack while writing this?
    # if ($class eq "singlet") {
    # send the sequence, the quality, and the name
    # $trimmed_sequence = $self->{o_trim}->trim_singlet
    #  ($self->{'contigs'}->{$contigname}->{consensus},
    #   $self->{'contigs'}->{$contigname}->{'quality'},$name);
    # }
    return $self->{'contigs'}->{$contigname}->{'sequence_trimmed'};
}

=head2 set_dash_present_in_sequence_name("yes")

 Title   : set_dash_present_in_sequence_name("yes")
 Usage   : $o_consed->set_dash_present_in_sequence_name("yes");
 Function: Deprecated. Part of an uncompleted thought. ("Oooh! Shiny!")
 Returns : Nothing.
 Args    : "yes" to set {dash_present_in_sequence_name} to 1
 Notes   : 

=cut

sub set_dash_present_in_sequence_name {
    my ($self,$dash_present) = @_;
    if ($dash_present eq "yes") {
	$self->{'dash_present_in_sequence_name'} = 1;
    }
    else {
	$self->{'dash_present_in_sequence_name'} = 0;
    }
    return $self->{'dash_present_in_sequence_name'};
} # end set_dash_present_in_sequence_name

=head2 set_doublets()

 Title   : set_doublets()
 Usage   : $o_consed->set_doublets();
 Function: Find pairs that have similar names and mark them as doublets
	   and set the {name}.
 Returns : 0 or 1.
 Args    : None.

A complicated subroutine that iterates over the
Bio::Tools::Alignment::Consed looking for contigs of 2. If the forward
and reverse designator are removed from each of the reads in
{'member_array'} and the remaining reads are the same, {name} is set
to that name and the contig's class is set as "doublet".  If any of
those cases fail the contig is marked as a "pair".

=cut

#' make my emacs happy

sub set_doublets {
    my ($self) = @_;
    # set the designators in the Bio::Tools::Alignment::Trim object

    $self->{'o_trim'}->set_designators($self->{'reverse_designator'},
				       $self->{'forward_designator'});
    foreach my $key_contig (sort keys %{$self->{'contigs'}}) {

	# if there is a member array (why would there not be? This should be a die()able offence
	# but for now I will leave it
	if ($self->{'contigs'}->{$key_contig}->{'member_array'}) {
	    # if there are two reads in this contig 
	    # i am pretty sure that this is wrong but i am keeping it for reference
	    # if (@{$self->{'contigs'}->{$key_contig}->{'member_array'}} == 2 || !$self->{'contigs'}->{$key_contig}->{'class'}) {
	    # <seconds later>
	    # <nod> WRONG. Was I on crack?
	    if (@{$self->{'contigs'}->{$key_contig}->{'member_array'}} == 2) {
		$self->{'contigs'}->{$key_contig}->{'num_members'} = 2;
		$self->debug("\tThere are 2 members! Looking for the contig name...\n");
		my $name = _get_contig_name($self,$self->{'contigs'}->{$key_contig}->{'member_array'});
		$self->debug("The name is $name\n") if defined $name;
		if ($name) {
		    $self->{'contigs'}->{$key_contig}->{'name'} = $name;
		    $self->{'contigs'}->{$key_contig}->{'class'} = "doublet";
		} else {
		    $self->debug("$key_contig is a pair.\n");
		    $self->{'contigs'}->{$key_contig}->{'class'} = "pair";
		}
	    }
            # this is all fair and good but what about singlets?
            # they have one reads in the member_array but certainly are not singletons
	    elsif (@{$self->{'contigs'}->{$key_contig}->{'member_array'}} == 1) {
		# set the name to be the name of the read
		$self->{'contigs'}->{$key_contig}->{name} = @{$self->{'contigs'}->{$key_contig}->{'member_array'}}[0];
		# set the number of members to be one
		$self->{'contigs'}->{$key_contig}->{num_members} = 1;
		# if this was a singlet, it would already belong to the class "singlet"
		# so leave it alone
		# if it is not a singlet, it is a singleton! lablel it appropriately
		unless ($self->{'contigs'}->{$key_contig}->{'class'}) {
		    $self->{'contigs'}->{$key_contig}->{'class'} = "singleton";
		}
	    }
            # set the multiplet characteristics
	    elsif (@{$self->{'contigs'}->{$key_contig}->{'member_array'}} >= 3) {
		$self->{'contigs'}->{$key_contig}->{'num_members'} = @{$self->{'contigs'}->{$key_contig}->{'member_array'}};
		$self->{'contigs'}->{$key_contig}->{'class'} = "multiplet";
	    }
	    $self->{'contigs'}->{$key_contig}->{'num_members'} = @{$self->{'contigs'}->{$key_contig}->{'member_array'}};

	}
    }
    $self->{'doublets_set'} = "done";
    return 0;
}				# end set_doublets

=head2 set_singlets

 Title   : set_singlets
 Usage   : $o_consed->set_singlets();
 Function: Read in a singlets file and place them into the
	   Bio::Tools::Alignment::Consed object.
 Returns : Nothing.
 Args    : A scalar to turn on verbose parsing of the singlets file.
 Notes   : 

=cut

sub set_singlets {
    # parse out the contents of the singlets file
    my ($self) = @_;
    $self->debug("Bio::Tools::Alignment::Consed Adding singlets to the contig hash...\n"); 
    my $full_filename = $self->{'filename'};
    $self->debug("Bio::Tools::Alignment::Consed::set_singlets: \$full_filename is $full_filename\n");
    $full_filename =~ s#\\#\/#g if $^O =~ m/mswin/i;
    $full_filename =~ m/(.*\/)(.*ace.*)$/; 			       
    my ($base_path,$filename) = ($1,$2);
    $self->debug("Bio::Tools::Alignment::Consed::set_singlets: singlets filename is $filename and \$base_path is $base_path\n");
    $filename =~ m/(.*)ace.*$/;
    my $singletsfile = $base_path.$1."singlets";
    $self->debug("\$singletsfile is $singletsfile\n");
     if (!-f $singletsfile) {
          # there is no singlets file.
          $self->{'singlets_set'} = "done";
          return;
     }
	$self->debug("$singletsfile is indeed a file. Trying to open it...\n");
    my $singlets_fh = Bio::Root::IO->new(-file => $singletsfile);
    my ($sequence,$name,$count);
    while ($_ = $singlets_fh->_readline()) {
	chomp $_;
	if (/\>/) {
	    if ($name && $sequence) {
		$self->debug("Adding $name with sequence $sequence to hash...\n");
		push @{$self->{'contigs'}->{$name}->{'member_array'}},$name;
		$self->{'contigs'}->{$name}->{'consensus'} = $sequence;
		$self->{'contigs'}->{$name}->{'name'} = $name;
		$self->{'contigs'}->{$name}->{"singlet"} = 1;
		$self->{'contigs'}->{$name}->{'class'} = "singlet";
	    }
	    $sequence = $name = undef;
	    $count++;
	    m/^\>(.*)\s\sCHROMAT/;
	    $name = $1;
	    if (!$name) {
		m/\>(\S+)\s/;
		$name = $1;
	    }
	}
	else { $sequence .= $_; }	
    }
    if ($name && $sequence) {
	$self->debug("Pushing the last of the singlets ($name)\n");
	@{$self->{'contigs'}->{$name}->{'member_array'}} = $name;
	$self->{'contigs'}->{$name}->{'consensus'} = $sequence;
	$self->{'contigs'}->{$name}->{'name'} = $name;
	$self->{'contigs'}->{$name}->{"singlet"} = 1;
	$self->{'contigs'}->{$name}->{'class'} = "singlet";
    }
    $self->debug("Bio::Tools::Alignment::Consed::set_singlets: Done adding singlets to the singlets hash.\n");
    $self->{'singlets_set'} = "done";
    return 0;
}				# end sub set_singlets

=head2 get_singlets()

 Title   : get_singlets()
 Usage   : $o_consed->get_singlets();
 Function: Return the keynames of the singlets.
 Returns : An array containing the keynames of all 
           Bio::Tools::Alignment::Consed sequences in the class "singlet".
 Args    : None.
 Notes   : 

=cut

sub get_singlets {
    # returns an array of singlet names
    # singlets have "singlet"=1 in the hash
    my $self = shift;
    if (!$self->{singlets_set}) {
	$self->debug("You need to set the singlets before you get them. Doing that now.");
	$self->set_singlets();
    }	

    my (@singlets,@array);
    foreach my $key (sort keys %{$self->{'contigs'}}) {
	# @array = @{$Consed::contigs{$key}->{'member_array'}};
	# somethimes a user will try to get a list of singlets before the classes for the rest of the
	# contigs has been set (see t/test.t for how I figured this out. <bah>			
	# so either way, just return class=singlets
	if (!$self->{'contigs'}->{$key}->{'class'}) {
	    # print("$key has no class. why?\n");
	}
	elsif ($self->{'contigs'}->{$key}->{'class'} eq "singlet") {
	    push @singlets,$key;
	}
    }
    return @singlets;
}

=head2 set_quality_by_name($name,$quality)

 Title   : set_quality_by_name($name,$quality)
 Usage   : $o_consed->set_quality_by_name($name,$quality);
 Function: Deprecated. Make the contig with {name} have {'quality'} $quality.
           Probably used for testing.
 Returns : Nothing.
 Args    : The name of a contig and a scalar for its quality.
 Notes   : Deprecated.

=cut

sub set_quality_by_name {
    # this is likely deprecated
    my ($self,$name,$quality) = shift;
    my $return;
    foreach (sort keys %{$self->{'contigs'}}) {
	if ($self->{'contigs'} eq "$name" || $self->{'contigs'}->{'name'} eq "$name") {
	    $self->{'contigs'}->{'quality'} = $quality;
	    $return=1;
	}
    }
    if ($return) { return "0"; } else { return "1"; }
}				# end set quality by name

=head2 set_singlet_quality()

 Title   : set_singlet_quality()
 Usage   : $o_consed->set_singlet_quality();
 Function: For each singlet, go to the appropriate file in phd_dir and read
           in the phred quality for that read and place it into {'quality'}
 Returns : 0 or 1.
 Args    : None.
 Notes   : This is the next subroutine that will receive substantial revision
           in the next little while. It really should eval the creation of
           Bio::Tools::Alignment::Phred objects and go from there. 

=cut

sub set_singlet_quality {
    my $self = shift;
    my $full_filename = $self->{'filename'};
    $full_filename =~ s#\\#\/#g if $^O =~ m/mswin/i;
    $full_filename =~ m/(.*\/)(.*)ace.*$/;
    my ($base_path,$filename) = ($1,"$2"."qual");
    my $singletsfile = $base_path.$filename;
    if (-f $singletsfile) {
	# print("$singletsfile is indeed a file. Trying to open it...\n");
    }
    else {
	$self->warn("$singletsfile is not a file. Sorry.\n");
	return;
    }
	my $singlets_fh = Bio::Root::IO->new(-file => $singletsfile);
	my ($sequence,$name,$count);
    my ($identity,$line,$quality,@qline);
    while ($line = $singlets_fh->_readline()) {
	chomp $line;
	if ($line =~ /^\>/) {
	    $quality = undef;
	    $line =~ m/\>(\S*)\s/;
	    $identity = $1;
	}
	else {
	    if ($self->{'contigs'}->{$identity}) {
		$self->{'contigs'}->{$identity}->{'quality'} .= "$line ";
	    }
	}

    }
    return 0;
}

=head2 set_contig_quality()

 Title   : set_contig_quality()
 Usage   : $o_consed->set_contig_quality();
 Function: Deprecated.
 Returns : Deprecated.
 Args    : Deprecated.
 Notes   : Deprecated. Really. Trust me.

=cut

sub set_contig_quality {
    # note: contigs _include_ singletons but _not_ singlets
    my ($self) = shift;
          # the unexpected results I am referring to here are a doubling of quality values.
          # the profanity I uttered on discovering this reminded me of the simpsons:
          # Ned Flanders: "That is the loudest profanity I have ever heard!"
     $self->warn("set_contig_quality is deprecated and will likely produce unexpected results");
    my $full_filename = $self->{'filename'};
    # Run_SRC3700_2000-08-01_73+74.fasta.screen.contigs.qual
    # from Consed.pm
    $full_filename =~ s#\\#\/#g if $^O =~ m/mswin/i;
    $full_filename =~ m/(.*\/)(.*)ace.*$/;
    my ($base_path,$filename) = ($1,"$2"."contigs.qual");
    my $singletsfile = $base_path.$filename;
    if (-f $singletsfile) {
	# print("$singletsfile is indeed a file. Trying to open it...\n");
    }
    else {
	$self->warn("Bio::Tools::Alignment::Consed::set_contig_quality $singletsfile is not a file. Sorry.\n");
	return;
    }
    my $contig_quality_fh = Bio::Root::IO->new(-file => $singletsfile);

    my ($sequence,$name,$count,$identity,$line,$quality);
    while ($line = $contig_quality_fh->_readline()) {
	chomp $line;
	if ($line =~ /^\>/) {
	    $quality = undef;
	    $line =~ m/\>.*Contig(\d+)\s/;
	    $identity = $1;
	}
	else {
	    if ($self->{'contigs'}->{$identity} ) {
		$self->{'contigs'}->{$identity}->{'quality'} .= " $line";
	    }
	}
    }
}				# end set_contig_quality

=head2 get_multiplets()

 Title   : get_multiplets()
 Usage   : $o_consed->get_multiplets();
 Function: Return the keynames of the multiplets.
 Returns : Returns an array containing the keynames of all 
           Bio::Tools::Alignment::Consed sequences in the class "multiplet".
 Args    : None.
 Notes   : 

=cut

sub get_multiplets {
	    # returns an array of multiplet names
	    # multiplets have # members > 2
    my $self = shift;
    my (@multiplets,@array);
    foreach my $key (sort keys %{$self->{'contigs'}}) {
	if ($self->{'contigs'}->{$key}->{'class'}) {
	    if ($self->{'contigs'}->{$key}->{'class'} eq "multiplet") {
		push @multiplets,$key;
	    }
	}
    }
    return @multiplets;
}

=head2 get_all_members()

  Title   : get_all_members()
  Usage   : @all_members = $o_consed->get_all_members();
  Function: Return a list of all of the read names in the 
            Bio::Tools::Alignment::Consed object.
  Returns : An array containing all of the elements in all of the
            {'member_array'}s.
  Args    : None.
  Notes   : 

=cut

sub get_all_members {
    my $self = shift;
    my @members;
    foreach my $key (sort keys %{$self->{'contigs'}}) {
	if ($key =~ /^singlet/) {
	    push @members,$self->{'contigs'}->{$key}->{'member_array'}[0];
	}
	elsif ($self->{'contigs'}->{$key}->{'member_array'}) {
	    push @members,@{$self->{'contigs'}->{$key}->{'member_array'}};
	}
	# else {
	#	print("Bio::Tools::Alignment::Consed: $key is _not_ an array. Pushing $self->{'contigs'}->{$key}->{'member_array'} onto \@members\n");
	#	push @members,$self->{'contigs'}->{$key}->{'member_array'};
	# }
    }
    return @members;
}

=head2 sum_lets($total_only)

 Title   : sum_lets($total_only)
 Usage   : $statistics = $o_consed->sum_lets($total_only);
 Function: Provide numbers for how many sequences were accounted for in the
           Bio::Tools::Alignment::Consed object.
 Returns : If a scalar is present, returns the total number of
           sequences accounted for in all classes. If no scalar passed
           then returns a string that looks like this:
           Singt/singn/doub/pair/mult/total : 2,0,1(2),0(0),0(0),4
           This example means the following: There were 1 singlets.
           There were 0 singletons.  There were 1 doublets for a total
           of 2 sequences in this class.  There were 0 pairs for a
           total of 0 sequences in this class.  There were 0
           multiplets for a total of 0 sequences in this class.  There
           were a total of 4 sequences accounted for in the
           Bio::Tools::Alignment::Consed object.   
 Args : A scalar is optional to change the way the numbers are returned.  
 Notes:

=cut

sub sum_lets {
    my ($self,$total_only) = @_;
    my ($count,$count_multiplets,$multiplet_count);
    my $singlets = &get_singlets($self); $count += $singlets;
    my $doublets = &get_doublets($self); $count += ($doublets * 2);
    my $pairs = &get_pairs($self); $count += ($pairs * 2);
    my $singletons = &get_singletons($self); $count += $singletons;
    my @multiplets = &get_multiplets($self);
    $count_multiplets = @multiplets;
    my $return_string;
    foreach (@multiplets) {
	my $number_members = $self->{'contigs'}->{$_}->{num_members};	
	$multiplet_count += $number_members;
    }
    if ($multiplet_count) {
	$count += $multiplet_count;
    }
    foreach (qw(multiplet_count singlets doublets pairs singletons
                multiplets count_multiplets)) {
	no strict 'refs';	# renege for the block
	if (!${$_}) {
            ${$_} = 0;
	}
    }
    if (!$multiplet_count) { $multiplet_count = 0; }
    if ($total_only) {
        return $count;
    }
    $return_string = "Singt/singn/doub/pair/mult/total : ".
        "$singlets,$singletons,$doublets(".
         ($doublets*2)."),$pairs(".($pairs*2).
        "),$count_multiplets($multiplet_count),$count";
    return $return_string;
}

=head2 write_stats()

 Title   : write_stats()
 Usage   : $o_consed->write_stats();
 Function: Write a file called "statistics" containing numbers similar to
	   those provided in sum_lets().
 Returns : Nothing. Write a file in $o_consed->{path} containing something
	   like this:

           0,0,50(100),0(0),0(0),100

           Where the numbers provided are in the format described in the
	   documentation for sum_lets().
 Args    : None.
 Notes   : This might break platform independence, I do not know.

See L<sum_lets()|sum_lets>

=cut

sub write_stats {
    # worry about platform dependence here?
    # oh shucksdarn.
    my $self = shift;
    my $stats_filename = $self->{'path'}."statistics";
    my $statistics_raw = $self->sum_lets;
    my ($statsfilecontents) = $statistics_raw =~ s/.*\ \:\ //g;
    umask 0001;
    my $fh = Bio::Root::IO->new(-file=>"$stats_filename");
    # open my $STATSFILE, '>', $stats_filename or print "Could not write the statsfile: $!\n");
    $fh->_print("$statsfilecontents");
    # close $STATSFILE;
    $fh->close();
}

=head2 get_singletons()

 Title   : get_singletons()
 Usage   : @singletons = $o_consed->get_singletons();
 Function: Return the keynames of the singletons.
 Returns : Returns an array containing the keynames of all
	   Bio::Tools::Alignment::Consed sequences in the class "singleton".
 Args    : None.
 Notes   : 

=cut

sub get_singletons {
		# returns an array of singleton names
		# singletons are contigs with one member (see consed documentation)
	my $self = shift;
	my (@singletons,@array);
	foreach my $key (sort keys %{$self->{'contigs'}}) {
		if ($self->{'contigs'}->{$key}->{'class'}) {
		    # print ("$key class: $self->{'contigs'}->{$key}->{'class'}\n");
		}
		else {
		    # print("$key belongs to no class. why?\n");
		}
		if ($self->{'contigs'}->{$key}->{'member_array'}) {
			@array = @{$self->{'contigs'}->{$key}->{'member_array'}};
		}
		my $num_array_elem = @array;
		if ($num_array_elem == 1 && $self->{'contigs'}->{$key}->{'class'} && $self->{'contigs'}->{$key}->{'class'} eq "singleton") { push @singletons,$key; }
	}
	return @singletons;
}

=head2 get_pairs()

 Title   : get_pairs()
 Usage   : @pairs = $o_consed->get_pairs();
 Function: Return the keynames of the pairs.
 Returns : Returns an array containing the keynames of all
           Bio::Tools::Alignment::Consed sequences in the class "pair".
 Args    : None.
 Notes   : 

=cut

sub get_pairs {
    # returns an array of pair contig names
    # a pair is a contig of two where the names do not match
    my $self = shift;
    my (@pairs,@array);
    foreach my $key (sort keys %{$self->{'contigs'}}) {
        if ($self->{'contigs'}->{$key}->{'member_array'}) {
            if (@{$self->{'contigs'}->{$key}->{'member_array'}} == 2 &&
                $self->{'contigs'}->{$key}->{'class'} eq "pair") {
                push @pairs,$key;
            }
        }
    }
    return @pairs;
}

=head2 get_name($contig_keyname)

 Title   : get_name($contig_keyname)
 Usage   : $name = $o_consed->get_name($contig_keyname);
 Function: Return the {name} for $contig_keyname.
 Returns : A string. ({name})
 Args    : A contig keyname.
 Notes   : 

=cut

sub get_name {
    my ($self,$contig) = @_;
    return $self->{'contigs'}->{$contig}->{'name'};
}

=head2 _get_contig_name(\@array_containing_reads)

 Title   : _get_contig_name(\@array_containing_reads)
 Usage   : $o_consed->_get_contig_name(\@array_containing_reads);
 Function: The logic for the set_doublets subroutine.
 Returns : The name for this contig.
 Args    : A reference to an array containing read names.
 Notes   : Depends on reverse_designator. Be sure this is set the way you
	   intend.

=cut

sub _get_contig_name {
    my ($self,$r_array) = @_;
    my @contig_members = @$r_array;
    my @name_nodir;
    foreach (@contig_members) {
        # how can I distinguish the clone name from the direction label?
        # look for $Consed::reverse_designator and $Consed::forward_designator
        # what if you do not find _any_ of those?
        my $forward_designator = $self->{'forward_designator'} || "f";
        my $reverse_designator = $self->{'reverse_designator'} || "r";
        my $any_hits = /(.+)($forward_designator.*)/ || /(.+)($reverse_designator.*)/||/(.+)(_.+)/;
        my $name = $1;
        my $suffix = $2;
        if ($name) {
            # print("\t\$name is $name ");
        }
        if ($suffix) {
            # print("and \$suffix is $suffix.\n");
        }
                                # Jee, I hope we get a naming convention soon
        if ($suffix) {
            if ($suffix =~ /^$forward_designator/ || $suffix =~ /^$reverse_designator/) {
                push @name_nodir,$name;
            }
				# bugwatch here! should this be unnested?
            else {
                push @name_nodir,"$name$suffix";
            }
        }
    }
    # print("\@name_nodir: @name_nodir\n");
    my $mismatch = 0;
    for (my $counter=0; $counter<@name_nodir;$counter++) {
        next if ($name_nodir[0] eq $name_nodir[$counter]);
        $mismatch = 1;
    }
    if ($mismatch == 0) {
        # print("\tYou have a cohesive contig named $name_nodir[0].\n\n");
        return $name_nodir[0];
    } else {
        # print("\tYou have mixed names in this contig.\n\n");
    }
}                               # end _get_contig_name

=head2 get_doublets()

 Title   : get_doublets()
 Usage   : @doublets = $o_consed->get_doublets();
 Function: Return the keynames of the doublets.
 Returns : Returns an array containing the keynames of all
           Bio::Tools::Alignment::Consed sequences in the class "doublet".
 Args    : None.
 Notes   : 

=cut

sub get_doublets {
    my $self = shift;
    if (!$self->{doublets_set}) {
        $self->warn("You need to set the doublets before you can get them. Doing that now.");
        $self->set_doublets();
    }
    my @doublets;
    foreach (sort keys %{$self->{'contigs'}}) {
        if ($self->{'contigs'}->{$_}->{name} && $self->{'contigs'}->{$_}->{'class'} eq "doublet") {
            push @doublets,$_;
        }
    }
    return @doublets;
}                               # end get_doublets

=head2 dump_hash()

 Title   : dump_hash()
 Usage   : $o_consed->dump_hash();
 Function: Use dumpvar.pl to dump out the Bio::Tools::Alignment::Consed
           object to STDOUT.
 Returns : Nothing.
 Args    : None.
 Notes   : I used this a lot in debugging.

=cut

sub dump_hash {
    my $self = shift;
    my $dumper = Dumpvalue->new();
    $self->debug( "Bio::Tools::Alignment::Consed::dump_hash - ".
                  "The following is the contents of the contig hash...\n");
    $dumper->dumpValue($self->{'contigs'});
}

=head2 dump_hash_compact()

 Title   : dump_hash_compact()
 Usage   : $o_consed->dump_hash_compact();
 Function: Dump out the Bio::Tools::Alignment::Consed object in a compact way.
 Returns : Nothing.
 Args    : Nothing.
 Notes   : Cleaner then dumpValue(), dumpHash(). I used this a lot in
           debugging.

=cut

sub dump_hash_compact {
    no strict 'refs';           # renege for the block
    my ($self,$sequence) = @_;
    # get the classes
    my @singlets = $self->get_singlets();
    my @singletons = $self->get_singletons();
    my @doublets = $self->get_doublets();
    my @pairs = $self->get_pairs();
    my @multiplets = $self->get_multiplets();
    print("Name\tClass\tMembers\tQuality?\n");
    foreach (@singlets) {
        my @members = $self->get_members($_);
        print($self->get_name($_)."\tsinglets\t".(join',',@members)."\t");
        if ($self->{'contigs'}->{$_}->{'quality'}) {
            print("qualities found here\n");
        } else {
            print("no qualities found here\n");
        }

    }
    foreach (@singletons) {
        my @members = $self->get_members($_);
        print($self->get_name($_)."\tsingletons\t".(join',',@members)."\t");
        if ($self->{'contigs'}->{$_}->{'quality'}) {
            print("qualities found here\n");
        } else {
            print("no qualities found here\n");
        }
    }
    foreach my $pair (@pairs) {
        my @members = $self->get_members($pair);
        my $name;
        if (!$self->get_name($pair)) {
            $name = "BLANK";
        } else {
            $name = $self->get_name($pair);
        }
        print("$name\tpairs\t".(join',',@members)."\n");
    }
    foreach (@doublets) {
        my @members = $self->get_members_by_name($_);
        print("$_\tdoublets\t".(join',',@members)."\t");
        my $contig_number = &get_contig_number_by_name($self,$_);
        if ($self->{'contigs'}->{$contig_number}->{'quality'}) {
            print("qualities found here\n");
        } else {
            print("no qualities found here\n");
        }
        # print($_."\tdoublets\t".(join',',@members)."\n");
    }
    foreach (@multiplets) {
        my @members = $self->get_members($_);
        print("Contig $_"."\tmultiplets\t".(join',',@members)."\n");
    }
}                               # end dump_hash_compact

=head2 get_phreds()

 Title   : get_phreds()
 Usage   : @phreds = $o_consed->get_phreds();
 Function: For each doublet in the Bio::Tools::Alignment::Consed hash, go
           and get the phreds for the top and bottom reads. Place them into
           {top_phreds} and {bottom_phreds}.
 Returns : Nothing.
 Args    : Nothing.

Requires parse_phd() and reverse_and_complement(). I realize that it
would be much more elegant to pull qualities as required but there
were certain "features" in the acefile that required a bit more
detailed work be done to get the qualities for certain parts of the
consensus sequence. In order to make _sure_ that this was done
properly I wrote things to do all steps and then I used dump_hash()
and checked each one to ensure expected behavior. I have never changed
this, so there you are.

=cut

sub get_phreds {
    # this subroutine is the target of a rewrite to use the Bio::Tools::Alignment::Phred object.
    my $self = shift;    
    my $current_contig;
    foreach $current_contig (sort keys %{$self->{'contigs'}}) {	
	if ($self->{'contigs'}->{$current_contig}->{'class'} eq "doublet") {
	    $self->debug("$current_contig is a doublet. Going to parse_phd for top($self->{'contigs'}->{$current_contig}->{'top_name'}) and bottom($self->{'contigs'}->{$current_contig}->{'bottom_name'})\n");
	    my $r_phreds_top = &parse_phd($self,$self->{'contigs'}->{$current_contig}->{'top_name'});
	    my $r_phreds_bottom = &parse_phd($self,$self->{'contigs'}->{$current_contig}->{'bottom_name'});
	    if ($self->{'contigs'}->{$current_contig}->{'top_complement'} eq "C") {
		# print("Reversing and complementing...\n");
		$r_phreds_top = &reverse_and_complement($r_phreds_top);
	    }
	    if ($self->{'contigs'}->{$current_contig}->{'bottom_complement'} eq "C") {
		$r_phreds_bottom = &reverse_and_complement($r_phreds_bottom);
	    }
	    $self->{'contigs'}->{$current_contig}->{'top_phreds'} = $r_phreds_top;
	    $self->{'contigs'}->{$current_contig}->{'bottom_phreds'} = $r_phreds_bottom;
	}
    }
}

=head2 parse_phd($read_name)

 Title   : parse_phd($read_name)
 Usage   : $o_consed->parse_phd($read_name);
 Function: Suck in the contents of a .phd file.
 Returns : A reference to an array containing the quality values for the read.
 Args    : The name of a read.
 Notes   : This is a significantly weak subroutine because it was always
	   intended that these functions, along with the functions provided by
	   get_phreds() be put into the Bio::SeqIO:phd module. This is done
           now but the Bio::Tools::Alignment::Consed module has not be
           rewritten to reflect this change.

See L<Bio::SeqIO::phd> for more information.

=cut

sub parse_phd {
    my ($self,$sequence_name) = @_;
    $self->debug("Parsing phd for $sequence_name\n");
    my $in_dna = 0;
    my $base_number = 0;
    my (@bases,@current_line);
    # print("parse_phd: $sequence_name\n");
    my $fh = Bio::Root::IO->new
        (-file=>"$self->{path}/../phd_dir/$sequence_name.phd.1");
    while ($fh->_readline()) {
	# print("Reading a line from a phredfile!\n");
	chomp;
	if (/^BEGIN_DNA/) { $in_dna = 1; next}
	if (/^END_DNA/) { last; }
	if (!$in_dna) { next; }
	push(@bases,$_);
    }
    return \@bases;
}

=head2 reverse_and_complement(\@source)

 Title   : reverse_and_complement(\@source)
 Usage   : $reference_to_array = $o_consed->reverse_and_complement(\@source);
 Function: A stub for the recursive routine reverse_recurse().
 Returns : A reference to a reversed and complemented array of phred data.
 Args    : A reference to an array of phred data.
 Notes   : 

=cut

sub reverse_and_complement {
    my $r_source = shift;
    my $r_destination;
    $r_destination = &reverse_recurse($r_source,$r_destination);
    return $r_destination;
}

=head2 reverse_recurse($r_source,$r_destination)

 Title   : reverse_recurse(\@source,\@destination)
 Usage   : $o_consed->reverse_recurse(\@source,\@destination);
 Function: A recursive routine to reverse and complement an array of
           phred data.
 Returns : A reference to an array containing reversed phred data.
 Args    : A reference to a source array and a reverence to a destination
	   array.

Recursion is kewl, but this sub should likely be _reverse_recurse.

=cut


sub reverse_recurse($$) {
    my ($r_source,my $r_destination) = @_;
    if (!@$r_source) {
        return $r_destination;
    }
    $_=pop(@$r_source);
    s/c/g/ || s/g/c/ || s/a/t/ || s/t/a/;
    push(@$r_destination,$_);
    &reverse_recurse($r_source,$r_destination);
}

=head2 show_missing_sequence()

 Title   : show_missing_sequence();
 Usage   : $o_consed->show_missing_sequence();
 Function: Used by set_trim_points_doublets() to fill in quality values where
	   consed (phrap?) set them to 0 at the beginning and/or end of the
	   consensus sequences.
 Returns : Nothing.
 Args    : None.

Acts on doublets only. Really very somewhat quite ugly. A disgusting
kludge. I<insert pride here> It was written stepwise with no real plan
because it was not really evident why consed (phrap?)  was doing this.

=cut

sub show_missing_sequence() {

    # decide which sequence should not have been clipped at consensus
    # position = 0

    my $self = shift;
    &get_phreds($self);
    my ($current_contig,@qualities);
    foreach $current_contig (sort keys %{$self->{'contigs'}}) {
	if ($self->{'contigs'}->{$current_contig}->{'class'} eq "doublet") {
	    my $number_leading_xs = 0;
	    my $number_trailing_xs = 0;
	    my $measurer = $self->{'contigs'}->{$current_contig}->{'quality'};
	    while ($measurer =~ s/^\ 0\ /\ /) {
		$number_leading_xs++;
	    }
	    while ($measurer =~ s/\ 0(\s*)$/$1/) {
		$number_trailing_xs++;
	    }
	    @qualities = split(' ',$self->{'contigs'}->{$current_contig}->{'quality'});
	    my $in_initial_zeros = 0;
	    for (my $count=0;$count<scalar(@qualities); $count++) {
		if ($qualities[$count] == 0) {
		    my ($quality,$top_phred_position,$bottom_phred_position,$top_phred_data,$bottom_phred_data);
		    # print("The quality of the consensus at ".($count+1)." is zero. Retrieving the real quality value.\n");
		    # how do I know which strand to get these quality values from????
		    # boggle
		    my $top_quality_here = $self->{'contigs'}->{$current_contig}->{'top_phreds'}->[0-$self->{'contigs'}->{$current_contig}->{'top_start'}+$count+1];
		    my $bottom_quality_here = $self->{'contigs'}->{$current_contig}->{'bottom_phreds'}->[1-$self->{'contigs'}->{$current_contig}->{'bottom_start'}+$count];
		    if (!$bottom_quality_here || (1-$self->{'contigs'}->{$current_contig}->{'bottom_start'}+$count)<0) {
			$bottom_quality_here = "not found";
		    }
		    if (!$top_quality_here) {
			$top_quality_here = "not found";
		    }
		    # print("Looking for quals at position $count of $current_contig: top position ".(0-$self->{'contigs'}->{$current_contig}->{top_start}+$count)." ($self->{'contigs'}->{$current_contig}->{top_name}) $top_quality_here , bottom position ".(1-$self->{'contigs'}->{$current_contig}->{bottom_start}+$count)." ($self->{'contigs'}->{$current_contig}->{bottom_name}) $bottom_quality_here\n"); 
		    if ($count<$number_leading_xs) {
			# print("$count is less then $number_leading_xs so I will get the quality from the top strand\n");
			# print("retrieved quality is ".$self->{'contigs'}->{$current_contig}->{top_phreds}[0-$self->{'contigs'}->{$current_contig}->{top_start}+$count+1]."\n");
			my $quality = $top_quality_here;
			$quality =~ /\S+\s(\d+)\s+/;
			$quality = $1;
			# print("retrieved quality for leading zero $count is $quality\n");
			# t 9 9226
			$qualities[$count] = $quality;
		    } else {
			# this part is tricky
			# if the contig is like this
			#      cccccccccccccccc
			# ffffffffffffffffff
			#          rrrrrrrrrrrrrrrrr
			# then take the quality value for the trailing zeros in the cons. seq from the r
			#
			# but if the contig is like this
			#      cccccccccccccccccc
			#      ffffffffffffffffffffffffffffffff
			# rrrrrrrrrrrrrrrrrrrrrrrxxxxxxxxr
			#                      ^^^
			# then any zeros that fall in the positions (^) must be decided whether the quality
			# is the qual from the f or r strand. I will use the greater number
			# does a similar situation exist for the leading zeros? i dunno
			#
			# print("$count is greater then $number_leading_xs so I will get the quality from the bottom strand\n");
			# print("retrieved quality is ".$contigs->{$current_contig}->{top_phreds}[0-$contigs->{$current_contig}->{top_start}+$count+1]."\n");
			# my ($quality,$top_phred_position,$bottom_phred_position,$top_phred_data,$bottom_phred_data);
			if ($bottom_quality_here eq "not found") {
			    # $top_phred_position = 1-$contigs->{$current_contig}->{bottom_start}+$count;
			    # print("Going to get quality from here: $top_phred_position of the top.\n");
			    # my $temp_quality - $contigs->{$current_contig}->{top_phreds}
			    # $quality = $contigs->{$current_contig}->{top_phreds}[$top_phred_position];
			    $top_quality_here =~ /\w+\s(\d+)\s/;
			    $quality = $1;
			} elsif ($top_quality_here eq "not found") {
			    # $bottom_phred_position = 1+$contigs->{$current_contig}->{bottom_start}+$count;
			    # print("Going to get quality from here: $bottom_phred_position of the bottom.\n");
			    # $quality = $contigs->{$current_contig}->{bottom_phreds}[$bottom_phred_position];
			    # print("Additional: no top quality but bottom is $quality\n");
			    $bottom_quality_here =~ /\w+\s(\d+)\s/;
			    $quality = $1;
			} else {
			    # print("Oh jeepers, there are 2 qualities to choose from at this position.\n");
			    # print("Going to compare these phred qualities: top: #$top_quality_here# bottom: #$bottom_quality_here#\n");
			    # now you have to compare them
			    # my $top_quality_phred = $contigs->{$current_contig}->{top_phreds}[$top_phred_position];
			    # #t 40 875#
			    # print("regexing #$top_quality_here#... ");
			    $top_quality_here =~ /\w\ (\d+)\s/;
			    my $top_quality = $1;
			    # print("$top_quality\nregexing #$bottom_quality_here#... ");
			    $bottom_quality_here =~ /\w\ (\d+)\s/;
			    my $bottom_quality = $1;
			    # print("$bottom_quality\n");
			    # print("top_quality: $top_quality bottom quality: $bottom_quality\n");
			    if ($bottom_quality > $top_quality) {
				# print("Chose to take the bottom quality: $bottom_quality\n");
				$quality = $bottom_quality;
			    } else {
				# print("Chose to take the top quality: $top_quality\n");
				$quality = $top_quality;
			    }
			}
			if (!$quality) {
			    # print("Warning: no quality value for $current_contig, position $count!\n");
			    # print("Additional data: top quality phred: $top_quality_here\n");
			    # print("Additional data: bottom quality phred: $bottom_quality_here\n");
			} else {
			    $qualities[$count] = $quality;
			}
		    }						
		}

	    }
	    unless (!@qualities) {
		$self->{'contigs'}->{$current_contig}->{'quality'} = join(" ",@qualities);
	    }
	    $self->{'contigs'}->{$current_contig}->{'bottom_phreds'} = undef;
	    $self->{'contigs'}->{$current_contig}->{'top_phreds'} = undef;
	    my $count = 1;
	}			# end foreach key
    }
}


1;
