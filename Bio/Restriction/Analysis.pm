#------------------------------------------------------------------
#
# BioPerl module Bio::Restriction::Analysis
#
# Cared for by Rob Edwards <redwards@utmem.edu>
#
# You may distribute this module under the same terms as perl itself
#------------------------------------------------------------------

## POD Documentation:

=head1 NAME

Bio::Restriction::Analysis - cutting sequences with restriction
enzymes

=head1 SYNOPSIS

  # analyze a DNA sequence for restriction enzymes
  use Bio::Restriction::Analysis;
  use Bio::PrimarySeq;
  use Data::Dumper;

  # get a DNA sequence from somewhere
  my $seq=new Bio::PrimarySeq
      (-seq =>'AGCTTAATTCATTAGCTCTGACTGCAACGGGCAATATGTCTC'.
       'TGTGTGGATTAAAAAAAGAGTGAGCTTCTGATAGCAGC',
       -primary_id => 'synopsis',
       -molecule => 'dna');

  # now start an analysis.
  # this is using the default set of enzymes
  my $ra=Bio::Restriction::Analysis->new(-seq=>$seq);

  # find unique cutters. This returns a
  # Bio::Restriction::EnzymeCollection object
  my $enzymes=$ra->unique_cutters;
  print "Unique cutters: ", join (', ', 
      map {$_->name} $enzymes->each_enzyme), "\n";

  # AluI is one them. Where does it cut?
  # This is will return an array of the sequence strings

  my $enz = 'AluI';
  my @frags=$ra->fragments($enz);
  # how big are the fragments?
  print "AluI fragment lengths: ", join(' & ', map {length $_} @frags), "\n";

  # You can also bypass fragments and call sizes directly:
  # to see all the fragment sizes
  print "All sizes: ", join " ", $ra->sizes($enz), "\n";
  # to see all the fragment sizes sorted
  print "All sizes, sorted ", join (" ", $ra->sizes($enz, 0, 1)), "\n";


  # how many times does each enzyme cut
  my $cuts=$ra->cuts_by_enzyme('BamHI');
  print "BamHI cuts $cuts times\n";

  # How many enzymes do not cut at all?
  print "There are ", scalar $ra->zero_cutters->each_enzyme,
        " enzymes that don't cut\n";

  # what about enzymes that cut twice?
  my $two_cutters=$ra->cutters(2);
  print join (" ", map {$_->name} $two_cutters->each_enzyme),
      " cut the sequence twice\n";

  # what are all the enzymes that cut, and how often do they cut
  printf "\n%-10s%s\n", 'Enzyme', 'Number of Cuts';
  my $all_cutters=$ra->cutters;
  map {
      printf "%-10s%s\n", $_->name, $ra->cuts_by_enzyme($_->name)
  } $all_cutters->each_enzyme;

  # Finally, we can interact the restriction enzyme object by
  # retrieving it from the collection object see the docs for
  # Bio::Restriction::Enzyme.pm
  my $enzobj=$enzymes->get_enzyme($enz);


=head1 DESCRIPTION

Bio::Restriction::Analysis describes the results of cutting a DNA
sequence with restriction enzymes.

To use this module you can pass a sequence object and optionally a
Bio::Restriction::EnzymeCollection that contains the enzyme(s) to cut the
sequences with. There is a default set of enzymes that will be loaded
if you don't pass in a Bio::Restriction::EnzymeCollection.

To cut a sequence, set up a Restriction::Analysis object with a sequence
like this:

  use Bio::Restriction::Analysis;
  my $ra=Bio::Restriction::Analysis->new(-seq=>$seqobj);

or

  my $ra=Bio::Restriction::Analysis->new
      (-seq=>$seqobj, -enzymes=>$enzs);

Then, to get the fragments for a particular enzyme use this:

  @fragments=$ra->fragments('EcoRI');

Note that the naming of restriction enzymes is that the last numbers
are usually Roman numbers (I, II, III, etc). You may want to use
something like this:

  # get a reference to an array of unique (single) cutters
  $singles = $re->unique_cutters;
  foreach my $enz ($singles->each_enzyme) {
      @fragments=$re->fragments($enz);
      ... do something here ...
  }

Note that if your sequence is circular, the first and last fragment
will be joined so that they are the appropriate length and sequence
for further analysis. However, this will change the start of the
sequence!

=head1 COMMENTS

I am trying to make this backwards compatible with
Bio::Tools::Restriction::Analysis Undoubtedly some things will break,
but we can fix things as we progress.....!

=head1 FEEDBACK

=head2 Mailing Lists 

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists. Your participation is much appreciated.

    bioperl-l@bioperl.org             - General discussion
    http://bioperl.org/MailList.shtml - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution. Bug reports can be submitted via email
or the web:

     bioperl-bugs@bio.perl.org
     http://bugzilla.bioperl.org/

=head1 AUTHOR

Rob Edwards, redwards@utmem.edu, 
Steve Chervitz, sac@bioperl.org

=head1 CONTRIBUTORS

Heikki Lehvaslaiho, heikki@ebi.ac.uk

=head1 COPYRIGHT

Copyright (c) 2003 Rob Edwards.  Some of this work is Copyright (c)
1997-2002 Steve A. Chervitz. All Rights Reserved.

This module is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=head1 SEE ALSO

L<Bio::Restriction::Enzyme>, 
L<Bio::Restriction::EnzymeCollection>

=head1 APPENDIX

Methods beginning with a leading underscore are considered private and
are intended for internal use by this module. They are not considered
part of the public interface and are described here for documentation
purposes only.

=cut

package Bio::Restriction::Analysis;
use Bio::Restriction::EnzymeCollection;
use Bio::Root::Root;
use strict;
use Data::Dumper;

use vars qw (@ISA);
@ISA = qw(Bio::Root::Root);

=head1 new

 Title     : new
 Function  : Initializes the restriction enzyme object
 Returns   : The Restriction::Analysis object 
 Arguments : 

	     $re_anal->new(-seq=$seqobj, 
                 -enzymes=>Restriction::EnzymeCollection object)
	     -seq requires a Bio::PrimarySeq object
	     -enzyme_collection is optional.
              If ommitted it will use the default set of enzymes

This is the place to start. Pass in a sequence, and you will be able
to get the fragments back out.  Several other things are available
like the number of zero cutters or single cutters.

=cut

sub new {
    my($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($seq,$enzymes) =
        $self->_rearrange([qw(
                              SEQ
                              ENZYMES
                             )], @args);

    $seq && $self->seq($seq);

    $enzymes ?  $self->enzymes($enzymes)
        :  ($self->{'_enzymes'} = Bio::Restriction::EnzymeCollection->new );

    # keep track of status
    $self->{'_cut'} = 0;

    return $self;

}

=head2 Methods to set parameters

=cut

=head2 seq

 Title    : seq
 Usage    : $ranalysis->seq($newval);
 Function : get/set method for the  sequence to be cut
 Example  : $re->seq($seq);
 Returns  : value of seq
 Args     : A Bio::PrimarySeqI dna object (optional)

=cut

sub seq {
     my $self = shift;
     if (@_) {
         my $seq = shift;
         $self->throw('Need a sequence object ['. ref $seq.  ']')
             unless $seq->isa('Bio::PrimarySeqI');
         $self->throw('Need a DNA sequence object ['. $seq->alphabet.  ']')
             unless $seq->alphabet eq 'dna';

         $self->{'_seq'} = $seq;
         $self->{'_cut'} = 0;
     }
     return $self->{'_seq'};
}

=head2 enzymes

 Title    : enzymes
 Usage    : $re->enzymes($newval)
 Function : gets/Set the restriction enzyme enzymes
 Example  : $re->enzymes('EcoRI')
 Returns  : reference to the collection
 Args     : an array of Bio::Restriction::EnzymeCollection and/or
            Bio::Restriction::Enzyme objects


The default object for this method is
Bio::Restriction::EnzymeCollection.  However, you can also pass it a
list of Bio::Restriction::Enzyme objects - even mixed with Collection
objects.  They will all be stored into one collection.

=cut

sub enzymes {
     my $self = shift;
     if (@_) {
         $self->{'_enzymes'} = Bio::Restriction::EnzymeCollection->new (-empty => 1)
             unless $self->{'_enzymes'};
         $self->{'_enzymes'}->enzymes(@_);
         $self->{'_cut'} = 0;
     }
     return $self->{'_enzymes'};
}


=head2 Perform the analysis

=cut

=head2 cut

 Title    : cut
 Usage    : $re->cut()
 Function : Cut the sequence with the enzymes
 Example  : $re->cut(); $re->cut('original');
 Returns  : $self
 Args     : 'original', optional

An explicit cut method is needed to pass arguments to it. The default
behaviour is to cut with each enzymes in the list separately. The
(unimplemented) option is to perform the cut with all the enzymes at
the same time: "double digest".

If you want to use the default setting there is no need to call cut
directly. Every method in the class that needs output checks the
object's internal status and recalculates the cuts if needed.

While the module is under construction, passing the string 'original'
as an argument, will run the code from older
Bio::Tools::RestrictionEnzyme module. Otherwise, the newer, untested
code will be executed. You have been warned.

=cut

sub cut {
    my ($self, $opt) = @_;

    $self->throw("A sequence must be supplied")
        unless $self->seq;

    $self->{maximum_cuts} = 0;

    $self->{'_number_of_cuts_by_enzyme'} = {};
    $self->{_number_of_cuts_by_cuts} = {};
    $self->{'_fragments'} = {};

    if ($opt and $opt eq 'original') {
        $self->_cuts;
    } else {
        $self->_new_cuts;
    }
    $self->{'_cut'} = 1;
    return $self;
}


=head2 Query the results of the analysis

=cut

=head2 fragments

  Title    : fragments
  Function : Retrieve the fragments that we cut
  Returns  : At the moment just returns a reference to an array of
             sequences, but I'd like this to return a collection of
             Bio::Seq objects, perhaps as a stream?  Returns 0 if not
             defined.
  Arguments: An enzyme name to retrieve the fragments for

For example this code will retrieve the fragments for all enzymes that
cut your sequence

  my $all_cutters = $analysis->cutters;
  foreach my $enz ($$all_cutters->each_enzyme}) {
      @fragments=$analysis->fragments($enz);
  }

=cut

sub fragments {
    my ($self, $enz) = @_;
    $self->cut unless $self->{'_cut'};
    $self->throw('no enzyme selected to get fragments for')
        unless $enz;
    return @{$self->{'_fragments'}->{$enz}};
}

=head2 sizes

  Title    : sizes
  Function : Retrieves an array with the sizes of the fragments
  Returns  : 

             Array that has the sizes of the fragments ordered from
             largest to smallest like they would appear in a gel.

  Arguments: 

             An enzyme name to retrieve the sizes for is required and
             if the optional second entry is set the sizes will be in
             kilobases to the nearest 0.1 kb, else it will be in
             bp. If the optional third entry is set the results will
             be sorted.

This is designed to make it easy to see what fragments you should get
on a gel!

You should be able to do these:

  # to see all the fragment sizes,
  print join "\n", @{$re->sizes($enz)}, "\n";
  # to see all the fragment sizes sorted
  print join "\n", @{$re->sizes($enz, 0, 1)}, "\n";
  # to see all the fragment sizes in kb sorted
  print join "\n", @{$re->sizes($enz, 1, 1)}, "\n";

=cut

sub sizes {
    my ($self, $enz, $kb, $sort) = @_;
    $self->throw('no enzyme selected to get fragments for')
        unless $enz;
    $self->cut unless $self->{'_cut'};
    my @frag;
    if (defined $self->{'_fragments'}->{$enz}) {
        foreach my $frag (@{$self->{'_fragments'}->{$enz}}) {
            my $len=length($frag);
            if ($kb) {
                $len=(int($len/100))/10;
            }
            push @frag, $len;
        }
        if ($sort) {
            @frag = sort {$b <=> $a} @frag;
        }
    }
    return @frag;
}

=head2 How many times does enzymes X cut?

=cut

=head2 cuts_by_enzyme

 Title     : cuts_by_enzyme
 Function  : Return the number of cuts for an enzyme
 Returns   : An integer with the number of times each enzyme cuts.
             Returns 0 if doesn't cut or undef if not defined
 Arguments : An enzyme name string


=cut

sub cuts_by_enzyme {
    my ($self, $enz)=@_;

    $self->throw("Need an enzyme name")
        unless defined $enz;
    $self->cut unless $self->{'_cut'};
    return $self->{'_number_of_cuts_by_enzyme'}->{$enz};
}

=head2 Which enzymes cut the sequence N times?

=cut

=head2 cutters

 Title     : cutters
 Function  : Find enzymes that cut a given number of times
 Returns   : a Bio::Restriction::EnzymeCollection
 Arguments : 1. exact time or lower limit,
                non-negative integer, optional
             2. upper limit, non-negative integer,
                larger or equalthan first, optional


If no argumets are given, the method returns all enzymes that do cut
the sequence. The argument zero, '0', is same as method
zero_cutters().  The argument one, '1', corresponds to unique_cutters.
If either of the limits is larger than number of cuts any enzyme cuts the
sequence, the that limit is automagically lowered. The method max_cuts()
gives the largest number of cuts.

See Also : L<unique_cutters|unique_cutters>,
L<zero_cutters|zero_cutters>, L<max_cuts|max_cuts>

=cut

sub cutters {
    my ($self, $a, $z) = @_;

    $self->cut unless $self->{'_cut'};

    my ($start, $end);
    if (defined $a) {
        $self->throw("Need a non-zero integer [$a]")
            unless $a =~ /^[+]?\d+$/;
        $start = $a;
    } else {
        $start = 1;
    }
    $start = $self->{'maximum_cuts'} if $start > $self->{'maximum_cuts'};

    if (defined $z) {
        $self->throw("Need a non-zero integer no smaller than start [0]")
            unless $z =~ /^[+]?\d+$/ and $z >= $a;
        $end = $z;
    }
    elsif (defined $a) {
        $end = $start;
    } else {
        $end = $self->{'maximum_cuts'};
    }
    $end = $self->{'maximum_cuts'} if $end > $self->{'maximum_cuts'};
    my $set = new Bio::Restriction::EnzymeCollection(-empty => 1);

    for (my $i=$start; $i<=$end; $i++) {
        $set->enzymes( @{$self->{_number_of_cuts_by_cuts}->{$i}} )
            if defined $self->{_number_of_cuts_by_cuts}->{$i};
    }

    return $set;
}


=head2 unique_cutters

 Title     : unique_cutters
 Function  : A special case if cutters() where enzymes only cut once
 Returns   : a Bio::Restriction::EnzymeCollection
 Arguments : -


See also:  L<cutters>, L<zero_cutters>

=cut

sub unique_cutters {
    shift->cutters(1);
}

=head2 zero_cutters

 Title     : zero_cutters
 Function  : A special case if cutters() where enzymes don't cut the sequence
 Returns   : a Bio::Restriction::EnzymeCollection
 Arguments : -

See also:  L<cutters>, L<unique_cutters>

=cut

sub zero_cutters {
    shift->cutters(0);
}

=head2 max_cuts

 Title     : max_cuts
 Function  : Find enzymes that cut a given number of times
 Returns   : a Bio::Restriction::EnzymeCollection
 Arguments : exact time or lower limit, non-negative integer, optional
             upper limit, non-negative integer, optional

This is not a very practical method, but if you are curious...

=cut

sub max_cuts { return shift->{maximum_cuts} }







=head2 Internal methods

=cut

=head2 _cuts

 Title     : _cuts
 Function  : Figures out which enzymes we know about and cuts the seequence
 Returns   : The Restriction::Analysis object 
 Arguments : This is an internal method that cuts the sequence and whatnot.
 Comments  : 

=cut


sub _cuts {
    my $self = shift;

    my $target_seq = uc $self->{_seq}->seq; # I have been burned on this before :)

    # first, find out all the enzymes that we have
    foreach my $enz ($self->{'_enzymes'}->each_enzyme) {
        # this part of the code came from Steve's module....originally
        my $cuts_after = $enz->cut;
        my ($site_3prime_seq, $site_5prime_seq);
        my $reSeq = $enz->seq;
        if ($cuts_after == 0) {
            $site_3prime_seq = '';
            $site_5prime_seq = $reSeq->seq();
        } elsif ($cuts_after == $reSeq->length) {
            $site_3prime_seq = $reSeq->seq();
            $site_5prime_seq = '';
        } elsif ($cuts_after > 0 && $cuts_after < $reSeq->length) {
            # it cuts within the sequence
            $site_3prime_seq = $reSeq->subseq(1, $enz->cut);
            $site_5prime_seq = $reSeq->subseq($enz->cut+1, $reSeq->length);
        } else {
            # man, we are screwed here, I am coming back to this
            if ($self->verbose) {
                $self->warn("Cutting outside the sequence is not implemented yet");
            }
            next;
        }

        my(@re_frags);
        my $seq = uc $self->_expanded_string($enz->string);

        if (!$enz->palindromic) {
            my $revseq = $self->_expanded_string($reSeq->revcom);
            $seq .= '|'.uc($revseq);
        }


        @re_frags = split(/$seq/i, $target_seq);
        ## Re-attach the split recognition site back to the frags
        ## since perl zapped them in the split() call.
        my($i);
        my $numFrags = scalar @re_frags;
        for ($i=0; $i<$numFrags; $i++) {
            $i < $#re_frags  and $re_frags[$i] = $re_frags[$i].$site_3prime_seq;
            $i > 0           and $re_frags[$i] = $site_5prime_seq.$re_frags[$i];
        }

        # get the number of times this enzyme cuts the sequence
        my $number_of_cuts = $#re_frags;
        $self->{_number_of_cuts_by_enzyme}->{$enz->name}=$number_of_cuts;
        @{$self->{$enz}->{fragments}} = @re_frags;
        push (@{$self->{_number_of_cuts_by_cuts}->{$number_of_cuts}}, $enz);
        if ($number_of_cuts > $self->{maximum_cuts}) {
            $self->{maximum_cuts}=$number_of_cuts;
        }
    }
}



=head2 _new_cuts

 Title     : _new_cuts
 Function  : Figures out which enzymes we know about and cuts the sequence.
 Returns   : Nothing.
 Arguments : None.
 Comments  : An internal method. This is the development part of rewriting _cuts().
             This will figure out where the sequence should be cut, and provide the 
	     appropriate fragments.

I have left _cuts() in for now, so I can compare this rewritten
version with the original to see that I have done it right!

=cut



sub _new_cuts {
    my $self = shift;

    my $target_seq=uc $self->{'_seq'}->seq; # I have been burned on this before :)


    # first, find out all the enzymes that we have
    foreach my $enz ($self->{'_enzymes'}->each_enzyme) {
        # get the cut site. This will have a ^ and N's in it.
        my $site=$enz->site;
        # split it into the two fragments for the sequence before and after.

        # BIG assumption - if the enzyme doesn't have a cut site, we're just guessing.
        # we should really throw a warning here
        unless ($site =~ /\^/) {
            $site .= "^";
            $self->warn("No cut site is known for $enz. \n".
                        "Therefore we are guessing that it is $site\n")
                if $self->verbose;
        }

        # the default values just stop an error from an undefined
        # string. But they don't affect the split.
        my ($beforeseq, $afterseq)= ('.', '.');
        # now perlify the sequences into a regexp
        ($beforeseq, $afterseq)=split /\^/, $self->_expanded_string($site);

        # now cut the sequence
        my @cuts = split /($beforeseq)($afterseq)/i, $target_seq;
        # now the array has extra elements --- the before and after!
        # we have:
        # element 0 sequence
        # element 1 3' end
        # element 2 5' end of next sequence
        # element 3 sequence
        # ....

        # we need to loop through the array and add the ends to the
        # appropriate parts of the sequence
        my $i=0;
        my @re_frags;

        if ($#cuts) {           # there is >1 element
            while ($i<$#cuts) {
                my $joinedseq;
                # the first sequence is a special case
                if ($i == 0) {
                    $joinedseq=$cuts[$i].$cuts[$i+1];
                } else {
                    $joinedseq=$cuts[$i-1].$cuts[$i].$cuts[$i+1];
                }
                push @re_frags, $joinedseq;
                $i+=3;
            }


            # for the last cut, we want to join the last two elements only
            # because of the last $i+=3, we should now be at the last position
            my $joinedseq;

            # there is a rare and unusual case where if the enzyme
            # cuts on the last base, it doesn't leave a sequence in
            # $i!

            unless ($cuts[$i]) {
                $joinedseq=$cuts[$i-1];
            } else {
                $joinedseq=$cuts[$i-1].$cuts[$i];
            }
            push @re_frags, $joinedseq;


            # now find out if the sequence is circular. If so, join
            # the first and last fragments

            if ($self->{'_seq'}->is_circular) {
                my $first = shift @re_frags;
                my $last = pop @re_frags;

                unshift (@re_frags, $last.$first);
            }

        } else {
            @re_frags=@cuts;
        }                       # the sequence was not cut.

        # get the number of times this enzyme cuts the sequence
        my $number_of_cuts = $#re_frags;
        $self->{_number_of_cuts_by_enzyme}->{$enz->name}=$number_of_cuts;
        @{$self->{'_fragments'}->{$enz->name}} = @re_frags;
        push (@{$self->{_number_of_cuts_by_cuts}->{$number_of_cuts}}, $enz);
        if ($number_of_cuts > $self->{maximum_cuts}) {
            $self->{maximum_cuts}=$number_of_cuts;
        }
    }
}


=head2 _expanded_string

 Title     : _expanded_string
 Function  : Expand nucleotide ambiguity codes to their representative letters
 Returns   : The full length string
 Arguments : The string to be expanded.

Stolen from the original RestrictionEnzyme.pm

=cut


sub _expanded_string {
    my ($self, $str) = @_;

    $str =~ s/N|X/\./g;
    $str =~ s/R/\[AG\]/g;
    $str =~ s/Y/\[CT\]/g;
    $str =~ s/S/\[GC\]/g;
    $str =~ s/W/\[AT\]/g;
    $str =~ s/M/\[AC\]/g;
    $str =~ s/K/\[TG\]/g;
    $str =~ s/B/\[CGT\]/g;
    $str =~ s/D/\[AGT\]/g;
    $str =~ s/H/\[ACT\]/g;
    $str =~ s/V/\[ACG\]/g;

    return $str;
}


1;
