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
      map {$_->name} $enzymes->unique_cutters), "\n";

  # AluI is one them. Where does it cut?
  # This is will return an array of the sequence strings

  my $enz = 'AluI';
  my @frags=$ra->fragments($enz);
  # how big are the fragments?
  print "AluI fragment lengths: ", join(' & ', map {length $_} @frags), "\n";

  # You can also bypass fragments and call sizes directly:
  # to see all the fragment sizes
  print "All sizes: ", join " ", $ra->sizes($enz), "\n";
  # to see all the fragment sizes sorted by size like on a gel
  print "All sizes, sorted ", join (" ", $ra->sizes($enz, 0, 1)), "\n";

  # how many times does each enzyme cut
  my $cuts=$ra->cuts_by_enzyme('BamHI');
  print "BamHI cuts $cuts times\n";

  # How many enzymes do not cut at all?
  print "There are ", scalar $ra->zero_cutters->each_enzyme,
        " enzymes that do not cut\n";

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
if you do not pass in a Bio::Restriction::EnzymeCollection.

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
for further analysis. This fragment will also be checked for cuts
by the enzyme(s).  However, this will change the start of the
sequence!

There are two separate algorithms used depending on whether your
enzyme has ambiguity. The non-ambiguous algoritm is a lot faster,
and if you are using very large sequences you should try and use
this algorithm. If you have a large sequence (e.g. genome) and 
want to use ambgiuous enzymes you may want to make seperate
Bio::Restriction::Enzyme objects for each of the possible
alternatives and make sure that you don't set is_ambiguous!

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
	     -enzymes is optional.
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
    
    $self->{maximum_cuts} = 0;

    $self->{'_number_of_cuts_by_enzyme'} = {};
    $self->{'_number_of_cuts_by_cuts'} = {};
    $self->{'_fragments'} = {};
    $self->{'_cut_positions'} = {}; # cut position is the real position 

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
 Example  : $re->cut(); $re->cut('single'); or $re->cut('multiple', $enzymecollection);
 Returns  : $self
 Args     : 'single' (optional), 'multiple' with enzyme collection.

An explicit cut method is needed to pass arguments to it. 

There are two varieties of cut. Single is the default, and need
not be explicitly called. This cuts the sequence with each
enzyme separately.

Multiple cuts a sequence with more than one enzyme. You must pass
it a Bio::Restriction::EnzymeCollection object of the set
of enzymes that you want to use in the double digest. The results
will be stored as an enzyme named "multiple_digest", so you can
use all the retrieval methods to get the data.

If you want to use the default setting there is no need to call cut
directly. Every method in the class that needs output checks the
object's internal status and recalculates the cuts if needed.

Note: cut doesn't now re-initialize everything before figuring
out cuts. This is so that you can do multiple digests, or add more
data or whatever. You'll have to use new to reset everything.

See also the comments in above about ambiguous and non-ambiguous
sequences.

=cut

sub cut {
    my ($self, $opt, $ec) = @_;

    # for the moment I have left this as a separate routine so
    # the user calls cuts rather than _cuts. This also initializes
    # some stuff we need to use.
  
    $self->throw("A sequence must be supplied")
        unless $self->seq;

    if (uc($opt) eq "MULTIPLE") {
      $self->throw("You must supply a separate enzyme collection for multiple digests") unless $ec;
      $self->_multiple_cuts($ec); # multiple digests
    } else {$self->_cuts}

    $self->{'_cut'} = 1;
    return $self;
}


=head2 Query the results of the analysis

=cut

=head2 positions

  Title    : positions
  Function : Retrieve the positions that an enzyme cuts at
  Returns  : An array of the positions that an enzyme cuts at
  Arguments: An enzyme name to retrieve the positions for
  Comments : The cut occurs after the base specified.
             Returns an empty array if the enzyme doesn't cut

=cut

sub positions {
    my ($self, $enz) = @_;
    $self->cut unless $self->{'_cut'};
    $self->throw('no enzyme selected to get positions for')
        unless $enz;
    
    if ($self->{'_cut_positions'}->{$enz}) {
       return @{$self->{'_cut_positions'}->{$enz}};
    }
    else {
       return ();
    }
}

=head2 fragments

  Title    : fragments
  Function : Retrieve the fragments that we cut
  Returns  : An array of the fragments retrieved. 
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
    
    # get the default cut positions
    my $cp= $self->{'_cut_positions'}->{$enz};
    
    my @fragments; 
    # if the sequence was not cut, $cp is not defined, 
    # and we can just return the whole sequence
    unless (defined $$cp[0]) {
       push (@fragments, $self->{'_seq'}->seq);
       return @fragments;
    }
    
    # make sure the list is in order. It should be for
    # most everything except multiple cuts.
    @$cp=sort {$a <=> $b} @$cp;

    # add each of the fragments to an array
    # we want to add the first and last fragments separately
    # add first fragment
    push @fragments, $self->{'_seq'}->subseq(1, ${$cp}[0]);
    # add all the intermediate fragments
    for (my $i=0; $i<$#$cp; $i++) {
      push (@fragments, $self->{'_seq'}->subseq(${$cp}[$i]+1, ${$cp}[$i+1]));
    }
    # add the last fragment. Remember this is not normally in $cp
    push @fragments, $self->{'_seq'}->subseq(${$cp}[$#$cp]+1, $self->{'_seq'}->length);	

    # if the sequence is circular we want to add the last and first
    # sequences together. UNLESS we have a site at the end of the sequence
    if ($self->{'_seq'}->is_circular && ${$cp}[$#$cp] != $self->{'_seq'}->length) {
     my ($first, $last)=(shift @fragments, pop @fragments);
     unshift @fragments, $last.$first;
    }

    $self->{'_fragments'}->{$enz}=\@fragments;
    # now we can return either a reference to an array or an array of references!
    return @fragments;
}

=head2 sizes

  Title    : sizes
  Function : Retrieves an array with the sizes of the fragments
  Returns  : Array that has the sizes of the fragments ordered from 
             largest to smallest like they would appear in a gel.
  Arguments: An enzyme name to retrieve the sizes for is required and
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
    $self->fragments($enz) unless (defined $self->{'_fragments'}->{$enz});
    
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
 Function  : Find the most number of cuts
 Returns   : The number of times the enzyme that cuts most cuts.
 Arguments : None

This is not a very practical method, but if you are curious...

=cut

sub max_cuts { return shift->{maximum_cuts} }







=head2 Internal methods

=cut

=head2 _cuts

 Title     : _cuts
 Function  : Figures out which enzymes we know about and cuts the sequence.
 Returns   : Nothing.
 Arguments : None.
 Comments  : An internal method. This will figure out where the sequence 
             should be cut, and provide the appropriate results.

=cut



sub _cuts {
    my $self = shift;

    my $target_seq=uc $self->{'_seq'}->seq; # I have been burned on this before :)


    # first, find out all the enzymes that we have
    foreach my $enz ($self->{'_enzymes'}->each_enzyme) {
        my ($beforeseq, $afterseq)=$self->_decon_enzyme($enz);
        # cut the sequence
	# this returns the number of cuts per enzyme

        # if the enzyme is ambiguous we need to use a regexp to find the cut site
	# otherwise we can use index (much faster)
	
	my $cut_positions;
	if ($enz->is_ambiguous) {
	   $cut_positions= $self->_ambig_cuts($beforeseq, $afterseq, $target_seq, $enz);
        } else {
	   $cut_positions= $self->_nonambig_cuts($beforeseq, $afterseq, $target_seq, $enz);
	}
	
	# deal with is_circular
	my $more_cut_positions=$self->_circular($beforeseq, $afterseq, $enz) if ($self->{'_seq'}->is_circular);

	# now deal with MultiSite and MultiCut enzymes
	# note I decided to do it this way because if you just cut with all the 
	# enzymes you'll get the default sites back, but I'll also add methods to retrieve
	# alternate data for multi enzymes

	# we do everything here that we did above, but without the explanation!!
        my @all_cuts;
	if ($enz->isa("Bio::Restriction::Enzyme::MultiSite") || $enz->isa("Bio::Restriction::Enzyme::MultiCut")) {
	     foreach my $other ($enz->others) {
    	          my ($beforeseq, $afterseq)=$self->_decon_enzyme($other);
        	  my $more_results=[];
	          if ($other->is_ambiguous) {
	              $more_results= $self->_ambig_cuts($beforeseq, $afterseq, $target_seq, $enz);
                  } else {
        	      $more_results= $self->_nonambig_cuts($beforeseq, $afterseq, $target_seq, $enz);
	          }
          
	          # deal with is_circular
	          my $more_more=[];
                  $more_more=$self->_circular($beforeseq, $afterseq, $enz) if ($self->{'_seq'}->is_circular);
                  my $non_pal_others=[];
                  unless ($other->is_palindromic) {
                       $non_pal_others=$self->_non_pal_enz($target_seq, $other);
                  }
	          push (@all_cuts, @$more_more, @$more_results, @$non_pal_others);
             }
	}
 
        # we need to deal with non-palindromic enzymes separately
        my $non_pal=[];
	unless ($enz->is_palindromic) {
	     $non_pal=$self->_non_pal_enz($target_seq, $enz);
         }

        # I pulled this out as a separate routine as I was using it a couple of times, then decided I only need it once.
        $self->{'_cut_positions'}->{$enz->name}=$self->_decon_cuts(@$cut_positions, @$more_cut_positions, @all_cuts, @$non_pal);
        my $number_of_cuts=scalar @{$self->{'_cut_positions'}->{$enz->name}};
    
        # now just store the number of cuts
	$self->{_number_of_cuts_by_enzyme}->{$enz->name}=$number_of_cuts;
        push (@{$self->{_number_of_cuts_by_cuts}->{$number_of_cuts}}, $enz);
        if ($number_of_cuts > $self->{maximum_cuts}) {
            $self->{maximum_cuts}=$number_of_cuts;
        }

    }
}

=head2 _non_pal_enz
  Title    : _non_pal_enz
  Function : Analyses non_palindromic enzymes for cuts in both ways
  Returns  : A reference to an array of cut positions
  Arguments: The sequence to check and the enzyme object

=cut

sub _non_pal_enz {
    my ($self, $target_seq, $enz) =@_;
    # add support for non-palindromic sequences
    # the enzyme is not the same forwards and backwards
    my $enzseq=$enz->revcom_site;
    unless ($enzseq =~ /\^/) {$enzseq .= "^"}
    my ($beforeseq, $afterseq)=('.', '.');
    ($beforeseq, $afterseq)=split /\^/, $enzseq;
    my $rc_cut=$enz->complementary_cut;
    
    # complementary cut is the position on the forward strand
    # correct for reverse strand - I think this is right
    my @all_cuts; # more cuts that we find
    my $more_results=[];
    if ($enz->is_ambiguous) {
          $more_results= $self->_ambig_cuts($beforeseq, $afterseq, $target_seq, $enz);
    } else {
          $more_results= $self->_nonambig_cuts($beforeseq, $afterseq, $target_seq, $enz);
    }
            
    # deal with is_circular
    my $more_more=[];
    $more_more=$self->_circular($beforeseq, $afterseq, $enz) if ($self->{'_seq'}->is_circular);
    push (@all_cuts, @$more_more, @$more_results);
    return \@all_cuts;
} 

=head2 _decon_enzyme

 Title     : _decon_enzyme
 Function  : An internal method to deconvolute and check the enzyme
 Returns   : The sequence before the cut and the sequence after the cut
 Arguments : A Bio::Restriction::Enzyme object

=cut

sub _decon_enzyme {
    my ($self, $enz)=@_;
    # get the cut site. This will have a ^ and N's in it.
    my $site=$enz->site;
    # split it into the two fragments for the sequence before and after.

    # BIG assumption - if the enzyme doesn't have a cut site, we're just guessing.
    # we should really throw a warning here
    unless ($site =~ /\^/) {
        $site = "^$site";
        $self->warn("No cut site is known for $enz. \n".
                    "Therefore we are guessing that it is $site\n")
            if $self->verbose;
    }

    # the default values just stop an error from an undefined
    # string. But they don't affect the split.
    my ($beforeseq, $afterseq)= ('.', '.');
    # now perlify the sequences into a regexp
    ($beforeseq, $afterseq)=split /\^/, $self->_expanded_string($site);
    return ($beforeseq, $afterseq);
}


=head2 _decon_cuts

 Title     : _decon_cuts
 Function  : An internal method to deconvolute the cut sites
 Returns   : The deconvoluted list
 Arguments : A list of cuts sites

=cut

sub _decon_cuts {
    my $self=shift;
    # there is a problem with circular that it can find a site that we have already found
    # but actually, the solution is simple. Because at the moment _cut_positions are the
    # actual locations along the sequence, they can't be redundant

    my %unique_cuts=map {$_=>1} @_;
    my @cuts;
    if (scalar keys %unique_cuts) {
	@cuts=sort {$a <=> $b} keys %unique_cuts;
    }
    return \@cuts;
}

=head2 _ambig_cuts

 Title     : _ambig_cuts
 Function  : An internal method to localize the cuts in the sequence
 Returns   : A reference to an array of cut positions
 Arguments : The separated enzyme site, the target sequence, and the enzyme object
 Comments  : This is a slow implementation but works for ambiguous sequences.
             Whenever possible, _nonambig_cuts should be used as it is a lot faster.

=cut

sub _ambig_cuts {
    my ($self, $beforeseq, $afterseq, $target_seq, $enz) = @_;
    
    # cut the sequence. This is done with split so we can use
    # regexp. Index would be faster!
    
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

    # I don't think we want the last fragment in. It is messing up the _circular
    # part of things. So I deleted this part of the code :)

    } else {
            # if we don't cut, leave the array empty
            #@re_frags=@cuts;
	    return [];
    } # the sequence was not cut.

    # now @re_frags has the fragments of all the sequences
    # but some people want to have this return the lengths
    # of the fragments.
    
    # in theory the actual cut sites should be the length
    # of the fragments in @re_frags
    
    # note, that now this is the only data that we are saving. We
    # will have to go back add regenerate re_frags. The reason is
    # that we can use this in _circular easier

    my @cut_positions = map {length($_)} @re_frags;

    # the cut positions are right now the lengths of the sequence, but
    # we need to add them all onto each other

    for (my $i=1; $i<=$#cut_positions; $i++) {
     $cut_positions[$i]+=$cut_positions[$i-1];
    }
   
    # in one of those oddities in life, 2 fragments mean an enzyme cut once
    # so $#re_frags is the number of cuts
    return \@cut_positions;
}


=head2 _nonambig_cuts

 Title     : _nonambig_cuts
 Function  : Figures out which enzymes we know about and cuts the sequence.
 Returns   : Nothing.
 Arguments : The separated enzyme site, the target sequence, and the enzyme object
 Comments  : An internal method. This will figure out where the sequence 
             should be cut, and provide the appropriate results.
	     This is a much faster implementation because it doesn't use a regexp,
	     but it can not deal with ambiguous sequences

=cut

sub _nonambig_cuts {
    my ($self, $beforeseq, $afterseq, $target_seq, $enz) = @_;
    
    my $index_posn=index($target_seq, $beforeseq.$afterseq);
    return [] if ($index_posn == -1); # there is no match to the sequence
    
    # there is at least one cut site
    my @cuts;
    while ($index_posn > -1) {
	  push (@cuts, $index_posn+length($beforeseq));
	  $index_posn=index($target_seq, $beforeseq.$afterseq, $index_posn+1);
    }
 
    return \@cuts;
}


=head2 _mulitple_cuts

 Title     : _multiple_cuts
 Function  : Figures out multiple digests
 Returns   : An array of the cut sites for multiply digested DNA
 Arguments : A Bio::Restriction::EnzymeCollection object
 Comments  : Double digests is one subset of this, but you can use
             as many enzymes as you want.

=cut

sub _multiple_cuts {
    my ($self, $ec)=@_;
    $self->cut unless $self->{'_cut'};
   
    # now that we are using positions rather than fragments
    # this is really easy
    my @cuts;
    foreach my $enz ($ec->each_enzyme) { 
       push @cuts, @{$self->{'_cut_positions'}->{$enz->name}} if defined $self->{'_cut_positions'}->{$enz->name};
    }
    @{$self->{'_cut_positions'}->{'multiple_digest'}}=sort {$a <=> $b} @cuts;

    my $number_of_cuts;
    
    $number_of_cuts=scalar @{$self->{'_cut_positions'}->{'multiple_digest'}};
    $self->{_number_of_cuts_by_enzyme}->{'multiple_digest'}=$number_of_cuts;
    push (@{$self->{_number_of_cuts_by_cuts}->{$number_of_cuts}}, 'multiple_digest');
    if ($number_of_cuts > $self->{maximum_cuts}) {
        $self->{maximum_cuts}=$number_of_cuts;
    }
}


=head2 _circular

 Title     : _circular
 Function  : Deals with circular sequences
 Returns   : Nothing.
 Arguments : None.
 Comments  : Help?!

There are two problems with circular sequences.

1. When you cut a sequence and rejoin fragments you could generate new cut sites.
2. There could be a cut site at the end of the sequence.

I think these may be the same problem, and so we're working on #2 first!

=cut

sub _circular {
    my ($self, $beforeseq, $afterseq, $enz) = @_;
    my $target_seq=uc $self->{'_seq'}->seq; # I have been burned on this before :)
   
    # the approach I am taking is to find out the longest enzyme in the collection
    # (I'll have to add a new function in enzyme collection for this)
    # and then add more than that sequence from the end of the sequence to the start
    # of the sequence, and map the new cut sites for each of the enzymes.

    # The cut sites that we are interested in must be within the length of the 
    # enzyme sequence from the start or the end.

    my $longest_enz=$self->{'_enzymes'}->longest_cutter;
    my $longest_cut=$longest_enz->recognition_length;
    # this is an error that I don't want to deal with at the moment
    $self->throw("Crap. The longest recognition site ($longest_cut) is longer than the".
      " length of the sequence") if ($longest_cut > $self->{'_seq'}->length);
   
   # newseq is just the last part of the sequence and the first part of the sequence
   # we don't want to go through and check the whole sequence again
   
   my ($first, $last)=(substr($target_seq, 0, $longest_cut),substr($target_seq, -$longest_cut));
   my $newseq=$last.$first;
   
   # now find the cut sites
   # if the enzyme is ambiguous we need to use a regexp to find the cut site
   # otherwise we can use index (much faster)
   my $cut_positions;
   if ($enz->is_ambiguous) {
      $cut_positions= $self->_ambig_cuts($beforeseq, $afterseq, $newseq, $enz);
   } else {
      $cut_positions=$self->_nonambig_cuts($beforeseq, $afterseq, $newseq, $enz);
   }
	
   
   return [] if (!$cut_positions); # the enzyme doesn't cut in the new fragment - likely to be default
   

   # now we are going to add things to _cut_positions
   # in this shema it doesn't matter if the site is there twice - 
   # we will take care of that later. Because we are using position
   # rather than frag or anything else, we can just
   # remove duplicates.
   my @circ_cuts;
   foreach my $cut (@$cut_positions) {
    if ($cut == length($last)) {
     # the cut is actually at position 0, but we're going to call this the
     # length of the sequence so we don't confuse no cuts with a 0 cut
     push (@circ_cuts, $self->{'_seq'}->length);
    }
    elsif ($cut < length($last)) {
     # the cut is before the end of the sequence
     # there is VERY likely to be an off by one error here
     push (@circ_cuts, $self->{'_seq'}->length - (length($last) - $cut));
    }
    else {
     # the cut is at the start of the sequence (position >=1)
     # there is VERY likely to be an off by one error here
     # note, we put this at the beginning of the array rather than the end!
     unshift (@circ_cuts, $cut-length($last));
    }
   }
   return \@circ_cuts;
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
