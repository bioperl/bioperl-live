# module Bio::PopGen::HtSNP.pm
# cared by Pedro M. Gomez-Fabre <pgf18872-at-gsk-dot-com>
#
#

=head1 NAME

Bio::PopGen::HtSNP.pm- Select htSNP from a haplotype set

=head1 SYNOPSIS

    use Bio::PopGen::HtSNP;

    my $obj = Bio::PopGen::HtSNP->new($hap,$snp,$pop);

=head1 DESCRIPTION

Select the minimal set of SNP that contains the full information about
the haplotype without redundancies.

Take as input the followin values:

=over 4

=item - the haplotype block (array of array).

=item - the snp id (array).

=item - family information and frequency (array of array).

=back

The final haplotype is generated in a numerical format and the SNP's
sets can be retrieve from the module.

B<considerations:>


- If you force to include a family with indetermination, the SNP's
with indetermination will be removed from the analysis, so consider
before to place your data set what do you really want to do.

- If two families have the same information (identical haplotype), one
of them will be removed and the removed files will be stored classify
as removed.

- Only are accepted for calculation A, C, G, T and - (as deletion) and
their combinations. Any other value as n or ? will be considered as
degenerations due to lack of information.

=head2 RATIONALE

On a haplotype set is expected that some of the SNP and their
variations contribute in the same way to the haplotype. Eliminating
redundancies will produce a minimal set of SNP's that can be used as
input for a taging selection process. On the process SNP's with the
same variation are clustered on the same group.

The idea is that because the tagging haplotype process is
exponential. All redundant information we could eliminate on the
tagging process will help to find a quick result.

=head2 CONSTRUCTORS

  my $obj = Bio::PopGen::HtSNP->new
    (-haplotype_block => \@haplotype_patterns,
     -snp_ids         => \@snp_ids,
     -pattern_freq    => \@pattern_name_and_freq);

where  $hap, $snp and $pop are in the format:

  my $hap = [
             'acgt',
             'agtc',
             'cgtc'
            ];                     # haplotype patterns' id

  my $snp = [qw/s1 s2 s3 s4/];     # snps' Id's

  my $pop = [
             [qw/ uno    0.20/],
             [qw/ dos    0.20/],
             [qw/ tres   0.15/],
            ];                     # haplotype_pattern_id    Frequency

=head2 OBJECT METHODS

    See Below for more detailed summaries.


=head1 DETAILS

=head2 How the process is working with one example

Let's begin with one general example of the code.

Input haplotype:

  acgtcca-t
  cggtagtgc
  cccccgtgc
  cgctcgtgc

The first thing to to is to B<split the haplotype> into characters.

  a       c       g       t       c       c       a       -       t
  c       g       g       t       a       g       t       g       c
  c       c       c       c       c       g       t       g       c
  c       g       c       t       c       g       t       g       c

Now we have to B<convert> the haplotype to B<Upercase>. This
will produce the same SNP if we have input a or A.

  A       C       G       T       C       C       A       -       T
  C       G       G       T       A       G       T       G       C
  C       C       C       C       C       G       T       G       C
  C       G       C       T       C       G       T       G       C

The program admit as values any combination of ACTG and - (deletions).
The haplotype is B<converted to number>, considering the first variation
as zero and the alternate value as 1 (see expanded description below).

  0       0       0       0       0       0       0       0       0
  1       1       0       0       1       1       1       1       1
  1       0       1       1       0       1       1       1       1
  1       1       1       0       0       1       1       1       1

Once we have the haplotype converted to numbers we have to generate the
snp type information for the haplotype.


B<SNP code = SUM ( value * multiplicity ^ position );>

    where:
      SUM is the sum of the values for the SNP
      value is the SNP number code (0 [generally for the mayor allele],
                                    1 [for the minor allele].
      position is the position on the block.

For this example the code is:

  0       0       0       0       0       0       0       0       0
  1       1       0       0       1       1       1       1       1
  1       0       1       1       0       1       1       1       1
  1       1       1       0       0       1       1       1       1
 ------------------------------------------------------------------
  14      10      12      4       2       14      14      14      14

  14 = 0*2^0 + 1*2^1 + 1*2^2 + 1*2^3
  12 = 0*2^0 + 1*2^1 + 0*2^2 + 1*2^3
  ....

Once we have the families classify. We will B<take> just the SNP's B<not
redundant>.

  14      10      12      4       2

This information will be B<passed to the tag module> is you want to tag
the htSNP.

Whatever it happens to one SNPs of a class will happen to a SNP of
the same class. Therefore you don't need to scan redundancies

=head2 Working with fuzzy data.

This module is designed to work with fuzzy data. As the source of the
haplotype is diverse. The program assume that some haplotypes can be
generated using different values. If there is any indetermination (? or n)
or any other degenerated value or invalid. The program will take away
This SNP and will leave that for a further analysis.

On a complex situation:

  a       c       g       t       ?       c       a       c       t
  a       c       g       t       ?       c       a       -       t
  c       g       ?       t       a       g       ?       g       c
  c       a       c       t       c       g       t       g       c
  c       g       c       t       c       g       t       g       c
  c       g       g       t       a       g       ?       g       c
  a       c       ?       t       ?       c       a       c       t

On this haplotype everything is happening. We have a multialelic variance.
We have indeterminations. We have deletions and we have even one SNP
which is not a real SNP.

The buiding process will be the same on this situation.

Convert the haplotype to uppercase.

  A       C       G       T       ?       C       A       C       T
  A       C       G       T       ?       C       A       -       T
  C       G       ?       T       A       G       ?       G       C
  C       A       C       T       C       G       T       G       C
  C       G       C       T       C       G       T       G       C
  C       G       G       T       A       G       ?       G       C
  A       C       ?       T       ?       C       A       C       T

All columns that present indeterminations will be removed from the analysis
on this Step.

hapotype after remove columns:

  A       C       T       C       C       T
  A       C       T       C       -       T
  C       G       T       G       G       C
  C       A       T       G       G       C
  C       G       T       G       G       C
  C       G       T       G       G       C
  A       C       T       C       C       T

All changes made on the haplotype matrix, will be also made on the SNP list.

  snp_id_1 snp_id_2 snp_id_4 snp_id_6 snp_id_8 snp_id_9

now the SNP that is not one SNP will be removed from the analysis.
SNP with Id snp_id_4 (the one with all T's).


because of the removing. Some of the families will become the same and will
be clustered. A posteriori analysis will diference these families.
but because of the indetermination can not be distinguish.

  A       C       C       C       T
  A       C       C       -       T
  C       G       G       G       C
  C       A       G       G       C
  C       G       G       G       C
  C       G       G       G       C
  A       C       C       C       T

The result of the mergering will go like:

  A       C       C       C       T
  A       C       C       -       T
  C       G       G       G       C
  C       A       G       G       C

Once again the changes made on the families and we merge the frequency (I<to be
implemented>)

Before to convert the haplotype into numbers we consider how many variations
we have on the set. On this case the variations are 3.

The control code will use on this situation base three as mutiplicity

  0       0       0       0       0
  0       0       0       1       0
  1       1       1       2       1
  1       2       1       2       1
 -----------------------------------
  36      63      36      75      36

And the minimal set for this combination is

  0       0       0
  0       0       1
  1       1       2
  1       2       2

B<NOTE:> this second example is a remote example an on normal conditions. This
conditions makes no sense, but as the haplotypes, can come from many sources
we have to be ready for all kind of combinations.


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

=head1 AUTHOR - Pedro M. Gomez-Fabre

Email pgf18872-at-gsk-dot-com


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::PopGen::HtSNP;
use Data::Dumper;
use Storable qw(dclone);

use vars qw ();
use strict;


use base qw(Bio::Root::Root);

my $USAGE = 'Usage:

    Bio::PopGen::HtSNP->new(-haplotype_block -ids -pattern_freq)

';

=head2 new

 Title   : new
 Function: constructor of the class.
 Usage   : $obj-> Bio::PopGen::HtSNP->new(-haplotype_block
                                          -snp_ids
                                          -pattern_freq)
 Returns : self hash
 Args    : input haplotype (array of array)
           snp_ids         (array)
           pop_freq        (array of array)
 Status  : public

=cut

sub new {
    my($class, @args) = @_;

    my $self = $class->SUPER::new(@args);
    my ($haplotype_block,
        $snp_ids,
        $pattern_freq    ) = $self->_rearrange([qw(HAPLOTYPE_BLOCK 
                                                   SNP_IDS
                                                   PATTERN_FREQ)],@args);

    if ($haplotype_block){
        $self->haplotype_block($haplotype_block);
    }
    else{
        $self->throw("Haplotype block has not been defined.
                      \n$USAGE");
    }
    if ($snp_ids){
        $self->snp_ids($snp_ids);
    }
    else{
        $self->throw("Array with ids has not been defined.
                      \n$USAGE");
    }
    if ($pattern_freq){
        $self->pattern_freq($pattern_freq);
    }
    else{
        $self->throw("Array with pattern id and frequency has not been defined.
                      \n$USAGE");
    }

    # if the input values are not well formed complained and exit.
    _check_input($self);

    _do_it($self);

    return $self;
}

=head2 haplotype_block 

 Title   : haplotype_block 
 Usage   : my $haplotype_block = $HtSNP->haplotype_block();
 Function: Get the haplotype block for a haplotype tagging selection
 Returns : reference of array 
 Args    : reference of array with haplotype pattern 


=cut

sub haplotype_block{
    my ($self) =shift;
    return $self->{'_haplotype_block'} = shift if @_;
    return $self->{'_haplotype_block'};
}

=head2 snp_ids 

 Title   : snp_ids 
 Usage   : my $snp_ids = $HtSNP->$snp_ids();
 Function: Get the ids for a haplotype tagging selection
 Returns : reference of array
 Args    : reference of array with SNP ids


=cut

sub snp_ids{
    my ($self) =shift;
    return $self->{'_snp_ids'} = shift if @_;
    return $self->{'_snp_ids'};
}


=head2 pattern_freq

 Title   : pattern_freq
 Usage   : my $pattern_freq = $HtSNP->pattern_freq();
 Function: Get the pattern id and frequency  for a haplotype
           tagging selection
 Returns : reference of array
 Args    : reference of array with SNP ids

=cut

sub pattern_freq{
    my ($self) =shift;
    return $self->{'_pattern_freq'} = shift if @_;
    return $self->{'_pattern_freq'};
}

=head2 _check_input

 Title   : _check_input
 Usage   : _check_input($self)
 Function: check for errors on the input
 Returns : self hash
 Args    : self
 Status  : internal

=cut

#------------------------
sub _check_input{
#------------------------

    my $self = shift;

    _haplotype_length_error($self);
    _population_error($self);

}

=head2 _haplotype_length_error

 Title   : _haplotype_length_error
 Usage   : _haplotype_length_error($self)
 Function: check if the haplotype length is the same that the one on the
           SNP id list. If not break and exit
 Returns : self hash
 Args    : self
 Status  : internal

=cut


#------------------------
sub _haplotype_length_error{
#------------------------

    my $self = shift;

    my $input_block = $self->haplotype_block();
    my $snp_ids     = $self->snp_ids();


    #############################
    # define error list
    #############################
    my $different_haplotype_length = 0;

    ##############################
    # get parameters used to find
    # the errors
    ##############################

    my $snp_number         = scalar @$snp_ids;
    my $number_of_families = scalar @$input_block;
    my $h                  = 0; # haplotype position


    ############################
    # haplotype length
    #
    # if the length differs from the number of ids
    ############################

    for ($h=0; $h<$#$input_block+1 ; $h++){
        if (length $input_block->[$h]  != $snp_number){
            $different_haplotype_length = 1;
            last;
        }
    }

    # haploytypes does not have the same length
    if ($different_haplotype_length){
       $self->throw("The number of snp ids is $snp_number and ".
            "the length of the family (". ($h+1) .") [".
             $input_block->[$h]."] is ".
             length $input_block->[$h], "\n");
    }
}

=head2 _population_error


 Title   : _population_error
 Usage   : _population_error($self)
 Function: use input_block and pop_freq test if the number of elements
           match. If doesn't break and quit.
 Returns : self hash
 Args    : self
 Status  : internal

=cut


#------------------------
sub _population_error{
#------------------------

    my $self = shift;

    my $input_block = $self->haplotype_block();
    my $pop_freq    = $self->pattern_freq();

    #############################
    # define error list
    #############################
    my $pop_freq_elements_error    = 0;  # matrix bad formed

    ##############################
    # get parameters used to find
    # the errors
    ##############################
    my $number_of_families = scalar @$input_block;

    my $pf         = 0; # number of elements on population frequency
    my $frequency  = 0; # population frequency
    my $p_f_length = 0;

    # check if the pop_freq array is well formed and if the number
    # of elements fit with the number of families

    #############################
    # check population frequency
    #
    # - population frequency matrix need to be well formed
    # - get the frequency
    # - calculate number of families on pop_freq
    #############################

    for  ($pf=0; $pf<$#$pop_freq+1; $pf++){
        $frequency += $pop_freq->[$pf]->[1];

        if ( scalar @{$pop_freq->[$pf]} !=2){
            $p_f_length = scalar @{$pop_freq->[$pf]};
            $pop_freq_elements_error = 1;
            last;
        }
    }

    ###########################
    ## error processing
    ###########################


    # The frequency shouldn't be greater than 1
    if ($frequency >1) {
        $self->warn("The frequency for this set is $frequency (greater than 1)\n");
    }

    # the haplotype matix is not well formed
    if ($pop_freq_elements_error){
        $self->throw("the frequency matrix is not well formed\n".
             "\nThe number of elements for pattern ".($pf+1)." is ".
             "$p_f_length\n".
             "It should be 2 for pattern \"@{$pop_freq->[$pf]}\"\n".
             "\nFormat should be:\n".
             "haplotype_id\t frequency\n"
            );
    }

    # the size does not fit on pop_freq array
    #  with the one in haplotype (input_block)
    if ($pf != $number_of_families) {
        $self->throw("The number of patterns on frequency array ($pf)\n".
             "does not fit with the number of haplotype patterns on \n". 
             "haplotype array ($number_of_families)\n");
    }
}

=head2 _do_it


 Title   : _do_it
 Usage   : _do_it($self)
 Function: Process the input generating the results.
 Returns : self hash
 Args    : self
 Status  : internal

=cut

#------------------------
sub _do_it{
#------------------------

    my $self = shift;

    # first we are goinf to define here all variables we are going to use
    $self -> {'w_hap'}          = [];
    $self -> {'w_pop_freq'}     = dclone ( $self ->pattern_freq() );
    $self -> {'deg_pattern'}    = {};
    $self -> {'snp_type'}       = {};  # type of snp on the set. see below
    $self -> {'alleles_number'} = 0;   # number of variations (biallelic,...)
    $self -> {'snp_type_code'}  = [];
    $self -> {'ht_type'}        = [];  # store the snp type used on the htSet
    $self -> {'split_hap'}      = [];
    $self -> {'snp_and_code'}   = [];


    # we classify the SNP under snp_type
    $self->{snp_type}->{useful_snp} = dclone ( $self ->snp_ids() );
    $self->{snp_type}->{deg_snp}    = []; # deg snp
    $self->{snp_type}->{silent_snp} = []; # not a real snp

    # split the haplotype
    _split_haplo ($self);

    # first we convert to upper case the haplotype
    # to make A the same as a for comparison
    _to_upper_case( $self -> {w_hap} );

    #######################################################
    # check if any SNP has indetermination. If any SNP has
    # indetermination this value will be removed.
    #######################################################
    _remove_deg ( $self );

    #######################################################
    # depending of the families you use some SNPs can be
    # silent. This silent SNP's are not used on the
    # creation of tags and has to be skipped from the
    # analysis.
    #######################################################
    _rem_silent_snp ( $self );

    #######################################################
    # for the remaining SNP's we have to check if two
    # families have the same value. If this is true, the families
    # will produce the same result and therefore we will not find
    # any pattern. So, the redundant families need to be take
    # away from the analysis. But also considered for a further
    # run.
    #
    # When we talk about a normal haplotype blocks this situation
    # makes no sense but if we remove one of the snp because the
    # degeneration two families can became the same.
    # these families may be analised on a second round
    #######################################################

    _find_deg_pattern ( $self );

    #################################################################
    # if the pattern list length is different to the lenght of the w_hap
    # we can tell that tow columns have been considered as the same one
    # and therefore we have to start to remove the values.
    # remove all columns with degeneration
    #
    # For this calculation we don't use the pattern frequency.
    # All patterns are the same, This selection makes
    # sense when you have different frequency.
    #
    # Note: on this version we don't classify the haplotype by frequency
    # but if you need to do it. This is the place to do it!!!!
    #
    # In reality you don't need to sort the values because you will remove
    # the values according to their values.
    #
    # But as comes from a hash, the order could be different and as a
    # consequence the code generate on every run of the same set could
    # differ. That is not important. In fact, does not matter but could
    # confuse people.
    #################################################################

    my @tmp =sort { $a <=> $b}
         keys %{$self -> {deg_pattern}}; # just count the families

    # if the size of the list is different to the size of the degenerated
    # family. There is degeneration. And the redundancies will be
    # removed.
    if($#tmp != $#{$self -> { w_hap } } ){
        _keep_these_patterns($self->{w_hap}, \@tmp);
        _keep_these_patterns($self->{w_pop_freq}, \@tmp);
    }

    #################################################################
    # the steps made before about removing snp and cluster families
    # are just needed pre-process the haplotype before.
    #
    # Now is when the fun starts.
    #
    #
    # once we have the this minimal matrix, we have to calculate the
    # max multipliticy for the values. The max number of alleles found
    # on the set. A normal haplotype is biallelic but we can not
    # reject multiple variations.
    ##################################################################

    _alleles_number ( $self );

    ##################################################################
    # Now we have to convert the haplotype into number
    #
    # A       C       C       -       T
    # C       A       G       G       C
    # A       C       C       C       T
    # C       G       G       G       C
    #
    # one haplotype like this transformed into number produce this result
    #
    # 0       0       0       0       0
    # 1       1       1       1       1
    # 0       0       0       2       0
    # 1       2       1       1       1
    #
    ##################################################################

    _convert_to_numbers( $self );

    ###################################################################
    # The next step is to calculate the type of the SNP.
    # This process is made based on the position of the SNP, the value
    # and its multiplicity.
    ###################################################################

    _snp_type_code( $self );

    ###################################################################
    # now we have all information we need to calculate the haplotype
    # tagging SNP htSNP
    ###################################################################

    _htSNP( $self );

    ###################################################################
    # patch:
    #
    # all SNP have a code. but if the SNP is not used this code must
    # be zero in case of silent SNP. This looks not to informative
    # because all the information is already there. But this method
    # compile the full set.
    ###################################################################

    _snp_and_code_summary( $self );
}

=head2 input_block

 Title   : input_block
 Usage   : $obj->input_block()
 Function: returns input block
 Returns : reference to array of array
 Args    : none
 Status  : public

=cut

#------------------------
sub input_block{
#------------------------

    my $self = shift;
    return $self -> {input_block};
}

=head2 hap_length

 Title   : hap_length
 Usage   : $obj->hap_length()
 Function: get numbers of SNP on the haplotype
 Returns : scalar
 Args    : none
 Status  : public

=cut

#------------------------
sub hap_length{
#------------------------

    my $self = shift;
    return scalar @{$self -> {'_snp_ids'}};
}


=head2 pop_freq

 Title   : pop_freq
 Usage   : $obj->pop_freq()
 Function: returns population frequency
 Returns : reference to array
 Args    : none
 Status  : public

=cut

#------------------------
sub pop_freq{
#------------------------

    my $self = shift;
    return $self -> {pop_freq}
}


=head2 deg_snp


 Title   : deg_snp
 Usage   : $obj->deg_snp()
 Function: returns snp_removes due to indetermination on their values
 Returns : reference to array
 Args    : none
 Status  : public

=cut

#------------------------
sub deg_snp{
#------------------------
    my $self = shift;
    return $self -> {snp_type} ->{deg_snp};
}


=head2 snp_type


 Title   : snp_type
 Usage   : $obj->snp_type()
 Function: returns hash with SNP type
 Returns : reference to hash
 Args    : none
 Status  : public

=cut

#------------------------
sub snp_type{
#------------------------
    my $self = shift;
    return $self -> {snp_type};
}


=head2 silent_snp


 Title   : silent_snp
 Usage   : $obj->silent_snp()
 Function: some SNP's are silent (not contibuting to the haplotype)
           and are not considering for this analysis
 Returns : reference to a array
 Args    : none
 Status  : public

=cut

#------------------------
sub silent_snp{
#------------------------
    my $self = shift;
    return $self -> {snp_type} ->{silent_snp};
}


=head2 useful_snp


 Title   : useful_snp
 Usage   : $obj->useful_snp()
 Function: returns list of SNP's that are can be used as htSNP. Some
           of them can produce the same information. But this is
           not considered here.
 Returns : reference to a array
 Args    : none
 Status  : public

=cut

#------------------------
sub useful_snp{
#------------------------
    my $self = shift;
    return $self -> {snp_type} ->{useful_snp};
}


=head2 ht_type


 Title   : ht_type
 Usage   : $obj->ht_type()
 Function: every useful SNP has a numeric code dependending of its
           value and position. For a better description see
           description of the module.
 Returns : reference to a array
 Args    : none
 Status  : public

=cut

#------------------------
sub ht_type{
#------------------------
    my $self = shift;
    return $self -> {ht_type};
}
=head2 ht_set


 Title   : ht_set
 Usage   : $obj->ht_set()
 Function: returns the minimal haplotype in numerical format. This
           haplotype contains the maximal information about the
           haplotype variations but with no redundancies. It's the
           minimal set that describes the haplotype.
 Returns : reference to an array of arrays
 Args    : none
 Status  : public

=cut

#------------------------
sub ht_set{
#------------------------
    my $self = shift;
    return $self -> {w_hap};
}

=head2 snp_type_code


 Title   : snp_type_code
 Usage   : $obj->snp_type_code()
 Function: returns the numeric code of the SNPs that need to be
           tagged that correspond to the SNP's considered in ht_set.
 Returns : reference to an array
 Args    : none
 Status  : public

=cut

#------------------------
sub snp_type_code{
#------------------------
    my $self = shift;
    return $self -> {snp_type_code};
}

=head2 snp_and_code


 Title   : snp_and_code
 Usage   : $obj->snp_and_code()
 Function: Returns the full list of SNP's and the code associate to
           them. If the SNP belongs to the group useful_snp it keep
           this code. If the SNP is silent the code is 0. And if the
           SNP is degenerated the code is -1.
 Returns : reference to an array of array
 Args    : none
 Status  : public

=cut

#------------------------
sub snp_and_code{
#------------------------
    my $self = shift;
    return $self -> {'snp_and_code'};
}

=head2 deg_pattern


 Title   : deg_pattern
 Usage   : $obj->deg_pattern()
 Function: Returns the a list with the degenerated haplotype.
           Sometimes due to degeneration some haplotypes looks
           the same and if we don't remove them it won't find
           any tag.
 Returns : reference to a hash of array
 Args    : none
 Status  : public

=cut

#------------------------
sub deg_pattern{
#------------------------
    my $self = shift;

    return $self -> {'deg_pattern'};
}

=head2 split_hap


 Title   : split_hap
 Usage   : $obj->split_hap()
 Function: simple representation of the haplotype base by base
           Same information that input haplotype but base based.
 Returns : reference to an array of array
 Args    : none
 Status  : public

=cut

#------------------------
sub split_hap{
#------------------------
    my $self = shift;
    return $self -> {'split_hap'};
}

=head2 _split_haplo

 Title   : _split_haplo
 Usage   : _split_haplo($self)
 Function: Take a haplotype and split it into bases
 Returns : self
 Args    : none
 Status  : internal

=cut

#------------------------
sub _split_haplo {
#------------------------
    my $self = shift;

    my $in  = $self ->{'_haplotype_block'};
    my $out = $self ->{'w_hap'};

    # split every haplotype and store the result into $out
    foreach (@$in){
        push @$out, [split (//,$_)];
    }

    $self -> {'split_hap'} = dclone ($out);
}

# internal method to convert the haplotype to uppercase


=head2 _to_upper_case


 Title   : _to_upper_case
 Usage   : _to_upper_case()
 Function: make SNP or in-dels Upper case
 Returns : self
 Args    : an AoA ref
 Status  : private

=cut

#------------------------
sub _to_upper_case {
#------------------------
    my ($arr) =@_;

    foreach my $aref (@$arr){
        foreach my $value (@{$aref} ){
            $value = uc $value;
        }
    }
}


=head2 _remove_deg


 Title   : _remove_deg
 Usage   : _remove_deg()
 Function: when have a indetermination or strange value this SNP
           is removed
 Returns : haplotype family set and degeneration list
 Args    : ref to an AoA and a ref to an array
 Status  : internal

=cut

#------------------------
sub _remove_deg {
#------------------------
    my $self = shift;

    my $hap         = $self->{w_hap};
    my $snp         = $self->{snp_type}->{useful_snp};
    my $deg_snp     = $self->{snp_type}->{deg_snp};

    my $rem = [];  # take the position of the array to be removed

    # first we work on the columns we have void values
    $rem = _find_indet($hap,$rem);  # find degenerated columns

    if (@$rem){

        # remove column on haplotype
        _remove_col($hap,$rem); # remove list

        # now remove the values from SNP id
        _remove_snp_id($snp,$deg_snp,$rem); # remove list
    }
}


=head2 _rem_silent_snp


 Title   : _rem_silent_snp
 Usage   : _rem_silent_snp()
 Function: there is the remote possibilty that one SNP won't be a
           real SNP on this situation we have to remove this SNP,
           otherwise the program won't find any tag
 Returns : nonthing
 Args    : ref to an AoA and a ref to an array
 Status  : internal

=cut

#------------------------
sub _rem_silent_snp {
#------------------------
    my $self = shift;

    my $hap         = $self->{w_hap};
    my $snp         = $self->{snp_type}->{useful_snp};
    my $silent_snp  = $self->{snp_type}->{silent_snp};

    my $rem = [];   # store the positions to be removed

    #find columns with no variation on the SNP, Real snp?
    $rem = _find_silent_snps($hap);

    if (@$rem){

        # remove column on haplotype
        _remove_col($hap,$rem);

        # remove the values from SNP id
        _remove_snp_id($snp,$silent_snp,$rem);
    }
}


=head2 _find_silent_snps


 Title   : _find_silent_snps
 Usage   :
 Function: list of snps that are not SNPs. All values for that
           SNPs on the set is the same one. Look stupid but can
           happend and if this happend you will not find any tag
 Returns : nothing
 Args    :
 Status  :

=cut

#------------------------
sub _find_silent_snps{
#------------------------
    my ($arr)=@_;

    my $list =[]; # no snp list;

    # determine the number of snp by the length of the first row.
    # we assume that the matrix is squared.
    my $colsn= @{$arr->[0]};

    for (my $i=0;$i<$colsn;$i++){
        my $different =0;  # check degeneration

        for my $r (1..$#$arr){
            if($arr->[0][$i] ne $arr->[$r][$i]){
                $different =1;
                last;
            }
        }

        if(!$different){
            push (@$list, $i);
        }
    }

    return $list;
}


=head2 _find_indet


 Title   : _find_indet
 Usage   :
 Function: find column (SNP) with invalid or degenerated values
           and store this values into the second parameter supplied.
 Returns : nothing
 Args    : ref to AoA and ref to an array
 Status  : internal

=cut

#------------------------
sub _find_indet{
#------------------------
    my ($arr, $list)=@_;

    foreach my $i(0..$#$arr){
        foreach my $j(0..$#{$arr->[$i]}){
            unless ($arr->[$i][$j] =~ /[ACTG-]/){
                if ($#$list<0){
                    push(@$list,$j);
                }
                else{
                    my $found =0;   # check if already exist the value
                    foreach my $k(0..$#$list){
                        $found =1 if ($list->[$k] eq $j);
                        last if ($found);
                    }
                    if(!$found){
                        push(@$list,$j);
                    }
                }
            }
        }
    }

    @$list = sort { $a <=> $b} @$list;

    return $list;
}

=head2 _remove_col

 Title   : _remove_col
 Usage   :
 Function: remove columns contained on the second array from
           the first arr
 Returns : nothing
 Args    : array of array reference and array reference
 Status  : internal

=cut

#------------------------
sub _remove_col{
#------------------------
    my ($arr,$rem)=@_;

    foreach my $col (reverse @$rem){
        splice @$_, $col, 1 for @$arr;
    }
}


=head2 _remove_snp_id

 Title   : _remove_snp_id
 Usage   :
 Function: remove columns contained on the second array from
           the first arr
 Returns : nothing
 Args    : array of array reference and array reference
 Status  : internal

=cut

#------------------------
sub _remove_snp_id{
#------------------------
    my ($arr,$removed,$rem_list)=@_;

    push @$removed, splice @$arr, $_, 1 foreach reverse @$rem_list;
}


=head2 _find_deg_pattern

 Title   : _find_deg_pattern
 Usage   :
 Function: create a list with the degenerated patterns
 Returns : @array
 Args    : a ref to AoA
 Status  : public

=cut

#------------------------
sub _find_deg_pattern{
#------------------------
    my $self  = shift;

    my $arr   = $self ->{w_hap};          # the working haplotype
    my $list  = $self ->{'deg_pattern'};  # degenerated patterns 

    # we have to check all elements
    foreach my $i(0..$#$arr){
        # is the element has not been used create a key
        unless  ( _is_on_hash ($list,\$i) ) {
            $list->{$i}=[$i];
        };

        foreach my $j($i+1..$#$arr){
            my $comp = compare_arrays($arr->[$i],$arr->[$j]);

            if($comp){
                # as we have no elements we push this into the list
                # check for the first element
                my $key = _key_for_value($list,\$i);

                push (@{$list->{$key}},$j);

                last;
            }
        }
    }

}

#------------------------
sub _key_for_value{
#------------------------
    my($hash,$value)=@_;

    foreach my $key (keys %$hash){
        if( _is_there(\@{$hash->{$key}},$value)){
            return $key;
        }
    }
}

#------------------------
sub _is_on_hash{
#------------------------
    my($hash,$value)=@_;

    foreach my $key (keys %$hash){
        if( _is_there(\@{$hash->{$key}},$value)){
            return 1;
        }
    }
}

#------------------------
sub _is_there{
#------------------------

    my($arr,$value)=@_;

    foreach my $el (@$arr){
        if ($el eq $$value){
            return 1;
        }
    }
}


=head2 _keep_these_patterns


 Title   : _keep_these_patterns
 Usage   :
 Function: this is a basic approach, take a LoL and a list,
           keep just the columns included on the list
 Returns : nothing
 Args    : an AoA and an array
 Status  : public

=cut

#------------------------
sub _keep_these_patterns{
#------------------------
    my ($arr,$list)=@_;

    # by now we just take one of the repetitions but you can weight
    # the values by frequency

    my @outValues=();

    foreach my $k (@$list){
        push @outValues, $arr->[$k];
    }

    #make arr to hold the new values
    @$arr= @{dclone(\@outValues)};

}


=head2 compare_arrays


 Title   : compare_arrays
 Usage   :
 Function: take two arrays and compare their values
 Returns : 1 if the two values are the same
           0 if the values are different
 Args    : an AoA and an array
 Status  : public

=cut

#------------------------
sub compare_arrays {
#------------------------
    my ($first, $second) = @_;
    return 0 unless @$first == @$second;
    for (my $i = 0; $i < @$first; $i++) {
        return 0 if $first->[$i] ne $second->[$i];
    }
    return 1;
}


=head2 _convert_to_numbers


 Title   : _convert_to_numbers
 Usage   : _convert_to_numbers()
 Function: tranform the haplotype into numbers. before to do that
           we have to consider the variation on the set.
 Returns : nonthing
 Args    : ref to an AoA and a ref to an array
 Status  : internal

=cut

#------------------------
sub _convert_to_numbers{
#------------------------
    my $self = shift;

    my $hap_ref = $self->{w_hap};
    my $mm      = $self->{alleles_number};

    # the first element is considered as zero. The first modification
    # is consider as one and so on.

    my $length = @{ @$hap_ref[0]};    #length of the haplotype

    for (my $c = 0; $c<$length;$c++){

        my @al=();

        for my $r (0..$#$hap_ref){

            push @al,$hap_ref->[$r][$c]
                unless _is_there(\@al,\$hap_ref->[$r][$c]);

            $hap_ref->[$r][$c] = get_position(\@al,\$hap_ref->[$r][$c]);
        }
    }
}


=head2 _snp_type_code


 Title   : _snp_type_code
 Usage   :
 Function:
           we have to create the snp type code for each version.
           The way the snp type is created is the following:

           we take the number value for every SNP and do the
           following calculation

           let be a SNP set as follow:

           0    0
           1    1
           1    2

           and multiplicity 3
           on this case the situation is:

           sum (value * multiplicity ^ position) for each SNP

           0 * 3 ^ 0 + 1 * 3 ^ 1 + 1 * 3 ^ 2 = 12
           0 * 3 ^ 0 + 1 * 3 ^ 1 + 2 * 3 ^ 2 = 21
 Returns : nothing
 Args    : $self
 Status  : private

=cut

#------------------------
sub _snp_type_code{
#------------------------
    my $self = shift;

    my $hap = $self->{w_hap};
    my $arr = $self->{snp_type_code};
    my $al  = $self->{alleles_number};

    my $length = @{ $hap->[0]};    #length of the haplotype

    for (my $c=0; $c<$length; $c++){
        for my $r (0..$#$hap){
            $arr->[$c] += $hap->[$r][$c] * $al ** $r;
        }
    }
}

#################################################
# return the position of an element in one array
# The element is always present on the array
#################################################

#------------------------
sub get_position{
#------------------------

    my($array, $value)=@_;

    for my $i(0..$#$array) {
        if ($array->[$i] eq $$value){
            return $i;
        }
    }

}


=head2 _alleles_number


 Title   : _alleles_number
 Usage   :
 Function: calculate the max number of alleles for a haplotype and
           if the number. For each SNP the number is stored and the
           max number of alleles for a SNP on the set is returned
 Returns : max number of alleles (a scalar storing a number)
 Args    : ref to AoA
 Status  : public

=cut

#------------------------
sub _alleles_number{
#------------------------

    my $self = shift;

    my $hap_ref = $self ->{w_hap};          # working haplotype

    my $length = @{ @$hap_ref[0]};    # length of the haplotype

    for (my $c = 0; $c<$length;$c++){

        my %alleles=();

        for my $r (0..$#$hap_ref){
            $alleles{ $hap_ref->[$r][$c] } =1; # new key for every new snp
        }

        # if the number of alleles for this column is
        # greater than before set $m value as allele number
        if ($self->{alleles_number} < keys %alleles) {
            $self->{alleles_number} = keys %alleles;
        }
    }
}


=head2 _htSNP


 Title   : _htSNP
 Usage   : _htSNP()
 Function: calculate the minimal set that contains all information of the
           haplotype.
 Returns : nonthing
 Args    : ref to an AoA and a ref to an array
 Status  : internal

=cut

#------------------------
sub _htSNP{
#------------------------
    my $self = shift;

    my $hap           = $self->{'w_hap'};
    my $type          = $self->{'snp_type_code'};
    my $set           = $self->{'ht_type'};
    my $out           = [];     # store the minimal set

    my $nc=0;        # new column for the output values

    # pass for every value of the snp_type_code
    for my $c (0..$#$type){

        my $exist =0;

        # every new value (not present) is pushed into set
        if ( ! _is_there( $set,\$type->[$c] ) ){
            push @$set, $type->[$c];

            $exist =1;

            for my $r(0..$#$hap){
                #save value of the snp for every SNP
                $out->[$r][$nc]= $hap->[$r][$c];
            }
        }

        if ($exist){ $nc++ };
    }

    @$hap = @{dclone $out};
}

=head2 _snp_and_code_summary

 Title   : _snp_and_code_summary
 Usage   : _snp_and_code_summary()
 Function: compile on a list all SNP and the code for each. This
           information can be also obtained combining snp_type and
           snp_type_code but on these results the information about
           the rest of SNP's are not compiled as table.

           0 will be silent SNPs
           -1 are degenerated SNPs
           and the rest of positive values are the code for useful SNP

 Returns : nonthing
 Args    : ref to an AoA and a ref to an array
 Status  : internal

=cut

#------------------------
sub _snp_and_code_summary{
#------------------------
    my $self = shift;

    my $snp_type_code = $self->{'snp_type_code'};
    my $useful_snp    = $self->{'snp_type'}->{'useful_snp'};
    my $silent_snp    = $self->{'snp_type'}->{'silent_snp'};
    my $deg_snp       = $self->{'snp_type'}->{'deg_snp'};
    my $snp_ids       = $self->snp_ids();
    my $snp_and_code  = $self->{'snp_and_code'};

    # walk all SNP's and generate code for each

    # do a practical thing. Consider all snp silent
    foreach my $i (0..$#$snp_ids){

        # assign zero to silent
        my $value=0;

        # active SNPs
        foreach my $j (0..$#$useful_snp){
            if ($snp_ids->[$i] eq $useful_snp->[$j]){
                $value = $snp_type_code->[$j];
                last;
            }
        }

        # assign -1 to degenerated
        foreach my $j (0..$#$deg_snp){
            if ($snp_ids->[$i] eq $deg_snp->[$j]){
                $value = -1;
                last;
            }
        }

        push @$snp_and_code, [$snp_ids->[$i], $value];

    }
}


1;


