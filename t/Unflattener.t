# -*-Perl-*- mode (to keep my emacs happy)
# $Id$

use strict;
use vars qw($DEBUG $TESTCOUNT);
BEGIN {     
    eval { require Test; };
    if( $@ ) {
	use lib 't';
    }
    use Test;
    $TESTCOUNT = 6;
    plan tests => $TESTCOUNT;
}

use Bio::Seq;
use Bio::SeqIO;
use Bio::Root::IO;
use Bio::SeqFeature::Tools::Unflattener;

ok(1);

my $verbosity = -1;   # Set to -1 for release version, so warnings aren't printed

my ($seq, @sfs);
my $unflattener = Bio::SeqFeature::Tools::Unflattener->new;


if (1) {
    my @path = ("t","data","AE003644_Adh-genomic.gb");
    # allow cmd line override
    if (@ARGV) {
	@path = (shift @ARGV);
    }
    $seq = getseq(@path);
    
    ok ($seq->accession_number, 'AE003644');
    my @topsfs = $seq->get_SeqFeatures;
    printf "TOP:%d\n", scalar(@topsfs);
    write_hier(@topsfs);
    
    # UNFLATTEN
    @sfs = $unflattener->unflatten_seq(-seq=>$seq,
                                     -group_tag=>'locus_tag');
    print "\n\nPOST PROCESSING:\n";
    write_hier(@sfs);
    printf "PROCESSED:%d\n", scalar(@sfs);
    ok(@sfs == 21);
}

# now try again, using a custom subroutine to link together features
$seq = getseq("t","data","AE003644_Adh-genomic.gb");
@sfs = $unflattener->unflatten_seq(-seq=>$seq,
                                 -group_tag=>'locus_tag',
                                 -resolver_method=>sub {
                                     my $self = shift;
                                     my ($sf, @candidate_container_sfs) = @_;
                                     if ($sf->has_tag('note')) {
                                         my @notes = $sf->get_tag_values('note');
                                         my @trnames = map {/from transcript\s+(.*)/;
                                                            $1} @notes;
                                         @trnames = grep {$_} @trnames;
                                         my $trname;
                                         if (@trnames == 0) {
                                             $self->throw("UNRESOLVABLE");
                                         }
                                         elsif (@trnames == 1) {
                                             $trname = $trnames[0];
                                         }
                                         else {
                                             $self->throw("AMBIGUOUS: @trnames");
                                         }
                                         my @container_sfs =
                                           grep {
                                               my ($product) =
                                                 $_->has_tag('product') ?
                                                   $_->get_tag_values('product') :
                                                     ('');
                                               $product eq $trname;
                                           } @candidate_container_sfs;
                                         if (@container_sfs == 0) {
                                             $self->throw("UNRESOLVABLE");
                                         }
                                         elsif (@container_sfs == 1) {
                                             # we got it!
                                             return $container_sfs[0];
                                         }
                                         else {
                                             $self->throw("AMBIGUOUS");
                                         }
                                         
                                     }
                                 });
$unflattener->feature_from_splitloc(-seq=>$seq);
print "\n\nPOST PROCESSING:\n";
write_hier(@sfs);
printf "PROCESSED2:%d\n", scalar(@sfs);
ok(@sfs == 21);

# try again; different sequence
# this is an E-Coli seq with no mRNA features;
# we just want to link all features directly with gene

$seq = getseq("t","data","D10483.gbk");

# UNFLATTEN
@sfs = $unflattener->unflatten_seq(-seq=>$seq,
				   -partonomy=>{'*'=>'gene'},
                                );
print "\n\nPOST PROCESSING:\n";
write_hier(@sfs);
printf "PROCESSED:%d\n", scalar(@sfs);
ok(@sfs == 98);


# try again; different sequence

$seq = getseq("t","data","X98338_Adh-mRNA.gb");

# UNFLATTEN
@sfs = $unflattener->unflatten_seq(-seq=>$seq,
                                 -partonomy=>{'*'=>'gene'},
                                );
                                 
print "\n\nPOST PROCESSING:\n";
write_hier(@sfs);
printf "PROCESSED:%d\n", scalar(@sfs);
ok(@sfs == 7);




sub write_hier {
    my @sfs = @_;
    _write_hier(0, @sfs);
}
sub _write_hier {
    my $indent = shift;
    my @sfs = @_;
    foreach my $sf (@sfs) {
        my $label;
        if ($sf->has_tag('product')) {
            ($label) = $sf->get_tag_values('product');
        }
        printf "%s%s $label\n", '  ' x $indent, $sf->primary_tag;
        my @sub_sfs = $sf->sub_SeqFeature;
        _write_hier($indent+1, @sub_sfs);
    }
}

sub getseq {
    my @path = @_;
    my $seqio =
      Bio::SeqIO->new('-file'=> Bio::Root::IO->catfile(
                                                       @path
                                                      ), 
                      '-format' => 'GenBank');
    $seqio->verbose($verbosity);

    my $seq = $seqio->next_seq();
    return $seq;
}
