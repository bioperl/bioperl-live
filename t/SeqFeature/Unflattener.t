# -*-Perl-*- Test Harness script for Bioperl
# $Id$

use strict;

BEGIN { 
    use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 10);
	
	use_ok('Bio::SeqIO');
	use_ok('Bio::SeqFeature::Tools::Unflattener');
}

my $verbosity = test_debug();

my ($seq, @sfs);
my $unflattener = Bio::SeqFeature::Tools::Unflattener->new();

if (1) {
    my @path = ("ribosome-slippage.gb");
    $seq = getseq(@path);
    
    my @topsfs = $seq->get_SeqFeatures;
    if( $verbosity > 0 ) {
	warn sprintf "TOP:%d\n", scalar(@topsfs);
	write_hier(@topsfs);
    }
    
    # UNFLATTEN
    $unflattener->verbose($verbosity);
    @sfs = $unflattener->unflatten_seq(-seq=>$seq,
				       -use_magic=>1);
    if( $verbosity > 0 ) {
	warn "\n\nPOST PROCESSING:\n";
	write_hier(@sfs);
	warn sprintf "PROCESSED:%d\n", scalar(@sfs);
    }
    is(@sfs, 3);
}

if (1) {
    my @path = ("AE003644_Adh-genomic.gb");
    $seq = getseq(@path);
    
    is ($seq->accession_number, 'AE003644');
    my @topsfs = $seq->get_SeqFeatures;
    if( $verbosity > 0 ) {
	warn sprintf "TOP:%d\n", scalar(@topsfs);
	write_hier(@topsfs);
    }
    
    # UNFLATTEN
    $unflattener->verbose($verbosity);
    @sfs = $unflattener->unflatten_seq(-seq=>$seq,
				       -group_tag=>'locus_tag');
    if( $verbosity > 0 ) {
	warn "\n\nPOST PROCESSING:\n";
	write_hier(@sfs);
	warn sprintf "PROCESSED:%d\n", scalar(@sfs);
    }
    is(@sfs, 21);
}

# now try again, using a custom subroutine to link together features
$seq = getseq("AE003644_Adh-genomic.gb");
@sfs = $unflattener->unflatten_seq
    (-seq=>$seq,
     -group_tag=>'locus_tag',
     -resolver_method => 
     sub {
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
		 return ($container_sfs[0]=>0);
	     }
	     else {
		 $self->throw("AMBIGUOUS");
	     }
	 }
     });
$unflattener->feature_from_splitloc(-seq=>$seq);
if( $verbosity > 0 ) {
    warn "\n\nPOST PROCESSING:\n";
    write_hier(@sfs);
    warn sprintf "PROCESSED2:%d\n", scalar(@sfs);
}
is(@sfs, 21);

# try again; different sequence
# this is an E-Coli seq with no mRNA features;
# we just want to link all features directly with gene

$seq = getseq("D10483.gbk");

# UNFLATTEN
@sfs = $unflattener->unflatten_seq(-seq=>$seq,
				   -partonomy=>{'*'=>'gene'},
                                );
if( $verbosity > 0 ) {
    warn "\n\nPOST PROCESSING:\n";
    write_hier(@sfs);
    warn sprintf "PROCESSED:%d\n", scalar(@sfs);
}
is(@sfs, 98);

# this sequence has no locus_tag or or gene tags
$seq = getseq("AY763288.gb");

# UNFLATTEN
@sfs = $unflattener->unflatten_seq(-seq=>$seq,
				   -use_magic=>1
                                  );
if( $verbosity > 0 ) {
    warn "\n\nPOST PROCESSING:\n";
    write_hier(@sfs);
    warn sprintf "PROCESSED:%d\n", scalar(@sfs);
}
is(@sfs, 3);


# try again; different sequence - dicistronic gene, mRNA record

$seq = getseq("X98338_Adh-mRNA.gb");

# UNFLATTEN
@sfs = $unflattener->unflatten_seq(-seq=>$seq,
                                 -partonomy=>{'*'=>'gene'},
                                );
if( $verbosity > 0 ) {                                 
    warn "\n\nPOST PROCESSING:\n";
    write_hier(@sfs);
    warn sprintf "PROCESSED:%d\n", scalar(@sfs);
}
is(@sfs, 7);

# try again; this sequence has no CDSs but rRNA present

$seq = getseq("no_cds_example.gb");

# UNFLATTEN
@sfs = $unflattener->unflatten_seq(-seq=>$seq,
                                 use_magic=>1
                                );
if( $verbosity > 0 ) {
    warn "\n\nPOST PROCESSING:\n";
    write_hier(@sfs);
    warn sprintf "PROCESSED:%d\n", scalar(@sfs);
}

my @all_sfs = $seq->get_all_SeqFeatures;

my @exons = grep { $_-> primary_tag eq 'exon' }  @all_sfs ; 

is(@exons, 2);



sub write_hier {
    my @sfs = @_;
    _write_hier(0, @sfs);
}

sub _write_hier {
    my $indent = shift;
    my @sfs = @_;
    foreach my $sf (@sfs) {
        my $label = '?';
        if ($sf->has_tag('product')) {
            ($label) = $sf->get_tag_values('product');
        }
        warn sprintf "%s%s $label\n", '  ' x $indent, $sf->primary_tag;
        my @sub_sfs = $sf->sub_SeqFeature;
        _write_hier($indent+1, @sub_sfs);
    }
}

sub getseq {
    my @path = @_;
    my $seqio =
      Bio::SeqIO->new('-file'=> test_input_file(@path), 
                      '-format' => 'GenBank');
    $seqio->verbose($verbosity);

    my $seq = $seqio->next_seq();
    return $seq;
}
