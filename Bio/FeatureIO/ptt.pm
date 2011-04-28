=pod

=head1 NAME

Bio::FeatureIO::ptt - read/write features in PTT format

=head1 SYNOPSIS

 # read features 
 my $fin = Bio::FeatureIO->new(-file=>'genes.ptt', -format=>'ptt');
 my @cds;
 while (my $f = $fin->next_feature) {
   push @cds, $f if $f->strand > 0;
 }

 # write features (NOT IMPLEMENTED)
 my $fout = Bio::FeatureIO->new(-fh=>\*STDOUT, -format=>'ptt');
 for my $f (@cds) {
   $fout->write_feature($f);
 }

=head1 DESCRIPTION

The PTT file format is a table of protein features. 
It is used mainly by NCBI who produce PTT files for 
all their published genomes found in L<ftp://ftp.ncbi.nih.gov/genomes/>.
It has the following format:

=over 4

=item Line 1

Description of sequence to which the features belong
 eg. "Leptospira interrogans chromosome II, complete sequence - 0..358943"

It is usually equivalent to the DEFINITION line of a Genbank file,
with the length of the sequence appended. It is unclear why "0" is 
used as a starting range, it should be "1".

=item Line 2

Number of feature lines in the table
 eg. "367 proteins"

=item Line 3

Column headers, tab separated
 eg. "Location  Strand  Length  PID Gene  Synonym Code  COG Product"

 Location : "begin..end" span of feature
 Strand   : "+" or "-"
 Length   : number of amino acids excluding the stop codon
 PID      : analogous to Genbank /db_xref="GI:xxxxxxxxx"
 Gene     : analogous to Genbank /gene="xxxx"
 Synonym  : analogous to Genbank /locus_tag="xxxx"
 Synonym  : analogous to Genbank /locus_tag="xxxx"
 COG      : CDD COG code with COG letter categories appended
 Product  : analogous to Genbank /product="xxxx"

=item Line 4 onwards

Feature lines, nine columns, tab separated, "-" used for empty fields
 eg. "2491..3423  + 310 24217063  metF  LB002 - COG0685E  5,10-methylenetetrahydrofolate reductase"


=back

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

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR - Torsten Seemann

Email torsten.seemann AT infotech.monash.edu.au

=head1 CONTRIBUTORS

Based on bed.pm and gff.pm by Allen Day.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::FeatureIO::ptt;

use strict;
use base qw(Bio::FeatureIO);
use Bio::SeqFeature::Generic;

# map tab-separated column number to field name
our %NAME_OF = (
  0 => 'Location',
  1 => 'Strand',
  2 => 'Length', 
  3 => 'PID', 
  4 => 'Gene',  
  5 => 'Synonym',
  6 => 'Code',  
  7 => 'COG', 
  8 => 'Product',
);
our $NUM_COL = 9;

=head2 _initialize

 Title   : _initialize
 Function: Reading? parses the header of the input
           Writing? 

=cut

sub _initialize {
  my($self,%arg) = @_;

  $self->SUPER::_initialize(%arg);

  if ($self->mode eq 'r') {
    # Line 1
    my $desc = $self->_readline();
    chomp $desc;
    $self->description($desc);
    # Line 2
    my $line = $self->_readline();
    $line =~ m/^(\d+) proteins/ or $self->throw("Invalid protein count");
    $self->protein_count($1);
    # Line 3
    $self->_readline();
  }
}

=head2 next_feature

 Title   : next_feature
 Usage   : $io->next_feature()
 Function: read the next feature from the PTT file
 Example : 
 Args    : 
 Returns : Bio::SeqFeatureI object

=cut

sub next_feature {
  my $self = shift;
  $self->mode eq 'r' || return; # returns if can't read next_feature when we're in write mode
  
  my $line = $self->_readline() or return; # returns if end of file, no more features?
  chomp $line;
  my @col = split m/\t/, $line;
  @col==$NUM_COL or $self->throw("Too many columns for PTT line");

  $col[0] =~ m/(\d+)\.\.(\d+)/ or $self->throw("Invalid location (column 1)");
  my $feat = Bio::SeqFeature::Generic->new(-start=>$1, -end=>$2, -primary=>'CDS');
  $col[1] =~ m/^([+-])$/ or $self->throw("Invalid strand (column 2)");
  $feat->strand($1 eq '+' ? +1 : -1);
  for my $i (2 .. $NUM_COL-1) {
    $feat->add_tag_value($NAME_OF{$i}, $col[$i]) if $col[$i] ne '-';
  }
  return $feat;
}

=head2 write_feature (NOT IMPLEMENTED)

 Title   : write_feature
 Usage   : $io->write_feature($feature)
 Function: write a Bio::SeqFeatureI object in PTT format
 Example : 
 Args    : Bio::SeqFeatureI object
 Returns : 

=cut

sub write_feature {
  shift->throw_not_implemented;
}

=head2 description

 Title   : description
 Usage   : $obj->description($newval)
 Function: set/get the PTT file description for/from line one
 Example : 
 Returns : value of description (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub description {
  my $self = shift;
  return $self->{'description'} = shift if @_;
  return $self->{'description'};
}

=head2 protein_count

 Title   : protein_count
 Usage   : $obj->protein_count($newval)
 Function: set/get the PTT protein count for/from line two
 Example : 
 Args    : on set, new value (a scalar or undef, optional)
 Returns : value of protein_count (a scalar)

=cut

sub protein_count {
  my $self = shift;
  return $self->{'protein_count'} = shift if @_;
  return $self->{'protein_count'};
}

1;
