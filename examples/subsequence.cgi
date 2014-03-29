#!/usr/bin/perl

# see http://zfish.nichd.nih.gov/tools/subsequence.cgi

# uncomment and modify the next two lines 
#  if your perl is in a nonstandard directory
#use lib '/disk3/local/lib/perl5/site_perl';
#use lib '/disk3/local/lib/perl5/';

use CGI qw/:standard :html3/;
use Bio::DB::GenBank;
use File::Temp;
use FileHandle;

print header,
  start_html(-title => 'find subsequence of large GenBank entries',-author => 'Jonathan_Epstein\@nih.gov');
print_form()    unless param;
print_results() if param;

sub print_results {
  $gb = new Bio::DB::GenBank;
  $accession = param('accession');
  eval {
    $seq = $gb->get_Seq_by_acc($accession); # Accession Number
  };
  if ($@) {
    print "***ERROR: accession $accession not found***\n";
    return;
  }
  $segment_start = param('start');
  $segment_end = param('length_or_end_value');
  $segment_end = $segment_start+$segment_end-1 if param('length_or_end_choice') eq 'Length';
  if ($segment_end<$segment_start || $segment_start<0) {
    print "***ERROR: invalid segment start and end values:$segment_start,$segment_end***\n";
    return;
  }
  $len = $seq->length();
  if ($segment_end>$len) {
    print "***ERROR: maximum length $len exceeded***\n";
    return;
  }
  $subseq = $seq->subseq ($segment_start,$segment_end);
  
  $name = "subsequence of $accession";
  $strand = "+";
  $strand = "-" if (param('reverse'));
  
  # For some reason, there seems to be a problem if you use the file
  # handle provided by File::Temp.  Similarly, there's a problem if you
  # pass a filename to BioPerl below rather than a file handle.  However,
  # constructing our own file handle and then passing it to BioPerl works
  # fine.
  (undef, $filename) = File::Temp::tempfile();
  $fh = new FileHandle "> $filename";
  $seqoutlong = Bio::SeqIO->new( '-format' => 'Fasta',-fh => $fh);
  $seqobj = Bio::PrimarySeq->new ( -seq => $subseq,
				   -id  => $name . "[length:$len]:" . $segment_start . "-" . $segment_end . "(" . $strand . "strand)",
				   -moltype => 'dna'
				 );
  $seqobj = $seqobj->revcom if ($strand ne "+");
  $seqoutlong->write_seq($seqobj);
  $fh->close;
  undef $fh;
  
  # Now we parse the FASTA file which was just generated, and perform
  # some simple conversions to HTML.   
  open my $TEMPORARY, '<', $filename or die "Could not read temporary file '$filename': $!\n";
  print "<tt>\n";
  while (<$TEMPORARY>) {
    print $_;
    print "<br>\n";
  }
  close $TEMPORARY;
  print "</tt>\n";
  unlink $filename;
}

sub print_form {
  print p("This web page permits you to extract a short subsequence of DNA from a large GenBank entry.  This is especially useful in an era of huge \"contigs\" of genomic DNA, where you only want to extract a few hundred base pairs for subsequent analysis.\n");
  
  print p,"This program also illustrates the power of ",a({-href => 'http://www.BioPerl.org/'}, "BioPerl"), ", a powerful set of tools for molecular biology analysis.  The ", a({-href => 'subsequence.pl.txt'}, "source code"), " for this program is less than 90 lines long.\n";
  
  print p,"You must specify the GenBank accession number along with a start position.  You may specify either the length of the subsequence you wish to extract or, equivalently, the endpoint.\n";
  
  print "The sequence may be reverse-complemented if you wish, e.g., the reverse complement of <font color=green>ATCGC</font> is <font color=yellow>GCGAT</font>.\n";
  
  print p,"To test this web page, try accession NT_004002, start 50000, length 400.\n";
  
  print start_form,table(
			 Tr(td("Enter your GenBank accession"),td(textfield(-name => 'accession',-size => 20))),
			 Tr(td("Start position"),td(textfield(-name => 'start',-size => 10))),
			 Tr(td("Specify length or end position"), td(radio_group (-name => 'length_or_end_choice',-values => [Length, End], default => Length))),
			 Tr(td("Length or end position"), td(textfield (-name => length_or_end_value,-size => 20))),
			 Tr(td("Reverse complement?"), td(checkbox (-name => 'reverse')))),
    submit ("Find my subsequence");
  
  print hr(),"Credits: Jonathan Epstein (Jonathan_Epstein\@nih.gov)";

}
