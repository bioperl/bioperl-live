#!/usr/bin/env perl

use strict;

while (<>) {
	# Process STDIN (the interpro2go file from geneontology.org
	# For each line, output to STDOUT the following:
	#	sourceaccession (IPR),targetaccession (GO)
	# where dbid is the human-readable non-unique IPR reference,
	# dbaccession is the unique IPR reference, and goaccession is the
	# unique GO accession.

	# Interpro lines look like this:
	# InterPro:IPR000018 P2Y4 purinoceptor > GO:purinergic nucleotide receptor activity, G-protein coupled ; GO:0045028

	/^InterPro:/ or next;
	s/^InterPro:(IPR\d+).*>.*;\s+(GO:\d+)/$1,$2/;
	print;
}
