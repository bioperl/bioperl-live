#!/usr/bin/env perl

use strict;
use URI;
use LWP::UserAgent;
use HTTP::Request;
use XML::LibXML;
use File::Temp qw(tempfile tempdir);

use Bio::DB::BioDB;
use Bio::Seq;
use Bio::Annotation::OntologyTerm;
use Bio::Ontology::Ontology;
use Bio::Ontology::Term;
use Bio::Ontology::RelationshipType;
use Bio::Ontology::Relationship;
use Bio::SeqIO;

# Open the config file.
my $root = XML::LibXML->new()->parse_file('conf.xml')->getDocumentElement() or die("Could not parse XML config file");
my $node;

#
# create the DBAdaptorI for our database
#
$node = @{$root->findnodes('database')}[0];
my $db = Bio::DB::BioDB->new(-database => "biosql",
                             -driver   => $node->findvalue('driver'),
                             -host     => $node->findvalue('host'),
                             -dbname   => $node->findvalue('dbname'),
                             -user     => $node->findvalue('username'),
                             -pass     => $node->findvalue('password')
                             );

#
# create the UserAgent to do our web fetches
#
my $ua = LWP::UserAgent->new();
$node = @{$root->findnodes('httpproxy')}[0];
$ua->proxy(['http','https'],$node->findvalue('url')) if $node;

print "Running parsers...\n";

# Loop through the parser children of parsers
for my $parser (@{$root->findnodes('parsers')}[0]->findnodes('parser')) {
	my $parsertype = $parser->findvalue('@type');
	my $sourcenamespace = $parser->findvalue('sourcenamespace');
	my $targetnamespace = $parser->findvalue('targetnamespace');
	my $relnamespace = "Term Importer Map";
	my $script = $parser->findvalue('../@scriptdir').'/'.$parser->findvalue('script');
	print "Parser ".$parser->findvalue('@name')." (type ".$parsertype.") lives at ".$script."\n";
	print "(links ".$sourcenamespace." to ".$targetnamespace.")\n";

	if (! ($parsertype eq 'TERM2TERM' || $parsertype eq 'SEQ2TERM')) {
		print("!!! ".$parsertype." is not a valid parser type. Skipping.");
		next;
	}

	my $uri = new URI;
	$uri->scheme($parser->findvalue('server/@protocol'));
	$uri->host($parser->findvalue('server/host'));
	$uri->path($parser->findvalue('server/filename'));
	print "Downloading data from: ".$uri->as_string()."\n";

	my $res = $ua->get($uri->as_string());
	die("Could not download data!") unless $res->is_success();

	my ($tmpfh,$tmpfilename) = tempfile();
	print $tmpfh $res->content();
	close($tmpfh);
	print "...downloaded to ".$tmpfilename.". Parser starting...\n";

	# Create the ontology references (for TERM2TERM) or database references (for SEQ2TERM)
	my $sourceontology = Bio::Ontology::Ontology->new(-name=>$sourcenamespace);
	my $targetontology;
	my $relontology;
	my $IS_A;
	if ($parsertype eq 'TERM2TERM') {
		$targetontology = Bio::Ontology::Ontology->new(-name=>$targetnamespace);
		$relontology = Bio::Ontology::Ontology->new(-name=>$relnamespace);
		$IS_A = Bio::Ontology::RelationshipType->get_instance( "IS_A" );
		$IS_A->ontology($relontology);
	}

	# Run the parser
	my ($tmpoutfh,$tmpoutfilename) = tempfile();
	system($script.' <'.$tmpfilename.' >'.$tmpoutfilename);

	# Read it's output - must be source (subject), target (object) CSV format with no quotes
	open(my $parsefh,$tmpoutfilename) or die("Could not read parser output!");
	my $counter = 0;
	while (<$parsefh>) {
		chomp; # Remove trailing newline
		my ($sourceid,$targetid) = split /,/,$_;

		my $sourceid_accession;
		my $sourceid_version;
		if ($parsertype eq 'SEQ2TERM') {
			($sourceid_accession,$sourceid_version) = split /\./,$sourceid;
		}

		#
	  	# look up source
		#
		my $sourceobject = ($parsertype eq 'TERM2TERM' ? Bio::Ontology::Term->new(-identifier=>$sourceid,-ontology=>$sourceontology) : Bio::Seq->new(-accession_number=>$sourceid_accession, -version=>$sourceid_version, -namespace=>$sourcenamespace));
            	my $sourceadp = $db->get_object_adaptor($sourceobject);
            	$sourceobject = $sourceadp->find_by_unique_key($sourceobject,-flat_only=>1);
		if (! $sourceobject) {
			print "...could not find source object ".$sourceid." - skipping...\n";
			next;
		}


		if ($parsertype eq 'TERM2TERM') {
			# Create a term_relation

			#
			# look up target
			#
			my $targetobject = Bio::Ontology::Term->new(-identifier=>$targetid,-ontology=>$targetontology);
            		my $targetadp = $db->get_object_adaptor($targetobject);
            		$targetobject = $targetadp->find_by_unique_key($targetobject);
			if (! $targetobject) {
				print "...could not find target object ".$targetid." (for source ".$sourceid.") - skipping...\n";
				next;
			}
			my $rel = $db->create_persistent(Bio::Ontology::Relationship->new(-subject_term=>$sourceobject,
								   -predicate_term=>$IS_A,
								   -object_term=>$targetobject,
								   -ontology=>$relontology
								  ));
            		$rel->create();
            		$rel->commit();
		} elsif ($parsertype eq 'SEQ2TERM') {
			# Create a bioentry_qualifier_value

   			my $ann = $db->create_persistent(Bio::Annotation::OntologyTerm->new(-identifier=>$targetid,-ontology=>$targetontology));
			$sourceobject->annotation()->add_Annotation($targetnamespace,$ann);
			$sourceobject->store();
			$sourceadp->commit();
		}

		$counter = $counter + 1;
		if ($counter % 1000 == 0) {
			print "...".$counter." records processed...\n";
		}	
	}
	close($parsefh);
	unlink($tmpfilename);
	unlink($tmpoutfilename);

	print "...parser finished with ".$counter." records.\n";
}

print "All parsers run.\n";

# All Done. 
