#!/usr/bin/perl -w

# $Id$

# This script generates a HTML or text FAQ from the XML version.
#
# based on twig_faq.pl by Michael Rodriguez
# Modified for Bioperl FAQ by Heikki Lehvaslaiho

use strict;
use XML::Twig;
use Text::Wrap;
use Getopt::Long;

my $text = '';
my $html = ''; 
GetOptions('text' => \$text, 'html' => \$html);
$text || $html || die "Usage: faq.pl (-t[ext] | -h[tml) <xml_file>";
$html and $text = ''; #preference for HTML output

my $file= shift || die "Usage faq.pl (-t[ext] | -h[tml) <xml_file>";

# those 2 globals are grabbed from the xml header and used to
# generate the html header
my $G_TITLE;   # title
my $G_VERSION; # version

my $SNB = -1;     # section index
my $QNB =  0;     # question index

# text mode globals
my $SEP = "-" x 75;
$SEP .= "\n";
my @QA;



#our ToC will contain two levels: header and question
my $TOC;       # the element containing the ToC
my $TOC_HEADER;# the element containing the latest ToC header

my $t= new XML::Twig( TwigHandlers =>
                      { title        => \&title,
                        version      => \&version,
                        header       => \&header,
                        credits      => \&section,
                        overview     => \&section,
			section      => \&erase,
			section_name => \&block,
                        'q'          => \&erase,
                        question     => \&question,
                        answer       => \&answer,
                        code         => \&code,
                        copyright    => \&copyright,
                      },
		      pretty_print   => 'indented',
		      keep_spaces_in => ['pre'],
		      empty_tags     => 'html',
                     );

$t->parsefile( $file);


my $faq= $t->root;
$faq->set_gi( 'html');
$faq->insert( 'body');

# create the html header
my $title= new XML::Twig::Elt( 'title', "$G_TITLE - version $G_VERSION");
my $html_header= new XML::Twig::Elt( 'header', $title);

$html_header->paste( $faq);

$faq->print if $html;

if ($text) {
    foreach (@QA) {
	print;
    }
    print "\n";
}

exit;

# end of main ##########################################################

sub title { 
    my( $t, $title)= @_;
    $title->set_gi( 'h1');      # title -> h1
    $G_TITLE= $title->text;     # store the title
}

sub version { 
    my( $t, $version)= @_;
    $version->set_gi( 'h2');    # version -> h2
    $G_VERSION= $version->text; # store the version
    $version->prefix( 'Version ');
}

sub author {
    my( $t, $author)= @_; }

sub header { 
    my( $t, $header)= @_;

    $header->set_gi( 'center'); # center the whole header

    # add the appropriate text around the authors
    my @author= $header->children( 'author');

    if ($text) {
	print "\n$G_TITLE\n-----------\nv. $G_VERSION\n\n";
	last unless @author;
	print "This FAQ maintained by:\n";

	foreach my $author( @author) {
	    print "* ", $author->text, "\n";
	}
	print "\n";
    }

    my $first_author= shift @author;
    if( $first_author) { 
	my $by= new XML::Twig::Elt( '#PCDATA', "Maintained by ");
        $by->paste( 'before', $first_author);
        $first_author->set_gi( 'b');      # author -> b
    }
    my $last_author= pop @author;
    if( $last_author) {
	my $and= new XML::Twig::Elt( '#PCDATA', " and ");
        $and->paste( 'before', $last_author);
        $last_author->set_gi( 'b');      # author -> b
    }
    foreach my $author( @author) {
	my $comma= new XML::Twig::Elt( '#PCDATA', ", ");
        $comma->paste( 'before', $author);
        $author->set_gi( 'b');      # author -> b
    }
}

sub section  { 
    my( $t, $section)= @_;


    if ($text) {
	print "\n$SEP\n", ucfirst($section->gi), "\n\n$SEP", ;
	foreach my $para ($section->children) {
	    print "\n\n", wrap("        ", "        ", $para->text. "\n");
	}

    }


    # add an <hr> 
    my $hr= new XML::Twig::Elt( 'hr');
    $hr->paste( $section);
    # add the title
    my $section_title= ucfirst( $section->gi);
    my $title= new XML::Twig::Elt( 'h3', $section_title);
    $title->paste( 'after', $hr);
    $section->set_gi( 'div');

}

sub erase  {
    my( $t, $elt)= @_;
    $elt->erase;
}


sub block {
    my( $t, $block)= @_;

    # add an <hr />
    my $hr= new XML::Twig::Elt(  'hr');
    $hr->paste( 'before', $block);

    # add the section number
    $SNB++;
    my $target= "$SNB";

    my $blocktext = $block->text;

    my $section=  new XML::Twig::Elt( 'h2', "$target. $blocktext" );
    $section->paste( 'after', $hr);

    # add a target element (so the ToC can link there)
    my $a= new XML::Twig::Elt( 'a');
    $a->set_att( name => $target);
    $a->paste( 'after', $section);


    # create the ToC if need be
    unless( $TOC) {

	if ($text) {
	    print "\n$SEP\n", "Contents" , "\n\n$SEP", ;
	}

	my $q= $block->parent;
        $TOC= new XML::Twig::Elt( 'blockquote');
        $TOC->paste( 'before', $q);
        my $toc_title=  new XML::Twig::Elt( 'h3', 'Contents');
        $toc_title->paste( 'before', $TOC);
	my $hr= new XML::Twig::Elt ( 'hr');
	$hr->paste( 'before', $toc_title);
    }

    my $entry= parse XML::Twig::Elt
	("<strong><a href=\"#$target\">$target. $blocktext</a></strong>");
    $entry->paste( 'last_child', $TOC);

    $TOC_HEADER= new XML::Twig::Elt( 'ul');
    $TOC_HEADER->paste( 'last_child', $TOC);

    # create the toc entry
    #generate_top_toc_entry( $block->id, $section->text);

    if ($text) {
	push @QA, "\n\n", $SEP;
	push @QA, "\n", $section->text, "\n\n";
	push @QA, $SEP, "\n";
	print "\n",$section->text, "\n\n";
    }

    # reset question counter for this section
    $QNB=0;

    # this is not a HTML tag
    $block->set_text('');
    $block->erase;

}

sub question { 
    my( $t, $question)= @_;

    # add an <hr> 
    my $hr= new XML::Twig::Elt( 'p');
    $hr->paste( 'before', $question);

    # store teh original text for headers
    my $question_text = $question->text;

    # add the question number
    $QNB++;
    my $target= "Q$SNB.$QNB";
    my $number= new XML::Twig::Elt( '#PCDATA', "$target: ");
    $number->paste( $question);

    # add a target element (so the ToC can link there)
    my $a= new XML::Twig::Elt( 'a');
    $a->set_att( name => $target);
    $a->paste( 'before', $question);

    # create the ToC if need be
    unless( $TOC) {
	my $q= $question->parent->parent;
        $TOC= new XML::Twig::Elt( 'ul');
        $TOC->paste( 'before', $q);
        my $toc_title=  new XML::Twig::Elt( 'h3', 'Content');
        $toc_title->paste( 'before', $TOC);
	my $hr= new XML::Twig::Elt ( 'hr');
	$hr->paste( 'before', $toc_title);
	$TOC_HEADER = $TOC; # only one level toc
    }
    # create the toc entry
    generate_toc_entry( $target, "<strong>$target: </strong>". $question_text);

    if ($text) {

	push @QA, "\n\n", wrap("  ", "        ", $question->text. "\n");
	print wrap("  ", "        ", $question->text. "\n");
    }

    $question->set_gi( 'h3');
}

sub answer { 
    my( $t, $answer)= @_;

    if ($text) {
	if ($answer->first_child && $answer->first_child->gi eq 'p') {
	    foreach my $para ($answer->children) {
		push @QA, "\n", wrap("        ", "        ", $para->text. "\n");
	    }
	} else {
	    push @QA, "\n", wrap("        ", "        ", $answer->text. "\n");
	}
    }


    # replace the answer by a 'p'$answer->first_child
    # unless the first child is already a 'p'
    my $first_child= $answer->first_child;
    if( $first_child->gi eq 'p') {
	$answer->erase;
    } else {
	$answer->set_gi( 'p');
    }
}

sub code {
    my( $t, $code)= @_;
    my $table = new XML::Twig::Elt( 'table');
    $table->set_att( bgcolor => "light grey", border => 0, cellspacing => 0, cellpadding => 10);
    my $tr = new XML::Twig::Elt( 'tr');
    $tr->paste($table);
    my $td = new XML::Twig::Elt( 'td');
    $td->paste($tr);
    my $font = new XML::Twig::Elt( 'font');
    $font->set_att(color => "blue");
    $font->paste($td);
    $table->paste('before', $code);
    $code->move($font);

    $code->set_gi( 'pre');

}


sub copyright {
    my( $t, $copyright)= @_;

    if ($text) {
	push @QA, "\n$SEP", wrap("", "", $copyright->text. "\n");
    }

    # add an <hr />
    my $hr= new XML::Twig::Elt( 'hr');
    $hr->paste( 'before', $copyright);

    $copyright->set_gi( 'p');
}



sub generate_toc_entry {
    my( $target, $text)= @_;
    my $entry= parse XML::Twig::Elt(  "<li><a href=\"#$target\">$text</a></li>");
    $entry->paste( 'last_child', $TOC_HEADER);
}

