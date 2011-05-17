#!/usr/bin/perl
# gSequence - Protein Sequence Control Panel
# by Lorenz Pollsk
#
# this is work in progress! use this only for testing

use Gtk;
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::SeqStats;
use Bio::SeqFeature::Generic;
use Bio::Index::Abstract;
use Bio::DB::GenBank;
use Bio::DB::GenPept;

init Gtk;

# constant
my $false = 0;
my $true = 1;

# global widgets
my ($main_notebook,@main_label,@seq_edit);
my $about_dialog;
my ($import_dialog,$import_entry,@import_buttons,$import_from);
my ($description_window,$description_edit);
my ($comment_window,$comment_edit,$current_comment,$comment_frame);
my ($seqstats_window,$seqstats_edit);
my ($dblink_window,@dblink_entry,$current_dblink,$dblink_clist,$dblink_handler_id);
my ($ref_window,@ref_entry,$current_ref,$ref_clist,$ref_handler_id);
my ($feature_window,@feature_entry,$current_feature_item,@feature_spinner,
    $feature_handler_id,$feature_tree);
my ($pref_window,@pref_entry);

# global file data
my @seq;
my @filename;
my @modified;
my @locked; # locked sequence for editing ?
my $current;

# menu
my @menu_items = ( { path        => '/_File',
		     type        => '<Branch>' },
		   { path        => '/File/_New',
		     accelerator => '<control>N',
		     callback    => \&new },
		   { path        => '/File/_Open SwissProt',
		     accelerator => '<control>O',
		     callback    => \&open_dialog },
		   { path        => '/File/_Save SwissProt',
		     accelerator => '<control>S',
		     callback    => \&save },
		   { path        => '/File/Save _As...',
		     callback    => \&saveas_dialog },
		   { path        => '/File/Close',
		     callback    => \&close },
		   { path        => '/File/sep1',
		     type        => '<Separator>' },
		   { path        => '/File/_Import from...',
		     type        => '<Branch>' },
                   { path        => '/File/Import from.../Remote DB',
		     type        => '<Branch>' },
                   { path        => '/File/Import from.../Remote DB/AceDB',
		     callback    => sub { &seq_import("ace"); } },
                   { path        => '/File/Import from.../Remote DB/GenPept',
		     callback    => sub { &seq_import("genpept"); } },
                   { path        => '/File/Import from.../Flat File Index',
		     type        => '<Branch>' },
                   { path        => '/File/Import from.../Flat File Index/Fasta',
		     callback    => sub { &seq_import("fasta"); } },
                   { path        => '/File/Import from.../Flat File Index/SwissProt',
		     callback    => sub { &seq_import("swissprot"); } },
                   { path        => '/File/Import from.../Flat File Index/SwissPfam',
		     callback    => sub { &seq_import("swisspfam"); } },
		   { path        => '/File/_Export to...' },
		   { path        => '/File/sep2',
		     type        => '<Separator>' },
		   { path        => '/File/_Quit',
		     callback    => sub { Gtk->exit( 0 ); } },

		   { path        => '/_Edit',
		     type        => '<Branch>' },
		   { path        => '/Edit/C_ut',
		     callback    => sub { $seq_edit[$current]->cut_clipboard(); },
		     accelerator => '<control>X' },
		   { path        => '/Edit/_Copy',
		     callback    => sub { $seq_edit[$current]->copy_clipboard(); },
		     accelerator => '<control>C' },
		   { path        => '/Edit/_Paste',
		     callback    => sub { $seq_edit[$current]->paste_clipboard(); },
		     accelerator => '<control>V' },
		   { path        => '/Edit/Select All',
		     callback    => sub { $seq_edit[$current]->select_region(0,-1); } },

		   { path        => '/_Specs',
		     type        => '<Branch>' },
		   { path        => '/Specs/_Sequence Stats',
		     callback    => sub {&update_seqstats_window(1);} },
		   { path        => '/Specs/sep1',
		     type        => '<Separator>' },
		   { path        => '/Specs/_Description',
		     callback    => sub {&update_description_window(1);} },
		   { path        => '/Specs/_Comments',
		     callback    => sub {&update_comment_window(1);} },
		   { path        => '/Specs/_DB Links',
		     callback    => sub {&update_dblink_window(1);} },
		   { path        => '/Specs/_References',
		     callback    => sub {&update_reference_window(1);} },
		   { path        => '/Specs/sep2',
		     type        => '<Separator>' },
		   { path        => '/Specs/_Features',
		     callback    => sub {&update_feature_window(1);} },

		   { path        => '/_Tools',
		     type        => '<Branch>' },
		   { path        => '/Tools/Code Table' },
		   { path        => '/Tools/sep1',
		     type        => '<Separator>' },
		   { path        => '/Tools/local Blast' },
		   { path        => '/Tools/local HMMER' },
		   { path        => '/Tools/hmmpfam' },
		   { path        => '/Tools/web Blast' },

		   { path        => '/_Options',
		     type        => '<Branch>' },
		   { path        => '/Options/_Preferences',
		     callback    => sub {&update_pref_window(1);} },

		   { path        => '/_Help',
		     type        => '<LastBranch>' },
		   { path        => '/Help/Help' },
		   { path        => '/Help/_About...',
		     callback    => sub { $about_dialog->show_all();} } );

### main

$current = 0;
&init_windows();
main Gtk;
exit( 0 );


### Subroutines

sub init_windows
{
    &init_main_window();
    &init_about_dialog();
    &init_import_dialog();
    &init_seqstats_window();
    &init_description_window();
    &init_comment_window();
    &init_dblink_window();
    &init_reference_window();
    &init_feature_window();
    &init_pref_window();
}

sub init_main_window
{
    # toplevel window
    my $window;
    $window = new Gtk::Window( 'toplevel' );
    $window->signal_connect( 'destroy', sub { Gtk->exit( 0 ); } );
    $window->set_title( "gSequence" );
    $window->set_usize( 600, 400 );

    # vertical box containing menu and text editor widget
    my $main_vbox;
    $main_vbox = new Gtk::VBox( $false, 1 );
    $main_vbox->border_width( 1 );
    $window->add( $main_vbox );

    # handlebox for menubar
    my $handlebox;
    $handlebox = new Gtk::HandleBox();
    $main_vbox->pack_start( $handlebox, $false, $true, 0 );

    # menubar
    my $menubar;
    $menubar = get_menu( $window );
    $handlebox->add( $menubar );

    # text widget
    $seq_edit[$current] = new Gtk::Text( undef, undef );
    $seq_edit[$current]->set_editable( $true );

    # vertical scrollbar for text widget
    my $scrollbar;
    $scrollbar = new Gtk::VScrollbar( $seq_edit[$current]->vadj );

    # horizontal box containing text widget and scrollbar
    my $seq_edit_hbox;
    $seq_edit_hbox = new Gtk::HBox( $false, 1 );
    $seq_edit_hbox->border_width( 1 );
    $seq_edit_hbox->pack_start( $seq_edit[$current], $true, $true, 0);
    $seq_edit_hbox->pack_end( $scrollbar, $false, $true, 0);

    $main_notebook = new Gtk::Notebook();
    $main_notebook->set_tab_pos( 'top' );

    $main_vbox->pack_end( $main_notebook, $true, $true, 0);

    # show everything
    $window->show_all();

    $main_notebook->signal_connect_after("switch-page",
       sub{ #$seq[$current]->seq($seq_edit[$current]->get_chars(0,-1)) 
	    #   if (defined($seq[$current]));
	    $current = $main_notebook->get_current_page(); 
	    &update_seq_data(); } );
}

sub get_menu
{
    my ( $window ) = @_;
    
    my $menubar;
    my $item_factory;
    my $accel_group;

    $accel_group = new Gtk::AccelGroup();

    # This function initializes the item factory.
    # Param 1: The type of menu - can be 'Gtk::MenuBar', 'Gtk::Menu',
    #          or 'Gtk::OptionMenu'.
    # Param 2: The path of the menu.
    # Param 3: The accelerator group.  The item factory sets up
    #          the accelerator table while generating menus.
    $item_factory = new Gtk::ItemFactory( 'Gtk::MenuBar',
					  '<main>',
					  $accel_group );

    # This function generates the menu items. Pass the item factory,
    # the number of items in the array, the array itself, and any
    # callback data for the the menu items.
    $item_factory->create_items( @menu_items );

    # Attach the new accelerator group to the window.
    $window->add_accel_group( $accel_group );

    # Finally, return the actual menu bar created by the item factory.
    #*menubar = gtk_item_factory_get_widget (item_factory, "&lt;main>");
    return ( $item_factory->get_widget( '<main>' ) );
}

sub new_seq_page
{
    my ($seq) = shift;
    my $curr;

    push @seq,$seq;
    $curr = @seq - 1;
    $main_label[$curr] = new Gtk::Label($seq[$curr]->id())
	if (defined($seq[$curr]));
    $main_label[$curr] = new Gtk::Label("<New>")
	if (!defined($seq[$curr]));

    # text widget
    $seq_edit[$curr] = new Gtk::Text( undef, undef );
    $seq_edit[$curr]->set_editable( $true );

    # vertical scrollbar for text widget
    my $scrollbar;
    $scrollbar = new Gtk::VScrollbar( $seq_edit[$curr]->vadj );

    # horizontal box containing text widget and scrollbar
    my $seq_edit_hbox;
    $seq_edit_hbox = new Gtk::HBox( $false, 1 );
    $seq_edit_hbox->border_width( 1 );
    $seq_edit_hbox->pack_start( $seq_edit[$curr], $true, $true, 0);
    $seq_edit_hbox->pack_end( $scrollbar, $false, $true, 0);

    $main_notebook->append_page( $seq_edit_hbox, $main_label[$curr] );
    $main_notebook->show_all();
    $main_notebook->set_page(-1);
}

sub seq_fetch
{
    my ($server,$port,$dir,$db); # read from preferences
    my ($dbobj);

    return if (!defined($import_from) || !($import_from));

    $dbobj = Bio::DB::GenPept->new() if ($import_from eq "genpept");
    $dbobj = Bio::DB::Ace->new(-host=>$server,-port=>$port)
	if ($import_from eq "ace");
    $dbobj = Bio::Index::Abstract->new("$dir/$db") 
	if ($import_from eq "fasta") ||
	   ($import_from eq "swissprot") ||
	   ($import_from eq "swisspfam");

    if( $import_buttons[0]->get_active() ) {
	&new_seq_page($dbobj->get_Seq_by_id($import_entry->get_text()));
    } else {
	&new_seq_page($dbobj->get_Seq_by_acc($import_entry->get_text()));
    }
}

sub seq_import
{
    ($import_from) = @_;
    my %names = ( "ace" => "AceDB",
			 "genpept" => "GenPept DB",
			 "fasta" => "Fasta Flat File",
			 "swissprot" => "SwissProt Flat File",
			 "swisspfam" => "SwissPfam Flat File"
		     );
    $import_dialog->set_title("Import from ".$names{$import_from});
    $import_entry->set_text("");
    $import_dialog->show_all();
}

sub init_import_dialog
{
    $import_dialog = new Gtk::Dialog();
    $import_dialog->border_width(5);

    # create the first button and add it to a box
    my $button = new Gtk::RadioButton( "Fetch by ID" );
    $import_dialog->vbox->pack_start($button,$false,$false,2);
         
    # create the second button and add it to a box
    $button = new Gtk::RadioButton( "Fetch by ACCESSION", $button );
    $import_dialog->vbox->pack_start($button,$false,$false,2);
    @import_buttons = $button->group();

    $import_entry = new Gtk::Entry();
    my $frame = new Gtk::Frame("Enter here:");
    $frame->add($import_entry);
    $import_dialog->vbox->pack_start( $frame, $true, $true, 5);
    
    my $bbox = new Gtk::HButtonBox();
    $bbox->set_layout("end");
    
    $button = new Gtk::Button( "OK" );
    $bbox->add( $button );
    $button->signal_connect("clicked",
           # OK button handler
           sub{ $import_dialog->hide();
		&seq_fetch();
	      });

    $button = new Gtk::Button( "Cancel" );
    $bbox->add( $button );
    $button->signal_connect("clicked",
           # close button handler
           sub{ $import_dialog->hide();
	      });

    $import_dialog->action_area->pack_start( $bbox, $true, $true, 0 );

    $import_dialog->signal_connect_after( "delete_event",
           # window delete handler
           sub{ $import_dialog->hide();
                return &Gtk::true;
	      });
}

sub open_dialog
{
    # Create a new file selection widget
    my $open_dialog = new Gtk::FileSelection( "Open File..." );
    # Connect the ok_button to open_ok_sel function
    $open_dialog->ok_button->signal_connect( "clicked",
					     \&ok_open_dialog,
					     $open_dialog );
    # Connect the cancel_button to destroy the widget
    $open_dialog->cancel_button->signal_connect( "clicked",
						 sub { $open_dialog->destroy(); } );
    $open_dialog->show();
}

# Get the selected filename
sub ok_open_dialog
  {
    my ( $widget, $file_selection ) = @_;
    push @filename, $file_selection->get_filename();

    $widget->parent->parent->parent->destroy();

    my $in = Bio::SeqIO->new(-file => $filename[-1] , '-format' => 'swiss');

    &new_seq_page($in->next_seq());
}

sub update_seq_data
{    
    $main_label[$current]->set_text($seq[$current]->id) if (defined($seq[$current]));
    $main_label[$current]->set_text("<New>") if (!defined($seq[$current]));

    $seq_edit[$current]->freeze();
    $seq_edit[$current]->delete_text(0,-1);
    $seq_edit[$current]->insert(undef,undef,undef,$seq[$current]->seq()) if (defined($seq[$current]));
    $seq_edit[$current]->thaw();

    &update_comment_window();
    &update_description_window();
    &update_seqstats_window();
    &update_dblink_window();
    &update_reference_window();
    &update_feature_window();
}

sub new
{
    &new_seq_page(undef);
}

sub close
{
}

sub save
{
    if (!defined($filename[$current])||!$filename[$current])
    {
	&saveas_dialog;
	return;
    }
    my $out = Bio::SeqIO->new(-file => ">$filename[$current]" , '-format' => 'swiss');
    $out->write_seq($seq[$current]);
}

sub saveas_dialog
{
    # Create a new file selection widget
    my $saveas_dialog = new Gtk::FileSelection( "Save As..." );
    # Connect the ok_button to saveas_ok_sel function
    $saveas_dialog->ok_button->signal_connect( "clicked",
					       \&ok_saveas_dialog,
					       $saveas_dialog );
    # Connect the cancel_button to destroy the widget
    $saveas_dialog->cancel_button->signal_connect( "clicked",
						   sub { $saveas_dialog->destroy(); } );
    $saveas_dialog->show();
}

# Get the selected filename and print it to the console
sub ok_saveas_dialog
  {
    my ( $widget, $file_selection ) = @_;
    my $filename = $file_selection->get_filename();
    $widget->parent->parent->parent->destroy();
    $filename[$current] = $filename;
    my $out = Bio::SeqIO->new(-file => ">$filename[$current]" , '-format' => 'swiss');
    $out->write_seq($seq[$current]);
  }

sub init_comment_window
{
    $current_comment = 0;
    
    $comment_window = new Gtk::Dialog();
    $comment_window->set_default_size(650,300);
    $comment_window->set_policy($false,$true,$false);
    $comment_window->set_title("Comments");
    $comment_window->border_width(5);
    
    # frame
    $comment_frame = new Gtk::Frame( "Comment[".$current_comment."]" );
    
    # text widget
    $comment_edit = new Gtk::Text( undef, undef );
    $comment_edit->set_editable( $true );
    $comment_edit->set_word_wrap( $true );
	
    # vertical scrollbar for text widget
    my $scrollbar;
    $scrollbar = new Gtk::VScrollbar( $comment_edit->vadj );
	
    # horizontal box containing text widget and scrollbar
    my $hbox;
    $hbox = new Gtk::HBox( $false, 1 );
    $hbox->border_width( 1 );
    $hbox->pack_start( $comment_edit, $true, $true, 0);
    $hbox->pack_end( $scrollbar, $false, $true, 0);
    $comment_frame->add($hbox);
    $comment_window->vbox->pack_start( $comment_frame, $true, $true, 5);

    my $bbox = new Gtk::HBox( $false, 5 );
    $bbox->border_width(10);
    my $arrow = new Gtk::Arrow('right','out');
    my $button = new Gtk::Button();
    $button->add($arrow);
    $bbox->pack_end( $button, $false, $false, 0);
    $button->signal_connect
	( "clicked", 
	  # next comment button handler
	  sub { return if !defined($seq[$current]);
		 &store_current_comment;
		$current_comment++ 
		    if ($current_comment <((scalar $seq[$current]->annotation->each_Comment)-1));
		&update_comment_window;
	    } );

    $arrow = new Gtk::Arrow('left','out');
    $button = new Gtk::Button();
    $button->add($arrow);
    $bbox->pack_end( $button, $false, $false, 0);
    $button->signal_connect( "clicked", 
           # prev comment button handler
	   sub { return if !defined($seq[$current]);
		 &store_current_comment;
		 $current_comment-- 
		     if ($current_comment > 0);
		 &update_comment_window;
	       } );

    $button = new Gtk::Button("Add");
    $bbox->pack_start( $button, $false, $false, 0);
    $button->signal_connect( "clicked",
           # add comment button handler 
           sub { return if !defined($seq[$current]);
		 &store_current_comment;
		 my $comment = new Bio::Annotation::Comment;
		 $comment->text("");
		 $seq[$current]->annotation->add_Comment( $comment );
		 $current_comment = $seq[$current]->annotation->each_Comment - 1;
		 &update_comment_window;
	       } );

     $button = new Gtk::Button("Delete");
     $bbox->pack_start( $button, $false, $false, 0);
     $button->signal_connect( "clicked", 
           # delete comment button handler
           sub { return if !defined($seq[$current]); 
		 $seq[$current]->annotation->remove_Comment( $current_comment );
		 $current_comment = $current_comment - 1 
		     if ($current_comment > 0);
		 &update_comment_window;		 
	       } );

     $comment_window->vbox->pack_end( $bbox, $false, $false, 0);

     $bbox = new Gtk::HButtonBox();
     $bbox->set_layout("end");
	
     $button = new Gtk::Button( "Close" );
     $bbox->add( $button );
     $button->signal_connect("clicked",
           # close button handler
           sub{ $comment_window->hide();
		&store_current_comment;
	      });

     $comment_window->action_area->pack_start( $bbox, $true, $true, 0 );
     $comment_window->signal_connect_after( "delete_event",
           # window delete handler
           sub{ $comment_window->hide();
		&store_current_comment;
                return &Gtk::true;
	      });
}

sub store_current_comment
{
    (($seq[$current]->annotation->each_Comment)[$current_comment])->
	text($comment_edit->get_chars(0,-1) )
	    if ((defined($seq[$current])) && ($seq[$current]->annotation->each_Comment));
}

sub update_comment_window
{
    my ($show_me) = @_;
    $comment_frame->set_label("Comment[".$current_comment."]");
    # insert comment text
    $comment_edit->freeze();
    $comment_edit->delete_text(0,-1);
    if (defined($seq[$current]))
    {
	my @comment = $seq[$current]->annotation->each_Comment;
	$comment_edit->insert(undef,undef,undef, $comment[$current_comment]->text)
	    if (@comment);
    }
    $comment_edit->thaw();
    
    $comment_window->show_all() if (defined($show_me));
}

sub init_description_window
{
    $description_window = new Gtk::Dialog();
    $description_window->set_default_size(620,250);
    $description_window->border_width(5);
    $description_window->set_title("Description");
    
    # frame
    my $description_frame = new Gtk::Frame( "Description" );
    
    # text widget
    $description_edit = new Gtk::Text( undef, undef );
    $description_edit->set_editable( $true );
    $description_edit->set_word_wrap( $true );	
    
    # vertical scrollbar for text widget
    my $scrollbar;
    $scrollbar = new Gtk::VScrollbar( $description_edit->vadj );
    
    # horizontal box containing text widget and scrollbar
    my $hbox;
    $hbox = new Gtk::HBox( $false, 1 );
    $hbox->border_width( 1 );
    $hbox->pack_start( $description_edit, $true, $true, 0);
    $hbox->pack_end( $scrollbar, $false, $true, 0);
    $description_frame->add($hbox);
    $description_window->vbox->pack_start( $description_frame, $true, $true, 5);
    
    my $bbox = new Gtk::HButtonBox();
    $bbox->set_layout("end");
    
    my $button = new Gtk::Button( "Close" );
    $bbox->add( $button );
    $button->signal_connect("clicked",
           # close button handler
           sub{ $description_window->hide();
		$seq[$current]->desc($description_edit->get_chars(0,-1))
		    if $description_edit->get_chars(0,-1);
	      });

    $description_window->action_area->pack_start( $bbox, $true, $true, 0 );
    $description_window->signal_connect_after( "delete_event",
           # window delete handler
           sub{ $description_window->hide();
		$seq[$current]->desc($description_edit->get_chars(0,-1))
		    if $description_edit->get_chars(0,-1);
                return &Gtk::true;
	      });
}

sub update_description_window
{
    my ($show_me) = @_;
    $description_edit->freeze();
    $description_edit->delete_text(0,-1);
    $description_edit->insert(undef,undef,undef,$seq[$current]->desc)
	if defined($seq[$current]) && defined($seq[$current]->desc);
    $description_edit->thaw();
    
    $description_window->show_all() if (defined($show_me));
}

sub init_seqstats_window
{
    $seqstats_window = new Gtk::Dialog();
    $seqstats_window->border_width(5);
    $seqstats_window->set_default_size(100,250);
    $seqstats_window->set_title("Sequence Statistics");

    # frame
    my $seqstats_frame = new Gtk::Frame( "Sequence Statistics" );
    
    # text widget
    $seqstats_edit = new Gtk::Text( undef, undef );
    $seqstats_edit->set_editable( $false );
    $seqstats_edit->set_word_wrap( $true );
    
    # vertical scrollbar for text widget
    my $scrollbar;
    $scrollbar = new Gtk::VScrollbar( $seqstats_edit->vadj );
    
    # horizontal box containing text widget and scrollbar
    my $hbox;
    $hbox = new Gtk::HBox( $false, 1 );
    $hbox->border_width( 1 );
    $hbox->pack_start( $seqstats_edit, $true, $true, 0);
    $hbox->pack_end( $scrollbar, $false, $true, 0);
    $seqstats_frame->add($hbox);
    $seqstats_window->vbox->pack_start( $seqstats_frame, $true, $true, 5);
    
    my $bbox = new Gtk::HButtonBox();
    $bbox->set_layout("end");
    
    my $button = new Gtk::Button( "Close" );
    $bbox->add( $button );
    $button->signal_connect("clicked",
       # close button handler
       sub{ $seqstats_window->hide();
	  });
    
    $seqstats_window->action_area->pack_start( $bbox, $true, $true, 0 );
    $seqstats_window->signal_connect_after( "delete_event",
       # window delete handler
       sub{ $seqstats_window->hide();
	    return &Gtk::true;
	  });
}

sub update_seqstats_window
{
    my ($show_me) = @_;
    my ($data,$weight,$count_hash,$percent);

    $seqstats_edit->freeze();
    $seqstats_edit->delete_text(0,-1);
    if (defined($seq[$current]))
    {
	$data = $seq[$current]->id."\n\n";
	$weight = Bio::Tools::SeqStats->get_mol_wt($seq[$current]->primary_seq);
	if ($$weight[0] == $$weight[1]) {
	    $data .= "Molecular weight of sequence equals to ".$$weight[0]."\n\n";
	} else {
	    $data .= "Molecular weight of sequence is greater than ";
	    $data .= $$weight[0]." and less than ".$$weight[1]."\n\n";
	}
	$count_hash = Bio::Tools::SeqStats->count_monomers($seq[$current]->primary_seq);
	$data .= "Amino Acids:\n";
	foreach (sort keys %$count_hash)
	{
	    $percent = sprintf "%.1f",
	    (($$count_hash{$_} / $seq[$current]->length)*100);
	    $data .= "${_}: ".$$count_hash{$_}." (${percent}%) \n"
	    }
	$seqstats_edit->insert(undef,undef,undef,$data)
	}
    $seqstats_edit->thaw();
    
    $seqstats_window->show_all() if (defined($show_me));
}

sub init_dblink_window
{
    $current_dblink = 0;
    
    $dblink_window = new Gtk::Dialog();
    $dblink_window->set_default_size(500,400);
    $dblink_window->set_policy($true,$true,$false);
    $dblink_window->set_title("Database Links");
    $dblink_window->border_width(5);
    
    # Create a scrolled window to pack the CList widget into
    my $scrolled_window = new Gtk::ScrolledWindow( undef, undef );
    $dblink_window->vbox->pack_start( $scrolled_window, $true, $true, 0 );
    $scrolled_window->set_policy( 'automatic', 'always' );

    # Create the CList. For this example we use 2 columns
    $dblink_clist = new_with_titles Gtk::CList( "Primary Id","Database" );

    # When a selection is made, we want to know about it. The callback
    # used is selection_made, and its code can be found further down
    $dblink_handler_id = $dblink_clist->signal_connect( "select_row", 
       sub{ return if (!defined($seq[$current]));
	    my ( $clist, $row ) = @_;
	    &store_current_dblink;
	    $current_dblink = $row;
	    &update_dblink_window;
          } );

    # It isn't necessary to shadow the border, but it looks nice :)
    $dblink_clist->set_shadow_type( 'out' );

    # What however is important, is that we set the column widths as
    # they will never be right otherwise. Note that the columns are
    # numbered from 0 and up (to 1 in this case).
    $dblink_clist->set_column_width( 0, 150 );

    # Add the CList widget to the vertical box
    $scrolled_window->add( $dblink_clist );

    my $bbox = new Gtk::HBox( $false, 5 );
    $bbox->border_width(10);
    my $arrow = new Gtk::Arrow('down','out');
    my $button = new Gtk::Button();
    $button->add($arrow);
    $bbox->pack_end( $button, $false, $false, 0);
    $button->signal_connect
	( "clicked", 
	  # next dblink button handler
	  sub { return if (!defined($seq[$current]));
		&store_current_dblink;
		$current_dblink++
		    if ($current_dblink <((scalar $seq[$current]->annotation->each_DBLink)-1));
		&update_dblink_window;
	      } );

    $arrow = new Gtk::Arrow('up','out');
    $button = new Gtk::Button();
    $button->add($arrow);
    $bbox->pack_end( $button, $false, $false, 0);
    $button->signal_connect( "clicked", 
           # prev comment button handler
	   sub { return if (!defined($seq[$current]));
		 &store_current_dblink;
		 $current_dblink--
		     if ($current_dblink > 0);
		 &update_dblink_window;
	       } );

    $button = new Gtk::Button("Add");
    $bbox->pack_start( $button, $false, $false, 0);
    $button->signal_connect( "clicked",
           # add comment button handler 
           sub { return if (!defined($seq[$current]));
		 &store_current_dblink;
		 my $dblink = new Bio::Annotation::DBLink;
		 $dblink->primary_id("<New>");
		 $seq[$current]->annotation->add_DBLink( $dblink );
		 $current_dblink = $seq[$current]->annotation->each_DBLink - 1;
		 $dblink_clist->append("","");
		 &update_dblink_window;
	       } );

     $button = new Gtk::Button("Delete");
     $bbox->pack_start( $button, $false, $false, 0);
     $button->signal_connect( "clicked", 
           # delete comment button handler
           sub { return if !defined($seq[$current]); 
		 $seq[$current]->annotation->remove_DBLink( $current_dblink );
		 $dblink_clist->remove($current_dblink);
		 $current_dblink-- if ($current_dblink > 0);
		 &update_dblink_window;
	       } );

     $dblink_window->vbox->pack_start( $bbox, $false, $false, 0);

    # horizontal box containing primary_id & optional_id entries
    my $hbox;
    $hbox = new Gtk::HBox( $true, 10 );
    $hbox->border_width( 1 );

    # text entries
    $dblink_entry[0] = new Gtk::Entry();
    my $frame = new Gtk::Frame("primary id");
    $frame->add($dblink_entry[0]);
    $hbox->pack_start( $frame, $true, $true, 0);

    $dblink_entry[1] = new Gtk::Entry();
    $frame = new Gtk::Frame("optional id");
    $frame->add($dblink_entry[1]);
    $hbox->pack_end( $frame, $true, $true, 0);

    $dblink_window->vbox->pack_start( $hbox, $false, $false, 5);

    $dblink_entry[2] = new Gtk::Entry();
    $frame = new Gtk::Frame("Database");
    $frame->add($dblink_entry[2]);
    $dblink_window->vbox->pack_start( $frame, $false, $false, 5);

    $dblink_entry[3] = new Gtk::Entry();
    $frame = new Gtk::Frame("Comment");
    $frame->add($dblink_entry[3]);
    $dblink_window->vbox->pack_end( $frame, $false, $false, 5);

     $bbox = new Gtk::HButtonBox();
     $bbox->set_layout("end");
	
     $button = new Gtk::Button( "Close" );
     $bbox->add( $button );
     $button->signal_connect("clicked",
           # close button handler
           sub{ $dblink_window->hide();
		&store_current_dblink;
	      });

     $dblink_window->action_area->pack_start( $bbox, $true, $true, 0 );
     $dblink_window->signal_connect_after( "delete_event",
           # window delete handler
           sub{ $dblink_window->hide();
		&store_current_dblink;
                return &Gtk::true;
	      });
}

sub store_current_dblink
{
    if ((defined($seq[$current])) && ($seq[$current]->annotation->each_DBLink))
    {
	(($seq[$current]->annotation->each_DBLink)[$current_dblink])->
	    primary_id($dblink_entry[0]->get_chars(0,-1) );		
	(($seq[$current]->annotation->each_DBLink)[$current_dblink])->
	    optional_id($dblink_entry[1]->get_chars(0,-1) );		
	(($seq[$current]->annotation->each_DBLink)[$current_dblink])->
	    database($dblink_entry[2]->get_chars(0,-1) );		
	(($seq[$current]->annotation->each_DBLink)[$current_dblink])->
	    comment($dblink_entry[3]->get_chars(0,-1) );		
    }
}

sub update_dblink_window
{
    my ($show_me) = @_;
    $dblink_window->show_all() if (defined($show_me));

    $dblink_clist->freeze();
    if (!defined($seq[$current]))
    {
	$dblink_clist->clear();
	$dblink_clist->thaw();
	foreach (@dblink_entry) { $_->set_text(""); }
	return;
    }
    my @dblinks = $seq[$current]->annotation->each_DBLink;
    # reset clist if rows are different to links
    if ($dblink_clist->rows != @dblinks) {
	$dblink_clist->clear();
	foreach (@dblinks) { $dblink_clist->append("",""); }
    }
    # redraw references
    for(my $i=0;$i<@dblinks;$i++)
    {
	$dblink_clist->set_text($i,0,$dblinks[$i]->primary_id);
	$dblink_clist->set_text($i,1,$dblinks[$i]->database);
    }
    # redraw text widgets
    foreach (@dblink_entry) { $_->set_text(""); }
    if (@dblinks)
    {
	$dblink_entry[0]->set_text($dblinks[$current_dblink]->primary_id);
	$dblink_entry[1]->set_text($dblinks[$current_dblink]->optional_id);
	$dblink_entry[2]->set_text($dblinks[$current_dblink]->database);
	$dblink_entry[3]->set_text($dblinks[$current_dblink]->comment);
    }

    $dblink_clist->moveto($current_dblink,0,0.3,0.0)
	if ($dblink_clist->row_is_visible($current_dblink) ne 'full');
    $dblink_clist->signal_handler_block($dblink_handler_id);
    $dblink_clist->select_row($current_dblink,0);
    $dblink_clist->signal_handler_unblock($dblink_handler_id);
    Gtk::CList::set_focus_row($dblink_clist,$current_dblink);
    $dblink_clist->thaw();
}

sub init_reference_window
{
    $current_ref = 0;
    
    $ref_window = new Gtk::Dialog();
    $ref_window->set_default_size(620,500);
    $ref_window->set_policy($true,$true,$false);
    $ref_window->set_title("References");
    $ref_window->border_width(5);
    
    # Create a scrolled window to pack the CList widget into
    my $scrolled_window = new Gtk::ScrolledWindow( undef, undef );
    $ref_window->vbox->pack_start( $scrolled_window, $true, $true, 0 );
    $scrolled_window->set_policy( 'automatic', 'always' );

    # Create the CList. For this example we use 2 columns
    $ref_clist = new_with_titles Gtk::CList( "Medline","Title","Authors" );

    # When a selection is made, we want to know about it. The callback
    # used is selection_made, and its code can be found further down
    $ref_handler_id = $ref_clist->signal_connect( "select_row", 
       sub{ return if (!defined($seq[$current]));
	    my ( $clist, $row ) = @_;
	    &store_current_reference;
	    $current_ref = $row;
	    &update_reference_window;
          } );

    # It isn't necessary to shadow the border, but it looks nice :)
    $ref_clist->set_shadow_type( 'out' );

    # What however is important, is that we set the column widths as
    # they will never be right otherwise. Note that the columns are
    # numbered from 0 and up (to 1 in this case).
    $ref_clist->set_column_width( 0, 70 );
    $ref_clist->set_column_width( 1, 350 );
    $ref_clist->set_column_width( 2, 300 );

    # Add the CList widget to the vertical box
    $scrolled_window->add( $ref_clist );

    my $bbox = new Gtk::HBox( $false, 5 );
    $bbox->border_width(10);
    my $arrow = new Gtk::Arrow('down','out');
    my $button = new Gtk::Button();
    $button->add($arrow);
    $bbox->pack_end( $button, $false, $false, 0);
    $button->signal_connect
	( "clicked", 
	  # next ref button handler
	  sub { return if (!defined($seq[$current]));
		&store_current_reference;
		$current_ref++
		    if ($current_ref <((scalar $seq[$current]->annotation->each_Reference)-1));
		&update_reference_window;
	    } );

    $arrow = new Gtk::Arrow('up','out');
    $button = new Gtk::Button();
    $button->add($arrow);
    $bbox->pack_end( $button, $false, $false, 0);
    $button->signal_connect( "clicked", 
           # prev comment button handler
	   sub { return if (!defined($seq[$current]));
		 &store_current_reference;
		 $current_ref--
		     if ($current_ref > 0);
		 &update_reference_window;
	       } );

    $button = new Gtk::Button("Add");
    $bbox->pack_start( $button, $false, $false, 0);
    $button->signal_connect( "clicked",
           # add comment button handler 
           sub { return if (!defined($seq[$current]));
		 &store_current_reference;
		 my $ref = new Bio::Annotation::Reference;
		 $ref->medline("<New>");
		 $seq[$current]->annotation->add_Reference( $ref );
		 $ref_clist->append("","","");
		 $current_ref = ($seq[$current]->annotation->each_Reference)-1;
		 &update_reference_window;
	       } );

     $button = new Gtk::Button("Delete");
     $bbox->pack_start( $button, $false, $false, 0);
     $button->signal_connect( "clicked", 
           # delete comment button handler
           sub { return if !defined($seq[$current]); 
		 $seq[$current]->annotation->remove_Reference( $current_ref );
		 $ref_clist->remove($current_ref);
		 $current_ref-- if ($current_ref > 0);
		 &update_reference_window;
	       } );

     $ref_window->vbox->pack_start( $bbox, $false, $false, 0);

    # horizontal box containing primary_id & optional_id entries
    my $hbox;
    $hbox = new Gtk::HBox( $true, 10 );
    $hbox->border_width( 1 );

    # text entries
    $ref_entry[0] = new Gtk::Entry();
    my $frame = new Gtk::Frame("Title");
    $frame->add($ref_entry[0]);
    $ref_window->vbox->pack_start( $frame, $false, $false, 5);

    $ref_entry[1] = new Gtk::Entry();
    $frame = new Gtk::Frame("Authors");
    $frame->add($ref_entry[1]);
    $ref_window->vbox->pack_start( $frame, $false, $false, 5);

    # horizontal box
    $hbox = new Gtk::HBox( $true, 10 );
    $hbox->border_width( 1 );

    # text entries
    $ref_entry[2] = new Gtk::Entry();
    $frame = new Gtk::Frame("Comment");
    $frame->add($ref_entry[2]);
    $hbox->pack_start( $frame, $true, $true, 0);

    $ref_entry[3] = new Gtk::Entry();
    $frame = new Gtk::Frame("Location");
    $frame->add($ref_entry[3]);
    $hbox->pack_end( $frame, $true, $true, 0);

    $ref_window->vbox->pack_start( $hbox, $false, $false, 5);

    # horizontal box
    $hbox = new Gtk::HBox( $false, 10 );
    $hbox->border_width( 1 );

    # text entries
    $ref_entry[4] = new Gtk::Entry();
    $frame = new Gtk::Frame("Medline");
    $frame->add($ref_entry[4]);
    $hbox->pack_start( $frame, $false, $false, 0);

#    $ref_entry[5] = new Gtk::Entry();
#    $frame = new Gtk::Frame("Start");
#    $frame->add($ref_entry[5]);
#    $hbox->pack_start( $frame, $false, $false, 0);

    $ref_entry[5] = new Gtk::Entry();
    $frame = new Gtk::Frame("Reference Position");
    $frame->add($ref_entry[5]);
    $hbox->pack_end( $frame, $true, $true, 0);

    $ref_window->vbox->pack_start( $hbox, $false, $false, 5);


     $bbox = new Gtk::HButtonBox();
     $bbox->set_layout("end");
	
     $button = new Gtk::Button( "Close" );
     $bbox->add( $button );
     $button->signal_connect("clicked",
           # close button handler
           sub{ $ref_window->hide();
		&store_current_reference;
	      });

     $ref_window->action_area->pack_start( $bbox, $true, $true, 0 );
     $ref_window->signal_connect_after( "delete_event",
           # window delete handler
           sub{ $ref_window->hide();
		&store_current_reference;
                return &Gtk::true;
	      });
}

sub store_current_reference
{
    if ((defined($seq[$current])) && ($seq[$current]->annotation->each_Reference))
    {
	(($seq[$current]->annotation->each_Reference)[$current_ref])->
	    title($ref_entry[0]->get_chars(0,-1) );		
	(($seq[$current]->annotation->each_Reference)[$current_ref])->
	    authors($ref_entry[1]->get_chars(0,-1) );		
	(($seq[$current]->annotation->each_Reference)[$current_ref])->
	    comment($ref_entry[2]->get_chars(0,-1) );		
	(($seq[$current]->annotation->each_Reference)[$current_ref])->
	    location($ref_entry[3]->get_chars(0,-1) );		
	(($seq[$current]->annotation->each_Reference)[$current_ref])->
	    medline($ref_entry[4]->get_chars(0,-1) );		
#	(($seq[$current]->annotation->each_Reference)[$current_ref])->
#	    start($ref_entry[5]->get_chars(0,-1) );		
	(($seq[$current]->annotation->each_Reference)[$current_ref])->
	    rp($ref_entry[5]->get_chars(0,-1) );		
    }
}

sub update_reference_window
{
    my ($show_me) = @_;
    $ref_window->show_all() if (defined($show_me));

    $ref_clist->freeze();
    if (!defined($seq[$current]))
    {
	$ref_clist->clear();
	$ref_clist->thaw();
	foreach (@ref_entry) { $_->set_text(""); }
	return;
    }
    my @refs = $seq[$current]->annotation->each_Reference;
    # reset clist if rows are different to references
    if ($ref_clist->rows != @refs) {
	$ref_clist->clear();
	foreach (@refs) { $ref_clist->append("","",""); }
    }
    # redraw references
    for(my $i=0;$i<@refs;$i++)
    {
	$ref_clist->set_text($i,0,$refs[$i]->medline)
	  if ($refs[$i]->medline);
	$ref_clist->set_text($i,1,$refs[$i]->title)
	  if ($refs[$i]->title);
	$ref_clist->set_text($i,2,$refs[$i]->authors)
	  if ($refs[$i]->authors);
    }
    # redraw text widgets
    foreach (@ref_entry) { $_->set_text(""); }
    if (@refs) {
	$ref_entry[0]->set_text($refs[$current_ref]->title);
	$ref_entry[1]->set_text($refs[$current_ref]->authors);
	$ref_entry[2]->set_text($refs[$current_ref]->comment);
	$ref_entry[3]->set_text($refs[$current_ref]->location);
	$ref_entry[4]->set_text($refs[$current_ref]->medline);
#	$ref_entry[5]->set_text($refs[$current_ref]->start);
	$ref_entry[5]->set_text($refs[$current_ref]->rp);
    }

    $ref_clist->moveto($current_ref,0,0.3,0.0)
	if ($ref_clist->row_is_visible($current_ref) ne 'full');
    $ref_clist->signal_handler_block($ref_handler_id);
    $ref_clist->select_row($current_ref,0);
    $ref_clist->signal_handler_unblock($ref_handler_id);
    Gtk::CList::set_focus_row($ref_clist,$current_ref);
    $ref_clist->thaw();
}


sub init_about_dialog {
 	my ($window,$bg,$tbox,$vbox,$hbox,$sep,$butbox,$button,$pixmap);
 	$about_dialog = new Gtk::Window("dialog");
 	$about_dialog->set_title("About gSequence");
 	$about_dialog->signal_connect_after("destroy" => 
					    sub { $about_dialog->hide; 
					          return &Gtk::true; });
 	$about_dialog->set_default_size('350','350');
 	$about_dialog->set_policy(1,1,0);
 	$window = $about_dialog->window;
 	$bg = $about_dialog->style->bg('normal');
 	$vbox= new Gtk::VBox(0,0);
 	$about_dialog->add($vbox);
    	$tbox = new Gtk::Label("\ngSequence\nAuthor: Lorenz Pollak\n\n
gSequence is cool! :-)\n(this text is to be written...)
\n");
 	$vbox->pack_start($tbox,1,1,1);
   
    	$hbox = new Gtk::HBox(0,0);   
    	$vbox->pack_start($hbox,0,0,0);
 	$sep = new Gtk::HSeparator;
   	$sep->set_usize(-1,5);
    	$vbox->pack_start($sep,0,1,0);
    
     	$butbox = new Gtk::HButtonBox;
     	$butbox->set_usize(-1,32);
     	$vbox->pack_start($butbox, 0,1,0); 
 	$button = new_with_label Gtk::Button("OK");
     	$button->set_usize(50,-1);
     	$button->signal_connect('clicked', sub { $about_dialog->hide; });
	$button->can_default(1);
	$button->grab_default;
	$butbox->add($button);
     
  return 1;
}    

sub init_feature_window
{
    $current_feature_item = 0;
    
    $feature_window = new Gtk::Dialog();
    $feature_window->set_default_size(500,400);
    $feature_window->set_policy($true,$true,$false);
    $feature_window->set_title("Sequence Features");
    $feature_window->border_width(5);

    my $pane = new Gtk::HPaned();
    $feature_window->vbox->pack_start( $pane, $true, $true, 0);
    $pane->set_handle_size( 10 );
    $pane->set_gutter_size( 8 );

    # Create a VBox for the Entry and Tree Scrolled Window
    my $vbox = new Gtk::VBox( $false, 0 );
    $pane->add1( $vbox );

    # Create a ScrolledWindow for the tree
    my $tree_scrolled_win = new Gtk::ScrolledWindow( undef, undef );
    $tree_scrolled_win->set_usize( 150, 400 );
    $vbox->pack_start( $tree_scrolled_win, $true, $true, 0 );
    $tree_scrolled_win->set_policy( 'automatic', 'automatic' );

    #my $list_scrolled_win = new Gtk::ScrolledWindow( undef, undef );
    #$list_scrolled_win->set_policy( 'automatic', 'automatic' );
    $vbox = new Gtk::VBox( $false, 0 );
    $pane->add2( $vbox );

    # add stuff to the vbox
    # text entries
    my $hbox = new Gtk::HBox( $true, 10 );

    $feature_entry[0] = new Gtk::Entry();
    my $frame = new Gtk::Frame("Primary Tag");
    $frame->add($feature_entry[0]);
    $hbox->pack_start( $frame, $true, $true, 0);

    $feature_entry[1] = new Gtk::Entry();
    $frame = new Gtk::Frame("Source Tag");
    $frame->add($feature_entry[1]);
    $hbox->pack_end( $frame, $true, $true, 0);

    $vbox->pack_start( $hbox, $false, $false, 5);

    $hbox = new Gtk::HBox( $true, 10 );

    my $adj = new Gtk::Adjustment( 0, 0, 0, 0, 0, 0 );
    $feature_spinner[0] = new Gtk::SpinButton( $adj, 0.0, 0 );
    $feature_spinner[0]->signal_connect( "changed", \&select_feature_region);
    $frame = new Gtk::Frame("Start");
    $frame->add($feature_spinner[0]);
    $hbox->pack_start( $frame, $true, $true, 0);

    $adj = new Gtk::Adjustment( 0, 0, 0, 0, 0, 0 );
    $feature_spinner[1] = new Gtk::SpinButton( $adj, 0.0, 0 );
    $feature_spinner[1]->signal_connect( "changed", \&select_feature_region);
    $frame = new Gtk::Frame("End");
    $frame->add($feature_spinner[1]);
    $hbox->pack_start( $frame, $true, $true, 0);

    $frame = new Gtk::Frame("Strand");
    $hbox->pack_start( $frame, $true, $true, 0);
    $frame = new Gtk::Frame("Score");
    $hbox->pack_start( $frame, $true, $true, 0);

    $vbox->pack_start( $hbox, $false, $false, 5);

    $feature_entry[2] = new Gtk::Entry();
    $frame = new Gtk::Frame("Description");
    $frame->add($feature_entry[2]);

    $vbox->pack_start( $frame, $false, $false, 5);

    my $bbox = new Gtk::HBox( $false, 5 );
    $bbox->border_width(10);
    my $button = new Gtk::Button("Add");
    $bbox->pack_start( $button, $false, $false, 0);
    $button->signal_connect( "clicked",
           # add comment button handler 
           sub { return if (!defined($seq[$current]));
		 &store_current_feature if ($current_feature_item);
		 my $feature = new Bio::SeqFeature::Generic;
		 $feature->primary_tag("<New>");
		 $seq[$current]->add_SeqFeature( $feature );
		 my $item_new = new_with_label Gtk::TreeItem( "<New>" );
		 $item_new->set_user_data( $feature );
		 $item_new->signal_connect( 'select', \&select_feature_item );
		 $current_feature_item->parent->append( $item_new )
		     if ($current_feature_item);
		 $feature_tree->append( $item_new ) if (!$current_feature_item);
		 $item_new->show();
		 $current_feature_item->deselect()
		     if ($current_feature_item);
		 $item_new->select();
	       } );
    $button = new Gtk::Button("Add Subfeature");
    $bbox->pack_start( $button, $false, $false, 0);
    $button->signal_connect( "clicked",
           # add comment button handler 
           sub { return if (!defined($seq[$current])||!$current_feature_item);
		 &store_current_feature;
		 my $feature = new Bio::SeqFeature::Generic;
		 $feature->primary_tag("<New>");
		 $feature->start($current_feature_item->get_user_data->start);
		 $feature->end($current_feature_item->get_user_data->end);
		 $current_feature_item->get_user_data->add_sub_SeqFeature( $feature );
		 my $new_subtree = new Gtk::Tree();
		 $current_feature_item->set_subtree( $new_subtree );
		 my $item_new = new_with_label Gtk::TreeItem( "<New>" );
		 $item_new->set_user_data( $feature );
		 $item_new->signal_connect( 'select', \&select_feature_item );
		 $new_subtree->append( $item_new );
		 $item_new->show();
		 $current_feature_item->deselect();
		 $current_feature_item->expand();
		 $item_new->select();
	       } );
     $button = new Gtk::Button("Delete");
     $bbox->pack_start( $button, $false, $false, 0);
     $button->signal_connect( "clicked", 
           # delete comment button handler
           sub { return if (!$current_feature_item); 
		 &store_current_feature;
		 my $flist = $seq[$current]->{_as_feat};
		 my $pos;
		 for(my $i=0;$i<@$flist;$i++) {
		   $pos=$i if $$flist[$i]==$current_feature_item->get_user_data();
		 }
		 splice @$flist, $pos, 1;
		 $seq[$current]->{_as_feat} = $flist;
		 $current_feature_item->parent->remove_item($current_feature_item);
		 $current_feature_item=0;
	       } );

     $vbox->pack_end( $bbox, $false, $false, 0);

    # Create root tree
    $feature_tree = new Gtk::Tree();
    $tree_scrolled_win->add_with_viewport( $feature_tree );
    $feature_tree->set_selection_mode( 'single' );
    $feature_tree->set_view_mode( 'item' );

    $bbox = new Gtk::HButtonBox();
    $bbox->set_layout("end");
    
    $button = new Gtk::Button( "Close" );
    $bbox->add( $button );
    $button->signal_connect("clicked",
			    # close button handler
			    sub{ $feature_window->hide();
				 &store_current_feature;
			     });
    
    $feature_window->action_area->pack_start( $bbox, $true, $true, 0 );
    $feature_window->signal_connect_after( "delete_event",
					   # window delete handler
					   sub{ $feature_window->hide();
						&store_current_feature;
						return &Gtk::true;
					    });
}

# Callback for expanding tree
sub expand_feature_tree
  {
    my ( $item, $subtree ) = @_;
    my ($feature,$subfeature,$item_new,$new_subtree);
    $feature = $item->get_user_data();

    foreach $subfeature ($feature->sub_SeqFeature)
      {
	  $item_new = new_with_label Gtk::TreeItem( $subfeature->primary_tag );
	  $item_new->set_user_data( $subfeature );
	  $item_new->signal_connect( 'select', \&select_feature_item );
	  $subtree->append( $item_new );
	  $item_new->show();
	  
	  if ( $subfeature->sub_SeqFeature )
	  {
	      $new_subtree = new Gtk::Tree();
	      $item_new->set_subtree( $new_subtree );
	      $item_new->signal_connect( 'expand',
					 \&expand_feature_tree,
					 $new_subtree );
	      $item_new->signal_connect( 'collapse', \&collapse_feature_tree );
	  }
	  $item_new->expand();
      }
  }


# Callback for collapsing tree
sub collapse_feature_tree
  {
    my ( $item ) = @_;

    my $subtree = new Gtk::Tree();

    $item->remove_subtree();
    $item->set_subtree( $subtree );
    $item->signal_connect( 'expand', \&expand_feature_tree, $subtree );
  }


sub store_current_feature
{
  if ((defined($seq[$current])) && ($seq[$current]->top_SeqFeatures) && ($current_feature_item))
  {
    my $current_feature = $current_feature_item->get_user_data();
    $current_feature->primary_tag( $feature_entry[0]->get_chars(0,-1) );		
    $current_feature->source_tag( $feature_entry[1]->get_chars(0,-1) );		
    if ($current_feature->has_tag("description"))
    {
      $current_feature->remove_tag("description");
      $current_feature->add_tag_value("description",
				      $feature_entry[2]->get_chars(0,-1));
    }
    $current_feature->start($feature_spinner[0]->get_value_as_int());
    $current_feature->end($feature_spinner[1]->get_value_as_int());
    # set tree item
    ($current_feature_item->children)[0]->set($current_feature->primary_tag);
  }
}

sub select_feature_item
{
    my ($widget) = @_;
    &store_current_feature;
    $current_feature_item->deselect()
      if $current_feature_item;
    $current_feature_item = $widget;
    &update_feature_paned2;
}

sub update_feature_paned2
{
  $feature_entry[0]->set_text("");
  $feature_entry[1]->set_text("");
  $feature_entry[2]->set_text("");

  return if (!defined($seq[$current])||(!$current_feature_item));
  my $current_feature = $current_feature_item->get_user_data();
  $feature_entry[0]->set_text($current_feature->primary_tag);
  $feature_entry[1]->set_text($current_feature->source_tag)
    if (defined($current_feature->source_tag));
  $feature_entry[2]->set_text(($current_feature->each_tag_value("description"))[0])
    if ($current_feature->has_tag("description"));
  my $adj = new Gtk::Adjustment($current_feature->start,
				0,
				$seq[$current]->length-1,
				1,
				1,
				0
			       );
  $feature_spinner[0]->set_adjustment($adj);
  $feature_spinner[0]->set_value($current_feature->start);
  $feature_spinner[0]->show_all();
  $adj = new Gtk::Adjustment($current_feature->end,
			     0,
			     $seq[$current]->length-1,
			     1,
			     1,
			     0
			    );
  $feature_spinner[1]->set_adjustment($adj);
  $feature_spinner[1]->set_value($current_feature->end);
  $feature_spinner[1]->show_all();
}

sub select_feature_region
{
  $seq_edit[$current]->freeze;
  $seq_edit[$current]->select_region($feature_spinner[0]->get_value_as_int(),
			   $feature_spinner[1]->get_value_as_int()+1);
  $seq_edit[$current]->thaw;
}

sub update_feature_window
{
    my ($show_me) = @_;
    $feature_window->show_all() if (defined($show_me));

    $feature_tree->clear_items(0,-1);
    if (!defined($seq[$current]))
    {
	&update_feature_paned2;
	return;
    }

    my ($item_new,$new_subtree);
    foreach ($seq[$current]->top_SeqFeatures)
      {
	  $item_new = new_with_label Gtk::TreeItem( $_->primary_tag );
	  $item_new->set_user_data( $_ );
	  $item_new->signal_connect( 'select', \&select_feature_item );
	  $feature_tree->append( $item_new );
	  if ( $_->sub_SeqFeature )
	  {
	      $new_subtree = new Gtk::Tree();
	      $item_new->set_subtree( $new_subtree );
	      $item_new->signal_connect( 'expand',
					 \&expand_feature_tree,
					 $new_subtree );
	      $item_new->signal_connect( 'collapse', \&collapse_feature_tree );
	  }
	  $item_new->expand();
      }
    $feature_tree->select_item($current_feature_item) 
      if $current_feature_item;
    $feature_tree->show_all();

    &update_feature_paned2;
}

sub store_prefs
{
}

sub update_pref_window
{
  $pref_window->show_all();
}

sub init_pref_window
{
  $pref_window = new Gtk::Dialog();
  $pref_window->set_default_size(500,400);
  $pref_window->set_policy($true,$true,$false);
  $pref_window->border_width( 5 );

  # Create a new notebook, place the position of the tabs
  my $notebook = new Gtk::Notebook();
  $pref_window->vbox->pack_start( $notebook, $true, $true, 0);
  $notebook->set_tab_pos( 'top' );

  my $main_vbox = new Gtk::VBox($false,10);

  my $label = new Gtk::Label( "Import Options" );
  my $frame = new Gtk::Frame("Flat File Indexes");
  my $vbox = new Gtk::VBox($false,10);
  $frame->add($vbox);
  $main_vbox->pack_start($frame,$false,$false,10);

  $notebook->append_page( $main_vbox, $label );

  my $hbox = new Gtk::HBox($false,0);

  $pref_entry[0] = new Gtk::Entry();
  $frame = new Gtk::Frame("Indexes Directory");
  $frame->add($pref_entry[0]);
  $hbox->pack_start( $frame, $true, $false, 0);

  $pref_entry[1] = new Gtk::Entry();
  $frame = new Gtk::Frame("Index Type");
  $frame->add($pref_entry[1]);
  $hbox->pack_start( $frame, $false, $false, 0);

  $vbox->pack_start( $hbox, $false, $false, 0);

  $pref_entry[2] = new Gtk::Entry();
  $frame = new Gtk::Frame("Fasta Index Name");
  $frame->add($pref_entry[2]);
  $vbox->pack_start( $frame, $false, $false, 0);

  $pref_entry[3] = new Gtk::Entry();
  $frame = new Gtk::Frame("SwissProt Index Name");
  $frame->add($pref_entry[3]);
  $vbox->pack_start( $frame, $false, $false, 0);

  $pref_entry[4] = new Gtk::Entry();
  $frame = new Gtk::Frame("SwissPfam Index Name");
  $frame->add($pref_entry[4]);
  $vbox->pack_start( $frame, $false, $false, 0);

  $frame = new Gtk::Frame("Remote DBs");
  $hbox = new Gtk::HBox($false,10);
  $frame->add($hbox);
  $main_vbox->pack_start($frame,$false,$false,10);

  $pref_entry[5] = new Gtk::Entry();
  $frame = new Gtk::Frame("AceDB host");
  $frame->add($pref_entry[5]);
  $hbox->pack_start( $frame, $true, $false, 0);

  $pref_entry[6] = new Gtk::Entry();
  $frame = new Gtk::Frame("AceDB port");
  $frame->add($pref_entry[6]);
  $hbox->pack_start( $frame, $false, $false, 0);

  $notebook->set_page( 0 );

  my $bbox = new Gtk::HButtonBox();
  $bbox->set_layout("end");

  my $button = new Gtk::Button( "Save" );
  $bbox->add( $button );
  $button->signal_connect("clicked",
			  # close button handler
			  sub{ $pref_window->hide();
			       &store_prefs();
			     });
  
  $button = new Gtk::Button( "Close" );
  $bbox->add( $button );
  $button->signal_connect("clicked",
			  # close button handler
			  sub{ $pref_window->hide();
			     });
  
  $pref_window->action_area->pack_start( $bbox, $true, $true, 0 );
  $pref_window->signal_connect_after( "delete_event",
					 # window delete handler
					 sub{ $pref_window->hide();
					      return &Gtk::true;
					    });
}
