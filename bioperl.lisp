
;; Perl mode set up

(assoc "\\.pl$" auto-mode-alist)
(setq auto-mode-alist (cons '("\\.pl$" . perl-mode) auto-mode-alist))

(assoc "\\.pm$" auto-mode-alist)
(setq auto-mode-alist (cons '("\\.pm$" . perl-mode) auto-mode-alist))

(defun perl-insert-start ()
  "Places #!..perl at the start of the script"
  (interactive)
  (goto-char (point-min))
  (insert "#!/usr/local/bin/perl\n"))



(defun bioperl-object-start (perl-object-name perl-caretaker-name caretaker-email)
  "Places standard bioperl object notation headers and footers"
  (interactive "sName of Object: \nsName of caretaker: \nsEmail: ")
  (insert "\n#\n# BioPerl module for " perl-object-name "\n#\n# Cared for by " perl-caretaker-name " <" caretaker-email ">\n#\n# Copyright " perl-caretaker-name "\n#\n# You may distribute this module under the same terms as perl itself\n\n")
  (insert "# POD documentation - main docs before the code\n\n")
  (insert "=head1 NAME\n\n" perl-object-name " - DESCRIPTION of Object\n\n")
  (insert "=head1 SYNOPSIS\n\nGive standard usage here\n\n")
  (insert "=head1 DESCRIPTION\n\nDescribe the object here\n\n")
  (insert "=head1 FEEDBACK\n\n=head2 Mailing Lists\n\n")
  (insert "User feedback is an integral part of the evolution of this\nand other Bioperl modules. Send your comments and suggestions preferably\n to one of the Bioperl mailing lists.\nYour participation is much appreciated.\n\n")
  (insert "  bioperl-l@bioperl.org          - General discussion\n  bioperl-guts-l@bioperl.org     - Technically-oriented discussion\n  http://bioperl.org/MailList.shtml             - About the mailing lists\n\n")
  (insert "=head2 Reporting Bugs\n\nReport bugs to the Bioperl bug tracking system to help us keep track\n the bugs and their resolution.\n Bug reports can be submitted via email or the web:\n\n")
  (insert "  bioperl-bugs@bioperl.org\n  http://bioperl.org/bioperl-bugs/\n\n")
  (insert "=head1 AUTHOR - " perl-caretaker-name "\n\nEmail " caretaker-email "\n\nDescribe contact details here\n\n")
  (insert "=head1 APPENDIX\n\nThe rest of the documentation details each of the object methods. Internal methods are usually preceded with a _\n\n=cut\n\n")
  (insert "\n# Let the code begin...\n\n")
  (insert "\npackage " perl-object-name ";\n")
  (insert "use vars qw($AUTOLOAD @ISA);\n")
  (insert "use strict;\n")
  (insert "\n# Object preamble - inherits from Bio::Root::Object\n")
  (insert "\nuse Bio::Root::Object;\n\n")
  (insert "\nuse AutoLoader;\n@ISA = qw(Bio::Root::Object Exporter);\n@EXPORT_OK = qw();\n")
  (insert "# new() is inherited from Bio::Root::Object\n\n")
  (insert "# _initialize is where the heavy stuff will happen when new is called\n\n")
  (insert "sub _initialize {\n  my($self,@args) = @_;\n\n  my $make = $self->SUPER::_initialize;\n\n# set stuff in self from @args\n return $make; # success - we hope!\n}\n")
  )

(defun bioperl-method (method-name)
  "puts in a bioperl method complete with pod boiler-plate"
  (interactive "smethod name:")
  (insert "=head2 " method-name "\n\n Title   : " method-name "\n Usage   :\n Function:\n Example :\n Returns : \n Args    :\n\n\n=cut\n\n")
  (insert "sub " method-name "{\n   my ($self,@args) = @_;\n")
  (save-excursion 
    (insert "\n\n}\n"))
  )


(defun bioperl-getset (field-name)
  "puts in a bioperl method for a get/set method complete with pod boiler-plate"
  (interactive "sfield name:")
  (insert "=head2 " field-name "\n\n Title   : " field-name "\n Usage   : $obj->" field-name "($newval)\n Function: \n Example : \n Returns : value of " field-name "\n Args    : newvalue (optional)\n\n\n=cut\n\n")
  (insert "sub " field-name "{\n   my ($obj,$value) = @_;\n   if( defined $value) {\n      $obj->{'" field-name "'} = $value;\n    }\n    return $obj->{'" field-name "'};\n")
  (insert "\n}\n"))

(setq perl-mode-hook 
      '(lambda ()
	 (define-key perl-mode-map "\C-c\C-h" 'perl-insert-start)
	 (define-key perl-mode-map "\C-c\C-b" 'bioperl-object-start)
	 (define-key perl-mode-map "\C-c\C-v" 'bioperl-getset)
	 (define-key perl-mode-map "\C-c\C-b" 'bioperl-method)
	 (define-key perl-mode-map "\C-c\C-z" 'compile)
	 (define-key perl-mode-map [menu-bar] (make-sparse-keymap))
	 (define-key perl-mode-map [menu-bar p]
	   (cons "BioPerl" (make-sparse-keymap "BioPerl")))
	 (define-key perl-mode-map [menu-bar p perl-script-start]
	   '("Insert script template" . perl-script-start))
	 (define-key perl-mode-map [menu-bar p bioperl-object-start]
	   '("bioperl object template" . bioperl-object-start))
	 (define-key perl-mode-map [menu-bar p bioperl-getset]
	   '("bioperl field func" . bioperl-getset))
	 (define-key perl-mode-map [menu-bar p bioperl-method]
	   '("bioperl method" . bioperl-method))
	 ))











