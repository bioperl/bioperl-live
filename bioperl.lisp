

(defun bioperl-method (method-name)
  "puts in a bioperl method complete with pod bioler-plate"
  (interactive "smethod name:")
  (insert "=head2 " method-name "\n\n Title   : " method-name "\n Usage   :\n Function:\n Example :\n Returns : \n Args    :\n\n\n=cut\n\n")
  (insert "sub " method-name "{\n   my ($self,@args) = @_;\n")
  (save-excursion 
    (insert "\n\n}\n"))
  )

(defun bioperl-object-start (perl-object-name perl-caretaker-name caretaker-email)
  "Places standard bioperl object notation headers and footers"
  (interactive "sName of Object: \nsName of caretaker: \nsEmail: ")
  (insert "\n#\n# BioPerl module for " perl-object-name "\n#\n# Cared for by " perl-caretaker-name " <" caretaker-email ">\n#\n# Copyright " perl-caretaker-name "\n#\n# You may distribute this module under the same terms as perl itself\n\n")
  (insert "# POD documentation - main docs before the code\n\n")
  (insert "=head1 NAME\n\n" perl-object-name " - DESCRIPTION of Object\n\n")
  (insert "=head1 SYNOPSIS\n\nGive standard usage here\n\n")
  (insert "=head1 DESCRIPTION\n\nDescribe the object here\n\n")
  (insert "=head1 CONTACT\n\nDescribe contact details here\n\n")
  (insert "=head1 APPENDIX\n\nThe rest of the documentation details each of the object methods. Internal methods are usually preceded with a _\n\n=cut\n\n")
  (insert "\n# Let the code begin...\n\n")
  (insert "\npackage " perl-object-name ";\n")
  (insert "use vars qw(@ISA);\n")
  (insert "use strict;\n")
  (insert "\n# Object preamble - inheriets from Bio::Root::Object\n")
  (insert "\nuse Bio::Root::Object;\n\n")
  (insert "\n@ISA = qw(Bio::Root::Object Exporter);\n")
  (insert "# new() is inherited from Bio::Root::Object\n\n")
  (insert "# _initialize is where the heavy stuff will happen when new is called\n\n")
  (insert "sub _initialize {\n  my($self,@args) = @_;\n\n  my $make = $self->SUPER::_initialize;\n\n# set stuff in self from @args\n return $make; # success - we hope!\n}\n")
  )




