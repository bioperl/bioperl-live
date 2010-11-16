package Bio::Tree::LabeledNodeI;
use strict;

use base qw(Bio::Root::RootI);

=head2 id

 Title   : id
 Usage   : $node->id($new_id_string);
 Function: Get/set the identifier / label for the node.
 Returns : Value of the node's identifier
 Args    : String, new value for the node's ID

=cut

sub id { shift->throw_not_implemented() }
sub name { shift->id(@_) }
sub label { shift->id(@_) }

1;
