# Let the code begin...
package Bio::Expression::FeatureGroup::FeatureGroupMas50;

use strict;

use base qw(Bio::Expression::FeatureGroup);
use vars qw($DEBUG);

use Class::MakeMethods::Emulator::MethodMaker
  get_set => [qw(
  
  probe_set_name stat_pairs stat_pairs_used
  signal detection detection_p_value
  stat_common_pairs signal_log_ratio
  signal_log_ratio_low
  signal_log_ratio_high change change_p_value
  positive negative pairs pairs_used
  pairs_inavg pos_fraction log_avg
  pos_neg avg_diff abs_call inc dec
  inc_ratio dec_ratio pos_change
  neg_change inc_dec dpos_dneg_ratio
  log_avg_ratio_change diff_call
  avg_diff_change b_a fold_change
  sort_score		 

  )],
;

1;
