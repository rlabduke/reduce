import libtbx
from libtbx.utils import warn_if_unexpected_md5_hexdigest
import os

if (self.env.is_ready_for_build()
    and self.env.dist_path("iotbx", default=None) is not None):
  if (not hasattr(libtbx, "manual_date_stamp")):
    path = "include/iotbx/pdb/hybrid_36_c.c"
  else:
    path = "pdb/hybrid_36_c.c"
  warn_if_unexpected_md5_hexdigest(
    path=self.env.under_dist(module_name="iotbx", path=path),
    expected_md5_hexdigests=[
      "e96f85e8ab863dd241893941c427a5aa", # SVN revision 6312
    ],
    hints=[
      "  Files to review:",
      "    iotbx/pdb/hybrid_36_c.c",
      "    reduce_src/hybrid_36_c.c"])
