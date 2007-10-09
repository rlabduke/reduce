from libtbx.utils import warn_if_unexpected_md5_hexdigest
import os

if (self.env.is_ready_for_build()
    and self.env.dist_path("iotbx", default=None) is not None):
  warn_if_unexpected_md5_hexdigest(
    path=self.env.under_dist(
      module_name="iotbx", path="include/iotbx/pdb/hybrid_36_c.c"),
    expected_md5_hexdigests=[
      "e96f85e8ab863dd241893941c427a5aa", # SVN revision 6312
    ],
    hints=[
      "  Files to review:",
      "    iotbx/include/iotbx/pdb/hybrid_36_c.c",
      "    reduce_src/hybrid_36_c.c"])
