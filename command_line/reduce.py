# LIBTBX_SET_DISPATCHER_NAME phenix.reduce
# LIBTBX_SET_DISPATCHER_NAME mmtbx.reduce
# LIBTBX_SET_DISPATCHER_NAME molprobity.reduce
import libtbx.load_env
import os, sys, subprocess

def raise_sorry(message, file_name):
  from libtbx.str_utils import show_string
  from libtbx.utils import Sorry
  raise Sorry(message % show_string(file_name))

def run(args):
  reduce_exe = libtbx.env.under_build("reduce/exe/reduce")
  if (os.name == "nt"):
    reduce_exe += ".exe"
  if (not os.path.isfile(reduce_exe)):
    raise_sorry("Missing phenix.reduce executable: %s", reduce_exe)
  db = libtbx.env.under_dist("reduce", "reduce_wwPDB_het_dict.txt")
  if (not os.path.isfile(db)):
    raise_sorry("Missing phenix.reduce database: %s", db)
  cmd = [reduce_exe, "-DB", db, "-ALLALT"]
  subprocess.call(cmd + args)

if (__name__ == "__main__"):
  run(sys.argv[1:])
