# LIBTBX_SET_DISPATCHER_NAME phenix.trim_pdb

import libtbx.phil.command_line
from libtbx import easy_run
from libtbx.utils import Sorry
import os.path
import sys

master_phil = libtbx.phil.parse("""
trim_pdb
  .caption = This utility will run phenix.reduce to remove all hydrogen \
    atoms.  We strongly recommend against depositing files manipulated in \
    this manner to the PDB, as it removes information required to reproduce \
    the published R-factors.  However, at intermediate stages of refinement \
    it may be necessary to strip hydrogens due to program incompatibilities.
  .style = auto_align caption_img:icons/custom/phenix.pdbtools.png
{
  pdb_in = None
    .type = path
    .short_caption = Input PDB file
    .style = bold file_type:pdb
  pdb_out = None
    .type = path
    .short_caption = Output file
    .style = bold file_type:pdb new_file
}
""")

def run (args=(), params=None, out=sys.stdout) :
  assert (len(args) > 0) or (params is not None)
  if (params is None) :
    interpreter = libtbx.phil.command_line.argument_interpreter(
      master_phil=master_phil,
      home_scope="")
    arg_phil = []
    for arg in args :
      if os.path.isfile(arg) :
        arg_phil.append(interpreter.process(arg="pdb_in=\"%s\"" %
          os.path.abspath(arg)))
      else :
        try :
          arg_phil.append(interpreter.process(arg=arg))
        except RuntimeError, e :
          raise Sorry("Parser error at '%s': %s" % (arg, str(e)))
    params = master_phil.fetch(sources=arg_phil).extract()
  validate_params(params)
  pdb_in = params.trim_pdb.pdb_in
  pdb_out = params.trim_pdb.pdb_out
  if (pdb_out is None) :
    pdb_out = os.path.join(os.getcwd(),
        os.path.splitext(os.path.basename(pdb_in))[0] + "_no_h.pdb")
  if os.path.exists(pdb_out) :
    os.remove(pdb_out)
  easy_run.call("phenix.reduce -Trim %s > %s" % (pdb_in, pdb_out))
  if os.path.exists(pdb_out) :
    print "Wrote %s" % pdb_out
    return pdb_out
  else :
    return None

def validate_params (params) :
  if (params.trim_pdb.pdb_in is None) :
    raise Sorry("No input file specified!")
  return True

if (__name__ == "__main__") :
  run(sys.argv[1:])
