from libtbx import easy_run
from libtbx.utils import Usage, file_size
import random
import sys, os

def run(args):
  if (len(args) != 1):
    raise Usage("python run_tests.py directory_with_pdb_files")
  pdb_dir = args[0]
  pdbs = os.listdir(pdb_dir)
  pdbs.sort()
  i_pdb = 0
  random.seed(0)
  while True:
    i_pdb += random.randrange(1000)
    if (i_pdb > len(pdbs)):
      break
    pdb = os.path.join(pdb_dir, pdbs[i_pdb])
    if (not os.path.isfile(pdb)): continue
    if (file_size(file_name=pdb) > 1000000): continue
    cmd = "reduce/exe/reduce " + pdb
    print cmd
    results = easy_run.fully_buffered(command=cmd)
    results.show_stderr()

if (__name__ == "__main__"):
  run(sys.argv[1:])
