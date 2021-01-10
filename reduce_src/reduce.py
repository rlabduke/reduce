##################################################################################
# This program is a rewrite of a subset of the main.cpp Reduce program into
# Python that makes use of the shared library that is produced using Boost.Python.
# It is a test program to validate that the Python wrapping worked.
#

import sys
import mmtbx_reduce_ext as reduce

def RunReduce(input, hetdatabase):
  ret = 0

  models = reduce.inputModels(input)
  for m in models:
    if reduce.getRemoveATOMHydrogens() or reduce.getRemoveOtherHydrogens():
      reduce.dropHydrogens(m, reduce.getRemoveATOMHydrogens(), reduce.getRemoveOtherHydrogens())

    UseSEGIDasChain = reduce.checkSEGIDs(m)

    dotBucket = reduce.DotSphManager(reduce.getVdwDotDensity())

    xyz = reduce.AtomPositions(2000, reduce.getDoOnlyAltA(), reduce.getUseXplorNames(),
            reduce.getUseOldNames(), reduce.getBackBoneModel(),
            reduce.getNBondCutoff(), reduce.getMinRegHBgap(),
            reduce.getMinChargedHBgap(),
            reduce.getBadBumpGapCut(), dotBucket,
            reduce.getProbeRadius(),
            reduce.getPenaltyMagnitude(), reduce.getOccupancyCutoff(),
            reduce.getVerbose(), reduce.getShowOrientScore(),
            reduce.getShowCliqueTicks()
          )

    infoPtr = m.begin()

    reduce.scanAndGroupRecords(m, xyz, infoPtr)

    adjNotes = reduce.StringVector()

    tally = reduce.getTally()
    tally._num_adj = 0
    reduce.setTally(tally)

    reduce.reduceList(hetdatabase, m, xyz, adjNotes)

    if not reduce.getStopBeforeOptimizing():
      ret = reduce.optimize(xyz, adjNotes)

      if reduce.getOKtoAdjust() and xyz.numChanges() > 0:
        xyz.describeChanges(m, infoPtr, adjNotes)

  reduce.outputRecords_all(models)

  return ret

if __name__ == '__main__':

  #==============================================================
  # Parse command-line arguments.  The 0th argument is the name
  # of the script. There can be -TRIM or -FLIP as an optional argument and
  # then the name of a PDB file to read.
  realParams = 0
  fileName = ""
  for i in range(1,len(sys.argv)):
    if sys.argv[i] == "-TRIM":
      reduce.setRemoveATOMHydrogens(True)
      reduce.setRemoveOtherHydrogens(True)
      reduce.setAddWaterHydrogens(False)
      reduce.setAddOtherHydrogens(False)
      reduce.setStopBeforeOptimizing(True)
    elif sys.argv[i] == "-FLIP":
      reduce.setBuildHisHydrogens(True)
      reduce.setSaveOHetcHydrogens(True)
      reduce.setRotExistingOH(True)
      reduce.setDemandFlipAllHNQs(True)
    else:
      fileName = sys.argv[i]

  if len(fileName) == 0:
    print('Usage:', sys.argv[0], "[-TRIM] PDB_FILE_NAME")
    sys.exit(-2)

  with open(fileName) as f:
    input = f.read()

  # @todo Point this at the general location.
  hetdatabase = reduce.CTab("/usr/local/reduce_wwPDB_het_dict.txt")

  sys.exit(RunReduce(input, hetdatabase))
