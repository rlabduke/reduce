##################################################################################
# This program is a rewrite of a subset of the main.cpp Reduce program into
# Python that makes use of the python.so (or python.dll) shared library that
# is produced by Boost.Python.
#

# @todo Add a main function and parse the command line
# @todo Make Python 2 and 3 compatible

import sys
import reduce

fileName = sys.argv[1]

with open(fileName) as f:
  input = f.read()

models = reduce.inputModels(input)
# print("Found "+str(len(models))+" models:")
for m in models:
  # print(" Model size "+str(m.size()))
  if reduce.getRemoveATOMHydrogens() or reduce.getRemoveOtherHydrogens():
    reduce.dropHydrogens(m, reduce.getRemoveATOMHydrogens(), reduce.getRemoveOtherHydrogens())

  UseSEGIDasChain = reduce.checkSEGIDs(m)
  # print(" UseSEGIDasChain = "+str(UseSEGIDasChain))

  # @todo Point this at the general location.
  hetdatabase = reduce.CTab("/usr/local/reduce_wwPDB_het_dict.txt")

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

sys.exit(ret)