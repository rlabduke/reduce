#include <boost/python.hpp>
#include "reduce.h"

using namespace boost::python;

// Defines a getter and setter for a named global variable of the specified type
#define getSet(type, name) type get ## name(void) { return name; } void set ## name(type val) { name = val; }

// Create getters and setters for the global variables that control behavior
getSet(bool,Verbose)
getSet(bool,KeepConnections)
getSet(bool,StandardizeRHBondLengths)
getSet(bool,ProcessConnHydOnHets)
getSet(bool,BuildHisHydrogens)
getSet(bool,SaveOHetcHydrogens)
getSet(bool,UseXplorNames)
getSet(bool,UseOldNames)
getSet(bool,BackBoneModel)
getSet(bool,DemandRotAllMethyls)
getSet(bool,RotExistingOH)
getSet(bool,NeutralTermini)
getSet(bool,DemandRotNH3)
getSet(bool,DemandRotExisting)
getSet(bool,DemandFlipAllHNQs)
getSet(bool,DoOnlyAltA)
getSet(bool,OKProcessMetMe)
getSet(bool,OKtoAdjust)
getSet(bool,ShowCliqueTicks)
getSet(bool,ShowOrientScore)
getSet(bool,StringInput)
getSet(bool,ShowCharges)
getSet(bool,UseNuclearDistances)
getSet(bool,UseSEGIDasChain)
getSet(bool,ProcessedFirstModel)
getSet(bool,RenameFlip)
getSet(bool,GenerateFinalFlip)

getSet(int,MaxAromRingDih)
getSet(int,NBondCutoff)
getSet(int,ExhaustiveLimit)

getSet(float,ProbeRadius)
getSet(float,VdwDotDensity)
getSet(float,OccupancyCutoff)
getSet(float,WaterBcutoff)
getSet(float,WaterOCCcutoff)
getSet(float,PenaltyMagnitude)
getSet(float,MinRegHBgap)
getSet(float,MinChargedHBgap)
getSet(float,BadBumpGapCut)
getSet(float,NonMetalBumpBias)
getSet(float,MetalBumpBias)
getSet(float,GapWidth)

getSet(std::string,OFile)
getSet(bool,UseSEGIDtoChainMap)
getSet(bool,StopBeforeOptimizing)
getSet(bool,AddWaterHydrogens)
getSet(bool,AddOtherHydrogens)
getSet(bool,RemoveATOMHydrogens)
getSet(bool,RemoveOtherHydrogens)

BOOST_PYTHON_MODULE(reduce)
{

  // Export the const char * variables, which will be read only.
  scope().attr("versionString") = versionString;
  scope().attr("shortVersion") = shortVersion;
  scope().attr("referenceString") = referenceString;
  scope().attr("electronicReference") = electronicReference;

  // Exports a getter and setter for a named global variable of the specified type
  #define exportGetSet(name) def("get" #name , get ## name, "Get value of " #name " global variable");\
                             def("set" #name , set ## name, "Set value of " #name " global variable");

  // Export getters and setters for the global variables that control behavior
  exportGetSet(Verbose);
  exportGetSet(KeepConnections);
  exportGetSet(StandardizeRHBondLengths)
  exportGetSet(ProcessConnHydOnHets)
  exportGetSet(BuildHisHydrogens)
	exportGetSet(SaveOHetcHydrogens)
	exportGetSet(UseXplorNames)
	exportGetSet(UseOldNames)
	exportGetSet(BackBoneModel)
	exportGetSet(DemandRotAllMethyls)
	exportGetSet(RotExistingOH)
	exportGetSet(NeutralTermini)
	exportGetSet(DemandRotNH3)
	exportGetSet(DemandRotExisting)
	exportGetSet(DemandFlipAllHNQs)
	exportGetSet(DoOnlyAltA)
	exportGetSet(OKProcessMetMe)
	exportGetSet(OKtoAdjust)
	exportGetSet(ShowCliqueTicks)
	exportGetSet(ShowOrientScore)
	exportGetSet(StringInput)
	exportGetSet(ShowCharges)
	exportGetSet(UseNuclearDistances)
	exportGetSet(UseSEGIDasChain)
	exportGetSet(ProcessedFirstModel)
	exportGetSet(RenameFlip)
	exportGetSet(GenerateFinalFlip)

	exportGetSet(MaxAromRingDih)
	exportGetSet(NBondCutoff)
	exportGetSet(ExhaustiveLimit)

	exportGetSet(ProbeRadius)
	exportGetSet(VdwDotDensity)
	exportGetSet(OccupancyCutoff)
	exportGetSet(WaterBcutoff)
	exportGetSet(WaterOCCcutoff)
	exportGetSet(PenaltyMagnitude)
	exportGetSet(MinRegHBgap)
	exportGetSet(MinChargedHBgap)
	exportGetSet(BadBumpGapCut)
	exportGetSet(NonMetalBumpBias)
	exportGetSet(MetalBumpBias)
	exportGetSet(GapWidth)

	exportGetSet(OFile)
	exportGetSet(UseSEGIDtoChainMap)
	exportGetSet(StopBeforeOptimizing)
	exportGetSet(AddWaterHydrogens)
	exportGetSet(AddOtherHydrogens)
	exportGetSet(RemoveATOMHydrogens)
	exportGetSet(RemoveOtherHydrogens)

  // Export the functions that will be called
  /// @todo

}

