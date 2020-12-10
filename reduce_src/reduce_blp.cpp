#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "reduce.h"

using namespace boost::python;

// Macro defines a getter and setter for a named global variable of the specified type
#define getSet(type, name) type get ## name(void) { return name; } void set ## name(type val) { name = val; }

// Create getters and setters for the global variables that control behavior
getSet(SummaryStats,Tally)

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

// Macro defines a getter and setter for a named global variable of the specified type
#define getSetMethod(class, name, type) type (class::*get ## name)() const = &class::name; void (class::*set ## name)(type) = &class::name;
#define getSetMethodTakeConst(class, name, type) type (class::*get ## name)() const = &class::name; void (class::*set ## name)(const type &) = &class::name;
#define getSetMethodReturnTakeConst(class, name, type) const type & (class::*get ## name)() const = &class::name; void (class::*set ## name)(const type &) = &class::name;

// Function pointers needed to handle overloaded functions
/*
getSetMethod(PDBrec, x, Coord);
getSetMethod(PDBrec, y, Coord);
getSetMethod(PDBrec, z, Coord);
getSetMethodTakeConst(PDBrec, loc, Point3d);
getSetMethodReturnTakeConst(PDBrec, elem, ElementInfo);
*/

// Describe and name compound classes that we need access to.
typedef std::list< std::shared_ptr<PDBrec> > ModelList;
typedef std::vector< ModelList > ModelsVector;

BOOST_PYTHON_MODULE(reduce)
{
  // Export the class objects that Python will need access to,
  // along with their shared-pointer and vector types.
  class_<SummaryStats>("SummaryStats", init<>())
    .def_readonly("_H_found", &SummaryStats::_H_found)
    .def_readonly("_H_HET_found", &SummaryStats::_H_HET_found)
    .def_readonly("_H_removed", &SummaryStats::_H_removed)
    .def_readonly("_H_HET_removed", &SummaryStats::_H_HET_removed)
    .def_readonly("_H_added", &SummaryStats::_H_added)
    .def_readonly("_H_HET_added", &SummaryStats::_H_HET_added)
    .def_readonly("_H_standardized", &SummaryStats::_H_standardized)
    .def_readonly("_H_HET_standardized", &SummaryStats::_H_HET_standardized)
    .def_readonly("_num_atoms", &SummaryStats::_num_atoms)
    .def_readonly("_conect", &SummaryStats::_conect)
    .def_readonly("_num_adj", &SummaryStats::_num_adj)
    .def_readonly("_num_renamed", &SummaryStats::_num_renamed)
  ;
  
/*
  class_<PDBrec, std::shared_ptr<PDBrec> >("PDBrec", init<>())
    .def(init<const PDB &>())
    .def(init<const PDBrec &>())
    .def("clone", &PDBrec::clone)
    .def(self == self)
    .add_property("x", getx, setx)
    .add_property("y", gety, sety)
    .add_property("z", getz, setz)
    .add_property("loc", getloc, setloc)
    .def("type", &PDBrec::type)
  ;
*/

  class_<ModelList>("ModelList")
    .def("size", &ModelList::size)
  ;
  class_<ModelsVector>("ModelsVector")
    .def(vector_indexing_suite<ModelsVector>() )
  ;
  /// @todo We need to get this working in a way that can be copied for Python
  //class_<CTab>("CTab", init<const std::string &, int>());
  class_<DotSphManager>("DotSphManager", init<float>());

  /// @todo We need to get this working in a way that can be copied for Python
/*
  class_<AtomPositions>("AtomPositions", init<int, bool, bool, bool, bool, int,
      float, float,
      float,
      DotSphManager &, float,
      float, float,
      bool, bool,
      bool, std::ostream &>());
*/

/*
  class_< std::vector< int > >("list_list_ptr_PDBrec")
    .def(vector_indexing_suite< std::vector< int >() )
  ;
*/


  // Export the global functions
  def("inputModels", inputModels, "Read the PDB records from the specified input stream.");
  def("dropHydrogens", dropHydrogens, "Drop hydrogens from a model in place.");
  def("optimize", optimize, "Optimize all of the records passed in in place.");
  def("checkSEGIDs", checkSEGIDs, "Check the list of PDB records to see if we should use segment ID as chain.");
  def("scanAndGroupRecords", scanAndGroupRecords);

  // Export the const char * variables, which will be read only.
  scope().attr("versionString") = versionString;
  scope().attr("shortVersion") = shortVersion;
  scope().attr("referenceString") = referenceString;
  scope().attr("electronicReference") = electronicReference;

  // Macro exports a getter and setter for a named global variable of the specified type
  #define exportGetSet(name) def("get" #name , get ## name, "Get value of " #name " global variable");\
                             def("set" #name , set ## name, "Set value of " #name " global variable");

  // Export getters and setters for the global variables that control behavior
  exportGetSet(Tally)

  exportGetSet(Verbose)
  exportGetSet(KeepConnections)
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
}

