// Enable functions with up to 20 parameters to be called.  Default of 15 is insufficient
#define BOOST_PYTHON_MAX_ARITY 20
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

// Describe and name compound classes that we need access to.
typedef std::list< std::shared_ptr<PDBrec> > ModelList;
typedef std::vector< ModelList > ModelsVector;
typedef std::list< std::shared_ptr<PDBrec> >::iterator ModelIterator;
typedef std::vector<std::string> StringVector;

// Desribe overloaded functions
BOOST_PYTHON_FUNCTION_OVERLOADS(outputRecords_all_overloads, outputRecords_all, 1, 2);

BOOST_PYTHON_MODULE(mmtbx_reduce_ext)
{
  // Export the class objects that Python will need access to,
  // along with their shared-pointer and vector types.
  class_<SummaryStats>("SummaryStats", init<>())
    .def_readwrite("_H_found", &SummaryStats::_H_found)
    .def_readwrite("_H_HET_found", &SummaryStats::_H_HET_found)
    .def_readwrite("_H_removed", &SummaryStats::_H_removed)
    .def_readwrite("_H_HET_removed", &SummaryStats::_H_HET_removed)
    .def_readwrite("_H_added", &SummaryStats::_H_added)
    .def_readwrite("_H_HET_added", &SummaryStats::_H_HET_added)
    .def_readwrite("_H_standardized", &SummaryStats::_H_standardized)
    .def_readwrite("_H_HET_standardized", &SummaryStats::_H_HET_standardized)
    .def_readwrite("_num_atoms", &SummaryStats::_num_atoms)
    .def_readwrite("_conect", &SummaryStats::_conect)
    .def_readwrite("_num_adj", &SummaryStats::_num_adj)
    .def_readwrite("_num_renamed", &SummaryStats::_num_renamed)
  ;
  
  // Select from among overloaded methods
  Coord (PDBrec::*getx)(void) const = &PDBrec::x;
  void (PDBrec::*setx)(Coord) = &PDBrec::x;
  Coord (PDBrec::*gety)(void) const = &PDBrec::y;
  void (PDBrec::*sety)(Coord) = &PDBrec::y;
  Coord (PDBrec::*getz)(void) const = &PDBrec::z;
  void (PDBrec::*setz)(Coord) = &PDBrec::z;
  //Point3D (PDBrec::*getloc)(void) const = &PDBrec::loc;
  //void (PDBrec::*setloc)(Point3D) = &PDBrec::loc;
  class_<PDBrec, std::shared_ptr<PDBrec> >("PDBrec", init<>())
    .def(init<const PDB &>())
    .def(init<const PDBrec &>())
    .def("clone", &PDBrec::clone)
    .def(self == self)
    .add_property("x", getx, setx)
    .add_property("y", gety, sety)
    .add_property("z", getz, setz)
    //.add_property("loc", getloc, setloc)
    .def("type", &PDBrec::type)
  ;

  // Select from among overloaded methods
  ModelIterator (ModelList::*mlbegin)(void) = &ModelList::begin;
  ModelIterator (ModelList::*mlend)(void) = &ModelList::end;
  class_<ModelList>("ModelList")
    .def("size", &ModelList::size)
    .def("begin", mlbegin)
    .def("end", mlend)
    .def("__iter__", iterator<ModelList>())
  ;
  class_<ModelIterator>("ModelIterator");

  class_<ModelsVector>("ModelsVector")
    .def(vector_indexing_suite<ModelsVector>() )
  ;
  
  class_<CTab>("CTab", init<const std::string &>());

  class_<DotSphManager>("DotSphManager", init<float>());

  class_<AtomPositions>("AtomPositions", init<int, bool, bool, bool, bool, int,
      float, float,
      float,
      DotSphManager &, float,
      float, float,
      bool, bool,
      bool>())
    .def("numChanges", &AtomPositions::numChanges)
    .def("describeChanges", &AtomPositions::describeChanges)
  ;

  class_<StringVector>("StringVector")
    .def(vector_indexing_suite<StringVector>() )
  ;

  // Export the global functions
  def("inputModels", inputModels, "Read the PDB records from the specified input stream.");
  def("dropHydrogens", dropHydrogens, "Drop hydrogens from a model in place.");
  def("optimize", optimize, "Optimize all of the records passed in in place.");
  def("checkSEGIDs", checkSEGIDs, "Check the list of PDB records to see if we should use segment ID as chain.");
  def("scanAndGroupRecords", scanAndGroupRecords);
  def("reduceList", reduceList);
  def("optimize", optimize);
  def("outputRecords_all", outputRecords_all, outputRecords_all_overloads());
  def("outputRecords_all_string", outputRecords_all_string);

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

