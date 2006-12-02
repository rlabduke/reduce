// name: StdResH.h
// author: J. Michael Word
// date written: 2/8/96 with many later modifications
// purpose: describes the connectivies, angles
//          and distances required to place
//          hydrogens on each standard residue

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#if defined(_MSC_VER)
#pragma warning(disable:4786) 
#pragma warning(disable:4800) 
#endif

#ifndef STDRESH_H
#define STDRESH_H 1

#include "AtomConn.h"
#include <list>
#include <map>
#include "ElementInfo.h"

// types of hydrogens-
// 1: HXR3 - requires just 4 atom centers
// 2: H2XR2- three atoms and an angle
// 3: H3XR - three atoms, an angle and dihedral
// 4: HXR2 - three atoms and a fudge factor (planar interpolating vectors)
// 5: HXR2 - three atoms and a fraction (planar based on angles)
// 6: HXY  - (linear) just two atoms
// 7: HX   - (can't create, can only adjust length) only one atom

class HydrogenPlanTable;
class StdResXtraInfo;

class StdResH {
public:
   StdResH(const char* rn, const char* ex) : _name(rn), _exclude(ex) {}
   ~StdResH() {}

   void addPlan(int type, const char* elem, const char *hname,
                const char *c1, const char *c2, const char *c3, const char *c4,
                double dist, double ang1, double ang2, int flags);

   const std::string& name() const { return _name; }
   const std::string& exclude() const { return _exclude; }

   // Be careful, std::list must be copied.
   std::list<atomPlacementPlan*> plans() { return _plans; }

   static const HydrogenPlanTable& HydPlanTbl();
   static const StdResXtraInfo&    ResXtraInfo();

private:
   static const HydrogenPlanTable *TheHydPlanTbl;    // class shared resource
   static const StdResXtraInfo    *TheResXtraInfo;   // class shared resource

   std::string                   _name;    // residue name
   std::string                   _exclude; // list of incompatible residues (for n-term, etc.)
   std::list<atomPlacementPlan*> _plans;
};

// -----------------------------------------------------------------------------
// implementation class for a table of Hydrogen building records by Residue type
// -----------------------------------------------------------------------------
class HydrogenPlanTable {
public:
   HydrogenPlanTable();
   ~HydrogenPlanTable() {
	   for (std::map<std::string, StdResH*>::const_iterator i = _restbl.begin(); i != _restbl.end(); ++i)
		   delete i->second;
   }

   StdResH* get(const std::string& resName) const {
	   std::map<std::string, StdResH*>::const_iterator i = _restbl.find(resName);
	   if (i != _restbl.end())
		   return i->second;
	   else
		   return NULL;
//	   return _restbl.get(resName); 
   }

private:
   HydrogenPlanTable(const HydrogenPlanTable&);           //can't copy
   HydrogenPlanTable& operator=(const HydrogenPlanTable&);//can't assign

   std::map<std::string, StdResH*> _restbl;

   struct addPlan_args {
     int type; const char* elem; const char *hname;
     const char *c1; const char *c2; const char *c3; const char *c4;
     double dist; double ang1; double ang2; int flags;
   };

   void
   insertStdResH(
     const char* rn,
     const char* ex,
     const addPlan_args* a);
};

// -----------------------------------------------------------------------------
// supplemental information about the standard residues (currently limited)
// -----------------------------------------------------------------------------
class StdResXtraInfo {
public:
   StdResXtraInfo();
   ~StdResXtraInfo() {}

   bool atomHasAttrib(const std::string& resName, const std::string& atomName, int attr) const {
	   std::map<std::string, int>::const_iterator i = _atomAttributes.find(resName + ":" + atomName);
	   if (i != _atomAttributes.end())
		   return i->second & attr;
	   else 
		   return FALSE;
//      int *p = _atomAttributes.get(makeKey(resName, atomName));
//      return (p == NULL) ? FALSE : ((*p) & attr);
   }

   bool match(const std::string& s1, const std::string& s2) const {
	   std::map<std::string, int>::const_iterator i = _atomAttributes.find(s1 + ":" + s2);
	   if (i != _atomAttributes.end())
		   return i->second != 0;
	   else 
		   return FALSE;
//      int *p = _atomAttributes.get(makeKey(s1, s2));
//      return (p == NULL) ? FALSE : ((*p) != 0);
   }

private:
   StdResXtraInfo(const StdResXtraInfo&);           //can't copy
   StdResXtraInfo& operator=(const StdResXtraInfo&);//can't assign

   std::string makeKey(const std::string& resName, const std::string& atomName) const {
      return resName + ":" + atomName;
   }

   std::map<std::string, int> _atomAttributes;
};
#endif
