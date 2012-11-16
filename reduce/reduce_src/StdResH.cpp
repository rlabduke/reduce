// name: StdResH.C
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

// types of hydrogens-
// 1: HXR3 - requires just 4 atom centers
// 2: H2XR2- three atoms and an angle
// 3: H3XR - three atoms, an angle and dihedral
// 4: HXR2 - three atoms and a fudge factor (planar interpolating vectors)
// 5: HXR2 - three atoms and a fraction (planar based on angles)
// 6: HXY  - (linear) just two atoms
// 7: HX   - (can't create, can only adjust length) only one atom

// "-C", "-O", "+N" refer to the prev and next residue

// the following expected conversions are done in other modules
// resname => key conversion ":HOH:DOD:H2O:D2O:WAT:TIP:SOL:MTO:" => "HOH"
// atoname => key conversion "?D??" => "?H??"

#include "StdResH.h"

const HydrogenPlanTable * StdResH::TheHydPlanTbl  = NULL;
const StdResXtraInfo    * StdResH::TheResXtraInfo = NULL;

// build TheHydPlanTbl and TheResXtraInfo only once, when first used

const HydrogenPlanTable& StdResH::HydPlanTbl()  {
   if (TheHydPlanTbl == NULL) { TheHydPlanTbl = new HydrogenPlanTable; }
   return *TheHydPlanTbl;
}
const StdResXtraInfo& StdResH::ResXtraInfo() {
   if (TheResXtraInfo == NULL) { TheResXtraInfo = new StdResXtraInfo; }
   return *TheResXtraInfo;
}

void StdResH::addPlan(int type, const char* elem, const char *hname,
					  const char *c1, const char *c2, const char *c3, const char *c4,
					  double dist, double ang1, double ang2, int flags) {
	
	ElementInfo *e = ElementInfo::StdElemTbl().element(elem);

	AtomConn names(hname, 0);
	names.addConn(c1);
	if (type <= 6) { names.addConn(c2); }
	if (type <= 5) { names.addConn(c3); }
	if (type == 1) { names.addConn(c4); }

	atomPlacementPlan* p = new atomPlacementPlan(type, *e, names, dist, ang1, ang2, flags);

	_plans.push_front(p);
//	_plans.insert(p);
}

void
HydrogenPlanTable::insertStdResH(
  const char* rn,
  const char* ex,
  const addPlan_args* a)
{
  StdResH* theRes = new StdResH(rn, ex);
  _restbl.insert(std::make_pair(theRes->name(), theRes));
  for(;a->elem;a++) {
    theRes->addPlan(
      a->type, a->elem, a->hname,
      a->c1, a->c2, a->c3, a->c4,
      a->dist, a->ang1, a->ang2, a->flags);
  }
}

HydrogenPlanTable::HydrogenPlanTable() {
//--------------------------------------------------------------------------
// the order of plans within each residue is the reverse of the order we want
// them to be generated because they will go into a sequence.
//--------------------------------------------------------------------------
  {
    static const addPlan_args args[] = {
      {5, "Hpol", " H", " N", " CA", "- C", "", 1.0,   0.0, 0.48, BONDBUMPFLAG},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("amide", "PRO", args); // *NON* N-terminal mc amide
  }
  {
    static const addPlan_args args[] = { //break-amide added 070925 by RMI
      {3, "Hpol", " H",   " N", " CA", " C", "", 1.0, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND|NH3FLAG},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("break-amide", "PRO", args); // N-terminal mc amide
  }
  {
    static const addPlan_args args[] = { //nt-amide updated 070703 by JJH
      {3, "Hpol", " HT3", " N", " CA", " C", "", 1.0, 109.5,  60.0, XPLORNAME|ROTATEONDEMAND|NH3FLAG},
      {3, "Hpol", " HT2", " N", " CA", " C", "", 1.0, 109.5, -60.0, XPLORNAME|ROTATEONDEMAND|NH3FLAG},
      {3, "Hpol", " HT1", " N", " CA", " C", "", 1.0, 109.5, 180.0, XPLORNAME|ROTATEONDEMAND|NH3FLAG},
      {3, "Hpol", " H3",   " N", " CA", " C", "", 1.0, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND|NH3FLAG},
      {3, "Hpol", " H2",   " N", " CA", " C", "", 1.0, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND|NH3FLAG},
      {3, "Hpol", " H1",   " N", " CA", " C", "", 1.0, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND|NH3FLAG},      
      {3, "Hpol", "3H",   " N", " CA", " C", "", 1.0, 109.5,  60.0, USEOLDNAMES|ROTATEONDEMAND|NH3FLAG},
      {3, "Hpol", "2H",   " N", " CA", " C", "", 1.0, 109.5, -60.0, USEOLDNAMES|ROTATEONDEMAND|NH3FLAG},
      {3, "Hpol", "1H",   " N", " CA", " C", "", 1.0, 109.5, 180.0, USEOLDNAMES|ROTATEONDEMAND|NH3FLAG},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("nt-amide", "PRO", args); // N-terminal mc amide
  }
  {
    static const addPlan_args args[] = { //nt-pro updated 070703 by JJH
      {2, "Hpol", " HT2", " N", " CD", " CA", "", 1.0, 126.5,   0.0, XPLORNAME|BONDBUMPFLAG},
      {2, "Hpol", " HT1", " N", " CD", " CA", "", 1.0,-126.5,   0.0, XPLORNAME|BONDBUMPFLAG},
      {2, "Hpol", " H3",   " N", " CD", " CA", "", 1.0, 126.5,   0.0, USENEWNAMES|BONDBUMPFLAG},
      {2, "Hpol", " H2",   " N", " CD", " CA", "", 1.0,-126.5,   0.0, USENEWNAMES|BONDBUMPFLAG},
      {2, "Hpol", "2H",   " N", " CD", " CA", "", 1.0, 126.5,   0.0, USEOLDNAMES|BONDBUMPFLAG},
      {2, "Hpol", "1H",   " N", " CD", " CA", "", 1.0,-126.5,   0.0, USEOLDNAMES|BONDBUMPFLAG},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("nt-pro", "", args); // N-terminal prolines
  }
  {
    static const addPlan_args args[] = {  // added plan for placing both hydrogens on a backbone only model RMI 070713
      {1, "H", " HA", " CA", " N", " C", " CB", 1.1,   0.0,   0.0, STRICTALTFLAG},
      {1, "H", " HA", " CA", " N", " C", " CB", 1.1,   0.0,   0.0, NOTBBMODEL|STRICTALTFLAG},
//      {2, "H", "2HA", " CA", " N", " C", "", 1.1, 126.5,   0.0, BACKBONEMODEL|USEOLDNAMES},
//      {2, "H", "1HA", " CA", " N", " C", "", 1.1,-126.5,   0.0, BACKBONEMODEL|USEOLDNAMES},
//      {2, "H", " HA2", " CA", " N", " C", "", 1.1, 126.5,   0.0, BACKBONEMODEL|USENEWNAMES},
//      {2, "H", " HA3", " CA", " N", " C", "", 1.1, 126.5,   0.0, BACKBONEMODEL|USENEWNAMES},
      {2, "H", " HA", " CA", " N", " C", "", 1.1,-126.5,   0.0, BACKBONEMODEL},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("alpha", "GLY", args); // mainchain alpha proton
  }
  {
    static const addPlan_args args[] = { //ribose phosphate backbone updated 070703 by JJH
      {3, "Hpol", "HO5'", " O5'", " C5'", " C4'", "",     1.0, 109.5, 180.0, USENEWNAMES|UNSUREDROPFLAG|ROTATEFLAG|IFNOPO4},
      {3, "Hpol", "HO3'", " O3'", " C3'", " C4'", "",     1.0, 109.5, 180.0, USENEWNAMES|UNSUREDROPFLAG|ROTATEFLAG|IFNOPO4},
      {3, "Hpol", "HO2'", " O2'", " C2'", " C3'", "",     1.0, 109.5, 180.0, USENEWNAMES|UNSUREDROPFLAG|ROTATEFLAG},
      {3, "Hpol", "HO5'", " O5*", " C5*", " C4*", "",     1.0, 109.5, 180.0, USENEWNAMES|UNSUREDROPFLAG|ROTATEFLAG|IFNOPO4},
      {3, "Hpol", "HO3'", " O3*", " C3*", " C4*", "",     1.0, 109.5, 180.0, USENEWNAMES|UNSUREDROPFLAG|ROTATEFLAG|IFNOPO4},
      {3, "Hpol", "HO2'", " O2*", " C2*", " C3*", "",     1.0, 109.5, 180.0, USENEWNAMES|UNSUREDROPFLAG|ROTATEFLAG},
      {1, "H",    " H2'", " C2'", " C3'", " C1'", " O2'", 1.1,   0.0,   0.0, USENEWNAMES|O2PRIMEFLAG},
      {2, "H",    "H2''", " C2'", " C3'", " C1'", "",     1.1, 126.5,   0.0, USENEWNAMES|NOO2PRIMEFLAG},
      {2, "H",    " H2'", " C2'", " C3'", " C1'", "",     1.1,-126.5,   0.0, USENEWNAMES|NOO2PRIMEFLAG},
      {1, "H",    " H2'", " C2*", " C3*", " C1*", " O2*", 1.1,   0.0,   0.0, USENEWNAMES|O2PRIMEFLAG},
      {2, "H",    "H2''", " C2*", " C3*", " C1*", "",     1.1, 126.5,   0.0, USENEWNAMES|NOO2PRIMEFLAG},
      {2, "H",    " H2'", " C2*", " C3*", " C1*", "",     1.1,-126.5,   0.0, USENEWNAMES|NOO2PRIMEFLAG},
      {1, "H",    " H3'", " C3'", " C4'", " C2'", " O3'", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H4'", " C4'", " C5'", " C3'", " O4'", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H3'", " C3*", " C4*", " C2*", " O3*", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H4'", " C4*", " C5*", " C3*", " O4*", 1.1,   0.0,   0.0, USENEWNAMES},
      {2, "H",    "H5''", " C5'", " C4'", " O5'", "",     1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H",    " H5'", " C5'", " C4'", " O5'", "",     1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H",    "H5''", " C5*", " C4*", " O5*", "",     1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H",    " H5'", " C5*", " C4*", " O5*", "",     1.1,-126.5,   0.0, USENEWNAMES},

      {3, "Hpol", "2HO*", " O2*", " C2*", " C3*", "",     1.0, 109.5, 180.0, USEOLDNAMES|XPLORNAME|UNSUREDROPFLAG|ROTATEFLAG},
      {1, "H",    " H2*", " C2*", " C3*", " C1*", " O2*", 1.1,   0.0,   0.0, USEOLDNAMES|XPLORNAME|O2PRIMEFLAG},
      {2, "H",    "2H2*", " C2*", " C3*", " C1*", "",     1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME|NOO2PRIMEFLAG},
      {2, "H",    "1H2*", " C2*", " C3*", " C1*", "",     1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME|NOO2PRIMEFLAG},
      {1, "H",    " H3*", " C3*", " C4*", " C2*", " O3*", 1.1,   0.0,   0.0, USEOLDNAMES|XPLORNAME},
      {1, "H",    " H4*", " C4*", " C5*", " C3*", " O4*", 1.1,   0.0,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H",    "2H5*", " C5*", " C4*", " O5*", "",     1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H",    "1H5*", " C5*", " C4*", " O5*", "",     1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {3, "Hpol", " H5T", " O5'", " C5'", " C4'", "",     1.0, 109.5, 180.0, USEOLDNAMES|XPLORNAME|UNSUREDROPFLAG|ROTATEFLAG|IFNOPO4},
      {3, "Hpol", " H3T", " O3'", " C3'", " C4'", "",     1.0, 109.5, 180.0, USEOLDNAMES|XPLORNAME|UNSUREDROPFLAG|ROTATEFLAG|IFNOPO4},

      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("ribose phosphate backbone", "", args); // phosphate chain (two different naming schemes)
  }
//--------------------------------------------------------------------------
  {
    static const addPlan_args args[] = { //LYS updated 070702 by JJH
      {3, "Hpol", " HZ3", " NZ", " CE", " CD", "", 1.0, 109.5,  60.0, XPLORNAME|ROTATEONDEMAND|NH3FLAG},
      {3, "Hpol", " HZ2", " NZ", " CE", " CD", "", 1.0, 109.5, -60.0, XPLORNAME|ROTATEONDEMAND|NH3FLAG},
      {3, "Hpol", " HZ1", " NZ", " CE", " CD", "", 1.0, 109.5, 180.0, XPLORNAME|ROTATEONDEMAND|NH3FLAG},
      {3, "Hpol", " HZ3",  " NZ", " CE", " CD", "", 1.0, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND|NH3FLAG},
      {3, "Hpol", " HZ2",  " NZ", " CE", " CD", "", 1.0, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND|NH3FLAG},
      {3, "Hpol", " HZ1",  " NZ", " CE", " CD", "", 1.0, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND|NH3FLAG},
      {3, "Hpol", "3HZ",  " NZ", " CE", " CD", "", 1.0, 109.5,  60.0, USEOLDNAMES|ROTATEONDEMAND|NH3FLAG},
      {3, "Hpol", "2HZ",  " NZ", " CE", " CD", "", 1.0, 109.5, -60.0, USEOLDNAMES|ROTATEONDEMAND|NH3FLAG},
      {3, "Hpol", "1HZ",  " NZ", " CE", " CD", "", 1.0, 109.5, 180.0, USEOLDNAMES|ROTATEONDEMAND|NH3FLAG},
      {2, "H", " HE3",  " CE", " CD", " NZ", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HE2",  " CE", " CD", " NZ", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", " HD3",  " CD", " CG", " CE", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HD2",  " CD", " CG", " CE", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", " HG3",  " CG", " CB", " CD", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HG2",  " CG", " CB", " CD", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", " HB3",  " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2",  " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},      
      {2, "H", "2HE",  " CE", " CD", " NZ", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HE",  " CE", " CD", " NZ", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "2HD",  " CD", " CG", " CE", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HD",  " CD", " CG", " CE", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "2HG",  " CG", " CB", " CD", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HG",  " CG", " CB", " CD", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "2HB",  " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB",  " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("LYS", "", args);
  }
  {
    static const addPlan_args args[] = { //GLY updated 070702 by JJH
      {2, "H", " HA3", " CA", " N", " C", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HA2", " CA", " N", " C", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", "2HA", " CA", " N", " C", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HA", " CA", " N", " C", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("GLY", "", args);
  }
  {
    static const addPlan_args args[] = { //GLU updated 070702 by JJH
      {2, "H", " HG3", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HG2", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", " HB3", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", "2HG", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HG", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("GLU", "", args);
  }
  {
    static const addPlan_args args[] = { //THR updated 070702 by JJH
      {3, "H", "HG23", " CG2", " CB", " CA", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HG22", " CG2", " CB", " CA", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HG21", " CG2", " CB", " CA", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "3HG2", " CG2", " CB", " CA", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2HG2", " CG2", " CB", " CA", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1HG2", " CG2", " CB", " CA", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "Hpol", " HG1", " OG1", " CB", " CA", "", 1.0, 109.5, 180.0, UNSUREDROPFLAG|ROTATEFLAG},
      {1, "H", " HB", " CB", " CA", " OG1", " CG2", 1.1,   0.0,   0.0, 0},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("THR", "", args);
  }
  {
    static const addPlan_args args[] = { //ALA updated 070702 by JJH
      {3, "H", " HB3", " CB", " CA", " N", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", " HB2", " CB", " CA", " N", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", " HB1", " CB", " CA", " N", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "3HB", " CB", " CA", " N", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2HB", " CB", " CA", " N", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1HB", " CB", " CA", " N", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("ALA", "", args);
  }
  {
    static const addPlan_args args[] = { //PHE updated 070702 by JJH
      {4, "Har", " HZ",  " CZ",  " CE1", " CE2", "", 1.1,   0.0,   0.0, 0},
      {4, "Har", " HE2", " CE2", " CZ",  " CD2", "", 1.1,   0.0,   0.0, 0},
      {4, "Har", " HE1", " CE1", " CD1", " CZ",  "", 1.1,   0.0,   0.0, 0},
      {4, "Har", " HD2", " CD2", " CE2", " CG",  "", 1.1,   0.0,   0.0, 0},
      {4, "Har", " HD1", " CD1", " CG",  " CE1", "", 1.1,   0.0,   0.0, 0},
      {2, "H", " HB3",  " CB",  " CA",  " CG",  "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2",  " CB",  " CA",  " CG",  "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", "2HB",  " CB",  " CA",  " CG",  "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB",  " CB",  " CA",  " CG",  "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("PHE", "", args);
  }
  {
    static const addPlan_args args[] = { //ARG updated 070702 by JJH
      {3, "Hpol", "HH22", " NH2", " CZ", " NE", "", 1.0, 120.0, 180.0, XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "HH21", " NH2", " CZ", " NE", "", 1.0, 120.0,   0.0, XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "HH12", " NH1", " CZ", " NE", "", 1.0, 120.0, 180.0, XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "HH11", " NH1", " CZ", " NE", "", 1.0, 120.0,   0.0, XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "HH22", " NH2", " CZ", " NE", "", 1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", "HH21", " NH2", " CZ", " NE", "", 1.0, 120.0,   0.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", "HH12", " NH1", " CZ", " NE", "", 1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", "HH11", " NH1", " CZ", " NE", "", 1.0, 120.0,   0.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", "2HH2", " NH2", " CZ", " NE", "", 1.0, 120.0, 180.0, USEOLDNAMES|BONDBUMPFLAG},
      {3, "Hpol", "1HH2", " NH2", " CZ", " NE", "", 1.0, 120.0,   0.0, USEOLDNAMES|BONDBUMPFLAG},
      {3, "Hpol", "2HH1", " NH1", " CZ", " NE", "", 1.0, 120.0, 180.0, USEOLDNAMES|BONDBUMPFLAG},
      {3, "Hpol", "1HH1", " NH1", " CZ", " NE", "", 1.0, 120.0,   0.0, USEOLDNAMES|BONDBUMPFLAG},
      {4, "Hpol", " HE", " NE", " CD", " CZ", "", 1.0,   0.0,   0.0, 0},
      {2, "H", " HD3", " CD", " CG", " NE", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HD2", " CD", " CG", " NE", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", " HG3", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HG2", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", " HB3", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", "2HD", " CD", " CG", " NE", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HD", " CD", " CG", " NE", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "2HG", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HG", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("ARG", "", args);
  }
  {
    static const addPlan_args args[] = { //HIS updated 070703 by JJH
      {4, "Ha+p", " HE2", " NE2", " CE1", " CD2", "", 1.0,   0.0,   0.0, XTRAFLAG|BONDBUMPFLAG},
      {4, "Har", " HE1", " CE1", " ND1", " NE2", "", 1.1,   0.0,   0.0, 0},
      {4, "Har", " HD2", " CD2", " NE2", " CG", "", 1.1,   0.0,   0.0, 0},
      {4, "Ha+p", " HD1", " ND1", " CG", " CE1", "", 1.0,   0.0,   0.0, XTRAFLAG|BONDBUMPFLAG},
      {2, "H", " HB3", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("HIS", "", args);
  }
  {
    static const addPlan_args args[] = { //MET updated 070703 by JJH
      {3, "H", " HE3", " CE", " SD", " CG", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEFLAG},
      {3, "H", " HE2", " CE", " SD", " CG", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEFLAG},
      {3, "H", " HE1", " CE", " SD", " CG", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEFLAG},
      {2, "H", " HG3", " CG", " CB", " SD", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HG2", " CG", " CB", " SD", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", " HB3", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {3, "H", "3HE", " CE", " SD", " CG", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEFLAG},
      {3, "H", "2HE", " CE", " SD", " CG", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEFLAG},
      {3, "H", "1HE", " CE", " SD", " CG", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEFLAG},
      {2, "H", "2HG", " CG", " CB", " SD", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HG", " CG", " CB", " SD", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("MET", "", args);
  }
  {
    static const addPlan_args args[] = { //ASP updated 070703 by JJH
      {2, "H", " HB3", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("ASP", "", args);
  }
  {
    static const addPlan_args args[] = { //SER updated 070703 by JJH reverse handedness of SER HB RMI 070731
      {3, "Hpol", " HG", " OG", " CB", " CA", "", 1.0, 109.5, 180.0, UNSUREDROPFLAG|ROTATEFLAG},
      {2, "H", " HB3", " CB", " CA", " OG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2", " CB", " CA", " OG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", "2HB", " CB", " CA", " OG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " OG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("SER", "", args);
  }
  {
    static const addPlan_args args[] = { //ASN updated 070703 by JJH
      {3, "Hpol", "HD22", " ND2", " CG", " OD1", "", 1.0, 120.0, 180.0, XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "HD21", " ND2", " CG", " OD1", "", 1.0, 120.0,   0.0, XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "HD22", " ND2", " CG", " OD1", "", 1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", "HD21", " ND2", " CG", " OD1", "", 1.0, 120.0,   0.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", "2HD2", " ND2", " CG", " OD1", "", 1.0, 120.0, 180.0, USEOLDNAMES|BONDBUMPFLAG},
      {3, "Hpol", "1HD2", " ND2", " CG", " OD1", "", 1.0, 120.0,   0.0, USEOLDNAMES|BONDBUMPFLAG},
      {2, "H", " HB3", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("ASN", "", args);
  }
  {
    static const addPlan_args args[] = { //TYR updated 070703 by JJH
      {3, "Hpol", " HH", " OH", " CZ", " CE1", "", 1.0, 109.5, 180.0, UNSUREDROPFLAG|ROTATEFLAG},
      {4, "Har", " HE2", " CE2", " CZ", " CD2", "", 1.1,   0.0,   0.0, 0},
      {4, "Har", " HE1", " CE1", " CD1", " CZ", "", 1.1,   0.0,   0.0, 0},
      {4, "Har", " HD2", " CD2", " CE2", " CG", "", 1.1,   0.0,   0.0, 0},
      {4, "Har", " HD1", " CD1", " CG", " CE1", "", 1.1,   0.0,   0.0, 0},
      {2, "H", " HB3", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("TYR", "", args);
  }
  {
    static const addPlan_args args[] = { //CYS updated 070703 by JJH
      {3, "Hpol", " HG", " SG", " CB", " CA", "", 1.3, 109.5,  180.0,UNSUREDROPFLAG|ROTATEFLAG},
      {2, "H", " HB3", " CB", " CA", " SG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2", " CB", " CA", " SG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", "2HB", " CB", " CA", " SG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " SG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("CYS", "", args);
  }
  {
    static const addPlan_args args[] = { //GLN updated 070703 by JJH
      {3, "Hpol", "HE22", " NE2", " CD", " OE1", "", 1.0, 120.0, 180.0, XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "HE21", " NE2", " CD", " OE1", "", 1.0, 120.0,   0.0, XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "HE22", " NE2", " CD", " OE1", "", 1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", "HE21", " NE2", " CD", " OE1", "", 1.0, 120.0,   0.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", "2HE2", " NE2", " CD", " OE1", "", 1.0, 120.0, 180.0, USEOLDNAMES|BONDBUMPFLAG},
      {3, "Hpol", "1HE2", " NE2", " CD", " OE1", "", 1.0, 120.0,   0.0, USEOLDNAMES|BONDBUMPFLAG},
      {2, "H", " HG3", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HG2", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", " HB3", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", "2HG", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HG", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("GLN", "", args);
  }
  {
    static const addPlan_args args[] = { //LEU updated 070703 by JJH
      {3, "H", "HD23", " CD2", " CG", " CB", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HD22", " CD2", " CG", " CB", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HD21", " CD2", " CG", " CB", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HD13", " CD1", " CG", " CB", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HD12", " CD1", " CG", " CB", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HD11", " CD1", " CG", " CB", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "3HD2", " CD2", " CG", " CB", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2HD2", " CD2", " CG", " CB", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1HD2", " CD2", " CG", " CB", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "3HD1", " CD1", " CG", " CB", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2HD1", " CD1", " CG", " CB", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1HD1", " CD1", " CG", " CB", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {1, "H", " HG", " CG", " CB", " CD1", " CD2", 1.1,   0.0,   0.0, 0},
      {2, "H", " HB3", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("LEU", "", args);
  }
  {
    static const addPlan_args args[] = { //PRO updated 070703 by JJH
      {2, "H", " HD3", " CD", " CG", " N", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HD2", " CD", " CG", " N", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", " HG3", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HG2", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", " HB3", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", "2HD", " CD", " CG", " N", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HD", " CD", " CG", " N", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "2HG", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HG", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("PRO", "", args);
  }
  {
    static const addPlan_args args[] = { //VAL updated 070703 by JJH
      {3, "H", "HG23", " CG2", " CB", " CA", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HG22", " CG2", " CB", " CA", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HG21", " CG2", " CB", " CA", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HG13", " CG1", " CB", " CA", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HG12", " CG1", " CB", " CA", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HG11", " CG1", " CB", " CA", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "3HG2", " CG2", " CB", " CA", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2HG2", " CG2", " CB", " CA", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1HG2", " CG2", " CB", " CA", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "3HG1", " CG1", " CB", " CA", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2HG1", " CG1", " CB", " CA", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1HG1", " CG1", " CB", " CA", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {1, "H", " HB", " CB", " CA", " CG1", " CG2", 1.1,   0.0,   0.0, 0},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("VAL", "", args);
  }
  {
    static const addPlan_args args[] = {  //ILE updated 070703 by JJH
      {3, "H", "HD13", " CD1", " CG1", " CB", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HD12", " CD1", " CG1", " CB", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HD11", " CD1", " CG1", " CB", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HG23", " CG2", " CB", " CA", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HG22", " CG2", " CB", " CA", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HG21", " CG2", " CB", " CA", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {2, "H", "HG13", " CG1", " CB", " CD1", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", "HG12", " CG1", " CB", " CD1", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {3, "H", "3HD1", " CD1", " CG1", " CB", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2HD1", " CD1", " CG1", " CB", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1HD1", " CD1", " CG1", " CB", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "3HG2", " CG2", " CB", " CA", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2HG2", " CG2", " CB", " CA", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1HG2", " CG2", " CB", " CA", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {2, "H", "2HG1", " CG1", " CB", " CD1", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HG1", " CG1", " CB", " CD1", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {1, "H", " HB", " CB", " CA", " CG1", " CG2", 1.1,   0.0,   0.0, 0},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("ILE", "", args);
  }
  {
    static const addPlan_args args[] = { //TRP updated 070703 by JJH
      {4, "Har", " HH2", " CH2", " CZ3", " CZ2", "", 1.1,   0.0,   0.0, 0},
      {4, "Har", " HZ3", " CZ3", " CE3", " CH2", "", 1.1,   0.0,   0.0, 0},
      {4, "Har", " HZ2", " CZ2", " CH2", " CE2", "", 1.1,   0.0,   0.0, 0},
      {4, "Har", " HE3", " CE3", " CD2", " CZ3", "", 1.1,   0.0,   0.0, 0},
      {4, "Ha+p", " HE1", " NE1", " CE2", " CD1", "", 1.0,   0.0,   0.0, 0},
      {4, "Har", " HD1", " CD1", " NE1", " CG", "", 1.1,   0.0,   0.0, 0},
      {2, "H", " HB3", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("TRP", "", args);
  }
//--------------------------------------------------------------------------
  {
    static const addPlan_args args[] = {
      {4, "Har", " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0},
      {4, "Har", " H5",  " C5",  " C4",  " C6",  "",    1.1,   0.0,   0.0, 0},
      {4, "Ha+p"," H3",  " N3",  " C4",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {1, "H",   " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USEOLDNAMES},
      {1, "H",   " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, XPLORNAME},
      {1, "H",   " H1'",  " C1*",  " O4*",  " C2*", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",   " H1'",  " C1'",  " O4'",  " C2'", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("  U", "", args);
  }
#ifdef LEFT_JUSTIFY_NUC_RES_OK
  {
    static const addPlan_args args[] = {
      {4, "Har", " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0},
      {4, "Har", " H5",  " C5",  " C4",  " C6",  "",    1.1,   0.0,   0.0, 0},
      {4, "Ha+p"," H3",  " N3",  " C4",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {1, "H",   " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USEOLDNAMES},
      {1, "H",   " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, XPLORNAME},
      {1, "H",   " H1'",  " C1*",  " O4*",  " C2*", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",   " H1'",  " C1'",  " O4'",  " C2'", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("U", "", args);
  }
#endif
  {
    static const addPlan_args args[] = {
      {4, "Har", " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0},
      {4, "Har", " H5",  " C5",  " C4",  " C6",  "",    1.1,   0.0,   0.0, 0},
      {4, "Ha+p"," H3",  " N3",  " C4",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {1, "H",   " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USEOLDNAMES},
      {1, "H",   " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, XPLORNAME},
      {1, "H",   " H1'",  " C1*",  " O4*",  " C2*", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",   " H1'",  " C1'",  " O4'",  " C2'", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("URA", "", args);
  }
//--------------------------------------------------------------------------
  // note C5M is an alternative name for C5A  --- C7 is the remediated name
  {
    static const addPlan_args args[] = {
      {4, "Har",  " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0},
      {3, "H",    "3H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "2H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "1H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {7, "H",    " H53", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {7, "H",    " H52", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {7, "H",    " H51", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {3, "H",    "3H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "2H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "1H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "3H5A", " C7", " C5",  " C4",  "",    1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "2H5A", " C7", " C5",  " C4",  "",    1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "1H5A", " C7", " C5",  " C4",  "",    1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {7, "H",    " H53", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {7, "H",    " H52", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {7, "H",    " H51", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {3, "H",    " H73", " C5M", " C5",  " C4",  "",    1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H72", " C5M", " C5",  " C4",  "",    1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H71", " C5M", " C5",  " C4",  "",    1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H73", " C5A", " C5",  " C4",  "",    1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H72", " C5A", " C5",  " C4",  "",    1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H71", " C5A", " C5",  " C4",  "",    1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H73", " C7", " C5",  " C4",  "",    1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H72", " C7", " C5",  " C4",  "",    1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H71", " C7", " C5",  " C4",  "",    1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {4, "Ha+p", " H3",  " N3",  " C4",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {1, "H",    " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USEOLDNAMES|XPLORNAME},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, XPLORNAME},
      {1, "H",    " H1'", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("  T", "", args);
  }
#ifdef LEFT_JUSTIFY_NUC_RES_OK
  {
    static const addPlan_args args[] = {
      {4, "Har",  " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0},
      {3, "H",    "3H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "2H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "1H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {7, "H",    " H53", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {7, "H",    " H52", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {7, "H",    " H51", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {3, "H",    "3H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "2H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "1H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "3H5A", " C7", " C5",  " C4",  "",    1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "2H5A", " C7", " C5",  " C4",  "",    1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "1H5A", " C7", " C5",  " C4",  "",    1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {7, "H",    " H53", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {7, "H",    " H52", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {7, "H",    " H51", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {3, "H",    " H73", " C5M", " C5",  " C4",  "",    1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H72", " C5M", " C5",  " C4",  "",    1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H71", " C5M", " C5",  " C4",  "",    1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H73", " C5A", " C5",  " C4",  "",    1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H72", " C5A", " C5",  " C4",  "",    1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H71", " C5A", " C5",  " C4",  "",    1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H73", " C7", " C5",  " C4",  "",    1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H72", " C7", " C5",  " C4",  "",    1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H71", " C7", " C5",  " C4",  "",    1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {4, "Ha+p", " H3",  " N3",  " C4",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {1, "H",    " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USEOLDNAMES|XPLORNAME},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, XPLORNAME}, 
      {1, "H",    " H1'", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("T", "", args);
  }
#endif
  {
    static const addPlan_args args[] = {
      {4, "Har",  " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0},
      {3, "H",    "3H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "2H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "1H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {7, "H",    " H53", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {7, "H",    " H52", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {7, "H",    " H51", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {3, "H",    "3H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "2H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "1H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "3H5A", " C7", " C5",  " C4",  "",    1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "2H5A", " C7", " C5",  " C4",  "",    1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "1H5A", " C7", " C5",  " C4",  "",    1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {7, "H",    " H53", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {7, "H",    " H52", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {7, "H",    " H51", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {3, "H",    " H73", " C5M", " C5",  " C4",  "",    1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H72", " C5M", " C5",  " C4",  "",    1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H71", " C5M", " C5",  " C4",  "",    1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H73", " C5A", " C5",  " C4",  "",    1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H72", " C5A", " C5",  " C4",  "",    1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H71", " C5A", " C5",  " C4",  "",    1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H73", " C7", " C5",  " C4",  "",    1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H72", " C7", " C5",  " C4",  "",    1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H71", " C7", " C5",  " C4",  "",    1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {4, "Ha+p", " H3",  " N3",  " C4",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {1, "H",    " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USEOLDNAMES|XPLORNAME},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, XPLORNAME},      
      {1, "H",    " H1'", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("THY", "", args);
  }
//--------------------------------------------------------------------------
  {
    static const addPlan_args args[] = {
      {4, "Har",  " H2",  " C2",  " N1",  " N3",  "",    1.1,   0.0,   0.0, 0},
      {3, "Hpol", "2H6",  " N6",  " C6",  " C5",  "",    1.0,-120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "1H6",  " N6",  " C6",  " C5",  "",    1.0, 120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {7, "Hpol", " H62", " N6",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {7, "Hpol", " H61", " N6",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {3, "Hpol", " H62",  " N6",  " C6",  " C5",  "",    1.0,-120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", " H61",  " N6",  " C6",  " C5",  "",    1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {4, "Har",  " H8",  " C8",  " N7",  " N9",  "",    1.1,   0.0,   0.0, 0},
      {1, "H",    " H1*", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, USEOLDNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, XPLORNAME},
      {1, "H",    " H1'",  " C1*",  " O4*",  " C2*",  " N9", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H1'",  " C1'",  " O4'",  " C2'",  " N9", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("  A", "", args);
  }
#ifdef LEFT_JUSTIFY_NUC_RES_OK
  {
    static const addPlan_args args[] = {
      {4, "Har",  " H2",  " C2",  " N1",  " N3",  "",    1.1,   0.0,   0.0, 0},
      {3, "Hpol", "2H6",  " N6",  " C6",  " C5",  "",    1.0,-120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "1H6",  " N6",  " C6",  " C5",  "",    1.0, 120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {7, "Hpol", " H62", " N6",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {7, "Hpol", " H61", " N6",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {3, "Hpol", " H62",  " N6",  " C6",  " C5",  "",    1.0,-120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", " H61",  " N6",  " C6",  " C5",  "",    1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {4, "Har",  " H8",  " C8",  " N7",  " N9",  "",    1.1,   0.0,   0.0, 0},
      {1, "H",    " H1*", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, USEOLDNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, XPLORNAME},
      {1, "H",    " H1'",  " C1*",  " O4*",  " C2*",  " N9", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H1'",  " C1'",  " O4'",  " C2'",  " N9", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("A", "", args);
  }
#endif
  {
    static const addPlan_args args[] = {
      {3, "Hpol", " H62",  " N6",  " C6",  " C5",  "",    1.0,-120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", " H61",  " N6",  " C6",  " C5",  "",    1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {1, "H",    " H1'",  " C1*",  " O4*",  " C2*",  " N9", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H1'",  " C1'",  " O4'",  " C2'",  " N9", 1.1,   0.0,   0.0, USENEWNAMES},

      {4, "Har",  " H2",  " C2",  " N1",  " N3",  "",    1.1,   0.0,   0.0, 0},
      {3, "Hpol", "2H6",  " N6",  " C6",  " C5",  "",    1.0,-120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "1H6",  " N6",  " C6",  " C5",  "",    1.0, 120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {7, "Hpol", " H62", " N6",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {7, "Hpol", " H61", " N6",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {4, "Har",  " H8",  " C8",  " N7",  " N9",  "",    1.1,   0.0,   0.0, 0},
      {1, "H",    " H1*", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, USEOLDNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, XPLORNAME}, 
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("ADE", "", args);
  }
//--------------------------------------------------------------------------
  {
    static const addPlan_args args[] = {
      {4, "Har",  " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0},
      {4, "Har",  " H5",  " C5",  " C4",  " C6",  "",    1.1,   0.0,   0.0, 0},
      {3, "Hpol", "2H4",  " N4",  " C4",  " N3",  "",    1.0,-120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "1H4",  " N4",  " C4",  " N3",  "",    1.0, 120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {7, "Hpol", " H42", " N4",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {7, "Hpol", " H41", " N4",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {3, "Hpol", " H42",  " N4",  " C4",  " N3",  "",    1.0,-120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", " H41",  " N4",  " C4",  " N3",  "",    1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {1, "H",    " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USEOLDNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, XPLORNAME},
      {1, "H",    " H1'", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("  C", "", args);
  }
#ifdef LEFT_JUSTIFY_NUC_RES_OK
  {
    static const addPlan_args args[] = {
      {4, "Har",  " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0},
      {4, "Har",  " H5",  " C5",  " C4",  " C6",  "",    1.1,   0.0,   0.0, 0},
      {3, "Hpol", "2H4",  " N4",  " C4",  " N3",  "",    1.0,-120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "1H4",  " N4",  " C4",  " N3",  "",    1.0, 120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {7, "Hpol", " H42", " N4",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {7, "Hpol", " H41", " N4",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {3, "Hpol", " H42",  " N4",  " C4",  " N3",  "",    1.0,-120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", " H41",  " N4",  " C4",  " N3",  "",    1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {1, "H",    " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USEOLDNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, XPLORNAME},
      {1, "H",    " H1'", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0} 
    };
    insertStdResH("C", "", args);
  }
#endif
  {
    static const addPlan_args args[] = {
      {4, "Har",  " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0},
      {4, "Har",  " H5",  " C5",  " C4",  " C6",  "",    1.1,   0.0,   0.0, 0},
      {3, "Hpol", "2H4",  " N4",  " C4",  " N3",  "",    1.0,-120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "1H4",  " N4",  " C4",  " N3",  "",    1.0, 120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {7, "Hpol", " H42", " N4",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {7, "Hpol", " H41", " N4",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {3, "Hpol", " H42",  " N4",  " C4",  " N3",  "",    1.0,-120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", " H41",  " N4",  " C4",  " N3",  "",    1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {1, "H",    " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USEOLDNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, XPLORNAME},
      {1, "H",    " H1'", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("CYT", "", args);
  }
//--------------------------------------------------------------------------
  {
    static const addPlan_args args[] = {
      {3, "Hpol", "2H2",  " N2",  " C2",  " N1",  "",    1.0,-120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "1H2",  " N2",  " C2",  " N1",  "",    1.0, 120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {7, "Hpol", " H22", " N2",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {7, "Hpol", " H21", " N2",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {3, "Hpol", " H22",  " N2",  " C2",  " N1",  "",    1.0,-120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", " H21",  " N2",  " C2",  " N1",  "",    1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {4, "Ha+p", " H1",  " N1",  " C6",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {4, "Har",  " H8",  " C8",  " N7",  " N9",  "",    1.1,   0.0,   0.0, 0},
      {1, "H",    " H1*", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, USEOLDNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, XPLORNAME},
      {1, "H",    " H1'", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("  G", "", args);
  }
#ifdef LEFT_JUSTIFY_NUC_RES_OK
  {
    static const addPlan_args args[] = {
      {3, "Hpol", "2H2",  " N2",  " C2",  " N1",  "",    1.0,-120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "1H2",  " N2",  " C2",  " N1",  "",    1.0, 120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {7, "Hpol", " H22", " N2",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {7, "Hpol", " H21", " N2",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {3, "Hpol", " H22",  " N2",  " C2",  " N1",  "",    1.0,-120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", " H21",  " N2",  " C2",  " N1",  "",    1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {4, "Ha+p", " H1",  " N1",  " C6",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {4, "Har",  " H8",  " C8",  " N7",  " N9",  "",    1.1,   0.0,   0.0, 0},
      {1, "H",    " H1*", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, USEOLDNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, XPLORNAME},
      {1, "H",    " H1'", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("G", "", args);
  }
#endif
  {
    static const addPlan_args args[] = {
      {3, "Hpol", "2H2",  " N2",  " C2",  " N1",  "",    1.0,-120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "1H2",  " N2",  " C2",  " N1",  "",    1.0, 120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {7, "Hpol", " H22", " N2",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},

      {3, "Hpol", " H21",  " N2",  " C2",  " N1",  "",    1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {4, "Ha+p", " H1",  " N1",  " C6",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {4, "Har",  " H8",  " C8",  " N7",  " N9",  "",    1.1,   0.0,   0.0, 0},
      {1, "H",    " H1*", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, USEOLDNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, XPLORNAME},
      {1, "H",    " H1'", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("GUA", "", args);
  }
//--------------------------------------------------------------------------
  { 
    static const addPlan_args args[] = {
      {4, "Har",  " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0},
      {3, "H",    "3H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "2H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "1H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {7, "H",    " H53", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {7, "H",    " H52", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {7, "H",    " H51", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {3, "H",    "3H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "2H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "1H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "3H5A", " C7", " C5",  " C4",  "",    1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "2H5A", " C7", " C5",  " C4",  "",    1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H",    "1H5A", " C7", " C5",  " C4",  "",    1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {7, "H",    " H53", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {7, "H",    " H52", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {7, "H",    " H51", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND},
      {3, "H",    " H73", " C5M", " C5",  " C4",  "",    1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H72", " C5M", " C5",  " C4",  "",    1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H71", " C5M", " C5",  " C4",  "",    1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H73", " C5A", " C5",  " C4",  "",    1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H72", " C5A", " C5",  " C4",  "",    1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H71", " C5A", " C5",  " C4",  "",    1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H73", " C7", " C5",  " C4",  "",    1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H72", " C7", " C5",  " C4",  "",    1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H",    " H71", " C7", " C5",  " C4",  "",    1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {4, "Ha+p", " H3",  " N3",  " C4",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {1, "H",    " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USEOLDNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, XPLORNAME},
      {1, "H",    " H1'", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH(" DT", "", args);
  }
  {
    static const addPlan_args args[] = {
      {4, "Har",  " H2",  " C2",  " N1",  " N3",  "",    1.1,   0.0,   0.0, 0},
      {3, "Hpol", "2H6",  " N6",  " C6",  " C5",  "",    1.0,-120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "1H6",  " N6",  " C6",  " C5",  "",    1.0, 120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {7, "Hpol", " H62", " N6",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {7, "Hpol", " H61", " N6",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {3, "Hpol", " H62",  " N6",  " C6",  " C5",  "",    1.0,-120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", " H61",  " N6",  " C6",  " C5",  "",    1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {4, "Har",  " H8",  " C8",  " N7",  " N9",  "",    1.1,   0.0,   0.0, 0},
      {1, "H",    " H1*", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, USEOLDNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, XPLORNAME},
      {1, "H",    " H1'",  " C1*",  " O4*",  " C2*",  " N9", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H1'",  " C1'",  " O4'",  " C2'",  " N9", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH(" DA", "", args);
  }
  {
    static const addPlan_args args[] = {
      {4, "Har",  " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0},
      {4, "Har",  " H5",  " C5",  " C4",  " C6",  "",    1.1,   0.0,   0.0, 0},
      {3, "Hpol", "2H4",  " N4",  " C4",  " N3",  "",    1.0,-120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "1H4",  " N4",  " C4",  " N3",  "",    1.0, 120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {7, "Hpol", " H42", " N4",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {7, "Hpol", " H41", " N4",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {3, "Hpol", " H42",  " N4",  " C4",  " N3",  "",    1.0,-120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", " H41",  " N4",  " C4",  " N3",  "",    1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {1, "H",    " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USEOLDNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, XPLORNAME},
      {1, "H",    " H1'", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH(" DC", "", args);
  }
  {
    static const addPlan_args args[] = {
      {3, "Hpol", "2H2",  " N2",  " C2",  " N1",  "",    1.0,-120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", "1H2",  " N2",  " C2",  " N1",  "",    1.0, 120.0, 180.0, USEOLDNAMES|XPLORNAME|BONDBUMPFLAG},
      {7, "Hpol", " H22", " N2",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {7, "Hpol", " H21", " N2",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {3, "Hpol", " H22",  " N2",  " C2",  " N1",  "",    1.0,-120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", " H21",  " N2",  " C2",  " N1",  "",    1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {4, "Ha+p", " H1",  " N1",  " C6",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG},
      {4, "Har",  " H8",  " C8",  " N7",  " N9",  "",    1.1,   0.0,   0.0, 0},
      {1, "H",    " H1*", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, USEOLDNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, XPLORNAME},
      {1, "H",    " H1'", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, USENEWNAMES},
      {1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH(" DG", "", args);
  }

//--------------------------------------------------------------------------
 {
    static const addPlan_args args[] = {
      {3, "H", "3HB2", " CB2", " CA", " N", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2HB2", " CB2", " CA", " N", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1HB2", " CB2", " CA", " N", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "3HB1", " CB1", " CA", " N", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2HB1", " CB1", " CA", " N", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1HB1", " CB1", " CA", " N", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "HB23", " CB2", " CA", " N", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HB22", " CB2", " CA", " N", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HB21", " CB2", " CA", " N", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HB13", " CB1", " CA", " N", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HB12", " CB1", " CA", " N", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HB11", " CB1", " CA", " N", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("AIB", "", args);
  }
  {
    static const addPlan_args args[] = {
    //  {3, "H", "3HG", " CG", " CB", " CA", "", 1.1, 109.5,  60.0, ROTATEONDEMAND}, //check old name.
      {3, "H", " HE2", " CG", " CB", " CA", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2HG", " CG", " CB", " CA", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1HG", " CG", " CB", " CA", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {3, "H", " HE2", " CG", " CB", " CA", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", " HG2", " CG", " CB", " CA", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", " HG1", " CG", " CB", " CA", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {2, "H", " HB2", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB1", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("ABU", "", args);
  }
  {
    static const addPlan_args args[] = {
  //    {3, "H", "3HH3", " CH3", " C", " O", "", 1.1, 109.5,   0.0, ROTATEONDEMAND},
  //    {3, "H", "2HH3", " CH3", " C", " O", "", 1.1, 109.5,-120.0, ROTATEONDEMAND},
  //    {3, "H", "1HH3", " CH3", " C", " O", "", 1.1, 109.5, 120.0, ROTATEONDEMAND},
     {3, "H", "3H", " CH3", " C", " O", "", 1.1, 109.5,   0.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
     {3, "H", "2H", " CH3", " C", " O", "", 1.1, 109.5,-120.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
     {3, "H", "1H", " CH3", " C", " O", "", 1.1, 109.5, 120.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
     {3, "H", " H3", " CH3", " C", " O", "", 1.1, 109.5,   0.0, USENEWNAMES|ROTATEONDEMAND},
     {3, "H", " H2", " CH3", " C", " O", "", 1.1, 109.5,-120.0, USENEWNAMES|ROTATEONDEMAND},
     {3, "H", " H1", " CH3", " C", " O", "", 1.1, 109.5, 120.0, USENEWNAMES|ROTATEONDEMAND}, 
     {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("ACE", "", args);
  }
  {
    static const addPlan_args args[] = {
      {2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", " HB2", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB1", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("ASX", "", args);
  }
  {
    static const addPlan_args args[] = {
      {2, "H", "2HG", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HG", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", " HG2", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HG1", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB1", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("GLX", "", args);
  }
  {
    static const addPlan_args args[] = {
      {3, "H", "3HE", " CE", "SED", " CG", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2HE", " CE", "SED", " CG", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1HE", " CE", "SED", " CG", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {2, "H", "2HG", " CG", " CB", "SED", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HG", " CG", " CB", "SED", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {3, "H", " HE3", " CE", "SED", " CG", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", " HE2", " CE", "SED", " CG", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", " HE1", " CE", "SED", " CG", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", " HE3", " CE", "SE", " CG", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", " HE2", " CE", "SE", " CG", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", " HE1", " CE", "SE", " CG", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {2, "H", " HG3", " CG", " CB", "SED", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HG2", " CG", " CB", "SED", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", " HG3", " CG", " CB", "SE", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HG2", " CG", " CB", "SE", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", " HB3", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("MSE", "", args);
  }
  {
    static const addPlan_args args[] = {
      {2, "H", "2HG", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HG", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", " HG3", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HG2", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {2, "H", " HB3", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, USENEWNAMES},
      {2, "H", " HB2", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("PCA", "", args);
  }
  {
    static const addPlan_args args[] = {
      {3, "Hpol", " H2", " N", "- C", "- O", "", 1.0, 120.0,   0.0, XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", " HN1", " N", "- C", "- O", "", 1.0, 120.0, 180.0, XPLORNAME|BONDBUMPFLAG},
      {3, "Hpol", " H2", " N", "- C", "- O", "", 1.0, 120.0,   0.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", " HN1", " N", "- C", "- O", "", 1.0, 120.0, 180.0, USENEWNAMES|BONDBUMPFLAG},
      {3, "Hpol", "2HN",  " N", "- C", "- O", "", 1.0, 120.0,   0.0, USEOLDNAMES|BONDBUMPFLAG},
      {3, "Hpol", "1HN",  " N", "- C", "- O", "", 1.0, 120.0, 180.0, USEOLDNAMES|BONDBUMPFLAG},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("NH2", "", args);
  }
  {
    static const addPlan_args args[] = {
  //    {3, "H", "3HH3", " CH3", " N",   "- C", "", 1.1, 109.5,  60.0, ROTATEONDEMAND},
  //    {3, "H", "2HH3", " CH3", " N",   "- C", "", 1.1, 109.5, -60.0, ROTATEONDEMAND},
  //    {3, "H", "1HH3", " CH3", " N",   "- C", "", 1.1, 109.5, 180.0, ROTATEONDEMAND},
  //    {4, "Hpol", " H",   " N",   " CH3", "- C", "", 1.0,   0.0,   0.02,0},
      {3, "H", "3H", " CH3", " N",   "- C", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2H", " CH3", " N",   "- C", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1H", " CH3", " N",   "- C", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {4, "Hpol", "2HN",   " N",   " CH3", "- C", "", 1.0,   0.0,   0.02,USEOLDNAMES|XPLORNAME},
//??  {4, "Hpol", "1HN",   " N",   " CH3", "- C", "", 1.0,   0.0,   0.02,USEOLDNAMES|XPLORNAME},
      {3, "H", " H3", " CH3", " N",   "- C", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", " H2", " CH3", " N",   "- C", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", " H1", " CH3", " N",   "- C", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {4, "Hpol", " HN2",   " N",   " CH3", "- C", "", 1.0,   0.0,   0.02, USENEWNAMES},
//??  {4, "Hpol", " HN1",   " N",   " CH3", "- C", "", 1.0,   0.0,   0.02, USEOLDNAMES|XPLORNAME},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("NME", "", args);
  }
  {
    static const addPlan_args args[] = {
  //    {4, "H", " H", " C", " O", "+ N", "", 1.1,   0.0,   0.0, 0},
      {4, "H", "2H", " C", " O", "+ N", "", 1.1,   0.0,   0.0, USEOLDNAMES|XPLORNAME},
      {4, "H", "1H", " C", " O", "+ N", "", 1.1,   0.0,   0.0, USEOLDNAMES|XPLORNAME},
      {4, "H", " H2", " C", " O", "+ N", "", 1.1,   0.0,   0.0, USENEWNAMES},
      {4, "H", " H1", " C", " O", "+ N", "", 1.1,   0.0,   0.0, USENEWNAMES},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("FOR", "", args);
  }
  {
    static const addPlan_args args[] = {
      {2, "H", "2HAD", " CAD", " C3D", " CBD", "", 1.1,-126.5,  0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HAD", " CAD", " C3D", " CBD", "", 1.1, 126.5,  0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "2HBD", " CBD", " CAD", " CGD", "", 1.1,-126.5,  0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HBD", " CBD", " CAD", " CGD", "", 1.1, 126.5,  0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "2HAA", " CAA", " C2A", " CBA", "", 1.1,-126.5,  0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HAA", " CAA", " C2A", " CBA", "", 1.1, 126.5,  0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "2HBA", " CBA", " CAA", " CGA", "", 1.1,-126.5,  0.0, USEOLDNAMES|XPLORNAME},
      {2, "H", "1HBA", " CBA", " CAA", " CGA", "", 1.1, 126.5,  0.0, USEOLDNAMES|XPLORNAME},
      {3, "H", "2HBC", " CBC", " CAC", " C3C", "", 1.1, 120.0,   0.0, USEOLDNAMES|XPLORNAME},
      {3, "H", "1HBC", " CBC", " CAC", " C3C", "", 1.1, 120.0, 180.0, USEOLDNAMES|XPLORNAME},
      {4, "H", " HAC", " CAC", " CBC", " C3C", "", 1.1,   0.0,   0.0, 0},
      {3, "H", "2HBB", " CBB", " CAB", " C3B", "", 1.1, 120.0,   0.0, USEOLDNAMES|XPLORNAME},
      {3, "H", "1HBB", " CBB", " CAB", " C3B", "", 1.1, 120.0, 180.0, USEOLDNAMES|XPLORNAME},
      {4, "H", " HAB", " CAB", " CBB", " C3B", "", 1.1,   0.0,   0.0, 0},
      {3, "H", "3HMD", " CMD", " C2D", " C1D", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2HMD", " CMD", " C2D", " C1D", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1HMD", " CMD", " C2D", " C1D", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "3HMC", " CMC", " C2C", " C1C", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2HMC", " CMC", " C2C", " C1C", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1HMC", " CMC", " C2C", " C1C", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "3HMB", " CMB", " C2B", " C1B", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2HMB", " CMB", " C2B", " C1B", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1HMB", " CMB", " C2B", " C1B", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "3HMA", " CMA", " C3A", " C4A", "", 1.1, 109.5, -60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "2HMA", " CMA", " C3A", " C4A", "", 1.1, 109.5,  60.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {3, "H", "1HMA", " CMA", " C3A", " C4A", "", 1.1, 109.5, 180.0, USEOLDNAMES|XPLORNAME|ROTATEONDEMAND},
      {4, "Har", " HHD", " CHD", " C1D", " C4C", "", 1.1,   0.0,   0.0, USENEWNAMES},
      {4, "Har", " HHC", " CHC", " C1C", " C4B", "", 1.1,   0.0,   0.0, USENEWNAMES},
      {4, "Har", " HHB", " CHB", " C1B", " C4A", "", 1.1,   0.0,   0.0, USENEWNAMES},
      {4, "Har", " HHA", " CHA", " C1A", " C4D", "", 1.1,   0.0,   0.0, USENEWNAMES},
      {2, "H", "HAD2", " CAD", " C3D", " CBD", "", 1.1,-126.5,  0.0, USENEWNAMES},
      {2, "H", "HAD1", " CAD", " C3D", " CBD", "", 1.1, 126.5,  0.0, USENEWNAMES},
      {2, "H", "HBD2", " CBD", " CAD", " CGD", "", 1.1,-126.5,  0.0, USENEWNAMES},
      {2, "H", "HBD1", " CBD", " CAD", " CGD", "", 1.1, 126.5,  0.0, USENEWNAMES},
      {2, "H", "HAA2", " CAA", " C2A", " CBA", "", 1.1,-126.5,  0.0, USENEWNAMES},
      {2, "H", "HAA1", " CAA", " C2A", " CBA", "", 1.1, 126.5,  0.0, USENEWNAMES},
      {2, "H", "HBA2", " CBA", " CAA", " CGA", "", 1.1,-126.5,  0.0, USENEWNAMES},
      {2, "H", "HBA1", " CBA", " CAA", " CGA", "", 1.1, 126.5,  0.0, USENEWNAMES},
      {3, "H", "HBC2", " CBC", " CAC", " C3C", "", 1.1, 120.0,   0.0, USENEWNAMES},
      {3, "H", "HBC1", " CBC", " CAC", " C3C", "", 1.1, 120.0, 180.0, USENEWNAMES},
      {3, "H", "HBB2", " CBB", " CAB", " C3B", "", 1.1, 120.0,   0.0, USENEWNAMES},
      {3, "H", "HBB1", " CBB", " CAB", " C3B", "", 1.1, 120.0, 180.0, USENEWNAMES},
      {3, "H", "HMD3", " CMD", " C2D", " C1D", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HMD2", " CMD", " C2D", " C1D", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HMD1", " CMD", " C2D", " C1D", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HMC3", " CMC", " C2C", " C1C", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HMC2", " CMC", " C2C", " C1C", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HMC1", " CMC", " C2C", " C1C", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HMB3", " CMB", " C2B", " C1B", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HMB2", " CMB", " C2B", " C1B", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HMB1", " CMB", " C2B", " C1B", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HMA3", " CMA", " C3A", " C4A", "", 1.1, 109.5, -60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HMA2", " CMA", " C3A", " C4A", "", 1.1, 109.5,  60.0, USENEWNAMES|ROTATEONDEMAND},
      {3, "H", "HMA1", " CMA", " C3A", " C4A", "", 1.1, 109.5, 180.0, USENEWNAMES|ROTATEONDEMAND},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("HEM", "", args);
  }
//--------------------------------------------------------------------------
  {
    static const addPlan_args args[] = {
      {7, "Hpol", " H2", " OD2", "", "", "", 1.0,   0.0,   0.0, USENEWNAMES|UNSUREDROPFLAG},
      {7, "Hpol", " H1", " OD2", "", "", "", 1.0,   0.0,   0.0, USENEWNAMES|UNSUREDROPFLAG},
      {7, "Hpol", " H2", " OH2", "", "", "", 1.0,   0.0,   0.0, USENEWNAMES|UNSUREDROPFLAG},
      {7, "Hpol", " H1", " OH2", "", "", "", 1.0,   0.0,   0.0, USENEWNAMES|UNSUREDROPFLAG},
      {7, "Hpol", " H2", " OW",  "", "", "", 1.0,   0.0,   0.0, USENEWNAMES|UNSUREDROPFLAG},
      {7, "Hpol", " H1", " OW",  "", "", "", 1.0,   0.0,   0.0, USENEWNAMES|UNSUREDROPFLAG},
      {7, "Hpol", " H2", " O",   "", "", "", 1.0,   0.0,   0.0, USENEWNAMES|UNSUREDROPFLAG},
      {7, "Hpol", " H1", " O",   "", "", "", 1.0,   0.0,   0.0, USENEWNAMES|UNSUREDROPFLAG},
      // alternative naming of hydrogens
      {7, "Hpol", "2H", " OD2", "", "", "", 1.0,   0.0,   0.0, USEOLDNAMES|XPLORNAME|UNSUREDROPFLAG},
      {7, "Hpol", "1H", " OD2", "", "", "", 1.0,   0.0,   0.0, USEOLDNAMES|XPLORNAME|UNSUREDROPFLAG},
      {7, "Hpol", "2H", " OH2", "", "", "", 1.0,   0.0,   0.0, USEOLDNAMES|XPLORNAME|UNSUREDROPFLAG},
      {7, "Hpol", "1H", " OH2", "", "", "", 1.0,   0.0,   0.0, USEOLDNAMES|XPLORNAME|UNSUREDROPFLAG},
      {7, "Hpol", "2H", " OW",  "", "", "", 1.0,   0.0,   0.0, USEOLDNAMES|XPLORNAME|UNSUREDROPFLAG},
      {7, "Hpol", "1H", " OW",  "", "", "", 1.0,   0.0,   0.0, USEOLDNAMES|XPLORNAME|UNSUREDROPFLAG},
      {7, "Hpol", "2H", " O",   "", "", "", 1.0,   0.0,   0.0, USEOLDNAMES|XPLORNAME|UNSUREDROPFLAG},
      {7, "Hpol", "1H", " O",   "", "", "", 1.0,   0.0,   0.0, USEOLDNAMES|XPLORNAME|UNSUREDROPFLAG},
      {0,0,0,0,0,0,0,0,0,0,0}
    };
    insertStdResH("HOH", "", args);
  }
}

// initialize xtra info about the standard residues
StdResXtraInfo::StdResXtraInfo() {
  struct raw_info_t {
    const char* resName;
    const char* atomName;
    int flag;
  };
  static const raw_info_t raw_info[] = {
    {"GLY","AminoAcid", 1},
    {"ALA","AminoAcid", 1},
    {"VAL","AminoAcid", 1},
    {"PHE","AminoAcid", 1},
    {"PRO","AminoAcid", 1},
    {"MET","AminoAcid", 1},
    {"ILE","AminoAcid", 1},
    {"LEU","AminoAcid", 1},
    {"ASP","AminoAcid", 1},
    {"GLU","AminoAcid", 1},
    {"LYS","AminoAcid", 1},
    {"ARG","AminoAcid", 1},
    {"SER","AminoAcid", 1},
    {"THR","AminoAcid", 1},
    {"TYR","AminoAcid", 1},
    {"HIS","AminoAcid", 1},
    {"CYS","AminoAcid", 1},
    {"ASN","AminoAcid", 1},
    {"GLN","AminoAcid", 1},
    {"TRP","AminoAcid", 1},
    {"ASX","AminoAcid", 1},
    {"GLX","AminoAcid", 1},
    {"ABU","AminoAcid", 1},
    {"AIB","AminoAcid", 1},
    {"ABU","AminoAcid", 1},
    {"MSE","AminoAcid", 1},
    {"PCA","AminoAcid", 1},

    {"  U","NucleicAcid", 1},
#ifdef LEFT_JUSTIFY_NUC_RES_OK
    {"U",  "NucleicAcid", 1},
#endif
    {"URA","NucleicAcid", 1},
    {"  T","NucleicAcid", 1},
#ifdef LEFT_JUSTIFY_NUC_RES_OK
    {"T",  "NucleicAcid", 1},
#endif
    {"THY","NucleicAcid", 1},
    {"  A","NucleicAcid", 1},
#ifdef LEFT_JUSTIFY_NUC_RES_OK
    {"A",  "NucleicAcid", 1},
#endif
    {"ADE","NucleicAcid", 1},
    {"  C","NucleicAcid", 1},
#ifdef LEFT_JUSTIFY_NUC_RES_OK
    {"C",  "NucleicAcid", 1},
#endif
    {"CYT","NucleicAcid", 1},
    {"  G","NucleicAcid", 1},
#ifdef LEFT_JUSTIFY_NUC_RES_OK
    {"G",  "NucleicAcid", 1},
#endif
    {"GUA","NucleicAcid", 1},

//    {"CTP","NucleicAcid", 1},
//    {"CDP","NucleicAcid", 1},
//    {"CMP","NucleicAcid", 1},
//    {"GDP","NucleicAcid", 1},
//    {"GMP","NucleicAcid", 1},
//    {"ATP","NucleicAcid", 1},
//    {"ADP","NucleicAcid", 1},
//    {"AMP","NucleicAcid", 1},
//    {"TTP","NucleicAcid", 1},
//    {"TDP","NucleicAcid", 1},
//    {"TMP","NucleicAcid", 1},
//    {"UTP","NucleicAcid", 1},
//    {"UDP","NucleicAcid", 1},
//    {"UMP","NucleicAcid", 1},
//    {"GSP","NucleicAcid", 1},
    {" DA","NucleicAcid", 1},
    {" DT","NucleicAcid", 1},
    {" DC","NucleicAcid", 1},
    {" DG","NucleicAcid", 1},

// (reduced VDW radii indicated for C=O carbon atoms
//  in-addition to AminoAcid mainchain COs

    {"ASP", " CG", ISACOFLAG},
    {"ASN", " CG", ISACOFLAG},
    {"ASX", " CG", ISACOFLAG},
    {"GLU", " CD", ISACOFLAG},
    {"GLN", " CD", ISACOFLAG},
    {"GLX", " CD", ISACOFLAG},

// (setting the HBACCEPTORFLAG & AROMATICFLAG for certain carbons & nitrogens

     // HB status of HIS may depend on the protonation state
    {"HIS", " ND1",HBACCEPTORFLAG|AROMATICFLAG},
    {"HIS", " NE2",HBACCEPTORFLAG|AROMATICFLAG},

    {"  A", " N1", HBACCEPTORFLAG|AROMATICFLAG},
    {"  A", " N3", HBACCEPTORFLAG|AROMATICFLAG},
    {"  A", " N7", HBACCEPTORFLAG|AROMATICFLAG},
    {"  C", " N3", HBACCEPTORFLAG|AROMATICFLAG},
    {"  G", " N3", HBACCEPTORFLAG|AROMATICFLAG},
    {"  G", " N7", HBACCEPTORFLAG|AROMATICFLAG},
#ifdef LEFT_JUSTIFY_NUC_RES_OK
    {"A",   " N1", HBACCEPTORFLAG|AROMATICFLAG},
    {"A",   " N3", HBACCEPTORFLAG|AROMATICFLAG},
    {"A",   " N7", HBACCEPTORFLAG|AROMATICFLAG},
    {"C",   " N3", HBACCEPTORFLAG|AROMATICFLAG},
    {"G",   " N3", HBACCEPTORFLAG|AROMATICFLAG},
    {"G",   " N7", HBACCEPTORFLAG|AROMATICFLAG},
#endif
    {"ADE", " N1", HBACCEPTORFLAG|AROMATICFLAG},
    {"ADE", " N3", HBACCEPTORFLAG|AROMATICFLAG},
    {"ADE", " N7", HBACCEPTORFLAG|AROMATICFLAG},
    {"CYT", " N3", HBACCEPTORFLAG|AROMATICFLAG},
    {"GUA", " N3", HBACCEPTORFLAG|AROMATICFLAG},
    {"GUA", " N7", HBACCEPTORFLAG|AROMATICFLAG},

// new DNA names rmi070719
    {" DA", " N1", HBACCEPTORFLAG|AROMATICFLAG},
    {" DA", " N3", HBACCEPTORFLAG|AROMATICFLAG},
    {" DA", " N7", HBACCEPTORFLAG|AROMATICFLAG},
    {" DC", " N3", HBACCEPTORFLAG|AROMATICFLAG},
    {" DG", " N3", HBACCEPTORFLAG|AROMATICFLAG},
    {" DG", " N7", HBACCEPTORFLAG|AROMATICFLAG},

// other HEM atom atributes lower down in file
    {"HEM", " N A", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " N B", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " N C", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " N D", HBACCEPTORFLAG|AROMATICFLAG},

// always assumed to be charged...
    {"AminoAcid",   " NT",  POSCHARGEFLAG},
    {"AminoAcid",   "1H"  , POSCHARGEFLAG},
    {"AminoAcid",   " H1"  , POSCHARGEFLAG},
    {"AminoAcid",   "2H"  , POSCHARGEFLAG},
    {"AminoAcid",   " H2"  , POSCHARGEFLAG},
    {"AminoAcid",   "3H"  , POSCHARGEFLAG},
    {"AminoAcid",   " H3"  , POSCHARGEFLAG},
    {"AminoAcid",   " HT1", POSCHARGEFLAG},
    {"AminoAcid",   " HT2", POSCHARGEFLAG},
    {"AminoAcid",   " HT3", POSCHARGEFLAG},
    {"AminoAcid",   "1D"  , POSCHARGEFLAG},
    {"AminoAcid",   " D1"  , POSCHARGEFLAG},
    {"AminoAcid",   "2D"  , POSCHARGEFLAG},
    {"AminoAcid",   " D2"  , POSCHARGEFLAG},
    {"AminoAcid",   "3D"  , POSCHARGEFLAG},
    {"AminoAcid",   " D3"  , POSCHARGEFLAG},
    {"AminoAcid",   " DT1", POSCHARGEFLAG},
    {"AminoAcid",   " DT2", POSCHARGEFLAG},
    {"AminoAcid",   " DT3", POSCHARGEFLAG},

    {"AminoAcid",   " OXT", NEGCHARGEFLAG},
    {"AminoAcid",   "1OXT", NEGCHARGEFLAG},
    {"AminoAcid",   "2OXT", NEGCHARGEFLAG},

    {"NucleicAcid", " P"  , NEGCHARGEFLAG},
    {"NucleicAcid", " O1P", NEGCHARGEFLAG},
    {"NucleicAcid", " OP1", NEGCHARGEFLAG},
    {"NucleicAcid", " O2P", NEGCHARGEFLAG},
    {"NucleicAcid", " OP2", NEGCHARGEFLAG},
    {"NucleicAcid", " PA" , NEGCHARGEFLAG},
    {"NucleicAcid", " PB" , NEGCHARGEFLAG},
    {"NucleicAcid", " PG" , NEGCHARGEFLAG},
    {"NucleicAcid", " O1A", NEGCHARGEFLAG},
    {"NucleicAcid", " O2A", NEGCHARGEFLAG},
    {"NucleicAcid", " O3A", NEGCHARGEFLAG},
    {"NucleicAcid", " O1B", NEGCHARGEFLAG},
    {"NucleicAcid", " O2B", NEGCHARGEFLAG},
    {"NucleicAcid", " O3B", NEGCHARGEFLAG},
    {"NucleicAcid", " O1G", NEGCHARGEFLAG},
    {"NucleicAcid", " O2G", NEGCHARGEFLAG},
    {"NucleicAcid", " O3G", NEGCHARGEFLAG},
    {"NucleicAcid", " S1G", NEGCHARGEFLAG},

    {"ASP","ChargedResidue", 1},
    {"GLU","ChargedResidue", 1},
    {"LYS","ChargedResidue", 1},
    {"ARG","ChargedResidue", 1},
    {"HEM","ChargedResidue", 1},

    {"ASP", " OD1",  NEGCHARGEFLAG},
    {"ASP", " OD2",  NEGCHARGEFLAG},
    {"GLU", " OE1",  NEGCHARGEFLAG},
    {"GLU", " OE2",  NEGCHARGEFLAG},

    {"LYS", " NZ" ,  POSCHARGEFLAG},
    {"LYS", "1HZ" ,  POSCHARGEFLAG},
    {"LYS", "2HZ" ,  POSCHARGEFLAG},
    {"LYS", "3HZ" ,  POSCHARGEFLAG},
    {"LYS", " HZ1",  POSCHARGEFLAG},
    {"LYS", " HZ2",  POSCHARGEFLAG},
    {"LYS", " HZ3",  POSCHARGEFLAG},
    {"LYS", "1DZ" ,  POSCHARGEFLAG},
    {"LYS", "2DZ" ,  POSCHARGEFLAG},
    {"LYS", "3DZ" ,  POSCHARGEFLAG},
    {"LYS", " DZ1",  POSCHARGEFLAG},
    {"LYS", " DZ2",  POSCHARGEFLAG},
    {"LYS", " DZ3",  POSCHARGEFLAG},

    {"ARG", " NE" ,  POSCHARGEFLAG},
    {"ARG", " NH1",  POSCHARGEFLAG},
    {"ARG", " NH2",  POSCHARGEFLAG},
    {"ARG", " HE" ,  POSCHARGEFLAG},
    {"ARG", "1HH1",  POSCHARGEFLAG},
    {"ARG", "2HH1",  POSCHARGEFLAG},
    {"ARG", "1HH2",  POSCHARGEFLAG},
    {"ARG", "2HH2",  POSCHARGEFLAG},
    {"ARG", "HH11",  POSCHARGEFLAG},
    {"ARG", "HH12",  POSCHARGEFLAG},
    {"ARG", "HH21",  POSCHARGEFLAG},
    {"ARG", "HH22",  POSCHARGEFLAG},
    {"ARG", "1DH1",  POSCHARGEFLAG},
    {"ARG", "2DH1",  POSCHARGEFLAG},
    {"ARG", "1DH2",  POSCHARGEFLAG},
    {"ARG", "2DH2",  POSCHARGEFLAG},
    {"ARG", "DH11",  POSCHARGEFLAG},
    {"ARG", "DH12",  POSCHARGEFLAG},
    {"ARG", "DH21",  POSCHARGEFLAG},
    {"ARG", "DH22",  POSCHARGEFLAG},

    {"HEM", " O1A",  NEGCHARGEFLAG},
    {"HEM", " O2A",  NEGCHARGEFLAG},
    {"HEM", " O1D",  NEGCHARGEFLAG},
    {"HEM", " O2D",  NEGCHARGEFLAG},
#ifdef AROMATICS_ACCEPT_HBONDS
// -----------------------------------------------------------
// here we treat the aromatic Pi-bonds as hydrogen bond acceptors
// -----------------------------------------------------------
    {"HEM", " C1A", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " C2A", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " C3A", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " C4A", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " C1B", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " C2B", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " C3B", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " C4B", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " C1C", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " C2C", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " C3C", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " C4C", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " C1D", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " C2D", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " C3D", HBACCEPTORFLAG|AROMATICFLAG},
    {"HEM", " C4D", HBACCEPTORFLAG|AROMATICFLAG},

    {"PHE", " CZ",  HBACCEPTORFLAG|AROMATICFLAG},
    {"PHE", " CE2", HBACCEPTORFLAG|AROMATICFLAG},
    {"PHE", " CE1", HBACCEPTORFLAG|AROMATICFLAG},
    {"PHE", " CD2", HBACCEPTORFLAG|AROMATICFLAG},
    {"PHE", " CD1", HBACCEPTORFLAG|AROMATICFLAG},
    {"PHE", " CG",  HBACCEPTORFLAG|AROMATICFLAG},

    {"TYR", " CZ",  HBACCEPTORFLAG|AROMATICFLAG},
    {"TYR", " CE2", HBACCEPTORFLAG|AROMATICFLAG},
    {"TYR", " CE1", HBACCEPTORFLAG|AROMATICFLAG},
    {"TYR", " CD2", HBACCEPTORFLAG|AROMATICFLAG},
    {"TYR", " CD1", HBACCEPTORFLAG|AROMATICFLAG},
    {"TYR", " CG",  HBACCEPTORFLAG|AROMATICFLAG},

    //{"HIS", " CD2", HBACCEPTORFLAG|AROMATICFLAG},
    //{"HIS", " CE1", HBACCEPTORFLAG|AROMATICFLAG},
    //{"HIS", " CG",  HBACCEPTORFLAG|AROMATICFLAG},
    {"HIS", " CD2", 0},
    {"HIS", " CE1", 0},
    {"HIS", " CG",  0},

    {"TRP", " CH2", HBACCEPTORFLAG|AROMATICFLAG},
    {"TRP", " CZ3", HBACCEPTORFLAG|AROMATICFLAG},
    {"TRP", " CZ2", HBACCEPTORFLAG|AROMATICFLAG},
    {"TRP", " CE3", HBACCEPTORFLAG|AROMATICFLAG},
    {"TRP", " CE2", HBACCEPTORFLAG|AROMATICFLAG},
    {"TRP", " NE1", HBACCEPTORFLAG|AROMATICFLAG},
    {"TRP", " CD2", HBACCEPTORFLAG|AROMATICFLAG},
    {"TRP", " CD1", HBACCEPTORFLAG|AROMATICFLAG},
    {"TRP", " CG",  HBACCEPTORFLAG|AROMATICFLAG},

// -----------------------------------------------------------
// pick up the parts of the bases not included above
// *** warning *** the bases below are missing the atoms
//                 which are listed above (atoms not duplicated
// -----------------------------------------------------------
    {"  A", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  A", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  A", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  A", " C6",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  A", " C8",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  A", " N9",  HBACCEPTORFLAG|AROMATICFLAG},

    {"  C", " N1",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  C", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  C", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  C", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  C", " C6",  HBACCEPTORFLAG|AROMATICFLAG},

    {"  G", " N1",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  G", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  G", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  G", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  G", " C6",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  G", " C8",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  G", " N9",  HBACCEPTORFLAG|AROMATICFLAG},

    {"  T", " N1",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  T", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  T", " N3",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  T", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  T", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  T", " C6",  HBACCEPTORFLAG|AROMATICFLAG},

    {"  U", " N1",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  U", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  U", " N3",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  U", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  U", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {"  U", " C6",  HBACCEPTORFLAG|AROMATICFLAG},

#ifdef LEFT_JUSTIFY_NUC_RES_OK
    {"A", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {"A", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {"A", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {"A", " C6",  HBACCEPTORFLAG|AROMATICFLAG},
    {"A", " C8",  HBACCEPTORFLAG|AROMATICFLAG},
    {"A", " N9",  HBACCEPTORFLAG|AROMATICFLAG},

    {"C", " N1",  HBACCEPTORFLAG|AROMATICFLAG},
    {"C", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {"C", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {"C", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {"C", " C6",  HBACCEPTORFLAG|AROMATICFLAG},

    {"G", " N1",  HBACCEPTORFLAG|AROMATICFLAG},
    {"G", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {"G", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {"G", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {"G", " C6",  HBACCEPTORFLAG|AROMATICFLAG},
    {"G", " C8",  HBACCEPTORFLAG|AROMATICFLAG},
    {"G", " N9",  HBACCEPTORFLAG|AROMATICFLAG},

    {"T", " N1",  HBACCEPTORFLAG|AROMATICFLAG},
    {"T", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {"T", " N3",  HBACCEPTORFLAG|AROMATICFLAG},
    {"T", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {"T", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {"T", " C6",  HBACCEPTORFLAG|AROMATICFLAG},

    {"U", " N1",  HBACCEPTORFLAG|AROMATICFLAG},
    {"U", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {"U", " N3",  HBACCEPTORFLAG|AROMATICFLAG},
    {"U", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {"U", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {"U", " C6",  HBACCEPTORFLAG|AROMATICFLAG},
#endif

    {"ADE", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {"ADE", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {"ADE", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {"ADE", " C6",  HBACCEPTORFLAG|AROMATICFLAG},
    {"ADE", " C8",  HBACCEPTORFLAG|AROMATICFLAG},
    {"ADE", " N9",  HBACCEPTORFLAG|AROMATICFLAG},

    {"CYT", " N1",  HBACCEPTORFLAG|AROMATICFLAG},
    {"CYT", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {"CYT", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {"CYT", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {"CYT", " C6",  HBACCEPTORFLAG|AROMATICFLAG},

    {"GUA", " N1",  HBACCEPTORFLAG|AROMATICFLAG},
    {"GUA", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {"GUA", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {"GUA", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {"GUA", " C6",  HBACCEPTORFLAG|AROMATICFLAG},
    {"GUA", " C8",  HBACCEPTORFLAG|AROMATICFLAG},
    {"GUA", " N9",  HBACCEPTORFLAG|AROMATICFLAG},

    {"THY", " N1",  HBACCEPTORFLAG|AROMATICFLAG},
    {"THY", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {"THY", " N3",  HBACCEPTORFLAG|AROMATICFLAG},
    {"THY", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {"THY", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {"THY", " C6",  HBACCEPTORFLAG|AROMATICFLAG},

    {"URA", " N1",  HBACCEPTORFLAG|AROMATICFLAG},
    {"URA", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {"URA", " N3",  HBACCEPTORFLAG|AROMATICFLAG},
    {"URA", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {"URA", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {"URA", " C6",  HBACCEPTORFLAG|AROMATICFLAG},
#endif

    {" DA", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DA", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DA", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DA", " C6",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DA", " C8",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DA", " N9",  HBACCEPTORFLAG|AROMATICFLAG},

    {" DC", " N1",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DC", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DC", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DC", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DC", " C6",  HBACCEPTORFLAG|AROMATICFLAG},

    {" DG", " N1",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DG", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DG", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DG", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DG", " C6",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DG", " C8",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DG", " N9",  HBACCEPTORFLAG|AROMATICFLAG},

    {" DT", " N1",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DT", " C2",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DT", " N3",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DT", " C4",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DT", " C5",  HBACCEPTORFLAG|AROMATICFLAG},
    {" DT", " C6",  HBACCEPTORFLAG|AROMATICFLAG},

    {0, 0, 0}
  };
  for(const raw_info_t* r=raw_info;r->resName;r++) {
    _atomAttributes.insert(
      std::make_pair(
        makeKey(r->resName,r->atomName),
        r->flag));
  }
}
