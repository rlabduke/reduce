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

 HydrogenPlanTable::HydrogenPlanTable() {
   StdResH *theRes = NULL;
//--------------------------------------------------------------------------
// the order of plans within each residue is the reverse of the order we want
// them to be generated because they will go into a sequence.
//--------------------------------------------------------------------------
   theRes = new StdResH("amide", "PRO"); // *NON* N-terminal mc amide
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(5, "Hpol", " H", " N", " CA", "- C", "", 1.0,   0.0, 0.48, BONDBUMPFLAG);

   theRes = new StdResH("nt-amide", "PRO"); // N-terminal mc amide
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "Hpol", " HT3", " N", " CA", " C", "", 1.0, 109.5,  60.0, XPLORNAME|ROTATEONDEMAND|NH3FLAG);
   theRes->addPlan(3, "Hpol", " HT2", " N", " CA", " C", "", 1.0, 109.5, -60.0, XPLORNAME|ROTATEONDEMAND|NH3FLAG);
   theRes->addPlan(3, "Hpol", " HT1", " N", " CA", " C", "", 1.0, 109.5, 180.0, XPLORNAME|ROTATEONDEMAND|NH3FLAG);
   theRes->addPlan(3, "Hpol", "3H",   " N", " CA", " C", "", 1.0, 109.5,  60.0, NOTXPLORNAME|ROTATEONDEMAND|NH3FLAG);
   theRes->addPlan(3, "Hpol", "2H",   " N", " CA", " C", "", 1.0, 109.5, -60.0, NOTXPLORNAME|ROTATEONDEMAND|NH3FLAG);
   theRes->addPlan(3, "Hpol", "1H",   " N", " CA", " C", "", 1.0, 109.5, 180.0, NOTXPLORNAME|ROTATEONDEMAND|NH3FLAG);

   theRes = new StdResH("nt-pro", ""); // N-terminal prolines
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(2, "Hpol", " HT2", " N", " CD", " CA", "", 1.0, 126.5,   0.0, XPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(2, "Hpol", " HT1", " N", " CD", " CA", "", 1.0,-126.5,   0.0, XPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(2, "Hpol", "2H",   " N", " CD", " CA", "", 1.0, 126.5,   0.0, NOTXPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(2, "Hpol", "1H",   " N", " CD", " CA", "", 1.0,-126.5,   0.0, NOTXPLORNAME|BONDBUMPFLAG);

   theRes = new StdResH("alpha", "GLY"); // mainchain alpha proton
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(1, "H", " HA", " CA", " N", " C", " CB", 1.1,   0.0,   0.0, STRICTALTFLAG);

   theRes = new StdResH("ribose phosphate backbone", ""); // phosphate chain (two different naming schemes)
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "Hpol", "5HO*", " O5*", " C5*", " C4*", "",     1.0, 109.5, 180.0, UNSUREDROPFLAG|ROTATEFLAG|IFNOPO4);
   theRes->addPlan(3, "Hpol", "3HO*", " O3*", " C3*", " C4*", "",     1.0, 109.5, 180.0, UNSUREDROPFLAG|ROTATEFLAG|IFNOPO4);
   theRes->addPlan(3, "Hpol", "2HO*", " O2*", " C2*", " C3*", "",     1.0, 109.5, 180.0, UNSUREDROPFLAG|ROTATEFLAG);
   theRes->addPlan(1, "H",    "1H2*", " C2*", " C3*", " C1*", " O2*", 1.1,   0.0,   0.0, O2PRIMEFLAG);
   theRes->addPlan(2, "H",    "2H2*", " C2*", " C3*", " C1*", "",     1.1, 126.5,   0.0, NOO2PRIMEFLAG);
   theRes->addPlan(2, "H",    "1H2*", " C2*", " C3*", " C1*", "",     1.1,-126.5,   0.0, NOO2PRIMEFLAG);
   theRes->addPlan(1, "H",    " H3*", " C3*", " C4*", " C2*", " O3*", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H4*", " C4*", " C5*", " C3*", " O4*", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(2, "H",    "2H5*", " C5*", " C4*", " O5*", "",     1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H",    "1H5*", " C5*", " C4*", " O5*", "",     1.1,-126.5,   0.0, 0);
   theRes->addPlan(3, "Hpol", " H5T", " O5'", " C5'", " C4'", "",     1.0, 109.5, 180.0, UNSUREDROPFLAG|ROTATEFLAG|IFNOPO4);
   theRes->addPlan(3, "Hpol", " H3T", " O3'", " C3'", " C4'", "",     1.0, 109.5, 180.0, UNSUREDROPFLAG|ROTATEFLAG|IFNOPO4);
   theRes->addPlan(3, "Hpol", " H2'", " O2'", " C2'", " C3'", "",     1.0, 109.5, 180.0, UNSUREDROPFLAG|ROTATEFLAG);
   theRes->addPlan(1, "H",    "H2''", " C2'", " C3'", " C1'", " O2'", 1.1,   0.0,   0.0, O2PRIMEFLAG);
   theRes->addPlan(2, "H",    " H2'", " C2'", " C3'", " C1'", "",     1.1, 126.5,   0.0, NOO2PRIMEFLAG);
   theRes->addPlan(2, "H",    "H2''", " C2'", " C3'", " C1'", "",     1.1,-126.5,   0.0, NOO2PRIMEFLAG);
   theRes->addPlan(1, "H",    " H3'", " C3'", " C4'", " C2'", " O3'", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H4'", " C4'", " C5'", " C3'", " O4'", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(2, "H",    " H5'", " C5'", " C4'", " O5'", "",     1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H",    "H5''", " C5'", " C4'", " O5'", "",     1.1,-126.5,   0.0, 0);
//--------------------------------------------------------------------------

   theRes = new StdResH("LYS", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "Hpol", " HZ3", " NZ", " CE", " CD", "", 1.0, 109.5,  60.0, XPLORNAME|ROTATEONDEMAND|NH3FLAG);
   theRes->addPlan(3, "Hpol", " HZ2", " NZ", " CE", " CD", "", 1.0, 109.5, -60.0, XPLORNAME|ROTATEONDEMAND|NH3FLAG);
   theRes->addPlan(3, "Hpol", " HZ1", " NZ", " CE", " CD", "", 1.0, 109.5, 180.0, XPLORNAME|ROTATEONDEMAND|NH3FLAG);
   theRes->addPlan(3, "Hpol", "3HZ",  " NZ", " CE", " CD", "", 1.0, 109.5,  60.0, NOTXPLORNAME|ROTATEONDEMAND|NH3FLAG);
   theRes->addPlan(3, "Hpol", "2HZ",  " NZ", " CE", " CD", "", 1.0, 109.5, -60.0, NOTXPLORNAME|ROTATEONDEMAND|NH3FLAG);
   theRes->addPlan(3, "Hpol", "1HZ",  " NZ", " CE", " CD", "", 1.0, 109.5, 180.0, NOTXPLORNAME|ROTATEONDEMAND|NH3FLAG);
   theRes->addPlan(2, "H", "2HE",  " CE", " CD", " NZ", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HE",  " CE", " CD", " NZ", "", 1.1,-126.5,   0.0, 0);
   theRes->addPlan(2, "H", "2HD",  " CD", " CG", " CE", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HD",  " CD", " CG", " CE", "", 1.1,-126.5,   0.0, 0);
   theRes->addPlan(2, "H", "2HG",  " CG", " CB", " CD", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HG",  " CG", " CB", " CD", "", 1.1,-126.5,   0.0, 0);
   theRes->addPlan(2, "H", "2HB",  " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB",  " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("GLY", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(2, "H", "2HA", " CA", " N", " C", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HA", " CA", " N", " C", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("GLU", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(2, "H", "2HG", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HG", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, 0);
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("THR", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "H", "3HG2", " CG2", " CB", " CA", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HG2", " CG2", " CB", " CA", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HG2", " CG2", " CB", " CA", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(3, "Hpol", " HG1", " OG1", " CB", " CA", "", 1.0, 109.5, 180.0, UNSUREDROPFLAG|ROTATEFLAG);
   theRes->addPlan(1, "H", " HB", " CB", " CA", " OG1", " CG2", 1.1,   0.0,   0.0, 0);

   theRes = new StdResH("ALA", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "H", "3HB", " CB", " CA", " N", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HB", " CB", " CA", " N", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HB", " CB", " CA", " N", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);

   theRes = new StdResH("PHE", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(4, "Har", " HZ",  " CZ",  " CE1", " CE2", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " HE2", " CE2", " CZ",  " CD2", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " HE1", " CE1", " CD1", " CZ",  "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " HD2", " CD2", " CE2", " CG",  "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " HD1", " CD1", " CG",  " CE1", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(2, "H", "2HB",  " CB",  " CA",  " CG",  "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB",  " CB",  " CA",  " CG",  "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("ARG", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "Hpol", "HH22", " NH2", " CZ", " NE", "", 1.0, 120.0, 180.0, XPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "HH21", " NH2", " CZ", " NE", "", 1.0, 120.0,   0.0, XPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "HH12", " NH1", " CZ", " NE", "", 1.0, 120.0, 180.0, XPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "HH11", " NH1", " CZ", " NE", "", 1.0, 120.0,   0.0, XPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "2HH2", " NH2", " CZ", " NE", "", 1.0, 120.0, 180.0, NOTXPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "1HH2", " NH2", " CZ", " NE", "", 1.0, 120.0,   0.0, NOTXPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "2HH1", " NH1", " CZ", " NE", "", 1.0, 120.0, 180.0, NOTXPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "1HH1", " NH1", " CZ", " NE", "", 1.0, 120.0,   0.0, NOTXPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(4, "Hpol", " HE", " NE", " CD", " CZ", "", 1.0,   0.0,   0.0, 0);
   theRes->addPlan(2, "H", "2HD", " CD", " CG", " NE", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HD", " CD", " CG", " NE", "", 1.1,-126.5,   0.0, 0);
   theRes->addPlan(2, "H", "2HG", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HG", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, 0);
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("HIS", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(4, "Ha+p", " HE2", " NE2", " CE1", " CD2", "", 1.0,   0.0,   0.0, XTRAFLAG|BONDBUMPFLAG);
   theRes->addPlan(4, "Har", " HE1", " CE1", " ND1", " NE2", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " HD2", " CD2", " NE2", " CG", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Ha+p", " HD1", " ND1", " CG", " CE1", "", 1.0,   0.0,   0.0, XTRAFLAG|BONDBUMPFLAG);
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("MET", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "H", "3HE", " CE", " SD", " CG", "", 1.1, 109.5,  60.0, ROTATEFLAG);
   theRes->addPlan(3, "H", "2HE", " CE", " SD", " CG", "", 1.1, 109.5, -60.0, ROTATEFLAG);
   theRes->addPlan(3, "H", "1HE", " CE", " SD", " CG", "", 1.1, 109.5, 180.0, ROTATEFLAG);
   theRes->addPlan(2, "H", "2HG", " CG", " CB", " SD", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HG", " CG", " CB", " SD", "", 1.1,-126.5,   0.0, 0);
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("ASP", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("SER", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "Hpol", " HG", " OG", " CB", " CA", "", 1.0, 109.5, 180.0, UNSUREDROPFLAG|ROTATEFLAG);
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " OG", "", 1.1,-126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " OG", "", 1.1, 126.5,   0.0, 0);

   theRes = new StdResH("ASN", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "Hpol", "HD22", " ND2", " CG", " OD1", "", 1.0, 120.0, 180.0, XPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "HD21", " ND2", " CG", " OD1", "", 1.0, 120.0,   0.0, XPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "2HD2", " ND2", " CG", " OD1", "", 1.0, 120.0, 180.0, NOTXPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "1HD2", " ND2", " CG", " OD1", "", 1.0, 120.0,   0.0, NOTXPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("TYR", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "Hpol", " HH", " OH", " CZ", " CE1", "", 1.0, 109.5, 180.0, UNSUREDROPFLAG|ROTATEFLAG);
   theRes->addPlan(4, "Har", " HE2", " CE2", " CZ", " CD2", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " HE1", " CE1", " CD1", " CZ", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " HD2", " CD2", " CE2", " CG", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " HD1", " CD1", " CG", " CE1", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("CYS", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "Hpol", " HG", " SG", " CB", " CA", "", 1.3, 109.5,  180.0,UNSUREDROPFLAG|ROTATEFLAG);
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " SG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " SG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("GLN", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "Hpol", "HE22", " NE2", " CD", " OE1", "", 1.0, 120.0, 180.0, XPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "HE21", " NE2", " CD", " OE1", "", 1.0, 120.0,   0.0, XPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "2HE2", " NE2", " CD", " OE1", "", 1.0, 120.0, 180.0, NOTXPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "1HE2", " NE2", " CD", " OE1", "", 1.0, 120.0,   0.0, NOTXPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(2, "H", "2HG", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HG", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, 0);
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("LEU", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "H", "3HD2", " CD2", " CG", " CB", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HD2", " CD2", " CG", " CB", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HD2", " CD2", " CG", " CB", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "3HD1", " CD1", " CG", " CB", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HD1", " CD1", " CG", " CB", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HD1", " CD1", " CG", " CB", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(1, "H", " HG", " CG", " CB", " CD1", " CD2", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("PRO", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(2, "H", "2HD", " CD", " CG", " N", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HD", " CD", " CG", " N", "", 1.1,-126.5,   0.0, 0);
   theRes->addPlan(2, "H", "2HG", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HG", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, 0);
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("VAL", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "H", "3HG2", " CG2", " CB", " CA", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HG2", " CG2", " CB", " CA", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HG2", " CG2", " CB", " CA", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "3HG1", " CG1", " CB", " CA", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HG1", " CG1", " CB", " CA", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HG1", " CG1", " CB", " CA", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(1, "H", " HB", " CB", " CA", " CG1", " CG2", 1.1,   0.0,   0.0, 0);

   theRes = new StdResH("ILE", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "H", "3HD1", " CD1", " CG1", " CB", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HD1", " CD1", " CG1", " CB", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HD1", " CD1", " CG1", " CB", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "3HG2", " CG2", " CB", " CA", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HG2", " CG2", " CB", " CA", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HG2", " CG2", " CB", " CA", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(2, "H", "2HG1", " CG1", " CB", " CD1", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HG1", " CG1", " CB", " CD1", "", 1.1,-126.5,   0.0, 0);
   theRes->addPlan(1, "H", " HB", " CB", " CA", " CG1", " CG2", 1.1,   0.0,   0.0, 0);

   theRes = new StdResH("TRP", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(4, "Har", " HH2", " CH2", " CZ3", " CZ2", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " HZ3", " CZ3", " CE3", " CH2", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " HZ2", " CZ2", " CH2", " CE2", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " HE3", " CE3", " CD2", " CZ3", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Ha+p", " HE1", " NE1", " CE2", " CD1", "", 1.0,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " HD1", " CD1", " NE1", " CG", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

//--------------------------------------------------------------------------

   theRes = new StdResH("  U", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(4, "Har", " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " H5",  " C5",  " C4",  " C6",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Ha+p"," H3",  " N3",  " C4",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(1, "H",   " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",   " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, 0);

#ifdef LEFT_JUSTIFY_NUC_RES_OK
   theRes = new StdResH("U", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(4, "Har", " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " H5",  " C5",  " C4",  " C6",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Ha+p"," H3",  " N3",  " C4",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(1, "H",   " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",   " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, 0);
#endif

   theRes = new StdResH("URA", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(4, "Har", " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " H5",  " C5",  " C4",  " C6",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Ha+p"," H3",  " N3",  " C4",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(1, "H",   " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",   " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, 0);

//--------------------------------------------------------------------------

   // note C5M is an alternative name for C5A

   theRes = new StdResH("  T", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(4, "Har",  " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(3, "H",    "3H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H",    "2H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H",    "1H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H53", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H52", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H51", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H",    "3H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H",    "2H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H",    "1H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H53", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H52", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H51", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(4, "Ha+p", " H3",  " N3",  " C4",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(1, "H",    " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, 0);

#ifdef LEFT_JUSTIFY_NUC_RES_OK
   theRes = new StdResH("T", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(4, "Har",  " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(3, "H",    "3H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H",    "2H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H",    "1H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H53", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H52", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H51", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H",    "3H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H",    "2H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H",    "1H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H53", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H52", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H51", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(4, "Ha+p", " H3",  " N3",  " C4",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(1, "H",    " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, 0);
#endif

   theRes = new StdResH("THY", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(4, "Har",  " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(3, "H",    "3H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H",    "2H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H",    "1H5M", " C5M", " C5",  " C4",  "",    1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H53", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H52", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H51", " C5M", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H",    "3H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H",    "2H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H",    "1H5A", " C5A", " C5",  " C4",  "",    1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H53", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H52", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(7, "H",    " H51", " C5A", "",     "",     "",    1.1,   0.0,   0.0, ROTATEONDEMAND);
   theRes->addPlan(4, "Ha+p", " H3",  " N3",  " C4",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(1, "H",    " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, 0);

//--------------------------------------------------------------------------

   theRes = new StdResH("  A", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(4, "Har",  " H2",  " C2",  " N1",  " N3",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(3, "Hpol", "2H6",  " N6",  " C6",  " C5",  "",    1.0,-120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "1H6",  " N6",  " C6",  " C5",  "",    1.0, 120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H62", " N6",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H61", " N6",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(4, "Har",  " H8",  " C8",  " N7",  " N9",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1*", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, 0);

#ifdef LEFT_JUSTIFY_NUC_RES_OK
   theRes = new StdResH("A", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(4, "Har",  " H2",  " C2",  " N1",  " N3",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(3, "Hpol", "2H6",  " N6",  " C6",  " C5",  "",    1.0,-120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "1H6",  " N6",  " C6",  " C5",  "",    1.0, 120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H62", " N6",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H61", " N6",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(4, "Har",  " H8",  " C8",  " N7",  " N9",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1*", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, 0);
#endif

   theRes = new StdResH("ADE", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(4, "Har",  " H2",  " C2",  " N1",  " N3",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(3, "Hpol", "2H6",  " N6",  " C6",  " C5",  "",    1.0,-120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "1H6",  " N6",  " C6",  " C5",  "",    1.0, 120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H62", " N6",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H61", " N6",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(4, "Har",  " H8",  " C8",  " N7",  " N9",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1*", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, 0);

//--------------------------------------------------------------------------

   theRes = new StdResH("  C", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(4, "Har",  " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har",  " H5",  " C5",  " C4",  " C6",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(3, "Hpol", "2H4",  " N4",  " C4",  " N3",  "",    1.0,-120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "1H4",  " N4",  " C4",  " N3",  "",    1.0, 120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H42", " N4",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H41", " N4",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(1, "H",    " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, 0);

#ifdef LEFT_JUSTIFY_NUC_RES_OK
   theRes = new StdResH("C", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(4, "Har",  " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har",  " H5",  " C5",  " C4",  " C6",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(3, "Hpol", "2H4",  " N4",  " C4",  " N3",  "",    1.0,-120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "1H4",  " N4",  " C4",  " N3",  "",    1.0, 120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H42", " N4",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H41", " N4",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(1, "H",    " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, 0);
#endif

   theRes = new StdResH("CYT", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(4, "Har",  " H6",  " C6",  " C5",  " N1",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har",  " H5",  " C5",  " C4",  " C6",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(3, "Hpol", "2H4",  " N4",  " C4",  " N3",  "",    1.0,-120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "1H4",  " N4",  " C4",  " N3",  "",    1.0, 120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H42", " N4",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H41", " N4",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(1, "H",    " H1*", " C1*", " O4*", " C2*", " N1", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1'", " C1'", " O4'", " C2'", " N1", 1.1,   0.0,   0.0, 0);

//--------------------------------------------------------------------------

   theRes = new StdResH("  G", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "Hpol", "2H2",  " N2",  " C2",  " N1",  "",    1.0,-120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "1H2",  " N2",  " C2",  " N1",  "",    1.0, 120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H22", " N2",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H21", " N2",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(4, "Ha+p", " H1",  " N1",  " C6",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(4, "Har",  " H8",  " C8",  " N7",  " N9",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1*", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, 0);

#ifdef LEFT_JUSTIFY_NUC_RES_OK
   theRes = new StdResH("G", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "Hpol", "2H2",  " N2",  " C2",  " N1",  "",    1.0,-120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "1H2",  " N2",  " C2",  " N1",  "",    1.0, 120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H22", " N2",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H21", " N2",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(4, "Ha+p", " H1",  " N1",  " C6",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(4, "Har",  " H8",  " C8",  " N7",  " N9",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1*", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, 0);
#endif

   theRes = new StdResH("GUA", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "Hpol", "2H2",  " N2",  " C2",  " N1",  "",    1.0,-120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "1H2",  " N2",  " C2",  " N1",  "",    1.0, 120.0, 180.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H22", " N2",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(7, "Hpol", " H21", " N2",  "",     "",     "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(4, "Ha+p", " H1",  " N1",  " C6",  " C2",  "",    1.0,   0.0,   0.0, BONDBUMPFLAG);
   theRes->addPlan(4, "Har",  " H8",  " C8",  " N7",  " N9",  "",    1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1*", " C1*", " O4*", " C2*", " N9", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(1, "H",    " H1'", " C1'", " O4'", " C2'", " N9", 1.1,   0.0,   0.0, 0);

//--------------------------------------------------------------------------

   theRes = new StdResH("AIB", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "H", "3HB2", " CB2", " CA", " N", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HB2", " CB2", " CA", " N", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HB2", " CB2", " CA", " N", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "3HB1", " CB1", " CA", " N", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HB1", " CB1", " CA", " N", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HB1", " CB1", " CA", " N", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);

   theRes = new StdResH("ABU", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "H", "3HG", " CG", " CB", " CA", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HG", " CG", " CB", " CA", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HG", " CG", " CB", " CA", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("ACE", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "H", "3HH3", " CH3", " C", " O", "", 1.1, 109.5,   0.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HH3", " CH3", " C", " O", "", 1.1, 109.5,-120.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HH3", " CH3", " C", " O", "", 1.1, 109.5, 120.0, ROTATEONDEMAND);

   theRes = new StdResH("ASX", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("GLX", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(2, "H", "2HG", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HG", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, 0);
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("MSE", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "H", "3HE", " CE", "SED", " CG", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HE", " CE", "SED", " CG", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HE", " CE", "SED", " CG", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(2, "H", "2HG", " CG", " CB", "SED", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HG", " CG", " CB", "SED", "", 1.1,-126.5,   0.0, 0);
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("PCA", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(2, "H", "2HG", " CG", " CB", " CD", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HG", " CG", " CB", " CD", "", 1.1,-126.5,   0.0, 0);
   theRes->addPlan(2, "H", "2HB", " CB", " CA", " CG", "", 1.1, 126.5,   0.0, 0);
   theRes->addPlan(2, "H", "1HB", " CB", " CA", " CG", "", 1.1,-126.5,   0.0, 0);

   theRes = new StdResH("NH2", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "Hpol", " HN2", " N", "- C", "- O", "", 1.0, 120.0,   0.0, XPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", " HN1", " N", "- C", "- O", "", 1.0, 120.0, 180.0, XPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "2HN",  " N", "- C", "- O", "", 1.0, 120.0,   0.0, NOTXPLORNAME|BONDBUMPFLAG);
   theRes->addPlan(3, "Hpol", "1HN",  " N", "- C", "- O", "", 1.0, 120.0, 180.0, NOTXPLORNAME|BONDBUMPFLAG);

   theRes = new StdResH("NME", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(3, "H", "3HH3", " CH3", " N",   "- C", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HH3", " CH3", " N",   "- C", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HH3", " CH3", " N",   "- C", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(4, "Hpol", " H",   " N",   " CH3", "- C", "", 1.0,   0.0,   0.02,0);

   theRes = new StdResH("FOR", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(4, "H", " H", " C", " O", "+ N", "", 1.1,   0.0,   0.0, 0);

   theRes = new StdResH("HEM", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(2, "H", "2HAD", " CAD", " C3D", " CBD", "", 1.1,-126.5,  0.0, 0);
   theRes->addPlan(2, "H", "1HAD", " CAD", " C3D", " CBD", "", 1.1, 126.5,  0.0, 0);
   theRes->addPlan(2, "H", "2HBD", " CBD", " CAD", " CGD", "", 1.1,-126.5,  0.0, 0);
   theRes->addPlan(2, "H", "1HBD", " CBD", " CAD", " CGD", "", 1.1, 126.5,  0.0, 0);
   theRes->addPlan(2, "H", "2HAA", " CAA", " C2A", " CBA", "", 1.1,-126.5,  0.0, 0);
   theRes->addPlan(2, "H", "1HAA", " CAA", " C2A", " CBA", "", 1.1, 126.5,  0.0, 0);
   theRes->addPlan(2, "H", "2HBA", " CBA", " CAA", " CGA", "", 1.1,-126.5,  0.0, 0);
   theRes->addPlan(2, "H", "1HBA", " CBA", " CAA", " CGA", "", 1.1, 126.5,  0.0, 0);
   theRes->addPlan(3, "H", "2HBC", " CBC", " CAC", " C3C", "", 1.1, 120.0,   0.0, 0);
   theRes->addPlan(3, "H", "1HBC", " CBC", " CAC", " C3C", "", 1.1, 120.0, 180.0, 0);
   theRes->addPlan(4, "H", " HAC", " CAC", " CBC", " C3C", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(3, "H", "2HBB", " CBB", " CAB", " C3B", "", 1.1, 120.0,   0.0, 0);
   theRes->addPlan(3, "H", "1HBB", " CBB", " CAB", " C3B", "", 1.1, 120.0, 180.0, 0);
   theRes->addPlan(4, "H", " HAB", " CAB", " CBB", " C3B", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(3, "H", "3HMD", " CMD", " C2D", " C1D", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HMD", " CMD", " C2D", " C1D", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HMD", " CMD", " C2D", " C1D", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "3HMC", " CMC", " C2C", " C1C", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HMC", " CMC", " C2C", " C1C", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HMC", " CMC", " C2C", " C1C", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "3HMB", " CMB", " C2B", " C1B", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HMB", " CMB", " C2B", " C1B", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HMB", " CMB", " C2B", " C1B", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "3HMA", " CMA", " C3A", " C4A", "", 1.1, 109.5, -60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "2HMA", " CMA", " C3A", " C4A", "", 1.1, 109.5,  60.0, ROTATEONDEMAND);
   theRes->addPlan(3, "H", "1HMA", " CMA", " C3A", " C4A", "", 1.1, 109.5, 180.0, ROTATEONDEMAND);
   theRes->addPlan(4, "Har", " HHD", " CHD", " C1D", " C4C", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " HHC", " CHC", " C1C", " C4B", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " HHB", " CHB", " C1B", " C4A", "", 1.1,   0.0,   0.0, 0);
   theRes->addPlan(4, "Har", " HHA", " CHA", " C1A", " C4D", "", 1.1,   0.0,   0.0, 0);

//--------------------------------------------------------------------------

   theRes = new StdResH("HOH", "");
   _restbl.insert(std::make_pair(theRes->name(), theRes));
   theRes->addPlan(7, "Hpol", " H2", " OD2", "", "", "", 1.0,   0.0,   0.0, UNSUREDROPFLAG);
   theRes->addPlan(7, "Hpol", " H1", " OD2", "", "", "", 1.0,   0.0,   0.0, UNSUREDROPFLAG);
   theRes->addPlan(7, "Hpol", " H2", " OH2", "", "", "", 1.0,   0.0,   0.0, UNSUREDROPFLAG);
   theRes->addPlan(7, "Hpol", " H1", " OH2", "", "", "", 1.0,   0.0,   0.0, UNSUREDROPFLAG);
   theRes->addPlan(7, "Hpol", " H2", " OW",  "", "", "", 1.0,   0.0,   0.0, UNSUREDROPFLAG);
   theRes->addPlan(7, "Hpol", " H1", " OW",  "", "", "", 1.0,   0.0,   0.0, UNSUREDROPFLAG);
   theRes->addPlan(7, "Hpol", " H2", " O",   "", "", "", 1.0,   0.0,   0.0, UNSUREDROPFLAG);
   theRes->addPlan(7, "Hpol", " H1", " O",   "", "", "", 1.0,   0.0,   0.0, UNSUREDROPFLAG);
   // alternative naming of hydrogens
   theRes->addPlan(7, "Hpol", "2H", " OD2", "", "", "", 1.0,   0.0,   0.0, UNSUREDROPFLAG);
   theRes->addPlan(7, "Hpol", "1H", " OD2", "", "", "", 1.0,   0.0,   0.0, UNSUREDROPFLAG);
   theRes->addPlan(7, "Hpol", "2H", " OH2", "", "", "", 1.0,   0.0,   0.0, UNSUREDROPFLAG);
   theRes->addPlan(7, "Hpol", "1H", " OH2", "", "", "", 1.0,   0.0,   0.0, UNSUREDROPFLAG);
   theRes->addPlan(7, "Hpol", "2H", " OW",  "", "", "", 1.0,   0.0,   0.0, UNSUREDROPFLAG);
   theRes->addPlan(7, "Hpol", "1H", " OW",  "", "", "", 1.0,   0.0,   0.0, UNSUREDROPFLAG);
   theRes->addPlan(7, "Hpol", "2H", " O",   "", "", "", 1.0,   0.0,   0.0, UNSUREDROPFLAG);
   theRes->addPlan(7, "Hpol", "1H", " O",   "", "", "", 1.0,   0.0,   0.0, UNSUREDROPFLAG);
}

// initialize xtra info about the standard residues
StdResXtraInfo::StdResXtraInfo() {
   _atomAttributes.insert(std::make_pair(makeKey("GLY","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("ALA","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("VAL","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("PHE","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("PRO","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("MET","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("ILE","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("LEU","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("ASP","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("GLU","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("LYS","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("ARG","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("SER","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("THR","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("TYR","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("HIS","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("CYS","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("ASN","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("GLN","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("TRP","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("ASX","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("GLX","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("ABU","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("AIB","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("ABU","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("MSE","AminoAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("PCA","AminoAcid"), 1));

   _atomAttributes.insert(std::make_pair(makeKey("  U","NucleicAcid"), 1));
#ifdef LEFT_JUSTIFY_NUC_RES_OK
   _atomAttributes.insert(std::make_pair(makeKey("U",  "NucleicAcid"), 1));
#endif
   _atomAttributes.insert(std::make_pair(makeKey("URA","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("  T","NucleicAcid"), 1));
#ifdef LEFT_JUSTIFY_NUC_RES_OK
   _atomAttributes.insert(std::make_pair(makeKey("T",  "NucleicAcid"), 1));
#endif
   _atomAttributes.insert(std::make_pair(makeKey("THY","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("  A","NucleicAcid"), 1));
#ifdef LEFT_JUSTIFY_NUC_RES_OK
   _atomAttributes.insert(std::make_pair(makeKey("A",  "NucleicAcid"), 1));
#endif
   _atomAttributes.insert(std::make_pair(makeKey("ADE","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("  C","NucleicAcid"), 1));
#ifdef LEFT_JUSTIFY_NUC_RES_OK
   _atomAttributes.insert(std::make_pair(makeKey("C",  "NucleicAcid"), 1));
#endif
   _atomAttributes.insert(std::make_pair(makeKey("CYT","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("  G","NucleicAcid"), 1));
#ifdef LEFT_JUSTIFY_NUC_RES_OK
   _atomAttributes.insert(std::make_pair(makeKey("G",  "NucleicAcid"), 1));
#endif
   _atomAttributes.insert(std::make_pair(makeKey("GUA","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("CTP","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("CDP","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("CMP","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("GDP","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("GMP","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("ATP","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("ADP","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("AMP","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("TTP","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("TDP","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("TMP","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("UTP","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("UDP","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("UMP","NucleicAcid"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("GSP","NucleicAcid"), 1));

// (reduced VDW radii indicated for C=O carbon atoms
//  in-addition to AminoAcid mainchain COs)

   _atomAttributes.insert(std::make_pair(makeKey("ASP", " CG"), ISACOFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ASN", " CG"), ISACOFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ASX", " CG"), ISACOFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("GLU", " CD"), ISACOFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("GLN", " CD"), ISACOFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("GLX", " CD"), ISACOFLAG));

// (setting the HBACCEPTORFLAG & AROMATICFLAG for certain carbons & nitrogens)

   // HB status of HIS may depend on the protonation state
   _atomAttributes.insert(std::make_pair(makeKey("HIS", " ND1"),HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HIS", " NE2"),HBACCEPTORFLAG|AROMATICFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("  A", " N1"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  A", " N3"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  A", " N7"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  C", " N3"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  G", " N3"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  G", " N7"), HBACCEPTORFLAG|AROMATICFLAG));
#ifdef LEFT_JUSTIFY_NUC_RES_OK
   _atomAttributes.insert(std::make_pair(makeKey("A",   " N1"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("A",   " N3"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("A",   " N7"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("C",   " N3"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("G",   " N3"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("G",   " N7"), HBACCEPTORFLAG|AROMATICFLAG));
#endif
   _atomAttributes.insert(std::make_pair(makeKey("ADE", " N1"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ADE", " N3"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ADE", " N7"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("CYT", " N3"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("GUA", " N3"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("GUA", " N7"), HBACCEPTORFLAG|AROMATICFLAG));

// other HEM atom atributes lower down in file
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " N A"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " N B"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " N C"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " N D"), HBACCEPTORFLAG|AROMATICFLAG));

// always assumed to be charged...
   _atomAttributes.insert(std::make_pair(makeKey("AminoAcid",   " NT"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("AminoAcid",   "1H"  ), POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("AminoAcid",   "2H"  ), POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("AminoAcid",   "3H"  ), POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("AminoAcid",   " HT1"), POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("AminoAcid",   " HT2"), POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("AminoAcid",   " HT3"), POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("AminoAcid",   "1D"  ), POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("AminoAcid",   "2D"  ), POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("AminoAcid",   "3D"  ), POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("AminoAcid",   " DT1"), POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("AminoAcid",   " DT2"), POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("AminoAcid",   " DT3"), POSCHARGEFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("AminoAcid",   " OXT"), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("AminoAcid",   "1OXT"), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("AminoAcid",   "2OXT"), NEGCHARGEFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("NucleicAcid", " P"  ), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("NucleicAcid", " O1P"), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("NucleicAcid", " O2P"), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("NucleicAcid", " PA" ), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("NucleicAcid", " PB" ), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("NucleicAcid", " PG" ), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("NucleicAcid", " O1A"), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("NucleicAcid", " O2A"), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("NucleicAcid", " O3A"), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("NucleicAcid", " O1B"), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("NucleicAcid", " O2B"), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("NucleicAcid", " O3B"), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("NucleicAcid", " O1G"), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("NucleicAcid", " O2G"), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("NucleicAcid", " O3G"), NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("NucleicAcid", " S1G"), NEGCHARGEFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("ASP","ChargedResidue"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("GLU","ChargedResidue"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("LYS","ChargedResidue"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("ARG","ChargedResidue"), 1));
   _atomAttributes.insert(std::make_pair(makeKey("HEM","ChargedResidue"), 1));

   _atomAttributes.insert(std::make_pair(makeKey("ASP", " OD1"),  NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ASP", " OD2"),  NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("GLU", " OE1"),  NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("GLU", " OE2"),  NEGCHARGEFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("LYS", " NZ" ),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("LYS", "1HZ" ),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("LYS", "2HZ" ),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("LYS", "3HZ" ),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("LYS", " HZ1"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("LYS", " HZ2"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("LYS", " HZ3"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("LYS", "1DZ" ),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("LYS", "2DZ" ),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("LYS", "3DZ" ),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("LYS", " DZ1"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("LYS", " DZ2"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("LYS", " DZ3"),  POSCHARGEFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("ARG", " NE" ),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", " NH1"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", " NH2"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", " HE" ),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", "1HH1"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", "2HH1"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", "1HH2"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", "2HH2"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", "HH11"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", "HH12"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", "HH21"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", "HH22"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", "1DH1"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", "2DH1"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", "1DH2"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", "2DH2"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", "DH11"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", "DH12"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", "DH21"),  POSCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ARG", "DH22"),  POSCHARGEFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("HEM", " O1A"),  NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " O2A"),  NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " O1D"),  NEGCHARGEFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " O2D"),  NEGCHARGEFLAG));
#ifdef AROMATICS_ACCEPT_HBONDS
// -----------------------------------------------------------
// here we treat the aromatic Pi-bonds as hydrogen bond acceptors
// -----------------------------------------------------------
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " C1A"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " C2A"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " C3A"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " C4A"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " C1B"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " C2B"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " C3B"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " C4B"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " C1C"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " C2C"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " C3C"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " C4C"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " C1D"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " C2D"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " C3D"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HEM", " C4D"), HBACCEPTORFLAG|AROMATICFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("PHE", " CZ"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("PHE", " CE2"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("PHE", " CE1"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("PHE", " CD2"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("PHE", " CD1"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("PHE", " CG"),  HBACCEPTORFLAG|AROMATICFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("TYR", " CZ"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("TYR", " CE2"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("TYR", " CE1"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("TYR", " CD2"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("TYR", " CD1"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("TYR", " CG"),  HBACCEPTORFLAG|AROMATICFLAG));

   //_atomAttributes.insert(std::make_pair(makeKey("HIS", " CD2"), HBACCEPTORFLAG|AROMATICFLAG));
   //_atomAttributes.insert(std::make_pair(makeKey("HIS", " CE1"), HBACCEPTORFLAG|AROMATICFLAG));
   //_atomAttributes.insert(std::make_pair(makeKey("HIS", " CG"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HIS", " CD2"), AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HIS", " CE1"), AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("HIS", " CG"),  AROMATICFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("TRP", " CH2"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("TRP", " CZ3"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("TRP", " CZ2"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("TRP", " CE3"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("TRP", " CE2"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("TRP", " NE1"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("TRP", " CD2"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("TRP", " CD1"), HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("TRP", " CG"),  HBACCEPTORFLAG|AROMATICFLAG));

// -----------------------------------------------------------
// pick up the parts of the bases not included above
// *** warning *** the bases below are missing the atoms
//                 which are listed above (atoms not duplicated)
// -----------------------------------------------------------
   _atomAttributes.insert(std::make_pair(makeKey("  A", " C2"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  A", " C4"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  A", " C5"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  A", " C6"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  A", " C8"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  A", " N9"),  HBACCEPTORFLAG|AROMATICFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("  C", " N1"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  C", " C2"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  C", " C4"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  C", " C5"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  C", " C6"),  HBACCEPTORFLAG|AROMATICFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("  G", " N1"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  G", " C2"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  G", " C4"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  G", " C5"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  G", " C6"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  G", " C8"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  G", " N9"),  HBACCEPTORFLAG|AROMATICFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("  T", " N1"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  T", " C2"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  T", " N3"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  T", " C4"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  T", " C5"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  T", " C6"),  HBACCEPTORFLAG|AROMATICFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("  U", " N1"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  U", " C2"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  U", " N3"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  U", " C4"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  U", " C5"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("  U", " C6"),  HBACCEPTORFLAG|AROMATICFLAG));

#ifdef LEFT_JUSTIFY_NUC_RES_OK
   _atomAttributes.insert(std::make_pair(makeKey("A", " C2"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("A", " C4"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("A", " C5"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("A", " C6"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("A", " C8"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("A", " N9"),  HBACCEPTORFLAG|AROMATICFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("C", " N1"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("C", " C2"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("C", " C4"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("C", " C5"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("C", " C6"),  HBACCEPTORFLAG|AROMATICFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("G", " N1"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("G", " C2"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("G", " C4"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("G", " C5"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("G", " C6"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("G", " C8"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("G", " N9"),  HBACCEPTORFLAG|AROMATICFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("T", " N1"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("T", " C2"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("T", " N3"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("T", " C4"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("T", " C5"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("T", " C6"),  HBACCEPTORFLAG|AROMATICFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("U", " N1"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("U", " C2"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("U", " N3"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("U", " C4"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("U", " C5"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("U", " C6"),  HBACCEPTORFLAG|AROMATICFLAG));
#endif

   _atomAttributes.insert(std::make_pair(makeKey("ADE", " C2"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ADE", " C4"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ADE", " C5"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ADE", " C6"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ADE", " C8"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("ADE", " N9"),  HBACCEPTORFLAG|AROMATICFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("CYT", " N1"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("CYT", " C2"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("CYT", " C4"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("CYT", " C5"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("CYT", " C6"),  HBACCEPTORFLAG|AROMATICFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("GUA", " N1"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("GUA", " C2"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("GUA", " C4"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("GUA", " C5"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("GUA", " C6"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("GUA", " C8"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("GUA", " N9"),  HBACCEPTORFLAG|AROMATICFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("THY", " N1"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("THY", " C2"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("THY", " N3"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("THY", " C4"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("THY", " C5"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("THY", " C6"),  HBACCEPTORFLAG|AROMATICFLAG));

   _atomAttributes.insert(std::make_pair(makeKey("URA", " N1"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("URA", " C2"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("URA", " N3"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("URA", " C4"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("URA", " C5"),  HBACCEPTORFLAG|AROMATICFLAG));
   _atomAttributes.insert(std::make_pair(makeKey("URA", " C6"),  HBACCEPTORFLAG|AROMATICFLAG));
#endif
}
