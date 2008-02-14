// name: ElementInfo.C
// author: J. Michael Word     date written: 6/12/97
// purpose: determine element from PDB atom name

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
#pragma warning(disable:4305) 
#endif

#include <iostream>
using std::endl;
using std::cerr;

#ifdef OLD_STD_HDRS
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#else
#include <cstring>
#include <cctype>
#include <cstdio>
using std::strstr;
using std::strcpy;
using std::toupper;
using std::sprintf;
#endif

#include <stdexcept>

#include "ElementInfo.h"
#include "StdResH.h"
//#include "PDBrec.h"

const StandardElementTable * ElementInfo::TheStdElemTbl = NULL;

// build TheStdElemTbl only once, when first used

const StandardElementTable& ElementInfo::StdElemTbl() {
   if (TheStdElemTbl == NULL) { TheStdElemTbl = new StandardElementTable; }
   return *TheStdElemTbl;
}

int basicChargeState(const char* atomname, const char* resname,
   int posFlag, int negFlag, ElementInfo *e) {
   int chargeState = 0;

   // histadine charge state not explicitly set here

   if (e && e->hasProp(METALIC_ATOM)) {
      chargeState |= posFlag;
   }
   else if (atomname && resname) {
      if (StdResH::ResXtraInfo().match(resname, "ChargedResidue")) {
	 if (StdResH::ResXtraInfo().atomHasAttrib(resname, atomname,
	    NEGCHARGEFLAG)) { chargeState |= negFlag; }
	 if (StdResH::ResXtraInfo().atomHasAttrib(resname, atomname,
	    POSCHARGEFLAG)) { chargeState |= posFlag; }
      }

      if (chargeState == 0) {
	 if (StdResH::ResXtraInfo().match(resname, "AminoAcid")) {
	    if (StdResH::ResXtraInfo().atomHasAttrib("AminoAcid", atomname,
	       NEGCHARGEFLAG)) { chargeState |= negFlag; }
	    if (StdResH::ResXtraInfo().atomHasAttrib("AminoAcid", atomname,
	       POSCHARGEFLAG)) { chargeState |= posFlag; }
	 }
	 else if (StdResH::ResXtraInfo().match(resname, "NucleicAcid")) {
	    if (StdResH::ResXtraInfo().atomHasAttrib("NucleicAcid", atomname,
	       NEGCHARGEFLAG)) { chargeState |= negFlag; }
	    if (StdResH::ResXtraInfo().atomHasAttrib("NucleicAcid", atomname,
	       POSCHARGEFLAG)) { chargeState |= posFlag; }
	 }
      }
   }
#ifdef DEBUGCHARGESTATE
if (chargeState != 0) {
   cerr << "DEBUG: ChargeState(" << resname << "," << atomname << ") "
   << (((chargeState & posFlag)!=0) ? "+" : "")
   << (((chargeState & negFlag)!=0) ? "-" : "")
   << endl;
}
#endif

   return chargeState;
}

bool fixupAmbigAtomName(char* atomname, const char* resname, char* segID) {
   bool rc = FALSE;
   if (atomname && resname
      && ::strstr(":ASN:asn:GLN:gln:", resname)
      && ::strstr(": AD1: ad1: AD2: ad2: AE1: ae1: AE2: ae2:", atomname)) {
      if (atomname[1] == 'A') {
	 if (atomname[3] == '1') {
	    atomname[1] = 'O';
	    ::strcpy(segID, "A>O");
	    rc = TRUE;
	 }
	 else if (atomname[3] == '2') {
	    atomname[1] = 'N';
	    ::strcpy(segID, "A>N");
	    rc = TRUE;
	 }
      }
      else if (atomname[1] == 'a') {
	 if (atomname[3] == '1') {
	    atomname[1] = 'o';
	    ::strcpy(segID, "a>o");
	    rc = TRUE;
	 }
	 else if (atomname[3] == '2') {
	    atomname[1] = 'n';
	    ::strcpy(segID, "a>n");
	    rc = TRUE;
	 }
      }
   }
   return rc;
}

bool fixAtomName(const char* atomname, const char* resname,int position) {
   char resn[6];
   ::sprintf(resn, ":%-3.3s:", resname);
   char buf[5] = "    ";
   int i;
   for (i = 0; i < 4; i++) { // uppercase the input
      if (atomname[i] == '\0') { break; }
#ifdef CHARFUNCMACROS
      buf[i] = toupper(atomname[i]);
#else
      buf[i] = ::toupper(atomname[i]);
#endif
   }
   buf[i] = '\0';
        switch(buf[position]) {
           case 'E': return strstr(HE_RESNAMES, resn) != NULL;
           case 'F': return strstr(HF_RESNAMES, resn) != NULL; 
           case 'G': return strstr(HG_RESNAMES, resn) != NULL;
           case 'O': return strstr(HO_RESNAMES, resn) != NULL;
           case 'S': return strstr(HS_RESNAMES, resn) != NULL;
        } 
   throw std::runtime_error("Internal Error: fixAtomName() " __FILE__);
}

ElementInfo::ElementInfoRep::ElementInfoRep(int atno,
		  const char* name, const char* fullName,
		  float eRad, float iRad, float covRad,
		  const char* color, int  flags) {

   _count = 1;

   _atno = atno;
   ::strcpy(_name = new char[::strlen(name)+1], name);
   ::strcpy(_fullName = new char[::strlen(fullName)+1], fullName);
   _eRad = eRad;
   _iRad = iRad;
   _covRad = covRad;
   ::strcpy(_color = new char[::strlen(color)+1], color);
   _flags = flags;
}

ElementInfo::ElementInfoRep::~ElementInfoRep() {
   delete[] _name;
   delete[] _fullName;
   delete[] _color;
}

void StandardElementTable::LayoutTable() {

   _explMaxRad = _implMaxRad = _covMaxRad = 0.0;

// For non-metals, explicit VDW radii from
// Gavezzotti, J. Am. Chem. Soc. (1983) 105, 5220-5225.
// or, if unavailable,
// Bondi, J. Phys. Chem. (1964), V68, N3, 441-451.
// Covalent and ionic radii from
// Advanced Inorganic Chemistry, Cotton & Wilkinson, 1962, p93.

//    atno                          explRad implRad covRad mageColors flags
insert( 0, "?",  "unknown",            1.00, 0.00, 0.00, "magenta", 0);
insert( 0, "ignore", "ignore",         0.00, 0.00, 0.00, "magenta", IGNORE);

insert( 1, "H",  "hydrogen",           1.17, 0.00, 0.30, "grey",   0);
insert( 1, "Har","hydrogen(aromatic)", 1.00, 0.00, 0.30, "grey",   0);
insert( 1, "Hpol","hydrogen(polar)",   1.00, 0.00, 0.30, "grey",   DONOR_ATOM);
insert( 1, "Ha+p",
         "hydrogen(aromatic&polar)",   1.00, 0.00, 0.30, "grey",   DONOR_ATOM);
	 
insert( 1, "HOd",
         "hydrogen(omnidirectional)",  1.00, 0.00, 0.30, "grey",   DONOR_ATOM|HB_ONLY_DUMMY);

insert( 6, "C",  "carbon",             1.75, 1.90, 0.77, "white",  0);
insert( 6, "Car","carbon(aromatic)",   1.75, 1.90, 0.77, "white",  ACCEPTOR_ATOM);
insert( 6, "C=O","carbon(carbonyl)",   1.65, 1.80, 0.77, "white",  0); //** seems to help
insert( 7, "N",  "nitrogen",           1.55, 1.70, 0.70, "sky",    0);
insert( 7, "Nacc","nitrogen(acceptor)",1.55, 1.70, 0.70, "sky",    ACCEPTOR_ATOM);
insert( 8, "O",  "oxygen",             1.40, 1.50, 0.66, "red",    ACCEPTOR_ATOM);
insert(15, "P",  "phosphorus",         1.80, 1.80, 1.10, "pink",   0);
insert(16, "S",  "sulfur",             1.80, 1.90, 1.04, "yellow", ACCEPTOR_ATOM);
insert(33, "As", "arsnic",             2.00, 2.10, 1.21, "grey",   0);
insert(34, "Se", "selenium",           1.90, 2.00, 1.17, "green",  0);

insert( 9, "F",  "fluorine",           1.30, 1.30, 0.58, "green",  ACCEPTOR_ATOM);
insert(17, "Cl", "chlorine",           1.77, 1.77, 0.99, "green",  ACCEPTOR_ATOM);
insert(35, "Br", "bromine",            1.95, 1.95, 1.14, "brown",  ACCEPTOR_ATOM);
insert(53, "I",  "iodine",             2.10, 2.10, 1.33, "brown",  ACCEPTOR_ATOM);

 // for most common metals we use Pauling's ionic radii
 // "covalent radii" = ionic + 0.74 (i.e., oxygenVDW(1.4) - oxygenCov(0.66))
 // because the ionic radii are usually calculated from Oxygen-Metal distance
insert( 3, "Li", "lithium",            0.60, 0.60, 1.34, "grey", METALIC_ATOM);
insert(11, "Na", "sodium",             0.95, 0.95, 1.69, "grey", METALIC_ATOM);
insert(13, "Al", "aluminum",           0.50, 0.50, 1.24, "grey", METALIC_ATOM);
insert(19, "K",  "potassium",          1.33, 1.33, 2.07, "grey", METALIC_ATOM);
insert(12, "Mg", "magnesium",          0.65, 0.65, 1.39, "grey", METALIC_ATOM);
insert(20, "Ca", "calcium",            0.99, 0.99, 1.73, "grey", METALIC_ATOM);
insert(25, "Mn", "manganese",          0.80, 0.80, 1.54, "grey", METALIC_ATOM);
insert(26, "Fe", "iron",               0.74, 0.74, 1.48, "grey", METALIC_ATOM);
insert(27, "Co", "cobolt",             0.70, 0.70, 1.44, "blue", METALIC_ATOM);
insert(28, "Ni", "nickel",             0.66, 0.66, 1.40, "grey", METALIC_ATOM);
insert(29, "Cu", "copper",             0.72, 0.72, 1.46,"orange",METALIC_ATOM);
insert(30, "Zn", "zinc",               0.71, 0.71, 1.45, "grey", METALIC_ATOM);
insert(37, "Rb", "rubidium",           1.48, 1.48, 2.22, "grey", METALIC_ATOM);
insert(38, "Sr", "strontium",          1.10, 1.10, 1.84, "grey", METALIC_ATOM);
insert(42, "Mo", "molybdenum",         0.93, 0.93, 1.67, "grey", METALIC_ATOM);
insert(47, "Ag", "silver",             1.26, 1.26, 2.00, "white",METALIC_ATOM);
insert(48, "Cd", "cadmium",            0.91, 0.91, 1.65, "grey", METALIC_ATOM);
insert(49, "In", "indium",             0.81, 0.81, 1.55, "grey", METALIC_ATOM);
insert(55, "Cs", "cesium",             1.69, 1.69, 2.43, "grey", METALIC_ATOM);
insert(56, "Ba", "barium",             1.29, 1.29, 2.03, "grey", METALIC_ATOM);
insert(79, "Au", "gold",               1.10, 1.10, 1.84, "gold", METALIC_ATOM);
insert(80, "Hg", "mercury",            1.00, 1.00, 1.74, "grey", METALIC_ATOM);
insert(81, "Tl", "thallium",           1.44, 1.44, 2.18, "grey", METALIC_ATOM);
insert(82, "Pb", "lead",               0.84, 0.84, 1.58, "grey", METALIC_ATOM);

// for other metals we use Shannon's ionic radii
// Acta Crystallogr. (1975) A32, pg751.
insert(23, "V",  "vanadium",           0.79, 0.79, 1.53, "grey", METALIC_ATOM);
insert(24, "Cr", "chromium",           0.73, 0.73, 1.47, "grey", METALIC_ATOM);
insert(52, "Te", "tellurium",          0.97, 0.97, 1.71, "grey", METALIC_ATOM);
insert(62, "Sm", "samarium",           1.08, 1.08, 1.82, "grey", METALIC_ATOM);
insert(64, "Gd", "gadolinium",         1.05, 1.05, 1.79, "grey", METALIC_ATOM);
insert(70, "Yb", "ytterbium",          1.14, 1.14, 1.88, "grey", METALIC_ATOM);
insert(74, "W",  "tungsten",           0.66, 0.66, 1.40, "grey", METALIC_ATOM);
insert(78, "Pt", "platinum",           0.63, 0.63, 1.37, "grey", METALIC_ATOM);
insert(92, "U",  "unanium",            1.03, 1.03, 1.77, "grey", METALIC_ATOM);

// Cotton & Wilkinson and also-
// L.E. Sutton (ed.) in Table of interatomic distances and configuration in molecules
// and ions, Supplement 1956-1959, Special publication No. 18, Chemical Society,
// London, UK, 1965 (as listed in web-elements by Mark Winter)
//                   http://www.shef.ac.uk/chemistry/web-elements

insert( 2, "He",  "helium",            1.60, 1.60, 0.00, "sky",             0);
insert( 4, "Be",  "beryllium",         0.31, 0.31, 0.90, "grey", METALIC_ATOM);
insert( 5, "B",   "boron",             0.20, 0.20, 0.86, "grey",            0);
insert(10, "Ne",  "neon",              1.60, 1.60, 0.00, "pink",            0);
insert(14, "Si",  "silicon",           2.10, 2.10, 1.17, "grey", METALIC_ATOM);
insert(18, "Ar",  "argon",             1.89, 1.89, 0.00, "orange",          0);
insert(21, "Sc",  "scandium",          0.68, 0.68, 0.44, "grey", METALIC_ATOM);
insert(22, "Ti",  "titanium",          0.75, 0.75, 1.49, "grey", METALIC_ATOM);
insert(31, "Ga",  "gallium",           0.53, 0.53, 1.27, "grey", METALIC_ATOM);
insert(32, "Ge",  "germanium",         0.60, 0.60, 1.34, "grey", METALIC_ATOM);
insert(36, "Kr",  "krypton",           2.01, 2.01, 1.15, "greentint",       0);
insert(39, "Y",   "yttrium",           0.90, 0.90, 1.64, "grey", METALIC_ATOM);
insert(40, "Zr",  "zirconium",         0.77, 0.77, 1.51, "grey", METALIC_ATOM);
insert(50, "Sn",  "tin",               0.71, 0.71, 1.45, "grey", METALIC_ATOM);
insert(51, "Sb",  "antimony",          2.20, 2.20, 1.41, "grey", METALIC_ATOM);
insert(54, "Xe",  "xenon",             2.18, 2.18, 1.28, "magenta",         0);
insert(57, "La",  "lanthanum",         1.03, 1.03, 1.77, "grey", METALIC_ATOM);
insert(58, "Ce",  "cerium",            0.87, 0.87, 1.61, "grey", METALIC_ATOM);
insert(87, "Fr",  "francium",          1.94, 1.94, 2.68, "grey", METALIC_ATOM);
insert(88, "Ra",  "radium",            1.62, 1.62, 2.36, "grey", METALIC_ATOM);
insert(90, "Th",  "thorium",           1.08, 1.08, 1.82, "grey", METALIC_ATOM);

// finally, we have a set of elements where the radii are unknown
// so we use estimates and extrapolations based on web-elements data
insert(41, "Nb",  "niobium",           0.86, 0.86, 1.40, "grey", METALIC_ATOM);
insert(43, "Tc",  "technetium",        0.71, 0.71, 1.25, "grey", METALIC_ATOM);
insert(44, "Ru",  "ruthenium",         0.82, 0.82, 1.36, "grey", METALIC_ATOM);
insert(45, "Rh",  "rhodium",           0.76, 1.76, 1.30, "grey", METALIC_ATOM);
insert(46, "Pd",  "palladium",         1.05, 1.05, 1.59, "grey", METALIC_ATOM);
insert(59, "Pr",  "praseodymium",      1.11, 1.11, 1.65, "grey", METALIC_ATOM);
insert(60, "Nd",  "neodymium",         1.10, 1.10, 1.64, "grey", METALIC_ATOM);
insert(61, "Pm",  "promethium",        1.15, 1.15, 1.89, "grey", METALIC_ATOM);
insert(63, "Eu",  "europium",          1.31, 1.31, 1.85, "grey", METALIC_ATOM);
insert(65, "Tb",  "terbium",           1.05, 1.05, 1.59, "grey", METALIC_ATOM);
insert(66, "Dy",  "dysprosium",        1.05, 1.05, 1.59, "grey", METALIC_ATOM);
insert(67, "Ho",  "holmium",           1.04, 1.04, 1.58, "grey", METALIC_ATOM);
insert(68, "Er",  "erbium",            1.03, 1.03, 1.57, "grey", METALIC_ATOM);
insert(69, "Tm",  "thulium",           1.02, 1.02, 1.56, "grey", METALIC_ATOM);
insert(71, "Lu",  "lutetium",          1.02, 1.02, 1.56, "grey", METALIC_ATOM);
insert(72, "Hf",  "hafnium",           0.85, 0.85, 1.46, "grey", METALIC_ATOM);
insert(73, "Ta",  "tantalum",          0.86, 0.86, 1.40, "grey", METALIC_ATOM);
insert(75, "Re",  "rhenium",           0.77, 0.77, 1.31, "grey", METALIC_ATOM);
insert(76, "Os",  "osmium",            0.78, 0.78, 1.32, "grey", METALIC_ATOM);
insert(77, "Ir",  "iridium",           0.80, 0.80, 1.34, "grey", METALIC_ATOM);
insert(83, "Bi",  "bismuth",           1.17, 1.17, 1.71, "grey", METALIC_ATOM);
insert(84, "Po",  "polonium",          0.99, 0.99, 1.53, "grey", METALIC_ATOM);
insert(85, "At",  "astatine",          0.91, 0.91, 1.45, "grey", METALIC_ATOM);
insert(86, "Rn",  "radon",             2.50, 2.50, 1.25, "pinktint", 0);
insert(89, "Ac",  "actinium",          1.30, 1.30, 2.00, "grey", METALIC_ATOM);
insert(91, "Pa",  "protoactinium",     1.10, 1.10, 1.85, "grey", METALIC_ATOM);
insert(93, "Np",  "neptunium",         1.00, 1.00, 1.72, "grey", METALIC_ATOM);
insert(94, "Pu",  "plutonium",         1.00, 1.00, 1.67, "grey", METALIC_ATOM);
insert(95, "Am",  "americium",         1.00, 1.00, 1.63, "grey", METALIC_ATOM);
insert(96, "Cm",  "curium",            1.00, 1.00, 1.60, "grey", METALIC_ATOM);
insert(97, "Bk",  "berkelium",         1.00, 1.00, 1.58, "grey", METALIC_ATOM);
insert(98, "Cf",  "californium",       1.00, 1.00, 1.57, "grey", METALIC_ATOM);
insert(99, "Es",  "einsteinium",       1.00, 1.00, 1.56, "grey", METALIC_ATOM);
insert(100,"Fm",  "fermium",           1.00, 1.00, 1.55, "grey", METALIC_ATOM);
insert(101,"Md",  "mendelevium",       1.00, 1.00, 1.55, "grey", METALIC_ATOM);
insert(102,"No",  "nobelium",          1.00, 1.00, 1.55, "grey", METALIC_ATOM);
}

StandardElementTable::~StandardElementTable() {
	for (std::map<std::string, ElementInfo*>::const_iterator i = _index.begin(); i != _index.end(); ++i)
		delete i->second;
}

bool StandardElementTable::insert(int atno,
		  const char* name, const char* fullName,
		  float eRad, float iRad, float covRad,
		  const char* color, int  flags) {

   if (  eRad > _explMaxRad){ _explMaxRad =   eRad; }
   if (  iRad > _implMaxRad){ _implMaxRad =   iRad; }
   if (covRad >  _covMaxRad){  _covMaxRad = covRad; }

   ElementInfo *ele = new ElementInfo(atno, name, fullName,
			   eRad, iRad, covRad, color, flags);

   _index.insert(std::make_pair(ele->atomName(), ele));
   return TRUE;
}

ElementInfo* StandardElementTable::lookupPDBatom(const char* name, const char* resname) const {
   const char *elementName = NULL;
   bool emitWarning = FALSE;
   char buf[5] = "    ";
   int i;
   for (i = 0; i < 4; i++) { // uppercase the input
      if (name[i] == '\0') { break; }
#ifdef CHARFUNCMACROS
      buf[i] = toupper(name[i]);
#else
      buf[i] = ::toupper(name[i]);
#endif
   }
   buf[i] = '\0';

   switch(buf[0]) {
   case '*': case '\'': case '"': case '`': case '_':
   case '+': case '-':  case ' ':
   case '0': case '1': case '2': case '3': case '4':
   case '5': case '6': case '7': case '8': case '9':
      switch(buf[1]) {
      case 'A':
	 switch(buf[3]) {
	 case '1': elementName = "O"; emitWarning = TRUE; break;
	 case '2': elementName = "N"; emitWarning = TRUE; break;
	 } break;
      case 'B': elementName = "B"; break;
      case 'C': elementName = "C"; break;
      case 'D': elementName = "H"; break;
      case 'F': elementName = "F"; break;
      case 'H': elementName = "H"; break;
	switch(buf[2]) {
     	case 'E': elementName = fixAtomName(name,resname,2) ? "He" : "H"; break;
     	case 'F': elementName = fixAtomName(name,resname,2) ? "Hf" : "H"; break;
	case 'G': elementName = fixAtomName(name,resname,2) ? "Hg" : "H"; break;
     	case 'O': elementName = fixAtomName(name,resname,2) ? "Ho" : "H"; break;
     	case 'S': elementName = fixAtomName(name,resname,2) ? "Hs" : "H"; break;
	default :elementName = "H"; break; 
	} break; 
      case 'I': elementName = "I"; break;
      case 'K': elementName = "K"; break;
      case 'N': elementName = "N"; break;
      case 'O': elementName = "O"; break;
      case 'P': elementName = "P"; break;
      case 'S': elementName = "S"; break;
      case 'U': elementName = "U"; break;
      case 'V': elementName = "V"; break;
      case 'W': elementName = "W"; break;
      case 'Y': elementName = "Y"; break;
      } break;
   case 'A':
      switch(buf[1]) {
      case 'C': elementName = "C";  emitWarning = TRUE;break;//nonstd!
      case 'G': elementName = "Ag"; break;
      case 'H': elementName = "H";  emitWarning = TRUE;break;
      case 'L': elementName = "Al"; break;
      case 'M': elementName = "Am"; break;
      case 'N': elementName = "N";  emitWarning = TRUE;break;
      case 'O': elementName = "O";  emitWarning = TRUE;break;
      case 'P': elementName = "P";  emitWarning = TRUE;break;
      case 'R': elementName = "Ar"; break;
      case 'S': elementName = "As"; break;
      case 'T': elementName = "At"; break;
      case 'U': elementName = "Au"; break;
      } break;
   case 'B':
      switch(buf[1]) {
      case 'A': elementName = "Ba"; break;
      case 'E': elementName = "Be"; break;
      case 'I': elementName = "Bi"; break;
      case 'K': elementName = "Bk"; break;
      case 'R': elementName = "Br"; break;
      } break;
   case 'C': 
      switch(buf[1]) {
      case 'A': elementName = "Ca"; break;
      case 'C': elementName = "C";  emitWarning = TRUE;break;
      case 'D': elementName = "Cd"; break;
      case 'E': elementName = "Ce"; break;
      case 'F': elementName = "Cf"; break;
      case 'H': elementName = "H";  emitWarning = TRUE;break;
      case 'L': elementName = "Cl"; break;
      case 'M': elementName = "Cm"; break;
      case 'N': elementName = "N";  emitWarning = TRUE;break;
      case 'O': elementName = "Co"; emitWarning = TRUE;break;
      case 'P': elementName = "P";  emitWarning = TRUE;break;
      case 'R': elementName = "Cr"; break;
      case 'S': elementName = "Cs"; break;
      case 'U': elementName = "Cu"; break;
      default:  elementName = "C";  emitWarning = TRUE;break;
      } break;
   case 'D':
      switch(buf[1]) {
      case 'Y': elementName = "Dy"; break;
      case 'C': elementName = "C";  emitWarning = TRUE;break;
      case 'H': elementName = "H";  emitWarning = TRUE;break;
      case 'N': elementName = "N";  emitWarning = TRUE;break;
      case 'O': elementName = "O";  emitWarning = TRUE;break;
      case 'P': elementName = "P";  emitWarning = TRUE;break;
      default:  elementName = "H";  emitWarning = TRUE;break;
      } break;
   case 'E':
      switch(buf[1]) {
      case 'R': elementName = "Er"; break;
      case 'S': elementName = "Es"; break;
      case 'U': elementName = "Eu"; break;
      case 'C': elementName = "C";  emitWarning = TRUE;break;
      case 'H': elementName = "H";  emitWarning = TRUE;break;
      case 'N': elementName = "N";  emitWarning = TRUE;break;
      case 'O': elementName = "O";  emitWarning = TRUE;break;
      case 'P': elementName = "P";  emitWarning = TRUE;break;
      } break;
   case 'F':
      switch(buf[1]) {
      case 'E': elementName = "Fe"; break;
      case 'M': elementName = "Fm"; break;
      case 'R': elementName = "Fr"; break;
      case 'C': elementName = "C";  emitWarning = TRUE;break;
      case 'H': elementName = "H";  emitWarning = TRUE;break;
      case 'N': elementName = "N";  emitWarning = TRUE;break;
      case 'O': elementName = "O";  emitWarning = TRUE;break;
      case 'P': elementName = "P";  emitWarning = TRUE;break;
      } break;
   case 'G':
      switch(buf[1]) {
      case 'A': elementName = "Ga"; break;
      case 'D': elementName = "Gd"; break;
      case 'E': elementName = "Ge"; break;
      case 'C': elementName = "C";  emitWarning = TRUE;break;
      case 'H': elementName = "H";  emitWarning = TRUE;break;
      case 'N': elementName = "N";  emitWarning = TRUE;break;
      case 'O': elementName = "O";  emitWarning = TRUE;break;
      case 'P': elementName = "P";  emitWarning = TRUE;break;
      } break;
   case 'H':
      switch(buf[1]) {
      case 'E': elementName = fixAtomName(name,resname,1) ? "He" : "H"; break;
      case 'F': elementName = fixAtomName(name,resname,1) ? "Hf" : "H"; break;
      case 'G': elementName = fixAtomName(name,resname,1) ? "Hg" : "H"; break;
      case 'O': elementName = fixAtomName(name,resname,1) ? "Ho" : "H"; break;
      case 'S': elementName = fixAtomName(name,resname,1) ? "Hs" : "H"; break;
      case 'H': elementName = "H"; break; // hopefully to get rid of warnings for pdb v3 names
      case 'D': elementName = "H"; break; // hopefully to get rid of warnings for pdb v3 names
//      case 'E': elementName = "He"; emitWarning = TRUE;break;
//      case 'F': elementName = "Hf"; emitWarning = TRUE;break;
//      case 'G': elementName = "Hg"; emitWarning = TRUE;break;//xplor Hgamma?
//      case 'O': elementName = "Ho"; emitWarning = TRUE;break;
      default:  elementName = "H"; emitWarning = TRUE;break;
      } break;
   case 'I':
      switch(buf[1]) {
      case 'N': elementName = "In"; break;
      case 'R': elementName = "Ir"; break;
      } break;
   case 'K':
      if (buf[1] == 'R') elementName = "Kr"; break;
   case 'L':
      switch(buf[1]) {
      case 'A': elementName = "La"; break;
      case 'I': elementName = "Li"; break;
      case 'U': elementName = "Lu"; break;
      } break;
   case 'M':
      switch(buf[1]) {
      case 'D': elementName = "Md"; break;
      case 'G': elementName = "Mg"; break;
      case 'N': elementName = "Mn"; break;
      case 'O': elementName = "Mo"; break;
      } break;
   case 'N':
      switch(buf[1]) {
      case 'A': elementName = "Na"; emitWarning = TRUE;break;
      case 'B': elementName = "Nb"; emitWarning = TRUE;break;
      case 'C': elementName = "C";  emitWarning = TRUE;break;
      case 'D': elementName = "Nd"; emitWarning = TRUE;break;
      case 'E': elementName = "Ne"; emitWarning = TRUE;break;
      case 'H': elementName = "H";  emitWarning = TRUE;break;
      case 'I': elementName = "Ni"; break;
      case 'N': elementName = "N";  emitWarning = TRUE;break;
      case 'O': elementName = "O";  emitWarning = TRUE;break;//nonstd!
      case 'P': elementName = "P";  emitWarning = TRUE;break;//nonstd!
      case 'S': elementName = "S";  emitWarning = TRUE;break;
      default:  elementName = "N";  emitWarning = TRUE;break;
      } break;
   case 'O':
      switch(buf[1]) {
      case 'S': elementName = "Os"; break;
      default:  elementName = "O";  emitWarning = TRUE;break;
      } break;
   case 'P':
      switch(buf[1]) {
      case 'A': elementName = "Pa"; emitWarning = TRUE;break;
      case 'B': elementName = "Pb"; emitWarning = TRUE;break;
      case 'D': elementName = "Pd"; emitWarning = TRUE;break;
      case 'M': elementName = "Pm"; break;
      case 'O': elementName = "Po"; break;
      case 'R': elementName = "Pr"; break;
      case 'T': elementName = "Pt"; break;
      case 'U': elementName = "Pu"; break;
      default:  elementName = "P";  emitWarning = TRUE;break;
      } break;
   case 'R':
      switch(buf[1]) {
      case 'A': elementName = "Ra"; break;
      case 'B': elementName = "Rb"; break;
      case 'E': elementName = "Re"; break;
      case 'H': elementName = "Rh"; break;
      case 'N': elementName = "Rn"; break;
      case 'U': elementName = "Ru"; break;
      } break;
   case 'S':
      switch(buf[1]) {
      case 'B': elementName = "Sb"; emitWarning = TRUE;break;
      case 'C': elementName = "Sc"; break;
      case 'E': elementName = "Se"; emitWarning = TRUE;break;
      case 'I': elementName = "Si"; break;
      case 'M': elementName = "Sm"; break;
      case 'N': elementName = "Sn"; break;
      case 'R': elementName = "Sr"; break;
      default:  elementName = "S";  emitWarning = TRUE;break;
      } break;
   case 'T':
      switch(buf[1]) {
      case 'A': elementName = "Ta"; break;
      case 'B': elementName = "Tb"; break;
      case 'C': elementName = "Tc"; break;
      case 'E': elementName = "Te"; break;
      case 'H': elementName = "Th"; break;
      case 'I': elementName = "Ti"; break;
      case 'L': elementName = "Tl"; break;
      case 'M': elementName = "Tm"; break;
      } break;
   case 'X':
      if (buf[1] == 'E') elementName = "Xe"; break;
   case 'Y':
      if (buf[1] == 'B') elementName = "Yb"; break;
   case 'Z':
      switch(buf[1]) {
      case 'N': elementName = "Zn"; break;
      case 'R': elementName = "Zr"; break;
      } break;
   default: break;
   }

   if (elementName == NULL) { emitWarning = TRUE;
      elementName = "C"; // default

      switch(buf[1]) { // punt on names
      case 'H': case 'D': elementName = "H"; break;
      case 'C': elementName = "C"; break;
      case 'N': elementName = "N"; break;
      case 'O': elementName = "O"; break;
      case 'P': elementName = "P"; break;
      case 'S': elementName = "S"; break;

      case 'I': elementName = "I"; break;
      case 'K': elementName = "K"; break;
      case 'V': elementName = "V"; break;
      case 'W': elementName = "W"; break;
      case 'U': elementName = "U"; break;

      case 'A':
	 switch(buf[2]) {
	 case 'G': elementName = "Ag"; break;
	 case 'L': elementName = "Al"; break;
	 case 'S': elementName = "As"; break;
	 case 'U': elementName = "Au"; break;
	 } break;
      case 'F':
	 if (buf[2] == 'E') elementName = "Fe"; break;
      case 'G':
	 if (buf[2] == 'D') elementName = "Gd"; break;
      case 'L':
	 if (buf[2] == 'I') elementName = "Li"; break;
      case 'M':
	 switch(buf[2]) {
	 case 'G': elementName = "Mg"; break;
	 case 'N': elementName = "Mn"; break;
	 case 'O': elementName = "Mo"; break;
	 } break;
      case 'Z':
	 if (buf[2] == 'N') elementName = "Zn"; break;
/* --------- default if we fall through... ---------------------*/
      default: elementName = "C"; break;
      }
   }

   ElementInfo* item = _index.find(elementName)->second;

   if (item && emitWarning) {
      cerr << "WARNING: atom " << name << " from " << resname << " will be treated as "
		  << item->fullName() << endl;
   }
   if (!item) {
      cerr << "WARNING: atom " << name << " from " << resname << " not recognized or unknown"
           << endl;
   }

   return item;
}

ElementInfo* StandardElementTable::element(const char *elementName) const {
	return _index.find(elementName)->second;
}

