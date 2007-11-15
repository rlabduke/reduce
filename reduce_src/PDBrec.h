// name: PDBrec.h
// author: J. Michael Word
// date written: 8/1/97
// purpose: Interface for PDBrec
//
// rmi 10-5-07 modified to allow assignment of Hybrid36 atom numbers to PDB records
//
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
#pragma warning(disable:4800) 
#endif

#ifndef PDBREC_H
#define PDBREC_H 1

#ifdef OLD_STD_HDRS
#include <string.h>
#else
#include <cstring>
#endif
#include "ElementInfo.h"
#include <iostream>
#include <map>
#include "pdb++.h"
#include "Point3d.h"
#include "utility.h"
#include "AtomDescr.h"
#include "hybrid_36_c.h"

// flags pertaining to this particular record
const int   NameModifiedFlag = (1 << 1);
const int NegativeChargeFlag = (1 << 2);
const int PositiveChargeFlag = (1 << 3);

#define WATER_RESNAMES \
 ":HOH:DOD:H2O:D2O:WAT:TIP:SOL:MTO:hoh:dod:h2o:d2o:wat:tip:sol:mto:"

/*#define HE_RESNAME \
// currently there are none RMI 070711

#define HF_RESNAME \
 ":PHF:HF3:HF5:"

#define HG_RESNAMES \
 ": HG:HG2:HGB:HGC:HGI:MAC:MBO:MMC:PHG:PMB:AAS:AMS:BE7:CMH:EMC:EMT:"

#define HO_RESNAME \
 ": HO:HO3:"

#define HS_RESNAME \
// currently there are none RMI 070711*/

class PDBrec {
private:
   static bool _MappingSEGIDtoChains;
   static std::map<std::string, char> _SEGtoChainMap;

   // internal representation of PDBrec which allows them to be shared
   struct PDBrecRep {
      void update(const PDB& r) {
	 _r = r;
	 _recMods = 0;

	 if (r.type() == PDB::ATOM || r.type() == PDB::HETATM) {
	    if (fixupAmbigAtomName(_r.atom.name,
				    _r.atom.residue.name, _r.atom.annotation)) {
	       _recMods |= NameModifiedFlag;
	    }
	    ElementInfo *e = ElementInfo::StdElemTbl().lookupPDBatom(_r.atom.name, _r.atom.residue.name);
	    if (e) { _e = *e; }

	    _recMods |= basicChargeState(_r.atom.name, _r.atom.residue.name,
			   PositiveChargeFlag, NegativeChargeFlag, e);
	 }
	 else if (r.type() == PDB::SIGATM || r.type() == PDB::ANISOU
	       || r.type() == PDB::SIGUIJ) {
	    if (fixupAmbigAtomName(_r.anisou.name,
				    _r.anisou.residue.name, _r.anisou.annotation)) {
	       _recMods |= NameModifiedFlag;
	    }
	    ElementInfo *e = ElementInfo::StdElemTbl().lookupPDBatom(_r.anisou.name, _r.atom.residue.name);
	    if (e) { _e = *e; }
	 }
	 _mark = 0;
	 _ok   = TRUE;
      }
      PDB         _r;
      ElementInfo _e;
      int         _mark;     // utility flag
      bool        _ok;       // is this a valid record?
      int         _recMods;  // adjustments to the record
};

public:
	PDBrec() : _rep(new PDBrecRep()), _i(new int(1)) {}
   PDBrec(const PDB& r) : _rep(new PDBrecRep()), _i(new int(1)) { _rep->update(r); }
  ~PDBrec() {
	  if (--*_i == 0) {
		  delete _rep;
		  delete _i;
	  }
  }

   PDBrec(const PDBrec& a): _rep(a._rep), _i(a._i) {++*_i;}
   PDBrec& operator=(const PDBrec& a) {
	   ++*a._i;
	   if (--*_i == 0) {
		   delete _i;
		   delete _rep;
	   }
	   _i = a._i;
	   _rep = a._rep;
	   return *this; 
   }

   void clone(PDBrec* p, bool setmark = FALSE);

   bool operator==(const PDBrec& a) const { return _rep->_r == a._rep->_r; }

   Coord x() const { return _rep->_r.atom.xyz[0]; }
   Coord y() const { return _rep->_r.atom.xyz[1]; }
   Coord z() const { return _rep->_r.atom.xyz[2]; }
   Point3d loc() const {
      return Point3d(_rep->_r.atom.xyz[0],
		     _rep->_r.atom.xyz[1],
		     _rep->_r.atom.xyz[2]);
   }
   void loc(const Point3d& p) {
      _rep->_r.atom.xyz[0] = p.x();
      _rep->_r.atom.xyz[1] = p.y();
      _rep->_r.atom.xyz[2] = p.z();
   }

   PDB::RecordType type() const { return _rep->_r.type(); }

   const ElementInfo& elem() const { return _rep->_e;  }
   void elem(const ElementInfo& e) { _rep->_e = e; }

   int mark() const { return _rep->_mark;  }
   void mark(int k) { _rep->_mark = k; }

   bool valid() const { return _rep->_ok == TRUE;  }
   void invalidateRecord() {
      _rep->_ok = FALSE;
      _rep->_r.type(PDB::UNKNOWN);
   }
   void partiallyInvalidateRecord() {
      _rep->_ok = FALSE;
   }
   void revalidateRecord() {
      _rep->_ok = TRUE;
   }
   bool atomNameModified() const { return _rep->_recMods & NameModifiedFlag; }

   int modelNum() const { return _rep->_r.model.num; }

   // if atom or hetatm:

   const char* atomname() const { return _rep->_r.atom.name;               }
   const char* Hy36Num()  const { return _rep->_r.atom.serialNum;          }
   int  atomno()          const { 
      int n; 
      const char* errmsg = hy36decode(5,  _rep->_r.atom.serialNum, 5, &n);
      return n; 
      //if (errmsg) throw std::runtime_error(errmsg); 
   }
   const char*  resname() const { return _rep->_r.atom.residue.name;       }
   char insCode()         const { return _rep->_r.atom.residue.insertCode; }
   const char* chain()    const { return _rep->_r.atom.residue.chainId;    }
   char one_char_chain()  const { return _rep->_r.atom.residue.chainId[1]; }
   char alt()             const { return _rep->_r.atom.altLoc;             }
   const char* Hy36resno()const { return _rep->_r.atom.residue.seqNum;     }
   int  resno()           const { 
      int n;
      const char* errmsg = hy36decode(4,  _rep->_r.atom.residue.seqNum, 4, &n);
      return n;
      //if (errmsg) throw std::runtime_error(errmsg);
   }
   float occupancy()      const { return _rep->_r.atom.occupancy;          }
   float tempFactor()     const { return _rep->_r.atom.tempFactor;         }
   const char*segidLabel() const{ return _rep->_r.atom.segID;              }
   const char*elemLabel()  const{ return _rep->_r.atom.element;            }
   const char*chargeLabel()const{ return _rep->_r.atom.charge;             }
   const char*annotation() const{ return _rep->_r.atom.annotation;         }

	AtomDescr getAtomDescr() const {
           int n;
           const char* errmsg = hy36decode(4, (*_rep)._r.atom.residue.seqNum, 4, &n);
           //if (errmsg) throw std::runtime_error(errmsg);
           AtomDescr thedesc( Point3d( (*_rep)._r.atom.xyz[0],(*_rep)._r.atom.xyz[1],(*_rep)._r.atom.xyz[2]), n, vdwRad()); return thedesc;
        }

   void atomname(const char* s)    {
      strncpy(_rep->_r.atom.name, s, 4);
      _rep->_r.atom.name[4] = '\0';
   }

   void atomno(int n)              { 
      const char* errmsg = hy36decode(5,  _rep->_r.atom.serialNum, 5, &n);
      //if (errmsg) throw std::runtime_error(errmsg);
   }
   void Hy36Num(int n)         { 
      char result[6]; 
      const char* errmsg = hy36encode(5, n, result); 
      strncpy(_rep->_r.atom.serialNum, result,6); 
      //_rep->_r.atom.serialNum[6] = '\0'; // RMI i don't know why, but this seems to need to be comment out
      //if (errmsg) throw std::runtime_error(errmsg); 
   }
   void alt(char a)                { _rep->_r.atom.altLoc = a;         }
   void occupancy(float o)         { _rep->_r.atom.occupancy = o;      }
   void tempFactor(float t)        { _rep->_r.atom.tempFactor = t;     }
   void  segidLabel(const char* s);
   void   elemLabel(const char* e) {
      strncpy(_rep->_r.atom.element, e, 2);
      _rep->_r.atom.element[2] = '\0';
   }
   void chargeLabel(const char* c) {
      strncpy(_rep->_r.atom.charge, c, 2);
      _rep->_r.atom.charge[2] = '\0';
   }
   void  annotation(const char* a) {
      strncpy(_rep->_r.atom.annotation,a,10);
      _rep->_r.atom.annotation[10] = '\0';
   }
   void terAtomno(int n)           { 
      const char* errmsg = hy36decode(5,  _rep->_r.atom.serialNum, 5, &n);
      //if (errmsg) throw std::runtime_error(errmsg);
   }
   void x(Coord xval) { _rep->_r.atom.xyz[0] = xval; }
   void y(Coord yval) { _rep->_r.atom.xyz[1] = yval; }
   void z(Coord zval) { _rep->_r.atom.xyz[2] = zval; }

   void getConect(int cvec[]) const; // for connection records
   void setConect(int cvec[]);

   bool hasProp(int p) const { return _rep->_e.hasProp(p);   }
   float      vdwRad() const { return _rep->_e.explRad();    }
   float     implRad() const { return _rep->_e.implRad();    }
   float      covRad() const { return _rep->_e.covRad();     }
   bool   isHydrogen() const { return _rep->_e.isHydrogen(); }

   bool isWater() const;
  
   bool isHe() const; 

   bool isHf() const;

   bool isHg() const;

   bool isHo() const;

   bool isHs() const;

   bool isCharged() const {
      return _rep->_recMods & (PositiveChargeFlag|NegativeChargeFlag);
   }
   bool isNegative() const {
      return _rep->_recMods & NegativeChargeFlag;
   }
   bool isPositive() const {
      return _rep->_recMods & PositiveChargeFlag;
   }
   void setNegative()    { _rep->_recMods |= NegativeChargeFlag; }
   void setPositive()    { _rep->_recMods |= PositiveChargeFlag; }
   void setNotNegative() { _rep->_recMods ^= NegativeChargeFlag; }
   void setNotPositive() { _rep->_recMods ^= PositiveChargeFlag; }

   friend std::ostream& operator << (std::ostream& s, const PDBrec& r);

   void MapSEGIDtoChain(); // updates chain assignment

   static bool MappingSEGIDtoChains() { return _MappingSEGIDtoChains; }
   static int InstallMapOfSEGIDstoChains(const std::string m);
   static void SEGIDtoChain(const char *seg, char c, char* dest);
   static void DumpSEGIDtoChainMap(std::ostream& s, const char *t);
   std::string recName() const;
   std::string stdFormatString() const { return recName(); };

protected:
   const PDB& r() const { return _rep->_r;     } // read-only access
         PDB& w()       { return _rep->_r; } // writable access

private:
   static std::string FormatSegToChainKey(const char *seg); // utility

   PDBrecRep* _rep;
   int* _i;
};

inline std::ostream& operator << (std::ostream& s, const PDBrec& pdbrec) {
   return s << (const PDB&)(pdbrec.r());
}
#endif
