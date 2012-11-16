// name: ElementInfo.h
// author: J. Michael Word     date written: 6/12/97
// purpose: define atom properties

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
#endif

#ifndef ELEMENTINFO_H
#define ELEMENTINFO_H 1

#ifdef OLD_STD_HDRS
#include <string.h>
#else
#include <cstring>
using std::strncpy;
using std::strcmp;
using std::strlen;
#endif
#include <string>
#include <map>

#define METALIC_ATOM    (1 <<  0)
#define DONOR_ATOM      (1 <<  1)
#define ACCEPTOR_ATOM   (1 <<  2)
#define HB_ONLY_DUMMY   (1 <<  3)
#define IGNORE          (1 <<  4)

#define HE_RESNAMES \
 "currently there are none RMI 070711"

#define HF_RESNAMES \
 ":PHF:HF3:HF5:"

#define HG_RESNAMES \
 ": HG:HG2:HGB:HGC:HGI:MAC:MBO:MMC:PHG:PMB:AAS:AMS:BE7:CMH:EMC:EMT:"

#define HO_RESNAMES \
 ": HO:HO3:"

#define HS_RESNAMES \
 "currently there are none RMI 070711"


class StandardElementTable;

// -----------------------------------------
// data associated with a given element type
// -----------------------------------------
class ElementInfo {
private:
   // internal representation of atom types
   // which allows them to be shared

   class ElementInfoRep {
   friend class ElementInfo;
   public:
      ElementInfoRep(int atno, const char* name, const char* fullName,
		  float eRad, float iRad, float covRad,
		  const char* color, int  flags);
      ~ElementInfoRep();

   private:
      int operator==(const ElementInfoRep& r) const {
	 return ::strcmp(_name, r._name) == 0;
      }
      int operator==(const char *s) const {
	 return ::strcmp(_name, s) == 0;
      }

      int _count; // reference count

      int    _atno;  // atomic number
      char*  _name;  // atom name
      char*  _fullName;  // long atom name
      float  _eRad;  // (explicit H) VDW radius
      float  _iRad;  // (implicit H) VDW radius
      float  _covRad;// covalent radius
      char*  _color; // dot color
      int    _flags; // element features
   };

public:
   ElementInfo() {
      _rep = new ElementInfoRep(0,"?","unknown", 0.0,0.0,0.0,"magenta",0);
   }
   ElementInfo(int atno,
		  const char* name, const char* fullName,
		  float eRad, float iRad, float covRad,
		  const char* color, int  flags) {
      _rep = new ElementInfoRep(atno, name, fullName,
			   eRad, iRad, covRad, color, flags);
   }
   virtual ~ElementInfo() { if (--_rep->_count <= 0) delete _rep; }

   // copy constructor
   ElementInfo(const ElementInfo& a) { _rep = a._rep; _rep->_count++; }
   // assignment (shallow copy)
   ElementInfo& operator=(const ElementInfo& a) {
      a._rep->_count++;
      if (--_rep->_count <= 0) delete _rep;
      _rep = a._rep;
      return *this;
   }
   // deep copy
   ElementInfo duplicate() {
      ElementInfo retval(atno(), atomName(), fullName(),
			explRad(), implRad(), covRad(), color(), _rep->_flags);
      return retval;
   }

   // equivalence
   int operator==(const ElementInfo& a) const {
      return operator==(a._rep->_name);
   }
   int operator==(const char *s) const {
      return ::strcmp(_rep->_name, s) == 0;
   }
   int operator!=(const ElementInfo& a) const {
      return operator!=(a._rep->_name);
   }
   int operator!=(const char *s) const {
      return ::strcmp(_rep->_name, s) != 0;
   }

   // access methods
   int   atno()     const { return _rep->_atno; }
   char* atomName() const { return _rep->_name; }
   char* fullName() const { return _rep->_fullName; }
   float explRad()  const { return _rep->_eRad; }
   float implRad()  const { return _rep->_iRad; }
   float covRad()   const { return _rep->_covRad; }
   char* color()    const { return _rep->_color; }

   int   hasProp(int p) const { return _rep->_flags & p; }
   int   isHydrogen()   const { return atno() == 1; }

   static const StandardElementTable& StdElemTbl();

private:
   static const StandardElementTable *TheStdElemTbl; // class shared resource

   ElementInfoRep *_rep;
};

// ----------------------------------------------
// abstract interface for a table of the elements
// ----------------------------------------------
class ElementTable {
protected:
         ElementTable() {}
public:
virtual ~ElementTable() {};

virtual ElementInfo* lookupPDBatom(const char* name, const char* resname) const = 0;

   float maxExplicitRadius() const { return _explMaxRad; }
   float maxImplicitRadius() const { return _implMaxRad; }
   float maxCovalentRadius() const { return _covMaxRad; }

virtual int size() const { return 0; }
protected:

   float _explMaxRad;
   float _implMaxRad;
   float  _covMaxRad;
};

// ------------------------------------------------
// implementation class for a table of the elements
// ------------------------------------------------
class StandardElementTable : public ElementTable {
public:
         StandardElementTable() { LayoutTable(); }
   virtual ~StandardElementTable();

   virtual ElementInfo* lookupPDBatom(const char* name, const char* resname) const;

   virtual int size() const { return _index.size(); }

   ElementInfo* element(const char *elementName) const;

private:
   void LayoutTable(); //define table layout
   StandardElementTable(const StandardElementTable&);           //can't copy
   StandardElementTable& operator=(const StandardElementTable&);//can't assign

   bool insert(int atno,
		  const char* name, const char* fullName,
		  float eRad, float iRad, float covRad,
		  const char* color, int  flags);

   std::map<std::string, ElementInfo*> _index;
};

int basicChargeState(const char* atomname, const char* resname,
	 int posFlag, int negFlag, ElementInfo *e);
bool fixupAmbigAtomName(char* atomname, const char* resname, char* segID);
#endif
