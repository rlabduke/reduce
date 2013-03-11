// name: ResBlk.h
// author: J. Michael Word
// date written: 7/15/97
// purpose: The records for a single residue and assoc. info

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#ifndef RESBLK_H
#define RESBLK_H 1

#include "PDBrec.h"
#include <list>
#include <iostream>

bool sameres(const PDBrec &r1, const PDBrec &rn);

// ----------------
//  A residue block
// ----------------
class ResBlk {
public:
	ResBlk(std::list<PDBrec*>& rlst, std::list<PDBrec*>::iterator& lit);
  ~ResBlk() {}

  bool valid(std::list<PDBrec*>& records) const { return (_insertPt != records.end()) && _resAtoms.size() > 0; }

   unsigned long size() const { return _resAtoms.size(); }

   const PDBrec& firstRec() const { return _firstRec; }

   void addPrevRec(PDBrec* rec) {
	   _resAtoms.insert(std::make_pair("-" + fixupAtomKey(*rec), rec));
   }
   void addNextRec(PDBrec* rec) {
	   _resAtoms.insert(std::make_pair("+" + fixupAtomKey(*rec), rec));
   }

   void get(const std::string& atomname, PDBrec* sp) {
	   std::multimap<std::string, PDBrec*>::const_iterator it = _resAtoms.lower_bound(atomname);
	   if (it != _resAtoms.end())
		   sp = it->second;
	   else
		   sp = NULL;
   }
   void get(const std::string& atomname, std::list<PDBrec*>& sp) {
	   std::multimap<std::string, PDBrec*>::const_iterator it2 = _resAtoms.upper_bound(atomname);
	   for (std::multimap<std::string, PDBrec*>::const_iterator it = _resAtoms.lower_bound(atomname); it != it2; ++it) {
//	   cerr << atomname << endl;
//	   for (std::multimap<std::string, PDBrec*>::const_iterator it = _resAtoms.lower_bound(atomname); it != _resAtoms.end() && it->first == atomname; ++it) {
		   sp.push_front(it->second);
	   }
   }

   bool contains(const std::string& atomname) const {
      return _resAtoms.find(atomname) != _resAtoms.end();
   }

   bool insertNewRec(std::list<PDBrec*>& records, PDBrec* rec) { // add a new atom to the list
	   if (_insertPt != records.end())
		   ++_insertPt;
	   _insertPt = records.insert(_insertPt, rec);
	return TRUE;
   }

   std::multimap<std::string, PDBrec*> atomIt() const {
      return _resAtoms;
   }
private:
   ResBlk(const ResBlk&);            // can't copy
   ResBlk& operator=(const ResBlk&); // can't assign

   std::string fixupAtomKey(const PDBrec &r) const;

   PDBrec                    _firstRec; // the first atom or het record
   std::multimap<std::string, PDBrec*> _resAtoms; // all of the atoms by atomname
   std::list<PDBrec*>::iterator _insertPt; // where we put new atom records
};

#endif
