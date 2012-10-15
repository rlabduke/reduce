// name: CTab.h
// author: J. Michael Word
// date written: 7/15/97
// purpose: Interface for PDB atom connection table

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

#ifndef CTAB_H
#define CTAB_H 1

#ifdef OLD_STD_HDRS
#include <stdio.h>
#else
#include <cstdio>
using std::FILE;
using std::fclose;
using std::ftell;
using std::fseek;
#endif

#include "AtomConn.h"
#include <list>
#include <map>

// -----------------------------------------
//  atom connections for a given residue
// -----------------------------------------
class ResConn {
public:
	ResConn(const char* resname, int sz_est) {}
	virtual ~ResConn() {
	   for (std::map<std::string, AtomConn*>::const_iterator i = _atomConn.begin(); i != _atomConn.end(); ++i)
		   delete i->second;
	}

   bool put(AtomConn *value) {
	   _atomConn.insert(std::make_pair(value->name(), value));
	   return TRUE;
   }
   AtomConn* get(const std::string &atomname) const {
	   std::map<std::string, AtomConn*>::const_iterator i = _atomConn.find(atomname);
	   if (i != _atomConn.end())
		   return i->second;
	   else
		   return NULL;
   }
   atomPlacementPlan* planHplacement(const std::string &atomname,const char* resname) const;

   // Be careful, std::list must be copied.
   std::list<atomPlacementPlan*> genHplans(const char* resname);

   const char*  resname() const { return _resname;       }

private:
   ResConn(const ResConn&);            // can't copy
   ResConn& operator=(const ResConn&); // can't assign
   std::map<std::string, AtomConn*> _atomConn;
   const char*       _resname;      // name of this residue
};

// -----------------------------------------
//  a specific location in an open file
// -----------------------------------------
class FileLoc {
public:
   FileLoc(FILE *stream, int xtra) : _value(xtra) {
      _location = ::ftell(stream);
   }
   virtual ~FileLoc() {}

   bool relocate(FILE *stream) const{
      return ::fseek(stream, _location, SEEK_SET) == 0;
   }
   int value() { return _value; }
private:

   long _location; // position in file (bytes from the beginning)
   int  _value;    // xtra info. (used to store table size)
};

// -------------------------------------------------------------
//  find residue connection tables in a connection database file
// -------------------------------------------------------------
class CTab {
public:
	CTab(const std::string& dbfile, int szest);
	virtual ~CTab() { 
		if (_fp) ::fclose(_fp); 
		for (std::map<std::string, FileLoc*>::const_iterator i = _filedict.begin(); i != _filedict.end(); ++i)
			delete i->second;
		for (std::map<std::string, ResConn*>::const_iterator j = _rescache.begin(); j != _rescache.end(); ++j)
			delete j->second;
	}

   bool ok() const { return _fp != NULL; }

   ResConn* findTable(const std::string &resname); // note: updates a cache...

   int numConn(const std::string &resname, const std::string &atomname); // note: updates a cache...

private:
   CTab(const CTab&);            // can't copy
   CTab& operator=(const CTab&); // can't assign

   FILE *                 _fp;
   std::map<std::string, FileLoc*>  _filedict;
   std::map<std::string, ResConn*>  _rescache;

   static const int DBbufsz;
};

#endif
