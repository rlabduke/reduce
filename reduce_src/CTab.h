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
#include <stack>
#include <memory>

extern bool UseNuclearDistances; //defined in reduce.cpp JJH

// -----------------------------------------
//  atom connections for a given residue
// -----------------------------------------
class ResConn {
public:
	ResConn() {}

   bool put(std::shared_ptr<AtomConn> value) {
	   _atomConn.insert(std::make_pair(value->name(), value));
	   return TRUE;
   }
   std::shared_ptr<AtomConn> get(const std::string &atomname) const {
     std::shared_ptr<AtomConn> ret;
	   std::map<std::string, std::shared_ptr<AtomConn> >::const_iterator i = _atomConn.find(atomname);
	   if (i != _atomConn.end())
		   ret = i->second;
     return ret;
   }
   // To find out there's a ring bonded to given H atom of methyl group (aromatic ring candidate) - Aram 07/09/12
   std::list<std::string> findRingBondedToMethyl(const std::string &atomname,const char* resname) const;

   std::shared_ptr<atomPlacementPlan> planHplacement(const std::string &atomname,const char* resname) const;

   // Be careful, std::list must be copied.
   std::list<std::shared_ptr<atomPlacementPlan> > genHplans(const char* resname);

private:
   //ResConn(const ResConn&);            // can't copy
   //ResConn& operator=(const ResConn&); // can't assign
   std::map<std::string, std::shared_ptr<AtomConn> > _atomConn;
};

// -------------------------------------------------------------
//  find residue connection tables in a connection database file
// -------------------------------------------------------------
class CTab {
public:
  /// @brief Read residue connection tables from specified file
  /// @param [in] dbfile Name of the file to read from
	CTab(const std::string& dbFileName);

   bool ok() const { return m_rescache.size() != 0; }

   std::shared_ptr<ResConn> findTable(const std::string &resname) const;

   int numConn(const std::string &resname, const std::string &atomname) const;

private:
   //CTab(const CTab&);            // can't copy
   //CTab& operator=(const CTab&); // can't assign

   std::map<std::string, std::shared_ptr<ResConn> >  m_rescache;
};

#endif
