// name: CliqueList.h
// author: J. Michael Word
// date written: 2/7/98
// purpose: Sequences containing collections of motion pointers
//          which refer to interacting residues, along with singletons

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#ifndef CLIQUELIST_H
#define CLIQUELIST_H 1

#include <iostream>
#include <list>
#include <string>
#include <vector>

class Mover;
typedef Mover* MoverPtr;

class CliqueList {
public:
   CliqueList() {}
  ~CliqueList() {}

   CliqueList(const CliqueList& cl) : _cliques(cl._cliques),
                                      _singles(cl._singles) {}

   CliqueList& operator=(const CliqueList& cl) {
      if (this != &cl) {
	 _cliques = cl._cliques; _singles = cl._singles;
      }
      return *this;
   }

   void insertClique(const std::list<MoverPtr>& lst) { _cliques.push_front(lst); }
   void insertSingleton(const MoverPtr    ptr) { _singles.push_front(ptr); }

   void sortSingletonsByDescr();

   const std::list< std::list<MoverPtr> >& cliques() { return _cliques; }
   const std::list<MoverPtr>& singles() { return _singles; }

   int numCliques()    { return _cliques.size(); }
   int numSingletons() { return _singles.size(); }

   void describe(std::ostream& os) const;
   void formatSingles(std::vector<std::string>& cliqueNotes) const;
   void formatClique(std::vector<std::string>& cliqueNotes, int c) const;
private:
	std::list< std::list<MoverPtr> > _cliques;
	std::list<MoverPtr>        _singles;
};
#endif
