// name: CliqueList.C
// author: J. Michael Word
// date written: 2/7/98
// purpose: Implementation for sequence of interacting residues
//          and singletons

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

#include <iostream>
using std::endl;

#ifdef OLD_STD_HDRS
#include <stdio.h>
#else
#include <cstdio>
using std::sprintf;
#endif

#include "CliqueList.h"
#include "Mover.h"

extern bool GenerateFinalFlip; // SJ - defined in reduce.cpp for checking if flips are being generated for scoring or final PDB file

void CliqueList::describe(std::ostream& os) const {
   int i=0, j=0;
   os << " Singles(size " << _singles.size() << ")";
   for(std::list<MoverPtr>::const_iterator sgl = _singles.begin(); sgl != _singles.end(); ++sgl) {
      if (++i > 4) { i = 0; os << endl << "  "; }
      os << ":" << (*sgl)->descr();
   }
   os << endl;

   i=0;
   for(std::list< std::list<MoverPtr> >::const_iterator ss = _cliques.begin(); ss != _cliques.end(); ++ss) {
      j = 0;
      os << " Set " << ++i << " (size " << (*ss).size() << ")";
	  std::list<MoverPtr> c_list = *ss;
      for(std::list<MoverPtr>::iterator c = c_list.begin(); c != c_list.end(); ++c) {
	 if (++j > 4) { j = 0; os << endl << "  "; }
	 os << ":" << (*c)->descr();
      }
      os << endl;
   }
}

void CliqueList::sortSingletonsByDescr() {
   mpSeqSort(_singles);
}

void CliqueList::formatSingles(std::vector<std::string>& cliqueNotes, AtomPositions& xyz) const { // SJ - 09/25/2015  added the last argument as that is need for doing the final flip

   for(std::list<MoverPtr>::const_iterator s = _singles.begin(); s != _singles.end(); ++s) {

      if ((*s)->canFlip() && GenerateFinalFlip) {
        int orientation = (*s)->orientation();
        (*s)->setOrientation(orientation,xyz);
      }

      cliqueNotes.push_back( (*s)->formatComment("Single ") );
   }
}

void CliqueList::formatClique(std::vector<std::string>& cliqueNotes, int c, AtomPositions& xyz) const { // SJ - 09/25/2015  added the last argument as that is need for doing the final flip

	std::list<MoverPtr> s_list;
	char buf[200];
	int i = 0;

	if (c >= 0 && _cliques.size() > c) {
		i = 0;
		for(std::list< std::list<MoverPtr> >::const_iterator ss = _cliques.begin(); ss != _cliques.end(); ++ss) {
			if (i++ == c) { s_list = *ss; s_list.sort( mpLess() ); break; } //apl adding sort before output begins
		}
	}
	if (s_list.size() > 0) {
		i = 0;
		for(std::list<MoverPtr>::iterator s = s_list.begin(); s != s_list.end(); ++s) {

      if ((*s)->canFlip() && GenerateFinalFlip) {
        int orientation = (*s)->orientation();
        (*s)->setOrientation(orientation,xyz);
      }

      ::sprintf(buf, "Set%2d.%d", c+1, i+1);
      cliqueNotes.push_back( (*s)->formatComment(buf) );

			i++;
		}
	}
}
