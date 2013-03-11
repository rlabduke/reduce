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

void CliqueList::formatSingles(std::vector<std::string>& cliqueNotes) const {
   char buf[200];

   for(std::list<MoverPtr>::const_iterator s = _singles.begin(); s != _singles.end(); ++s) {

      if ((*s)->type() == Mover::ROTATE_METHYL) {
	 ::sprintf(buf, "USER  MOD Single :%s:%s:sc=%8.3g%c  (180deg=%.3g%s)",
	    (*s)->descr().c_str(), (*s)->describeOrientation().c_str(),
	    (*s)->bestScore(),(((*s)->bestHasBadBump())?'!':' '),
	    (*s)->initScore(),(((*s)->initHasBadBump())?"!":""));
      }
      else if ((*s)->canFlip()) {
	 char classcode = ' ';
	 if ((*s)->flipStateHasBadBump(0) && (*s)->flipStateHasBadBump(1)) {
	    classcode = 'C';
	 }
	 else if (abs((*s)->flipMaxScore(0) - (*s)->flipMaxScore(1)) <= 0.5) {
	    classcode = 'X';
	 }
	 else if ((*s)->flipState() != 0) {
	    classcode = 'F';
	 }
	 else {
	    classcode = 'K';
	 }
	 ::sprintf(buf, "USER  MOD Single :%s:%s:sc=%8.3g%c %c(o=%.2g%s,f=%.2g%s)",
	    (*s)->descr().c_str(), (*s)->describeOrientation().c_str(),
	    (*s)->bestScore(),(((*s)->bestHasBadBump())?'!':' '), classcode,
	    (*s)->flipMaxScore(0),(((*s)->flipStateHasBadBump(0))?"!":""),
	    (*s)->flipMaxScore(1),(((*s)->flipStateHasBadBump(1))?"!":""));
      }
      else {
	 ::sprintf(buf, "USER  MOD Single :%s:%s:sc=%8.3g%c",
	    (*s)->descr().c_str(), (*s)->describeOrientation().c_str(),
	    (*s)->bestScore(), (((*s)->bestHasBadBump())?'!':' '));
      }
      cliqueNotes.push_back(buf);
   }
}

void CliqueList::formatClique(std::vector<std::string>& cliqueNotes, int c) const {
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

			if ((*s)->type() == Mover::ROTATE_METHYL) {
				::sprintf(buf, "USER  MOD Set%2d.%d:%s:%s:sc=%8.3g%c  (180deg=%.3g%s)",
					c+1, i+1,
					(*s)->descr().c_str(), (*s)->describeOrientation().c_str(),
					(*s)->bestScore(),(((*s)->bestHasBadBump())?'!':' '),
					(*s)->initScore(),(((*s)->initHasBadBump())?"!":""));
			}
			else if ((*s)->canFlip()) {
				char classcode = ' ';
				if ((*s)->flipStateHasBadBump(0) && (*s)->flipStateHasBadBump(1)) {
					classcode = 'C';
				}
				else if (abs((*s)->flipMaxScore(0) - (*s)->flipMaxScore(1)) <= 0.5) {
					classcode = 'X';
				}
				else if ((*s)->flipState() != 0) {
					classcode = 'F';
				}
				else {
					classcode = 'K';
				}
				::sprintf(buf, "USER  MOD Set%2d.%d:%s:%s:sc=%8.3g%c %c(o=%.2g%s,f=%.2g%s)",
					c+1, i+1,
					(*s)->descr().c_str(), (*s)->describeOrientation().c_str(),
					(*s)->bestScore(),(((*s)->bestHasBadBump())?'!':' '), classcode,
					(*s)->flipMaxScore(0),(((*s)->flipStateHasBadBump(0))?"!":""),
					(*s)->flipMaxScore(1),(((*s)->flipStateHasBadBump(1))?"!":""));
			}
			else {
				::sprintf(buf, "USER  MOD Set%2d.%d:%s:%s:sc=%8.3g%c",
					c+1, i+1,
					(*s)->descr().c_str(), (*s)->describeOrientation().c_str(),
					(*s)->bestScore(), (((*s)->bestHasBadBump())?'!':' '));
			}
			cliqueNotes.push_back(buf);
			i++;
		}
	}
}
