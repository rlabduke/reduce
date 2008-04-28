// name: ResBlk.C
// author: J. Michael Word
// date written: 8/12/97
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

#include "ResBlk.h"

#ifdef OLD_STD_HDRS
#include <ctype.h>
#include <string.h>
#else
#include <cctype>
#include <cstring>
using std::toupper;
#endif

// are these two records in the same residue?
bool sameres(const PDBrec &r1, const PDBrec &rn) {
   return (::strcmp(r1.resname(), rn.resname()) == 0
	    && r1.resno() == rn.resno()
	    && r1.insCode() == rn.insCode()
	    && strcmp(r1.chain(), rn.chain()) == 0 );
}

// build a residue block by advancing in the list (updates the input iterator)
ResBlk::ResBlk(std::list<PDBrec*>& rlst, std::list<PDBrec*>::iterator& lit) {
	_insertPt = lit;

	while (_insertPt != rlst.end()) {
		PDBrec* r1 = *_insertPt;

		++_insertPt;

		if (r1->type() == PDB::ATOM
			|| r1->type() == PDB::HETATM) { // find the first atom

			_firstRec = *r1;

			_resAtoms.insert(std::make_pair(fixupAtomKey(*r1), r1));

			while (_insertPt != rlst.end()) {  // now find the rest

				PDBrec* rn = *_insertPt;

				if (rn->type() == PDB::ATOM
					|| rn->type() == PDB::HETATM) {

					if (! sameres(*r1, *rn)) { break; }

					_resAtoms.insert(std::make_pair(fixupAtomKey(*rn), rn));
				}
				else if (rn->type() == PDB::SIGATM
					|| rn->type() == PDB::ANISOU
					|| rn->type() == PDB::SIGUIJ) {

					if (! sameres(*r1, *rn)) { break; }
				}
				else { break; } // residue ends at a non-(atom, etc.)

				++_insertPt;
			}
			break;
		}
	}
	lit = _insertPt;
//	lit.sync(_insertPt); // *tell the outside world where to look next*

	// need to keep a pointer to the end of the res.
	//apl 2007/03/07.  March iterator backwards even if the list end has been reached.
	//if (_insertPt != rlst.end()) {--_insertPt;}
	--_insertPt;
//	if (_insertPt) { _insertPt--;     }
//	else           { _insertPt.end(); }
}

// creates the name string we use as a key for an atom based on atomname
std::string ResBlk::fixupAtomKey(const PDBrec &r) const {
   char buf[5];
   const char *p = r.atomname();
   int i = 0;

   for(i=0; i < 4 && p[i] != '\0'; i++) { // copy out up to four characters
      buf[i] = toupper(p[i]);             // and make sure they are uppercase
   }
   buf[i] = '\0';

   if( r.isHydrogen() ) { // try to change deuterium to hydrogen
      // this will miss some deuterated files in xplor format

      if (buf[0] == 'D' && buf[1] != 'D' && buf[1] != 'H') {
	 buf[0] = 'H';
      }
      else if (buf[0] != 'D' && buf[0] != 'H' && buf[1] == 'D') {
	 buf[1] = 'H';
      }
      // fix deuterium bug -  080428 by JJH
      // for version 3.0 Hydrogens/Deuteriums with 4 character names
     else if (buf[0] != ' ' && buf[1] != ' ' && buf[2] != ' ' && buf[3] != ' ' && buf[0] == 'D'){
         buf[0] = 'H';
     }
   }

   if((::strcmp(r.resname(), "NH2" ) == 0)
   && (::strcmp("1H", buf) == 0 || ::strcmp("2H", buf) == 0)) {
      buf[2] = 'N';	      // use the standard names 1HN,2HN
      buf[3] = '\0';
   }
   else if (::strcmp("2OXT", buf) == 0) {
      buf[0] = ' ';
      buf[1] = 'O';	      // store under the standard name O
      buf[2] = '\0';
   }
   else if (::strcmp(" NT", buf) == 0) {
      buf[0] = ' ';
      buf[1] = 'N';	      // store under the standard name N
      buf[2] = '\0';
   }
   else if ((buf[1] == 'A') &&
         ( ((::strcmp(r.resname(), "ASN" ) == 0) && (buf[2] == 'D'))
      ||   ((::strcmp(r.resname(), "GLN" ) == 0) && (buf[2] == 'E')) )) {

	 // arbitrary decision about ambiguous assignments O_1 & N_2

         if      (buf[3] == '1') { buf[1] = 'O'; }
         else if (buf[3] == '2') { buf[1] = 'N'; }
   }

   return buf;
}
