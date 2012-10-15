// name: CTab.C
// author: J. Michael Word
// date written: 7/15/97
// purpose: Implementation for PDB atom connection table

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
#pragma warning(disable:4800) 
#endif

#include <iostream>
using std::cerr;
using std::endl;

#ifdef OLD_STD_HDRS
#include <stdio.h>
#else
#include <cstdio>
using std::fopen;
using std::fgets;
#endif

#include "CTab.h"
#include "ElementInfo.h"
#include "utility.h"

atomPlacementPlan* ResConn::planHplacement(const std::string &atomname, const char* resname) const {
   int nn = 0, nh = 0, whichH = 0, type = 0, flags = 0;
   float dist = 0.0, ang1 = 0.0, ang2 = 0.0;

   ElementInfo *e = ElementInfo::StdElemTbl().lookupPDBatom(atomname.c_str(), resname);
   if (e && e->isHydrogen()) {
      AtomConn *hc = get(atomname);
      if (hc && hc->num_conn() > 0)  {
		  AtomConn names(atomname, hc->order());
	 std::string x1name = hc->conn(0);
	 names.addConn(x1name); // H atom connections will be heavy atoms
	 AtomConn *x1c = get(x1name);
	 if (x1c && x1c->num_conn() > 0)  {
	    for (int i = 0; i < x1c->num_conn(); i++) {
			std::string xc = x1c->conn(i);
	       ElementInfo *xe = ElementInfo::StdElemTbl().lookupPDBatom(xc.c_str(), resname);
	       if (xe && xe->isHydrogen()) {
		  ++nh;
		  if (xc == names.name()) { whichH = nh; }
	       }
	       else {
		  names.addConn(xc);
		  nn++;
	       }
	    }

	    ElementInfo *x1e = ElementInfo::StdElemTbl().lookupPDBatom(x1name.c_str(), resname);
	    dist = 1.1; // default X-H distance
	    bool polarH = FALSE;
	    if (x1e) {
	       if (x1e->atno()==7 || x1e->atno()==8) { // N & O
		  dist = 1.0;
		  polarH = TRUE;
	       }
	       if (x1e->atno()==16) { // S
		  dist = 1.3;
		  polarH = TRUE;
	       }
	    }

	    type = 0;
	    if (nn == 3) { type = 1; }
	    else if (nn == 2) {
	       if (x1c->num_conn() == 4) {
		  type = 2;
		  ang1 = (whichH == 1) ? 126.5 : -126.5;
	       }
	       else if (x1c->num_conn() == 3) { type = 4; }
	    }
	    else if (nn == 1) {
		  AtomConn *x2c = get(names.conn(1));
		  if (x2c && x2c->num_conn() > 0)  {
		     for (int k = 0; k < x2c->num_conn(); k++) {
			std::string xxc = x2c->conn(k);
			if (xxc != x1name) {
			   ElementInfo *xxe = ElementInfo::StdElemTbl().lookupPDBatom(xxc.c_str(), resname);
			   if (!xxe || !(xxe->isHydrogen())) {
				   names.addConn(xxc);
			      if (x1c->num_conn() == 4) {
				 type = 3; 
				 ang1 = 109.5;
				 ang2 = (whichH == 1) ? 180.0
				      : (whichH == 2) ? 60.0 : -60.0;
				 std::string x2name = names.conn(1);
				 ElementInfo *x2e = ElementInfo::StdElemTbl().lookupPDBatom(x2name.c_str(), resname);
				 if (x2e->atno()==16) { // S, so this is like a MET methyl
				    flags |= ROTATEFLAG;
				 }
				 else {
				    flags |= ROTATEONDEMAND;
				    if (x1e->atno()==7) { // NH3
				       flags |= NH3FLAG;
				    }
				 }
			      }
			      else if (x1c->num_conn() == 3) {
				 type = 3; 
				 ang1 = 120.0;
				 ang2 = (whichH == 1) ? 0.0 : 180.0;
			      }
			      else if (x1c->num_conn() == 2) {
				 if (x1e && (x1e->atno()==6)) {// Carbon
				    type = 6; // linear (triple-bond)
				 }
				 else {
				    type = 3; 
				    ang1 = 109.5;
				    ang2 = 180.0;
				    if (polarH) {
				       flags |= UNSUREDROPFLAG|ROTATEFLAG;
				    }
			         }
			      }
			      break;
			   }
			}
		     }
		  }
		  if (type == 0) { type = 7; } // no atom to make dihedral - just adjust length
		  // ** in this case, the rotamer is unclear
	    }
	    else if (nn == 0) { type = 7; } // water, etc. (only one heavy atom) - just adjust length

	    if (type != 0) {
	       ElementInfo *ederrived = NULL;

	       // Here we make the assumption that a type 4 proton is
	       // aromatic knowing full well this is sometimes incorrect.

	       if (type == 4 && polarH) {
		  ederrived = ElementInfo::StdElemTbl().element("Ha+p");
	       }
	       else if (type == 4) {
		  ederrived = ElementInfo::StdElemTbl().element("Har");
	       }
	       else if (polarH) {
		  ederrived = ElementInfo::StdElemTbl().element("Hpol");
	       }

	       if (ederrived) { e = ederrived; }

	       switch(type) {
	       case 1: names.limitConnections(4); break;
	       case 2: names.limitConnections(3); break;
	       case 3: names.limitConnections(3); break;
	       case 4: names.limitConnections(3); break;
	       case 5: names.limitConnections(3); break;
	       case 6: names.limitConnections(2); break;
	       case 7: names.limitConnections(1); break;
	       }

	       flags |= BONDBUMPFLAG;

	       return new atomPlacementPlan(type, *e, names, dist,
	                                        ang1, ang2, flags);
	    }
	    else { cerr << "ERROR ResConn::connect(" << atomname
	                << "): unresolved" << endl; }
	 }
      }
   }
   return NULL;
}

std::list<atomPlacementPlan*> ResConn::genHplans(const char* resname) {
	std::list<atomPlacementPlan*> plans;
	std::map<std::string, AtomConn*>::const_iterator i = _atomConn.begin();
	
	while (i != _atomConn.end()) {
		atomPlacementPlan *p = planHplacement(i->first, resname);
		if (p) {
			plans.push_front(p); // store a copy of this part of the plan
		}
		++i;
	}
//	return sort(plans);
	return plans;
}

// (excessive) record size for reading connection database file
const int CTab::DBbufsz = 500;

CTab::CTab(const std::string& dbfile, int sz) {
	char buf[DBbufsz+1], resname[4];
	int n;

	_fp = ::fopen(dbfile.c_str(), "r");

	if (_fp == NULL) {
		cerr << "ERROR CTab(" << dbfile << "): could not open" << endl;
		return;
	}

	while (::fgets(buf, DBbufsz, _fp)) {
		if (::strncasecmp(buf, "RESIDUE", 7) == 0) {
			if (0 > column_sscanf(buf, "%10 %3s %6d", resname, &n)) {
				cerr << "ERROR CTab(" << dbfile << ", " << sz
					<< "): scan error" << endl;
			}
			else {
				FileLoc *loc = new FileLoc(_fp, n);
				_filedict.insert(std::make_pair(std::string(resname), loc));
			}
		}
	}
}

// not const because it modifies the cache
ResConn* CTab::findTable(const std::string &resname) {
	char buf[DBbufsz+1];
	char an[10][5];     // holds ten strings of length 4
	int m = 0, n = 0;
        //char temp[4];

	std::map<std::string, ResConn*>::iterator iter1 = _rescache.find(resname);
	if (iter1 != _rescache.end())
		return iter1->second;
	//   ResConn *currTbl = _rescache.get(resname);

	std::map<std::string, FileLoc*>::iterator iter2 = _filedict.find(resname);
	FileLoc *loc;
	if (iter2 != _filedict.end())
		loc = iter2->second;
	else
		loc = NULL;
	if (loc == NULL || _fp == NULL || loc->relocate(_fp) == FALSE)
		return NULL;

	ResConn *currTbl = new ResConn(resname.c_str(), loc->value());
	_rescache.insert(std::make_pair(resname, currTbl));
//	_rescache.put(resname, currTbl);

	while (::fgets(buf, DBbufsz, _fp)) {
		if (::strncasecmp(buf, "CONECT", 6) == 0) {
			if (0 > column_sscanf(buf, "%11 %4s %4d%4s %4s %4s %4s %4s %4s %4s %4s %4s",
				an[0], &n,  an[1], an[2], an[3], an[4],
				an[5], an[6], an[7], an[8], an[9])) {
				cerr << "ERROR CTab::findTable(" << resname
					<< "): scan error" << endl;
				break;
			}
			else {
				AtomConn *connectedAtoms = new AtomConn(an[0],++m);
				currTbl->put(connectedAtoms);
				if (n > 9) { n = 9; } // overflow
				for (int i=1; i <= n; i++) {
					connectedAtoms->addConn(an[i]);
				}
			}
		}
		else if (::strncasecmp(buf, "END", 3) == 0) {
			break;
		}
		else if (::strncasecmp(buf, "RESIDUE", 7) == 0) {
			break; // starting next residue
		}
	}
	return currTbl;
}

// not const because it modifies the cache

int CTab::numConn(const std::string &resname, const std::string &atomname) {
   ResConn *tbl = findTable(resname);
   if (tbl) {
      AtomConn *hc = tbl->get(atomname);
      if (hc) { return hc->num_conn(); }
   }
   return 0;
}
