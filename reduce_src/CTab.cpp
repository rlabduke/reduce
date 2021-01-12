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

std::list<std::string> ResConn::findRingBondedToMethyl(const std::string &atomname,const char* resname) const {
	int nh = 0;

	ElementInfo *e = ElementInfo::StdElemTbl().lookupPDBatom(atomname.c_str(), resname);
	if (e && e->isHydrogen()) { // start only from hydrogens

		std::stack<std::shared_ptr<AtomConn> > stkAtoms;
		std::stack<int> stkAtomsDepth;
		std::vector<bool> visited(_atomConn.size());
		std::list<std::string> L;

		// initialize all nodes to false
		for (int i=0; i < visited.size(); i++)
			visited.at(i)=false;

		std::shared_ptr<AtomConn> hc = get(atomname);
		if (hc && hc->num_conn() > 0)  {

			std::string x1name = hc->conn(0); // heavy atom
			std::shared_ptr<AtomConn> x1c = get(x1name);
			int x1cDepth = 0;

			ElementInfo *x1e = ElementInfo::StdElemTbl().lookupPDBatom(x1name.c_str(), resname);
			if (x1c && x1e->atno()==6) { // Test to see if this hydrogen is in methyl group
				std::string xc = "";
				for (int i = x1c->num_conn()-1; i >= 0; i--) {
					ElementInfo *xe = ElementInfo::StdElemTbl().lookupPDBatom(x1c->conn(i).c_str(), resname);

					if (xe && xe->isHydrogen()) { // Count number of hydrogen atoms
						++nh;
					}

					if (xe && !xe->isHydrogen()) { // heavy atom
						xc = x1c->conn(i);
					}
				}

				//std::cout << std::endl << "nh: " << nh << " xc: " << xc;

				if (nh == 3) {
					visited[x1c->order()]=true;
					std::shared_ptr<AtomConn> x2c = get(xc);
					stkAtoms.push(x2c);
					stkAtomsDepth.push(x1cDepth+1);
					x1name = xc;
				}
			}

			while (!stkAtoms.empty()) {
				x1c = stkAtoms.top();
				x1cDepth = stkAtomsDepth.top();
				stkAtoms.pop();
				stkAtomsDepth.pop();

				if ( x1c && x1c->num_conn() > 1    // Atom exists and has more than 1 connected atoms (other than it's parent)
					 && x1cDepth < 7               // only report 5 or 6 member rings
					 && x1cDepth > 4 )  {          // exclude 3 member rings - JJH 130326

					// std::cout << std::endl << "here: " << x1c->order() << x1c->name() << "-" << x1cDepth;
					visited[x1c->order()]=true;
					if (L.size() > x1cDepth-1)
						L.resize(x1cDepth-1);
					std::string xparent = L.empty() ? "NONE" : L.back();
					L.push_back(x1c->name());

					for (int i = x1c->num_conn()-1; i >= 0; i--) {
						std::string xc = x1c->conn(i);
						if (xparent != xc) {
							if (xc == x1name) { // Found cycle
								return L;
							}

							std::shared_ptr<AtomConn> x2c = get(xc);
							if (x2c && !visited[x2c->order()]) {
								// std::cout << "(x2c: " << x2c->order() << x2c->name();
								stkAtoms.push(x2c);
								stkAtomsDepth.push(x1cDepth+1);
							}
						}
					}
				}
			}
		}
	}

	std::list<std::string> emptyList;
	return emptyList;
}

std::shared_ptr<atomPlacementPlan> ResConn::planHplacement(const std::string &atomname,
                                           const char* resname) const {
   int nn = 0, nh = 0, whichH = 0, type = 0, flags = 0;
   float dist = 0.0, ang1 = 0.0, ang2 = 0.0;

   ElementInfo *e = ElementInfo::StdElemTbl().lookupPDBatom(atomname.c_str(), resname);
   if (e && e->isHydrogen()) {
	   std::shared_ptr<AtomConn> hc = get(atomname);
      if (hc && hc->num_conn() > 0)  {
		  AtomConn names(atomname, hc->order());
	 std::string x1name = hc->conn(0);
	 names.addConn(x1name); // H atom connections will be heavy atoms
	 std::shared_ptr<AtomConn> x1c = get(x1name);
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
	    // default X-H distance
	    if (UseNuclearDistances) {
	      dist = 1.09;
	    }
	    else{
	      dist = 0.97;
	    }
	    bool polarH = FALSE;
	    if (x1e) {
	      if (x1e->atno()==7) { // N
	        if (UseNuclearDistances) {
		      dist = 1.02; // nuclear
		    }
		    else {
		      dist = 0.86; // electron cloud
		    }
	        polarH = TRUE;
	      }
	      else if (x1e->atno()==8) { // O
		    if (UseNuclearDistances) {
		      dist = 0.98; // nuclear
		    }
		    else {
		      dist = 0.84; // electron cloud
		    }
	        polarH = TRUE;
	      }
	      else if (x1e->atno()==16) { // S
		    if (UseNuclearDistances) {
		      dist = 1.3; // nuclear
		    }
		    else {
		      dist = 1.2; // electron cloud
		    }
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
			std::shared_ptr<AtomConn> x2c = get(names.conn(1));
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

	       return std::make_shared<atomPlacementPlan>(type, *e, names, dist,
	                                        ang1, ang2, flags);
	    }
	    else { cerr << "ERROR ResConn::connect(" << atomname
	                << "): unresolved" << endl; }
	 }
      }
   }
   std::shared_ptr<atomPlacementPlan> ret;
   return ret;
}

std::list<std::shared_ptr<atomPlacementPlan> > ResConn::genHplans(const char* resname) {
	std::list<std::shared_ptr<atomPlacementPlan> > plans;
	std::map<std::string, std::shared_ptr<AtomConn> >::const_iterator i = _atomConn.begin();

	while (i != _atomConn.end()) {
		std::shared_ptr<atomPlacementPlan> p = planHplacement(i->first,resname);
		if (p) {
			plans.push_front(p); // store a copy of this part of the plan
		}
		++i;
	}
//	return sort(plans);
	return plans;
}


CTab::CTab(const std::string& dbFileName)
{
  const unsigned DBbufsz = 500;
  FILE *                 fp;

	char buf[DBbufsz+1], resname[4];
	int n;

	fp = ::fopen(dbFileName.c_str(), "r");
	if (fp == NULL) {
		cerr << "ERROR CTab(" << dbFileName << "): could not open" << endl;
		return;
	}

	// Parse the entries in the file.  Each residue starts with the
	// keyword RESIDUE at the start of its line.  There are a number
	// of CONECT lines associated with a residue and then END, or a
	// new RESIDUE, or the end of the file finished the residue.
	// There are other lines (HET HETSYN FORMUL, blank) that should be
	// ignored.
	bool inResidue = false;
	std::shared_ptr<ResConn> curResidue;
	std::string curName;
	while (::fgets(buf, DBbufsz, fp)) {
    
		// Take action based on the type of line we find.
		if (strncasecmp_cp(buf, "RESIDUE", 7) == 0) {
			// Whenever we find a RESIDUE line, we finish any existing
			// residue.
			if (curResidue) {
				m_rescache.insert(std::make_pair(curName, curResidue));
				curResidue.reset();
			}
			if (0 > column_sscanf(buf, "%10 %3s %6d", resname, &n)) {
				cerr << "ERROR CTab(" << dbFileName
					<< ", ): scan error on residue line: " << buf
					<< " (skipping residue)" << endl;
			} else {
				// Whenever we find a valid RESIDUE line, we start a new
				// residue.
				curName = resname;
				curResidue = std::make_shared<ResConn>();
			}

		} else if (strncasecmp_cp(buf, "END", 3) == 0) {
			// Finish any existing residue.
			if (curResidue) {
				m_rescache.insert(std::make_pair(curName, curResidue));
				curResidue.reset();
			}

		} else if (strncasecmp_cp(buf, "CONECT", 6) == 0) {
			char an[10][5];     // holds ten strings of length 4
			int m = 0, n = 0;

			// Add the line into the current residue
			if (0 > column_sscanf(buf, "%11 %4s %4d%4s %4s %4s %4s %4s %4s %4s %4s %4s",
					an[0], &n,  an[1], an[2], an[3], an[4],
					an[5], an[6], an[7], an[8], an[9])) {
				cerr << "ERROR CTab(" << dbFileName << "): Error in residue " << curName
					<< ": scan error in line: " << buf << endl;
			} else {
				if (curResidue) {
					std::shared_ptr<AtomConn> connectedAtoms = std::make_shared<AtomConn>(an[0],++m);
					curResidue->put(connectedAtoms);
					if (n > 9) { n = 9; } // overflow
					for (int i=1; i <= n; i++) {
						connectedAtoms->addConn(an[i]);
					}
				} else {
					cerr << "ERROR CTab(" << dbFileName << "): Error CONECT outide of residue"
						<< ": scan error in line: " << buf << endl;
				}
			}

		} else {
			// Ignore all other lines
		}

	}
  
	// Insert the last residue, if there was one.
	if (curResidue) {
		m_rescache.insert(std::make_pair(curName, curResidue));
		curResidue.reset();
	}

  ::fclose(fp);
}

std::shared_ptr<ResConn> CTab::findTable(const std::string &resname) const {

  std::shared_ptr<ResConn> ret;
	std::map<std::string, std::shared_ptr<ResConn> >::const_iterator iter = m_rescache.find(resname);
	if (iter != m_rescache.end()) {
		ret = iter->second;
  }

  return ret;
}

int CTab::numConn(const std::string &resname, const std::string &atomname) const {
  std::shared_ptr<ResConn> tbl = findTable(resname);
  if (tbl) {
    std::shared_ptr<AtomConn> hc = tbl->get(atomname);
    if (hc) {
      return hc->num_conn();
    }
  }
  return 0;
}
