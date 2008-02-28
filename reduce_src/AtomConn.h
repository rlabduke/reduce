// name: AtomConn.h
// author: J. Michael Word
// date written: 7/15/97
// purpose: Connected atoms and proton placement plans

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#ifndef ATOMCONN_H
#define ATOMCONN_H 1

#include <vector>
#include "Point3d.h"
#include "ElementInfo.h"
#include "utility.h"

// -----------------
//  atom plan flags
// -----------------

//NOTE:  STRICTALTFLAG is used only to prevent alpha H from worrying about sidechain alternate conformations
#define BONDBUMPFLAG    (1<<0)
#define STRICTALTFLAG   (1<<1)
#define   O2PRIMEFLAG   (1<<2)
#define NOO2PRIMEFLAG   (1<<3)
#define XTRAFLAG        (1<<4)
#define UNSUREDROPFLAG  (1<<5)
#define ROTATEFLAG      (1<<6)
#define ROTATEONDEMAND  (1<<7)
#define AROMATICFLAG    (1<<8)
#define HBACCEPTORFLAG  (1<<9)
#define HBDONORFLAG    (1<<10)
#define IFNOPO4        (1<<11)
#define NH3FLAG        (1<<12)
#define ISACOFLAG      (1<<13)

#define NEGCHARGEFLAG  (1<<14)
#define POSCHARGEFLAG  (1<<15)

// the following two flags are used to mark polar hydrogens where
// XPLOR has a different naming convention

#define NOTXPLORNAME    (1<<16)
#define    XPLORNAME    (1<<17)

// the following two flags are used to mark new vs. old naming conventions
#define USEOLDNAMES     (1<<18)
#define USENEWNAMES     (1<<19)

// the following two flags are used to allow two hydrogens 
// to be built on Calpha of backbone models
#define    NOTBBMODEL	(1<<20)
#define BACKBONEMODEL	(1<<21)

// -----------------------------------------
//  an atom and connected atoms
// -----------------------------------------
class AtomConn {
public:
	AtomConn(const std::string& name, int ord) : _name(name), _order(ord) {_emptyString = "";}
  virtual ~AtomConn() {}

   AtomConn(const AtomConn& a)
	   : _name(a._name), _order(a._order) { _neighbor = a._neighbor; _emptyString = "";}
   AtomConn& operator=(const AtomConn& a);

   const std::string& name() const { return _name; }
   int num_conn() const { return _neighbor.size(); }
   int order() const { return _order; }

   void addConn(const std::string& c) { _neighbor.push_back(c); }

   void limitConnections(int n);
  
   // return the "i"th connection [0..num_conn-1] for this atom
   const std::string& conn(int i) const;
   
   // relational operator allow sequences of AtomConns to be sorted
   bool operator<(const AtomConn& a) const { return _order < a._order; }

private:

   std::string       _name;      // name of this atom
   std::string       _emptyString;
   int          _order;     // number used when sorting AtomConns
   std::vector<std::string> _neighbor;  // list of connected atom names
};

// ---------------------------
//  new atom construction map
// ---------------------------
class atomPlacementPlan : public AtomConn {
public:
   atomPlacementPlan(int t, ElementInfo& e, const AtomConn &c,
			 float d, float a1, float a2, int f)
     : AtomConn(c), _type(t), _elem(e),
                    _dist(d), _ang1(a1), _ang2(a2), _flags(f) {}

   virtual ~atomPlacementPlan() {}

   atomPlacementPlan(const atomPlacementPlan& a);
   atomPlacementPlan& operator=(const atomPlacementPlan& a);

           int type() const { return _type;  }
   const ElementInfo& elem() const { return _elem;  }
         float dist() const { return _dist;  }
int hasFeature(int f) const { return _flags & f; }

bool placeH(const std::vector<Point3d>& loc, Point3d& hpos) const;

   static Point3d calcLoc(int locType,
                      const Point3d& a1, const Point3d& a2,
                      const Point3d& a3, const Point3d& a4,
                      float len, float ang, float xtra);

   static Point3d type1position(
                      const Point3d& a1pos, const Point3d& a2pos,
                      const Point3d& a3pos, const Point3d& a4pos,
                      float bondlen);
   static Point3d type2position(
                      const Point3d& a1pos, const Point3d& a2pos,
                      const Point3d& a3pos,
                      float bondlen, float angle, float fudge);
   static Point3d type3position(
                      const Point3d& a1pos, const Point3d& a2pos,
                      const Point3d& a3pos,
                      float bondlen, float theta, float phi);
   static Point3d type4position(
                      const Point3d& a1pos, const Point3d& a2pos,
                      const Point3d& a3pos,
                      float bondlen, float fudge);
   static Point3d type5position(
                      const Point3d& a1pos, const Point3d& a2pos,
                      const Point3d& a3pos,
                      float bondlen, float fract);
   static Point3d type6position(
                      const Point3d& a1pos, const Point3d& a2pos,
                      float bondlen);
private:
   int         _type;  // connection type
   ElementInfo _elem;  // type of atom
   float       _dist;  // proton to heavy atom distance
   float       _ang1;  // angle one
   float       _ang2;  // angle two
   int         _flags; // flags
};

#endif
