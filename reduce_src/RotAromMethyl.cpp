// name: RotAromMethyl.C
// author: J. Michael Word, modified by Aram Han
// date written: 2/7/98, modified 8/13/12
// purpose: Implementation for RotAromMethyl

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#ifdef OLD_STD_HDRS
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#else
#include <cstdio>
#include <cctype>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
using std::cerr;
using std::endl;
using std::strcpy;
using std::exit;
#include <sstream>
using std::cout;
#endif

#include "RotAromMethyl.h"
#include "AtomPositions.h"

#define START_ANGLE 150.0
#define ROUGH_STEP  180

RotAromMethyl::RotAromMethyl(const Point3d& a, const Point3d& b,
                     const double ang, const PDBrec& heavyAtom)
   : RotMethyl(a, b, ang, heavyAtom)
{
}

int RotAromMethyl::numOrientations(SearchStrategy ss) const
{
	return 2;
}

bool RotAromMethyl::setOrientation(int oi, float delta, AtomPositions &xyz,
                                                     SearchStrategy ss)
{
	const double oldTheta = angle();
	const double    theta = orientationAngle(oi, Mover::LOW_RES) + delta;

	if (abs(theta-oldTheta) > 0.1) {
		for(std::list< std::shared_ptr<PDBrec> >::const_iterator hydlist = _rot.begin(); hydlist != _rot.end(); ++hydlist) {
			setHydAngle(**hydlist, oldTheta, theta, xyz);
		}
		angle(theta);
	}
	rememberOrientation(oi, ss);
	return TRUE;
}
