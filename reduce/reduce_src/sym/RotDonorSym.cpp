// name: RotDonorSym.cpp
// author: J. Michael Word
// date written: 2/7/98
// purpose: Symmetry-aware implementation for RotDonor

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
#include <stdlib.h>
#else
#include <cstdio>
#include <cctype>
#include <cmath>
#include <cstdlib>
using std::exit;
#endif

#include "../RotDonor.h"
#include "../AtomPositions.h"

// for isHBDonorOrAcceptorFlipped
#include "../FlipMemo.h"

#define START_ANGLE 180.0
#define ROUGH_STEP   10
#define FINE_STEP     1

#define OVERLAPCUTOFF_1    -0.001
#define OVERLAPCUTOFF_3    -0.7
#define OVERLAPCUTOFF_5    -0.9

void RotDonor::finalize(int nBondCutoff, bool useXplorNames, bool useOldNames, bool bbModel, 
                        AtomPositions &xyz, DotSphManager& dotBucket) {
	int i = 0;

	const double origAngle = angle();

	if (isComplete()) {

		// pre-build lists of bonded atoms

		const double approxNbondDistLimit = 3.0 + 0.5*nBondCutoff; // just a rule of thumb
		i = 0;
		_rot.push_front(&_heavyAtom);
		for (std::list<PDBrec*>::const_iterator alst = _rot.begin(); alst != _rot.end(); ++alst) {
			PDBrec* thisAtom = *alst;
			std::list<PDBrec*> *temp = new std::list<PDBrec*>();
			std::list< std::pair<PDBrec*, Point3d> > neighborList = xyz.neighbors(thisAtom->loc(), thisAtom->covRad(),
				approxNbondDistLimit);
			bondedList(thisAtom, neighborList, nBondCutoff, _rot, temp);
			_bnded.push_back(temp);
		}
		_rot.pop_front();
		// -------------------------------------------------------
		// determine all the angles we will consider

		// first identify nearby acceptors

		std::list<double> angs;

		const ElementInfo& elemHpol = * ElementInfo::StdElemTbl().element("Hpol");

		// O-H & N-H bond len == 1.0, S-H bond len == 1.3
		const double XHbondlen = (_heavyAtom.elem().atno() == 16) ? 1.3 : 1.0;

		std::list< std::pair<PDBrec*, Point3d> > nearby_list = xyz.neighbors(_heavyAtom.loc(), 0.001, 4.0);
		PDBrec* rec = NULL;
		for(std::list< std::pair<PDBrec*,Point3d> >::const_iterator nearby = nearby_list.begin(); nearby != nearby_list.end(); ++nearby) {
			rec = nearby->first;
			//rec = *nearby;
			if (rec->hasProp(ACCEPTOR_ATOM)
				|| FlipMemo::isHBDonorOrAcceptorFlipped(*rec, useXplorNames, useOldNames, bbModel)) {

				double HBoverlap = distance2(_heavyAtom.loc(), nearby->second)
					- ( elemHpol.explRad() + rec->vdwRad()
					+ XHbondlen );

				if (HBoverlap < OVERLAPCUTOFF_1
					&& ! diffAltLoc(_heavyAtom, *rec)) {

					const double ang = origAngle +
						dihedral(nearby->second, _p1, _p2, (*(_rot.begin()))->loc());
//						dihedral((*nearby).loc(), _p1, _p2, (*_rot).loc());

					angs.push_front(clampAngle(ang));
					if (HBoverlap < OVERLAPCUTOFF_3) {
						angs.push_front(clampAngle(ang - 1.5*ROUGH_STEP));
						angs.push_front(clampAngle(ang + 1.5*ROUGH_STEP));
					}
					if (HBoverlap < OVERLAPCUTOFF_5) {
						angs.push_front(clampAngle(ang - 3.0*ROUGH_STEP));
						angs.push_front(clampAngle(ang + 3.0*ROUGH_STEP));
					}
				}

			}
		}

		// -------------------------------------------------------
		// next, look for least-bumping angle not near an acceptor

		double leastBumpAngle = START_ANGLE;
		double bestScore = LowestMoverScore;

		for (i=0; i < (360/ROUGH_STEP); i++) {
			const double oldTheta = angle();
			const double theta = START_ANGLE +
				((i+1)>>1) * ((i&1) ? 1.0 : -1.0) * ROUGH_STEP;
			bool tryThisAngle = TRUE;
			for(std::list<double>::iterator chi = angs.begin(); chi != angs.end(); ++chi) {
				const double delta = abs(clampAngle(clampAngle(theta) - (*chi)));
				if (delta < ROUGH_STEP) {
					tryThisAngle = FALSE;
					break;
				}
			}

			if (tryThisAngle) {
				for(std::list<PDBrec*>::const_iterator hydlist = _rot.begin(); hydlist != _rot.end(); ++hydlist) {
					setHydAngle(**hydlist, oldTheta, theta, xyz);
				}
				angle(theta);
				const double bs = bumpsThisAngle(xyz, dotBucket);
				if (bs > bestScore) { bestScore = bs; leastBumpAngle = theta; }
			}
		}

		// re-position to least bump angle
		for(std::list<PDBrec*>::const_iterator hydlist = _rot.begin(); hydlist != _rot.end(); ++hydlist) {
			setHydAngle(**hydlist, angle(), leastBumpAngle, xyz);
		}
		angle(leastBumpAngle);
		angs.push_front(clampAngle(leastBumpAngle));

		// -------------------------------------------------------
		// merge both sets of angles into the orientation angle array

		_nori = angs.size();
		_oang.resize(_nori);
		int j = 0;
		angs.sort();
		angs.reverse();
		for(std::list<double>::iterator chi = angs.begin(); chi != angs.end(); ++chi) {
			if (j > 0) {
				const double delta = abs(_oang[j-1] - (*chi));
				if ((delta >  (0.5*ROUGH_STEP))
					&& (*chi  >=((0.5*ROUGH_STEP)-180.0)) ) {
					_oang[j++] = *chi;
				}
			}
			else { _oang[j++] = *chi; }
		}
		_nori = j; // <-- the actual length after the merge operation
		_oang.resize(_nori);
	}
}

int RotDonor::makebumpers(NeighborList<BumperPoint*>& sym_bblks, 
						  int rn, float& maxVDWrad) {
	int an = 0;
	const double dtheta = 10.0; // fineness of rotation angle scan
	const double scanAngle = 180.0;
	BumperPoint* bp;
	for(std::list<PDBrec*>::const_iterator it = _rot.begin(); it != _rot.end(); ++it) {
		PDBrec* a = *it;
		for (double theta = -scanAngle; theta < scanAngle; theta += dtheta) {
			Point3d p(a->loc().rotate(theta, _p2, _p1));
			bp = new BumperPoint(p, rn, an++, a->vdwRad());
			sym_bblks.insert(bp);
			if (a->vdwRad() > maxVDWrad) { maxVDWrad = a->vdwRad(); }
		}
	}
	return an;
}
