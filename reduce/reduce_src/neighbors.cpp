// name: neighbors.C
// author: J. Michael Word
// date written: 8/21/97
// purpose: generic 3d range searching in a MultiDict

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#include "neighbors.h"

#include <map>

// requires class L to be copyable and to answer
// the message loc() with a Point3d

template <class L>
std::list<L> neighbors(const Point3d& p, Coord minrange, Coord maxrange,
					   const std::multimap<LocBlk, L>& locdict) {
	const Coord span = LocBlk::chunksz(1.0);

	const Point3d r(maxrange, maxrange, maxrange);

	const Coord minsq = minrange * minrange;
	const Coord maxsq = maxrange * maxrange;

	Point3d qm, qp;
	LocBlk::chunk(qm, p-r);
	LocBlk::chunk(qp, p+r);

	const Coord x0 = LocBlk::chunksz(qm.x());
	const Coord y0 = LocBlk::chunksz(qm.y());
	const Coord z0 = LocBlk::chunksz(qm.z());
	const Coord x1 = LocBlk::chunksz(qp.x()) + LocBlkTolerance;
	const Coord y1 = LocBlk::chunksz(qp.y()) + LocBlkTolerance;
	const Coord z1 = LocBlk::chunksz(qp.z()) + LocBlkTolerance;

	typename std::list< L > inrange;
	typename std::multimap<LocBlk, L >::const_iterator it, it2;
	
	for (Coord qx=x0; qx <= x1; qx+=span) {
		for (Coord qy=y0; qy <= y1; qy+=span) {
			for (Coord qz=z0; qz <= z1; qz+=span) {
				LocBlk k(qx + LocBlkTolerance,
					qy + LocBlkTolerance,
					qz + LocBlkTolerance);

				it2 = locdict.upper_bound(k);
				for (it = locdict.lower_bound(k); it != it2; ++it) {
					const Coord dsq = distanceSquared((it->second)->loc(), p);
					if (dsq >= minsq && dsq <= maxsq) {
						inrange.push_front(it->second);
					}
				}
			} // for z
		} // for y
	} // for x

	return inrange;
}
