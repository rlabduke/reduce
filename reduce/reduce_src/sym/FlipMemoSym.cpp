// name: FlipMemoSym.cpp
// author: J. Michael Word
// date written: 2/7/98
// purpose: Symmetry-aware implementation for FlipMemo

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
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#else
#include <cctype>
#include <cmath>
#include <cstdlib>
using std::exit;
using std::toupper;
#endif

#include "../FlipMemo.h"
#include "../AtomConn.h"
#include "../AtomPositions.h"

extern bool UseNuclearDistances; //defined in reduce.cpp JJH

int FlipMemo::makebumpers(NeighborList<BumperPoint*>& sym_bblks, int rn, float& maxVDWrad) {
	int i = 0, an = 0;
	BumperPoint* bp;
	if (_isComplete) {
		for (i = 0; i < _resFlip[_resType].numBmpr; i++) { // regular
			const int f1 = i + 1;
			bp = new BumperPoint(_origLoc[f1], rn, an++, _wrkAtom[f1].vdwRad());
			sym_bblks.insert(bp);

			if (_wrkAtom[f1].vdwRad() > maxVDWrad) { maxVDWrad = _wrkAtom[f1].vdwRad(); }
		}
		for (i = 0; i < _resFlip[_resType].numPP; i++) { // flipped
			const int f2 = _resFlip[_resType].numPnts - i;
			bp = new BumperPoint(_origLoc[f2], rn, an++, _wrkAtom[f2].vdwRad());
			sym_bblks.insert(bp);

			if (_wrkAtom[f2].vdwRad() > maxVDWrad) { maxVDWrad = _wrkAtom[f2].vdwRad(); }
		}
	}
	return an;
}
