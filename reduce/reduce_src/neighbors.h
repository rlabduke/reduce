// name: neighbors.h
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

#ifndef NEIGHBORS_H
#define NEIGHBORS_H 1

#include "LocBlk.h"

// requires class L to answer the message loc() with a const Point3d&
// and is copyable

// returns a sequence of Ls within range of p from locdict

template <class L>
std::list<L> neighbors(const Point3d& p, Coord minrange, Coord maxrange,
		     const std::multimap<LocBlk, L>& locdict);

#ifdef INCTEMPLATEDEFNS
#include "neighbors.cpp"
#endif

#endif
