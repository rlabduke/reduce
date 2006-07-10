// name: Hash.h
// author: J. Michael Word
// date written: 8/7/97
// purpose: Interface file for a dictionary class
//          which maps (sequences of) copies of objects to keys.

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#ifndef HASH_H
#define HASH_H 1

unsigned long hash( const char *key, unsigned long M );
unsigned long hash( long key, unsigned long M );
unsigned long hash( char key, unsigned long M );
unsigned long genHashM( const unsigned long sz );

#endif
