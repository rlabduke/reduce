// name: Hash.C
// author: J. Michael Word
// date written: 8/7/97
// purpose: Implementation for Hash

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#include "Hash.h"

// -----------------------------------------------------
//  these are the main functions we use to hash our keys
// -----------------------------------------------------
#ifdef USEOLDSTRINGHASH
unsigned long hash(const char *key, unsigned long M) {
   unsigned long h = 0;

// this was what we were using initially

   for (int i = 0; key[i] != '\0'; i++) {
      h = (64*h + key[i]) % M;
   }
   return h;
}
#else
unsigned long hash(const char *key, unsigned long M) {
   unsigned long h = 0, g = 0;

// Hash function from PJ Weinberger's compiler as described 
// in Aho, Sethi & Ullman, Compilers, Principles, Techniques & Tools
// 1986, Addison Wesley, pg 436

   for (int i = 0; key[i] != '\0'; i++) {
      h = (h<<4) + key[i];
      if ((g = h & 0xf0000000)) {
	 h ^= g >> 24;
	 h ^= g;
      }
   }
   return h % M;
}
#endif

unsigned long hash(long key, unsigned long M) {
   unsigned long h = (key < 0) ? -key : key;

   return h % M;
}
unsigned long hash(char key, unsigned long M) {
   unsigned long h = key;

   return h % M;
}

// ----------------------------------------------------
//  return a prime number to be used as hash table size
// ----------------------------------------------------
unsigned long genHashM( const unsigned long sz ) {
   struct Plist { unsigned long prime; float fract; };
   static Plist plist[] = { {3, 1.0}, {7, 1.0}, {47, 1.0}, {137, 1.0},
  {257, 1.0}, {577, 1.0}, {1117, 1.0}, {1987, 1.1},
  {2287, 1.2}, {2423, 1.3}, {2617, 1.4},
  {2741, 1.5}, {3079, 1.6}, {4057, 1.7}, {6133, 1.8},
  {7919, 2.0}, {9733, 3.0}, {11657, 5.0}, {13499, 17.0},
  {29443, 10.0}, {63809, 15.0}, {81799, 20.0}, {104729, 99.0} };

   int max = (sizeof(plist)/sizeof(Plist));
   unsigned long M = plist[ max - 1 ].prime;

   for (int i=0; i < max; i++) {
      if ( sz < (plist[i].prime*plist[i].fract) ) {
         M = plist[ i ].prime;
         break;
      }
   }
   return M;
}
