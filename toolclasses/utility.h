/* name: utility.h              */
/* author: J. Michael Word      */
/* date written: 2/26/96        */
/* purpose: utility functions   */

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#ifndef UTILITY_H
#define UTILITY_H 1

#include <string>
#include <algorithm>

#ifndef BOOLPREDEFINED
typedef int bool;
#endif

#ifndef PIPREDEFINED
const double PI = 3.14159265358979323846264;
#endif

#if defined(__DECCXX_VER)
#include <fstream>
#elif defined(__APPLE_CC__) && __APPLE_CC__ != 1 && __APPLE_CC__ <= 1671
// pass
#else
#ifndef TFPREDEFINED
const bool TRUE  = 1;
const bool FALSE = 0;
#endif
#endif

#ifndef ABSPREDEFINED
template <class T>
inline T abs(const T& x) {
   return (x < 0) ? -x : x;
}
#endif
#ifndef MIN3PREDEFINED
template <class T>
inline T min(const T& a, const T& b, const T& c) {
   return (a < b) ? ((a < c) ? a : c) : ((b < c) ? b : c);
}
#endif
#ifndef MAX3PREDEFINED
template <class T>
inline T max(const T& a, const T& b, const T& c) {
   return (a < b) ? ((b < c) ? c : b) : ((a < c) ? c : b);
}
#endif
#ifndef MEDIAN3PREDEFINED
template <class T>
inline const T& median(const T& a, const T& b, const T& c) {
   if (a < b) {
           if (b < c) { return b; }
      else if (a < c) { return c; }
      else            { return a; }
   }
   else if (a < c) { return a; }
   else if (b < c) { return c; }
   else            { return b; }
}
#endif
#ifndef SWAPPREDEFINED
// swaps two objects
template <class T>
inline void swap2(T& a, T& b) { T temp = a; a = b; b = temp; }
#endif

void note(const char *message);
void warn(const char *message);
void errmsg(const char *message);
void halt(const char *message);

int trimStr(char *str);

void copyChars(char *to, const char *from, int n);

#if defined(__DECCXX_VER) || defined(_MSC_VER)
#define NEEDSTRCASECMP
#endif
#ifdef NEEDSTRCASECMP
int strncasecmp(const char *buf, const char *pat, int sz);
int strcasecmp(const char *buf, const char *pat);
#endif

double clampAngle(double a, int min=-180);

int compArgStr(const char *str, const char *arg, int n);

int parseInteger(const char *str, int start, int len);
double parseReal(const char *str, int start, int len);

// column/fixed field formatted io routines from pdb++ classes
int column_sscanf(const char *, const char *, ...);
int column_sprintf(char *, const char *, ...);

std::string toUppercase(const char* a);

struct DeleteObject {
	template<class T>
	void operator() (const T* ptr) const {
		delete ptr;
	}
};

#endif

