// Stringclass.C
// Author: J. Michael Word
// Date Written: 4/92
// Purpose: Methods for Stringclass manipulaton

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#include "Stringclass.h"

#ifdef OLD_STD_HDRS
#include <ctype.h>
#include <stdio.h>
#else
#include <cctype>
#include <cstdio>
#endif

// size of buffers used when converting numbers to strings
const int NumBufMaxLen = 20;

// make a Stringclass from a char
Stringclass::Stringclass(const char ch) {
   char b[2];
   b[0] = ch; b[1] = '\0';
   rep = new StringRep(b);
}

// make a Stringclass from an int
Stringclass::Stringclass(const int n) {
   char b[NumBufMaxLen];
   ::sprintf(b, "%d", n);
   rep = new StringRep(b);
}

// make a Stringclass from an int
Stringclass::Stringclass(const int n, const char *fmt) {
   char b[NumBufMaxLen];
   ::sprintf(b, fmt, n);
   rep = new StringRep(b);
}

// make a Stringclass from an double
Stringclass::Stringclass(const double f) {
   char b[NumBufMaxLen];
   ::sprintf(b, "%g", f);
   rep = new StringRep(b);
}

// make a Stringclass from an double
Stringclass::Stringclass(const double f, const char *fmt) {
   char b[NumBufMaxLen];
   ::sprintf(b, fmt, f);
   rep = new StringRep(b);
}

// convert string to an integer -- format has default: %d
int Stringclass::asint(const char *fmt) const {
   int retval = 0;
   if (::sscanf(rep->rep, fmt, &retval) != 1) return 0;
   else return retval;
}

// convert string to a float -- format has default: %g
float Stringclass::asfloat(const char *fmt) const {
   float retval = 0.0;
   if (::sscanf(rep->rep, fmt, &retval) != 1) return 0.0;
   else return retval;
}

// generate a new lower case Stringclass
Stringclass Stringclass::asLowercase() {
   int i;
   const int len = length();
   char *buf = new char[len + 1];
   for (i = 0; i < len; i++) {
#ifdef CHARFUNCMACROS
      buf[i] = tolower(rep->rep[i]);
#else
      buf[i] = ::tolower(rep->rep[i]);
#endif
   }
   buf[i] = '\0';
   Stringclass retval( &buf );
   return retval;
}

// generate a new upper case Stringclass
Stringclass Stringclass::asUppercase() {
   int i;
   const int len = length();
   char *buf = new char[len + 1];
   for (i = 0; i < len; i++) {
#ifdef CHARFUNCMACROS
      buf[i] = toupper(rep->rep[i]);
#else
      buf[i] = ::toupper(rep->rep[i]);
#endif
   }
   buf[i] = '\0';
   Stringclass retval( &buf );
   return retval;
}

// substring: Stringclass s("abcdef"); s(3,2) == Stringclass("de")
Stringclass Stringclass::operator()(int fr, int sz) const {
   const int len = length();
   if (fr >= len || fr < 0) { fr = sz = 0; }
   // neg size means remainder
   if ((fr+sz) > len || sz < 0) { sz = length() - fr; }

   char *buf = new char[sz + 1];
   if (sz) {
      ::strncpy(buf, rep->rep + fr, sz);
      buf[sz] = '\0';
   }
   else { buf[0] = '\0'; }

   Stringclass retval( &buf );
   return retval;
}

// substring: return portion after the last character ch
Stringclass Stringclass::afterLast(const char ch) const {
   const char *p = ::strrchr(rep->rep, ch);
   return Stringclass((p == NULL) ? "" : p+1);
}

// combine:  Str + "abc"
Stringclass operator+(const Stringclass& S, const char* s) {
   char *buf = new char[S.length() + ::strlen(s) + 1];
   ::strcpy(buf, (const char *)S);
   ::strcat(buf, s);
   Stringclass retval( &buf );
   return retval;
}

// combine:  "abc" + Str
Stringclass operator+(const char* s, const Stringclass& S) {
   char *buf = new char[::strlen(s) + S.length() + 1];
   ::strcpy(buf, s);
   ::strcat(buf, (const char *)S);
   Stringclass retval( &buf );
   return retval;
}

// combine:  StrA + StrB
Stringclass Stringclass::operator+(const Stringclass& s) const {
   char *buf = new char[s.length() + length() + 1];
   ::strcpy(buf, rep->rep);
   ::strcat(buf, s.rep->rep);
   Stringclass retval( &buf );
   return retval;
}

// stream output
std::ostream& operator<<(std::ostream& os, const Stringclass& s) {
   os << (const char *)s;
   return os;
}
