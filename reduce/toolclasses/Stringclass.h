// Stringclass.h
// Author: J. Michael Word
// Date Written: 4/92
// Purpose: Class for Stringclass manipulaton

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#ifndef __CLASS_STRINGCLASS_H__
#define __CLASS_STRINGCLASS_H__

#include <iostream>
#ifdef OLD_STD_HDRS
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#else
#include <cstdio>
#include <cstring>
#include <cctype>
using std::strcpy;
using std::strlen;
using std::strcat;
using std::strcmp;
using std::strchr;
using std::strstr;
using std::strrchr;
using std::strncpy;
using std::sprintf;
using std::sscanf;
using std::tolower;
using std::toupper;
#endif

#ifndef BOOLPREDEFINED
typedef int bool;
#endif

class Stringclass {
private:
   // internal representation of a character array
   // which allows them to be shared

   class StringRep {
   friend class Stringclass;
   public:
      StringRep(const char *s) {
         ::strcpy(rep = new char[::strlen(s)+1], s);
         count = 1;
      }
      ~StringRep() { delete[] rep; }
   private:
      // special StringRep constructor used when
      // Stringclass makes the buffer
      StringRep(char ** const r) { rep = *r; *r = 0; count = 1; }

      char *rep;
      int count;
   };

   // special Stringclass constructor used to optimize concatination
   Stringclass(char ** const r) { rep = new StringRep(r); }

public:
   // make...
   Stringclass() { rep = new StringRep(""); }
   Stringclass(const Stringclass& s) { rep = s.rep; rep->count++; }
   Stringclass(const char *s) { rep = new StringRep(s); }
   Stringclass(const char);
   Stringclass(const int);
   Stringclass(const double);
   Stringclass(const int, const char *fmt);
   Stringclass(const double, const char *fmt);
   // break...
   ~Stringclass() { if (--rep->count <= 0) delete rep; }

   // feature access
   bool empty() const { return (rep->rep)[0] == '\0'; }
   int length() const { return ::strlen(rep->rep); }
   const char *array() const { return rep->rep; }
   int asint(const char *fmt="%d") const;
   float asfloat(const char *fmt="%g") const;

   // conversion
   operator const char *() const { return rep->rep; }

   Stringclass asLowercase();
   Stringclass asUppercase();

   // shallow copy
   Stringclass& operator=(const Stringclass& s) {
      s.rep->count++;
      if (--rep->count <= 0) delete rep;
      rep = s.rep;
      return *this;
   }
   // deep copy
   Stringclass duplicate() const { Stringclass retval(rep->rep); return retval; }

   // equivalence
   int operator==(const Stringclass& s) const {
      return operator==(s.rep->rep);
   }
   int operator==(const char *s) const {
      return ::strcmp(rep->rep, s) == 0;
   }
   int operator!=(const Stringclass& s) const {
      return operator!=(s.rep->rep);
   }
   int operator!=(const char *s) const {
      return ::strcmp(rep->rep, s) != 0;
   }

   // relationship
   int operator<(const Stringclass& s) const {
      return operator<(s.rep->rep);
   }
   int operator<(const char *s) const {
      return ::strcmp(rep->rep, s) < 0;
   }
   int operator>(const Stringclass& s) const {
      return operator>(s.rep->rep);
   }
   int operator>(const char *s) const {
      return ::strcmp(rep->rep, s) > 0;
   }

   // element access
   char operator[](const int i) const { return rep->rep[i]; };

   int operator[](const char ch) const { // find character
      const char *p = ::strchr(rep->rep, ch);
      return (p == NULL) ? (-1) : (p - rep->rep);
   };

   int operator[](const char *s) const { // find substring
      const char *p = ::strstr(rep->rep, s);
      return (p == NULL) ? (-1) : (p - rep->rep);
   };
   int operator[](const Stringclass& s) const {
      return operator[](s.rep->rep);
   };
   int contains(const char *s) const { // is s a substring?
      return(::strstr(rep->rep, s) != NULL);
   };
   int contains(const Stringclass& s) const {
      return(::strstr(rep->rep, s.rep->rep) != NULL);
   };

   // substring
   Stringclass operator()(int, int=-1) const;
   Stringclass afterLast(const char ch) const;

   // concatenation
   Stringclass operator+(const Stringclass&) const;
   friend Stringclass operator+(const char*, const Stringclass&);
   friend Stringclass operator+(const Stringclass&, const char*);

   int last(const char ch) const { // find last character
      const char *p = ::strrchr(rep->rep, ch);
      return (p == NULL) ? (-1) : (p - rep->rep);
   };

private:
   StringRep *rep;
};

std::ostream& operator<<(std::ostream&, const Stringclass&);
Stringclass operator+(const Stringclass& S, const char* s);
Stringclass operator+(const char* s, const Stringclass& S);
#endif
