/* name: utility.C              */
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

#include <stdio.h>
#include <string>
#include <ctype.h>
#include <stdlib.h>
#include "utility.h"
#include <string.h>

void note(const char *message) {
	fprintf(stderr, "%s\n", message);
}

void warn(const char *message) {
	fprintf(stderr, "WARNING: %s\n", message);
}

void errmsg(const char *message) {
	fprintf(stderr, "ERROR: %s\n", message);
}

void halt(const char *message) {
	fprintf(stderr, "ERROR: %s\n", message);
	exit(1);
}

int trimStr(char *s) {
   int nrem = 0;

   for (int i = strlen(s); i > 0 && isspace(s[i-1]); i--) {
      s[i-1] = '\0'; nrem++;
   }
   return nrem;
}

void copyChars(char *to, const char *from, int n) {
	int i;

	for(i=0; i<n; i++) { to[i] = from[i]; }
}

#ifdef NEEDSTRCASECMP
int strncasecmp(const char *buf, const char *pat, int sz) {
	int rc = 0;
	for(int i=0; i < sz; i++) {
		if (tolower(buf[i]) != tolower(pat[i])) { rc = 1; break; }
		else if (buf[i] == '\0') { break; }
	}
	return rc;
}
int strcasecmp(const char *buf, const char *pat) {
	int rc = 0;
	for(int i=0; buf[i] && pat[i]; i++) {
		if (tolower(buf[i]) != tolower(pat[i])) { rc = 1; break; }
	}
	return rc;
}
#endif

double clampAngle(double a, int min) {
   if (min == 0) {
      while (a >= 360.0) { a -= 360.0; }
      while (a <    0.0) { a += 360.0; }
   }
   else {
      while (a > 360.0+min) { a -= 360.0; }
      while (a <=  0.0+min) { a += 360.0; }
   }
   return a;
}

int compArgStr(const char *str, const char *arg, int min) {
	int i, max;
	char s, a;

	if (!str || !arg) return 0;

	max = strlen(arg);

	for(i=0; i<max; i++) {
		s = toupper(str[i]);
		a = toupper(arg[i]);

		if (i >= min && (s == '\0' || s == '.' || isdigit(s))) {
			break; /* good ending point */
		}
		else if (s != a) {
			i = 0; /* failed to match */
			break;
		}
	}

	return i;
}

int parseInteger(const char *str, int start, int len) {
	register int value = 0;
	register char ch;
	int neg = 0, inside = 0;

	if (!str || start < 0) { return 0; }
	str += start;

	while((len-- > 0) && *str) {
		ch = *str++;
		if ((ch >='0') && (ch <= '9')) {
			value = (10*value) + (ch - '0');
			inside = 1;
		}
		else if (ch == '+' && !inside) {
			inside = 1;
		}
		else if (ch == '-' && !inside) {
			neg = 1;
			inside = 1;
		}
		else if (isspace(ch) && !inside) { /* nothing */ }
		else break; /* end of integer */
	}
	return (neg?-value:value);
}

double parseReal(const char *str, int start, int len) {
   double value = 0.0, scale = 1.0, expscale = 1.0, expfact = 10.0;
   int expval = 0;
   register char ch;
   int inside = 0, infract = 0, inexp = 0, insn = 0, esn = 0;

   if (!str || start < 0) { return 0; }
   str += start;

   while((len-- > 0) && *str) {
      ch = *str++;
      if (inexp) {
	 if ((ch >='0') && (ch <= '9')) {
	    expval = (10*expval) + (ch - '0');
	    esn = 1;
	 }
	 else if (ch == '+' && !esn) {
	    esn = 1;
	 }
	 else if (ch == '-' && !esn) {
	    expfact = 0.1;
	    esn = 1;
	 }
	 else break; /* end of real */
      }
      else if ((ch >='0') && (ch <= '9')) {
	 value = (10.0*value) + (ch - '0');
	 if (infract) { scale *= 0.1; }
	 inside = 1;
      }
      else if (ch == '+' && !inside && !insn) {
	 insn = 1;
      }
      else if (ch == '-' && !inside && !insn) {
	 scale = -1.0;
	 insn = 1;
      }
      else if (ch == '.' && !infract) {
	 inside = infract = 1;
      }
      else if ((ch == 'e' || ch == 'E') && !inexp) {
	 if (!inside) { value = 1.0; }
	 inexp = inside = 1;
      }
      else if (isspace(ch) && !inside && !insn) { /* nothing */ }
      else break; /* end of real */
   }
   if (expval) {
      for(;expval; expval--) { expscale *= expfact; }
   }
   return value*scale*expscale;
}

std::string
toUppercase(const char* a) {
	int i = 0;
	while (a[i++] != '\0')
		;
	const int len = i;
	char *buf = new char[len + 1];
	for (i = 0; i < len; i++) {
#ifdef CHARFUNCMACROS
		buf[i] = toupper(a[i]);
#else
		buf[i] = ::toupper(a[i]);
#endif
	}
	buf[i] = '\0';
	std::string retval( buf );
	delete [] buf;
	return retval;
}
