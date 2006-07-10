//
//	Copyright (c) 1989,1992 The Regents of the University of California.
//	All rights reserved.
//
//	Redistribution and use in source and binary forms are permitted
//	provided that the above copyright notice and this paragraph are
//	duplicated in all such forms and that any documentation,
//	advertising materials, and other materials related to such
//	distribution and use acknowledge that the software was developed
//	by the University of California, San Francisco.  The name of the
//	University may not be used to endorse or promote products derived
//	from this software without specific prior written permission.
//	THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
//	IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
//	WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//

#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include "utility.h"

	// scratch must be big enough to hold the largest number
static char	scratch[256];

#define	OVERFLOW_CHAR	'*'

static char	*outint(int, int, int, char, char, int, char *, char);
static char	*outunsigned(unsigned int, int, char, int, char *);
static char	*outstr(char *, int, int, char, int, char *);
static char	*outfloat(double, int, int, char, int, char *);
static char	*outexp(double, int, int, char, int, char *);
static char	*e_out(int, char *);

int column_sprintf(char *outbuf, const char *fmt, ...) {
	va_list		argv;
	char		*p;
	const char	*f;
	int		field1, field2;
	char		c, fill_char;
	int		inum;
	unsigned 	unum;
	double		fnum;
	int		left_justify;

	va_start(argv, fmt);
	f = fmt;
	p = outbuf;
	while (*f) {
		if (*f == '%') {
			f++;
			if (*f == '-')
				left_justify = 1, f++;
			else
				left_justify = 0;

			if (*f == '0')
				fill_char = '0', f++;
			else
				fill_char = ' ';

			if (isdigit(*f)) {
				field1 = *f++ - '0';
				while (isdigit(*f))
					field1 = field1 * 10 + *f++ - '0';
			}
			else
				field1 = -1;

			if (*f == '.') {
				f++;
				field2 = 0;
				while (isdigit(*f))
					field2 = field2 * 10 + *f++ - '0';
			}
			else
				field2 = -1;

			if (*f == 'l' || *f == 'h')
				f++;

			while (isspace(*f))
				f++;
			switch (*f) {
			  case 'c':
				c = (char) va_arg(argv, int);
				if (c == '\0')
					c = ' ';
				if (left_justify)
					*p++ = c;
				while (--field1 > 0)
					*p++ = fill_char;
				if (!left_justify)
					*p++ = c;
				break;
			  case 'd':
			  case 'D':
				inum = va_arg(argv, int);
				p = outint(inum, field1, 10, fill_char, 'a',
					left_justify, p, (*f == 'D') ? ' ':'0');
				break;
			  case 'e':
				fnum = va_arg(argv, double);
				if (field2 < 0)
					field2 = 6;
				p = outexp(fnum, field1, field2, fill_char,
					left_justify, p);
				break;
			  case 'f':
				fnum = va_arg(argv, double);
				if (field2 < 0)
					field2 = 6;
				p = outfloat(fnum, field1, field2, fill_char,
					left_justify, p);
				break;
			  case 'o':
				inum = va_arg(argv, int);
				p = outint(inum, field1, 8, fill_char, 'a',
					left_justify, p, '0');
				break;
			  case 's':
				p = outstr(va_arg(argv, char *), field1, field2,
					fill_char, left_justify, p);
				break;
			  case 'u':
				unum = va_arg(argv, unsigned);
				p = outunsigned(unum, field1, fill_char,
					left_justify, p);
				break;
			  case 'x':
				inum = va_arg(argv, int);
				p = outint(inum, field1, 16, fill_char, 'a',
					left_justify, p, '0');
				break;
			  case 'X':
				inum = va_arg(argv, int);
				p = outint(inum, field1, 16, fill_char, 'A',
					left_justify, p, '0');
				break;
			  default:
				if (left_justify)
					*p++ = *f;
				while (--field1 > 0)
					*p++ = fill_char;
				if (!left_justify)
					*p++ = *f;
				break;
			}
			f++;
		}
		else if (*f == '\\') {		/* Special character */
			switch (*++f) {
			  case 'n':
				*p++ = '\n';
				break;
			  case 'r':
				*p++ = '\r';
				break;
			  case 'b':
				*p++ = '\b';
				break;
			  case 't':
				*p++ = '\t';
				break;
			  case 'f':
				*p++ = '\f';
				break;
			  case '0': case '1': case '2': case '3':
			  case '4': case '5': case '6': case '7':
				inum = *f++ - '0';
				if (*f >= '0' && *f <= '7') {
					inum = inum * 8 + *f++ - '0';
					if (*f >= '0' && *f <= '7')
						inum = inum * 8 + *f++ - '0';
				}
				f--;
				*p++ = (char) inum;
				break;
			  default:
				*p++ = *f;
			}
			f++;
		}
		else				/* Normal character */
			*p++ = *f++;
	}
	*p = '\0';
	va_end(argv);
	return p - outbuf;
}

static char *
e_out(int width, char *where)
{
	while (width-- > 0)
		*where++ = OVERFLOW_CHAR;
	return where;
}

static char *
outint(int value, int width, int radix, char fill_char, char hex,
					int left_justify, char *p, char zero)
{
	char	*s;
	int	n;
	int	negative;

	if (value < 0)
		negative = 1, value = -value, width--;
	else
		negative = 0;
	s = scratch;
	if (value)
		do {
			n = value % radix;
			*s++ = n < 10 ? '0' + n : hex + n - 10;
			value /= radix;
		} while (value);
	else
		*s++ = zero;
	n = s - scratch;
	if (width != -1 && n > width)
		return e_out(width + negative, p);

	if (negative && fill_char == '0')
		*p++ = '-';
	if (!left_justify)
		while (width-- > n)
			*p++ = fill_char;
	if (negative && fill_char == ' ')
		*p++ = '-';
	while (--s >= scratch)
		*p++ = *s;
	if (left_justify)
		while (width-- > n)
			*p++ = fill_char;
	return p;
}

static char *
outunsigned(unsigned int value, int width, char fill_char, int left_justify,
									char *p)
{
	char	*s;
	int	n;

	s = scratch;
	while (value) {
		*s++ = value % 10 + '0';
		value /= 10;
	}
	n = s - scratch;
	if (n == 0)
		*s++ = '0', n = 1;
	if (width != -1 && n > width)
		return e_out(width, p);

	if (!left_justify)
		while (width-- > n)
			*p++ = fill_char;
	while (--s >= scratch)
		*p++ = *s;
	if (left_justify)
		while (width-- > n)
			*p++ = fill_char;
	return p;
}

static char *
outstr(char *s, int width, int maxstr, char fill_char, int left_justify, char *p)
{
	int	len;

	len = strlen(s);
	if (maxstr >= 0 && len > maxstr)
		len = maxstr;
	if (width != -1 && len > width)
		return e_out(width, p);

	if (!left_justify)
		while (width-- > len)
			*p++ = fill_char;
	else
		width -= len;
	while (len--)
		*p++ = *s++;
	if (left_justify)
		while (width-- > 0)
			*p++ = fill_char;
	return p;
}

static char *
outfloat(double value, int width, int nplace, char fill_char, int left_justify,
									char *p)
{
	int	i, intval;
	char	*place, *to, *from;
	int	negative;

	negative = value < 0.0 ? 1 : 0;
		
	if (negative)
		value = -value;

	for (i = 0; i < nplace; i++)
		value *= 10.0;

	intval = (int) (value + 0.5);

	if (width == -1)
		width = nplace + 4;		/* TODO: fix */
	else if (nplace + (nplace == 0 ? 1 : 2) > width)
		return e_out(width, p);

	for (place = p + width - 1; place >= p + width - nplace; place--) {
		*place = '0' + intval % 10;
		intval /= 10;
	}

	if (nplace > 0)
		*place-- = '.';

	if (intval == 0)
		*place-- = '0';

	for (; place >= p; place--) {
		if (intval == 0)
			break;
		else {
			*place = '0' + intval % 10;
			intval /= 10;
		}
	}

	if (intval != 0)
		return e_out(width, p);

	if (place < p && negative)
		return e_out(width, p);

	if (left_justify) {
		for (from = place + 1, to = (negative ? p + 1 : p);
						from < p + width; from++, to++)
			*to = *from;
		for (; to < p + width; to++)
			*to = fill_char;
		if (negative)
			*p = '-';
	} else {
		for (to = place; to >= p; to--)
			*to = fill_char;
		if (negative)
			if (fill_char == ' ')
				*place = '-';
			else
				*p = '-';
	}

	return p + width;
}

static char *
outexp(double value, int width, int nplace, char fill_char, int left_justify,
									char *p)
{
	int	n;
	char	*s;
	int	negative;
	double	fraction;

	if (value < 0)
		negative = 1, value = -value, width--;
	else
		negative = 0;

	n = 0;
	while (value > 10)
		n++, value /= 10;
	if (value)
		while (value < 1)
			n--, value *= 10;

	s = scratch;
	if (n < 0) {
		n = -n;
		*s++ = n % 10 + '0';
		*s++ = n / 10 + '0';
		*s++ = '-';
	}
	else {
		*s++ = n % 10 + '0';
		*s++ = n / 10 + '0';
		*s++ = '+';
	}
	*s = 'e';

	s = scratch + nplace + 4;	/* 4 == strlen("e+00") */
	fraction = value - (int) value;
	for (n = 0; n < nplace; n++) {
		fraction *= 10.0;
		*--s = '0' + (int) fraction;
		fraction -= (int) fraction;
	}

	s = scratch + nplace + 4;
	if (nplace)
		*s++ = '.';
	n = (int) value;
	if (n)
		*s++ = n % 10 + '0';
	else
		*s++ = '0';
	n = s - scratch;
	if (width != -1 && n > width)
		return e_out(width + negative, p);

	if (negative && fill_char == '0')
		*p++ = '-';
	if (!left_justify)
		while (width-- > n)
			*p++ = fill_char;
	if (negative && fill_char == ' ')
		*p++ = '-';
	while (--s >= scratch)
		*p++ = *s;
	if (left_justify)
		while (width-- > n)
			*p++ = fill_char;
	return p;
}
