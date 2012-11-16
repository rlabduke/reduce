//
// with Modifications by J. Michael Word, Duke University
//
//	Copyright (c) 1992 The Regents of the University of California.
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
//	subroutine for reading PDB format files
//

# include	"pdb++.h"
extern "C" {
#ifdef OLD_STD_HDRS
#include <ctype.h>
#include <string.h>
#else
#include <cctype>
#include <cstring>
using std::islower;
using std::toupper;
#endif
}
#if defined(__DECCXX_VER) || defined(_MSC_VER)
#define NEEDSTRCASECMP
#endif
#ifdef NEEDSTRCASECMP
int strncasecmp(const char *buf, const char *pat, int sz);
int strcasecmp(const char *buf, const char *pat);
#endif

int PDB::pdbrunInputVersion = PDB::PDBRUNVersion;	// just in case
int PDB::pdbrunOutputVersion = PDB::PDBRUNVersion;	// just in case

PDB::GfxType
PDB::getGfxType(const char *buf)
{
	switch (buf[0]) {
	case 'L': case 'l':
		if (strcasecmp(buf + 1, "INE-LOOP") == 0)
			return GFX_LINE_LOOP;
		if (strcasecmp(buf + 1, "INE-STRIP") == 0)
			return GFX_LINE_STRIP;
		if (strcasecmp(buf + 1, "INES") == 0)
			return GFX_LINES;
		break;
	case 'M': case 'm':
		if (strcasecmp(buf + 1, "ARKERS") == 0)
			return GFX_MARKERS;
		break;
	case 'P': case 'p':
		if (strcasecmp(buf + 1, "OINTS") == 0)
			return GFX_POINTS;
		if (strcasecmp(buf + 1, "OLYGON") == 0)
			return GFX_POLYGON;
		break;
	case 'Q': case 'q':
		if (strcasecmp(buf + 1, "UAD-STRIP") == 0)
			return GFX_QUAD_STRIP;
		if (strcasecmp(buf + 1, "UADS") == 0)
			return GFX_QUADS;
		break;
	case 'T': case 't':
		if (strcasecmp(buf + 1, "RIANGLE-FAN") == 0)
			return GFX_TRIANGLE_FAN;
		if (strcasecmp(buf + 1, "RIANGLE-STRIP") == 0)
			return GFX_TRIANGLE_STRIP;
		if (strcasecmp(buf + 1, "RIANGLES") == 0)
			return GFX_TRIANGLES;
		break;
	}
	return GFX_UNKNOWN;
}

const char *
PDB::gfxChars(GfxType gt)
{
	switch (gt) {
	default:		return "UNKNOWN";
	case GFX_POINTS:	return "POINTS";
	case GFX_MARKERS:	return "MARKERS";
	case GFX_LINES:		return "LINES";
	case GFX_LINE_STRIP:	return "LINE-STRIP";
	case GFX_LINE_LOOP:	return "LINE-LOOP";
	case GFX_TRIANGLES:	return "TRIANGLES";
	case GFX_TRIANGLE_STRIP:	return "TRIANGLE-STRIP";
	case GFX_TRIANGLE_FAN:	return "TRIANGLE-FAN";
	case GFX_QUADS:		return "QUADS";
	case GFX_QUAD_STRIP:	return "QUAD-STRIP";
	case GFX_POLYGON:	return "POLYGON";
	}
}

static PDB::RecordType
pdbrun5Type(const char *buf)
{
	switch (buf[0]) {
	case 'A': case 'a':
		if (strncasecmp(buf + 1, "NGLE ", 5) == 0)
			return PDB::USER_ANGLE;
		if (strncasecmp(buf + 1, "TPOS ", 5) == 0)
			return PDB::USER_ATPOS;
		break;
	case 'B': case 'b':
		if (strncasecmp(buf + 1, "GCOLOR ", 7) == 0)
			return PDB::USER_BGCOLOR;
		break;
	case 'C': case 'c':
		if (strncasecmp(buf + 1, "HAIN ", 5) == 0)
			return PDB::USER_CHAIN;
		if (strncasecmp(buf + 1, "NAME ", 5) == 0)
			return PDB::USER_CNAME;
		if (strncasecmp(buf + 1, "OLOR ", 5) == 0)
			return PDB::USER_COLOR;
		break;
	case 'D': case 'd':
		if (strncasecmp(buf + 1, "ISTANCE ", 8) == 0)
			return PDB::USER_DISTANCE;
		break;
	case 'E': case 'e':
		if (strncasecmp(buf + 1, "NDOBJ ", 6) == 0)
			return PDB::USER_ENDOBJ;
		if (strncasecmp(buf + 1, "YEPOS ", 6) == 0)
			return PDB::USER_EYEPOS;
		break;
	case 'F': case 'f':
		if (strncasecmp(buf + 1, "ILE ", 4) == 0)
			return PDB::USER_FILE;
		if (strncasecmp(buf + 1, "OCUS ", 5) == 0)
			return PDB::USER_FOCUS;
		break;
	case 'G': case 'g':
		if (strncasecmp(buf + 1, "FX ", 3) != 0)
			break;
		if (strncasecmp(buf + 4, "COLOR ", 6) == 0)
			return PDB::USER_GFX_COLOR;
		if (strncasecmp(buf + 4, "DRAW ", 5) == 0)
			return PDB::USER_GFX_DRAW;
		if (strncasecmp(buf + 4, "FONT ", 5) == 0)
			return PDB::USER_GFX_FONT;
		if (strncasecmp(buf + 4, "LABEL ", 6) == 0)
			return PDB::USER_GFX_LABEL;
		if (strncasecmp(buf + 4, "MARKER ", 7) == 0)
			return PDB::USER_GFX_MARKER;
		if (strncasecmp(buf + 4, "MOVE ", 5) == 0)
			return PDB::USER_GFX_MOVE;
		if (strncasecmp(buf + 4, "POINT ", 6) == 0)
			return PDB::USER_GFX_POINT;
		break;
	case 'O': case 'o':
		if (strncasecmp(buf + 1, "BJECT ", 6) == 0)
			return PDB::USER_OBJECT;
		break;
	case 'P': case 'p':
		if (strncasecmp(buf + 1, "DBRUN ", 6) == 0)
			return PDB::USER_PDBRUN;
		break;
	case 'R': case 'r':
		if (strncasecmp(buf + 1, "ADIUS ", 6) == 0)
			return PDB::USER_RADIUS;
		break;
	case 'V': case 'v':
		if (strncasecmp(buf + 1, "IEWPORT ", 8) == 0)
			return PDB::USER_VIEWPORT;
		break;
	case 'W': case 'w':
		if (strncasecmp(buf + 1, "INDOW ", 6) == 0)
			return PDB::USER_WINDOW;
		break;
	}
	return PDB::USER;
}

static PDB::RecordType
pdbrun6Type(const char *buf)
{
	switch (buf[0]) {
	case 'A': case 'a':
		if (strncasecmp(buf + 1, "NGLE ", 5) == 0)
			return PDB::USER_ANGLE;
		if (strncasecmp(buf + 1, "TPOS ", 5) == 0)
			return PDB::USER_ATPOS;
		break;
	case 'B': case 'b':
		if (strncasecmp(buf + 1, "GCOLOR ", 7) == 0)
			return PDB::USER_BGCOLOR;
		break;
	case 'C': case 'c':
		if (strncasecmp(buf + 1, "HAIN ", 5) == 0)
			return PDB::USER_CHAIN;
		if (strncasecmp(buf + 1, "NAME ", 5) == 0)
			return PDB::USER_CNAME;
		if (strncasecmp(buf + 1, "OLOR ", 5) == 0)
			return PDB::USER_COLOR;
		break;
	case 'D': case 'd':
		if (strncasecmp(buf + 1, "ISTANCE ", 8) == 0)
			return PDB::USER_DISTANCE;
		break;
	case 'E': case 'e':
		if (strncasecmp(buf + 1, "NDOBJ", 5) == 0
		&& (buf[6] == '\0' || buf[6] == '\n' || buf[6] == '\r' || buf[6] == ' '))
			return PDB::USER_ENDOBJ;
		if (strncasecmp(buf + 1, "YEPOS ", 6) == 0)
			return PDB::USER_EYEPOS;
		break;
	case 'F': case 'f':
		if (strncasecmp(buf + 1, "ILE ", 4) == 0)
			return PDB::USER_FILE;
		if (strncasecmp(buf + 1, "OCUS ", 5) == 0)
			return PDB::USER_FOCUS;
		break;
	case 'G': case 'g':
		if (buf[1] != 'F' || buf[2] != 'X' || buf[3] != ' ')
			break;
		switch (buf[4]) {
		case 'B': case 'b':
			if (strncasecmp(buf + 5, "EGIN ", 5) == 0)
				return PDB::USER_GFX_BEGIN;
			break;
		case 'C': case 'c':
			if (strncasecmp(buf + 5, "OLOR ", 5) == 0)
				return PDB::USER_GFX_COLOR;
			break;
		case 'E': case 'e':
			if (buf[5] == 'N' && buf[6] == 'D'
			&& (buf[7] == '\0' || buf[7] == '\n' || buf[7] == '\r' || buf[7] == ' '))
				return PDB::USER_GFX_END;
			break;
		case 'F': case 'f':
			if (strncasecmp(buf + 5, "ONT ", 4) == 0)
				return PDB::USER_GFX_FONT;
			break;
		case 'L': case 'l':
			if (strncasecmp(buf + 5, "ABEL ", 5) == 0)
				return PDB::USER_GFX_LABEL;
			break;
		case 'N': case 'n':
			if (strncasecmp(buf + 5, "ORMAL ", 6) == 0)
				return PDB::USER_GFX_NORMAL;
			break;
		case 'T': case 't':
			if (strncasecmp(buf + 5, "EXTPOS ", 7) == 0)
				return PDB::USER_GFX_TEXTPOS;
			break;
		case 'V': case 'v':
			if (strncasecmp(buf + 5, "ERTEX ", 6) == 0)
				return PDB::USER_GFX_VERTEX;
			break;
		}
		break;
	case 'M': case 'm':
		if (strncasecmp(buf + 1, "ARK ", 4) == 0)
			return PDB::USER_MARK;
		if (strncasecmp(buf + 1, "ARKNAME ", 6) == 0)
			return PDB::USER_MARKNAME;
		break;
	case 'O': case 'o':
		if (strncasecmp(buf + 1, "BJECT", 5) == 0
		&& (buf[6] == '\0' || buf[6] == '\n' || buf[6] == '\r' || buf[6] == ' '))
			return PDB::USER_OBJECT;
		break;
	case 'P': case 'p':
		if (strncasecmp(buf + 1, "DBRUN ", 6) == 0)
			return PDB::USER_PDBRUN;
		break;
	case 'R': case 'r':
		if (strncasecmp(buf + 1, "ADIUS ", 6) == 0)
			return PDB::USER_RADIUS;
		break;
	case 'V': case 'v':
		if (strncasecmp(buf + 1, "IEWPORT ", 6) == 0)
			return PDB::USER_VIEWPORT;
		break;
	case 'W': case 'w':
		if (strncasecmp(buf + 1, "INDOW ", 6) == 0)
			return PDB::USER_WINDOW;
		break;
	}
	return PDB::USER;
}

PDB::RecordType
PDB::getType(const char *buf)
{
	char	rt[4];		// PDB record type
	int	i;

	for (i = 0; buf[i] != '\0' && buf[i] != '\n' && buf[i] != '\r' && i < 4; i += 1) {
		if (islower(buf[i]))
			rt[i] = toupper(buf[i]);
		else
			rt[i] = buf[i];
	}
	if (i < 4)
		for (; i < 4; i += 1)
			rt[i] = ' ';

	switch (rt[0]) {

	case 'A':
		switch (rt[1]) {
		case 'G':
			if (rt[2] == 'R' && rt[3] == 'D')
				return AGRDES;
			if (rt[2] == 'G' && rt[3] == 'R')
				return AGGRGT;
			break;
		case 'N':
			if (rt[2] == 'I' && rt[3] == 'S')
				return ANISOU;
			break;
		case 'T':
			if (rt[2] == 'O' && rt[3] == 'M')
				return ATOM;
			break;
		case 'U':
			if (rt[2] == 'T' && rt[3] == 'H')
				return AUTHOR;
			break;
		}
		break;

	case 'C':
		switch (rt[1]) {
		case 'M':
			if (rt[2] == 'P' && rt[3] == 'D')
				return CMPDES;
			if (rt[2] == 'P' && rt[3] == 'O')
				return CMPONT;
			break;
		case 'O':
			if (rt[2] == 'M' && rt[3] == 'P')
				return COMPND;
			if (rt[2] == 'N' && rt[3] == 'E')
				return CONECT;
			break;
		case 'R':
			if (rt[2] == 'Y' && rt[3] == 'S')
				return CRYST1;
			break;
		}
		break;

	case 'E':
		switch (rt[1]) {
		case 'N':
			if (rt[2] == 'D' && rt[3] == ' ')
				return END;
			if (rt[2] == 'D' && rt[3] == 'M')
				return ENDMDL;
			break;
		case 'X':
			if (rt[2] == 'P' && rt[3] == 'D')
				return EXPDTA;
			break;
		}
		break;

	case 'F':
		switch (rt[1]) {
		case 'T':
			if (rt[2] == 'N' && rt[3] == 'O')
				return FTNOTE;
			break;
		case 'O':
			if (rt[2] == 'R' && rt[3] == 'M')
				return FORMUL;
			break;
		}
		break;

	case 'H':
		if (rt[1] != 'E')
			break;
		if (rt[2] == 'T' && rt[3] == 'A')
			return HETATM;
		if (rt[2] == 'A' && rt[3] == 'D')
			return HEADER;
		if (rt[2] == 'T' && rt[3] == ' ')
			return HET;
		if (rt[2] == 'L' && rt[3] == 'I')
			return HELIX;
		break;

	case 'J':
		if (rt[1] == 'R' && rt[2] == 'N' && rt[3] == 'L')
			return JRNL;
		break;

	case 'M':
		switch (rt[1]) {
		case 'A':
			if (rt[2] == 'S' && rt[3] == 'T')
				return MASTER;
			break;
		case 'O':
			if (rt[2] == 'D' && rt[3] == 'E')
				return MODEL;
			break;
		case 'T':
			if (rt[2] == 'R' && rt[3] == 'I')
				return MTRIX;
			if (rt[2] == 'X' && rt[3] == 'D')
				return MTXDES;
			break;
		}
		break;

	case 'O':
		switch (rt[1]) {
		case 'B':
			if (rt[2] == 'S' && rt[3] == 'L')
				return OBSLTE;
			break;
		case 'R':
			if (rt[2] == 'I' && rt[3] == 'G')
				return ORIGX;
			break;
		}
		break;

	case 'R':
		if (rt[1] != 'E')
			break;
		if (rt[2] == 'M' && rt[3] == 'A')
			return REMARK;
		if (rt[2] == 'V' && rt[3] == 'D')
			return REVDAT;
		break;

	case 'S':
		switch (rt[1]) {

		case 'C':
			if (rt[2] == 'A' && rt[3] == 'L')
				return SCALE;
			break;

		case 'E':
			if (rt[2] == 'Q' && rt[3] == 'R')
				return SEQRES;
			break;

		case 'H':
			if (rt[2] == 'E' && rt[3] == 'E')
				return SHEET;
			break;

		case 'I':
			if (rt[2] == 'T' && rt[3] == 'E')
				return SITE;
			if (rt[2] == 'G' && rt[3] == 'A')
				return SIGATM;
			if (rt[2] == 'G' && rt[3] == 'U')
				return SIGUIJ;
			break;

		case 'O':
			if (rt[2] == 'U' && rt[3] == 'R')
				return SOURCE;
			break;

		case 'P':
			if (rt[2] == 'R' && rt[3] == 'S')
				return SPRSDE;
			break;

		case 'S':
			if (rt[2] == 'B' && rt[3] == 'O')
				return SSBOND;
			break;

		case 'Y':
			if (rt[2] == 'M' && rt[3] == 'D')
				return SYMDES;
			if (rt[2] == 'M' && rt[3] == 'O')
				return SYMOP;
			break;
		}
		break;

	case 'T':
		switch (rt[1]) {
		case 'E':
			if (rt[2] == 'R' && rt[3] == ' ')
				return TER;
			break;
		case 'R':
			if (rt[2] == 'N' && rt[3] == 'S')
				return TRNSFM;
			break;
		case 'U':
			if (rt[2] == 'R' && rt[3] == 'N')
				return TURN;
			break;
		case 'V':
			if (rt[2] == 'E' && rt[3] == 'C')
				return TVECT;
			break;
		}
		break;

	case 'U':
		if (rt[1] == 'S' && rt[2] == 'E' && rt[3] == 'R')
			switch (pdbrunInputVersion) {
			case 1: case 2: case 3: case 4: case 5:
				return pdbrun5Type(buf + 6);
			case 6:
				return pdbrun6Type(buf + 6);
			default:
				if (strncasecmp(buf + 6, "PDBRUN ", 7) == 0)
					return USER_PDBRUN;
				return USER;
			}
		break;
	}
	return UNKNOWN;
}
