/*
 * with Modifications by J. Michael Word, Duke University
 *
 *	Copyright (c) 1993 The Regents of the University of California.
 *	All rights reserved.
 *
 *	Redistribution and use in source and binary forms are permitted
 *	provided that the above copyright notice and this paragraph are
 *	duplicated in all such forms and that any documentation,
 *	advertising materials, and other materials related to such
 *	distribution and use acknowledge that the software was developed
 *	by the University of California, San Francisco.  The name of the
 *	University may not be used to endorse or promote products derived
 *	from this software without specific prior written permission.
 *	THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 *	IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 *	WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
/* extra field included in ANISOU, ATOM and related records */
/* UNKNOWN */	"%-80.80s",
/* ANISOU */	"%6 %5d %4s%c%4s%c%4d%c %7d%7d%7d%7d%7d%7d%2 %4s%2s%2s %10s", /* SIGUIJ */   /**non-standard**/
/* ATOM */	"%6 %5d %4s%c%4s%c%4d%c%3 %8f%8f%8f%6R%6f%6 %4s%2s%2s %10s", /* HETATM, SIGATM */   /**non-standard**/
/* AUTHOR */	"%9 %c%60s%10s",		/* COMPND, EXPDTA, JRNL, SOURCE */
/* COMPND */	"%9 %c%60s%10s",			/* AUTHOR */
/* CONECT */	"%6 %5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%9 %10s",
/* CRYST1 */	"%6 %9f%9f%9f%7f%7f%7f %11s%4d%10s",
/* END */	"%70 %10s",
/* FORMUL */	"%8 %2d  %4s%2d%c%51s%10s",
/* FTNOTE */	"%7 %3d %59s%10s",	/*  REMARK, SYMDES, MTXDES, CMPDES, AGRDES */
/* HEADER */	"%10 %40s%9s  %c%4s%4 %10s",
/* HELIX */	"%7 %3d %3s %4s%c %4d%c %4s%c %4d%c%2d%30s%10s",
/* HET */	"%7 %4s %c%4d%c  %5d%5 %40s%10s",
/* HETATM */	"%6 %5d %4s%c%4s%c%4d%c%3 %8f%8f%8f%6f%6f%6 %4s%2s%2s %10s",	/* ATOM */   /**non-standard**/
/* JRNL */	"%9 %c%60s%10s",			/* AUTHOR */
/* MASTER */	"%10 %5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%10s",
/* MTRIX */	"%5 %d %3d%10f%10f%10f%5 %10f   %2d%10 %10s",
/* OBSLTE */	"%8 %2d %9s %4s%6 %4s %4s %4s %4s %4s %4s %4s %4s%10s",
/* ORIGX */	"%5 %d%4 %10f%10f%10f%5 %10f%15 %10s",	/* SCALE */
/* REMARK */	"%7 %3d %59s%10s",			/* FTNOTE */
/* REVDAT */	"%7 %3d%2d %9s %7s %c%7 %31s%10s",
/* SCALE */	"%5 %d%4 %10f%10f%10f%5 %10f%15 %10s",
/* SEQRES */	"%6 %4d %c %4d  %4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%4s%10s",
/* SHEET */	"%6 %4d %3s%2d %4s%c%4d%c %4s%c%4d%c%2d %4s%4s%c%4d%c %4s%4s%c%4d%c%10s",
/* SIGATM */	"%6 %5d %4s%c%4s%c%4d%c%3 %8f%8f%8f%6f%6f%6 %4s%2s%2s %10s",	/* ATOM */   /**non-standard**/
/* SIGUIJ */	"%6 %5d %4s%c%4s%c%4d%c %7d%7d%7d%7d%7d%7d%2 %4s%2s%2s %10s", /* ANISOU */   /**non-standard**/
/* SITE */	"%7 %3d %3s %2d %4s%c%4d%c %4s%c%4d%c %4s%c%4d%c %4s%c%4d%c%9 %10s",
/* SOURCE */	"%9 %c%60s%10s",			/* AUTHOR */
/* SPRSDE */	"%8 %2d %9s %4s%6 %4s %4s %4s %4s %4s %4s %4s %4s%10s",
/* SSBOND */	"%7 %3d %4s%c %4d%c   %4s%c %4d%c%4 %30s%10s",
/* TER */	"%6 %5d%6 %4s%c%4d%c%43 %10s",
/* TURN */	"%7 %3d %3s %4s%c%4d%c %4s%c%4d%c%4 %30s%10s",
/* TVECT */	"%7 %3d%10f%10f%10f%30s%10s",
/* USER */	"%4 %2s%64s%10s",
/* MODEL */	"%9 %5d%56 %10s",
/* ENDMDL */	"%70 %10s",
/* EXPDTA */	"%9 %c%60s%10s",			/* AUTHOR */
/* SYMDES */	"%7 %3d %59s%10s",			/* FTNOTE */
/* SYMOP */	"%5 %d %3d%10f%10f%10f%5 %10f%15 %10s",
/* MTXDES */	"%7 %3d %59s%10s",			/* FTNOTE */
/* CMPDES */	"%7 %3d %59s%10s",			/* FTNOTE */
/* CMPONT */	"%7 %3d %4s%c %4d%c %4s%c %4d%c%36 %10s",
/* TRNSFM */	"%7 %3d %3d %3d%52 %10s",
/* AGRDES */	"%7 %3d %59s%10s",			/* FTNOTE */
/* AGGRGT */	"%7 %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d%10s",
