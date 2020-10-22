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
"%s", /* UNKNOWN */
/* extra field included in ANISOU, ATOM and related records */
/* changed %5d %4s%c to %-5s %-3s%-2s to force three character resnames and two character chainIds */
"ANISOU%-5s %-4s%c%-3s%-2s%4s%c %7d%7d%7d%7d%7d%7d  %-4s%2s%-2s %-10s",	/* SIGUIJ */  /**non-standard**/
"ATOM  %-5s %-4s%c%-3s%-2s%4s%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%-2s %-10s", /* HETATM, SIGATM */  /**non-standard**/
"AUTHOR   %c%-60s%s",			/* COMPND, EXPDTA, JRNL, SOURCE */
"COMPND   %c%-60s%s",					/* AUTHOR */
"CONECT%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-5s         %s",
"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d%s",
"END                                                                   %s",
"FORMUL  %2D  %-4s%2D%c%-51s%s",
"FTNOTE %3D %-59s%s",					/* REMARK */
"HEADER    %-40s%-11s%c%-4s    %s",
"HELIX  %3D %3s %-3s%-2s %4s%c %-3s%-2s %4s%c%2D%-30s%s",
"HET    %-3s %-2s%4s%c  %5d     %-40s%s",
"HETATM%5s %-4s%c%-3s%-2s%4s%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%-2s %-10s", /**non-standard**/
"JRNL     %c%-60s%s",					/* AUTHOR */
"MASTER    %5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%s",
"MTRIX%1d %3d%10.6f%10.6f%10.6f     %10.5f   %2D          %s",
"OBSLTE  %2D %-9s %-10s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-4s%s",
"ORIGX%1d    %10.6f%10.6f%10.6f     %10.5f               %s",		/* SCALE */
"REMARK %3D %-59s%s",
"REVDAT %3D%2D %-9s %-7s %c       %-31s%s",
"SCALE%1d    %10.6f%10.6f%10.6f     %10.5f               %s",		/* ORIGX */
"SEQRES%4d%-2s %4d  %-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s%-4s%s",
"SHEET %4D %3s%2d %3s%-2s%4s%c %-3s%-2s%4s%c%2d %-4s%-3s%-2s%4s%c %-4s%-3s%-2s%4s%c%s",
"SIGATM%-5s %-4s%c%-3s%-2s%4s%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%-2s %-10s", /**non-standard**/
"SIGUIJ%-5s %-4s%c%-3s%2s%4s%c %7d%7d%7d%7d%7d%7d  %-4s%2s%-2s %-10s",	/* ANISOU */  /**non-standard**/
"SITE   %3d %3s %2d %-3s%-2s%4s%c %-3s%-2s%4s%c %-3s%-2s%4s%c %-3s%-2s%4s%c         %s",
"SOURCE   %c%-60s%s",					/* AUTHOR */
"SPRSDE  %2D %-9s %-10s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-4s%s",
"SSBOND %3D %-3s%-2s %4s%c   %-3s%-2s %4s%c    %-30s%s",
"TER   %-5s      %-3s%-2s%4s%c                                           %s",
"TURN   %3D %3s %-3s%-2s%4s%c %-3s%-2s%4s%c    %-30s%s",
"TVECT  %3D%10.5f%10.5f%10.5f%-30s%s",
"USER%-2s%-95s",
"MODEL    %5d                                                        %s",
"ENDMDL                                                                %s",
"EXPDTA   %c%-60s%s",					/* AUTHOR */
"SYMDES %3d %59s %s",					/* FTNOTE */
"SYMOP%1d %3d%10.6f%10.6f%10.6f     %10.5f               %s",
"MTXDES %3d %59s %s",					/* FTNOTE */
"CMPDES %3d %59s %s",					/* FTNOTE */
"CMPONT %3d %-3s%-2s %4s%c %-3s%-2s %4s%c                                    %s",
"TRNSFM %3d %3d %3d                                                    %s",
"AGRDES %3d %59s %s",					/* FTNOTE */
"AGGRGT %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d%s",
