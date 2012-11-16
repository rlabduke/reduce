//
// with Modifications by J. Michael Word, Duke University
//
//	Copyright (c) 1989, 1992 The Regents of the University of California.
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
#ifdef OLD_STD_HDRS
#include <string.h>
#else
#include <cstring>
using std::memset;
#endif

# ifndef NULL
# define	NULL		0
# endif

# ifndef HASSSCANFEXTERN
extern "C" int	sscanf(const char *, const char *, ...);
# endif

//
//	for each pdb record type there is a format reading in the
//	record values and for printing them out.
//
//	The actual format of a line written, is the print format
//	followed by blank padding to 72 characters, followed by
//	8 characters of file and line information.
//

static char const * const pdbRecordFormat[PDB::NUM_TYPES] = {
#include "read_format.i"
};

static char const * const pdbrun5[PDB::NUM_USER] = {
#include "pdbrun5_read.i"
};

static char const * const pdbrun6[PDB::NUM_USER] = {
#include "pdbrun6_read.i"
};

#include <stdio.h>
PDB::PDB(const char *buf)
{
	initialize_everything();

//	char dummy='\0'; // rmi added for 3 character resname 
//                          removed when 2 character chains were added
	const char	*fmt;
	Sheet		*sh;
	Residue		*sha0, *sha1;

	// convert pdb record to C structure

	memset(&rType, 0, sizeof(RecordType));
	memset(&unknown, 0, sizeof(Unknown));
	rType = getType(buf);
	if (rType < USER_PDBRUN)
		fmt = pdbRecordFormat[rType];
	else if (pdbrunInputVersion < 6)
		fmt = pdbrun5[rType - USER_PDBRUN];
	else
		fmt = pdbrun6[rType - USER_PDBRUN];
	switch (rType) {

	default:
	case UNKNOWN:
unknown:
		rType = UNKNOWN;		// in case of goto
		fmt = pdbRecordFormat[rType];
		(void) ::sprintf(unknown.junk, fmt, buf);
		break;

	case AGGRGT:
		if (0 > sscanf(buf, fmt, &aggrgt.serialNum,
				&aggrgt.numComponents,
				&aggrgt.cmpontSerialNums[0],
				&aggrgt.cmpontSerialNums[1],
				&aggrgt.cmpontSerialNums[2],
				&aggrgt.cmpontSerialNums[3],
				&aggrgt.cmpontSerialNums[4],
				&aggrgt.cmpontSerialNums[5],
				&aggrgt.cmpontSerialNums[6],
				&aggrgt.cmpontSerialNums[7],
				&aggrgt.cmpontSerialNums[8],
				&aggrgt.cmpontSerialNums[9],
				&aggrgt.cmpontSerialNums[10],
				&aggrgt.cmpontSerialNums[11],
				&aggrgt.cmpontSerialNums[12],
				&aggrgt.cmpontSerialNums[13],endOfLineInfo))
			goto unknown;
		break;

	case AGRDES:
	case CMPDES:
	case FTNOTE:
	case MTXDES:
	case REMARK:
	case SYMDES:
		if (0 > sscanf(buf, fmt, &agrdes.num, agrdes.text,endOfLineInfo))
			goto unknown;
		break;

	case ANISOU:
	case SIGUIJ:
		if (0 > sscanf(buf, fmt, &anisou.serialNum, anisou.name,
				&anisou.altLoc, anisou.residue.name, 
				&anisou.residue.chainId,
				&anisou.residue.seqNum,
				&anisou.residue.insertCode,
				&anisou.u[0], &anisou.u[1], &anisou.u[2],
				&anisou.u[3], &anisou.u[4], &anisou.u[5],
				anisou.segID, anisou.element, anisou.charge,
				anisou.annotation)) // *annotaton* is non-standard
			goto unknown;
		break;

	case ATOM:
	case HETATM:
	case SIGATM:
		if (0 > sscanf(buf, fmt, &atom.serialNum, atom.name,
				&atom.altLoc, atom.residue.name,  
				&atom.residue.chainId, &atom.residue.seqNum,
				&atom.residue.insertCode, &atom.xyz[0],
				&atom.xyz[1], &atom.xyz[2], &atom.occupancy,
				&atom.tempFactor, atom.segID, atom.element, atom.charge,
				atom.annotation)) // *annotaton* is non-standard
			goto unknown;
		break; 

	case AUTHOR:
	case COMPND:
	case EXPDTA:
	case JRNL:
	case SOURCE:
		if (0 > sscanf(buf, fmt, &author.continuation, author.data,endOfLineInfo))
			goto unknown;
		break;

	case CONECT:
		if (0 > sscanf(buf, fmt, &conect.serialNum,
				&conect.covalent[0], &conect.covalent[1],
				&conect.covalent[2], &conect.covalent[3],
				&conect.bonds[0].hydrogen[0],
				&conect.bonds[0].hydrogen[1],
				&conect.bonds[0].salt,
				&conect.bonds[1].hydrogen[0],
				&conect.bonds[1].hydrogen[1],
				&conect.bonds[1].salt,endOfLineInfo))
			goto unknown;
		break;

	case CMPONT:
		if (0 > sscanf(buf, fmt, &cmpont.seqNum,
				cmpont.residues[0].name,
				&cmpont.residues[0].chainId,
				&cmpont.residues[0].seqNum,
				&cmpont.residues[0].insertCode,
				cmpont.residues[1].name,
				&cmpont.residues[1].chainId,
				&cmpont.residues[1].seqNum,
				&cmpont.residues[1].insertCode,endOfLineInfo))
			goto unknown;
		break;

	case CRYST1:
		if (0 > sscanf(buf, fmt, &cryst1.a, &cryst1.b, &cryst1.c,
				&cryst1.alpha, &cryst1.beta, &cryst1.gamma,
				cryst1.spaceGrp, &cryst1.z,endOfLineInfo))
			goto unknown;
		break;

	case END:
	case ENDMDL:
		sscanf(buf, fmt, endOfLineInfo);
		break;

	case FORMUL:
		if (0 > sscanf(buf, fmt, &formul.component, formul.hetId,
				&formul.continuation, &formul.exclude,
				formul.formula,endOfLineInfo))
			goto unknown;
		break;

	case HEADER:
		if (0 > sscanf(buf, fmt, header.classification,
				header.timestamp, &header.type, header.id,endOfLineInfo))
			goto unknown;
		break;

	case HELIX:
		if (0 > sscanf(buf, fmt, &helix.serialNum, helix.id,
				helix.residues[0].name,
				&helix.residues[0].chainId,
				&helix.residues[0].seqNum,
				&helix.residues[0].insertCode,
				helix.residues[1].name,
				&helix.residues[1].chainId,
				&helix.residues[1].seqNum,
				&helix.residues[1].insertCode,
				&helix.type, helix.comment,endOfLineInfo))
			goto unknown;
		break;

	case HET:
		if (0 > sscanf(buf, fmt, het.hetGrp.name,
				&het.hetGrp.chainId, &het.hetGrp.seqNum,
				&het.hetGrp.insertCode, &het.numAtoms,
				het.text,endOfLineInfo))
			goto unknown;
		break;

	case MASTER:
		if (0 > sscanf(buf, fmt, &master.numRemark, &master.numFtnote,
				&master.numHet, &master.numHelix,
				&master.numSheet, &master.numTurn,
				&master.numSite, &master.numTransform,
				&master.numCoordinate, &master.numTer,
				&master.numConect, &master.numSeqres,endOfLineInfo))
			goto unknown;
		break;

	case MODEL:
		if (0 > sscanf(buf, fmt, &model.num,endOfLineInfo))
			goto unknown;
		break;

	case MTRIX:
		if (0 > sscanf(buf, fmt, &mtrix.rowNum, &mtrix.serialNum,
				&mtrix.m1, &mtrix.m2, &mtrix.m3, &mtrix.v,
				&mtrix.given,endOfLineInfo))
			goto unknown;
		break;

	case OBSLTE:
		if (0 > sscanf(buf, fmt, &obslte.continuation, obslte.timestamp,
				obslte.oldId, obslte.idMap[0],
				obslte.idMap[1], obslte.idMap[2],
				obslte.idMap[3], obslte.idMap[4],
				obslte.idMap[2], obslte.idMap[6],
				obslte.idMap[7],endOfLineInfo))
			goto unknown;
		break;

	case ORIGX:
		if (0 > sscanf(buf, fmt, &origx.rowNum, &origx.o1, &origx.o2,
				&origx.o3, &origx.t,endOfLineInfo))
			goto unknown;
		break;

	case REVDAT:
		if (0 > sscanf(buf, fmt, &revdat.modification,
				&revdat.continuation, revdat.timestamp,
				revdat.id, &revdat.modType,
				revdat.corrections,endOfLineInfo))
			goto unknown;
		break;

	case SCALE:
		if (0 > sscanf(buf, fmt, &scale.rowNum, &scale.s1, &scale.s2,
				&scale.s3, &scale.u,endOfLineInfo))
			goto unknown;
		break;

	case SEQRES:
		if (0 > sscanf(buf, fmt, &seqres.serialNum, &seqres.chainId,
				&seqres.count, seqres.names[0], seqres.names[1],
				seqres.names[2], seqres.names[3],
				seqres.names[4], seqres.names[5],
				seqres.names[6], seqres.names[7],
				seqres.names[8], seqres.names[9],
				seqres.names[10], seqres.names[11],
				seqres.names[12],endOfLineInfo))
			goto unknown;
		break;

	case SHEET:
		sh = &sheet;
		sha0 = &sh->atoms[0].residue;
		sha1 = &sh->atoms[1].residue;
		if (0 > sscanf(buf, fmt, &sh->strandNum, sh->id, &sh->count,
				sh->residues[0].name, &sh->residues[0].chainId,
				&sh->residues[0].seqNum,
				&sh->residues[0].insertCode,
				sh->residues[1].name, &sh->residues[1].chainId,
				&sh->residues[1].seqNum,
				&sh->residues[1].insertCode, &sh->sense,
				sh->atoms[0].name, sha0->name, &sha0->chainId,
				&sha0->seqNum, &sha0->insertCode,
				sh->atoms[1].name, sha1->name, &sha1->chainId,
				&sha1->seqNum, &sha1->insertCode,endOfLineInfo))
			goto unknown;
		break;

	case SITE:
		if (0 > sscanf(buf, fmt, &site.seqNum, site.id, &site.count,
				site.residues[0].name,
				&site.residues[0].chainId,
				&site.residues[0].seqNum,
				&site.residues[0].insertCode,
				site.residues[1].name,
				&site.residues[1].chainId,
				&site.residues[1].seqNum,
				&site.residues[1].insertCode,
				site.residues[2].name,
				&site.residues[2].chainId,
				&site.residues[2].seqNum,
				&site.residues[2].insertCode,
				site.residues[3].name,
				&site.residues[3].chainId,
				&site.residues[3].seqNum,
				&site.residues[3].insertCode,endOfLineInfo))
			goto unknown;
		break;

	case SPRSDE:
		if (0 > sscanf(buf, fmt, &sprsde.continuation,
				sprsde.timestamp, sprsde.id,
				sprsde.supersede[0], sprsde.supersede[1],
				sprsde.supersede[2], sprsde.supersede[3],
				sprsde.supersede[4], sprsde.supersede[5],
				sprsde.supersede[6], sprsde.supersede[7],endOfLineInfo))
			goto unknown;
		break;

	case SSBOND:
		if (0 > sscanf(buf, fmt, &ssbond.seqNum,
				ssbond.residues[0].name,
				&ssbond.residues[0].chainId,
				&ssbond.residues[0].seqNum,
				&ssbond.residues[0].insertCode,
				ssbond.residues[1].name,
				&ssbond.residues[1].chainId,
				&ssbond.residues[1].seqNum,
				&ssbond.residues[1].insertCode,
				ssbond.comment,endOfLineInfo))
			goto unknown;
		break;

	case SYMOP:
		if (0 > sscanf(buf, fmt, &symop.rowNum, &symop.serialNum,
				&symop.s1, &symop.s2, &symop.s3, &symop.t,endOfLineInfo))
			goto unknown;
		break;

	case TER:
		if (0 > sscanf(buf, fmt, &ter.serialNum, ter.residue.name,
				&ter.residue.chainId, &ter.residue.seqNum,
				&ter.residue.insertCode,endOfLineInfo))
			goto unknown;
		break;

	case TRNSFM:
		if (0 > sscanf(buf, fmt, &trnsfm.resultSerialNum,
				&trnsfm.applySerialNum,
				&trnsfm.sourceSerialNum,endOfLineInfo))
			goto unknown;
		break;

	case TURN:
		if (0 > sscanf(buf, fmt, &turn.seqNum, turn.id,
				turn.residues[0].name,
				&turn.residues[0].chainId,
				&turn.residues[0].seqNum,
				&turn.residues[0].insertCode,
				turn.residues[1].name,
				&turn.residues[1].chainId,
				&turn.residues[1].seqNum,
				&turn.residues[1].insertCode, turn.comment,endOfLineInfo))
			goto unknown;
		break;

	case TVECT:
		if (0 > sscanf(buf, fmt, &tvect.serialNum, &tvect.t1,
				&tvect.t2, &tvect.t3, tvect.comment,endOfLineInfo))
			goto unknown;
		break;

user:
		rType = USER;
		fmt = pdbRecordFormat[rType];
	case USER:
		if (0 > sscanf(buf, fmt, user.subtype, user.text,endOfLineInfo))
			goto unknown;
		break;

	case USER_PDBRUN:
		if (0 > sscanf(buf, fmt, &userPdbrun.version))
			goto user;
		pdbrunInputVersion = userPdbrun.version;
		break;

	case USER_EYEPOS:
		if (0 > sscanf(buf, fmt, &userEyePos.xyz[0],
				&userEyePos.xyz[1], &userEyePos.xyz[2]))
			goto user;
		break;

	case USER_ATPOS:
		if (0 > sscanf(buf, fmt, &userAtPos.xyz[0],
				&userAtPos.xyz[1], &userAtPos.xyz[2]))
			goto user;
		break;

	case USER_WINDOW:
		if (0 > sscanf(buf, fmt, &userWindow.left, &userWindow.right,
				&userWindow.bottom, &userWindow.top,
				&userWindow.hither, &userWindow.yon))
			goto user;
		break;

	case USER_FOCUS:
		if (0 > sscanf(buf, fmt, &userFocus.focus))
			goto user;
		break;

	case USER_VIEWPORT:
		if (0 > sscanf(buf, fmt, &userViewport.xmin, &userViewport.xmax,
				&userViewport.ymin, &userViewport.ymax))
			goto user;
		break;

	case USER_BGCOLOR:
		if (pdbrunInputVersion < 6) {
			if (0 > ::sscanf(buf, fmt, &userBgColor.rgb[0],
					&userBgColor.rgb[1],
					&userBgColor.rgb[2]))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userBgColor.rgb[0],
				&userBgColor.rgb[1], &userBgColor.rgb[2]))
			goto user;
		break;

	case USER_ANGLE:
		if (pdbrunInputVersion < 6) {
			if (0 > sscanf(buf, fmt, &userAngle.which,
					&userAngle.atom0, &userAngle.atom1,
					&userAngle.atom2, &userAngle.atom3,
					&userAngle.angle))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userAngle.atom0,
				&userAngle.atom1, &userAngle.atom2,
				&userAngle.atom3, &userAngle.angle))
			goto user;
		break;

	case USER_DISTANCE:
		if (pdbrunInputVersion < 6) {
			if (0 > sscanf(buf, fmt, &userDistance.which,
					&userDistance.atom0,
					&userDistance.atom1,
					&userDistance.distance))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userDistance.atom0,
				&userDistance.atom1, &userDistance.distance))
			goto user;
		break;

	case USER_FILE:
		if (pdbrunInputVersion < 6) {
			if (0 > sscanf(buf, fmt, userFile.filename))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userFile.model,
							userFile.filename))
			goto user;
		break;

	case USER_MARKNAME:
		if (pdbrunInputVersion < 6
		|| 0 > sscanf(buf, fmt, userMarkname.markname))
			goto user;
		break;

	case USER_MARK:
		if (pdbrunInputVersion < 6
		|| 0 > sscanf(buf, fmt, userMark.markname))
			goto user;
		break;

	case USER_CNAME:
		if (pdbrunInputVersion < 6) {
			if (0 > ::sscanf(buf, fmt, userCName.name,
					&userCName.rgb[0], &userCName.rgb[1],
					&userCName.rgb[2]))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userCName.rgb[0],
				&userCName.rgb[1], &userCName.rgb[2],
				userCName.name))
			goto user;
		break;

	case USER_COLOR:
		if (pdbrunInputVersion < 6) {
			if (0 > ::sscanf(buf, fmt, userColor.spec,
					&userColor.rgb[0], &userColor.rgb[1],
					&userColor.rgb[2]))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userColor.rgb[0],
				&userColor.rgb[1], &userColor.rgb[2],
				userColor.spec))
			goto user;
		break;

	case USER_RADIUS:
		if (0 > sscanf(buf, fmt, &userRadius.radius))
			goto user;
		break;

	case USER_OBJECT:
		if (pdbrunInputVersion < 6) {
			if (0 > sscanf(buf, fmt, &userObject.model))
				goto user;
		}
		break;

	case USER_ENDOBJ:
		if (pdbrunInputVersion < 6) {
			if (0 > sscanf(buf, fmt, &userEndObj.model))
				goto user;
		}
		break;

	case USER_CHAIN:
		if (pdbrunInputVersion < 6) {
			if (0 > ::sscanf(buf, fmt, &userChain.atom0,
					&userChain.atom1))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userChain.atom0,
				&userChain.atom1))
			goto user;
		break;

	case USER_GFX_BEGIN:
		if (pdbrunInputVersion < 6
		|| 0 > sscanf(buf, fmt, userGfxBegin.unknown))
			goto user;
		userGfxBegin.primitive = getGfxType(userGfxBegin.unknown);
		break;

	case USER_GFX_END:
		if (pdbrunInputVersion < 6)
			goto user;
		break;

	case USER_GFX_COLOR:
		if (pdbrunInputVersion < 6) {
			if (0 > ::sscanf(buf, fmt, userGfxColor.spec,
					&userGfxColor.rgb[0],
					&userGfxColor.rgb[1],
					&userGfxColor.rgb[2]))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userGfxColor.rgb[0],
				&userGfxColor.rgb[1], &userGfxColor.rgb[2],
				userGfxColor.spec))
			goto user;
		break;

	case USER_GFX_NORMAL:
		if (pdbrunInputVersion < 6
		|| 0 > sscanf(buf, fmt, &userGfxNormal.xyz[0],
				&userGfxNormal.xyz[1],
				&userGfxNormal.xyz[2]))
			goto user;
		break;

	case USER_GFX_VERTEX:
		if (pdbrunInputVersion < 6
		|| 0 > sscanf(buf, fmt, &userGfxVertex.xyz[0],
				&userGfxVertex.xyz[1],
				&userGfxVertex.xyz[2]))
			goto user;
		break;

	case USER_GFX_FONT:
		if (pdbrunInputVersion < 6) {
			if (0 > ::sscanf(buf, fmt, userGfxFont.name,
					&userGfxFont.size))
				goto user;
		} else if (0 > sscanf(buf, fmt, &userGfxFont.size,
				userGfxFont.name))
			goto user;
		break;

	case USER_GFX_TEXTPOS:
		if (pdbrunInputVersion < 6
		|| 0 > sscanf(buf, fmt, &userGfxTextPos.xyz[0],
				&userGfxTextPos.xyz[1], &userGfxTextPos.xyz[2]))
			goto user;
		break;

	case USER_GFX_LABEL:
		if (pdbrunInputVersion < 6) {
			if (0 > ::sscanf(buf, fmt, &userGfxLabel.xyz[0],
					&userGfxLabel.xyz[1],
					&userGfxLabel.xyz[2],
					userGfxLabel.text))
				goto user;
		} else if (0 > sscanf(buf, fmt, userGfxLabel.text))
			goto user;
		// TODO: process text?
		break;

	case USER_GFX_MOVE:
		if (pdbrunInputVersion >= 6
		|| 0 > sscanf(buf, fmt, &userGfxMove.xyz[0],
				&userGfxMove.xyz[1], &userGfxMove.xyz[2]))
			goto user;
		break;

	case USER_GFX_DRAW:
		if (pdbrunInputVersion >= 6
		|| 0 > sscanf(buf, fmt, &userGfxDraw.xyz[0],
				&userGfxDraw.xyz[1], &userGfxDraw.xyz[2]))
			goto user;
		break;

	case USER_GFX_MARKER:
		if (pdbrunInputVersion >= 6
		|| 0 > sscanf(buf, fmt, &userGfxMarker.xyz[0],
				&userGfxMarker.xyz[1],
				&userGfxMarker.xyz[2]))
			goto user;
		break;

	case USER_GFX_POINT:
		if (pdbrunInputVersion >= 6
		|| 0 > sscanf(buf, fmt, &userGfxPoint.xyz[0],
				&userGfxPoint.xyz[1], &userGfxPoint.xyz[2]))
			goto user;
		break;
	}
}
