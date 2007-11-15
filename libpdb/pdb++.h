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
//	Based on Brookhaven National Laboratory Protein Data Bank, Feb 1992
//
//	C structure declarations
//

#ifndef PDB_H
#define	PDB_H

#include <iostream>

#define PDB_WORKAROUND 1

#ifndef BOOLPREDEFINED
typedef int bool;
#endif

// JMW  8/6/98:  BufLen has been expanded from 81 to 101 to include a non-standard
//               .annotaton field in Atom and Aniosu type records so that the SegID
//               elem and charge fields could be free from unneccessary modification
// RMI  10/3/07  Hy36Num added to allow usage of the Hybrid 36 atom numbering scheme
// RMI 10/10/07  Modified to allow two character chain ids with chainId[3] 
 
class PDB {
public:
#ifdef PDB_WORKAROUND
	enum { BufLen = 101, PDBRUNVersion = 6 };
#else
	static const int BufLen = 101;		// PDB record length (80 + 20 + '\0')
	static const int PDBRUNVersion = 6;	// Best version generated
#endif
	typedef char	Date[10];
        typedef char    Hy36Num[6];             // hybrid 36 atom number
        typedef char    Hy36seqNum[5];          // hybrid 36 residue number
	typedef char	AName[5];		// atom name - NO2*
	typedef char	RName[5];		// residue name - ALA
	typedef char	PName[5];		// pdb name - 9lyz
	typedef char	Id[4];			// generic short id field
	typedef double	Real;			// size of floating point

	struct Residue {			// residue info
		RName	   name;
		char	   chainId[3];
		Hy36seqNum seqNum;
		char	   insertCode;
	};

	// graphics primitive types
	enum GfxType {
		GFX_UNKNOWN, GFX_POINTS, GFX_MARKERS, GFX_LINES,
		GFX_LINE_STRIP, GFX_LINE_LOOP, GFX_TRIANGLES,
		GFX_TRIANGLE_STRIP, GFX_TRIANGLE_FAN, GFX_QUADS,
		GFX_QUAD_STRIP, GFX_POLYGON
	};

	//
	//	structures declarations for each record type
	//

	struct Unknown {
		char	junk[BufLen];
	};
	struct Aggrgt {
		int	serialNum;
		int	numComponents;
		int	cmpontSerialNums[14];
		char	reclabel[8];
	};
	struct Agrdes {
		int	num;
		char	text[60];
		char	reclabel[8];
	};
	struct Anisou {
		Hy36Num	serialNum;
		AName	name;
		char	altLoc;
		Residue	residue;
		int	u[6];
		char	segID[5];
		char	element[3];
		char	charge[3];
		char	annotation[11]; // *annotaton* is non-standard
	};
	struct Atom {
		Hy36Num	serialNum;
		AName	name;
		char	altLoc;
		Residue	residue;
		Real	xyz[3];
		Real	occupancy, tempFactor;
		char	segID[5];
		char	element[3];
		char	charge[3];
		char	annotation[11];  // *annotaton* is non-standard
	};
	struct Author {
		char	data[61];
		char	continuation;
		char	reclabel[8];
	};
	typedef Agrdes	Cmpdes;
	struct Cmpont {
		Hy36seqNum	seqNum;
		Residue	residues[2];
		char	reclabel[8];
	};
	typedef Author	Compnd;
	struct Conect {
		Hy36Num	serialNum;
		char	covalent[4][6];
		struct {
			char	hydrogen[2][6];
			char	salt[6];
		}	bonds[2];
		char	reclabel[8];
	};
	struct Cryst1 {
		Real	a, b, c;
		Real	alpha, beta, gamma;
		char	spaceGrp[12];
		int	z;
		char	reclabel[8];
	};
	struct End {
		char	reclabel[8];
	};
	struct Endmdl {
		char	reclabel[8];
	};
	typedef Author	Expdta;
	struct Formul {
		int	component;
		RName	hetId;
		int	continuation;
		char	exclude;	// * to exclude
		char	formula[52];
		char	reclabel[8];
	};
	typedef Agrdes	Ftnote;
	struct Header {
		char	classification[41];
		Date	timestamp;
		PName	id;
		char	type;
		char	reclabel[8];
	};
	struct Helix {
		int	serialNum;
		Id	id;
		Residue	residues[2];
		int	type;
		char	comment[31];
		char	reclabel[8];
	};
	struct Het {
		Residue	hetGrp;
		int	numAtoms;
		char	text[41];
		char	reclabel[8];
	};
	typedef Atom	Hetatm;
	typedef Author	Jrnl;
	struct Master {
		int	numRemark;
		int	numFtnote;
		int	numHet;
		int	numHelix;
		int	numSheet;
		int	numTurn;
		int	numSite;
		int	numTransform;
		int	numCoordinate;
		int	numTer;
		int	numConect;
		int	numSeqres;
		char	reclabel[8];
	};
	struct Model {
		int	num;
		char	reclabel[8];
	};
	struct Mtrix {
		int	rowNum;
		int	serialNum;
		Real	m1, m2, m3, v;
		int	given;
		char	reclabel[8];
	};
	typedef Agrdes	Mtxdes;
	struct Obslte {
		int	continuation;
		Date	timestamp;
		PName	oldId;
		PName	idMap[8];
		char	reclabel[8];
	};
	struct Origx {
		int	rowNum;
		Real	o1, o2, o3, t;
		char	reclabel[8];
	};
	typedef Ftnote	Remark;
	struct Revdat {
		int	modification;
		int	continuation;
		Date	timestamp;
		char	id[8];
		char	modType;
		char	corrections[31];
		char	reclabel[8];
	};
	struct Scale {
		int	rowNum;
		Real	s1, s2, s3, u;
		char	reclabel[8];
	};
	struct Seqres {
		int	serialNum;
		char	chainId[3];
		int	count;
		RName	names[13];
		char	reclabel[8];
	};
	struct Sheet {
		int	strandNum;
		Id	id;
		int	count;
		Residue	residues[2];
		int	sense;
		struct {
			AName	name;
			Residue	residue;
		}		atoms[2];
		char	reclabel[8];
	};
	typedef Atom	Sigatm;
	typedef Anisou	Siguij;
	struct Site {
		Hy36seqNum	seqNum;
		Id	id;
		int	count;
		Residue	residues[4];
		char	reclabel[8];
	};
	typedef Author	Source;
	struct Sprsde {
		int	continuation;
		Date	timestamp;
		PName	id;
		PName	supersede[8];
		char	reclabel[8];
	};
	struct Ssbond {
		Hy36seqNum	seqNum;
		Residue	residues[2];
		char	comment[31];
		char	reclabel[8];
	};
	typedef Agrdes	Symdes;
	struct Symop {
		int	rowNum;
		int	serialNum;
		Real	s1, s2, s3, t;
		char	reclabel[8];
	};
	struct Ter {
		Hy36Num	serialNum;
		Residue	residue;
		char	reclabel[8];
	};
	struct Trnsfm {
		int	resultSerialNum;
		int	applySerialNum;
		int	sourceSerialNum;
		char	reclabel[8];
	};
	struct Turn {
		Hy36seqNum	seqNum;
		Id	id;
		Residue	residues[2];
		char	comment[31];
		char	reclabel[8];
	};
	struct Tvect {
		int	serialNum;
		Real	t1, t2, t3;
		char	comment[31];
		char	reclabel[8];
	};
	struct User {
		char	subtype[3];
		char	text[67];
		char	reclabel[8];
	};
	struct UserPdbrun {
		int	version;
	};
	struct UserEyePos {
		Real	xyz[3];
	};
	typedef UserEyePos UserAtPos;
	struct UserWindow {
		Real	left, right, bottom, top, hither, yon;
	};
	struct UserFocus {
		Real	focus;
	};
	struct UserViewport {
		Real	xmin, xmax, ymin, ymax;
	};
	struct UserBgColor {
		Real	rgb[3];
	};
	struct UserAngle {
		int	atom0, atom1, atom2, atom3;
		Real	angle;
		int	which;			// version 5 -- obsolete
	};
	struct UserDistance {
		int	atom0, atom1;
		Real	distance;
		int	which;			// version 5 -- obsolete
	};
	struct UserFile {
		char	filename[62];
		int	model;			// not in version 5
	};
	struct UserMarkname {			// not in version 5
		char	markname[58];
	};
	typedef UserMarkname UserMark;		// not in version 5
	struct UserCName {
		Real	rgb[3];
		char	name[39];
	};
	struct UserColor {
		Real	rgb[3];
		char	spec[39];
	};
	struct UserRadius {
		Real	radius;
	};
	struct UserObject {
		int	model;			// version 5 -- obsolete
	};
	struct UserEndObj {
		int	model;			// version 5 -- obsolete
	};
	struct UserChain {
		int	atom0, atom1;
	};
	struct UserGfxBegin {			// not in version 5
		GfxType	primitive;
		char	unknown[33];
	};
	// no structure for USER  GFX END
	typedef UserColor UserGfxColor;
	struct UserGfxNormal {			// not in version 5
		Real	xyz[3];
	};
	typedef UserGfxNormal UserGfxVertex;
	struct UserGfxFont {
		int	size;
		char	name[54];
	};
	struct UserGfxTextPos {			// not in version 5
		Real	xyz[3];
	};
	struct UserGfxLabel {
		Real	xyz[3];			// version 5 -- obsolete
		char	text[57];		// 27 in version 5
	};
	struct UserGfxMove {			// version 5 -- obsolete
		Real	xyz[3];
	};
	typedef UserGfxMove UserGfxDraw;	// version 5 -- obsolete
	typedef UserGfxMove UserGfxMarker;	// version 5 -- obsolete
	typedef UserGfxMove UserGfxPoint;	// version 5 -- obsolete

	enum RecordType { UNKNOWN, ANISOU, ATOM, AUTHOR, COMPND, CONECT, CRYST1,
		END, FORMUL, FTNOTE, HEADER, HELIX, HET, HETATM, JRNL, MASTER,
		MTRIX, OBSLTE, ORIGX, REMARK, REVDAT, SCALE, SEQRES, SHEET,
		SIGATM, SIGUIJ, SITE, SOURCE, SPRSDE, SSBOND, TER, TURN, TVECT,
		USER, MODEL, ENDMDL, EXPDTA, SYMDES, SYMOP, MTXDES, CMPDES,
		CMPONT, TRNSFM, AGRDES, AGGRGT,
		USER_PDBRUN, USER_EYEPOS, USER_ATPOS, USER_WINDOW, USER_FOCUS,
		USER_VIEWPORT, USER_BGCOLOR, USER_ANGLE, USER_DISTANCE,
		USER_FILE, USER_MARKNAME, USER_MARK, USER_CNAME, USER_COLOR,
		USER_RADIUS, USER_OBJECT, USER_ENDOBJ, USER_CHAIN,
		USER_GFX_BEGIN, USER_GFX_END, USER_GFX_COLOR, USER_GFX_NORMAL,
		USER_GFX_VERTEX, USER_GFX_FONT, USER_GFX_TEXTPOS,
		USER_GFX_LABEL,
		USER_GFX_MOVE, USER_GFX_DRAW, USER_GFX_MARKER, USER_GFX_POINT
	};
#ifdef PDB_WORKAROUND
	enum { NUM_TYPES = AGGRGT + 1,
		NUM_USER = USER_GFX_POINT - USER_PDBRUN + 1 };
#else
	static const int	NUM_TYPES = AGGRGT + 1;
	static const int	NUM_USER = USER_GFX_POINT - USER_PDBRUN + 1;
#endif
private:
	RecordType	rType;
	static int	pdbrunInputVersion, pdbrunOutputVersion;
	static int	byteCmp(const PDB &l, const PDB &r);
	
	void initialize_everything();
public:
	union {
		Unknown	unknown;
		Aggrgt	aggrgt;
		Agrdes	agrdes;
		Anisou	anisou;
		Atom	atom;
		Author	author;
		Cmpdes	cmpdes;
		Cmpont	cmpont;
		Compnd	compnd;
		Conect	conect;
		Cryst1	cryst1;
		End     end;    // added
		Endmdl  endmdl; // added
		Expdta	expdta;
		Formul	formul;
		Ftnote	ftnote;
		Header	header;
		Helix	helix;
		Het	het;
		Hetatm	hetatm;
		Jrnl	jrnl;
		Master	master;
		Model	model;
		Mtrix	mtrix;
		Mtxdes	mtxdes;
		Obslte	obslte;
		Origx	origx;
		Remark	remark;
		Revdat	revdat;
		Scale	scale;
		Seqres	seqres;
		Sheet	sheet;
		Sigatm	sigatm;
		Siguij	siguij;
		Site	site;
		Source	source;
		Sprsde	sprsde;
		Ssbond	ssbond;
		Symdes	symdes;
		Symop	symop;
		Ter	ter;
		Trnsfm	trnsfm;
		Turn	turn;
		Tvect	tvect;
		User	user;
		UserPdbrun	userPdbrun;
		UserEyePos	userEyePos;
		UserAtPos	userAtPos;
		UserWindow	userWindow;
		UserFocus	userFocus;
		UserViewport	userViewport;
		UserBgColor	userBgColor;
		UserAngle	userAngle;
		UserDistance	userDistance;
		UserFile	userFile;
		UserMarkname	userMarkname;
		UserMark	userMark;
		UserCName	userCName;
		UserColor	userColor;
		UserRadius	userRadius;
		UserObject	userObject;
		UserEndObj	userEndObj;
		UserChain	userChain;
		UserGfxBegin	userGfxBegin;
		// no UserGfxEnd structure
		UserGfxColor	userGfxColor;
		UserGfxNormal	userGfxNormal;
		UserGfxVertex	userGfxVertex;
		UserGfxFont	userGfxFont;
		UserGfxTextPos	userGfxTextPos;
		UserGfxLabel	userGfxLabel;
		UserGfxMove	userGfxMove;
		UserGfxDraw	userGfxDraw;
		UserGfxMarker	userGfxMarker;
		UserGfxPoint	userGfxPoint;
	};
	char endOfLineInfo[9]; // added to keep identifying info

	PDB();
	PDB(RecordType t);
	PDB(const char *buf);
	
	RecordType	type() const { return rType; }
	void		type(RecordType t);
	const char	*chars() const;
	static int	PdbrunInputVersion() { return pdbrunInputVersion; }
	static int	PdbrunOutputVersion() { return pdbrunOutputVersion; }
	static void	PdbrunInputVersion(int v) { pdbrunInputVersion = v; }
	static void	PdbrunOutputVersion(int v) { pdbrunOutputVersion = v; }
	static RecordType
			getType(const char *buf);
	static GfxType	getGfxType(const char *buf);
	static const char
			*gfxChars(GfxType gt);
	static int	sscanf(const char *, const char *, ...);
	static int	sprintf(char *, const char *, ...);

	inline bool operator==(const PDB &r) const {
				if (rType != r.rType)
					return 0;
				return byteCmp(*this, r) == 0;
			}
	inline bool operator!=(const PDB &r) const {
				if (rType != r.rType)
					return 1;
				return byteCmp(*this, r) != 0;
			}

	friend std::istream	&operator>>(std::istream &s, PDB &p);
};

inline std::ostream &
operator<<(std::ostream &s, const PDB &p)
{
	s << p.chars();
	return s;
}

# endif // PDB_H
