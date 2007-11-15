// name: PDBrec.C
// author: J. Michael Word
// date written: 8/1/97
// purpose: Implementation for PDBrec

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#include "PDBrec.h"

#include <iostream>
using std::endl;

#ifdef OLD_STD_HDRS
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#else
#include <cstdio>
#include <cctype>
#include <cstring>
using std::toupper;
using std::strstr;
using std::strtok;
#endif

bool PDBrec::_MappingSEGIDtoChains = FALSE;
std::map<std::string, char> PDBrec::_SEGtoChainMap;

void PDBrec::clone(PDBrec* p, bool setmark) {
   p->_rep->_r          = _rep->_r;
   p->_rep->_e          = _rep->_e;
   if (setmark) {
//      p->_rep->_mark    = _rep->_mark;
//      p->_rep->_ok      = _rep->_ok;
//      p->_rep->_recMods = _rep->_recMods;
   }
   else {
      p->_rep->_mark    = 0;
      p->_rep->_ok      = TRUE;
      p->_rep->_recMods = 0;
   }
}

void PDBrec::getConect(int cvec[]) const {
   int n;
   const char* errmsg = hy36decode(5,  _rep->_r.conect.serialNum, 5, &n);
   //if (errmsg) throw std::runtime_error(errmsg);
   cvec[0] = n; 
   errmsg = hy36decode(5,  _rep->_r.conect.covalent[0], 5, &n);
   //if (errmsg) throw std::runtime_error(errmsg);
   cvec[1] = n;
   errmsg = hy36decode(5,  _rep->_r.conect.covalent[1], 5, &n);
   //if (errmsg) throw std::runtime_error(errmsg);
   cvec[2] = n;
   errmsg = hy36decode(5,  _rep->_r.conect.covalent[2], 5, &n);
   //if (errmsg) throw std::runtime_error(errmsg);
   cvec[3] = n;
   errmsg = hy36decode(5,  _rep->_r.conect.covalent[3], 5, &n);
   //if (errmsg) throw std::runtime_error(errmsg);
   cvec[4] = n;
   errmsg = hy36decode(5,  _rep->_r.conect.bonds[0].hydrogen[0], 5, &n);
   //if (errmsg) throw std::runtime_error(errmsg);
   cvec[5] = n;
   errmsg = hy36decode(5,  _rep->_r.conect.bonds[0].hydrogen[1], 5, &n);
   //if (errmsg) throw std::runtime_error(errmsg);
   cvec[6] = n;
   errmsg = hy36decode(5,  _rep->_r.conect.bonds[0].salt, 5, &n);
   //if (errmsg) throw std::runtime_error(errmsg);
   cvec[7] = n;
   errmsg = hy36decode(5,  _rep->_r.conect.bonds[1].hydrogen[0], 5, &n);
   //if (errmsg) throw std::runtime_error(errmsg);
   cvec[8] = n;
   errmsg = hy36decode(5,  _rep->_r.conect.bonds[1].hydrogen[1], 5, &n);
   //if (errmsg) throw std::runtime_error(errmsg);
   cvec[9] = n;
   errmsg = hy36decode(5,  _rep->_r.conect.bonds[1].salt, 5, &n);
   //if (errmsg) throw std::runtime_error(errmsg);
   cvec[10] = n;

//   cvec[0] = _rep->_r.conect.serialNum;
//   cvec[1] = _rep->_r.conect.covalent[0];
//   cvec[2] = _rep->_r.conect.covalent[1];
//   cvec[3] = _rep->_r.conect.covalent[2];
//   cvec[4] = _rep->_r.conect.covalent[3];
//   cvec[5] = _rep->_r.conect.bonds[0].hydrogen[0];
//   cvec[6] = _rep->_r.conect.bonds[0].hydrogen[1];
//   cvec[7] = _rep->_r.conect.bonds[0].salt;
//   cvec[8] = _rep->_r.conect.bonds[1].hydrogen[0];
//   cvec[9] = _rep->_r.conect.bonds[1].hydrogen[1];
//   cvec[10]= _rep->_r.conect.bonds[1].salt;
}

void PDBrec::setConect(int cvec[]) {
   char Hy36_cvec[11][6]; 
   const char* errmsg = hy36encode(5,  cvec[0], Hy36_cvec[0]);
   Hy36_cvec[0][5]='\0';
   strncpy(_rep->_r.conect.serialNum, Hy36_cvec[0], 6); 
   _rep->_r.conect.serialNum[5]='0'; 
   //if (errmsg) throw std::runtime_error(errmsg);

   errmsg = hy36encode(5, cvec[1], Hy36_cvec[1]);
   Hy36_cvec[1][5]='\0';
   strncpy(_rep->_r.conect.covalent[0], Hy36_cvec[1],6);
   _rep->_r.conect.covalent[0][5]='\0';
   //if (errmsg) throw std::runtime_error(errmsg);
 
   errmsg = hy36encode(5, cvec[2], Hy36_cvec[2]);
   Hy36_cvec[2][5]='\0';
   strncpy(_rep->_r.conect.covalent[1], Hy36_cvec[2],6);
   _rep->_r.conect.covalent[1][5]='\0';
   //if (errmsg) throw std::runtime_error(errmsg);

   errmsg = hy36encode(5, cvec[3], Hy36_cvec[3]);
   Hy36_cvec[3][5]='\0';
   strncpy(_rep->_r.conect.covalent[2], Hy36_cvec[3],6);
   _rep->_r.conect.covalent[2][5]='\0';
   //if (errmsg) throw std::runtime_error(errmsg);

   errmsg = hy36encode(5, cvec[4], Hy36_cvec[4]);
   Hy36_cvec[4][5]='\0';
   strncpy(_rep->_r.conect.covalent[3], Hy36_cvec[4],6);
   _rep->_r.conect.covalent[3][5]='\0';
   //if (errmsg) throw std::runtime_error(errmsg);

   errmsg = hy36encode(5, cvec[5], Hy36_cvec[5]);
   Hy36_cvec[5][5]='\0';
   strncpy(_rep->_r.conect.bonds[0].hydrogen[0], Hy36_cvec[5],6); 
   _rep->_r.conect.bonds[0].hydrogen[0][5]='\0';
   //if (errmsg) throw std::runtime_error(errmsg);

   errmsg = hy36encode(5, cvec[6], Hy36_cvec[6]);
   Hy36_cvec[6][5]='\0';
   strncpy(_rep->_r.conect.bonds[0].hydrogen[1], Hy36_cvec[6],6);
   _rep->_r.conect.bonds[0].hydrogen[1][5]='\0';
   //if (errmsg) throw std::runtime_error(errmsg);

   errmsg = hy36encode(5, cvec[7], Hy36_cvec[7]);
   Hy36_cvec[7][5]='\0';
   strncpy(_rep->_r.conect.bonds[0].salt, Hy36_cvec[7],6);
   _rep->_r.conect.bonds[0].salt[5]='\0';
   //if (errmsg) throw std::runtime_error(errmsg);

   errmsg = hy36encode(5, cvec[8], Hy36_cvec[8]);
   Hy36_cvec[8][5]='\0';
   strncpy(_rep->_r.conect.bonds[1].hydrogen[0], Hy36_cvec[8],6);
   _rep->_r.conect.bonds[1].hydrogen[0][5]='\0';
   //if (errmsg) throw std::runtime_error(errmsg);

   errmsg = hy36encode(5, cvec[9], Hy36_cvec[9]);
   Hy36_cvec[9][5]='\0';
   strncpy(_rep->_r.conect.bonds[1].hydrogen[1], Hy36_cvec[9],6);
   _rep->_r.conect.bonds[1].hydrogen[1][5]='\0';
   //if (errmsg) throw std::runtime_error(errmsg);

   errmsg = hy36encode(5, cvec[10], Hy36_cvec[10]);
   Hy36_cvec[10][5]='\0';
   strncpy(_rep->_r.conect.bonds[1].salt, Hy36_cvec[10],6);
   _rep->_r.conect.bonds[1].salt[5]='\0';
   //if (errmsg) throw std::runtime_error(errmsg);

//   _rep->_r.conect.serialNum            = cvec[0]; 
//   _rep->_r.conect.covalent[0]          = cvec[1];
//   _rep->_r.conect.covalent[1]          = cvec[2];
//   _rep->_r.conect.covalent[2]          = cvec[3];
//   _rep->_r.conect.covalent[3]          = cvec[4];
//   _rep->_r.conect.bonds[0].hydrogen[0] = cvec[5];
//   _rep->_r.conect.bonds[0].hydrogen[1] = cvec[6];
//   _rep->_r.conect.bonds[0].salt        = cvec[7];
//   _rep->_r.conect.bonds[1].hydrogen[0] = cvec[8];
//   _rep->_r.conect.bonds[1].hydrogen[1] = cvec[9];
//   _rep->_r.conect.bonds[1].salt        = cvec[10];
}

std::string PDBrec::recName() const {
	char fmtbuf[30];
	::sprintf(fmtbuf, "%-2.2s%4d%c%-3.3s%-4.4s%c",
		chain(), resno(), insCode(),
		resname(), atomname(), alt());
	return fmtbuf;
}

bool PDBrec::isWater() const {
   char fmtbuf[10];
   ::sprintf(fmtbuf, ":%-3.3s:", resname());
   return strstr(WATER_RESNAMES, fmtbuf) != NULL;
}

/*bool PDBrec::isHe() const {
   char resn[6];
   ::sprintf(resn, ":%-3.3s:", resname());
   return strstr(HE_RESNAMES, resn) != NULL;
}

bool PDBrec::isHf() const {
   char resn[6];
   ::sprintf(resn, ":%-3.3s:", resname());
   return strstr(HF_RESNAMES, resn) != NULL;
}

bool PDBrec::isHg() const {
   char resn[6];
   ::sprintf(resn, ":%-3.3s:", resname());
   return strstr(HG_RESNAMES, resn) != NULL;
}

bool PDBrec::isHo() const {
   char resn[6];
   ::sprintf(resn, ":%-3.3s:", resname());
   return strstr(HO_RESNAMES, resn) != NULL;
}

bool PDBrec::isHs() const {
   char resn[6];
   ::sprintf(resn, ":%-3.3s:", resname());
   return strstr(HS_RESNAMES, resn) != NULL;
}*/

#define SEGMAPBUFFSZ 500

int PDBrec::InstallMapOfSEGIDstoChains(const std::string map) {
   int mapentries = 0;
   if (!map.empty()) {
      char SEGmapBuffer[SEGMAPBUFFSZ+1];
      std::string sid;
      char   cid = ' ';
      strncpy(SEGmapBuffer, map.c_str(), SEGMAPBUFFSZ);
      SEGmapBuffer[SEGMAPBUFFSZ] = '\0';
      char *t = strtok(SEGmapBuffer, ","); // 1st segid
      while (t) {
         sid = FormatSegToChainKey(t);
	 t = strtok(NULL, ","); // chainid
	 if (t) {
	    cid = (*t && (*t != '_')) ? toupper(t[0])
	                              : ' ';
//	    if (*sid && cid) {
	    if (!sid.empty() && cid) {
			_SEGtoChainMap.insert(std::make_pair(sid, cid));
	       mapentries++;
	    }
	    t = strtok(NULL, ","); // next segid
         }
      }
      if (mapentries) { _MappingSEGIDtoChains = TRUE; }
   }
   return mapentries;
}

void PDBrec::segidLabel(const char* s) {
   strncpy(_rep->_r.atom.segID, s, 4);
   _rep->_r.atom.segID[4] = '\0';
   if (_MappingSEGIDtoChains) {
      SEGIDtoChain(s, one_char_chain(), _rep->_r.atom.residue.chainId);
   }
}

void PDBrec::MapSEGIDtoChain() {
   if (_MappingSEGIDtoChains) {
      SEGIDtoChain(
         segidLabel(), one_char_chain(), _rep->_r.atom.residue.chainId);
   }
}

void PDBrec::SEGIDtoChain(const char *seg, char cdefault, char* dest) {
   dest[0] = dest[1] = ' '; 
   dest[2] = '\0'; 
   if (_MappingSEGIDtoChains && seg) {
      std::string s = FormatSegToChainKey(seg);
	  std::map<std::string, char>::const_iterator i = _SEGtoChainMap.find(s);
		if (i != _SEGtoChainMap.end())
			//return i->second;
                        dest[1] = i->second;
                        return;
   }
	//return cdefault;
        dest[1] = cdefault;
        return;
}

std::string PDBrec::FormatSegToChainKey(const char *seg) {
   char buf[5];
   buf[0] = buf[1] = buf[2] = buf[3] = ' ';
   buf[4] = '\0';
   if (seg) {
      int i = 0;
      for(char *s = (char *)seg; *s && i < 4; s++) {
        buf[i++] = (*s == '_') ? ' ' : toupper(*s);
      }
   }
   return buf;
}

void PDBrec::DumpSEGIDtoChainMap(std::ostream& s, const char *t) {
   if (_MappingSEGIDtoChains) {
		for (std::map<std::string, char>::const_iterator i = _SEGtoChainMap.begin();
		i != _SEGtoChainMap.end(); ++i)
			s << t << "segID \"" << i->first
			<< "\" to chainID \"" << i->second << "\"" << endl;
   }
}
