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

#ifdef OLD_STD_HDRS
#include <stdio.h>
#include <ctype.h>
#else
#include <cstdio>
#include <cctype>
using std::endl;
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
   cvec[0] = _rep->_r.conect.serialNum;
   cvec[1] = _rep->_r.conect.covalent[0];
   cvec[2] = _rep->_r.conect.covalent[1];
   cvec[3] = _rep->_r.conect.covalent[2];
   cvec[4] = _rep->_r.conect.covalent[3];
   cvec[5] = _rep->_r.conect.bonds[0].hydrogen[0];
   cvec[6] = _rep->_r.conect.bonds[0].hydrogen[1];
   cvec[7] = _rep->_r.conect.bonds[0].salt;
   cvec[8] = _rep->_r.conect.bonds[1].hydrogen[0];
   cvec[9] = _rep->_r.conect.bonds[1].hydrogen[1];
   cvec[10]= _rep->_r.conect.bonds[1].salt;
}

void PDBrec::setConect(int cvec[]) {
   _rep->_r.conect.serialNum            = cvec[0];
   _rep->_r.conect.covalent[0]          = cvec[1];
   _rep->_r.conect.covalent[1]          = cvec[2];
   _rep->_r.conect.covalent[2]          = cvec[3];
   _rep->_r.conect.covalent[3]          = cvec[4];
   _rep->_r.conect.bonds[0].hydrogen[0] = cvec[5];
   _rep->_r.conect.bonds[0].hydrogen[1] = cvec[6];
   _rep->_r.conect.bonds[0].salt        = cvec[7];
   _rep->_r.conect.bonds[1].hydrogen[0] = cvec[8];
   _rep->_r.conect.bonds[1].hydrogen[1] = cvec[9];
   _rep->_r.conect.bonds[1].salt        = cvec[10];
}

std::string PDBrec::recName() const {
	char fmtbuf[30];
	::sprintf(fmtbuf, "%c%4d%c%-3.3s%-4.4s%c",
		chain(), resno(), insCode(),
		resname(), atomname(), alt());
	return fmtbuf;
}

bool PDBrec::isWater() const {
   char fmtbuf[10];
   ::sprintf(fmtbuf, ":%-3.3s:", resname());
   return strstr(WATER_RESNAMES, fmtbuf) != NULL;
}

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
      _rep->_r.atom.residue.chainId =
      SEGIDtoChain(s, chain());
   }
}

void PDBrec::MapSEGIDtoChain() {
   if (_MappingSEGIDtoChains) {
      _rep->_r.atom.residue.chainId =
         SEGIDtoChain(segidLabel(), chain());
   }
}

char PDBrec::SEGIDtoChain(const char *seg, char cdefault) {
   if (_MappingSEGIDtoChains && seg) {
      std::string s = FormatSegToChainKey(seg);
	  std::map<std::string, char>::const_iterator i = _SEGtoChainMap.find(s);
		if (i != _SEGtoChainMap.end())
			return i->second;
   }
	return cdefault;
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

void PDBrec::DumpSEGIDtoChainMap(ostream& s, const char *t) {
   if (_MappingSEGIDtoChains) {
		for (std::map<std::string, char>::const_iterator i = _SEGtoChainMap.begin();
		i != _SEGtoChainMap.end(); ++i)
			s << t << "segID \"" << i->first
			<< "\" to chainID \"" << i->second << "\"" << endl;
   }
}
