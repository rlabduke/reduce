// name: Dict.C
// author: J. Michael Word
// date written: 8/7/97
// purpose: Implementation for Dict

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#include "Dict.h"

// ---------------------------------------------------
//  make an empty dictionary for an approx. sz entries
// ---------------------------------------------------
template <class K, class T>
Dict<K,T>::Dict(const unsigned long sz) {

   _M = genHashM(sz);
   _n = 0;

   _Heads = new DictLink<K,T> *[_M];

   for (unsigned long i = 0; i < _M; i++) {
      _Heads[i] = NULL;
   }
}

// ---------------------
//  destroy a dictionary
// ---------------------
template <class K, class T>
Dict<K,T>::~Dict() {

   reset();	     // drop all entries
   delete [] _Heads;
}

// -----------------------------------------------
//  associate a value with a key in the dictionary
// -----------------------------------------------
template <class K, class T>
bool Dict<K,T>::put(const K& key, const T& value) {
   unsigned long head = hash(key, _M);
   DictLink<K,T> *prev = NULL, *curr = _Heads[head];

   while (curr && (curr->key() < key)) {
      prev = curr;
      curr = curr->next();
   }

   if (curr && (curr->key() == key)) {	       // keys match
      curr->value(value);
      return TRUE;
   }
   else {			 // no previous entry with this key
      DictLink<K,T> *x = new DictLink<K,T>(key, value);
      if (prev==NULL) { x->next(_Heads[head]); _Heads[head] = x;  }
      else            { x->next(prev->next(x)); }
      _n++;
      return FALSE;
   }
}

// ------------------------------------------------
// Associate a value* with a key in the dictionary.
// The memory pointed to is owned by and freed by
// the dict.
// ------------------------------------------------
template <class K, class T>
bool Dict<K,T>::put(const K& key, T* vp) {
   unsigned long head = hash(key, _M);
   DictLink<K,T> *prev = NULL, *curr = _Heads[head];

   while (curr && (curr->key() < key)) {
      prev = curr;
      curr = curr->next();
   }

   if (curr && (curr->key() == key)) {	       // keys match
      curr->value(vp);
      return TRUE;
   }
   else {			 // no previous entry with this key
      DictLink<K,T> *x = new DictLink<K,T>(key, vp);
      if (prev==NULL) { x->next(_Heads[head]); _Heads[head] = x;  }
      else            { x->next(prev->next(x)); }
      _n++;
      return FALSE;
   }
}

// --------------------------------------------------------
//  look up an entry by key and return the associated value
// --------------------------------------------------------
template <class K, class T>
bool Dict<K,T>::get(const K &key, T& retval) const {
   DictLink<K,T> *curr = _Heads[hash(key, _M)];

   while (curr && (curr->key() < key)) {
      curr = curr->next();
   }

   if (curr && (curr->key() == key) && curr->value()) {
	 retval = *(curr->value());
	 return TRUE;
   }
   return FALSE;
}

// --------------------------------------------------------
//  look up an entry by key and return the value pointer
// --------------------------------------------------------
template <class K, class T>
T* Dict<K,T>::get(const K &key) const {
   DictLink<K,T> *curr = _Heads[hash(key, _M)];

   while (curr && (curr->key() < key)) {
      curr = curr->next();
   }

   if (curr && (curr->key() == key)) {
      return curr->value();
   }
   else {
      return NULL;
   }
}

// --------------------------------------------------------
//  look up an entry by key and return if it exists
// --------------------------------------------------------
template <class K, class T>
bool Dict<K,T>::exists(const K &key) const {
   DictLink<K,T> *curr = _Heads[hash(key, _M)];

   while (curr && (curr->key() < key)) {
      curr = curr->next();
   }

   return (curr && (curr->key() == key));
}
// ---------------------------------------
//  remove one of our key/sequence entries
// ---------------------------------------
template <class K, class T>
bool Dict<K,T>::drop(const K &key) {
   unsigned long head = hash(key, _M);
   DictLink<K,T> *prev = NULL, *curr = _Heads[head];

   while (curr && (curr->key() < key)) {
      prev = curr;
      curr = curr->next();
   }

   if (curr && (curr->key() == key)) {		// keys match
      if (prev==NULL) { _Heads[head] = curr->next(NULL);  }
      else            {     prev->next(curr->next(NULL)); }
      delete curr;
      _n--;
      return TRUE;
   }
   else { return FALSE; }
}

// ----------------------
//  drop all dict entries
// ----------------------
template <class K, class T>
void Dict<K,T>::reset() {
   DictLink<K,T> *curr = NULL, *next = NULL;

   for (unsigned long i = 0; i < _M; i++) {
      for (curr = _Heads[i]; curr; curr = next) {
	 next = curr->next(NULL);
	 delete curr;
      }
   }

   _n = 0;
}

// -----------------------
//  display dict structure
// -----------------------
template <class K, class T>
void Dict<K,T>::dumpStructure(ostream& os) const {
   os << endl << _n << " keys in [" << _M << "] (";
   for (unsigned long i = 0; i < _M; i++) {
      int clump = 1;
      for (DictLink<K,T> *curr = _Heads[i]; curr; curr = curr->next()) {
	 os << '#';
	 if (curr->next()) { clump++; }
      }
      if (clump > 3 ) { os << '(' << i << ')'; }
      os << ' ';
   }
   os << ')' << endl;
}

// ----------------------------------------
//  routine to find and return the next key
// ----------------------------------------
template <class K, class T>
bool DictIterator<K,T>::next(K &nextkey) {
   if (_i > _dict._M) { return 0; }

   if (_p != NULL) { _p = _p->next(); }

   while (_p == NULL) {
      if (++_i > _dict._M) { return FALSE; }

      _p = _dict._Heads[_i-1];
   }
   nextkey = _p->key(); // return next key

   return TRUE;
}
