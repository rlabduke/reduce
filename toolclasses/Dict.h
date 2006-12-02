// name: Dict.h
// author: J. Michael Word
// date written: 8/7/97
// purpose: Interface file for a dictionary class
//          which maps (single) objects to keys.

// KEYS: a key type needs the following properties:
// it must be copyable (a const copy is made)
// and it should be assignable for the iterator to work
// it must have valid == and < operators
// and it must have a hash function like:
//     unsigned long hash(const KEYTYPE& k, unsigned long M);

// VALUES: a value type needs to be copyable
// or else a ptr can be passed (in this case the dict will free).

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#ifndef DICT_H
#define DICT_H 1

#include "Hash.h"
#include "utility.h"

#include <iostream>

// simple single-link list of key/value-seq pairs

template <class K, class T> class DictLink {
public:
   DictLink()                  : _next(NULL), _value(NULL) {}
   DictLink(const K& k, const T& v)
                               : _next(NULL), _key(k), _value(new T(v)) {}
   DictLink(const K& k, T* vp) : _next(NULL), _key(k), _value(vp) {}

   ~DictLink() { delete _value; delete _next; _next = NULL; }

   const K& key() const { return _key;   } // access methods
       T* value() const { return _value; }

   DictLink<K,T>* next() const { return _next; }

   void value(const T& v) {
      if (_value) {delete _value;};
      _value = new T(v);
   }
   void value(T* vp)      {
      if (_value) {delete _value;};
      _value = vp;
   }

   DictLink<K,T> * next(DictLink<K,T> *nxt) { // update link
      DictLink<K,T> * prev = _next;
      _next = nxt;
      return prev;  //NOTE: returns the old next pointer!
   }
private:
   DictLink(const DictLink&);            // can't copy
   DictLink& operator=(const DictLink&); // can't assign
   
   DictLink<K, T> *_next;
   const K         _key;
   T*	           _value;
};

template <class K, class T> class DictIterator;

// -----------------------------------------
//  dictionary maps value (pointers) to keys
// -----------------------------------------
template <class K, class T> class Dict {
friend class DictIterator<K,T>;
public:
   Dict(const unsigned long sz=100);
   virtual ~Dict();

   bool put(const K& key, const T& value);
   bool put(const K& key, T* value);

   bool get(const K &key, T& retval) const;
   T*   get(const K &key) const;

   bool exists(const K& key) const;

   bool drop(const K& key);
   void reset();
   
   unsigned long size() const { return _n; } // num of keys

   void dumpStructure(std::ostream& os) const;
private:
   Dict(const Dict&);            // can't copy
   Dict& operator=(const Dict&); // can't assign
   
   unsigned long _M;       // table size (a prime number)
   unsigned long _n;       // current number of entries
   DictLink<K,T> ** _Heads; // array of entry pointers
};

// -----------------------------------------
//  iterate through the keys in a dictionary
// -----------------------------------------
template <class K, class T> class DictIterator {
public:
   DictIterator(const Dict<K,T> &d) : _dict(d), _i(0), _p(NULL) {}
virtual ~DictIterator() {}

   bool next(K& nextkey);

   void reset() { _i = 0; _p = NULL; }
private:
   const Dict<K,T> &_dict;// the dictionary we are looking at
   unsigned long    _i;   // which head were we on (1 .. _M)
   DictLink<K,T>*   _p;   // _p->next() is the next entry to consider
};

#ifdef INCTEMPLATEDEFNS
#include "Dict.cpp"
#endif

#endif
