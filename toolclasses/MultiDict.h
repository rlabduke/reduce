// name: MultiDict.h
// author: J. Michael Word
// date written: 8/7/97
// purpose: Interface file for a dictionary class
//          which maps (multiple) objects to keys.

// KEYS: a key type needs the following properties:
// it must be copyable (a const copy is made)
// and it should be assignable for the iterator to work
// it must have valid == and < operators
// and it must have a hash function like:
//     unsigned long hash(const KEYTYPE& k, unsigned long M);

// VALUES: a value type needs only to be copyable
// a sequence of values is maintained for each key

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#ifndef MULTIDICT_H
#define MULTIDICT_H 1

#include "Hash.h"
#include "Seq.h"
#include "utility.h"

#include <iostream>

// simple single-link list of key/value-seq pairs

template <class K, class T> class MultiLink {
public:
   MultiLink()           : _next(NULL)          {}
   MultiLink(const K& k) : _next(NULL), _key(k) {}

   ~MultiLink() { delete _next; _next = NULL; }

   const K& key() const { return _key;   } // access methods
   Seq<T> value() const { return _value; }

   MultiLink<K,T>* next() const { return _next; }

   void insert(const T& v) { _value.insert(v); }
   int count() const { return _value.len(); }

   MultiLink<K,T> * next(MultiLink<K,T> *nxt) { // update link
      MultiLink<K,T> * prev = _next;
      _next = nxt;
      return prev;  //NOTE: returns the old next pointer!
   }
private:
   MultiLink(const MultiLink&);            // can't copy
   MultiLink& operator=(const MultiLink&); // can't assign
   
   MultiLink<K, T> *_next;
   const K          _key;
   Seq<T>	    _value;
};

template <class K, class T> class MultiDictIterator;

// -----------------------------------------
//  dictionary maps value (pointers) to keys
// -----------------------------------------
template <class K, class T> class MultiDict {
friend class MultiDictIterator<K,T>;
public:
   MultiDict(const unsigned long sz=100);
   virtual ~MultiDict();

   int put(const K& key, const T& value);
   Seq<T> get(const K& key) const;
   bool exists(const K& key) const;

   bool drop(const K& key);
   void reset();
   
   unsigned long size() const { return _n; } // num of keys

   void dumpStructure(std::ostream& os) const;
private:
   MultiDict(const MultiDict&);            // can't copy
   MultiDict& operator=(const MultiDict&); // can't assign
   
   unsigned long _M;       // table size (a prime number)
   unsigned long _n;       // current number of entries
   MultiLink<K,T> ** _Heads; // array of entry pointers
};

// -----------------------------------------
//  iterate through the keys in a dictionary
// -----------------------------------------
template <class K, class T> class MultiDictIterator {
public:
   MultiDictIterator(const MultiDict<K,T> &d)
				 : _dict(d), _i(0), _p(NULL) {}
virtual ~MultiDictIterator() {}

   int next(K& nextkey);

   void reset() { _i = 0; _p = NULL; }
private:
   const MultiDict<K,T> &_dict;// the dictionary we are looking at
   unsigned long        _i;   // which head were we on (1 .. _M)
   MultiLink<K,T>*        _p;   // _p->next() is the next entry to consider
};

#ifdef INCTEMPLATEDEFNS
#include "MultiDict.cpp"
#endif

#endif
