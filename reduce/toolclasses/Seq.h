// Name: Seq.h
// Author: J. Michael Word
// Date Written: 10/24/96
// Purpose: Interface for Seq and Seq_item classes
//          which are used to manage a sequence of objects.
//          From: Ruminations on C++, Koenig and Moo, Addison Wesley, 1996
#ifndef SEQ_H
#define SEQ_H 1

#ifdef OLD_STD_HDRS
#include <assert.h>
#else
#include <cassert>
#endif

#ifndef BOOLPREDEFINED
typedef int bool;
#endif
#ifdef BRACKETOPERPARMS
#define PARM_TYPE_T <T>
#else
#define PARM_TYPE_T
#endif

template <class T> class Seq;

template <class T> class Seq_item {
   friend class Seq<T>;

public:
   const T& data() const { return _data; }
private:
   int _use;
   const T _data;
   Seq_item* _next;

   friend int  operator== PARM_TYPE_T(const Seq_item&, const Seq_item&);
   friend int  operator!= PARM_TYPE_T(const Seq_item&, const Seq_item&);

   Seq_item(const T& t, Seq_item* s): _use(1), _data(t),
                                      _next(s) { if (s) s->_use++; }
   Seq_item(const T& t): _use(1), _data(t), _next(0) { }
};

template <class T> class Seq {
public:
   Seq(): _item(0), _len(0) { }
   Seq(const T& t, const Seq& x): _item(new Seq_item<T>(t, x._item)),
                                  _len(x._len + 1) { }
   Seq(const Seq& s): _item(s._item), _len(s._len) { if (_item) _item->_use++; }
   ~Seq() { destroy(_item); }
   Seq& operator=(const Seq&);
   Seq& operator++();
   Seq  operator++(int);
   const T& operator*() const {
      assert(_item != 0);
      return _item->_data;
   }
   T hd() const {
      assert(_item != 0);
      return _item->_data;
   }
   Seq tl() const {
      assert(_item != 0);
      return Seq<T>(_item->_next, _len-1);
   }
   unsigned len() const { return _len; }
   operator bool() const { return _item != 0; }
   int     empty() const { return _item == 0; }
   Seq& insert(const T& t) {
      if (_item) _item->_use--;                   // added!!
      _item = new Seq_item<T>(t, _item);
      _len++;
      return *this;
   }
   Seq& rev();
   Seq& operator+=(Seq);

   friend int operator== PARM_TYPE_T(const Seq&, const Seq&);
   friend int operator!= PARM_TYPE_T(const Seq&, const Seq&);
   friend Seq operator + PARM_TYPE_T(const Seq&, const Seq&);

   Seq_item<T>* start() const { return _item; } // JSS: iterator that doesn't ref count (I hope) 

   // next() has been added to make == and other friend functions workable
	// JSS: next() made public for iterator
   static Seq_item<T>* next(Seq_item<T>* it) { return (it?it->_next:0); }
private:
   Seq(Seq_item<T>* s, unsigned l): _item(s), _len(l) { if (s) s->_use++; }

   void destroy(Seq_item<T>*);
   Seq_item<T>* owntail();

   Seq_item<T>* _item;
   unsigned _len;
};

template <class T>
int operator==(const Seq_item<T>& op1, const Seq_item<T>& op2);

template <class T>
int operator!=(const Seq_item<T>& op1, const Seq_item<T>& op2);

template <class T>
int operator +(const Seq_item<T>& op1, const Seq_item<T>& op2);

template <class T>
inline Seq<T> cons(const T& t, const Seq<T>& s) { return Seq<T>(t, s); }

template <class T> int length(Seq<T>); // for debugging only

template <class T> Seq<T> sort(const Seq<T>& x);

template <class T> Seq<T> merge(Seq<T>, Seq<T>);
template <class T> void split(Seq<T>, Seq<T>&, Seq<T>&);
template <class T> Seq<T> mergesort(Seq<T>&);

#ifdef INCTEMPLATEDEFNS
#include "Seq.cpp"
#endif

#endif
