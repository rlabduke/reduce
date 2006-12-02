// Name: Seq.C
// Author: J. Michael Word
// Date Written: 10/24/96
// Purpose: Implementation for Seq and Seq_item classes
//          which are used to manage a sequence of objects.
//          From: Ruminations on C++, Koenig and Moo, Addison Wesley, 1996
#include "Seq.h"

template <class T>
int operator==(const Seq_item<T>& op1, const Seq_item<T>& op2) {
   return (op1._data == op2._data);
}

template <class T>
int operator!=(const Seq_item<T>& op1, const Seq_item<T>& op2) {
   return !(op1._data == op2._data);
}
template <class T>
Seq<T>& Seq<T>::operator=(const Seq<T>& s) {
   if (s._item) { s._item->_use++; }
   destroy(_item);
   _item = s._item;
   _len  = s._len;
   return *this;
}

template <class T>
Seq<T>& Seq<T>::operator++() {
   if (_item) {
      Seq_item<T>* p = _item->_next;
      if (--_item->_use == 0) { delete _item; }
      else if (p) { p->_use++; }            // <--- modification
      _item = p;
      _len--;
   }
   return *this;
}

template <class T>
Seq<T> Seq<T>::operator++(int) {
   Seq<T> ret(*this);
   if (_item) {
      --_item->_use;
      _item = _item->_next;
      if (_item) { _item->_use++; }
      _len--;
   }
   return ret;
}

template <class T>
Seq<T>& Seq<T>::rev() {
   if (_item) {
      Seq_item<T>* k = owntail();
      Seq_item<T>* curr = _item;
      Seq_item<T>* behind = 0;
      do {
         Seq_item<T>* ahead = curr->_next;
         curr->_next = behind;
         behind = curr;
         curr = ahead;
      } while(curr);
      _item = k;
   }
   return *this;
}

template <class T>
Seq<T>& Seq<T>::operator+=(Seq<T> s) {
   if (s) {
      Seq_item<T>* final = owntail();
      if (s._item) { s._item->_use++; }
      if (final) { final->_next = s._item; }
      else       { _item        = s._item; }
      _len  += s._len;
   }
   return *this;
}

template <class T>
int operator==(const Seq<T>& op1, const Seq<T>& op2) {
   if (op1._len != op2._len) { return 0; }
   Seq_item<T>* p = op1._item;
   Seq_item<T>* q = op2._item;
   while(p != q) {
      assert(p != 0 && q != 0);
      if (*p != *q) { return 0; }
      p = Seq<T>::next(p);        // <-- fixed bug
      q = Seq<T>::next(q);
   }
   return 1;
}

template <class T>
int operator!=(const Seq<T>& op1, const Seq<T>& op2) {
   return !(op1 == op2);
}

template <class T>
Seq<T> operator+(const Seq<T>& op1, const Seq<T>& op2) {
   Seq<T> r(op1);
   return r += op2;
}

template <class T>
void Seq<T>::destroy(Seq_item<T>* it) { // NOTE: does not set item or len
   while(it && --it->_use == 0) {
      Seq_item<T>* next = it->_next;
      delete it;
      it = next;
   }
}

template <class T>
Seq_item<T>* Seq<T>::owntail() {
   if (_item == 0) { return 0; }
   Seq_item<T>* i = _item;
   Seq_item<T>** p = &_item;
   while(i->_use == 1) {
      if (i->_next == 0) { return i; }
      p = &i->_next;
      i = i->_next;
   }
   *p = new Seq_item<T>(i->_data);
   --i->_use;
   i = i->_next;
   Seq_item<T>* j = *p;
   while(i) {
      j->_next = new Seq_item<T>(i->_data);
      i = i->_next;
      j = j->_next;
   }
   return j;
}


// for debugging only
template <class T>
int length(Seq<T> s) {
   int n = 0;
   while (s++) { n++; }
   return n;
}

template <class T>
Seq<T> sort(const Seq<T>& x) {
   Seq<T> xx = x;
   return mergesort(xx);
}

// function templates for sorting and merging
template <class T>
Seq<T> merge(Seq<T> x, Seq<T> y) {
   Seq<T> r;
   while (x && y) {
      if (x.hd() < y.hd()) {
         r.insert(x.hd());
         x++;
      } else {
         r.insert(y.hd());
         y++;
      }
   }
   while (x) {
      r.insert(x.hd());
      x++;
   }
   while (y) {
      r.insert(y.hd());
      y++;
   }
   return r.rev();
}

template <class T>
void split(Seq<T> x, Seq<T>& y, Seq<T>& z) {
   while(x) {
      y.insert(x.hd());
      if (++x) {
         z.insert(x.hd());
         x++;
      }
   }
}

template <class T>
Seq<T> mergesort(Seq<T>& x) {
   if (!x || !x.tl()) { return x; }
   Seq<T> p, q;
   split(x, p, q);
   return merge(mergesort(p), mergesort(q));
}
