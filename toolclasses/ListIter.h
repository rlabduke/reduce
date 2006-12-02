// Name: ListIter.h
// Author: J. Michael Word
// Date Written: 8/1/97
// Purpose: List iterator interface (includes non-const list iterator)
// Reference: based on (with extensive changes),
//     An Introduction to Object-Oriented Programming in C++
//     Graham Seed, Springer-Verlag, 1996.

#ifndef LISTITER_H
#define LISTITER_H 1

#include "Iter.h"
#include "List.h"
#ifndef BOOLPREDEFINED
typedef int bool;
#endif

template <class T>
class ListIter : public Iterator<T> {
public:
   ListIter(const List<T>& l) : _node(l._first), _list(l) {}
   ListIter(const ListIter<T>& li)
			      : _node(li._node), _list(li._list) {}

   void sync(const ListIter<T>& li) { _node = li._node; }

   // does current entry exist?
   bool valid()    const { return _node != NULL; }
   operator bool() const { return _node != NULL; }

   int size() const { return _list.size(); }

   // data access and modification
   const T& operator *() const {
      assert(_node);
      return _list.linkData(_node);
   }

   void reset() { _node = _list._first; } // go to the first entry
   void begin() { _node = _list._first; }
   void end()   { _node = _list._last;  } // go to last entry

   bool next(T& e);   // return current entry and advance
   bool prev(T& e);   // return current entry and retreat

   bool operator++();    // advance to next entry (prefix)
   bool operator++(int); // advance to next entry (postfix)
   bool operator--();    // retreat to prev entry (prefix)
   bool operator--(int); // retreat to prev entry (postfix)
protected:
   ListIter<T>& operator=(const ListIter<T>&); // can't assign

   DblLnkLstNode<T>* _node;  // current node
   const List<T>&    _list;  // the list
};

// NonConst-Iterator which modifies with non-const lists
template <class T>
class NonConstListIter : public ListIter<T> {
public:
   NonConstListIter(List<T>& l) : ListIter<T>(l), _NClist(l) {}
   NonConstListIter(const NonConstListIter<T>& li)
			      : ListIter<T>(li), _NClist(li._NClist) {}

   T& data() const {
      assert(_node);
      return _NClist.linkData(_node);
   }

   void update(const T& e) const {
      assert(_node);
      _NClist.linkData(_node) = e;
   }

   bool drop(); // remove the current element

   bool insertBefore(const T& e);
   bool insertAfter (const T& e);

protected:
   NonConstListIter<T>& operator=(const NonConstListIter<T>&);

   List<T>& _NClist;  // con-constant reference to the list
};

#ifdef INCTEMPLATEDEFNS
#include "ListIter.cpp"
#endif

#endif
