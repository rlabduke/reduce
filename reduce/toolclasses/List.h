// Name: List.h
// Author: J. Michael Word
// Date Written: 8/1/97
// Purpose: doubly linked list interface
// Reference: based on (with extensive changes),
//     An Introduction to Object-Oriented Programming in C++
//     Graham Seed, Springer-Verlag, 1996.

#ifndef LIST_H
#define LIST_H 1

#include <iostream>
#ifdef OLD_STD_HDRS
#include <limits.h>
#include <assert.h>
#else
#include <climits>
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

template <class T> class List;
template <class T> class ListIter;
template <class T> class NonConstListIter;

template <class T>
class DblLnkLstNode {
   friend class List<T>;
private:
   DblLnkLstNode() : _data(), _prev(NULL), _next(NULL) {}
   DblLnkLstNode(const T& d,
                 DblLnkLstNode* p=NULL,
                 DblLnkLstNode* n=NULL) : _data(d), _prev(p), _next(n) {}

   DblLnkLstNode(const DblLnkLstNode&);              // can't copy
   DblLnkLstNode& operator = (const DblLnkLstNode&); // can't assign
	// JSS: deleting connected nodes is the responsiblity the List
public:

private:
   T                 _data; // stored data value
   DblLnkLstNode<T> *_prev;
   DblLnkLstNode<T> *_next;
};

template <class T>
class List {
   friend class ListIter<T>;
   friend class NonConstListIter<T>;
public:
   enum { START = 0, END = INT_MAX };

   List() : _num_elem(0), _first(NULL), _last(NULL) {}
   List(const List& l) { dupList(l); }
  ~List() { 
    while (_first != NULL) {
      DblLnkLstNode<T> *entry = _first;
      _first = _first->_next;
      delete entry;
    }
    _num_elem = 0 ; 
    _first = _last = NULL; 
  }

   int  size() const { return _num_elem; }
   bool empty() const { return _num_elem == 0; }

   int  search(const T& obj) const;

   void insertBefore(const T& obj, int loc = START);
   void prepend(const T& obj);
   void append(const T& obj);

   bool drop(int loc); // remove the entry at loc [0..size()-1]

   void sort();
   void rev(); // reverse

   const T& operator*() const { assert(_num_elem > 0); return _first->_data; }

   List& operator = (const List& l); // assignment

   const T& operator [] (int index) const;
         T& operator [] (int index);   // data is accessible

   bool operator == (const List& l) const; // comparisons
   bool operator != (const List& l) const;
   bool operator <  (const List& l) const;

   List& operator += (const List& l); // concatenate

   friend std::ostream& operator << PARM_TYPE_T (std::ostream& s, const List<T>& l);
   friend std::ostream& outputRecords(std::ostream& os, const List<T>& l);

#ifndef LISTFRIENDFIX
private:
#endif

   int _num_elem; // num nodes

   DblLnkLstNode<T> *_first;
   DblLnkLstNode<T> *_last;

   void dupList(const List& l);  // used in copy constructor and assignment
   bool drop(DblLnkLstNode<T>*); // unlink and delete a node
   void insertBefore(const T& obj, DblLnkLstNode<T>*);
   void insertAfter (const T& obj, DblLnkLstNode<T>*);

   DblLnkLstNode<T> * mergesort(DblLnkLstNode<T> *x); // sort utilities
   DblLnkLstNode<T> * merge(DblLnkLstNode<T> *x, DblLnkLstNode<T> *y);
   void split(DblLnkLstNode<T> *x, DblLnkLstNode<T> **y, DblLnkLstNode<T> **z);
   DblLnkLstNode<T> * addToHead(DblLnkLstNode<T> **toList, DblLnkLstNode<T> **fromList);

   // private friend access to link data
   T& linkData(DblLnkLstNode<T> *ptr) const {
      return ptr->_data;
   }
   DblLnkLstNode<T> * linkNext(DblLnkLstNode<T> *ptr) const {
      return ptr->_next;
   }
   DblLnkLstNode<T> * linkPrev(DblLnkLstNode<T> *ptr) const {
      return ptr->_prev;
   }
};

template <class T> std::ostream& operator << (std::ostream& s, const List<T>& l);

#ifdef INCTEMPLATEDEFNS
#include "List.cpp"
#endif

#endif
