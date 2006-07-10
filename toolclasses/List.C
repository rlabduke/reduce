// Name: List.C
// Author: J. Michael Word
// Date Written: 8/1/97
// Purpose: doubly linked list implementation
// Reference: based on (with extensive changes),
//     An Introduction to Object-Oriented Programming in C++
//     Graham Seed, Springer-Verlag, 1996.

#include "List.h"

// returns list index of a given object; else -1
template <class T>
int List<T>::search(const T& obj) const {
   int i = 0;
   for (DblLnkLstNode<T>* ptr(_first); ptr; ptr=ptr->_next) {
      if (ptr->_data == obj) { return(i); }
      i++;
   }
   return(-1);
}

// inserts an object into the list before specified the location
template <class T>
void List<T>::insertBefore(const T& obj, int loc) {

   if (_first == NULL) { // empty list
      _first = _last = new DblLnkLstNode<T>(obj);
      _num_elem = 1;
   }
   else {
      if (loc == List<T>::START || loc == 0) {
	 insertBefore(obj, _first);
      }
      else if (loc == List<T>::END || loc >= _num_elem) { 
	 insertAfter(obj, _last);
      }
      else {
	 if (loc >= 0 || loc < _num_elem) {
	    DblLnkLstNode<T>* ptr(_first);

	    for (int i=0; i < loc; i++) {
	       ptr = ptr->_next;
	    }
	    insertBefore(obj, ptr);
	 }
      }
   }
}

// places an object at head of list
template <class T>
void List<T>::prepend(const T& obj) {
   insertBefore(obj, List<T>::START);
}

// appends an object to list
template <class T>
void List<T>::append(const T& obj){
   insertBefore(obj, List<T>::END);
}

// removes an object from a given position [0 .. size()-1]
template <class T>
bool List<T>::drop(int loc) {
   bool successFlag = 0;
   if (_first) { // i.e. non empty
           if (loc == List<T>::START
	    || loc == 0)            { successFlag = drop(_first); }
      else if (loc == List<T>::END
            || loc == _num_elem-1 ) { successFlag = drop(_last);  }
      else {
	 if (loc >= 0 || loc < _num_elem) {
	    DblLnkLstNode<T>* ptr(_first);

	    for (int i=0; i < loc; i++) {
	       ptr = ptr->_next;
	    }
	    successFlag = drop(ptr);
	 }
      }
   }
   return successFlag;
}

// sorts a List object in ascending order
template <class T>
void List<T>::sort () {
   DblLnkLstNode<T>* p = mergesort(_first); // uses just the _next links
   _first = p;

   if (p) {    // now build back the _prev links
      while (p->_next) {
	 p->_next->_prev = p;
	 p = p->_next;
      }
   }
   _last = p;
}

template <class T>
void List<T>::rev() {
   DblLnkLstNode<T> *p = NULL;
   while (_first) { addToHead(&p, &_first); }

   _first = p;

   if (p) {    // now build back the _prev links
      while (p->_next) {
	 p->_next->_prev = p;
	 p = p->_next;
      }
   }
   _last = p;
}

// unlinks and deletes a node
template <class T>
bool List<T>::drop(DblLnkLstNode<T> *ptr) {

   if (ptr) {
      if (ptr->_next) { ptr->_next->_prev = ptr->_prev; }
      if (ptr->_prev) { ptr->_prev->_next = ptr->_next; }
      if (_first == ptr) { _first = ptr->_next; }
      if (_last  == ptr) { _last  = ptr->_prev; }
      _num_elem--;
      ptr->_next = NULL;
      ptr->_prev = NULL;
      delete ptr;
      return 1;
   }
   return 0;
}

template <class T>
void List<T>::insertBefore(const T& obj, DblLnkLstNode<T>*ptr) {
   if (ptr) {
      DblLnkLstNode<T>* x = new DblLnkLstNode<T>(obj);
      x->_next = ptr;
      if (ptr->_prev) {
	 ptr->_prev->_next = x;
	 x->_prev = ptr->_prev;
      }
      ptr->_prev = x;
      if (_first == ptr) { _first = x; }
      _num_elem++;
   }
}

template <class T>
void List<T>::insertAfter(const T& obj, DblLnkLstNode<T>*ptr) {
   if (ptr) {
      DblLnkLstNode<T>* x = new DblLnkLstNode<T>(obj);
      x->_prev = ptr;
      if (ptr->_next) {
	 ptr->_next->_prev = x;
	 x->_next = ptr->_next;
      }
      ptr->_next = x;
      if (_last == ptr) { _last = x; }
      _num_elem++;
   }
}

// assignment operator
template <class T>
List<T>& List<T>::operator = (const List& l) {

   if (this != &l) {
      delete _first;   // delete existing nodes...

      _first = NULL;   // tidy up for good measure
      _last  = NULL;
      _num_elem = 0;

      dupList(l);  // build new list from old
   }
   return *this;
}

// subscript operator - const
template <class T>
const T& List<T>::operator [] (int index) const {
   if      (index == List<T>::START) { index = 0; }
   else if (index == List<T>::END)   { index = _num_elem - 1; }
   assert(index >= 0 && index < _num_elem) ;

   DblLnkLstNode<T>* loc(_first);

   for (int i=0; i < index; i++) { loc = loc->_next; }
   return loc->_data;
}

// subscript operator - non-const
template <class T>
T& List<T>::operator [] (int index) {
   if      (index == List<T>::START) { index = 0; }
   else if (index == List<T>::END)   { index = _num_elem - 1; }
   assert(index >= 0 && index < _num_elem) ;

   DblLnkLstNode<T>* loc(_first);

   for (int i=0; i < index; i++) { loc = loc->_next; }
   return loc->_data;
}

// equality operator
template <class T>
bool List<T>::operator == (const List& l) const {
   if (_num_elem != l._num_elem) { return FALSE; }

   DblLnkLstNode<T>* ptr(_first);      // initialise pointers to each
   DblLnkLstNode<T>* l_ptr(l._first);

   for (int i=0; i<_num_elem; i++, ptr=ptr->_next, l_ptr=l_ptr->_next) {
      if (! (ptr->_data == l_ptr->_data)) {
	 return FALSE;
      }
   }
   return TRUE;
}

// inequality operator
template <class T>
bool    List<T>::operator != (const List& l) const {
   return !(*this == l);
}

// less-than operator
template <class T>
bool List<T>::operator < (const List& l) const {
   if (_num_elem != l._num_elem) { return (_num_elem < l._num_elem); }

   DblLnkLstNode<T>* ptr(_first);      // initialise pointers to each
   DblLnkLstNode<T>* l_ptr(l._first);

   for (int i=0; i<_num_elem; i++, ptr=ptr->_next, l_ptr=l_ptr->_next) {
      if (! (ptr->_data == l_ptr->_data)) {
	 return (ptr->_data < l_ptr->_data);
      }
   }
   return FALSE;
}

// concatenate
template <class T>
List<T>& List<T>::operator+=(const List& l) {
   int n = l.size();
   DblLnkLstNode<T>* l_ptr(l._first);

   for (int i=0; i < n; i++) {
      insertAfter(l_ptr->_data, _last);
      l_ptr=l_ptr->_next;
   }
   return *this;
}

// output operator
template <class T>
ostream& operator << (ostream& s, const List<T>& l) {
   DblLnkLstNode<T>* ptr(l._first);  // point to first node in list

   s << "(" ;
   for (int i=0; i < l._num_elem; i++, ptr=l.linkNext(ptr)) {
      s << l.linkData(ptr) << ((i != l._num_elem - 1) ? ", " : "");
   }
   s << ")";
   return s;
}

// private method used in copy constructor and assignment operator
template <class T>
void List<T>::dupList(const List& l) {
   _num_elem = 0;

   DblLnkLstNode<T>* old_node(l._first);  // old list
   DblLnkLstNode<T>* new_node(NULL);      // new list

   if (old_node) {            // copy first node
      new_node = new DblLnkLstNode<T>(old_node->_data);
      _num_elem++;
      old_node = old_node->_next;
      _first = new_node;
   }
   else { _first = NULL; }
                              // copy remaining nodes
   while (old_node) {
      new_node->_next = new DblLnkLstNode<T>(old_node->_data, new_node);
      _num_elem++;
      old_node = old_node->_next;
      new_node = new_node->_next;
   }
   _last = new_node;
}

// sort work routines
// they manipulate singly linked lists
// the full sort routine will build the backward links

// remove singly linked item from head of one list to head of another
template <class T>
DblLnkLstNode<T>* List<T>::addToHead(DblLnkLstNode<T> **toList, DblLnkLstNode<T> **fromList) {
   if (*fromList) {
      DblLnkLstNode<T> *next = (*fromList)->_next;
      (*fromList)->_next = (*toList);
      (*toList) = (*fromList);
      (*fromList) = next;
   }
   return (*toList);
}

// merge two sorted singly linked lists into one sorted list
template <class T>
DblLnkLstNode<T>* List<T>::merge(DblLnkLstNode<T> *x, DblLnkLstNode<T> *y) {
   DblLnkLstNode<T> *lst=NULL;
   while (x && y) {
      if (x->_data < y->_data) { addToHead(&lst, &x); }
      else { addToHead(&lst, &y); }
   }
   while (x) { addToHead(&lst, &x); }
   while (y) { addToHead(&lst, &y); }

   DblLnkLstNode<T> *revLst = NULL;	// now we reverse the list
   while (lst) { addToHead(&revLst, &lst); }
   return revLst;
}

// split singly linked list evenly into two lists
template <class T>
void List<T>::split(DblLnkLstNode<T> *x, DblLnkLstNode<T> **y, DblLnkLstNode<T> **z) {
   while(x) {
      addToHead(y, &x);
      if (x) { addToHead(z, &x); }
   }
}

// sort a singly linked list
template <class T>
DblLnkLstNode<T>* List<T>::mergesort(DblLnkLstNode<T> *x) {
   DblLnkLstNode<T> *p=NULL, *q=NULL;

   if (x == NULL || x->_next == NULL) { return x; }
   split(x, &p, &q);
   return merge(mergesort(p), mergesort(q));
}
