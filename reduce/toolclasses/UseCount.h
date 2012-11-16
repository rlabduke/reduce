// Name: UseCount.h
// Author: J. Michael Word
// Date Written: 10/25/96
// Purpose: Interface for useage Count class
//          From: Ruminations on C++, Koenig and Moo, Addison Wesley, 1996

#ifndef USECOUNT_H
#define USECOUNT_H 1

#ifndef BOOLPREDEFINED
typedef int bool;
#endif

class UseCount {
public:
   UseCount(): _p(new int(1)) {}
   UseCount(const UseCount& u): _p(u._p) { ++*_p; }
   ~UseCount() { if (--*_p == 0) delete _p; }
   bool only() { return *_p == 1; }
   bool reattach(const UseCount&);
   bool makeonly();
private:
   UseCount& operator=(const UseCount&); // unimplemented to prevent copy
   int* _p;
};

#endif
