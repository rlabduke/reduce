// Name: UseCount.h
// Author: J. Michael Word
// Date Written: 10/25/96
// Purpose: Implementation for useage Count class
//          From: Ruminations on C++, Koenig and Moo, Addison Wesley, 1996
#include "UseCount.h"

bool UseCount::reattach(const UseCount& u) {
   bool rc = 0;
   ++*u._p;
   if (--*_p == 0) {  // <-- fixed bug in book
      delete _p;
      rc = 1;
   }
   _p = u._p;
   return rc;
}

bool UseCount::makeonly() {
   if (*_p == 1) { return 0; }
   --*_p;
   _p = new int(1);
   return 1;
}
