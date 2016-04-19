//Name: DisjointSets.h
//Written By: J. Michael Word
//DateWritten: 11/23/97
//Purpose: manage union/find data structure on a set of integer identifiers
//Based On: Algorithms, Data Structures, and Problem Solving With C++
//          by Mark Allen Weiss, page 740.

#ifndef DISJOINTSETS_H
#define DISJOINTSETS_H 1

#include <iostream>

class DisjointSets {
public:
   DisjointSets(int sz);
  ~DisjointSets();
   DisjointSets(const DisjointSets& ds);            // reference copy
   DisjointSets& operator=(const DisjointSets& ds); //reference assign

   void connect(int x, int y) {
      if (x != y) { rootUnion(findRoot(x), findRoot(y)); }
   }
   int operator[](int x) { return findRoot(x); }
   int size() const { return *_size; }
   int size(int x);             // size of set containing x
   int numGroups(int singletons=0) const;

   int related(int x, int y);   // are x and y related?
   int singleton(int x) const;  // is x unrelated to others?
   int id(int x) const;         // is x an identifier?

   int** subsets(); // create an array of subset vectors

   // workhorse routines
   void rootUnion(int root1, int root2);
   int findRoot(int x);
   
   // debugging
   void dumpArray(std::ostream& os) const;
private:
   int *_array;
   int *_refcnt;
   int *_size;

   // utilities used in subsets
   void hsortFlat(int vec[], int len) const;
   void siftdownFlat(int vec[], int l, int u) const;
   int flatRoot(int x) const;
};

// this finds the root assuming a pre-flattened array
inline int DisjointSets::flatRoot(int x) const {
   return (_array[x] < 0) ? x : _array[x];
}

void freeDJsubsets(int**); // clean up after one of our subset arrays
#endif
