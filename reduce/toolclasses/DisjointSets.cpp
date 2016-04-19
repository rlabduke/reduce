//Name: DisjointSets.C
//Written By: J. Michael Word
//DateWritten: 11/23/97
//Purpose: manage union/find data structure on a set of integer identifiers
//Based On: Algorithms, Data Structures, and Problem Solving With C++
//          by Mark Allen Weiss, page 740.

#include "DisjointSets.h"

DisjointSets::DisjointSets(int sz) {
   _size   = new int; *_size   = sz;
   _refcnt = new int; *_refcnt =  1;
   _array  = new int [sz];
   for(int i=0; i < sz; i++) {
      _array[i] = -1; // if > 0 contains parent else neg max depth
   }
}

DisjointSets::~DisjointSets() {
   if (--*_refcnt <= 0) {
      delete [] _array;
      delete    _refcnt;
      delete    _size;
   }
}

DisjointSets::DisjointSets(const DisjointSets& ds) {
   (*ds._refcnt)++;
   _array  = ds._array;
   _refcnt = ds._refcnt;
   _size   = ds._size;
}

DisjointSets& DisjointSets::operator=(const DisjointSets& ds) {
   (*ds._refcnt)++;
   if (--*_refcnt <= 0) {
      delete [] _array;
      delete    _refcnt;
      delete    _size;
   }
   _array  = ds._array;
   _refcnt = ds._refcnt;
   _size   = ds._size;
   return *this;
}

// size of a subgroup is encoded in root as negative number
int DisjointSets::size(int x) {
   if ((x >= 0) && (x < *_size)) {
      return(-_array[findRoot(x)]); // keeping neg set size
   }
   return 0; // invalid index
}

int DisjointSets::related(int x, int y) {
   if ((x >= 0) && (x < *_size)
    && (y >= 0) && (y < *_size)) {
      return(findRoot(x) == findRoot(y));
   }
   return 0; // invalid index
}

int DisjointSets::singleton(int x) const {
   if ((x >= 0) && (x < *_size)) {
      return(_array[x] == -1);
   }
   return 0; // invalid index
}

int DisjointSets::id(int x) const {
   if ((x >= 0) && (x < *_size)) {
      return(_array[x] < 0);
   }
   return 0; // invalid index
}

// Array values at root nodes are negative numbers counting
// the size of the tree.

void DisjointSets::rootUnion(int root1, int root2) {
   if ((root1 != root2)
    && (root1 >= 0) && (root1 < *_size)
    && (root2 >= 0) && (root2 < *_size)) {
      if (_array[root2] < _array[root1]) {
	 _array[root2] += _array[root1];   // update size
	 _array[root1] = root2; // make root2 the new root node
      }
      else {
	 _array[root1] += _array[root2];   // update size
	 _array[root2] = root1; // make root1 the new root node
      }
   }
}

int DisjointSets::findRoot(int x) {
   int r = -1;  // return given for invalid index
   if((x >= 0) && (x < *_size)) {
      r = x;
      while(_array[r] >= 0) { r = _array[r]; } // find head of tree
      while(x != r && _array[x] != r) {
	 int t = _array[x];
	 _array[x] = r;    // optimize by compressing tree
	 x = t;
      }
   }
   return r;
}

int DisjointSets::numGroups(int singletons) const {
   int count = 0;
   for(int i=0; i < *_size; i++) {
      if (_array[i] < 0
       && (singletons || _array[i] < -1)) { count++; }
   }
   return count;
}

// build a new table listing the elements of each subset
// entry zero of each subset array is the number of entries
// entries in ss[0] are the subset identifiers
// (   subsets are indexed from 1 to ss[ 0][0])
// (each subset is indexed from 1 to ss[ss][0])
// a call to freeDJsubsets(ss) will free the memory for the table

int** DisjointSets::subsets() {
   int i=0, j=0, k=0;

// first we count the subsets with more than one member
// this step will also make sure findRoot flattens each tree
   int *indx  = new int [*_size];
   int numSS = 0, numSSEntries = 0;
   for(i=0; i < *_size; i++) {
      if (_array[findRoot(i)] < -1) {
	 indx[numSSEntries++] = i;
      }
      if (_array[i] < -1) { numSS++; }
   }

   hsortFlat(indx, numSSEntries);

   int **subsetList = new int *[numSS+1];

   subsetList[0]    = new int[numSS+1];
   subsetList[0][0] = numSS;
   j = 0;
   for(i=0; i < numSSEntries; i++) {
      if (_array[indx[i]] < -1
        && j < numSS) { // sanity check for j
	 subsetList[0][++j] = indx[i];
      }
   }
   k = 0; // offset into the index array
   for(j=1; j <= numSS; j++) {
      const int grpsz = DisjointSets::size(subsetList[0][j]);
      subsetList[j] = new int[grpsz+1];
      subsetList[j][0] = grpsz;
      for(i=1; i <= grpsz; i++) {
	 if (k < numSSEntries) { // sanity check for k
	    subsetList[j][i] = indx[k++];
	 }
      }
   }
   delete [] indx;
   return subsetList;
}

// free up the array of arrays data structure created by subsets
// (including the zeroing out of indexing information)
void freeDJsubsets(int** ss) {
   for(int j=1; j <= ss[0][0]; j++) {
      ss[j][0] = 0;
      delete [] ss[j];
      ss[j] = NULL;
   }
   ss[0][0] = 0;
   delete [] ss[0];
   ss[0] = NULL;
   delete [] ss;
}

// for debugging
void DisjointSets::dumpArray(std::ostream& os) const {
   os << "{ ";
   for(int i=0; i < *_size; i++) {
      os << _array[i] << " ";
   }
   os << "}";
}

inline void swapInts(int& a, int& b) { int temp = a; a = b; b = temp; }

// heapsort
// vec is an array of indexes to (pre-flattened) subset data
void DisjointSets::hsortFlat(int vec[], int len) const {
   // (internal working indexes are 1 based)

   int i;                                 // build heap
   for(i = len >> 1; i >= 1; i--) { // starting in the center
      siftdownFlat(vec, i, len);
   }
   for(i = len; i >= 2; i--) { // put max values at the end
      swapInts(vec[0], vec[i-1]);
      siftdownFlat(vec, 1, i-1);
   }
}

// siftdown is a heap building utility func.
void DisjointSets::siftdownFlat(int vec[], int l, int u) const {
//    (internal working indexes are 1 based)
//  pre-condition: maximal heap(l+1, u)
// post-condition: maximal heap(l,   u)

   int i = l;
   for(;;) { // maxheap(l, u) except between i and its children
      int c = i << 1; // double
      if (c > u) { break; }
      if (c+1 <= u
       && flatRoot(vec[c-1]) < flatRoot(vec[c])) { c++; }
      if (! (flatRoot(vec[i-1]) < flatRoot(vec[c-1]))) {
	 break;
      }
      swapInts(vec[c-1], vec[i-1]);
      i = c;
   }
}
