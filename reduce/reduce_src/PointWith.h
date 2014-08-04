//Much of the symmetry code requires an additional point to be carried around
//with each atom. This class simplifies that process and is usable by both
//the sym and nosym versions.

#ifndef _POINTWITH_H_
#define _POINTWITH_H_

#include "Point3d.h"

template <class L>
class PointWith {
private:
     Point3d _point;
     L _atom;

public:
     PointWith() {}
     PointWith(L atom) {
          _atom = atom;
          _point = atom->loc();
     }

     PointWith(L atom, Point3d point) {
          _atom = atom;
          _point = point;
     }

     Point3d getPoint() const {
          return _point;
     }

     L getAtom() const {
          return _atom;
     }

     Point3d& getPoint() {
          return _point;
     }

     L& getAtom() {
          return _atom;
     }
};

#endif
