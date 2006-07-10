
#include "AtomDescr.h"

ostream& operator<<(ostream& os, const AtomDescr& A)
{
		os << "[" << A.getAtomPos() << ", " << A.getAtomResNum() << ", " << A.getAtRadius() << ", "<< A.getOriginalAtomPtr() << "]";
		return os;
}
