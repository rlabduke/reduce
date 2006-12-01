
#include "AtomDescr.h"

std::ostream& operator<<(std::ostream& os, const AtomDescr& A)
{
		os << "[" << A.getAtomPos() << ", " << A.getAtomResNum() << ", " << A.getAtRadius() << ", "<< A.getOriginalAtomPtr() << "]";
		return os;
}
