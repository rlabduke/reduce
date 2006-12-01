
#ifndef  _HACK_SORT_CLASS_FOR_PDBRECS_H_
#define  _HACK_SORT_CLASS_FOR_PDBRECS_H_

//Lightweight class to represent an atom 
//Meant to be sorted, compared, and easily allocated.
//put in an interface to this class to answer my questions easily
//and with a lot less memory required.

#include "Point3d.h"
#include <iostream>

class PDBrec;

class AtomDescr
{
public:
	AtomDescr() : _atomPos(0,0,0), _atomResNo(0), _atRadius(0), _atRadiusSqr(0), _originalAtom( 0 ) {};
	AtomDescr(
		Point3d atomPos,
		int resNum,
		float atRadius
	) : 
		_atomPos(atomPos), _atomResNo(resNum),
		_atRadius(atRadius),
		_atRadiusSqr( atRadius*atRadius),
		_originalAtom( 0 )
	{};
	
	~AtomDescr() {};
	
	AtomDescr(const AtomDescr& rhs) {
		_atomPos = rhs._atomPos;
		_atomResNo = rhs._atomResNo;
		_atRadius = rhs._atRadius; 
		_atRadiusSqr = rhs._atRadiusSqr;
		_originalAtom = rhs._originalAtom;
	};
	
	void setAtomPos(Point3d sn) {_atomPos = sn; return;}
	
	void setAtomResNum(int rn) {_atomResNo = rn; return;};
	
	void setAtRadius(float atRadius) {_atRadius = atRadius; _atRadiusSqr = atRadius*atRadius; return;};
	Point3d  getAtomPos() const {return _atomPos;};
	int      getAtomResNum() const {return _atomResNo;};
	float   getAtRadius() const {return _atRadius;} ;
	
	bool operator==(const AtomDescr& rhs) const
	{
		return ((_atomResNo == rhs._atomResNo) && 
			(_atRadius == rhs._atRadius) && 
			(distanceSquared(_atomPos, rhs._atomPos) < .000001 ) );
	};

	bool operator!=(const AtomDescr& rhs) const
        {
          return !operator==(rhs);
        }
	
	bool operator < (const AtomDescr& rhs) const
	{
		if (_atomResNo == rhs._atomResNo)
		{
			if (_atRadius == rhs._atRadius)
			{  
				if (_atomPos.x() == rhs._atomPos.x() )
				{
					if (_atomPos.y() == rhs._atomPos.y() )
					{
						if (_atomPos.z() == rhs._atomPos.z() )
						{
							return false;
						}
						else
						{
							return _atomPos.z() < rhs._atomPos.z();
						}
					}
					else
					{
						return _atomPos.y() < rhs._atomPos.y();
					}
				}
				
				else
				{
					return _atomPos.x() < rhs._atomPos.x();
				}
			}
			else
			{
				return _atRadius < rhs._atRadius;
			}
		}
		else
		{
			return _atomResNo < rhs._atomResNo;
		}
	}
	bool intersects(const AtomDescr& a) const
	{
		return (distanceSquared(a._atomPos, _atomPos) < (a._atRadius + _atRadius) * (a._atRadius + _atRadius));
	}
	
	void setOriginalAtomPtr(PDBrec* originalAtom) {_originalAtom = originalAtom;}
	
	PDBrec* getOriginalAtomPtr() const {return _originalAtom;}
	
	bool containsPoint( Point3d const & p ) const{ return (distanceSquared( p , _atomPos) < _atRadiusSqr );}
private:
	Point3d  _atomPos;
	int  		_atomResNo;
	float   _atRadius;
	float _atRadiusSqr;
	PDBrec* _originalAtom;
};

std::ostream& operator<<(std::ostream& os, const AtomDescr& A);

#endif //_HACK_SORT_CLASS_FOR_PDBRECS_H_

