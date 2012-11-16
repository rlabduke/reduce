#include <cstdio>
#include <list>
#include <cctbx/sgtbx/direct_space_asu/proto/direct_space_asu.h>
#include <boost/shared_ptr.hpp> // necessary for pair_tables.h
#include <cctbx/crystal/pair_tables.h>
#include "LocBlk.h"
#include <utility>
//#include "PDBrec.h"
typedef scitbx::vec3<double> double3;
//typedef PDBrec* L;

class NeighborList{
private:
	cctbx::crystal::pair_sym_table _pair_sym_table;
	Coord _max_range;
	bool _init;
//	scitbx::af::const_ref< double3 > _sites_frac;
	scitbx::af::shared< double3 > sites_frac;
	cctbx::uctbx::unit_cell _unit_cell;

public:

	NeighborList();
	NeighborList(Coord max_range);

	template <class L> 
	int init(const std::vector<L> atoms, scitbx::af::double6 cell, char* space_grp);

	std::list< std::pair<int,Point3d> > get_neighbors(int index, Coord min_range, Coord max_range=-1.0f);

	int size();
};
