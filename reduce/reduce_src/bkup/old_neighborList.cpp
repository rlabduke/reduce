#include "neighborList.h"

namespace cc {
  using namespace cctbx::sgtbx;
  using namespace cctbx::uctbx;
  using namespace cctbx::sgtbx::asu;
  using namespace cctbx::crystal::direct_space_asu;
  using cctbx::crystal::pair_asu_table;
  using cctbx::crystal::pair_sym_table;
  using cctbx::crystal::pair_sym_dict;
  using cctbx::crystal::pair_sym_ops;
};


NeighborList::NeighborList(){_init=false;}

NeighborList::NeighborList(Coord max_range){
	_max_range=max_range;
	_init=false;
}

template <class L>
int NeighborList::init(const std::vector<L> atoms, scitbx::af::double6 cell, char* space_grp){

	_init=true;

	fprintf(stderr,"Received space group type %s\n",space_grp);
	_unit_cell=cc::unit_cell(cell);
	cctbx::sgtbx::space_group spaceGrp=cctbx::sgtbx::space_group(cctbx::sgtbx::space_group_symbols(space_grp));

	// Add sites to _sites_frac

	for(typename std::vector<L >::const_iterator itr=atoms.begin(); itr != atoms.end(); itr++){
		Point3d site=(*itr)->loc();
		sites_frac.push_back( double3(site.x(), site.y(), site.z()) );
	}

	scitbx::af::const_ref< double3 > _sites_frac=sites_frac.const_ref();

	// Construct pair_sym_table
	cc::direct_space_asu metric_free_asu( spaceGrp.type() );
	cc::float_asu<> float_asu = metric_free_asu.as_float_asu(_unit_cell, 1.0E-6);
	boost::shared_ptr< cc::asu_mappings<> > asu_mappings(
	new cc::asu_mappings<> (spaceGrp, float_asu, _max_range) );
	asu_mappings->process_sites_frac(_sites_frac); /// change 0.5 to 0
	cc::pair_asu_table<> pair_asu_table(asu_mappings);
	pair_asu_table.add_all_pairs(_max_range);
	_pair_sym_table = pair_asu_table.extract_pair_sym_table();

	return _pair_sym_table.size();
}

std::list< std::pair<int,Point3d> > NeighborList::get_neighbors(int index, Coord min_range, Coord max_range){
	std::list< std::pair<int, Point3d> > neighbors;

	if (!_init)
		return neighbors;

	if (max_range < 0.0)
		max_range = _max_range;

	scitbx::af::const_ref< double3 > _sites_frac=sites_frac.const_ref();

	const cctbx::fractional<> frac_i = _sites_frac[index];
	const cc::pair_sym_dict pair_sym_dict = _pair_sym_table[index];
	for(cc::pair_sym_dict::const_iterator pair=pair_sym_dict.begin(); pair!=pair_sym_dict.end(); ++pair ){
		const int j = pair->first;
		const cctbx::fractional<> frac_j = _sites_frac[j];
		const cc::pair_sym_ops sym_ops = pair->second;
		for(cc::pair_sym_ops::const_iterator sym_op = sym_ops.begin();sym_op!=sym_ops.end(); ++sym_op){
			const cctbx::fractional<> frac_ji = (*sym_op) * frac_j;
			Coord dist=_unit_cell.distance(frac_i, frac_ji);
			if((dist>min_range)&& (dist < max_range)){
				neighbors.push_back(std::make_pair(j,Point3d(frac_ji[0],frac_ji[1],frac_ji[2])));
			}
		}
	}
	return neighbors;
}

int NeighborList::size(){return _pair_sym_table.size();}
