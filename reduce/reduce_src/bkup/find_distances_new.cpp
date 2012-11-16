#include <cstdio>

#include <cctbx/sgtbx/direct_space_asu/proto/direct_space_asu.h>
#include <boost/shared_ptr.hpp> // necessary for pair_tables.h
#include <cctbx/crystal/point_neighbors.h>
#include <cctbx/crystal/neighbors_simple.h>

// @@@@@@@@
// This file is a translation of find_distances_using_cpp_objects.py into C++
// @@@@@@@@

namespace cc {
  using namespace cctbx::sgtbx;
  using namespace cctbx::uctbx;
  using namespace cctbx::sgtbx::asu;
  using namespace cctbx::crystal::direct_space_asu;
  using namespace cctbx::crystal::neighbors;
  typedef scitbx::vec3<double> double3;
};

void find_distances(
  const cc::unit_cell &unit_cell,
  const cc::space_group &space_group,
  const scitbx::af::const_ref< cc::double3 > &sites_frac_par,
  double distance_cutoff)
{
  scitbx::af::shared< cc::double3 > n_sites_frac=unit_cell.orthogonalize(sites_frac_par);
  scitbx::af::shared< cc::double3 > old_sites_frac=unit_cell.fractionalize(n_sites_frac.const_ref());
  scitbx::af::const_ref< cc::double3 > sites_frac=old_sites_frac.const_ref();
  cc::direct_space_asu metric_free_asu( space_group.type() );
  cc::float_asu<> float_asu = metric_free_asu.as_float_asu(unit_cell, 1.0E-6);
  boost::shared_ptr< cc::asu_mappings<> > asu_mappings(
    new cc::asu_mappings<> (space_group, float_asu, distance_cutoff) );
  asu_mappings->process_sites_frac(sites_frac, 0.5);

  cc::point_neighbors<> ngbr_struct(
          asu_mappings,
          distance_cutoff*(1+1.0E-6));

  std::vector< cc::neighborSymmetryIndexAndDistanceSq<> > ngbr_list;
  std::vector< cc::neighborSymmetryIndexAndDistanceSq<> >::iterator itr;

  typedef cc::simple_pair_generator<> base_t;
  typedef base_t::asu_mappings_t asu_mappings_t;

  scitbx::af::const_ref<asu_mappings_t::array_of_mappings_for_one_site>
  const& mappings = asu_mappings->mappings_const_ref();

  for(unsigned i=0; i<2; ++i) // af::shared  array
  {
  //  cc::double3 p(sites_frac[i]);
   cctbx::crystal::direct_space_asu::asu_mapping_index mi;
   
   mi.i_seq=i;
   mi.i_sym=0;
   ngbr_list=ngbr_struct.neighborsOfPoint(mappings[i][0].mapped_site());
   cctbx::fractional<> i_pos=sites_frac[i];
//   std::printf("%d pos %f %f %f\n",i,i_pos[0],i_pos[1],i_pos[2]);
    for(itr=ngbr_list.begin(); itr!=ngbr_list.end(); itr++){
	cc::rt_mx rot_i=asu_mappings->get_rt_mx(i,0);
	cc::rt_mx rot_j=asu_mappings->get_rt_mx(itr->ngbr_seq,itr->ngbr_sym);
	cc::rt_mx rot_ji=rot_i.inverse().multiply(rot_j);
	const cctbx::fractional<> frac_j = rot_ji*sites_frac[itr->ngbr_seq];	
	cctbx::cartesian<> ngbr_pos=mappings[itr->ngbr_seq][itr->ngbr_sym].mapped_site();
	double dist=unit_cell.distance(i_pos,frac_j);
	std::printf("%d%d%1.3f ngbr:%f %f %f sym: %d %-20s %f\n",i,itr->ngbr_seq,sqrt(itr->dist_sq),ngbr_pos[0],ngbr_pos[1],ngbr_pos[2],
		itr->ngbr_sym,rot_ji.as_xyz().c_str(),dist);
    }
  }

  ngbr_struct.reposition(0,sites_frac[1],0.5);
  ngbr_struct.reposition(1,sites_frac[0],0.5);
  for(unsigned i=0; i<2; ++i) // af::shared  array
  {
  //  cc::double3 p(sites_frac[i]);
   cctbx::crystal::direct_space_asu::asu_mapping_index mi;
   
   mi.i_seq=i;
   mi.i_sym=0;
   ngbr_list=ngbr_struct.neighborsOfPoint(mappings[i][0].mapped_site());
   cctbx::fractional<> i_pos=sites_frac[1-i];
//   std::printf("%d pos %f %f %f\n",i,i_pos[0],i_pos[1],i_pos[2]);
    for(itr=ngbr_list.begin(); itr!=ngbr_list.end(); itr++){
	cc::rt_mx rot_i=asu_mappings->get_rt_mx(i,0);
	cc::rt_mx rot_j=asu_mappings->get_rt_mx(itr->ngbr_seq,itr->ngbr_sym);
	cc::rt_mx rot_ji=rot_i.inverse().multiply(rot_j);
	const cctbx::fractional<> frac_j = rot_ji*sites_frac[1-itr->ngbr_seq];	
	cctbx::cartesian<> ngbr_pos=mappings[itr->ngbr_seq][itr->ngbr_sym].mapped_site();
	double dist=unit_cell.distance(i_pos,frac_j);
	std::fprintf(stderr,"%d%d%1.3f ngbr:%f %f %f sym: %d %-20s %f\n",1-i,1-itr->ngbr_seq,sqrt(itr->dist_sq),ngbr_pos[0],ngbr_pos[1],ngbr_pos[2],
		itr->ngbr_sym,rot_ji.as_xyz().c_str(),dist);
    }
  }

}

int main(int argc, char *[])
{
  CCTBX_ASSERT( argc==1 );
  scitbx::af::double6 cell(5.01, 5.01, 5.47, 90.0, 90.0, 120.0);
  cc::unit_cell unit_cell(cell);
  cc::space_group space_group(cctbx::sgtbx::space_group_symbols("P6222"));
  scitbx::af::shared< cc::double3 > sites_frac;
  sites_frac.push_back( cc::double3(0.5, 0.5, 1.0/3.0) );
  sites_frac.push_back( cc::double3(0.197, -0.197, 0.83333) );
  find_distances(unit_cell, space_group, sites_frac.const_ref(), 5.0);
  return 0;
}

