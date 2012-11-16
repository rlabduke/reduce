
#include <cstdio>
#include <fstream>
#include <cctbx/sgtbx/direct_space_asu/proto/direct_space_asu.h>
#include <boost/shared_ptr.hpp> // necessary for pair_tables.h
#include <cctbx/crystal/pair_tables.h>

// @@@@@@@@
// This file is a translation of find_distances_using_cpp_objects.py into C++
// @@@@@@@@

namespace cc {
  using namespace cctbx::sgtbx;
  using namespace cctbx::uctbx;
  using namespace cctbx::sgtbx::asu;
  using namespace cctbx::crystal::direct_space_asu;
  using cctbx::crystal::pair_asu_table;
  using cctbx::crystal::pair_sym_table;
  using cctbx::crystal::pair_sym_dict;
  using cctbx::crystal::pair_sym_ops;
  typedef scitbx::vec3<double> double3;
};

void find_distances(
  const cc::unit_cell &unit_cell,
  const cc::space_group &space_group,
  const scitbx::af::const_ref< cc::double3 > &sites_frac,
  double distance_cutoff)
{
  // Get the stuff needed for construcing mappings
  cc::direct_space_asu metric_free_asu( space_group.type() );
  cc::float_asu<> float_asu = metric_free_asu.as_float_asu(unit_cell, 1.0E-6);

  // Create the asu mappings... this is a container that maps sites to asymmetric units
  boost::shared_ptr< cc::asu_mappings<> > asu_mappings(
    new cc::asu_mappings<> (space_group, float_asu, distance_cutoff) );
  asu_mappings->process_sites_frac(sites_frac);//, 0.5);    


  // Managed table of pair interactions based on direct_space_asu::asu_mappings
  // ... does this mean anything to you?
  //Uses neighbors::fast_pair_generator to add all pairs with 
  //    distances <= distance_cutoff*(1+epsilon).
  //All symmetrically equivalent pairs are automatically generated.
  cc::pair_asu_table<> pair_asu_table(asu_mappings);   
  pair_asu_table.add_all_pairs(distance_cutoff);



  // Table of pair interactions indexed by i_seq.
  cc::pair_sym_table pair_sym_table = pair_asu_table.extract_pair_sym_table();



  // stopped here trying to figure out wtf is going on with this code
  for(unsigned i=0; i<pair_sym_table.size(); ++i) // af::shared  array
  {
    std::printf("i: %d\n", i);
    const cctbx::fractional<> frac_i = sites_frac[i];
    const cc::pair_sym_dict pair_sym_dict = pair_sym_table[i];
    for(cc::pair_sym_dict::const_iterator pair=pair_sym_dict.begin();
      pair!=pair_sym_dict.end(); ++pair ) // std::map
    {
      const int j = pair->first;
      std::printf("  j: %d\n", j);
      const cctbx::fractional<> frac_j = sites_frac[j];
      const cc::pair_sym_ops sym_ops = pair->second;
      for(cc::pair_sym_ops::const_iterator sym_op = sym_ops.begin();
        sym_op!=sym_ops.end(); ++sym_op) // std::vector
      {
        const cctbx::fractional<> frac_ji = (*sym_op) * frac_j;
          std::printf("    %-20s %8.3f\n", sym_op->as_xyz().c_str(),
            unit_cell.distance(frac_i, frac_ji));
      }
    }
  }
   std::cout<<"Pair sym table has size "<<pair_sym_table.size()<<std::endl;
}

int main(int argc, char *[])
{
  CCTBX_ASSERT( argc==1 );
 // scitbx::af::double6 cell(5.01, 5.01, 5.47, 90.0, 90.0, 120.0);
  scitbx::af::double6 cell(360.0, 360.0, 360.0, 90.0, 90.0, 90.0);
  cc::unit_cell unit_cell(cell);
  cc::space_group space_group(cctbx::sgtbx::space_group_symbols("F 41 3 2"));
  //cc::space_group space_group(cctbx::sgtbx::space_group_symbols("P6222"));
  scitbx::af::shared< cc::double3 > sites_frac;

  std::ifstream fptr("1ZFN");

  double x,y,z;
  int n=0;
  while (fptr){	
	fptr>>x;
	fptr>>y;
	fptr>>z;
//	std::cout<<x<<" "<<y<<" "<<z<<std::endl;
  	sites_frac.push_back( cc::double3(x,y,z) );
	n++;
  }
//  sites_frac.push_back( cc::double3(0.5, 0.5, 1.0/3.0) );
 // sites_frac.push_back( cc::double3(0.197, -0.197, 0.83333) );
  find_distances(unit_cell, space_group, sites_frac.const_ref(), 5.0);
  std::cout<<"Number of atoms: "<<n<<std::endl;
  return 0;
}


