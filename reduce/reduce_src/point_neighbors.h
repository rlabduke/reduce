#ifndef CCTBX_CRYSTAL_POINT_NEIGHBORS_H
#define CCTBX_CRYSTAL_POINT_NEIGHBORS_H

#include <cctbx/crystal/direct_space_asu.h>
#include <scitbx/cubicles.h>

namespace cctbx { namespace crystal { namespace neighbors {

template <typename FloatType=double>
struct neighborSymmetryIndexAndDistanceSq{
	unsigned ngbr_seq;
	unsigned ngbr_sym;
	bool diff_copy;

	FloatType dist_sq;
};


  template <typename FloatType=double, typename IntShiftType=int>
  class point_neighbors{
      typedef direct_space_asu::asu_mappings<FloatType, IntShiftType> asu_mappings_t;

      FloatType epsilon_;
      typedef std::vector<direct_space_asu::asu_mapping_index> box_content_t;
      scitbx::cubicles<box_content_t, FloatType> cubicles_;
      scitbx::vec3<unsigned> n_boxes_;
      boost::shared_ptr<asu_mappings_t> asu_mappings_owner_;
      asu_mappings_t* asu_mappings_;
      unsigned mappings_size_;
      FloatType distance_cutoff_sq_;
	
      public:

	point_neighbors(){};

	point_neighbors(boost::shared_ptr< direct_space_asu::asu_mappings<FloatType, IntShiftType> > const& asu_mappings,FloatType const& epsilon=1.e-6)
	:
	epsilon_(epsilon)
	{
	CCTBX_ASSERT(epsilon > 0);
	CCTBX_ASSERT(epsilon < 0.01);
	this->asu_mappings_owner_ = asu_mappings;
	this->asu_mappings_ = asu_mappings.get();
	}

	point_neighbors(cctbx::sgtbx::space_group spaceGrp, cctbx::uctbx::unit_cell cell, FloatType max_range, FloatType const& epsilon=1.e-6)
        :
        epsilon_(epsilon)
	{
		CCTBX_ASSERT(epsilon > 0);
		CCTBX_ASSERT(epsilon < 0.01);
		cctbx::sgtbx::asu::direct_space_asu metric_free_asu( spaceGrp.type() );
                cctbx::crystal::direct_space_asu::float_asu<> float_asu = metric_free_asu.as_float_asu(cell, 1.0E-6);
                asu_mappings_owner_=boost::shared_ptr<asu_mappings_t>(new cctbx::crystal::direct_space_asu::asu_mappings<>(spaceGrp, float_asu, max_range));
		asu_mappings_=asu_mappings_owner_.get();
	}

	void init_cubicles(FloatType distance_cutoff, FloatType const& min_cubicle_edge=5){
		this->distance_cutoff_sq_ = distance_cutoff*distance_cutoff;
		direct_space_asu::float_asu<FloatType> asu_buffer(asu_mappings_->asu().add_buffer(distance_cutoff));
		cubicles_=scitbx::cubicles<box_content_t, FloatType>(asu_buffer.box_min(true),
          							     asu_buffer.box_max(true)-asu_buffer.box_min(true),
								     std::max(distance_cutoff, min_cubicle_edge),
          							     epsilon_);
		/*cubicles_=scitbx::cubicles<box_content_t, FloatType>(asu_mappings_->mapped_sites_min(),
          							     asu_mappings_->mapped_sites_span(),
								     std::max(distance_cutoff, min_cubicle_edge),
          							     epsilon_);*/
		n_boxes_=cubicles_.ref.accessor();
		af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
		  const& mappings = this->asu_mappings_->mappings_const_ref();
		direct_space_asu::asu_mapping_index mi;

		for(mi.i_seq=0;mi.i_seq<mappings.size();mi.i_seq++) {
		  for(mi.i_sym=0; mi.i_sym<mappings[mi.i_seq].size(); mi.i_sym++) {
		    std::size_t i1d_cub = cubicles_.ref.accessor()(
		      cubicles_.i_cubicle(mappings[mi.i_seq][mi.i_sym].mapped_site()));
	  	      cubicles_.ref[i1d_cub].push_back(mi);
	  	  }
		}
	}

	void insert(fractional<FloatType> const& original_site, bool cubiclized){
		asu_mappings_->process(original_site);
		af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
		  const& mappings = this->asu_mappings_->mappings_const_ref();
		insert_into_cubicles(mappings.size()-1,cubiclized);
	}

	void insert_many_sites(af::const_ref<scitbx::vec3<FloatType> > const& original_sites,  FloatType const& min_distance_sym_equiv=0.5){
		asu_mappings_->process_sites_frac(original_sites,min_distance_sym_equiv);
		af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
		  const& mappings = this->asu_mappings_->mappings_const_ref();
		direct_space_asu::asu_mapping_index mi;

		for(mi.i_seq=0;mi.i_seq<mappings.size();mi.i_seq++) {
		  for(mi.i_sym=0; mi.i_sym<mappings[mi.i_seq].size(); mi.i_sym++) {
		    std::size_t i1d_cub = cubicles_.ref.accessor()(
		      cubicles_.i_cubicle(mappings[mi.i_seq][mi.i_sym].mapped_site()));
		    cubicles_.ref[i1d_cub].push_back(mi);
		  }
		}
	}

	void reposition(unsigned seq,fractional<FloatType> const& original_site,FloatType const& min_distance_sym_equiv=0.5){
		delete_from_cubicles(seq);
		asu_mappings_owner_->replace(seq, original_site,min_distance_sym_equiv);
		insert_into_cubicles(seq,true);
	}

	void insert_into_cubicles(int seq, bool cubiclized){
	
	  af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
	  const& mappings = this->asu_mappings_->mappings_const_ref();
	  direct_space_asu::asu_mapping_index mi;
	  mi.i_seq=seq;

	  if (cubiclized){		
	    for(mi.i_sym=0; mi.i_sym<mappings[mi.i_seq].size(); mi.i_sym++) {
	      std::size_t i1d_cub = cubicles_.ref.accessor()(
	        cubicles_.i_cubicle(mappings[mi.i_seq][mi.i_sym].mapped_site()));
	      cubicles_.ref[i1d_cub].push_back(mi);
	    }
	  }
	}

	void delete_from_cubicles(int seq){
	  
		af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
		const& mappings = this->asu_mappings_->mappings_const_ref();

		for(unsigned i_sym=0; i_sym<mappings[seq].size(); i_sym++) {
			std::size_t i1d_cub = cubicles_.ref.accessor()(cubicles_.i_cubicle(mappings[seq][i_sym].mapped_site()));
			for (box_content_t::iterator itr=cubicles_.ref[i1d_cub].begin();itr!=cubicles_.ref[i1d_cub].end();){
				if (itr->i_seq == seq){
					itr=cubicles_.ref[i1d_cub].erase(itr);
				}
				else
					itr++;
			}
	  	}
	}

	void print(int seq, std::ostream &os){
		af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
		const& mappings = this->asu_mappings_->mappings_const_ref();

		for(unsigned i_sym=0; i_sym<mappings[seq].size(); i_sym++) {
			cartesian<FloatType> pt=mappings[seq][i_sym].mapped_site();
			os<<pt[0]<<" "<<pt[1]<<" "<<pt[2]<<std::endl;
		}
        }

	cctbx::cartesian<FloatType> get_coordinates(neighborSymmetryIndexAndDistanceSq<FloatType> pair){
		af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
		const& mappings = this->asu_mappings_->mappings_const_ref();
		return mappings[pair.ngbr_seq][pair.ngbr_sym].mapped_site();
	}
	
	unsigned get_symmetry_index(neighborSymmetryIndexAndDistanceSq<FloatType> pair){
		af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
		const& mappings = this->asu_mappings_->mappings_const_ref();
		return mappings[pair.ngbr_seq][pair.ngbr_sym].i_sym_op();
	}
	
	af::tiny<float, 12> get_rt_mx(neighborSymmetryIndexAndDistanceSq<FloatType> pair){
//		af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
//		const& mappings = this->asu_mappings_->mappings_const_ref();
		return asu_mappings_->get_rt_mx(pair.ngbr_seq,pair.ngbr_sym).as_float_array();
	}
	
	int get_sym_size(int seq){
		af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
		const& mappings = this->asu_mappings_->mappings_const_ref();
		return mappings[seq].size();
	}
	
	cartesian<FloatType> map_to_asu(cartesian<FloatType> query, cctbx::sgtbx::rt_mx* transformation_matrix=NULL){
		direct_space_asu::float_asu<FloatType> asu=asu_mappings_->asu();
		fractional<FloatType> qfrac=asu.unit_cell().fractionalize(query);
		
		asu_mappings_->process(qfrac);//,mindistequiv

		af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
		const& mappings = this->asu_mappings_->mappings_const_ref();

		std::size_t sz=mappings.size();
		cartesian<FloatType> p=mappings[sz-1][0].mapped_site();

		if (transformation_matrix != NULL)
			*transformation_matrix=asu_mappings_->get_rt_mx(sz-1,0);

		asu_mappings_->discard_last();

		return p;
	}

	int total_size(){
		af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
		const& mappings = this->asu_mappings_->mappings_const_ref();
		scitbx::vec3<FloatType> space_min=asu_mappings_->mapped_sites_min();
          	scitbx::vec3<FloatType> space_span=asu_mappings_->mapped_sites_span();
		std::cerr<<"Min: "<<space_min[0]<<" "<<space_min[1]<<" "<<space_min[2]<<std::endl; 
		std::cerr<<"Span: "<<space_span[0]<<" "<<space_span[1]<<" "<<space_span[2]<<std::endl; 
		int retval =0;
		for(int i=0; i<mappings.size(); i++){
			for (int j=0; j<mappings[i].size(); j++){
				cctbx::cartesian<FloatType> current=mappings[i][j].mapped_site();
			}
			retval+=mappings[i].size();
		}
		scitbx::vec3<unsigned> ind;
		const box_content_t* box_ngbr;
		std::cerr<<"Cubicle sizes: "<<std::endl;
		for(ind[0]=0;ind[0]<n_boxes_[0];ind[0]++)
			for(ind[1]=0;ind[1]<n_boxes_[1];ind[1]++)
				for(ind[2]=0;ind[2]<n_boxes_[2];ind[2]++){
					box_ngbr=&cubicles_.ref(ind);
					std::cerr<<box_ngbr->size()<<" ";
				}

		std::cerr<<std::endl;

		std::cerr<<n_boxes_[0]<<" "<<n_boxes_[1]<<" "<<n_boxes_[2]<<std::endl;
		return retval;
	}

	int get_cell_contents(cartesian<FloatType> query, int* &contents){
		cartesian<FloatType> p=map_to_asu(query);
		scitbx::vec3<unsigned> box_indices=cubicles_.i_cubicle(p);

		const box_content_t* box_ngbr=&cubicles_.ref(box_indices);
		typename box_content_t::const_iterator box_itr;
		int retval=box_ngbr->size();
		contents=new int[retval];
		int i=0;
		for (box_itr=box_ngbr->begin(),i=0; box_itr != box_ngbr->end(); box_itr++, i++)
			contents[i]=box_itr->i_seq;
		return retval;
	}

	std::vector< neighborSymmetryIndexAndDistanceSq<FloatType> > neighborsOfPoint(cartesian<FloatType> query, cctbx::sgtbx::rt_mx* transformation_matrix, bool check_diff_copy=false)
	{
		cartesian<FloatType> p=map_to_asu(query,transformation_matrix);
		scitbx::vec3<unsigned> box_indices=cubicles_.i_cubicle(p);
		scitbx::vec3<unsigned> ind_max,ind_min,ind;
		af::const_ref<typename asu_mappings_t::array_of_mappings_for_one_site>
		const& mappings = this->asu_mappings_->mappings_const_ref();

		for (int i=0; i<3; i++){
			ind_min[i]= ( box_indices[i] == 0 ? 0 : box_indices[i]-1);
			ind_max[i]= ( box_indices[i] == n_boxes_[i]-1 ? n_boxes_[i]-1 : box_indices[i]+1);
		}

		int ngbr_size=0;
      		const box_content_t* box_ngbr;
		typename box_content_t::const_iterator box_itr;

		neighborSymmetryIndexAndDistanceSq<FloatType> elem;
		std::vector< neighborSymmetryIndexAndDistanceSq<FloatType> > ngbr_list;

		for (ind[0]=ind_min[0]; ind[0]<=ind_max[0]; ind[0]++){
		for (ind[1]=ind_min[1]; ind[1]<=ind_max[1]; ind[1]++){
		for (ind[2]=ind_min[2]; ind[2]<=ind_max[2]; ind[2]++){
			box_ngbr=&cubicles_.ref(ind);
			for (box_itr=box_ngbr->begin(); box_itr != box_ngbr->end(); box_itr++){
				cartesian<FloatType> ngbr_pos=mappings[box_itr->i_seq][box_itr->i_sym].mapped_site();
				elem.dist_sq=(ngbr_pos-p).length_sq();

				if (elem.dist_sq > this->distance_cutoff_sq_) continue;
				
				elem.ngbr_seq=box_itr->i_seq;
				elem.ngbr_sym=box_itr->i_sym;
				
				if (check_diff_copy){
					af::tiny<float, 12> ngbr_mat=asu_mappings_->get_rt_mx(elem.ngbr_seq,elem.ngbr_sym).as_float_array();
					af::tiny<float, 12> query_mat=transformation_matrix->as_float_array();

					double mat_diff_sum=0;
					for (int i=0; i<12; i++){
						double mat_i_diff=ngbr_mat[i]-query_mat[i];
						mat_diff_sum+=mat_i_diff*mat_i_diff;
					}
					elem.diff_copy=(mat_diff_sum > 1e-6);
				}
				
				ngbr_list.push_back(elem);
			}
		}}};

		return ngbr_list;
	}

	std::vector< neighborSymmetryIndexAndDistanceSq<FloatType> > neighborsOfPoint(FloatType x, FloatType y, FloatType z){
		cartesian<FloatType> p(x,y,z);
		return neighborsOfPoint(p);
	}

   };

}}} // namespace cctbx::crystal::neighbors

#endif // CCTBX_CRYSTAL_POINT_NEIGHBORS_H
