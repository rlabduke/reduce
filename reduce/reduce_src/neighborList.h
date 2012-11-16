#include <cstdio>
#include <list>
#include <cctbx/sgtbx/direct_space_asu/proto/direct_space_asu.h>
#include <boost/shared_ptr.hpp> // necessary for pair_tables.h
#include "point_neighbors.h"
#include <cctbx/crystal/neighbors_simple.h>
#include <utility>

#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H 1

typedef scitbx::vec3<double> double3;

namespace cc {
  using namespace cctbx::sgtbx;
  using namespace cctbx::uctbx;
  using namespace cctbx::sgtbx::asu;
  using namespace cctbx::crystal::direct_space_asu;
  using namespace cctbx::crystal::neighbors;
  typedef scitbx::vec3<double> double3;
};

static double rt_mx_dist2(scitbx::af::tiny<float, 12> m1,scitbx::af::tiny<float, 12> m2){
	double dist2=0;
	for (int i=0; i<12; i++){
		double diff=m1[i]-m2[i];
		dist2 += (diff*diff);
	}
	return dist2;
}

static int get_chain_index(int index){
	if (index < 26)
		return 65+index;
	else if (index < 36)
		return 22+index;
	else
		return 61+index;
}

template <class L>
class NeighborList{

	private:
	cc::point_neighbors<> ngbr_struct;
	Coord _max_range;
	cc::unit_cell _unit_cell;
	bool _cubiclized;
//	bool _nosym;
//	Coord _offset[3];
//	Point3d _offset_p3d;

	public:

	std::vector<L> _atom_list;

	NeighborList(){
		_max_range=-1.0f;
		_cubiclized=false;
//		_nosym=false;
	}

	NeighborList(Coord max_range){
		_max_range=max_range;
		_cubiclized=false;
//		_nosym=false;
	}

	void clear_list(){
		for(typename std::vector<L >::const_iterator itr=_atom_list.begin(); itr != _atom_list.end(); itr++)
			delete	*itr;
	}

  	void get_unit_cell_params(double uparams[]){
		uparams[0]=_unit_cell.parameters()[0];
		uparams[1]=_unit_cell.parameters()[1];
		uparams[2]=_unit_cell.parameters()[2];
	}

	NeighborList(scitbx::af::double6 cell, char* space_grp, Coord distance_cutoff, bool NoSym=false, Coord bbox[]=NULL){//, Coord box_center[]=NULL){
		fprintf(stderr,"Received space group type %s\n",space_grp);
		_max_range=distance_cutoff;
//		_nosym=false;
		cctbx::sgtbx::space_group spaceGrp;
		if (!NoSym){
			try {
				spaceGrp=cctbx::sgtbx::space_group(cctbx::sgtbx::space_group_symbols(space_grp));
			}
			catch (cctbx::error e){
				NoSym=true;
				fprintf(stderr,"Unrecognized space group\n");
			}
		}
		if (NoSym) {
//			_nosym=true;
			strcpy(space_grp,"P 1");
			spaceGrp=cctbx::sgtbx::space_group(cctbx::sgtbx::space_group_symbols(space_grp));
			cell[0]=bbox[0]*10;
			cell[1]=bbox[1]*10;
			cell[2]=bbox[2]*10;
//			_offset[0]=cell[0]/2-box_center[0];
//			_offset[1]=cell[1]/2-box_center[1];
//			_offset[2]=cell[2]/2-box_center[2];
//			_offset_p3d=Point3d(_offset[0],_offset[1],_offset[2]);
			
		}
		cc::direct_space_asu metric_free_asu( spaceGrp.type() );
		_unit_cell=cc::unit_cell(cell);
		cc::float_asu<> float_asu = metric_free_asu.as_float_asu(_unit_cell, 1.0E-6);
		boost::shared_ptr< cc::asu_mappings<> > asu_mappings(
		  new cc::asu_mappings<> (spaceGrp, float_asu, _max_range) );
		ngbr_struct=cc::point_neighbors<>(asu_mappings);
		_cubiclized=false;
	}


	void init(std::vector<L> atoms){
		_atom_list=atoms;
		// Add sites to _sites_frac
		scitbx::af::shared< double3 > sites;
		for(typename std::vector<L >::const_iterator itr=atoms.begin(); itr != atoms.end(); itr++){
			Point3d site=(*itr)->loc();
//			if (_nosym)		// ************NOSYMCHANGES*************
//				site=site+_offset_p3d;
			sites.push_back( double3(site.x(), site.y(), site.z()) );
		}
		
		scitbx::af::shared< double3 > sites_frac=_unit_cell.fractionalize(sites.const_ref());
		ngbr_struct.insert_many_sites(sites_frac.const_ref(), 0.5);

	}

	void insert(L atom){
		atom->set_index(_atom_list.size());
		_atom_list.push_back(atom);
		Point3d site=atom->loc();
//		if (_nosym)		// ************NOSYMCHANGES*************
//			site=site+_offset_p3d;
		cctbx::cartesian<> pt(site.x(),site.y(),site.z());
		cctbx::fractional<> pt_frac=_unit_cell.fractionalize(pt);
		ngbr_struct.insert(pt_frac,_cubiclized);
	}

	void cubiclize(){
		ngbr_struct.init_cubicles(_max_range*(1+1.0E-6));
		_cubiclized=true;
	}

	std::list< std::pair<L,Point3d> > get_neighbors(Point3d p, Coord min_range, Coord max_range=-1.0f){
		std::list< std::pair<L, Point3d> > neighbors;
//		if (_nosym)		// ************NOSYMCHANGES*************
//			p=p+_offset_p3d;
		cctbx::cartesian<> query(p.x(),p.y(),p.z());

		if (!_cubiclized){
			std::cerr<<"***********DATA not cubiclized**********"<<std::endl;
			return neighbors;
		}

		if (max_range < 0.0)
			max_range = _max_range;

		cctbx::sgtbx::rt_mx transformation_matrix;
		std::vector< cc::neighborSymmetryIndexAndDistanceSq<> > ngbrs_cctbx= ngbr_struct.neighborsOfPoint(query, &transformation_matrix);
		std::vector< cc::neighborSymmetryIndexAndDistanceSq<> >::const_iterator itr;
		for(itr=ngbrs_cctbx.begin(); itr!=ngbrs_cctbx.end(); ++itr ){
			cctbx::cartesian<> ngbr_pos=ngbr_struct.get_coordinates(*itr);

			cctbx::fractional<> ngbr_trans=_unit_cell.orthogonalize(transformation_matrix.inverse()*(_unit_cell.fractionalize(ngbr_pos)));
			Coord dist=sqrt(itr->dist_sq);
//			std::cerr<<"ngbrdists "<<itr->dist_sq<<" "<<(ngbr_trans-query).length_sq()<<std::endl;
			if((dist>min_range)&& (dist < max_range)){
//				if (_nosym)		// ************NOSYMCHANGES*************
//				neighbors.push_back(std::make_pair(_atom_list[(int)itr->ngbr_seq], Point3d(ngbr_trans[0],ngbr_trans[1],ngbr_trans[2])-_offset_p3d));
//				else
				neighbors.push_back(std::make_pair(_atom_list[(int)itr->ngbr_seq], Point3d(ngbr_trans[0],ngbr_trans[1],ngbr_trans[2])));
			}
		}
		return neighbors;
	}

	void reposition(int seq, Point3d p, L atom){
//		if (_nosym)		// ************NOSYMCHANGES*************
//			p=p+_offset_p3d;
		cctbx::cartesian<> query(p.x(),p.y(),p.z());
		cctbx::fractional<> q=_unit_cell.fractionalize(query);
//		delete _atom_list[seq];	// Uncomment when _xyzblocks is taken out.
		_atom_list[seq]=atom;
		ngbr_struct.reposition((unsigned)seq, q);
		_atom_list[seq]->set_index(seq); // shouldn't have to do this explicitly
	}

	int find(L atom){
		for (int i=0; i<_atom_list.size(); i++)
			if (*atom == *(_atom_list[i]))
				return i;
		return -1;
	}

	int find(L atom, Point3d p){
		int retval=-1;
		int *cell_contents=NULL;
//		if (_nosym)		// ************NOSYMCHANGES*************
//			p=p+_offset_p3d;
		cctbx::cartesian<> q(p.x(),p.y(),p.z());
		int num=ngbr_struct.get_cell_contents(q,cell_contents);
		for (int i=0; i<num; i++){
			if (*atom == *(_atom_list[cell_contents[i]])){
				retval=cell_contents[i];
				break;
			}
		}
		if (retval < 0){
			std::cerr<<"Candidates for: "<<(*atom)<<std::endl;
			for (int i=0; i< num; i++)
				std::cerr<<*(_atom_list[cell_contents[i]])<<std::endl;
			std::cerr<<std::endl;
		}
		delete cell_contents;
		return retval;
	}

	int total_size(){
		return ngbr_struct.total_size();
	}

	void print(int seq, std::ostream &os){
		ngbr_struct.print(seq,os);
	}
	
	// Debugging code: Might not work for BumperPoints and changes the loc  information in PDB
	void print_all(std::ostream& out, std::ostream& pdbout){
		int sz=_atom_list.size();
		int *sym_size=new int[sz];
		std::list< scitbx::af::tiny<float,12> > rt_mxs;

		// Get maximum number of symmetric copies for each atom
		for(int i=0; i<sz; i++)
			sym_size[i]=ngbr_struct.get_sym_size(i);
		cc::neighborSymmetryIndexAndDistanceSq<> pair;

		cctbx::cartesian<> sym_site;
		int sym_index;
		for (pair.ngbr_seq=0; pair.ngbr_seq<sz; pair.ngbr_seq++){
			Point3d site=_atom_list[pair.ngbr_seq]->loc();
			for (pair.ngbr_sym=0; pair.ngbr_sym<sym_size[pair.ngbr_seq]; pair.ngbr_sym++){
				sym_site=ngbr_struct.get_coordinates(pair);

				scitbx::af::tiny<float, 12> rt_mx=ngbr_struct.get_rt_mx(pair);
				sym_index=0;
				std::list< scitbx::af::tiny<float,12> >::iterator itr=rt_mxs.begin();
				bool not_found=true;
				while (itr != rt_mxs.end()){
					if (rt_mx_dist2(*itr, rt_mx) < 1e-6){
						not_found=false;
						break;
					}
					sym_index++;
					itr++;
				}

				if (not_found)
					rt_mxs.push_back(rt_mx);

				Point3d sym_coords(sym_site[0],sym_site[1],sym_site[2]);

				//**** Changing coordinates before printing
				out<<sym_coords.x()<<" "<<sym_coords.y()<<" "<<sym_coords.z()<<" "<<_atom_list[pair.ngbr_seq]->vdwRad()<<" "<<sym_index<<" ";
				out<<_atom_list[pair.ngbr_seq]->resname()<<_atom_list[pair.ngbr_seq]->resno()<<" ";
				out<<_atom_list[pair.ngbr_seq]->atomname()<<std::endl;
				_atom_list[pair.ngbr_seq]->loc(sym_coords);
				_atom_list[pair.ngbr_seq]->set_chain_id(get_chain_index(sym_index));
				pdbout<<*_atom_list[pair.ngbr_seq]<<std::endl;
			}
		}
		rt_mxs.clear();
		delete sym_size;
	}

	void print_all_sym_ngbrs(float cutoff, std::ostream& pdbout){
		int sz=_atom_list.size();
		Point3d* atom_centers;
		atom_centers=new Point3d[sz];

		// Get maximum number of symmetric copies for each atom
		for(int i=0; i<sz; i++)
			atom_centers[i]=_atom_list[i]->loc();

		cctbx::cartesian<> sym_site,ngbr_trans;
//		cctbx::fractional<> ngbr_trans;
		cctbx::sgtbx::rt_mx transformation_matrix;
		int sym_index;
		for (int i=0; i<sz; i++){
			cctbx::cartesian<> curr_site(atom_centers[i].x(),atom_centers[i].y(),atom_centers[i].z());
			std::vector< cc::neighborSymmetryIndexAndDistanceSq<> > ngbrs_cctbx= ngbr_struct.neighborsOfPoint(curr_site, &transformation_matrix, true);
			std::vector< cc::neighborSymmetryIndexAndDistanceSq<> >::const_iterator itr;
			for(itr=ngbrs_cctbx.begin(); itr!=ngbrs_cctbx.end(); ++itr ){
				if (!itr->diff_copy)
					continue;

				if (sqrt(itr->dist_sq) > _atom_list[i]->vdwRad()+_atom_list[itr->ngbr_seq]->vdwRad()+cutoff)
					continue;
				
				sym_site=ngbr_struct.get_coordinates(*itr);
				ngbr_trans=_unit_cell.orthogonalize(transformation_matrix.inverse()*(_unit_cell.fractionalize(sym_site)));
				Point3d sym_coords(ngbr_trans[0],ngbr_trans[1],ngbr_trans[2]);

				//**** Changing coordinates before printing
				_atom_list[itr->ngbr_seq]->loc(sym_coords);
				pdbout<<*_atom_list[itr->ngbr_seq]<<std::endl;
			}
		}
		delete[] atom_centers;
	}
};
#endif
