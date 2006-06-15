#include "GraphToHoldScores.h"
#include "AtomPositions.h"
#include "Mover.h"
#include "PDBrec.h"

void sort_three(int a, int b, int c, int & x, int & y, int & z)
{
	if ( a < b )
	{
		if ( b < c )
		{
			x = a;
			y = b;
			z = c;
		}
		else 
		{
			z = b;
			if ( a < c )
			{
				x = a;
				y = c;
			}
			else
			{
				x = c;
				y = a;
			}
		}
	}
	else
	{
		if ( c < b )
		{				
			x = c;
			y = b;
			z = a;
		}
		else
		{
			x = b;
			if ( a < c )
			{
				y = a;
				z = c;
			}
			else
			{
				y = c;
				z = a;
			}
		}
	}
}


template < class T >
void set_intersect_apl( 
	std::list< T > const & list1, 
	std::list< T > const & list2,
	std::list< T > & intersection
)
{
	//std::cerr << "clearing" << std::endl;
	intersection.clear();
	//std::cerr << "done clearing" << std::endl;
	typename std::list< T >::const_iterator iter1 = list1.begin();
	typename std::list< T >::const_iterator iter2 = list2.begin();
	while ( iter1 != list1.end() && iter2 != list2.end() )
	{
		//std::cerr << "iter1: " << *iter1 << " iter2: " << *iter2 << std::endl;
		if ( *iter1 == *iter2 )
		{
			//std::cerr << "set_intersect_apl: " << *iter1 << std::endl;
			intersection.push_back( *iter1 );
			++iter1;
			++iter2;
		}
		else if ( *iter1 < *iter2 )
		{
			++iter1;
		}
		else
		{
			++iter2;
		}
	}
	//std::cerr << "leaving set_intersect_apl" << std::endl;
}


template <class T >
void set_difference_apl(
	std::list< T > & sortedList1,
	std::list< T > const & sortedList2 
)
{
	typename std::list< T >::iterator list1iter = sortedList1.begin();
	typename std::list< T >::const_iterator list2iter = sortedList2.begin();
	while( list1iter != sortedList1.end() && list2iter != sortedList2.end()  )
	{
		//std::cerr << "set diff: iter1: " << *list1iter << " iter2: " << *list2iter << std::endl;

		if ( *list1iter == *list2iter )
		{
			//std::cerr << "set_diff drop: " << *list1iter << std::endl;			
			typename std::list< T >::iterator list1iter_next = list1iter;
			++list1iter_next;
			sortedList1.erase( list1iter );
			list1iter = list1iter_next;
			//++list2iter;
		}
		else if ( *list1iter < *list2iter )
		{
			++list1iter;
		}
		else {
			++list2iter;
		}
	}
}

template <class T >
void set_difference_apl(
	std::list< T > & sortedList1,
	std::vector< T > const & sortedList2 
)
{
	typename std::list< T >::iterator list1iter = sortedList1.begin();
	typename std::vector< T >::const_iterator list2iter = sortedList2.begin();
	while( list1iter != sortedList1.end() && list2iter != sortedList2.end()  )
	{
		//std::cerr << "set diff vect: iter1: " << *list1iter << " iter2: " << *list2iter << std::endl;

		if ( *list1iter == *list2iter )
		{
			//std::cerr << "set_diff vect drop: " << *list1iter << std::endl;			
			typename std::list< T >::iterator list1iter_next = list1iter;
			++list1iter_next;
			sortedList1.erase( list1iter );
			list1iter = list1iter_next;
			//++list2iter;
		}
		else if ( *list1iter < *list2iter )
		{
			++list1iter;
		}
		else {
			++list2iter;
		}
	}
}

//----------------------------------------------------------------------------//
//------------------------------ Vertex Class --------------------------------//
//----------------------------------------------------------------------------//


Vertex_ths::Vertex_ths(
	GraphToHoldScores * owner,
	AtomPositions * xyz,
	int index,
	int num_states 
)
:
	owner_( owner ),
	xyz_( xyz ),
	index_( index ),
	mover_( 0 ),
	//num_atoms_( -1 ),
	num_states_( num_states ),
	num_incident_d2edges_( 0 ),
	num_incident_d3edges_( 0 ),
	any_high_order_overlap_( false ),
	scores_( new float[ num_states ] ),
	penalties_( num_states, 0.0 ),
	statesEnabled_( num_states, true ),
	enabledStates2OriginalStates_( num_states, 0 )
{
	for (int ii = 0; ii < num_states_; ++ii)
	{
		scores_[ ii ] = 0;
		enabledStates2OriginalStates_[ ii ] = ii;
	}	
}

Vertex_ths::~Vertex_ths()
{
	delete [] scores_; scores_ = 0;
}

void Vertex_ths::setMover( Mover* mover )
{
	mover_ = mover;
}		

void
Vertex_ths::obtainAtomsFromMover()
{
  	std::list< AtomDescr > atoms = mover_->getAtDescOfAllPos( xyz_->getMaxFoundVDWRad() );
	atoms.sort();
	atoms.unique();
  	atoms_.resize( atoms.size() );
  	std::copy( atoms.begin(), atoms.end(), atoms_.begin() );
	
	bg_ats_overlapping_.resize( getNumAtoms() );
	bg_ats_scoring_overlapping_.resize( getNumAtoms() );
	atom_in_high_order_overlap_.resize( getNumAtoms() );
	std::fill( atom_in_high_order_overlap_.begin(), atom_in_high_order_overlap_.end(), false );	
}

int Vertex_ths::getNumStates() const
{
	return num_states_;
}

int Vertex_ths::getNumAtoms() const
{
	return atoms_.size();
}

void Vertex_ths::getAtoms( std::vector< AtomDescr > & atoms ) const
{
	assert( atoms_.size() == atoms.size() );
	std::copy( atoms_.begin(), atoms_.end(), atoms.begin() );
}


void
Vertex_ths::getOverlappingBackgroundAtoms()
{
	for( int ii = 0; ii < getNumAtoms(); ++ii )
	{	
		//std::cerr << "getOverlappingBackgroundAtoms[ " << index_ << "]:  find bumping for atom: " << atoms_[ ii ]<< std::endl;
		std::list<PDBrec*> all_bumping;
		xyz_->CollectBumping( atoms_[ ii ], all_bumping);
		
		//for( std::list< PDBrec* >::iterator iter = all_bumping.begin(); iter != all_bumping.end(); ++iter)
		//{
		//	std::cerr << atoms_[ ii ] << " bumping ALL atoms: " << (*iter)->getAtomDescr() << std::endl;
		//}
		
		std::list<PDBrec*> unboundBumping = all_bumping; //copy
		//mover_->dropBondedFromBumpingListForPDBrec( unboundBumping, atoms_[ ii ].getOriginalAtomPtr(), xyz_->getNBondCutoff() );
		 
		// for( std::list< PDBrec* >::iterator iter = unboundBumping.begin(); iter != unboundBumping.end(); ++iter)
		// {
		// 	std::cerr << atoms_[ ii ] << " bumping unbound atoms: " << (*iter)->getAtomDescr() << std::endl;
		// }
		
		for (std::list< PDBrec* >::iterator iter = all_bumping.begin();
			iter != all_bumping.end(); ++iter)
		{
			bg_ats_overlapping_[ ii ].push_back( (*iter)->getAtomDescr() );
		}
		bg_ats_overlapping_[ ii ].sort();
		bg_ats_overlapping_[ ii ].unique();
		
		
		for (std::list< PDBrec* >::iterator iter = unboundBumping.begin();
			iter != unboundBumping.end(); ++iter)
		{
			bg_ats_scoring_overlapping_[ ii ].push_back( (*iter)->getAtomDescr() );
		}
		bg_ats_scoring_overlapping_[ ii ].sort();
		bg_ats_scoring_overlapping_[ ii ].unique();
		
		//remove clique atoms from bumping list
		set_difference_apl( bg_ats_overlapping_[ ii ], atoms_ );
		set_difference_apl( bg_ats_scoring_overlapping_[ ii ], atoms_ );

		/*
		for( std::list< AtomDescr >::iterator iter = bg_ats_overlapping_[ ii ].begin(); 
			iter != bg_ats_overlapping_[ ii ].end(); ++iter)
		{
			std::cerr << atoms_[ ii ] << " bumping ALL atoms: " << *iter << std::endl;
		}		
		for( std::list< AtomDescr >::iterator iter = bg_ats_scoring_overlapping_[ ii ].begin(); 
			iter != bg_ats_scoring_overlapping_[ ii ].end(); ++iter)
		{
			std::cerr << atoms_[ ii ] << " bumping Scoring atoms: " << *iter << std::endl;
		}
		/**/
		
	}
}

void
Vertex_ths::dropMoverAtomsFromBackgroundAtomLists( Vertex_ths * other )
{
	for (int ii = 0; ii < getNumAtoms(); ++ii)
	{
		set_difference_apl( bg_ats_overlapping_[ ii ], other->atoms_ );
		set_difference_apl( bg_ats_scoring_overlapping_[ ii ], other->atoms_ );
	}
	for (int ii = 0; ii < other->getNumAtoms(); ++ii)
	{
		set_difference_apl( other->bg_ats_overlapping_[ ii ], atoms_ );
		set_difference_apl( other->bg_ats_scoring_overlapping_[ ii ], atoms_ );
	}
	//std::cerr << "After set differences: " << index_ << " and " << other->index_ << std::endl;
	//for( int ii = 0; ii < getNumAtoms(); ++ii )
	//{	
	//	for ( std::list< AtomDescr >::iterator iter = bg_ats_scoring_overlapping_[ ii ].begin();
	//		iter != bg_ats_scoring_overlapping_[ ii ].end(); ++iter )
	//	{
	//		std::cerr << atoms_[ ii ] << " bumping: " << *iter << std::endl;
	//	}
	//}
	//for (int ii = 0; ii < other->getNumAtoms(); ++ii)
	//{
	//	for ( std::list< AtomDescr >::iterator iter = other->bg_ats_scoring_overlapping_[ ii ].begin();
	//		iter != other->bg_ats_scoring_overlapping_[ ii ].end(); ++iter )
	//	{
	//		std::cerr << other->atoms_[ ii ] << " bumping: " << *iter << std::endl;
	//	}
	//}
	
}

bool
Vertex_ths::anyMoverOverlap( Vertex_ths * other ) const
{
	for (int ii = 0; ii < getNumAtoms(); ++ii )
	{
		for (int jj = 0; jj < other->getNumAtoms(); ++jj)
		{
			if ( atoms_[ ii ].intersects( other->atoms_[ jj ] ) )
			{
				return true;
			}
		}
	}
	return false;
}

void
Vertex_ths::detectBackgroundAtomsInThreeWayOverlap()
{
	std::vector< std::list< AtomDescr > > backgroundAtomsToScoreOnDeg2Edges( getNumAtoms() );
	for (std::list< DegreeTwoEdge_ths* >::iterator iter = deg2edges_.begin();
		iter != deg2edges_.end(); ++iter )
	{
		for ( int ii = 0; ii < getNumAtoms(); ++ii )
		{
			std::list< AtomDescr > temp = (*iter)->getScoringBackgroundAtomsOverlappingAtom( index_, ii );
			backgroundAtomsToScoreOnDeg2Edges[ ii ].splice(
				backgroundAtomsToScoreOnDeg2Edges[ ii ].end(),
				temp );
		}
	}
	
	for ( int ii = 0; ii < getNumAtoms(); ++ii )
	{
		backgroundAtomsToScoreOnDeg2Edges[ ii ].sort();
	}
	
	for (int ii = 0; ii < getNumAtoms(); ++ii )
	{

		std::list< AtomDescr >::iterator iter = backgroundAtomsToScoreOnDeg2Edges[ ii ].begin();
		while( iter != backgroundAtomsToScoreOnDeg2Edges[ ii ].end() )
		{
			//std::cerr << "detectBackgroundAtomsInThreeWayOverlap: " << index_ << " atom " << atoms_[ ii ] << " " << *iter << std::endl;

			std::list< AtomDescr >::iterator iternext = iter;
			++iternext;
			if ( iternext != backgroundAtomsToScoreOnDeg2Edges[ ii ].end() && *iter == *iternext )
			{
				++iternext;
				if (*iter == *iternext )
				{
					boot( ii );
				}
				else
				{
					//std::cerr << "Atom in psuedo 3-way overlap -- mover: " << index_ << " " << atoms_[ ii ] << " bg: " << *iter << std::endl;
					noteBackgroundAtomInPsuedo3WOForAtom( ii, *iter );
				}
			}
			++iter;
		}
	}
}

void
Vertex_ths::detectThreeWayOverlap( Vertex_ths * other1, Vertex_ths * other2 )
{
	assert( index_ < other1->index_ && other1->index_ < other2->index_ );
	for (int ii = 0; ii < getNumAtoms(); ++ii)
	{
		for (int jj = 0; jj < other1->getNumAtoms(); ++jj )
		{
			if ( atoms_[ ii ].intersects( other1->atoms_[ jj ] ) )
			{
				for (int kk = 0; kk < other2->getNumAtoms(); ++kk )
				{
					if ( ! atoms_[ ii ].intersects( other2->atoms_[ kk ] ) ) continue;
					if ( ! other1->atoms_[ jj ].intersects( other2->atoms_[ kk ] ) ) continue;
					
					DegreeThreeEdge_ths * d3e = owner_->getDegree3Edge( index_, other1->index_, other2->index_ );
					d3e->addMoverAtomsToScoreForAtom( index_, ii, other1->atoms_[ jj ] );
					d3e->addMoverAtomsToScoreForAtom( index_, ii, other2->atoms_[ kk ] );
					d3e->addMoverAtomsToScoreForAtom( other1->index_, jj, atoms_[ ii ] );
					d3e->addMoverAtomsToScoreForAtom( other1->index_, jj, other2->atoms_[ kk ] );
					d3e->addMoverAtomsToScoreForAtom( other2->index_, kk, atoms_[ ii ] );
					d3e->addMoverAtomsToScoreForAtom( other2->index_, kk, other1->atoms_[ jj ] );
					
					DegreeTwoEdge_ths * d2e = owner_->getDegree2Edge( index_, other1->index_ );
					d2e->addOverlappingBackgroundAtomsForPairToD3E( ii, jj, other2->atoms_[ kk ], d3e );
					
					d2e = owner_->getDegree2Edge( index_, other2->index_ );
					d2e->addOverlappingBackgroundAtomsForPairToD3E( ii, kk, other1->atoms_[ jj ], d3e );
					
					d2e = owner_->getDegree2Edge( other1->index_, other2->index_ );
					d2e->addOverlappingBackgroundAtomsForPairToD3E( jj, kk, atoms_[ ii ], d3e );
				}
			}
		}
	}
}

//atom in high order overlap -- either with other movers or with a background atom 
void
Vertex_ths::boot(
	int atom
)
{
	
	//if (! any_high_order_overlap_ ) std::cerr << "Detected high order overlap for vertex " << index_ << std::endl;	
	any_high_order_overlap_ = true;
	
	//return; //short circuit for debugging / testing how neccessary high order overlap is.
	
	atom_in_high_order_overlap_[ atom ] = true;
	std::list< int > movers_to_force_in_clique;
	for (std::list< DegreeTwoEdge_ths * >::iterator iter = deg2edges_.begin();
		iter != deg2edges_.end(); ++iter )
	{
		if ( (*iter)->atomHasOverlap( index_, atom ) )
		{
			movers_to_force_in_clique.push_back( (*iter)->getOtherNodeIndex( index_ ) );
		}
	}
	movers_to_force_in_clique.push_back( index_ );
	movers_to_force_in_clique.sort();
	owner_->forceClique( movers_to_force_in_clique );
}

std::list< DegreeTwoEdge_ths * >::iterator 
Vertex_ths::addEdge(  DegreeTwoEdge_ths * edge )
{
	deg2edges_.push_front( edge );
	return deg2edges_.begin();
}

std::list< DegreeThreeEdge_ths * >::iterator 
Vertex_ths::addEdge( DegreeThreeEdge_ths * edge)
{
	deg3edges_.push_front( edge );
	return deg3edges_.begin();
}

void Vertex_ths::dropEdge( std::list< DegreeTwoEdge_ths * >::iterator iter )
{
	deg2edges_.erase( iter );
}

void Vertex_ths::dropEdge( std::list< DegreeThreeEdge_ths * >::iterator iter )
{
	deg3edges_.erase( iter );
}

std::list< AtomDescr >
Vertex_ths::getAtomsOverlappingBothMoverAndBGAtom
(
	AtomDescr mover_atom,
	AtomDescr bg_atom
) const
{
	std::list< AtomDescr > overlaps_both;
	for (int ii = 0; ii < getNumAtoms(); ++ii )
	{
		if ( atoms_[ ii ].intersects( mover_atom ) && atoms_[ii].intersects( bg_atom ) )
		{
			overlaps_both.push_back( atoms_[ ii ] );
		}
	}
	return overlaps_both;
}

std::list< AtomDescr > const & 
Vertex_ths::getAllBumpingBackgroundAtoms( int atom ) const
{
	return bg_ats_overlapping_[ atom ];
}

std::list< AtomDescr > const &
Vertex_ths::getBumpingScoringBackgroundAtoms( int atom ) const
{
	return bg_ats_scoring_overlapping_[ atom ];
}


void Vertex_ths::informIncidentEdgesAboutBootedAtoms()
{
	if (! any_high_order_overlap_ ) return;
	for (std::list< DegreeTwoEdge_ths * >::iterator iter = deg2edges_.begin(); 
		iter != deg2edges_.end(); ++iter)
	{
		(*iter)->noteBootedAtoms( index_, atom_in_high_order_overlap_ );
	}
	for (std::list< DegreeThreeEdge_ths * >::iterator iter = deg3edges_.begin(); 
		iter != deg3edges_.end(); ++iter)
	{
		(*iter)->noteBootedAtoms( index_, atom_in_high_order_overlap_ );
	}
}

void Vertex_ths::removeBGAtomsForAtom(
	int atom,
	std::list< AtomDescr > const & bg_atoms_scored_on_edge
)
{
	set_difference_apl( bg_ats_scoring_overlapping_[ atom ], bg_atoms_scored_on_edge );
}

void Vertex_ths::finalizeScoreVectors()
{
	//std::cerr << "entering vertex_ths::finalizeScoreVectors for " << index_ << " with # atoms: " << getNumAtoms() << std::endl;
	int count_atoms_to_score = 0;
	std::vector< bool > score_atom( getNumAtoms(), false );
	for (int ii = 0; ii < getNumAtoms(); ++ii )
	{
		if ( ! atom_in_high_order_overlap_[ ii ]
			&& ! bg_ats_scoring_overlapping_[ ii ].empty() )
		{
			//std::cerr << "Score atom " << ii << " on "  << index_ << std::endl;
			score_atom[ ii ] = true;
			++count_atoms_to_score;
		}
		else
		{
			//std::cerr << "Do not score atom " << ii << " on "  << index_ << std::endl;
		}
	}
	atoms_to_score_as_vertex_.resize( count_atoms_to_score );
	int count = 0;
	for (int ii = 0; ii < getNumAtoms(); ++ii)
	{
		if ( ! score_atom[ ii ] )
		{
			bg_ats_scoring_overlapping_[ ii ].clear(); 
			continue;
		}

		/*
		std::cerr << "Atoms to score on vertex " << index_ << " for atom " << atoms_[ ii ] << std::endl;
		for( std::list< AtomDescr >::iterator iter = bg_ats_scoring_overlapping_[ ii ].begin(); 
			iter != bg_ats_scoring_overlapping_[ ii ].end(); ++iter )
		{
			std::cerr << "Atom: " << *iter << std::endl;
		}/**/
				
		atoms_to_score_as_vertex_[ count ].first = atoms_[ ii ];
		atoms_to_score_as_vertex_[ count ].second.resize( bg_ats_scoring_overlapping_[ ii ].size() );
		std::copy(
			bg_ats_scoring_overlapping_[ ii ].begin(),
			bg_ats_scoring_overlapping_[ ii ].end(),
			atoms_to_score_as_vertex_[ count ].second.begin() );
		
		bg_ats_scoring_overlapping_[ ii ].clear();
		++count;
	}

}

void
Vertex_ths::score()
{
	//std::cerr << "Scoring vertex: " << index_ << std::endl;
	mover_->initOrientation( *xyz_ );
	double penalty;
	for( int ii = 0; ii < num_states_; ++ii )
	{
		scores_[ ii ] = xyz_->determineScoreForMover( mover_, atoms_to_score_as_vertex_, penalty );
		penalties_[ ii ] = penalty;
		//std::cerr << ii << " " << scores_[ ii ]<< std::endl;
		mover_->nextOrientation( *xyz_ );
	}
	//std::cerr << "Finished scoring vertex: " << index_ << std::endl;
}

//Methods for state disabling
void
Vertex_ths::enableAllStates()
{
	std::fill( statesEnabled_.begin(), statesEnabled_.end(), true );
	numEnabledStates_ = num_states_;
	for (int ii = 0; ii < num_states_; ++ii)
	{
		enabledStates2OriginalStates_[ ii ] = ii;
	}
}

void
Vertex_ths::disableState( int state )
{
	assert( 0 <= state && state < num_states_);
	statesEnabled_[ state ] = false;
}

void
Vertex_ths::setStateDisablingComplete()
{
	numEnabledStates_ = 0;
	for (int ii = 0; ii < num_states_; ++ii )
	{
		if ( statesEnabled_[ ii ] )
		{
			enabledStates2OriginalStates_[ numEnabledStates_ ] = ii;
			++numEnabledStates_;
		}
	}
}

int
Vertex_ths::convertEnabledStateToRegularState( int enabled_state ) const
{
	assert( enabled_state >= 0 && enabled_state < numEnabledStates_ );
	//std::cerr << "convertingES2RegS[ " << enabled_state << "]--> " << enabledStates2OriginalStates_[ enabled_state ] << " on vertex " << index_ << std::endl;
	return enabledStates2OriginalStates_[ enabled_state ];
}

//Methods for interaction graph initialization
int Vertex_ths::getNumEnabledStates() const
{
	return numEnabledStates_;
}

float Vertex_ths::getScoreForEnabledState( int enabled_state ) const
{
	return scores_[ convertEnabledStateToRegularState( enabled_state ) ];
}

Mover* Vertex_ths::getMover() const
{
	return mover_;
}

bool Vertex_ths::getHasAnyHighOrderOverlap() const
{
  return any_high_order_overlap_;
  //return false;
}

std::list< AtomDescr > Vertex_ths::getAtomsInHighOrderOverlap() const
{
	std::list< AtomDescr > highOrderOverlapAtoms;
	for (int ii = 0; ii < getNumAtoms(); ++ii)
	{
		if ( atom_in_high_order_overlap_[ ii ])
		{
			highOrderOverlapAtoms.push_back( atoms_[ ii ] );
		}
	}
	return highOrderOverlapAtoms;
}

void
Vertex_ths::getPenalties( std::vector< float > & penalties ) const
{
	std::copy( penalties_.begin(), penalties_.end(), penalties.begin());
}

void Vertex_ths::assignState( int state )
{
	assert( state >= 0 && state <= num_states_ );
	assignedState_ = state;
}

float Vertex_ths::getScoreForAssignedState() const
{
	return scores_[ assignedState_ ];
}

int Vertex_ths::getAssignedState() const
{
	return assignedState_;
}

void
Vertex_ths::noteBackgroundAtomInPsuedo3WOForAtom( int mover_atom, AtomDescr bgatom )
{
	int fn_other = -1;
	int sn_other = -1;
	//std::list< AtomDescr > other_overlapping_atoms[ 2 ];
	for (std::list< DegreeTwoEdge_ths * >::iterator iter = deg2edges_.begin(); iter != deg2edges_.end(); ++iter )
	{
		if ( (*iter)->bgAtomOverlapsMovingAtomOnEdge( index_, mover_atom, bgatom ) )
		{
			Vertex_ths * other = (*iter)->getOtherNode( index_ );
			if ( fn_other == -1)
			{
				fn_other = (*iter)->getOtherNodeIndex( index_ );
				//other_overlapping_atoms[ 0 ] = other->getAtomsOverlappingBothMoverAndBGAtom( atoms_[ mover_atom ], bgatom );
			}
			else
			{
				sn_other = (*iter)->getOtherNodeIndex( index_ );
				//other_overlapping_atoms[ 1 ] = other->getAtomsOverlappingBothMoverAndBGAtom( atoms_[ mover_atom ], bgatom );
				break;
			}
		}
	}
	assert( sn_other != -1 );
	
	int fn( -1 ), sn( -1 ), tn( -1 );
	sort_three( index_, fn_other, sn_other, fn, sn, tn );

	DegreeThreeEdge_ths * d3e = owner_->getDegree3Edge( fn, sn, tn );
	d3e->addBackgroundAtomToScoreForMoverAtom( index_, mover_atom, bgatom );
	//for (int ii = 0; ii < 2; ++ii)
	//{
	//	for (std::list< AtomDescr >::iterator iter = other_overlapping_atoms[ ii ].begin();
	//		iter != other_overlapping_atoms[ ii ].end(); ++iter)
	//	{
	//		d3e->addMoverAtomsToScoreForAtom( index_, mover_atom, *iter );
	//	}
	//}
}

//----------------------------------------------------------------------------//
//----------------------------- Degree Two Edge ------------------------------//
//----------------------------------------------------------------------------//


DegreeTwoEdge_ths::DegreeTwoEdge_ths(
	GraphToHoldScores * owner,
	AtomPositions * xyz,
	int vertex1,
	int vertex2
)
:
	owner_( owner ),
	xyz_( xyz )
{
	vertex_indices_[ 0 ] = vertex1;
	vertex_indices_[ 1 ] = vertex2;
	for (int ii = 0; ii < 2; ++ii)
	{
		vertex_ptrs_[ ii ] = owner_->getVertexPtr( vertex_indices_[ ii ] );
		pos_in_verts_edge_list_[ ii ] = vertex_ptrs_[ ii ]->addEdge( this );
		num_states_[ ii ] = vertex_ptrs_[ ii ]->getNumStates();
		num_atoms_[ ii ] = vertex_ptrs_[ ii ]->getNumAtoms();
		movers_[ ii ] = vertex_ptrs_[ ii ]->getMover();
		bg_ats_scoring_for_atom_[ ii ].resize( num_atoms_[ ii ] );
		bg_ats_all_for_atom_[ ii ].resize( num_atoms_[ ii ] );
		mover_atoms_atom_overlaps_[ ii ].resize( num_atoms_[ ii ] );
		bg_atoms_overlapping_pair_[ ii ].resize( num_atoms_[ ii ] );
		atoms_[ ii ].resize( num_atoms_[ ii ] );
		vertex_ptrs_[ ii ]->getAtoms( atoms_[ ii ] );
	}
	scores_ = new float * [ num_states_[ 0 ] ];
	for (int ii = 0; ii < num_states_[ 0 ]; ++ii )
	{
		scores_[ ii ] = new float[ num_states_[ 1 ] ];
		for (int jj = 0; jj < num_states_[ 1 ]; ++jj )
		{
			scores_[ ii ][ jj ] = 0;
		}
	}
}

DegreeTwoEdge_ths::~DegreeTwoEdge_ths()
{
	if ( scores_ != 0 )
	{
		for (int ii = 0; ii < num_states_[ 0 ]; ++ii )
		{
			delete [] scores_[ ii ]; scores_[ ii ] = 0;
		}
		delete [] scores_; scores_ = 0;
	}
	vertex_ptrs_[ 0 ]->dropEdge( pos_in_verts_edge_list_[ 0 ] );
	vertex_ptrs_[ 1 ]->dropEdge( pos_in_verts_edge_list_[ 1 ] );
	owner_->dropEdge( pos_in_owners_edge_list_ );

}

void
DegreeTwoEdge_ths::setPosInOwnersEdgeList( std::list< DegreeTwoEdge_ths * >::iterator iter )
{
	pos_in_owners_edge_list_ = iter;
}

bool DegreeTwoEdge_ths::sameEdge( int vertex1, int vertex2 ) const
{
	return ( vertex1 == vertex_indices_[ 0 ] && vertex2 == vertex_indices_[ 1 ] );
}

Vertex_ths * DegreeTwoEdge_ths::getOtherNode( int vertex_index ) const
{
	int other = whichVertex( vertex_index ) == 0 ? 1 : 0;
	return vertex_ptrs_[ other ];
}

int DegreeTwoEdge_ths::getFirstNodeIndex() const
{
	return vertex_indices_[ 0 ];
}

int DegreeTwoEdge_ths::getSecondNodeIndex() const
{
	return vertex_indices_[ 1 ];
}

int DegreeTwoEdge_ths::getOtherNodeIndex( int vertex_index ) const
{
	int other = whichVertex( vertex_index ) == 0 ? 1 : 0;
	return vertex_indices_[ other ];
}

bool DegreeTwoEdge_ths::incidentUpon( int vertex_index ) const
{
	return vertex_indices_[ 0 ] == vertex_index || vertex_indices_[ 1 ] == vertex_index;
}

void
DegreeTwoEdge_ths::detectBackgroundAtomsToScoreForEdge()
{
	//std::cerr << "detectBackgroundAtomsToScoreForEdge(); " << vertex_indices_[ 0 ] << " " << vertex_indices_[ 1 ] << std::endl;
	//std::cerr << num_atoms_[ 0 ] << " and " << num_atoms_[ 1 ] << " atoms" << std::endl;
	//std::cerr << atoms_[0].size() << " and " << atoms_[1].size() << std::endl;
	for (int ii = 0; ii < num_atoms_[ 0 ]; ++ii )
	{
		//std::cerr << "ii atom: " << atoms_[0][ii] << std::endl;
		for (int jj = 0; jj < num_atoms_[ 1 ]; ++jj )
		{
			//std::cerr << "jj atom: " << atoms_[1][jj] << std::endl;
			//std::cerr << "ii: " << ii << " jj: " << jj << std::endl;
			if ( atoms_[ 0 ][ ii ].intersects( atoms_[ 1 ][ jj ] ) )
			{
				//std::cerr << "overlap!" << std::endl;
				mover_atoms_atom_overlaps_[ 0 ][ ii ].push_back( atoms_[ 1 ][ jj ]);
				mover_atoms_atom_overlaps_[ 1 ][ jj ].push_back( atoms_[ 0 ][ ii ]);
				
				std::list< AtomDescr > const & ii_overlapping_bg = vertex_ptrs_[ 0 ]->getAllBumpingBackgroundAtoms( ii );
				std::list< AtomDescr > const & jj_overlapping_bg = vertex_ptrs_[ 1 ]->getAllBumpingBackgroundAtoms( jj );

				// for (std::list< AtomDescr >::const_iterator iter = ii_overlapping_bg.begin();
				// 	iter != ii_overlapping_bg.end(); ++iter )
				// {
				// 	std::cerr << "ii_overlapping_bg: " << *iter << " for " << ii << " and " << jj << std::endl;
				// }
				// for (std::list< AtomDescr >::const_iterator iter = jj_overlapping_bg.begin();
				// 	iter != jj_overlapping_bg.end(); ++iter )
				// {
				// 	std::cerr << "jj_overlapping_bg: " << *iter << " for " << ii << " and " << jj << std::endl;
				// }

				//std::cerr << "about to call set_intersect_apl" << std::endl;				
				std::list< AtomDescr > bg_intersecting_both;	
				set_intersect_apl( ii_overlapping_bg, jj_overlapping_bg, bg_intersecting_both );
				
				// for (std::list< AtomDescr >::iterator iter = bg_intersecting_both.begin();
				// 	iter != bg_intersecting_both.end(); ++iter )
				// {
				// 	std::cerr << "bg_intersecting_both: " << *iter << " for " << ii << " and " << jj << std::endl;
				// }
								
				//bg_atoms_overlapping_pair_[ 0 ][ ii ].push_back( bg_intersecting_both );
				//bg_atoms_overlapping_pair_[ 1 ][ jj ].push_back( bg_intersecting_both );

				std::list< AtomDescr > const & ii_scoring_overlapping_bg = vertex_ptrs_[ 0 ]->getBumpingScoringBackgroundAtoms( ii );

				std::list< AtomDescr > ii_bg_score_on_edge;
				set_intersect_apl( ii_scoring_overlapping_bg, bg_intersecting_both, ii_bg_score_on_edge );
				bg_atoms_overlapping_pair_[ 0 ][ ii ].push_back( ii_bg_score_on_edge );
				//std::cerr << "push back 1 " << std::endl;
				for (std::list< AtomDescr >::iterator iter = ii_bg_score_on_edge.begin();
					iter != ii_bg_score_on_edge.end(); ++iter )
				{
					bg_ats_scoring_for_atom_[ 0 ][ ii ].push_back( *iter );
				}
				for (std::list< AtomDescr >::iterator iter = bg_intersecting_both.begin();
					iter != bg_intersecting_both.end(); ++iter )
				{
					bg_ats_all_for_atom_[ 0 ][ ii ].push_back( *iter );
				}

				// for (std::list< AtomDescr >::iterator iter = bg_intersecting_both.begin();
				// 	iter != bg_intersecting_both.end(); ++iter )
				// {
				// 	std::cerr << "bg_intersecting_both: " << *iter << " for " << ii << " and " << jj << std::endl;
				// }
				//std::cerr << "push back 2 " << std::endl;
				
				std::list< AtomDescr > const & jj_scoring_overlapping_bg = vertex_ptrs_[ 1 ]->getBumpingScoringBackgroundAtoms( jj );
				std::list< AtomDescr > jj_bg_score_on_edge;
				set_intersect_apl( jj_scoring_overlapping_bg, bg_intersecting_both, jj_bg_score_on_edge );
				bg_atoms_overlapping_pair_[ 1 ][ jj ].push_back( jj_bg_score_on_edge );

				for (std::list< AtomDescr >::iterator iter = jj_bg_score_on_edge.begin();
					iter != jj_bg_score_on_edge.end(); ++iter )
				{
					bg_ats_scoring_for_atom_[ 1 ][ jj ].push_back( *iter );
				}
				for (std::list< AtomDescr >::iterator iter = bg_intersecting_both.begin();
					iter != bg_intersecting_both.end(); ++iter )
				{
					bg_ats_all_for_atom_[ 1 ][ jj ].push_back( *iter );
				}

				//std::cerr << "going to next iteration" << std::endl;
			}
		}
	}
	
	for (int ii = 0; ii < 2; ++ii )
	{
		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj )
		{
			bg_ats_scoring_for_atom_[ ii ][ jj ].sort();
			bg_ats_scoring_for_atom_[ ii ][ jj ].unique();
			bg_ats_all_for_atom_[ ii ][ jj ].sort();
			bg_ats_all_for_atom_[ ii ][ jj ].unique();
			//for (std::list< AtomDescr >::iterator iter = bg_ats_scoring_for_atom_[ ii ][ jj ].begin();
			//	iter != bg_ats_scoring_for_atom_[ ii ][ jj ].end(); ++iter )
			//{
			//	std::cerr << "bg_ats_scoring_for_atom_[ ii ][ jj ]: " << *iter << " for mover " << vertex_indices_[ ii ] << " atom " << atoms_[ii][jj] << std::endl;
			//}
			//for (std::list< AtomDescr >::iterator iter = mover_atoms_atom_overlaps_[ ii ][ jj ].begin();
			//	iter != mover_atoms_atom_overlaps_[ ii ][ jj ].end(); ++iter )
			//{
			//	std::cerr << "mover_atom_atom_overlaps_[ ii ][ jj ]: " << *iter << " for mover " << vertex_indices_[ ii ]  << " atom " << atoms_[ii][jj] << std::endl;
			//}
		}
	}
}

void
DegreeTwoEdge_ths::addOverlappingBackgroundAtomsForPairToD3E(
	int atom_on_vertex1,
	int atom_on_vertex2,
	AtomDescr atom_on_vertex3,
	DegreeThreeEdge_ths * d3e
)
{
	std::list< std::list< AtomDescr > >::iterator bg_list_iter1 =
		bg_atoms_overlapping_pair_[ 0 ][ atom_on_vertex1 ].begin();
		
	for (std::list< AtomDescr >::iterator iter = mover_atoms_atom_overlaps_[ 0 ][ atom_on_vertex1 ].begin();
		iter != 	mover_atoms_atom_overlaps_[ 0 ][ atom_on_vertex1 ].end(); ++iter )
	{
		if ( *iter == atoms_[ 1 ][ atom_on_vertex2 ] )
		{
			for( std::list< AtomDescr >::iterator bg_atom_iter = bg_list_iter1->begin();
				bg_atom_iter != bg_list_iter1->end(); ++bg_atom_iter )
			{
				if ( atom_on_vertex3.intersects( *bg_atom_iter ) )
				{
					d3e->addBackgroundAtomToScoreForMoverAtom( vertex_indices_[ 0 ], atom_on_vertex1, *bg_atom_iter );
				}
				//d3e->addBackgroundAtomToScoreForMoverAtom( vertex_indices_[ 1 ], atom_on_vertex2, *bg_atom_iter );
			}
			break;
		}
		++bg_list_iter1;
	}
	
	std::list< std::list< AtomDescr > >::iterator bg_list_iter2 =
	bg_atoms_overlapping_pair_[ 1 ][ atom_on_vertex2 ].begin();
		
	for (std::list< AtomDescr >::iterator iter = mover_atoms_atom_overlaps_[ 1 ][ atom_on_vertex2 ].begin();
		iter != 	mover_atoms_atom_overlaps_[ 1 ][ atom_on_vertex2 ].end(); ++iter )
	{
		if ( *iter == atoms_[ 0 ][ atom_on_vertex1 ] )
		{
			for( std::list< AtomDescr >::iterator bg_atom_iter = bg_list_iter2->begin();
				bg_atom_iter != bg_list_iter2->end(); ++bg_atom_iter )
			{
				if ( atom_on_vertex3.intersects( *bg_atom_iter ) )
				{ 
					d3e->addBackgroundAtomToScoreForMoverAtom( vertex_indices_[ 1 ], atom_on_vertex2, *bg_atom_iter );
				}
			}
			break;
		}
		++bg_list_iter2;
	}
}

std::list< AtomDescr >
DegreeTwoEdge_ths::getScoringBackgroundAtomsOverlappingAtom( int vertex_index, int atom )
{
	return bg_ats_scoring_for_atom_[ whichVertex( vertex_index ) ][ atom ];
}

bool
DegreeTwoEdge_ths::bgAtomOverlapsMovingAtomOnEdge( int vertex_index, int atom, AtomDescr bgatom ) const
{
	return ( std::find( bg_ats_all_for_atom_[ whichVertex( vertex_index ) ][ atom ].begin(),
		bg_ats_all_for_atom_[ whichVertex( vertex_index ) ][ atom ].end(),
		bgatom)  != bg_ats_all_for_atom_[ whichVertex( vertex_index ) ][ atom ].end() );
}

bool
DegreeTwoEdge_ths::atomHasOverlap( int vertex_index, int atom ) const
{
	return ( ! bg_ats_scoring_for_atom_[ whichVertex( vertex_index ) ][ atom ].empty() 
		|| ! mover_atoms_atom_overlaps_[ whichVertex( vertex_index ) ][ atom ].empty() );
}

void
DegreeTwoEdge_ths::noteBootedAtoms(
	int vertex_index,
	std::vector< bool  > const & booted_atoms
)
{
	int whichvertex = whichVertex( vertex_index );
	for ( int ii = 0; ii < num_atoms_[ whichvertex ]; ++ii )
	{
		if ( booted_atoms[ ii ] )
		{
			bg_ats_scoring_for_atom_[ whichvertex ][ ii ].clear();
			bg_ats_all_for_atom_[ whichvertex ][ ii ].clear();
			mover_atoms_atom_overlaps_[ whichvertex ][ ii ].clear();
		}
	}
}

void DegreeTwoEdge_ths::tellVerticesToNotScoreBGAtomsScoredOnD2Edge() const
{
	for (int ii = 0; ii < 2; ++ii)
	{
		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj )
		{
			//std::cerr << "BG Atoms scored on edge for node " << vertex_indices_[ ii ] << " for atom " << atoms_[ ii ][ jj ] << std::endl;
			//for ( std::list< AtomDescr >::const_iterator iter = bg_ats_scoring_for_atom_[ ii ][ jj ].begin();
			//	iter != bg_ats_scoring_for_atom_[ ii ][ jj ].end(); ++iter)
			//{
			//	std::cerr << "atom: " << *iter << std::endl;
			//}
			vertex_ptrs_[ ii ]->removeBGAtomsForAtom( jj, bg_ats_scoring_for_atom_[ ii ][ jj ]);
		}
	}
}

void
DegreeTwoEdge_ths::dropBGAtomsOverlappingAtom(
	int vertex_index,
	int atom,
	std::list< AtomDescr > const & bgats
)
{
	set_difference_apl( 
		bg_ats_scoring_for_atom_[ whichVertex( vertex_index ) ][ atom ],
		bgats
	);	
}

void
DegreeTwoEdge_ths::dropMoverAtomsOverlappingAtom(
	int vertex_index,
	int atom,
	std::list< AtomDescr > const & moverats
)
{
	set_difference_apl( 
		mover_atoms_atom_overlaps_[ whichVertex( vertex_index ) ][ atom ],
		moverats
	);	
}

void DegreeTwoEdge_ths::finalizeScoreVectors()
{
	for (int ii = 0; ii < 2; ++ii )
	{
		int count_atoms_to_score = 0;
		std::vector< bool > score_atom( num_atoms_[ ii ], false );
		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj )
		{
			if ( ! bg_ats_scoring_for_atom_[ ii ][ jj ].empty()
				|| ! mover_atoms_atom_overlaps_[ ii ][ jj ].empty() )
			{
				score_atom[ jj ] = true;
				++count_atoms_to_score;
			}
		}
	
		atoms_to_score_[ ii ].resize( count_atoms_to_score );
		int count = 0;
		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj)
		{
			if ( ! score_atom[ jj ] ) continue;
			std::list< AtomDescr > all_atoms_to_score;
			all_atoms_to_score.splice(
				all_atoms_to_score.end(),
				bg_ats_scoring_for_atom_[ ii ][ jj ] );
			all_atoms_to_score.splice(
				all_atoms_to_score.end(),
				mover_atoms_atom_overlaps_[ ii ][ jj ]
			);
			all_atoms_to_score.sort();
			
			/*
			std::cerr << "Atoms To score on edge (" << vertex_indices_[ 0 ] << ", " << vertex_indices_[ 1 ];
			std::cerr << " for v: " << vertex_indices_[ ii ] << " for at: " << atoms_[ ii ][ jj ] << std::endl;
			for ( std::list< AtomDescr >::iterator iter = all_atoms_to_score.begin();
				iter != all_atoms_to_score.end(); ++iter)
			{
				std::cerr << "atom: " << *iter << std::endl;
			}
			/**/
				
			atoms_to_score_[ ii ][ count ].first = atoms_[ ii ][ jj ];
			atoms_to_score_[ ii ][ count ].second.resize( all_atoms_to_score.size() );
			std::copy(
				all_atoms_to_score.begin(),
				all_atoms_to_score.end(),
				atoms_to_score_[ ii ][ count ].second.begin() );

			++count;
		}
	}
}

void
DegreeTwoEdge_ths::score()
{
	//std::cerr << "Scoring edge: " << vertex_indices_[ 0 ] << " " << vertex_indices_[ 1 ] << std::endl;
	if ( nothingToScore() ) return;
	movers_[ 0 ]->initOrientation( *xyz_ );
	double temp;
	for( int ii = 0; ii < num_states_[ 0 ]; ++ii )
	{
		//std::cerr << "ii: " << ii << std::endl;
		movers_[ 1 ]->initOrientation( *xyz_ );
		for (int jj = 0; jj < num_states_[ 1 ]; ++jj )
		{
			//std::cerr << "jj: " << jj << std::endl;
			for (int kk = 0; kk < 2; ++kk )
			{
				if ( nothingToScore( kk ) ) continue;
				scores_[ ii ][ jj ] += xyz_->determineScoreForMover( movers_[ kk ], atoms_to_score_[ kk ], temp );
			}
			//std::cerr << ii << " " << jj << " " << scores_[ ii ][ jj ] << std::endl;
			movers_[ 1 ]->nextOrientation( *xyz_ );
		}
		movers_[ 0 ]->nextOrientation( *xyz_ );
	}
	//std::cerr << "Finished scoring edge: " << vertex_indices_[ 0 ] << " " << vertex_indices_[ 1 ] << std::endl;
}

float
DegreeTwoEdge_ths::getScoreForEnabledStates(
	int enabled_state1,
	int enabled_state2 ) const
{
	int states[] = {enabled_state1, enabled_state2};
	int original_states[2];
	for (int ii = 0; ii < 2; ++ii)
	{
		original_states[ ii ] = vertex_ptrs_[ ii ]->convertEnabledStateToRegularState( states[ ii ] );
	}
	return scores_[ original_states[ 0 ] ][ original_states[ 1 ] ];
}

float DegreeTwoEdge_ths::getScoreForAssignedState() const
{
	return scores_[ vertex_ptrs_[ 0 ]->getAssignedState() ][ vertex_ptrs_[ 1 ]->getAssignedState() ];
}


int
DegreeTwoEdge_ths::whichVertex( int vertex_index ) const
{
	for (int ii = 0; ii < 2; ++ii )
	{
		if (vertex_indices_[ ii ] == vertex_index )
		{
			return ii;
		}
	}
	std::cerr << "CRITICAL ERROR IN whichVertex( " << vertex_index << ") called on edge [" <<
	std::cerr << vertex_indices_[ 0 ] << ", " << vertex_indices_[ 1 ] << "]" << std::endl;
	assert(false);
	exit(1);
}

bool DegreeTwoEdge_ths::nothingToScore() const
{
	//std::cerr << "nothingToScore?" << std::endl;
	for (int ii = 0; ii < 2; ++ii )
	{
		if ( ! atoms_to_score_[ ii ].empty() ) return false;
	}
	//std::cerr << "Indeed, nothing to score" << std::endl;
	return true;
}

bool DegreeTwoEdge_ths::nothingToScore( int vertex ) const
{
	return atoms_to_score_[ vertex ].empty();
}

//----------------------------------------------------------------------------//
//------------------------- Degree Three Edge Class --------------------------//
//----------------------------------------------------------------------------//

DegreeThreeEdge_ths::DegreeThreeEdge_ths(
	GraphToHoldScores * owner,
	AtomPositions * xyz,
	int vertex1,
	int vertex2,
	int vertex3
)
:
	owner_( owner ),
	xyz_( xyz )
{
	vertex_indices_[ 0 ] = vertex1;
	vertex_indices_[ 1 ] = vertex2;
	vertex_indices_[ 2 ] = vertex3;
	//std::cerr << "Creating degree three edge";
	//for (int ii = 0; ii < 3; ++ii)
	//{
	//	std::cerr << " " << vertex_indices_[ ii ];
	//}
	//std::cerr << std::endl;
	
	for (int ii = 0; ii < 3; ++ii)
	{
		vertex_ptrs_[ ii ] = owner_->getVertexPtr( vertex_indices_[ ii ] );
		pos_in_verts_edge_list_[ ii ] = vertex_ptrs_[ ii ]->addEdge( this );
		num_states_[ ii ] = vertex_ptrs_[ ii ]->getNumStates();
		num_atoms_[ ii ] = vertex_ptrs_[ ii ]->getNumAtoms();
		bg_ats_for_atom_[ ii ].resize( num_atoms_[ ii ] );
		mover_atoms_atom_overlaps_[ ii ].resize( num_atoms_[ ii ] );
		movers_[ ii ] = vertex_ptrs_[ ii ]->getMover();
		atoms_[ ii ].resize( num_atoms_[ ii ] );
		vertex_ptrs_[ ii ]->getAtoms( atoms_[ ii ] );
	}
	d2e_[ 0 ] = owner_->getDegree2Edge( vertex1, vertex2 );
	d2e_[ 1 ] = owner_->getDegree2Edge( vertex1, vertex3 );
	d2e_[ 2 ] = owner_->getDegree2Edge( vertex2, vertex3 );
	
	scores_table_ = new float**[ num_states_[ 0 ] ];	
	for (int ii = 0; ii < num_states_[ 0 ]; ++ii)
	{
		scores_table_[ ii ] = new float*[ num_states_[ 1 ] ];
		for (int jj = 0; jj < num_states_[ 1 ]; ++jj )
		{
			scores_table_[ ii ][ jj ] = new float[ num_states_[ 2 ] ];
			for (int kk = 0; kk < num_states_[ 2 ]; ++kk )
			{
				scores_table_[ ii ][ jj ][ kk ] = 0;
			}
		}
	}
}	

DegreeThreeEdge_ths::~DegreeThreeEdge_ths()
{
	if (scores_table_ != 0 ) 
	{	
		for (int ii = 0; ii < num_states_[ 0 ]; ++ii)
		{
			if ( scores_table_[ ii ] == 0 ) continue;
			for (int jj = 0; jj < num_states_[ 1 ]; ++jj )
			{
				delete [] scores_table_[ ii ][ jj ];scores_table_[ ii ][ jj ] = 0;
			}
			delete [] scores_table_[ ii ];scores_table_[ ii ] = 0;
		}
		delete [] scores_table_;scores_table_ = 0;
	}
	vertex_ptrs_[ 0 ]->dropEdge( pos_in_verts_edge_list_[ 0 ] );
	vertex_ptrs_[ 1 ]->dropEdge( pos_in_verts_edge_list_[ 1 ] );
	vertex_ptrs_[ 2 ]->dropEdge( pos_in_verts_edge_list_[ 2 ] );
	owner_->dropEdge( pos_in_owners_edge_list_ );
}

void
DegreeThreeEdge_ths::setPosInOwnersEdgeList( std::list< DegreeThreeEdge_ths * >::iterator iter )
{
	pos_in_owners_edge_list_ = iter;
}

bool DegreeThreeEdge_ths::sameEdge( int vertex1, int vertex2, int vertex3 ) const
{
	return ( vertex1 == vertex_indices_[ 0 ] && vertex2 == vertex_indices_[ 1 ] && vertex3 == vertex_indices_[ 2 ] );
}

int DegreeThreeEdge_ths::getFirstNodeIndex() const
{
	return vertex_indices_[ 0 ];
}

int DegreeThreeEdge_ths::getSecondNodeIndex() const
{
	return vertex_indices_[ 1 ];
}

int DegreeThreeEdge_ths::getThirdNodeIndex() const
{
	return vertex_indices_[ 2 ];
}


void
DegreeThreeEdge_ths::addBackgroundAtomToScoreForMoverAtom(
	int vertex,
	int atom,
	AtomDescr atom_overlapping
)
{
	int which_of_three = whichVertex( vertex );
	bg_ats_for_atom_[ which_of_three ][ atom ].push_back( atom_overlapping );
	bg_ats_for_atom_[ which_of_three ][ atom ].sort();
	bg_ats_for_atom_[ which_of_three ][ atom ].unique();
}

void
DegreeThreeEdge_ths::addMoverAtomsToScoreForAtom(
	int vertex,
	int atom,
	AtomDescr atom_overlapping
)
{
	int which_of_three = whichVertex( vertex );
	mover_atoms_atom_overlaps_[ which_of_three ][ atom ].push_back( atom_overlapping );
	mover_atoms_atom_overlaps_[ which_of_three ][ atom ].sort();
	mover_atoms_atom_overlaps_[ which_of_three ][ atom ].unique();
}

bool
DegreeThreeEdge_ths::sharesTwoVertices( DegreeThreeEdge_ths const * other ) const
{
	if ( vertex_indices_[ 0 ] == other->vertex_indices_[ 0 ] )
	{
		if (vertex_indices_[ 1 ] == other->vertex_indices_[ 1 ] ||
			vertex_indices_[ 1 ] == other->vertex_indices_[ 2 ] ||
			vertex_indices_[ 2 ] == other->vertex_indices_[ 2 ] ||
			vertex_indices_[ 2 ] == other->vertex_indices_[ 1 ])
		{
			return true;
		}
	}
	else if (vertex_indices_[ 0 ] == other->vertex_indices_[ 1 ] )
	{
		if ( vertex_indices_[ 1 ] == other->vertex_indices_[ 2 ] ||
			vertex_indices_[ 2 ] == other->vertex_indices_[ 2 ] )
		{
			return true;
		}
	}
	else if (vertex_indices_[ 1 ] == other->vertex_indices_[ 1 ]  &&
		vertex_indices_[ 2 ] == other->vertex_indices_[ 2 ] )
	{
		return true;
	}
	
	return false;
	
}

void
DegreeThreeEdge_ths::detect_shared_atom_pair( DegreeThreeEdge_ths * other ) const
{
	if ( ! sharesTwoVertices( other ) ) return;
	int this_fn_shared = -1;
	int this_sn_shared = -1;
	int other_fn_shared = -1;
	int other_sn_shared = -1;
	
	//std::cerr << "detect_shared_atom_pair comparing";
	//for (int ii = 0; ii < 3; ++ii)
	//{
	//	std::cerr << " " << vertex_indices_[ ii ];
	//}
	//std::cerr << " and";
	//for (int ii = 0; ii < 3; ++ii)
	//{
	//	std::cerr << " " << other->vertex_indices_[ ii ];
	//}
	//std::cerr << std::endl;
	
	if ( vertex_indices_[ 0 ] == other->vertex_indices_[ 0 ] )
	{
		this_fn_shared = 0;
		other_fn_shared = 0;
		if (vertex_indices_[ 1 ] == other->vertex_indices_[ 1 ] )
		{
			this_sn_shared = 1;
			other_sn_shared = 1;
		}
		else if ( vertex_indices_[ 1 ] == other->vertex_indices_[ 2 ] )
		{
			this_sn_shared = 1;
			other_sn_shared = 2;
		}
		else if ( vertex_indices_[ 2 ] == other->vertex_indices_[ 2 ] ) 
		{
			this_sn_shared = 2;
			other_sn_shared = 2;
		}
		else if ( vertex_indices_[ 2 ] == other->vertex_indices_[ 1 ] )
		{
			this_sn_shared = 2;
			other_sn_shared = 1;
		}
	}
	else if (vertex_indices_[ 0 ] == other->vertex_indices_[ 1 ] )
	{
		this_fn_shared = 0;
		other_fn_shared = 1;
		if ( vertex_indices_[ 1 ] == other->vertex_indices_[ 2 ] )
		{
			this_sn_shared = 1;
			other_sn_shared = 2;
		}
		else if (vertex_indices_[ 2 ] == other->vertex_indices_[ 2 ] )
		{
			this_sn_shared = 2;
			other_sn_shared = 2;
		}
	}
	else if (vertex_indices_[ 1 ] == other->vertex_indices_[ 1 ] &&
		vertex_indices_[ 2 ] == other->vertex_indices_[ 2 ] )
	{
		this_fn_shared = 1; 
		other_fn_shared = 1;
		this_sn_shared = 2;
		other_sn_shared = 2;
	}
	else
	{
		std::cerr << "CRITICAL ERROR IN detect_shared_atom_pair: between edges [" << 
		std::cerr << vertex_indices_[ 0 ] << ", ";
		std::cerr << vertex_indices_[ 1 ] << ", ";
		std::cerr << vertex_indices_[ 2 ] << "] and [";
		std::cerr << other->vertex_indices_[ 0 ] << ", ";
		std::cerr << other->vertex_indices_[ 1 ] << ", ";
		std::cerr << other->vertex_indices_[ 2 ] << "]" << std::endl;
	}
	// std::cerr << "found:";
	// std::cerr << " this_fn_shared " << this_fn_shared;
	// std::cerr << " this_sn_shared " << this_sn_shared;
	// std::cerr << " other_fn_shared " << other_fn_shared;
	// std::cerr << " other_sn_shared " << other_sn_shared << std::endl;
	
	compareSharedAtomPairs( other, this_fn_shared, other_fn_shared );
	compareSharedAtomPairs( other, this_sn_shared, other_sn_shared );
}

void DegreeThreeEdge_ths::detectBackgroundAtomScoredTwiceForAtom
( 
	DegreeThreeEdge_ths const * other
) const
{
	for (int ii = 0; ii < 3; ++ii )
	{
		for (int jj = 0; jj < 3; ++jj )
		{
			if ( vertex_indices_[ ii ] == other->vertex_indices_[ jj ] )
			{
				compareSharedBackgroundAtoms( other, ii, jj );
			}
		}
	}
	
}

void 
DegreeThreeEdge_ths::tellD2EdgesToNotScoreBGandMoverAtomsScoredOnD3Edge() const
{
	for (int ii = 0; ii < 3; ++ii )
	{
		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj )
		{
			if ( bg_ats_for_atom_[ ii ][ jj ].empty() && mover_atoms_atom_overlaps_[ ii ][ jj ].empty() )
			{
				continue;
			}
			for (int kk = 0; kk < 3; ++kk )
			{
				if ( ! d2e_[ kk ]->incidentUpon( vertex_indices_[ ii ] ) )
				{ 
					continue;
				}
				d2e_[ kk ]->dropBGAtomsOverlappingAtom(
					vertex_indices_[ ii ],
					jj,
					bg_ats_for_atom_[ ii ][ jj ]
				);
				d2e_[ kk ]->dropMoverAtomsOverlappingAtom(
					vertex_indices_[ ii ],
					jj,
					mover_atoms_atom_overlaps_[ ii ][ jj ]
				);
			}
		}
	}
}

void
DegreeThreeEdge_ths::noteBootedAtoms(
	int vertex_index,
	std::vector< bool  > const & booted_atoms
)
{
	int whichvertex = whichVertex( vertex_index );
	for ( int ii = 0; ii < num_atoms_[ whichvertex ]; ++ii )
	{
		if ( booted_atoms[ ii ] )
		{
			bg_ats_for_atom_[ whichvertex ][ ii ].clear();
			mover_atoms_atom_overlaps_[ whichvertex ][ ii ].clear();
		}
	}
}

void DegreeThreeEdge_ths::finalizeScoreVectors()
{
	for (int ii = 0; ii < 3; ++ii )
	{
		int count_atoms_to_score = 0;
		std::vector< bool > score_atom( num_atoms_[ ii ], false );
		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj )
		{
			if ( ! bg_ats_for_atom_[ ii ][ jj ].empty()
				|| ! mover_atoms_atom_overlaps_[ ii ][ jj ].empty() )
			{
				score_atom[ jj ] = true;
				++count_atoms_to_score;
			}
		}
	
		atoms_to_score_[ ii ].resize( count_atoms_to_score );
		int count = 0;
		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj)
		{
			if ( ! score_atom[ jj ] ) continue;
			std::list< AtomDescr > all_atoms_to_score;
			all_atoms_to_score.splice(
				all_atoms_to_score.end(),
				bg_ats_for_atom_[ ii ][ jj ] );
			all_atoms_to_score.splice(
				all_atoms_to_score.end(),
				mover_atoms_atom_overlaps_[ ii ][ jj ]
			);
			all_atoms_to_score.sort();
				
			atoms_to_score_[ ii ][ count ].first = atoms_[ ii ][ jj ];
			atoms_to_score_[ ii ][ count ].second.resize( all_atoms_to_score.size() );
			std::copy(
				all_atoms_to_score.begin(),
				all_atoms_to_score.end(),
				atoms_to_score_[ ii ][ count ].second.begin() );
			
			/*
			for ( int kk = 0; kk < atoms_to_score_[ ii ][ count ].second.size(); ++kk)
			{
				std::cerr << "ats2score: v" << vertex_indices_[ ii ] << " on d3e";
				for (int ll = 0; ll < 3; ++ll)
				{
					std::cerr << " " << vertex_indices_[ ll ];
				}
				std::cerr << " atom " << atoms_to_score_[ ii ][ count ].first;
				std::cerr << " with " <<  (atoms_to_score_[ ii ][ count ].second)[ kk ] << std::endl;
			}
			/**/


			++count;
		}
	}

}

void
DegreeThreeEdge_ths::score()
{
	if ( nothingToScore() ) return;
	//std::cerr << "Scoring degree 3 hyperedge: " << vertex_indices_[ 0 ] << " ";
	//std::cerr << vertex_indices_[ 1 ] << " ";
	//std::cerr << vertex_indices_[ 2 ] << std::endl;
	movers_[ 0 ]->initOrientation( *xyz_ );
	double temp;
	for( int ii = 0; ii < num_states_[ 0 ]; ++ii )
	{
		movers_[ 1 ]->initOrientation( *xyz_ );
		for (int jj = 0; jj < num_states_[ 1 ]; ++jj )
		{
			movers_[ 2 ]->initOrientation( *xyz_ );
			for (int kk = 0; kk < num_states_[ 2 ]; ++kk )
			{
				for (int ll = 0; ll < 3; ++ll )
				{
					if ( nothingToScore( ll ) ) continue;
					scores_table_[ ii ][ jj ][ kk ] += xyz_->determineScoreForMover( movers_[ ll ], atoms_to_score_[ ll ], temp );
				}
				//std::cerr << "d3e score: " << ii << " " << jj << " " << kk << " = " << scores_table_[ ii ][ jj ][ kk ] << std::endl;
				movers_[ 2 ]->nextOrientation( *xyz_ );
			}
			movers_[ 1 ]->nextOrientation( *xyz_ );
		}
		movers_[ 0 ]->nextOrientation( *xyz_ );
	}
}

float
DegreeThreeEdge_ths::getScoreForEnabledStates(
	int enabled_state1,
	int enabled_state2,
	int enabled_state3) const
{
	int states[] = {enabled_state1, enabled_state2, enabled_state3};
	int original_states[3];
	for (int ii = 0; ii < 3; ++ii)
	{
		original_states[ ii ] = vertex_ptrs_[ ii ]->convertEnabledStateToRegularState( states[ ii ] );
	}
	return scores_table_[ original_states[ 0 ]][ original_states[ 1 ] ][ original_states[ 2 ]];
}


float DegreeThreeEdge_ths::getScoreForAssignedState() const
{
	return scores_table_[ vertex_ptrs_[ 0 ]->getAssignedState() ][ vertex_ptrs_[ 1 ]->getAssignedState() ][ vertex_ptrs_[ 2 ]->getAssignedState() ];
}

int
DegreeThreeEdge_ths::whichVertex( int vertex ) const
{
	for (int ii = 0; ii < 3; ++ii)
	{
		if ( vertex_indices_[ ii ] == vertex ) return ii;
	}
	std::cerr << "CRITICAL ERROR: did not find sought vertex: " << vertex << " out of [ ";
	std::cerr << vertex_indices_[ 0 ] << ", " << vertex_indices_[ 1 ] << ", " << vertex_indices_[ 2 ];
	std::cerr << "]" << std::endl;
	assert( false );
	exit(1);
}

void
DegreeThreeEdge_ths::compareSharedAtomPairs(
	DegreeThreeEdge_ths const  *  other,
	int vertex_this,
	int vertex_other
) const
{
	//std::cerr << "compareSharedAtomPairs: " << mover_atoms_atom_overlaps_[ vertex_this ].size() <<  " ";
	//std::cerr << other->mover_atoms_atom_overlaps_[ vertex_other ].size() << std::endl;
	//std::cerr << "comparing vertices: " << vertex_indices_[ vertex_this ]  << " and ";
	//std::cerr << other->vertex_indices_[ vertex_other ] << std::endl;
	for (int ii = 0; ii < num_atoms_[ vertex_this ]; ++ii )
	{
	
		std::list< AtomDescr > atomsScoredOnBothEdges;
		set_intersect_apl( mover_atoms_atom_overlaps_[ vertex_this ][ ii ],
			other->mover_atoms_atom_overlaps_[ vertex_other ][ ii ],
			atomsScoredOnBothEdges );
		if ( ! atomsScoredOnBothEdges.empty() )
		{
			//Discovered atoms that would be doubly counted
			vertex_ptrs_[ vertex_this ]->boot( ii );
		}
	}
}

void
DegreeThreeEdge_ths::compareSharedBackgroundAtoms
(
	DegreeThreeEdge_ths const  *  other,
	int vertex_this,
	int vertex_other
) const
{
	for (int ii = 0; ii < num_atoms_[ vertex_this ]; ++ii )
	{
		std::list< AtomDescr > bgatomsScoredOnBothEdges;
		set_intersect_apl( bg_ats_for_atom_[ vertex_this ][ ii ],
			other->bg_ats_for_atom_[ vertex_other ][ ii ],
			bgatomsScoredOnBothEdges );
		if ( ! bgatomsScoredOnBothEdges.empty() )
		{
			//Discovered atom / bg atom pair that would be doubly counted
			vertex_ptrs_[ vertex_this ]->boot( ii );
		}
	}	
}

bool DegreeThreeEdge_ths::nothingToScore() const
{
	for (int ii = 0; ii < 3; ++ii )
	{
		if ( ! atoms_to_score_[ ii ].empty() ) return false;
	}
	return true;
}

bool DegreeThreeEdge_ths::nothingToScore( int vertex ) const
{
	return atoms_to_score_[ vertex ].empty();
}

GraphToHoldScores::GraphToHoldScores( 
	AtomPositions * xyz,
	std::vector< int > const & num_states,
	std::vector< Mover * > const & movers )
:
	xyz_( xyz ),
	num_vertices_( movers.size() ),
	num_deg2edges_( 0 ),
	num_deg3edges_( 0 ),
	vertices_( num_vertices_, static_cast< Vertex_ths * > ( 0 ) )
{
	//std::cerr << "allocateConnectivityTable();" << std::endl;
	allocateConnectivityTable();	
	//std::cerr << "instantiateVertices( num_states, movers );" << std::endl;
	instantiateVertices( num_states, movers );
	//std::cerr << "obtainAtomsFromMovers();" << std::endl;
	obtainAtomsFromMovers();
	//std::cerr << "getOverlappingBackgroundAtoms();" << std::endl;
	getOverlappingBackgroundAtoms();
	//std::cerr << "dropMoverAtomsFromBackgroundAtomLists();" << std::endl;
	dropMoverAtomsFromBackgroundAtomLists();
	//std::cerr << "addDegreeTwoEdges();" << std::endl;
	addDegreeTwoEdges();
	//std::cerr << "promoteBackgroundAtomsToDegreeTwoEdges();" << std::endl;
	promoteBackgroundAtomsToDegreeTwoEdges();
	//std::cerr << "addDegreeThreeEdges();" << std::endl;
	addDegreeThreeEdges();
	///std::cerr << "detectFourWayInteractions();" << std::endl;
	detectFourWayInteractions();
	//std::cerr << "cascadeScoringInfo();" << std::endl;
	cascadeScoringInfo();
	//std::cerr << "finalizeScoreVectors();" << std::endl;
	finalizeScoreVectors();
	if (xyz_->outputNotice() ) std::cerr << "score hypergraph" << std::endl;
	score();
};

GraphToHoldScores::~GraphToHoldScores()
{
	if ( connectivity_ != 0 )
	{
		for (int ii = 0; ii < num_vertices_; ++ii )
		{
			delete [] connectivity_[ ii ];
		}
		delete [] connectivity_; connectivity_ = 0;
	}
	
	std::list< DegreeTwoEdge_ths * >::iterator iter2 = deg2edges_.begin();
	while ( iter2 != deg2edges_.end() )
	{
		std::list< DegreeTwoEdge_ths * >::iterator nextiter2 = iter2;
		++nextiter2;
		delete (*iter2);
		iter2 = nextiter2;
	}

	std::list< DegreeThreeEdge_ths * >::iterator iter3 = deg3edges_.begin();
	while ( iter3 != deg3edges_.end() )
	{
		std::list< DegreeThreeEdge_ths  * >::iterator nextiter3 = iter3;
		++nextiter3;
		delete (*iter3);
		iter3 = nextiter3;
	}
	for (int ii = 0; ii < num_vertices_; ++ii )
	{
		delete vertices_[ ii ];
	}
}

int
GraphToHoldScores::getNumNodes() const
{
	return num_vertices_;
}

Vertex_ths *
GraphToHoldScores::getVertexPtr( int vertex ) const
{
	nodeInRange( vertex );
	return vertices_[ vertex ];
}

DegreeTwoEdge_ths*
GraphToHoldScores::getDegree2Edge( int fn, int sn)
{
	nodeInRange( fn );
	nodeInRange( sn );
	assert( fn < sn );
	DegreeTwoEdge_ths * edge = findDeg2Edge( fn, sn );
	if (edge == 0 ) edge = addDegree2Edge( fn, sn );
	return edge;
}

DegreeThreeEdge_ths*
GraphToHoldScores::getDegree3Edge( int fn, int sn, int tn )
{
	nodeInRange( fn );
	nodeInRange( sn );
	nodeInRange( tn );
	assert( fn < sn && sn < tn );
	DegreeThreeEdge_ths * edge = findDeg3Edge(fn, sn, tn );
	if (edge == 0 ) edge = addDegree3Edge( fn, sn, tn );
	return edge;
}

void GraphToHoldScores::forceClique( std::list< int > const & movers_in_clique )
{
	std::list< int >::const_iterator iter_outer = movers_in_clique.begin();
	while( iter_outer != movers_in_clique.end() )
	{
		std::list< int >::const_iterator iter_inner = iter_outer;
		++iter_inner;
		while (iter_inner != movers_in_clique.end() )
		{
			assert( *iter_outer < *iter_inner );
			if ( ! connectivity_[ *iter_outer ][ *iter_inner ] )
				addDegree2Edge( *iter_outer, *iter_inner );
			++iter_inner;
		}
		++iter_outer;
	}
}

void GraphToHoldScores::dropEdge( std::list< DegreeTwoEdge_ths * >::iterator iter )
{
	deg2edges_.erase( iter );
	--num_deg2edges_;
}

void GraphToHoldScores::dropEdge( std::list< DegreeThreeEdge_ths * >::iterator iter )
{
	deg3edges_.erase( iter );
	--num_deg3edges_;
}

void GraphToHoldScores::getPenalties(
	std::vector< std::vector< float > > & penalties, 
	std::list< float > & allPenalties )
{
	for ( int ii = 0; ii < num_vertices_; ++ii)
	{
		vertices_[ ii ]->getPenalties( penalties[ ii ] );
		for ( int jj = 0; jj < penalties[ ii ].size(); ++jj)
		{
			allPenalties.push_back( penalties[ ii ][ jj ] );
		}
	}
	allPenalties.sort();
	allPenalties.unique();
}


//Methods for state disabling
void GraphToHoldScores::setAllStatesEnabled()
{
	for (int ii = 0; ii < getNumNodes(); ++ii )
	{
		vertices_[ ii ]->enableAllStates(); 
	}
}

void GraphToHoldScores::disableStateOnNode( int vertex, int state )
{
	nodeInRange( vertex );
	vertices_[ vertex ]->disableState( state );
}

void GraphToHoldScores::setStateDisablingCompleteForNow()
{
	for (int ii = 0; ii < getNumNodes(); ++ii )
	{
		vertices_[ ii ]->setStateDisablingComplete(); 
	}
}

int 
GraphToHoldScores::convertEnabledState2OriginalStateEnumerationOnNode( 
	int vertex, 
	int enabled_state ) const
{
	nodeInRange( vertex );
	return vertices_[ vertex ]->convertEnabledStateToRegularState( enabled_state );
}

//Methods for initializing interaction graph
AtomPositions * GraphToHoldScores::getAtomPositionsPointer() const
{
	return xyz_;
}

int GraphToHoldScores::getNumStatesForNode(int vertex_index ) const
{
	nodeInRange( vertex_index );
	return vertices_[ vertex_index ]->getNumStates();
}

int GraphToHoldScores::getNumEnabledStatesForNode( int vertex_index ) const
{
	nodeInRange( vertex_index );
	return vertices_[ vertex_index ]->getNumEnabledStates();
}

float GraphToHoldScores::getNodeScoreForState( int vertex_id, int state) const
{
	nodeInRange( vertex_id );
	return vertices_[ vertex_id ]->getScoreForEnabledState( state );
}

bool GraphToHoldScores::getNodeHasAnyHighOrderOverlap( int vertex_id ) const
{
	nodeInRange( vertex_id );
	return vertices_[ vertex_id ]->getHasAnyHighOrderOverlap();
}

std::list< AtomDescr > GraphToHoldScores::getAtomsInHighOrderOverlapForNode( int vertex_id ) const
{
	nodeInRange( vertex_id );
	return vertices_[ vertex_id ]->getAtomsInHighOrderOverlap();
}

Mover* GraphToHoldScores::getMoverForNode( int vertex_id ) const
{
	nodeInRange( vertex_id );
	return vertices_[ vertex_id ]->getMover();
}

void GraphToHoldScores::setD2EIteratorAtBegining()
{
	d2eiter_ = deg2edges_.begin();
}

void GraphToHoldScores::incrementD2EIterator()
{
	++d2eiter_;
}

bool GraphToHoldScores::getD2EIteratorAtEnd() const
{
	return d2eiter_ == deg2edges_.end();
}

int GraphToHoldScores::getFirstIndexFocusedD2E() const
{
	return (*d2eiter_)->getFirstNodeIndex();
}

int GraphToHoldScores::getSecondIndexFocusedD2E() const
{
	return (*d2eiter_)->getSecondNodeIndex();
}

float GraphToHoldScores::getScoreForFocusedD2E( int state1, int state2 ) const
{
	return (*d2eiter_)->getScoreForEnabledStates( state1, state2 );
}

void GraphToHoldScores::setD3EIteratorAtBegining()
{
	d3eiter_ = deg3edges_.begin();
}

void GraphToHoldScores::incrementD3EIterator()
{
	++d3eiter_;
}

bool GraphToHoldScores::getD3EIteratorAtEnd() const
{
	return d3eiter_ == deg3edges_.end();
}

int GraphToHoldScores::getFirstIndexFocusedD3E() const
{
	return (*d3eiter_)->getFirstNodeIndex();
}

int GraphToHoldScores::getSecondIndexFocusedD3E() const
{
	return (*d3eiter_)->getSecondNodeIndex();
}

int GraphToHoldScores::getThirdIndexFocusedD3E() const
{
	return (*d3eiter_)->getThirdNodeIndex();
}

float
GraphToHoldScores::getScoreForFocusedD3E( int state1, int state2, int state3) const
{
	return (*d3eiter_)->getScoreForEnabledStates( state1, state2, state3 );
}


float
GraphToHoldScores::getScoreForStateAssignment(
	std::vector< int > const & states
)
{
	for (int ii = 0; ii < num_vertices_; ++ii)
	{
		vertices_[ ii ]->assignState( states[ ii ] );
	}
	
	float score = 0;
	for (int ii = 0; ii < num_vertices_; ++ii )
	{
		score += vertices_[ ii ]->getScoreForAssignedState();
	}
	for (std::list< DegreeTwoEdge_ths* >::const_iterator iter = deg2edges_.begin();
		iter != deg2edges_.end(); ++iter)
	{
		score += (*iter)->getScoreForAssignedState();
	}
	for (std::list< DegreeThreeEdge_ths* >::const_iterator iter = deg3edges_.begin();
		iter != deg3edges_.end(); ++iter)
	{
		score += (*iter)->getScoreForAssignedState();
	}
	return score;
}

//Private methods
void
GraphToHoldScores::allocateConnectivityTable()
{
	connectivity_ = new bool* [ num_vertices_];
	for (int ii = 0; ii < num_vertices_; ++ii)
	{
		connectivity_[ ii ] = new bool [ num_vertices_ ];
		for (int jj = 0; jj < num_vertices_; ++jj )
		{
			connectivity_[ ii ][ jj ] = false;
		}
	}
}

void
GraphToHoldScores::instantiateVertices
( 
	std::vector< int > const & num_states,
	std::vector< Mover * > const & movers
)
{
	for (int ii = 0; ii < num_vertices_; ++ii )
	{
		vertices_[ ii ] = new Vertex_ths( this, xyz_, ii, num_states[ ii ]);
		vertices_[ ii ]->setMover( movers[ ii ] );
	}
}

void
GraphToHoldScores::obtainAtomsFromMovers()
{
	for (int ii = 0; ii < getNumNodes(); ++ii)
	{
		vertices_[ ii ]->obtainAtomsFromMover();
	}
}

void
GraphToHoldScores::getOverlappingBackgroundAtoms()
{
	for (int ii = 0; ii < getNumNodes(); ++ii)
	{
		vertices_[ ii ]->getOverlappingBackgroundAtoms();
	}
}

void
GraphToHoldScores::dropMoverAtomsFromBackgroundAtomLists()
{
	for (int ii = 0; ii < getNumNodes(); ++ii)
	{
		for (int jj = ii+1; jj < getNumNodes(); ++jj)
		{
			vertices_[ ii ]->dropMoverAtomsFromBackgroundAtomLists( vertices_[ jj ] );
		}
	}
}

void
GraphToHoldScores::addDegreeTwoEdges()
{
	for (int ii = 0; ii < getNumNodes(); ++ii)
	{
		for (int jj = ii+1; jj < getNumNodes(); ++jj)
		{
			if (vertices_[ ii ]->anyMoverOverlap( vertices_[ jj ] ) )
			{
				addDegree2Edge( ii, jj );
			}
		}
	}	
}

void
GraphToHoldScores::promoteBackgroundAtomsToDegreeTwoEdges()
{
	for ( std::list< DegreeTwoEdge_ths * >::iterator iter = deg2edges_.begin();
		iter != deg2edges_.end(); ++iter )
	{
		(*iter)->detectBackgroundAtomsToScoreForEdge();
	}
}

void
GraphToHoldScores::addDegreeThreeEdges()
{
	for (int ii = 0; ii < getNumNodes(); ++ii)
	{
		vertices_[ ii ]->detectBackgroundAtomsInThreeWayOverlap();
	}

	if ( getNumNodes() < 3 ) return;
	for (int ii = 0; ii < getNumNodes(); ++ii)
	{
		for (int jj = ii+1; jj < getNumNodes(); ++jj)
		{
			if ( ! connectivity_[ ii ][ jj ] ) continue;
			for (int kk = jj+1; kk < getNumNodes(); ++kk)
			{
				if ( ! connectivity_[ ii ][ kk ] ) continue;
				if ( ! connectivity_[ jj ][ kk ] ) continue;
				
				vertices_[ ii ]->detectThreeWayOverlap( vertices_[ jj ], vertices_[ kk ]);
			}
		}
	}
}

void
GraphToHoldScores::detectFourWayInteractions()
{
	//compare all degree3 hyperedges
	std::list< DegreeThreeEdge_ths * >::iterator iter_outer = deg3edges_.begin();
	while ( iter_outer != deg3edges_.end() )
	{
		std::list< DegreeThreeEdge_ths * >::iterator iter_inner = iter_outer;
		++iter_inner;
		while ( iter_inner != deg3edges_.end() )
		{
			(*iter_outer)->detect_shared_atom_pair( *iter_inner );
			(*iter_outer)->detectBackgroundAtomScoredTwiceForAtom( *iter_inner );
			++iter_inner;
		}
		++iter_outer;
	}
}

void
GraphToHoldScores::cascadeScoringInfo()
{
	for ( int ii = 0; ii < getNumNodes(); ++ii )
	{
		vertices_[ ii ]->informIncidentEdgesAboutBootedAtoms();
	}
	
	for( std::list< DegreeTwoEdge_ths * >::iterator iter = deg2edges_.begin();
		iter != deg2edges_.end(); ++iter )
	{
		(*iter)->tellVerticesToNotScoreBGAtomsScoredOnD2Edge();
	}
	
	for( std::list< DegreeThreeEdge_ths * >::iterator iter = deg3edges_.begin();
		iter != deg3edges_.end(); ++iter )
	{
		(*iter)->tellD2EdgesToNotScoreBGandMoverAtomsScoredOnD3Edge();
	}

}

void
GraphToHoldScores::finalizeScoreVectors()
{
	//std::cerr << "Entering finalize score vectors" << std::endl;
	for ( int ii = 0; ii < getNumNodes(); ++ii )
	{
		vertices_[ ii ]->finalizeScoreVectors();
	}

	for( std::list< DegreeTwoEdge_ths * >::iterator iter = deg2edges_.begin();
		iter != deg2edges_.end(); ++iter )
	{
		(*iter)->finalizeScoreVectors();
	}

	for( std::list< DegreeThreeEdge_ths * >::iterator iter = deg3edges_.begin();
		iter != deg3edges_.end(); ++iter )
	{
		(*iter)->finalizeScoreVectors();
	}
	
}

void
GraphToHoldScores::score()
{
	for ( int ii = 0; ii < getNumNodes(); ++ii )
	{
		vertices_[ ii ]->score();
	}
	
	for( std::list< DegreeTwoEdge_ths * >::iterator iter = deg2edges_.begin();
		iter != deg2edges_.end(); ++iter )
	{
		(*iter)->score();
	}
	
	for( std::list< DegreeThreeEdge_ths * >::iterator iter = deg3edges_.begin();
		iter != deg3edges_.end(); ++iter )
	{
		(*iter)->score();
	}
}

void
GraphToHoldScores::nodeInRange( int node ) const
{	assert( 0 <= node && node < getNumNodes() ); }

DegreeTwoEdge_ths*
GraphToHoldScores::findDeg2Edge( int fn, int sn ) const
{
	for (std::list< DegreeTwoEdge_ths * >::const_iterator iter = deg2edges_.begin();
		iter != deg2edges_.end(); ++iter )
	{
		if ( (*iter)->sameEdge( fn, sn ) )
			return *iter;
	}
	return 0;
}

DegreeTwoEdge_ths*
GraphToHoldScores::addDegree2Edge( int node1, int node2 )
{
	nodeInRange( node1 );
	nodeInRange( node2 );
	
	DegreeTwoEdge_ths * new_edge = new DegreeTwoEdge_ths( this, xyz_, node1, node2 );
	deg2edges_.push_front( new_edge );
	new_edge->setPosInOwnersEdgeList( deg2edges_.begin() );
	
	connectivity_[ node1 ][ node2 ] = true;
	connectivity_[ node2 ][ node1 ] = true;
	++num_deg2edges_;
	
	return new_edge;	
}

DegreeThreeEdge_ths*
GraphToHoldScores::findDeg3Edge( int fn, int sn, int tn ) const
{
	for (std::list< DegreeThreeEdge_ths * >::const_iterator iter = deg3edges_.begin();
		iter != deg3edges_.end(); ++iter )
	{
		if ( (*iter)->sameEdge( fn, sn, tn ) )
			return *iter;
	}
	return 0;
}

DegreeThreeEdge_ths*
GraphToHoldScores::addDegree3Edge( int fn, int sn, int tn )
{
	//std::cerr << "Adding degree 3 edge: " << fn << ", " << sn << ", " << tn << std::endl;
	DegreeThreeEdge_ths * new_edge = new DegreeThreeEdge_ths( this, xyz_, fn, sn, tn );
	deg3edges_.push_front( new_edge );
	new_edge->setPosInOwnersEdgeList( deg3edges_.begin() );
	++num_deg3edges_;
	return new_edge;
}

