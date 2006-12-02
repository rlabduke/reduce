#include "GraphToHoldScores.h"
#include "AtomPositions.h"
#include "DotSph.h"
#include "Mover.h"
#include "PDBrec.h"
#include "Point3d.h"

#ifdef OLD_STD_HDRS
#include <assert.h>
#else
#include <cassert>
using std::exit;
#endif

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


void sort_four( int a, int b, int c, int d, int & w, int & x, int & y, int & z)
{
	int four_ints[ 4 ];
	int three_ints[ 3 ];
	sort_three( a, b, c, three_ints[ 0 ], three_ints[ 1 ], three_ints[ 2 ]);
	int offset = 0;
	for (int ii = 0; ii < 3; ++ii)
	{
		if ( offset == 0 && d < three_ints[ ii ] )
		{
			four_ints[ ii ] = d;
			offset = 1;
		}
		four_ints[ ii + offset ] = three_ints[ ii ];
	}
	w = four_ints[ 0 ];
	x = four_ints[ 1 ];
	y = four_ints[ 2 ];
	z = four_ints[ 3 ];
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

template < class T >
void set_merge_apl(
	std::list< T > & sortedList1,
	std::list< T > const & sortedList2 
)
{
	//std::cerr << "Begin set merge" << std::endl;
	typename std::list< T >::iterator list1iter = sortedList1.begin();
	typename std::list< T >::const_iterator list2iter = sortedList2.begin();
	while( list1iter != sortedList1.end() && list2iter != sortedList2.end()  )
	{
		//std::cerr << "set merge: iter1: " << *list1iter << " iter2: " << *list2iter << std::endl;

		if ( *list1iter == *list2iter )
		{
			//std::cerr << "merge equal: " << *list1iter << std::endl;			
			++list1iter;
			++list2iter;
		}
		else if ( *list1iter < *list2iter )
		{
			++list1iter;
		}
		else {
			//std::cerr << "merge insert: " << *list2iter << std::endl;			
			sortedList1.insert( list1iter, *list2iter );
			++list2iter;
		}
	}
	while ( list2iter != sortedList2.end() )
	{
		//std::cerr << "merge tail insert: " << *list2iter << std::endl;
		sortedList1.insert( list1iter, *list2iter );
		++list2iter;
	}
	//std::cerr << "end set merge" << std::endl;

}

DotsForAtom::DotsForAtom() : 
	theAtom_(),
	dotSphereManager_( 0 ),
	dotsOn_( 0 ),
	indsOfDotsOn_( 0 ),
	indsCurrent_( false ),
	numDots_( 0 ),
	numDotsOn_( 0 )
{}

DotsForAtom::~DotsForAtom()
{
	delete [] dotsOn_;
	delete [] indsOfDotsOn_;
}

DotsForAtom::DotsForAtom( DotsForAtom const & rhs ) :
	theAtom_( rhs.theAtom_ ),
	dotSphereManager_( rhs.dotSphereManager_  ),
	dotsOn_( new bool[ rhs.numDots_ ]  ),
	indsOfDotsOn_( new int[ rhs.numDots_ ] ),
	indsCurrent_( rhs.indsCurrent_ ),
	numDots_( rhs.numDots_  ),
	numDotsOn_( rhs.numDotsOn_ )
{
	for (int ii = 0; ii < numDots_; ++ii)
	{
		dotsOn_[ ii ] = rhs.dotsOn_[ ii ];
		indsOfDotsOn_[ ii ] = rhs.indsOfDotsOn_[ ii ];
	}
}

DotsForAtom const & 
DotsForAtom::operator = ( DotsForAtom const & rhs )
{
	theAtom_ = rhs.theAtom_;
	dotSphereManager_ = rhs.dotSphereManager_;
	delete [] dotsOn_; dotsOn_ = new bool[ rhs.numDots_ ];
	delete [] indsOfDotsOn_; indsOfDotsOn_ = new int[ rhs.numDots_ ];
	indsCurrent_ = rhs.indsCurrent_;
	numDots_ = rhs.numDots_;
	numDotsOn_ = rhs.numDotsOn_;

	for (int ii = 0; ii < numDots_; ++ii)
	{
		dotsOn_[ ii ] = rhs.dotsOn_[ ii ];
		indsOfDotsOn_[ ii ] = rhs.indsOfDotsOn_[ ii ];
	}
	return *this;
}
	
void DotsForAtom::setDotManager( DotSphManager * dotSphereManager )
{
	assert( dotSphereManager_ == 0 );
	dotSphereManager_ = dotSphereManager;
}

void DotsForAtom::setAtom( AtomDescr theAtom )
{
	theAtom_ = theAtom;
	DotSph const & dotSphere = dotSphereManager_->fetch( theAtom_.getAtRadius() );
	numDots_ = dotSphere.count();
	delete [] dotsOn_; dotsOn_ = new bool[ numDots_ ];
	delete [] indsOfDotsOn_; indsOfDotsOn_ = new int[ numDots_ ];
	for (int ii = 0; ii < numDots_; ++ii )
	{
		dotsOn_[ ii ] = false;
		indsOfDotsOn_[ ii ] = -1;
	}
	indsCurrent_ = true;
	numDotsOn_ = 0;
}

void DotsForAtom::turnOnAllDots()
{
	for (int ii = 0; ii < numDots_; ++ii)
	{
		dotsOn_[ ii ] = true;
		indsOfDotsOn_[ ii ] = ii;
	}
	indsCurrent_ = true;
	numDotsOn_ = numDots_;
}

void DotsForAtom::turnOffAllDots()
{
	if ( numDotsOn_ == 0 ) return;
	indsCurrent_ = true;
	if (numDotsOn_ == 0 ) return;
	for (int ii = 0; ii < numDots_; ++ii)
	{
		dotsOn_[ ii ] = false;
		indsOfDotsOn_[ ii ] = -1;
	}
	numDotsOn_ = 0;
}

void DotsForAtom::findContained( std::list< AtomDescr > const & atoms )
{
	//std::cerr << "Entering findContained" << std::endl;
	turnOffAllDots();
	DotSph const & dotSphere = dotSphereManager_->fetch( theAtom_.getAtRadius() );
	for (int ii = 0; ii < numDots_; ++ii)
	{
		Point3d const p = theAtom_.getAtomPos() + dotSphere[ ii ];
		for (std::list< AtomDescr >::const_iterator iter = atoms.begin(); 
			iter != atoms.end(); ++iter )
		{
			if ( iter->containsPoint( p ) )
			{
				dotsOn_[ ii ] = true;
				indsOfDotsOn_[ numDotsOn_ ] = ii;
				++numDotsOn_;
				break;
			}
		}
	}
	//std::cerr << "Leaving findContained" << std::endl;
}


bool DotsForAtom::dotOn( int dot ) const
{
	return dotsOn_[ dot ];
}

int DotsForAtom::numDotsOn() const
{
	return numDotsOn_;
}

bool DotsForAtom::anyDotsShared( DotsForAtom const & other ) const
{
	assert( numDots_ == other.numDots_ );
	if ( numDotsOn_ == 0 || other.numDotsOn_ == 0 ) return false;

	updateIndsOfDotsOn();
	other.updateIndsOfDotsOn();
		
	int iter_this = 0;
	int iter_other = 0;
	while ( iter_this != numDotsOn_ && iter_other != other.numDotsOn_ )
	{
		//assert( 0 <= indsOfDotsOn_[ iter_this ] && indsOfDotsOn_[ iter_this ] < numDots_ );
		//assert( 0 <= other.indsOfDotsOn_[ iter_other ] && other.indsOfDotsOn_[ iter_other ] < numDots_ );
		if ( indsOfDotsOn_[ iter_this ] == other.indsOfDotsOn_[ iter_other ])
		{
			return true;
		}
		else if (indsOfDotsOn_[ iter_this ] < other.indsOfDotsOn_[ iter_other ])
		{
			++iter_this;
		}
		else
		{
			++iter_other;
		}
	}
	return false;
}

void DotsForAtom::intersect( DotsForAtom const & other, DotsForAtom & intersection ) const
{
	assert( numDots_ == other.numDots_ && numDots_ == intersection.numDots_ );
	
	if ( numDotsOn_ == 0 || other.numDotsOn_ == 0 ) return;

	updateIndsOfDotsOn();
	other.updateIndsOfDotsOn();
	intersection.turnOffAllDots();
		
	int iter_this = 0;
	int iter_other = 0;
	while ( iter_this != numDotsOn_ && iter_other != other.numDotsOn_ )
	{
		//assert( 0 <= indsOfDotsOn_[ iter_this ] && indsOfDotsOn_[ iter_this ] < numDots_ );
		//assert( 0 <= other.indsOfDotsOn_[ iter_other ] && other.indsOfDotsOn_[ iter_other ] < numDots_ );

		if ( indsOfDotsOn_[ iter_this ] == other.indsOfDotsOn_[ iter_other ])
		{
			intersection.dotsOn_[ indsOfDotsOn_[ iter_this] ] = true;
			intersection.indsOfDotsOn_[ intersection.numDotsOn_ ] = indsOfDotsOn_[ iter_this ];
			++intersection.numDotsOn_;
			++iter_this;
			++iter_other;
		}
		else if (indsOfDotsOn_[ iter_this ] < other.indsOfDotsOn_[ iter_other ])
		{
			++iter_this;
		}
		else
		{
			++iter_other;
		}
	}
}

void DotsForAtom::subtractDotsOn( DotsForAtom const & other )
{
	assert( numDots_ == other.numDots_ );
	other.updateIndsOfDotsOn();
	for (int ii = 0; ii < other.numDotsOn_; ++ii)
	{
		if ( dotsOn_[ other.indsOfDotsOn_[ ii] ] )
		{
			dotsOn_[ other.indsOfDotsOn_[ ii] ] = false;
			--numDotsOn_;
			indsCurrent_ = false;
		}
	}
}

void DotsForAtom::addDots( DotsForAtom const & other )
{
	assert( numDots_ == other.numDots_ );
	other.updateIndsOfDotsOn();
	for (int ii = 0; ii < other.numDotsOn_; ++ii)
	{
		if ( ! dotsOn_[ other.indsOfDotsOn_[ ii] ] )
		{
			dotsOn_[ other.indsOfDotsOn_[ ii] ] = true;
			++numDotsOn_;
			indsCurrent_ = false;
		}
	}
}

std::ostream & 
operator << (std::ostream & os, DotsForAtom const & dots )
{
	dots.print( os );
	return os;
}

void
DotsForAtom::print( std::ostream & os ) const
{
	updateIndsOfDotsOn();
	os << "Dots: " << std::endl;
	for (int ii = 0; ii < numDotsOn_; ++ii )
	{
		os << indsOfDotsOn_[ ii ] << " ";
		if ( ii % 20 == 19 )
		{
			os << std::endl;
		}
	}
}

void DotsForAtom::updateIndsOfDotsOn() const
{
	if ( indsCurrent_ ) return;
	int countOnDotsSeen = 0;
	for (int ii = 0; ii < numDots_; ++ii)
	{
		if ( dotsOn_[ ii ] )
		{
			indsOfDotsOn_[ countOnDotsSeen ] = ii;
			++countOnDotsSeen;
		}
	}
	indsCurrent_ = true;
}

//----------------------------------------------------------------------------//
//------------------------------ Vertex Class --------------------------------//
//----------------------------------------------------------------------------//


Vertex_ths::Vertex_ths(
	GraphToHoldScores * owner,
	AtomPositions * xyz,
	DotSphManager * dotSphereManager,
	int index,
	int num_states 
)
:
	owner_( owner ),
	xyz_( xyz ),
	dotSphereManager_( dotSphereManager ),
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
  	//std::cerr << "Mover " << index_ << " with " << atoms_.size() << " atoms" << std::endl;
  	std::copy( atoms.begin(), atoms.end(), atoms_.begin() );
	dotsForAtom_.resize( atoms_.size() );
	for (int ii = 0; ii < atoms_.size(); ++ii)
	{
		dotsForAtom_[ ii ].setDotManager( dotSphereManager_ );
		dotsForAtom_[ ii ].setAtom( atoms_[ ii ] );
		dotsForAtom_[ ii ].turnOnAllDots();
	}
	
	atom_in_high_order_overlap_.resize( getNumAtoms() );
	std::fill( atom_in_high_order_overlap_.begin(), atom_in_high_order_overlap_.end(), false );
	dotsInHOO_.resize( getNumAtoms() );
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
Vertex_ths::detectThreeWayOverlap( Vertex_ths * other1, Vertex_ths * other2 )
{
	assert( index_ < other1->index_ && other1->index_ < other2->index_ );
	//std::cerr << "Entering Vertex_ths::detectThreeWayOverlap with: " << index_;
	//std::cerr << " " << other1->index_;
	//std::cerr << " " << other2->index_ << std::endl;

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
					
					//std::cerr << "3-way mover overlap: ";
					//std::cerr << index_ << " ";
					//std::cerr << other1->index_ << " ";
					//std::cerr << other2->index_ << " ";
					//std::cerr << atoms_[ ii ] << " ";
					//std::cerr << other1->atoms_[ jj ] << " ";
					//std::cerr << other2->atoms_[ kk ] << std::endl;
					
					DegreeThreeEdge_ths * d3e = owner_->getDegree3Edge( index_, other1->index_, other2->index_ );
					d3e->setMoverAtomHas3Way( index_, ii );
					d3e->setMoverAtomHas3Way( other1->index_, jj );
					d3e->setMoverAtomHas3Way( other2->index_, kk );

				}
			}
		}
	}
	//std::cerr << "Leaving Vertex_ths::detectThreeWayOverlap with: " << index_;
	//std::cerr << " " << other1->index_;
	//std::cerr << " " << other2->index_ << std::endl;
		
}

//atom in high order overlap -- either with other movers or with a background atom 
void
Vertex_ths::boot(
	int atom,
	DotsForAtom const & dotsToBoot,
	std::list< int > const & moversInHOO
)
{
	if (! any_high_order_overlap_ && xyz_->outputNotice() )
	{
		std::cerr << "Detected high order overlap for vertex " << index_ << std::endl;	
	}
	any_high_order_overlap_ = true;
	
	//return; //short circuit for debugging / testing how neccessary high order overlap is.
	
	atom_in_high_order_overlap_[ atom ] = true;
	dotsForAtom_[ atom ].addDots( dotsToBoot );
	owner_->forceClique( moversInHOO );
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

std::list< DegreeFourEdge_ths * >::iterator 
Vertex_ths::addEdge( DegreeFourEdge_ths * edge )
{
	deg4edges_.push_front( edge );
	return deg4edges_.begin();
}


void Vertex_ths::dropEdge( std::list< DegreeTwoEdge_ths * >::iterator iter )
{
	deg2edges_.erase( iter );
}

void Vertex_ths::dropEdge( std::list< DegreeThreeEdge_ths * >::iterator iter )
{
	deg3edges_.erase( iter );
}

void Vertex_ths::dropEdge( std::list< DegreeFourEdge_ths * >::iterator iter )
{
	deg4edges_.erase( iter );
}

void Vertex_ths::informIncidentEdgesAboutBootedAtoms()
{
	if (! any_high_order_overlap_ ) return;
	
	for (int ii = 0; ii < getNumAtoms(); ++ii)
	{
		if ( atom_in_high_order_overlap_[ ii ] )
		{
			dotsForAtom_[ ii ].subtractDotsOn( dotsInHOO_[ ii ] );
		}
	}
	
	for (std::list< DegreeTwoEdge_ths * >::iterator iter = deg2edges_.begin(); 
		iter != deg2edges_.end(); ++iter)
	{
		(*iter)->noteBootedAtoms( index_, atom_in_high_order_overlap_, dotsInHOO_);
	}
	for (std::list< DegreeThreeEdge_ths * >::iterator iter = deg3edges_.begin(); 
		iter != deg3edges_.end(); ++iter)
	{
		(*iter)->noteBootedAtoms( index_, atom_in_high_order_overlap_, dotsInHOO_ );
	}
	for (std::list< DegreeFourEdge_ths * >::iterator iter = deg4edges_.begin(); 
		iter != deg4edges_.end(); ++iter)
	{
		(*iter)->noteBootedAtoms( index_, atom_in_high_order_overlap_, dotsInHOO_ );
	}
}

void Vertex_ths::removeDotsScoredOnD2Edge(
	int atom,
	DotsForAtom const & dots_scored_on_edge )
{
	dotsForAtom_[ atom ].subtractDotsOn( dots_scored_on_edge );
}

void Vertex_ths::finalizeScoreVectors()
{
	//std::cerr << "entering vertex_ths::finalizeScoreVectors for " << index_ << " with # atoms: " << getNumAtoms() << std::endl;

	int count_atoms_to_score = 0;
	std::vector< bool > score_atom( getNumAtoms(), false );
	for (int ii = 0; ii < getNumAtoms(); ++ii )
	{
		if ( dotsForAtom_[ ii ].numDotsOn() != 0  )
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
	dotsToScoreOnAtoms_.resize( count_atoms_to_score );
	int count = 0;
	for (int ii = 0; ii < getNumAtoms(); ++ii)
	{
		if ( ! score_atom[ ii ] )
		{
			continue;
		}
				
		dotsToScoreOnAtoms_[ count ].first = atoms_[ ii ];
		dotsToScoreOnAtoms_[ count ].second = & dotsForAtom_[ ii ];
		
		//std::cerr << "Score: " << atoms_[ ii ] << " " << dotsForAtom_[ ii ] << std::endl;
		++count;
	}
	
	if ( ! any_high_order_overlap_ ) return;
	int numAtomsInHOO = 0;
	for (int ii = 0; ii < getNumAtoms(); ++ii )
	{
		if ( atom_in_high_order_overlap_[ ii ] )
		{
			++numAtomsInHOO;
		}
	}
	dotsToScoreInHOO_.resize( numAtomsInHOO );
	numAtomsInHOO = 0;
	for (int ii = 0; ii < getNumAtoms(); ++ii )
	{
		if ( atom_in_high_order_overlap_[ ii ] )
		{
			dotsToScoreInHOO_[ numAtomsInHOO ].first = atoms_[ ii ];
			dotsToScoreInHOO_[ numAtomsInHOO ].second = & ( dotsInHOO_[ ii ] );
			++numAtomsInHOO;
		}
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
		scores_[ ii ] = xyz_->determineScoreForMover( mover_, dotsToScoreOnAtoms_, penalty );
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
}

std::vector< std::pair< AtomDescr, DotsForAtom * > >
Vertex_ths::getAtomsInHighOrderOverlap() const
{
	return dotsToScoreInHOO_;
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

//----------------------------------------------------------------------------//
//----------------------------- Degree Two Edge ------------------------------//
//----------------------------------------------------------------------------//


DegreeTwoEdge_ths::DegreeTwoEdge_ths(
	GraphToHoldScores * owner,
	AtomPositions * xyz,
	DotSphManager * dotSphereManager,
	int vertex1,
	int vertex2
)
:
	owner_( owner ),
	xyz_( xyz ),
	dotSphereManager_( dotSphereManager )
{
	vertex_indices_[ 0 ] = vertex1;
	vertex_indices_[ 1 ] = vertex2;
	//std::cerr << "Creating degree-2 edge: " << vertex1 << " " << vertex2 << std::endl;
	for (int ii = 0; ii < 2; ++ii)
	{
		vertex_ptrs_[ ii ] = owner_->getVertexPtr( vertex_indices_[ ii ] );
		pos_in_verts_edge_list_[ ii ] = vertex_ptrs_[ ii ]->addEdge( this );
		num_states_[ ii ] = vertex_ptrs_[ ii ]->getNumStates();
		num_atoms_[ ii ] = vertex_ptrs_[ ii ]->getNumAtoms();
		movers_[ ii ] = vertex_ptrs_[ ii ]->getMover();
		mover_atoms_atom_overlaps_[ ii ].resize( num_atoms_[ ii ] );
		atoms_[ ii ].resize( num_atoms_[ ii ] );
		vertex_ptrs_[ ii ]->getAtoms( atoms_[ ii ] );
		dotsForAtom_[ ii ].resize( atoms_[ ii ].size() );
		//std::cerr << ii << " with " << atoms_[ ii ].size() << " atoms " << std::endl;
		for (int jj = 0; jj < atoms_[ ii ].size(); ++jj)
		{
			dotsForAtom_[ ii ][ jj ].setDotManager( dotSphereManager_ );
			dotsForAtom_[ ii ][ jj ].setAtom( atoms_[ ii ][ jj ] );
		}

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
	return (vertex_indices_[ 0 ] == vertex_index || vertex_indices_[ 1 ] == vertex_index );
}

void
DegreeTwoEdge_ths::detectDotsToScoreOnEdge()
{
	//std::cerr << "detectDotsToScoreOnEdge(); " << vertex_indices_[ 0 ] << " " << vertex_indices_[ 1 ] << std::endl;
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
			}
		}
	}
	
	for (int ii = 0; ii < 2; ++ii )
	{
		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj )
		{
			//std::cerr << "FindContained: " << ii << " " << jj << std::endl;
			dotsForAtom_[ ii ][ jj ].findContained( mover_atoms_atom_overlaps_[ ii ][ jj ] );
		}
	}
}

DotsForAtom const & 
DegreeTwoEdge_ths::getDotsForAtom( int vertex, int atom )
{
	return dotsForAtom_[ whichVertex( vertex ) ][ atom ];
}

void
DegreeTwoEdge_ths::noteBootedAtoms(
	int vertex_index,
	std::vector< bool  > const & booted_atoms,
	std::vector< DotsForAtom > const & dotsInHOO
)
{
	int whichvertex = whichVertex( vertex_index );
	for ( int ii = 0; ii < num_atoms_[ whichvertex ]; ++ii )
	{
		if ( booted_atoms[ ii ] )
		{
			dotsForAtom_[ whichvertex ][ ii ].subtractDotsOn( dotsInHOO[ ii ] );
		}
	}
}

void DegreeTwoEdge_ths::tellVerticesToNotScoreDotsScoredOnD2Edge() const
{
	for (int ii = 0; ii < 2; ++ii)
	{
		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj )
		{
			vertex_ptrs_[ ii ]->removeDotsScoredOnD2Edge( jj, dotsForAtom_[ ii ][ jj ]);
		}
	}
}

void DegreeTwoEdge_ths::dropDotsScoredOnD3Edge(
	int vertex,
	int atom,
	DotsForAtom const & dots )
{
	dotsForAtom_[ whichVertex( vertex )][ atom ].subtractDotsOn( dots );
}


void DegreeTwoEdge_ths::finalizeScoreVectors()
{
	//std::cerr << "Entering D2Edge::finalizeScoreVectors() for " << vertex_indices_[ 0 ] << " " << vertex_indices_[ 1 ] << std::endl;
	for (int ii = 0; ii < 2; ++ii )
	{
		int count_atoms_to_score = 0;
		std::vector< bool > score_atom( num_atoms_[ ii ], false );
		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj )
		{
			if ( dotsForAtom_[ ii ][ jj ].numDotsOn() != 0)
			{
				score_atom[ jj ] = true;
				++count_atoms_to_score;
			}
		}
	
		dotsToScoreOnAtoms_[ ii ].resize( count_atoms_to_score );
		int count = 0;
		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj)
		{
			if ( ! score_atom[ jj ] ) continue;
			dotsToScoreOnAtoms_[ ii ][ count ].first = atoms_[ ii ][ jj ];
			dotsToScoreOnAtoms_[ ii ][ count ].second = & dotsForAtom_[ ii ][ jj ];
			++count;
			//std::cerr << "Score: " << atoms_[ ii ][ jj ] << " " << dotsForAtom_[ ii ][ jj ] << std::endl;
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
				scores_[ ii ][ jj ] += xyz_->determineScoreForMover( movers_[ kk ], dotsToScoreOnAtoms_[ kk ], temp );
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
	std::cerr << "CRITICAL ERROR IN whichVertex(" << vertex_index << ") called on edge [" <<
	std::cerr << vertex_indices_[ 0 ] << ", " << vertex_indices_[ 1 ] << "]" << std::endl;
	assert(false);
	exit(1);
        return 0; // to avoid warnings
}

bool DegreeTwoEdge_ths::nothingToScore() const
{
	//std::cerr << "nothingToScore?" << std::endl;
	for (int ii = 0; ii < 2; ++ii )
	{
		if ( ! dotsToScoreOnAtoms_[ ii ].empty() ) return false;
	}
	//std::cerr << "Indeed, nothing to score" << std::endl;
	return true;
}

bool DegreeTwoEdge_ths::nothingToScore( int vertex ) const
{
	return dotsToScoreOnAtoms_[ vertex ].empty();
}

//----------------------------------------------------------------------------//
//------------------------- Degree Three Edge Class --------------------------//
//----------------------------------------------------------------------------//

DegreeThreeEdge_ths::DegreeThreeEdge_ths(
	GraphToHoldScores * owner,
	AtomPositions * xyz,
	DotSphManager * dotSphereManager,
	int vertex1,
	int vertex2,
	int vertex3
)
:
	owner_( owner ),
	xyz_( xyz ),
	dotSphereManager_( dotSphereManager )
{
	vertex_indices_[ 0 ] = vertex1;
	vertex_indices_[ 1 ] = vertex2;
	vertex_indices_[ 2 ] = vertex3;
	//std::cerr << "Creating degree-3 edge";
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
		movers_[ ii ] = vertex_ptrs_[ ii ]->getMover();
		atoms_[ ii ].resize( num_atoms_[ ii ] );
		vertex_ptrs_[ ii ]->getAtoms( atoms_[ ii ] );
		atomHasThreeWay_[ ii ].resize( num_atoms_[ ii ]);
		std::fill( atomHasThreeWay_[ ii ].begin(), atomHasThreeWay_[ ii ].end(), false );

		dotsForAtom_[ ii ].resize( atoms_[ ii ].size() );
		for (int jj = 0; jj < atoms_[ ii ].size(); ++jj)
		{
			dotsForAtom_[ ii ][ jj ].setDotManager( dotSphereManager_ );
			dotsForAtom_[ ii ][ jj ].setAtom( atoms_[ ii ][ jj ] );
		}
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

DotsForAtom const & 
DegreeThreeEdge_ths::getDotsForAtom( int vertex, int atom )
{
	return dotsForAtom_[ whichVertex( vertex ) ][ atom ];
}

void
DegreeThreeEdge_ths::setMoverAtomHas3Way(
	int vertex,
	int atom
)
{
	atomHasThreeWay_[ whichVertex( vertex ) ][ atom ] = true;
}

void
DegreeThreeEdge_ths::detectDotsToScoreOnEdge()
{
	//std::cerr << "Entering D3Edge::detectDotsToScoreOnEdge: ";
	//std::cerr << vertex_indices_[ 0 ] << " ";
	//std::cerr << vertex_indices_[ 1 ] << " ";
	//std::cerr << vertex_indices_[ 2 ] << std::endl;
	for (int ii = 0; ii < 3; ++ii )
	{
		//std::cerr << "ii: " << ii << std::endl;
		DegreeTwoEdge_ths * d2edgesIncidentOnII[ 2 ] = {0,0};
		int count = 0;
		for (int jj = 0; jj < 3; ++jj )
		{
			if ( d2e_[ jj ]->incidentUpon( vertex_indices_[ ii ] ))
			{
				d2edgesIncidentOnII[ count ] = d2e_[ jj ];
				++count;
			}
		}
		assert( count == 2 );
		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj)
		{
			//std::cerr << "jj: " << jj << std::endl;
			if ( ! atomHasThreeWay_[ ii ][ jj ] )
			{
				continue;
			}
			
			//std::cerr << "ThreeWay jj: " << jj << std::endl;
			d2edgesIncidentOnII[ 0 ]->getDotsForAtom( vertex_indices_[ ii ], jj ).intersect( 
				d2edgesIncidentOnII[ 1 ]->getDotsForAtom( vertex_indices_[ ii ], jj ),
				dotsForAtom_[ ii ][ jj ] );
		}
	}
	//std::cerr << "Leaving D3Edge::detectDotsToScoreOnEdge: ";
	//std::cerr << vertex_indices_[ 0 ] << " ";
	//std::cerr << vertex_indices_[ 1 ] << " ";
	//std::cerr << vertex_indices_[ 2 ] << std::endl;

}


void
DegreeThreeEdge_ths::detectDotsDoublyCounted( DegreeThreeEdge_ths const * other ) const
{
	int countSharedVertices = 0;
	for (int ii = 0; ii < 3; ++ii)
	{
		for (int jj = 0; jj < 3; ++jj)
		{
			if ( vertex_indices_[ ii ] == vertex_indices_[ jj ] )	
			{
				++countSharedVertices;
			}
		}
	}
	
	if ( countSharedVertices != 2 ) return;
	
	for (int ii = 0; ii < 3; ++ii)
	{
		for (int jj = 0; jj < 3; ++jj)
		{
			if ( vertex_indices_[ ii ] == vertex_indices_[ jj ] )	
			{
				lookForSharedDots( other, ii, jj );
			}
		}
	}
}

void 
DegreeThreeEdge_ths::tellD2EdgesToNotScoreDotsScoredOnD3Edge() const
{
	for (int ii = 0; ii < 3; ++ii )
	{
		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj )
		{
			if ( dotsForAtom_[ ii ][ jj ].numDotsOn() == 0 )
			{
				continue;
			}
			for (int kk = 0; kk < 3; ++kk )
			{
				if ( ! d2e_[ kk ]->incidentUpon( vertex_indices_[ ii ] ) )
				{ 
					continue;
				}
				d2e_[ kk ]->dropDotsScoredOnD3Edge(
					vertex_indices_[ ii ],
					jj,
					dotsForAtom_[ ii ][ jj ]
				);
			}
		}
	}
}

void 
DegreeThreeEdge_ths::dropDotsScoredOnD4Edge(
	int vertex_index,
	int atom,
	DotsForAtom const & dotsOnD4Edge
)
{
	dotsForAtom_[ whichVertex( vertex_index )][ atom ].subtractDotsOn( dotsOnD4Edge );
}

void
DegreeThreeEdge_ths::noteBootedAtoms(
	int vertex_index,
	std::vector< bool  > const & booted_atoms,
	std::vector< DotsForAtom > const & dotsInHOO
)
{
	int whichvertex = whichVertex( vertex_index );
	for ( int ii = 0; ii < num_atoms_[ whichvertex ]; ++ii )
	{
		if ( booted_atoms[ ii ] )
		{
			dotsForAtom_[ whichvertex ][ ii ].subtractDotsOn( dotsInHOO[ ii ] );
		}
	}
}

void DegreeThreeEdge_ths::finalizeScoreVectors()
{
	//std::cerr << "Entering D3Edge::finalizeScoreVectors() for " << vertex_indices_[ 0 ] << " " << vertex_indices_[ 1 ] << " " << vertex_indices_[ 2 ] << std::endl;
	for (int ii = 0; ii < 3; ++ii )
	{
		//std::cerr << "ii: " << ii << std::endl;
		int count_atoms_to_score = 0;
		std::vector< bool > score_atom( num_atoms_[ ii ], false );
		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj )
		{
			//std::cerr << "jj: " << jj << std::endl;
			if ( dotsForAtom_[ ii ][ jj ].numDotsOn() != 0  )
			{
				score_atom[ jj ] = true;
				++count_atoms_to_score;
			}
		}
	
		//std::cerr << "count_atoms_to_score " << count_atoms_to_score  << std::endl;
		if (count_atoms_to_score == 0) continue;
		
		
		dotsToScoreOnAtoms_[ ii ].resize( count_atoms_to_score );
		int count = 0;
		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj)
		{
			//std::cerr << "jj: " << jj << std::endl;
			if ( ! score_atom[ jj ] ) continue;
			//std::cerr << "Score jj: " << jj << " count=" << count << " count_atoms_to_score "<< count_atoms_to_score << std::endl;
			
			dotsToScoreOnAtoms_[ ii ][ count ].first = atoms_[ ii ][ jj ];
			dotsToScoreOnAtoms_[ ii ][ count ].second = &( dotsForAtom_[ ii ][ jj ] );
			//std::cerr << "Score: " << atoms_[ ii ][ jj ] << " " << dotsForAtom_[ ii ][ jj ] << std::endl;
			++count;
		}
	}
	//std::cerr << "Leaving D3Edge::finalizeScoreVectors() for " << vertex_indices_[ 0 ] << " " << vertex_indices_[ 1 ] << " " << vertex_indices_[ 2 ] << std::endl;
	
}

void
DegreeThreeEdge_ths::score()
{

	//std::cerr << "Scoring degree 3 hyperedge: " << vertex_indices_[ 0 ] << " ";
	//std::cerr << vertex_indices_[ 1 ] << " ";
	//std::cerr << vertex_indices_[ 2 ] << " with num states: ";
	//std::cerr << num_states_[ 0 ] << " ";
	//std::cerr << num_states_[ 1 ] << " ";
	//std::cerr << num_states_[ 2 ] << std::endl;

	if ( nothingToScore() ) return;

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
					scores_table_[ ii ][ jj ][ kk ] += 
						xyz_->determineScoreForMover(
						movers_[ ll ], 
						dotsToScoreOnAtoms_[ ll ],
						temp );
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
	int enabled_state3
) const
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
        return 0; // to avoid warnings
}

void
DegreeThreeEdge_ths::lookForSharedDots( 
	DegreeThreeEdge_ths const * other,
	int vertex_this,
	int vertex_other
) const
{
	std::list< int > verts;
	for (int jj = 0; jj < 3; ++jj)
	{
		verts.push_back( vertex_indices_[ jj ] );
		verts.push_back( other->vertex_indices_[ jj ] );
	}
	verts.sort();
	verts.unique();
	assert( verts.size() == 4 );
	
	int v[ 4 ];
	int count = 0;
	for (std::list< int >::iterator iter = verts.begin(); iter != verts.end(); ++iter )
	{
		v[ count ] = *iter;
		++count;
	}
	
	for ( int ii = 0; ii < num_atoms_[ vertex_this ]; ++ii )
	{
		if ( dotsForAtom_[ vertex_this ][ ii ].anyDotsShared(
			other->dotsForAtom_[ vertex_other ][ ii ] ) )
		{
			DegreeFourEdge_ths * d4e = owner_->getDegree4Edge( v[ 0 ], v[ 1 ], v[ 2 ], v[ 3 ] );
			d4e->setAtomHas4WayOverlap( vertex_indices_[ vertex_this ], ii );
		}
	}
}

bool DegreeThreeEdge_ths::nothingToScore() const
{
	for (int ii = 0; ii < 3; ++ii )
	{
		if ( ! dotsToScoreOnAtoms_[ ii ].empty() ) return false;
	}
	return true;
}

bool DegreeThreeEdge_ths::nothingToScore( int vertex ) const
{
	return dotsToScoreOnAtoms_[ vertex ].empty();
}

//----------------------------------------------------------------------------//
//------------------------- Degree Four Edge Class --------------------------//
//----------------------------------------------------------------------------//

DegreeFourEdge_ths::DegreeFourEdge_ths(
	GraphToHoldScores * owner,
	AtomPositions * xyz,
	DotSphManager * dotSphereManager,
	int vertex1,
	int vertex2,
	int vertex3,
	int vertex4
)
:
	owner_( owner ),
	xyz_( xyz ),
	dotSphereManager_( dotSphereManager )
{
	vertex_indices_[ 0 ] = vertex1;
	vertex_indices_[ 1 ] = vertex2;
	vertex_indices_[ 2 ] = vertex3;
	vertex_indices_[ 3 ] = vertex4;
	
	//if ( xyz_->outputNotice() ) 
	//{
	//	std::cerr << "Creating degree four edge";
	//	for (int ii = 0; ii < 4; ++ii)
	//	{
	//		std::cerr << " " << vertex_indices_[ ii ];
	//	}
	//	std::cerr << std::endl;
	//}
	
	for (int ii = 0; ii < 4; ++ii)
	{
		vertex_ptrs_[ ii ] = owner_->getVertexPtr( vertex_indices_[ ii ] );
		pos_in_verts_edge_list_[ ii ] = vertex_ptrs_[ ii ]->addEdge( this );
		num_states_[ ii ] = vertex_ptrs_[ ii ]->getNumStates();
		num_atoms_[ ii ] = vertex_ptrs_[ ii ]->getNumAtoms();
		movers_[ ii ] = vertex_ptrs_[ ii ]->getMover();
		atoms_[ ii ].resize( num_atoms_[ ii ] );
		vertex_ptrs_[ ii ]->getAtoms( atoms_[ ii ] );
		atomHasFourWay_[ ii ].resize( num_atoms_[ ii ]);
		std::fill( atomHasFourWay_[ ii ].begin(), atomHasFourWay_[ ii ].end(), false );
		
		dotsForAtom_[ ii ].resize( atoms_[ ii ].size() );
		for (int jj = 0; jj < atoms_[ ii ].size(); ++jj)
		{
			dotsForAtom_[ ii ][ jj ].setDotManager( dotSphereManager_ );
			dotsForAtom_[ ii ][ jj ].setAtom( atoms_[ ii ][ jj ] );
		}
	}
	
	scores_table_ = new float***[ num_states_[ 0 ] ];	
	for (int ii = 0; ii < num_states_[ 0 ]; ++ii)
	{
		scores_table_[ ii ] = new float**[ num_states_[ 1 ] ];
		for (int jj = 0; jj < num_states_[ 1 ]; ++jj )
		{
			scores_table_[ ii ][ jj ] = new float* [ num_states_[ 2 ] ];
			for (int kk = 0; kk < num_states_[ 2 ]; ++kk )
			{
				scores_table_[ ii ][ jj ][ kk ] = new float[ num_states_[ 3 ]];
				for (int ll = 0; ll < num_states_[ 3 ]; ++ll )
				{
					scores_table_[ ii ][ jj ][ kk ][ ll ] = 0;
				}
			}
		}
	}
}	

DegreeFourEdge_ths::~DegreeFourEdge_ths()
{
	if (scores_table_ != 0 ) 
	{	
		for (int ii = 0; ii < num_states_[ 0 ]; ++ii)
		{
			if ( scores_table_[ ii ] == 0 ) continue;
			for (int jj = 0; jj < num_states_[ 1 ]; ++jj )
			{
				for (int kk = 0; kk < num_states_[ 2 ]; ++kk)
				{
					delete [] scores_table_[ ii ][ jj ][ kk ]; scores_table_[ ii ][ jj ][ kk ] = 0;
				}
				delete [] scores_table_[ ii ][ jj ];scores_table_[ ii ][ jj ] = 0;
			}
			delete [] scores_table_[ ii ];scores_table_[ ii ] = 0;
		}
		delete [] scores_table_;scores_table_ = 0;
	}
	for (int ii = 0; ii < 4; ++ii)
	{
		vertex_ptrs_[ ii ]->dropEdge( pos_in_verts_edge_list_[ ii ] );
	}
	owner_->dropEdge( pos_in_owners_edge_list_ );
}

void
DegreeFourEdge_ths::setPosInOwnersEdgeList( std::list< DegreeFourEdge_ths * >::iterator iter )
{
	pos_in_owners_edge_list_ = iter;
}

bool DegreeFourEdge_ths::sameEdge( int vertex1, int vertex2, int vertex3, int vertex4 ) const
{
	return ( vertex1 == vertex_indices_[ 0 ] 
		&& vertex2 == vertex_indices_[ 1 ] 
		&& vertex3 == vertex_indices_[ 2 ] 
		&& vertex4 == vertex_indices_[ 3 ]);
}

int DegreeFourEdge_ths::getFirstNodeIndex() const
{
	return vertex_indices_[ 0 ];
}

int DegreeFourEdge_ths::getSecondNodeIndex() const
{
	return vertex_indices_[ 1 ];
}

int DegreeFourEdge_ths::getThirdNodeIndex() const
{
	return vertex_indices_[ 2 ];
}

int DegreeFourEdge_ths::getFourthNodeIndex() const
{
	return vertex_indices_[ 3 ];
}

void DegreeFourEdge_ths::setAtomHas4WayOverlap(
	int vertex,
	int atom
)
{
	atomHasFourWay_[ whichVertex( vertex )][ atom ] = true;
}

void DegreeFourEdge_ths::detectDotsToScoreOnEdge()
{
	for (int ii = 0; ii < 4; ++ii)
	{
		int others[ 3 ];
		int count = 0;
		for (int jj = 0; jj < 4; ++jj)
		{
			if (ii == jj ) continue;
			others[ count ] = jj;
			++count;
		}
		
		int groupOfThreeA[ 3 ] = { -1, -1, -1};
		int groupOfThreeB[ 3 ] = { -1, -1, -1};
		
		int offsetA = 0;
		for (int jj = 0; jj < 2; ++jj)
		{
			if ( ii < others[ jj ] )
			{
				groupOfThreeA[ jj ] = ii;
				offsetA = 1;
			}
			groupOfThreeA[ jj + offsetA ] = others[ jj ];
		}
		int offsetB = 0;
		for (int jj = 1; jj < 3; ++jj)
		{
			if ( ii < others[ jj ] )
			{
				groupOfThreeB[ jj - 1 ] = ii;
				offsetB = 1;
			}
			groupOfThreeB[ jj - 1 + offsetB ] = others[ jj ];
		}
		//std::cerr << "D4E: ii: " << ii << " groupA: " << groupOfThreeA[ 0 ] << " " << groupOfThreeA[ 1 ] << " " << groupOfThreeA[ 2 ] << std::endl;
		//std::cerr << "D4E: ii: " << ii << " groupB: " << groupOfThreeB[ 0 ] << " " << groupOfThreeB[ 1 ] << " " << groupOfThreeB[ 2 ] << std::endl;
		
		DegreeThreeEdge_ths * d3eA = owner_->getDegree3Edge( groupOfThreeA[ 0 ], groupOfThreeA[ 1 ], groupOfThreeA[ 2 ]);
		DegreeThreeEdge_ths * d3eB = owner_->getDegree3Edge( groupOfThreeB[ 0 ], groupOfThreeB[ 1 ], groupOfThreeB[ 2 ]);

		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj)
		{
			if ( ! atomHasFourWay_[ ii ][ jj ] ) continue;
			{
				d3eA->getDotsForAtom( vertex_indices_[ ii ], jj ).intersect(
					d3eB->getDotsForAtom( vertex_indices_[ ii ], jj ),
					dotsForAtom_[ ii ][ jj ] );
			}
		}
	}	
}

void
DegreeFourEdge_ths::detectDotsDoublyCounted( DegreeFourEdge_ths const * other ) const
{
	int countSharedVertices = 0;
	for (int ii = 0; ii < 4; ++ii)
	{
		for (int jj = 0; jj < 4; ++jj)
		{
			if ( vertex_indices_[ ii ] == vertex_indices_[ jj ] )	
			{
				++countSharedVertices;
			}
		}
	}
	
	if ( countSharedVertices != 3 ) return;
	
	for (int ii = 0; ii < 4; ++ii)
	{
		for (int jj = 0; jj < 4; ++jj)
		{
			if ( vertex_indices_[ ii ] == vertex_indices_[ jj ] )	
			{
				lookForSharedDots( other, ii, jj );
			}
		}
	}
}


void 
DegreeFourEdge_ths::tellD3EdgesToNotScoreDotsScoredOnD4Edge() const
{
	//std::cerr << "tellD3EdgesToNotScoreBGandMoverAtomsScoredOnD4Edge: [";
	//for (int nn = 0; nn < 4; ++nn)
	//{
	//	std::cerr<< " " << vertex_indices_[ nn ];
	//}
	//std::cerr << std::endl;

	for (int ii = 0; ii < 4; ++ii )
	{
		int subset[ 3 ];
		int count = 0;
		for (int jj = 0; jj < 4; ++jj)
		{
			if (ii == jj ) continue;
			subset[ count ] = jj;
			++count;
		}
		
		DegreeThreeEdge_ths * d3e = owner_->getDegree3Edge( 
			vertex_indices_[ subset[ 0 ]],
			vertex_indices_[ subset[ 1 ]],
			vertex_indices_[ subset[ 2 ]]);
		
		for (int jj = 0; jj < 3; ++jj)
		{
			for (int kk = 0; kk < num_atoms_[ subset[ jj ] ]; ++kk)
			{
				if ( ! atomHasFourWay_[ subset[ jj ]][ kk ] ) continue;
				d3e->dropDotsScoredOnD4Edge(
					vertex_indices_[ subset[ jj ] ],
					kk,
					dotsForAtom_[ subset[ jj ]][ kk ] );
			}
		}
	}		
}			

void
DegreeFourEdge_ths::noteBootedAtoms(
	int vertex_index,
	std::vector< bool  > const & booted_atoms,
	std::vector< DotsForAtom > const & dotsInHOO
)
{
	int whichvertex = whichVertex( vertex_index );
	for ( int ii = 0; ii < num_atoms_[ whichvertex ]; ++ii )
	{
		if ( booted_atoms[ ii ] )
		{
			dotsForAtom_[ whichvertex ][ ii ].subtractDotsOn( dotsInHOO[ ii ] );
		}
	}
}

void DegreeFourEdge_ths::finalizeScoreVectors()
{
	//std::cerr << "Entering D2Edge::finalizeScoreVectors() for " << vertex_indices_[ 0 ] << " " << vertex_indices_[ 1 ]<< " " << vertex_indices_[ 2 ]<< " " << vertex_indices_[ 3 ] << std::endl;
	for (int ii = 0; ii < 4; ++ii )
	{
		int count_atoms_to_score = 0;
		std::vector< bool > score_atom( num_atoms_[ ii ], false );
		for (int jj = 0; jj < num_atoms_[ ii ]; ++jj )
		{
			if ( dotsForAtom_[ ii ][ jj ].numDotsOn() != 0 )
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
			dotsToScoreOnAtoms_[ ii ][ count ].first = atoms_[ ii ][ jj ];
			dotsToScoreOnAtoms_[ ii ][ count ].second = &( dotsForAtom_[ ii ][ jj ] );
			++count;
		}
	}
	
}

void
DegreeFourEdge_ths::score()
{
	if ( nothingToScore() ) return;
	//std::cerr << "Scoring degree 4 hyperedge: " << vertex_indices_[ 0 ] << " ";
	//std::cerr << vertex_indices_[ 1 ] << " ";
	//std::cerr << vertex_indices_[ 2 ] << " ";
	//std::cerr << vertex_indices_[ 3 ] << std::endl;
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
				for (int ll = 0; ll < num_states_[ 3 ]; ++ll )
				{
					movers_[ 3 ]->initOrientation( *xyz_ );
					for ( int mm = 0; mm < 4; ++mm)
					{
						if ( nothingToScore( mm ) ) continue;
						scores_table_[ ii ][ jj ][ kk ][ ll ] += 
							xyz_->determineScoreForMover(
							movers_[ mm ], 
							dotsToScoreOnAtoms_[ mm ],
							temp );

					}
					//std::cerr << "d4e score: " << ii << " " << jj << " " << kk << " " << ll;
					//std::cerr<< " = " << scores_table_[ ii ][ jj ][ kk ][ ll ] << std::endl;

					movers_[ 3 ]->nextOrientation( *xyz_ );
				}
				movers_[ 2 ]->nextOrientation( *xyz_ );
			}
			movers_[ 1 ]->nextOrientation( *xyz_ );
		}
		movers_[ 0 ]->nextOrientation( *xyz_ );
	}
}

float
DegreeFourEdge_ths::getScoreForEnabledStates(
	int enabled_state1,
	int enabled_state2,
	int enabled_state3,
	int enabled_state4) const
{
	int states[] = {enabled_state1, enabled_state2, enabled_state3, enabled_state4};
	int original_states[4];
	for (int ii = 0; ii < 4; ++ii)
	{
		original_states[ ii ] = vertex_ptrs_[ ii ]->convertEnabledStateToRegularState( states[ ii ] );
	}
	return scores_table_[ original_states[ 0 ]][ original_states[ 1 ] ][ original_states[ 2 ]][ original_states[ 3 ]];
}


float DegreeFourEdge_ths::getScoreForAssignedState() const
{
	return scores_table_[ vertex_ptrs_[ 0 ]->getAssignedState() ]
		[ vertex_ptrs_[ 1 ]->getAssignedState() ]
		[ vertex_ptrs_[ 2 ]->getAssignedState() ]
		[ vertex_ptrs_[ 3 ]->getAssignedState() ];
}

int
DegreeFourEdge_ths::whichVertex( int vertex ) const
{
	for (int ii = 0; ii < 4; ++ii)
	{
		if ( vertex_indices_[ ii ] == vertex ) return ii;
	}
	std::cerr << "DegreeFourEdge_ths: CRITICAL ERROR. Did not find sought vertex: " << vertex << " out of [ ";
	std::cerr << vertex_indices_[ 0 ] << ", " << vertex_indices_[ 1 ] << ", ";
	std::cerr << vertex_indices_[ 2 ] << ", " << vertex_indices_[ 3 ];
	std::cerr << "]" << std::endl;
	assert( false );
	exit(1);
        return 0; // to avoid warnings
}

void
DegreeFourEdge_ths::lookForSharedDots( 
	DegreeFourEdge_ths const * other,
	int vertex_this,
	int vertex_other
) const
{
	std::list< int > verts;
	for (int jj = 0; jj < 4; ++jj)
	{
		verts.push_back( vertex_indices_[ jj ] );
		verts.push_back( other->vertex_indices_[ jj ] );
	}
	verts.sort();
	verts.unique();
		
	for ( int ii = 0; ii < num_atoms_[ vertex_this ]; ++ii )
	{
		if ( dotsForAtom_[ vertex_this ][ ii ].anyDotsShared(
			other->dotsForAtom_[ vertex_other ][ ii ] ) )
		{
			DotsForAtom dotsInHOO;
			dotsInHOO.setDotManager( dotSphereManager_ );
			dotsInHOO.setAtom( atoms_[ vertex_this ][ ii ] );

			dotsForAtom_[ vertex_this ][ ii ].intersect(
				other->dotsForAtom_[ vertex_other ][ ii ],
				dotsInHOO );
			
			vertex_ptrs_[ vertex_this ]->boot( ii, dotsInHOO, verts );
		}
	}
}

bool DegreeFourEdge_ths::nothingToScore() const
{
	for (int ii = 0; ii < 4; ++ii )
	{
		if ( ! dotsToScoreOnAtoms_[ ii ].empty() ) return false;
	}
	return true;
}

bool DegreeFourEdge_ths::nothingToScore( int vertex ) const
{
	return dotsToScoreOnAtoms_[ vertex ].empty();
}

//----------------------------------------------------------------------------//
//------------------------------ Graph Class --------------------------------//
//----------------------------------------------------------------------------//


GraphToHoldScores::GraphToHoldScores( 
	AtomPositions * xyz,
	DotSphManager * dotSphereManager,
	std::vector< int > const & num_states,
	std::vector< Mover * > const & movers )
:
	xyz_( xyz ),
	dotSphereManager_( dotSphereManager ),
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
	//std::cerr << "addDegreeTwoEdges();" << std::endl;
	addDegreeTwoEdges();
	//std::cerr << "addDegreeThreeEdges();" << std::endl;
	addDegreeThreeEdges();
	//std::cerr << "detectFourWayInteractions();" << std::endl;
	detectFourWayInteractions();
	//std::cerr << "detectFiveWayInteractions();" << std::endl;
	detectFiveWayInteractions();
	//std::cerr << "cascadeScoringInfo();" << std::endl;
	cascadeScoringInfo();
	//std::cerr << "finalizeScoreVectors();" << std::endl;
	finalizeScoreVectors();
	if (xyz_->outputNotice() ) std::cerr << " Computing dot scores" << std::endl;
	score();
}

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
	std::list< DegreeFourEdge_ths * >::iterator iter4 = deg4edges_.begin();
	while ( iter4 != deg4edges_.end() )
	{
		std::list< DegreeFourEdge_ths  * >::iterator nextiter4 = iter4;
		++nextiter4;
		delete (*iter4);
		iter4 = nextiter4;
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

DegreeFourEdge_ths* 
GraphToHoldScores::getDegree4Edge( int fn, int sn, int tn, int fthn )
{
	nodeInRange( fn );
	nodeInRange( sn );
	nodeInRange( tn );
	nodeInRange( fthn );
	assert( fn < sn && sn < tn && tn < fthn );
	DegreeFourEdge_ths * edge = findDeg4Edge(fn, sn, tn, fthn );
	if (edge == 0 ) edge = addDegree4Edge( fn, sn, tn, fthn );
	return edge;
}
	
bool GraphToHoldScores::degree2EdgeExists( int fn, int sn ) const
{	return (findDeg2Edge( fn, sn ) != 0 ); }

bool GraphToHoldScores::degree3EdgeExists( int fn, int sn, int tn ) const
{	return (findDeg3Edge( fn, sn, tn ) != 0 ); }

bool GraphToHoldScores::degree4EdgeExists( int fn, int sn, int tn, int fthn ) const
{	return ( findDeg4Edge( fn, sn, tn, fthn ) != 0 ); }

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

void GraphToHoldScores::dropEdge( std::list< DegreeFourEdge_ths * >::iterator iter )
{
	deg4edges_.erase( iter );
	--num_deg4edges_;
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

bool GraphToHoldScores::anyNodeWithAllStatesDisabled() const
{

	for (int ii = 0; ii < getNumNodes(); ++ii )
	{
		if ( vertices_[ ii ]->getNumEnabledStates() == 0 )
		{
			return true;
		}
	}
	return false;
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

std::vector< std::pair< AtomDescr, DotsForAtom * > > 
GraphToHoldScores::getAtomsInHighOrderOverlapForNode( int vertex_id ) const
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
        // conversion to const_iterator works around Tru64 cxx issue
        return std::list< DegreeTwoEdge_ths * >::const_iterator(d2eiter_)
            == deg2edges_.end();
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
        // conversion to const_iterator works around Tru64 cxx issue
        return std::list< DegreeThreeEdge_ths * >::const_iterator(d3eiter_)
	    == deg3edges_.end();
}

void GraphToHoldScores::setD4EIteratorAtBegining()
{
	d4eiter_ = deg4edges_.begin();
}

void GraphToHoldScores::incrementD4EIterator()
{
	++d4eiter_;
}

bool GraphToHoldScores::getD4EIteratorAtEnd() const
{
        // conversion to const_iterator works around Tru64 cxx issue
        return std::list< DegreeFourEdge_ths * >::const_iterator(d4eiter_)
            == deg4edges_.end();
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

int GraphToHoldScores::getFirstIndexFocusedD4E() const
{
	return (*d4eiter_)->getFirstNodeIndex();
}

int GraphToHoldScores::getSecondIndexFocusedD4E() const
{
	return (*d4eiter_)->getSecondNodeIndex();
}

int GraphToHoldScores::getThirdIndexFocusedD4E() const
{
	return (*d4eiter_)->getThirdNodeIndex();
}

int GraphToHoldScores::getFourthIndexFocusedD4E() const
{
	return (*d4eiter_)->getFourthNodeIndex();
}

float
GraphToHoldScores::getScoreForFocusedD4E( int state1, int state2, int state3, int state4) const
{
	return (*d4eiter_)->getScoreForEnabledStates( state1, state2, state3, state4 );
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
	for (std::list< DegreeFourEdge_ths* >::const_iterator iter = deg4edges_.begin();
		iter != deg4edges_.end(); ++iter)
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
		vertices_[ ii ] = new Vertex_ths( 
			this, xyz_,  dotSphereManager_,
			ii, num_states[ ii ]);
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
	
	for (std::list< DegreeTwoEdge_ths* >::iterator iter = deg2edges_.begin();
		iter != deg2edges_.end(); ++iter)
	{
		(*iter)->detectDotsToScoreOnEdge();
	}
}
void
GraphToHoldScores::addDegreeThreeEdges()
{

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
	
	for (std::list< DegreeThreeEdge_ths * >::iterator iter = deg3edges_.begin();
		iter != deg3edges_.end(); ++iter )
	{
		(*iter)->detectDotsToScoreOnEdge();
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
			(*iter_outer)->detectDotsDoublyCounted( *iter_inner );
			++iter_inner;
		}
		++iter_outer;
	}

	for (std::list< DegreeFourEdge_ths * >::iterator iter = deg4edges_.begin();
		iter != deg4edges_.end(); ++iter )
	{
		(*iter)->detectDotsToScoreOnEdge();
	}
}

void GraphToHoldScores::detectFiveWayInteractions()
{
	//compare all degree4 hyperedges
	std::list< DegreeFourEdge_ths * >::iterator iter_outer = deg4edges_.begin();
	while ( iter_outer != deg4edges_.end() )
	{
		std::list< DegreeFourEdge_ths * >::iterator iter_inner = iter_outer;
		++iter_inner;
		while ( iter_inner != deg4edges_.end() )
		{
			(*iter_outer)->detectDotsDoublyCounted( *iter_inner );
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
		(*iter)->tellVerticesToNotScoreDotsScoredOnD2Edge();
	}
	
	for( std::list< DegreeThreeEdge_ths * >::iterator iter = deg3edges_.begin();
		iter != deg3edges_.end(); ++iter )
	{
		(*iter)->tellD2EdgesToNotScoreDotsScoredOnD3Edge();
	}

	for( std::list< DegreeFourEdge_ths * >::iterator iter = deg4edges_.begin();
		iter != deg4edges_.end(); ++iter )
	{
		(*iter)->tellD3EdgesToNotScoreDotsScoredOnD4Edge();
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
	
	for( std::list< DegreeFourEdge_ths * >::iterator iter = deg4edges_.begin();
		iter != deg4edges_.end(); ++iter )
	{
		(*iter)->finalizeScoreVectors();
	}
	//std::cerr << "Leaving finalize score vectors" << std::endl;
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

	for( std::list< DegreeFourEdge_ths * >::iterator iter = deg4edges_.begin();
		iter != deg4edges_.end(); ++iter )
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
	
	DegreeTwoEdge_ths * new_edge = new DegreeTwoEdge_ths( this, xyz_,  dotSphereManager_, node1, node2 );
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
	DegreeThreeEdge_ths * new_edge = new DegreeThreeEdge_ths( 
		this, xyz_,  dotSphereManager_,
		fn, sn, tn );
	
	deg3edges_.push_front( new_edge );
	new_edge->setPosInOwnersEdgeList( deg3edges_.begin() );
	++num_deg3edges_;
	return new_edge;
}

DegreeFourEdge_ths*
GraphToHoldScores::findDeg4Edge( int fn, int sn, int tn, int fthn ) const
{
	for (std::list< DegreeFourEdge_ths * >::const_iterator iter = deg4edges_.begin();
		iter != deg4edges_.end(); ++iter )
	{
		if ( (*iter)->sameEdge( fn, sn, tn, fthn ) )
			return *iter;
	}
	return 0;
}

DegreeFourEdge_ths*
GraphToHoldScores::addDegree4Edge( int fn, int sn, int tn, int fthn )
{
	//std::cerr << "Adding degree 4 edge: " << fn << ", " << sn << ", " << tn << ", " << fthn << std::endl;
	DegreeFourEdge_ths * new_edge = new DegreeFourEdge_ths( 
		this, xyz_, dotSphereManager_,
		fn, sn, tn, fthn );
	
	deg4edges_.push_front( new_edge );
	new_edge->setPosInOwnersEdgeList( deg4edges_.begin() );
	++num_deg4edges_;
	return new_edge;
}
