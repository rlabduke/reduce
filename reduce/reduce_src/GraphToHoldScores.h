#ifndef GRAPH_TO_HOLD_SCORES_H
#define GRAPH_TO_HOLD_SCORES_H

#include <list>
#include <vector>

class Mover;
class AtomPositions;
class DotSphManager;

#include "AtomDescr.h"

void sort_three(int, int, int, int &, int &, int & );

template < class T >
void set_intersect_apl( std::list< T > const & list1, std::list< T > const & list2, std::list< T > & intersection );

template < class T >
void set_difference_apl( std::list< T > & list1, std::list< T > const & list2 );

template < class T >
void set_difference_apl( std::list< T > & list1, std::vector< T > const & list2 );

class DotsForAtom;
class GraphToHoldScores;
class Vertex_ths;
class DegreeTwoEdge_ths;
class DegreeThreeEdge_ths;
class DegreeFourEdge_ths;

class DotsForAtom
{
public:
	DotsForAtom();
	~DotsForAtom();
	DotsForAtom( DotsForAtom const & rhs );
	DotsForAtom const & operator = ( DotsForAtom const & rhs );
	
	void setDotManager( DotSphManager * dotSphereManager );
	void setAtom( AtomDescr theAtom );
	void findContained( std::list< AtomDescr > const & atoms );
	void turnOnAllDots();
	void turnOffAllDots();
	bool dotOn( int ) const;
	int numDotsOn() const;

	bool anyDotsShared( DotsForAtom const & other ) const;
	void intersect( DotsForAtom const & other, DotsForAtom & intersection ) const;
	void subtractDotsOn( DotsForAtom const & other );
	void addDots( DotsForAtom const & other );

	void print(std::ostream & os ) const;
private:
	void updateIndsOfDotsOn() const;
	
	AtomDescr theAtom_;
	DotSphManager * dotSphereManager_;
	bool * dotsOn_;
	mutable int * indsOfDotsOn_;
	mutable bool indsCurrent_;
	int numDots_;
	int numDotsOn_;

};

std::ostream & operator << (std::ostream & os, DotsForAtom const & );

class Vertex_ths
{
public:
	Vertex_ths( 
		GraphToHoldScores * owner, 
		AtomPositions * xyz, 
		DotSphManager * dotSphereManager,
		int index,
		int num_states
	);
	~Vertex_ths();
	
	void setMover( Mover* mover );		
	void obtainAtomsFromMover();
	int getNumStates() const;
	int getNumAtoms() const;
	void getAtoms( std::vector< AtomDescr > & atoms ) const;
	//void getOverlappingBackgroundAtoms();
	//void dropMoverAtomsFromBackgroundAtomLists( Vertex_ths * other );
	bool anyMoverOverlap( Vertex_ths * other ) const;
	//void detectBackgroundAtomsInThreeWayOverlap();
	void detectThreeWayOverlap( Vertex_ths * other1, Vertex_ths * other_2 );
	void boot( int atom, DotsForAtom const & dotsToBoot, std::list< int > const & moversInHOO);

	std::list< DegreeTwoEdge_ths * >::iterator addEdge( DegreeTwoEdge_ths * );
	std::list< DegreeThreeEdge_ths * >::iterator addEdge( DegreeThreeEdge_ths * );
	std::list< DegreeFourEdge_ths * >::iterator addEdge( DegreeFourEdge_ths * );
		
	void dropEdge( std::list< DegreeTwoEdge_ths * >::iterator iter );
	void dropEdge( std::list< DegreeThreeEdge_ths * >::iterator iter );
	void dropEdge( std::list< DegreeFourEdge_ths * >::iterator iter );

	void informIncidentEdgesAboutBootedAtoms();
	void removeDotsScoredOnD2Edge( int atom, DotsForAtom const & dots_scored_on_edge );
	void finalizeScoreVectors();
	void score();

	//Methods for state disabling
	void enableAllStates();
	void disableState( int state );
	void setStateDisablingComplete();
	int convertEnabledStateToRegularState( int enabled_state ) const;

	//Methods for interaction graph initialization
	int getNumEnabledStates() const;
	float getScoreForEnabledState( int enabled_state ) const;
	Mover* getMover() const;
	bool getHasAnyHighOrderOverlap() const;
	
	std::vector< std::pair< AtomDescr, DotsForAtom * > >
	getAtomsInHighOrderOverlap() const;
	
	void getPenalties( std::vector< float > & penalties ) const;
	
	void assignState( int state );
	float getScoreForAssignedState() const;
	int getAssignedState() const;
private:
	//void noteBackgroundAtomInPsuedo3WOForAtom( int mover_atom, AtomDescr bgatom );
	//void noteBackgroundAtomInPsuedo4WOForAtom( int mover_atom, AtomDescr bgatom );

	GraphToHoldScores * owner_;
	AtomPositions * xyz_;
	DotSphManager * dotSphereManager_;

	int const index_;

	std::list< DegreeTwoEdge_ths* > deg2edges_;
	std::list< DegreeThreeEdge_ths* > deg3edges_;
	std::list< DegreeFourEdge_ths* > deg4edges_;
	
	Mover* mover_;
	//int num_atoms_;
	int num_states_;
	int num_incident_d2edges_;
	int num_incident_d3edges_;
	int num_incident_d4edges_;
	
	std::vector< AtomDescr > atoms_;
	std::vector< DotsForAtom > dotsForAtom_;
	std::vector< std::pair< AtomDescr, DotsForAtom * > > dotsToScoreOnAtoms_;

	bool any_high_order_overlap_;
	std::vector< bool > atom_in_high_order_overlap_;
	std::vector< DotsForAtom > dotsInHOO_;
	std::vector< std::pair< AtomDescr, DotsForAtom * > > dotsToScoreInHOO_;
	
	float* scores_;
	std::vector< float > penalties_;
	
	std::vector< bool > statesEnabled_;
	std::vector< int > enabledStates2OriginalStates_;
	int numEnabledStates_;
	
	int assignedState_; //testing purposes only

	Vertex_ths( Vertex_ths const & );
	Vertex_ths();
};

class DegreeTwoEdge_ths
{
public:
	DegreeTwoEdge_ths( 
		GraphToHoldScores * owner,
		AtomPositions * xyz,
		DotSphManager * dotSphereManager,
		int index1,
		int index2
	);
	~DegreeTwoEdge_ths();
	void setPosInOwnersEdgeList( std::list< DegreeTwoEdge_ths * >::iterator iter );

	bool sameEdge( int vertex1, int vertex2 ) const;
	Vertex_ths * getOtherNode( int vertex_index ) const;
	int getFirstNodeIndex() const;
	int getSecondNodeIndex() const;
	int getOtherNodeIndex( int vertex_index ) const;
	bool incidentUpon( int vertex_index ) const;
	void detectDotsToScoreOnEdge();
	DotsForAtom const & getDotsForAtom( int vertex, int atom );
	void noteBootedAtoms(
		int vertex_index,
		std::vector< bool  > const & booted_atoms,
		std::vector< DotsForAtom > const & dotsInHOO
	);
	void tellVerticesToNotScoreDotsScoredOnD2Edge() const;
	void dropDotsScoredOnD3Edge(
		int vertex,
		int atom,
		DotsForAtom const & dots );
	
	
	void finalizeScoreVectors();
	void score();
	
	//methods for interaction graph initialization
	float getScoreForEnabledStates( int enabled_state1, int enabled_state2 ) const;
	
	float getScoreForAssignedState() const;
	
private:
	int whichVertex( int index ) const;
	bool nothingToScore() const;
	bool nothingToScore( int vertex ) const; 
	
	GraphToHoldScores * owner_;
	AtomPositions * xyz_;
	DotSphManager * dotSphereManager_;

	int vertex_indices_[2];
	Vertex_ths* vertex_ptrs_[2];
	int num_states_[2];
	int num_atoms_[2];
	Mover* movers_[2];
	std::list< DegreeTwoEdge_ths* >::iterator pos_in_verts_edge_list_[2];
	std::list< DegreeTwoEdge_ths* >::iterator pos_in_owners_edge_list_;
	
	std::vector< AtomDescr > atoms_[ 2 ];
	std::vector< DotsForAtom > dotsForAtom_[ 2 ];
	std::vector< std::list< AtomDescr > > mover_atoms_atom_overlaps_[ 2 ];

	std::vector< std::pair< AtomDescr, DotsForAtom * > > dotsToScoreOnAtoms_[ 2 ];
	//std::vector< std::pair< AtomDescr, std::vector< AtomDescr > > > atoms_to_score_[ 2 ];
		
	float** scores_;
};

class DegreeThreeEdge_ths
{
public:
	DegreeThreeEdge_ths( 
		GraphToHoldScores * owner,
		AtomPositions * xyz,
		DotSphManager * dotSphereManager,
		int index1,
		int index2,
		int index3
	);
	~DegreeThreeEdge_ths();
	
	void setPosInOwnersEdgeList( std::list< DegreeThreeEdge_ths * >::iterator iter );
	bool sameEdge( int vertex1, int vertex2, int vertex3 ) const;
	int getFirstNodeIndex() const;
	int getSecondNodeIndex() const;
	int getThirdNodeIndex() const;	
	
	DotsForAtom const & 
	getDotsForAtom( int vertex, int atom );
	
	void
	setMoverAtomHas3Way(
		int vertex,
		int atom
	);
	
	void detectDotsToScoreOnEdge();
	void detectDotsDoublyCounted( DegreeThreeEdge_ths const * other ) const;

	void tellD2EdgesToNotScoreDotsScoredOnD3Edge() const;
	void dropDotsScoredOnD4Edge(
		int vertex_index,
		int atom,
		DotsForAtom const & dotsOnD4Edge
	);
	void noteBootedAtoms(
		int vertex_index,
		std::vector< bool  > const & booted_atoms,
		std::vector< DotsForAtom > const & dotsInHOO
	);	
	void finalizeScoreVectors();
	void score();

	//methods for interaction graph initialization
	float getScoreForEnabledStates( 
		int enabled_state1, 
		int enabled_state2, 
		int enabled_state3 ) const;

	float getScoreForAssignedState() const;

private:
	int whichVertex( int vertex ) const;
	
	void
	lookForSharedDots( 
		DegreeThreeEdge_ths const * other,
		int vertex_this,
		int vertex_other
	) const;

	bool nothingToScore() const;
	bool nothingToScore( int vertex ) const; 
	
	GraphToHoldScores * owner_;
	AtomPositions * xyz_;
	DotSphManager * dotSphereManager_;

	int vertex_indices_[3];
	Vertex_ths* vertex_ptrs_[3];
	int num_states_[3];
	int num_atoms_[3];
	Mover* movers_[3];
	DegreeTwoEdge_ths * d2e_[ 3 ]; // 0 == 0-1 edge, 1 == 0-2 edge, 3 == 1-2 edge
	
	std::list< DegreeThreeEdge_ths* >::iterator pos_in_verts_edge_list_[ 3 ];
	std::list< DegreeThreeEdge_ths* >::iterator pos_in_owners_edge_list_;
	
	std::vector< AtomDescr > atoms_[ 3 ];
	std::vector< bool > atomHasThreeWay_[ 3 ];
	std::vector< DotsForAtom > dotsForAtom_[ 3 ];
	std::vector< std::pair< AtomDescr, DotsForAtom * > > dotsToScoreOnAtoms_[ 3 ];

		
	float*** scores_table_;

};

class DegreeFourEdge_ths
{
public:
	DegreeFourEdge_ths( 
		GraphToHoldScores * owner, 
		AtomPositions * xyz, 
		DotSphManager * dotSphereManager, 
		int index1, 
		int index2, 
		int index3, 
		int index4
	);
	~DegreeFourEdge_ths();
	
	void setPosInOwnersEdgeList( std::list< DegreeFourEdge_ths * >::iterator iter );
	bool sameEdge( int vertex1, int vertex2, int vertex3, int vertex4 ) const;
	int getFirstNodeIndex() const;
	int getSecondNodeIndex() const;
	int getThirdNodeIndex() const;
	int getFourthNodeIndex() const;
	
	void
	setAtomHas4WayOverlap(
		int vertex,
		int atom
	);

	void detectDotsToScoreOnEdge();
	void detectDotsDoublyCounted( DegreeFourEdge_ths const * other ) const;
	void tellD3EdgesToNotScoreDotsScoredOnD4Edge() const;

	void noteBootedAtoms(
		int vertex_index,
		std::vector< bool  > const & booted_atoms,
		std::vector< DotsForAtom > const & dotsInHOO
	);
	void finalizeScoreVectors();
	void score();

	//methods for interaction graph initialization
	float getScoreForEnabledStates( 
		int enabled_state1, 
		int enabled_state2, 
		int enabled_state3,
		int enabled_state4 ) const;

	float getScoreForAssignedState() const;

private:
	int whichVertex( int vertex ) const;

	void
	lookForSharedDots( 
		DegreeFourEdge_ths const * other,
		int vertex_this,
		int vertex_other
	) const;
	
	bool nothingToScore() const;
	bool nothingToScore( int vertex ) const; 
	
	GraphToHoldScores * owner_;
	AtomPositions * xyz_;
	DotSphManager * dotSphereManager_;

	int vertex_indices_[4];
	Vertex_ths* vertex_ptrs_[4];
	int num_states_[4];
	int num_atoms_[4];
	Mover* movers_[4];
	//DegreeTwoEdge_ths * d2e_[ 4 ][ 4 ];
	
	std::list< DegreeFourEdge_ths* >::iterator pos_in_verts_edge_list_[ 4 ];
	std::list< DegreeFourEdge_ths* >::iterator pos_in_owners_edge_list_;
	
	std::vector< AtomDescr > atoms_[ 4 ];
	std::vector< DotsForAtom > dotsForAtom_[ 4 ];
	std::vector< bool > atomHasFourWay_[ 4 ];
	std::vector< std::pair< AtomDescr, DotsForAtom * > > dotsToScoreOnAtoms_[ 4 ];

	std::vector< std::pair< AtomDescr, std::vector< AtomDescr > > >  atoms_to_score_[ 4 ];
		
	float**** scores_table_;

};

class GraphToHoldScores
{
public:
	GraphToHoldScores( 
		AtomPositions * xyz,
		DotSphManager * dotSphereManager, 
		std::vector< int > const & num_states,
		std::vector< Mover * > const & movers
	);
	~GraphToHoldScores();

	int getNumNodes() const;
	Vertex_ths * getVertexPtr( int vertex_index ) const;
	DegreeTwoEdge_ths* getDegree2Edge( int fn, int sn );
	DegreeThreeEdge_ths* getDegree3Edge( int fn, int sn, int tn );
	DegreeFourEdge_ths* getDegree4Edge( int fn, int sn, int tn, int fthn );
	
	bool degree2EdgeExists( int fn, int sn ) const;
	bool degree3EdgeExists( int fn, int sn, int tn ) const;
	bool degree4EdgeExists( int fn, int sn, int tn, int fthn ) const;
	
	void forceClique( std::list< int > const & movers_in_clique );
	void dropEdge( std::list< DegreeTwoEdge_ths * >::iterator iter );
	void dropEdge( std::list< DegreeThreeEdge_ths * >::iterator iter );
	void dropEdge( std::list< DegreeFourEdge_ths * >::iterator iter );

	void getPenalties( std::vector< std::vector< float > > & penalties, std::list< float > & allPenalties );
	
	//Methods for state disabling
	void setAllStatesEnabled();
	void disableStateOnNode( int vertex, int state );
	void setStateDisablingCompleteForNow();
	bool anyNodeWithAllStatesDisabled() const;
	int convertEnabledState2OriginalStateEnumerationOnNode( int vertex, int enabled_state ) const;

	//Methods for initializing interaction graph
	AtomPositions * getAtomPositionsPointer() const;
	int getNumStatesForNode(int vertex_index ) const;
	int getNumEnabledStatesForNode( int vertex_index ) const;
	float getNodeScoreForState( int vertex_id, int state) const;
	bool getNodeHasAnyHighOrderOverlap( int vertex_id ) const;
	std::vector< std::pair< AtomDescr, DotsForAtom * > > getAtomsInHighOrderOverlapForNode( int vertex_id ) const;
	Mover* getMoverForNode( int vertex_id ) const;
	
	void setD2EIteratorAtBegining();
	void incrementD2EIterator();
	bool getD2EIteratorAtEnd() const;
	int getFirstIndexFocusedD2E() const;
	int getSecondIndexFocusedD2E() const;
	float getScoreForFocusedD2E( int state1, int state2 ) const;
	
	void setD3EIteratorAtBegining();
	void incrementD3EIterator();
	bool getD3EIteratorAtEnd() const;

	void setD4EIteratorAtBegining();
	void incrementD4EIterator();
	bool getD4EIteratorAtEnd() const;

	int getFirstIndexFocusedD3E() const;
	int getSecondIndexFocusedD3E() const;
	int getThirdIndexFocusedD3E() const;
	float getScoreForFocusedD3E( int state1, int state2, int state3) const;

	int getFirstIndexFocusedD4E() const;
	int getSecondIndexFocusedD4E() const;
	int getThirdIndexFocusedD4E() const;
	int getFourthIndexFocusedD4E() const;
	float getScoreForFocusedD4E( int state1, int state2, int state3, int state4) const;
	
	//testing purposes only
	float getScoreForStateAssignment( std::vector< int > const & );
	
private:
	void allocateConnectivityTable();
	void instantiateVertices( 
		std::vector< int > const & num_states,
		std::vector< Mover * > const & movers
	);
	void obtainAtomsFromMovers();
	void addDegreeTwoEdges();
	void addDegreeThreeEdges();
	void detectFourWayInteractions();
	void detectFiveWayInteractions();
	void cascadeScoringInfo();
	void finalizeScoreVectors();
	void score();
	
	void nodeInRange( int node ) const;
	DegreeTwoEdge_ths* findDeg2Edge( int fn, int sn ) const;
	DegreeTwoEdge_ths* addDegree2Edge( int node1, int node2 );
	DegreeThreeEdge_ths* findDeg3Edge( int fn, int sn, int tn ) const;
	DegreeThreeEdge_ths* addDegree3Edge( int fn, int sn, int tn );
	DegreeFourEdge_ths* findDeg4Edge( int fn, int sn, int tn, int fthn ) const;
	DegreeFourEdge_ths* addDegree4Edge( int fn, int sn, int tn, int fthn );

	AtomPositions * xyz_;
	DotSphManager * dotSphereManager_;
	int num_vertices_;
	int num_deg2edges_;
	int num_deg3edges_;
	int num_deg4edges_;
	
	std::vector< Vertex_ths * > vertices_;
	std::list< DegreeTwoEdge_ths * > deg2edges_;
	std::list< DegreeThreeEdge_ths * > deg3edges_;
	std::list< DegreeFourEdge_ths * > deg4edges_;
	
	std::list< DegreeTwoEdge_ths * >::iterator d2eiter_;
	std::list< DegreeThreeEdge_ths * >::iterator d3eiter_;
	std::list< DegreeFourEdge_ths * >::iterator d4eiter_;
	
	bool ** connectivity_;
	
	//no default or copy constructors
	GraphToHoldScores();
	GraphToHoldScores( GraphToHoldScores const & );
};

#endif
