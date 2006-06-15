#ifndef GRAPH_TO_HOLD_SCORES_H
#define GRAPH_TO_HOLD_SCORES_H

#include <list>
#include <vector>

class Mover;
class AtomPositions;
#include "AtomDescr.h"

void sort_three(int, int, int, int &, int &, int & );

template < class T >
void set_intersect_apl( std::list< T > const & list1, std::list< T > const & list2, std::list< T > & intersection );

template < class T >
void set_difference_apl( std::list< T > & list1, std::list< T > const & list2 );

template < class T >
void set_difference_apl( std::list< T > & list1, std::vector< T > const & list2 );

class GraphToHoldScores;
class Vertex_ths;
class DegreeTwoEdge_ths;
class DegreeThreeEdge_ths;

class Vertex_ths
{
public:
	Vertex_ths( GraphToHoldScores * owner, AtomPositions * xyz, int index, int num_states );
	~Vertex_ths();
	
	void setMover( Mover* mover );		
	void obtainAtomsFromMover();
	int getNumStates() const;
	int getNumAtoms() const;
	void getAtoms( std::vector< AtomDescr > & atoms ) const;
	void getOverlappingBackgroundAtoms();
	void dropMoverAtomsFromBackgroundAtomLists( Vertex_ths * other );
	bool anyMoverOverlap( Vertex_ths * other ) const;
	void detectBackgroundAtomsInThreeWayOverlap();
	void detectThreeWayOverlap( Vertex_ths * other1, Vertex_ths * other_2 );
	void boot( int atom );

	std::list< DegreeTwoEdge_ths * >::iterator addEdge( DegreeTwoEdge_ths * );
	std::list< DegreeThreeEdge_ths * >::iterator addEdge( DegreeThreeEdge_ths * );
	void dropEdge( std::list< DegreeTwoEdge_ths * >::iterator iter );
	void dropEdge( std::list< DegreeThreeEdge_ths * >::iterator iter );
	std::list< AtomDescr > getAtomsOverlappingBothMoverAndBGAtom( AtomDescr mover_atom, AtomDescr bg_atom ) const;

	std::list< AtomDescr > const & getAllBumpingBackgroundAtoms( int atom ) const;
	std::list< AtomDescr > const & getBumpingScoringBackgroundAtoms( int atom ) const;
	void informIncidentEdgesAboutBootedAtoms();
	void removeBGAtomsForAtom( int atom, std::list< AtomDescr > const & bg_atoms_scored_on_edge );
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
	std::list< AtomDescr > getAtomsInHighOrderOverlap() const;
	void getPenalties( std::vector< float > & penalties ) const;
	
	void assignState( int state );
	float getScoreForAssignedState() const;
	int getAssignedState() const;
private:
	void noteBackgroundAtomInPsuedo3WOForAtom( int mover_atom, AtomDescr bgatom );

	GraphToHoldScores * owner_;
	AtomPositions * xyz_;

	int const index_;

	std::list< DegreeTwoEdge_ths* > deg2edges_;
	std::list< DegreeThreeEdge_ths* > deg3edges_;
	
	Mover* mover_;
	//int num_atoms_;
	int num_states_;
	int num_incident_d2edges_;
	int num_incident_d3edges_;
	
	std::vector< AtomDescr > atoms_;
	std::vector< std::list< AtomDescr > > bg_ats_overlapping_;
	std::vector< std::list< AtomDescr > > bg_ats_scoring_overlapping_;
	
	std::vector< std::pair< AtomDescr, std::vector< AtomDescr > > > atoms_to_score_as_vertex_;

	bool any_high_order_overlap_;
	std::vector< bool > atom_in_high_order_overlap_;
	
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
	DegreeTwoEdge_ths( GraphToHoldScores * owner, AtomPositions * xyz, int index1, int index2);
	~DegreeTwoEdge_ths();
	void setPosInOwnersEdgeList( std::list< DegreeTwoEdge_ths * >::iterator iter );

	bool sameEdge( int vertex1, int vertex2 ) const;
	Vertex_ths * getOtherNode( int vertex_index ) const;
	int getFirstNodeIndex() const;
	int getSecondNodeIndex() const;
	int getOtherNodeIndex( int vertex_index ) const;
	bool incidentUpon( int vertex_index ) const;
	void detectBackgroundAtomsToScoreForEdge();
	void addOverlappingBackgroundAtomsForPairToD3E(
		int atom_on_vertex1,
		int atom_on_vertex2,
		AtomDescr atom_on_vertex3,
		DegreeThreeEdge_ths * d3e
	);
	std::list< AtomDescr > getScoringBackgroundAtomsOverlappingAtom( int vertex_index, int atom );
	bool bgAtomOverlapsMovingAtomOnEdge( int vertex_index, int atom, AtomDescr bgatom ) const;
	bool atomHasOverlap( int vertex_index, int atom ) const;
	void noteBootedAtoms( int vertex_index, std::vector< bool  > const & booted_atoms );
	void tellVerticesToNotScoreBGAtomsScoredOnD2Edge() const;
	void dropBGAtomsOverlappingAtom(
		int vertex_index,
		int atom,
		std::list< AtomDescr > const & bgats
	);
	void dropMoverAtomsOverlappingAtom(
		int vertex_index,
		int atom,
		std::list< AtomDescr > const & moverats
	);
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

	int vertex_indices_[2];
	Vertex_ths* vertex_ptrs_[2];
	int num_states_[2];
	int num_atoms_[2];
	Mover* movers_[2];
	std::list< DegreeTwoEdge_ths* >::iterator pos_in_verts_edge_list_[2];
	std::list< DegreeTwoEdge_ths* >::iterator pos_in_owners_edge_list_;
	
	std::vector< AtomDescr > atoms_[ 2 ];
	std::vector< std::list< AtomDescr > > bg_ats_scoring_for_atom_[ 2 ];
	std::vector< std::list< AtomDescr > > bg_ats_all_for_atom_[ 2 ];
	std::vector< std::list< AtomDescr > > mover_atoms_atom_overlaps_[ 2 ];
	std::vector< std::list< std::list< AtomDescr > > > bg_atoms_overlapping_pair_[ 2 ];

	std::vector< std::pair< AtomDescr, std::vector< AtomDescr > > > atoms_to_score_[ 2 ];
		
	float** scores_;
};

class DegreeThreeEdge_ths
{
public:
	DegreeThreeEdge_ths( GraphToHoldScores * owner, AtomPositions * xyz, int index1, int index2, int index3);
	~DegreeThreeEdge_ths();
	
	void setPosInOwnersEdgeList( std::list< DegreeThreeEdge_ths * >::iterator iter );
	bool sameEdge( int vertex1, int vertex2, int vertex3 ) const;
	int getFirstNodeIndex() const;
	int getSecondNodeIndex() const;
	int getThirdNodeIndex() const;	
	void addBackgroundAtomToScoreForMoverAtom(
		int vertex,
		int atom,
		AtomDescr atom_overlapping
	);
	
	void
	addMoverAtomsToScoreForAtom(
		int vertex,
		int atom,
		AtomDescr atom_overlapping
	);

	bool
	sharesTwoVertices( DegreeThreeEdge_ths const * other ) const;

	void detect_shared_atom_pair( DegreeThreeEdge_ths * other ) const;
	void detectBackgroundAtomScoredTwiceForAtom( DegreeThreeEdge_ths const * other ) const;
	void tellD2EdgesToNotScoreBGandMoverAtomsScoredOnD3Edge() const;
	void noteBootedAtoms( int vertex_index, std::vector< bool  > const & booted_atoms );
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
	
	void compareSharedAtomPairs(
		DegreeThreeEdge_ths const  *  other,
		int vertex_this,
		int vertex_other
	) const;
	
	void
	compareSharedBackgroundAtoms(
		DegreeThreeEdge_ths const  *  other,
		int vertex_this,
		int vertex_other
	) const;

	bool nothingToScore() const;
	bool nothingToScore( int vertex ) const; 
	
	GraphToHoldScores * owner_;
	AtomPositions * xyz_;

	int vertex_indices_[3];
	Vertex_ths* vertex_ptrs_[3];
	int num_states_[3];
	int num_atoms_[3];
	Mover* movers_[3];
	DegreeTwoEdge_ths * d2e_[ 3 ]; // 0 == 0-1 edge, 1 == 0-2 edge, 3 == 1-2 edge
	
	std::list< DegreeThreeEdge_ths* >::iterator pos_in_verts_edge_list_[ 3 ];
	std::list< DegreeThreeEdge_ths* >::iterator pos_in_owners_edge_list_;
	
	std::vector< AtomDescr > atoms_[ 3 ];
	std::vector< std::list< AtomDescr > > bg_ats_for_atom_[ 3 ];
	std::vector< std::list< AtomDescr > > mover_atoms_atom_overlaps_[ 3 ];

	std::vector< std::pair< AtomDescr, std::vector< AtomDescr > > >  atoms_to_score_[ 3 ];
		
	float*** scores_table_;

};

class GraphToHoldScores
{
public:
	GraphToHoldScores( 
		AtomPositions * xyz,
		std::vector< int > const & num_states,
		std::vector< Mover * > const & movers
	);
	~GraphToHoldScores();

	int getNumNodes() const;
	Vertex_ths * getVertexPtr( int vertex_index ) const;
	DegreeTwoEdge_ths* getDegree2Edge( int fn, int sn );
	DegreeThreeEdge_ths* getDegree3Edge( int fn, int sn, int tn );
	void forceClique( std::list< int > const & movers_in_clique );
	void dropEdge( std::list< DegreeTwoEdge_ths * >::iterator iter );
	void dropEdge( std::list< DegreeThreeEdge_ths * >::iterator iter );

	void getPenalties( std::vector< std::vector< float > > & penalties, std::list< float > & allPenalties );
	
	//Methods for state disabling
	void setAllStatesEnabled();
	void disableStateOnNode( int vertex, int state );
	void setStateDisablingCompleteForNow();
	int convertEnabledState2OriginalStateEnumerationOnNode( int vertex, int enabled_state ) const;

	//Methods for initializing interaction graph
	AtomPositions * getAtomPositionsPointer() const;
	int getNumStatesForNode(int vertex_index ) const;
	int getNumEnabledStatesForNode( int vertex_index ) const;
	float getNodeScoreForState( int vertex_id, int state) const;
	bool getNodeHasAnyHighOrderOverlap( int vertex_id ) const;
	std::list< AtomDescr > getAtomsInHighOrderOverlapForNode( int vertex_id ) const;
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
	
	int getFirstIndexFocusedD3E() const;
	int getSecondIndexFocusedD3E() const;
	int getThirdIndexFocusedD3E() const;
	float getScoreForFocusedD3E( int state1, int state2, int state3) const;
	
	//testing purposes only
	float getScoreForStateAssignment( std::vector< int > const & );
	
private:
	void allocateConnectivityTable();
	void instantiateVertices( 
		std::vector< int > const & num_states,
		std::vector< Mover * > const & movers
	);
	void obtainAtomsFromMovers();
	void getOverlappingBackgroundAtoms();
	void dropMoverAtomsFromBackgroundAtomLists();
	void addDegreeTwoEdges();
	void promoteBackgroundAtomsToDegreeTwoEdges();
	void addDegreeThreeEdges();
	void detectFourWayInteractions();
	void cascadeScoringInfo();
	void finalizeScoreVectors();
	void score();
	
	void nodeInRange( int node ) const;
	DegreeTwoEdge_ths* findDeg2Edge( int fn, int sn ) const;
	DegreeTwoEdge_ths* addDegree2Edge( int node1, int node2 );
	DegreeThreeEdge_ths* findDeg3Edge( int fn, int sn, int tn ) const;
	DegreeThreeEdge_ths* addDegree3Edge( int fn, int sn, int tn );

	AtomPositions * xyz_;
	int num_vertices_;
	int num_deg2edges_;
	int num_deg3edges_;
	
	std::vector< Vertex_ths * > vertices_;
	std::list< DegreeTwoEdge_ths * > deg2edges_;
	std::list< DegreeThreeEdge_ths * > deg3edges_;
	
	std::list< DegreeTwoEdge_ths * >::iterator d2eiter_;
	std::list< DegreeThreeEdge_ths * >::iterator d3eiter_;
	
	bool ** connectivity_;
	
	//no default or copy constructors
	GraphToHoldScores();
	GraphToHoldScores( GraphToHoldScores const & );
};

#endif
