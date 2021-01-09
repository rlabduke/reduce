// name: MoveableNode.h
// author: Andrew Leaver-Fay
// date written: 4/15/06
// purpose: Interface for Interaction Graph Classes

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************


#ifndef  MoveableNode_h
#define  MoveableNode_h

#include <iostream>

#include "AtomDescr.h"
#include "GraphToHoldScores.h"

class AtomPositions;
class Mover;

#include <vector>
#include <list>

class OptimizedNodeStateVector;	//Keeps optimal node state for node eliminated through tree reduction
class OptimizedNodeStateMatrix; //Keeps optimal node state for node eliminated through cycle reduction
class NodeScoreVector;			//The one dimensional array that keeps track of the "internal" score for a node in each of its states
class EdgeScoreMatrix;			//The two dimensional array that keeps track of the interaction dependent score for a pair of nodes
class NodeAndStatePair;			//A glorified struct that holds a node's index and it's state.
class NetworkStateVector;		//The planned interface between the Richardson's code and the NodeAndEdgeManager for reporting the optimal network state
class MoveableNode;				//Class representing nodes in the interaction graph
class EdgeBetweenMoveableNodes;	//Class representing edges in the interaction graph
class Degree3Hyperedge;
class Degree4Hyperedge;
class NodeAndEdgeManager;		//Singleton entity meant to serve as the main interface for this algorithm and code that uses it
class QueueOfToBeContractedNodes;//Singleton queue of nodes that have been marked as ready to be contracted
class DummyNetworkDescriptorClass; //Testing class used to describe some arbitrary network to the NodeAndEdgeManager

typedef int NodeState;

//Nodes are indexed from 0 through all parts of the interface.
//NodeStates are indexed from 0 as well.
//In accessing information about edges through the indices of their nodes, the order is smaller node index first, larger node index second
//		- e.g. if I am asking for the score of an edge between nodes indexed by 5 and 10, where the NodeState of node 5 is 3 and the NodeState
//		of node 10 is 1 then the call is getEdgeScore(3,1).
//A sentinal value of -1 is used for all cases of uninitialized values in indexing nodes or node states.

class OptimizedNodeStateVector
{
public:
	OptimizedNodeStateVector();
	OptimizedNodeStateVector(NodeState maxNodeStates);
	OptimizedNodeStateVector(int depNodeIndex, NodeState maxNodeStates);
	~OptimizedNodeStateVector();
	OptimizedNodeStateVector(const OptimizedNodeStateVector& rhs);
	void setMaxNodeStates(NodeState maxNodeStates);
	NodeState getOptimumNodeState(NodeState dependentNodeState);
	void setOptimizedNodeState(NodeState dependentNodeState, NodeState ownerNodeState);
private:
	std::vector<NodeState >	_theNSVector;
	int					_dependentNodeIndex;
};

class OptimizedNodeStateMatrix
{
public:
	OptimizedNodeStateMatrix();
	OptimizedNodeStateMatrix(NodeState firstNodeMaxStates, NodeState secondNodeMaxStates);
	OptimizedNodeStateMatrix(int firstNodeIndex, int secondNodeIndex, NodeState firstNodeMaxStates, NodeState secondNodeMaxStates);
	~OptimizedNodeStateMatrix();
	OptimizedNodeStateMatrix(const OptimizedNodeStateMatrix& rhs); //shallow copy
	void DeepCopy(const OptimizedNodeStateMatrix& toBeCopied);

	NodeState getOptimumNodeState(NodeState firstNodeState, NodeState secondNodeState);
	void setOptimumNodeState(NodeState firstNodeState, NodeState secondNodeState, NodeState ownerNodeState);

private:
	void				deallocateMatrix();

	NodeState**			_theNSMatrix;
	int					_1stNodeIndex;
	int					_2ndNodeIndex;
	int					_1stNodeMaxStates;
	int					_2ndNodeMaxStates;

};

class OptimizedNodeStateMatrix3
{
public:
	OptimizedNodeStateMatrix3();
	OptimizedNodeStateMatrix3(NodeState firstNodeMaxStates, NodeState secondNodeMaxStates, NodeState thirdNodeMaxStates);
	OptimizedNodeStateMatrix3(int firstNodeIndex, int secondNodeIndex, int thirdNodeIndex, NodeState firstNodeMaxStates, NodeState secondNodeMaxStates, NodeState thirdNodeMaxStates);
	~OptimizedNodeStateMatrix3();

	NodeState getOptimumNodeState(NodeState firstNodeState, NodeState secondNodeState, NodeState thirdNodeState);
	void setOptimumNodeState(NodeState firstNodeState, NodeState secondNodeState, NodeState thirdNodeState, NodeState ownerNodeState);

private:
	void				deallocateMatrix();

	NodeState***		_theNSMatrix;
	int					_1stNodeIndex;
	int					_2ndNodeIndex;
	int					_3rdNodeIndex;
	int					_1stNodeMaxStates;
	int					_2ndNodeMaxStates;
	int					_3rdNodeMaxStates;

};


class NodeScoreVector
{
public:
	NodeScoreVector();
	NodeScoreVector(NodeState ownerMaxStates);
	~NodeScoreVector();
	NodeScoreVector(const NodeScoreVector& rhs);

	double getNodeScore(NodeState ownersState);
	void   setNodeScore(NodeState ownersState, double theScore);
	void   addToNodeScore(NodeState ownersState, double theScore);
private:
	std::vector<double >		_theNScrVect;
};

class EdgeScoreMatrix
{
public:
	EdgeScoreMatrix();
	EdgeScoreMatrix(NodeState firstNodeMaxStates, NodeState secondNodeMaxStates);
	~EdgeScoreMatrix();
	EdgeScoreMatrix(const EdgeScoreMatrix& rhs);			//shallow copy - shouldn't be called in any case
	void DeepCopy(const EdgeScoreMatrix& toBeCopied);		//implemented, but never used

	double getEdgeScore(NodeState firstNodeState, NodeState secondNodeState);
	void   setEdgeScore(NodeState firstNodeState, NodeState secondNodeState, double theScore);
	void   addEdgeScore(NodeState firstNodeState, NodeState secondNodeState, double theScore);

private:
	void				deallocateMatrix();
	double**			_theEScrMatrix;
	int					_1stNodeIndex;
	int					_2ndNodeIndex;
	NodeState			_1stNodeMaxStates;
	NodeState			_2ndNodeMaxStates;
};

class NodeAndStatePair
{
public:
	NodeAndStatePair();
	NodeAndStatePair(int index, NodeState state);
	~NodeAndStatePair();
	NodeAndStatePair(const NodeAndStatePair& rhs);
	int		  getNodeIndex();
	NodeState getNodeState();
	void	  setNodeIndex(int);
	void	  setNodeState(NodeState);

private:
	int		  _nodeIndex;
	NodeState _nodeState;
};

class MoveableNode
{
public:
	MoveableNode();
	~MoveableNode();
	//MoveableNode(const MoveableNode& rhs);		//unimplemented since shallow copy is sufficient - should not be called in any case
	MoveableNode(NodeState maxStates);				

	void				setNodeStateScore(NodeState theNodeState, double theScore);
	double				getNodeStateScore(NodeState theNodeState); 
	void 				eliminate();
	void 				beNotifiedDependencyEliminated(int i);
	NodeState			getNumberOfPossibleStates();
	void 				setNumberOfPossibleStates(NodeState maxStates);
	int  				getNodeIndex();
	void				setNodeIndex(int i);
	bool 				getIsEliminated();
	void 				setIsEliminated(bool eliminationState);
	void 				combineNodeScoreWithOptEdgeAndNodeScore(NodeScoreVector& optNodeAndEdgeScoreVect);
	void 				addEdge(EdgeBetweenMoveableNodes* newEdge);
	//double				getOptimalNetworkScore();			//To be called iff size(_edgeList) == 0
	NodeState			getNodeStateInOptimalNetworkConf(std::vector<NodeState>& theOptNetworkState);
	int					getNumberOfEdges();
	std::list<EdgeBetweenMoveableNodes* >*  getEdgeList();

	bool getHasAnyDegree3Hyperedges() const;
	bool hasD3Edge( int fn, int sn, int tn );
	Degree3Hyperedge * getD3Edge( int fn, int sn, int tn );
	bool getHasAnyDegree4Hyperedges() const;
	bool getHasAnyHighOrderOverlap() const { return _hasAnyHighOrderOverlap;}
	std::list< Degree3Hyperedge * > &  getDegree3HyperedgeList();
	std::list< Degree4Hyperedge * > &  getDegree4HyperedgeList();
	
	void addDegree3Hyperedge( Degree3Hyperedge * edge );
	void addDegree4Hyperedge( Degree4Hyperedge * edge );
	void deleteDegreeThreeHyperedge( Degree3Hyperedge * edge );
	void deleteDegreeFourHyperedge( Degree4Hyperedge * edge );
	void setMover( Mover* mover ) {_mover = mover;}
	Mover* getMover() {return _mover;}
	
	void setAtomsInHighOrderOverlap( 
		std::vector< std::pair< AtomDescr, DotsForAtom * > > const & atsInHOO
	)
	{ 	
		_hasAnyHighOrderOverlap = true; 
		_atomsInHighOrderOverlap = atsInHOO;
	}
   
   std::vector< std::pair< AtomDescr, DotsForAtom * > > getAtomsInHighOrderOverlap() {return _atomsInHighOrderOverlap;}
private:
	void 				eliminateThroughTreeReduction();
	void 				eliminateThroughCycleReduction();
	void				eliminateThroughS3Reduction();
	void				eliminateSingleton();
	void				eliminateOneStateVertex();
	
	int									_index;
	NodeState							_maxNodeStates;
	std::list<EdgeBetweenMoveableNodes* >	_edgeList;
	NodeScoreVector*					_theNodeScoreVect;	
	bool								_isEliminated;
	int									_eliminationOrder;	//0 = not eliminated, 1 = s_1 elimination, 2 = s_2 elimination, 3 = singleton elimination, 4 = s_3 elimination
	OptimizedNodeStateVector*			_optimalStateVector;
	NodeState _optimalState; //for singleton elimination
	double _optimalScore;
	int									_nodeIndexDetStateFromTreeReduction;
	
	OptimizedNodeStateMatrix*			_optimalStateMatrix;
	int									_firstNodeIndexDetStateFromCycleReduction;
	int									_secondNodeIndexDetStateFromCycleReduction;
	
	OptimizedNodeStateMatrix3* _optimalStateMatrix3;
	int _nodesDeterminingOptimalState[ 3 ];

	bool	_hasAnyDegree3Hyperedges;
	std::list< Degree3Hyperedge * > _degree3hyperedges;

	bool	_hasAnyDegree4Hyperedges;
	std::list< Degree4Hyperedge * > _degree4hyperedges;

	bool _hasAnyHighOrderOverlap;
	std::vector< std::pair< AtomDescr, DotsForAtom * > > _atomsInHighOrderOverlap;
	Mover* _mover;
};

class EdgeBetweenMoveableNodes 
{
public:
	EdgeBetweenMoveableNodes(NodeState totStatesFirstNode, NodeState totStatesSecondNode);
	~EdgeBetweenMoveableNodes();
	//EdgeBetweenMoveableNodes(const EdgeBetweenMoveableNodes& rhs);	//Shallow copy is desired if any copy is called at all
	
	void   setEdgeScore(NodeState firstNodeState, NodeState secondNodeState, double theScore);
	void   addEdgeScore( NodeState firstNodeState, NodeState secondNodeState, double theScore);
	double getEdgeScore(NodeState firstNodeState, NodeState secondNodeState);
	double getEdgeScore(int Node1sIndex, NodeState Node1sState, int Node2sIndex, NodeState Node2sState);
	void   setFirstNodePtr(MoveableNode* firstNode);
	void   setSecondNodePtr(MoveableNode* secondNode);
	void   setFirstNodeIndex(int firstNodeIndex);
	void   setSecondNodeIndex(int secondNodeIndex);
	int    getFirstNodeIndex();
	int    getSecondNodeIndex();
	bool   getIsEliminated();
	void   setIsEliminated(bool isEliminated = true);
	int	   getIndexOfOtherNode(int indexOfOneNode);
	MoveableNode* getPointerToOtherNode(int indexOfOneNode);
	
private:
	EdgeScoreMatrix* _theEdgeScoreMatrix;
	bool			 _isEliminated;
	int				 _firstNodeIndex;
	int			     _secondNodeIndex;
	MoveableNode*    _firstNode;
	MoveableNode*    _secondNode;

};

class Degree3Hyperedge
{
public:
	
	Degree3Hyperedge( int vertex1, int vertex2, int vertex3);
	~Degree3Hyperedge();
		
	bool sameEdge(int fn, int sn, int tn ) const;
	int getVertexIndex( int ind );
	void setScore( int state1, int state2, int state3, float energy);
	void setVertexOrder( int vert1, int vert2, int vert3); 
	void setNaturalVertexOrder() const;
	float getScore( int state1, int state2, int state3) const;
	float getScoreGivenSetVertexOrdering( int state_vert1, int state_vert2, int state_vert3) const;
	void addToScoreGivenSetVertexOrdering( int state_vert1, int state_vert2, int state_vert3, float score); 
private:
	int getIndex( int state1, int state2, int state3) const;

	int _nodeIndices[ 3 ];
	MoveableNode * _moveableNode[ 3 ];
	int _numStates[ 3 ];
	int _total_state_combos;
	mutable int _indexMultipliers [ 3 ];
	mutable bool _inNaturalVertexOrder;
	float * _scores;

	Degree3Hyperedge();
	Degree3Hyperedge( Degree3Hyperedge const & rhs );

};

class Degree4Hyperedge
{
public:
	
	Degree4Hyperedge( int vertex1, int vertex2, int vertex3,int vertex4);
	~Degree4Hyperedge();
		
	int getVertexIndex( int ind );
	void setScore( int state1, int state2, int state3, int state4, float energy);
	void setVertexOrder( int vert1, int vert2, int vert3, int vert4); 
	void setNaturalVertexOrder() const;
	float getScore( int state1, int state2, int state3, int state4) const;
	float getScoreGivenSetVertexOrdering( int state_vert1, int state_vert2, int state_vert3, int state_vert4) const;

private:
	int getIndex( int state1, int state2, int state3, int state4) const;

	int _nodeIndices[ 4 ];
	MoveableNode * _moveableNode[ 4 ];
	int _numStates[ 4 ];
	int _total_state_combos;
	mutable int _indexMultipliers [ 4 ];
	mutable bool _inNaturalVertexOrder;
	float * _scores;

	Degree4Hyperedge();
	Degree4Hyperedge( Degree4Hyperedge const & rhs );

};



//Initialize the NodeAndEdgeManager with a clique
//Tell it to determine the optimal network configuration
//then ask it either for the optimalNetworkConfiguration or the ScoreOfOptimalNetworkConfiguration or both.
//Once computeOptimalNetworkConfiguration is called, the original interaction graph is changed - both in scores and connectivity

class NodeAndEdgeManager
{
private:
	NodeAndEdgeManager();				//Singleton
	static bool							_instanceFlag;
	static NodeAndEdgeManager*			_theNaEManager;

	void eliminateQueuedNodes();
	void initiateTreeReduction();		//Trees are reduced, the disappearing nodes and edges in a tree being contracted are said to be eliminated
	void initiateCycleReduction();	
	void initiateSafeS3Reduction();
	void initiateUnsafeS3Reduction(); //better than brute force
	
	bool matchesAP86_CPrime( int node, int neighbors[3] );
	bool matchesAP86_CDoublePrime( int node, int neighbors[3] );
	
	bool bruteForceIrreducibleSubgraph();
	void computeOptimalNetworkConfigurationAfterOptimization();
	void writeOutOptimalNetworkState();
	
	std::vector<MoveableNode* >				_theMNVector;
	std::list<EdgeBetweenMoveableNodes* >	_theEBMNList;
	int									_numNodes;
	int									_numEdges;
	
	std::list<int> _PartialOrderStateDeterminationStack; //Evaluate Final Network configuration in reverse order of addition
	// this exists as a partial order of node dependencies - the state of the top node index in the stack
	// depended only on the states of those nodes previously removed from the stack if any.
	// Finding the order of all members of a stack destroys it - I only need to go through it once
	// to read out the partial order - but I may as well use a list to represent the stack so I can avoid
	// destorying it if maybe I want to use it again.

	int _numUnEliminatedNodes;	//Initialized to _numNodes
	bool _cycleReductionBegun;	//Don't add nodes to the cycleReduction queue until all initial tree reduction has been done - not because this is necessary, but because it is expedient
	bool _safeS3ReductionBegun;
	bool _unsafeS3ReductionBegun;
	bool _optimalSolutionFound;  //Keep track of whether work needs to be done first to report the optimal network state
	std::vector<NodeState > _optimalNetworkStateVector;  //Compute this once and keep it if needed more than once
	double _optimalNetworkScore;  

	int _indFirstNodeLastEdgeLookup;  //Used by the cycle reduction algorithm.  I currently do not have a way to get the edge between two nodes in constant time - this saves previous work for reuse.
	int _indSecondNodeLastEdgeLookup; 
	EdgeBetweenMoveableNodes* _ptrToEdgeOfLastEdgeLookup;
	AtomPositions * _xyz;
	
	bool _timeLimitExists;
	double _timeLimit;
	
	bool** _connectivity;
	std::vector< int > _vertexDegrees;
	double _effortCounter;
public:
	~NodeAndEdgeManager();
	static NodeAndEdgeManager* getInstance();
	static bool betterThan(double first, double second);
	//void NotifyNodeOfDependencyElimination(int eliminatedNodeIndex, int dependentNodeIndex);
	void BeNotifiedOfEliminatedNode(int indexContractedNode);
	bool computeOptimalNetworkConfiguration();
	double getScoreOfOptimalNetworkConfiguration();
	void InitializeNetwork(GraphToHoldScores & gths);
	void setTimeLimit( double time_limit_in_seconds );
	void clear();																			
	bool cycleReductionHasBegun() const;
	bool safeS3ReductionHasBegun() const;
	bool considerNodeForSafeS3Reduction( int node );
	bool existsEdgeBetween(int firstNodeIndex, int secondNodeIndex);
	EdgeBetweenMoveableNodes* getEdgeBetweenMoveableNodes(int firstNodeIndex, int secondNodeIndex);
	EdgeBetweenMoveableNodes* addNewEdge(int firstNodeIndex, int secondNodeIndex);
	void haveReportedOptimalNetworkScoreForSingleton(double optScore);
	bool d3edgeExistsBetween( int fn, int sn, int tn );
	Degree3Hyperedge * getD3Edge( int fn, int sn, int tn );
	
	void bruteForceOriginalGraph();
	std::vector<NodeState> const & getOptimalNetworkState(){return _optimalNetworkStateVector;};	
	int getNumStatesForNode( int node ) const;
	MoveableNode* getMoveableNode( int node ) const { return _theMNVector[ node ]; }
	void noteEffort( double effort );
};

class QueueOfToBeEliminatedNodes
{
public:
	static QueueOfToBeEliminatedNodes* getInstance();
	~QueueOfToBeEliminatedNodes();
	void addNodeIndexForTreeReduction(int nodeIndex);
	void addNodeIndexForCycleReduction(int nodeIndex);
	void addNodeIndexForSafeS3Reduction( int nodeIndex );
	void addNodeIndexForUnsafeS3Reduction( int nodeIndex );
	bool isEmpty();
	int  nextNodeForReduction();

private:
	QueueOfToBeEliminatedNodes();		//Singleton
	static bool							_instanceFlag;
	static QueueOfToBeEliminatedNodes*	_theQueue;
	std::list<int>						_QForTreeRed;
	std::list<int>						_QForCycleRed;
	std::list<int>						_QForSafeS3Red;
	std::list<int>						_QForUnsafeS3Red;
};

class DummyNetworkDescriptorClass
{
public:
	DummyNetworkDescriptorClass(int numNodes, int numEdges, AtomPositions * xyz);
	~DummyNetworkDescriptorClass();
	DummyNetworkDescriptorClass(const DummyNetworkDescriptorClass& rhs);

	void setMaxStatesForNode(int nodeIndex, NodeState maxNodeStates);
	void setEdgeBetweenNodes(int EdgeNumber, int firstNodeIndex, int secondNodeIndex);
	void setScoreForNode(int nodeIndex, NodeState firstNodeState, double theScore);
	void setScoreForEdge(int EdgeNumber, NodeState firstNodeState, NodeState secondNodeState, double theScore);

	int getNumberOfNodes() const;
	int getNumberOfEdges() const;
	NodeState getMaxStatesForNode(int nodeIndex) const;
	int getFirstNodeIndexForEdge(int edgeIndex) const;
	int getSecondNodeIndexForEdge(int edgeIndex) const;
	double getNodeScoreForState(int nodeIndex, NodeState theState) const;
	double getEdgeScoreForEdge(int edgeNumber, NodeState firstState, NodeState secondState) const;

	int addDegreeThreeHyperedge();
	void setNodesDeg3HypererdgeIncidentUpon( int hyperedge, int first_node, int second_node, int third_node );
	void setMoverForD3H( int hyperedge, int node_index, Mover * mover );
	void setAtomsInThreeWayOverlap( int hyperedge, int node_index, std::list<AtomDescr> & atomlist );
	void addAtomsToThreeWayOverlap( int hyperedge, int node_index, std::list<AtomDescr> & atomlist );
	void setAllAtomsForMover( int hyperedge, int node_index, std::list<AtomDescr> & atomlist );

	AtomPositions * getAtomPositionsPointer() const;
	int getNumHyperedges() const;
	int getVertexIdForD3H(int edge, int vertexID ) const;
	Mover* getMoverForD3H( int edge, int vertexID ) const;
	std::list< AtomDescr > const & getAtomsIn3WayOverlapForD3H( int edge, int vertexID ) const;
		
	void setAllStatesEnabled();
	void disableStateOnNode( int vertex, int state );
	void setStateDisablingCompleteForNow();
	int convertEnabledState2OriginalStateEnumerationOnNode( int vertex, int enabled_state ) const;

	bool getEdgeExists(int vertex1, int vertex2) const;
	int addEdge();
	
	bool getDegree3HyperedgeExists(int v1, int v2, int v3 ) const;
	bool nodesOnD3H( int v1, int v2 ) const;
	int getNodeOnD3HWithVertex( int vertex ) const;
	int getDegree3HyperedgeIndex(int v1, int v2, int v3 ) const;
	void finalizeD3HEs();
	
private:
	
	void putAtomsIn3WOverlapOnAllAppropriateHEs();
	void forceCliquesForD3HEsThatShareAtoms();
	bool HEsShareIncidentNode( int const he1, int const he2 );
	bool shareAtomsBetweenHEs( int const he1, int const he2 );
	void forceCliqueForD3HEPairIfShareAtoms(int const d3h1, int const d3h2 );
	void forceCliqueForD3HEPair( int const d3h1, int const d3h2 );
	
	AtomPositions * _xyz;
	int _numNodes;
	int _numEdges;
	std::vector<int> _firstNodeIndexForEdgeArray;
	std::vector<int> _secondNodeIndexForEdgeArray;
	NodeState* _nodeMaxStatesArray;
	double** _nodeScoreArrayOfArrays;
	std::vector< double** > _edgeScoreArrayOfMatrices;
	int _numberDegree3Hyperedges;
	bool _alreadyFinalizedD3HEs;
	std::vector< std::vector< int > > _verticesForDegree3Hyperedge;
	std::vector< std::vector< Mover* > > _moversForDegree3Hyperedge;
	std::vector< std::vector< std::list< AtomDescr > > > _atomsInThreeWayOverlap;
	std::vector< std::vector< std::list< AtomDescr > > > _allAtomsOnMover;

	bool _anyStatesDisabled;
	std::vector< int > _numEnabledStatesForVertex;
	std::vector< std::vector< bool > > _stateEnabled;
	std::vector< std::vector< int > > _enabled2Original;
};

#endif
