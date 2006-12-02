// name: MoveableNode.C
// author: Andrew Leaver-Fay
// date written: 4/15/06
// purpose: Implementation for Interaction Graph Classes

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************
#if defined(_MSC_VER)
#pragma warning(disable:4786) 
#endif

#include "Mover.h"
#include "AtomPositions.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
using std::ifstream;
using std::cout;
using std::cerr;
using std::endl;
using std::ios_base;

#ifdef OLD_STD_HDRS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <time.h>
#else
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cassert>
#include <ctime>
using std::exit;
using std::clock;
using std::clock_t;
#endif

#include "MoveableNode.h"

//---------------------------------------------------------------------------------------------------------------

OptimizedNodeStateVector::OptimizedNodeStateVector() { _dependentNodeIndex = -1;}
OptimizedNodeStateVector::OptimizedNodeStateVector(NodeState maxNodeStates) {_dependentNodeIndex = -1; _theNSVector.assign((int) maxNodeStates, (NodeState) - 1);}
OptimizedNodeStateVector::OptimizedNodeStateVector(int depNodeIndex, NodeState maxNodeStates) : _dependentNodeIndex(depNodeIndex) {_theNSVector.assign((int) maxNodeStates, (NodeState) -1);}
OptimizedNodeStateVector::~OptimizedNodeStateVector() {}
OptimizedNodeStateVector::OptimizedNodeStateVector(const OptimizedNodeStateVector& rhs) {_theNSVector = rhs._theNSVector; _dependentNodeIndex = rhs._dependentNodeIndex;}
void OptimizedNodeStateVector::setMaxNodeStates(NodeState maxNodeStates)
{
	_theNSVector.resize((int) maxNodeStates);

	for (int i=0;i<maxNodeStates;i++)
		_theNSVector[i] = (NodeState) -1;

	return;
}
NodeState OptimizedNodeStateVector::getOptimumNodeState(NodeState dependentNodeState)
{
	if (((int) dependentNodeState < _theNSVector.size() ) && ((int) dependentNodeState >= 0))
	{
		return _theNSVector[(int) dependentNodeState];
	}
	else
	{
		std::cerr << "Error Accessing Negative or Out of Bounds Index for OptimizedNodeStateVector. -1 returned: " << dependentNodeState << endl;
		return -1;
	}
}
void OptimizedNodeStateVector::setOptimizedNodeState(NodeState dependentNodeState, NodeState ownerNodeState)
{
	if (((int) dependentNodeState < _theNSVector.size() ) && ((int) dependentNodeState >= 0))
	{
		_theNSVector[(int) dependentNodeState] = ownerNodeState;
		return;
	}
	else
	{
		std::cerr << "Error in Assigning with Negative or Out of Bounds Index for OptimizedNodeStateVector. No assignment Performed: " << dependentNodeState << endl;
		return;
	}
}


//---------------------------------------------------------------------------------------------------------------

OptimizedNodeStateMatrix::OptimizedNodeStateMatrix() : _theNSMatrix(NULL), _1stNodeIndex(-1), _2ndNodeIndex(-1), _1stNodeMaxStates((NodeState) -1), _2ndNodeMaxStates((NodeState) -1) {}
OptimizedNodeStateMatrix::OptimizedNodeStateMatrix(NodeState firstNodeMaxStates, NodeState secondNodeMaxStates) : _1stNodeIndex(-1), _2ndNodeIndex(-1), _1stNodeMaxStates(firstNodeMaxStates), _2ndNodeMaxStates(secondNodeMaxStates)
{
	_theNSMatrix = NULL;
	if ((_1stNodeMaxStates > 0) && (_2ndNodeMaxStates > 0))
	{
		_theNSMatrix = new NodeState*[_1stNodeMaxStates];
		for (int i = 0; i < _1stNodeMaxStates; i++)
		{
			_theNSMatrix[i] = new NodeState[_2ndNodeMaxStates];
			for (int j = 0; j < _2ndNodeMaxStates; j++)
			{
				_theNSMatrix[i][j] = (NodeState) -1;
			}
		}
	}
}
OptimizedNodeStateMatrix::OptimizedNodeStateMatrix(int firstNodeIndex, int secondNodeIndex, NodeState firstNodeMaxStates, NodeState secondNodeMaxStates) : _1stNodeIndex(firstNodeIndex), _2ndNodeIndex(secondNodeIndex), _1stNodeMaxStates(firstNodeMaxStates), _2ndNodeMaxStates(secondNodeMaxStates)
{
	_theNSMatrix = NULL;
	if ((_1stNodeMaxStates > 0) && (_2ndNodeMaxStates > 0))
	{
		_theNSMatrix = new NodeState*[_1stNodeMaxStates];
		for (int i = 0; i < _1stNodeMaxStates; i++)
		{
			_theNSMatrix[i] = new NodeState[_2ndNodeMaxStates];
			for (int j = 0; j < _2ndNodeMaxStates; j++)
			{
				_theNSMatrix[i][j] = (NodeState) -1;
			}
		}
	}
}
OptimizedNodeStateMatrix::~OptimizedNodeStateMatrix()
{
	deallocateMatrix();
	_1stNodeMaxStates = (NodeState) -1;
	_2ndNodeMaxStates = (NodeState) -1;
}
OptimizedNodeStateMatrix::OptimizedNodeStateMatrix(const OptimizedNodeStateMatrix& rhs)
{
	deallocateMatrix();

	_theNSMatrix = rhs._theNSMatrix;
	_1stNodeMaxStates = rhs._1stNodeMaxStates;
	_2ndNodeMaxStates = rhs._2ndNodeMaxStates;
	_1stNodeIndex     = rhs._1stNodeIndex;
	_2ndNodeIndex     = rhs._2ndNodeIndex;
	//TEST 
	//cerr << "Copy constructor for OptimizedNodeStateMatrix" << endl;
}
void OptimizedNodeStateMatrix::DeepCopy(const OptimizedNodeStateMatrix& toBeCopied)
{
	_1stNodeMaxStates = toBeCopied._1stNodeMaxStates;
	_2ndNodeMaxStates = toBeCopied._2ndNodeMaxStates;
	_1stNodeIndex     = toBeCopied._1stNodeIndex;
	_2ndNodeIndex     = toBeCopied._2ndNodeIndex;

	deallocateMatrix();
	if ((_1stNodeMaxStates > 0) && (_2ndNodeMaxStates > 0))
	{
		_theNSMatrix = new NodeState*[_1stNodeMaxStates];
		for (int i = 0; i < _1stNodeMaxStates; i++)
		{
			_theNSMatrix[i] = new NodeState[_2ndNodeMaxStates];
			for (int j = 0; j < _2ndNodeMaxStates; j++)
			{
				_theNSMatrix[i][j] = toBeCopied._theNSMatrix[i][j];
			}
		}
	}
	return;
}
void OptimizedNodeStateMatrix::deallocateMatrix()
{
	if (_theNSMatrix != NULL)
	{
		for (int i = 0; i < _1stNodeMaxStates; i++)
		{
			if (_theNSMatrix[i] != NULL)
			{
				delete [] _theNSMatrix[i];
				_theNSMatrix[i] = NULL;
			}
		}
		delete [] _theNSMatrix;
		_theNSMatrix = NULL;
	}
	return;
}

NodeState OptimizedNodeStateMatrix::getOptimumNodeState(NodeState firstNodeState, NodeState secondNodeState)
{
	if (((int) firstNodeState >= 0) && (firstNodeState < _1stNodeMaxStates))
		if (((int) secondNodeState >= 0) && (secondNodeState < _2ndNodeMaxStates))
			return _theNSMatrix[(int) firstNodeState][(int) secondNodeState];
		else
			cerr << "Second node index negative or exceeding secondNodeMaxStates: " << secondNodeState << endl;
	else
		cerr << "First node index negative or exceeding firstNodeMaxStates: " << firstNodeState << endl;
	return (NodeState) -1;
}
void OptimizedNodeStateMatrix::setOptimumNodeState(NodeState firstNodeState, NodeState secondNodeState, NodeState ownerNodeState)
{
	if (((int) firstNodeState >= 0) && (firstNodeState < _1stNodeMaxStates))
		if (((int) secondNodeState >= 0) && (secondNodeState < _2ndNodeMaxStates))
			_theNSMatrix[(int) firstNodeState][(int) secondNodeState] = ownerNodeState;
		else
			cerr << "Second node index negative or exceeding secondNodeMaxStates: " << secondNodeState << endl;
	else
		cerr << "First node index negative or exceeding firstNodeMaxStates: " << firstNodeState << endl;
	return;
}
//---------------------------------------------------------------------------------------------------------------

OptimizedNodeStateMatrix3::OptimizedNodeStateMatrix3() : _theNSMatrix(NULL), _1stNodeIndex(-1), _2ndNodeIndex(-1), _1stNodeMaxStates((NodeState) -1), _2ndNodeMaxStates((NodeState) -1) {}
OptimizedNodeStateMatrix3::OptimizedNodeStateMatrix3(
	NodeState firstNodeMaxStates, 
	NodeState secondNodeMaxStates,
	NodeState thirdNodeMaxStates) 
: 
	_1stNodeIndex(-1), _2ndNodeIndex(-1), _3rdNodeIndex( -1 ),
	_1stNodeMaxStates(firstNodeMaxStates), 
	_2ndNodeMaxStates(secondNodeMaxStates),
	_3rdNodeMaxStates(thirdNodeMaxStates)
{
	_theNSMatrix = NULL;
	if ((_1stNodeMaxStates > 0) && (_2ndNodeMaxStates > 0) && (_3rdNodeMaxStates > 0))
	{
		_theNSMatrix = new NodeState**[_1stNodeMaxStates];
		for (int i = 0; i < _1stNodeMaxStates; i++)
		{
			_theNSMatrix[i] = new NodeState*[_2ndNodeMaxStates];
			for (int j = 0; j < _2ndNodeMaxStates; j++)
			{
				_theNSMatrix[i][j] = new NodeState[_3rdNodeMaxStates];
				for (int k = 0; k < _3rdNodeMaxStates; ++k)
				{
					_theNSMatrix[i][j][k] = (NodeState) -1;
				}
			}
		}
	}
}
OptimizedNodeStateMatrix3::OptimizedNodeStateMatrix3(
	int firstNodeIndex,
	int secondNodeIndex,
	int thirdNodeIndex,
	NodeState firstNodeMaxStates,
	NodeState secondNodeMaxStates,
	NodeState thirdNodeMaxStates
)
	: 
	_1stNodeIndex(firstNodeIndex),
	_2ndNodeIndex(secondNodeIndex),
	_3rdNodeIndex(thirdNodeIndex),
	_1stNodeMaxStates(firstNodeMaxStates), 
	_2ndNodeMaxStates(secondNodeMaxStates),
	_3rdNodeMaxStates(thirdNodeMaxStates)
{
	_theNSMatrix = NULL;
	if ((_1stNodeMaxStates > 0) && (_2ndNodeMaxStates > 0) && (_3rdNodeMaxStates > 0))
	{
		_theNSMatrix = new NodeState**[_1stNodeMaxStates];
		for (int i = 0; i < _1stNodeMaxStates; i++)
		{
			_theNSMatrix[i] = new NodeState*[_2ndNodeMaxStates];
			for (int j = 0; j < _2ndNodeMaxStates; j++)
			{
				_theNSMatrix[i][j] = new NodeState[_3rdNodeMaxStates];
				for (int k = 0; k < _3rdNodeMaxStates; ++k)
				{
					_theNSMatrix[i][j][k] = (NodeState) -1;
				}
			}
		}
	}
}

OptimizedNodeStateMatrix3::~OptimizedNodeStateMatrix3()
{
	deallocateMatrix();
	_1stNodeMaxStates = (NodeState) -1;
	_2ndNodeMaxStates = (NodeState) -1;
	_3rdNodeMaxStates = (NodeState) -1;
}

void OptimizedNodeStateMatrix3::deallocateMatrix()
{
	if (_theNSMatrix != NULL)
	{
		for (int i = 0; i < _1stNodeMaxStates; i++)
		{
			if (_theNSMatrix[i] != NULL)
			{
				for ( int j = 0; j < _2ndNodeMaxStates; j++)
				{
					delete [] _theNSMatrix[i][j]; _theNSMatrix[i][j] = 0;
				}
				_theNSMatrix[i] = 0;
			}
		}
		delete [] _theNSMatrix;_theNSMatrix = NULL;
	}
	return;
}

NodeState OptimizedNodeStateMatrix3::getOptimumNodeState(
	NodeState firstNodeState,
	NodeState secondNodeState,
	NodeState thirdNodeState)
{
	if (((int) firstNodeState >= 0) && (firstNodeState < _1stNodeMaxStates))
	{
		if (((int) secondNodeState >= 0) && (secondNodeState < _2ndNodeMaxStates))
		{
			if (((int) thirdNodeState >= 0) && (thirdNodeState < _3rdNodeMaxStates))
			{
				return _theNSMatrix[(int) firstNodeState][(int) secondNodeState][(int) thirdNodeState];
			}
			else
			{
				cerr << "Third node index negative or exceeding secondNodeMaxStates: " << secondNodeState << endl;
			}
		}
		else
		{
			cerr << "Second node index negative or exceeding secondNodeMaxStates: " << secondNodeState << endl;
		}
	}
	else
	{
		cerr << "First node index negative or exceeding firstNodeMaxStates: " << firstNodeState << endl;
	}
	return (NodeState) -1;
}

void OptimizedNodeStateMatrix3::setOptimumNodeState(
	NodeState firstNodeState,
	NodeState secondNodeState,
	NodeState thirdNodeState,
	NodeState ownerNodeState
)
{
	if (((int) firstNodeState >= 0) && (firstNodeState < _1stNodeMaxStates))
	{
		if (((int) secondNodeState >= 0) && (secondNodeState < _2ndNodeMaxStates))
		{
			if (((int) thirdNodeState >= 0) && (thirdNodeState < _3rdNodeMaxStates))
			{
				_theNSMatrix[(int) firstNodeState][(int) secondNodeState][(int) thirdNodeState] = ownerNodeState;
			}
			else
			{
				cerr << "Third node index negative or exceeding secondNodeMaxStates: " << secondNodeState << endl;
			}
		}
		else
		{
			cerr << "Second node index negative or exceeding secondNodeMaxStates: " << secondNodeState << endl;
		}
	}
	else
	{
		cerr << "First node index negative or exceeding firstNodeMaxStates: " << firstNodeState << endl;
	}
	return;
}

//---------------------------------------------------------------------------------------------------------------

NodeScoreVector::NodeScoreVector() : _theNScrVect(0) {}
NodeScoreVector::NodeScoreVector(NodeState ownerMaxStates) {_theNScrVect.assign((std::vector<double>::size_type)ownerMaxStates, 0.);}
NodeScoreVector::~NodeScoreVector() {}
//NodeScoreVector::NodeScoreVector(const NodeScoreVector& rhs) //unimplemented
double NodeScoreVector::getNodeScore(NodeState ownersState)
{
	if (((int) ownersState >= 0) && ((int) ownersState < _theNScrVect.size()))
		return _theNScrVect[(int) ownersState];
	else
	{
		cerr << "Out of bounds index into NodeScoreVector (neg or >= MaxStates).  Zero returned: " << ownersState << endl;
		return 0;
	}

}

void   NodeScoreVector::setNodeScore(NodeState ownersState, double theScore)
{
	if (((int) ownersState >= 0) && ((int) ownersState < _theNScrVect.size()))
	{
		_theNScrVect[(int) ownersState] = theScore;
		return;
	}
	else
	{
		cerr << "Out of bounds index into NodeScoreVector (neg or >= MaxStates).  Zero returned: " << ownersState << endl;
		return;
	}
}

void   NodeScoreVector::addToNodeScore(NodeState ownersState, double theScore)
{
	if (((int) ownersState >= 0) && ((int) ownersState < _theNScrVect.size()))
	{
		_theNScrVect[(int) ownersState] += theScore;
		return;
	}
	else
	{
		cerr << "Out of bounds index into NodeScoreVector (neg or >= MaxStates).  Zero returned: " << ownersState << endl;
		return;
	}
}

//---------------------------------------------------------------------------------------------------------------

EdgeScoreMatrix::EdgeScoreMatrix() : _theEScrMatrix(NULL), _1stNodeMaxStates((NodeState) -1), _2ndNodeMaxStates((NodeState) -1) {}
EdgeScoreMatrix::EdgeScoreMatrix(NodeState firstNodeMaxStates, NodeState secondNodeMaxStates) : _1stNodeMaxStates(firstNodeMaxStates),
	 _2ndNodeMaxStates(secondNodeMaxStates) 
{
	if ((_1stNodeMaxStates > 0) && (_2ndNodeMaxStates > 0))
	{
		_theEScrMatrix = new double*[_1stNodeMaxStates];
		for (int i = 0; i < _1stNodeMaxStates; i++)
		{
			_theEScrMatrix[i] = new double[_2ndNodeMaxStates];
			for (int j = 0; j < _2ndNodeMaxStates; j++)
				_theEScrMatrix[i][j] = 0;
		}
	}
}
EdgeScoreMatrix::~EdgeScoreMatrix()
{
	deallocateMatrix();
	_1stNodeMaxStates = -1;
	_2ndNodeMaxStates = -1;
}
EdgeScoreMatrix::EdgeScoreMatrix(const EdgeScoreMatrix& rhs)
{
	deallocateMatrix();
	_theEScrMatrix = rhs._theEScrMatrix;
	_1stNodeMaxStates = rhs._1stNodeMaxStates;
	_2ndNodeMaxStates = rhs._2ndNodeMaxStates;
	_1stNodeIndex     = rhs._1stNodeIndex;
	_2ndNodeIndex     = rhs._2ndNodeIndex;
	//TEST 
	//cerr << "Copy constructor for EdgeScoreMatrix" << endl;
}
void EdgeScoreMatrix::DeepCopy(const EdgeScoreMatrix& toBeCopied)
{
	deallocateMatrix();
	_1stNodeMaxStates = toBeCopied._1stNodeMaxStates;
	_2ndNodeMaxStates = toBeCopied._2ndNodeMaxStates;
	_1stNodeIndex     = toBeCopied._1stNodeIndex;
	_2ndNodeIndex     = toBeCopied._2ndNodeIndex;
	if ((_1stNodeMaxStates > 0) && (_2ndNodeMaxStates > 0))
	{
		_theEScrMatrix = new double*[_1stNodeMaxStates];
		for (int i = 0; i < _1stNodeMaxStates; i++)
		{
			_theEScrMatrix[i] = new double[_2ndNodeMaxStates];
			for (int j = 0; j < _2ndNodeMaxStates; j++)
				_theEScrMatrix[i][j] = toBeCopied._theEScrMatrix[i][j];
		}
	}

	return;
}

double EdgeScoreMatrix::getEdgeScore(NodeState firstNodeState, NodeState secondNodeState)
{
	if (((int) firstNodeState >= 0) && (firstNodeState < _1stNodeMaxStates))
		if (((int) secondNodeState >= 0) && (secondNodeState < _2ndNodeMaxStates))
			return _theEScrMatrix[(int) firstNodeState][(int) secondNodeState];
		else
			cerr << "EdgeScoreMatrix: Second node index negative or exceeding secondNodeMaxStates: " << secondNodeState << endl;
	else
		cerr << "EdgeScoreMatrix: First node index negative or exceeding firstNodeMaxStates: " << firstNodeState << endl;
	return 0;
}
void   EdgeScoreMatrix::setEdgeScore(NodeState firstNodeState, NodeState secondNodeState, double theScore)
{
	if (((int) firstNodeState >= 0) && (firstNodeState < _1stNodeMaxStates))
		if (((int) secondNodeState >= 0) && (secondNodeState < _2ndNodeMaxStates))
			_theEScrMatrix[(int) firstNodeState][(int) secondNodeState] = theScore;
		else
			cerr << "EdgeScoreMatrix: Second node index negative or exceeding secondNodeMaxStates: " << secondNodeState << endl;
	else
		cerr << "EdgeScoreMatrix: First node index negative or exceeding firstNodeMaxStates: " << firstNodeState << endl;
	return;
}

void   EdgeScoreMatrix::addEdgeScore(NodeState firstNodeState, NodeState secondNodeState, double theScore)
{
	if (((int) firstNodeState >= 0) && (firstNodeState < _1stNodeMaxStates))
		if (((int) secondNodeState >= 0) && (secondNodeState < _2ndNodeMaxStates))
			_theEScrMatrix[(int) firstNodeState][(int) secondNodeState] += theScore;
		else
			cerr << "EdgeScoreMatrix: Second node index negative or exceeding secondNodeMaxStates: " << secondNodeState << endl;
	else
		cerr << "EdgeScoreMatrix: First node index negative or exceeding firstNodeMaxStates: " << firstNodeState << endl;
	return;
}


void EdgeScoreMatrix::deallocateMatrix()
{
	if (_theEScrMatrix != NULL)
	{
		for (int i = 0; i < _1stNodeMaxStates; i++)
		{
			if (_theEScrMatrix[i] != NULL)
			{
				delete [] _theEScrMatrix[i];
				_theEScrMatrix[i] = NULL;
			}
		}
		delete [] _theEScrMatrix;
		_theEScrMatrix = NULL;
	}
	return;
}
//---------------------------------------------------------------------------------------------------------------

NodeAndStatePair::NodeAndStatePair() : _nodeIndex(-1), _nodeState((NodeState) -1) {}
NodeAndStatePair::NodeAndStatePair(int index, NodeState state) : _nodeIndex(index), _nodeState(state) {}
//NodeAndStatePair::~NodeAndStatePair();
//NodeAndStatePair::NodeAndStatePair(const NodeAndStatePair& rhs);
int		  NodeAndStatePair::getNodeIndex() {return _nodeIndex;}
NodeState NodeAndStatePair::getNodeState() {return _nodeState;}
void	  NodeAndStatePair::setNodeIndex(int index) {_nodeIndex = index; return;}
void	  NodeAndStatePair::setNodeState(NodeState state) {_nodeState = state; return;}

//---------------------------------------------------------------------------------------------------------------

MoveableNode::MoveableNode() : _index(-1), _maxNodeStates((NodeState) -1), _theNodeScoreVect(NULL),
_isEliminated(false), _eliminationOrder(0), _optimalStateVector(NULL), _nodeIndexDetStateFromTreeReduction(-1),
_optimalStateMatrix(NULL), _firstNodeIndexDetStateFromCycleReduction(-1), _secondNodeIndexDetStateFromCycleReduction(-1), 
_hasAnyDegree3Hyperedges( false ), _hasAnyHighOrderOverlap( false ), _mover( NULL )
{}
MoveableNode::~MoveableNode()
{
	if (_theNodeScoreVect != NULL)
	{
		delete _theNodeScoreVect;
		_theNodeScoreVect = NULL;
	}
	if (_optimalStateVector != NULL)
	{
		delete _optimalStateVector;
		_optimalStateVector = NULL;
	}
	if (_optimalStateMatrix != NULL)
	{
		delete _optimalStateMatrix;
		_optimalStateMatrix = NULL;
	}
	if ( _degree3hyperedges.size() != 0 )
	{
		for (std::list< Degree3Hyperedge * >::iterator iter = _degree3hyperedges.begin(); 
			iter != _degree3hyperedges.end(); )
		{
			std::list< Degree3Hyperedge * >::iterator iternext = iter;
			++iternext;
			delete (*iter); iter = iternext;
		}
	}
	if ( _degree4hyperedges.size() != 0 )
	{
		for (std::list< Degree4Hyperedge * >::iterator iter = _degree4hyperedges.begin(); 
			iter != _degree4hyperedges.end(); )
		{
			std::list< Degree4Hyperedge * >::iterator iternext = iter;
			++iternext;
			delete (*iter); iter = iternext;
		}
	}
}

//MoveableNode::MoveableNode(const MoveableNode& rhs);
MoveableNode::MoveableNode(NodeState maxStates) : 
	_index(-1),
	_maxNodeStates(maxStates),
        _theNodeScoreVect(NULL),
	_isEliminated(false), _eliminationOrder(0), _optimalStateVector(NULL),
	_nodeIndexDetStateFromTreeReduction(-1),
	_optimalStateMatrix(NULL), _firstNodeIndexDetStateFromCycleReduction(-1),
	_secondNodeIndexDetStateFromCycleReduction(-1),
        _hasAnyDegree3Hyperedges(false),_hasAnyDegree4Hyperedges(false),
	_hasAnyHighOrderOverlap( false ), _mover( NULL )
{
	_theNodeScoreVect = new NodeScoreVector(_maxNodeStates);
}

void				MoveableNode::setNodeStateScore(NodeState theNodeState, double theScore) 
{
	_theNodeScoreVect->setNodeScore(theNodeState,theScore);
	return;
}
double				MoveableNode::getNodeStateScore(NodeState theNodeState) {return _theNodeScoreVect->getNodeScore(theNodeState);}
void 				MoveableNode::eliminate()
{
	//This general purpose function will decide which of the two elimination functions should be called
	//and then will be responsible for notifying the NodeAndEdgeManager, it's edge(s), and it's neighbor(s)
	//of the fact that it has been eliminated.  A neighbor of an eliminated node must remove the edge that
	//connected the two from it's edge list, and then decide whether or not it should add itself to the 
	//queue of to-be-eliminated nodes.  The edges need to be marked as eliminated so that they can be ignored
	//by the brute-force algorithm that might require being called later.  The NodeAndEdgeManager will add the
	//index of this node to the top of the partialOrderDependency stack - the optimal state of this node can
	//be determined using only the states of the nodes that still remain in the graph at the time this function is called.
	//This function will also store the elimination order for this node, which will later be used to determine
	//the optimal state of this node from the states of other nodes already determined to be optimal.

	//Since a node may be put into the elimination queue for cycle reduction, and then have it's connectivity change so that
	//it becomes part of a tree (and was eliminated by tree reduction) we need to first check to see that the node still requires
	//our attention - if not, we go on.
	if (_isEliminated) return;

	//TEST
	//cerr << "Begining Elimination of Node: " << _index << endl;

	NodeAndEdgeManager* theNaEManager = NodeAndEdgeManager::getInstance();

	if ( _maxNodeStates == 1 )
	{
		eliminateOneStateVertex();
		_eliminationOrder = 3; //same as singleton elimination
	}
	else if (_edgeList.size() == 1)
	{
		eliminateThroughTreeReduction();
		_eliminationOrder = 1;
	}
	else if (_edgeList.size() == 2)
	{
		eliminateThroughCycleReduction();
		_eliminationOrder = 2;

	}
	else if (_edgeList.size() == 3)
	{
		eliminateThroughS3Reduction();
		_eliminationOrder = 4;
	}
	else if (_edgeList.size() == 0)
	{
		//Entire graph has been eliminated - this is the last node, "it's" optimum score
		//represents the optimum score for the entirety of the graph.  From the state of
		//this node in it's optimum configuration, the states of the rest of the nodes in the
		//graph can be inferred.
		
		//Compute nothing as of now.  Notify theNaEManager that this node, too, has been eliminated

		//TEST
		//cerr << "Last Node in reduced graph being eliminated " << endl;

		eliminateSingleton();
		_eliminationOrder = 3; 
	}
	else {
		return;
	}

	
	for (std::list<EdgeBetweenMoveableNodes*>::iterator iter = _edgeList.begin();
		iter != _edgeList.end(); ++iter)
	{
	
		(*iter)->setIsEliminated(true);
		MoveableNode* otherNodePtr = (*iter)->getPointerToOtherNode(_index);
		otherNodePtr->beNotifiedDependencyEliminated(_index);
	}
	_edgeList.clear();
	_isEliminated = true;
	theNaEManager->BeNotifiedOfEliminatedNode(_index);

	//cerr << "Finished Elimination of Node: " << _index << endl;
	return;

}
void 				MoveableNode::eliminateThroughTreeReduction()
{
	//For the sake of this description assume the following: 
	//This node is node "i", we know it has a single edge, so call the node it is connected to node "j"
	//There are n_i possible states for node i and there are n_j possible states for node j.
	//This function will calculate the optimal state for node i in all of node j's possible states,
	//and store this state in the newly allocated _optimalNodeStateVector, an array with n_j entires.
	//It will take the optimal score for each of the n_j states of node j and add them to the 
	//node j's _nodeScoreVector.  Other bookkeeping is taken care of by eliminate().

	//"Optimal" is defined within a static function of the NodeAndEdgeManager.  For energies
	// Optimal(a,b) (a if a < b; b o.w.).  For use with 'reduce', Optimal(a,b) (a if a > b; b o.w.)

	//Tree elimination is done in O(n_j*n_i) time.

	//cerr << "Begining Tree Reduction for Node: " << _index << endl;
	if ( _hasAnyDegree3Hyperedges )
	{
		std::cerr << "Critical Error: degree three hyperedge on node with only one neighbor; d3he requires node have two neighbors! " << std::endl;
	}

	std::list<EdgeBetweenMoveableNodes*>::iterator theCommonEdge = _edgeList.begin();
	//TEST
	//cerr << "Common Edge: " << (*theCommonEdge) << endl;
	MoveableNode* otherNodePtr = (*theCommonEdge)->getPointerToOtherNode(_index);
	//TEST
	//cerr << "otherNodePtr = " << otherNodePtr << endl;
	int otherNodeIndex = otherNodePtr->getNodeIndex();
	NodeState numStatesOtherNode = otherNodePtr->getNumberOfPossibleStates();
	NodeScoreVector optScoreVector = NodeScoreVector( numStatesOtherNode);

	_optimalStateVector = new OptimizedNodeStateVector( numStatesOtherNode);
	_nodeIndexDetStateFromTreeReduction = otherNodeIndex;

	NodeAndEdgeManager* theNaEManager = NodeAndEdgeManager::getInstance();
	theNaEManager->noteEffort( numStatesOtherNode * _maxNodeStates );

	double bestScore, holdScore;
	NodeState bestState;
	for (NodeState nsIndOther = (NodeState) 0; nsIndOther < numStatesOtherNode; nsIndOther++)
	{
		bestScore = _theNodeScoreVect->getNodeScore((NodeState) 0) + (*theCommonEdge)->getEdgeScore(_index,(NodeState) 0, otherNodeIndex, nsIndOther);
		bestState = 0;
		for (NodeState nsIndThis = (NodeState) 1; nsIndThis < _maxNodeStates; nsIndThis++)
		{
			holdScore = _theNodeScoreVect->getNodeScore(nsIndThis) + (*theCommonEdge)->getEdgeScore(_index, nsIndThis, otherNodeIndex, nsIndOther);
			if (NodeAndEdgeManager::betterThan(holdScore, bestScore))
			{
				bestScore = holdScore;
				bestState = nsIndThis;
			}
		}
		optScoreVector.setNodeScore(nsIndOther, bestScore);
		_optimalStateVector->setOptimizedNodeState(nsIndOther, bestState);
	}
	otherNodePtr->combineNodeScoreWithOptEdgeAndNodeScore(optScoreVector);

	//TEST
	//cerr << "Leaving eliminateThroughTreeReduction();" << endl;
	
	return;
}

void 				MoveableNode::eliminateThroughCycleReduction()
{
	//For the sake of this description, assume the following:
	//This is node "a", we already know it has two edges to two distinct nodes, call them "b" and "c"
	//Node a has n_a possible states.  Node b has n_b possible states.  Node c has n_c possible states.
	//Node b has e_b edges, one of which connects to a.  Node c has e_c edges, one of which connects to a.
	//This function will calculate the optimal state for node a for all n_b * n_c configurations of states for nodes b and c.
	//It will create an edge between nodes b and c if one does not exist with scores coming from the _nodeScoreVertex
	//of node a and the edges that connect nodes a and b and nodes a and c.  If an edge between b and c already exists
	//then it will add the optimal scores to the edgeScoreMatrix

	//"Optimal" is defined within a static function of the NodeAndEdgeManager.  For energies
	// Optimal(a,b) (a if a < b; b o.w.).  For use with 'reduce', Optimal(a,b) (a if a > b; b o.w.)

	//Cycle reduction runs in O(n_a*n_b*n_c + min(e_b, e_c) ) time since NodeAndEdgeManager::existsEdgeBetween() runs in min(e_b, e_c) time.

	//cerr << "Begining Cycle Reduction for node: " << _index << endl;

	std::list<EdgeBetweenMoveableNodes*>::iterator edgeIter = _edgeList.begin();
	EdgeBetweenMoveableNodes* firstEdgePtr  = *edgeIter;
	edgeIter++;
	EdgeBetweenMoveableNodes* secondEdgePtr = *edgeIter;
	int indexFirstOtherNode  = firstEdgePtr->getIndexOfOtherNode(_index);
	int indexSecondOtherNode = secondEdgePtr->getIndexOfOtherNode(_index);

	int tempIndex;
	EdgeBetweenMoveableNodes* tempEdgePtr;
	
	if (indexFirstOtherNode > indexSecondOtherNode) //Swap so that the first edge connects this node with the node of the smaller index
	{
		tempIndex            = indexFirstOtherNode;
		indexFirstOtherNode  = indexSecondOtherNode;
		indexSecondOtherNode = tempIndex;
		tempEdgePtr          = firstEdgePtr;
		firstEdgePtr         = secondEdgePtr;
		secondEdgePtr        = tempEdgePtr;
	}

	//int nodeIndices[ 3 ] = { indexFirstOtherNode, indexSecondOtherNode, _index };

	MoveableNode* firstOtherNodePtr  = firstEdgePtr->getPointerToOtherNode(_index);
	MoveableNode* secondOtherNodePtr = secondEdgePtr->getPointerToOtherNode(_index);

	NodeState numStatesFirstOther  = firstOtherNodePtr->getNumberOfPossibleStates();
	NodeState numStatesSecondOther = secondOtherNodePtr->getNumberOfPossibleStates(); 

	NodeState nsIndFirstOther, nsIndSecondOther, nsIndThis;

	double bestScore;
	NodeState bestState;
	_optimalStateMatrix = new OptimizedNodeStateMatrix(indexFirstOtherNode, indexSecondOtherNode, numStatesFirstOther, numStatesSecondOther);
	_firstNodeIndexDetStateFromCycleReduction  = indexFirstOtherNode;
	_secondNodeIndexDetStateFromCycleReduction = indexSecondOtherNode;

	double holdScore;

	NodeAndEdgeManager* theNaEManager = NodeAndEdgeManager::getInstance();
	if (theNaEManager->existsEdgeBetween(indexFirstOtherNode,indexSecondOtherNode))
	{
		Degree3Hyperedge * hyperedge = 0;
		if ( _hasAnyDegree3Hyperedges )
		{
			//std::cerr << "Node Has degree 3 hyperedge! " << std::endl;
			//std::cerr << "Node 1: " << nodeIndices[ 0 ];
			//std::cerr << " Node 2: " << nodeIndices[ 1 ];
			//std::cerr << " node 3: " << nodeIndices[ 2 ];
	
			hyperedge = *(_degree3hyperedges.begin());
			hyperedge->setVertexOrder( indexFirstOtherNode, indexSecondOtherNode, _index );
		}

		//TEST
		//cerr << "Edge found between the neighbors of " << _index << endl;
		EdgeBetweenMoveableNodes* edgeBtwnOtherNodes = theNaEManager->getEdgeBetweenMoveableNodes(indexFirstOtherNode,indexSecondOtherNode);
		for (nsIndFirstOther = 0; nsIndFirstOther < numStatesFirstOther; nsIndFirstOther++)
		{
			//if ( _hasAnyDegree3Hyperedges ) {hyperedge->resetMover( nodeIndices[ 1 ] );}

			for (nsIndSecondOther = 0; nsIndSecondOther < numStatesSecondOther; nsIndSecondOther++)
			{
			//if ( _hasAnyDegree3Hyperedges ){  hyperedge->resetMover( nodeIndices[ 2 ] );}

				
				bestScore = _theNodeScoreVect->getNodeScore((NodeState) 0) + firstEdgePtr->getEdgeScore(_index,(NodeState) 0, indexFirstOtherNode,nsIndFirstOther) + secondEdgePtr->getEdgeScore(_index,(NodeState) 0, indexSecondOtherNode, nsIndSecondOther);
				if ( _hasAnyDegree3Hyperedges ) {
					float threeWayDotScore = hyperedge->getScoreGivenSetVertexOrdering(nsIndFirstOther, nsIndSecondOther, 0); 
					//hyperedge->incrementMover( nodeIndices[ 2 ] ); 
					//std::cerr << nsIndFirstOther << " " << nsIndSecondOther << " : " << 0 << " = " << bestScore << " " << threeWayDotScore << " " << bestScore + threeWayDotScore << std::endl;
					bestScore += threeWayDotScore;
				}
				
				bestState = (NodeState) 0;
				for (nsIndThis = 1; nsIndThis < _maxNodeStates; nsIndThis++)
				{
					holdScore = _theNodeScoreVect->getNodeScore(nsIndThis) + firstEdgePtr->getEdgeScore(_index, nsIndThis, indexFirstOtherNode, nsIndFirstOther) + secondEdgePtr->getEdgeScore(_index, nsIndThis, indexSecondOtherNode, nsIndSecondOther);
					if ( _hasAnyDegree3Hyperedges )
					{
						//std::cerr << nodeIndices[ 2 ] << " -- " << nsIndThis << std::endl;
						float threeWayDotScore = hyperedge->getScoreGivenSetVertexOrdering(nsIndFirstOther, nsIndSecondOther, nsIndThis); 
						
						//hyperedge->incrementMover( nodeIndices[ 2 ] );
						//std::cerr << nsIndFirstOther << " " << nsIndSecondOther << " : " << nsIndThis << " = " << holdScore << " " << threeWayDotScore << " " << holdScore + threeWayDotScore << std::endl;
						holdScore += threeWayDotScore;
					}
					if (NodeAndEdgeManager::betterThan(holdScore, bestScore))
					{
						bestScore = holdScore;
						bestState = nsIndThis;
					}
				}
				_optimalStateMatrix->setOptimumNodeState(nsIndFirstOther,nsIndSecondOther,bestState);
				holdScore = edgeBtwnOtherNodes->getEdgeScore(nsIndFirstOther,nsIndSecondOther);
				edgeBtwnOtherNodes->setEdgeScore(nsIndFirstOther,nsIndSecondOther, (holdScore + bestScore));

			}
		}
		if (_hasAnyDegree3Hyperedges) {delete hyperedge; hyperedge = 0;}

	} 
	else
	{
		//TEST
		//cerr << "Edge NOT found between the neighbors of " << _index << endl;

		EdgeBetweenMoveableNodes* newEdgePtr = theNaEManager->addNewEdge(indexFirstOtherNode, indexSecondOtherNode);
		for (nsIndFirstOther = 0; nsIndFirstOther < numStatesFirstOther; nsIndFirstOther++)
		{

			for (nsIndSecondOther = 0; nsIndSecondOther < numStatesSecondOther; nsIndSecondOther++)
			{
				bestScore = _theNodeScoreVect->getNodeScore((NodeState) 0) + firstEdgePtr->getEdgeScore(_index,(NodeState) 0, indexFirstOtherNode,nsIndFirstOther) + secondEdgePtr->getEdgeScore(_index,(NodeState) 0, indexSecondOtherNode, nsIndSecondOther);
				bestState = (NodeState) 0;
				for (nsIndThis = 1; nsIndThis < _maxNodeStates; nsIndThis++)
				{
					holdScore = _theNodeScoreVect->getNodeScore(nsIndThis) + firstEdgePtr->getEdgeScore(_index, nsIndThis, indexFirstOtherNode, nsIndFirstOther) + secondEdgePtr->getEdgeScore(_index, nsIndThis, indexSecondOtherNode, nsIndSecondOther);
					if (NodeAndEdgeManager::betterThan(holdScore, bestScore))
					{
						bestScore = holdScore;
						bestState = nsIndThis;
					}
				}
				_optimalStateMatrix->setOptimumNodeState(nsIndFirstOther,nsIndSecondOther,bestState);
				newEdgePtr->setEdgeScore(nsIndFirstOther,nsIndSecondOther, bestScore);
			}
		}
	}
	theNaEManager->noteEffort( numStatesFirstOther * numStatesSecondOther * _maxNodeStates );
	return;
}

void
MoveableNode::eliminateThroughS3Reduction()
{
	//std::cerr << "eliminateThroughS3Reduction: " << _index << std::endl;
	assert( _edgeList.size() == 3);
	int nodes[ 3 ];
	int nodes_sorted[ 3 ];
	int count = 0;
	for ( std::list< EdgeBetweenMoveableNodes* >::iterator iter = _edgeList.begin();
		iter != _edgeList.end(); ++iter )
	{
		nodes[ count ] = (*iter)->getIndexOfOtherNode( _index );
		++count;
	}
	
	sort_three( nodes[ 0 ], nodes[ 1 ], nodes[ 2 ], nodes_sorted[ 0 ], nodes_sorted[ 1 ], nodes_sorted[ 2 ]);
	//std::cerr << "Neigbors: " << nodes_sorted[ 0 ] << " " << nodes_sorted[ 1 ] << " " <<  nodes_sorted[ 2 ] << std::endl;
	NodeAndEdgeManager* theNaEManager = NodeAndEdgeManager::getInstance();
	for (int ii = 0; ii < 3; ++ii)
	{
		for (int jj = ii+1; jj < 3; ++jj)
		{
			if ( ! theNaEManager->existsEdgeBetween( nodes_sorted[ ii ], nodes_sorted[ jj ]) )
			{
				theNaEManager->addNewEdge( nodes_sorted[ ii ], nodes_sorted[ jj ] );
			}
		}
	}
	
	Degree3Hyperedge * neighbors_d3he =  
		theNaEManager->d3edgeExistsBetween( nodes_sorted[ 0 ], nodes_sorted[ 1 ], nodes_sorted[ 2 ] )?
		theNaEManager->getD3Edge( nodes_sorted[ 0 ], nodes_sorted[ 1 ], nodes_sorted[ 2 ] ) :
		new Degree3Hyperedge( nodes_sorted[ 0 ], nodes_sorted[ 1 ], nodes_sorted[ 2 ] );

	int numStates[ 3 ] = {-1, -1, -1 };
	for (int ii = 0; ii < 3; ++ii)
	{
		numStates[ ii ] = theNaEManager->getMoveableNode( nodes_sorted[ ii ] )->getNumberOfPossibleStates();
	}
	theNaEManager->noteEffort( numStates[ 0 ] * numStates[1] * numStates[2] * _maxNodeStates );
	
	//std::cerr << "numStates:  " << numStates[ 0 ] << " " << numStates[1] << " " <<  numStates[2] << std::endl;
	
	//collect info for single degree-4 edge
	Degree4Hyperedge* degree4hyperedge (0);
	if ( ! _degree4hyperedges.empty() )
	{
		//std::cerr << "Degree-4 Hyperedge present in reduction" << std::endl;
		assert(_degree4hyperedges.size() == 1);
		std::list< Degree4Hyperedge * >::iterator iter = _degree4hyperedges.begin();
		degree4hyperedge = *iter;
		degree4hyperedge->setVertexOrder(
			nodes_sorted[ 0 ],
			nodes_sorted[ 1 ],
			nodes_sorted[ 2 ],
			_index );
	}
	
	//collect info for degree-3 edges
	count = 0;
	int const numd3edges = _degree3hyperedges.size();
	std::vector< Degree3Hyperedge* > d3hes( numd3edges );
	std::vector< std::pair< int, int > > d3whichVert( numd3edges );
	for ( std::list< Degree3Hyperedge* >::iterator iter = _degree3hyperedges.begin(); 
		iter != _degree3hyperedges.end(); ++iter )
	{
		d3hes[ count ] = *iter;
		int count2 = 0;
		for (int ii = 0; ii < 3; ++ii)
		{
			if ( (*iter)->getVertexIndex( ii ) != _index )
			{
				int whichOfThreeNeighbors = -1;
				for (int jj = 0; jj < 3; ++jj)
				{
					if ((*iter)->getVertexIndex( ii ) == nodes_sorted[ jj ]  )
					{
						whichOfThreeNeighbors = jj;
					}
				}
				if (count2 == 0 )
				{
					d3whichVert[ count ].first = whichOfThreeNeighbors;
				}
				else
				{
					d3whichVert[ count ].second = whichOfThreeNeighbors;
				}
				nodes[ count2 ] = (*iter)->getVertexIndex( ii );
				++count2;
			}
		}
		nodes[ 2 ] = _index;
		(*iter)->setVertexOrder( nodes[ 0 ], nodes[ 1 ], nodes[ 2 ] );
		++count;
	}
	
	//collect info for degree-2 hyperedges
	EdgeBetweenMoveableNodes * edges[ 3 ];
	int otherNodeForEdge[ 3 ];
	int whichOtherNode[ 3 ];
	count = 0;
	for ( std::list< EdgeBetweenMoveableNodes* >::iterator iter = _edgeList.begin();
		iter != _edgeList.end(); ++iter )
	{
		edges[ count ] = (*iter);
		otherNodeForEdge[ count ] = (*iter)->getIndexOfOtherNode( _index );
		for (int ii = 0; ii < 3; ++ii)
		{
			if ( otherNodeForEdge[ count ] == nodes_sorted[ ii ] )
			{
				whichOtherNode[ count ] = ii;
			}
		}
		++count;
	}
	
	_optimalStateMatrix3 = new OptimizedNodeStateMatrix3( numStates[ 0 ], numStates[ 1 ], numStates[ 2 ] );
	int stateAssignment[ 4 ];
	for (int ii = 0; ii < numStates[ 0 ]; ++ii)
	{
		stateAssignment[ 0 ] = ii;
		for (int jj = 0; jj < numStates[ 1 ]; ++jj )
		{
			stateAssignment[ 1 ] = jj;
			for (int kk = 0; kk < numStates[ 2 ]; ++kk )
			{
				stateAssignment[ 2 ] = kk;
				float best_score = -1234;
				int best_state = -1;
				
				for (int ll = 0; ll < _maxNodeStates; ++ll )
				{
					stateAssignment[ 3 ] = ll;
					
					float score = _theNodeScoreVect->getNodeScore( ll ) ;
					for (int mm = 0; mm < 3; ++mm)
					{
						score += edges[ mm ]->getEdgeScore( _index, ll, 
							otherNodeForEdge[ mm ], stateAssignment[ whichOtherNode[ mm ] ] );
					}
					for (int mm = 0; mm < numd3edges; ++mm)
					{
						score += d3hes[mm]->getScoreGivenSetVertexOrdering(
							stateAssignment[ d3whichVert[mm].first],
							stateAssignment[ d3whichVert[mm].second], 
							ll );
					}
					if ( degree4hyperedge )
					{
						score += degree4hyperedge->getScoreGivenSetVertexOrdering( ii, jj, kk, ll );
					}
					
					if (ll == 0 || score > best_score )
					{
						best_score = score;
						best_state = ll;
					}
				}
				_optimalStateMatrix3->setOptimumNodeState( ii, jj, kk, best_state );
				neighbors_d3he->addToScoreGivenSetVertexOrdering( ii, jj, kk, best_score );
			}
		}
	}
	
	for (int ii = 0; ii < 3; ++ii)
	{
		_nodesDeterminingOptimalState[ ii ] = nodes_sorted[ ii ];
	}
	
	if (degree4hyperedge) delete degree4hyperedge;
	for ( int ii = 0; ii < numd3edges; ++ii )
	{
		delete d3hes[ ii ];
	}
	//for (int ii = 0; ii < 3; ++ii)
	//{
	//	delete edges[ii];
	//}
				
}

void
MoveableNode::eliminateSingleton()
{
	//Find the optimal state for this node
	
	//std::cerr << "Entering Eliminate Singelton" << std::endl;

	_eliminationOrder = 3;
	NodeAndEdgeManager* theNaEManager = NodeAndEdgeManager::getInstance();
	_optimalState = (NodeState) 0;
	_optimalScore = _theNodeScoreVect->getNodeScore(_optimalState);
	for (NodeState nsindex = (NodeState) 1; nsindex < _maxNodeStates; nsindex++)
	{
		if (NodeAndEdgeManager::betterThan(_theNodeScoreVect->getNodeScore(nsindex), _optimalScore))
		{
			_optimalState = nsindex;
			_optimalScore = _theNodeScoreVect->getNodeScore(nsindex);
		}
	}

	theNaEManager->haveReportedOptimalNetworkScoreForSingleton(_optimalScore);
	theNaEManager->noteEffort( _maxNodeStates );
	//std::cerr << "Exiting Eliminate Singelton" << std::endl;
}

void
MoveableNode::eliminateOneStateVertex()
{
	assert( _maxNodeStates == 1 );
	NodeAndEdgeManager* theNaEManager = NodeAndEdgeManager::getInstance();

	for (std::list< EdgeBetweenMoveableNodes* >::iterator iter = _edgeList.begin();
		iter != _edgeList.end(); ++iter)
	{
		MoveableNode* other = (*iter)->getPointerToOtherNode( _index );
		//std::cerr << "Eliminate one state vertes -- d2edge: other: " << other->_index << std::endl;
		for (int ii = 0; ii < other->_maxNodeStates; ++ii)
		{
			
			other->_theNodeScoreVect->addToNodeScore( ii, (*iter)->getEdgeScore( _index, 0, other->_index, ii ));
		}
		theNaEManager->noteEffort( other->_maxNodeStates );
	}
	
	for (std::list< Degree3Hyperedge* >::iterator iter = _degree3hyperedges.begin();
		iter != _degree3hyperedges.end(); ++iter )
	{
		int otherIndices[ 2 ] = {-1, -1};
		int countSeen = 0;
		for (int ii = 0; ii < 3; ++ii)
		{
			if ( (*iter)->getVertexIndex( ii ) != _index )
			{
				otherIndices[ countSeen ] = (*iter)->getVertexIndex( ii );
				++countSeen;
			}
		}
		//std::cerr << "Eliminate one state vertes -- d3edge: others " << otherIndices[ 0 ] << " " << otherIndices[ 1 ] << std::endl;
		EdgeBetweenMoveableNodes* edge = theNaEManager->getEdgeBetweenMoveableNodes( otherIndices[ 0 ], otherIndices[ 1 ]);
		int otherNumStates[ 2 ];
		for (int ii = 0; ii < 2; ++ii)
		{
			otherNumStates[ ii ] = theNaEManager->getNumStatesForNode( otherIndices[ ii ] );
		}
		(*iter)->setVertexOrder( otherIndices[ 0 ], otherIndices[ 1 ], _index );
		for (int ii = 0; ii < otherNumStates[ 0 ]; ++ii)
		{
			for (int jj = 0; jj < otherNumStates[ 1 ]; ++jj )
			{
				edge->addEdgeScore( ii, jj, (*iter)->getScoreGivenSetVertexOrdering( ii, jj, 0 ) );
			}
		}
		theNaEManager->noteEffort( otherNumStates[ 0 ] * otherNumStates[ 1 ] );
	}
	
	for (std::list< Degree4Hyperedge* >::iterator iter = _degree4hyperedges.begin();
		iter != _degree4hyperedges.end(); ++iter )
	{
		int otherIndices[ 3 ] = {-1, -1, -1};
		int countSeen = 0;
		for (int ii = 0; ii < 4; ++ii)
		{
			if ( (*iter)->getVertexIndex( ii ) != _index )
			{
				otherIndices[ countSeen ] = (*iter)->getVertexIndex( ii );
				++countSeen;
			}
		}
		//std::cerr << "Eliminate one state vertes -- d4edge: others " << otherIndices[ 0 ] << " " << otherIndices[ 1 ] << " " << otherIndices[ 2 ] << std::endl;
		Degree3Hyperedge * edge = theNaEManager->getD3Edge( otherIndices[ 0 ], otherIndices[ 1 ], otherIndices[ 2 ]);
		edge->setNaturalVertexOrder();
		int otherNumStates[ 3 ];
		for (int ii = 0; ii < 3; ++ii)
		{
			otherNumStates[ ii ] = theNaEManager->getNumStatesForNode( otherIndices[ ii ] );
		}
		(*iter)->setVertexOrder( otherIndices[ 0 ], otherIndices[ 1 ], otherIndices[ 2 ], _index );
		for (int ii = 0; ii < otherNumStates[ 0 ]; ++ii)
		{
			for (int jj = 0; jj < otherNumStates[ 1 ]; ++jj )
			{
				for (int kk = 0; kk < otherNumStates[ 2 ]; ++kk)
				{
					edge->addToScoreGivenSetVertexOrdering( ii, jj, kk, (*iter)->getScoreGivenSetVertexOrdering( ii, jj, kk,  0 ) );
				}
			}
		}
		theNaEManager->noteEffort( otherNumStates[ 0 ] * otherNumStates[ 1 ] * otherNumStates[ 2 ]);
	}
	
	for (std::list< Degree3Hyperedge* >::iterator iter = _degree3hyperedges.begin();
		iter != _degree3hyperedges.end(); )
	{
		std::list< Degree3Hyperedge* >::iterator nextiter = iter;
		++nextiter;
		delete *iter; iter = nextiter;
	}
	for (std::list< Degree4Hyperedge* >::iterator iter = _degree4hyperedges.begin();
		iter != _degree4hyperedges.end(); )
	{
		std::list< Degree4Hyperedge* >::iterator nextiter = iter;
		++nextiter;
		delete *iter; iter = nextiter;
	}

	_optimalState = 0;
	_optimalScore = _theNodeScoreVect->getNodeScore( 0 );
	theNaEManager->haveReportedOptimalNetworkScoreForSingleton( _optimalScore );
	//std::cerr << "Exiting Eliminate One State Vertex" << std::endl;
}


void 				MoveableNode::beNotifiedDependencyEliminated(int i)
{
	//Find the edge that connects this node to node i and drop it from the _edgeList.
	//If the dropping of that edge means that this node can now be eliminated, add this
	//edge to the elimination queue.

	//Calling this node "node j", and claiming node j has e_j edges currently attached to it,
	//This function runs in O(e_j) time.

	//TEST
	//cerr << "beNotifiedDependencyEliminated being called on node " << _index << " by node " << i << endl;

	NodeAndEdgeManager* theNaEManager = NodeAndEdgeManager::getInstance();

	std::list<EdgeBetweenMoveableNodes*>::iterator iter = _edgeList.begin();
	int otherNodeIndex;
	bool notFound = true;
	while ((notFound) && (iter != _edgeList.end()))
	{
		//TEST
		//cerr << "Edge in list " << (*iter) << endl;
		otherNodeIndex = (*iter)->getIndexOfOtherNode(_index);
		if (otherNodeIndex == i)
		{
			std::list<EdgeBetweenMoveableNodes*>::iterator iter2 = iter;
			++iter;
			_edgeList.erase(iter2);
//			iter.drop();
			notFound = false;
		}
		else
			iter++;
	}
	if (notFound)
	{
		std::cerr << "Edge between this node and the other node (" << _index << " and " << i << " ) not found\n.";
		std::cerr << "Critical Error: Graph incorrectly initialized!" << std::endl;
		exit(1);
	}

	if (_edgeList.size() == 1)
	{
		QueueOfToBeEliminatedNodes* q = QueueOfToBeEliminatedNodes::getInstance();
		q->addNodeIndexForTreeReduction(_index);
	}
	else if ((theNaEManager->cycleReductionHasBegun()) && (_edgeList.size() == 2))
	{
		QueueOfToBeEliminatedNodes* q = QueueOfToBeEliminatedNodes::getInstance();
		q->addNodeIndexForCycleReduction(_index);
	}
	else if ( theNaEManager->safeS3ReductionHasBegun() &&
		(_edgeList.size() == 3) &&
		theNaEManager->considerNodeForSafeS3Reduction( _index ))
	{
		QueueOfToBeEliminatedNodes* q = QueueOfToBeEliminatedNodes::getInstance();
		q->addNodeIndexForCycleReduction(_index);
	}
	return;
}
NodeState			MoveableNode::getNumberOfPossibleStates() {return _maxNodeStates;}
void 				MoveableNode::setNumberOfPossibleStates(NodeState maxStates) {
	_maxNodeStates = maxStates; 
	if (_theNodeScoreVect != NULL) 
		delete _theNodeScoreVect;
	_theNodeScoreVect = new NodeScoreVector(_maxNodeStates);
	return;
}
int  				MoveableNode::getNodeIndex() {return _index;}
void				MoveableNode::setNodeIndex(int i) {_index = i; return;}
bool 				MoveableNode::getIsEliminated() {return _isEliminated;}
void 				MoveableNode::setIsEliminated(bool eliminationState) {_isEliminated = eliminationState; return;}
void 				MoveableNode::combineNodeScoreWithOptEdgeAndNodeScore(NodeScoreVector& optNodeAndEdgeScoreVect)
{
	for (NodeState stateInd = (NodeState) 0; stateInd < _maxNodeStates; stateInd++)
		_theNodeScoreVect->addToNodeScore(stateInd, optNodeAndEdgeScoreVect.getNodeScore(stateInd));
	return;
}
void 				MoveableNode::addEdge(EdgeBetweenMoveableNodes* newEdge) {_edgeList.push_back(newEdge); return;}
NodeState			MoveableNode::getNodeStateInOptimalNetworkConf(std::vector<NodeState>& theOptNetworkState)
{
	NodeState firstNodesState;
	NodeState secondNodesState;

	if (_eliminationOrder == 0)
	{
		std::cerr << "Node not eliminated, internal error in calling of getNodeStateInOptimalNetworkConf.  Node index: " << _index << endl;
		return (NodeState) -1;
	}
	else if (_eliminationOrder == 1)
	{
		//Eliminated through tree reduction
		firstNodesState = theOptNetworkState[_nodeIndexDetStateFromTreeReduction];
		return _optimalStateVector->getOptimumNodeState(firstNodesState);
	}
	else if (_eliminationOrder == 2)
	{
		//Eliminated through cycle reduction
		firstNodesState  = theOptNetworkState[_firstNodeIndexDetStateFromCycleReduction];
		secondNodesState = theOptNetworkState[_secondNodeIndexDetStateFromCycleReduction];
		return _optimalStateMatrix->getOptimumNodeState(firstNodesState,secondNodesState);
	}
	else if (_eliminationOrder == 3)
	{
		return _optimalState;
	}
	else if (_eliminationOrder == 4 )
	{
		return _optimalStateMatrix3->getOptimumNodeState( 
			theOptNetworkState[ _nodesDeterminingOptimalState[ 0 ]],
			theOptNetworkState[ _nodesDeterminingOptimalState[ 1 ]],
			theOptNetworkState[ _nodesDeterminingOptimalState[ 2 ]]);
	}
	else 
	{
		std::cerr << "CRITICAL ERROR: Unrecognized elimination order in getNodeStateInOptimalNetworkConf: " << _eliminationOrder << endl;
		return -1;
	}
}


int					MoveableNode::getNumberOfEdges() {return _edgeList.size();}

std::list<EdgeBetweenMoveableNodes* >* MoveableNode::getEdgeList() {return &_edgeList;}

bool MoveableNode::getHasAnyDegree3Hyperedges() const
{
	return _hasAnyDegree3Hyperedges;
}
bool MoveableNode::hasD3Edge( int fn, int sn, int tn )
{
	for (std::list< Degree3Hyperedge* >::iterator iter = _degree3hyperedges.begin();
		iter != _degree3hyperedges.end(); ++iter )
	{
		if ( (*iter)->sameEdge( fn, sn, tn ) )
		{
			return true;
		}
	}
	return false;
}

Degree3Hyperedge * 
MoveableNode::getD3Edge( int fn, int sn, int tn )
{
	for (std::list< Degree3Hyperedge* >::iterator iter = _degree3hyperedges.begin();
		iter != _degree3hyperedges.end(); ++iter )
	{
		if ( (*iter)->sameEdge( fn, sn, tn ) )
		{
			return *iter;
		}
	}
	return 0;
}

bool MoveableNode::getHasAnyDegree4Hyperedges() const
{
	return _hasAnyDegree4Hyperedges;
}

std::list< Degree3Hyperedge * > &  MoveableNode::getDegree3HyperedgeList()
{
	return _degree3hyperedges;
}

std::list< Degree4Hyperedge * > &  MoveableNode::getDegree4HyperedgeList()
{
	return _degree4hyperedges;
}

void MoveableNode::addDegree3Hyperedge( Degree3Hyperedge * edge )
{
	_hasAnyDegree3Hyperedges = true;
	_degree3hyperedges.push_back( edge );
}

void MoveableNode::addDegree4Hyperedge( Degree4Hyperedge * edge )
{
	_hasAnyDegree4Hyperedges = true;
	_degree4hyperedges.push_back( edge );
}

void MoveableNode::deleteDegreeThreeHyperedge( Degree3Hyperedge * edge )
{
	std::list< Degree3Hyperedge* >::iterator iter = std::find(
		_degree3hyperedges.begin(), 
		_degree3hyperedges.end(),
		edge);
	if (iter != _degree3hyperedges.end() )
	{
		//std::cerr << "Erasing iterator: " << *iter << " " << edge << " from node " << _index << std::endl;
		_degree3hyperedges.erase( iter );
	}
	else
	{
		std::cerr << "Error in deleting degree three hyperedge from graph!" << std::endl;
	}
	_hasAnyDegree3Hyperedges = _degree3hyperedges.size() != 0;
	
}

void MoveableNode::deleteDegreeFourHyperedge( Degree4Hyperedge * edge )
{
	std::list< Degree4Hyperedge* >::iterator iter = std::find(
		_degree4hyperedges.begin(), 
		_degree4hyperedges.end(),
		edge);
	if (iter != _degree4hyperedges.end() )
	{
		//std::cerr << "Erasing iterator: " << *iter << " " << edge << " from node " << _index << std::endl;
		_degree4hyperedges.erase( iter );
	}
	else
	{
		std::cerr << "Error in deleting degree four hyperedge from graph!" << std::endl;
	}
	_hasAnyDegree4Hyperedges = _degree4hyperedges.size() != 0;
	
}
//---------------------------------------------------------------------------------------------------------------

EdgeBetweenMoveableNodes::EdgeBetweenMoveableNodes(NodeState totStatesFirstNode, NodeState totStatesSecondNode) :
_isEliminated(false), _firstNodeIndex(-1), _secondNodeIndex(-1),
_firstNode(NULL), _secondNode(NULL)
{
	_theEdgeScoreMatrix = new EdgeScoreMatrix(totStatesFirstNode, totStatesSecondNode);
}
EdgeBetweenMoveableNodes::~EdgeBetweenMoveableNodes() { delete _theEdgeScoreMatrix; _theEdgeScoreMatrix = NULL;}
	
void
EdgeBetweenMoveableNodes::setEdgeScore(
	NodeState firstNodeState,
	NodeState secondNodeState,
	double theScore
)
{
	_theEdgeScoreMatrix->setEdgeScore(firstNodeState, secondNodeState, theScore);
	return;
}

void
EdgeBetweenMoveableNodes::addEdgeScore(
	NodeState firstNodeState,
	NodeState secondNodeState,
	double theScore
)
{
	_theEdgeScoreMatrix->addEdgeScore(firstNodeState, secondNodeState, theScore); 
	return;
}

double EdgeBetweenMoveableNodes::getEdgeScore(NodeState firstNodeState, NodeState secondNodeState) {return _theEdgeScoreMatrix->getEdgeScore(firstNodeState,secondNodeState);}
double EdgeBetweenMoveableNodes::getEdgeScore(int Node1sIndex, NodeState Node1sState, int Node2sIndex, NodeState Node2sState)
{
	//The convention, mentioned earlier, is that the state of the node with the smaller index goes first 
	//in the call to _theEdgeScoreMatrix->getEdgeScore.  The addition of this overloaded call to 
	//getEdgeScore allows cleaner code within the cycle elimination routine in MoveableNode.

	if (Node1sIndex < Node2sIndex)
		return _theEdgeScoreMatrix->getEdgeScore(Node1sState,Node2sState);
	else
		return _theEdgeScoreMatrix->getEdgeScore(Node2sState,Node1sState);
}
void   EdgeBetweenMoveableNodes::setFirstNodePtr(MoveableNode* firstNode) {_firstNode = firstNode; return;}
void   EdgeBetweenMoveableNodes::setSecondNodePtr(MoveableNode* secondNode) {_secondNode = secondNode; return;}
void   EdgeBetweenMoveableNodes::setFirstNodeIndex(int firstNodeIndex) {_firstNodeIndex = firstNodeIndex; return;}
void   EdgeBetweenMoveableNodes::setSecondNodeIndex(int secondNodeIndex) {_secondNodeIndex = secondNodeIndex; return;}
int    EdgeBetweenMoveableNodes::getFirstNodeIndex() {return _firstNodeIndex;}
int    EdgeBetweenMoveableNodes::getSecondNodeIndex() {return _secondNodeIndex;}
bool   EdgeBetweenMoveableNodes::getIsEliminated() {return _isEliminated;}
void   EdgeBetweenMoveableNodes::setIsEliminated(bool isEliminated) {_isEliminated = isEliminated; return;}
int	   EdgeBetweenMoveableNodes::getIndexOfOtherNode(int indexOfOneNode)
{
	//TEST
	//cerr << "First Node Index: " << _firstNodeIndex << ".  Second Node Index " << _secondNodeIndex << endl;

	if (_firstNodeIndex == indexOfOneNode)
		return _secondNodeIndex;
	else if (_secondNodeIndex == indexOfOneNode)
		return _firstNodeIndex;
	else
	{
		std::cerr << "getIndexOfOtherNode() called on edge(" << _firstNodeIndex << ", " << _secondNodeIndex << ") lacking the calling node: " << indexOfOneNode << endl;
		return -1;
	}
}

MoveableNode* EdgeBetweenMoveableNodes::getPointerToOtherNode(int indexOfOneNode)
{
	if (_firstNodeIndex == indexOfOneNode)
		return _secondNode;
	else if (_secondNodeIndex == indexOfOneNode)
		return _firstNode;
	else
	{
		std::cerr << "getPointerToOtherNode() called on edge(" << _firstNodeIndex << ", " << _secondNodeIndex << ") lacking the calling node:" << indexOfOneNode << endl;
		return NULL;
	}
}

//------------------- Degree 3 Hyperedge ---------------------------

Degree3Hyperedge::Degree3Hyperedge( int vert1, int vert2, int vert3 )
	: 
	_inNaturalVertexOrder( false )
{
	//std::cerr << "Creating degree 3 hyperedge: " << vert1 << " " << vert2 << " " << vert3 << std::endl;

	_nodeIndices[ 0 ] = vert1;
	_nodeIndices[ 1 ] = vert2;
	_nodeIndices[ 2 ] = vert3;
	NodeAndEdgeManager * theNaEM = NodeAndEdgeManager::getInstance();
	for (int ii = 0; ii < 3; ++ii )
	{
		_moveableNode[ ii ] = theNaEM->getMoveableNode( _nodeIndices[ ii ] );
		_moveableNode[ ii ]->addDegree3Hyperedge( this );
		_numStates[ ii ] = _moveableNode[ ii ]->getNumberOfPossibleStates( );
		//std::cerr << _numStates[ ii ] << " ";
	}
	//std::cerr << std::endl;
	setNaturalVertexOrder();
	_total_state_combos = _indexMultipliers[ 0 ] * _numStates[ 0 ];
	//std::cerr << "_total_state_combos " << _total_state_combos << std::endl;
	_scores = new float[ _total_state_combos ];
	for (int ii = 0; ii < _total_state_combos; ++ii )
	{
		_scores[ ii ] = 0;
	}
	
		
}

Degree3Hyperedge::~Degree3Hyperedge()
{
	delete [] _scores;
	for (int ii = 0; ii < 3; ++ii)
	{
		_moveableNode[ ii ]->deleteDegreeThreeHyperedge( this );
	}	
}

bool Degree3Hyperedge::sameEdge(int fn, int sn, int tn ) const
{
	return ( _nodeIndices[ 0 ] == fn &&
		_nodeIndices[ 1 ] == sn &&
		_nodeIndices[ 2 ] == tn );
}

int Degree3Hyperedge::getVertexIndex( int ind )
{
	assert( ind >= 0 && ind < 3 );
	return _nodeIndices[ ind ];
}

void Degree3Hyperedge::setScore( int state1, int state2, int state3, float energy)
{
	setNaturalVertexOrder();
	_scores[ getIndex( state1, state2, state3) ] = energy;
}

void Degree3Hyperedge::setVertexOrder( int vert1, int vert2, int vert3)
{
	_inNaturalVertexOrder = true;
	//std::cerr << "D3HE: " << _nodeIndices[ 0 ] << " " << _nodeIndices[ 1 ] << " " << _nodeIndices[ 2 ] << std::endl;
	//std::cerr << "setVertexOrder: " << vert1 << " " << vert2 << " " << vert3 << " _total_state_combos= " << _total_state_combos << std::endl;
	int threeVerts[] = {vert1, vert2, vert3};
	for (int ii = 0; ii < 3; ++ii )
	{
		_indexMultipliers[ ii ] = 1;
		for (int jj = 2; jj >= 0; --jj )
		{
			if ( threeVerts[ ii ] != _nodeIndices[ jj ] )
			{
				//std::cerr << "Degree3Hyperedge::setVertexOrder _indexMultipliers[ "<< ii <<" ] *= _numStates[ "<< jj << "]; _total_state_combos=" << _total_state_combos << std::endl;
				_indexMultipliers[ ii ] *= _numStates[ jj ];
			}
			else
			{
				_inNaturalVertexOrder &= (ii == jj);
				break;
			}
		}
	}
	
}

void Degree3Hyperedge::setNaturalVertexOrder() const
{
	if ( _inNaturalVertexOrder ) return;
	_indexMultipliers[ 2 ] = 1;
	_indexMultipliers[ 1 ] = _numStates[ 2 ];
	_indexMultipliers[ 0 ] = _numStates[ 2 ] * _numStates[ 1 ];
	_inNaturalVertexOrder = true;
}

float Degree3Hyperedge::getScore( int state1, int state2, int state3) const
{
	setNaturalVertexOrder();
	return _scores[ getIndex( state1, state2, state3) ];
}

float Degree3Hyperedge::getScoreGivenSetVertexOrdering( int state_vert1, int state_vert2, int state_vert3) const
{
	return _scores[ getIndex( state_vert1, state_vert2, state_vert3) ];
}

void Degree3Hyperedge::addToScoreGivenSetVertexOrdering( 
	int state_vert1, 
	int state_vert2, 
	int state_vert3,
	float score)
{
	_scores[ getIndex( state_vert1, state_vert2, state_vert3) ] += score;
}

int Degree3Hyperedge::getIndex( int state1, int state2, int state3) const
{
	int index = state1 * _indexMultipliers[ 0 ] + state2 * _indexMultipliers[ 1 ] + state3 * _indexMultipliers[ 2 ];
	assert( 0 <= index && index < _total_state_combos );
	return index;
}
//------------------- Degree 4 Hyperedge ---------------------------

Degree4Hyperedge::Degree4Hyperedge( int vert1, int vert2, int vert3, int vert4 )
	: 
	_inNaturalVertexOrder( false )
{
	_nodeIndices[ 0 ] = vert1;
	_nodeIndices[ 1 ] = vert2;
	_nodeIndices[ 2 ] = vert3;
	_nodeIndices[ 3 ] = vert4;
	NodeAndEdgeManager * theNaEM = NodeAndEdgeManager::getInstance();
	for (int ii = 0; ii < 4; ++ii )
	{
		_moveableNode[ ii ] = theNaEM->getMoveableNode( _nodeIndices[ ii ] );
		_moveableNode[ ii ]->addDegree4Hyperedge( this );
		_numStates[ ii ] = _moveableNode[ ii ]->getNumberOfPossibleStates( );
	}
	setNaturalVertexOrder();
	_total_state_combos = _indexMultipliers[ 0 ] * _numStates[ 0 ];
	_scores = new float[ _total_state_combos ];
	for (int ii = 0; ii < _total_state_combos; ++ii )
	{
		_scores[ ii ] = 0;
	}
	
	//std::cerr << "Creating degree 4 hyperedge!" << std::endl;	
}

Degree4Hyperedge::~Degree4Hyperedge()
{
	delete [] _scores;
	for (int ii = 0; ii < 4; ++ii)
	{
		_moveableNode[ ii ]->deleteDegreeFourHyperedge( this );
	}	
}

int Degree4Hyperedge::getVertexIndex( int ind )
{
	assert( ind >= 0 && ind < 4 );
	return _nodeIndices[ ind ];
}

void Degree4Hyperedge::setScore( int state1, int state2, int state3, int state4, float energy)
{
	setNaturalVertexOrder();
	_scores[ getIndex( state1, state2, state3, state4) ] = energy;
}

void Degree4Hyperedge::setVertexOrder( int vert1, int vert2, int vert3, int vert4 )
{
	_inNaturalVertexOrder = true;

	int fourVerts[] = {vert1, vert2, vert3, vert4};
	for (int ii = 0; ii < 4; ++ii )
	{
		_indexMultipliers[ ii ] = 1;
		for (int jj = 3; jj >= 0; --jj )
		{
			if ( fourVerts[ ii ] != _nodeIndices[ jj ] )
			{
				_indexMultipliers[ ii ] *= _numStates[ jj ];
			}
			else
			{
				_inNaturalVertexOrder &= (ii == jj);
				break;
			}
		}
	}
	
}

void Degree4Hyperedge::setNaturalVertexOrder() const
{
	if ( _inNaturalVertexOrder ) return;
	_indexMultipliers[ 3 ] = 1;
	_indexMultipliers[ 2 ] = _numStates[ 3 ];
	_indexMultipliers[ 1 ] = _indexMultipliers[ 2 ] * _numStates[ 2 ];
	_indexMultipliers[ 0 ] = _indexMultipliers[ 1 ] * _numStates[ 1 ];
	_inNaturalVertexOrder = true;
}

float Degree4Hyperedge::getScore( int state1, int state2, int state3, int state4) const
{
	setNaturalVertexOrder();
	return _scores[ getIndex( state1, state2, state3, state4) ];
}

float Degree4Hyperedge::getScoreGivenSetVertexOrdering( int state_vert1, int state_vert2, int state_vert3, int state_vert4) const
{
	return _scores[ getIndex( state_vert1, state_vert2, state_vert3, state_vert4) ];
}

int Degree4Hyperedge::getIndex( int state1, int state2, int state3, int state4) const
{
	int index = state1 * _indexMultipliers[ 0 ] + 
		state2 * _indexMultipliers[ 1 ] + 
		state3 * _indexMultipliers[ 2 ] +
		state4 * _indexMultipliers[ 3 ];
	assert( 0 <= index && index < _total_state_combos );
	return index;
}


//---------------------------------------------------------------------------------------------------------------

NodeAndEdgeManager* NodeAndEdgeManager::_theNaEManager = NULL;
bool NodeAndEdgeManager::_instanceFlag = false;

bool NodeAndEdgeManager::betterThan(double first, double second) {return (first > second);} // Semantics of "optimal" localized to the NaEManager.

NodeAndEdgeManager::NodeAndEdgeManager() 
: 
	_theMNVector(0), _numNodes(-1), _numEdges(-1), _numUnEliminatedNodes(-1),
	_cycleReductionBegun(false), _optimalSolutionFound(false), _optimalNetworkScore(0),
	_indFirstNodeLastEdgeLookup(-1), _indSecondNodeLastEdgeLookup(-1), 
	_ptrToEdgeOfLastEdgeLookup(NULL), _xyz( 0 ), _timeLimitExists( false ),
	_timeLimit( -1 ), _connectivity( 0 ), _effortCounter( 0 )
	
{ 
	//TEST
	//cerr << "Node and Edge Manager Constructor Succeeded" << endl;
}
void NodeAndEdgeManager::initiateTreeReduction()
{
	//Touch every node once.  Put all leaf nodes in the elimination queue.

	//TEST
	//cerr << "Initiating Tree Reduction" << endl;

	QueueOfToBeEliminatedNodes* the2BElimQ = QueueOfToBeEliminatedNodes::getInstance();
	for (int i=0; i<_numNodes; i++)
	{
		if (_theMNVector[i]->getNumberOfPossibleStates() == 1)
		{
			_theMNVector[i]->eliminate();
		}
	}
	for (int i=0; i<_numNodes; i++)
	{
		if (_theMNVector[i]->getNumberOfEdges() == 0)
		{
			_theMNVector[i]->eliminate();
		}
		else if (_theMNVector[i]->getNumberOfEdges() == 1)
		{
			the2BElimQ->addNodeIndexForTreeReduction(i);
		}
	}

	return;
}
void NodeAndEdgeManager::initiateCycleReduction()
{
	//Touch every node once.  Put all nodes that have only two edges into the elimination queue.

	//TEST
	//cerr << "Initiating cycle reduction" << endl;
	QueueOfToBeEliminatedNodes* the2BElimQ = QueueOfToBeEliminatedNodes::getInstance();
	for (int i=0;i<_numNodes;i++)
	{
		if (_theMNVector[i]->getNumberOfEdges() == 2)
			the2BElimQ->addNodeIndexForCycleReduction(i);
	}
	_cycleReductionBegun = true;
	return;
}

void NodeAndEdgeManager::initiateSafeS3Reduction()
{
	QueueOfToBeEliminatedNodes* the2BElimQ = QueueOfToBeEliminatedNodes::getInstance();
	for (int i=0;i<_numNodes;i++)
	{
		if ( _vertexDegrees[ i ] == 3 && considerNodeForSafeS3Reduction( i ) )
		{
			the2BElimQ->addNodeIndexForSafeS3Reduction(i);
		}
	}
	_safeS3ReductionBegun = true;
	return;
}

void NodeAndEdgeManager::initiateUnsafeS3Reduction()
{
	QueueOfToBeEliminatedNodes* the2BElimQ = QueueOfToBeEliminatedNodes::getInstance();
	for (int i=0;i<_numNodes;i++)
	{
		if ( _vertexDegrees[ i ] == 3 )
		{
			the2BElimQ->addNodeIndexForUnsafeS3Reduction(i);
		}
	}
	_unsafeS3ReductionBegun = true;
	return;
}

bool NodeAndEdgeManager::considerNodeForSafeS3Reduction( int node )
{
	if ( _vertexDegrees[ node ] != 3 ) return false;

	int neighbors[ 3 ];
	int countNeighbors = 0;
	for ( int ii = 0; ii < _numNodes; ++ii)
	{
		if ( _connectivity[ node ][ ii ] )
		{
			neighbors[ countNeighbors ] = ii;
			++countNeighbors;
		}
	}
	
	//arnborg & proskurowski, 86
	bool foundAdjacentNeighbors = false;
	for (int ii = 0; ii < 3; ++ii)
	{
		for (int jj = ii+1; jj < 3; ++jj )
		{
			if ( _connectivity[ neighbors[ ii ]][ neighbors[ jj ]] )
			{
				foundAdjacentNeighbors = true;
				break;
			}
		}
		if ( foundAdjacentNeighbors ) break;
	}
	
	if ( foundAdjacentNeighbors ) return true;

	//look for two graphs, C' and C'' from AP86	
	if ( matchesAP86_CPrime( node, neighbors ) ) return true;
	if ( matchesAP86_CDoublePrime( node, neighbors ) ) return true;
	return false;
}

bool NodeAndEdgeManager::matchesAP86_CPrime( int node, int neighbors[3] )
{
	//C', ask if node is isomorphic to vertex v (or u or w but not vertex x ) from figure 3
	int vertexX = -1;
	for (int ii = 0; ii < 3; ++ii)
	{
		if ( vertexX == -1 )
		{
			if (_vertexDegrees[ neighbors[ ii ]] == 3 )
			{
				vertexX = ii;
			}
			else if ( _vertexDegrees[ neighbors[ ii ]] <= 3 )
			{
				return false;
			}
		}
		else if ( _vertexDegrees[ neighbors[ ii ]] <= 3 )
		{
			return false;
		}
	}

	// vertexX from figure has degree3; other two vertices have degree no less than 4
	int xneighbors[ 2 ];
	int countXneighbors = 0;
	for (int ii = 0; ii < _numNodes; ++ii )
	{
		if (ii == node ) continue;
		if ( _connectivity[ neighbors[ vertexX ]][ ii ] )
		{
			xneighbors[ countXneighbors ] = ii;
			++countXneighbors;
		}
	}
	//other two neighbors of node must be adjacent to x's neighbors
	for (int ii = 0; ii < 3; ++ii)
	{
		if ( ii == vertexX ) continue;
		for ( int jj = 0; jj < 2; ++jj )
		{
			if ( ! _connectivity[ neighbors[ ii ]][ xneighbors[ jj ]] )
				return false;
		}
	}
	return true;
}

bool NodeAndEdgeManager::matchesAP86_CDoublePrime( int node, int neighbors[3] )
{
	//C'', ask if node is isomorphic to vertex v (or u) from figure 3
	
	//look for a neighbor of node's first neighbor that is adjacent to node's
	//other two neighbors
	for (int ii = 0; ii < _numNodes; ++ii)
	{
		if ( ii == node ) continue;
		if ( _connectivity[ neighbors[ 0 ]][ ii ] )
		{
			if ( _connectivity[ ii ][ neighbors[ 1 ] ] &&
				_connectivity[ ii ][ neighbors[ 2 ] ])
			{
				return true;
			}
		}
	}
	return false;
}

	

bool NodeAndEdgeManager::bruteForceIrreducibleSubgraph()
{
	//remainingNodesNewInd[i] = -1 if node i has been eliminated - otherwise, it's k-1, where 
	//there are k-1 uneliminated nodes with a smaller index.
	//(I use "k-1" in both places above to emphasise the face that I'm indexing by zero).
	//remainingNodesOldInd[k] = i where i is the (k+1)st smallest-indexed-uneliminated node.
	//I iterate across all possible states of the irreducible network in lexographical order, storing
	//the best network configuration encountered and it's associated score.
	//Having concluded, this sets the _optimalSolutionFound flag to true and saves the optimal network
	//state and scores in their respective variables.

	//For an irreducible subgraph of k nodes with the ith node having n_i configurations, and e edges
	//this algorithm runs in O( PI(i=1,k, n_i) * (k + e)) time.

	//This brute forcing will also work if no optimization has been done.

	//returns true if the projected time spent in brute force enumeration exceeds the requested limit
	//and computation is abandoned, returns false otherwise.

	//TEST
	if (_xyz->outputNotice() ) 
	{
		std::cerr << " Beginning begining brute force enumeration on graph with ";
		std::cerr << _numUnEliminatedNodes << " nodes" << endl;
	}

	std::vector<int> remainingNodesNewInd;
	remainingNodesNewInd.assign(_numNodes, -1);
	std::vector<int> remainingNodesOldInd;
	remainingNodesOldInd.assign(_numUnEliminatedNodes,-1);
	std::vector<MoveableNode*> remNodePtrVect;
	remNodePtrVect.assign(_numUnEliminatedNodes, static_cast<MoveableNode*> (0));
	std::vector<NodeState> remNodeMaxStates;
	remNodeMaxStates.assign(_numUnEliminatedNodes,(NodeState) -1);
	std::list<EdgeBetweenMoveableNodes*> remainingEdges;

	int i;
	int countUnEliminated = 0;
	for (i = 0; i < _numNodes; i++)
	{
		if (!_theMNVector[i]->getIsEliminated())
		{
			if (countUnEliminated >= _numUnEliminatedNodes)
			{
				std::cerr << "Error in counting numUnEliminatedNodes - iterated past end of remaining nodes vector. " << endl;
			}
			else
			{
				remainingNodesNewInd[i] = countUnEliminated;
				remainingNodesOldInd[countUnEliminated] = i;
				remNodePtrVect[countUnEliminated] = _theMNVector[i];
			}
			countUnEliminated++;

		}
	}
	if (countUnEliminated != _numUnEliminatedNodes)
	{
		std::cerr << "Error in counting numUnEliminatedNodes - found " << countUnEliminated << " while expecting " << _numUnEliminatedNodes << "." << endl;
	}

	std::list<EdgeBetweenMoveableNodes*>::iterator edgePtrIter = _theEBMNList.begin();
	int firstNodeInd, secondNodeInd;
	while(edgePtrIter != _theEBMNList.end())
	{
		if (! (*edgePtrIter)->getIsEliminated())
		{
			remainingEdges.push_back(*edgePtrIter);
			firstNodeInd  = (*edgePtrIter)->getFirstNodeIndex();
			secondNodeInd = (*edgePtrIter)->getSecondNodeIndex();
			(*edgePtrIter)->setFirstNodeIndex(remainingNodesNewInd[firstNodeInd]);
			(*edgePtrIter)->setSecondNodeIndex(remainingNodesNewInd[secondNodeInd]);

		}
		edgePtrIter++;
	}

	double runningProduct = 1;
	for (i = 0; i < _numUnEliminatedNodes; i++)
	{
		remNodeMaxStates[i] = _theMNVector[remainingNodesOldInd[i]]->getNumberOfPossibleStates();
		runningProduct *= remNodeMaxStates[i];
	}
	noteEffort( runningProduct );

	if (_xyz->outputNotice() ) 
	{
		std::cerr << " Brute Force Enumeration: ";
		std::cerr << runningProduct << " network states requiring examination." << std::endl;
	}

	int one_percent_complete = static_cast< int > (runningProduct / 100 );
	int first_stopping_point = one_percent_complete < 100000 ? one_percent_complete : 100000;
	int count_cycles = 0;
	int percent_complete = 0;
	double time_for_one_percent = 0;
	
	std::vector<NodeState> currNetworkState;
	currNetworkState.assign(_numUnEliminatedNodes, (NodeState) 0);
	std::vector<NodeState> bestNetworkState;
	bestNetworkState.assign(_numUnEliminatedNodes, (NodeState) 0);

	std::list<EdgeBetweenMoveableNodes*>::iterator edgePtrIterNew = remainingEdges.begin();
	
	std::list<Degree3Hyperedge*> degree3Hyperedges;
	std::list<Degree4Hyperedge*> degree4Hyperedges;

	bool anyHighOrderOverlap = false;	
	std::vector< bool > nodeHasHighOrderOverlap( _numUnEliminatedNodes, false );
	std::vector< Mover* > nodesCorrespondingMover( _numUnEliminatedNodes, static_cast< Mover*> (0) );
	std::vector< std::vector< std::pair< AtomDescr, DotsForAtom * > > > atomsInHighOrderOverlap( _numUnEliminatedNodes );

	clock_t start_time = clock();
	for (int ii = 0; ii < _numUnEliminatedNodes; ++ii)
	{
		nodesCorrespondingMover[ ii ] = remNodePtrVect[ ii ]->getMover();
		if ( remNodePtrVect[ ii ]->getHasAnyHighOrderOverlap() )
		{
			//int ii_node_index = remainingNodesOldInd[ ii ];
			//std::cerr << "Node " << ii_node_index << " has high order overlap" << std::endl;
		
			anyHighOrderOverlap = true;
			nodeHasHighOrderOverlap[ ii ] = true;
			atomsInHighOrderOverlap[ ii ] = remNodePtrVect[ ii ]->getAtomsInHighOrderOverlap();
		}
	}
	
	if ( anyHighOrderOverlap )
	{
		for (int ii = 0; ii < _numUnEliminatedNodes; ++ii )
		{
			//std::cerr << "Initializing mover " << ii << std::endl;
			nodesCorrespondingMover[ ii ]->initOrientation( *_xyz );
		}
	}
	
	for (int ii = 0; ii < _numUnEliminatedNodes; ++ii )
	{
		if ( remNodePtrVect[ ii ]->getHasAnyDegree3Hyperedges() )
		{
			//std::cerr << "Collecting degree3 hyperedge for mover " << ii << std::endl;
			degree3Hyperedges.splice( degree3Hyperedges.end(), remNodePtrVect[ ii ]->getDegree3HyperedgeList() );
		}
		if ( remNodePtrVect[ ii ]->getHasAnyDegree4Hyperedges() )
		{
			//std::cerr << "Collecting degree3 hyperedge for mover " << ii << std::endl;
			degree4Hyperedges.splice( degree4Hyperedges.end(), remNodePtrVect[ ii ]->getDegree4HyperedgeList() );
		}
	}
	//std::cerr << "All d3e's collected" << std::endl;
	degree3Hyperedges.sort();
	degree3Hyperedges.unique();
	std::vector< Degree3Hyperedge* > d3he( degree3Hyperedges.size() );
	std::copy( degree3Hyperedges.begin(), degree3Hyperedges.end(), d3he.begin() );
	std::vector< std::vector< int > > d3he_verts( d3he.size() );
	for (int ii = 0; ii < d3he.size(); ++ii )
	{
		//std::cerr << "Initializing d3e: " << ii << " " << d3he[ ii ] << "..."; 
		d3he[ ii ]->setNaturalVertexOrder();
		d3he_verts[ ii ].resize( 3 );
		for (int jj = 0; jj < 3; ++jj )
		{
			d3he_verts[ ii ][ jj ] = remainingNodesNewInd[ d3he[ ii ]->getVertexIndex( jj ) ];
		}
		//std::cerr << "...done" << std::endl;
	}

	degree4Hyperedges.sort();
	degree4Hyperedges.unique();
	std::vector< Degree4Hyperedge* > d4he( degree4Hyperedges.size() );
	std::copy( degree4Hyperedges.begin(), degree4Hyperedges.end(), d4he.begin() );
	std::vector< std::vector< int > > d4he_verts( d4he.size() );
	for (int ii = 0; ii < d4he.size(); ++ii )
	{
		//std::cerr << "Initializing d4e: " << ii << " " << d4he[ ii ] << "..."; 
		d4he[ ii ]->setNaturalVertexOrder();
		d4he_verts[ ii ].resize( 4 );
		for (int jj = 0; jj < 4; ++jj )
		{
			d4he_verts[ ii ][ jj ] = remainingNodesNewInd[ d4he[ ii ]->getVertexIndex( jj ) ];
		}
		//std::cerr << "...done" << std::endl;
	}

	
	std::vector< EdgeBetweenMoveableNodes* > d2he( remainingEdges.size() );
	std::copy( remainingEdges.begin(), remainingEdges.end(), d2he.begin() );
	std::vector< std::pair< int, int > > d2he_verts( remainingEdges.size() );
	std::vector< float > d2he_scores( remainingEdges.size() );
	for (int ii = 0; ii < d2he.size(); ++ii )
	{
		//d2he_verts[ ii ].resize( 2 );
		d2he_verts[ ii ].first = d2he[ ii ]->getFirstNodeIndex();
		d2he_verts[ ii ].second = d2he[ ii ]->getSecondNodeIndex();
	}
	
	int nodeIndex = 0;
	bool notDone = true;
	bool firstNetworkState = true;

	double runningScore;
	double bestScore = -1;
	//std::cerr << "Beginning enumeration..." << std::endl;
	int smallestIncremented = -1;
	double penalty;
	while (notDone)
	{
		//evaluate score of network in current state
		nodeIndex = 0;
//		edgePtrIterNew.reset();
		runningScore = 0;
		while (nodeIndex < _numUnEliminatedNodes)
		{
			//TEST
			//cerr << "Scoring... Node Index = " << nodeIndex << endl;
			runningScore += remNodePtrVect[nodeIndex]->getNodeStateScore(currNetworkState[nodeIndex]);
			if ( nodeHasHighOrderOverlap[ nodeIndex ] )
			{
				runningScore += _xyz->determineScoreForMover(
					nodesCorrespondingMover[ nodeIndex ],
					atomsInHighOrderOverlap[ nodeIndex ],
					penalty );
			}
					
			nodeIndex++;
		}
		for (int ii = 0; ii < d2he.size(); ++ii)
		{
			if ( d2he_verts[ ii ].second < smallestIncremented)
			{
				runningScore += d2he_scores[ ii ];
			}
			else{
				runningScore += d2he_scores[ ii ] = d2he[ ii ]->getEdgeScore(
					currNetworkState[d2he_verts[ ii ].first],
					currNetworkState[d2he_verts[ ii ].second]);
			}
		}
		for ( int ii = 0; ii < d3he.size(); ++ii )
		{
			runningScore += d3he[ ii ]->getScore(
				currNetworkState[ d3he_verts[ ii ][ 0 ]],
				currNetworkState[ d3he_verts[ ii ][ 1 ]],
				currNetworkState[ d3he_verts[ ii ][ 2 ]]);
		}
		
		for ( int ii = 0; ii < d4he.size(); ++ii )
		{
			runningScore += d4he[ ii ]->getScore(
				currNetworkState[ d4he_verts[ ii ][ 0 ]],
				currNetworkState[ d4he_verts[ ii ][ 1 ]],
				currNetworkState[ d4he_verts[ ii ][ 2 ]],
				currNetworkState[ d4he_verts[ ii ][ 3 ]]);
		}

	
		//save the network state if it's the best score found so far
		if (firstNetworkState)
		{
			firstNetworkState = false;
			bestScore = runningScore;
		}
		else
		{
			if (NodeAndEdgeManager::betterThan(runningScore, bestScore))
			{
				bestScore = runningScore;
				for (i = 0; i < _numUnEliminatedNodes; i++)
					bestNetworkState[i] = currNetworkState[i];
			}
		}

		//increment to the next network state, if there is no next state, set notDone to false
		//using lexographical order
		for (i=(_numUnEliminatedNodes - 1); i >= 0; i--)
		{
			if (currNetworkState[i] < (remNodeMaxStates[i] - 1))
			{
				currNetworkState[i]++;
				if ( anyHighOrderOverlap )
				{
					nodesCorrespondingMover[ i ]->nextOrientation( *_xyz );
				}
				smallestIncremented = i;
				break;
			}
			else
			{
				if (i != 0)
				{
					currNetworkState[i] = 0;
					if ( anyHighOrderOverlap )
					{
						nodesCorrespondingMover[ i ]->initOrientation( *_xyz );
					}
				}
				else
				{
					notDone = false;
				}

			}
		}
		
		++count_cycles;

		if ( _timeLimitExists && count_cycles == first_stopping_point )
		{
			clock_t stoptime = clock();
			double time_remaining = ((double) stoptime - start_time ) * 
				( ((double ) (runningProduct )/ first_stopping_point) - 1 ) / CLOCKS_PER_SEC;
			if ( time_remaining > _timeLimit )
			{
				//projected running time exceeds the time limit that brute force has been given
				//quit
				int hours_remaining = static_cast< int > ( time_remaining ) / 3600;
				int minutes_remaining = static_cast< int >  ( time_remaining ) / 60  - 60 * hours_remaining;
				int seconds_remaining = static_cast< int > ( time_remaining ) - 60 * minutes_remaining - 3600 * hours_remaining;
				
				if ( _xyz->outputNotice() )
				{
					std::cerr << " Projected time remaining for this single optimization: " << hours_remaining << "h ";
					std::cerr << minutes_remaining << "m " << seconds_remaining << "s" << std::endl;
				}

				return true;
			}
			else
			{
				//projected running time comes in under the limit; keep processing
				_timeLimitExists = false;
			}
		}
		
		if ( count_cycles == one_percent_complete ) { 
			++percent_complete;
			count_cycles = 0;
			if ( percent_complete == 1 )
			{
				clock_t stoptime = clock();
				time_for_one_percent = ( static_cast< double > ( stoptime - start_time ) ) / CLOCKS_PER_SEC;
			}
			if (percent_complete  == 10 )
			{
				clock_t stoptime = clock();
				time_for_one_percent = ( static_cast< double > ( stoptime - start_time ) ) / (10 * CLOCKS_PER_SEC);
			}
			
			if ( percent_complete <= 10 || (percent_complete % 10 == 0 ) )
			{
				double time_remaining = (100 - percent_complete) * time_for_one_percent;
				int hours_remaining = static_cast< int > ( time_remaining ) / 3600;
				int minutes_remaining = static_cast< int >  ( time_remaining ) / 60  - 60 * hours_remaining;
				int seconds_remaining = static_cast< int > ( time_remaining ) - 60 * minutes_remaining - 3600 * hours_remaining;
				if ( _xyz->outputNotice() )
				{
					std::cerr << " Brute force: percent complete= " << percent_complete << "; time remaining= " << hours_remaining << "h ";
					std::cerr << minutes_remaining << "m " << seconds_remaining << "s" << std::endl;
				}
			}
			
		}
	} // end while

	//Now save what you've found
	_optimalNetworkScore += bestScore;
	if ( _xyz->outputNotice() ) std::cerr << " Exhaustive search ended. best network score = " << _optimalNetworkScore << std::endl;

	for (i=0; i < _numUnEliminatedNodes; i++)
	{
		//TEST
		//cerr << "Saving optimal solution. i = " << i ;
		//cerr << ".  remainingNodeOldInd[ " << i << "] = " << remainingNodesOldInd[i];
		//cerr << ".  bestNetworkState[ " << i << "] = " << bestNetworkState[i] << endl;
		_optimalNetworkStateVector[remainingNodesOldInd[i]] = bestNetworkState[i];
	}

	_optimalSolutionFound = true;
	return false;
}

NodeAndEdgeManager::~NodeAndEdgeManager()
{
	clear();
	_instanceFlag = false;
	_theNaEManager = NULL;
}
NodeAndEdgeManager* NodeAndEdgeManager::getInstance()
{
	if (!_instanceFlag)
	{
		_theNaEManager = new NodeAndEdgeManager();
		_instanceFlag = true;
	}
	return _theNaEManager;
}
//void NodeAndEdgeManager::NotifyNodeOfDependencyElimination(int eliminatedNodeIndex, int dependentNodeIndex);
void NodeAndEdgeManager::BeNotifiedOfEliminatedNode(int indexContractedNode) 
{
	_PartialOrderStateDeterminationStack.push_front(indexContractedNode);
	_numUnEliminatedNodes--;
	_vertexDegrees[ indexContractedNode ] = 0;
	for (int ii = 0; ii < _numNodes; ++ii )
	{
		if ( _connectivity[ indexContractedNode ][ ii ] )
		{
			--_vertexDegrees[ ii ];
			_connectivity[ indexContractedNode ][ ii ] = false;
			_connectivity[ ii ][ indexContractedNode ]= false;
		}
	}
	
	return;
}

void NodeAndEdgeManager::eliminateQueuedNodes()
{
	QueueOfToBeEliminatedNodes* the2BElimQ = QueueOfToBeEliminatedNodes::getInstance();
	int indexOfNextNode2BElim = -1;
	while (!the2BElimQ->isEmpty())
	{
		indexOfNextNode2BElim = the2BElimQ->nextNodeForReduction();
		_theMNVector[indexOfNextNode2BElim]->eliminate();
	}
}

bool NodeAndEdgeManager::computeOptimalNetworkConfiguration()
{
	//Initialize must have already been called.

	//First, trim all trees.
	//Then start reducing cycles, reducing all trees that are formed along the way.
	//If that hasn't finished off the network, then brute force the remaining subgraph.
	//Finally, determine the states of all nodes in the network.

	if ( _xyz->outputNotice() ) std::cerr << " Beginning Optimization" << endl;

	initiateTreeReduction();
	eliminateQueuedNodes();
	
	initiateCycleReduction();
	eliminateQueuedNodes();
	
	initiateSafeS3Reduction();
	eliminateQueuedNodes();
	
	initiateUnsafeS3Reduction();
	eliminateQueuedNodes();
	
	if (_numUnEliminatedNodes == 0)
	{
		//TEST
		if ( _xyz->outputNotice() ) std::cerr << " Dynamic programming succeeded to fully optimize hypergraph" << endl;
		_optimalSolutionFound = true;
		computeOptimalNetworkConfigurationAfterOptimization();
	}
	else
	{
		bool abandoned = bruteForceIrreducibleSubgraph();
		if ( abandoned ) return true;
		
		computeOptimalNetworkConfigurationAfterOptimization();
	}

	//std::cerr << "Optimal network state computed with effort of " << _effortCounter << std::endl;

	//writeOutOptimalNetworkState();
	return false;
}
void NodeAndEdgeManager::computeOptimalNetworkConfigurationAfterOptimization()
{
	// The stack which holds the partial order giving the sequence in which node states
	// can be determined for the eliminated nodes should be iterated across.  
	// The optimal network score has already been computed.

	// With k nodes in the network, this is done in O(k) time.

	std::list<int>::const_iterator partOrderStackIter = _PartialOrderStateDeterminationStack.begin();
	NodeState currNodeBestNodeState;
	while (partOrderStackIter != _PartialOrderStateDeterminationStack.end())
	{
		currNodeBestNodeState = _theMNVector[*partOrderStackIter]->getNodeStateInOptimalNetworkConf(_optimalNetworkStateVector);
		_optimalNetworkStateVector[*partOrderStackIter] = currNodeBestNodeState;
		++partOrderStackIter;
	}
	return;
}
void NodeAndEdgeManager::InitializeNetwork(GraphToHoldScores & gths)
{

	//TEST
	//cerr << "Entering InitializeNetwork(gths)" << endl;

	//Go ahead and make sure that things within the gths are legit to the extent allowable.
	int i, numNodes, firstNodeIndex, secondNodeIndex;
	NodeState j, k, firstNodeMaxStates, secondNodeMaxStates;
	numNodes = gths.getNumNodes();	

	clear();
	_numNodes = numNodes;
	_theMNVector.resize(_numNodes);
	_optimalNetworkStateVector.resize(_numNodes);
	_numUnEliminatedNodes = _numNodes;
	
	_connectivity = new bool*[ _numNodes ];
	for (int ii = 0; ii < _numNodes; ++ii )
	{
		_connectivity[ ii ] = new bool[ _numNodes ];
		for (int jj = 0; jj < _numNodes; ++jj )
		{
			_connectivity[ ii ][ jj ] = false;
		}
	}
	_vertexDegrees.resize( _numNodes );
	std::fill( _vertexDegrees.begin(), _vertexDegrees.end(), 0 );
	

	for (i = 0; i < _numNodes; i++)
	{
		firstNodeMaxStates = gths.getNumEnabledStatesForNode(i);
		if ( firstNodeMaxStates == 0 )
		{
			std::cerr << "Critical error in initialization: vertex " << i << " with 0 states." << std::endl; 
			return;
		}
		_theMNVector[i] = new MoveableNode(firstNodeMaxStates);
		_theMNVector[i]->setNodeIndex(i);
		for (j = (NodeState) 0; j < firstNodeMaxStates; j++)
		{
			_theMNVector[i]->setNodeStateScore(j, gths.getNodeScoreForState(i,j));
		}
		if ( gths.getNodeHasAnyHighOrderOverlap( i ) )
		{
			if ( _xyz->outputNotice() )
			{
				std::cerr << "InitializeNetwork:: Node " << i << " with high order overlap " << std::endl;
			}
			_theMNVector[i]->setAtomsInHighOrderOverlap( gths.getAtomsInHighOrderOverlapForNode( i ) );
		}
		_theMNVector[i]->setMover( gths.getMoverForNode( i ) );
	}

	EdgeBetweenMoveableNodes* newEdge = NULL;
	gths.setD2EIteratorAtBegining();
	while ( ! gths.getD2EIteratorAtEnd() )
	{
		firstNodeIndex = gths.getFirstIndexFocusedD2E();
		secondNodeIndex = gths.getSecondIndexFocusedD2E();
		firstNodeMaxStates = _theMNVector[firstNodeIndex]->getNumberOfPossibleStates();
		secondNodeMaxStates = _theMNVector[secondNodeIndex]->getNumberOfPossibleStates();
		newEdge = new EdgeBetweenMoveableNodes(firstNodeMaxStates, secondNodeMaxStates);
		_theEBMNList.push_back(newEdge);
		
		_connectivity[ firstNodeIndex ][ secondNodeIndex ] = true;
		_connectivity[ secondNodeIndex ][ firstNodeIndex ] = true;
		++_vertexDegrees[ firstNodeIndex ];
		++_vertexDegrees[ secondNodeIndex ];
		
		newEdge->setFirstNodeIndex(firstNodeIndex);
		newEdge->setSecondNodeIndex(secondNodeIndex);
		newEdge->setFirstNodePtr(_theMNVector[firstNodeIndex]);
		newEdge->setSecondNodePtr(_theMNVector[secondNodeIndex]);
		_theMNVector[firstNodeIndex]->addEdge(newEdge);
		_theMNVector[secondNodeIndex]->addEdge(newEdge);
		for (j = (NodeState) 0; j < firstNodeMaxStates; j++)
		{
			for (k = (NodeState) 0; k< secondNodeMaxStates; k++)
			{
				newEdge->setEdgeScore(j, k, gths.getScoreForFocusedD2E( j, k));
			}
		}
		gths.incrementD2EIterator();
	}

	//std::cerr << "Add Degree-3 Hyperedges?  How many: " << gths.getNumHyperedges() << endl;
	_xyz = gths.getAtomPositionsPointer();
	gths.setD3EIteratorAtBegining();
	while ( ! gths.getD3EIteratorAtEnd() )
	{
		int fn = gths.getFirstIndexFocusedD3E();
		int sn = gths.getSecondIndexFocusedD3E();
		int tn = gths.getThirdIndexFocusedD3E();
		Degree3Hyperedge * d3h = new Degree3Hyperedge( fn, sn, tn );
		firstNodeMaxStates = _theMNVector[fn]->getNumberOfPossibleStates();
		secondNodeMaxStates = _theMNVector[sn]->getNumberOfPossibleStates();
		int thirdNodeMaxStates = _theMNVector[ tn ]->getNumberOfPossibleStates();
		for (j = (NodeState) 0; j < firstNodeMaxStates; j++)
		{
			for (k = (NodeState) 0; k< secondNodeMaxStates; k++)
			{
				for (int l = (NodeState ) 0; l< thirdNodeMaxStates; ++l)
				{
					d3h->setScore(j, k, l, gths.getScoreForFocusedD3E( j, k, l));
				}
			}
		}
		gths.incrementD3EIterator();
	}
	
	gths.setD4EIteratorAtBegining();
	while ( ! gths.getD4EIteratorAtEnd() )
	{
		int fn = gths.getFirstIndexFocusedD4E();
		int sn = gths.getSecondIndexFocusedD4E();
		int tn = gths.getThirdIndexFocusedD4E();
		int fthn = gths.getFourthIndexFocusedD4E();
		Degree4Hyperedge * d4h = new Degree4Hyperedge( fn, sn, tn, fthn );
		firstNodeMaxStates = _theMNVector[fn]->getNumberOfPossibleStates();
		secondNodeMaxStates = _theMNVector[sn]->getNumberOfPossibleStates();
		int thirdNodeMaxStates = _theMNVector[ tn ]->getNumberOfPossibleStates();
		int fourthNodeMaxStates = _theMNVector[ fthn ]->getNumberOfPossibleStates();
		for (j = (NodeState) 0; j < firstNodeMaxStates; j++)
		{
			for (k = (NodeState) 0; k< secondNodeMaxStates; k++)
			{
				for (int l = (NodeState ) 0; l< thirdNodeMaxStates; ++l)
				{
					for (int m = (NodeState ) 0; m< fourthNodeMaxStates; ++m)
					{
						d4h->setScore(j, k, l, m, gths.getScoreForFocusedD4E( j, k, l, m));
					}
				}
			}
		}
		gths.incrementD4EIterator();
	}


	//TEST
	//std::cerr << "Network Initialized From GraphToHoldScores" << endl;
	
	return;
}

void NodeAndEdgeManager::setTimeLimit( double time_limit_in_seconds )
{
	_timeLimitExists = true;
	_timeLimit = time_limit_in_seconds;
}

double NodeAndEdgeManager::getScoreOfOptimalNetworkConfiguration()
{
	if (_optimalSolutionFound)
		return _optimalNetworkScore;
	else
	{
		computeOptimalNetworkConfiguration();
		return _optimalNetworkScore;
	}
	
}
void NodeAndEdgeManager::clear()
{
	//Deallocate everything.

	//TEST
	//std::cerr << "Entering Clear().  _numNodes = " << _numNodes << endl;

	for (int ii = 0; ii < _theMNVector.size(); ++ii)
	{
		if (_theMNVector[ii] != NULL)
		{	
			delete (_theMNVector[ii]); _theMNVector[ii] = NULL;
		}
	}
	_theMNVector.resize(0);

	std::list<EdgeBetweenMoveableNodes*>::iterator edgeIter = _theEBMNList.begin();
	while (edgeIter != _theEBMNList.end())
	{
		if ((*edgeIter) != NULL)
			delete (*edgeIter);		
		edgeIter++;
	}
	_theEBMNList.clear();	

	if ( _connectivity )
	{
		for (int ii = 0; ii < _numNodes; ++ii)
		{
			delete [] _connectivity[ ii ]; _connectivity[ ii ] = 0;
		}
		delete [] _connectivity; _connectivity = 0;
	}

	_PartialOrderStateDeterminationStack.clear();

	_optimalNetworkStateVector.resize(0);
	_ptrToEdgeOfLastEdgeLookup = NULL;
	_indFirstNodeLastEdgeLookup = -1;
	_indSecondNodeLastEdgeLookup = -1;
	_numNodes = 0;
	_numEdges = 0;
	_numUnEliminatedNodes = 0;
	_cycleReductionBegun = false;
	_safeS3ReductionBegun = false;
	_unsafeS3ReductionBegun = false;
	_optimalSolutionFound = false;
	_optimalNetworkScore = 0;
	_timeLimitExists = false;
	_timeLimit = -1.0;
	_effortCounter = 0;

	//TEST
	//cerr << "Leaving Clear()" << endl;

	return;
}
bool NodeAndEdgeManager::cycleReductionHasBegun() const
{return _cycleReductionBegun;}

bool NodeAndEdgeManager::safeS3ReductionHasBegun() const
{return _safeS3ReductionBegun;}

bool NodeAndEdgeManager::existsEdgeBetween(int firstNodeIndex, int secondNodeIndex)
{
	//Call these two nodes b and c.  Node b has e_b edges, node c has e_c edges.
	//This runs in O(min(e_b, e_c)) time.

	if ((firstNodeIndex == _indFirstNodeLastEdgeLookup) && (secondNodeIndex == _indSecondNodeLastEdgeLookup))
	{
		return (_ptrToEdgeOfLastEdgeLookup != NULL);
	}

	_indFirstNodeLastEdgeLookup  = firstNodeIndex;
	_indSecondNodeLastEdgeLookup = secondNodeIndex;
	_ptrToEdgeOfLastEdgeLookup = NULL;

	std::list<EdgeBetweenMoveableNodes* >* firstEdgeListPtr  = _theMNVector[firstNodeIndex]->getEdgeList();
	std::list<EdgeBetweenMoveableNodes* >* secondEdgeListPtr = _theMNVector[secondNodeIndex]->getEdgeList();
	int searchForNodeIndex      = -1;
	int nodeOfSearechedEdgeList = -1;

	if (firstEdgeListPtr->size() < secondEdgeListPtr->size())
	{
		std::list<EdgeBetweenMoveableNodes*>::iterator iter = (*firstEdgeListPtr).begin();
		searchForNodeIndex = secondNodeIndex;
		nodeOfSearechedEdgeList = firstNodeIndex;
		while (iter != (*firstEdgeListPtr).end())			//Duplicated code since iter can't be declared outside of the if statement (unless I make it a pointer)
		{
			if ((*iter)->getIndexOfOtherNode(nodeOfSearechedEdgeList) == searchForNodeIndex)
			{
				_ptrToEdgeOfLastEdgeLookup = *iter;
				return true;
			}
			iter++;
		}
	}
	else
	{
		std::list<EdgeBetweenMoveableNodes*>::iterator iter = (*secondEdgeListPtr).begin();
		searchForNodeIndex = firstNodeIndex;
		nodeOfSearechedEdgeList = secondNodeIndex;
		while (iter != (*secondEdgeListPtr).end())			//Duplicated code since iter can't be declared outside of the if statement (unless I make it a pointer)
		{
			if ((*iter)->getIndexOfOtherNode(nodeOfSearechedEdgeList) == searchForNodeIndex)
			{
				_ptrToEdgeOfLastEdgeLookup = *iter;
				return true;
			}
			iter++;
		}
	}
	
	return false;
}

void NodeAndEdgeManager::haveReportedOptimalNetworkScoreForSingleton(double optScore)
{
	_optimalNetworkScore += optScore;
	_optimalSolutionFound = true;
	return;
}

bool NodeAndEdgeManager::d3edgeExistsBetween( int fn, int sn, int tn )
{
	return _theMNVector[ fn ]->hasD3Edge( fn, sn, tn );
}
	
Degree3Hyperedge * NodeAndEdgeManager::getD3Edge( int fn, int sn, int tn )
{
	return _theMNVector[ fn ]->getD3Edge( fn, sn, tn );
}

EdgeBetweenMoveableNodes* NodeAndEdgeManager::getEdgeBetweenMoveableNodes(int firstNodeIndex, int secondNodeIndex)
{
	if ((firstNodeIndex == _indFirstNodeLastEdgeLookup) && (secondNodeIndex == _indSecondNodeLastEdgeLookup))
	{
		return _ptrToEdgeOfLastEdgeLookup;
	}
	else
	{
		if (existsEdgeBetween(firstNodeIndex,secondNodeIndex))
			return _ptrToEdgeOfLastEdgeLookup;
		else
			return NULL;
	}
}

EdgeBetweenMoveableNodes* NodeAndEdgeManager::addNewEdge(int firstNodeIndex, int secondNodeIndex)
{
	//This is called by MoveableNode::eliminateNodeThroughCycleReduction() if nodes b and c, connected
	//to the node being eliminated, node a, do not share an edge.  The edge is new'ed here, put in the
	//NodeAndEdgeManager _theEBMNList and the _edgeLists of both nodes b and c.  A pointer to this edge
	//is returned so that the elimination function can fill out the edgeScoreMatrix.

	//We've changed the graph - so the previous existsEdgeBetween() query is no longer valid.
	//Store the new edge instead.
	
	if (firstNodeIndex > secondNodeIndex)
	{
		std::cerr << "Critical Error: firstNodeIndex greater than secondNodeIndex in call to addNewEdge! " << endl;
	}

	MoveableNode* firstNodePtr  = _theMNVector[firstNodeIndex];
	MoveableNode* secondNodePtr = _theMNVector[secondNodeIndex];
	NodeState firstNodeMaxStates  = firstNodePtr->getNumberOfPossibleStates();
	NodeState secondNodeMaxStates = secondNodePtr->getNumberOfPossibleStates();
	EdgeBetweenMoveableNodes* newEdge = new EdgeBetweenMoveableNodes(firstNodeMaxStates,secondNodeMaxStates);
	_theEBMNList.push_back(newEdge);
	newEdge->setFirstNodeIndex(firstNodeIndex);
	newEdge->setSecondNodeIndex(secondNodeIndex);
	newEdge->setFirstNodePtr(firstNodePtr);
	newEdge->setSecondNodePtr(secondNodePtr);
	firstNodePtr->addEdge(newEdge);
	secondNodePtr->addEdge(newEdge);

	_connectivity[ firstNodeIndex ][ secondNodeIndex ] = true;
	_connectivity[ secondNodeIndex ][ firstNodeIndex ] = true;
	++_vertexDegrees[ firstNodeIndex ];
	++_vertexDegrees[ secondNodeIndex ];

	_indFirstNodeLastEdgeLookup  = firstNodeIndex;
	_indSecondNodeLastEdgeLookup = secondNodeIndex;
	_ptrToEdgeOfLastEdgeLookup   = newEdge;

	return newEdge;

}

void NodeAndEdgeManager::bruteForceOriginalGraph()
{
	//the brute force algorithm works reguardless of whether the graph is irreducible.
	bruteForceIrreducibleSubgraph();
	//writeOutOptimalNetworkState();
	return;
}

int 
NodeAndEdgeManager::getNumStatesForNode( int node ) const
{
	return _theMNVector[ node ]->getNumberOfPossibleStates();
}

void NodeAndEdgeManager::writeOutOptimalNetworkState()
{
	std::cerr << "_________________________________________" << std::endl;
	std::cerr << "Outputting optimal network configuration:" << std::endl;
	for (int i = 0; i < _numNodes; i++)
		std::cerr << "Node Index: " << i << "; Node State: " << (int) _optimalNetworkStateVector[i] << std::endl;
	std::cerr << "Optimal network score: " << _optimalNetworkScore << std::endl;
	std::cerr << "_________________________________________" << std::endl;
}

void NodeAndEdgeManager::noteEffort( double effort )
{
  //std::cerr << "Noting effort: " << effort << " to add to running tally: " << _effortCounter << std::endl;
  _effortCounter += effort;
}

//---------------------------------------------------------------------------------------------------------------

QueueOfToBeEliminatedNodes* QueueOfToBeEliminatedNodes::_theQueue = NULL;
bool QueueOfToBeEliminatedNodes::_instanceFlag = false;

QueueOfToBeEliminatedNodes* QueueOfToBeEliminatedNodes::getInstance()
{
	if (!_instanceFlag)
	{
		_theQueue = new QueueOfToBeEliminatedNodes();
		_instanceFlag = true;
	}
	return _theQueue;
}
QueueOfToBeEliminatedNodes::~QueueOfToBeEliminatedNodes()
{
	_instanceFlag = false;
	_theQueue = NULL;
}
void QueueOfToBeEliminatedNodes::addNodeIndexForTreeReduction(int nodeIndex)
{
	_QForTreeRed.push_back(nodeIndex);
	return;
}
void QueueOfToBeEliminatedNodes::addNodeIndexForCycleReduction(int nodeIndex)
{
	_QForCycleRed.push_back(nodeIndex);
	return;
}

void QueueOfToBeEliminatedNodes::addNodeIndexForSafeS3Reduction( int nodeIndex )
{
	_QForSafeS3Red.push_back( nodeIndex );
	return;
}

void QueueOfToBeEliminatedNodes::addNodeIndexForUnsafeS3Reduction( int nodeIndex )
{
	//sort by state space size; eliminate nodes with the greatest # states, first
	NodeAndEdgeManager * theNaEM = NodeAndEdgeManager::getInstance();
	int newNodeNumStates = theNaEM->getNumStatesForNode( nodeIndex );
	
	_QForUnsafeS3Red.push_front( nodeIndex );
	std::list< int >::iterator newNode = _QForUnsafeS3Red.begin();
	std::list< int >::iterator nextNode = newNode;
	++nextNode;
	
	while ( nextNode != _QForUnsafeS3Red.end() )
	{
		if ( theNaEM->getNumStatesForNode( *nextNode ) > newNodeNumStates )
		{
			*newNode = *nextNode;
			*nextNode = nodeIndex;
			newNode = nextNode;
			++nextNode;
		}
		else{
			break;
		}
		
	}
	return;
}


bool QueueOfToBeEliminatedNodes::isEmpty()
{
	return (_QForTreeRed.empty() &&
		_QForCycleRed.empty() &&
		_QForSafeS3Red.empty() &&
		_QForUnsafeS3Red.empty() );
}

int  QueueOfToBeEliminatedNodes::nextNodeForReduction() 
{	
	//Part of the algorithm is hidden here: always reduce a tree node if that option is available
	// otherwise, go and reduce a cycle node.  Some nodes may be placed in both queues.  This requires
	// the ::eliminate() function to check first if the node it is acting on hasn't already been eliminated.
	// since a node can be placed in each queue at most once, a linear time bound is preserved.

	int holdIndex;
	if (!_QForTreeRed.empty()) 
	{
		std::list<int>::iterator iter = _QForTreeRed.begin();
		holdIndex = *iter;
		_QForTreeRed.pop_front();
		return holdIndex;
	}
	else if (!_QForCycleRed.empty() )
	{
		std::list<int>::iterator iter = _QForCycleRed.begin();
		holdIndex = *iter;
		_QForCycleRed.pop_front();
		return holdIndex;
	}
	else if (!_QForSafeS3Red.empty() )
	{
		std::list<int>::iterator iter = _QForSafeS3Red.begin();
		holdIndex = *iter;
		_QForSafeS3Red.pop_front();
		return holdIndex;
	}
	else if (!_QForUnsafeS3Red.empty() )
	{
		std::list<int>::iterator iter = _QForUnsafeS3Red.begin();
		holdIndex = *iter;
		_QForUnsafeS3Red.pop_front();
		return holdIndex;
	}
	else
	{
		return -1;
	}
}
QueueOfToBeEliminatedNodes::QueueOfToBeEliminatedNodes()  {}
