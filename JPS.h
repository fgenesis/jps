#ifndef JUMP_POINT_SEARCH_H
#define JUMP_POINT_SEARCH_H

// Jump point search implementation
// Based on the paper http://users.cecs.anu.edu.au/~dharabor/data/papers/harabor-grastien-aaai11.pdf
// Jumper (https://github.com/Yonaba/Jumper) served as reference for this implementation.

#include <algorithm>
#include <vector>
#include <map>
#include <cmath>

#ifdef _DEBUG
#include <cassert>
#define JPS_ASSERT(cond) assert(cond)
#else
#define JPS_ASSERT(cond)
#endif

namespace JPS {

typedef std::vector<std::pair<unsigned, unsigned> > PathVector;

namespace Internal {
class Node
{
public:
	Node() {}
	Node(unsigned x, unsigned y) : parent(0), x(x), y(y), f(0), g(0), flags(0) {}
	unsigned x, y;
	unsigned f, g;
	const Node *parent;

	inline void setOpen() { flags |= 1; }
	inline void setClosed() { flags |= 2; }
	inline unsigned char isOpen() const { return flags & 1; }
	inline unsigned char isClosed() const { return flags & 2; }
	inline void clearState() { flags = 0; f = 0; g = 0; }

	inline bool operator==(const Node& o)
	{
		return x == o.x && y == o.y;
	}
private:
	unsigned char flags;
};
} // end namespace Internal

namespace Heuristic
{
	inline unsigned Manhattan(const Internal::Node *a, const Internal::Node *b)
	{
		return abs(int(a->x - b->x)) + abs(int(a->y - b->y));
	}

	inline unsigned Euclidean(const Internal::Node *a, const Internal::Node *b)
	{
		float fx = float(int(a->x - b->x));
		float fy = float(int(a->y - b->y));
		return unsigned(int(sqrtf(fx*fx + fy*fy)));
	}
} // end namespace heuristic

namespace Internal {

typedef std::vector<Node*> NodeVector;

class OpenList
{
public:
	inline void push(Node *node)
	{
		nodes.push_back(node);
		std::push_heap(nodes.begin(), nodes.end(), _compare);
	}
	inline Node *pop()
	{
		std::pop_heap(nodes.begin(), nodes.end(), _compare);
		Node *node = nodes.back();
		nodes.pop_back();
		return node;
	}
	inline bool empty() const
	{
		return nodes.empty();
	}
	inline void fixup(const Node *item)
	{
		std::make_heap(nodes.begin(), nodes.end(), _compare);

		// FIXME: rebuild heap, but only what's necessary
		/*size_t i = 0;
		for( ; i < nodes.size(); ++i)
			if(nodes[i] == item)
			{
				std::make_heap(nodes.begin() + i, nodes.end(), _compare);
				return;
			}*/
	}

protected:
	static inline bool _compare(const Node *a, const Node *b)
	{
		return a->f > b->f;
	}
	NodeVector nodes;
};

template <typename GRID> class Searcher
{
public:
	Searcher(const GRID& g)
		: grid(g)
	{}

	bool findPath(std::vector<std::pair<unsigned, unsigned> >& path, const Node& start, const Node& end, bool detail);
	Node *getNode(unsigned x, unsigned y);


private:
	const GRID& grid;
	Node *startNode;
	Node *endNode;
	OpenList open;
	NodeVector wrkmem;

	typedef std::map<std::pair<unsigned, unsigned>, Node> NodeGrid;
	NodeGrid nodegrid;

	Node * _addNodeWrk(unsigned x, unsigned y);
	void identifySuccessors(const Node *n);
	void findNeighbors(const Node *n);
	Node *jump(Node *n, const Node *parent);
};

template <typename GRID> inline Node *Searcher<GRID>::getNode(unsigned x, unsigned y)
{
	if(grid(x, y))
	{
		NodeGrid::iterator it = nodegrid.find(std::make_pair(x, y));
		if(it == nodegrid.end())
		{
			NodeGrid::iterator ins = nodegrid.insert(it, std::make_pair(std::make_pair(x, y), Node(x, y)));
			return &ins->second;
		}
		return &it->second;
	}
	return 0;
}

template <typename GRID> inline Node *Searcher<GRID>::_addNodeWrk(unsigned x, unsigned y)
{
	Node *n = getNode(x, y);
	if(n)
		wrkmem.push_back(n);
	return n;
}

template <typename GRID> Node *Searcher<GRID>::jump(Node *n, const Node *parent)
{
	if(!n)
		return 0;
	const unsigned x = n->x;
	const unsigned y = n->y;
	if(!grid(x, y))
		return 0;
	if(n == endNode)
		return n;

	int dx = x - parent->x;
	int dy = y - parent->y;

	if(dx && dy)
	{
		if( (grid(x-dx, y+dy) && !grid(x-dx, y)) || (grid(x+dx, y-dy) && !grid(x, y-dy)) )
			return n;
	}
	else if(dx)
	{
		if( (grid(x+dx, y+1) && !grid(x, y+1)) || (grid(x+dx, y-1) && !grid(x, y-1)) )
			return n;
	}
	else if(dy)
	{
		if( (grid(x+1, y+dy) && !grid(x+1, y)) || (grid(x-1, y+dy) && !grid(x-1, y)) )
			return n;
	}

	if(dx && dy)
	{
		if(jump(getNode(x+dx, y), n))
			return n;
		if(jump(getNode(x, y+dy), n))
			return n;
	}

	// TODO: get rid of this recursion
	if(grid(x+dx, y) || grid(x, y+dy))
		return jump(getNode(x+dx, y+dy), n);

	return 0;
}

template <typename GRID> void Searcher<GRID>::findNeighbors(const Node *n)
{
	const unsigned x = n->x;
	const unsigned y = n->y;
	wrkmem.clear();

#define JPS_CHECKGRID(dx, dy) (grid(x+(dx), y+(dy)))
#define JPS_ADDNODE(dx, dy) _addNodeWrk(x + (dx), y + (dy))
#define JPS_ADDNODE_NT(dx, dy) do { if(grid(x+(dx),y) || grid(x,y+(dy))) JPS_ADDNODE(dx, dy); } while(0)

	if(!n->parent)
	{
		// straight moves
		JPS_ADDNODE(-1, 0);
		JPS_ADDNODE(0, -1);
		JPS_ADDNODE(0, 1);
		JPS_ADDNODE(1, 0);

		// diagonal moves + prevent tunneling
		JPS_ADDNODE_NT(-1, -1);
		JPS_ADDNODE_NT(-1, 1);
		JPS_ADDNODE_NT(1, -1);
		JPS_ADDNODE_NT(1, 1);

		return;
	}

	// jump directions (both -1, 0, or 1)
	int dx = x - n->parent->x;
	dx /= std::max(abs(dx), 1);
	int dy = y - n->parent->y;
	dy /= std::max(abs(dy), 1);

	if(dx && dy)
	{
		// diagonal
		// natural neighbors
		bool walkX = !!JPS_ADDNODE(dx, 0);
		bool walkY = !!JPS_ADDNODE(0, dy);
		if(walkX || walkY)
			JPS_ADDNODE(dx, dy);

		// forced neighbors
		if(walkY && !JPS_CHECKGRID(-dx,0))
			JPS_ADDNODE(-dx, dy);

		if(walkX && !JPS_CHECKGRID(0,-dy))
			JPS_ADDNODE(dx, -dy);
	}
	else if(dx)
	{
		// along X axis
		JPS_ADDNODE(dx, 0);
		// Forced neighbors are up and down ahead along X
		if(!JPS_CHECKGRID(0, 1))
			JPS_ADDNODE(dx, 1);
		if(!JPS_CHECKGRID(0,-1))
			JPS_ADDNODE(dx,-1);
	}
	else if(dy)
	{
		// along Y axis
		JPS_ADDNODE(0, dy);
		// Forced neighbors are left and right ahead along Y
		if(!JPS_CHECKGRID(1, 0))
			JPS_ADDNODE(1, dy);
		if(!JPS_CHECKGRID(-1, 0))
			JPS_ADDNODE(-1,dy);
	}
#undef JPS_ADDNODE
#undef JPS_ADDNODE_NT
#undef JPS_CHECKGRID
}

template <typename GRID> void Searcher<GRID>::identifySuccessors(const Node *n)
{
	findNeighbors(n);
	for(NodeVector::reverse_iterator it = wrkmem.rbegin(); it != wrkmem.rend(); ++it)
	{
		Node *nb = *it;
		Node *jn = jump(nb, n);
		if(jn && !jn->isClosed())
		{
			unsigned extraG = Heuristic::Euclidean(jn, n);
			unsigned newG = n->g + extraG;
			if(!jn->isOpen() || newG < jn->g)
			{
				jn->g = newG;
				jn->f = jn->g + Heuristic::Manhattan(jn, endNode);
				jn->parent = n;
				if(!jn->isOpen())
				{
					open.push(jn);
					jn->setOpen();
				}
				else
					open.fixup(jn);
			}
		}
	}
	wrkmem.clear();
}

template <typename GRID> bool Searcher<GRID>::findPath(std::vector<std::pair<unsigned, unsigned> >& path, const Node& start, const Node& end, bool detail)
{
	for(NodeGrid::iterator it = nodegrid.begin(); it != nodegrid.end(); ++it)
		it->second.clearState();

	this->startNode = getNode(start.x, start.y);
	this->endNode = getNode(end.x, end.y);

	open.push(startNode);
	do
	{
		Node *n = open.pop();
		n->setClosed();
		if(n == endNode)
		{
			size_t pos = path.size();
			if(detail)
			{
				const Node *next = n;
				const Node *prev = n->parent;
				do
				{
					const unsigned x = next->x, y = next->y;
					int dx = int(prev->x - x);
					int dy = int(prev->y - y);
					JPS_ASSERT(!dx || !dy || abs(dx) == abs(dy)); // known to be straight diagonal
					const int steps = std::max(abs(dx), abs(dy));
					dx /= std::max(abs(dx), 1);
					dy /= std::max(abs(dy), 1);
					int dxa = 0, dya = 0;
					for(int i = 0; i < steps; ++i)
					{
						path.push_back(std::make_pair(x+dxa, y+dya));
						dxa += dx;
						dya += dy;
					}
					next = prev;
					prev = prev->parent;
				}
				while (prev);
				path.push_back(std::make_pair(next->x, next->y));
			}
			else
			{
				const Node *next = n;
				do
				{
					path.push_back(std::make_pair(next->x, next->y));
					next = next->parent;
				}
				while (next);
			}
			std::reverse(path.begin() + pos, path.end());
#ifdef _DEBUG
			std::cout << "Nodes generated: " << nodegrid.size() << std::endl;
#endif
			return true;
		}
		identifySuccessors(n);
	}
	while (!open.empty());
#ifdef _DEBUG
	std::cout << "Nodes generated: " << nodegrid.size() << std::endl;
#endif
	return false;
}

} // end namespace Internal

// GRID: expected to overload operator()(x, y), return true if position is walkable, false if not.
template <typename GRID> bool findPath(PathVector& path, const GRID& grid, unsigned startx, unsigned starty, unsigned endx, unsigned endy, bool detail = false)
{
	const Internal::Node startNode(startx, starty);
	const Internal::Node endNode(endx, endy);
	Internal::Searcher<GRID> search(grid);
	return search.findPath(path, startNode, endNode, detail);
}

} // end namespace JPS


#endif
