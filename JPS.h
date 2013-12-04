#ifndef JUMP_POINT_SEARCH_H
#define JUMP_POINT_SEARCH_H

// Public domain Jump Point Search implementation by False.Genesis
// Please keep the following source information intact when you use this file in your own projects:
// This file originates from: https://github.com/fgenesis/jps
// Based on the paper http://users.cecs.anu.edu.au/~dharabor/data/papers/harabor-grastien-aaai11.pdf
// by Daniel Harabor & Alban Grastien.
// Jumper (https://github.com/Yonaba/Jumper) and PathFinding.js (https://github.com/qiao/PathFinding.js)
// served as reference for this implementation.

// Usage:
/*
// Define a class that overloads `operator()(x, y) const`, returning a value that can be treated as boolean.
// You are responsible for bounds checking!
// You want your operator() to be as fast as possible, as it will be called a LOT.

struct MyGrid
{
	inline bool operator()(unsigned x, unsigned y) const
	{
		if(x < width && y < height) // Unsigned will wrap if < 0
			... return true if terrain at (x, y) is walkable.
	}
	unsigned width, height;
};

// Then you can retrieve a path:

MyGrid grid;
// ... set grid with, height, and whatever
bool detailed = false; // set this to true if you want a detailed single-step path
                       // (e.g. if you plan to further mangle the path yourself)
JPS::PathVector path; // The resulting path will go here.


// Single-call interface:
bool found = JPS::findPath(path, grid, startx, starty, endx, endy, detailed);


// Alternatively, if you want more control:

JPS::Searcher<MyGrid> search(grid);
while(true)
{
	// ..stuff happening ...

	// build path incrementally from waypoints
	JPS::Position a, b, c, d; // some waypoints
	search.findPath(path, a, b);
	search.findPath(path, b, c);
	search.findPath(path, c, d);

	if(search.findPath(path2, JPS::Pos(startx, starty), JPS::Pos(endx, endy), detailed))
	{
		// ...handle failure...
	}
	// ... more stuff happening ...

	// At convenient times, you can clean up accumulated nodes to reclaim memory.
	// This is never necessary, but performance will drop if too many cached nodes exist.
	if(mapWasReloaded)
		search.freeMemory();
}

// Further remarks can be found at the bottom of this file.
*/


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

// If this is defined, compare all jumps against recursive reference implementation (only if _DEBUG is defined)
//#define JPS_VERIFY

namespace JPS {

struct Position
{
	unsigned x, y;

	inline bool operator==(const Position& p) const
	{
		return x == p.x && y == p.y;
	}
	inline bool operator!=(const Position& p) const
	{
		return x != p.x || y != p.y;
	}

	// for sorting
	inline bool operator<(const Position& p) const
	{
		return y < p.y || (y == p.y && x < p.x);
	}
};

typedef std::vector<Position> PathVector;

// ctor function to keep Position a real POD struct.
inline static Position Pos(unsigned x, unsigned y)
{
	Position p;
	p.x = x;
	p.y = y;
	return p;
}

static const Position npos = Pos(-1, -1);

namespace Internal {
class Node
{
public:
	Node(unsigned x, unsigned y) : pos(Pos(x, y)), f(0), g(0), parent(0), flags(0) {}
	const Position pos;
	unsigned f, g;
	const Node *parent;

	inline void setOpen() { flags |= 1; }
	inline void setClosed() { flags |= 2; }
	inline unsigned char isOpen() const { return flags & 1; }
	inline unsigned char isClosed() const { return flags & 2; }
	inline void clearState() { f = 0; g = 0, parent = 0; flags = 0; }

private:
	unsigned char flags;

	bool operator==(const Node& o); // not implemented, nodes should not be compared
};
} // end namespace Internal

namespace Heuristic
{
	inline unsigned Manhattan(const Internal::Node *a, const Internal::Node *b)
	{
		return abs(int(a->pos.x - b->pos.x)) + abs(int(a->pos.y - b->pos.y));
	}

	inline unsigned Euclidean(const Internal::Node *a, const Internal::Node *b)
	{
		float fx = float(int(a->pos.x - b->pos.x));
		float fy = float(int(a->pos.y - b->pos.y));
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
		JPS_ASSERT(node);
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
	{
		wrkmem.reserve(8);
	}

	void freeMemory();

	bool findPath(PathVector& path, const Position& start, const Position& end, bool detail = false);

private:

	typedef std::map<Position, Node> NodeGrid;

	const GRID& grid;
	Node *endNode;
	OpenList open;
	NodeVector wrkmem;

	NodeGrid nodegrid;

	Node *getNode(unsigned x, unsigned y);
	Node * _addNodeWrk(unsigned x, unsigned y);
	void identifySuccessors(const Node *n);
	void findNeighbors(const Node *n);
	Node *jump(Node *n, const Node *parent);
	Position jumpP(Position p, const Position& src);
	Position jumpD(Position p, int dx, int dy);
	Position jumpX(Position p, int dx);
	Position jumpY(Position p, int dy);
#ifdef JPS_VERIFY
	Position jumpPRec(const Position& p, const Position& src);
#endif
};

template <typename GRID> inline Node *Searcher<GRID>::getNode(unsigned x, unsigned y)
{
	if(grid(x, y))
	{
		const Position p = Pos(x, y);
		NodeGrid::iterator it = nodegrid.find(p);
		if(it == nodegrid.end())
		{
			NodeGrid::iterator ins = nodegrid.insert(it, std::make_pair(p, Node(x, y)));
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
	Position p = jumpP(n->pos, parent->pos);
#ifdef JPS_VERIFY
	JPS_ASSERT(p == jumpPRec(n->pos, parent->pos));
#endif
	return p != npos ? getNode(p.x, p.y) : 0;
}

template <typename GRID> Position Searcher<GRID>::jumpP(Position p, const Position& src)
{
	JPS_ASSERT(grid(p.x, p.y));

	int dx = int(p.x - src.x);
	int dy = int(p.y - src.y);
	JPS_ASSERT(dx || dy);

	if(dx && dy)
		return jumpD(p, dx, dy);
	else if(dx)
		return jumpX(p, dx);
	else if(dy)
		return jumpY(p, dy);

	// not reached
	JPS_ASSERT(false);
	return npos;
}

template <typename GRID> Position Searcher<GRID>::jumpD(Position p, int dx, int dy)
{
	JPS_ASSERT(grid(p.x, p.y));
	JPS_ASSERT(dx && dy);

	const Position& endpos = endNode->pos;

	while(true)
	{
		if(p == endpos)
			return p;

		const unsigned x = p.x;
		const unsigned y = p.y;

		if( (grid(x-dx, y+dy) && !grid(x-dx, y)) || (grid(x+dx, y-dy) && !grid(x, y-dy)) )
			return p;

		const bool gdx = grid(x+dx, y);
		const bool gdy = grid(x, y+dy);

		if(gdx && jumpX(Pos(x+dx, y), dx) != npos)
			return p;

		if(gdy && jumpY(Pos(x, y+dy), dy) != npos)
			return p;

		if((gdx || gdy) && grid(x+dx, y+dy))
		{
			p.x += dx;
			p.y += dy;
		}
		else
			break;
	}

	return npos;
}

template <typename GRID> Position Searcher<GRID>::jumpX(Position p, int dx)
{
	JPS_ASSERT(dx);
	JPS_ASSERT(grid(p.x, p.y));

	const unsigned y = p.y;
	const Position& endpos = endNode->pos;
	while(true)
	{
		if(p == endpos)
			return p;

		const unsigned x = p.x;

		if((grid(x+dx, y+1) && !grid(x, y+1) ) || (grid(x+dx, y-1)  && !grid(x, y-1)))
			return p;

		if(grid(x+dx, y))
			p.x += dx;
		else
			break;
	}

	return npos;
}

template <typename GRID> Position Searcher<GRID>::jumpY(Position p, int dy)
{
	JPS_ASSERT(dy);
	JPS_ASSERT(grid(p.x, p.y));

	const unsigned x = p.x;
	const Position& endpos = endNode->pos;
	while(true)
	{
		if(p == endpos)
			return p;

		const unsigned y = p.y;

		if(((grid(x+1, y+dy) && !grid(x+1, y)) || (grid(x-1, y+dy)  && !grid(x-1, y))))
			return p;

		if(grid(x, y+dy))
			p.y += dy;
		else
			break;
	}

	return npos;
}

#ifdef JPS_VERIFY
// Recursive reference implementation -- for comparison only
template <typename GRID> Position Searcher<GRID>::jumpPRec(const Position& p, const Position& src)
{
	unsigned x = p.x;
	unsigned y = p.y;
	if(!grid(x, y))
		return npos;
	if(p == endNode->pos)
		return p;

	int dx = int(x - src.x);
	int dy = int(y - src.y);
	JPS_ASSERT(dx || dy);

	if(dx && dy)
	{
		if( (grid(x-dx, y+dy) && !grid(x-dx, y)) || (grid(x+dx, y-dy) && !grid(x, y-dy)) )
			return p;
	}
	else if(dx)
	{
		if( (grid(x+dx, y+1) && !grid(x, y+1)) || (grid(x+dx, y-1) && !grid(x, y-1)) )
			return p;
	}
	else if(dy)
	{
		if( (grid(x+1, y+dy) && !grid(x+1, y)) || (grid(x-1, y+dy) && !grid(x-1, y)) )
			return p;
	}

	if(dx && dy)
	{
		if(jumpPRec(Pos(x+dx, y), p) != npos)
			return p;
		if(jumpPRec(Pos(x, y+dy), p) != npos)
			return p;
	}

	if(grid(x+dx, y) || grid(x, y+dy))
		return jumpPRec(Pos(x+dx, y+dy), p);

	return npos;
}
#endif

template <typename GRID> void Searcher<GRID>::findNeighbors(const Node *n)
{
	const unsigned x = n->pos.x;
	const unsigned y = n->pos.y;
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
	int dx = int(x - n->parent->pos.x);
	dx /= std::max(abs(dx), 1);
	int dy = int(y - n->parent->pos.y);
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
		JPS_ASSERT(jn != n);
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

template <typename GRID> bool Searcher<GRID>::findPath(PathVector& path, const Position& start, const Position& end, bool detail /* = false */)
{
	for(NodeGrid::iterator it = nodegrid.begin(); it != nodegrid.end(); ++it)
		it->second.clearState();

	if(start == end)
	{
		if(grid(end.x, end.y)) // There is only a path if this single position is walkable.
		{
			path.push_back(end);
			return true;
		}
		return false;
	}

	// If start or end point are obstructed, don't even start
	{
		Node *startNode = getNode(start.x, start.y);
		if(!startNode)
			return false;
		if(!(endNode = getNode(end.x, end.y)))
			return false;
		open.push(startNode);
	}

	do
	{
		Node *n = open.pop();
		n->setClosed();
		if(n == endNode)
		{
			size_t offset = path.size();
			if(detail)
			{
				const Node *next = n;
				const Node *prev = n->parent;
				bool first = true;
				do
				{
					const unsigned x = next->pos.x, y = next->pos.y;
					int dx = int(prev->pos.x - x);
					int dy = int(prev->pos.y - y);
					JPS_ASSERT(!dx || !dy || abs(dx) == abs(dy)); // known to be straight, if diagonal
					const int steps = std::max(abs(dx), abs(dy));
					dx /= std::max(abs(dx), 1);
					dy /= std::max(abs(dy), 1);
					int dxa = 0, dya = 0;
					for(int i = 0; i < steps; ++i)
					{
						if(first)
							first = false;
						else
							path.push_back(Pos(x+dxa, y+dya));
						dxa += dx;
						dya += dy;
					}
					next = prev;
					prev = prev->parent;
				}
				while (prev);
				path.push_back(next->pos);
			}
			else
			{
				const Node *next = n;
				do
				{
					path.push_back(next->pos);
					next = next->parent;
				}
				while (next->parent);
			}
			std::reverse(path.begin() + offset, path.end());
			return true;
		}
		identifySuccessors(n);
	}
	while (!open.empty());
	return false;
}

template<typename GRID> void Searcher<GRID>::freeMemory()
{
	NodeGrid v;
	nodegrid.swap(v);
	// other containers known to be empty.
}

} // end namespace Internal

using Internal::Searcher;

// path: If the function returns true, the path is stored in this vector.
//       The path does NOT contain the starting position.
//       The vector does not have to be empty. The function does not clear it;
//       instead, the new path positions are appended at the end.
//       This allows building a path incrementally.
// grid: expected to overload operator()(x, y), return true if position is walkable, false if not.
// detail: If true, create exhaustive step-by-step path.
//         If false, only return waypoints. Waypoints are guaranteed to be on a straight line (vertically, horizontally, or diagonally),
//         and there is no obstruction between waypoints.
template <typename GRID> bool findPath(PathVector& path, const GRID& grid, unsigned startx, unsigned starty, unsigned endx, unsigned endy, bool detail = false)
{
	Searcher<GRID> search(grid);
	return search.findPath(path, Pos(startx, starty), Pos(endx, endy), detail);
}



} // end namespace JPS


#endif
