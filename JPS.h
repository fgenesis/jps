#ifndef JUMP_POINT_SEARCH_H
#define JUMP_POINT_SEARCH_H

// Public domain Jump Point Search implementation by False.Genesis
// Please keep the following source information intact when you use this file in your own projects:
// This file originates from: https://github.com/fgenesis/jps
// Based on the paper http://users.cecs.anu.edu.au/~dharabor/data/papers/harabor-grastien-aaai11.pdf
// by Daniel Harabor & Alban Grastien.
// Jumper (https://github.com/Yonaba/Jumper) and PathFinding.js (https://github.com/qiao/PathFinding.js)
// served as reference for this implementation.

// ====== COMPILE CONFIG ======

// If this is defined, compare all jumps against recursive reference implementation (only if _DEBUG is defined)
//#define JPS_VERIFY

// NYI
//#define JPS_USE_HASHMAP

// ============================

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
// ... set grid width, height, and whatever
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

#ifdef JPS_USE_HASHMAP
#include "JPSUtilHashmap.h"
#endif

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

	inline bool isValid() const { return x != unsigned(-1); }
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

namespace Internal {

static const Position npos = Pos(-1, -1);

class Node
{
public:
	Node(const Position& p) : f(0), g(0), pos(p), parent(0), flags(0) {}
	unsigned f, g;
	const Position pos;
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
	inline void clear()
	{
		nodes.clear();
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

	void freeMemory();

	bool findPath(PathVector& path, const Position& start, const Position& end, bool detail = false);

private:

#ifdef JPS_USE_HASHMAP
	typedef Util::HashMap<Position, Node> NodeGrid;
#else
	typedef std::map<Position, Node> NodeGrid;
#endif

	const GRID& grid;
	Node *endNode;
	OpenList open;

	NodeGrid nodegrid;

	Node *getNode(const Position& pos);
	void identifySuccessors(const Node *n);
	unsigned findNeighbors(const Node *n, Position *wptr) const;
	Position jumpP(const Position& p, const Position& src) const;
	Position jumpD(Position p, int dx, int dy) const;
	Position jumpX(Position p, int dx) const;
	Position jumpY(Position p, int dy) const;
	void generatePath(PathVector& path, bool detail) const;
#ifdef JPS_VERIFY
	Position jumpPRec(const Position& p, const Position& src) const;
#endif
};

template <typename GRID> inline Node *Searcher<GRID>::getNode(const Position& pos)
{
	JPS_ASSERT(grid(pos.x, pos.y));
	return &nodegrid.insert(std::make_pair(pos, Node(pos))).first->second;

	/*if(grid(x, y))
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
	return 0;*/
}

template <typename GRID> Position Searcher<GRID>::jumpP(const Position &p, const Position& src) const
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

template <typename GRID> Position Searcher<GRID>::jumpD(Position p, int dx, int dy) const
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

		if(gdx && jumpX(Pos(x+dx, y), dx).isValid())
			return p;

		if(gdy && jumpY(Pos(x, y+dy), dy).isValid())
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

template <typename GRID> inline Position Searcher<GRID>::jumpX(Position p, int dx) const
{
	JPS_ASSERT(dx);
	JPS_ASSERT(grid(p.x, p.y));

	const unsigned y = p.y;
	const Position& endpos = endNode->pos;

	unsigned a = ~((!!grid(p.x, y+1)) | ((!!grid(p.x, y-1)) << 1));

	while(true)
	{
		const unsigned xx = p.x + dx;
		const unsigned b = (!!grid(xx, y+1)) | ((!!grid(xx, y-1)) << 1);

		if((b & a) || p == endpos)
			return p;
		if(!grid(xx, y))
			return npos;

		p.x += dx;
		a = ~b;
	}
}

template <typename GRID> inline Position Searcher<GRID>::jumpY(Position p, int dy) const
{
	JPS_ASSERT(dy);
	JPS_ASSERT(grid(p.x, p.y));

	const unsigned x = p.x;
	const Position& endpos = endNode->pos;

	unsigned a = ~((!!grid(x+1, p.y)) | ((!!grid(x-1, p.y)) << 1));

	while(true)
	{
		const unsigned yy = p.y + dy;
		const unsigned b = (!!grid(x+1, yy)) | ((!!grid(x-1, yy)) << 1);

		if((a & b) || p == endpos)
			return p;
		if(!grid(x, yy))
			return npos;

		p.y += dy;
		a = ~b;
	}
}

#ifdef JPS_VERIFY
// Recursive reference implementation -- for comparison only
template <typename GRID> Position Searcher<GRID>::jumpPRec(const Position& p, const Position& src) const
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
		if(jumpPRec(Pos(x+dx, y), p).isValid())
			return p;
		if(jumpPRec(Pos(x, y+dy), p).isValid())
			return p;
	}

	if(grid(x+dx, y) || grid(x, y+dy))
		return jumpPRec(Pos(x+dx, y+dy), p);

	return npos;
}
#endif

template <typename GRID> unsigned Searcher<GRID>::findNeighbors(const Node *n, Position *wptr) const
{
	Position *w = wptr;
	const unsigned x = n->pos.x;
	const unsigned y = n->pos.y;

#define JPS_CHECKGRID(dx, dy) (grid(x+(dx), y+(dy)))
#define JPS_ADDPOS(dx, dy) 	do { if(grid(x+(dx), y+(dy))) { *w++ = Pos(x+(dx), y+(dy)); }} while(0)
#define JPS_ADDPOS_NT(dx, dy) do { if(grid(x+(dx),y) || grid(x,y+(dy))) JPS_ADDPOS(dx, dy); } while(0)

	if(!n->parent)
	{
		// straight moves
		JPS_ADDPOS(-1, 0);
		JPS_ADDPOS(0, -1);
		JPS_ADDPOS(0, 1);
		JPS_ADDPOS(1, 0);

		// diagonal moves + prevent tunneling
		JPS_ADDPOS_NT(-1, -1);
		JPS_ADDPOS_NT(-1, 1);
		JPS_ADDPOS_NT(1, -1);
		JPS_ADDPOS_NT(1, 1);

		return unsigned(w - wptr);
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
		bool walkX = false;
		bool walkY = false;
		if((walkX = grid(x+dx, y)))
			*w++ = Pos(x+dx, y);
		if((walkY = grid(x, y+dy)))
			*w++ = Pos(x, y+dy);

		if(walkX || walkY)
			JPS_ADDPOS(dx, dy);

		// forced neighbors
		if(walkY && !JPS_CHECKGRID(-dx,0))
			JPS_ADDPOS(-dx, dy);

		if(walkX && !JPS_CHECKGRID(0,-dy))
			JPS_ADDPOS(dx, -dy);
	}
	else if(dx)
	{
		// along X axis
		JPS_ADDPOS(dx, 0);
		// Forced neighbors are up and down ahead along X
		if(!JPS_CHECKGRID(0, 1))
			JPS_ADDPOS(dx, 1);
		if(!JPS_CHECKGRID(0,-1))
			JPS_ADDPOS(dx,-1);
	}
	else if(dy)
	{
		// along Y axis
		JPS_ADDPOS(0, dy);
		// Forced neighbors are left and right ahead along Y
		if(!JPS_CHECKGRID(1, 0))
			JPS_ADDPOS(1, dy);
		if(!JPS_CHECKGRID(-1, 0))
			JPS_ADDPOS(-1,dy);
	}
#undef JPS_ADDPOS
#undef JPS_ADDPOS_NT
#undef JPS_CHECKGRID

	return unsigned(w - wptr);
}

template <typename GRID> void Searcher<GRID>::identifySuccessors(const Node *n)
{
	Position buf[8];
	const int num = findNeighbors(n, &buf[0]);
	for(int i = num-1; i >= 0; --i)
	{
		// Invariant: A node is only a valid neighbor if the corresponding grid position is walkable (asserted in jumpP)
		Position jp = jumpP(buf[i], n->pos);
#ifdef JPS_VERIFY
		JPS_ASSERT(jp == jumpPRec(buf[i], n->pos));
#endif
		if(!jp.isValid())
			continue;

		// Now that the grid position is definitely a valid jump point, we have to create the actual node.
		Node *jn = getNode(jp);
		JPS_ASSERT(jn && jn != n);
		if(!jn->isClosed())
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
}

template <typename GRID> void Searcher<GRID>::generatePath(PathVector& path, bool detail) const
{
	size_t offset = path.size();
	if(detail)
	{
		const Node *next = endNode;
		const Node *prev = endNode->parent;
		do
		{
			const unsigned x = next->pos.x, y = next->pos.y;
			int dx = int(prev->pos.x - x);
			int dy = int(prev->pos.y - y);
			JPS_ASSERT(!dx || !dy || abs(dx) == abs(dy)); // known to be straight, if diagonal
			const int steps = std::max(abs(dx), abs(dy)); // leave out starting position
			dx /= std::max(abs(dx), 1);
			dy /= std::max(abs(dy), 1);
			int dxa = 0, dya = 0;
			for(int i = 0; i < steps; ++i)
			{
				path.push_back(Pos(x+dxa, y+dya));
				dxa += dx;
				dya += dy;
			}
			next = prev;
			prev = prev->parent;
		}
		while (prev);
	}
	else
	{
		const Node *next = endNode;
		do
		{
			path.push_back(next->pos);
			next = next->parent;
		}
		while (next->parent);
	}
	std::reverse(path.begin() + offset, path.end());
}

template <typename GRID> bool Searcher<GRID>::findPath(PathVector& path, const Position& start, const Position& end, bool detail /* = false */)
{
	for(NodeGrid::iterator it = nodegrid.begin(); it != nodegrid.end(); ++it)
		it->second.clearState();

	if(start == end)
	{
		// There is only a path if this single position is walkable.
		// But since the starting position is omitted, there is nothing to do here.
		return grid(end.x, end.y);
	}

	// If start or end point are obstructed, don't even start
	if(!grid(start.x, start.y) || !grid(end.x, end.y))
		return false;

	open.push(getNode(start));
	endNode = getNode(end);
	JPS_ASSERT(endNode);

	do
	{
		Node *n = open.pop();
		n->setClosed();
		if(n == endNode)
		{
			open.clear();
			generatePath(path, detail);
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
//       The path does NOT contain the starting position, i.e. if start and end are the same, the resulting path has no elements.
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
