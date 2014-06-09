#ifndef JUMP_POINT_SEARCH_H
#define JUMP_POINT_SEARCH_H

// Public domain Jump Point Search implementation by False.Genesis
// Very fast pathfinding for uniform cost grids.
// Supports incremental pathfinding.

// Please keep the following source information intact when you use this file in your own projects:
// This file originates from: https://github.com/fgenesis/jps
// Based on the paper http://users.cecs.anu.edu.au/~dharabor/data/papers/harabor-grastien-aaai11.pdf
// by Daniel Harabor & Alban Grastien.
// Jumper (https://github.com/Yonaba/Jumper) and PathFinding.js (https://github.com/qiao/PathFinding.js)
// served as reference for this implementation.
// If you use this, attribution would be nice, but is not necessary.

// ====== COMPILE CONFIG ======

// If this is defined, compare all jumps against recursive reference implementation (only if _DEBUG is defined)
#define JPS_VERIFY

// If this is defined, search nodes will cache jump points for each possible direction.
// Test your use cases whether jump point caching gives a speedup for you or not.
// If you only use the single-call interface, caching will not be of any use so leave this commented out.
// Search nodes will use 4 + 8*5*sizeof(void*) bytes more memory with caching enabled.
// IMPORTANT: If caching is on and the grid changes, you need to call Searcher::flushCache() before starting another search!!
#define JPS_CACHE

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
unsigned step = 0; // set this to 1 if you want a detailed single-step path
                   // (e.g. if you plan to further mangle the path yourself),
                   // or any other higher value to output every Nth position.
JPS::PathVector path; // The resulting path will go here.


// Single-call interface:
bool found = JPS::findPath(path, grid, startx, starty, endx, endy, step);


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

	if(!search.findPath(path2, JPS::Pos(startx, starty), JPS::Pos(endx, endy), step))
	{
		// ...handle failure...
	}
	// ... more stuff happening ...

	// At convenient times, you can clean up accumulated nodes to reclaim memory.
	// This is never necessary, but performance will drop if too many cached nodes exist.
	if(mapWasReloaded)
		search.flushCache();
}

// Further remarks about the super lazy single-call function can be found at the bottom of this file.

// -------------------------------
// --- Incremental pathfinding ---
// -------------------------------

First, call findPathInit(Position start, Position end). Don't forget to check the return value,
as it may return NO_PATH if one or both of the points are obstructed,
or FOUND_PATH if the points are equal and not obstructed.
If it returns NEED_MORE_STEPS then the next part can start.

Repeatedly call findPathStep(int limit) until it returns NO_PATH or FOUND_PATH.
For consistency, you will want to ensure that the grid does not change between subsequent calls;
if the grid changes, parts of the path may go through a now obstructed area.
If limit is 0, it will perform the pathfinding in one go. Values > 0 abort the search
as soon as possible after the number of steps was exceeded, returning NEED_MORE_STEPS.
Use getStepsDone() after some test runs to find a good value for the limit.

Lastly, generate the actual path points from a successful run via findPathFinish(PathVector& path, unsigned step = 0).
Like described above, path points are appended, and granularity can be adjusted with the step parameter.
Returns false if the pathfinding did not finish or generating the path failed.

*/


#include <stdlib.h>
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

enum Result
{
	NO_PATH = 0,
	FOUND_PATH = 1,
	NEED_MORE_STEPS = 2
};

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

struct JumpResult
{
	JumpResult(const Position& p, bool hitEnd) : p(p), hitEnd(hitEnd) {}
	inline bool operator==(const JumpResult& o) const { return p == o.p && hitEnd == o.hitEnd; }
	inline bool operator!=(const JumpResult& o) const { return !(*this == o); }
	Position p;
	bool hitEnd;
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
static const JumpResult njmp(npos, false);

enum
{
	MAX_DIRECTIONS = 8,
	MAX_NEIGHBORS = 5 // max. possible directions if the search is moving
};

class Node
{
public:
	Node(const Position& p) : f(0), g(0), pos(p), parent(0)
	{
		for(unsigned i = 0; i < (sizeof(flags) / sizeof(flags[0])); ++i)
			flags[i] = 0;
		// jumpPoints array can stay uninitialized
	}

	unsigned f, g;
	const Position pos;
	const Node *parent;

	inline void     setOpen()        {        flags[0] |= (1 << 31); }
	inline unsigned isOpen() const   { return flags[0] &  (1 << 31); }
	inline void     setClosed()      {        flags[0] |= (1 << 30); }
	inline unsigned isClosed() const { return flags[0] &  (1 << 30); }

	inline void clearState()
	{
		f = 0;
		g = 0;
		parent = NULL;
		// clear open and closed bits(known jump points stay intact)
		flags[0] &= ~(3 << 30);
	}


private:
	bool operator==(const Node& o); // not implemented, nodes should not be compared

	// 0 = unused, O = open?, C = closed?
	// a-h: cache bits for all directions
	// [bits  0-31]: OC00 dddc ccbb baaa
	// [bits 32-63]: 0000 hhhg ggff feee
#ifndef JPS_CACHE
	unsigned flags[1];
#else
	unsigned flags[2];
public:

	inline void setNumJumpPoints(int dx, int dy, unsigned n)
	{
		JPS_ASSERT(n <= MAX_NEIGHBORS);
		const unsigned both = dx & dy & 1;
		const unsigned tristart = bitpos(dx, dy, both);
		const unsigned mask = 7 << tristart;
		unsigned f = flags[both];
		f &= ~mask; // mask out old number of jump points if present
		f |= (n+1) << tristart; // set num
		flags[both] = f;
	}

	inline Node **getJumpPointsForDirection(int dx, int dy, int *num)
	{
		const unsigned both = dx & dy & 1;
		const unsigned tri = bitidx(dx, dy, both);
		JPS_ASSERT(tri <= 3);
		*num = ((flags[both] >> (3*tri)) & 7) - 1; // will result in -1 if number of jump points is unknown
		return &jumpPoints[(both << 2) | tri][0];
	}

private:

	// use 2 groups of tuples (dx, dy):
	// both set: (-1, -1), (-1, 1), (1, -1), (1, 1)
	// one set:  (0, -1), (0, 1), (-1, 0), (1, 0)
	// This function will assign each 4 possibilities in each group a unique index [0..3].
	// Note: the bit magic below was determined experimentally. There might be a faster way to do this.
	inline static unsigned bitidx(unsigned ux, unsigned uy, unsigned both)
	{
		const unsigned dxs = ux >> 31;
		const unsigned a = ((uy&1) ^ (uy>>31)) << 1;
		const unsigned b = (ux&1) + dxs;
		const unsigned c = (~both) & dxs;
		return (a ^ b) + c; // this is now in [0..3]
	}

	inline static unsigned bitpos(unsigned ux, unsigned uy, unsigned both)
	{
		return 3 * bitidx(ux, uy, both); // We have 3-bit blocks
	}

	/*inline static unsigned bitdir(unsigned ux, unsigned uy, unsigned both)
	{
		return (both << 2) | bitidx(ux, uy, ux & uy & 1); // [0..7]
	}*/

	Node *jumpPoints[MAX_DIRECTIONS][MAX_NEIGHBORS]; // cached jump points for each direction

#endif
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
typedef std::map<Position, Node> NodeGrid;

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
	}

protected:
	static inline bool _compare(const Node *a, const Node *b)
	{
		return a->f > b->f;
	}
	NodeVector nodes;
};

// Base class for thin template to avoid template bloat
class SearcherBase
{
public:
	SearcherBase() : endNode(NULL), skip(1), stepsRemain(0), stepsDone(0) {}

	void flushCache();

	bool generatePath(PathVector& path, unsigned step) const;
	Node *getNode(const Position& pos);
	void followSuccessors(Node *n, unsigned njump, Node **jumpPoints);

	Node *endNode;
	int skip;
	int stepsRemain;
	size_t stepsDone;
	OpenList open;

	NodeGrid nodegrid;
};

template <typename GRID> class Searcher
{
public:
	Searcher(const GRID& g)
		: grid(g)
	{}

	// single-call
	bool findPath(PathVector& path, Position start, Position end, unsigned step = 0);

	// incremental pathfinding
	Result findPathInit(Position start, Position end);
	Result findPathStep(int limit);
	bool findPathFinish(PathVector& path, unsigned step = 0);

	// misc
	inline void flushCache() { base.flushCache(); }
	inline void setSkip(int s) { base.skip = std::max(1, s); }
	inline size_t getStepsDone() const { return base.stepsDone; }
	inline size_t getNodesExpanded() const { return base.nodegrid.size(); }

private:

	const GRID& grid; // intentionally very first member
	SearcherBase base;

	void identifySuccessors(Node *n);
	bool identifyJumpPoints(Node *n, Node **jumpPoints, int *njump);
	unsigned findNeighbors(const Node *n, Position *wptr) const;
	unsigned findNeighborsNoParent(const Node *n, Position *wptr) const;
	JumpResult jumpP(const Position& p, const Position& src);
	JumpResult jumpD(Position p, int dx, int dy);
	JumpResult jumpX(Position p, int dx);
	JumpResult jumpY(Position p, int dy);

#ifdef JPS_VERIFY
	JumpResult jumpPRec(const Position& p, const Position& src) const;
#endif
};

// ------------- Start of thin base methods ---------------------

void SearcherBase::flushCache()
{
	nodegrid.clear();
	endNode = NULL;
	open.clear();
	// other containers known to be empty.
}

inline Node *SearcherBase::getNode(const Position& pos)
{
	return &nodegrid.insert(std::make_pair(pos, Node(pos))).first->second;
}

bool SearcherBase::generatePath(PathVector& path, unsigned step) const
{
	if(!endNode)
		return false;
	size_t offset = path.size();
	if(step)
	{
		const Node *next = endNode;
		const Node *prev = endNode->parent;
		if(!prev)
			return false;
		do
		{
			const unsigned x = next->pos.x, y = next->pos.y;
			int dx = int(prev->pos.x - x);
			int dy = int(prev->pos.y - y);
			JPS_ASSERT(!dx || !dy || abs(dx) == abs(dy)); // known to be straight, if diagonal
			const int steps = std::max(abs(dx), abs(dy));
			dx /= std::max(abs(dx), 1);
			dy /= std::max(abs(dy), 1);
			dx *= int(step);
			dy *= int(step);
			int dxa = 0, dya = 0;
			for(int i = 0; i < steps; i += step)
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
		if(!next->parent)
			return false;
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

void SearcherBase::followSuccessors(Node *n, unsigned njump, Node **jumpPoints)
{
	for(unsigned i = 0; i < njump; ++i)
	{
		Node *jn = jumpPoints[i];
		if(!jn->isClosed())
		{
			const unsigned extraG = Heuristic::Euclidean(jn, n);
			const unsigned newG = n->g + extraG;
			if(!jn->isOpen() || newG < jn->g)
			{
				jn->g = newG;
				jn->f = newG + Heuristic::Manhattan(jn, endNode);
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

// ---------------- Start of template code ----------------

template <typename GRID> JumpResult Searcher<GRID>::jumpP(const Position &p, const Position& src)
{
	JPS_ASSERT(grid(p.x, p.y));

	int dx = int(p.x - src.x);
	int dy = int(p.y - src.y);
	JPS_ASSERT(dx || dy);
	printf("jumpP(%u, %u), dx = %d, dy = %d\n", p.x, p.y, dx, dy);

	if(dx && dy)
		return jumpD(p, dx, dy);
	else if(dx)
		return jumpX(p, dx);
	else if(dy)
		return jumpY(p, dy);

	// not reached
	JPS_ASSERT(false);
	return njmp;
}

template <typename GRID> JumpResult Searcher<GRID>::jumpD(Position p, int dx, int dy)
{
	printf("jumpD(%u, %u), dx = %d, dy = %d\n", p.x, p.y, dx, dy);
	JPS_ASSERT(grid(p.x, p.y));
	JPS_ASSERT(dx && dy);

	const Position endpos = base.endNode->pos;
	int steps = 0;
	bool hitEnd = false;

	while(true)
	{
		if(p == endpos)
		{
			hitEnd = true;
			break;
		}

		++steps;
		const unsigned x = p.x;
		const unsigned y = p.y;

		if( (grid(x-dx, y+dy) && !grid(x-dx, y)) || (grid(x+dx, y-dy) && !grid(x, y-dy)) )
			break;

		const bool gdx = grid(x+dx, y);
		const bool gdy = grid(x, y+dy);

		if(gdx)
		{
			JumpResult jr = jumpX(Pos(x+dx, y), dx);
			if(jr.p.isValid())
			{
				hitEnd = jr.hitEnd;
				break;
			}
		}

		if(gdy)
		{
			JumpResult jr = jumpY(Pos(x, y+dy), dy);
			if(jr.p.isValid())
			{
				hitEnd = jr.hitEnd;
				break;
			}
		}

		if((gdx || gdy) && grid(x+dx, y+dy))
		{
			p.x += dx;
			p.y += dy;
		}
		else
		{
			p = npos;
			break;
		}
	}
	base.stepsDone += (unsigned)steps;
	base.stepsRemain -= steps;
	printf(" <- (%u, %u)\n", p.x, p.y);
	return JumpResult(p, hitEnd);
}

template <typename GRID> inline JumpResult Searcher<GRID>::jumpX(Position p, int dx)
{
	printf("jumpX(%u, %u), dx = %d\n", p.x, p.y, dx);
	bool v = p.x == 344 && p.y == 376 && dx == 1;
	JPS_ASSERT(dx);
	JPS_ASSERT(grid(p.x, p.y));

	const unsigned y = p.y;
	const Position endpos = base.endNode->pos;
	const unsigned skip = base.skip;
	int steps = 0;
	bool hitEnd = false;

	unsigned a = ~((!!grid(p.x, y+skip)) | ((!!grid(p.x, y-skip)) << 1));

	while(true)
	{
		if(p == endpos)
		{
			hitEnd = true;
			break;
		}
		const unsigned xx = p.x + dx;
		const unsigned b = (!!grid(xx, y+skip)) | ((!!grid(xx, y-skip)) << 1);
		if(v)
			printf("# xx = %u, a = %u, b = %u, endpos = (%u, %u), steps = %u\n", xx, a, b, endpos.x, endpos.y, steps);

		if(b & a)
			break;
		if(!grid(xx, y))
		{
			p = npos;
			break;
		}

		p.x = xx;
		a = ~b;
		++steps;
	}

	base.stepsDone += (unsigned)steps;
	base.stepsRemain -= steps;
	printf(" <- (%u, %u)\n", p.x, p.y);
	return JumpResult(p, hitEnd);
}

template <typename GRID> inline JumpResult Searcher<GRID>::jumpY(Position p, int dy)
{
	printf("jumpY(%u, %u), dy = %d\n", p.x, p.y, dy);
	JPS_ASSERT(dy);
	JPS_ASSERT(grid(p.x, p.y));

	const unsigned x = p.x;
	const Position endpos = base.endNode->pos;
	const unsigned skip = base.skip;
	int steps = 0;
	bool hitEnd = false;

	unsigned a = ~((!!grid(x+skip, p.y)) | ((!!grid(x-skip, p.y)) << 1));

	while(true)
	{
		if(p == endpos)
		{
			hitEnd = true;
			break;
		}

		const unsigned yy = p.y + dy;
		const unsigned b = (!!grid(x+skip, yy)) | ((!!grid(x-skip, yy)) << 1);

		if(a & b)
			break;
		if(!grid(x, yy))
		{
			p = npos;
			break;
		}

		p.y = yy;
		a = ~b;
	}

	base.stepsDone += (unsigned)steps;
	base.stepsRemain -= steps;
	printf(" <- (%u, %u)\n", p.x, p.y);
	return JumpResult(p, hitEnd);
}

#ifdef JPS_VERIFY
// Recursive reference implementation -- for comparison only
template <typename GRID> JumpResult Searcher<GRID>::jumpPRec(const Position& p, const Position& src) const
{
	const unsigned x = p.x;
	const unsigned y = p.y;
	const int skip = base.skip;
	if(!grid(x, y))
		return JumpResult(npos, false);
	if(p == base.endNode->pos)
		return JumpResult(p, true);

	int dx = int(x - src.x);
	int dy = int(y - src.y);
	JPS_ASSERT(dx || dy);

	if(dx && dy)
	{
		if( (grid(x-dx, y+dy) && !grid(x-dx, y)) || (grid(x+dx, y-dy) && !grid(x, y-dy)) )
			return JumpResult(p, false);
	}
	else if(dx)
	{
		if( (grid(x+dx, y+skip) && !grid(x, y+skip)) || (grid(x+dx, y-skip) && !grid(x, y-skip)) )
			return JumpResult(p, false);
	}
	else if(dy)
	{
		if( (grid(x+skip, y+dy) && !grid(x+skip, y)) || (grid(x-skip, y+dy) && !grid(x-skip, y)) )
			return JumpResult(p, false);
	}

	if(dx && dy)
	{
		JumpResult jr = jumpPRec(Pos(x+dx, y), p);
		if(jr.p.isValid())
			return JumpResult(p, jr.hitEnd);
		jr = jumpPRec(Pos(x, y+dy), p);
		if(jr.p.isValid())
			return JumpResult(p, jr.hitEnd);
	}

	if(grid(x+dx, y) || grid(x, y+dy))
		return jumpPRec(Pos(x+dx, y+dy), p);

	return njmp;
}
#endif


#define JPS_CHECKGRID(dx, dy) (grid(x+(dx), y+(dy)))
#define JPS_ADDPOS(dx, dy) 	do { *w++ = Pos(x+(dx), y+(dy)); } while(0)
#define JPS_ADDPOS_CHECK(dx, dy) do { if(JPS_CHECKGRID(dx, dy)) JPS_ADDPOS(dx, dy); } while(0)
#define JPS_ADDPOS_NO_TUNNEL(dx, dy) do { if(grid(x+(dx),y) || grid(x,y+(dy))) JPS_ADDPOS_CHECK(dx, dy); } while(0)

template <typename GRID> unsigned Searcher<GRID>::findNeighborsNoParent(const Node *n, Position *wptr) const
{
	JPS_ASSERT(!n->parent);

	Position *w = wptr;
	const int skip = base.skip;
	const unsigned x = n->pos.x;
	const unsigned y = n->pos.y;

	// straight moves
	JPS_ADDPOS_CHECK(-skip, 0);
	JPS_ADDPOS_CHECK(0, -skip);
	JPS_ADDPOS_CHECK(0, skip);
	JPS_ADDPOS_CHECK(skip, 0);

	// diagonal moves + prevent tunneling
	JPS_ADDPOS_NO_TUNNEL(-skip, -skip);
	JPS_ADDPOS_NO_TUNNEL(-skip, skip);
	JPS_ADDPOS_NO_TUNNEL(skip, -skip);
	JPS_ADDPOS_NO_TUNNEL(skip, skip);

	return unsigned(w - wptr);
}

template <typename GRID> unsigned Searcher<GRID>::findNeighbors(const Node *n, Position *wptr) const
{
	JPS_ASSERT(n->parent);

	Position *w = wptr;
	const int skip = base.skip;
	const unsigned x = n->pos.x;
	const unsigned y = n->pos.y;

	// jump directions (both -1, 0, or 1)
	int dx = int(x - n->parent->pos.x);
	int dy = int(y - n->parent->pos.y);

	if(dx && dy)
	{
		dx = (dx / abs(dx)) * skip;
		dy = (dy / abs(dy)) * skip;

		// diagonal
		// natural neighbors
		bool walkX = false;
		bool walkY = false;
		if((walkX = grid(x+dx, y)))
			*w++ = Pos(x+dx, y);
		if((walkY = grid(x, y+dy)))
			*w++ = Pos(x, y+dy);

		if(walkX || walkY)
			JPS_ADDPOS_CHECK(dx, dy);

		// forced neighbors
		if(walkY && !JPS_CHECKGRID(-dx,0))
			JPS_ADDPOS_CHECK(-dx, dy);

		if(walkX && !JPS_CHECKGRID(0,-dy))
			JPS_ADDPOS_CHECK(dx, -dy);
	}
	else if(dx)
	{
		dx = (dx / abs(dx)) * skip;

		// along X axis
		if(JPS_CHECKGRID(dx, 0))
		{
			JPS_ADDPOS(dx, 0);

			// Forced neighbors (+ prevent tunneling)
			if(!JPS_CHECKGRID(0, skip))
				JPS_ADDPOS_CHECK(dx, skip);
			if(!JPS_CHECKGRID(0,-skip))
				JPS_ADDPOS_CHECK(dx,-skip);
		}
	}
	else if(dy)
	{
		dy = (dy / abs(dy)) * skip;

		// along Y axis
		if(JPS_CHECKGRID(0, dy))
		{
			JPS_ADDPOS(0, dy);

			// Forced neighbors (+ prevent tunneling)
			if(!JPS_CHECKGRID(skip, 0))
				JPS_ADDPOS_CHECK(skip, dy);
			if(!JPS_CHECKGRID(-skip, 0))
				JPS_ADDPOS_CHECK(-skip,dy);
		}
	}

	return unsigned(w - wptr);
}

#undef JPS_ADDPOS
#undef JPS_ADDPOS_CHECK
#undef JPS_ADDPOS_NO_TUNNEL
#undef JPS_CHECKGRID

template <typename GRID> bool Searcher<GRID>::identifyJumpPoints(Node *n, Node **jumpPoints, int *njumpP)
{
	int njump = 0;
	Position buf[MAX_DIRECTIONS];
	const int num = n->parent ? findNeighbors(n, &buf[0]) : findNeighborsNoParent(n, &buf[0]);
	bool hitEnd = false;
	for(int i = num-1; i >= 0; --i)
	{
		printf("Jump from (%u, %u) to (%u, %u)\n", n->pos.x, n->pos.y, buf[i].x, buf[i].y);

		// Invariant: A node is only a valid neighbor if the corresponding grid position is walkable (asserted in jumpP)
		JumpResult jr = jumpP(buf[i], n->pos);

#ifdef JPS_VERIFY
		{
			JumpResult jrcheck = jumpPRec(buf[i], n->pos);
			JPS_ASSERT(jr == jrcheck);
		}
#endif

		if(!jr.p.isValid())
			continue;

		hitEnd = hitEnd || jr.hitEnd;

		// Now that the grid position is definitely a valid jump point, we have to create the actual node.
		Node *jn = base.getNode(jr.p);
		JPS_ASSERT(jn && jn != n);
		printf("# Hit end: %d\n", jr.hitEnd);

		jumpPoints[njump++] = jn;
	}
	*njumpP = njump;
	printf("#### Hit end: %d\n", hitEnd);
	return hitEnd; // indicate not to cache jump point if traversal was aborted early because end was hit
}

template <typename GRID> void Searcher<GRID>::identifySuccessors(Node *n)
{
	Node *jumpPoints[MAX_DIRECTIONS];
	Node **jumpPointsUse = &jumpPoints[0];
	int njump;

#ifndef JPS_CACHE
	identifyJumpPoints(n, jumpPoints, &njump);
#else
	// Jump points for this node stay constant as long as the grid stays constant.
	// If the jump points are not known yet, collect them.
	// If they are known, the slow map lookup can be skipped.

	// No parent -- this is the initial node, dx == dy == 0.
	// Do not cache jump points in this case
	if(!n->parent)
		identifyJumpPoints(n, jumpPoints, &njump);
	else
	{
		int dx = int(n->pos.x - n->parent->pos.x);
		int dy = int(n->pos.y - n->parent->pos.y);
		JPS_ASSERT(dx || dy);
		JPS_ASSERT(!dx || !dy || abs(dx) == abs(dy));
		dx /= std::max(abs(dx), 1);
		dy /= std::max(abs(dy), 1);

		// TODO: only use cache if end node is not between n and end node

		Node **jumpPointsPtr = n->getJumpPointsForDirection(dx, dy, &njump);

		if(njump < 0)
		{
			// only cache jump points if end wasn't hit
			if(!identifyJumpPoints(n, jumpPoints, &njump))
			{
				JPS_ASSERT(njump <= MAX_NEIGHBORS);
				for(int i = 0; i < njump; ++i)
					jumpPointsPtr[i] = jumpPoints[i];
				n->setNumJumpPoints(dx, dy, njump);

				printf("[NEW]   JP for node %p (%u, %u), direction (%d, %d), parent %p: %d\n", n, n->pos.x, n->pos.y, dx, dy, n->parent, njump);
				for(int i = 0; i < njump; ++i)
					printf("           (%u, %u)\n", jumpPointsPtr[i]->pos.x, jumpPointsPtr[i]->pos.y);
			}
		}
		else
		{
			jumpPointsUse = jumpPointsPtr; // just use that pointer without copying contents

		#ifdef JPS_VERIFY
			int dummy, cc = 0;
			for(int y = -1; y <= 1; ++y)
				for(int x = -1; x <= 1; ++x)
					if(x || y)
						cc += jumpPointsPtr == n->getJumpPointsForDirection(x, y, &dummy);
			JPS_ASSERT(cc == 1);

			printf("[CHECK] JP for node %p (%u, %u), direction (%d, %d), parent %p: %d\n", n, n->pos.x, n->pos.y, dx, dy, n->parent, njump);
			JPS_ASSERT(njump <= MAX_NEIGHBORS);
			Node *jp[MAX_NEIGHBORS];
			int njj;
			bool hitEnd = identifyJumpPoints(n, &jp[0], &njj);
			(void)njj;
			JPS_ASSERT(!hitEnd); // should not be cached if end was hit
			JPS_ASSERT(njump == njj);
			for(int i = 0; i < njump; ++i)
			{
				printf("           (%u, %u)\n", jp[i]->pos.x, jp[i]->pos.y);
				JPS_ASSERT(jp[i] == jumpPointsPtr[i]); // array is uninitialized for i >= njump
			}
		#endif
		}
	}
#endif

	base.followSuccessors(n, njump, jumpPointsUse);
}

template <typename GRID> bool Searcher<GRID>::findPath(PathVector& path, Position start, Position end, unsigned step /* = 0 */)
{
	Result res = findPathInit(start, end);

	// If this is true, the resulting path is empty (findPathFinish() would fail, so this needs to be checked before)
	if(res == FOUND_PATH)
		return true;

	while(true)
	{
		switch(res)
		{
			case NEED_MORE_STEPS:
				res = findPathStep(0);
				break; // the switch

			case FOUND_PATH:
				return findPathFinish(path, step);

			case NO_PATH:
			default:
				return false;
		}
	}
}

template <typename GRID> Result Searcher<GRID>::findPathInit(Position start, Position end)
{
	for(NodeGrid::iterator it = base.nodegrid.begin(); it != base.nodegrid.end(); ++it)
		it->second.clearState();
	base.open.clear();
	base.endNode = NULL;
	base.stepsDone = 0;

	// If skip is > 1, make sure the points are aligned so that the search will always hit them
	const unsigned skip = base.skip;
	start.x = (start.x / skip) * skip;
	start.y = (start.y / skip) * skip;
	end.x = (end.x / skip) * skip;
	end.y = (end.y / skip) * skip;

	if(start == end)
	{
		// There is only a path if this single position is walkable.
		// But since the starting position is omitted, there is nothing to do here.
		return grid(end.x, end.y) ? FOUND_PATH : NO_PATH;
	}

	// If start or end point are obstructed, don't even start
	if(!grid(start.x, start.y) || !grid(end.x, end.y))
		return NO_PATH;

	base.open.push(base.getNode(start));
	base.endNode = base.getNode(end);
	JPS_ASSERT(base.endNode);

	return NEED_MORE_STEPS;
}

template <typename GRID> Result Searcher<GRID>::findPathStep(int limit)
{
	base.stepsRemain = limit;
	do
	{
		if(base.open.empty())
			return NO_PATH;
		Node *n = base.open.pop();
		n->setClosed();
		if(n == base.endNode)
			return FOUND_PATH;
		identifySuccessors(n);
	}
	while(base.stepsRemain >= 0);
	return NEED_MORE_STEPS;
}

template<typename GRID> bool Searcher<GRID>::findPathFinish(PathVector& path, unsigned step /* = 0 */)
{
	return base.generatePath(path, step);
}


} // end namespace Internal

using Internal::Searcher;

// Single-call convenience function
//
// path: If the function returns true, the path is stored in this vector.
//       The path does NOT contain the starting position, i.e. if start and end are the same,
//       the resulting path has no elements.
//       The vector does not have to be empty. The function does not clear it;
//       instead, the new path positions are appended at the end.
//       This allows building a path incrementally.
//
// grid: expected to overload operator()(x, y), return true if position is walkable, false if not.
//
// step: If 0, only return waypoints.
//       If 1, create exhaustive step-by-step path.
//       If N, put in one position for N blocks travelled, or when a waypoint is hit.
//       All returned points are guaranteed to be on a straight line (vertically, horizontally, or diagonally),
//       and there is no obstruction between any two consecutive points.
//       Note that this parameter does NOT influence the pathfinding in any way;
//       it only controls the coarseness of the output path.
//
// skip: If you know your map data well enough, this can be set to > 1 to speed up pathfinding even more.
//       Warning: Start and end positions will be rounded down to the nearest <skip>-aligned position,
//       so make sure to give appropriate positions so they do not end up in a wall.
//       This will also skip through walls if they are less than <skip> blocks thick at any reachable position.
template <typename GRID> bool findPath(PathVector& path, const GRID& grid, unsigned startx, unsigned starty, unsigned endx, unsigned endy,
                                       unsigned step = 0, int skip = 0, // optional params
                                       size_t *stepsDone = NULL, size_t *nodesExpanded = NULL // for information
                                       )
{
	Searcher<GRID> search(grid);
	search.setSkip(skip);
	bool found = search.findPath(path, Pos(startx, starty), Pos(endx, endy), step);
	if(stepsDone)
		*stepsDone = search.getStepsDone();
	if(nodesExpanded)
		*nodesExpanded = search.getNodesExpanded();
	return found;
}



} // end namespace JPS


#endif
