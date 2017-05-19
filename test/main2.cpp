#include "../JPS.h"
#include <iostream>
#include "ScenarioLoader.h"
#include <fstream>
#include <stdlib.h>
#include <assert.h>

// Testing material from http://www.movingai.com/benchmarks/

static void die(const char *msg)
{
	std::cerr << msg << std::endl;
	abort();
}

struct MapGrid
{
	MapGrid(const char *file)
	{
		std::ifstream in(file);
		if(!in)
			die(file);

		std::string s;
		std::getline(in, s);
		in >> s >> h;
		in >> s >> w;
		in >> s;

		while(in >> s)
			if(s.length() == w)
				lines.push_back(s);

		if(h != lines.size())
			die("Wrong number of lines");

		std::cout << "W: " << w << "; H: " << h << "; Total cells: " << (w*h) << std::endl;
	}

	bool operator()(unsigned x, unsigned y) const
	{
		if(x < w && y < h)
		{
			const char c = lines[y][x];
			switch(c)
			{
				case '.':
				case 'G':
				case 'S':
					return true;
			}
		}
		return false;
	}

	unsigned w, h;
	std::vector<std::string> lines;
};

static double pathcost(unsigned startx, unsigned starty, const JPS::PathVector& path)
{
	unsigned lastx = startx;
	unsigned lasty = starty;
	double accu = 0;
	assert(path.empty() || path[0] != JPS::Pos(startx, starty));
	for(size_t i = 0; i < path.size(); ++i)
	{
		unsigned x = path[i].x;
		unsigned y = path[i].y;

		int dx = int(x - lastx);
		int dy = int(y - lasty);

		accu += sqrt(double(dx*dx + dy*dy));
		
		lastx = x;
		lasty = y;
	}
	return accu;
}

double runScenario(const char *file)
{
	ScenarioLoader loader(file);
	if(!loader.GetNumExperiments())
		die(file);
	MapGrid grid(loader.GetNthExperiment(0).GetMapName());
	double sum = 0;
	JPS::PathVector path;
	for(int i = 0; i < loader.GetNumExperiments(); ++i)
	{
		const Experiment& ex = loader.GetNthExperiment(i);
		path.clear();
		size_t stepsDone, nodesExpanded;
		int runs = 0;

		// single-call
		//bool found = JPS::findPath(path, grid, ex.GetStartX(), ex.GetStartY(), ex.GetGoalX(), ex.GetGoalY(), 0, 0, true, &stepsDone, &nodesExpanded);

		// Testing incremental runs
		bool found = false;
		JPS::Searcher<MapGrid> search(grid);
		JPS::Result res = search.findPathInit(JPS::Pos(ex.GetStartX(), ex.GetStartY()), JPS::Pos(ex.GetGoalX(), ex.GetGoalY()));
		if(res == JPS::EMPTY_PATH)
			found = true;
		else
		{
			while(res == JPS::NEED_MORE_STEPS)
			{
				++runs;
				res = search.findPathStep(10000);
			}
			found = (res == JPS::FOUND_PATH) && search.findPathFinish(path, 0);
		}
		stepsDone = search.getStepsDone();
		nodesExpanded = search.getNodesExpanded();


		if(!found)
		{
			printf("#### [%s:%d] PATH NOT FOUND: (%d, %d) -> (%d, %d)\n",
				file, i, ex.GetStartX(), ex.GetStartY(), ex.GetGoalX(), ex.GetGoalY());
			die("Path not found!"); // all paths known to be valid, so this is bad
			continue;
		}

		// Starting position is NOT included in vector
		double cost = pathcost(ex.GetStartX(), ex.GetStartY(), path);

		//if(cost > ex.GetDistance()+0.5f)
			printf("[%s] [%s:%d] Path len: %.3f (%.3f); Diff: %.3f; Steps: %u; Nodes: %u; Runs: %u\n",
				(cost > ex.GetDistance()+0.5f ? "##" : "  "), file, i, cost, ex.GetDistance(),
				fabs(cost - ex.GetDistance()), (unsigned)stepsDone, (unsigned)nodesExpanded, runs);

		sum += cost;
	}
	return sum;
}

int main(int argc, char **argv)
{
	double sum = 0;
	for(int i = 1; i < argc; ++i)
		sum += runScenario(argv[i]);

	std::cout << "Total distance travelled: " << sum << std::endl;

	return 0;
}

