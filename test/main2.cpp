#include "../JPS.h"
#include <iostream>
#include "ScenarioLoader.h"
#include <fstream>
#include <stdlib.h>

// Testing material from http://www.movingai.com/benchmarks/

const char *scenarios[] =
{
	"maps/AR0011SR.map.scen",
	"maps/den011d.map.scen",
	"maps/den602d.map.scen",
	"maps/hrt201n.map.scen",
	NULL
};

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

		std::cout << "W: " << w << "; H: " << h << "; Total cells: " << (w*h) << std::endl;
		std::cout << "Lines: " << lines.size() << std::endl;
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

static float pathcost(const JPS::PathVector& path)
{
	unsigned lastx = path[0].x;
	unsigned lasty = path[0].y;
	float accu = 0;
	for(size_t i = 1; i < path.size(); ++i)
	{
		unsigned x = path[i].x;
		unsigned y = path[i].y;

		int dx = int(x - lastx);
		int dy = int(y - lasty);

		accu += sqrtf(float(dx*dx + dy*dy));
		
		lastx = x;
		lasty = y;
	}
	return accu;
}

void runScenario(const char *file)
{
	ScenarioLoader loader(file);
	if(!loader.GetNumExperiments())
		die(file);
	MapGrid grid(loader.GetNthExperiment(0).GetMapName());
	for(int i = 0; i < loader.GetNumExperiments(); ++i)
	{
		Experiment ex = loader.GetNthExperiment(i);
		JPS::PathVector path;
		bool found = JPS::findPath(path, grid, ex.GetStartX(), ex.GetStartY(), ex.GetGoalX(), ex.GetGoalY(), true);
		if(!found)
			die("Path not found");

		float cost = pathcost(path);

		printf("Path len: %.3f; Expected: %.3f\n", cost, ex.GetDistance());
		if(cost > ex.GetDistance()+0.5f)
			printf(" -- PATH TOO LONG\n");
	}
}

int main(int argc, char **argv)
{
	for(unsigned i = 0; scenarios[i]; ++i)
		runScenario(scenarios[i]);

	return 0;
}

