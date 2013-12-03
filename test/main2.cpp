#include "../JPS.h"
#include <iostream>
#include "ScenarioLoader.h"
#include <fstream>
#include <stdlib.h>


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
	unsigned lastx = path[0].first;
	unsigned lasty = path[0].second;
	float accu = 0;
	for(size_t i = 1; i < path.size(); ++i)
	{
		unsigned x = path[i].first;
		unsigned y = path[i].second;

		int dx = int(x - lastx);
		int dy = int(y - lasty);

		accu += sqrtf(float(dx*dx + dy*dy));
		
		lastx = x;
		lasty = y;
	}
	return accu;
}

int main(int argc, char **argv)
{
	ScenarioLoader loader("maps/AR0011SR.map.scen");
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

	return 0;
}
