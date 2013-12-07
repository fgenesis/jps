#include "../JPS.h"
#include <iostream>
#include <string>
#include <string.h>
#include <assert.h>

static const char *data[] =
{
	"##############################################",
	"#                           #                #",
	"#                   #       #         ####   #",
	"#                   #       #        #       #",
	"#                   #   2   #       #        #",
	"#                   #       #      #5        #",
	"#                   #      #     ######      #",
	"#      #############            #            #",
	"#                  3#           #            #",
	"#                   ############             #",
	"#               ####                         #",
	"#                                        #####",
	"#                 #                        1 #",
	"#                 #         4                #",
	"##############################################",
	NULL
};

struct MyGrid
{
	~MyGrid()
	{
		delete [] out;
	}

	MyGrid(const char *d[])
		: mapdata(d)
	{
		w = -1;
		h = 0;
		for(; mapdata[h]; ++h)
			w = std::min(w, strlen(mapdata[h]));

		out = new std::string[h];
		for(unsigned i = 0; i < h; ++i)
			out[i] = mapdata[i];

		std::cout << "W: " << w << "; H: " << h << "; Total cells: " << (w*h) << std::endl;
	}

	bool operator()(unsigned x, unsigned y) const
	{
		if(x < w && y < h)
		{
			if(mapdata[y][x] == '#')
			{
				out[y][x] = '@';
			}
			else
			{
				out[y][x] = '.';
				return true;
			}
		}
		return false;
	}

	unsigned w, h;
	const char **mapdata;
	std::string *out;
};


int main(int argc, char **argv)
{
	MyGrid grid(data);

	JPS::PathVector waypoints;
	for(char a = '1'; a <= '9'; ++a)
	{
		for(unsigned y = 0; y < grid.h; ++y)
		{
			const char *sp = strchr(data[y], a);
			if(sp)
			{
				waypoints.push_back(JPS::Pos(sp - data[y], y));
			}
		}
	}

	unsigned step = argc > 1 ? atoi(argv[1]) : 0;
	std::cout << "Calculating path with step " << step << std::endl;

	JPS::PathVector path;
	for(size_t i = 1; i < waypoints.size(); ++i)
	{
		bool found = JPS::findPath(path, grid, waypoints[i-1].x, waypoints[i-1].y, waypoints[i].x, waypoints[i].y, step);
		if(found)
		{
			assert(path[0] != waypoints[i-1]);
		}
		else
		{
			std::cout << "Path not found!" << std::endl;
			break;
		}
	}


#define PUT(x, y, v) (grid.out[(y)][(x)] = (v))

	unsigned c = 0;
	for(JPS::PathVector::iterator it = path.begin(); it != path.end(); ++it)
		PUT(it->x, it->y, (c++ % 26) + 'a');

	for(unsigned i = 0; i < grid.h; ++i)
		std::cout << grid.out[i] << std::endl;

	return 0;
}
