#include "../JPS.h"
#include <iostream>
#include <string>
#include <string.h>

/*
static const char *data[] =
{
	"############",
	"#   #      #",
	"# s #  ### #",
	"#   #  # e #",
	"#      #   #",
	"############",
	NULL
};
*/

static const char *data[] =
{
	"##############################################",
	"#                           #                #",
	"#                   #       #        #####   #",
	"#                   #       #       ##       #",
	"#                   #   e   #      ##        #",
	"#                   #       #     ##         #",
	"#                   #      ##    ##          #",
	"#      ##############           ##           #",
	"#                   #           #            #",
	"#                   #############            #",
	"#               #####                        #",
	"#                                            #",
	"#                 #                     s    #",
	"#                 #                          #",
	"##############################################",
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

	unsigned startx, starty, endx, endy;
	for(unsigned y = 0; y < grid.h; ++y)
	{
		const char *s = data[y];
		const char *sp = strchr(s, 's');
		if(sp)
		{
			startx = sp - s;
			starty = y;
		}
		const char *ep = strchr(s, 'e');
		if(ep)
		{
			endx = ep - s;
			endy = y;
		}
	}

#define PUT(x, y, v) (grid.out[(y)][(x)] = (v))

	JPS::PathVector path;
	bool found = JPS::findPath(path, grid, startx, starty, endx, endy, true);

	if(found)
	{
		unsigned c = 0;
		for(JPS::PathVector::iterator it = path.begin(); it != path.end(); ++it)
			PUT(it->first, it->second, (c++ % 26) + 'a');
	}

	for(unsigned i = 0; i < grid.h; ++i)
		std::cout << grid.out[i] << std::endl;

	if(!found)
		std::cout << "No path found!" << std::endl;

	return found ? 0 : 1;
}
