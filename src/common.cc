#include <cstring>
#include <libgen.h>

#include "common.h"


// Placed here because there's no common.cc at the moment
const std::string &sourcedir(const char *arg)
{
	static std::string dir;
	if(arg)
	{
		std::vector<char> v(arg, arg + std::strlen(arg) + 1);
		dir = dirname(&v[0]);
	}
	return dir;
}

