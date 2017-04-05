#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "recordlocker.cc"


void dostuff()
{
	std::string dir = "lockertestdir";
	mkdir(dir.c_str(), 0777);

	recordlocker locker(dir, 50);
	for(int gen = 1; gen < 50; ++gen)
	{
		while(auto r = locker.lock_record(gen))
		{
			std::ostringstream fname;
			fname << dir << "/" << r << ".out";

//			usleep(rand() % 100000);

			std::ofstream out(fname.str(), std::ios::app);
			char host[1024];
			gethostname(host, 1024);
			out << gen << "\t" << host << "\t" << getpid() << "\t" << locker.get_value() << "\n";
			out.close();

			locker.unlock_record(gen);

//			usleep(rand() % 10000);
		}
	}
}


int main()
{
	try
	{
		dostuff();
	}
	catch(std::string s)
	{
		std::cerr << s << "\n";
		return 1;
	}
}
